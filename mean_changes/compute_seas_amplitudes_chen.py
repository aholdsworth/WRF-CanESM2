import numpy as np
from netCDF4 import Dataset
import datetime
import scipy.stats
import os
import xarray as xr
import gc
import pandas as pd
# -------- CONFIGURATION --------
variable = 't'       # 't', 'pr', or 'wind'
period   = 'rcp85'
domain   = 'd03'

basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'

BASE_DATES = {
    "hist"  : datetime.datetime(1986, 1, 1),
    "rcp85" : datetime.datetime(2046, 1, 1),
    "rcp45" : datetime.datetime(2046, 1, 1),
}


if variable == "t":
    var      = 'T2'
    filename = 't_' + domain + '_daily'
elif variable == "pr":
    var      = 'pr'
    filename = 'pr_' + domain + '_daily'
elif variable == "wind":
    var      = 'wspd'
    filename = "wind_" + domain + '_daily_wspd'
else:
    raise ValueError("Unsupported variable.")

OUTPUT_DIR = '/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/MEANS/'

# ======================================================================
# HELPER FUNCTIONS
# ======================================================================

def load_data(path: str, varname: str, scenario: str):
    """
    Load a variable from a NetCDF file.

    Returns
    -------
    data      : np.ndarray  (time, y, x)  float32
    datetimes : list of datetime.datetime
    lat       : np.ndarray  (y, x) or (y,)  float32
    lon       : np.ndarray  (y, x) or (x,)  float32
    """
    print(f"    Opening {path} …")
    nc    = Dataset(path, "r")
    raw   = nc.variables[varname][:]
    times = nc.variables["time"][:]
    lat   = nc.variables["lat"][:]
    lon   = nc.variables["lon"][:]
    nc.close()

    # MaskedArray → plain float32
    if isinstance(raw, np.ma.MaskedArray):
        data = np.ma.filled(raw, fill_value=np.nan).astype(np.float32)
    else:
        data = raw.astype(np.float32)

    # Build datetime list from hours-since base date
    base_date = BASE_DATES[scenario]
    datetimes = [
        base_date + datetime.timedelta(hours=float(h)) for h in times
    ]

    # Strip masks from lat / lon
    if isinstance(lat, np.ma.MaskedArray):
        lat = np.ma.filled(lat, np.nan)
    if isinstance(lon, np.ma.MaskedArray):
        lon = np.ma.filled(lon, np.nan)

    lat = lat.astype(np.float32)
    lon = lon.astype(np.float32)

    print(f"    shape={data.shape}  dtype={data.dtype}  "
          f"NaNs={int(np.isnan(data).sum())}")
    return data, datetimes, lat, lon



def to_dataarray(data, datetimes, lat, lon, name):
    time_index = np.array(datetimes, dtype="datetime64[ns]")

    if lat.ndim == 1:
        da = xr.DataArray(
            data,
            dims   = ["time", "lat", "lon"],
            coords = {"time": time_index, "lat": lat, "lon": lon},
            name   = name,
        )
    else:
        da = xr.DataArray(
            data,
            dims   = ["time", "y", "x"],
            coords = {
                "time" : time_index,
                "lat"  : (["y", "x"], lat),
                "lon"  : (["y", "x"], lon),
            },
            name = name,
        )
    return da


def get_seasonal(da):
    """
    Compute mean seasonal cycle and derived diagnostics.

    Returns
    -------
    xr.Dataset with amplitude, day_of_max, day_of_min
    """
    # Drop Feb 29
    da = da.sel(time=~((da.time.dt.month == 2) & (da.time.dt.day == 29)))

    # Mean seasonal cycle  (365 day-of-year groups)
    seasonal_cycle = da.groupby("time.dayofyear").mean(dim="time")

    # Peak-to-peak amplitude
    amplitude = (
        seasonal_cycle.max(dim="dayofyear")
        - seasonal_cycle.min(dim="dayofyear")
    )

    # Day of max / min (1-indexed)
    day_of_max = seasonal_cycle.argmax(dim="dayofyear") + 1
    day_of_min = seasonal_cycle.argmin(dim="dayofyear") + 1

    return xr.Dataset({
        "amplitude"  : amplitude,
        "day_of_max" : day_of_max,
        "day_of_min" : day_of_min,
    })


def circular_day_diff(future, hist, cycle=365):
    """
    Shortest signed circular difference in day-of-year.
    Result in (−cycle/2, +cycle/2].
    """
    diff = future - hist
    diff = ((diff + cycle / 2) % cycle) - cycle / 2
    return diff


def weighted_masked_mean(data_2d, weights_2d, mask_2d):
    """
    Area-weighted mean of data_2d where mask_2d is True.

    Parameters
    ----------
    data_2d    : np.ndarray (y, x)
    weights_2d : np.ndarray (y, x)  – e.g. cos(lat)
    mask_2d    : np.ndarray (y, x)  – boolean

    Returns
    -------
    scalar float
    """
    masked_data    = np.where(mask_2d, data_2d,    np.nan)
    masked_weights = np.where(mask_2d, weights_2d, 0.0)

    numerator   = np.nansum(masked_data * masked_weights)
    denominator = np.nansum(masked_weights)

    return float(numerator / denominator) if denominator != 0 else np.nan
# ======================================================================
# LOAD DATA
# ======================================================================

hist_path = os.path.join(basepath, f"historical/variables_complete/{filename}.nc")
fut_path  = os.path.join(basepath, f"{period}/variables_complete/{filename}.nc")

data_hist, time_hist, lat, lon = load_data(hist_path, var, 'hist')
data_fut,  time_fut,  lat, lon = load_data(fut_path,  var,  period)

# ── Load land mask from geo_em file ──────────────────────────────────────────
geo_em_file = os.path.join(basepath, f"domain/geo_em.{domain}.nc")   # update path if needed

landmask = xr.open_dataset(geo_em_file)['LANDMASK'].squeeze()

# Rename WRF dimension names to match our data
if 'south_north' in landmask.dims:
    landmask = landmask.rename({'south_north': 'y'})
if 'west_east' in landmask.dims:
    landmask = landmask.rename({'west_east': 'x'})

land_mask_np  = (landmask.values == 1)   # True where land
ocean_mask_np = (landmask.values == 0)   # True where ocean

# Sanity check
assert land_mask_np.shape == data_hist.shape[1:], (
    f"Land mask shape {land_mask_np.shape} does not match "
    f"data spatial shape {data_hist.shape[1:]}"
)

# Pre-compute cosine-latitude weights (2-D)
if lat.ndim == 1:
    lat_2d = np.broadcast_to(lat[:, np.newaxis], land_mask_np.shape).astype(np.float32)
else:
    lat_2d = lat.astype(np.float32)

weights_2d = np.cos(np.deg2rad(lat_2d)).astype(np.float32)
del lat_2d
gc.collect()


# ── Wrap into xr.DataArray ────────────────────────────────────────────────────
print("  Building DataArrays …")
da_hist = to_dataarray(data_hist, time_hist, lat, lon, name=var)
da_scen = to_dataarray(data_fut,  time_fut,  lat, lon, name=var)

del data_hist, time_hist, data_fut, time_fut
gc.collect()

VARIABLES = {
    variable: {"hist": da_hist, period: da_scen},
}

# ======================================================================
# MAIN LOOP
# ======================================================================

for var_name, das in VARIABLES.items():
    print(f"\n{'='*60}")
    print(f"  Processing : {var_name.upper()}  |  scenario : {period}")
    print(f"{'='*60}")

    # ── a) Seasonal diagnostics ───────────────────────────────────────────
    print(f"  [{var_name}] Computing historical seasonal cycle …")
    hist_ds = get_seasonal(das["hist"])

    print(f"  [{var_name}] Computing {period} seasonal cycle …")
    scen_ds = get_seasonal(das[period])

    # ── b) Relative change in amplitude ──────────────────────────────────
    print(f"  [{var_name}] Computing relative change in amplitude …")
    hist_amp = hist_ds["amplitude"]
    scen_amp = scen_ds["amplitude"]

    with np.errstate(invalid="ignore", divide="ignore"):
        rel_amp_change = xr.where(
            hist_amp != 0,
            100.0 * (scen_amp - hist_amp) / hist_amp,
            np.nan,
        )
    rel_amp_change.name = "rel_amp_change"
    rel_amp_change.attrs.update({
        "long_name" : f"Relative change in seasonal amplitude ({var_name})",
        "units"     : "%",
        "note"      : f"100 * ({period} - hist) / hist",
    })

    # ── c) Change in day of maximum ───────────────────────────────────────
    print(f"  [{var_name}] Computing change in day of max / min …")
    delta_day_of_max = circular_day_diff(
        scen_ds["day_of_max"], hist_ds["day_of_max"]
    )
    delta_day_of_max.name = "delta_day_of_max"
    delta_day_of_max.attrs.update({
        "long_name" : f"Change in day of maximum ({var_name})",
        "units"     : "days",
        "note"      : "Positive = later in the year; range (−182, +182]",
    })

    # ── d) Change in day of minimum ───────────────────────────────────────
    delta_day_of_min = circular_day_diff(
        scen_ds["day_of_min"], hist_ds["day_of_min"]
    )
    delta_day_of_min.name = "delta_day_of_min"
    delta_day_of_min.attrs.update({
        "long_name" : f"Change in day of minimum ({var_name})",
        "units"     : "days",
        "note"      : "Positive = later in the year; range (−182, +182]",
    })

    # ── e) Attribute updates ──────────────────────────────────────────────
    hist_ds["amplitude"].attrs.update({
        "long_name" : f"Historical seasonal amplitude ({var_name})",
        "units"     : das["hist"].attrs.get("units", "unknown"),
    })
    scen_ds["amplitude"].attrs.update({
        "long_name" : f"{period} seasonal amplitude ({var_name})",
        "units"     : das["hist"].attrs.get("units", "unknown"),
    })

    # ── f) Build output Dataset ───────────────────────────────────────────
    ds_var = xr.Dataset({
        "rel_amp_change"       : rel_amp_change,
        "delta_day_of_max"     : delta_day_of_max,
        "delta_day_of_min"     : delta_day_of_min,
        "hist_amplitude"       : hist_ds["amplitude"],
        "hist_day_of_max"      : hist_ds["day_of_max"],
        "hist_day_of_min"      : hist_ds["day_of_min"],
        f"{period}_amplitude"  : scen_ds["amplitude"],
        f"{period}_day_of_max" : scen_ds["day_of_max"],
        f"{period}_day_of_min" : scen_ds["day_of_min"],
    })

    ds_var.attrs = {
        "description" : (
            f"Seasonal cycle diagnostics and {period} vs historical changes "
            f"for variable: {var_name}"
        ),
        "variable"    : var_name,
        "scenario"    : period,
        "method"      : (
            "Chen et al. (2019) peak-to-peak amplitude; "
            "circular phase shift"
        ),
        "created_by"  : "seasonal_change_loop.py",
    }

    # ── g) Global summary (all grid points) ──────────────────────────────
    print(f"  [{var_name}] rel_amp_change  : "
          f"mean = {float(rel_amp_change.mean(skipna=True)):+.2f} %")
    print(f"  [{var_name}] delta_day_of_max: "
          f"mean = {float(delta_day_of_max.mean(skipna=True)):+.2f} days")
    print(f"  [{var_name}] delta_day_of_min: "
          f"mean = {float(delta_day_of_min.mean(skipna=True)):+.2f} days")

    # ── h) Land / ocean weighted averages ────────────────────────────────
    print(f"  [{var_name}] Computing land / ocean weighted averages …")
# ── h) Land / ocean weighted averages ────────────────────────────────
    print(f"  [{var_name}] Computing land / ocean weighted averages …")

    # Grab units string for column names
    data_units = das["hist"].attrs.get("units", "unknown")

    # Extract 2-D numpy arrays
    pct_np      = rel_amp_change.values.astype(np.float32)
    dmax_np     = delta_day_of_max.values.astype(np.float32)
    dmin_np     = delta_day_of_min.values.astype(np.float32)
    hist_amp_np = hist_ds["amplitude"].values.astype(np.float32)
    scen_amp_np = scen_ds["amplitude"].values.astype(np.float32)

    results = {}
    for region, mask in [("Land", land_mask_np), ("Ocean", ocean_mask_np)]:
        results[region] = {
            "Relative Amplitude Change (%)"       : weighted_masked_mean(pct_np,      weights_2d, mask),
            "Delta Day of Maximum (days)"          : weighted_masked_mean(dmax_np,     weights_2d, mask),
            "Delta Day of Minimum (days)"          : weighted_masked_mean(dmin_np,     weights_2d, mask),
            f"Historical Amplitude ({data_units})" : weighted_masked_mean(hist_amp_np, weights_2d, mask),
            f"{period} Amplitude ({data_units})"   : weighted_masked_mean(scen_amp_np, weights_2d, mask),
        }

    # Build DataFrame and save
    df = pd.DataFrame(results).T
    df.index.name = "Region"
    df = df.round(4)

    print(f"\n{df.to_string()}\n")

    csv_path = f"{OUTPUT_DIR}{var_name}_seasonalAmp_landocean_avg_{period}.csv"
    df.to_csv(csv_path)
    print(f"  [{var_name}] Regional summary (CSV) → {csv_path}")

    # Free temporary arrays
    del pct_np, dmax_np, dmin_np, hist_amp_np, scen_amp_np, df
    gc.collect()

print(f"\nAll variables processed for scenario: {period}")
