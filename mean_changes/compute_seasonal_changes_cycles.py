import numpy as np
from netCDF4 import Dataset
import datetime
import scipy.stats
import os

import xarray as xr
import gc
# -------- CONFIGURATION --------
variable = 't'  # 't', 'pr', or 'wind'
period = 'rcp85'
domain = 'd03'

basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
BASE_DATES = {
    "hist"  : datetime.datetime(1986, 1, 1),
    "rcp85" : datetime.datetime(2046, 1, 1),
    "rcp45" : datetime.datetime(2046, 1, 1),
}
if variable == "t":
    var = 'T2'
    filename = 't_' + domain + '_daily'
elif variable == "pr":
    var = 'pr'
    filename = 'pr_' + domain + '_daily'
elif variable == "wind":
    var = 'wspd'
    filename = "wind_" + domain + '_daily_wspd'
else:
    raise ValueError("Unsupported variable.")

# -------- HELPER FUNCTIONS --------
def load_data(path: str, varname: str, scenario: str):
    """
    Load a variable from a NetCDF file using netCDF4.Dataset.

    Returns
    -------
    data      : np.ndarray  (time, lat, lon)  – MaskedArray → plain float32
    datetimes : list of datetime.datetime
    lat       : np.ndarray  1-D
    lon       : np.ndarray  1-D
    """
    print(f"    Opening {path} …")
    nc        = Dataset(path, "r")
    raw       = nc.variables[varname][:]          # MaskedArray
    times     = nc.variables["time"][:]
    lat       = nc.variables["lat"][:]
    lon       = nc.variables["lon"][:]
    nc.close()                                    # close immediately

    # ── convert MaskedArray → plain float32 ndarray ───────────────────────
    if isinstance(raw, np.ma.MaskedArray):
        data = np.ma.filled(raw, fill_value=np.nan).astype(np.float32)
    else:
        data = raw.astype(np.float32)

    # ── build datetime list from hours-since base date ────────────────────
    base_date = BASE_DATES[scenario]
    datetimes = [base_date + datetime.timedelta(hours=float(h)) for h in times]

    # ── also strip masks from lat / lon if present ────────────────────────
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
        # regular grid
        da = xr.DataArray(
            data,
            dims   = ["time", "lat", "lon"],
            coords = {"time": time_index, "lat": lat, "lon": lon},
            name   = name,
        )
    else:
        # curvilinear / 2-D lat-lon grid
        da = xr.DataArray(
            data,
            dims   = ["time", "y", "x"],
            coords = {
                "time" : time_index,
                "lat"  : (["y", "x"], lat),   # ← 2-D coord needs explicit dims
                "lon"  : (["y", "x"], lon),
            },
            name = name,
        )
    return da


# ── Helper: get_seasonal ──────────────────────────────────────────────────────
def get_seasonal(da):
    """
    Compute the mean seasonal cycle and derived diagnostics for a DataArray.

    Parameters
    ----------
    da : xr.DataArray  (time, lat, lon)

    Returns
    -------
    xr.Dataset with amplitude, day_of_max, day_of_min
    """
    # 1. Drop Feb 29
    da = da.sel(time=~((da.time.dt.month == 2) & (da.time.dt.day == 29)))

    # 2. Mean seasonal cycle (365 groups)
    seasonal_cycle = da.groupby("time.dayofyear").mean(dim="time")

    # 3. Amplitude = peak-to-peak of seasonal cycle
    amplitude = (
        seasonal_cycle.max(dim="dayofyear")
        - seasonal_cycle.min(dim="dayofyear")
    )

    # 4. Day of max / min (1-indexed)
    day_of_max = seasonal_cycle.argmax(dim="dayofyear") + 1
    day_of_min = seasonal_cycle.argmin(dim="dayofyear") + 1

    return xr.Dataset({
        "amplitude"  : amplitude,
        "day_of_max" : day_of_max,
        "day_of_min" : day_of_min,
    })


# ── Helper: circular day difference (−182 … +182) ────────────────────────────
def circular_day_diff(future, hist, period=365):
    """
    Compute the shortest signed difference between two day-of-year fields.
    Result is in (−period/2, +period/2].

    Parameters
    ----------
    future, hist : xr.DataArray  – day-of-year values
    period       : int           – length of the cycle (365)

    Returns
    -------
    xr.DataArray with values in (−182, +182]
    """
    diff = future - hist
    diff = ((diff + period / 2) % period) - period / 2
    return diff

# ── Configuration ─────────────────────────────────────────────────────────────
# -------- LOAD DATA --------
hist_path = os.path.join(basepath, f"historical/variables_complete/{filename}.nc")
fut_path = os.path.join(basepath, f"{period}/variables_complete/{filename}.nc")

data_hist, time_hist,lat,lon = load_data(hist_path, var, 'hist')
data_fut, time_fut,lat,lon = load_data(fut_path, var, period)

print(data_hist)
# ── 2. Wrap into xr.DataArray ─────────────────────────────────────
print("  Building DataArrays …")
da_hist = to_dataarray(data_hist, time_hist, lat, lon, name=var)
da_scen = to_dataarray(data_fut, time_fut, lat, lon, name=var)

# free raw numpy arrays immediately
del data_hist, time_hist, data_fut, time_fut
gc.collect()


VARIABLES = {
    variable   : {"hist": da_hist,   period: da_scen},   # replace with your DataArrays
}

OUTPUT_DIR = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/MEANS/'


# ── Main loop ─────────────────────────────────────────────────────────────────
for var_name, das in VARIABLES.items():
    print(f"\n{'='*60}")
    print(f"  Processing : {var_name.upper()}  |  scenario : {period}")
    print(f"{'='*60}")

    # ── a) Compute seasonal diagnostics ──────────────────────────────────
    print(f"  [{var_name}] Computing historical seasonal cycle …")
    hist_ds = get_seasonal(das["hist"])

    print(f"  [{var_name}] Computing {period} seasonal cycle …")
    scen_ds = get_seasonal(das[period])

    # ── b) Relative change in amplitude (%) ──────────────────────────────
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

    # ── c) Change in day of maximum (circular) ────────────────────────────
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

    # ── d) Change in day of minimum (circular) ────────────────────────────
    delta_day_of_min = circular_day_diff(
        scen_ds["day_of_min"], hist_ds["day_of_min"]
    )
    delta_day_of_min.name = "delta_day_of_min"
    delta_day_of_min.attrs.update({
        "long_name" : f"Change in day of minimum ({var_name})",
        "units"     : "days",
        "note"      : "Positive = later in the year; range (−182, +182]",
    })

    # ── e) Also save the raw hist / scenario amplitudes and days ─────────
    #     Useful for downstream analysis without re-running get_seasonal
    hist_ds["amplitude"].attrs.update({
        "long_name" : f"Historical seasonal amplitude ({var_name})",
        "units"     : das["hist"].attrs.get("units", "unknown"),
    })
    scen_ds["amplitude"].attrs.update({
        "long_name" : f"{period} seasonal amplitude ({var_name})",
        "units"     : das["hist"].attrs.get("units", "unknown"),
    })

    # ── f) Build per-variable Dataset ─────────────────────────────────────
    ds_var = xr.Dataset({
        "rel_amp_change"       : rel_amp_change,
        "delta_day_of_max"     : delta_day_of_max,
        "delta_day_of_min"     : delta_day_of_min,
        # raw diagnostics for reference
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
        "method"      : "Chen et al. (2019) peak-to-peak amplitude; circular phase shift",
        "created_by"  : "seasonal_change_loop.py",
    }

    # ── g) Print quick summary ────────────────────────────────────────────
    print(f"  [{var_name}] rel_amp_change  : "
          f"mean = {float(rel_amp_change.mean(skipna=True)):+.2f} %")
    print(f"  [{var_name}] delta_day_of_max: "
          f"mean = {float(delta_day_of_max.mean(skipna=True)):+.2f} days")
    print(f"  [{var_name}] delta_day_of_min: "
          f"mean = {float(delta_day_of_min.mean(skipna=True)):+.2f} days")

    # ── h) Save to NetCDF ─────────────────────────────────────────────────
    out_path = f"{OUTPUT_DIR}{var_name}_seasonalAmp_change_{period}.nc"
    print(f"  [{var_name}] Writing → {out_path} …")
    ds_var.to_netcdf(out_path)
    print(f"  [{var_name}] Saved.")

    # ── i) Free memory ────────────────────────────────────────────────────
    del hist_ds, scen_ds, hist_amp, scen_amp
    del rel_amp_change, delta_day_of_max, delta_day_of_min, ds_var
    gc.collect()

print(f"\nAll variables processed for scenario: {period}")



