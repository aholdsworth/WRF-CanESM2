import xarray as xr
import numpy as np
import os
import datetime
from dask import delayed, compute

# --------------------
# Settings
# --------------------
basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
out_path = '/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/percentiles/'

variables = {
    "tas": ["T2", [10, 90]],   # 2m temperature, 10th & 90th
    "pr":  ["pr", [25, 75]],       # precipitation, 75th
}
#periods = ["hist", "rcp45", "rcp85"]
periods = ["rcp85"]
domain = 'd03'
file_map = {'tas': f't_{domain}_daily', 'pr': f'pr_{domain}_daily'}

# --------------------
# Helper for time decoding
# --------------------
def decode_time(ds, period):
    if "time" in ds:
        try:
            ds = xr.decode_cf(ds)
        except Exception:
            base_date = datetime.datetime(1986, 1, 1) if period == "hist" else datetime.datetime(2046, 1, 1)
            ds["time"] = [base_date + datetime.timedelta(hours=int(h)) for h in ds["time"].values]
    return ds

def rolling_doy_percentile(da, p, window):
    """
    Compute rolling-window DOY percentile of DataArray da.
    da: xarray DataArray with time dimension.
    p: percentile (0-100).
    window: odd number of days in window.
    """
    doy = da['time'].dt.dayofyear
    daysinyear = 366
    half = (window-1)//2

    # Prepare empty list
    percs = []
    for i in range(1, daysinyear+1):
        # wrap-around selection 
        if i <= half:
            mask = (doy >= daysinyear-half+(i-1)) | (doy <= i+half)
        elif i >= daysinyear-half:
            j = daysinyear - i + 1  # similar to old code
            mask = (doy >= i-half) | (doy <= (3-j))
        else:
            mask = (doy >= i-half) & (doy <= i+half)

        subset = da.sel(time=mask)
        percs.append(np.nanpercentile(subset, p, axis=0))

    # Stack back into DataArray (dayofyear, y, x)
    perc_da = xr.DataArray(
        np.stack(percs, axis=0),
        coords={'dayofyear': np.arange(1, daysinyear+1),
                'lon': da['lon'], 'lat': da['lat']},
        dims=('dayofyear','lon','lat')
    )
    return perc_da

# --------------------
# Compute and save function
# --------------------
@delayed
def compute_and_save(infile, var, varname, period, p, window):
    ds = xr.open_dataset(infile, chunks={"time": 365})
    ds = decode_time(ds, period)
    da = ds[varname]

    per = rolling_doy_percentile(da, p, window)

    out_file = f"{out_path}{var}{p}p_{period}_roll{window}.nc"
    per.name = f"{var}_{p}p"
    per.to_netcdf(out_file)

    return f"Saved {out_file}"

# --------------------
# Main execution
# --------------------
def main():
    tasks = []
    for var, (varname, percs) in variables.items():
        filename = file_map[var]
        for period in periods:
            if period == "hist":
                infile = os.path.join(basepath, f"historical/variables_complete/{filename}.nc")
            else:
                infile = os.path.join(basepath, f"{period}/variables_complete/{filename}.nc")

            if not os.path.exists(infile):
                print(f"!! Skipping {infile} (not found)")
                continue

            print(f"Scheduling {infile} ...")
            for p in percs:
                # choose window size depending on variable if you like
                window = 5 if var=="tas" else 29
                tasks.append(compute_and_save(infile, var, varname, period, p, window))

    # Run tasks in parallel
    results = compute(*tasks, scheduler="processes")
    for r in results:
        print(r)

    print("All percentiles computed.")

if __name__ == "__main__":
    main()

