#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Count drizzle days per season for historical and future daily precipitation.
Uses dask for chunked/lazy loading.
"""

import xarray as xr

# --------------------------
# User settings
# --------------------------
period='rcp45'
basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
out_path = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/percentiles/'
hist_file = basepath+'/historical/variables_complete/pr_d03_daily.nc'
fut_file  = basepath+period+'/variables_complete/pr_d03_daily.nc'
varname   = 'pr'       # precip variable name in file
drizzle_threshold = 0.2  # mm/day threshold for drizzle
out_hist = out_path+'drizzle_days_hist'+'_'+str(drizzle_threshold)+'.nc'
out_fut  = out_path+'drizzle_days_'+period+'_'+str(drizzle_threshold)+'.nc'

# Set chunk sizes to fit your memory / cluster
# (time is most important, y and x optional)
chunks = {'time': 365, 'y': 200, 'x': 200}

# --------------------------
# Function
# --------------------------
def drizzle_day_count(da, drizzle_threshold=1.0):
    """
    Count number of drizzle days (< drizzle_threshold) per season.
    Returns a DataArray with dims ('season','y','x').
    """
    # Flag drizzle days lazily (boolean mask)
    drizzle = da < drizzle_threshold

    # Add season coordinate
    drizzle = drizzle.assign_coords(season=drizzle['time'].dt.season)

    # Sum over time dimension per season (lazy with dask)
    counts = drizzle.groupby('season').sum(dim='time')

    return counts

# --------------------------
# Historical
# --------------------------
print("Loading historical precip …")
ds_hist = xr.open_dataset(hist_file, chunks=chunks)
pr_hist = ds_hist[varname]

print("Computing drizzle-day counts for historical …")
drizzle_hist = drizzle_day_count(pr_hist, drizzle_threshold)

# Name and save to NetCDF (this triggers the dask computation)
drizzle_hist.name = f"drizzle_days_{drizzle_threshold:g}mm"
drizzle_hist.to_netcdf(out_hist, compute=True)
print(f"Saved {out_hist}")

# --------------------------
# Future
# --------------------------
print("Loading future precip …")
ds_fut = xr.open_dataset(fut_file, chunks=chunks)
pr_fut = ds_fut[varname]

print("Computing drizzle-day counts for future …")
drizzle_fut = drizzle_day_count(pr_fut, drizzle_threshold)

drizzle_fut.name = f"drizzle_days_{drizzle_threshold:g}mm"
drizzle_fut.to_netcdf(out_fut, compute=True)
print(f"Saved {out_fut}")

print("All done.")

