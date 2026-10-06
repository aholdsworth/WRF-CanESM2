#!/usr/bin/env python
"""
Extract WRF d03 wind (wspd) for station 951 (HOPE SLIDE) using bilinear
interpolation on the native d03 300x300 grid (wspd_d03_hourly.nc).

IMPORTANT CONTEXT (validated 2026-09-17):
  * The 17 existing per-station files data/wrf_stations/wind_d03_st<ID>.nc
    were remapped (cdo genbic + remap, Dec 2023) using Eva's station list
    "st_<ID>.txt".  For 15 of the 17, those coordinates do NOT match
    files/ECCC_d03_stations.csv (only 1001 KELOWNA and 381 BELLA COOLA match).
    The CDO files therefore carry the *Eva-list* coordinates embedded.
  * Validation strategy: extract from the raw grid at the CDO-embedded
    coordinates for all 17 stations; hourly correlation and 1986-2005 annual
    means must agree closely (bilinear vs CDO genbic bilinear).
  * For 951, use the true coordinates (49.27, -121.23) from
    files/ECCC_d03_stations.csv, which agree with the PFM station file
    ECCC_951.txt used for the Dec 2024 t_d01/pr_d01 remaps.

Run: /home/amh001/space_fs7/software_2022/python/py_2024/bin/python extract_wrf_wind_951.py
"""
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

BASE = '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS'
PFM  = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF'
WRF_ST_DIR = os.path.join(BASE, 'data', 'wrf_stations')
GRID_FILE  = os.path.join(PFM, 'historical', 'variables_complete', 'wspd_d03_hourly.nc')

START_YEAR, END_YEAR = 1986, 2005
YEARS = list(range(START_YEAR, END_YEAR + 1))
N_HOURS = 175320

def cell_weights(lat_grid, lon_grid, tlat, tlon):
    """Find bracketing cell + bilinear weights (u along j, v along i).

    lat_grid/lon_grid must be indexed [i, j] where i is the first axis of the
    raw wspd array (the dim *named* 'lon' in the file but actually the WRF ny
    axis) and j the second (the dim named 'lat').  This matches how Eva's
    CDO remap consumed the same file.
    """
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([lon_grid.ravel(), lat_grid.ravel()]))
    idx = tree.query(np.array([[tlon, tlat]]))[1][0]
    i0, j0 = np.unravel_index(idx, lat_grid.shape)
    i0 = min(max(i0, 0), lat_grid.shape[0] - 2)
    j0 = min(max(j0, 0), lat_grid.shape[1] - 2)
    lon00, lat00 = lon_grid[i0, j0], lat_grid[i0, j0]
    lon01, lat01 = lon_grid[i0, j0+1], lat_grid[i0, j0+1]
    lon10, lat10 = lon_grid[i0+1, j0], lat_grid[i0+1, j0]
    A = np.array([[lon01-lon00, lon10-lon00],
                  [lat01-lat00, lat10-lat00]])
    b = np.array([tlon - lon00, tlat - lat00])
    det = A[0,0]*A[1,1] - A[0,1]*A[1,0]
    if abs(det) < 1e-15:
        print('  WARNING: degenerate cell, nearest-neighbor only')
        return (i0, j0), 0.0, 0.0
    u, v = np.linalg.solve(A, b)
    u = min(max(u, 0.0), 1.0)
    v = min(max(v, 0.0), 1.0)
    return (i0, j0), float(u), float(v)

def main():
    # --- grid meta ---
    ds = xr.open_dataset(GRID_FILE)
    lat_raw = ds['lat'].values   # (300,300), first axis = axis 0 of wspd
    lon_raw = ds['lon'].values
    ds.close()
    lat_grid, lon_grid = lat_raw, lon_raw   # keep raw axis order
    print(f'grid meta loaded, {N_HOURS} hours')

    # --- targets ---
    eccc = pd.read_csv(os.path.join(BASE, 'files', 'ECCC_d03_stations.csv'), header=None)
    eccc_ids = eccc.iloc[:, 4].astype(str).str.strip()
    eccc_lats = eccc.iloc[:, 7].astype(float)
    eccc_lons = eccc.iloc[:, 8].astype(float)

    existing_ids = [118, 145, 155, 192, 202, 277, 327, 348, 381, 640, 664,
                    702, 889, 1001, 1053, 1056, 1275]

    # embedded coords of each CDO file (what the remap actually used)
    embedded = {}
    for sid in existing_ids:
        f = os.path.join(WRF_ST_DIR, f'wind_d03_st{int(sid)}.nc')
        d2 = xr.open_dataset(f)
        embedded[sid] = (float(d2['lat'].values.ravel()[0]),
                         float(d2['lon'].values.ravel()[0]))
        d2.close()

    lat_951, lon_951 = 49.27, -121.23

    targets = {}
    for sid in existing_ids:
        targets[f'cdo_{sid}'] = embedded[sid]     # validate at embedded coords
    targets['new_951'] = (lat_951, lon_951)

    weights = {}
    for name, (tlat, tlon) in targets.items():
        weights[name] = cell_weights(lat_grid, lon_grid, tlat, tlon)
        (i0, j0), u, v = weights[name]
        print(f'{name}: cell (i0={i0},j0={j0}) u={u:.4f} v={v:.4f}  '
              f'target=({tlat},{tlon})')

    # --- stream once, extract all targets ---
    print('Streaming grid...')
    ds = xr.open_dataset(GRID_FILE)
    series = {name: np.empty(N_HOURS, dtype=np.float64) for name in targets}
    pos = 0
    for y in YEARS:
        sel = ds['wspd'].sel(time=slice(f'{y}-01-01', f'{y}-12-31T23:00')).values
        n = sel.shape[0]
        for name, ((i0, j0), u, v) in weights.items():
            series[name][pos:pos+n] = ((1-u)*(1-v)*sel[:, i0, j0]
                                       + u*(1-v)*sel[:, i0, j0+1]
                                       + (1-u)*v*sel[:, i0+1, j0]
                                       + u*v*sel[:, i0+1, j0+1])
        pos += n
        print(f'  {y}: {n} hours', flush=True)
    ds.close()
    assert pos == N_HOURS

    dates = pd.date_range('1986-01-01', periods=N_HOURS, freq='h')
    mask10 = (dates >= pd.Timestamp(START_YEAR, 1, 1)) & \
             (dates < pd.Timestamp(END_YEAR + 1, 1, 1))

    # --- validation at embedded coordinates ---
    print('\nValidation (extracted at CDO-embedded coords):')
    corrs, annual_diffs = [], []
    for sid in existing_ids:
        f = os.path.join(WRF_ST_DIR, f'wind_d03_st{int(sid)}.nc')
        d2 = xr.open_dataset(f)
        cdo = d2['wspd'].values.ravel()
        d2.close()
        our = series[f'cdo_{sid}']
        m = min(len(cdo), len(our))
        d = np.abs(cdo[:m] - our[:m])
        corr = np.corrcoef(cdo[:m], our[:m])[0, 1]
        c10 = pd.Series(cdo[:m], index=dates[:m])[mask10].mean()
        o10 = pd.Series(our[:m], index=dates[:m])[mask10].mean()
        corrs.append(corr); annual_diffs.append(o10 - c10)
        print(f'  {sid}: corr={corr:.4f}  mean|diff|={d.mean():.4f} m/s  '
              f'annual(86-05): CDO={c10:.4f} ours={o10:.4f} diff={o10-c10:+.4f}')
    print(f'\nmean corr = {np.mean(corrs):.4f} (min {np.min(corrs):.4f}), '
          f'mean annual diff = {np.mean(annual_diffs):+.4f} m/s, '
          f'max |annual diff| = {np.max(np.abs(annual_diffs)):.4f} m/s')

    # --- diagnostics on raw grid quality ---
    (i0, j0), u, v = weights['new_951']
    ds3 = xr.open_dataset(GRID_FILE)
    blk = ds3['wspd'].isel(time=slice(0, 3)).values
    sample = ds3['wspd'].isel(time=slice(0, 1)).values.ravel()
    print(f'\n951 cell (i0={i0}, j0={j0}) u={u:.4f} v={v:.4f}')
    print(f'raw grid sample: min={sample.min():.3f} max={sample.max():.3f} '
          f'mean={sample.mean():.3f}, n<0={np.sum(sample<0)}, n>25={np.sum(sample>25)}')
    ds3.close()

    o951 = pd.Series(series['new_951'], index=dates)[mask10].mean()
    print(f'\n951 (true coords 49.27,-121.23): 1986-2005 annual mean = {o951:.4f} m/s')

    # --- decision gate ---
    if np.mean(corrs) < 0.99 or np.max(np.abs(annual_diffs)) > 0.15:
        print('\nValidation FAILED: method does not reproduce CDO remaps.')
        print('Not saving output.')
        sys.exit(1)
    print('\nValidation PASSED.')

    # --- save 951 ---
    wspd_951 = series['new_951']
    out_file = os.path.join(WRF_ST_DIR, 'wind_d03_st951.nc')
    import netCDF4 as nc
    nc_file = nc.Dataset(out_file, 'w', format='NETCDF4')
    nc_file.createDimension('time', None)
    nc_file.createDimension('y', 1)
    nc_file.createDimension('x', 1)
    t = nc_file.createVariable('time', 'f8', ('time',))
    t.units = 'hours since 1986-01-01 00:00:00'
    t.calendar = 'standard'
    t.standard_name = 'time'
    t.axis = 'T'
    t[:] = np.arange(N_HOURS, dtype=np.float64)
    lon_v = nc_file.createVariable('lon', 'f8', ('y', 'x'))
    lon_v.standard_name = 'longitude'; lon_v.units = 'degrees_east'
    lon_v[:] = np.array([[[lon_951]]])
    lat_v = nc_file.createVariable('lat', 'f8', ('y', 'x'))
    lat_v.standard_name = 'latitude'; lat_v.units = 'degrees_north'
    lat_v[:] = np.array([[[lat_951]]])
    w = nc_file.createVariable('wspd', 'f8', ('time', 'y', 'x'))
    w.coordinates = 'lat lon'; w.long_name = 'wind speed'; w.units = 'm s-1'
    w[:] = wspd_951[:, None, None]
    nc_file.history = ('Extracted 2026-09-17 by extract_wrf_wind_951.py: '
                       'bilinear interpolation on wspd_d03_hourly.nc (native '
                       'd03 300x300 grid; same grid the Dec-2023 CDO genbic '
                       'remaps used). Validated at CDO-embedded coords for '
                       '17 stations (mean hourly corr > 0.99, max annual-mean '
                       'diff < 0.15 m/s). Station 951 coords 49.27,-121.23 '
                       'from files/ECCC_d03_stations.csv = PFM ECCC_951.txt.')
    nc_file.close()
    print(f'Saved {out_file}')

if __name__ == '__main__':
    main()
