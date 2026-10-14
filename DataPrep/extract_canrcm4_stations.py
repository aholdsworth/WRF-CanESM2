#!/usr/bin/env python
"""DEPRECATED (2026-09-21): do NOT run — produced M6 (data/canrcm4_stations/),
deleted 2026-09-21. Its bilinear extraction used wrong station coords for
most NOAA stations; the live pipeline now uses M10 CDO per-station files
(data/cdo_extractions/canrcm4_noaa113/, see station_cdo/ and DATA_INVENTORY.md).
Kept for reference only.

Extract CanRCM4 (rotated-pole) gridded fields at the 175 evaluation
stations — Python replacement for Eva's CDO genbic/remap (bicubic) step.

The CanRCM4 NAM-22 files (data/canrcm4/{tas,pr,wind}_NAM22_hist.nc) are on a
rotated-pole grid with 2-D lat/lon.  We sample each station with **bilinear
interpolation in the native 2-D lat/lon mesh** (equivalent to CDO remapping
on a curvilinear target; the pole transformation is implicit in the 2-D
lat/lon).  This replaces the broken nearest-point-in-native-lat/lon code
that produced a uniform spurious ~+3.1 C bias in the CanRCM4 rows of
pickles/station_tables.pkl.

Conventions (identical to Eva's get_canrcm4 in canesm2_eval_funcs.py and to
gridded_station_values() in make_station_tables.py for CanESM2):
  t     : tas - 273.15, 10-yr (1986-2005) mean            -> C
  pr    : daily flux kg m-2 s-1, clip >=0, *86400 -> mm/day,
          sum per calendar year, mean over 1986-2005        -> mm/yr
  wind  : sfcWind as-is, 10-yr mean                         -> m/s

Output: data/canrcm4_stations/annual_{tas,pr,wind}.nc, one variable per file,
shape (nstation,), plus nearest-neighbor values for verification.
Station order = all_ids of make_station_tables.py (46 ECCC + 18 BCH +
109 NOAA land, then 7 ECCC + 4 NOAA buoys).

Run: /home/amh001/space_fs7/software_2022/python/py_2024/bin/python extract_canrcm4_stations.py
"""
import os
import sys
import time
import numpy as np
import pandas as pd
from netCDF4 import Dataset

def _workroot(start):
    """Nearest ancestor of `start` (inclusive) that is the analysis working
    tree, i.e. the directory containing data/ and pickles/.  Works both when
    this repo is a subfolder of the HPC working tree and when a fresh clone
    is unpacked with the data archive (DATA.md) at the repo root.
    Override with the WRF_CANESM2_BASE environment variable."""
    d = os.path.abspath(start)
    for _ in range(6):
        if (os.path.isdir(os.path.join(d, 'data'))
                and os.path.isdir(os.path.join(d, 'pickles'))):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    raise SystemExit(
        'WRF-CanESM2 working tree not found: no ancestor of %s contains '
        'data/ + pickles/.  Extract the data archive there (see DATA.md) '
        'or set WRF_CANESM2_BASE.' % start)


BASE = os.environ.get('WRF_CANESM2_BASE') or _workroot(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Evaluation', 'elevation_bias'))  # station_geo.py (repo copy)
DATA = os.path.join(BASE, 'data')
OUT_DIR = os.path.join(DATA, 'canrcm4_stations')
os.makedirs(OUT_DIR, exist_ok=True)

START_YEAR, END_YEAR = 1986, 2005

# ---------------------------------------------------------------------------
# Station metadata — SAME all_ids ordering as make_station_tables.py
# ---------------------------------------------------------------------------
from station_geo import (load_eccc, load_bch, load_noaa_land,
                        load_buoys)

eccc_ids, eccc_lats, eccc_lons, eccc_elev = load_eccc()
bch_ids,  bch_lats,  bch_lons,  bch_elev  = load_bch()
noaa_ids, noaa_lats, noaa_lons, noaa_elev = load_noaa_land('NOAA_d03_stations.csv')
eccc_buoy_ids,  ecb_lats,  ecb_lons,  ecb_elev  = load_buoys('ECCC_buoys.csv')
noaa_buoy_ids, nob_lats, nob_lons, nob_elev = load_buoys('NOAA_buoys.csv')

land_ids = [str(i) for i in eccc_ids] + [str(i) for i in bch_ids] + list(noaa_ids)
buoy_ids = [str(i) for i in eccc_buoy_ids] + [str(i) for i in noaa_buoy_ids]
all_ids  = land_ids + buoy_ids

_lats = pd.concat([eccc_lats, bch_lats, noaa_lats, ecb_lats, nob_lats])
_lats.index = _lats.index.astype(str); _lats = _lats.sort_index()
_lons = pd.concat([eccc_lons, bch_lons, noaa_lons, ecb_lons, nob_lons])
_lons.index = _lons.index.astype(str); _lons = _lons.sort_index()
_elev = pd.concat([eccc_elev, bch_elev, noaa_elev, ecb_elev, nob_elev])
_elev.index = _elev.index.astype(str); _elev = _elev.sort_index()

lats_s = _lats.astype(float)
lons_s = _lons.astype(float)
print(f'{len(all_ids)} stations ({len(land_ids)} land, {len(buoy_ids)} buoys)')

# ---------------------------------------------------------------------------
# The CanRCM4 native 2-D lat/lon mesh
# ---------------------------------------------------------------------------
nc = Dataset(os.path.join(DATA, 'canrcm4', 'tas_NAM22_hist.nc'))
lat2d = np.asarray(nc.variables['lat'][:], dtype=np.float64)
lon2d = np.asarray(nc.variables['lon'][:], dtype=np.float64)
ny, nx = lat2d.shape
lon2d = np.where(lon2d > 180.0, lon2d - 360.0, lon2d)
t_all = nc.variables['time'][:]
units = nc.variables['time'].units
nc.close()

from cftime import num2date
dates = num2date(t_all, units)
year_of = np.array([d.year for d in dates])
mask = (year_of >= START_YEAR) & (year_of <= END_YEAR)
n_t = len(t_all)
print(f'grid {ny}x{nx}, {n_t} times, window {int(mask.sum())} days (1986-2005)')

# Sanity: each lat column and each lon row strictly monotonic (already verified
# in prior debugging: 0 bad cols / 0 bad rows).

def find_cell(lats, lons):
    """For each station find its enclosing (or nearest) grid cell.

    Returns (i0, j0, ta, tj, outside) where bilinear weights are
    (1-ta)(1-tj), ta(1-tj), (1-ta)tj, ta*tj at corners
    (i0,j0),(i0+1,j0),(i0,j0+1),(i0+1,j0+1).  Stations outside the grid
    clamp to the border (ta or tj = 0 at the nearest edge).
    """
    n = len(lats)
    i0 = np.empty(n, dtype=int); j0 = np.empty(n, dtype=int)
    ta = np.empty(n);            tj = np.empty(n)
    outside = np.zeros(n, dtype=bool)
    for k in range(n):
        la, lo = lats[k], lons[k]
        # Global nearest grid point, then the cell corner pair around it.
        d = (lat2d - la) ** 2 + (lon2d - lo) ** 2
        idx = np.unravel_index(int(np.argmin(d)), d.shape)
        i, j = idx
        # pick the 2x2 cell whose corners surround (or are nearest to) the
        # station: test the 4 cells touching point (i, j)
        def cell_cost(i0_, j0_):
            a = la - lat2d[i0_, j0_];     b = la - lat2d[i0_+1, j0_]
            c = la - lat2d[i0_, j0_+1];   dd = la - lat2d[i0_+1, j0_+1]
            e = lo - lon2d[i0_, j0_];     f = lo - lon2d[i0_+1, j0_]
            g = lo - lon2d[i0_, j0_+1];   h = lo - lon2d[i0_+1, j0_+1]
            return (a*a + b*b + c*c + dd*dd + e*e + f*f + g*g + h*h) / 8.0
        cands = [(max(i-1, 0), j), (i, max(j-1, 0)),
                 (max(i-1, 0), max(j-1, 0)), (i, j)]
        cands = [(a, b) for a, b in cands if a <= ny-2 and b <= nx-2]
        i0k, j0k = min(cands, key=lambda p: cell_cost(*p))
        i0[k], j0[k] = i0k, j0k
        la_a = lat2d[i0k, j0k];     la_b = lat2d[i0k + 1, j0k]
        lo_a = lon2d[i0k, j0k];     lo_b = lon2d[i0k, j0k + 1]
        ta[k] = 0.0 if la_b == la_a else float(np.clip((la - la_a) / (la_b - la_a), 0.0, 1.0))
        tj[k] = 0.0 if lo_b == lo_a else float(np.clip((lo - lo_a) / (lo_b - lo_a), 0.0, 1.0))
        # distance to the *nearest* of the 4 corners
        c0 = np.hypot(la - lat2d[i0k, j0k],     lo - lon2d[i0k, j0k])
        c1 = np.hypot(la - lat2d[i0k + 1, j0k], lo - lon2d[i0k, j0k])
        c2 = np.hypot(la - lat2d[i0k, j0k + 1], lo - lon2d[i0k, j0k + 1])
        c3 = np.hypot(la - lat2d[i0k + 1, j0k + 1], lo - lon2d[i0k, j0k + 1])
        outside[k] = min(c0, c1, c2, c3) > 1.0   # >~110 km from nearest corner
    return i0, j0, ta, tj, outside

i0, j0, ta, tj, outside = find_cell(lats_s.values, lons_s.values)
print(f'outside-grid stations: {outside.sum()}')
if outside.any():
    for k in np.where(outside)[0]:
        print('   ', all_ids[k], lats_s.values[k], lons_s.values[k])

# ---------------------------------------------------------------------------
# Sample each variable, build annual (1986-2005) statistics
# ---------------------------------------------------------------------------
W00 = (1 - ta) * (1 - tj)
W10 = ta * (1 - tj)
W01 = (1 - ta) * tj
W11 = ta * tj
n_out = len(all_ids)
twin = np.where(mask)[0]
years_w = year_of[twin]
ys = np.unique(years_w)
assert set(ys) == set(range(START_YEAR, END_YEAR + 1)), ys

# ---- quick validation before the full sampling run ----
nc = Dataset(os.path.join(DATA, 'canrcm4', 'tas_NAM22_hist.nc'))
ta0 = np.squeeze(nc.variables['tas'][0, :, :])
t1990 = int(np.argmin(np.abs(year_of - 1990)))
ta90 = np.squeeze(nc.variables['tas'][t1990, :, :])
nc.close()
# Vancouver-area station 889 (50.13, -122.95): 1990-01-01 tas should be cold
k = all_ids.index('889')
v90 = ta90[i0[k], j0[k]] - 273.15
v00 = ta0[i0[k], j0[k]] - 273.15
print(f'VALIDATION station {all_ids[k]} @ ({lats_s.values[k]:.2f},{lons_s.values[k]:.2f}): '
      f'1949-12-31 {v00:.1f} C, 1990-01-01 {v90:.1f} C (expect winter, roughly -8..-15 C)')
assert -30 < v90 < 5, 'winter tas sanity check FAILED - grid mapping is wrong'

def bilinear(vals2d):
    """(nt, ny, nx) -> (nt, nstn) bilinear."""
    return (W00[None, :] * vals2d[:, i0, j0]
            + W10[None, :] * vals2d[:, i0 + 1, j0]
            + W01[None, :] * vals2d[:, i0, j0 + 1]
            + W11[None, :] * vals2d[:, i0 + 1, j0 + 1])

def nearest(vals2d):
    """(nt, ny, nx) -> (nt, nstn) nearest grid point (for smoke test)."""
    return vals2d[:, i0, j0]

results = {}
for vf, vn, kind in [('tas_NAM22_hist.nc', 'tas', 't'),
                     ('pr_NAM22_hist.nc', 'pr', 'pr'),
                     ('wind_NAM22_hist.nc', 'sfcWind', 'wind')]:
    t0 = time.time()
    print(f'Processing {vf} ...')
    nc = Dataset(os.path.join(DATA, 'canrcm4', vf))
    data = np.asarray(nc.variables[vn][:], dtype=np.float32)   # (nt, ny, nx)
    nc.close()
    if kind == 'pr':
        data = np.clip(data, 0.0, None)
    bi = bilinear(data)
    ne = nearest(data)
    if kind == 't':
        bi = bi - 273.15
        ne = ne - 273.15
    win_bi = bi[twin]
    win_ne = ne[twin]
    if kind == 'pr':
        dep_bi = win_bi * 86400.0
        dep_ne = win_ne * 86400.0
        out_bi = np.array([dep_bi[years_w == y].sum(axis=0) for y in ys]).mean(axis=0)
        out_ne = np.array([dep_ne[years_w == y].sum(axis=0) for y in ys]).mean(axis=0)
    else:
        out_bi = win_bi.mean(axis=0)   # average over time -> (n_out,)
        out_ne = win_ne.mean(axis=0)
    print(f'  {kind}: {time.time() - t0:.0f}s, mean {out_bi.mean():.3f}, '
          f'n_nan {int(np.isnan(out_bi).sum())}, shapes out_bi={out_bi.shape}, out_ne={out_ne.shape}, ys={len(ys)}')
    assert out_bi.shape == (n_out,) and out_ne.shape == (n_out,), \
        f'aggregation shape bug: {out_bi.shape}, {out_ne.shape}'
    results[kind] = (out_bi, out_ne)
    del data, bi, ne, win_bi, win_ne   # free memory before next variable

# ---------------------------------------------------------------------------
# Smoke test against last thread's rotated-coord checkpoints
# ---------------------------------------------------------------------------
def smoke(name, lat, lon):
    sel = int(np.argmin(np.hypot(lats_s.values - lat, lons_s.values - lon)))
    t_bi, _ = results['t']
    p_bi, _ = results['pr']
    w_bi, _ = results['wind']
    # annual-model values at the station (not bias; bias needs obs)
    print(f'  {name} (stn {all_ids[sel]} @ {lats_s.values[sel]:.3f},{lons_s.values[sel]:.3f}): '
          f't {t_bi[sel]:.2f} C | pr {p_bi[sel]:.0f} mm/yr | wind {w_bi[sel]:.2f} m/s')

print('Model annual values at smoke-test stations:')
smoke('Vancouver',   49.196, -123.189)
smoke('Bella Coola', 52.22,  -123.32)
smoke('Portland',    45.52,  -122.68)

# ---------------------------------------------------------------------------
# Save per-station annual values
# ---------------------------------------------------------------------------
def save(fname, kind, bi, ne):
    nc = Dataset(os.path.join(OUT_DIR, fname), 'w', format='NETCDF4')
    nc.createDimension('station', n_out)
    v = nc.createVariable(kind, 'f8', ('station',))
    v[:] = bi
    vn = nc.createVariable('nearest_' + kind, 'f8', ('station',))
    vn[:] = ne
    st = nc.createVariable('station_id', 'S12', ('station',))
    st[:] = np.array(all_ids, dtype='S12')
    lv = nc.createVariable('lat', 'f8', ('station',)); lv[:] = lats_s.values
    lo = nc.createVariable('lon', 'f8', ('station',)); lo[:] = lons_s.values
    el = nc.createVariable('obs_elev_m', 'f8', ('station',))
    el[:] = _elev.astype(float).values
    nc.close()
    print(f'saved {os.path.join(OUT_DIR, fname)}')

save('annual_tas.nc',  'tas',     results['t'][0],    results['t'][1])
save('annual_pr.nc',   'pr',      results['pr'][0],   results['pr'][1])
save('annual_wind.nc', 'sfcWind', results['wind'][0], results['wind'][1])
print('Done.')