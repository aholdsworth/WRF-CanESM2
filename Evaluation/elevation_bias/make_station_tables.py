#!/usr/bin/env python
"""Build the SINGLE canonical Rule-1 land station table (root producer).

This is now the ONLY root producer for the land verification tables and the
single obs code path for every downstream figure.  Observations are pooled
through seasonal_pooling.obs_pool (Rule 1: 90% annual completeness gate;
period reported only if >= MIN_YEARS=15 of 20 years survive) and the skill
model side is masked to the obs surviving years (seasonal_pooling.skill).
The legacy per-agency ungated obs loaders (_noaa_obs_annual, _eccc_elev_obs,
_bch_obs_annual, _eccc_hourly_wind) were removed -- this file no longer has
its own obs code path.

For each of the 5 datasets (WRF D03, D02, D01, CanESM2, CanRCM4) a DataFrame
is written with, per station:
    station, lat, lon, obs_elev_m,
    elevation_bias_m,      # model grid elevation - station elevation
    elev_method,
    temperature_bias_C     # annual (1986-2005) model T2 - obs T, over shared
                           #   surviving years (deg C)
    precipitation_bias_pct # 100*(model-pr - obs-pr)/obs-pr  (%)
    wind_bias_pct          # 100*(model-wind - obs-wind)/obs-wind  (%)
    per var in {t, pr, wind}:
        {var}_mae, {var}_rmse, {var}_r2, {var}_n_years
The existing bias column names / units are KEPT for consumer compatibility
(t in C; pr & wind in %).  A top-level '_version' tag is added.

Model data:
  * WRF d01/d02/d03 : data/wrf_canon/{t,pr,wind}_{d0X}_{AGENCY}_{ID}.nc
    (hourly per-station; t T2-273.15, pr cumulative->diff clip>=0, wind wspd)
  * CanESM2 / CanRCM4 : data/cdo_extractions/{canesm2_raw,canrcm4}_noaa113/
    per-station CDO daily files (M9/M10).

A station-variable cell is NaN when obs has < MIN_YEARS surviving years OR
the model file is absent.

make_skill_tables.py is now a seconds-long PROJECTION of this table into
skill_seasonal.pkl / skill_summary.pkl with the legacy column names, so
make_bias_figure.py / make_elev_bias_grid.py / make_abs_skill_table.py need
no edits.

Run:
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      /home/amh001/space_fs7/software_2022/python/py_2024/bin/python \
      make_station_tables.py
Output: <BASE>/pickles/station_tables.pkl  (dict of DataFrames + '_version')
"""
import os
import sys
import glob
import datetime
import warnings

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date
from scipy.interpolate import NearestNDInterpolator

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

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
PFM  = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF'
RUN  = 'historical'
START_YEAR, END_YEAR = 1986, 2005
TABLE_VERSION = 'rule1-canonical-land-v1-20260929'

sys.path.insert(0, BASE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # station_geo.py (repo copy)
from station_geo import (load_eccc, load_bch, load_noaa,
                        read_geo_em, read_orog,
                        weighted_model_elevs, _can_weight_file,
                        modeled_elevs)
import seasonal_pooling as sp

DOMAIN_DIR = os.path.join(BASE, 'plotting_files')
DATA = os.path.join(BASE, 'data')
WRF_ST_DIR  = os.path.join(DATA, 'wrf_canon')
OUT_PKL     = os.path.join(BASE, 'pickles', 'station_tables.pkl')
os.makedirs(os.path.dirname(OUT_PKL), exist_ok=True)

# ---------------------------------------------------------------------------
# Station metadata (canonical 173 land-only set, same as the 2026-09-29
# USW-fixed table: 46 ECCC + 18 BCH + 100 USC + 9 USW).  Buoys are scored via
# the validated buoy table pickles elsewhere and are NOT part of this table.
# ---------------------------------------------------------------------------
eccc_ids, eccc_lats, eccc_lons, eccc_elev = load_eccc()
bch_ids,  bch_lats,  bch_lons,  bch_elev  = load_bch()
noaa_ids, noaa_lats, noaa_lons, noaa_elev = load_noaa('NOAA_d03_stations.csv')

land_ids = ([str(i) for i in eccc_ids] + [str(i) for i in bch_ids]
            + [str(i) for i in noaa_ids])
buoy_ids = []
all_ids  = land_ids
buoy_set = set()

_lats = pd.concat([eccc_lats, bch_lats, noaa_lats])
_lats.index = _lats.index.astype(str); _lats = _lats.sort_index()
_lons = pd.concat([eccc_lons, bch_lons, noaa_lons])
_lons.index = _lons.index.astype(str); _lons = _lons.sort_index()
_elev = pd.concat([eccc_elev, bch_elev, noaa_elev])
_elev.index = _elev.index.astype(str); _elev = _elev.sort_index()

lats_s = _lats.astype(float).reindex([str(s) for s in all_ids])
lons_s = _lons.astype(float).reindex([str(s) for s in all_ids])
elev_s = _elev.astype(float).reindex([str(s) for s in all_ids])
assert lats_s.notna().all() and lons_s.notna().all() and elev_s.notna().all(), \
    'station metadata misalignment after reindex'


def agency_of(sid):
    """'E' (ECCC, numeric), 'N' (NOAA USC/USW), 'B' (BCH 3-letter)."""
    if sid.isdigit():
        return 'E'
    if sid.startswith(('USC', 'USW')):
        return 'N'
    return 'B'


# ---------------------------------------------------------------------------
# Topography & elevation bias (unchanged from the legacy table)
# ---------------------------------------------------------------------------
def make_interp(lat2d, lon2d, topo2d):
    pts = np.column_stack([np.ravel(lon2d), np.ravel(lat2d)])
    vals = np.ravel(topo2d)
    return NearestNDInterpolator(pts, vals)

lat_d03, lon_d03, topo_d03 = read_geo_em('geo_em.d03.nc')
lat_d02, lon_d02, topo_d02 = read_geo_em('geo_em.d02.nc')
lat_d01, lon_d01, topo_d01 = read_geo_em('geo_em.d01.nc')
lats_c2, lons_c2, topo_c2 = read_orog('orog_CanESM2.nc')
lats_r4, lons_r4, topo_r4 = read_orog('orog_CanRCM4.nc')
pts_all = np.column_stack([lons_s.values, lats_s.values])

elev_bias = {}
elev_method = {}

for name, (lat2d, lon2d, topo2d) in [('CanESM2', (lats_c2, lons_c2, topo_c2)),
                                     ('CanRCM4', (lats_r4, lons_r4, topo_r4))]:
    mtopo, _w = weighted_model_elevs(name, list(all_ids))
    warn_ids = [s for s in mtopo[mtopo.isna()].index.tolist() if s not in buoy_set]
    if warn_ids:
        nearest = modeled_elevs(lat2d, lon2d, topo2d, lats_s, lons_s, elev_s)
        for sid in warn_ids:
            print(f'  WARN elev {name}: {sid}: no weight file -> nearest-GP '
                  f'(model_elev={nearest.loc[sid, "elev"]:.1f})')
            mtopo.loc[sid] = nearest.loc[sid, 'elev']
    _land = [s for s in all_ids if s not in buoy_set]
    assert mtopo.loc[_land].notna().all(), f'{name}: unresolved land model elevation'
    _eb = mtopo.values - elev_s.values
    _eb[[s in buoy_set for s in all_ids]] = np.nan
    elev_bias[name] = pd.Series(_eb, index=[str(s) for s in all_ids])
    elev_method[name] = pd.Series(
        ['bicubic' if w is not None else 'nearest'
         for w in [_can_weight_file(name, s) for s in all_ids]],
        index=[str(s) for s in all_ids])

for name, dom, (clat, clon, ctopo) in [('D03', 'd03', (lat_d03, lon_d03, topo_d03)),
                                       ('D02', 'd02', (lat_d02, lon_d02, topo_d02)),
                                       ('D01', 'd01', (lat_d01, lon_d01, topo_d01))]:
    wmelev, wwarns = weighted_model_elevs(dom, list(all_ids))
    wmelev_orig = wmelev.copy()
    clat_c = 0.25 * (clat[:-1, :-1] + clat[1:, :-1] + clat[:-1, 1:] + clat[1:, 1:])
    clon_c = 0.25 * (clon[:-1, :-1] + clon[1:, :-1] + clon[:-1, 1:] + clon[1:, 1:])
    cpts = np.column_stack([clon_c.ravel(), clat_c.ravel()])
    cinterp = NearestNDInterpolator(cpts, ctopo.ravel())
    ctopo_all = cinterp(pts_all)
    for w in wwarns:
        print(f'  WARN elev {name}: {w}')
    for i, sid in enumerate(all_ids):
        if sid in buoy_set:
            continue
        if np.isnan(wmelev.iloc[i]):
            wmelev.iloc[i] = ctopo_all[i]
    assert all(not np.isnan(wmelev.iloc[i]) for i, s in enumerate(all_ids)
               if s not in buoy_set), f'{name}: unresolved land model elevation'
    _eb = wmelev.values - elev_s.values
    _eb[[s in buoy_set for s in all_ids]] = np.nan
    elev_bias[name] = pd.Series(_eb, index=[str(s) for s in all_ids])
    elev_method[name] = pd.Series(
        ['bicubic' if not np.isnan(orig) else 'corner'
         for orig in wmelev_orig.values], index=[str(s) for s in all_ids])

# ---------------------------------------------------------------------------
# Model data: WRF per-station hourly files (data/wrf_canon)
# ---------------------------------------------------------------------------
def wrf_st_files(var, dom):
    """station file map for WRF per-station files of one var/domain.
    data/wrf_canon naming: {var}_{dom}_{AGENCY}_{ID}.nc with AGENCY in
    {ECCC, BCH, NOAA}; keys are 'E'+ClimateID / 'B'+code / 'N'+USC*/USW*."""
    out = {}
    d = WRF_ST_DIR
    prefix = var + '_' + dom + '_'
    for f in os.listdir(d):
        if not f.endswith('.nc') or not f.startswith(prefix):
            continue
        rest = f[len(prefix):-3]
        if rest.startswith('ECCC_'):
            out['E' + rest[5:]] = os.path.join(d, f)
        elif rest.startswith('BCH_'):
            out['B' + rest[4:]] = os.path.join(d, f)
        elif rest.startswith('NOAA_'):
            out['N' + rest[5:]] = os.path.join(d, f)
    return out


def wrf_station_daily(path, var_name, kind):
    """Load one per-station WRF file -> DAILY Series (t: C, pr: mm/day,
    wind: m/s).  pr uses diff of the cumulative field, clipped >= 0."""
    nc = Dataset(path, mode='r')
    t = nc.variables['time'][:].astype(float)
    data = np.squeeze(nc.variables[var_name][:])
    nc.close()
    dates = pd.to_datetime(datetime.datetime(START_YEAR, 1, 1)) + pd.to_timedelta(t, unit='h')
    if kind == 't':
        h = pd.Series(data - 273.15, index=dates)
        return h.groupby(h.index.normalize()).mean()
    if kind == 'wind':
        h = pd.Series(data, index=dates)
        return h.groupby(h.index.normalize()).mean()
    # pr: cumulative -> increments, clip <0, drop the first (no prior value)
    inc = np.clip(np.diff(data), 0.0, None)
    d = pd.Series(inc, index=dates[1:].normalize())
    return d.groupby(d.index).sum()


def wrf_model_annuals(ids, var, dom):
    """dict sid -> Series(annual, index=years) for WRF per-station files."""
    files = wrf_st_files(var, dom)
    out = {}
    for sid in ids:
        key = agency_of(sid) + sid
        p = files.get(key)
        if p is None:
            continue
        try:
            daily = wrf_station_daily(p, {'t': 'T2', 'pr': 'pr', 'wind': 'wspd'}[var], var)
            out[sid] = sp.model_annual_from_series(daily, 'sum' if var == 'pr' else 'mean')
        except Exception as e:
            print(f'  WARN WRF {dom} {var} {sid}: {e}', flush=True)
    return out


# ---------------------------------------------------------------------------
# Model data: per-station CDO daily files (CanESM2 M9 / CanRCM4 M10)
# ---------------------------------------------------------------------------
def can_model_annuals(var, model):
    """dict sid -> Series(annual, index=years) from per-station CDO daily
    files in data/cdo_extractions/{canesm2_raw,canrcm4}_noaa113/."""
    vname = {'t': 'tas', 'pr': 'pr', 'wind': 'sfcWind'}[var]
    d = os.path.join(DATA, sp.M10_CanESM2_DIR if model == 'canesm2' else sp.M10_CanRCM4_DIR)
    lo, hi = sp.M10_EVAL_SLICE[model]
    out = {}
    for sid in land_ids:
        f = os.path.join(d, f'{vname}_{sp._m10_file_tag(sid)}.nc')
        if not os.path.exists(f):
            continue
        try:
            nc = Dataset(f, mode='r')
            v = nc.variables[vname][:].astype(float).ravel()
            t = nc.variables['time'][:].astype(float)
            units, cal = nc.variables['time'].units, nc.variables['time'].calendar
            nc.close()
            dates = num2date(t[lo:hi], units, calendar=cal)
            s = pd.Series(v[lo:hi], index=pd.DatetimeIndex([
                pd.Timestamp(int(dd.year), int(dd.month), int(dd.day)) for dd in dates]))
            if var == 't':
                s = s - 273.15
            elif var == 'pr':
                s = s * 86400.0   # daily-mean kg m-2 s-1 -> mm per day
            out[sid] = sp.model_annual_from_series(s, 'sum' if var == 'pr' else 'mean')
        except Exception as e:
            print(f'  WARN {model} {var} {sid}: {e}', flush=True)
    return out


# ---------------------------------------------------------------------------
# Step 1: obs side -- Rule-1 obs_pool for every land station & variable
# ---------------------------------------------------------------------------
print('Step 1: obs pooling (Rule 1, 90% annual gate, >=15 yrs)...', flush=True)
obs_pool_cache = {}     # (sid, kind) -> obs_pool dict or None
for i, sid in enumerate(all_ids):
    ag = agency_of(sid)
    for kind in ('t', 'pr', 'wind'):
        obs_pool_cache[(sid, kind)] = sp.obs_pool(sid, ag, kind, DATA)
    if (i + 1) % 50 == 0:
        print(f'  {i+1}/{len(all_ids)}', flush=True)


def combine(model_annual, obs_rec):
    """Like-for-like skill (model masked to obs surviving years).
    Returns dict or None when obs is absent / < MIN_YEARS years."""
    if model_annual is None or obs_rec is None:
        return None
    if obs_rec['n_years'] < sp.MIN_YEARS:
        return None
    b, mae, rmse, r2, n = sp.skill(model_annual, obs_rec['annual'])
    return dict(bias=b, mae=mae, rmse=rmse, r2=r2, n_years=n,
                obs_period=obs_rec['period'])


def make_table(name, get_model):
    """get_model: dict var -> dict sid -> Series(annual).  Returns a DataFrame
    with the bias columns + per-var mae/rmse/r2/n_years."""
    rows = {}
    for sid in all_ids:
        row = {}
        for var in ('t', 'pr', 'wind'):
            ma = get_model.get(var, {}).get(sid)
            obs_rec = obs_pool_cache.get((sid, var))
            res = combine(ma, obs_rec)
            if res is None:
                for k in ('bias', 'mae', 'rmse', 'r2', 'n_years', 'obs_period'):
                    row[f'{var}_{k}'] = np.nan
            else:
                row[f'{var}_bias'] = res['bias']
                row[f'{var}_mae'] = res['mae']
                row[f'{var}_rmse'] = res['rmse']
                row[f'{var}_r2'] = res['r2']
                row[f'{var}_n_years'] = res['n_years']
                row[f'{var}_obs_period'] = res['obs_period']
        rows[sid] = row
    df = pd.DataFrame.from_dict(rows, orient='index')
    df = df.reindex([str(s) for s in all_ids])
    df.insert(0, 'station', df.index)
    df['lat'] = lats_s.values
    df['lon'] = lons_s.values
    df['obs_elev_m'] = elev_s.values
    df['elevation_bias_m'] = elev_bias[name].values
    df['elev_method'] = elev_method[name].values
    # Canonical bias columns (kept names/units for consumer compatibility):
    # t in C; pr & wind as % (obs period > 0 guard).
    df['temperature_bias_C'] = df['t_bias'].values
    for var in ('pr', 'wind'):
        op = df[f'{var}_obs_period'].values
        mb = df[f'{var}_bias'].values
        df[f'{var}_bias'] = np.where(np.isnan(op) | np.isnan(mb) | (op <= 0),
                                     np.nan, 100.0 * mb / op)
        df[f'precipitation_bias_pct' if var == 'pr' else 'wind_bias_pct'] = \
            df[f'{var}_bias'].values
    return df


tables = {}

# --- WRF d01/d02/d03 ---
for name, dom in [('D03', 'd03'), ('D02', 'd02'), ('D01', 'd01')]:
    anns = {}
    for var in ('t', 'pr', 'wind'):
        print(f'{name} {var} (WRF)...', flush=True)
        anns[var] = wrf_model_annuals(all_ids, var, dom)
    tables[name] = make_table(name, anns)

# --- CanESM2 / CanRCM4 (per-station CDO daily) ---
for name, model in [('CanESM2', 'canesm2'), ('CanRCM4', 'canrcm4')]:
    anns = {}
    for var in ('t', 'pr', 'wind'):
        print(f'{name} {var} (CDO {model})...', flush=True)
        anns[var] = can_model_annuals(var, model)
    tables[name] = make_table(name, anns)

# ---------------------------------------------------------------------------
# Save & report
# ---------------------------------------------------------------------------
out = {'_version': TABLE_VERSION}
out.update(tables)
with open(OUT_PKL, 'wb') as f:
    import pickle as _pk
    _pk.dump(out, f)
print(f'\nSaved {OUT_PKL}')
for name, df in tables.items():
    print(f"\n=== {name} ===")
    print(df[['temperature_bias_C', 'precipitation_bias_pct', 'wind_bias_pct']]
          .notna().sum().to_string())
    print(df[['elevation_bias_m', 'temperature_bias_C', 'precipitation_bias_pct',
              'wind_bias_pct']].mean(numeric_only=True).round(3).to_string())
