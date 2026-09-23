#!/usr/bin/env python3
"""Build the season-first pooled bias & skill tables (convention of
seasonal_pooling_handoff_20260918.md, implemented in seasonal_pooling.py).

Does NOT touch pickles/station_tables.pkl (kept as the legacy catalogue).
Produces:
  pickles/skill_seasonal.pkl  : dict dataset -> DataFrame (one row per station)
      station, lat, lon, obs_elev_m, elevation_bias_m,
      per var in {t, pr, wind}:
        {var}_bias        : seasonal-pooled model-obs bias
                            (t in C; pr & wind in %)
        {var}_mae, {var}_rmse : over shared surviving years
        {var}_n_years    : n shared surviving years used
        {var}_obs_period : obs period value (season-first pooled)
        {var}_dropped    : '; ' joined "year(reason)" list of dropped obs years
  pickles/skill_summary.pkl / .csv : 15 rows (5 datasets x 3 vars):
      dataset, var, bias, mae, rmse, n_stations, n_years_min, n_years_max

Convention (see seasonal_pooling.py docstring):
  * obs gated per season at 90% completeness -> year -> period (>=15 yrs)
  * model (complete) restricted to the obs surviving years
  * bias = diff of annual means; MAE/RMSE = error-first over shared years

Run:  /home/amh001/space_fs7/software_2022/python/py_2024/bin/python make_skill_tables.py
"""
import os
import sys
import pickle
import datetime
import warnings

import numpy as np
import pandas as pd
from netCDF4 import Dataset

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Root of the analysis working tree (data/, pickles/). Defaults to the HPC
# working tree used for the analysis; override with ELEV_BIAS_BASE.
BASE = os.environ.get('ELEV_BIAS_BASE',
    '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS')
DATA = os.path.join(BASE, 'data')
sys.path.insert(0, BASE)

import seasonal_pooling as sp
from plot_elev_bias_revised import (load_eccc, load_bch, load_noaa_land,
                                    load_buoys, read_geo_em, read_orog,
                                    weighted_model_elevs, _can_weight_file,
                                    modeled_elevs)

OUT_SKILL = os.path.join(BASE, 'pickles', 'skill_seasonal.pkl')
OUT_SUM   = os.path.join(BASE, 'pickles', 'skill_summary.pkl')
OUT_CSV   = os.path.join(BASE, 'pickles', 'skill_summary.csv')
os.makedirs(os.path.join(BASE, 'pickles'), exist_ok=True)

WRF_ST_DIR  = os.path.join(DATA, 'wrf_stations')
RAW_DIR     = os.path.join(DATA, 'canesm2_raw')
# (M6 canrcm4_stations/ deleted 2026-09-21; CanRCM4 rows use M10 via sp.canrcm4_daily_series)

# ---------------------------------------------------------------------------
# Station metadata (same set/order as make_station_tables.py: 174 land + 11 buoys,
# 951 excluded)
# ---------------------------------------------------------------------------
eccc_ids, eccc_lats, eccc_lons, eccc_elev = load_eccc()
bch_ids,  bch_lats,  bch_lons,  bch_elev  = load_bch()
noaa_ids, noaa_lats, noaa_lons, noaa_elev = load_noaa_land('NOAA_d03_stations.csv')
eccc_buoy_ids,  ecb_lats,  ecb_lons,  ecb_elev  = load_buoys('ECCC_buoys.csv')
noaa_buoy_ids, nob_lats, nob_lons, nob_elev = load_buoys('NOAA_buoys.csv')

land_ids = [str(i) for i in eccc_ids] + [str(i) for i in bch_ids] + list(noaa_ids)
buoy_ids = [str(i) for i in eccc_buoy_ids] + [str(i) for i in noaa_buoy_ids]
all_ids  = land_ids + buoy_ids
buoy_set = set(buoy_ids)

_lats = pd.concat([eccc_lats, bch_lats, noaa_lats, ecb_lats, nob_lats])
_lats.index = _lats.index.astype(str); _lats = _lats.sort_index()
_lons = pd.concat([eccc_lons, bch_lons, noaa_lons, ecb_lons, nob_lons])
_lons.index = _lons.index.astype(str); _lons = _lons.sort_index()
_elev = pd.concat([eccc_elev, bch_elev, noaa_elev, ecb_elev, nob_elev])
_elev.index = _elev.index.astype(str); _elev = _elev.sort_index()
lats_s = _lats.astype(float).reindex([str(s) for s in all_ids])
lons_s = _lons.astype(float).reindex([str(s) for s in all_ids])
elev_s = _elev.astype(float).reindex([str(s) for s in all_ids])
assert lats_s.notna().all() and lons_s.notna().all() and elev_s.notna().all()

def agency_of(sid):
    if sid in buoy_ids:
        return None
    if sid.isdigit():
        return 'E'
    if sid.startswith('USC'):
        return 'N'
    return 'B'

# ---------------------------------------------------------------------------
# Topography / elevation bias (same logic as make_station_tables.py)
# ---------------------------------------------------------------------------
from scipy.interpolate import NearestNDInterpolator

def make_interp(lat2d, lon2d, topo2d):
    pts = np.column_stack([np.ravel(lon2d), np.ravel(lat2d)])
    return NearestNDInterpolator(pts, np.ravel(topo2d))

# WRF elevation bias uses the extraction-consistent CDO bicubic weights (the
# orography actually baked into the model data we compare).  Any station without
# a usable T-file weight file (e.g. buoys without T files) falls back to the
# corner nearest-grid-point lookup over HGT_M.  CanESM2/CanRCM4 use the same
# extraction-consistent method with the exact M9/M10 weight files (see
# plot_elev_bias_revised._can_weight_file), falling back to nearest-GP over orog
# only where no weight file exists.
lat_d03, lon_d03, topo_d03 = read_geo_em('geo_em.d03.nc')
lat_d02, lon_d02, topo_d02 = read_geo_em('geo_em.d02.nc')
lat_d01, lon_d01, topo_d01 = read_geo_em('geo_em.d01.nc')
lats_c2, lons_c2, topo_c2 = read_orog('orog_CanESM2.nc')
lats_r4, lons_r4, topo_r4 = read_orog('orog_CanRCM4.nc')
pts_all = np.column_stack([lons_s.values, lats_s.values])
elev_bias = {}

# CanESM2 & CanRCM4: weight-consistent model elevation (fallback: nearest-GP)
for name, (lat2d, lon2d, topo2d) in [('CanESM2', (lats_c2, lons_c2, topo_c2)),
                                     ('CanRCM4', (lats_r4, lons_r4, topo_r4))]:
    mtopo, _w = weighted_model_elevs(name, list(all_ids))
    warn_ids = [s for s in mtopo[mtopo.isna()].index if str(s) not in buoy_set]
    if warn_ids:
        nearest = modeled_elevs(lat2d, lon2d, topo2d)
        for sid in warn_ids:
            print(f'  WARN elev {name}: {sid}: no weight file -> nearest-GP '
                  f'(model_elev={nearest.loc[sid, "elev"]:.1f})', flush=True)
            mtopo.loc[sid] = nearest.loc[sid, 'elev']
    assert not mtopo[[str(s) not in buoy_set for s in mtopo.index]].isna().any(), \
        f'{name}: unresolved land model elevation'
    # Land-only elevation bias (buoys -> NaN): buoys have no meaningful
    # orography and would stack at bias=0 (the vertical-line artefact).
    _eb = mtopo.values - elev_s.values
    _eb[[s in buoy_set for s in all_ids]] = np.nan
    elev_bias[name] = pd.Series(_eb, index=[str(s) for s in all_ids])

# WRF D01/D02/D03: weight-based model elevation with corner fallback.
for name, dom, (clat, clon, ctopo) in [('D03', 'd03', (lat_d03, lon_d03, topo_d03)),
                                       ('D02', 'd02', (lat_d02, lon_d02, topo_d02)),
                                       ('D01', 'd01', (lat_d01, lon_d01, topo_d01))]:
    wmelev, wwarns = weighted_model_elevs(dom, list(all_ids))
    # Corner fallback for stations without a usable T-file weight file.
    # clat/clon are C-grid corners (ny+1 x nx+1) while ctopo (HGT_M) is the
    # A-grid (ny x nx); use cell centres (mean of the 4 surrounding corners)
    # so the points match ctopo's shape.
    clat_c = 0.25 * (clat[:-1, :-1] + clat[1:, :-1] + clat[:-1, 1:] + clat[1:, 1:])
    clon_c = 0.25 * (clon[:-1, :-1] + clon[1:, :-1] + clon[:-1, 1:] + clon[1:, 1:])
    cpts = np.column_stack([clon_c.ravel(), clat_c.ravel()])
    cinterp = NearestNDInterpolator(cpts, ctopo.ravel())
    ctopo_all = cinterp(pts_all)
    for w in wwarns:
        print(f'  WARN elev {name}: {w}', flush=True)
    for i, sid in enumerate(all_ids):
        if sid in buoy_set:
            continue  # buoys excluded from the elevation diagnostic
        if np.isnan(wmelev.iloc[i]):
            wmelev.iloc[i] = ctopo_all[i]
    assert all(not np.isnan(wmelev.iloc[i]) for i, s in enumerate(all_ids)
               if s not in buoy_set), f'{name}: unresolved land model elevation'
    # Land-only elevation bias (buoys -> NaN), same rationale as the Can* rows.
    _eb = wmelev.values - elev_s.values
    _eb[[s in buoy_set for s in all_ids]] = np.nan
    elev_bias[name] = pd.Series(_eb, index=[str(s) for s in all_ids])

# ---------------------------------------------------------------------------
# Step 1: obs side — season/gate/year-mask for every land station & variable
# ---------------------------------------------------------------------------
print('Step 1: obs seasonal pooling (gate + year mask)...', flush=True)
obs_pool_cache = {}    # (sid, kind) -> obs_pool dict or None
obs_drop_report = []   # (sid, agency, kind, dropped list)
for i, sid in enumerate(land_ids):
    ag = agency_of(sid)
    for kind in ('t', 'pr', 'wind'):
        r = sp.obs_pool(sid, ag, kind, DATA)
        obs_pool_cache[(sid, kind)] = r
        if r is None:
            obs_drop_report.append((sid, ag, kind, ['no obs data']))
        elif r['dropped']:
            obs_drop_report.append((sid, ag, kind, r['dropped']))
    if (i + 1) % 50 == 0:
        print(f'  {i+1}/{len(land_ids)}', flush=True)

# ---------------------------------------------------------------------------
# Step 2: WRF model annuals (per station, per domain, per var), masked to the
# obs surviving years (model data complete -> just restrict to those years)
# ---------------------------------------------------------------------------
def wrf_st_files(var, dom):
    out = {}
    d = WRF_ST_DIR
    BCH3 = {'ALU','BCK','BLN','CLO','CMU','CMX','COQ','CQM','CRU','DAI','DLU','ECL',
            'GOC','MIS','NTY','STA','WAH','WOL'}
    for f in os.listdir(d):
        if not f.endswith('.nc') or not f.startswith(var + '_' + dom + '_'):
            continue
        rest = f[len(var + '_' + dom + '_'):-3]
        if dom == 'd01':
            if rest.startswith('ECCC_'):
                out['E' + rest[5:]] = os.path.join(d, f)
            elif rest.startswith('NOAA_'):
                out['N' + rest[5:]] = os.path.join(d, f)
            elif rest.startswith('BCH_'):
                out['B' + rest[4:]] = os.path.join(d, f)
            elif rest in BCH3:
                out['B' + rest] = os.path.join(d, f)
            else:
                import sys
                print(f'WARN wrf_st_files: untagged d01 file skipped: {f}', file=sys.stderr)
        else:
            if rest.startswith('st'):
                sid2 = rest[2:]
                if sid2.startswith('USC') or sid2.startswith('USW'):
                    out['N' + sid2] = os.path.join(d, f)
                else:
                    out['E' + sid2] = os.path.join(d, f)
            elif rest.startswith('BCH_'):
                out['B' + rest[4:]] = os.path.join(d, f)
            else:
                import sys
                print(f'WARN wrf_st_files: untagged {dom} file skipped: {f}', file=sys.stderr)
    return out

def wrf_station_daily(path, var_name, kind):
    """Load one per-station WRF file -> DAILY Series (t: C, pr: mm/day,
    wind: m/s).  pr uses diff of the cumulative field, clipped >= 0."""
    nc = Dataset(path, mode='r')
    t = nc.variables['time'][:]
    data = np.squeeze(nc.variables[var_name][:])
    nc.close()
    dates = pd.to_datetime(datetime.datetime(sp.START_YEAR, 1, 1)) \
              + pd.to_timedelta(t, unit='h')
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
    """dict sid -> Series(annual values, index=years) for WRF per-station files."""
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
# Step 3: gridded CanESM2 model annuals (nearest grid point, daily)
# ---------------------------------------------------------------------------
def canesm2_model_annuals(var):
    """dict sid -> Series(annual, index=years) from per-station CDO daily files
    (M9, station_cdo/): data/cdo_extractions/canesm2_raw_noaa113/<var>_<tag>.nc,
    1850-2005 (56940 d), masked to 1986-2005 (days 49640..56939) before
    annual aggregation.  Rewired from M3 nearest-point on 2026-09-21."""
    from cftime import num2date
    vname = {'t': 'tas', 'pr': 'pr', 'wind': 'sfcWind'}[var]
    lo, hi = sp.M10_EVAL_SLICE['canesm2']
    out = {}
    for sid in land_ids:
        f = os.path.join(DATA, sp.M10_CanESM2_DIR, f'{vname}_{sp._m10_file_tag(sid)}.nc')
        if not os.path.exists(f):
            continue
        try:
            nc = Dataset(f, mode='r')
            v = nc.variables[vname][:].astype(float).ravel()
            t = nc.variables['time'][:]
            units, cal = nc.variables['time'].units, nc.variables['time'].calendar
            nc.close()
            dates = num2date(t[lo:hi], units, calendar=cal)
            s = pd.Series(v[lo:hi], index=pd.DatetimeIndex([
                pd.Timestamp(int(d.year), int(d.month), int(d.day)) for d in dates]))
            if var == 't':
                s = s - 273.15
            elif var == 'pr':
                s = s * 86400.0   # daily-mean kg m-2 s-1 -> mm per day
            out[sid] = sp.model_annual_from_series(s, 'sum' if var == 'pr' else 'mean')
        except Exception as e:
            print(f'  WARN CanESM2 {var} {sid}: {e}', flush=True)
    return out

# ---------------------------------------------------------------------------
# Step 4: CanRCM4 pre-extracted annuals (bilinear per-station annual values;
# model complete -> mean of annuals over surviving years is the period value)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Assemble per-dataset tables
# ---------------------------------------------------------------------------
def combine(model_annual, obs_rec):
    """Apply like-for-like year mask + skill. Returns dict of per-var results
    for one (dataset, station): {'bias','mae','rmse','n_years','obs_period'}.
    model_annual: Series(annual, index=year) or None.
    obs_rec: obs_pool dict or None."""
    if model_annual is None or obs_rec is None:
        return None
    b, mae, rmse, n = sp.skill(model_annual, obs_rec['annual'])
    return dict(bias=b, mae=mae, rmse=rmse, n_years=n,
                obs_period=obs_rec['period'])

def make_table(name, get_model):
    """get_model: dict var -> dict sid -> Series(annual)."""
    rows = {}
    for sid in all_ids:
        row = {}
        for var in ('t', 'pr', 'wind'):
            ma = None
            if sid in land_ids:
                ma = get_model.get(var, {}).get(sid)
            obs_rec = obs_pool_cache.get((sid, var))
            if obs_rec is not None and obs_rec['n_years'] < sp.MIN_YEARS:
                res = None
            else:
                res = combine(ma, obs_rec)
            if res is None:
                for k in ('bias', 'mae', 'rmse', 'n_years', 'obs_period'):
                    row[f'{var}_{k}'] = np.nan
                row[f'{var}_dropped'] = np.nan
            else:
                row[f'{var}_bias'] = res['bias']
                row[f'{var}_mae'] = res['mae']
                row[f'{var}_rmse'] = res['rmse']
                row[f'{var}_n_years'] = res['n_years']
                row[f'{var}_obs_period'] = res['obs_period']
                row[f'{var}_dropped'] = np.nan
        rows[sid] = row
    df = pd.DataFrame.from_dict(rows, orient='index')
    df = df.reindex([str(s) for s in all_ids])
    df.insert(0, 'station', df.index)
    df['lat'] = lats_s.values
    df['lon'] = lons_s.values
    df['obs_elev_m'] = elev_s.values
    df['elevation_bias_m'] = elev_bias[name].values
    # percentage bias for pr & wind (obs period > 0 guard); t stays in C
    for var in ('pr', 'wind'):
        op = df[f'{var}_obs_period'].values
        mb = df[f'{var}_bias'].values
        pct = np.where(np.isnan(op) | np.isnan(mb) | (op <= 0), np.nan,
                       100.0 * mb / op)
        df[f'{var}_bias'] = pct
    return df

tables = {}

# --- WRF d01/d02/d03 ---
print('Step 2: WRF model annuals (d03/d02/d01)...', flush=True)
for name, dom in [('D03', 'd03'), ('D02', 'd02'), ('D01', 'd01')]:
    anns = {}
    for var in ('t', 'pr', 'wind'):
        print(f'  {name} {var}...', flush=True)
        anns[var] = wrf_model_annuals(land_ids, var, dom)
    tables[name] = make_table(name, anns)

# --- CanESM2 ---
print('Step 3: CanESM2 gridded model annuals...', flush=True)
anns = {}
for var in ('t', 'pr', 'wind'):
    print(f'  CanESM2 {var}...', flush=True)
    try:
        anns[var] = canesm2_model_annuals(var)
    except Exception as e:
        print(f'  WARN CanESM2 {var}: {e}', flush=True)
        anns[var] = {}
tables['CanESM2'] = make_table('CanESM2', anns)

# --- CanRCM4: per-station CDO daily (M10, station_cdo/).  Now we have yearly
# model values -> full skill (MAE/RMSE) on the obs surviving years, like WRF.
print('Step 4: CanRCM4 CDO per-station (M10)...', flush=True)
anns = {}
for var in ('t', 'pr', 'wind'):
    out = {}
    for sid in land_ids:
        try:
            s = sp.canrcm4_daily_series(sid, var, DATA)
        except Exception:
            s = None
        if s is None or len(s.dropna()) == 0:
            continue
        out[sid] = sp.model_annual_from_series(s, 'sum' if var == 'pr' else 'mean')
    anns[var] = out
    print(f'  CanRCM4 {var}: {len(out)} stations', flush=True)
tables['CanRCM4'] = make_table('CanRCM4', anns)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
with open(OUT_SKILL, 'wb') as f:
    pickle.dump(tables, f)
print(f'\nSaved {OUT_SKILL}')

# ---------------------------------------------------------------------------
# Domain summary: 15 rows (5 datasets x 3 vars)
# ---------------------------------------------------------------------------
sum_rows = []
for name, df in tables.items():
    for var, unit in (('t', 'C'), ('pr', '%'), ('wind', '%')):
        sub = df[[f'{var}_bias', f'{var}_mae', f'{var}_rmse', f'{var}_n_years']].dropna(
            subset=[f'{var}_bias'])
        if len(sub):
            ny = sub[f'{var}_n_years']
            sum_rows.append({
                'dataset': name, 'var': var, 'unit': unit,
                'bias': sub[f'{var}_bias'].mean(),
                'mae': sub[f'{var}_mae'].mean(),
                'rmse': sub[f'{var}_rmse'].mean(),
                'n_stations': len(sub),
                'n_years_min': int(ny.min()), 'n_years_max': int(ny.max()),
            })
        else:
            sum_rows.append({'dataset': name, 'var': var, 'unit': unit,
                             'bias': np.nan, 'mae': np.nan, 'rmse': np.nan,
                             'n_stations': 0, 'n_years_min': np.nan,
                             'n_years_max': np.nan})
summary = pd.DataFrame(sum_rows,
                       columns=['dataset', 'var', 'unit', 'bias', 'mae', 'rmse',
                                'n_stations', 'n_years_min', 'n_years_max'])
with open(OUT_SUM, 'wb') as f:
    pickle.dump(summary, f)
summary.to_csv(OUT_CSV, index=False)
print(f'Saved {OUT_SUM} and {OUT_CSV}')

# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------
print('\n=== Domain summary (season-first pooled) ===')
print(summary.to_string(index=False,
                        formatters={'bias': '{:.3f}'.format,
                                    'mae': '{:.3f}'.format,
                                    'rmse': '{:.3f}'.format}))

print('\n=== Per-dataset non-NaN bias counts ===')
for name, df in tables.items():
    print(f"{name:9s} t={df['t_bias'].notna().sum():3d} "
          f"pr={df['pr_bias'].notna().sum():3d} "
          f"wind={df['wind_bias'].notna().sum():3d}")

print('\n=== Dropped station-seasons (obs side), first 40 ===')
n_dropped_stations = 0
for sid, ag, kind, dropped in obs_drop_report:
    if kind == 'wind' and ag == 'B':
        continue  # BCH has no wind obs at all (expected)
    if not dropped:
        continue
    n_dropped_stations += 1
    if n_dropped_stations <= 40:
        parts = []
        for it in dropped[:4]:
            try:
                y, r = it
                parts.append(f'{y}: {r}')
            except (TypeError, ValueError):
                parts.append(str(it))
        detail = '; '.join(parts)
        more = f' (+{len(dropped)-4} more)' if len(dropped) > 4 else ''
        print(f'  {sid} [{ag}] {kind}: {detail}{more}')
print(f'  ... total station-variable pairs with dropped years: {n_dropped_stations}')
