#!/usr/bin/env python3
"""
buoy_tas_bias.py — buoy near-surface air-temperature (tas) bias metrics,
3 model columns (2026-09-23).  Mirrors buoy_wind_bias.py exactly
(same 90% cascade, same apples-to-apples shared-year pairing, same pooling,
same metrics), with temperature substitutions:

  * obs  = "Temperature:Air" (column 4, degrees C) from the NOAA
    MB_<ID>_HM.mat files, restricted to 1986-01-01 .. 2005-12-31,
    NaN = missing (no exact-zero convention for temperature),
    raw ~15-min samples -> hourly mean of each hour's pair -> daily mean.
    Converted to K (+273.15) to match the model units.
    Sensor is the buoy met-package air-temperature probe (~1.5-3 m above
    the sea surface); model fields are 2 m (WRF T2) / near-surface
    (CanESM2/CanRCM4 tas).  No height scaling is applied to temperature
    (log-profile scaling is for momentum, not a conserved scalar).
  * model sources (per-buoy 1x1 files, all >= daily):
      WRF D01    data/wrf_stations/t_d01_{ECCC,NOAA}_buoy_<ID>.nc
                 (re-extracted 2026-09-23 from t_d01_hourly.nc, verified
                  embedded coords = truth; the original PFM pull had
                  rotated ECCC filenames - see extract_buoy_t_d01.sh)
      CanESM2    data/cdo_extractions/canesm2_raw_buoys/tas_buoy_<ID>.nc
                 (daily, 1850-2005, sliced to 1986-2005; K)
      CanRCM4    data/cdo_extractions/canrcm4_buoys/tas_buoy_<ID>.nc
                 (daily, 1986-2005; K)
      WRF D02/D03: no buoy t extractions yet -> columns are pending and
                  omitted from this run.
  * D01 75 km sub-grid placement caveat: buoy position up to ~0.35 deg
    from the nearest cell centre (documented, same as wind D01).

Outputs (per --tag suffix):
  files/buoy_tas_bias<tag>.csv          (per-buoy x model)
  files/buoy_tas_bias_pooled<tag>.csv   (pooled over the >= --min-years core)
  stdout: surviving-years report + main table + pooled row.

Usage:  python buoy_tas_bias.py [--min-years 8] [--year-gate 3|4]
                                [--tag _strict] [--obs-root DIR]

Environment: py_2024 only. No network. Read-only obs.
"""
import argparse
import os
import re
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seasonal_pooling as sp  # GATE_FRAC, SEASONS, YEARS, ...

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MATLAB_EPOCH = 719529.0   # datenum of 1970-01-01 00:00 (classic MATLAB epoch)
AIR_COL = 4               # "Temperature:Air", degrees C (col 0 time, 1 V, 2 U,
                          # 3 pressure, 4 air T, 5 water T, ...)

TRUTH = {  # .mat-header coords (lon, lat in deg), station_cdo/buoys/buoy_descs.csv
    '46131': (-124.99, 49.91), '46132': (-127.93, 49.74),
    '46134': (-123.5, 48.65),  '46146': (-123.73, 49.34),
    '46204': (-128.77, 51.38), '46206': (-126.0, 48.83),
    '46207': (-129.91, 50.88),
    '46005': (-131.09, 46.143), '46029': (-124.487, 46.163),
    '46041': (-124.739, 47.352),
    '46050': (-124.514667, 44.6775), '46087': (-124.726333, 48.493),
    '46088': (-123.164667, 48.333667), '46089': (-125.819167, 45.893333),
}
MODELS = ('WRF_d01', 'CanESM2', 'CanRCM4')
BASE = os.path.dirname(os.path.abspath(__file__))

# D01 agency map (which t_d01_<AGENCY>_buoy_<ID>.nc file each buoy uses)
D01_AGENCY = {
    '46131': 'ECCC', '46132': 'ECCC', '46134': 'ECCC', '46146': 'ECCC',
    '46204': 'ECCC', '46206': 'ECCC', '46207': 'ECCC',
    '46005': 'NOAA', '46029': 'NOAA', '46041': 'NOAA',
    '46050': 'NOAA', '46087': 'NOAA', '46088': 'NOAA', '46089': 'NOAA',
}


# ---------------------------------------------------------------------------
# Obs: .mat parsing (temperature)
# ---------------------------------------------------------------------------
def parse_mat(path):
    """Return (DatetimeIndex, air_temp degC) for one MB_<ID>_HM.mat."""
    m = loadmat(path, squeeze_me=True)
    vk = [k for k in m.keys() if 'VALUES' in k][0]
    v = m[vk]
    dnum = np.asarray(v[:, 0], dtype=float)
    # Snap to the hour the sample represents (15-min datenums land ~3.3 us
    # off the true hour; same snap as the wind pipeline).
    secs = (dnum - MATLAB_EPOCH) * 86400.0
    times = pd.to_datetime(np.round(secs / 3600.0) * 3600.0, unit='s')
    air = np.asarray(v[:, AIR_COL], dtype=float)
    return times, air


def obs_window(times, air):
    """Restrict to 1986-2005.  Returns (DatetimeIndex, degC series)."""
    mask = (times >= pd.Timestamp(1986, 1, 1)) & (times < pd.Timestamp(2006, 1, 1))
    return times[mask], pd.Series(air[mask], index=times[mask])


# ---------------------------------------------------------------------------
# Model series loaders (K, on a datetime index)
# ---------------------------------------------------------------------------
def _units_to_index(t):
    u = t.units
    m = re.match(r'(\w+)s? since (\d{4})-(\d{1,2})-(\d{1,2})', u)
    unit, y, mo, dd = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
    epoch = pd.Timestamp(y, mo, dd)
    return epoch + pd.to_timedelta(t[:], unit=unit)


def nc_hourly_t2(path):
    d = nc.Dataset(path)
    t = _units_to_index(d.variables['time'])
    v = d.variables['T2'][:, 0, 0]
    d.close()
    return pd.Series(v, index=t)


def nc_tas_daily(path, var='tas'):
    d = nc.Dataset(path)
    t = _units_to_index(d.variables['time'])
    fv = d.variables[var]._FillValue if hasattr(d.variables[var], '_FillValue') else None
    v = d.variables[var][:, 0, 0]
    d.close()
    if fv is not None:
        v = np.where(v == fv, np.nan, v)
    s = pd.Series(v, index=t)
    return s[(s.index >= pd.Timestamp(1986, 1, 1)) & (s.index < pd.Timestamp(2006, 1, 1))]


def wrf_d01_t2(buoy_id):
    """Locate the D01 T2 file whose EMBEDDED coords match truth (±0.05 deg).
    Returns (path, note) or (None, reason).  The 2026-09-23 re-extraction
    should make the by-name file correct, but we verify by coords anyway."""
    tlon, tlat = TRUTH[buoy_id]
    d0 = os.path.join(BASE, 'data/wrf_stations')
    cands = []
    for f in os.listdir(d0):
        m = re.match(r't_d01_(ECCC|NOAA)_buoy_(\w+)\.nc$', f)
        if not m:
            continue
        fp = os.path.join(d0, f)
        d = nc.Dataset(fp)
        la = float(d.variables['lat'][:].ravel()[0])
        lo = float(d.variables['lon'][:].ravel()[0])
        d.close()
        if abs(la - tlat) <= 0.05 and abs(lo - tlon) <= 0.05:
            cands.append((f, la, lo))
    if not cands:
        return None, 'no d01 file at true location'
    own = [c for c in cands if c[0].endswith('_' + buoy_id + '.nc')]
    if len(own):
        return os.path.join(d0, own[0][0]), 'embedded-coord match'
    pick = cands[0][0]
    others = ','.join(c[0] for c in cands if c[0] != pick)
    return os.path.join(d0, pick), (
        'embedded-coord match (rotated filename: %s holds %s; others: %s)'
        % (pick, buoy_id, others))


# ---------------------------------------------------------------------------
# Cascade: hourly -> daily -> month -> season -> year (identical to wind)
# ---------------------------------------------------------------------------
def daily_from_hourly(hourly, gate=sp.GATE_FRAC, per_day=24, max_gap_h=3):
    h = hourly.dropna()
    h = h[(h.index >= pd.Timestamp(1986, 1, 1)) & (h.index < pd.Timestamp(2006, 1, 1))]
    if h.empty:
        return pd.Series(dtype=float)
    h = h[~h.index.duplicated(keep='last')]
    g = h.groupby(h.index.floor('D'))
    days = g.mean()
    counts = g.count()
    ok = counts >= gate * per_day
    if per_day >= 2:
        day0 = h.index.floor('D')
        keep_days = day0[day0.isin(days.index[ok])].unique()
        if len(keep_days):
            grid = h.reindex(
                pd.DatetimeIndex(np.repeat(keep_days, per_day)) +
                pd.to_timedelta(np.tile(np.arange(per_day), len(keep_days)), unit='h'))
            miss = grid.isna().values
            run = np.zeros(len(miss), dtype=int)
            cur = 0
            for i in range(len(miss)):
                if i % per_day == 0:
                    cur = 0
                if miss[i]:
                    cur += 1
                    run[i] = cur
                else:
                    cur = 0
            maxrun = pd.Series(run, index=grid.index).groupby(
                grid.index.normalize()).max()
            bad = maxrun[maxrun > max_gap_h].index
            ok = ok.reindex(days.index, fill_value=False)
            ok.loc[bad.intersection(days.index)] = False
    return days[ok]


def cascade_to_annual(daily, gate=sp.GATE_FRAC, n_seasons=4):
    seas = pd.DataFrame(np.nan, index=list(sp.YEARS), columns=list(sp.SEASONS))
    seas.index.name = 'year'
    d = daily.dropna()
    if d.empty:
        return (pd.Series(np.nan, index=list(sp.YEARS)), seas,
                [(y, 'no days') for y in sp.YEARS])
    d = d.groupby(d.index.normalize()).mean()
    cal = d.groupby(d.index.to_period('M'))
    nhave = cal.count()
    ncal = {p: (p.to_timestamp() + pd.offsets.MonthEnd(0)).day
            for p in nhave.index}
    months = cal.mean()
    months = months[[p for p in months.index if nhave[p] >= gate * ncal[p]]]

    seas_of_m = pd.Series(months.index.month, index=months.index).map(sp.SEASON_OF_MONTH)
    seas_of_y = pd.Series(months.index.year, index=months.index)
    for y in sp.YEARS:
        for st in sp.SEASONS:
            msk = ((seas_of_y == y) & (seas_of_m == st))
            if int(msk.sum()) >= 2:   # 90% of 3 months -> at least 2 of 3
                seas.loc[y, st] = float(months.values[msk].mean())

    n_st = seas.notna().sum(axis=1)
    keep = n_st >= n_seasons
    annual = seas.mean(axis=1).where(keep)
    annual.name = 'annual'
    dropped = [(y, 'season(s) <90%% complete (%d/%d): %s' % (
                  int(n_st.loc[y]), len(sp.SEASONS), ','.join(
                      st for st in sp.SEASONS if np.isnan(seas.loc[y, st]))))
               for y in sp.YEARS if not keep.loc[y]]
    return annual, seas, dropped


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-years', type=int, default=8,
                    help='min shared surviving years for a buoy to enter the CORE pool')
    ap.add_argument('--year-gate', type=int, default=3, choices=(3, 4),
                    help='min passing seasons for a surviving year (3=relaxed, 4=strict)')
    ap.add_argument('--tag', default='', help='output filename suffix')
    ap.add_argument('--obs-root',
                    default='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/DATA/BuoyData')
    args = ap.parse_args()

    os.makedirs(os.path.join(BASE, 'files'), exist_ok=True)
    os.makedirs(os.path.join(BASE, 'pickles'), exist_ok=True)

    rows, surv, surv_obs = [], [], []
    t2_cache = {}

    def load_t2_cached(path):
        if path not in t2_cache:
            t2_cache[path] = nc_hourly_t2(path)
        return t2_cache[path]

    for bid in sorted(TRUTH):
        mat = os.path.join(args.obs_root, 'MB_%s_HM.mat' % bid)
        if not os.path.exists(mat):
            print('%s: no .mat — skipped' % bid, file=sys.stderr)
            continue
        times, air = parse_mat(mat)
        t_idx, air_obs = obs_window(times, air)
        if np.all(np.isnan(air_obs.values)):
            print('%s: all-NaN air T 1986-2005 — skipped' % bid, file=sys.stderr)
            continue

        # ~15-min samples -> hourly mean of the hour's pair (Eva's convention).
        # Convert to K NOW so obs and model share units end-to-end.
        obs_hourly = (air_obs + 273.15).resample('h').mean()
        obs_daily = daily_from_hourly(obs_hourly)
        obs_annual, obs_seas, obs_dropped = cascade_to_annual(
            obs_daily, n_seasons=args.year_gate)
        surv_years = [int(y) for y in sp.YEARS if not np.isnan(obs_annual.loc[y])]
        surv_obs.append(dict(buoy=bid, n_obs_years=len(surv_years),
                             obs_years=''.join(str(y)[-2:] for y in surv_years)))

        flags, series = {}, {}
        # WRF D01 (hourly T2, K) — matched by embedded coords
        path, note = wrf_d01_t2(bid)
        if path is None:
            series['WRF_d01'] = pd.Series(dtype=float)
            flags['WRF_d01'] = note
        else:
            series['WRF_d01'] = load_t2_cached(path)
            flags['WRF_d01'] = note
        # CanESM2 + CanRCM4 (daily tas, K)
        for model, sub in (('CanESM2', 'canesm2_raw_buoys'),
                           ('CanRCM4', 'canrcm4_buoys')):
            p = os.path.join(BASE, 'data/cdo_extractions', sub,
                             'tas_buoy_%s.nc' % bid)
            if os.path.exists(p):
                series[model] = nc_tas_daily(p)
                flags[model] = ''
            else:
                series[model] = pd.Series(dtype=float)
                flags[model] = 'file missing'

        for model in MODELS:
            s = series[model]
            if s.empty:
                rows.append(dict(buoy=bid, model=model, n_years=0, obs=np.nan,
                                 model_mean=np.nan, bias=np.nan, bias_pct=np.nan,
                                 mae=np.nan, rmse=np.nan, flag=flags[model]))
                surv.append(dict(buoy=bid, model=model, n_years=0, years='',
                                 flag=flags[model]))
                continue
            hourly = (s.index[1] - s.index[0]) < pd.Timedelta(hours=2)
            m_daily = daily_from_hourly(s) if hourly else daily_from_hourly(s, per_day=1)
            m_annual, m_seas, _ = cascade_to_annual(m_daily, n_seasons=args.year_gate)
            obs_a = obs_annual.reindex(surv_years)
            mod_a = m_annual.reindex(surv_years)
            ok = obs_a.notna() & mod_a.notna()
            n = int(ok.sum())
            if n >= 2:
                dm = (mod_a - obs_a)[ok]
                r = dict(buoy=bid, model=model, n_years=n,
                         obs=float(obs_a[ok].mean()),
                         model_mean=float(mod_a[ok].mean()),
                         bias=float(dm.mean()),
                         bias_pct=float(100.0 * dm.mean() / obs_a[ok].mean()),
                         mae=float(dm.abs().mean()),
                         rmse=float(np.sqrt((dm ** 2).mean())),
                         flag=flags[model])
            else:
                r = dict(buoy=bid, model=model, n_years=n,
                         obs=float(obs_a[ok].mean()) if n else np.nan,
                         model_mean=float(mod_a[ok].mean()) if n else np.nan,
                         bias=np.nan, bias_pct=np.nan, mae=np.nan, rmse=np.nan,
                         flag=flags[model])
            rows.append(r)
            surv.append(dict(buoy=bid, model=model, n_years=n,
                             years=''.join(str(y)[-2:] for y in surv_years) if n else '',
                             flag=flags[model]))

    tab = pd.DataFrame(rows).set_index(['buoy', 'model'])
    survdf = pd.DataFrame(surv).set_index(['buoy', 'model'])
    surv_obs_df = pd.DataFrame(surv_obs).set_index('buoy')

    pool_ids = sorted({r['buoy'] for r in rows if r['n_years'] >= args.min_years})
    pool_rows = []
    for model in MODELS:
        idx = [(b, model) for b in pool_ids if (b, model) in tab.index]
        if not idx:
            continue
        sub = tab.loc[idx]
        obs_pm, mod_pm = sub['obs'].mean(), sub['model_mean'].mean()
        pool_rows.append(dict(model=model, n_buoys=len(sub), obs=obs_pm,
                              model_mean=mod_pm, bias=mod_pm - obs_pm,
                              bias_pct=100.0 * (mod_pm - obs_pm) / obs_pm,
                              mae=sub['mae'].mean(), rmse=sub['rmse'].mean()))
    pool = pd.DataFrame(pool_rows).set_index('model') if pool_rows else pd.DataFrame()

    tag = args.tag
    tab.to_csv(os.path.join(BASE, 'files/buoy_tas_bias%s.csv' % tag))
    pool.to_csv(os.path.join(BASE, 'files/buoy_tas_bias_pooled%s.csv' % tag))
    with open(os.path.join(BASE, 'pickles/buoy_tas_tables%s.pkl' % tag), 'wb') as fh:
        import pickle
        pickle.dump(dict(tab=tab, pool=pool, min_years=args.min_years,
                         truth=TRUTH), fh)

    print('\nBuoy tas (near-surface air T) bias — 3 model columns — '
          'obs 1986-2005, K — year gate: >=%d/4 seasons\n' % args.year_gate)
    print('=== SURVIVING YEARS (headline; obs cascade, 90%% gate) ===')
    o = surv_obs_df.reset_index()
    print(o.to_string(index=False))
    print('\n=== PER BUOY x MODEL ===')
    print(tab.reset_index().to_string(index=False, float_format=lambda x: '%.2f' % x))
    print('\n=== POOLED (core: >= %d shared yrs: %s) ==='
          % (args.min_years, ', '.join(pool_ids)))
    print(pool.to_string(float_format=lambda x: '%.2f' % x))


if __name__ == '__main__':
    main()
