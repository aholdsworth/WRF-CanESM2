#!/usr/bin/env python3
"""
buoy_wind_bias.py — buoy wind-speed bias metrics, 5 model columns (2026-09-22).

Implements the locked decisions in buoy_wind_bias_handoff_20260922.md §2:
  * obs  = sqrt(V_north^2 + U_east^2) from the NOAA MB_<ID>_HM.mat files,
    restricted to 1986-01-01 .. 2005-12-31, exact-zero U/V = missing
    (Eva's convention), raw ~15-min samples -> hourly mean of each hour's
    pair. Height factor x ln(10/z0)/ln(5/z0) = 1.0684 applied to OBS only
    (5 m anemometer -> 10 m), per SalishSea/WRF/Historical2023/plot_WS.ipynb
    cell 18 getSpeed; model series are used natively at 10 m (no scaling).
  * z0 = 0.0002 (same z0 as SalishSea/WRF/Historical2023/WindEvaluation.ipynb).
  * comparison level = DAILY values (hourly WRF -> mean of surviving hours;
    daily CanESM2/CanRCM4 -> as-is).
  * pooling rules (2026-09-25; same as the land pipeline):
      Rule 1 (annual): mean of day values, year valid only if >= 90% of the
             year's expected days are present.  A day is excluded if >10% of
             its hours are missing or any gap > 3 h.
      Rule 2 (seasonal): a season-year is valid only if all 3 of its months
             are each >= 90% of their calendar days (meteorological DJF).
  * Apples-to-apples: model annuals are always computed on the SAME
    surviving-year set as the obs (shared-year intersection).
  * metrics per buoy x model: n_years, obs period mean, model period mean,
    bias % = 100*(m-o)/o, MAE, RMSE (on the annual values).
    pooled row = mean over the CORE set: buoys with >= --min-years (default 8)
    shared surviving years (handoff's 15 is unattainable for these buoys —
    obs records are 1989-2005 with ~10% wind NaN, so 15 clean years do not
    exist).

Model sources (per-buoy 1x1 files):
  WRF d01  data/wrf_stations/wind_d01_{ECCC,NOAA}_buoy_<ID>.nc
           !! rotated filenames: match by EMBEDDED coords (±0.05 deg),
              never by filename (handoff §4). 46131 has no file -> NaN.
              46207 is duplicated (46132-file + 46207-file): use the
              46207-file, flagged.
  WRF d02  data/wrf_stations/wind_d02_st<ID>.nc   (M11, 2026-09-22)
  WRF d03  data/wrf_stations/wind_d03_st<ID>.nc   (M11, 2026-09-22)
  CanESM2  data/cdo_extractions/canesm2_raw_buoys/sfcWind_buoy_<ID>.nc
           (daily, 365_day, days since 1850-1-1; mask to 1986-2005)
  CanRCM4  data/cdo_extractions/canrcm4_buoys/sfcWind_buoy_<ID>.nc
           (daily, 365_day, days since 1949-12-01)

Outputs (per --tag suffix, e.g. --tag _strict for the locked 4/4 run):
  pickles/buoy_wind_tables<tag>.pkl   (per-buoy x model annual/seasonal tables)
  files/buoy_wind_bias<tag>.csv       (main table)
  files/buoy_wind_bias_pooled<tag>.csv
  stdout: surviving-years report (headline), main table, pooled row.

Usage:  python buoy_wind_bias.py [--min-years 8] [--tag _strict]
                                [--tag _strict] [--obs-root DIR]

Environment: py_2024 only. No network. Read-only obs + PFM.
"""
import argparse
import glob
import os
import pickle
import re
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seasonal_pooling as sp  # ANNUAL_GATE, SEASONS, YEARS, MIN_YEARS, SEASON_OF_MONTH, ...

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
Z0 = 0.0002
H_OBS, H_MODEL = 5.0, 10.0
# Height factor: scale OBS (5 m anemometer) UP to model height (10 m), per
# SalishSea/WRF/Historical2023/plot_WS.ipynb cell 18 (getSpeed):
#   SS = (ln(10/z0)/ln(5/z0)) * hypot(U,V), labelled "wind speed at 10m".
# H5 = ln(10/z0)/ln(5/z0) = 1.0684 (obs 5m -> model 10m), applied to OBS.
H5 = np.log(H_MODEL / Z0) / np.log(H_OBS / Z0)   # = 1.0684, obs(5m) -> model(10m)
MATLAB_EPOCH = 719529.0   # datenum of 1970-01-01 00:00 (classic MATLAB epoch)

D03_BOX = (43.40286, 54.66032, -133.75034, -116.614975)  # la_lo,la_hi,lo_lo,lo_hi

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
MODELS = ('WRF_d01', 'WRF_d02', 'WRF_d03', 'CanESM2', 'CanRCM4')
BASE = os.path.dirname(os.path.abspath(__file__))


def _in_d03_box(lon, lat):
    la_lo, la_hi, lo_lo, lo_hi = D03_BOX
    return (la_lo <= lat <= la_hi) and (lo_lo <= lon <= lo_hi)


# ---------------------------------------------------------------------------
# Obs: .mat parsing
# ---------------------------------------------------------------------------
def parse_mat(path):
    """Return (DatetimeIndex, V, U) for one MB_<ID>_HM.mat.

    Time column 0 is a MATLAB datenum (days+fraction since 1970-01-01;
    classic MATLAB epoch datenum(1970,1,1) = 719529).  Verified against
    header START/END TIME (46207: 1989-10-18 10:00 -> 2016-10-10 22:00,
    matches to the second).  Data are ~2.5-min samples.
    Col 1 = northward V, col 2 = eastward U (m/s).
    """
    m = loadmat(path, squeeze_me=True)
    keys = list(m.keys())
    vk = [k for k in keys if 'VALUES' in k][0]
    v = m[vk]
    dnum = np.asarray(v[:, 0], dtype=float)
    # Snap to the hour the sample represents: 15-min datenums land ~3.3 us
    # off the true hour (23:59:59.999996662 vs 00:00:00.000003338), which
    # silently splits hours across resample bins.
    secs = (dnum - MATLAB_EPOCH) * 86400.0
    times = pd.to_datetime(np.round(secs / 3600.0) * 3600.0, unit='s')
    return times, np.asarray(v[:, 1], float), np.asarray(v[:, 2], float)


def obs_window(times, V, U):
    """Restrict to 1986-2005; exact-zero U or V = missing (Eva's convention).
    Returns (DatetimeIndex, speed m/s array)."""
    mask = (times >= pd.Timestamp(1986, 1, 1)) & (times < pd.Timestamp(2006, 1, 1))
    t = times[mask]
    Vv, Uu = V[mask].copy(), U[mask].copy()
    bad = (Vv == 0.0) | (Uu == 0.0) | np.isnan(Vv) | np.isnan(Uu)
    spd = np.full(len(t), np.nan)
    ok = ~bad
    spd[ok] = np.hypot(Vv[ok], Uu[ok])
    return pd.DatetimeIndex(t), spd


# ---------------------------------------------------------------------------
# Model series loaders (m/s at 10 m, on a datetime index)
# ---------------------------------------------------------------------------
def nc_hourly_series(path):
    d = nc.Dataset(path)
    t = d.variables['time']
    u = t.units
    tv = t[:]
    v = d.variables['wspd'][:, 0, 0]
    d.close()
    m = re.match(r'(\w+)s? since (\d{4})-(\d{1,2})-(\d{1,2})', u)
    unit, y, mo, dd = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
    epoch = pd.Timestamp(y, mo, dd)
    times = epoch + pd.to_timedelta(tv, unit=unit)
    return pd.Series(v, index=times)


def nc_daily_series(path, var='sfcWind'):
    d = nc.Dataset(path)
    t = d.variables['time']
    u = t.units
    tv = t[:]
    fill = d.variables[var]._FillValue if hasattr(d.variables[var], '_FillValue') else None
    v = d.variables[var][:, 0, 0]
    d.close()
    if fill is not None:
        v = np.where(v == fill, np.nan, v)
    m = re.match(r'(\w+)s? since (\d{4})-(\d{1,2})-(\d{1,2})', u)
    unit, y, mo, dd = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
    epoch = pd.Timestamp(y, mo, dd)
    times = epoch + pd.to_timedelta(tv, unit=unit)
    s = pd.Series(v, index=times)
    return s[(s.index >= pd.Timestamp(1986, 1, 1)) & (s.index < pd.Timestamp(2006, 1, 1))]


def wrf_d01_match(buoy_id):
    """Find the d01 wind file whose EMBEDDED coords match truth (±0.05 deg).
    Returns (basename, note) or (None, reason)."""
    tlon, tlat = TRUTH[buoy_id]
    cands = []
    for f in glob.glob(os.path.join(BASE, 'data/wrf_stations/wind_d01_*.nc')):
        d = nc.Dataset(f)
        la = float(d.variables['lat'][:].ravel()[0])
        lo = float(d.variables['lon'][:].ravel()[0])
        d.close()
        if abs(la - tlat) <= 0.05 and abs(lo - tlon) <= 0.05:
            cands.append((os.path.basename(f), la, lo))
    if not cands:
        return None, 'no d01 file at true location'
    if len(cands) > 1:
        own = [c for c in cands if c[0].endswith('_' + buoy_id + '.nc')]
        pick = own[0] if own else cands[0]
        others = ','.join(c[0] for c in cands if c != pick)
        return pick[0], ('DUPLICATE (%d files): using %s; others: %s'
                         % (len(cands), pick[0], others))
    c = cands[0]
    note = 'embedded-coord match'
    if not c[0].endswith('_' + buoy_id + '.nc'):
        note += ' (rotated filename: %s holds buoy %s location)' % (c[0], buoy_id)
    return c[0], note


# ---------------------------------------------------------------------------
# Cascade: hourly -> daily -> month -> season -> year (handoff §2)
# ---------------------------------------------------------------------------
def daily_from_hourly(hourly, gate=0.75, per_day=24):
    """Daily value = mean of the day's hours; day valid if >= gate*per_day
    hours are present. No intra-day gap filter (matches the land pipeline:
    seasonal_pooling annual_frame only gates on days-present).

    BUOY-ONLY gate (2026-09-26): default day-gate 0.75 (18 h of 24), deliberately
    LOOSER than the land pipeline's 0.90.  Buoys have many short days and the
    strict 90% day-gate was dropping ~2,200 days/buoy.  Local to the buoy
    scripts; seasonal_pooling.py (land) is untouched."""
    h = hourly.dropna()
    h = h[(h.index >= pd.Timestamp(1986, 1, 1)) & (h.index < pd.Timestamp(2006, 1, 1))]
    if h.empty:
        return pd.Series(dtype=float)
    h = h[~h.index.duplicated(keep='last')]
    g = h.groupby(h.index.floor('D'))
    days = g.mean()
    counts = g.count()
    ok = counts >= gate * per_day
    return days[ok]


def cascade_to_annual(daily, gate=0.80, min_months=10):
    """Rule 1 + Rule 2 pooling (2026-09-25).

    Rule 1 (annual): mean of the day values, year valid only if >= gate of
      the year's expected days are present.
    Rule 2 (seasonal): a season-year is valid only if all 3 of its months are
      each >= gate of their calendar days (a monthly cascade, so "present"
      means 90%-complete, not merely non-empty).  Meteorological DJF.
    Returns (annual, seas, dropped).
    """
    seas = pd.DataFrame(np.nan, index=list(sp.YEARS), columns=list(sp.SEASONS))
    seas.index.name = 'year'
    d = daily.dropna()
    if d.empty:
        return (pd.Series(np.nan, index=list(sp.YEARS)), seas,
                [(y, 'no days') for y in sp.YEARS])
    d = d.groupby(d.index.normalize()).mean()   # merge same-day duplicates

    # ---- Rule 1: annual = direct day->year, >= gate of days AND >= min_months ----
    annual = pd.Series(np.nan, index=list(sp.YEARS))
    import calendar as _cal
    yr = pd.Series(d.index.year, index=d.index)
    for y in sp.YEARS:
        sel = (yr == y)
        n_exp = 366 if _cal.isleap(y) else 365
        if int(sel.sum()) >= gate * n_exp:
            # month-count floor: require data in >= min_months of the 12 months
            if d.index[sel.values].month.nunique() < min_months:
                continue
            annual.loc[y] = float(d.values[sel.values].mean())
    annual.name = 'annual'

    # ---- Rule 2: seasonal = all 3 months present (each >= gate) ----
    cal = d.groupby(d.index.to_period('M'))
    nhave = cal.count()
    ncal = {p: (p.to_timestamp() + pd.offsets.MonthEnd(0)).day
            for p in nhave.index}
    months = cal.mean()
    months = months[[p for p in months.index if nhave[p] >= gate * ncal[p]]]

    mo = pd.Series(months.index.month, index=months.index)
    ay = pd.Series(months.index.year + (mo == 12), index=months.index)
    seas_name = pd.Series(mo, index=months.index).map(
        lambda m: sp.SEASON_OF_MONTH[m][0])
    for y in sp.YEARS:
        for st in sp.SEASONS:
            msk = (ay == y) & (seas_name == st)
            if not msk.any():
                continue
            if not set(sp.SEASON_MONTHS[st]) <= set(mo[msk].tolist()):
                continue            # a whole month of the season is absent
            seas.loc[y, st] = float(months.values[msk].mean())

    keep = annual.notna()
    dropped = [(y, 'year <%.0f%% complete (Rule 1)' % (100 * gate))
               for y in sp.YEARS if not keep.loc[y]]
    return annual, seas, dropped


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-years', type=int, default=8,
                    help='min shared surviving years for a buoy to enter the '
                         'CORE pool (default 8; handoff default was 15, which '
                         'is unattainable for these buoys)')
    ap.add_argument('--gate', type=float, default=0.80,
                    help='Rule-1 annual gate (fraction of the year days that must '
                         'be present; buoy default 0.80, relaxed from land 0.90)')
    ap.add_argument('--day-gate', type=float, default=0.75,
                    help="fraction of a day's 24 h that must be present for the "
                         "day to be valid (buoy default 0.75 = 18 h)")
    ap.add_argument('--min-months', type=int, default=10,
                    help='min distinct months with data for a year to be valid '
                         '(buoy default 10)')
    ap.add_argument('--tag', default='',
                    help='suffix added to output filenames (e.g. _strict)')
    ap.add_argument('--obs-root',
                    default=os.path.join(BASE, 'data', 'buoy_obs'),
                    help='dir with MB_<ID>_HM.mat (local copies in '
                         'data/buoy_obs; source of record: '
                         '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/DATA/BuoyData)')
    args = ap.parse_args()

    os.makedirs(os.path.join(BASE, 'pickles'), exist_ok=True)
    os.makedirs(os.path.join(BASE, 'files'), exist_ok=True)

    rows, surv, pk, surv_obs = [], [], {}, []
    d01_cache = {}  # basename -> hourly series (29 files, loaded once each... cache by path)

    def load_hourly_cached(path):
        if path not in d01_cache:
            d01_cache[path] = nc_hourly_series(path)
        return d01_cache[path]

    for bid in sorted(TRUTH):
        tlon, tlat = TRUTH[bid]
        if not _in_d03_box(tlon, tlat):
            print('%s: outside D03 box — skipped' % bid, file=sys.stderr)
            continue
        mat = os.path.join(args.obs_root, 'MB_%s_HM.mat' % bid)
        if not os.path.exists(mat):
            print('%s: no .mat — skipped' % bid, file=sys.stderr)
            continue
        times, V, U = parse_mat(mat)
        obs_idx, obs_spd = obs_window(times, V, U)
        if np.all(np.isnan(obs_spd)):
            print('%s: all-NaN wind 1986-2005 — skipped' % bid, file=sys.stderr)
            continue

        # Obs are ~15-min samples (NOT 2.5-min; handoff §4 "hourly" is wrong).
        # Hourly value = mean of that hour's raw pair (Eva's resample('H')
        # convention). The 24-h day gate in daily_from_hourly then applies to
        # genuine hours.
        obs_hourly = pd.Series(obs_spd * H5, index=obs_idx)  # 5m -> 10m
        obs_hourly = obs_hourly.resample('h').mean()
        obs_daily = daily_from_hourly(obs_hourly, gate=args.day_gate)
        obs_annual, obs_seas, obs_dropped = cascade_to_annual(obs_daily, gate=args.gate, min_months=args.min_months)
        surv_years = [int(y) for y in sp.YEARS if not np.isnan(obs_annual.loc[y])]
        surv_obs.append(dict(buoy=bid, n_obs_years=len(surv_years),
                             obs_years=''.join(str(y)[-2:] for y in surv_years)))

        flags, series = {}, {}
        fname, note = wrf_d01_match(bid)
        if fname is None:
            series['WRF_d01'] = pd.Series(dtype=float)
            flags['WRF_d01'] = note
        else:
            series['WRF_d01'] = load_hourly_cached(
                os.path.join(BASE, 'data/wrf_stations', fname))
            flags['WRF_d01'] = note
        for dom in ('d02', 'd03'):
            # 46005 is OFF the D03 grid (~0.97 deg west of the SW corner). The
            # nearest-cell weight fallback that CDO 1.9.9 `remap` applied for it
            # mis-mapped to the NE corner cell (verified: extracted value
            # matches raw cell (159,299) to 0.0000 m/s, not the intended
            # (160,0)), so the D03 value for 46005 is a wrong-cell artifact.
            # Drop 46005 from D03 only; keep it in D01 and D02 (user 2026-09-22).
            if bid == '46005' and dom == 'd03':
                series['WRF_d03'] = pd.Series(dtype=float)
                flags['WRF_d03'] = 'excluded: off D03 domain (wrong nearest cell)'
                continue
            p = os.path.join(BASE, 'data/wrf_stations/wind_%s_st%s.nc' % (dom, bid))
            if os.path.exists(p):
                series['WRF_%s' % dom] = nc_hourly_series(p)
                flags['WRF_%s' % dom] = ''
            else:
                series['WRF_%s' % dom] = pd.Series(dtype=float)
                flags['WRF_%s' % dom] = 'file missing'
        for model, sub in (('CanESM2', 'canesm2_raw_buoys'), ('CanRCM4', 'canrcm4_buoys')):
            p = os.path.join(BASE, 'data/cdo_extractions', sub,
                             'sfcWind_buoy_%s.nc' % bid)
            if os.path.exists(p):
                series[model] = nc_daily_series(p)
                flags[model] = ''
            else:
                series[model] = pd.Series(dtype=float)
                flags[model] = 'file missing'

        for model in MODELS:
            s = series[model]
            if s.empty:
                rows.append(dict(buoy=bid, model=model, n_years=0, obs=np.nan,
                                 model_mean=np.nan, bias_pct=np.nan, mae=np.nan,
                                 rmse=np.nan, flag=flags[model]))
                surv.append(dict(buoy=bid, model=model, n_years=0, years='',
                                 flag=flags[model]))
                continue
            hourly = (s.index[1] - s.index[0]) < pd.Timedelta(hours=2)
            m_daily = daily_from_hourly(s, per_day=1, gate=args.gate) if not hourly else daily_from_hourly(s, gate=args.day_gate)
            m_annual, m_seas, _ = cascade_to_annual(m_daily, gate=args.gate, min_months=args.min_months)
            obs_a = obs_annual.reindex(surv_years)
            mod_a = m_annual.reindex(surv_years)
            ok = obs_a.notna() & mod_a.notna()
            n = int(ok.sum())
            if n >= 2:  # metrics on shared surviving years; 2 = minimal
                # meaningful annual-mean sample (n=1 is a single year, no
                # period mean; reported with NaN metrics, n shown)
                dm = (mod_a - obs_a)[ok]
                r = dict(buoy=bid, model=model, n_years=n,
                         obs=float(obs_a[ok].mean()), model_mean=float(mod_a[ok].mean()),
                         bias_pct=float(100.0 * dm.mean() / obs_a[ok].mean()),
                         mae=float(dm.abs().mean()),
                         rmse=float(np.sqrt((dm ** 2).mean())), flag=flags[model])
            else:
                r = dict(buoy=bid, model=model, n_years=n,
                         obs=float(obs_a[ok].mean()) if n else np.nan,
                         model_mean=float(mod_a[ok].mean()) if n else np.nan,
                         bias_pct=np.nan, mae=np.nan, rmse=np.nan, flag=flags[model])
            rows.append(r)
            surv.append(dict(buoy=bid, model=model, n_years=n,
                             years=''.join(str(y)[-2:] for y in surv_years) if n else '',
                             flag=flags[model]))
            pk['%s|%s' % (bid, model)] = dict(
                annual_obs=obs_annual, annual_model=m_annual,
                seas_obs=obs_seas, seas_model=m_seas,
                obs_dropped=obs_dropped)

    tab = pd.DataFrame(rows).set_index(['buoy', 'model'])
    survdf = pd.DataFrame(surv).set_index(['buoy', 'model'])
    surv_obs_df = pd.DataFrame(surv_obs).set_index('buoy')

    # Pool = mean over the CORE set: buoys with >= min_years shared surviving
    # years (default 8; the handoff's 15 is unattainable for these buoys —
    # obs records are 1989-2005 with ~10% wind NaN, so 15 clean years do not
    # exist). Core-set pooling is the apples-to-apples headline row.
    # Pinned pooled buoy set = the buoys present in the WRF D03 column.
    # Pinning the set to D03 makes the pooled obs baseline identical across all
    # five model columns (previously D03 pooled over fewer buoys, e.g. 46005 is
    # off the D03 domain and only D03 drops it).  Each model column pools over
    # this same set; a column missing a buoy simply skips that buoy.
    core = {r['buoy'] for r in rows if r['n_years'] >= args.min_years}
    d03_ids = {r['buoy'] for r in rows if r['model'] == 'WRF_d03' and r['n_years'] >= 1}
    pool_ids = sorted(core & d03_ids)
    pool_rows = []
    for model in MODELS:
        idx = [(b, model) for b in pool_ids if (b, model) in tab.index and
               tab.loc[(b, model), 'n_years'] >= 1]
        if not idx:
            continue
        sub = tab.loc[idx]
        obs_pm, mod_pm = sub['obs'].mean(), sub['model_mean'].mean()
        pool_rows.append(dict(model=model, n_buoys=len(sub), obs=obs_pm,
                              model_mean=mod_pm,
                              bias_pct=100.0 * (mod_pm - obs_pm) / obs_pm,
                              mae=sub['mae'].mean(), rmse=sub['rmse'].mean()))
    pool = pd.DataFrame(pool_rows).set_index('model') if pool_rows else pd.DataFrame()

    tag = args.tag
    tab.to_csv(os.path.join(BASE, 'files/buoy_wind_bias%s.csv' % tag))
    pool.to_csv(os.path.join(BASE, 'files/buoy_wind_bias_pooled%s.csv' % tag))
    with open(os.path.join(BASE, 'pickles/buoy_wind_tables%s.pkl' % tag), 'wb') as fh:
        pickle.dump(dict(tab=tab, pool=pool, per_buoy=pk,
                         min_years=args.min_years, H5=float(H5), truth=TRUTH), fh)

    print('\nBuoy wind bias — 5 model columns — obs 1986-2005 — '
          'obs 5m->10m x%.4f (z0=%.5f) — Rule 1 annual (day>=%.0f%%, year>=%.0f%%, >=%d months)\n'
          % (H5, Z0, 100 * args.day_gate, 100 * args.gate, args.min_months))
    print('=== SURVIVING YEARS (headline; Rule 1, day>=%.0f%%, year>=%.0f%%, >=%d months) ==='
          % (100 * args.day_gate, 100 * args.gate, args.min_months))
    o = surv_obs_df.reset_index()
    print(o[['buoy', 'n_obs_years', 'obs_years']].rename(
        columns={'n_obs_years': 'n_years', 'obs_years': 'years'}).to_string(index=False))

    print('\n=== MAIN TABLE (metrics on shared surviving years; '
          'metrics computed if n>=2) ===')
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(tab.reset_index()[['buoy', 'model', 'n_years', 'obs', 'model_mean',
                                 'bias_pct', 'mae', 'rmse', 'flag']].to_string(index=False))

    print('\n=== POOLED (pinned buoy set = WRF D03 column: %s) ==='
          % ', '.join(pool_ids))
    if len(pool):
        print(pool.to_string())
    print('\nWrote files/buoy_wind_bias%s.csv, files/buoy_wind_bias_pooled%s.csv, '
          'pickles/buoy_wind_tables%s.pkl' % (tag, tag, tag))


if __name__ == '__main__':
    main()
