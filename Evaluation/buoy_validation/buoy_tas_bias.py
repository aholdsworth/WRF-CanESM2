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
      WRF D02    data/wrf_stations/t_d02_st<ID>.nc  (hourly T2, K; the 6
                 rotated Sep-15 PFM ECCC files were re-extracted 2026-09-23
                 with verified coords - see extract_buoy_t_parallel.sh)
      WRF D03    data/wrf_stations/t_d03_st<ID>.nc  (hourly T2, K; same
                 re-extraction; 46005 EXCLUDED - its nearest D03 cell is
                 1.16 deg away, an off-domain artifact)
  * D01 75 km sub-grid placement caveat: buoy position up to ~0.35 deg
    from the nearest cell centre (documented, same as wind D01).

Outputs (per --tag suffix):
  files/buoy_tas_bias<tag>.csv          (per-buoy x model)
  files/buoy_tas_bias_pooled<tag>.csv   (pooled over the >= --min-years core)
  stdout: surviving-years report + main table + pooled row.

Usage:  python buoy_tas_bias.py [--min-years 8] [--tag _strict] [--obs-root DIR]

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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'elevation_bias'))  # shared seasonal_pooling.py (single copy)
import seasonal_pooling as sp  # ANNUAL_GATE, SEASONS, YEARS, SEASON_OF_MONTH, ...

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
MODELS = ('WRF_d01', 'WRF_d02', 'WRF_d03', 'CanESM2', 'CanRCM4')

# D03 exclusion: 46005's nearest D03 cell is 1.156 deg away (off-domain
# artifact at the domain edge); all other buoys are within 0.018 deg.
D03_EXCLUDE = {'46005': 'off D03 domain (nearest cell 1.16 deg away)'}
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

# Fingerprint of the (day-gate, year-gate, min-months) cascade that a saved
# obs pickle was built with.  The obs daily series in
# pickles/buoy_tas_obs_daily.pkl is gate-dependent, so reading one built
# under different gates would silently give wrong surviving years.
def _obs_cache_fp(day_gate, gate, min_months):
    return {'day_gate': float(day_gate), 'gate': float(gate),
            'min_months': int(min_months)}

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


def wrf_nest_t2(buoy_id, dom, tol=0.05):
    """Locate the T2 file for nest `dom` (d01/d02/d03) whose EMBEDDED coords
    match the buoy truth within `tol`.  Returns (path, note) or (None, reason).
    D01 files are named t_d01_<AGENCY>_buoy_<ID>.nc; D02/D03 files are named
    t_<dom>_st<ID>.nc.  Matching is by embedded coords (not just filename)
    so a rotated/misnamed file can never enter the analysis silently."""
    tlon, tlat = TRUTH[buoy_id]
    d0 = os.path.join(BASE, 'data/wrf_stations')
    pat = (r't_d01_(ECCC|NOAA)_buoy_(\w+)\.nc$' if dom == 'd01'
           else r't_%s_st(\w+)\.nc$' % dom)
    cands = []
    for f in os.listdir(d0):
        m = re.match(pat, f)
        if not m:
            continue
        fp = os.path.join(d0, f)
        d = nc.Dataset(fp)
        la = float(d.variables['lat'][:].ravel()[0])
        lo = float(d.variables['lon'][:].ravel()[0])
        d.close()
        if abs(la - tlat) <= tol and abs(lo - tlon) <= tol:
            cands.append((f, la, lo))
    if not cands:
        return None, 'no %s file at true location' % dom
    own = [c for c in cands if c[0].endswith(buoy_id + '.nc')]
    if len(own):
        return os.path.join(d0, own[0][0]), 'embedded-coord match'
    pick = cands[0][0]
    others = ','.join(c[0] for c in cands if c[0] != pick)
    return os.path.join(d0, pick), (
        'embedded-coord match (rotated filename: %s holds %s; others: %s)'
        % (pick, buoy_id, others))


def wrf_d01_t2(buoy_id):
    return wrf_nest_t2(buoy_id, 'd01')


def wrf_d02_t2(buoy_id):
    return wrf_nest_t2(buoy_id, 'd02')


def wrf_d03_t2(buoy_id):
    # 3 km cell: use a looser tol so we don't reject a file that sits on the
    # right cell but whose centre is ~0.008 deg (half a 3 km cell) from truth.
    return wrf_nest_t2(buoy_id, 'd03', tol=0.2)


# ---------------------------------------------------------------------------
# Cascade: hourly -> daily -> month -> season -> year (identical to wind)
# ---------------------------------------------------------------------------
def daily_from_hourly(hourly, gate=0.75, per_day=24):
    """daily = mean of the hours present; day valid if >= gate*per_day hours.
    No intra-day gap filter (matches the land pipeline: seasonal_pooling
    annual_frame only gates on days-present, never on intra-day gaps).

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
    """BUOY-ONLY annual/seasonal pooling (2026-09-26, relaxed year gate).

    Rule 1 (annual, relaxed): mean of the day values.  A year is valid only if
      BOTH (a) >= gate of the year's expected days are present AND
          (b) the year has data in at least `min_months` of its 12 months.
      The month-count floor keeps the all-12-months intent even when the
      day-count gate is relaxed (a year with 85% complete days but only
      8 months of data is still dropped).
    Rule 2 (seasonal): a season-year is valid only if all 3 of its months are
      each >= gate of their calendar days (a monthly cascade, so "present"
      means gate-complete, not merely non-empty).  Meteorological DJF.

    Defaults are the relaxed buoy values (0.80 / 10 months), NOT the land
    pipeline's ANNUAL_GATE=0.90.  seasonal_pooling.py is left unchanged."""
    seas = pd.DataFrame(np.nan, index=list(sp.YEARS), columns=list(sp.SEASONS))
    seas.index.name = 'year'
    d = daily.dropna()
    if d.empty:
        return (pd.Series(np.nan, index=list(sp.YEARS)), seas,
                [(y, 'no days') for y in sp.YEARS])
    d = d.groupby(d.index.normalize()).mean()   # merge same-day duplicates

    # ---- Rule 1: annual = direct day->year, >= gate of days AND >= min_months ----
    annual = pd.Series(np.nan, index=list(sp.YEARS))
    yr = pd.Series(d.index.year, index=d.index)
    import calendar as _cal
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
    # analysis year: Dec -> following year (DJF(y)=Dec(y-1)+Jan-Feb(y))
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
                    help='min shared surviving years for a buoy to enter the CORE pool')
    ap.add_argument('--gate', type=float, default=0.80,
                    help='Rule-1 annual gate (fraction of the year days that must '
                         'be present; buoy default 0.80, relaxed from land 0.90)')
    ap.add_argument('--day-gate', type=float, default=0.75,
                    help="fraction of a day's 24 h that must be present for the "
                         "day to be valid (buoy default 0.75 = 18 h)")
    ap.add_argument('--min-months', type=int, default=10,
                    help='min distinct months with data for a year to be valid '
                         '(buoy default 10)')
    ap.add_argument('--tag', default='', help='output filename suffix')
    ap.add_argument('--season-pool', action='store_true',
                    help='pool surviving (year, season) pairs with equal '
                         'weight instead of pooling annual means (season '
                         'weighting variant; metrics on >=4 season-pairs)')
    ap.add_argument('--obs-root', default=None,
                    help='dir with MB_<ID>_HM.mat; if unset, the cached obs '
                         'DAILY series are loaded from '
                         'pickles/buoy_tas_obs_daily.pkl and the SAME '
                         'cascade_to_annual is re-run (build with '
                         '--dump-obs-daily). Caching daily, not seasonal '
                         'frames, keeps the two paths identical.')
    ap.add_argument('--dump-obs-daily', action='store_true',
                    help='write pickles/buoy_tas_obs_daily.pkl (obs daily '
                         'series per buoy, with a gate fingerprint) and exit')
    args = ap.parse_args()

    os.makedirs(os.path.join(BASE, 'files'), exist_ok=True)
    os.makedirs(os.path.join(BASE, 'pickles'), exist_ok=True)

    rows, surv, surv_obs = [], [], []
    _DUMP = {}
    t2_cache = {}

    def load_t2_cached(path):
        if path not in t2_cache:
            t2_cache[path] = nc_hourly_t2(path)
        return t2_cache[path]

    for bid in sorted(TRUTH):
        if args.obs_root is not None:
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
            obs_daily = daily_from_hourly(obs_hourly, gate=args.day_gate)
            obs_annual, obs_seas, obs_dropped = cascade_to_annual(obs_daily, gate=args.gate, min_months=args.min_months)
        else:
            # No --obs-root: reuse the cached obs DAILY series and re-run the
            # SAME cascade_to_annual, so this path is byte-identical to the
            # --obs-root path (no lossy season-first averaging).
            _op = os.path.join(BASE, 'pickles/buoy_tas_obs_daily.pkl')
            if not os.path.exists(_op):
                raise SystemExit('no --obs-root and no pickles/buoy_tas_obs_daily.pkl; '
                                 'run once with: python buoy_tas_bias.py '
                                 '--obs-root <dir> --dump-obs-daily')
            import pickle as _pk
            _loaded = _pk.load(open(_op, 'rb'))
            _fp = _loaded.get('fingerprint') if isinstance(_loaded, dict) and 'fingerprint' in _loaded else None
            _want = _obs_cache_fp(args.day_gate, args.gate, args.min_months)
            if _fp is None or _fp != _want:
                raise SystemExit(
                    'pickles/buoy_tas_obs_daily.pkl was built with gates %s '
                    'but this run uses %s. Re-dump it with: '
                    'python buoy_tas_bias.py --obs-root <dir> --dump-obs-daily '
                    '(and the same --day-gate/--gate/--min-months)'
                    % (_fp, _want))
            _daily = _loaded['frames']
            if bid not in _daily:
                print('%s: no obs daily — skipped' % bid, file=sys.stderr)
                continue
            obs_daily = _daily[bid].copy()
            obs_annual, obs_seas, obs_dropped = cascade_to_annual(
                obs_daily, gate=args.gate, min_months=args.min_months)
        if getattr(args, 'dump_obs_daily', False):
            _DUMP[bid] = obs_daily
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
        # WRF D02 / D03 (hourly T2, K) — matched by embedded coords
        for model, loader in (('WRF_d02', wrf_d02_t2), ('WRF_d03', wrf_d03_t2)):
            if model == 'WRF_d03' and bid in D03_EXCLUDE:
                series[model] = pd.Series(dtype=float)
                flags[model] = D03_EXCLUDE[bid]
                continue
            path, note = loader(bid)
            if path is None:
                series[model] = pd.Series(dtype=float)
                flags[model] = note
            else:
                series[model] = load_t2_cached(path)
                flags[model] = note
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
            m_daily = daily_from_hourly(s, gate=args.day_gate) if hourly else daily_from_hourly(s, per_day=1, gate=args.gate)
            m_annual, m_seas, _ = cascade_to_annual(m_daily, gate=args.gate, min_months=args.min_months)
            if getattr(args, 'season_pool', False):
                # SEASON-POOL variant: flat mean over surviving (year, season)
                # pairs within the *shared* years (obs surviving years) —
                # seasons weighted equally, 3/4-season years included, same
                # year core as the headline annual-pool path. Residuals are
                # on season means (matching the headline's mean-of-means
                # convention, so the two MAEs are comparable).
                ovals, mvals = [], []
                for st in sp.SEASONS:
                    o, m = obs_seas[st], m_seas[st]
                    o = o.reindex(surv_years)
                    m = m.reindex(surv_years)
                    ok_s = o.notna() & m.notna()
                    ovals.extend(o[ok_s].values)
                    mvals.extend(m[ok_s].values)
                ovals = np.array(ovals, dtype=float); mvals = np.array(mvals, dtype=float)
                n = len(ovals)
                if n >= 4:
                    dm = mvals - ovals
                    r = dict(buoy=bid, model=model, n_years=n,
                             obs=float(ovals.mean()),
                             model_mean=float(mvals.mean()),
                             bias=float(dm.mean()),
                             bias_pct=float(100.0 * dm.mean() / ovals.mean()),
                             mae=float(np.abs(dm).mean()),
                             rmse=float(np.sqrt((dm ** 2).mean())),
                             flag=flags[model])
                else:
                    r = dict(buoy=bid, model=model, n_years=n,
                             obs=float(ovals.mean()) if n else np.nan,
                             model_mean=float(mvals.mean()) if n else np.nan,
                             bias=np.nan, bias_pct=np.nan, mae=np.nan,
                             rmse=np.nan, flag=flags[model])
            else:
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

    if getattr(args, 'dump_obs_daily', False):
        import pickle as _pk
        with open(os.path.join(BASE, 'pickles/buoy_tas_obs_daily.pkl'), 'wb') as fh:
            _pk.dump({'fingerprint': _obs_cache_fp(args.day_gate, args.gate, args.min_months),
                      'frames': _DUMP}, fh)
        print('dumped obs daily series for %d buoys (gates %s) -> pickles/buoy_tas_obs_daily.pkl'
              % (len(_DUMP), _obs_cache_fp(args.day_gate, args.gate, args.min_months)))
        return
    tab = pd.DataFrame(rows).set_index(['buoy', 'model'])
    survdf = pd.DataFrame(surv).set_index(['buoy', 'model'])
    surv_obs_df = pd.DataFrame(surv_obs).set_index('buoy')

    # Core = obs-driven shared years (same set for all model columns),
    # identical for headline and season-pool runs.
    # Pinned pooled buoy set = the buoys present in the WRF D03 column.
    # Pinning the set to D03 makes the pooled obs baseline identical across all
    # five model columns (previously D03 pooled over fewer buoys, e.g. 46005 is
    # off the D03 domain and only D03 drops it).  Each model column pools over
    # this same set; a column missing a buoy simply skips that buoy.
    core = {b for b in surv_obs_df.index
            if surv_obs_df.loc[b, 'n_obs_years'] >= args.min_years}
    d03_ids = {b for b in core
               if (b, 'WRF_d03') in tab.index and tab.loc[(b, 'WRF_d03'), 'n_years'] >= 1}
    pool_ids = sorted(d03_ids)
    pool_rows = []
    for model in MODELS:
        idx = [(b, model) for b in pool_ids if (b, model) in tab.index and
               tab.loc[(b, model), 'n_years'] >= 1]
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

    print('\nBuoy tas (near-surface air T) bias — 5 model columns — '
          'obs 1986-2005, K — Rule 1 annual (day>=%.0f%%, year>=%.0f%% of days, >=%d months)\n'
          % (100 * args.day_gate, 100 * args.gate, args.min_months))
    print('=== SURVIVING YEARS (headline; Rule 1, day>=%.0f%%, year>=%.0f%%, >=%d months) ==='
          % (100 * args.day_gate, 100 * args.gate, args.min_months))
    o = surv_obs_df.reset_index()
    print(o.to_string(index=False))
    print('\n=== PER BUOY x MODEL ===')
    print(tab.reset_index().to_string(index=False, float_format=lambda x: '%.2f' % x))
    print('\n=== POOLED (pinned buoy set = WRF D03 column: %s) ==='
          % ', '.join(pool_ids))
    print(pool.to_string(float_format=lambda x: '%.2f' % x))


if __name__ == '__main__':
    main()
