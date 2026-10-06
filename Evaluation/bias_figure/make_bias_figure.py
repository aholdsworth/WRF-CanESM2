#!/usr/bin/env python
"""
make_bias_figure.py -- COMPUTE step for the bias figure (D03, 1986-2005).

Layout of the final figure (see bias_figure.ipynb, which is the PLOT step):
  3 columns (T / PR / wind) x 2 rows (bias map on top, model-vs-obs below).
  Land stations only (ECCC + BCH + NOAA; buoys excluded).

Workflow (compute first, then plot):
  1. This script reads all obs + WRF D03 per-station files once,
     computes everything, and saves:
       pickles/bias_figure.pkl          (all arrays needed by the notebook)
       pickles/bias_filter_comparison.csv  (gated vs raw annual bias per station)
  2. bias_figure.ipynb loads the pickle and draws the panels in seconds.
  Delete pickles/bias_figure.pkl to force a recompute.

Conventions (shared with seasonal_pooling.py / taylor_data.py):
  * GATED annuals: obs pooled day -> season (90% per-season completeness)
    -> year (year dropped if any season missing), period value needs
    >= 15 surviving years; model restricted to the same surviving years.
  * RAW annuals: simple mean/sum over all 20 years (taylor_annuals.pkl),
    no season gate.
  * pr bias in %, wind bias in %, t bias in degC.
"""
import os
import sys
import glob
import time
import datetime
import pickle
import warnings

import numpy as np
import pandas as pd
from multiprocessing import Pool
from netCDF4 import Dataset

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, 'data')
FILES = os.path.join(BASE, 'files')
PFM_WRF = os.path.join(DATA, 'wrf_canon')

sys.path.insert(0, os.path.join(BASE, '..', 'elevation_bias'))  # shared seasonal_pooling.py
import seasonal_pooling as sp
import taylor_data as td

START_YEAR, END_YEAR = sp.START_YEAR, sp.END_YEAR
YEARS = sp.YEARS
SEASONS = sp.SEASONS
N_YEARS = len(YEARS)
VARS = ('t', 'pr', 'wind')

OUT_PKL = os.path.join(BASE, 'pickles', 'bias_figure.pkl')
OUT_CSV = os.path.join(BASE, 'pickles', 'bias_filter_comparison.csv')
WIND_CACHE = os.path.join(BASE, 'pickles', '_wind_obs_cache.pkl')
CACHE_VERSION = 5   # 5 = 2026-09-29: canonical station_tables.pkl is now
                   # the single Rule-1 root producer (obs via sp.obs_pool,
                   # skill via sp.skill); land station sets per var = D03
                   # finite bias column: t 146 / pr 159 / wind 15 (sparse
                   # stations dropped); skill column projected from the
                   # table by make_skill_tables.py. Buoy obs+model only from
                   # the validated buoy pickles.

# ---------------------------------------------------------------------------
# Buoy station sets (D03, pinned to the WRF D03 column) - SAME core sets as
# the published buoy_validation summary (buoy_summary_table.tex, 2026-09-24):
#   wind: 7-buoy D03 core (>=8 shared surviving years in the D03 column) --
#         the same core as files/buoy_wind_bias_pooled.csv / the report
#         (2026-09-29: widened from 5 to 7 to match the published pooled CSV)
#   tas : 8 buoys (>=6 shared years; the tas core is these 8, D03-pinned)
# Their obs+model period means come straight from the validated
# buoy_{wind,tas}_tables.pkl 'WRF_d03' rows (no re-pooling here), so the
# bias figure's buoy points are exactly the published buoy-validation set.
# ---------------------------------------------------------------------------
WIND_BUOYS_D03 = ['46041', '46131', '46132', '46146',
                  '46204', '46206', '46207']
TAS_BUOYS_D03  = ['46041', '46050', '46131', '46132', '46146',
                  '46204', '46206', '46207']


def _buoy_truth_coords():
    """(lon, lat) per buoy id from the validated buoy tables (14 buoys)."""
    try:
        w = pickle.load(open(os.path.join(BASE, 'pickles',
                                          'buoy_wind_tables.pkl'), 'rb'))
        return w['truth']
    except Exception:
        t = pickle.load(open(os.path.join(BASE, 'pickles',
                                          'buoy_tas_tables.pkl'), 'rb'))
        return t['truth']


def _buoy_d03_row(table_pkl, buoy_id):
    """The 'WRF_d03' row of a buoy tables pickle for one buoy id, or None."""
    try:
        tab = pickle.load(open(os.path.join(BASE, 'pickles', table_pkl),
                               'rb'))['tab']
        if (buoy_id, 'WRF_d03') in tab.index:
            return tab.loc[(buoy_id, 'WRF_d03')]
    except Exception:
        pass
    return None

# Populated in main() before the wind Pool starts (forked workers inherit it).
_WIND_CACHE = None

# Pre-computed D03 buoy annuals (obs & model, from the validated buoy tables);
# keyed by (sid, var), values are 20-yr arrays aligned to YEARS.  Filled in
# main() before the Pool starts so forked workers inherit it (shared memory).
_BUOY_ANNUAL = {}
# Pre-computed D03 buoy SEASONAL means (obs & model), keyed (sid, var),
# values (4,) in SEASONS order, NaN where the season has no valid years.
# Only populated for variables whose validated table stores seasonal frames
# (wind does; tas does not, so its seasonal scatter points stay off).
_BUOY_SEAS = {}


def _load_buoy_annuals():
    """Read the validated buoy tables' WRF_d03 rows into 20-yr annual arrays.

    wind -> m/s (as stored) ; tas -> K -> C.  A year is kept only if the
    per-year model value is present in the D03 column; the obs year value is
    carried through (the published table's shared-year pairing).
    """
    out = {}
    # wind: per_buoy lives inside buoy_wind_tables.pkl (key '<sid>|WRF_d03')
    try:
        pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                           'buoy_wind_tables.pkl'), 'rb'))
        per = pk.get('per_buoy', {})
        for sid in WIND_BUOYS_D03:
            d = per.get('%s|WRF_d03' % sid)
            if not d:
                continue
            o = np.asarray(d['annual_obs'], float).reshape(-1)
            m = np.asarray(d['annual_model'], float).reshape(-1)
            if len(o) == N_YEARS and len(m) == N_YEARS:
                out[(sid, 'wind')] = (o, m)
            # seasonal frames (20 x 4) already Rule-2 gated in buoy_wind_bias
            so = np.asarray(d.get('seas_obs'), float).reshape(-1, 4)
            sm = np.asarray(d.get('seas_model'), float).reshape(-1, 4)
            if so.shape == (N_YEARS, 4) and sm.shape == (N_YEARS, 4):
                _BUOY_SEAS[(sid, 'wind')] = (
                    np.array([np.nanmean(so[:, j]) for j in range(4)]),
                    np.array([np.nanmean(sm[:, j]) for j in range(4)]))
    except Exception:
        pass
    # tas: per-buoy D03 annuals in buoy_tas_d03_per_buoy.pkl (K -> C)
    try:
        pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                           'buoy_tas_d03_per_buoy.pkl'), 'rb'))
        for sid in TAS_BUOYS_D03:
            d = pk.get(sid)
            if not d:
                continue
            o = np.asarray(d['annual_obs'], float).reshape(-1) - 273.15
            m = np.asarray(d['annual_model'], float).reshape(-1) - 273.15
            if len(o) == N_YEARS and len(m) == N_YEARS:
                out[(sid, 't')] = (o, m)
            # seasonal frames (4, in SEASONS order), K -> C
            so = np.asarray(d.get('seas_obs'), float).reshape(-1)
            sm = np.asarray(d.get('seas_model'), float).reshape(-1)
            if so.shape == (len(SEASONS),) and sm.shape == (len(SEASONS),):
                _BUOY_SEAS[(sid, 't')] = (so - 273.15, sm - 273.15)
    except Exception:
        pass
    return out


# ---------------------------------------------------------------------------
# ECCC hourly-wind obs cache (shared across pool workers via shared memory).
# Each station is 20 years of hourly CSV (~130-200 s to read from GPFS);
# reading it once in the main process and sharing the dict through the forked
# Pool's shared memory cuts ~35 min of wall time to ~4 min.
# ---------------------------------------------------------------------------
def _build_wind_cache(eche_ids):
    if os.path.exists(WIND_CACHE):
        try:
            with open(WIND_CACHE, 'rb') as f:
                c = pickle.load(f)
            if len(c) == len(eche_ids):
                print('wind obs cache: loaded %d stations from %s'
                      % (len(c), os.path.basename(WIND_CACHE)), flush=True)
                return c
        except Exception:
            pass
    cache = {}
    t_start = time.time()
    for i, sid in enumerate(eche_ids, 1):
        t0 = time.time()
        cache[sid] = sp.eccc_hourly_wind_series(sid, DATA)
        dt = time.time() - t0
        print('   wind cache %2d/%d  %s: %.0f s  (total %.0f s)' %
              (i, len(eche_ids), sid, dt, time.time() - t_start), flush=True)
    with open(WIND_CACHE, 'wb') as f:
        pickle.dump(cache, f)
    return cache


# ---------------------------------------------------------------------------
# Station metadata (land stations only, no buoys)
# ---------------------------------------------------------------------------
def land_coords(meta=None):
    """lat/lon/elev for the land stations, from the CANONICAL table
    pickles/station_tables.pkl (the same 173 land rows that define the
    per-variable station sets -- single source of truth).  (Buoys get
    their own figure.)  meta is unused; the signature is kept for
    call-site compatibility."""
    d = td._table_d03()
    land = d[['station', 'lat', 'lon', 'obs_elev_m']].copy()
    land['station'] = land['station'].astype(str)
    land = land.set_index('station')
    lat = land['lat'].astype(float)
    lon = land['lon'].astype(float)
    el  = land['obs_elev_m'].astype(float)
    lat.index = lat.index.astype(str)
    lon.index = lon.index.astype(str)
    el.index = el.index.astype(str)

    # D03 buoys (truth coords, sea level).  Union of the wind & tas buoy
    # sets; a buoy appears in the figure only for the variable(s) it has.
    truth = _buoy_truth_coords()
    buoy_ids = sorted(set(WIND_BUOYS_D03) | set(TAS_BUOYS_D03))
    b_lat = pd.Series({b: truth[b][1] for b in buoy_ids if b in truth}, dtype=float)
    b_lon = pd.Series({b: truth[b][0] for b in buoy_ids if b in truth}, dtype=float)
    b_el  = pd.Series({b: 0.0 for b in b_lat.index}, dtype=float)
    lat = pd.concat([lat, b_lat])
    lon = pd.concat([lon, b_lon])
    el  = pd.concat([el, b_el])

    lat = lat[~lat.index.duplicated()]
    lon = lon[~lon.index.duplicated()]
    el  = el[~el.index.duplicated()]
    return lat.sort_index(), lon.sort_index(), el.sort_index()


# ---------------------------------------------------------------------------
# Gated per-station annual model & obs (season-first pool, >= MIN_YEARS)
# ---------------------------------------------------------------------------
def _daily_series(sid, var):
    """Obs daily series for one station via taylor_data agency convention.
    Returns (series, freq) or (None, None)."""
    agency = td._agency_for(sid)
    if agency == 'ECCC':
        if var == 'wind':
            if _WIND_CACHE is not None and sid in _WIND_CACHE:
                s = _WIND_CACHE[sid]
            else:
                s = sp.eccc_hourly_wind_series(sid, DATA)
            return (s, 'H') if s is not None else (None, None)
        s = sp.eccc_daily_series(sid, var, DATA)
        return (s, 'D') if s is not None else (None, None)
    if agency == 'NOAA':
        s = sp.noaa_daily_series(sid, var, DATA)
        return (s, 'D') if s is not None else (None, None)
    s = sp.bch_daily_series(sid, var, DATA)
    return (s, 'D') if s is not None else (None, None)


def _station_worker(args):
    """One station/variable: gated annual obs & model + seasonal model/obs.

    Returns dict with keys: obs_ann, mod_ann, seas_mod (4,), seas_obs (4,),
    seas_bias (4,)  -- all NaN where the gates fail.
    """
    global _WIND_CACHE
    sid, var = args
    kind = 'sum' if var == 'pr' else 'mean'
    out = {'obs_ann': np.nan, 'mod_ann': np.nan,
           'seas_mod': np.full(len(SEASONS), np.nan),
           'seas_obs': np.full(len(SEASONS), np.nan),
           'seas_bias': np.full(len(SEASONS), np.nan)}
    # --- D03 buoy fast path: consume the validated buoy-table annuals ---
    # (buoys can't go through the daily season-gate; their shared-year
    #  pairing + metrics are already locked in buoy_validation).  Annuals
    #  for all buoy vars; seasonal frames where the validated table stores
    #  them (wind).  The seasonal scatter panel picks these up.
    if (sid, var) in _BUOY_ANNUAL:
        o_arr, m_arr = _BUOY_ANNUAL[(sid, var)]
        yrs = [y for y in YEARS if np.isfinite(m_arr[y - START_YEAR])]
        if len(yrs) >= 2:      # buoy minimum (handoff: n>=2 for metrics)
            ov = o_arr[[y - START_YEAR for y in yrs]]
            mv = m_arr[[y - START_YEAR for y in yrs]]
            out['obs_ann'] = float(np.nanmean(ov))
            out['mod_ann'] = float(np.nanmean(mv))
        if (sid, var) in _BUOY_SEAS:
            out['seas_obs'], out['seas_mod'] = _BUOY_SEAS[(sid, var)]
        return out
    try:
        series, freq = _daily_series(sid, var)
        if series is None or len(series.dropna()) == 0:
            return out
        pool = sp.obs_seasonal_pool(series, kind, freq=freq)
        if pool['n_years'] < sp.MIN_YEARS:
            return out
        p = td._wrf_file(var, 'd03', sid)
        if p is None:
            return out
        try:
            m_arr = td._wrf_annual_series(p, var)   # 20-annual array
        except Exception:
            return out
        m_ann = pd.Series(m_arr, index=YEARS)
        yrs = [y for y in pool['years'] if np.isfinite(m_ann.get(y, np.nan))]
        if len(yrs) < sp.MIN_YEARS:
            return out
        out['obs_ann'] = float(pool['period'])
        out['mod_ann'] = float(m_ann.loc[yrs].mean())
        try:
            m_daily = _wrf_daily(p, var)
            m_seas = sp.seasonal_year_frame(m_daily, kind, freq='D')
        except Exception:
            m_seas = None
        if m_seas is None:
            return out
        # Seasonal values use the Rule-2 3-months-present gate; the year must
        # also be a valid (Rule-1) annual year and have all 4 model seasons.
        yrs2 = [y for y in yrs
                if m_seas.loc[y].notna().all() and pool['year_mask'].loc[y]]
        for i, seas_name in enumerate(SEASONS):
            if len(yrs2) < sp.MIN_YEARS:
                break
            ov = pool['seas'].loc[yrs2, seas_name].values
            mv = m_seas.loc[yrs2, seas_name].values
            out['seas_obs'][i] = np.nanmean(ov)
            out['seas_mod'][i] = np.nanmean(mv)
            out['seas_bias'][i] = (np.nanmean(100.0 * (mv - ov) / ov)
                                   if var == 'pr' else np.nanmean(mv - ov))
    except Exception:
        pass
    return out


def _compute_var(var, meta, land_idx, n_workers=16):
    """Parallel per-station compute for one variable."""
    sids = [s for s in td.station_set(var, meta) if s in land_idx]
    tasks = [(s, var) for s in sids]
    with Pool(n_workers) as pool:
        res = pool.map(_station_worker, tasks, chunksize=4)
    keep = [True] * len(tasks)
    seen = {}
    for i, (s, _) in enumerate(tasks):       # keep FIRST entry per station id
        if s in seen:
            keep[i] = False
        else:
            seen[s] = i
    tasks, res = ([t for t, k in zip(tasks, keep) if k],
                  [r for r, k in zip(res, keep) if k])
    idx = pd.Index([s for s, _ in tasks])
    o_g = pd.Series([r['obs_ann'] for r in res], index=idx, dtype=float)
    m_g = pd.Series([r['mod_ann'] for r in res], index=idx, dtype=float)
    sm = pd.DataFrame([r['seas_mod'] for r in res], index=idx, columns=SEASONS)
    so = pd.DataFrame([r['seas_obs'] for r in res], index=idx, columns=SEASONS)
    sb = pd.DataFrame([r['seas_bias'] for r in res], index=idx, columns=SEASONS)
    return o_g, m_g, sm, so, sb


def _wrf_daily(path, kind):
    """WRF D03 per-station file -> daily Series (t: C, pr: mm/d, wind: m/s)."""
    nc = Dataset(path, mode='r')
    t = np.squeeze(nc.variables['time'][:]).astype(float)
    vname = {'t': 'T2', 'pr': 'pr', 'wind': 'wspd'}[kind]
    v = np.squeeze(nc.variables[vname][:]).astype(float)
    nc.close()
    dates = pd.to_datetime(datetime.datetime(START_YEAR, 1, 1)) \
              + pd.to_timedelta(t, unit='h')
    s = pd.Series(v, index=dates)
    if kind == 't':
        s = (s - 273.15).round(1)
    elif kind == 'pr':
        inc = np.clip(np.diff(s.values), 0.0, None)
        s = pd.Series(inc, index=dates[1:])
    h = s.groupby(s.index.normalize())
    return h.sum() if kind == 'pr' else h.mean()


# ---------------------------------------------------------------------------
def bias_of(obs, mod, var):
    return 100.0 * (mod - obs) / obs if var == 'pr' else mod - obs


def main():
    os.makedirs(os.path.dirname(OUT_PKL), exist_ok=True)

    meta = td.load_meta()
    LAT, LON, ELEV = land_coords(meta)
    land_idx = set(LAT.index)
    print(f'stations: {len(LAT)} (ECCC/BCH/NOAA land + D03 buoys)')

    # Load the validated D03 buoy annuals once (forked workers inherit it).
    global _BUOY_ANNUAL
    _BUOY_ANNUAL = _load_buoy_annuals()
    print(f'D03 buoy annuals loaded: {sorted(set(k[0] for k in _BUOY_ANNUAL))}')

    taylor = {var: td.build(var) for var in VARS}
    raw = {}   # var -> (obs_period, mod_period): scalar per station
    raw_ann = {}   # var -> (obs: Series-of-arrays, mod: Series-of-arrays)
    for var in VARS:
        o, models = taylor[var]
        m = models['D03']
        o_s = pd.Series({str(k): np.asarray(v, float) for k, v in o.items()})
        m_s = pd.Series({str(k): np.asarray(v, float) for k, v in m.items()})

        # The land Taylor series carries no buoys, so the map bias for the
        # tas buoys would otherwise be NaN and they'd be dropped from the
        # temperature map (the wind Taylor series happens to include the D03
        # wind buoys, so only tas was affected).  Add the D03 buoy annual
        # arrays to o_s/m_s so the raw/raw_ann slices below include them
        # (obs & model period means over the shared years, from _BUOY_ANNUAL).
        for (sid, bvar), (o_arr, m_arr) in _BUOY_ANNUAL.items():
            if bvar == var:
                o_s[str(sid)] = np.asarray(o_arr, float)
                m_s[str(sid)] = np.asarray(m_arr, float)

        common = o_s.index.intersection(m_s.index).intersection(land_idx)
        buoy_in = [str(sid) for (sid, bvar) in _BUOY_ANNUAL if bvar == var
                   and str(sid) in o_s.index]
        keep = list(common) + [b for b in buoy_in if b not in set(common)]
        raw_o = o_s.loc[keep].apply(np.nanmean)
        raw_m = m_s.loc[keep].apply(np.nanmean)
        raw_ann_o = o_s.loc[keep]
        raw_ann_m = m_s.loc[keep]
        raw[var] = (raw_o, raw_m)
        raw_ann[var] = (raw_ann_o, raw_ann_m)
        print(f'{var}: raw (taylor) annuals for {len(common)} land + '
              f'{len(buoy_in)} D03 buoys')

    out = {
        '_version': CACHE_VERSION,
        '_generated': datetime.datetime.utcnow().isoformat(timespec='seconds'),
        'stations': {
            'lat': LAT, 'lon': LON, 'elev': ELEV,
            'buoy_ids': sorted(set(WIND_BUOYS_D03) | set(TAS_BUOYS_D03)),
            'per_var': {var: list(td.station_set(var, meta)) for var in VARS},
        },
        'raw':      {var: (raw[var][0], raw[var][1]) for var in VARS},
        'raw_ann':  {var: (raw_ann[var][0], raw_ann[var][1]) for var in VARS},
        'gated':    {},
        'seasonal': {},
        'skill':    {},
    }

    cmp_frames = []
    for var in VARS:
        t0 = time.time()
        print(f'== {var}: gated annuals + seasonal values ==', flush=True)
        if var == 'wind':
            global _WIND_CACHE
            # ECCC-hourly obs cache over the canonical land wind set (table
            # rows; buoys come from the validated tables, not this cache).
            _WIND_CACHE = _build_wind_cache(td.station_ids_from_table('wind'))
        o_g, m_g, sm, so, sb = _compute_var(var, meta, land_idx)
        out['gated'][var] = (o_g, m_g)
        out['seasonal'][var] = {'mod': sm, 'obs': so, 'bias': sb}
        print(f'   ({time.time() - t0:.1f} s)', flush=True)

        # skill_seasonal.pkl annual bias (Eva's published numbers)
        skill = pickle.load(open(os.path.join(BASE, 'pickles',
                                              'skill_seasonal.pkl'), 'rb'))['D03']
        out['skill'][var] = skill[f'{var}_bias']

        # ---- two-filter comparison (gated vs raw annual bias) ----
        b_gated = bias_of(o_g, m_g, var)
        o_raw, m_raw = raw[var]
        b_raw = bias_of(o_raw.reindex(o_g.index),
                        m_raw.reindex(o_g.index), var)
        cmp_frames.append(pd.DataFrame({
            'var': var, 'station': o_g.index,
            'bias_gated': b_gated.values, 'bias_raw': b_raw.values,
            'obs_gated': o_g.values, 'mod_gated': m_g.values,
            'obs_raw': o_raw.reindex(o_g.index).values,
            'mod_raw': m_raw.reindex(o_g.index).values,
        }))
        both = pd.DataFrame({'var': var, 'station': o_g.index,
                             'bias_gated': b_gated.values,
                             'bias_raw': b_raw.values}).dropna()
        d = (both['bias_gated'] - both['bias_raw']).to_numpy()
        print(f'   gated: {int(b_gated.notna().sum())} stations | '
              f'matched n={len(both)}, mean|diff|={float(np.abs(d).mean()):.3f}, '
              f'max|diff|={float(np.abs(d).max()):.3f}')
        print(f'   seasonal values at '
              f'{int(sb.notna().any(axis=1).sum())} stations')

    pd.concat(cmp_frames, ignore_index=True).to_csv(OUT_CSV, index=False,
                                                      float_format='%.4f')
    with open(OUT_PKL, 'wb') as f:
        pickle.dump(out, f)
    print(f'\nsaved {OUT_PKL}')
    print(f'saved {OUT_CSV}')


if __name__ == '__main__':
    main()
