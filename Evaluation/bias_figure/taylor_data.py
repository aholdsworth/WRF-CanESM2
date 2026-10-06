#!/usr/bin/env python
"""
taylor_data.py -- per-variable, PER-STATION ANNUAL series for the publication
Taylor figure (T / PR / wind, 1986-2005).

CONVENTION (re-anchored 2026-09-26; see taylor_3panel_handoff.md):
  * fundamental unit is the PER-STATION ANNUAL value, 1986-2005:
        t    : annual MEAN (degC)     -> 20 values / station
        pr   : annual SUM   (mm)      -> 20 values / station
        wind : annual MEAN (m/s)      -> 20 values / station
  * the figure (taylor_figure.py) pools (station x year) samples over the
    years where BOTH obs and model are valid (Eva's convention; NO anomaly
    detrending).  Model years are complete, so this is exactly like-for-like.
  * obs year validity comes from the SAME loaders/gates as the skill & bias
    pipelines (seasonal_pooling.py is the single source of truth):
        land daily (ECCC t/pr, NOAA t/pr, BCH t/pr): sp.annual_frame
            Rule-1 90%-of-days gate
        land wind (ECCC hourly, NOAA AWND): RELAXED wind-obs convention
            (2026-09-29): annual mean over whatever hours/days exist that
            year (no 90% completeness gate).  Reason: ECCC 192/277 have
            only ~50-67% of hours each year (no complete year exists), so
            the Rule-1 gate would drop them entirely, while the canonical
            station_tables.pkl (ungated obs mean) carries them (+30.2 /+
            +14.2%).  Wind obs only -- t/pr still use the strict 90% gate.
            (freq='H' for ECCC hourly files; NOAA AWND is daily)
        buoys: the validated relaxed-gate (75/80/10m) pooled arrays from
            buoy_wind_bias.py / buoy_tas_bias.py (table pickles, below)
  * station sets (land) come from the CANONICAL land table
    pickles/station_tables.pkl (2026-09-29 USW-fixed, 173 land rows): the
    per-variable set is the D03 rows with a finite bias column, so
        t    : 156 land + 8 D03 tas buoys (164)
        pr   : 173 land
        wind : 21 land + 7 D03 wind buoys (28)
    Buoys are scored separately (below) and are NOT part of the table's
    land rows; the buoy obs/model arrays come only from the validated
    buoy pickles (see "BUOY MODEL SERIES").

This module returns, per variable:
    obs    : dict  station -> 1D array (20 annual values; NaN = invalid year)
    models : dict  label -> dict  station -> 1D array (20 annual values)

Data sources (all local, under this workspace):
  obs (land)  data/eccc_obs, data/eccc_hourly_obs, data/noaa_obs, data/bch_obs
              via seasonal_pooling series loaders (shared with skill/bias)
  obs (buoy)  pickles/buoy_wind_tables.pkl (per_buoy '<ID>|WRF_d03' annual_obs)
              pickles/buoy_tas_d03_per_buoy.pkl (annual_obs, K -> C)
  model WRF   data/wrf_canon (per-station hourly files; make_wrf_canon.py)
  model Can*  data/cdo_extractions/{canesm2_raw_noaa113, canrcm4_noaa113}
              (land) and data/cdo_extractions/{canesm2_raw_buoys,
              canrcm4_buoys} (buoys; validated by the buoy scripts)

BUOY MODEL SERIES (important): taken from the VALIDATED arrays -- the exact
ones used by the published buoy tables and the bias figure -- NOT re-parsed
from data/wrf_canon (whose 46xxx buoy files are CORRUPT D01-era artifacts:
d02 wind = 1-yr duplicated time axis; d03 wind = rotated embedded coords AND
values; d02/d03 t/pr same artifact class):
  * wind: buoy_wind_tables.pkl per_buoy['<ID>|WRF_d01|d02|d03' / 'CanESM2' /
          'CanRCM4']['annual_model'] (buoy_wind_bias.py, relaxed gate)
  * t   : D03 column from buoy_tas_d03_per_buoy.pkl (validated D03-pinned set;
          46050 has no D03 WRF file -> no D03 t point); CanESM2/CanRCM4 from
          the validated cdo_extractions/*_buoys daily files (buoy_tas_bias.py
          sources) with the same annual-mean pooling.
  The validated buoy WRF scripts used the good
  data/wrf_stations/wind_d{02,03}_st<ID>.nc re-extractions (verified 2026-09-26
  to reproduce the table arrays to 0.0000 m/s).

Run self-test:  python taylor_data.py
"""
import os, glob, datetime, warnings
import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

_BASE = os.path.dirname(os.path.abspath(__file__))
import sys as _sys
_sys.path.insert(0, os.path.join(_BASE, '..', 'elevation_bias'))  # shared seasonal_pooling.py
import seasonal_pooling as sp

BASE = _BASE
DATA = os.path.join(BASE, 'data')
FILES = os.path.join(BASE, 'files')
PFM_WRF = os.path.join(DATA, 'wrf_canon')
# Validated per-buoy Can-model extraction dirs (per-buoy daily 1986-2005 files,
# validated by the buoy scripts).  NOTE the cdo_extractions/ prefix -- the
# earlier raw path omitted it, so Can* buoy lookups silently returned empty.
CDO_CAN_BUOYS = {'canesm2_raw_buoys': os.path.join(DATA, 'cdo_extractions',
                                                   'canesm2_raw_buoys'),
                 'canrcm4_buoys': os.path.join(DATA, 'cdo_extractions',
                                               'canrcm4_buoys')}
# Canonical WRF per-station dir (built by make_wrf_canon.py).
# ONE naming convention: {var}_{dom}_{AGENCY}_{ID}.nc  (AGENCY = ECCC|NOAA|BCH)
# READ LAND STATIONS ONLY -- the 46xxx buoy files are corrupt D01-era
# artifacts (see module docstring); all buoy series come from the validated
# table pickles / cdo_extractions/*_buoys.

START_YEAR, END_YEAR = 1986, 2005
YEARS = list(range(START_YEAR, END_YEAR + 1))
N_YEARS = len(YEARS)

# Disk cache of the full annual tables. Delete to force a recompute.
CACHE_PKL = os.path.join(BASE, 'pickles', 'taylor_annuals.pkl')
CACHE_VERSION = 10         # 10 = 2026-09-30: add filtered D03 buoy points to
                           # CanESM2/CanRCM4 (wind 7 from the validated
                           # buoy_wind_tables.pkl, t 8 from the validated
                           # cdo_extractions/*_buoys tas files) so all 5 model
                           # columns share the same obs set.  9 = 2026-09-30:
                           # canon-fill -- station_tables.pkl
                           # rebuilt with the promoted D01/1001 model files so
                           # all 5 datasets are equal (t 146 / pr 159 /
                           # wind 16).  Raw rebuild re-derives membership from
                           # the table (station_ids_from_table) and now reads
                           # the 30 promoted wrf_canon files (18 BCH + 82 +
                           # 202 for D01 t/pr; wind 1001 d01/d02/d03; 5
                           # USW d01 wind).
                           # 8 = 2026-09-29: UNIFY -- re-indexed from the
                           # v7 cache onto the canonical Rule-1 land set
                           # (station_tables.pkl D03 finite bias column:
                           # t 146 / pr 159 / wind 15); sparse land stations
                           # dropped from every obs + model dict, arrays kept
                           # byte-identical.  A raw rebuild from this version
                           # reproduces the same membership automatically
                           # (station_ids_from_table reads the table).
                           # 7 = 2026-09-29: land station sets from the
                           # canonical station_tables.pkl (USW-fixed, 173
                           # land rows; per-var = D03 finite bias column:
                           # t 156 / pr 173 / wind 21); buoy ids excluded
                           # from ALL _wrf_annual calls (no buoy model value
                           # ever parsed from data/wrf_canon); t D01/D02 buoy
                           # model arrays added from the validated
                           # buoy_tas_d03_per_buoy.pkl (annual_model_d01/
                           # d02).  Buoys: 8 tas / 7 wind (validated pickles).
# Canonical land table: single source of truth for the land station set,
# coords and elevations shared by make_bias_figure.py.
TABLE_PKL = os.path.join(BASE, 'pickles', 'station_tables.pkl')
TABLE_VAR_COL = {'t': 'temperature_bias_C', 'pr': 'precipitation_bias_pct',
                 'wind': 'wind_bias_pct'}

_TABLE_CACHE = None
def _table_d03():
    """The D03 frame of pickles/station_tables.pkl (canonical 173 land rows)."""
    global _TABLE_CACHE
    if _TABLE_CACHE is None:
        import pickle
        _TABLE_CACHE = pickle.load(open(TABLE_PKL, 'rb'))['D03']
    return _TABLE_CACHE


def station_ids_from_table(var):
    """Land station ids (as str) for `var`: the D03 table rows with a finite
    bias column.  This is the canonical land set -- single source of truth
    shared by taylor_data and make_bias_figure."""
    d = _table_d03()
    return [str(s) for s in d.loc[d[TABLE_VAR_COL[var]].notna(), 'station']]


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
def load_meta():
    e = pd.read_csv(os.path.join(FILES, 'ECCC_d03_stations.csv'), header=None)
    eccc_ids = [str(i) for i in e.iloc[:, 4]]
    eche = pd.read_csv(os.path.join(FILES, 'ECCC_d03_stations_hourly.csv'), header=None)
    eche_ids = [str(i) for i in eche.iloc[:, 4]]
    b = pd.read_csv(os.path.join(FILES, 'BCH_d03_stations.csv'))
    bch_ids = [str(i) for i in b['STATION_NO']]

    def noaa(fn):
        d = pd.read_csv(os.path.join(FILES, fn))
        return [str(i) for i in d['STATION']]
    n_tas = noaa('NOAA_d03_stations_tas.csv')
    n_pr = noaa('NOAA_d03_stations.csv')
    n_w = noaa('NOAA_d03_stations_wind.csv')
    return dict(eccc_ids=eccc_ids, eche_ids=eche_ids, bch_ids=bch_ids,
                n_tas_ids=n_tas, n_pr_ids=n_pr, n_w_ids=n_w)


def station_set(var, meta=None):
    """Full station set for `var`: canonical land set from
    station_tables.pkl (see station_ids_from_table) + the D03 buoy set
    (obs/model from the validated buoy pickles only)."""
    import make_bias_figure as mbf
    land = station_ids_from_table(var)
    if var == 't':
        return land + mbf.TAS_BUOYS_D03
    if var == 'wind':
        return land + mbf.WIND_BUOYS_D03
    return land


# ---------------------------------------------------------------------------
# Observations -> per-station ANNUAL series
# Land: seasonal_pooling loaders + Rule-1 gate (identical to skill/bias).
# Buoys: validated pooled annual arrays (relaxed gate) from the table pickles.
# ---------------------------------------------------------------------------
def _annual_array(s, kind, freq):
    """sp.annual_frame (Rule 1, 90% gate) -> (20,) array aligned to YEARS,
    or None if no year survives.  freq='D' daily, 'H' hourly."""
    if s is None or len(s) == 0:
        return None
    ann = sp.annual_frame(s, kind, freq=freq)
    arr = np.full(N_YEARS, np.nan)
    for j, y in enumerate(YEARS):
        v = ann.loc[y]
        if pd.notna(v):
            arr[j] = float(v)
    return None if np.all(np.isnan(arr)) else arr


def _land_daily_annual(ids, agency_fn, var):
    """ECCC/NOAA/BCH DAILY obs -> per-station annual series via
    seasonal_pooling loaders + Rule-1 gate.  agency_fn(sid) -> daily Series."""
    kind = 'mean' if var == 't' else 'sum'
    out = {}
    for sid in ids:
        arr = _annual_array(agency_fn(str(sid)), kind, 'D')
        if arr is not None:
            out[sid] = arr
    return out


def _land_wind_hourly_annual(ids):
    """Wind obs -> per-station annual MEAN (m/s), RELAXED convention
    (2026-09-29): each year is the mean over the hours/days that EXIST
    (>=1 valid sample), with NO 90% completeness gate.  ECCC hourly CSVs
    via sp.eccc_hourly_wind_series (already m/s); NOAA AWND via
    sp.noaa_daily_series('wind') (m/s).  This matches the canonical
    station_tables.pkl wind obs (ungated mean), so the table's finite wind
    set (192/277, ~50-67% data completeness) can be plotted.  t/pr obs are
    unaffected and still use the strict Rule-1 90% gate."""
    out = {}
    for sid in ids:
        s = sp.eccc_hourly_wind_series(str(sid), DATA)
        if s is None or len(s) == 0:
            s = sp.noaa_daily_series(str(sid), 'wind', DATA)
        if s is None or len(s) == 0:
            continue
        s = s[(s.index >= pd.Timestamp(START_YEAR, 1, 1))
              & (s.index <  pd.Timestamp(END_YEAR + 1, 1, 1))]
        ann = s.groupby(s.index.year).mean()
        arr = np.full(N_YEARS, np.nan)
        for y, v in ann.items():
            if START_YEAR <= int(y) <= END_YEAR and pd.notna(v):
                arr[int(y) - START_YEAR] = float(v)
        if not np.all(np.isnan(arr)):
            out[str(sid)] = arr
    return out


# -- buoy obs: validated pooled annual arrays (relaxed gate, buoy scripts) --
_BUOY_WIND_OBS = None
def _buoy_wind_obs():
    global _BUOY_WIND_OBS
    if _BUOY_WIND_OBS is None:
        import pickle
        pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                           'buoy_wind_tables.pkl'), 'rb'))
        per = pk.get('per_buoy', {})
        out = {}
        for bid, d in per.items():
            if not bid.endswith('|WRF_d03'):
                continue
            o = np.asarray(d.get('annual_obs'), float).reshape(-1)
            if len(o) == N_YEARS and not np.all(np.isnan(o)):
                out[bid.split('|')[0]] = o
        _BUOY_WIND_OBS = out
    return _BUOY_WIND_OBS


_BUOY_TAS_OBS = None
def _buoy_tas_obs():
    global _BUOY_TAS_OBS
    if _BUOY_TAS_OBS is None:
        import pickle
        pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                           'buoy_tas_d03_per_buoy.pkl'), 'rb'))
        out = {}
        for bid, d in pk.items():
            o = np.asarray(d.get('annual_obs'), float).reshape(-1)
            if len(o) == N_YEARS and not np.all(np.isnan(o)):
                out[str(bid)] = o - 273.15      # K -> C
        _BUOY_TAS_OBS = out
    return _BUOY_TAS_OBS


def obs_annual(var, meta):
    """Per-station annual obs series for the variable's station set."""
    import make_bias_figure as mbf
    if var == 't':
        o = {}
        o.update(_land_daily_annual(meta['eccc_ids'],
                                    lambda s: sp.eccc_daily_series(s, 't', DATA), 't'))
        o.update(_land_daily_annual(meta['bch_ids'],
                                    lambda s: sp.bch_daily_series(s, 't', DATA), 't'))
        o.update(_land_daily_annual(meta['n_tas_ids'],
                                    lambda s: sp.noaa_daily_series(s, 't', DATA), 't'))
        ob = _buoy_tas_obs()
        for sid in mbf.TAS_BUOYS_D03:
            if sid in ob:
                o[sid] = ob[sid]
    elif var == 'pr':
        o = {}
        o.update(_land_daily_annual(meta['eccc_ids'],
                                    lambda s: sp.eccc_daily_series(s, 'pr', DATA), 'pr'))
        o.update(_land_daily_annual(meta['bch_ids'],
                                    lambda s: sp.bch_daily_series(s, 'pr', DATA), 'pr'))
        o.update(_land_daily_annual(meta['n_pr_ids'],
                                    lambda s: sp.noaa_daily_series(s, 'pr', DATA), 'pr'))
    elif var == 'wind':
        o = {}
        # canonical land wind set from station_tables.pkl (2026-09-29)
        o.update(_land_wind_hourly_annual(station_ids_from_table('wind')))
        ob = _buoy_wind_obs()
        for sid in mbf.WIND_BUOYS_D03:
            if sid in ob:
                o[sid] = ob[sid]
    return o


# ---------------------------------------------------------------------------
# Model data -> per-station ANNUAL series
# ---------------------------------------------------------------------------
def _agency_for(sid):
    # must match make_wrf_canon.canonical_name() exactly
    if sid.startswith('46'):
        return 'NOAA' if sid in ('46005', '46029', '46041', '46050',
                                 '46087', '46088', '46089') else 'ECCC'
    if sid in ('desw1', 'wpow1'):
        return 'NOAA'
    if sid.startswith('USC') or sid.startswith('USW'):
        return 'NOAA'
    if sid.isdigit():
        return 'ECCC'
    return 'BCH'


def _wrf_file(var, dom, sid):
    # canonical convention: {var}_{dom}_{AGENCY}_{ID}.nc
    p = os.path.join(PFM_WRF, '%s_%s_%s_%s.nc'
                     % (var, dom, _agency_for(sid), sid))
    return p if os.path.exists(p) else None


def _wrf_annual_series(path, kind):
    """Per-station ANNUAL series from a WRF per-station file.
    t: (T2-273.15).round(1) hourly -> annual mean.
    wind: wspd hourly -> annual mean.
    pr: cumulative field -> clip(diff,0) hourly increments -> drop Jan-1 00:00
         discontinuity -> annual SUM (mm).  (Eva: groupby('Y').sum().)"""
    nc = Dataset(path, mode='r')
    t = np.squeeze(nc.variables['time'][:]).astype(float)
    vname = {'t': 'T2', 'pr': 'pr', 'wind': 'wspd'}[kind]
    v = np.squeeze(nc.variables[vname][:]).astype(float)
    nc.close()
    dates = pd.to_datetime(datetime.datetime(START_YEAR, 1, 1)) + pd.to_timedelta(t, unit='h')
    s = pd.Series(v, index=dates)
    if kind == 't':
        s = (s - 273.15).round(1)
    elif kind == 'wind':
        pass
    else:
        inc = np.clip(np.diff(s.values), 0.0, None)
        s = pd.Series(inc, index=dates[1:])
        s = s[~((s.index.day == 1) & (s.index.month == 1) & (s.index.hour == 0))]
    if kind == 'pr':
        ann = s.groupby(s.index.year).sum()
    else:
        ann = s.groupby(s.index.year).mean()
    arr = np.full(N_YEARS, np.nan)
    for y, val in ann.items():
        if START_YEAR <= y <= END_YEAR:
            arr[y - START_YEAR] = val
    return arr


def _wrf_annual(ids, var, dom):
    """Land WRF annuals from data/wrf_canon ONLY.  Buoy ids (46*) are
    EXCLUDED for every variable: their validated model series are overlaid
    from the table pickles in model_annual(), and the 46xxx wrf_canon files
    are corrupt D01-era artifacts (see module docstring)."""
    out = {}
    for sid in ids:
        if str(sid).startswith('46'):
            continue
        p = _wrf_file(var, dom, sid)
        if p is None:
            continue
        try:
            arr = _wrf_annual_series(p, var)
        except Exception:
            continue
        if not np.all(np.isnan(arr)):
            out[sid] = arr
    return out


def _can_daily_annual(sid, var, subdir):
    """Per-station annual array from a cdo_extractions daily file (one of the
    validated M10/buoy extraction dirs).  Returns {sid: arr} or {}."""
    vname = {'t': 'tas', 'pr': 'pr', 'wind': 'sfcWind'}[var]
    cands = []
    if subdir in (sp.M10_CanESM2_DIR, sp.M10_CanRCM4_DIR):
        cands.append(os.path.join(DATA, subdir,
                                  '%s_%s.nc' % (vname, sp._m10_file_tag(str(sid)))))
    else:
        # validated per-buoy Can-model extraction dir (cdo_extractions/ prefix)
        cands.append(os.path.join(CDO_CAN_BUOYS[subdir],
                                  '%s_buoy_%s.nc' % (vname, sid)))
    for f in cands:
        if not os.path.exists(f):
            continue
        try:
            nc = Dataset(f, mode='r')
            v = nc.variables[vname][:].astype(float).ravel()
            t = nc.variables['time'][:]
            units, cal = nc.variables['time'].units, nc.variables['time'].calendar
            nc.close()
            if subdir == sp.M10_CanESM2_DIR:
                lo, hi = sp.M10_EVAL_SLICE['canesm2']
            elif subdir == sp.M10_CanRCM4_DIR:
                lo, hi = sp.M10_EVAL_SLICE['canrcm4']
            else:
                lo, hi = 0, len(v)          # *_buoys files already 1986-2005
            v, t = v[lo:hi], t[lo:hi]
            dates = num2date(t, units, calendar=cal)
            s = pd.Series(v, index=pd.DatetimeIndex(
                [pd.Timestamp(int(dd.year), int(dd.month), int(dd.day))
                 for dd in dates]))
            if var == 't':
                s = s - 273.15
            elif var == 'pr':
                s = s * 86400.0            # daily-mean kg m-2 s-1 -> mm/day
            ann = sp.model_annual_from_series(s, 'sum' if var == 'pr' else 'mean')
            arr = np.full(N_YEARS, np.nan)
            for y, val in ann.items():
                if START_YEAR <= int(y) <= END_YEAR:
                    arr[int(y) - START_YEAR] = val
            if not np.all(np.isnan(arr)):
                return {sid: arr}
        except Exception:
            continue
    return {}


def _canesm2_annual(var, ids):
    """Land-only CanESM2 annuals from the M10 noaa113 dir.  Buoy ids (46*) are
    EXCLUDED here and overlaid from the validated buoy values in model_annual()
    (mirrors _wrf_annual, which skips 46*)."""
    out = {}
    for sid in ids:
        if sid.startswith('46'):
            continue
        a = _can_daily_annual(sid, var, sp.M10_CanESM2_DIR)
        out.update(a)
    return out


def _canrcm4_annual(var, ids):
    """Land-only CanRCM4 annuals from the M10 noaa113 dir (buoys overlaid
    separately in model_annual())."""
    out = {}
    for sid in ids:
        if sid.startswith('46'):
            continue
        a = _can_daily_annual(sid, var, sp.M10_CanRCM4_DIR)
        out.update(a)
    return out


# -- buoy WRF models: VALIDATED per-buoy x per-model annual arrays ----------
_BUOY_WRF_MODEL_CACHE = None
def _buoy_wrf_model_annuals():
    """dict label -> dict sid -> (20,) annual mean, for D01/D02/D03, from
    buoy_wind_tables.pkl (wind) and the D03 tas per-buoy pickle.  These are
    the exact arrays behind the published buoy tables / bias figure."""
    global _BUOY_WRF_MODEL_CACHE
    if _BUOY_WRF_MODEL_CACHE is not None:
        return _BUOY_WRF_MODEL_CACHE
    import pickle
    out = {}
    pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                       'buoy_wind_tables.pkl'), 'rb'))
    lab_map = {'WRF_d01': 'D01', 'WRF_d02': 'D02', 'WRF_d03': 'D03'}
    for bid, d in pk.get('per_buoy', {}).items():
        sid, lab = bid.split('|')
        lab = lab_map.get(lab)
        if lab is None:
            continue
        a = np.asarray(d.get('annual_model'), float).reshape(-1)
        if len(a) == N_YEARS and not np.all(np.isnan(a)):
            out.setdefault(('wind', lab), {})[sid] = a
    pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                       'buoy_tas_d03_per_buoy.pkl'), 'rb'))
    for sid, d in pk.items():
        o = np.asarray(d.get('annual_obs'), float).reshape(-1)
        if len(o) != N_YEARS or np.all(np.isnan(o)):
            continue
        # D03 column (validated D03-pinned set)
        m = np.asarray(d.get('annual_model'), float).reshape(-1)
        if len(m) == N_YEARS and not np.all(np.isnan(m)):
            out.setdefault(('t', 'D03'), {})[str(sid)] = m - 273.15
        # D01/D02 columns (validated wrf_nest_t2 selection, cascaded and
        # masked to the same D03 obs-surviving years) -- never from
        # data/wrf_canon, whose 46xxx t files are corrupt artifacts.
        for dom, key in (('D01', 'annual_model_d01'), ('D02', 'annual_model_d02')):
            md = np.asarray(d.get(key), float).reshape(-1)
            if len(md) == N_YEARS and not np.all(np.isnan(md)):
                out.setdefault(('t', dom), {})[str(sid)] = md - 273.15
    _BUOY_WRF_MODEL_CACHE = out
    return out


# -- buoy Can models: VALIDATED per-buoy annual arrays (D03-pinned set) -------
_CAN_BUOY_MODEL_CACHE = None
def _can_buoy_model_annuals():
    """dict (var, label) -> dict sid -> (20,) annual, for CanESM2/CanRCM4, on
    the filtered D03 buoy set (make_bias_figure.WIND_BUOYS_D03 / TAS_BUOYS_D03).
    These are the validated arrays behind the published buoy tables / bias
    figure -- the SAME obs reference the WRF columns are scored against.
      * wind: buoy_wind_tables.pkl per_buoy['<ID>|<label>']['annual_model']
      * t   : validated cdo_extractions/<can>_buoys/tas_buoy_<ID>.nc (K->C),
              via _can_daily_annual (corrected cdo_extractions/ path).
    The validated per-buoy tas table (buoy_tas_d03_per_buoy.pkl) carries WRF
    arrays only, so the t Can* values come from the validated extraction files.
    """
    global _CAN_BUOY_MODEL_CACHE
    if _CAN_BUOY_MODEL_CACHE is not None:
        return _CAN_BUOY_MODEL_CACHE
    import pickle
    import make_bias_figure as mbf
    out = {}
    # wind: from the validated buoy wind table (per-buoy Can* annual_model)
    pk = pickle.load(open(os.path.join(BASE, 'pickles',
                                       'buoy_wind_tables.pkl'), 'rb'))
    per = pk.get('per_buoy', {})
    for label in ('CanESM2', 'CanRCM4'):
        for sid in mbf.WIND_BUOYS_D03:
            d = per.get('%s|%s' % (sid, label))
            if d is None:
                continue
            a = np.asarray(d.get('annual_model'), float).reshape(-1)
            if len(a) == N_YEARS and not np.all(np.isnan(a)):
                out.setdefault(('wind', label), {})[str(sid)] = a
    # t: from the validated per-buoy Can extraction files (K -> C)
    for label, subdir in (('CanESM2', 'canesm2_raw_buoys'),
                          ('CanRCM4', 'canrcm4_buoys')):
        for sid in mbf.TAS_BUOYS_D03:
            a = _can_daily_annual(str(sid), 't', subdir)
            out.setdefault(('t', label), {}).update(a)
    _CAN_BUOY_MODEL_CACHE = out
    return out


def model_annual(var, ids):
    out = {}
    for label, dom in [('D03', 'd03'), ('D02', 'd02'), ('D01', 'd01')]:
        out[label] = _wrf_annual(ids, var, dom)
    out['CanESM2'] = _canesm2_annual(var, ids)
    out['CanRCM4'] = _canrcm4_annual(var, ids)
    # overlay the validated buoy WRF series (buoy ids are NOT read from
    # data/wrf_canon -- those files are corrupt D01-era artifacts)
    for (v, label), bysid in _buoy_wrf_model_annuals().items():
        if v == var and label in out:
            out[label].update({s: a for s, a in bysid.items() if s in ids})
    # overlay the validated buoy Can-model series so CanESM2/CanRCM4 are scored
    # on the same obs set (land + filtered D03 buoys) as the WRF columns
    for (v, label), bysid in _can_buoy_model_annuals().items():
        if v == var and label in out:
            out[label].update({s: a for s, a in bysid.items() if s in ids})
    return out


# ---------------------------------------------------------------------------
# Assemble + cache
# ---------------------------------------------------------------------------
def build(var, meta=None, use_cache=True):
    """Return (obs: dict sid->array, models: dict label->dict sid->array).
    Disk-cached in pickles/taylor_annuals.pkl so the figure re-runs in seconds."""
    if use_cache and os.path.exists(CACHE_PKL):
        import pickle
        try:
            cache = pickle.load(open(CACHE_PKL, 'rb'))
            if cache.get('_version') == CACHE_VERSION and var in cache:
                return cache[var]
        except Exception:
            pass
    meta = meta or load_meta()
    ids = station_set(var, meta)
    obs = obs_annual(var, meta)
    models = model_annual(var, ids)
    result = (obs, models)
    if use_cache:
        import pickle
        cache = {'_version': CACHE_VERSION, var: result}
        # preserve other already-built vars if present and same version
        if os.path.exists(CACHE_PKL):
            try:
                old = pickle.load(open(CACHE_PKL, 'rb'))
                if old.get('_version') == CACHE_VERSION:
                    for k in old:
                        if k not in cache and k != '_version':
                            cache[k] = old[k]
            except Exception:
                pass
        pickle.dump(cache, open(CACHE_PKL, 'wb'))
    return result


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
def _period(s):
    arrs = [np.asarray(x, dtype=float) for x in s.values() if x is not None]
    if not arrs:
        return np.array([])
    a = np.concatenate(arrs)
    return a[~np.isnan(a)]

if __name__ == '__main__':
    meta = load_meta()
    for var in ('t', 'pr', 'wind'):
        print('===', var, '===')
        obs, models = build(var, meta)
        o = _period(obs)
        print('  obs  n_st=%3d  n_ann=%4d  mean=%9.3f  std=%8.3f' %
              (len(obs), len(o), o.mean(), o.std()))
        for k, s in models.items():
            a = _period(s)
            if len(a):
                print('  %-8s n_st=%3d  n_ann=%4d  mean=%9.3f  std=%8.3f' %
                      (k, len(s), len(a), a.mean(), a.std()))
            else:
                print('  %-8s n_st=  0  (none)' % k)
