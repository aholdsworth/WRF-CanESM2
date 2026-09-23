#!/usr/bin/env python3
"""
seasonal_pooling.py — season-first pooling convention for bias & skill metrics.

Convention (agreed 2026-09-18, see seasonal_pooling_handoff_20260918.md):

  Pooling hierarchy (obs side):  day -> season -> year -> period
  with equal weight at each level.  Seasons are DJF / MAM / JJA / SON.

  Completeness gate ("10% rule" per season):
    * a season value is NaN if < 90% of its expected days (or hours, for
      hourly wind) are present.  Expected days/hours are computed from the
      actual calendar for each year (leap years included), never hardcoded.
    * a year is NaN if any of its 4 seasons is NaN.
    * the period value (1986-2005) is the mean of the surviving annual
      values, reported only if >= MIN_YEARS (15) of the 20 years survive.

  Like-for-like comparison:
    * the surviving-year set is defined by the OBS series (the only series
      that can have gaps); the model series (complete) is restricted to the
      same years before averaging.
    * bias = mean(model_years) - mean(obs_years)
             (difference of means; identical to mean of differences when the
              year sets match)
    * MAE  = mean( |model_year - obs_year| )        (error first, then mean)
    * RMSE = sqrt( mean( (model_year - obs_year)^2 ) )
    * reported per station: bias, MAE, RMSE, n_years; all NaN if
      n_years < MIN_YEARS.

  Citations for the pooling hierarchy:
    * Jones & McAvoy (2010), "Guide to Climate Indices", WMO/TUBITAK,
      section 3.2.1 (seasonal and annual values as means of lower-level
      values).
    * WMO Guide to Climatological Practices, WMO No. 100.

This module is the single source of truth for the obs-side gate/mask and the
skill() aggregation; it is shared by make_skill_tables.py (pickle producer)
and mae_vs_eva.py (fast replication of the "OUR" side) so the two can never
drift apart.
"""
import calendar
import datetime
import glob
import os
import re

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
START_YEAR, END_YEAR = 1986, 2005
YEARS = tuple(range(START_YEAR, END_YEAR + 1))
SEASONS = ('DJF', 'MAM', 'JJA', 'SON')
MIN_YEARS = 15        # min surviving years for a period value / skill report
GATE_FRAC = 0.90      # a season needs >= 90% of its expected days/hours

SEASON_OF_MONTH = {1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
                   6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON',
                   11: 'SON', 12: 'DJF'}

# BCH stations that record p24 in tenths of mm in the wide CSVs
# (verified 2026-09: their /10 annual sums are climatologically sensible).
BCH_TENTHS_MM = {'BLN', 'DLU', 'MIS', 'NTY', 'WOL'}

# ---------------------------------------------------------------------------
# Calendar
# ---------------------------------------------------------------------------
def season_expected_days(year, season):
    """Expected number of days in one calendar season (leap-aware)."""
    if season == 'DJF':
        return 31 + (29 if calendar.isleap(year) else 28) + 31
    if season == 'MAM':
        return 31 + 30 + 31
    if season == 'JJA':
        return 30 + 31 + 31
    if season == 'SON':
        return 31 + 30 + 31
    raise ValueError(season)

def expected_days_table(years=YEARS):
    """dict (year, season) -> expected days; hours = days * 24."""
    return {(y, s): season_expected_days(y, s) for y in years for s in SEASONS}

# ---------------------------------------------------------------------------
# Season / year aggregation (core of the convention)
# ---------------------------------------------------------------------------
def seasonal_year_frame(daily, kind, freq='D'):
    """Pool a daily (or hourly) obs series to a (20 years x 4 seasons) frame.

    daily : pd.Series indexed by datetime, values numeric (NaN = missing).
            For wind this is the HOURLY series (already in m/s).
    kind  : 'mean' (t, wind) or 'sum' (pr).
    freq  : 'D' for daily data (expected count = days), 'H' for hourly
            (expected count = hours = days*24).

    Returns DataFrame index=years, columns=SEASONS; NaN where the season is
    < GATE_FRAC complete.  A season with zero observed values is NaN.
    """
    exp = expected_days_table()
    s = daily.dropna()
    if len(s) == 0:
        return pd.DataFrame(np.nan, index=list(YEARS), columns=list(SEASONS))
    s = s[(s.index >= pd.Timestamp(START_YEAR, 1, 1)) &
          (s.index < pd.Timestamp(END_YEAR + 1, 1, 1))]
    if len(s) == 0:
        return pd.DataFrame(np.nan, index=list(YEARS), columns=list(SEASONS))

    year = pd.Series(s.index.year, index=s.index)
    seas = pd.Series(s.index.month, index=s.index).map(SEASON_OF_MONTH)

    vals = np.full((len(YEARS), len(SEASONS)), np.nan)
    ymap = {y: i for i, y in enumerate(YEARS)}
    smap = {st: i for i, st in enumerate(SEASONS)}
    for y in YEARS:
        for st in SEASONS:
            m = ((year == y) & (seas == st)).values
            n_obs = int(m.sum())
            n_exp = exp[(y, st)] * (24 if freq == 'H' else 1)
            if n_obs < GATE_FRAC * n_exp:
                continue
            v = s.values[m].sum() if kind == 'sum' else s.values[m].mean()
            vals[ymap[y], smap[st]] = float(v)
    df = pd.DataFrame(vals, index=list(YEARS), columns=list(SEASONS))
    df.index.name = 'year'
    return df

def year_mask(seas_df):
    """Surviving-year mask: a year survives iff all 4 seasons are non-NaN."""
    return seas_df.notna().all(axis=1)

def annual_from_seasons(seas_df):
    """Annual value per year = mean of its 4 seasonal values (NaN if any
    season is NaN).  For pr the seasonal values are seasonal SUMS, so the
    annual is the mean of the 4 seasonal sums."""
    ann = seas_df.mean(axis=1).where(seas_df.notna().all(axis=1))
    ann.name = 'annual'
    return ann

def period_value(annual):
    """1986-2005 period value: mean of surviving annual values, NaN unless
    >= MIN_YEARS years survive."""
    a = annual.dropna()
    if len(a) < MIN_YEARS:
        return np.nan
    return float(a.mean())

def obs_seasonal_pool(daily, kind, freq='D'):
    """Full obs-side pipeline for one station/variable.

    Returns dict with:
      seas        : (20 x 4) seasonal frame (pre-year-pooling)
      annual      : Series (20,) annual values (NaN for dropped years)
      year_mask   : Series (20,) bool
      years       : list[int] surviving years
      n_years     : int
      period      : float | NaN  (period value of the OBS series)
      dropped     : list of (year, reason) for dropped years
    """
    seas = seasonal_year_frame(daily, kind, freq=freq)
    mask = year_mask(seas)
    annual = annual_from_seasons(seas)
    dropped = []
    for y in YEARS:
        if not mask.loc[y]:
            bad = [st for st in SEASONS if np.isnan(seas.loc[y, st])]
            dropped.append((y, 'season(s) <90% complete: ' + ','.join(bad)))
    return dict(seas=seas, annual=annual, year_mask=mask,
                years=[int(y) for y in YEARS if mask.loc[y]],
                n_years=int(mask.sum()),
                period=period_value(annual), dropped=dropped)

# ---------------------------------------------------------------------------
# Skill metrics (model vs obs over shared surviving years)
# ---------------------------------------------------------------------------
def skill(model_annual, obs_annual, min_years=MIN_YEARS):
    """bias / MAE / RMSE / n_years over the shared surviving years.

    model_annual : pd.Series or dict, annual model values (index/keys = year)
    obs_annual   : pd.Series, annual obs values (NaN = dropped year)

    The model series is assumed complete; the shared year set is the obs
    surviving years (obs_annual non-NaN).  Returns (bias, mae, rmse, n) or
    (NaN, NaN, NaN, n) if fewer than min_years shared years.
    """
    if isinstance(model_annual, dict):
        model_annual = pd.Series(model_annual)
    m = model_annual.reindex(list(YEARS)).astype(float)
    o = obs_annual.reindex(list(YEARS)).astype(float)
    ok = o.notna() & m.notna()
    if int(ok.sum()) < min_years:
        return np.nan, np.nan, np.nan, int(ok.sum())
    dm = (m - o)[ok]
    bias = float(dm.mean())
    mae = float(dm.abs().mean())
    rmse = float(np.sqrt((dm ** 2).mean()))
    return bias, mae, rmse, int(ok.sum())

def pct_bias(model_period, obs_period):
    """Percentage bias with the obs>0 guard (pr, wind)."""
    if model_period is None or obs_period is None:
        return np.nan
    if np.isnan(obs_period) or obs_period <= 0:
        return np.nan
    return 100.0 * (model_period - obs_period) / obs_period

# ---------------------------------------------------------------------------
# Model-side annual values (model data is complete; just restrict to years)
# ---------------------------------------------------------------------------
def model_annual_from_series(series, kind):
    """Annual model values from a raw series (hourly WRF, daily gridded).
    kind: 'sum' (pr increments) or 'mean' (t already converted to C, wind m/s)."""
    s = series[(series.index >= pd.Timestamp(START_YEAR, 1, 1)) &
               (series.index < pd.Timestamp(END_YEAR + 1, 1, 1))].dropna()
    if len(s) == 0:
        return pd.Series(dtype=float)
    yr = s.index.year
    if kind == 'sum':
        return s.groupby(yr).sum()
    return s.groupby(yr).mean()

# ---------------------------------------------------------------------------
# Observations: per-agency daily-series builders
# ---------------------------------------------------------------------------
def eccc_daily_series(sid, kind, base_dir):
    """ECCC daily CSVs -> Series indexed by date.
    kind: 't' -> 'Mean Temp (C)'; 'pr' -> 'Total Rain (mm)' (Eva's column)."""
    d = os.path.join(base_dir, 'eccc_obs', str(sid))
    if not os.path.isdir(d):
        return None
    col = 'Mean Temp (\u00b0C)' if kind == 't' else 'Total Rain (mm)'
    files = sorted(glob.glob(os.path.join(d, '*_P1D.csv')))
    if not files:
        return None
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, usecols=['Date/Time', col], encoding='utf-8-sig')
            dfs.append(df)
        except Exception:
            continue
    if not dfs:
        return None
    df = pd.concat(dfs, ignore_index=True)
    idx = pd.to_datetime(df['Date/Time'])
    s = pd.to_numeric(df[col], errors='coerce')
    s.index = idx
    return s

def eccc_hourly_wind_series(sid, base_dir):
    """ECCC hourly CSVs -> Series of 'Wind Spd (km/h)' in m/s, hourly index."""
    d = os.path.join(base_dir, 'eccc_hourly_obs', str(sid))
    if not os.path.isdir(d):
        return None
    files = sorted(glob.glob(os.path.join(d, '*_P1H.csv')))
    if not files:
        return None
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, usecols=['Date/Time (UTC)', 'Wind Spd (km/h)'],
                             encoding='utf-8-sig')
            dfs.append(df)
        except Exception:
            continue
    if not dfs:
        return None
    df = pd.concat(dfs, ignore_index=True)
    idx = pd.to_datetime(df['Date/Time (UTC)'])
    s = pd.to_numeric(df['Wind Spd (km/h)'], errors='coerce') / 3.6
    s.index = idx
    return s

_NOAA_CACHE = {}
def noaa_read(base_dir):
    """Lazy-load the 4 GHCN-d 5-yr block CSVs -> dict sid -> DataFrame."""
    if _NOAA_CACHE:
        return _NOAA_CACHE
    cols = ['STATION', 'DATE', 'TAVG', 'TMIN', 'TMAX', 'PRCP', 'AWND']
    files = [os.path.join(base_dir, 'noaa_obs', f'{y0}-{y0+4}.csv')
             for y0 in (1986, 1991, 1996, 2001)]
    if not all(os.path.exists(p) for p in files):
        return {}
    dfs = []
    for p in files:
        df = pd.read_csv(p, usecols=cols, low_memory=False)
        df['DATE'] = pd.to_datetime(df['DATE'])
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    out = {}
    for s, g in df.groupby('STATION'):
        g = g.copy()
        for c in ['TAVG', 'TMIN', 'TMAX', 'PRCP', 'AWND']:
            g[c] = pd.to_numeric(g[c], errors='coerce')
        out[str(s)] = g
    _NOAA_CACHE.update(out)
    return _NOAA_CACHE

def noaa_daily_series(sid, kind, base_dir):
    """NOAA GHCN-d -> daily Series.
    kind: 't'    -> TAVG where present, else (TMIN+TMAX)/2 (daily-level
                    fallback; pooling is what changed, not variable choice)
          'pr'   -> PRCP in mm (GHCN-D PRCP column is already 0.01 mm, i.e. mm)
          'wind' -> AWND in m/s (mph * 0.44704)
    """
    g = noaa_read(base_dir).get(str(sid))
    if g is None or len(g) == 0:
        return None
    idx = g['DATE']
    if kind == 't':
        s = g['TAVG'].fillna((g['TMIN'] + g['TMAX']) / 2.0)
    elif kind == 'pr':
        s = g['PRCP']
    elif kind == 'wind':
        s = g['AWND'] * 0.44704
    else:
        return None
    s = s.copy()
    s.index = idx
    return s

_BCH_CACHE = {}
def bch_read(base_dir):
    """Lazy-load the BCH wide CSVs -> dict kind -> DataFrame (date x station)."""
    if _BCH_CACHE:
        return _BCH_CACHE
    out = {}
    for kind, fn in (('tn', 'BCH_tnId.csv'), ('tx', 'BCH_txId.csv'),
                     ('p24', 'BCH_p24Id.csv')):
        p = os.path.join(base_dir, 'bch_obs', fn)
        if not os.path.exists(p):
            return {}
        d = pd.read_csv(p, index_col=0)
        d.index = pd.to_datetime(d.index).tz_localize(None)
        out[kind] = d
    _BCH_CACHE.update(out)
    return _BCH_CACHE

def bch_daily_series(sid, kind, base_dir):
    """BCH daily wide CSVs -> Series.
    kind: 't'  -> (tn+tx)/2 ; 'pr' -> p24 (tenths-mm stations /10)."""
    d = bch_read(base_dir)
    if not d:
        return None
    if kind == 't':
        if str(sid) not in d['tn'].columns or str(sid) not in d['tx'].columns:
            return None
        s = (d['tn'][str(sid)] + d['tx'][str(sid)]) / 2.0
    elif kind == 'pr':
        if str(sid) not in d['p24'].columns:
            return None
        s = d['p24'][str(sid)].copy()
        if sid in BCH_TENTHS_MM:
            s = s / 10.0
    else:
        return None   # no BCH wind
    return s

# ---------------------------------------------------------------------------
# Convenience: full obs pipeline for one station & variable
# ---------------------------------------------------------------------------
def obs_pool(sid, agency, kind, base_dir):
    """agency: 'E' (ECCC), 'N' (NOAA), 'B' (BCH).
    kind: 't', 'pr', 'wind'.
    Returns the obs_seasonal_pool() dict, or None if no data / no such var
    for the agency (e.g. BCH wind)."""
    if agency == 'E':
        if kind == 'wind':
            series, freq = eccc_hourly_wind_series(sid, base_dir), 'H'
        else:
            series, freq = eccc_daily_series(sid, kind, base_dir), 'D'
    elif agency == 'N':
        series, freq = noaa_daily_series(sid, kind, base_dir), 'D'
    elif agency == 'B':
        series, freq = bch_daily_series(sid, kind, base_dir), 'D'
    else:
        return None
    if series is None or len(series.dropna()) == 0:
        return None
    pk = 'sum' if kind == 'pr' else 'mean'
    return obs_seasonal_pool(series, pk, freq=freq)

# ---------------------------------------------------------------------------
# CanRCM4 / CanESM2 raw per-station CDO extractions (113-NOAA common set)
#
# These are the CDO bicubic per-station daily files built by station_cdo/
# (see DATA_INVENTORY.md M9/M10). They REPLACE the old bilinear stand-in
# annuals in data/canrcm4_stations/ (M6), which had wrong station coords and
# was deleted 2026-09-21.
#
#   CanRCM4  : data/cdo_extractions/canrcm4_noaa113/<var>_NOAA_<ID>.nc
#              1986-2005, 7300 days, 365_day, epoch 1949-12-01.
#   CanESM2  : data/cdo_extractions/canesm2_raw_noaa113/<var>_NOAA_<ID>.nc
#              1850-2005, 56940 days, 365_day, epoch 1850-1-1 (mask to
#              1986-2005 downstream = indices 49640..56939).
#
# var names in the files: tas [K], pr [kg m-2 s-1 daily-mean rate], sfcWind.
# ---------------------------------------------------------------------------
M10_CanRCM4_DIR  = 'cdo_extractions/canrcm4_noaa113'
M10_CanESM2_DIR  = 'cdo_extractions/canesm2_raw_noaa113'
M10_VARNAME = {'t': 'tas', 'pr': 'pr', 'wind': 'sfcWind'}
M10_N_DAYS = {'canrcm4': 7300, 'canesm2': 56940}
M10_EVAL_SLICE = {'canrcm4': (0, 7300), 'canesm2': (49640, 56940)}

# ECCC M9/M10 files are tagged with the WMO/numeric pipeline sid (e.g. 1016335),
# but the pipeline land_ids use the ClimateID (e.g. 82).  Resolve ClimateID ->
# WMO sid from files/ECCC_d03_stations.csv (col3 = WMO sid, col4 = ClimateID).
_ECCC_CID2SID = None
def _eccc_cid2sid():
    global _ECCC_CID2SID
    if _ECCC_CID2SID is None:
        import pandas as pd
        base = os.path.dirname(os.path.abspath(__file__))
        try:
            df = pd.read_csv(os.path.join(base, 'files', 'ECCC_d03_stations.csv'), header=None)
            _ECCC_CID2SID = dict(zip(df.iloc[:,4].astype(str), df.iloc[:,3].astype(str)))
        except Exception:
            _ECCC_CID2SID = {}
    return _ECCC_CID2SID


def _m10_file_tag(sid):
    """Agency tag for an M10 station file name.
    NOAA sids start with USC/USW; BCH sids are 3-letter codes (ALU, BCK, ...);
    everything else is ECCC.  For ECCC, a numeric ClimateID (e.g. 82) is
    resolved to its WMO/numeric pipeline sid (e.g. 1016335), which is how the
    M9/M10 files are tagged.  A sid already in WMO form (e.g. 1113581) or a
    mixed id (1101N65) passes through unchanged."""
    if sid.startswith(('USC', 'USW')):
        return f'NOAA_{sid}'
    if re.fullmatch(r'[A-Z]{3}', sid):
        return f'BCH_{sid}'
    m = _eccc_cid2sid()
    return f'ECCC_{m.get(sid, sid)}'


def _m10_file(base_dir, model_dir, vname, sid):
    return os.path.join(base_dir, model_dir, f'{vname}_{_m10_file_tag(sid)}.nc')


def canrcm4_daily_series(sid, kind, base_dir):
    """Daily CanRCM4 (M10) station series, index=1986-01-01..2005-12-30, as a
    pd.Series. t -> C, pr -> mm (daily), wind -> m/s. NaN if no file."""
    from netCDF4 import Dataset
    vname = M10_VARNAME[kind]
    f = _m10_file(base_dir, M10_CanRCM4_DIR, vname, sid)
    if not os.path.exists(f):
        return None
    nc = Dataset(f, mode='r')
    v = nc.variables[vname][:].astype(float).ravel()
    t = nc.variables['time'][:]
    units, cal = nc.variables['time'].units, nc.variables['time'].calendar
    nc.close()
    from cftime import num2date
    dates = num2date(t, units, calendar=cal)
    idx = pd.DatetimeIndex([pd.Timestamp(int(d.year), int(d.month), int(d.day))
                            for d in dates])
    s = pd.Series(v, index=idx)
    if kind == 't':
        s = s - 273.15
    elif kind == 'pr':
        s = s * 86400.0          # daily-mean kg m-2 s-1 -> mm per day
    return s


def cdo_station_value(sid, kind, base_dir, model='canrcm4'):
    """1986-2005 period value for one station from an M9/M10 CDO file.
    kind: 't' (C mean), 'pr' (mm/yr mean of yearly sums), 'wind' (m/s mean).
    Returns float or NaN if the file is absent."""
    vname = M10_VARNAME[kind]
    if model == 'canrcm4':
        d = os.path.join(base_dir, M10_CanRCM4_DIR)
    else:
        d = os.path.join(base_dir, M10_CanESM2_DIR)
    f = _m10_file(base_dir, d, vname, sid)
    if not os.path.exists(f):
        return np.nan
    from netCDF4 import Dataset
    nc = Dataset(f, mode='r')
    v = nc.variables[vname][:].astype(float).ravel()
    t = nc.variables['time'][:]
    units, cal = nc.variables['time'].units, nc.variables['time'].calendar
    nc.close()
    lo, hi = M10_EVAL_SLICE[model]
    v = v[lo:hi]
    from cftime import num2date
    dates = num2date(t[lo:hi], units, calendar=cal)
    years = np.array([d.year for d in dates])
    if kind == 't':
        return float(np.nanmean(v - 273.15))
    if kind == 'wind':
        return float(np.nanmean(v))
    # pr: mm per day -> sum per calendar year -> mean over years
    mmday = v * 86400.0
    ys = np.unique(years)
    out = np.array([mmday[years == y].sum() for y in ys])
    return float(np.nanmean(out))
