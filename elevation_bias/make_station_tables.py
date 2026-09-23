#!/usr/bin/env python
"""Build per-station bias tables for the 5x4 elevation-bias diagnostic grid.

For each of the 5 datasets (WRF D03, D02, D01, CanESM2, CanRCM4) produce a
DataFrame with columns:
    station, lat, lon, obs_elev_m,
    elevation_bias_m,     # model grid elevation - station elevation
    temperature_bias_C,   # annual mean (1986-2005) model T2/tas - obs T
    precipitation_bias_pct, # annual mean (1986-2005) 100*(model-pr - obs-pr)/obs-pr
    wind_bias_pct         # annual mean (1986-2005) 100*(model-wind - obs-wind)/obs-wind

Model data:
  * WRF d01/d02/d03 : hourly per-station files in
      <PFM>/<run>/variables_stations/{t,pr,wind}_{d0X}_*.nc
    (d01 ECCC/NOAA files: {t,pr}_d01_{AGENCY}_{ID}.nc)
  * CanESM2 raw     : gridded <PFM>/gridded_model_daily/CanESM2_raw/{tas,pr,wind}_hist.nc
    sampled at the nearest grid point of each station.
  * CanRCM4         : gridded <PFM>/gridded_model_daily/CanRCM4/*_NAM22_hist.nc
    sampled at the nearest grid point of each station.

Observations (historical period 1986-2005, 10-yr annual means):
  * ECCC: /gpfs/fs7/dfo/hpcmc/pfm/evg000/obs/1986_2005_stations/<ID>/
          en_climate_daily_BC_<ClimateID>_<year>_P1D.csv  (daily, 1986-2005)
          T  = mean of 'Mean Temp (°C)'
          PR = sum of 'Total Rain (mm)' per year (Eva's convention), mm/yr
          wind = not available in the daily files -> NaN
  * NOAA (USC...): NOT in mounted paths -> NaN (loader keeps the original
          /daily/NOAA/{1986-1990,...,2001-2005}.csv reading if the dir appears)
  * BCH (ALU, ...): NOT in mounted paths -> NaN (original BCH_tn/tx/p24 CSVs)

Elevation bias (model grid topography - station elevation) is recomputed here
with the same nearest-gridpoint logic as plot_elev_bias_revised.py.

Run:  /home/amh001/space_fs7/software_2022/python/py_2024/bin/python make_station_tables.py
Output: <BASE>/pickles/station_tables.pkl  (dict of DataFrames, one per dataset)
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

BASE = os.environ.get('ELEV_BIAS_BASE',
    '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS')
PFM  = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF'
EVOBS = '/gpfs/fs7/dfo/hpcmc/pfm/evg000'
RUN  = 'historical'
START_YEAR, END_YEAR = 1986, 2005

sys.path.insert(0, BASE)
from plot_elev_bias_revised import (load_eccc, load_bch, load_noaa_land,
                                    load_buoys, read_geo_em, read_orog,
                                    weighted_model_elevs, _can_weight_file,
                                    modeled_elevs)
import seasonal_pooling as sp

DOMAIN_DIR = os.path.join(BASE, 'plotting_files')
STATION_FILES_DIR = os.path.join(PFM, 'station_files')
DATA = os.path.join(BASE, 'data')
ECCC_OBS_DIR  = os.path.join(DATA, 'eccc_obs')
NOAA_OBS_DIR  = os.path.join(DATA, 'noaa_obs')   # NOT YET COPIED (absent from mounts) -> NaN
BCH_OBS_DIR   = os.path.join(DATA, 'bch_obs')     # NOT YET COPIED (absent from mounts) -> NaN
WRF_ST_DIR   = os.path.join(DATA, 'wrf_stations')
ECCC_HOURLY_DIR = os.path.join(DATA, 'eccc_hourly_obs')  # <ID>/*.csv, 1986-2005 hourly
RAW_DIR      = os.path.join(DATA, 'canesm2_raw')
RCM_DIR      = os.path.join(DATA, 'canrcm4')
# (M6 canrcm4_stations/ deleted 2026-09-21; CanRCM4 rows use M10 via sp.cdo_station_value)
OUT_PKL      = os.path.join(BASE, 'pickles', 'station_tables.pkl')
os.makedirs(os.path.dirname(OUT_PKL), exist_ok=True)

# ---------------------------------------------------------------------------
# Station metadata
# ---------------------------------------------------------------------------
eccc_ids, eccc_lats, eccc_lons, eccc_elev = load_eccc()
# 951 (HOPE SLIDE, BC) is now INCLUDED (2026-09-21): it has CDO per-station
# model data (M9/M10, pipeline sid 1113581) and is in the 177-station
# inventory.  The old EXCLUDED_STATIONS = {'951'} (consistency to-do,
# 2026-09-17, DATA_README.md Known issues #4) is removed.
EXCLUDED_STATIONS = set()
EXCLUDED_LAND = [s for s in EXCLUDED_STATIONS]
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

lats_s = _lats.astype(float)
lons_s = _lons.astype(float)
elev_s = _elev.astype(float)
# FIX: all arrays must be aligned to all_ids order (natural: ECCC, BCH, NOAA,
# buoys).  The raw concatenated series are string-sorted by index, which is a
# DIFFERENT order from all_ids -> reindex explicitly.
lats_s = lats_s.reindex([str(s) for s in all_ids])
lons_s = lons_s.reindex([str(s) for s in all_ids])
elev_s = elev_s.reindex([str(s) for s in all_ids])
assert lats_s.notna().all() and lons_s.notna().all() and elev_s.notna().all(), \
    'station metadata misalignment after reindex'

# ---------------------------------------------------------------------------
# Topography & elevation bias (same as plot_elev_bias_revised.py)
# ---------------------------------------------------------------------------
def topo_lats_lons(lat2d, lon2d):
    """Flatten a 2D (lat, lon) pair into a 1D array of (lon, lat) points."""
    return np.column_stack([np.ravel(lon2d), np.ravel(lat2d)])

def make_interp(lat2d, lon2d, topo2d):
    pts = topo_lats_lons(lat2d, lon2d)
    vals = np.ravel(topo2d)
    return NearestNDInterpolator(pts, vals)

# WRF elevation bias uses the extraction-consistent CDO bicubic weights (the
# orography actually baked into the model data we compare).  Any station without
# a usable T-file weight file (e.g. buoys without T files) falls back to the
# corner nearest-grid-point lookup over HGT_M.  CanESM2/CanRCM4 use the same
# extraction-consistent method with the exact M9/M10 weight files (see
# plot_elev_bias_revised._can_model_elevs), falling back to nearest-GP over orog
# only where no weight file exists; `elev_method` records which was used.
lat_d03, lon_d03, topo_d03 = read_geo_em('geo_em.d03.nc')
lat_d02, lon_d02, topo_d02 = read_geo_em('geo_em.d02.nc')
lat_d01, lon_d01, topo_d01 = read_geo_em('geo_em.d01.nc')
lats_c2, lons_c2, topo_c2 = read_orog('orog_CanESM2.nc')
lats_r4, lons_r4, topo_r4 = read_orog('orog_CanRCM4.nc')
pts_all = np.column_stack([lons_s.values, lats_s.values])

elev_bias = {}
elev_method = {}

# CanESM2 & CanRCM4: weight-consistent model elevation (fallback: nearest-GP)
for name, (lat2d, lon2d, topo2d) in [('CanESM2', (lats_c2, lons_c2, topo_c2)),
                                     ('CanRCM4', (lats_r4, lons_r4, topo_r4))]:
    mtopo, _w = weighted_model_elevs(name, list(all_ids))
    warn_ids = [s for s in mtopo[mtopo.isna()].index.tolist() if s not in buoy_set]
    if warn_ids:
        nearest = modeled_elevs(lat2d, lon2d, topo2d)
        for sid in warn_ids:
            print(f'  WARN elev {name}: {sid}: no weight file -> nearest-GP '
                  f'(model_elev={nearest.loc[sid, "elev"]:.1f})')
            mtopo.loc[sid] = nearest.loc[sid, 'elev']
    # Buoys stay NaN (excluded from the elevation diagnostic); assert land resolves.
    _land = [s for s in all_ids if s not in buoy_set]
    assert mtopo.loc[_land].notna().all(), f'{name}: unresolved land model elevation'
    # Elevation bias is LAND-ONLY: buoys have no meaningful orography
    # (obs_elev=0, and CanRCM4/CanESM2 orog marks their coastal cells as
    # ocean=0 m), so they were stacked at bias=0 -> the vertical-line
    # artefact.  Exclude them from the elevation diagnostic entirely.
    _eb = mtopo.values - elev_s.values
    _eb[[s in buoy_set for s in all_ids]] = np.nan
    elev_bias[name] = pd.Series(_eb, index=[str(s) for s in all_ids])
    elev_method[name] = pd.Series(
        ['bicubic' if w is not None else 'nearest'
         for w in [_can_weight_file(name, s) for s in all_ids]],
        index=[str(s) for s in all_ids])

# WRF D01/D02/D03: weight-based model elevation with corner fallback.
for name, dom, (clat, clon, ctopo) in [('D03', 'd03', (lat_d03, lon_d03, topo_d03)),
                                       ('D02', 'd02', (lat_d02, lon_d02, topo_d02)),
                                       ('D01', 'd01', (lat_d01, lon_d01, topo_d01))]:
    wmelev, wwarns = weighted_model_elevs(dom, list(all_ids))
    wmelev_orig = wmelev.copy()  # NaN where the weight lookup failed
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
        print(f'  WARN elev {name}: {w}')
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
    # WRF elev_method: 'bicubic' where a T-file weight was used, 'corner' for
    # the logged fallbacks (buoys / stations without T files).
    elev_method[name] = pd.Series(
        ['bicubic' if not np.isnan(orig) else 'corner'
         for orig in wmelev_orig.values], index=[str(s) for s in all_ids])

# ---------------------------------------------------------------------------
# Model data: WRF per-station hourly files
# ---------------------------------------------------------------------------
def wrf_st_files(var, dom):
    """Map station id -> file path for WRF per-station files of one var/domain."""
    out = {}
    d = WRF_ST_DIR
    # d02/d03 naming: {var}_{dom}_st{ID}.nc  (ID may be numeric ECCC or USxxx NOAA)
    # d01  ECCC/NOAA/BCH naming: {var}_d01_{AGENCY}_{ID}.nc
    #      BCH d01 files on PFM are ALSO untagged: {var}_d01_{ID}.nc (3 letters)
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
                sid = rest[2:]
                if sid.startswith('USC') or sid.startswith('USW'):
                    out['N' + sid] = os.path.join(d, f)
                else:
                    out['E' + sid] = os.path.join(d, f)
            elif rest.startswith('BCH_'):
                out['B' + rest[4:]] = os.path.join(d, f)
            else:
                import sys
                print(f'WARN wrf_st_files: untagged {dom} file skipped: {f}', file=sys.stderr)
    return out

def wrf_station_series(path, var_name, kind):
    """Load one per-station WRF file; return Series indexed by date (hours).

    kind: 't'  -> T2 - 273.15
          'pr' -> hourly increment (diff of cumulative pr, clip >=0)
          'wind' -> wspd as-is
    """
    nc = Dataset(path, mode='r')
    t = nc.variables['time'][:]
    data = np.squeeze(nc.variables[var_name][:])
    nc.close()
    dates = pd.to_datetime(
        datetime.datetime(START_YEAR, 1, 1) + pd.to_timedelta(t, unit='h'))
    if kind == 't':
        v = (data - 273.15)
    elif kind == 'pr':
        v = np.clip(np.diff(data), 0.0, None)
        dates = dates[1:]
    elif kind == 'wind':
        v = data
    return pd.Series(v, index=dates)

def annual_mean(series, kind):
    """10-year (1986-2005) annual value: mean for t/wind, sum for pr (mm/yr)."""
    s = series[(series.index >= datetime.datetime(START_YEAR, 1, 1)) &
               (series.index <  datetime.datetime(END_YEAR + 1, 1, 1))]
    if kind == 'pr':
        return s.sum()
    return s.mean()

def wrf_annual_means(agency_ids, tag_prefix, var, kind, dom):
    """Return dict station_id -> annual value for WRF files of one domain."""
    files = wrf_st_files(var, dom)
    out = {}
    for sid in agency_ids:
        key = tag_prefix + str(sid)
        p = files.get(key)
        if p is None:
            continue
        try:
            series = wrf_station_series(p, {'t': 'T2', 'pr': 'pr', 'wind': 'wspd'}[kind], kind)
            out[str(sid)] = annual_mean(series, kind)
        except Exception as e:
            print(f"  WARN {dom} {var} {sid}: {e}")
    return out

# ---------------------------------------------------------------------------
# Model data: gridded CanESM2 raw / CanRCM4 (nearest grid point)
# ---------------------------------------------------------------------------
def gridded_station_values(dirpath, varfile, var_name, kind):
    """Load one gridded file, return array (nstation) of 1986-2005 annual value.

    kind: 't' (K->C mean), 'pr' (daily flux kg m-2 s-1 -> mm/yr), 'wind' (m/s mean)
    Files are daily (365_day calendar); cftime is used to build the date mask.
    """
    from cftime import num2date
    nc = Dataset(os.path.join(dirpath, varfile), mode='r')
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    data = np.squeeze(nc.variables[var_name][:])
    t = nc.variables['time'][:]
    units = nc.variables['time'].units
    nc.close()
    if lat.ndim == 1:
        lons2d, lats2d = np.meshgrid(lon, lat)
    else:
        lons2d, lats2d = lon, lat
    dates = num2date(t, units)
    # NearestNDInterpolator wants values shaped (npts, k); data is (nt, npts)
    ip = NearestNDInterpolator(np.column_stack([lons2d.ravel(), lats2d.ravel()]),
                               data.reshape(len(t), -1).T)
    vals = ip(np.column_stack([lons_s.values, lats_s.values]))  # (nstn, nt)
    if kind == 't':
        vals = vals - 273.15
    mask = np.array([(d.year >= START_YEAR) and (d.year <= END_YEAR) for d in dates])
    period = vals[:, mask]   # (nstn, nperiod)
    if kind == 'pr':
        # daily mean flux -> depth mm; sum per calendar year, mean over years
        depths = period * 86400.0   # mm per day (1 kg m-2 = 1 mm)
        years = np.array([d.year for d in dates if (d.year >= START_YEAR) and (d.year <= END_YEAR)])
        ys = np.unique(years)
        out = np.array([depths[:, years == y].sum(axis=1) for y in ys])
        out = out.mean(axis=0)       # (nstn) mm/yr
    else:
        out = period.mean(axis=1)    # (nstn)
    return out

# ---------------------------------------------------------------------------
# Observations
#   ECCC: per-station daily CSVs (Eva's "en_climate_daily" format)
#         <ECCC_OBS_DIR>/<ID>/en_climate_daily_BC_<ClimateID>_<year>_P1D.csv
#   NOAA: <NOAA_OBS_DIR>/1986-1990.csv ... 2001-2005.csv (station-indexed,
#         DATE column, TMIN/TMAX/PRCP inches/AWND mph)  [not in current mount]
#   BCH : <BCH_OBS_DIR>/BCH_{tn,tx,p24}<ID>.csv          [not in current mount]
# ---------------------------------------------------------------------------
def _eccc_hourly_wind(sid):
    """1986-2005 ECCC hourly 'Wind Spd (km/h)' -> mean over all hourly values
    (km/h * 0.277778 -> m/s)."""
    d = os.path.join(ECCC_HOURLY_DIR, str(sid))
    if not os.path.isdir(d):
        return np.nan
    files = sorted(glob.glob(os.path.join(d, '*_P1H.csv')))
    if not files:
        return np.nan
    try:
        dfs = []
        for f in files:
            df = pd.read_csv(f, usecols=['Date/Time (UTC)', 'Wind Spd (km/h)'],
                             encoding='utf-8-sig')
            df = df.dropna(subset=['Wind Spd (km/h)'])
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        df['Date/Time (UTC)'] = pd.to_datetime(df['Date/Time (UTC)'])
        df = df[(df['Date/Time (UTC)'] >= datetime.datetime(START_YEAR, 1, 1)) &
                (df['Date/Time (UTC)'] <  datetime.datetime(END_YEAR + 1, 1, 1))]
        v = df['Wind Spd (km/h)'].astype(float)
        if len(v) == 0:
            return np.nan
        return float(v.mean() * 0.277778)   # km/h -> m/s
    except Exception as e:
        print(f"  WARN obs ECCC hourly wind {sid}: {e}")
    return np.nan

def _eccc_elev_obs(sid, kind):
    """1986-2005 ECCC obs: T = mean 'Mean Temp (C)'; PR = sum 'Total Rain (mm)'."""
    if kind == 'wind':
        return _eccc_hourly_wind(sid)
    if kind not in ('t', 'pr'):
        return np.nan
    d = os.path.join(ECCC_OBS_DIR, str(sid))
    if not os.path.isdir(d):
        return np.nan
    col = 'Mean Temp (°C)' if kind == 't' else 'Total Rain (mm)'
    files = sorted(glob.glob(os.path.join(d, f'*_{START_YEAR}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{START_YEAR + 1}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{START_YEAR + 2}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{START_YEAR + 3}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{START_YEAR + 4}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{END_YEAR - 4}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{END_YEAR - 3}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{END_YEAR - 2}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{END_YEAR - 1}_P1D.csv')) +
                   glob.glob(os.path.join(d, f'*_{END_YEAR}_P1D.csv')))
    if not files:
        return np.nan
    try:
        dfs = []
        for f in files:
            df = pd.read_csv(f, usecols=['Date/Time', col], encoding='utf-8-sig')
            df = df.dropna(subset=[col])
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        df['Date/Time'] = pd.to_datetime(df['Date/Time'])
        df = df[(df['Date/Time'] >= datetime.datetime(START_YEAR, 1, 1)) &
                (df['Date/Time'] <  datetime.datetime(END_YEAR + 1, 1, 1))]
        if len(df) == 0:
            return np.nan
        v = df[col].astype(float)
        if kind == 't':
            return float(v.mean())
        # pr: sum per calendar year, then mean over years -> mm/yr
        yearly = v.groupby(df['Date/Time'].dt.year).sum()
        return float(yearly.mean())
    except Exception as e:
        print(f"  WARN obs ECCC {sid} {kind}: {e}")
        return np.nan

_NOAA_CACHE = {}
def _noaa_read():
    """Lazy-load the 4 long-format GHCN-d 5-yr block CSVs once, returning
    dict station_id -> DataFrame(DATE, TAVG, TMIN, TMAX, PRCP, AWND),
    restricted to 1986-2005."""
    if _NOAA_CACHE:
        return _NOAA_CACHE
    cols = ['STATION', 'DATE', 'TAVG', 'TMIN', 'TMAX', 'PRCP', 'AWND']
    files = [os.path.join(NOAA_OBS_DIR, f'{y0}-{y0+4}.csv')
             for y0 in (1986, 1991, 1996, 2001)]
    if not all(os.path.exists(p) for p in files):
        return {}
    dfs = []
    for p in files:
        df = pd.read_csv(p, usecols=cols, low_memory=False)
        df['DATE'] = pd.to_datetime(df['DATE'])
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df = df[(df['DATE'] >= datetime.datetime(START_YEAR, 1, 1)) &
            (df['DATE'] <  datetime.datetime(END_YEAR + 1, 1, 1))]
    out = {}
    for s, g in df.groupby('STATION'):
        g = g.copy()
        for c in ['TAVG', 'TMIN', 'TMAX', 'PRCP', 'AWND']:
            g[c] = pd.to_numeric(g[c], errors='coerce')
        out[s] = g
    _NOAA_CACHE.update(out)
    return _NOAA_CACHE

def _noaa_obs_annual(sid, kind):
    """NOAA GHCN-d long-format CSVs: T = mean TAVG (C), PR = yearly-sum
    PRCP (inches in the GHCN-d CSV) x25.4 -> mm/yr, wind = mean AWND (mph) x0.44704 -> m/s.
    If TAVG is absent for the station, fall back to mean((TMIN+TMAX)/2).

    NOTE (2026-09-23 noaainch fix): the 2026-09-23 "prmm fix" wrongly dropped this
    x25.4 on the assumption the PRCP column was already mm. It is in INCHES
    (see DATA_INVENTORY.md O3, DATA_README.md, Missing_Obs_Handoff.md B1, and the
    data itself: Beaverton USC00350595 annual sum ~891 raw -> 22633 mm/yr x25.4,
    matching model pr ~20521 mm/yr -> ~-9% D03 bias; raw 891 gave a bogus +2203%).
    Restored here."""
    try:
        g = _noaa_read().get(sid)
        if g is None:
            return np.nan
        if kind == 't':
            s = g['TAVG'].dropna()
            if len(s) == 0:
                t = (g['TMIN'] + g['TMAX']) / 2.0
                s = t.dropna()
            return float(np.nanmean(s)) if len(s) else np.nan
        if kind == 'pr':
            s = g['PRCP'].dropna()
            if len(s) == 0:
                return np.nan
            yearly = s.groupby(g.loc[s.index, 'DATE'].dt.year).sum()
            return float(yearly.mean() * 25.4)   # PRCP in inches -> mm
        if kind == 'wind':
            s = g['AWND'].dropna()
            return float(np.nanmean(s) * 0.44704) if len(s) else np.nan
    except Exception as e:
        print(f"  WARN obs NOAA {sid} {kind}: {e}")
    return np.nan

def _bch_obs_annual(sid, kind):
    """BCH daily obs, 1986-2005 annual mean.

    Two layouts are supported:
      * Eva's per-station:  BCH_tn<sid>.csv / BCH_tx<sid>.csv / BCH_p24<sid>.csv
        (single station column, date index)
      * wide multi-station: BCH_tnId.csv / BCH_txId.csv / BCH_p24Id.csv
        (one column per station, date index)  [current data/bch_obs layout]
    """
    if kind == 't':
        pmin, pmax = f'BCH_tn{sid}.csv', f'BCH_tx{sid}.csv'
    elif kind == 'pr':
        pmin = pmax = f'BCH_p24{sid}.csv'
    else:
        return np.nan   # wind not in BCH daily files
    pmin, pmax = os.path.join(BCH_OBS_DIR, pmin), os.path.join(BCH_OBS_DIR, pmax)
    if not os.path.exists(pmin):
        # fall back to the wide multi-station files
        wide_min = 'BCH_tnId.csv' if kind == 't' else 'BCH_p24Id.csv'
        wide_max = 'BCH_txId.csv' if kind == 't' else 'BCH_p24Id.csv'
        pmin = os.path.join(BCH_OBS_DIR, wide_min)
        pmax = os.path.join(BCH_OBS_DIR, wide_max)
        if not os.path.exists(pmin):
            return np.nan
    # 5 stations (BLN, DLU, MIS, NTY, WOL) record p24 in tenths of mm in the
    # wide CSVs (verified: their /10 annual sums are climatologically sensible
    # while the raw values are ~10x every other station).
    TENTHS_MM = {'BLN', 'DLU', 'MIS', 'NTY', 'WOL'}
    try:
        dmin = pd.read_csv(pmin, index_col=0)
        dmin.index = pd.to_datetime(dmin.index).tz_localize(None)
        smin = dmin[sid].dropna()
        if kind == 'pr':
            if sid in TENTHS_MM:
                smin = smin / 10.0
            yearly = smin.groupby(smin.index.year).sum()
            return float(yearly.mean())
        dmax = pd.read_csv(pmax, index_col=0)
        dmax.index = pd.to_datetime(dmax.index).tz_localize(None)
        smax = dmax[sid].dropna()
        t = (smin + smax) / 2.0
        # annual mean of daily means over 1986-2005
        t = t[(t.index >= pd.Timestamp(1986, 1, 1)) & (t.index < pd.Timestamp(2006, 1, 1))]
        return float(t.groupby(t.index.year).mean().mean())
    except Exception as e:
        print(f"  WARN obs BCH {sid} {kind}: {e}")
        return np.nan

def obs_annual(sid, agency, kind):
    """agency: 'E' (ECCC), 'N' (NOAA), 'B' (BCH). Returns 20-yr annual value."""
    if agency == 'E':
        return _eccc_elev_obs(sid, kind)
    if agency == 'N':
        return _noaa_obs_annual(sid, kind)
    return _bch_obs_annual(sid, kind)

# Precompute obs for all stations (once)
print("Computing obs annual means (1986-2005)...")
obs_t, obs_pr, obs_w = {}, {}, {}
for i, sid in enumerate(all_ids):
    agency = 'E' if sid.isdigit() else ('N' if sid.startswith('USC') else
                                        'B' if not sid.startswith('US') else 'N')
    # ECCC ids are numeric; BCH ids are 3 letters; NOAA land are USCxxxx
    if sid.isdigit():
        agency = 'E'
    elif sid.startswith('USC'):
        agency = 'N'
    else:
        agency = 'B'  # BCH (or buoy handled separately)
    if sid in buoy_ids:
        continue  # no land obs txt for buoys
    obs_t[sid]  = obs_annual(sid, agency, 't')
    obs_pr[sid] = obs_annual(sid, agency, 'pr')
    obs_w[sid]  = obs_annual(sid, agency, 'wind')
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(all_ids)}")

# ---------------------------------------------------------------------------
# Assemble tables per dataset
# ---------------------------------------------------------------------------
# WRF station file availability
print("Finding WRF per-station files...")
wrf_files = {}
for dom in ['d01', 'd02', 'd03']:
    for var in ['t', 'pr', 'wind']:
        wrf_files[(dom, var)] = wrf_st_files(var, dom)

# Gridded CanESM2 / CanRCM4 values (per station, per var)
print("Sampling gridded CanESM2 / CanRCM4 at station locations...")
grid_vals = {}
# CanESM2: per-station CDO daily (M9, station_cdo/) -> 1986-2005 annual value.
# Rewired from M3 nearest-point (data/canesm2_raw, +3.4 C bias source) to M9
# on 2026-09-21: all 177 stations now have CDO bicubic per-station files.
for kind, vname in [('t', 'tas'), ('pr', 'pr'), ('wind', 'sfcWind')]:
    print(f"  CanESM2 {kind} (CDO per-station M9)...")
    try:
        grid_vals[('CanESM2', kind)] = np.array([
            sp.cdo_station_value(sid, kind, DATA, 'canesm2') for sid in land_ids])
    except Exception as e:
        print(f"  WARN CanESM2 {kind}: {e}")
        grid_vals[('CanESM2', kind)] = None
# CanRCM4: per-station CDO daily (M10, station_cdo/) -> 1986-2005 annual value.
# Aligned to land_ids. ECCC/BCH stations are not in the 113-NOAA common set
# (no M10 file) -> NaN (previously M6 bilinear stand-in, now deleted).
for kind, vname in [('t', 'tas'), ('pr', 'pr'), ('wind', 'sfcWind')]:
    print(f"  CanRCM4 {kind} (CDO per-station M10)...")
    try:
        grid_vals[('CanRCM4', kind)] = np.array([
            sp.cdo_station_value(sid, kind, DATA, 'canrcm4') for sid in land_ids])
    except Exception as e:
        print(f"  WARN CanRCM4 {kind}: {e}")
        grid_vals[('CanRCM4', kind)] = None

tables = {}
for name, dom in [('D03', 'd03'), ('D02', 'd02'), ('D01', 'd01')]:
    print(f"WRF {dom} t: computing...")
    t_b, pr_b, w_b = {}, {}, {}
    # temperature
    for sid in land_ids:
        key = ('E' if sid.isdigit() else ('N' if sid.startswith('USC') else 'B')) + sid
        p = wrf_files[(dom, 't')].get(key)
        if p is None:
            continue
        try:
            m = annual_mean(wrf_station_series(p, 'T2', 't'), 't')
            o = obs_t.get(sid)
            if not np.isnan(o):
                t_b[sid] = m - o
        except Exception as e:
            print(f"  WARN {name} t {sid}: {e}")
    for sid in buoy_ids:
        # buoy obs: not in obs dirs -> obs-based bias stays NaN
        pass
    tables[name] = pd.DataFrame({
        'station': all_ids,
        'lat': lats_s.values,
        'lon': lons_s.values,
        'obs_elev_m': elev_s.values,
        'elevation_bias_m': elev_bias[name].values,
        'elev_method': elev_method[name].values,
        'temperature_bias_C': [t_b.get(s, np.nan) for s in all_ids],
        'precipitation_bias_pct': [pr_b.get(s, np.nan) for s in all_ids],
        'wind_bias_pct': [w_b.get(s, np.nan) for s in all_ids],
    })

# (D02/D03 pr & wind, D01 all vars, and raw/RCM rows)
for name, dom in [('D03', 'd03'), ('D02', 'd02'), ('D01', 'd01')]:
    tdf = tables[name]
    pr_b, w_b = {}, {}
    for sid in land_ids:
        agency = 'E' if sid.isdigit() else ('N' if sid.startswith('USC') else 'B')
        key = agency + sid
        # pr
        p = wrf_files[(dom, 'pr')].get(key)
        if p is not None:
            try:
                m = annual_mean(wrf_station_series(p, 'pr', 'pr'), 'pr')
                o = obs_pr.get(sid)
                if not np.isnan(o) and o > 0:
                    pr_b[sid] = 100.0 * (m - o) / o
            except Exception as e:
                print(f"  WARN {name} pr {sid}: {e}")
        # wind
        p = wrf_files[(dom, 'wind')].get(key)
        if p is not None:
            try:
                m = annual_mean(wrf_station_series(p, 'wspd', 'wind'), 'wind')
                o = obs_w.get(sid)
                if not np.isnan(o) and o > 0:
                    w_b[sid] = 100.0 * (m - o) / o
            except Exception as e:
                print(f"  WARN {name} wind {sid}: {e}")
    tdf['precipitation_bias_pct'] = [pr_b.get(s, np.nan) for s in all_ids]
    tdf['wind_bias_pct'] = [w_b.get(s, np.nan) for s in all_ids]
    tables[name] = tdf
    print(f"WRF {dom} done: t={sum(1 for v in tdf['temperature_bias_C'] if not np.isnan(v))}, "
          f"pr={tdf['precipitation_bias_pct'].notna().sum()}, "
          f"wind={tdf['wind_bias_pct'].notna().sum()}")

# CanESM2 & CanRCM4 (gridded, nearest point, land stations only)
for src, kind_map in [('CanESM2', {'t': 't', 'pr': 'pr', 'wind': 'wind'}),
                      ('CanRCM4', {'t': 't', 'pr': 'pr', 'wind': 'wind'})]:
    print(f"{src}: assembling table...")
    t_b, pr_b, w_b = {}, {}, {}
    for j, sid in enumerate(land_ids):
        o_t  = obs_t.get(sid)
        o_pr = obs_pr.get(sid)
        o_w  = obs_w.get(sid)
        m_t  = grid_vals[(src, 't')][j] if grid_vals[(src, 't')] is not None else np.nan
        m_pr = grid_vals[(src, 'pr')][j] if grid_vals[(src, 'pr')] is not None else np.nan
        m_w  = grid_vals[(src, 'wind')][j] if grid_vals[(src, 'wind')] is not None else np.nan
        if not np.isnan(m_t) and not np.isnan(o_t):
            t_b[sid] = m_t - o_t
        if not np.isnan(m_pr) and not np.isnan(o_pr) and o_pr > 0:
            pr_b[sid] = 100.0 * (m_pr - o_pr) / o_pr
        if not np.isnan(m_w) and not np.isnan(o_w) and o_w > 0:
            w_b[sid] = 100.0 * (m_w - o_w) / o_w
    tables[src] = pd.DataFrame({
        'station': all_ids,
        'lat': lats_s.values,
        'lon': lons_s.values,
        'obs_elev_m': elev_s.values,
        'elevation_bias_m': elev_bias[src].values,
        'elev_method': elev_method[src].values,
        'temperature_bias_C': [t_b.get(s, np.nan) for s in all_ids],
        'precipitation_bias_pct': [pr_b.get(s, np.nan) for s in all_ids],
        'wind_bias_pct': [w_b.get(s, np.nan) for s in all_ids],
    })

# ---------------------------------------------------------------------------
# Save & report
# ---------------------------------------------------------------------------
with open(OUT_PKL, 'wb') as f:
    pd.to_pickle(tables, f)
print(f"\nSaved {OUT_PKL}")
for name, df in tables.items():
    print(f"\n=== {name} ===")
    print(df[['temperature_bias_C', 'precipitation_bias_pct', 'wind_bias_pct']]
          .notna().sum().to_string())
    print(df[['elevation_bias_m', 'temperature_bias_C', 'precipitation_bias_pct',
              'wind_bias_pct']].mean(numeric_only=True).round(2).to_string())
