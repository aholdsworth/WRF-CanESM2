#!/usr/bin/env python
"""
station_geo.py -- side-effect-free station-metadata / topography / model-elevation
loaders, extracted from plot_elev_bias_revised.py (2026-09-29).

This is a LEAF module: importing it runs NO side effects (no os.makedirs, no
print, no WRFDomainLib import, no station assembly).  It is the single source of
truth for:
  * per-agency station metadata loaders (load_eccc/bch/noaa/noaa_land/buoys)
  * domain topography readers (read_geo_em, read_orog)
  * the model-elevation machinery used by both the elevation-bias figure
    (plot_elev_bias_revised.py) and the root table producer
    (make_station_tables.py): wrf_t_files, _weight_file_from_history,
    _weighted_elev, weighted_model_elevs, get_nearest_gridpoint, modeled_elevs,
    _can_weight_file.

Consumers:
  * plot_elev_bias_revised.py   (elevation-bias figure)
  * make_station_tables.py      (root producer of station_tables.pkl)
  * make_elev_bias_grid.py      (read_geo_em, read_orog)
  * extract_canrcm4_stations.py (load_eccc/bch/noaa_land/buoys)

Note: modeled_elevs is stateless here (takes the elevation Series explicitly)
rather than reading plot_elev_bias_revised's module-level `elev`, so it can be
shared by the root producer.
"""

import os
import re
import json
import sys

import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
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


BASE        = os.environ.get('WRF_CANESM2_BASE') or _workroot(os.path.dirname(os.path.abspath(__file__)))
STATION_DIR = os.path.join(BASE, 'files')            # unzipped station lists
DOMAIN_DIR  = os.path.join(BASE, 'plotting_files')   # geo_em / orog / namelist

# Per-station WRF CDO-extracted files (t_d0X_*.nc); their `history` attribute
# records the exact CDO bicubic weight file used to remap each station.
WRF_ST_DIR  = os.path.join(BASE, 'data', 'wrf_stations')

# One obsolete weight-file directory: the 2022 D03 evg000 set was later merged
# into the sibling CanESM2_weight_files/ dir, so history paths pointing at the
# "_historical" dir are rewritten to the merged location.
_WGT_HIST = '/gpfs/fs7/dfo/hpcmc/pfm/evg000/verification/CanESM2_weight_files_historical/'
_WGT_NEW  = '/gpfs/fs7/dfo/hpcmc/pfm/evg000/verification/CanESM2_weight_files/'

# Can* (M9 CanESM2 / M10 CanRCM4) station CDO bicubic weight sets.
#   * pfm/.../CanESM2-WRF/weight_files/{CanESM2_raw,CanRCM4}/weights_<sid>.nc
#     are the exact files our M9/M10 extractions used (recorded in each
#     output's `history` attribute); NOAA stations are tagged by GHCND id
#     (USC00...).  They are read-only.
#   * weight_files_gen/{CanESM2_raw,CanRCM4}/weights_<sid>.nc are our
#     regenerated set (ECCC sid-tagged, BCH code-tagged); the Beaverton gate
#     confirmed them bit-identical to the pfm set where both exist.
_PFG_CANW = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/weight_files'
_GEN_CANW = os.path.join(BASE, 'weight_files_gen')

# ECCC pipeline sid (WMO) -> ClimateID, from the M9/M10 extraction map.
# Loaded lazily (once) on first use by _can_weight_file so that importing this
# module has no side effects and the JSON read is not forced on consumers that
# never touch Can* elevation.
_ECCC_SID2CID = None


def _eccc_sid2cid():
    global _ECCC_SID2CID
    if _ECCC_SID2CID is None:
        _map = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'station_cdo', 'eccc_bch_desc_map.json')
        if not os.path.exists(_map):
            _map = os.path.join(BASE, 'station_cdo', 'eccc_bch_desc_map.json')
        with open(_map) as _f:
            _ECCC_SID2CID = json.load(_f)['ECCC']
    return _ECCC_SID2CID


# ----------------------------------------------------------------------------
# Station metadata
# ----------------------------------------------------------------------------
def load_eccc():
    df = pd.read_csv(os.path.join(STATION_DIR, 'ECCC_d03_stations.csv'), header=None)
    ids  = list(df.iloc[:, 4])
    # 951 (HOPE SLIDE, BC) is now INCLUDED (2026-09-21): it has CDO-extracted
    # per-station model data (M9/M10, pipeline sid 1113581) and is in the
    # 177-station inventory.  The old EXCLUDE={'951'} (consistency to-do,
    # 2026-09-17, DATA_README.md Known issues #4) is removed.
    lats = df.iloc[:, 7].copy(); lats.index = [str(i) for i in df.iloc[:, 4]]
    lons = df.iloc[:, 8].copy(); lons.index = [str(i) for i in df.iloc[:, 4]]
    elev = df.iloc[:, 11].copy(); elev.index = [str(i) for i in df.iloc[:, 4]]
    lats = lats.loc[[str(i) for i in ids]]
    lons = lons.loc[[str(i) for i in ids]]
    elev = elev.loc[[str(i) for i in ids]]
    return ids, lats, lons, elev


def load_bch():
    df = pd.read_csv(os.path.join(STATION_DIR, 'BCH_d03_stations.csv'))
    ids  = list(df["STATION_NO"])
    lats = df['Y'].copy(); lats.index = ids
    lons = df['X'].copy(); lons.index = ids
    elev = df["ELEV"].copy(); elev.index = ids
    return ids, lats, lons, elev


def load_noaa(fname):
    df = pd.read_csv(os.path.join(STATION_DIR, fname))
    ids  = list(df.iloc[:, 0])
    lats = df.iloc[:, 2].copy(); lats.index = ids
    lons = df.iloc[:, 3].copy(); lons.index = ids
    elev = df["ELEVATION"].copy(); elev.index = ids
    return ids, lats, lons, elev


def load_noaa_land(fname):
    """Land-only NOAA station list (excludes buoys)."""
    df = pd.read_csv(os.path.join(STATION_DIR, fname))
    df = df[~df['STATION'].str.startswith('USW')]
    ids  = list(df['STATION'])
    lats = df['LATITUDE'].copy(); lats.index = ids
    lons = df['LONGITUDE'].copy(); lons.index = ids
    elev = df['ELEVATION'].copy(); elev.index = ids
    return ids, lats, lons, elev


def load_buoys(fname):
    df = pd.read_csv(os.path.join(STATION_DIR, fname))
    ids  = list(df["STATION_ID"])
    lats = df['Y'].copy(); lats.index = ids
    lons = df['X'].copy(); lons.index = ids
    elev = pd.Series(0.0, index=ids)   # buoys -> sea level
    return ids, lats, lons, elev


# ----------------------------------------------------------------------------
# Domain topography
# ----------------------------------------------------------------------------
def read_geo_em(fname):
    nc = Dataset(os.path.join(DOMAIN_DIR, fname), mode='r')
    lat  = np.squeeze(nc.variables['XLAT_C'][:])
    lon  = np.squeeze(nc.variables['XLONG_C'][:])
    topo = np.squeeze(nc.variables['HGT_M'][:])
    nc.close()
    return lat, lon, topo


def read_orog(fname, lon_subtract=360.0):
    """orog_CanESM2 has 1D lat/lon; orog_CanRCM4 has 2D lat/lon grids."""
    nc = Dataset(os.path.join(DOMAIN_DIR, fname), mode='r')
    lat  = np.squeeze(nc.variables['lat'][:])
    lon  = np.squeeze(nc.variables['lon'][:])
    topo = np.squeeze(nc.variables['orog'][:])
    nc.close()
    if lon.ndim == 1:
        lon = lon - lon_subtract
        lons2d, lats2d = np.meshgrid(lon, lat)
    else:
        lon  = lon - lon_subtract
        lons2d, lats2d = lon, lat
    return lats2d, lons2d, topo


# ----------------------------------------------------------------------------
# WRF per-station T-file -> CDO weight file -> model elevation
# ----------------------------------------------------------------------------
def wrf_t_files(dom):
    """Map station id -> path of the per-station WRF temperature file for a
    domain.  The T file's `history` attribute records the exact CDO bicubic
    weight file that was used to remap that station onto the model grid, so it
    is the provenance record for the model elevation.  Naming:
      d01:  {t}_d01_{ECCC|NOAA|BCH}_{ID}.nc  (BCH d01 also untagged 3-letter)
      d02/d03: {t}_{dom}_st{ID}.nc (ID = numeric ECCC or USxxx NOAA) + BCH_{ID}
    """
    out = {}
    d = WRF_ST_DIR
    BCH3 = {'ALU','BCK','BLN','CLO','CMU','CMX','COQ','CQM','CRU','DAI','DLU','ECL',
            'GOC','MIS','NTY','STA','WAH','WOL'}
    for f in os.listdir(d):
        if not f.endswith('.nc') or not f.startswith('t_' + dom + '_'):
            continue
        rest = f[len('t_' + dom + '_'):-3]
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
                print(f'WARN wrf_t_files: untagged d01 file skipped: {f}', file=sys.stderr)
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
                print(f'WARN wrf_t_files: untagged {dom} file skipped: {f}', file=sys.stderr)
    return out


_WEIGHT_RE = re.compile(r'remap,[^,]+,(\S+?\.nc)')


def _weight_file_from_history(path):
    """Return the CDO remap weight-file path recorded in a T file's history
    attribute, or None if not found.  Rewrites the obsolete 2022 D03
    "_historical" weight dir to its merged location."""
    try:
        nc = Dataset(path, mode='r')
        h = str(nc.getncattr('history'))
        nc.close()
    except Exception:
        return None
    m = _WEIGHT_RE.search(h)
    if not m:
        return None
    w = m.group(1)
    if _WGT_HIST in w:
        w = w.replace(_WGT_HIST, _WGT_NEW)
    return w


def _weighted_elev(wfile, hgt):
    """Bicubic-weighted model elevation for one station from a CDO weight file
    and the domain A-grid HGT_M.  `src_address` flat index -> (i, j) via
    divmod(flat, nx).  Ocean-masked src cells (HGT_M NaN) drop out of the
    sum.  Returns (elev, weight_sum) or None on any failure."""
    ny, nx = hgt.shape
    try:
        w = Dataset(wfile, mode='r')
        sa = [int(x) for x in w['src_address'][:]]
        rm = np.asarray(w['remap_matrix'][:], dtype=float)
        w.close()
    except Exception:
        return None
    if sa:
        if any(f < 0 or f >= ny * nx for f in sa):
            return None
    vals = np.array([hgt[divmod(f, nx)] for f in sa])
    if len(sa) == 0:
        return None
    if np.isnan(vals).any():
        # masked (ocean) source cell: drop that link from the sum
        good = ~np.isnan(vals)
        if not good.any():
            return None
        vals = vals[good]
        rm = rm[good]
    # The weight file's destination is a single station point
    # (dst_grid_dims=[1,1]); CDO pads the remap_matrix to (nsrc, 4) where
    # column 0 is the station and columns 1-3 are zero-weight placeholders.
    # Station value = sum_k rm[k, 0] * HGT[k]  (the station column).
    elev = float((rm[:, 0] * vals).sum())
    return elev, float(rm[:, 0].sum())


def weighted_model_elevs(dom, ids):
    """Per-station model elevation for a WRF domain using the *exact* CDO
    bicubic weights that were used to extract each station's T series.  This is
    the orography actually baked into the model data we compare, so it is the
    consistent reference for the elevation-bias diagnostic.

    Returns (series, warnings): `series` is a pd.Series indexed by station id
    (model elevation in m, NaN where unresolved); `warnings` is a list of
    human-readable strings for every station that could not be resolved via a
    weight file (caller should fall back to the corner lookup).
    """
    if dom in ('CanESM2', 'CanRCM4'):
        _, _, hgt = read_orog(f'orog_{dom}.nc')
        files = None  # weight file resolved per-station by _can_weight_file()
    else:
        n = int(dom[2])  # 'd01' -> 1
        _, _, hgt = read_geo_em(f'geo_em.d0{n}.nc')
        files = wrf_t_files(dom)
    out = {}
    warns = []
    for sid in ids:
        # Map the raw station id to the wrf_t_files key namespace
        # (E-prefixed ECCC / N-prefixed NOAA / B-prefixed BCH):
        #   numeric ECCC id             -> 'E'+id
        #   'USC.../USW...' NOAA id     -> 'N'+id
        #   'B'+buoy (ECCC/NOAA buoy)   -> no T file; NaN quietly (the real
        #                                  corner fallback is done by caller)
        #   bare 3-letter BCH id        -> 'B'+id
        #   already 'B'+code (module)   -> as-is
        if sid.startswith('B') and not sid.isdigit():
            k, is_buoy = sid, True
        elif sid.startswith('US'):
            k, is_buoy = 'N' + sid, False
        elif sid.isdigit():
            k, is_buoy = 'E' + sid, False
        else:
            k, is_buoy = 'B' + sid, False
        p = _can_weight_file(dom, sid) if files is None else files.get(k)
        if p is None:
            if files is None or not is_buoy:
                if files is None:
                    warns.append(f'{dom} {sid}: no weight file -> nearest-GP fallback')
                else:
                    warns.append(f'{dom} {sid}: no T file (t_{dom}_*) -> corner fallback')
            out[sid] = np.nan
            continue
        wf = p if files is None else _weight_file_from_history(p)
        if wf is None:
            warns.append(f'{dom} {sid}: no remap weight in history -> corner fallback')
            out[sid] = np.nan
            continue
        if not os.path.exists(wf):
            warns.append(f'{dom} {sid}: weight file missing {wf} -> corner fallback')
            out[sid] = np.nan
            continue
        res = _weighted_elev(wf, hgt)
        if res is None:
            warns.append(f'{dom} {sid}: weight eval failed ({wf}) -> corner fallback')
            out[sid] = np.nan
            continue
        out[sid] = res[0]
    idx = [str(s) for s in ids]
    series = pd.Series([out.get(s, np.nan) for s in idx], index=idx, name='model_elev')
    return series, warns


# ----------------------------------------------------------------------------
# Nearest-grid-point topography lookup
# ----------------------------------------------------------------------------
def get_nearest_gridpoint(target_lat, target_lon, lat, lon, topo):
    """Nearest grid point by Chebyshev (max) distance in lat/lon space.
    lat/lon are 2D grids matching topo's (lat, lon) index order."""
    abslat = np.abs(lat - target_lat)
    abslon = np.abs(lon - target_lon)
    d = np.maximum(abslon, abslat)
    x, y = np.where(d == d.min())
    return topo[x[0], y[0]]


def modeled_elevs(lat, lon, topo, lats, lons, elev):
    """Nearest-grid-point model elevation for each station in `elev`.

    `lats`/`lons`/`elev` are pd.Series indexed by station id (same index, same
    order).  For grid point i the station `elev.index[i]` at (lats.iloc[i],
    lons.iloc[i]) is looked up on the (lat, lon, topo) grid.  Returns a
    DataFrame with a single column 'elev' indexed by `elev.index`.

    Stateless (takes the station lat/lon/elev explicitly rather than reading a
    module-level list) so it can be shared by the root producer.
    """
    out = []
    for i in range(len(elev)):
        out.append(get_nearest_gridpoint(lats.iloc[i], lons.iloc[i], lat, lon, topo))
    return pd.DataFrame(out, columns=['elev'], index=elev.index)


def _can_weight_file(src, sid):
    """Weight file for one Can* station, or None if it does not exist."""
    if sid.startswith('US'):
        p = os.path.join(_PFG_CANW, src if src != 'CanESM2' else 'CanESM2_raw', f'weights_{sid}.nc')
        if os.path.exists(p):
            return p
        return None  # pfm set is NOAA-only; no local NOAA weights were generated
    if sid.isdigit():
        # table id is the ClimateID; local weights are pipeline-sid tagged
        wmo = [k for k, v in _eccc_sid2cid().items() if v == sid]
        if not wmo:
            return None
        p = os.path.join(_GEN_CANW, 'CanESM2_raw' if src == 'CanESM2' else 'CanRCM4', f'weights_{wmo[0]}.nc')
        return p if os.path.exists(p) else None
    # BCH 3-letter code (tagged identically in the local gen set)
    p = os.path.join(_GEN_CANW, 'CanESM2_raw' if src == 'CanESM2' else 'CanRCM4', f'weights_{sid}.nc')
    return p if os.path.exists(p) else None
