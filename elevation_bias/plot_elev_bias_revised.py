#!/usr/bin/env python
"""
Revised CanESM2-WRF elevation-bias figure script.

Derived from Eva's plot_elev_error.py / plotting_elev_bias.py.

Changes from the original:
  * Self-contained: no dependency on /Users/evagnegy/... (Eva's Mac) paths,
    canesm2_eval_funcs, or make_colorbar.
  * All paths point at this sandbox workspace (plotting_files/ + files/).
  * Reproducible at 11 pt fonts (FONT = 11).
  * Produces the observed-station-elevation map plus the model-minus-obs
    elevation-bias maps for D03 / D02 / D01 / CanESM2 / CanRCM4.

Run with the py_2024 environment:
  /home/amh001/space_fs7/software_2022/python/py_2024/bin/python plot_elev_bias_revised.py
"""

import os
import sys
import re
import json
import warnings

import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cf

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ----------------------------------------------------------------------------
# Paths & fonts
# ----------------------------------------------------------------------------
# Root of the analysis working tree (data/, pickles/, plotting_files/, files/).
# Defaults to the HPC working tree used for the analysis; override with
# ELEV_BIAS_BASE if you relocate the tree.
BASE   = os.environ.get('ELEV_BIAS_BASE',
    '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS')
STATION_DIR = os.path.join(BASE, 'files')          # unzipped station lists
DOMAIN_DIR  = os.path.join(BASE, 'plotting_files') # geo_em / orog / namelist
REPO_DIR    = os.path.join(BASE, 'WRF-CanESM2')    # WRFDomainLib
OUT_DIR     = os.path.join(BASE, 'figures', 'elevation_bias')
os.makedirs(OUT_DIR, exist_ok=True)

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
with open(os.path.join(BASE, 'station_cdo', 'eccc_bch_desc_map.json')) as _f:
    _ECCC_SID2CID = json.load(_f)['ECCC']

sys.path.insert(0, REPO_DIR)
import WRFDomainLib

FONT = 11  # requested font size (points)

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

eccc_ids, eccc_lats, eccc_lons, eccc_elev   = load_eccc()
bch_ids,  bch_lats,  bch_lons,  bch_elev    = load_bch()
noaa_ids, noaa_lats, noaa_lons, noaa_elev   = load_noaa('NOAA_d03_stations.csv')
noaa_buoy_ids, noaa_buoy_lats, noaa_buoy_lons, noaa_buoy_elev = load_buoys('NOAA_buoys.csv')
eccc_buoy_ids, eccc_buoy_lats, eccc_buoy_lons, eccc_buoy_elev = load_buoys('ECCC_buoys.csv')

# Combine all stations (same station set the original used for elevation bias).
lats = pd.concat([eccc_lats, bch_lats, noaa_lats, eccc_buoy_lats, noaa_buoy_lats])
lats.index = lats.index.astype(str)
lats = lats.sort_index()

lons = pd.concat([eccc_lons, bch_lons, noaa_lons, eccc_buoy_lons, noaa_buoy_lons])
lons.index = lons.index.astype(str)
lons = lons.sort_index()

elev = pd.concat([eccc_elev, bch_elev, noaa_elev, eccc_buoy_elev, noaa_buoy_elev])
elev.index = elev.index.astype(str)
elev = elev.sort_index()

# Land-only elevation diagnostic: the 11 buoys (7 ECCC + 4 NOAA) have no
# meaningful orography (obs_elev=0, and CanRCM4/CanESM2 orog marks their
# coastal cells as ocean=0 m), so they were being stacked at bias=0 -> the
# vertical-line artefact in the Can* rows.  Excluded from the elevation calc.
buoy_all_ids = {str(i) for i in list(eccc_buoy_ids) + list(noaa_buoy_ids)}
_lats_dropped = len(lats.index.intersection(buoy_all_ids))
lats = lats.drop(lats.index.intersection(buoy_all_ids))
lons = lons.drop(lons.index.intersection(buoy_all_ids))
elev = elev.drop(elev.index.intersection(buoy_all_ids))
print(f"Elevation diagnostic: {len(elev)} land stations (dropped {_lats_dropped} buoys)")

print(f"Loaded {len(elev)} land stations "
      f"(ECCC={len(eccc_ids)}, BCH={len(bch_ids)}, NOAA={len(noaa_ids)}, "
      f"ECCC-buoys={len(eccc_buoy_ids)}, NOAA-buoys={len(noaa_buoy_ids)})")

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

lat_d03, lon_d03, topo_d03 = read_geo_em('geo_em.d03.nc')
lat_d02, lon_d02, topo_d02 = read_geo_em('geo_em.d02.nc')
lat_d01, lon_d01, topo_d01 = read_geo_em('geo_em.d01.nc')
lats_canesm2, lons_canesm2, topo_canesm2 = read_orog('orog_CanESM2.nc')
lats_canrcm4, lons_canrcm4, topo_canrcm4 = read_orog('orog_CanRCM4.nc')

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

def modeled_elevs(lat, lon, topo):
    out = []
    for i in range(len(elev)):
        out.append(get_nearest_gridpoint(lats.iloc[i], lons.iloc[i], lat, lon, topo))
    return pd.DataFrame(out, columns=['elev'], index=elev.index)

print("Computing model elevation for each station...")
# WRF domains: use the extraction-consistent CDO bicubic weights (the orography
# actually baked into the model data we compare).  Falls back to the corner
# nearest-grid-point lookup for any station without a usable T-file weight file
# (e.g. buoys without T files).
def _wrf_model_elevs(dom, lat, lon, topo):
    s, warns = weighted_model_elevs(dom, list(elev.index))
    corner = modeled_elevs(lat, lon, topo)
    for sid in s.index:
        if np.isnan(s.loc[sid]):
            s.loc[sid] = corner.loc[sid, 'elev']
            print(f"  WARN {dom} {sid}: weight lookup fell back to corner "
                  f"(model_elev={s.loc[sid]:.1f})")
    return s

# Can* (M9 CanESM2 / M10 CanRCM4): same extraction-consistent method as WRF,
# using the exact CDO bicubic weight files each station's M9/M10 extraction was
# remapped with (recorded in the output `history`), applied to the same orog
# file the model data was remapped from:
#   NOAA  -> pfm weight_files/<set>/weights_<GHCND>.nc   (read-only, used as-is)
#   ECCC  -> weight_files_gen/<set>/weights_<pipeline sid>.nc (bit-validated)
#   BCH   -> weight_files_gen/<set>/weights_<code>.nc
# Stations with no weight file anywhere fall back to the nearest-grid-point
# lookup over orog (logged, and reported to the caller for an `elev_method`
# flag).  This is stricter than Eva's plot_elev_error.py, which uses pure
# nearest-GP for all rows.
def _can_weight_file(src, sid):
    """Weight file for one Can* station, or None if it does not exist."""
    if sid.startswith('US'):
        p = os.path.join(_PFG_CANW, src if src != 'CanESM2' else 'CanESM2_raw', f'weights_{sid}.nc')
        if os.path.exists(p):
            return p
        return None  # pfm set is NOAA-only; no local NOAA weights were generated
    if sid.isdigit():
        # table id is the ClimateID; local weights are pipeline-sid tagged
        wmo = [k for k, v in _ECCC_SID2CID.items() if v == sid]
        if not wmo:
            return None
        p = os.path.join(_GEN_CANW, 'CanESM2_raw' if src == 'CanESM2' else 'CanRCM4', f'weights_{wmo[0]}.nc')
        return p if os.path.exists(p) else None
    # BCH 3-letter code (tagged identically in the local gen set)
    p = os.path.join(_GEN_CANW, 'CanESM2_raw' if src == 'CanESM2' else 'CanRCM4', f'weights_{sid}.nc')
    return p if os.path.exists(p) else None

def _can_model_elevs(src, lat2d, lon2d, topo2d):
    """Per-station Can* model elevation from the extraction's own bicubic
    weights, with nearest-GP fallback for stations without a weight file.
    Returns (series, fallback_ids)."""
    s, warns = weighted_model_elevs(src, list(elev.index))
    nearest = modeled_elevs(lat2d, lon2d, topo2d)
    fallback = []
    for sid in s.index:
        if np.isnan(s.loc[sid]):
            s.loc[sid] = nearest.loc[sid, 'elev']
            fallback.append(sid)
            print(f"  WARN {src} {sid}: no weight file -> nearest-GP "
                  f"(model_elev={s.loc[sid]:.1f})")
    return s, fallback

me_d03 = _wrf_model_elevs('d03', lat_d03, lon_d03, topo_d03)
me_d02 = _wrf_model_elevs('d02', lat_d02, lon_d02, topo_d02)
me_d01 = _wrf_model_elevs('d01', lat_d01, lon_d01, topo_d01)
me_raw, raw_fallback = _can_model_elevs('CanESM2', lats_canesm2, lons_canesm2, topo_canesm2)
me_rcm, rcm_fallback = _can_model_elevs('CanRCM4', lats_canrcm4, lons_canrcm4, topo_canrcm4)

elev_bias_d03 = me_d03.reindex(elev.index) - elev
elev_bias_d02 = me_d02.reindex(elev.index) - elev
elev_bias_d01 = me_d01.reindex(elev.index) - elev
elev_bias_raw = me_raw.reindex(elev.index) - elev
elev_bias_rcm = me_rcm.reindex(elev.index) - elev

print("Elevation bias summary (model - obs, m):")
for name, s in [('D03', elev_bias_d03), ('D02', elev_bias_d02),
                ('D01', elev_bias_d01), ('CanESM2', elev_bias_raw),
                ('CanRCM4', elev_bias_rcm)]:
    print(f"  {name:10s} mean={s.mean():8.1f}  median={s.median():8.1f}  "
          f"RMS={np.sqrt((s**2).mean()):8.1f}  min={s.min():8.1f}  max={s.max():8.1f}")

# ----------------------------------------------------------------------------
def _make_figures():
    WPS_FILE = os.path.join(DOMAIN_DIR, 'namelist.wps.txt')
    # Map base (D03 Lambert domain, as in original plot_all_d03)
    wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = \
        WRFDomainLib.calc_wps_domain_info(WPS_FILE)
    corner_x3, corner_y3 = WRFDomainLib.reproject_corners(
        corner_lon_full[2, :], corner_lat_full[2, :], wpsproj, latlonproj)

    def new_domain_ax():
        fig = plt.figure(figsize=(8, 7), dpi=200)
        ax = fig.add_subplot(1, 1, 1, projection=wpsproj)

        cmap, vmin, vmax = 'terrain', 0, 3000
        ax.pcolormesh(lon_d02, lat_d02, topo_d02, cmap=cmap, vmin=vmin, vmax=vmax,
                      alpha=0.3, transform=ccrs.PlateCarree(), zorder=0)
        ax.pcolormesh(lon_d03, lat_d03, topo_d03, cmap=cmap, vmin=vmin, vmax=vmax,
                      alpha=0.4, transform=ccrs.PlateCarree(), zorder=0)

        # NOTE: cf.OCEAN / cf.BORDERS / cf.STATES require Natural Earth shapefiles
        # that are not present and cannot be downloaded (no network egress). The
        # D02+D03 topography pcolormesh already renders ocean as low terrain, so the
        # map remains complete without these cosmetic features.

        # D03 box
        rxf = corner_x3[0] / 65
        ryf = -corner_y3[0] / 13
        ax.add_patch(mpl.patches.Rectangle(
            (corner_x3[0] + rxf, corner_y3[0] + ryf),
            length_x[2], length_y[2],
            fill=None, lw=2, edgecolor='red', zorder=2))
        ax.text(-3700000, 700000, 'D03', va='top', ha='left',
                fontweight='bold', size=FONT * 1.4, color='red', zorder=2)

        ax.set_extent([-131, -119, 46, 52], crs=ccrs.PlateCarree())

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linestyle='--', alpha=1)
        gl.top_labels = True
        gl.right_labels = True
        gl.bottom_labels = False
        gl.left_labels = False
        gl.xlocator = mpl.ticker.FixedLocator(np.arange(-180, -49, 4))
        gl.ylocator = mpl.ticker.FixedLocator(np.arange(0, 81, 4))
        gl.label_style = {'inline': False}
        fig._elev_gl = gl   # resized on save
        return fig, ax

    def add_colorbar(fig, ax, cmap, vmin, vmax, label, ticks=None, extend='neither'):
        cbar_ax = fig.add_axes([0.15, 0.08, 0.73, 0.025])
        sm = cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(vmin, vmax))
        cb = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal', extend=extend,
                          ticks=ticks)
        cb.ax.tick_params(labelsize=FONT)
        cb.ax.set_xlabel(label, fontsize=FONT)
        cb.outline.set_linewidth(0.8)
        return cb

    def style_scatter(ax):
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontsize(FONT)

    def save(fig, title):
        out = os.path.join(OUT_DIR, title + '.png')
        # Gridline label artists are generated during the savefig draw; resize
        # them after the first render, then re-render so the size change is
        # captured.
        fig.savefig(out, dpi=300, bbox_inches='tight')
        gl = getattr(fig, '_elev_gl', None)
        if gl is not None:
            for lbl in gl.label_artists:
                lbl.set_fontsize(FONT)
            fig.savefig(out, dpi=300, bbox_inches='tight')
        print(f"  saved {out}")

    def plot_obs_elevations():
        fig, ax = new_domain_ax()
        ax.scatter(lons, lats, c=elev, s=35, cmap='terrain', vmin=0, vmax=3000,
                   transform=ccrs.PlateCarree(), edgecolor='k', linewidth=0.4,
                   zorder=3, marker='o')
        add_colorbar(fig, ax, 'terrain', 0, 3000, 'Elevation (m)',
                     ticks=np.arange(0, 3001, 500), extend='neither')
        style_scatter(ax)
        save(fig, 'obs_station_elevations')
        plt.close(fig)

    def plot_elev_bias(data, title, cmap, vmin, vmax, extend='both'):
        fig, ax = new_domain_ax()
        ax.scatter(lons, lats, c=data, s=35, cmap=cmap, vmin=vmin, vmax=vmax,
                   transform=ccrs.PlateCarree(), edgecolor='k', linewidth=0.4,
                   zorder=3, marker='o')
        add_colorbar(fig, ax, cmap, vmin, vmax, 'Elevation bias (m)',
                     ticks=np.arange(vmin, vmax + 1, 500), extend=extend)
        style_scatter(ax)
        save(fig, title)
        plt.close(fig)

    print("Plotting observed station elevations...")
    plot_obs_elevations()

    BIAS_VMIN, BIAS_VMAX = -1500, 1500
    bias_cmap = cm.get_cmap('PRGn', 24)
    print("Plotting elevation-bias maps...")
    plot_elev_bias(elev_bias_d03, 'CanESM2-WRF_D03_elev_bias', bias_cmap, BIAS_VMIN, BIAS_VMAX)
    plot_elev_bias(elev_bias_d02, 'CanESM2-WRF_D02_elev_bias', bias_cmap, BIAS_VMIN, BIAS_VMAX)
    plot_elev_bias(elev_bias_d01, 'CanESM2-WRF_D01_elev_bias', bias_cmap, BIAS_VMIN, BIAS_VMAX)
    plot_elev_bias(elev_bias_raw, 'CanESM2_elev_bias', bias_cmap, BIAS_VMIN, BIAS_VMAX)
    plot_elev_bias(elev_bias_rcm, 'CanRCM4_elev_bias', bias_cmap, BIAS_VMIN, BIAS_VMAX)
    print("Done.")

if __name__ == '__main__':
    _make_figures()
