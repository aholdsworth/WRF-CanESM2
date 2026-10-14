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


BASE   = os.environ.get('WRF_CANESM2_BASE') or _workroot(os.path.dirname(os.path.abspath(__file__)))
STATION_DIR = os.path.join(BASE, 'files')          # unzipped station lists
DOMAIN_DIR  = os.path.join(BASE, 'plotting_files') # geo_em / orog / namelist
REPO_DIR    = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))    # repo root (WRFDomainLib)
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
_map = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'station_cdo', 'eccc_bch_desc_map.json')
if not os.path.exists(_map):
    _map = os.path.join(BASE, 'station_cdo', 'eccc_bch_desc_map.json')
with open(_map) as _f:
    _ECCC_SID2CID = json.load(_f)['ECCC']

sys.path.insert(0, REPO_DIR)
import WRFDomainLib

FONT = 11  # requested font size (points)

# ----------------------------------------------------------------------------
# Station metadata -- now provided by the leaf module station_geo.py (single
# source of truth; imported here for back-compat so `from
# plot_elev_bias_revised import load_*` keeps working).
# ----------------------------------------------------------------------------
from station_geo import (load_eccc, load_bch, load_noaa, load_noaa_land,
                         load_buoys)


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
# Domain topography + model-elevation machinery now live in station_geo.py
# (single source of truth, shared with make_station_tables.py).  Imported for
# back-compat: `from plot_elev_bias_revised import read_geo_em, read_orog, ...`
# keeps working.
# ----------------------------------------------------------------------------
from station_geo import (read_geo_em, read_orog, wrf_t_files,
                         _weight_file_from_history, _weighted_elev,
                         weighted_model_elevs, get_nearest_gridpoint,
                         modeled_elevs, _can_weight_file)



if __name__ == '__main__':
    # --- Domain topography (recomputed here so importing this module has
    # no I/O side effects) ---
    lat_d03, lon_d03, topo_d03 = read_geo_em('geo_em.d03.nc')
    lat_d02, lon_d02, topo_d02 = read_geo_em('geo_em.d02.nc')
    lat_d01, lon_d01, topo_d01 = read_geo_em('geo_em.d01.nc')
    lats_canesm2, lons_canesm2, topo_canesm2 = read_orog('orog_CanESM2.nc')
    lats_canrcm4, lons_canrcm4, topo_canrcm4 = read_orog('orog_CanRCM4.nc')
    print("Computing model elevation for each station...")
    # WRF domains: use the extraction-consistent CDO bicubic weights (the orography
    # actually baked into the model data we compare).  Falls back to the corner
    # nearest-grid-point lookup for any station without a usable T-file weight file
    # (e.g. buoys without T files).
    def _wrf_model_elevs(dom, lat, lon, topo):
        s, warns = weighted_model_elevs(dom, list(elev.index))
        corner = modeled_elevs(lat, lon, topo, lats, lons, elev)
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
    def _can_model_elevs(src, lat2d, lon2d, topo2d):
        """Per-station Can* model elevation from the extraction's own bicubic
        weights, with nearest-GP fallback for stations without a weight file.
        Returns (series, fallback_ids)."""
        s, warns = weighted_model_elevs(src, list(elev.index))
        nearest = modeled_elevs(lat2d, lon2d, topo2d, lats, lons, elev)
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
    _make_figures()
