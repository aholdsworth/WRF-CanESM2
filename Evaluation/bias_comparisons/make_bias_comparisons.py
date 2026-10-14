#!/usr/bin/env python
"""Reproduce the bias_comparisons figure.

Layout: 4 rows (models, top -> bottom: D02, D01, CanRCM4, CanESM2) x
3 columns (variables): Temperature bias (deg C), Precipitation bias (%),
Wind speed bias (%).

Each panel: the model's own topography as background shading (terrain
colormap, 0..3000 m; each panel zoomed to its own topo coverage so the
frame edges are clean) with per-station bias dots coloured on the
column's diverging scale.  Colormaps match the local git repo's mean
plots (WRF-CanESM2/mean_changes/plot_means_*.ipynb):
  T    -> RdBu_r   (-4..+4 deg C)
  PR   -> BrBG     (-100..+100 %)
  wind -> PuOr     (-100..+100 %)
One shared colorbar per column.  Stations missing a variable are omitted.

Data sources (2026-09-29 rewire: aligned with the bias_figure gated set):
  * pickles/bias_figure.pkl      station set + obs reference: the Rule-1
    gated period values ('gated'[var][0], MIN_YEARS=15, incl. buoys) --
    exactly the stations shown in figures/bias/bias_figure.png (t 154,
    pr 159, wind 22 finite), plus their lat/lon from 'stations'.
  * data/wrf_canon/*.nc          row model period values for the D02/D01
    panels, same files & pooled conventions as make_station_tables.py
    (t/wind: period mean; pr: mean of annual sums, mm/yr).
  * seasonal_pooling.cdo_station_value  row model period values for the
    CanESM2 / CanRCM4 panels (M9/M10 CDO per-station files, same source
    as the station-table build).
  * pickles/buoy_{tas,wind}_tables.pkl  validated buoy obs & model
    period values ('tab'[(buoy, model)]: obs, model_mean) for the buoy
    dots (8 tas buoys, 7 wind buoys; no pr buoys -- matches Fig-1).
  * plotting_files/geo_em.dXX.nc WRF HGT_M topography (d01-d03)
  * plotting_files/orog_Can*.nc  CanESM2 / CanRCM4 orog topography

Dot value = row model period value - gated obs (t: deg C; pr/wind:
100*(m-o)/o).  Land rows are data-driven NaN where the row model file is
absent; the station SET and obs reference are always the gated Fig-1 set.

Output: figures/bias_comparisons.png  (does NOT touch the placeholder
bias_comparisons.png in the workspace root).
"""
import os
import sys
import pickle
import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.projections import register_projection
from netCDF4 import Dataset

register_projection(GeoAxes)

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
DOMAIN_DIR = os.path.join(BASE, 'plotting_files')
DATA = os.path.join(BASE, 'data')
WRF_ST_DIR = os.path.join(DATA, 'wrf_canon')
START_YEAR, END_YEAR = 1986, 2005

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'elevation_bias'))  # shared seasonal_pooling.py
import seasonal_pooling as sp
OUT_DIR = os.path.join(BASE, 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
OUT = os.path.join(OUT_DIR, 'bias_comparisons.png')


# --- Topography readers (same convention as plot_elev_bias_revised.py) -------
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


# --- Panel definitions -------------------------------------------------------
ROWS = [
    ('D02',     'D02',     'geo_em.d02.nc',   'wrf'),
    ('D01',     'D01',     'geo_em.d01.nc',   'wrf'),
    ('CanRCM4', 'CanRCM4', 'orog_CanRCM4.nc', 'orog'),
    ('CanESM2', 'CanESM2', 'orog_CanESM2.nc', 'orog'),
]
# column: (key, title, cmap, half-range)  -- cmaps from WRF-CanESM2 mean plots
COLS = [
    ('t',    'Temperature bias (°C)', 'RdBu_r', 4.0),
    ('pr',   'Precipitation bias (%)', 'BrBG',  100.0),
    ('wind', 'Wind speed bias (%)',   'PuOr',   100.0),
]

MAP_EXTENT = [-128.5, -121.5, 45.5, 52.5]   # same window as elev_bias_grid
TOPO_VMIN, TOPO_VMAX = 0, 3000
FONT = 8

# --- Load gated obs (Fig-1 station set + obs reference) -----------------------
with open(os.path.join(BASE, 'pickles', 'bias_figure.pkl'), 'rb') as f:
    bf = pickle.load(f)
bf_st = bf['stations']
# gated[var] = (gated_obs, gated_model) Series, indexed by station id string
gated_obs = {k: bf['gated'][k][0] for k in ('t', 'pr', 'wind')}
gated_ids = {k: [i for i in gated_obs[k].index if np.isfinite(gated_obs[k][i])]
             for k in gated_obs}
lat_s = bf_st['lat']
lon_s = bf_st['lon']

# --- Buoy obs/model period values (validated tables, never wrf_canon) --------
with open(os.path.join(BASE, 'pickles', 'buoy_tas_tables.pkl'), 'rb') as f:
    _bt = pickle.load(f)
with open(os.path.join(BASE, 'pickles', 'buoy_wind_tables.pkl'), 'rb') as f:
    _bw = pickle.load(f)
TAS_BUOY_TAB = _bt['tab']
WIND_BUOY_TAB = _bw['tab']
# row label -> model key inside the buoy tables
BUOY_MODEL_KEY = {'D02': 'WRF_d02', 'D01': 'WRF_d01',
                  'CanRCM4': 'CanRCM4', 'CanESM2': 'CanESM2'}
# Fig-1 buoy cores (make_bias_figure.py): 8 tas, 7 wind (no 46050 in wind)
WIND_BUOYS_D03 = ['46041', '46131', '46132', '46146', '46204', '46206', '46207']


def buoy_period(buoy, tab, row_label, field):
    """Period obs/model value for one buoy row; (obs, model) or (NaN, NaN)."""
    r = tab.loc[(str(buoy), BUOY_MODEL_KEY[row_label])]
    return float(r[field]), float(r['model_mean'])


# --- Row model period values (same files as make_station_tables.py) ----------
def _wrf_st_files(var, dom):
    """station_id -> file path for data/wrf_canon {var}_{dom}_{AGENCY}_{ID}.nc,
    keyed as 'E'+ClimateID / 'B'+code / 'N'+USC*/USW* (table convention)."""
    out = {}
    prefix = var + '_' + dom + '_'
    for f in os.listdir(WRF_ST_DIR):
        if not f.endswith('.nc') or not f.startswith(prefix):
            continue
        rest = f[len(prefix):-3]
        if rest.startswith('ECCC_'):
            out['E' + rest[5:]] = os.path.join(WRF_ST_DIR, f)
        elif rest.startswith('BCH_'):
            out['B' + rest[4:]] = os.path.join(WRF_ST_DIR, f)
        elif rest.startswith('NOAA_'):
            out['N' + rest[5:]] = os.path.join(WRF_ST_DIR, f)
    return out


def _wrf_station_series(path, var_name, kind):
    """One per-station WRF file -> Series (same convention as the table build)."""
    nc = Dataset(path, mode='r')
    t = nc.variables['time'][:]
    data = np.squeeze(nc.variables[var_name][:])
    nc.close()
    dates = pd.to_datetime(
        datetime.datetime(START_YEAR, 1, 1) + pd.to_timedelta(t, unit='h'))
    if kind == 't':
        v = data - 273.15
    elif kind == 'pr':
        v = np.clip(np.diff(data), 0.0, None)
        dates = dates[1:]
    elif kind == 'wind':
        v = data
    return pd.Series(v, index=dates)


def _annual_mean(series, kind):
    """Period value: mean for t/wind; mean of annual sums for pr (mm/yr)."""
    s = series[(series.index >= datetime.datetime(START_YEAR, 1, 1)) &
               (series.index < datetime.datetime(END_YEAR + 1, 1, 1))]
    if kind == 'pr':
        return float(s.groupby(s.index.year).sum().mean())
    return float(s.mean())


WRF_FILES = {}
for _dom in ('d01', 'd02'):
    for _var in ('t', 'pr', 'wind'):
        WRF_FILES[(_dom, _var)] = _wrf_st_files(_var, _dom)

# CanESM2/CanRCM4 row model values: same M9/M10 CDO files as the table build.
CAN_MODEL = {}
for _src, _mkey in (('CanESM2', 'canesm2'), ('CanRCM4', 'canrcm4')):
    CAN_MODEL[_src] = {}
    for _k in ('t', 'pr', 'wind'):
        CAN_MODEL[_src][_k] = {sid: sp.cdo_station_value(sid, _k, DATA, _mkey)
                               for sid in bf_st['lat'].index}

# --- Topographies, masked to the map extent ----------------------------------
def _wrf_topo(fname):
    """geo_em d0X: XLAT_C/XLONG_C are (ny+1, nx+1) corners; HGT_M centres.
    Pass 2D corners straight to pcolormesh (curvilinear grid)."""
    lat_c, lon_c, topo = read_geo_em(fname)
    lat_ctr = 0.25 * (lat_c[:-1, :-1] + lat_c[1:, :-1] + lat_c[:-1, 1:] + lat_c[1:, 1:])
    lon_ctr = 0.25 * (lon_c[:-1, :-1] + lon_c[1:, :-1] + lon_c[:-1, 1:] + lon_c[1:, 1:])
    m = (lat_ctr >= MAP_EXTENT[2] - 0.5) & (lat_ctr <= MAP_EXTENT[3] + 0.5) & \
        (lon_ctr >= MAP_EXTENT[0] + 0.5) & (lon_ctr <= MAP_EXTENT[1] - 0.5)
    return lon_c, lat_c, np.ma.masked_where(~m, topo)

def _orog_topo(fname):
    lats, lons, topo = read_orog(fname)
    m = (lats >= MAP_EXTENT[2] - 0.5) & (lats <= MAP_EXTENT[3] + 0.5) & \
        (lons >= MAP_EXTENT[0] + 0.5) & (lons <= MAP_EXTENT[1] - 0.5)
    return lons, lats, np.ma.masked_where(~m, topo)

cmap_topo = cm.get_cmap('terrain', 16)
col_cmaps = {key: cm.get_cmap(cname) for key, _, cname, _ in COLS}
norms = {key: mcolors.TwoSlopeNorm(vmin=-half, vcenter=0, vmax=half)
         for key, _, _, half in COLS}

# --- Figure layout: 4 x 3 grid of maps + one colorbar per column --------------
def row_model_value(row_label, sid, col_key):
    """Row model period value for one land station (NaN if no file)."""
    if row_label in ('D02', 'D01'):
        dom = {'D02': 'd02', 'D01': 'd01'}[row_label]
        agency = 'E' if sid.isdigit() else (
            'N' if sid.startswith(('USC', 'USW')) else 'B')
        p = WRF_FILES[(dom, col_key)].get(agency + sid)
        if p is None:
            return np.nan
        return _annual_mean(_wrf_station_series(
            p, {'t': 'T2', 'pr': 'pr', 'wind': 'wspd'}[col_key], col_key), col_key)
    return CAN_MODEL[row_label][col_key].get(sid, np.nan)


def panel_points(row_label, col_key):
    """(lons, lats, bias values) for one panel: gated land stations + buoys."""
    o = gated_obs[col_key]
    lons, lats, vals = [], [], []
    for sid in gated_ids[col_key]:
        if sid in bf_st['buoy_ids']:
            continue  # handled below via the validated buoy tables
        om = o[sid]
        if col_key == 't':
            m = row_model_value(row_label, sid, 't')
            if not np.isfinite(m):
                continue
            b = m - om
        else:
            m = row_model_value(row_label, sid, col_key)
            if not np.isfinite(m) or om <= 0:
                continue
            b = 100.0 * (m - om) / om
        lons.append(float(lon_s[sid]))
        lats.append(float(lat_s[sid]))
        vals.append(b)
    # buoy dots from the validated table pickles (obs & model)
    if col_key == 't':
        buoys, tab = list(bf_st['buoy_ids']), TAS_BUOY_TAB
    elif col_key == 'wind':
        buoys = [b for b in WIND_BUOYS_D03
                 if b in WIND_BUOY_TAB.index.get_level_values(0)]
        tab = WIND_BUOY_TAB
    else:  # pr: no buoy dots (matches Fig-1)
        buoys, tab = [], None
    for b in buoys:
        if b not in tab.index.get_level_values(0):
            continue
        ob, mb = buoy_period(b, tab, row_label, 'obs')
        if not (np.isfinite(ob) and np.isfinite(mb)):
            continue
        bb = mb - ob if col_key == 't' else 100.0 * (mb - ob) / ob
        lons.append(float(lon_s[b]))
        lats.append(float(lat_s[b]))
        vals.append(bb)
    return np.array(lons), np.array(lats), np.array(vals)


n_rows, n_cols = len(ROWS), len(COLS)
fig = plt.figure(figsize=(7.6, 6.6))
gs = gridspec.GridSpec(n_rows, n_cols, left=0.09, right=0.97,
                       top=0.965, bottom=0.16, wspace=0.22, hspace=0.34)

col_mesh = {}
for ri, (row_label, table_key, fname, kind) in enumerate(ROWS):
    tlons, tlats, ttopo = _wrf_topo(fname) if kind == 'wrf' else _orog_topo(fname)

    for ci, (col_key, col_title, _, half) in enumerate(COLS):
        lons, lats, var = panel_points(row_label, col_key)
        m = ~np.isnan(var)
        ax = fig.add_subplot(gs[ri, ci], projection=ccrs.PlateCarree())

        # panel tag
        ax.text(0.04, 0.94, chr(ord('a') + ri * n_cols + ci), transform=ax.transAxes,
                fontweight='bold', fontsize=FONT, ha='left', va='top', zorder=10)

        ax.pcolormesh(tlons, tlats, ttopo, cmap=cmap_topo,
                      vmin=TOPO_VMIN, vmax=TOPO_VMAX,
                      transform=ccrs.PlateCarree())
        # Zoom to the intersection of the map window and this model's own
        # topo coverage, so every panel frame is filled with terrain
        # (no rough masked edges).
        tlons_flat = tlons.compressed()
        tlats_flat = tlats.compressed()
        ext = [max(MAP_EXTENT[0], np.nanmin(tlons_flat)),
               min(MAP_EXTENT[1], np.nanmax(tlons_flat)),
               max(MAP_EXTENT[2], np.nanmin(tlats_flat)),
               min(MAP_EXTENT[3], np.nanmax(tlats_flat))]
        ax.set_extent(ext, crs=ccrs.PlateCarree())

        mesh = ax.scatter(lons[m], lats[m], c=var[m], s=8,
                          cmap=col_cmaps[col_key], norm=norms[col_key],
                          edgecolors='k', linewidths=0.25,
                          zorder=5, transform=ccrs.PlateCarree())
        col_mesh.setdefault(ci, mesh)

        if ci == 0:
            ax.text(-0.14, 0.5, row_label, transform=ax.transAxes,
                    fontsize=FONT + 1, ha='center', va='center', rotation='vertical')
        if ri == n_rows - 1:
            ax.set_title(col_title, fontsize=FONT)

        for spine in ax.spines.values():
            spine.set_linewidth(0.7)
            spine.set_color('black')

# --- One shared colorbar per column -------------------------------------------
for ci, (col_key, col_title, _, half) in enumerate(COLS):
    ax_pos = gs[n_rows - 1, ci].get_position(fig)
    x0 = ax_pos.x0 + 0.25 * ax_pos.width
    xw = 0.50 * ax_pos.width
    cax = fig.add_axes([x0, 0.045, xw, 0.012])
    cbar = fig.colorbar(col_mesh[ci], cax=cax, orientation='horizontal',
                        ticks=[-half, 0, half], extend='both')
    cbar.set_label(col_title, fontsize=FONT)
    cbar.outline.set_linewidth(0.7)

fig.savefig(OUT, dpi=200, bbox_inches='tight')
print('Saved', OUT)
