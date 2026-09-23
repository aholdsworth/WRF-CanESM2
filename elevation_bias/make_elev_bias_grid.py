#!/usr/bin/env python
"""Build the 5×4 elevation-bias diagnostic grid figure.

Layout (per user's manual description of Eva's figure):
  5 rows  : D3 (WRF d03), D15 (WRF d02), D75 (WRF d01), CanRCM4, CanESM2
  4 cols  : (1) model topography + stations (per-row, shows resolution diff)
            (2) T-bias scatter  (3) PR-bias scatter  (4) wind-bias scatter
  Scatter cols: elevation bias (m, -1500..1500) on the y-axis, the variable
  bias on the x-axis.  Y-tick labels only on the outermost (T) column,
  x-tick labels only where the panel has data (D75 PR/wind have none),
  x-axis labels only on the bottom row.
  Portrait, shared terrain colorbar under the map column, bold row labels,
  panel tags (a)–(t), OLS fit + R², origin crosshairs.

Run: /home/amh001/space_fs7/software_2022/python/py_2024/bin/python make_elev_bias_grid.py
Output: figures/elevation_bias/elevation_bias_grid.png
"""
import os
import sys
import numpy as np
import pandas as pd
import pickle
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.projections import register_projection
from scipy import stats

register_projection(GeoAxes)

warnings.filterwarnings("ignore")

BASE = os.environ.get('ELEV_BIAS_BASE',
    '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS')
OUT  = os.path.join(BASE, 'figures', 'elevation_bias', 'elevation_bias_grid.png')
os.makedirs(os.path.dirname(OUT), exist_ok=True)

FONT = 13
FONT_S = 11
plt.rcParams.update({
    'font.size': FONT,
    'axes.labelsize': FONT,
    'xtick.labelsize': FONT_S,
    'ytick.labelsize': FONT_S,
    'axes.titlesize': FONT,
    'legend.fontsize': FONT_S,
    'figure.dpi': 150,
})

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# SEASONAL=1 -> read the season-first pooled bias layer (pickles/skill_seasonal.pkl,
# produced by make_skill_tables.py); default (unset) -> legacy catalogue
# (pickles/station_tables.pkl, produced by make_station_tables.py).  The
# season-first layer carries the same bias column names plus mae/rmse/n_years.
SEASONAL = os.environ.get('SEASONAL', '0') == '1'
_PKL = 'skill_seasonal.pkl' if SEASONAL else 'station_tables.pkl'
print(f'Using bias layer: {_PKL}')
with open(os.path.join(BASE, 'pickles', _PKL), 'rb') as f:
    tables = pickle.load(f)

# Row definitions: (label, table_key, colour, topo_file)
# topo_file = ('wrf', geo_em) or ('orog', orog_file)
ROWS = [
    ('D3',      'D03',     '#1f77b4', ('wrf',  'geo_em.d03.nc')),   # blue
    ('D15',     'D02',     '#ff7f0e', ('wrf',  'geo_em.d02.nc')),   # orange
    ('D75',     'D01',     '#2ca02c', ('wrf',  'geo_em.d01.nc')),   # green (emphasised)
    ('CanRCM4', 'CanRCM4', '#9467bd', ('orog', 'orog_CanRCM4.nc')),  # purple
    ('CanESM2', 'CanESM2', '#d62728', ('orog', 'orog_CanESM2.nc')),  # red
]

# Column definitions
COLS = [
    ('elev',  'Model topography',            None,   None),
    ('t',     'Temperature bias (°C)',    -17.5,  12.5),
    ('pr',    'Precipitation bias (%)',     -400,  400),
    ('wind',  'Wind speed bias (%)',       -175.0,  575.0),
]

MAP_EXTENT = [-132, -115, 43, 55]

# ---------------------------------------------------------------------------
# Load per-row topography
# ---------------------------------------------------------------------------
sys.path.insert(0, BASE)
from plot_elev_bias_revised import read_geo_em, read_orog

DOMAIN_DIR = os.path.join(BASE, 'plotting_files')

def _wrf_topo(fname):
    """geo_em d0X: XLAT_C/XLONG_C are (ny+1, nx+1) cell corners; HGT_M is
    (ny, nx) cell centres. Build centre coordinates and mask off-domain."""
    lat_c, lon_c, topo = read_geo_em(fname)
    lat1d = (lat_c[:-1, 0] + lat_c[1:, 0]) / 2.0
    lon1d = (lon_c[0, :-1] + lon_c[0, 1:]) / 2.0
    lons, lats = np.meshgrid(lon1d, lat1d)
    m = (lats >= MAP_EXTENT[2] - 0.5) & (lats <= MAP_EXTENT[3] + 0.5) & \
        (lons >= MAP_EXTENT[0] + 0.5) & (lons <= MAP_EXTENT[1] - 0.5)
    return lons, lats, np.ma.masked_where(~m, topo)

def _orog_topo(fname):
    lats, lons, topo = read_orog(fname)
    m = (lats >= MAP_EXTENT[2] - 0.5) & (lats <= MAP_EXTENT[3] + 0.5) & \
        (lons >= MAP_EXTENT[0] + 0.5) & (lons <= MAP_EXTENT[1] - 0.5)
    return lons, lats, np.ma.masked_where(~m, topo)

topos = {}
for row_label, _, _, (kind, fname) in ROWS:
    topos[row_label] = _wrf_topo(fname) if kind == 'wrf' else _orog_topo(fname)

# Shared terrain scale across rows
_all_topo = np.concatenate([t[2].compressed() for t in topos.values()])
TOPO_VMIN, TOPO_VMAX = float(np.min(_all_topo)), float(np.max(_all_topo))
cmap_topo = cm.get_cmap('terrain', 16)

# Elevation-bias colour scale for station points on the map panels
# (PRGn: purple = negative bias, green = positive bias), matching the
# per-dataset bias maps in plot_elev_bias_revised.py.
BIAS_VMIN, BIAS_VMAX = -1500, 1500
cmap_bias = cm.get_cmap('PRGn', 24)

# ---------------------------------------------------------------------------
# Figure setup
# ---------------------------------------------------------------------------
n_rows, n_cols = 5, 4
fig = plt.figure(figsize=(16, 20))
gs = gridspec.GridSpec(
    n_rows, n_cols,
    left=0.06, right=0.95, top=0.96, bottom=0.06,
    wspace=0.35, hspace=0.45,
    width_ratios=[1.1, 1, 1, 1],
)

# Track which row/col get x-axis labels (bottom row)
panel_tag = 0
def next_tag():
    global panel_tag
    tag = chr(ord('a') + panel_tag)
    panel_tag += 1
    return tag

for ri, (row_label, table_key, row_color, _) in enumerate(ROWS):
    df = tables[table_key]
    elev = df['elevation_bias_m'].values
    def _col(new, old):
        """Season-first layer uses {t,pr,wind}_bias; legacy uses the long
        names.  Return whichever exists."""
        return df[new].values if new in df.columns else df[old].values
    t_bias = _col('t_bias', 'temperature_bias_C')
    pr_bias = _col('pr_bias', 'precipitation_bias_pct')
    wind_bias = _col('wind_bias', 'wind_bias_pct')
    lats = df['lat'].values
    lons = df['lon'].values

    is_d75 = (row_label == 'D75')

    # Which scatter columns have points in this row?  Needed to cull
    # interior tick labels only where panels are actually populated
    # (D75 has no PR or wind points).
    col_has_pts = {
        't':    int((~np.isnan(elev) & ~np.isnan(t_bias)).sum()),
        'pr':   int((~np.isnan(elev) & ~np.isnan(pr_bias)).sum()),
        'wind': int((~np.isnan(elev) & ~np.isnan(wind_bias)).sum()),
    }

    for ci, (col_key, col_label, xmin, xmax) in enumerate(COLS):
        if col_key == 'elev':
            ax = GeoAxes(fig, gs[ri, ci].get_position(fig),
                         projection=ccrs.PlateCarree())
            fig.add_axes(ax)
        else:
            ax = fig.add_subplot(gs[ri, ci])

        # --- Panel tag ---
        tag = next_tag()
        ax.set_title(f'({tag})', fontweight='bold', color='#8B0000',
                     fontsize=FONT, loc='left', pad=4)

        if col_key == 'elev':
            # --- Map panel: this row's model topography + station points ---
            ax.set_facecolor('#f0f0f0')
            tlons, tlats, ttopo = topos[row_label]
            im = ax.pcolormesh(
                tlons, tlats, ttopo,
                cmap=cmap_topo,
                vmin=TOPO_VMIN, vmax=TOPO_VMAX,
                transform=ccrs.PlateCarree(),
            )
            # Gridlines (no coastlines - offline)
            ax.gridlines(draw_labels=True, linewidth=0.3, color='grey',
                         alpha=0.5, linestyle=':')
            ax.set_extent(MAP_EXTENT, crs=ccrs.PlateCarree())

            # Station points coloured by this row's elevation bias
            # (PRGn: purple = low bias, green = high bias, ±1500 m).
            # All 175 stations in every row (elevation bias is defined
            # for all of them).
            sc = ax.scatter(
                lons, lats,
                c=elev,
                s=18, cmap=cmap_bias,
                norm=mcolors.TwoSlopeNorm(vmin=BIAS_VMIN, vcenter=0,
                                          vmax=BIAS_VMAX),
                edgecolors='k', linewidths=0.3, zorder=5,
                transform=ccrs.PlateCarree(),
            )
            ax.set_aspect('equal')

        else:
            # --- Scatter panels ---
            # y = elevation bias (m); x = the variable's bias
            if col_key == 't':
                var = t_bias
            elif col_key == 'pr':
                var = pr_bias
            else:
                var = wind_bias

            m = ~np.isnan(elev) & ~np.isnan(var)
            xs, ys = var[m], elev[m]
            n_pts = len(xs)  # == col_has_pts[col_key]
            print(f'  panel {row_label:8s} {col_key:4s}: n_pts={n_pts}', flush=True)

            ax.scatter(xs, ys, s=18, color=row_color,
                       edgecolors='k', linewidths=0.3, zorder=5, alpha=0.85)

            # Origin crosshairs
            ax.axhline(0, color='grey', linewidth=0.5, linestyle='--', zorder=1)
            ax.axvline(0, color='grey', linewidth=0.5, linestyle='--', zorder=1)

            # OLS fit (variable bias ~ elevation bias; R is the Pearson
            # correlation of the same pair, matching Eva's original
            # convention of reporting the correlation of the two axes)
            if n_pts >= 3:
                slope, intercept, r_value, p_value, std_err = stats.linregress(ys, xs)
                y_fit = np.linspace(ys.min(), ys.max(), 100)
                ax.plot(slope * y_fit + intercept, y_fit,
                        'k--', linewidth=1, zorder=4)
                r2 = r_value ** 2
                ax.text(0.95, 0.95, f'R² = {r2:.3f}',
                        transform=ax.transAxes, ha='right', va='top',
                        fontsize=FONT_S, style='italic',
                        bbox=dict(boxstyle='round,pad=0.2', fc='white',
                                  ec='grey', alpha=0.7))

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-1500, 1500)
            ax.set_yticks([-1500, -750, 0, 750, 1500])

            # Interior tick labels: keep only on the outermost edge of each
            # column block (T: left, wind: right; PR: bottom row only, since
            # D75 has no PR points).  Tick marks stay, labels are culled.
            if col_key == 't':
                if ri < n_rows - 1:
                    for lbl in ax.get_xticklabels():
                        lbl.set_visible(False)
                # y tick labels culled on the T column (inner left edge)
                for lbl in ax.get_yticklabels():
                    lbl.set_visible(False)
            elif col_key == 'wind':
                if ri < n_rows - 1:
                    for lbl in ax.get_xticklabels():
                        lbl.set_visible(False)
            else:  # pr
                # bottom row: labels only if this panel has data
                # (D75's PR panel is empty)
                if n_pts == 0 or ri < n_rows - 1:
                    for lbl in ax.get_xticklabels():
                        lbl.set_visible(False)

        # --- Axis labels ---
        if col_key != 'elev':
            # x-label: the variable's bias, bottom row only
            if ri == n_rows - 1:
                ax.set_xlabel(col_label, fontsize=FONT)
            # y-label: elevation bias, T column only (outermost scatter col)
            if col_key == 't':
                ax.set_ylabel('Elevation bias (m)', fontsize=FONT)
            # PR/wind: no y-label, but y tick labels kept

        if ri == n_rows - 1 and col_key == 'elev':
            ax.set_xlabel('Longitude (°W)', fontsize=FONT)

        # Row emphasis for D75
        if is_d75:
            for spine in ax.spines.values():
                spine.set_linewidth(2.5)
                spine.set_color('black')

# ---------------------------------------------------------------------------
# Row labels (bold, rotated, left of each map)
# ---------------------------------------------------------------------------
for ri, (row_label, table_key, row_color, _) in enumerate(ROWS):
    fig.text(0.015, 1 - (ri + 0.5) / n_rows * 0.9, row_label,
             fontsize=FONT + 2, fontweight='bold',
             rotation=90, va='center', ha='center',
             color=row_color)

# Colorbars under the map column: model topography (left half) and
# elevation bias (right half, for the coloured station points).
cbar_ax1 = fig.add_axes([0.06, 0.015, 0.10, 0.012])
cbar1 = mpl.colorbar.ColorbarBase(
    cbar_ax1, cmap=cmap_topo,
    norm=mpl.colors.Normalize(vmin=TOPO_VMIN, vmax=TOPO_VMAX),
    orientation='horizontal', ticks=[0, 500, 1000, 1500, 2000, 2500, 3000],
)
cbar1.set_label('Model topography (m)', fontsize=FONT_S)
cbar1.outline.set_linewidth(0.8)

cbar_ax2 = fig.add_axes([0.22, 0.015, 0.10, 0.012])
cbar2 = mpl.colorbar.ColorbarBase(
    cbar_ax2, cmap=cmap_bias,
    norm=mcolors.TwoSlopeNorm(vmin=BIAS_VMIN, vcenter=0, vmax=BIAS_VMAX),
    orientation='horizontal', ticks=[-1500, -750, 0, 750, 1500],
)
cbar2.set_label('Elevation bias (m)', fontsize=FONT_S)
cbar2.outline.set_linewidth(0.8)

fig.suptitle('Elevation-Bias Diagnostic Grid (1986–2005)',
             fontsize=FONT + 4, fontweight='bold', y=0.985)

# Debug report: tick-label visibility per scatter panel
if os.environ.get('GRID_DEBUG'):
    for ri in range(n_rows):
        for ci in range(1, n_cols):
            gp = gs[ri, ci].get_position(fig)
            best, best_d = None, 1e9
            for a in fig.axes:
                if isinstance(a, GeoAxes):
                    continue
                p = a.get_position()
                d = abs(p.x0 - gp.x0) + abs(p.y0 - gp.y0)
                if d < best_d:
                    best, best_d = a, d
            ax = best
            xv = sum(1 for l in ax.get_xticklabels() if l.get_visible())
            yv = sum(1 for l in ax.get_yticklabels() if l.get_visible())
            print(f'DEBUG {ROWS[ri][0]:8s} {COLS[ci][0]:4s} '
                  f'xlabels={xv} ylabels={yv} xlim={ax.get_xlim()} ylim={ax.get_ylim()}')

fig.savefig(OUT, bbox_inches='tight', dpi=150)
print(f'Saved {OUT}')
