#!/usr/bin/env python
"""
taylor_figure.py -- publication-quality 3-panel Taylor diagram
(T / PR / wind) for the CanESM2-WRF vs CanESM2 / CanRCM4 evaluation.

Convention (2026-09-26, re-anchored -- see taylor_3panel_handoff.md):
  * Each model/obs quantity is the per-station ANNUAL series 1986-2005
    (t/wind: annual mean; pr: annual sum in mm).
  * For each model, over the stations it shares with obs, we pool the
    (station x year) samples over the years where BOTH obs and model are
    valid.  NO anomaly detrending.
  * refstd / std / corr:
        refstd = std(obs pooled sample)   (union of all models' shared stations)
        std    = std(model pooled sample)
        corr   = Pearson(model, obs) on the pooled sample

Layout revision (spacing fixes):
  * correlation tick labels + "Correlation" title are drawn with annotate()
    (offsets in points) instead of axisartist's rotated labels, which
    collided with each other and with the axis title;
  * panels are square in inches and separated by a wider gap so std tick
    labels never touch the neighbouring panel;
  * RMS-contour labels share one diagonal in every panel (bug fix: the
    target angle used to leak from the first contour into the rest).

Run with the py_2024 environment:
  /home/amh001/space_fs7/software_2022/python/py_2024/bin/python taylor_figure.py
"""
import os, sys
import warnings
import numpy as np
from scipy.stats import pearsonr

import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF
from matplotlib.projections import PolarAxes

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------
# Paths & style
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


BASE = os.environ.get('WRF_CANESM2_BASE') or _workroot(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bias_figure'))  # taylor_data.py
OUT_DIR = os.path.join(BASE, 'figures', 'taylor')
os.makedirs(OUT_DIR, exist_ok=True)
OUT = os.path.join(OUT_DIR, 'taylor_3panel.png')

# Figure is built at its FINAL printed size (Wiley / IJoC: 80-180 mm wide),
# so the point sizes below are what the reader sees -- no rescaling.
FIG_W_MM = 180.0                       # full page width
FIG_W = FIG_W_MM / 25.4                # inches
FIG_H = 3.3
FONT = 9
FONT_S = 7.5
plt.rcParams.update({
    'font.size': FONT,
    'axes.labelsize': FONT_S,
    'xtick.labelsize': FONT_S,
    'ytick.labelsize': FONT_S,
    'axes.titlesize': FONT,
    'legend.fontsize': FONT_S,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype': 42,                # embed TrueType, editable text
    'ps.fonttype': 42,
    'savefig.dpi': 600,
})

MODELS = ['D03', 'D02', 'D01', 'CanESM2', 'CanRCM4']
LABELS = {'D03': 'CanESM2-WRF D3', 'D02': 'CanESM2-WRF D15',
          'D01': 'CanESM2-WRF D75', 'CanESM2': 'CanESM2', 'CanRCM4': 'CanRCM4'}
COLORS = {'D03': '#1f77b4', 'D02': '#ff7f0e', 'D01': '#2ca02c',
          'CanESM2': '#d62728', 'CanRCM4': '#9467bd'}
MARKERS = {'D03': 's', 'D02': '^', 'D01': 'o', 'CanESM2': 'P', 'CanRCM4': 'p'}

VARS = [
    ('t',    '(a) Temperature',    r'$\degree{}$C'),
    ('pr',   '(b) Precipitation',  'mm y$^{-1}$'),
    ('wind', '(c) Wind speed',     'm s$^{-1}$'),
]
VAR_UNITS = {key: unit for key, _, unit in VARS}

CORR_TICKS = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0)


# ----------------------------------------------------------------------------
# Taylor diagram
# ----------------------------------------------------------------------------
class TaylorDiagram(object):
    def __init__(self, refstd, fig, rect, label='Observations',
                 srange=(0, 1.6), max_std=None, unit='unit',
                 n_intervals=5, corr_ticks=CORR_TICKS):
        self.refstd = refstd
        self.corr_ticks = corr_ticks
        tr = PolarAxes.PolarTransform()

        rlocs = np.array(corr_ticks, dtype=float)
        self.tmax = np.pi / 2.0
        tlocs = np.arccos(rlocs)
        gl1 = GF.FixedLocator(tlocs)
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        self.smin = srange[0] * self.refstd
        self.smax = max_std if max_std is not None else srange[1] * self.refstd
        std_ticks = np.linspace(self.smin, self.smax, n_intervals + 1)
        gl2 = GF.FixedLocator(std_ticks)

        ghelper = FA.GridHelperCurveLinear(
            tr, extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1, grid_locator2=gl2)

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_axes(ax)

        # Arc: keep the ticks, but draw the labels ourselves (see
        # _draw_corr_labels) so nothing collides with the axis title.
        ax.axis['top'].set_axis_direction('bottom')
        ax.axis['top'].toggle(ticklabels=False, label=False)

        # Horizontal std axis (bottom of the figure)
        ax.axis['left'].set_axis_direction('bottom')
        ax.axis['left'].label.set_text('Standard deviation (%s)' % unit)
        ax.axis['left'].label.set_pad(2)

        # Vertical std axis (left of the figure)
        ax.axis['right'].set_axis_direction('top')
        ax.axis['right'].toggle(ticklabels=True)
        ax.axis['right'].major_ticklabels.set_axis_direction('left')
        if self.smin:
            ax.axis['bottom'].toggle(ticklabels=False, label=False)
        else:
            ax.axis['bottom'].set_visible(False)

        self._ax = ax
        self.ax = ax.get_aux_axes(tr)

        l, = self.ax.plot([0.001], self.refstd, 'k*', ls='', ms=13,
                          zorder=1000, clip_on=False, label='observations')
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', lw=1.0, zorder=1)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, **kwargs):
        c = np.clip(corrcoef, -1.0, 1.0)
        l, = self.ax.plot(np.arccos(c), stddev, **kwargs)
        self.samplePoints.append(l)
        return l

    def add_contours(self, levels, **kwargs):
        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax, 120),
                             np.linspace(0, self.tmax, 120))
        rms = np.sqrt(self.refstd**2 + rs**2 - 2 * self.refstd * rs * np.cos(ts))
        return self.ax.contour(ts, rs, rms, levels, **kwargs)

    def label_contours(self, ct, fmt='%.1f', fontsize=FONT_S, color='0.35',
                       corr=0.55):
        """Label every RMS contour along one diagonal of constant correlation.

        For a contour of RMS difference L, the radius r at angle a satisfies
            r^2 - 2 refstd cos(a) r + (refstd^2 - L^2) = 0.
        Contours with L < refstd never reach small-corr diagonals (they hug
        the reference star), so for those the angle is pushed just past the
        tangent point.  Each label uses the root closest to mid-radius.
        """
        r_mid = 0.5 * (self.smin + self.smax)
        placed = []
        for lvl in ct.levels:
            a = np.arccos(np.clip(corr, -1, 1))          # local, per level
            if lvl < self.refstd:
                a = max(a, np.arcsin(lvl / self.refstd) + 0.08)
            a = min(a, 0.97 * self.tmax)
            b = -2.0 * self.refstd * np.cos(a)
            disc = b * b - 4.0 * (self.refstd**2 - lvl**2)
            if disc < 0:
                continue
            roots = [(-b - np.sqrt(disc)) / 2.0, (-b + np.sqrt(disc)) / 2.0]
            roots = [r for r in roots
                     if self.smin + 1e-9 <= r <= self.smax - 1e-9]
            if not roots:
                continue
            rr = min(roots, key=lambda r: abs(r - r_mid))
            self.ax.text(a, rr, fmt % lvl, ha='center', va='center',
                         fontsize=fontsize, color=color, zorder=500,
                         bbox=dict(boxstyle='round,pad=0.12', fc='white',
                                   ec='none', alpha=0.8, lw=0))
            placed.append((lvl, a, rr))
        return placed

    def _draw_corr_labels(self, fontsize=FONT_S, pad_pt=5, title_pad_pt=40):
        """Upright correlation tick labels just outside the arc, plus the
        'Correlation' title rotated along the arc at 45 degrees.  Offsets are
        in points, so spacing is independent of the data scale."""
        for c in self.corr_ticks:
            a = np.arccos(c)
            dx, dy = np.cos(a), np.sin(a)
            ha = 'left' if dx > 0.5 else 'center'
            va = 'bottom' if dy > 0.5 else 'center'
            self.ax.annotate(str(float(c)), xy=(a, self.smax),
                             xytext=(pad_pt * dx, pad_pt * dy),
                             textcoords='offset points', ha=ha, va=va,
                             fontsize=fontsize, annotation_clip=False)
        a = np.pi / 4.0
        self.ax.annotate('Correlation', xy=(a, self.smax),
                         xytext=(title_pad_pt * np.cos(a),
                                 title_pad_pt * np.sin(a)),
                         textcoords='offset points', ha='center', va='center',
                         rotation=-45, rotation_mode='anchor',
                         fontsize=fontsize, annotation_clip=False)

    def put_ticks_out(self):
        self._ax.axis[:].major_ticks.set_tick_out(True)

    def _cosmetics(self):
        for loc in self._ax.axis.values():
            try:
                loc.major_ticklabels.set_fontsize(FONT_S)
                loc.major_ticklabels.set_pad(3)
            except Exception:
                pass
        self._ax.axis['left'].label.set_fontsize(FONT_S)
        self._draw_corr_labels()


# ----------------------------------------------------------------------------
# Stats
# ----------------------------------------------------------------------------
def _shared_samples(obs, model):
    """(obs_samples, model_samples) pooled over shared stations/years."""
    o_parts, m_parts = [], []
    for s in sorted(set(obs) & set(model)):
        o = np.asarray(obs[s], dtype=float).reshape(-1)
        m = np.asarray(model[s], dtype=float).reshape(-1)
        if o.size != m.size:
            continue
        mask = ~(np.isnan(o) | np.isnan(m))
        if mask.sum() < 2:
            continue
        o_parts.append(o[mask])
        m_parts.append(m[mask])
    if not o_parts:
        return np.array([]), np.array([])
    return np.concatenate(o_parts), np.concatenate(m_parts)


def taylor_stats(obs, models, label):
    """(refstd, std, corr, n_stations) for one model over its shared stations."""
    m = models.get(label, {})
    o, mm = _shared_samples(obs, m)
    n_st = len([s for s in (set(obs) & set(m))
                if np.asarray(obs[s]).size == np.asarray(m[s]).size])
    if len(o) < 5:
        return np.nan, np.nan, np.nan, n_st
    return float(np.std(o)), float(np.std(mm)), float(pearsonr(o, mm)[0]), n_st


def panel_refstd(obs, models):
    """refstd from OBS over the union of all models' shared station-years."""
    keep = set()
    for m in MODELS:
        keep |= set(obs) & set(models.get(m, {}))
    o_parts = []
    for s in sorted(keep):
        o = np.asarray(obs[s], dtype=float)
        o = o[~np.isnan(o)]
        if o.size:
            o_parts.append(o)
    if not o_parts:
        return np.nan
    return float(np.std(np.concatenate(o_parts)))


# ----------------------------------------------------------------------------
# Auto-scaling
# ----------------------------------------------------------------------------
def _nice_ceiling(x):
    """Round x UP to a 'nice' number (1/2/2.5/5/10 x a power of ten)."""
    if not np.isfinite(x) or x <= 0:
        return 1.0
    exp = np.floor(np.log10(x))
    frac = x / 10 ** exp
    for m in (1, 2, 2.5, 5, 10):
        if frac <= m:
            return m * 10 ** exp
    return 10 * 10 ** exp


def _round_axis_ceiling(x, n_intervals=5):
    """Smallest (1/2/5 x 10^k step) * n_intervals that covers x."""
    if not np.isfinite(x) or x <= 0:
        return n_intervals, 1.0
    e = int(np.floor(np.log10(x)))
    cands = sorted(m * 10 ** d for d in range(e - 2, e + 3) for m in (1, 2, 5))
    step = next((s for s in cands if s * n_intervals >= x), cands[-1])
    return step * n_intervals, step


def _auto_levels(max_std, n=4):
    step = _nice_ceiling(max_std / (n + 1))
    return [step * (i + 1) for i in range(n)]


def _draw_panel(fig, rect, data, var):
    obs, models = data[var]

    refstd = panel_refstd(obs, models)
    stats = {}
    for m in MODELS:
        _, std, corr, n = taylor_stats(obs, models, m)
        stats[m] = (std, corr, n)

    finite_stds = [std for std, _, _ in stats.values() if not np.isnan(std)]
    panel_max_raw = max([refstd] + finite_stds)
    max_std, tick_step = _round_axis_ceiling(panel_max_raw * 1.15, n_intervals=5)
    rms_levels = _auto_levels(max_std)
    fmt = '%.0f' if tick_step >= 1 else ('%.1f' if tick_step >= 0.1 else '%.2f')

    dia = TaylorDiagram(refstd, fig=fig, rect=rect, label='Observations',
                        max_std=max_std, unit=VAR_UNITS[var])
    dia.put_ticks_out()
    ct = dia.add_contours(rms_levels, colors='0.55', linewidths=0.9, zorder=0)
    dia.label_contours(ct, fmt=fmt, corr=0.55)   # same diagonal in all panels
    for m in MODELS:
        std, corr, n = stats[m]
        if not (np.isnan(std) and np.isnan(corr)):
            dia.add_sample(std, corr, marker=MARKERS[m], ms=7,
                           color=COLORS[m], mec=COLORS[m], mew=0.8,
                           zorder=3, label=LABELS[m])
    dia._cosmetics()
    return dia


def make_figure():
    import taylor_data as td

    meta = td.load_meta()
    data = {v: td.build(v, meta) for v in ('t', 'pr', 'wind')}

    # ------------------------------------------------------------------
    # Layout (figure fractions) at the final 180 mm width.  Panels are
    # square in inches; the legend is a single row beneath the panels so
    # the three diagrams can use the full width.  No bbox_inches='tight':
    # the saved size is exactly FIG_W x FIG_H.
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    x0 = 0.065
    panel_w = 0.235
    panel_gap = 0.095
    panel_h = panel_w * FIG_W / FIG_H          # square in inches
    panel_y = 0.27                             # tick labels + axis title + legend row

    rects = [[x0 + i * (panel_w + panel_gap), panel_y, panel_w, panel_h]
             for i in range(3)]

    for i, (var, title, unit) in enumerate(VARS):
        _draw_panel(fig, rects[i], data, var)
        fig.text(rects[i][0] + rects[i][2] / 2.0,
                 panel_y + panel_h + 0.135, title,
                 ha='center', va='bottom', fontsize=FONT)

    handles = [plt.Line2D([], [], marker='*', ls='', ms=12, color='black',
                          mec='black', mew=0.8, label='Observations')]
    handles += [plt.Line2D([], [], marker=MARKERS[m], ls='', ms=7,
                           color=COLORS[m], mec=COLORS[m], mew=0.8,
                           label=LABELS[m]) for m in MODELS]
    fig.legend(handles=handles, loc='lower center', ncol=len(handles),
               bbox_to_anchor=(0.5, 0.0), frameon=False, handlelength=1.0,
               handletextpad=0.4, columnspacing=1.2, fontsize=FONT_S)

    base = os.path.splitext(OUT)[0]
    fig.savefig(base + '.png', dpi=600)                    # raster, 600 dpi
    fig.savefig(base + '.pdf')                             # vector
    fig.savefig(base + '.tiff', dpi=600,
                pil_kwargs={'compression': 'tiff_lzw'})    # production TIFF
    print('saved:', base + '.{png,pdf,tiff}')
    plt.close(fig)


if __name__ == '__main__':
    make_figure()