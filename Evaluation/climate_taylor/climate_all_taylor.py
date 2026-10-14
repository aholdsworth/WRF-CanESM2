#!/usr/bin/env python3
"""climate_all_taylor.py -- 3-panel "climate" Taylor diagram
(T / PR / wind) on 20-yr (1986-2005) station climatologies,
land + ocean buoys.  Same figure as climate_all_taylor.ipynb.
"""
# Exactly the taylor_3panel figure, but each station is reduced to its
# 20-year annual-mean climatology before pooling (T / PR / wind panels).
import os, pickle, warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import taylor_figure as tf
warnings.filterwarnings('ignore')

# Path to the per-station x per-year obs/model annual arrays (v10).
# Regenerate in the main verification workspace (NOTEBOOKS/) via
#   python taylor_data.py
# which reads the raw WRF/CanESM2/CanRCM4 + obs data and writes the cache.
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
PKL = os.environ.get(
    'TAYLOR_ANNUALS_PKL',
    os.path.join(BASE, 'pickles', 'taylor_annuals.pkl'))
with open(PKL, 'rb') as f:
    TA = pickle.load(f)
print('taylor cache version:', TA['_version'])
MODELS = tf.MODELS
def climo(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return a.mean() if a.size else np.nan
# build (obs_climo, models_climo) per var over the FULL 3-panel station set
DATA = {}
for v in ('t', 'pr', 'wind'):
    obs, models = TA[v]
    o = {s: climo(obs[s]) for s in obs}
    m = {mm: {s: climo(models[mm][s]) for s in models[mm]} for mm in MODELS}
    DATA[v] = (o, m)
    print(v, 'obs stations =', len(o))


# Taylor inputs on climatology values (per model over its shared stations)
def climo_refstd(obs, models):
    keep = set()
    for m in MODELS:
        keep |= set(obs) & set(models.get(m, {}))
    vals = [obs[s] for s in sorted(keep) if np.isfinite(obs[s])]
    return float(np.std(vals))
def climo_stats(obs, models, label):
    m = models.get(label, {})
    sh = [s for s in sorted(set(obs) & set(m))
          if np.isfinite(obs[s]) and np.isfinite(m[s])]
    if len(sh) < 3:
        return np.nan, np.nan, np.nan, len(sh)
    o = np.array([obs[s] for s in sh]); mm = np.array([m[s] for s in sh])
    return float(np.std(mm)), float(np.corrcoef(o, mm)[0, 1]), None, len(sh)
print('%-4s %-8s %5s %9s %9s %7s' % ('var', 'model', 'n', 'std', 'refstd', 'r'))
for v in ('t', 'pr', 'wind'):
    o, m = DATA[v]
    ref = climo_refstd(o, m)
    for lab in MODELS:
        std, corr, _, n = climo_stats(o, m, lab)
        print('%-4s %-8s %5d %9.3f %9.3f %7.3f' % (v, lab, n, std, ref, corr))


# ---- Draw the 3-panel figure (mirrors tf.make_figure layout, climatology data) ----
FIG_W = tf.FIG_W_MM / 25.4
FIG_H = tf.FIG_H
def draw_panel(fig, rect, data, var, refstd, max_std, levels):
    obs, models = data
    stats = {}
    for mm in MODELS:
        std, corr, _, n = climo_stats(obs, models, mm)
        stats[mm] = (std, corr, n)
    dia = tf.TaylorDiagram(refstd, fig=fig, rect=rect, label='Observations',
                           max_std=max_std, unit=tf.VAR_UNITS[var])
    dia.put_ticks_out()
    ct = dia.add_contours(levels, colors='0.55', linewidths=0.9, zorder=0)
    # label_contours(corr=0.55) labels each level where it crosses the
    # shared corr=0.55 diagonal but silently skips levels whose contour
    # never reaches that diagonal.  In the wind panel (refstd 2.09,
    # axis to 2.5) even the top contour only spans corr in
    # [0.087, 0.996], so the shared diagonal misses EVERY level and the
    # panel is left unlabeled.  Fallback: label where each contour
    # crosses the mid-radius ring r = (smin+smax)/2 (guaranteed for any
    # level within [refstd - r_mid, refstd + r_mid]); a level too close
    # to the refstd star to intersect the ring (disc < 0 there) is
    # placed at the contour's own maximum-RMS point (corr = 0) instead.
    # White halo keeps the labels legible on the grey lines.
    if not dia.label_contours(ct, corr=0.55):
        # label_contours() re-clamps the angle with
        # a = max(a, arcsin(l/refstd) + 0.08), so a per-level
        # corr arg does not reliably place the label on the
        # mid-radius ring.  Place each label at the point where the
        # contour crosses r = r_mid (the root closest to mid-radius),
        # or, for a level too close to the refstd star to reach the
        # ring, at the contour's maximum-RMS point (corr = 0).
        r_mid = 0.5 * (dia.smin + dia.smax)
        # 1 decimal: the step is 0.5, and the integer format would
        # label the 0.5 contour "0" (ambiguous with no zero contour).
        fmt = '%.1f'
        for lvl in ct.levels:
            arg = (dia.refstd ** 2 + r_mid ** 2 - lvl ** 2) / (
                2.0 * dia.refstd * r_mid)
            a = np.arccos(np.clip(arg, -1.0, 1.0))
            b = -2.0 * dia.refstd * np.cos(a)
            disc = b * b - 4.0 * (dia.refstd ** 2 - lvl ** 2)
            if disc < 0:
                continue
            roots = [(-b - np.sqrt(disc)) / 2.0,
                     (-b + np.sqrt(disc)) / 2.0]
            roots = [r for r in roots
                     if dia.smin + 1e-9 <= r <= dia.smax - 1e-9]
            if not roots:
                continue
            rr = min(roots, key=lambda r: abs(r - r_mid))
            dia.ax.text(a, rr, fmt % lvl, ha='center', va='center',
                        fontsize=tf.FONT_S, color='0.35', zorder=500,
                        bbox=dict(boxstyle='round,pad=0.12', fc='white',
                                  ec='none', alpha=0.8, lw=0))
    for mm in MODELS:
        std, corr, n = stats[mm]
        if not (np.isnan(std) and np.isnan(corr)):
            dia.add_sample(std, corr, marker=tf.MARKERS[mm], ms=7,
                           color=tf.COLORS[mm], mec=tf.COLORS[mm], mew=0.8,
                           zorder=3, label=tf.LABELS[mm])
    dia._cosmetics()

fig = plt.figure(figsize=(FIG_W, FIG_H))
x0 = 0.065; panel_w = 0.235; panel_gap = 0.095
panel_h = panel_w * FIG_W / FIG_H
panel_y = 0.27
rects = [[x0 + i * (panel_w + panel_gap), panel_y, panel_w, panel_h] for i in range(3)]
for i, (var, title, unit) in enumerate(tf.VARS):
    o, m = DATA[var]
    refstd = climo_refstd(o, m)
    fin = [climo_stats(o, m, mm)[0] for mm in MODELS]
    fin = [s for s in fin if not np.isnan(s)]
    max_std, _ = tf._round_axis_ceiling(max([refstd] + fin) * 1.15, n_intervals=5)
    levels = tf._auto_levels(max_std)
    draw_panel(fig, rects[i], DATA[var], var, refstd, max_std, levels)
    fig.text(rects[i][0] + rects[i][2] / 2.0, panel_y + panel_h + 0.135, title,
             ha='center', va='bottom', fontsize=tf.FONT)
handles = [plt.Line2D([], [], marker='*', ls='', ms=12, color='black',
                      mec='black', mew=0.8, label='Observations (climatology)')]
handles += [plt.Line2D([], [], marker=tf.MARKERS[mm], ls='', ms=7,
                       color=tf.COLORS[mm], mec=tf.COLORS[mm], mew=0.8,
                       label=tf.LABELS[mm]) for mm in MODELS]
fig.legend(handles=handles, loc='lower center', ncol=len(handles),
           bbox_to_anchor=(0.5, 0.0), frameon=False, handlelength=1.0,
           handletextpad=0.4, columnspacing=1.2, fontsize=tf.FONT_S)
OUT_DIR = os.path.join(BASE, 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
OUTF = os.path.join(OUT_DIR, 'climate_all_taylor.png')
fig.savefig(OUTF, dpi=600)
fig.savefig(OUTF.replace('.png', '.pdf'))
print('saved', OUTF)
plt.close(fig)

