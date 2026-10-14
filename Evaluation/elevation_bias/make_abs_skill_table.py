#!/usr/bin/env python
"""
Reproducible builder for the absolute-unit skill table
(results/skill_table_bias_mae_rmse.tex + results/skill_abs_table.csv).

Why this script exists
----------------------
skill_summary.csv reports precipitation (pr) and wind (wind) errors as
percentages relative to each station's own reference annual. For
precipitation the per-station "annual" in the pooling pipeline is the mean of
the four seasonal SUMS, i.e. roughly a quarter of the true annual, so the %
denominator is unit-inconsistent and cannot be rescaled back to mm/day.
Wind errors in skill_summary.csv are already stored in m/s despite the '%'
unit label, but for a single consistent table we recompute everything here.

Method
------
- Land stations only (164 of the 175 in pickles/skill_seasonal.pkl; the 11
  buoy stations are excluded because the CDO M10 extractions are
  land-only). 1986-2005, sign convention model - obs.
- Temperature (T): PINNED from results/skill_summary.csv (validated
  reference; D03: bias -1.03, MAE 1.46, RMSE 1.64, n=121).
- Precipitation (P, mm/day) and wind speed (W, m/s): recomputed directly in
  absolute units from annual means (model - obs, per station, then averaged
  over stations).
  * WRF: data/wrf_stations/{t,pr,wind}_<dom>_*.nc; T2 -> C, pr daily
    increments clipped at 0 summed and /365, wspd mean.
  * CDO M10: data/cdo_extractions/canrcm4_noaa113 and
    canesm2_raw_noaa113; tas -> C, pr kg m-2 -> mm/day (*86400/365),
    sfcWind m/s.
- Observed annuals are read via seasonal_pooling.obs_pool and cached to
  results/_obs_annuals_abs_cache.pkl (obs wind reads ECCC hourly files, so
  the cache is important).

Outputs
-------
- results/skill_abs_table.csv          (machine-readable values)
- results/skill_table_bias_mae_rmse.tex (the LaTeX table)
"""
import datetime
import json
import os
import pickle
import sys
import time

import numpy as np
import pandas as pd
from cftime import num2date
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(BASE, 'data')
RESULTS = os.path.join(HERE, 'results')
OBS_CACHE = os.path.join(RESULTS, '_obs_annuals_abs_cache.pkl')
CSV_OUT = os.path.join(RESULTS, 'skill_abs_table.csv')
TEX_OUT = os.path.join(RESULTS, 'skill_table_bias_mae_rmse.tex')

sys.path.insert(0, HERE)
import seasonal_pooling as sp  # noqa: E402

DATASETS = ['D03', 'D02', 'D01', 'CanRCM4', 'CanESM2']
# Dataset keys: the WRF domains d03/d02/d01 are named D03/D02/D01 everywhere in
# this repo (pickles, skill_summary.csv) even though d02 is the 15 km nest and
# d01 the 75 km nest. Display names use the km convention.
LABEL = {
    'D03': 'CanESM2-WRF D3',
    'D02': 'CanESM2-WRF D15',
    'D01': 'CanESM2-WRF D75',
    'CanRCM4': 'CanRCM4',
    'CanESM2': 'CanESM2',
}
DOM = {'D03': 'd03', 'D02': 'd02', 'D01': 'd01'}
CDO_DIR = {'CanRCM4': 'canrcm4_noaa113', 'CanESM2': 'canesm2_raw_noaa113'}


def agency(sid):
    if sid.isdigit():
        return 'E'
    if sid.startswith(('USC', 'USW')):
        return 'N'
    if sid.isupper() and len(sid) == 3:
        return 'B'
    return None


def land_ids():
    skill = pickle.load(open(os.path.join(BASE, 'pickles', 'skill_seasonal.pkl'), 'rb'))
    return [str(s) for s in skill['D03']['station'] if agency(str(s)) in ('E', 'N', 'B')]


def m10_tag(sid, cid2wmo):
    if sid.startswith(('USC', 'USW')):
        return f'NOAA_{sid}'
    if sid.isupper() and len(sid) == 3:
        return f'BCH_{sid}'
    return f'ECCC_{cid2wmo.get(sid, sid)}'


def load_obs_annuals(ids):
    """Observed annual means per (station, var); pr in mm/day, t in C, wind in m/s."""
    try:
        with open(OBS_CACHE, 'rb') as fh:
            obs = pickle.load(fh)
        print(f'loaded obs cache: {len(obs)} entries')
        return obs
    except FileNotFoundError:
        pass
    obs = {}
    for i, sid in enumerate(ids):
        ag = agency(sid)
        for var in ('t', 'pr', 'wind'):
            r = sp.obs_pool(sid, ag, var, DATA)
            if r is None:
                continue
            ann = r['annual'].astype(float).copy()
            if var == 'pr':
                ann = ann / 365.0  # obs_pool returns mm/yr -> mm/day
            obs[(sid, var)] = ann
        if (i + 1) % 20 == 0:
            with open(OBS_CACHE, 'wb') as fh:
                pickle.dump(obs, fh)
            print(f'  obs computed for {i + 1}/{len(ids)} stations', flush=True)
    with open(OBS_CACHE, 'wb') as fh:
        pickle.dump(obs, fh)
    print(f'obs cache written: {len(obs)} entries')
    return obs


def wrf_files(prefix, dom):
    """Map station key (E/WMO or N/USC id or B/CCC) -> file path for one variable."""
    out = {}
    d = os.path.join(DATA, 'wrf_stations')
    bch3 = {'ALU', 'BCK', 'BLN', 'CLO', 'CMU', 'CMX', 'COQ', 'CQM', 'CRU',
            'DAI', 'DLU', 'ECL', 'GOC', 'MIS', 'NTY', 'STA', 'WAH', 'WOL'}
    for f in os.listdir(d):
        if not f.endswith('.nc') or not f.startswith(prefix + '_' + dom + '_'):
            continue
        rest = f[len(prefix + '_' + dom + '_'):-3]
        if dom == 'd01':
            if rest.startswith('ECCC_'):
                out['E' + rest[5:]] = os.path.join(d, f)
            elif rest.startswith('NOAA_'):
                out['N' + rest[5:]] = os.path.join(d, f)
            elif rest.startswith('BCH_'):
                out['B' + rest[4:]] = os.path.join(d, f)
            elif rest in bch3:
                out['B' + rest] = os.path.join(d, f)
        elif rest.startswith('st'):
            sid2 = rest[2:]
            out[('N' if sid2.startswith(('USC', 'USW')) else 'E') + sid2] = os.path.join(d, f)
    return out


def wrf_annual(path, var):
    """Annual means: T2 -> C; pr (running total) -> mm/day; wspd -> m/s."""
    nc = Dataset(path)
    t = nc.variables['time'][:]
    data = np.squeeze(nc.variables[var][:])
    nc.close()
    dates = pd.to_datetime(datetime.datetime(1986, 1, 1)) + pd.to_timedelta(t, unit='h')
    if var in ('T2', 'wspd'):
        h = pd.Series((data - 273.15) if var == 'T2' else data, index=dates)
        return h.groupby(h.index.year).mean()
    inc = np.clip(np.diff(data), 0.0, None)
    d = pd.Series(inc, index=dates[1:].normalize())
    return d.groupby(d.index.year).sum() / 365.0


def cdo_annual(dataset, var, sid, cid2wmo):
    vname = {'t': 'tas', 'pr': 'pr', 'wind': 'sfcWind'}[var]
    f = os.path.join(DATA, 'cdo_extractions', CDO_DIR[dataset],
                     f'{vname}_{m10_tag(sid, cid2wmo)}.nc')
    if not os.path.exists(f):
        return None
    nc = Dataset(f)
    v = nc.variables[vname][:].astype(float).ravel()
    tt = nc.variables['time'][:]
    u = nc.variables['time'].units
    c = nc.variables['time'].calendar
    nc.close()
    lo, hi = sp.M10_EVAL_SLICE[dataset]
    yrs = np.array([int(x.year) for x in num2date(tt[lo:hi], u, calendar=c)])
    vv = v[lo:hi]
    if var == 't':
        vv = vv - 273.15
    elif var == 'pr':
        vv = vv * 86400.0 / 365.0  # kg m-2 -> mm/day
    s = pd.Series(vv, index=yrs)
    return s.groupby(s.index).mean()


def skill_abs(model_ann, obs_ann):
    m = pd.Series(model_ann).astype(float).reindex(sp.YEARS)
    o = obs_ann.reindex(sp.YEARS).astype(float)
    sel = o.notna() & m.notna()
    if int(sel.sum()) < sp.MIN_YEARS:
        return None
    e = (m - o)[sel].values
    return dict(bias=float(np.mean(e)),
                mae=float(np.mean(np.abs(e))),
                rmse=float(np.sqrt(np.mean(e ** 2))),
                n=int(sel.sum()))


def mean_skill(rows):
    """(bias, mae, rmse, n_stations) averaged over station rows."""
    if not rows:
        return (np.nan, np.nan, np.nan, 0)
    return (float(np.mean([r['bias'] for r in rows])),
            float(np.mean([r['mae'] for r in rows])),
            float(np.mean([r['rmse'] for r in rows])),
            len(rows))


def t_from_summary():
    """Pinned T (C) rows from the validated results/skill_summary.csv."""
    df = pd.read_csv(os.path.join(RESULTS, 'skill_summary.csv'))
    t = {}
    for name in DATASETS:
        r = df[(df['dataset'] == name) & (df['var'] == 't')].iloc[0]
        t[name] = (float(r['bias']), float(r['mae']), float(r['rmse']), int(r['n_stations']))
    return t


def compute(dataset, ids, obs, cid2wmo):
    res = {'P': [], 'W': []}
    if dataset in DOM:
        dom = DOM[dataset]
        pf = wrf_files('pr', dom)
        wf = wrf_files('wspd', dom)
        for sid in ids:
            key = agency(sid) + sid
            for var, files, tag, out in (('pr', pf, 'pr', 'P'), ('wind', wf, 'wspd', 'W')):
                p = files.get(key)
                if p is None or (sid, var) not in obs:
                    continue
                try:
                    ma = wrf_annual(p, tag)
                except Exception:
                    continue
                r = skill_abs(ma, obs[(sid, var)])
                if r:
                    res[out].append(r)
    else:
        for sid in ids:
            for var, out in (('pr', 'P'), ('wind', 'W')):
                if (sid, var) not in obs:
                    continue
                try:
                    ma = cdo_annual(dataset, var, sid, cid2wmo)
                except Exception:
                    continue
                if ma is None:
                    continue
                r = skill_abs(ma, obs[(sid, var)])
                if r:
                    res[out].append(r)
    return {'P': mean_skill(res['P']), 'W': mean_skill(res['W'])}


def fmt_signed(x):
    return f'{x:+.2f}'


def fmt_val(x, signed=False):
    if signed:
        return fmt_signed(x)
    return f'{x:.2f}'


def main():
    t0 = time.time()
    os.makedirs(RESULTS, exist_ok=True)
    ids = land_ids()
    print(f'{len(ids)} land stations')

    eid = pd.read_csv(os.path.join(BASE, 'files', 'ECCC_d03_stations.csv'), header=None)
    cid2wmo = dict(zip(eid.iloc[:, 4].astype(str), eid.iloc[:, 3].astype(str)))

    obs = load_obs_annuals(ids)
    t_rows = t_from_summary()

    table = {}
    for name in DATASETS:
        pw = compute(name, ids, obs, cid2wmo)
        table[name] = {'T': t_rows[name], 'P': pw['P'], 'W': pw['W']}
        tb, tm, tr, tn = table[name]['T']
        pb, pm, pr_, pn = table[name]['P']
        wb, wm, wr, wn = table[name]['W']
        print(f'{LABEL[name]:16s} T[{tb:+.2f} {tm:.2f} {tr:.2f} n={tn}] '
              f'P[{pb:+.3f} {pm:.3f} {pr_:.3f} n={pn}] '
              f'W[{wb:+.3f} {wm:.3f} {wr:.3f} n={wn}]', flush=True)

    # ---- CSV
    rows = []
    for name in DATASETS:
        for var in ('T', 'P', 'W'):
            unit = {'T': 'C', 'P': 'mm/day', 'W': 'm/s'}[var]
            b, m, r, n = table[name][var]
            rows.append(dict(model=LABEL[name], var=var, unit=unit,
                             bias=b, mae=m, rmse=r, n_stations=n))
    pd.DataFrame(rows).to_csv(CSV_OUT, index=False)
    print(f'wrote {CSV_OUT}')

    # ---- LaTeX
    lines = [
        r'\begin{table}[!htbp]',
        r'\caption{\label{table:bias} Mean annual bias and error across all land '
        r'stations for',
        r'2\,m air temperature (T; $^\circ$C), precipitation (P; mm/day), and wind '
        r'speed (Wspd; m/s).',
        r'Bias $=$ model $-$ obs (signed); MAE and RMSE are the mean absolute and '
        r'root-mean-square error, each averaged over stations. Station counts: T, '
        r'120--121; P, 111--129; Wspd, 10 (wind-observed land stations only).}',
        r'\centering',
        r'\begin{tabular}{l ccc ccc ccc}',
        r'\hline',
        r' & \multicolumn{3}{c}{T ($^\circ$C)} & \multicolumn{3}{c}{P (mm/day)} & '
        r'\multicolumn{3}{c}{Wspd (m/s)} \\',
        r'\cline{2-4}\cline{5-7}\cline{8-10}',
        r'Model & Bias & MAE & RMSE & Bias & MAE & RMSE & Bias & MAE & RMSE \\',
        r'\hline',
    ]
    for name in DATASETS:
        tb, tm, tr, tn = table[name]['T']
        pb, pm, pr_, pn = table[name]['P']
        wb, wm, wr, wn = table[name]['W']
        lines.append(
            f'{LABEL[name]:16s} & {fmt_val(tb, True)} & {fmt_val(tm)} & {fmt_val(tr)}'
            f' & {fmt_val(pb, True)} & {fmt_val(pm)} & {fmt_val(pr_)}'
            f' & {fmt_val(wb, True)} & {fmt_val(wm)} & {fmt_val(wr)} \\\\')
    lines += [r'\hline', r'\end{tabular}', r'\end{table}', '']
    with open(TEX_OUT, 'w') as fh:
        fh.write('\n'.join(lines))
    print(f'wrote {TEX_OUT}')
    print(f'total {round(time.time() - t0, 1)} s')


if __name__ == '__main__':
    main()
