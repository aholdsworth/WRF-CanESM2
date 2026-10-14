#!/usr/bin/env python
"""Compatibility PROJECTION of the canonical Rule-1 land table.

Since 2026-09-29 (unify-obs-pipeline) `make_station_tables.py` is the single
root producer: it computes obs through `seasonal_pooling.obs_pool` (Rule 1:
90% annual gate, period only if >= MIN_YEARS=15 of 20 years survive) and
skill through `seasonal_pooling.skill` (model masked to the obs surviving
years).  This script no longer recomputes anything — it projects
`pickles/station_tables.pkl` into the legacy
`pickles/skill_seasonal.pkl` / `skill_summary.pkl` (and .csv) column names so
`make_bias_figure.py`, `make_elev_bias_grid.py` and `make_abs_skill_table.py`
run unchanged.

skill_seasonal.pkl : dict dataset -> DataFrame (one row per station), columns
    station, lat, lon, obs_elev_m, elevation_bias_m,
    per var in {t, pr, wind}:
        {var}_bias        (t in C; pr & wind in % — identical values to the
                           table's temperature_bias_C / precipitation_bias_pct
                           / wind_bias_pct)
        {var}_mae, {var}_rmse, {var}_r2, {var}_n_years, {var}_obs_period
    plus a top-level '_source' tag naming the upstream table version.

skill_summary.pkl / .csv : 15 rows (5 datasets x 3 vars):
    dataset, var, unit, bias, mae, rmse, r2, n_stations, n_years_min,
    n_years_max

Run:  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \\
        /home/amh001/space_fs7/software_2022/python/py_2024/bin/python \\
        make_skill_tables.py
"""
import os
import pickle
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

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
TABLE_PKL = os.path.join(BASE, 'pickles', 'station_tables.pkl')
OUT_SKILL = os.path.join(BASE, 'pickles', 'skill_seasonal.pkl')
OUT_SUM   = os.path.join(BASE, 'pickles', 'skill_summary.pkl')
OUT_CSV   = os.path.join(BASE, 'pickles', 'skill_summary.csv')

DATASETS = ['D03', 'D02', 'D01', 'CanESM2', 'CanRCM4']
# table column -> legacy skill column, per var
VAR_COLS = {'t': 'temperature_bias_C', 'pr': 'precipitation_bias_pct',
            'wind': 'wind_bias_pct'}


def project(name, df):
    """Project one dataset frame of station_tables.pkl to legacy columns."""
    out = pd.DataFrame(index=df.index)
    for col in ('station', 'lat', 'lon', 'obs_elev_m', 'elevation_bias_m'):
        out[col] = df[col].values
    for var in ('t', 'pr', 'wind'):
        out[f'{var}_bias'] = df[VAR_COLS[var]].values
        for k in ('mae', 'rmse', 'r2', 'n_years', 'obs_period'):
            out[f'{var}_{k}'] = df[f'{var}_{k}'].values
    # legacy skill frames carried no index name of their own
    out.index.name = None
    return out


def main():
    with open(TABLE_PKL, 'rb') as f:
        table = pickle.load(f)
    src_version = table.get('_version', 'unknown')
    if not all(k in table for k in DATASETS):
        sys.exit(f'ERROR: {TABLE_PKL} missing frames {DATASETS}; rebuild '
                 'make_station_tables.py first')

    tables = {name: project(name, table[name]) for name in DATASETS}

    with open(OUT_SKILL, 'wb') as f:
        pickle.dump(tables, f)
    print(f'Saved {OUT_SKILL} (projected from {TABLE_PKL}, '
          f'source version {src_version})')

    # ---- domain summary (same schema as before) ----
    sum_rows = []
    for name in DATASETS:
        df = tables[name]
        for var, unit in (('t', 'C'), ('pr', '%'), ('wind', '%')):
            sub = df[[f'{var}_bias', f'{var}_mae', f'{var}_rmse', f'{var}_r2',
                      f'{var}_n_years']].dropna(subset=[f'{var}_bias'])
            if len(sub):
                ny = sub[f'{var}_n_years']
                sum_rows.append({
                    'dataset': name, 'var': var, 'unit': unit,
                    'bias': sub[f'{var}_bias'].mean(),
                    'mae': sub[f'{var}_mae'].mean(),
                    'rmse': sub[f'{var}_rmse'].mean(),
                    'r2': sub[f'{var}_r2'].mean(),
                    'n_stations': len(sub),
                    'n_years_min': int(ny.min()), 'n_years_max': int(ny.max()),
                })
            else:
                sum_rows.append({'dataset': name, 'var': var, 'unit': unit,
                                 'bias': np.nan, 'mae': np.nan, 'rmse': np.nan,
                                 'r2': np.nan, 'n_stations': 0,
                                 'n_years_min': np.nan, 'n_years_max': np.nan})
    summary = pd.DataFrame(sum_rows,
                           columns=['dataset', 'var', 'unit', 'bias', 'mae',
                                    'rmse', 'r2', 'n_stations', 'n_years_min',
                                    'n_years_max'])
    with open(OUT_SUM, 'wb') as f:
        pickle.dump(summary, f)
    summary.to_csv(OUT_CSV, index=False)
    print(f'Saved {OUT_SUM} and {OUT_CSV}')

    # ---- console report ----
    print('\n=== Domain summary (Rule-1 canonical land table) ===')
    print(summary.to_string(index=False,
                            formatters={'bias': '{:.3f}'.format,
                                        'mae': '{:.3f}'.format,
                                        'rmse': '{:.3f}'.format,
                                        'r2': '{:.3f}'.format}))

    print('\n=== Per-dataset non-NaN bias counts ===')
    for name in DATASETS:
        df = tables[name]
        print(f"{name:9s} t={df['t_bias'].notna().sum():3d} "
              f"pr={df['pr_bias'].notna().sum():3d} "
              f"wind={df['wind_bias'].notna().sum():3d}")

    # ---- identity check: projection must equal the table ----
    bad = []
    for name in DATASETS:
        src = table[name]
        dst = tables[name]
        for var in ('t', 'pr', 'wind'):
            for k, tcol in (('bias', VAR_COLS[var]),
                            ('mae', f'{var}_mae'),
                            ('rmse', f'{var}_rmse'),
                            ('r2', f'{var}_r2'),
                            ('n_years', f'{var}_n_years'),
                            ('obs_period', f'{var}_obs_period')):
                a, b = src[tcol].values, dst[f'{var}_{k}'].values
                if not np.allclose(a, b, equal_nan=True):
                    bad.append((name, var, k))
    if bad:
        sys.exit(f'ERROR: projection identity check failed: {bad}')
    print('\nIdentity check passed: all projected columns == table columns.')


if __name__ == '__main__':
    main()
