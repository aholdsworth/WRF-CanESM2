#!/usr/bin/env python
"""
make_buoy_wind_obs.py -- extract the buoy wind OBS annual series (10 m wind
speed, m/s, 1986-2005) for the Taylor figure.

Source: pickles/buoy_wind_tables.pkl (the verified 2026-09-22
buoy_wind_bias.py run, relaxed year-gate).  per_buoy['<ID>|<model>']
['annual_obs'] holds the 20-yr observed annual means (identical across
models for a given buoy -- obs are model-independent).  The raw .mat files
(amh001/DATA/BuoyData) are not available in this sandbox, so we use the
verified pickle values directly rather than reconstructing hourly data.

Output: pickles/buoy_wind_obs_annuals.pkl -> {id: np.array(20) or None}
(None = buoy has no obs years 1986-2005 -> excluded from the Taylor wind
panel, e.g. 46005 / 46050 / 46087 / 46088 / 46089).

Run: /home/amh001/space_fs7/software_2022/python/py_2024/bin/python make_buoy_wind_obs.py
"""
import os
import pickle
import numpy as np

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
PKL_IN = os.path.join(BASE, 'pickles', 'buoy_wind_tables.pkl')
PKL_OUT = os.path.join(BASE, 'pickles', 'buoy_wind_obs_annuals.pkl')

# buoys with a verified D01 WRF wind file in data/wrf_canon/
D01_WIND_BUOYS = ['46029', '46041', '46131', '46132', '46134',
                  '46146', '46204', '46206', '46207']


def main():
    data = pickle.load(open(PKL_IN, 'rb'))
    per_buoy = data['per_buoy']
    out = {}
    for bid in sorted(set(k.split('|')[0] for k in per_buoy)):
        anns = [per_buoy[k]['annual_obs']
                for k in per_buoy
                if k.split('|')[0] == bid and 'annual_obs' in per_buoy[k]]
        if not anns:
            out[bid] = None
            continue
        # all model entries carry the same obs series; sanity-check equality
        a0 = np.asarray(anns[0], dtype=float)
        for a in anns[1:]:
            a = np.asarray(a, dtype=float)
            ok = np.allclose(a0, a, equal_nan=True) or np.all(
                np.isnan(a) & np.isnan(a0))
            if not ok:
                raise SystemExit('obs annuals disagree for buoy %s' % bid)
        out[bid] = a0
        n = int(np.isfinite(a0).sum())
        tag = ' (D01 wind)' if bid in D01_WIND_BUOYS else ''
        print('%s: n_obs_years=%2d  mean=%s%s'
              % (bid, n, ('%.3f' % np.nanmean(a0)) if n else '  n/a', tag))
    pickle.dump(out, open(PKL_OUT, 'wb'))
    print('\nsaved: %s (%d buoys, %d with obs years)'
          % (PKL_OUT, len(out), sum(1 for v in out.values()
                                    if v is not None and np.isfinite(v).any())))


if __name__ == '__main__':
    main()
