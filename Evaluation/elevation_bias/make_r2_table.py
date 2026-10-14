#!/usr/bin/env python3
"""make_r2_table.py v2 — compute Table-1 R^2 (both definitions) FAST.

Phase A (once): pool obs annuals (Rule-1 gate) for all 173 land x 3 vars and
cache to pickles/_obs_annual_cache.pkl  (the ~20-min part).
Phase B (seconds, cached): build each model's annuals and compute, per
(model, var):
  st_mean_r2 : mean of per-station corr(model_annual, obs_annual) over shared yrs
  pooled_r2  : corr over ALL pooled (station x year) pairs, squared
Writes pickles/r2_table.pkl = {model: {var: {st_mean_r2, st_n, pooled_r, pooled_r2, pooled_n}}}

Run with OPENBLAS_NUM_THREADS=1 (avoid 129-thread corrcoef spin):
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python make_r2_table.py
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python make_r2_table.py --cached
"""
import os, sys, pickle, datetime, warnings
import numpy as np, pandas as pd
from netCDF4 import Dataset
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
DATA = os.path.join(BASE, 'data')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # seasonal_pooling.py / station_geo.py (same dir)
import seasonal_pooling as sp
WRF_ST_DIR = os.path.join(DATA, 'wrf_canon')
CACHE = os.path.join(BASE, 'pickles', '_obs_annual_cache.pkl')

def _loaders():
    from station_geo import load_eccc, load_bch, load_noaa
    land = [str(i) for i in load_eccc()[0]] + [str(i) for i in load_bch()[0]] + list(load_noaa('NOAA_d03_stations.csv')[0])
    return land
def wrf_st_files(var, dom):
    out = {}
    d = WRF_ST_DIR
    BCH3 = {'ALU','BCK','BLN','CLO','CMU','CMX','COQ','CQM','CRU','DAI','DLU','ECL',
            'GOC','MIS','NTY','STA','WAH','WOL'}
    for f in os.listdir(d):
        if not f.endswith('.nc') or not f.startswith(var + '_' + dom + '_'):
            continue
        rest = f[len(var + '_' + dom + '_'):-3]
        # wrf_canon naming: {var}_{dom}_{AGENCY}_{ID}.nc (ECCC_ / BCH_ / NOAA_)
        if rest.startswith('ECCC_'):
            out['E' + rest[5:]] = os.path.join(d, f)
        elif rest.startswith('BCH_'):
            out['B' + rest[4:]] = os.path.join(d, f)
        elif rest.startswith('NOAA_'):
            out['N' + rest[5:]] = os.path.join(d, f)
        else:
            import sys
            print(f'WARN wrf_st_files: untagged {dom} file skipped: {f}', file=sys.stderr)
    return out
def wrf_station_daily(path, var_name, kind):
    """Load one per-station WRF file -> DAILY Series (t: C, pr: mm/day,
    wind: m/s).  pr uses diff of the cumulative field, clipped >= 0."""
    nc = Dataset(path, mode='r')
    t = nc.variables['time'][:]
    data = np.squeeze(nc.variables[var_name][:])
    nc.close()
    dates = pd.to_datetime(datetime.datetime(sp.START_YEAR, 1, 1)) \
              + pd.to_timedelta(t, unit='h')
    if kind == 't':
        h = pd.Series(data - 273.15, index=dates)
        return h.groupby(h.index.normalize()).mean()
    if kind == 'wind':
        h = pd.Series(data, index=dates)
        return h.groupby(h.index.normalize()).mean()
    # pr: cumulative -> increments, clip <0, drop the first (no prior value)
    inc = np.clip(np.diff(data), 0.0, None)
    d = pd.Series(inc, index=dates[1:].normalize())
    return d.groupby(d.index).sum()
def wrf_model_annuals(ids, var, dom):
    """dict sid -> Series(annual values, index=years) for WRF per-station files."""
    files = wrf_st_files(var, dom)
    out = {}
    for sid in ids:
        key = agency_of(sid) + sid
        p = files.get(key)
        if p is None:
            continue
        try:
            daily = wrf_station_daily(p, {'t': 'T2', 'pr': 'pr', 'wind': 'wspd'}[var], var)
            out[sid] = sp.model_annual_from_series(daily, 'sum' if var == 'pr' else 'mean')
        except Exception as e:
            print(f'  WARN WRF {dom} {var} {sid}: {e}', flush=True)
    return out

def canesm2_model_annuals(var):
    from cftime import num2date
    vname = {'t':'tas','pr':'pr','wind':'sfcWind'}[var]
    lo, hi = sp.M10_EVAL_SLICE['canesm2']; out = {}
    for sid in land_ids:
        f = os.path.join(DATA, sp.M10_CanESM2_DIR, f'{vname}_{sp._m10_file_tag(sid)}.nc')
        if not os.path.exists(f): continue
        try:
            nc = Dataset(f, mode='r')
            v = nc.variables[vname][:].astype(float).ravel(); t = nc.variables['time'][:]
            units, cal = nc.variables['time'].units, nc.variables['time'].calendar; nc.close()
            dates = num2date(t[lo:hi], units, calendar=cal)
            s = pd.Series(v[lo:hi], index=pd.DatetimeIndex([pd.Timestamp(int(d.year),int(d.month),int(d.day)) for d in dates]))
            if var=='t': s = s-273.15
            elif var=='pr': s = s*86400.0
            out[sid] = sp.model_annual_from_series(s, 'sum' if var=='pr' else 'mean')
        except Exception as e:
            print(f'  WARN CanESM2 {var} {sid}: {e}', flush=True)
    return out
def canrcm4_model_annuals(var):
    out = {}
    for sid in land_ids:
        try: s = sp.canrcm4_daily_series(sid, var, DATA)
        except Exception: s = None
        if s is None or len(s.dropna())==0: continue
        out[sid] = sp.model_annual_from_series(s, 'sum' if var=='pr' else 'mean')
    return out

land_ids = _loaders()
def agency_of(sid):
    if sid.isdigit(): return 'E'
    if sid.startswith(('USC','USW')): return 'N'
    return 'B'

def st_r2(model_annual, o):
    if model_annual is None: return None
    m = model_annual.reindex(list(sp.YEARS)).astype(float)
    ok = o.notna() & m.notna()
    if int(ok.sum()) < sp.MIN_YEARS: return None
    r = np.corrcoef(m[ok].values, o[ok].values)[0,1]
    return float(r) if np.isfinite(r) else None

def pooled(model_map, o):
    Mv, Ov = [], []
    for sid in land_ids:
        m = model_map.get(sid)
        if m is None: continue
        mm = m.reindex(list(sp.YEARS)).astype(float)
        ok = o[sid].notna() & mm.notna()
        if int(ok.sum()) < sp.MIN_YEARS: continue
        Mv.extend(mm[ok].values); Ov.extend(o[sid][ok].values)
    if len(Mv) < 3: return float('nan'), 0
    r = np.corrcoef(Mv, Ov)[0,1]
    return (float(r), float(r*r), len(Mv)) if np.isfinite(r) else (float('nan'), float('nan'), len(Mv))

def model_annuals(name):
    anns = {}
    if name in ('D03','D02','D01'):
        dom = {'D03':'d03','D02':'d02','D01':'d01'}[name]
        for var in ('t','pr','wind'):
            print(f'  {name} {var}...', flush=True); anns[var] = wrf_model_annuals(land_ids, var, dom)
    elif name == 'CanESM2':
        for var in ('t','pr','wind'):
            print(f'  CanESM2 {var}...', flush=True)
            try: anns[var] = canesm2_model_annuals(var)
            except Exception as e: print(f'  WARN {e}', flush=True); anns[var]={}
    elif name == 'CanRCM4':
        for var in ('t','pr','wind'):
            print(f'  CanRCM4 {var}...', flush=True); anns[var] = canrcm4_model_annuals(var)
    return anns

# ---- Phase A: obs annuals (cached) ----
if '--cached' in sys.argv and os.path.exists(CACHE):
    print('loading obs cache', flush=True)
    obs_ann = pickle.load(open(CACHE,'rb'))   # {var: {sid: Series}}
    obs_ny  = pickle.load(open(CACHE.replace('.pkl','_ny.pkl'),'rb')) if os.path.exists(CACHE.replace('.pkl','_ny.pkl')) else None
else:
    print('Phase A: obs pooling (Rule-1)...', flush=True)
    obs_ann, obs_ny = {}, {}
    for var in ('t','pr','wind'):
        obs_ann[var], obs_ny[var] = {}, {}
    for i, sid in enumerate(land_ids):
        ag = agency_of(sid)
        for var in ('t','pr','wind'):
            rec = sp.obs_pool(sid, ag, var, DATA)
            if rec is None:
                obs_ann[var][sid] = pd.Series(np.nan, index=list(sp.YEARS)); obs_ny[var][sid] = 0
            else:
                obs_ann[var][sid] = rec['annual'].reindex(list(sp.YEARS))
                obs_ny[var][sid] = int(rec['n_years'])
        if (i+1) % 40 == 0: print(f'  obs {i+1}/{len(land_ids)}', flush=True)
    pickle.dump(obs_ann, open(CACHE,'wb'))
    pickle.dump(obs_ny, open(CACHE.replace('.pkl','_ny.pkl'),'wb'))
    print('saved obs cache', flush=True)

def ok_years(sid, var):
    o = obs_ann[var][sid]
    ny = obs_ny[var][sid] if obs_ny is not None else int(o.notna().sum())
    return (o, ny >= sp.MIN_YEARS)

# ---- Phase B: models + r2 ----
results = {}
for name in ('D03','D02','D01','CanESM2','CanRCM4'):
    print(f'dataset {name}', flush=True)
    ma = model_annuals(name)
    per = {}
    for var in ('t','pr','wind'):
        rs = [r for r in (st_r2(ma.get(var,{}).get(sid), obs_ann[var][sid]) for sid in land_ids
                          if (obs_ny[var][sid] if obs_ny is not None else obs_ann[var][sid].notna().sum() >= sp.MIN_YEARS))
              if r is not None]
        pr, pr2, pn = pooled(ma.get(var,{}), obs_ann[var])
        per[var] = dict(st_mean_r2=(float(np.mean(rs)) if rs else float('nan')), st_n=len(rs),
                        pooled_r=pr, pooled_r2=pr2, pooled_n=pn)
    results[name] = per
    print('  ' + '  '.join(f"{v}: st_mean={per[v]['st_mean_r2']:.3f}(n={per[v]['st_n']}) pooled_r2={per[v]['pooled_r2']:.3f}(n={per[v]['pooled_n']})" for v in ('t','pr','wind')), flush=True)

pickle.dump(results, open(os.path.join(BASE,'pickles','r2_table.pkl'),'wb'))
print('Saved pickles/r2_table.pkl')
