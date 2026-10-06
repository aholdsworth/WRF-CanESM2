#!/usr/bin/env python
"""1986-2005 gridded D03 climatology from the CF-compliant WRF files.

Output: pickles/climo_d03_hist.pkl
  {var: dict(lat=2d, lon=2d, clim=2d)}  var in ('t','pr','wind')
  t    : daily-mean T2 -> 20-yr mean, degC
  pr   : daily-sum pr  -> mean daily total, mm/day
  wind : daily-mean wspd -> 20-yr mean, m/s
All three are simple time means of the daily values (uniform 1986-2005).
"""
import os, pickle, time
import numpy as np
import xarray as xr

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = '/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/WRF_FILES/cf_compliant/data/historical/'
OUT = os.path.join(BASE, 'pickles', 'climo_d03_hist.pkl')

FILES = {'t':    ('CanESM2-WRF_t2_d03_historical_daily_1986_2005_300x300.nc',    'T2'),
         'pr':   ('CanESM2-WRF_pr_d03_historical_daily_1986_2005_300x300.nc',    'pr'),
         'wind': ('CanESM2-WRF_wspd_d03_historical_daily_1986_2005_300x300.nc',  'wspd')}

out = {}
for var, (fname, vname) in FILES.items():
    t0 = time.time()
    ds = xr.open_dataset(os.path.join(SRC, fname))
    x = ds[vname]
    if var == 't':
        x = x - 273.15
    clim = np.asarray(x.mean(dim='time'), dtype=np.float32)
    lat, lon = (np.asarray(ds[k], dtype=np.float32) for k in ('lat', 'lon'))
    out[var] = dict(lat=lat, lon=lon, clim=np.ascontiguousarray(clim))
    ds.close()
    c = out[var]['clim']
    print('%-5s min=%.3f max=%.3f p50=%.3f   (%.0f s)' %
          (var, c.min(), c.max(), np.nanpercentile(c, 50), time.time() - t0))

with open(OUT, 'wb') as f:
    pickle.dump(dict(version=1, source=SRC, out=out), f, protocol=4)
print('wrote', OUT, '%.1f MB' % (os.path.getsize(OUT) / 1e6))
