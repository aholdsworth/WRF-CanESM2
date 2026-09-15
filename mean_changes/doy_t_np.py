import numpy as np, pandas as pd, os
from netCDF4 import Dataset
BASE='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/WRF_FILES/cf_compliant/data/'
OUT='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS/WRF-CanESM2/mean_changes/output_local/'
geo='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS/WRF-CanESM2/geo_em.d03.nc'
stem='CanESM2-WRF_t2_d03'

g=Dataset(geo); lm=np.asarray(g.variables['LANDMASK'][:]); g.close()
land=(lm==1); ocean=(lm==0)

def seasonal(path):
    nc=Dataset(path)
    t=np.ma.filled(np.ma.asarray(nc.variables['T2'][:],dtype=np.float32),np.nan)
    ts=np.asarray(nc.variables['time'][:])
    lat=np.asarray(nc.variables['lat'][:],dtype=np.float64)
    units=nc.variables['time'].units; cal=nc.variables['time'].calendar
    nc.close()
    import cftime
    times=cftime.num2date(ts, units, calendar=cal)
    keep=np.array([not (dt.month==2 and dt.day==29) for dt in times])
    doy=np.array([dt.timetuple().tm_yday for dt in times[keep]],dtype=np.int32)
    T=t[keep]
    D=(doy-1).astype(np.intp)
    S=np.full((365,300,300),np.nan)
    counts=np.bincount(D,minlength=366).astype(np.float64)
    tmp=np.full((366,300,300),0.0)
    np.add.at(tmp, D, T)
    S=tmp[:365]/counts[:365,None,None]
    imax=S.argmax(axis=0)+1; imin=S.argmin(axis=0)+1
    return imax,imin,lat

import datetime
hi_max,hi_min,lahist=seasonal(BASE+f'historical/{stem}_historical_daily_1986_2005_300x300.nc')
w=np.cos(np.deg2rad(lahist))
def circ(f,h,c=365):
    d=f-h; return ((d+c/2)%c)-c/2
rows=[]
for scen in ['rcp45','rcp85']:
    fm,fmmin,la=seasonal(BASE+f'{scen}/{stem}_{scen}_daily_2046_2065_300x300.nc')
    dmax, dmin=circ(fm,hi_max), circ(fmmin,hi_min)
    for region,m in [('Land',land),('Ocean',ocean)]:
        den=np.nansum(np.where(m,w,0.0))
        num=lambda a: np.nansum(np.where(m,a.astype(np.float64),np.nan)*w)/den
        rows.append(dict(scenario=scen,region=region,
            hist_day_max=num(hi_max),fut_day_max=num(fm),delta_day_max=num(dmax),
            hist_day_min=num(hi_min),fut_day_min=num(fmmin),delta_day_min=num(dmin)))
df=pd.DataFrame(rows)
df.to_csv(OUT+'doy_T2_land_ocean.csv',index=False)
print(df.to_string(index=False))
