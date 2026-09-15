import numpy as np, xarray as xr, pandas as pd, os
BASE='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/WRF_FILES/cf_compliant/data/'
OUT='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS/WRF-CanESM2/mean_changes/output_local/'
geo='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS/WRF-CanESM2/geo_em.d03.nc'
stem='CanESM2-WRF_t2_d03'

lm=xr.open_dataset(geo)['LANDMASK'].squeeze().values
land=(lm==1); ocean=(lm==0)
def circ(f,h,c=365):
    d=f-h; return ((d+c/2)%c)-c/2

h=xr.open_dataset(BASE+f'historical/{stem}_historical_daily_1986_2005_300x300.nc')['T2']
w=np.cos(np.deg2rad(h.lat.values)).astype(np.float64)
ch=h.groupby('time.dayofyear').mean(dim='time')
hmax,chmin=ch.argmax(dim='dayofyear')+1, ch.argmin(dim='dayofyear')+1
hamp=ch.max(dim='dayofyear')-ch.min(dim='dayofyear')
del h
rows=[]
for scen in ['rcp45','rcp85']:
    f=xr.open_dataset(BASE+f'{scen}/{stem}_{scen}_daily_2046_2065_300x300.nc')['T2']
    cf=f.groupby('time.dayofyear').mean(dim='time')
    fmax,fmin=cf.argmax(dim='dayofyear')+1, cf.argmin(dim='dayofyear')+1
    famp=cf.max(dim='dayofyear')-cf.min(dim='dayofyear')
    del f,cf
    dmax, dmin = circ(fmax,hmax), circ(fmin,chmin)
    for region,m in [('Land',land),('Ocean',ocean)]:
        den=np.nansum(np.where(m,w,0.0))
        num=lambda a: np.nansum(np.where(m,a.values,np.nan)*w)/den
        rows.append(dict(scenario=scen,region=region,
            hist_day_max=num(hmax), fut_day_max=num(fmax), delta_day_max=num(dmax),
            hist_day_min=num(chmin), fut_day_min=num(fmin), delta_day_min=num(dmin)))
df=pd.DataFrame(rows)
df.to_csv(OUT+'doy_T2_land_ocean.csv',index=False)
print(df.to_string(index=False))
