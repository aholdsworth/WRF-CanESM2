import numpy as np, xarray as xr, pandas as pd, os, time as t
BASE='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/WRF_FILES/cf_compliant/data/'
OUT='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS/WRF-CanESM2/mean_changes/output_local/'
geo='/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS/WRF-CanESM2/geo_em.d03.nc'
fname={'t':'CanESM2-WRF_t2_d03','pr':'CanESM2-WRF_pr_d03','wind':'CanESM2-WRF_wspd_d03'}
varname={'t':'T2','pr':'pr','wind':'wspd'}

lm = xr.open_dataset(geo)['LANDMASK'].squeeze()
if 'south_north' in lm.dims: lm = lm.rename({'south_north':'y','west_east':'x'})
lm = lm.values
land = (lm==1); ocean=(lm==0)

def circ(f,h,c=365):
    d=f-h; return ((d+c/2)%c)-c/2

rows=[]
for var,stem in fname.items():
    h=xr.open_dataset(BASE+f'historical/{stem}_historical_daily_1986_2005_300x300.nc')[varname[var]]
    lat=h.lat
    w=np.cos(np.deg2rad(lat.values)).astype(np.float64)
    ch = h.groupby('time.dayofyear').mean(dim='time')
    hmax,chmin = ch.argmax(dim='dayofyear')+1, ch.argmin(dim='dayofyear')+1
    hamp = ch.max(dim='dayofyear')-ch.min(dim='dayofyear')
    del h
    for scen in ['rcp45','rcp85']:
        f=xr.open_dataset(BASE+f'{scen}/{stem}_{scen}_daily_2046_2065_300x300.nc')[varname[var]]
        cf = f.groupby('time.dayofyear').mean(dim='time')
        fmax,fmin = cf.argmax(dim='dayofyear')+1, cf.argmin(dim='dayofyear')+1
        famp = cf.max(dim='dayofyear')-cf.min(dim='dayofyear')
        del f,cf
        dmax=circ(fmax,hmax); dmin=circ(fmin,chmin)
        relamp=100*(famp-hamp)/hamp
        for region,m in [('Land',land),('Ocean',ocean)]:
            num=lambda a: np.nansum(np.where(m,a,np.nan)*w)
            den=np.nansum(np.where(m,w,0.0))
            rows.append(dict(var=var,scenario=scen,region=region,
                hist_max=float(num(hmax.values)/den), fut_max=float(num(fmax.values)/den),
                dmax=float(num(dmax.values)/den),
                hist_min=float(num(chmin.values)/den), fut_min=float(num(fmin.values)/den),
                dmin=float(num(dmin.values)/den),
                rel_amp_change_pct=float(num(relamp.values)/den)))
df=pd.DataFrame(rows)
df.to_csv(OUT+'doy_chen_results.csv',index=False)
print(df.to_string(index=False))
