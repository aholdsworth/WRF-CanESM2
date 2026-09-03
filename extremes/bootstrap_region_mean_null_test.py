import os
import sys
import datetime
import pickle
from collections import defaultdict

import numpy as np
import xarray as xr
from netCDF4 import Dataset
from joblib import Parallel, delayed


# ------------------ Load data ------------------
def load_data(path, varname):
    nc = Dataset(path, 'r')
    data = nc.variables[varname][:]
    times = nc.variables['time'][:]
    base_date = datetime.datetime(1986, 1, 1) if 'historical' in path else datetime.datetime(2046, 1, 1)
    datetimes = [base_date + datetime.timedelta(hours=float(h)) for h in times]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    return data, datetimes, lat, lon


# ------------------ Shift seasons, keep yearly blocks ------------------
def shift_seasons_to_blocks(data, times):
    seasons_out = defaultdict(list)

    years = np.array([t.year for t in times])
    months = np.array([t.month for t in times])

    season_year = np.copy(years)
    season_year[months == 12] += 1

    unique_years = np.unique(season_year)

    for yr in unique_years:
        djf_idx = np.where((season_year == yr) & ((months == 12) | (months <= 2)))[0]
        mam_idx = np.where((years == yr) & (months >= 3) & (months <= 5))[0]
        jja_idx = np.where((years == yr) & (months >= 6) & (months <= 8))[0]
        son_idx = np.where((years == yr) & (months >= 9) & (months <= 11))[0]

        if len(djf_idx) > 0:
            seasons_out['DJF'].append(data[djf_idx])
        if len(mam_idx) > 0:
            seasons_out['MAM'].append(data[mam_idx])
        if len(jja_idx) > 0:
            seasons_out['JJA'].append(data[jja_idx])
        if len(son_idx) > 0:
            seasons_out['SON'].append(data[son_idx])

    return seasons_out


# ------------------ Precip top-quartile filtering ------------------
def filter_top_quartile(hist_data, target_data, percentile=75):
    thresh = np.nanpercentile(hist_data, percentile)
    filtered = np.where(target_data < thresh, np.nan, target_data)
    return filtered, thresh


# ------------------ Land/Ocean mask ------------------
def load_land_ocean_mask(geo_em_file):
    ds = xr.open_dataset(geo_em_file)
    landmask = ds['LANDMASK'].squeeze()

    rename_dict = {}
    if 'south_north' in landmask.dims:
        rename_dict['south_north'] = 'lat'
    if 'west_east' in landmask.dims:
        rename_dict['west_east'] = 'lon'
    if rename_dict:
        landmask = landmask.rename(rename_dict)

    landmask = landmask.values
    return (landmask == 1), (landmask == 0)


# ------------------ Core percentile-change calculator ------------------
def compute_change_field(hist_concat, fut_concat, perc, relative='no', minusmed='no', variable='t'):
    p0 = 75 if variable == 'pr' else 50

    perc_hist = np.nanpercentile(hist_concat, perc, axis=0)
    perc_fut = np.nanpercentile(fut_concat, perc, axis=0)

    if minusmed == 'yes':
        med_hist = np.nanpercentile(hist_concat, p0, axis=0)
        med_fut = np.nanpercentile(fut_concat, p0, axis=0)

        if perc > p0:
            perc_hist = perc_hist - med_hist
            perc_fut = perc_fut - med_fut
        else:
            perc_hist = med_hist - perc_hist
            perc_fut = med_fut - perc_fut

    if relative == 'yes':
        with np.errstate(divide='ignore', invalid='ignore'):
            delta = 100.0 * (perc_fut - perc_hist) / perc_hist
    else:
        delta = perc_fut - perc_hist

    return delta


# ------------------ One pooled-null replicate ------------------
def single_null_iteration(_, pooled_blocks, n_hist, n_fut, perc, relative, minusmed,
                          variable, land_flat, ocean_flat):
    rng = np.random.default_rng()

    n_total = len(pooled_blocks)

    # Pooled bootstrap under null: resample seasonal-years from the pooled set
    hist_idx = rng.choice(n_total, size=n_hist, replace=True)
    fut_idx = rng.choice(n_total, size=n_fut, replace=True)

    hist_concat = np.concatenate([pooled_blocks[i] for i in hist_idx], axis=0)
    fut_concat = np.concatenate([pooled_blocks[i] for i in fut_idx], axis=0)

    delta = compute_change_field(
        hist_concat,
        fut_concat,
        perc=perc,
        relative=relative,
        minusmed=minusmed,
        variable=variable
    )

    delta_flat = delta.reshape(-1)

    land_mean = np.nanmean(delta_flat[land_flat]) if land_flat is not None else np.nan
    ocean_mean = np.nanmean(delta_flat[ocean_flat]) if ocean_flat is not None else np.nan

    return land_mean, ocean_mean


# ------------------ Seasonal pooled-null wrapper ------------------
def pooled_null_bootstrap_region(
    season,
    hist_blocks,
    fut_blocks,
    land_mask,
    ocean_mask,
    variable,
    perc,
    relative='no',
    minusmed='no',
    n_boot=10,
    n_procs=8
):
    n_hist = len(hist_blocks)
    n_fut = len(fut_blocks)

    print(f'{season}: n_hist_years = {n_hist}, n_fut_years = {n_fut}')

    pooled_blocks = hist_blocks + fut_blocks

    land_flat = land_mask.reshape(-1) if land_mask is not None else None
    ocean_flat = ocean_mask.reshape(-1) if ocean_mask is not None else None

    # Observed statistic from the real split
    hist_concat_obs = np.concatenate(hist_blocks, axis=0)
    fut_concat_obs = np.concatenate(fut_blocks, axis=0)

    delta_obs = compute_change_field(
        hist_concat_obs,
        fut_concat_obs,
        perc=perc,
        relative=relative,
        minusmed=minusmed,
        variable=variable
    )
    delta_obs_flat = delta_obs.reshape(-1)

    land_obs = np.nanmean(delta_obs_flat[land_flat]) if land_flat is not None else np.nan
    ocean_obs = np.nanmean(delta_obs_flat[ocean_flat]) if ocean_flat is not None else np.nan

    # Null distribution
    results = Parallel(n_jobs=n_procs)(
        delayed(single_null_iteration)(
            None, pooled_blocks, n_hist, n_fut, perc, relative, minusmed,
            variable, land_flat, ocean_flat
        )
        for _ in range(n_boot)
    )

    land_null = np.array([r[0] for r in results], dtype=np.float64)
    ocean_null = np.array([r[1] for r in results], dtype=np.float64)

    # Two-sided p-values relative to null centered at no-change
    land_p = np.nanmean(np.abs(land_null) >= np.abs(land_obs))
    ocean_p = np.nanmean(np.abs(ocean_null) >= np.abs(ocean_obs))

    out = {
        'season': season,
        'land_obs': land_obs,
        'ocean_obs': ocean_obs,
        'land_null_mean': np.nanmean(land_null),
        'ocean_null_mean': np.nanmean(ocean_null),
        'land_null_sd': np.nanstd(land_null, ddof=1),
        'ocean_null_sd': np.nanstd(ocean_null, ddof=1),
        'land_null_ci': np.nanpercentile(land_null, [2.5, 97.5]),
        'ocean_null_ci': np.nanpercentile(ocean_null, [2.5, 97.5]),
        'land_p_value': land_p,
        'ocean_p_value': ocean_p,
        'land_null_samples': land_null,
        'ocean_null_samples': ocean_null,
    }

    return out



# ------------------ Main script ------------------
if __name__ == "__main__":
    print("The script has the name %s" % sys.argv[0])

    variable = sys.argv[1]
    relative = sys.argv[2]
    minusmed = sys.argv[3]
    print('varibale ', variable)
    period = 'rcp45'
    domain = 'd03'
    basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
    geo_em_file = f'{basepath}/domain/geo_em.{domain}.nc'

    var_map = {'tmin': 'T2', 't': 'T2', 'pr': 'pr', 'wind': 'wspd'}
    file_map = {
        'tmin': f't_{domain}_daily',
        't': f't_{domain}_daily',
        'pr': f'pr_{domain}_daily',
        'wind': f'wind_{domain}_daily_wspd'
    }

    perc = 5 if variable == 'tmin' else 95
    var = var_map[variable]
    filename = file_map[variable]

    hist_path = os.path.join(basepath, f"historical/variables_complete/{filename}.nc")
    fut_path = os.path.join(basepath, f"{period}/variables_complete/{filename}.nc")

    print('loading files')
    data_hist, time_hist, lat, lon = load_data(hist_path, var)
    data_fut, time_fut, lat, lon = load_data(fut_path, var)

    if variable == 'pr':
        data_hist, thresh = filter_top_quartile(data_hist, data_hist, percentile=75)
        data_fut = np.where(data_fut < thresh, np.nan, data_fut)

    print('building seasonal-year blocks')
    seasons_hist = shift_seasons_to_blocks(data_hist, time_hist)
    seasons_fut = shift_seasons_to_blocks(data_fut, time_fut)

    print('loading land/ocean masks')
    land_mask, ocean_mask = load_land_ocean_mask(geo_em_file)

    save_path = '/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/Extremes/'
    os.makedirs(save_path, exist_ok=True)

    tag = '_rel' if relative == 'yes' else ''
    tag2 = f'_med{75 if variable == "pr" else 50}' if minusmed == 'yes' else ''

#    seasons_to_run = ['DJF', 'MAM', 'JJA', 'SON']
    seasons_to_run = ['JJA', 'SON']
#    seasons_to_run = ['DJF', 'MAM']

    for season in seasons_to_run:
        out_file = os.path.join(
            save_path,
            f'{var}_{season}_{period}_region_null_boots_1000_perc_{perc}{tag}{tag2}.pkl'
        )

        # Skip if already done
        if os.path.exists(out_file):
            print(f"{season} already exists, skipping.")
            continue

        print(f'Processing {season}...')

        res = pooled_null_bootstrap_region(
            season=season,
            hist_blocks=seasons_hist[season],
            fut_blocks=seasons_fut[season],
            land_mask=land_mask,
            ocean_mask=ocean_mask,
            variable=variable,
            perc=perc,
            relative=relative,
            minusmed=minusmed,
            n_boot=1000,
            n_procs=8
        )

        with open(out_file, 'wb') as f:
            pickle.dump(res, f)

        print(f'Saved {season} to {out_file}')
