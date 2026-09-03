import os
import datetime
from collections import defaultdict
import numpy as np
from netCDF4 import Dataset
import pickle
import sys
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


# ------------------ Seasonal-year blocks ------------------
def shift_seasons_to_blocks(data, times):
    seasons_out = defaultdict(list)

    years = np.array([t.year for t in times])
    months = np.array([t.month for t in times])

    season_year = years.copy()
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



# ------------------ Single bootstrap iteration (BLOCK) ------------------
def single_iteration_blocks(_, pooled_blocks, n_hist, n_fut,
                            perc_delta, perc, relative, minusmed, variable):

    rng = np.random.default_rng()
    n_total = len(pooled_blocks)

    hist_idx = rng.choice(n_total, size=n_hist, replace=True)
    fut_idx  = rng.choice(n_total, size=n_fut, replace=True)

    hist_concat = np.concatenate([pooled_blocks[i] for i in hist_idx], axis=0)
    fut_concat  = np.concatenate([pooled_blocks[i] for i in fut_idx], axis=0)

    p0 = 75 if variable == 'pr' else 50


    if minusmed == 'yes':
        perc_hist = np.nanpercentile(hist_concat, perc, axis=0)
        med_hist  = np.nanpercentile(hist_concat, p0,   axis=0)
        perc_fut  = np.nanpercentile(fut_concat,  perc, axis=0)
        med_fut   = np.nanpercentile(fut_concat,  p0,   axis=0)
        if perc > p0:
            perc_hist = perc_hist - med_hist
            perc_fut  = perc_fut - med_fut
        else:
            perc_hist = med_hist - perc_hist
            perc_fut  = med_fut - perc_fut
    else:
        perc_hist = np.nanpercentile(hist_concat, perc, axis=0)
        perc_fut  = np.nanpercentile(fut_concat, perc, axis=0)
    if relative == "yes":
        with np.errstate(divide='ignore', invalid='ignore'):
            diff = 100. * (perc_fut - perc_hist) / perc_hist
    else:
        diff = perc_fut - perc_hist

    return diff >= perc_delta


# ------------------ Bootstrap function ------------------
def bootstrappin(season, hist_blocks, fut_blocks, varname,
                 iters=1000, perc=95, relative='no', minusmed='no', n_procs=32):

    print(f'{season}: n_hist_years = {len(hist_blocks)}, n_fut_years = {len(fut_blocks)}')

    n_hist = len(hist_blocks)
    n_fut  = len(fut_blocks)

    pooled_blocks = hist_blocks + fut_blocks

    p0 = 75 if varname == 'pr' else 50

    tag = '_rel' if relative == 'yes' else ''
    tag2 = f'_med{p0}' if minusmed == 'yes' else ''

    output_file = os.path.join(
        save_path,
        f'{varname}_{season}_{period}_boots_years_{iters}_perc_{perc}{tag}{tag2}.pkl'
    )

    # Observed statistic
    hist_concat = np.concatenate(hist_blocks, axis=0)
    fut_concat  = np.concatenate(fut_blocks, axis=0)
    if minusmed == 'yes':
        perc_hist = np.nanpercentile(hist_concat, perc, axis=0)
        med_hist  = np.nanpercentile(hist_concat, p0,   axis=0)
        perc_fut  = np.nanpercentile(fut_concat,  perc, axis=0)
        med_fut   = np.nanpercentile(fut_concat,  p0,   axis=0)

        if perc > p0:
            perc_hist = perc_hist - med_hist
            perc_fut  = perc_fut - med_fut
        else:
            perc_hist = med_hist - perc_hist
            perc_fut  = med_fut - perc_fut
    else:
        perc_hist = np.nanpercentile(hist_concat, perc, axis=0)
        perc_fut  = np.nanpercentile(fut_concat,  perc, axis=0)

    if relative == "yes":
        with np.errstate(divide='ignore', invalid='ignore'):
            perc_delta = 100. * (perc_fut - perc_hist) / perc_hist
    else:
        perc_delta = perc_fut - perc_hist

    # Bootstrap
    results = Parallel(n_jobs=n_procs)(
        delayed(single_iteration_blocks)(
            None, pooled_blocks, n_hist, n_fut,
            perc_delta, perc, relative, minusmed, varname
        )
        for _ in range(iters)
    )

    count = np.sum(results, axis=0)
    p_value = count / iters

    with open(output_file, 'wb') as f:
        pickle.dump({'perc_delta': perc_delta, 'p_value': p_value}, f)

    return p_value


# ------------------ Run all seasons ------------------
def run_bootstrap_all_seasons(varname, seasons_hist, seasons_fut,
                             save_path, n_boot=1000, perc=95,
                             relative='no', minusmed='no', n_procs=32):

    results = {}

#    for season in ['SON']:
    for season in ['JJA']:
    #for season in ['DJF', 'MAM', 'JJA', 'SON']:
        print(f"Bootstrapping {season}...")
        p_val = bootstrappin(
            season,
            seasons_hist[season],
            seasons_fut[season],
            varname,
            iters=n_boot,
            perc=perc,
            relative=relative,
            minusmed=minusmed,
            n_procs=n_procs
        )
        results[season] = p_val

    return results


# ------------------ Main ------------------
print("The script has the name %s" % (sys.argv[0]))

variable = sys.argv[1]
relative = sys.argv[2]
minusmed = sys.argv[3]

print('The script is running for ', variable, ' with rel', relative, ' and minusmed', minusmed)
period = 'rcp85'
domain = 'd03'
print('The scenario is ', period, ' and domain is ', domain)
basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'

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

print('loading file')
data_hist, time_hist, lat, lon = load_data(hist_path, var)
data_fut, time_fut, lat, lon = load_data(fut_path, var)


print('building seasonal-year blocks')
seasons_hist = shift_seasons_to_blocks(data_hist, time_hist)
seasons_fut = shift_seasons_to_blocks(data_fut, time_fut)

save_path = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/Extremes/'

print('running bootstrap')
results = run_bootstrap_all_seasons(
    var,
    seasons_hist,
    seasons_fut,
    save_path,
    n_boot=1000,
    perc=perc,
    relative=relative,
    minusmed=minusmed,
    n_procs=8
)
