import numpy as np
from netCDF4 import Dataset
import datetime
import scipy.stats
import os

# -------- CONFIGURATION --------
variable = 't'  # 't', 'pr', or 'wind'
period = 'rcp45'
domain = 'd03'
out_file = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/MEANS/seasonal_deltas_and_pvals_{variable}_{period}_{domain}.npz'

basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
#if period=='hist':
#    data_path = basepath+'historical/variables_complete/'
#else:
#    data_path = basepath+period+'/variables_complete/'

if variable == "t":
    var = 'T2'
    filename = 't_' + domain + '_daily'
elif variable == "pr":
    var = 'pr'
    filename = 'pr_' + domain + '_daily'
elif variable == "wind":
    var = 'wspd'
    filename = "wind_" + domain + '_daily_wspd'
else:
    raise ValueError("Unsupported variable.")

# -------- HELPER FUNCTIONS --------

def load_data(path, varname):
    nc = Dataset(path, 'r')
    data = nc.variables[varname][:]
    times = nc.variables['time'][:]
    base_date = datetime.datetime(1986, 1, 1) if 'hist' in path else datetime.datetime(2046, 1, 1)
    datetimes = [base_date + datetime.timedelta(hours=h) for h in times]
    
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    return data, datetimes, lat,lon

def group_by_season(data, time, season):
    season_months = {
        "DJF": [12, 1, 2],
        "MAM": [3, 4, 5],
        "JJA": [6, 7, 8],
        "SON": [9, 10, 11]
    }
    season_vals = {}
    for i, date in enumerate(time):
        if date.month in season_months[season]:
            year = date.year if date.month != 12 else date.year + 1
            if year not in season_vals:
                season_vals[year] = []
            season_vals[year].append(data[i])
    result = []
    for year, vals in sorted(season_vals.items()):
        stacked = np.stack(vals)
        result.append((year, stacked.mean(axis=0)))
    return result

def compute_seasonal_means(data, time):
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    results = {}
    for season in seasons:
        results[season] = group_by_season(data, time, season)
    return results

def compute_annual_means(data, time):
    annual_vals = {}
    for i, date in enumerate(time):
        year = date.year
        if year not in annual_vals:
            annual_vals[year] = []
        annual_vals[year].append(data[i])
    result = []
    for year, vals in sorted(annual_vals.items()):
        stacked = np.stack(vals)
        result.append((year, stacked.mean(axis=0)))
    return result

def separate_years(mean_list):
    return np.stack([x[1] for x in mean_list])

def seasonal_change_and_ttest(hist_means, fut_means, var_type='t'):
    hist_vals = separate_years(hist_means)
    fut_vals = separate_years(fut_means)
    if var_type == 't':
        delta = fut_vals.mean(axis=0) - hist_vals.mean(axis=0)
    else:
        delta = ((fut_vals.mean(axis=0) - hist_vals.mean(axis=0)) / hist_vals.mean(axis=0)) * 100
    ttest = scipy.stats.ttest_ind(hist_vals, fut_vals, axis=0)
    return delta, ttest.pvalue

# -------- LOAD DATA --------
hist_path = os.path.join(basepath, f"historical/variables_complete/{filename}.nc")
fut_path = os.path.join(basepath, f"{period}/variables_complete/{filename}.nc")

data_hist, time_hist,lat,lon = load_data(hist_path, var)
data_fut, time_fut,lat,lon = load_data(fut_path, var)

# -------- COMPUTE MEANS --------
seasonal_hist = compute_seasonal_means(data_hist, time_hist)
seasonal_fut = compute_seasonal_means(data_fut, time_fut)

annual_hist = compute_annual_means(data_hist, time_hist)
annual_fut = compute_annual_means(data_fut, time_fut)

# -------- COMPUTE DELTAS + T-TESTS --------
seasons = ['DJF', 'MAM', 'JJA', 'SON']
seasonal_deltas = {}
seasonal_pvals = {}

for s in seasons:
    delta, pval = seasonal_change_and_ttest(seasonal_hist[s], seasonal_fut[s], var_type=variable)
    seasonal_deltas[s] = delta
    seasonal_pvals[s] = pval

delta_ANN, pval_ANN = seasonal_change_and_ttest(annual_hist, annual_fut, var_type=variable)

# -------- SAVE RESULTS --------
np.savez_compressed(out_file,
    **{f"{s}_delta": seasonal_deltas[s] for s in seasons},
    **{f"lat": lat},
    **{f"lon": lon},
    **{f"{s}_pval": seasonal_pvals[s] for s in seasons},
    ANN_delta=delta_ANN, ANN_pval=pval_ANN
)

print(f"Saved results to {out_file}")

