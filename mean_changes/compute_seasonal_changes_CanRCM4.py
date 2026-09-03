import numpy as np
from netCDF4 import Dataset
import datetime
import scipy.stats
import os
from netCDF4 import Dataset, num2date

# -------- CONFIGURATION --------
variable = 't'  # 't', 'pr', or 'wind'
period = 'rcp85'
domain = 'CanRCM4'
out_file = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/MEANS/seasonal_deltas_and_pvals_{variable}_{period}_{domain}.npz'

basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/gridded_model_daily/CanRCM4/'
if variable=='wind':
   var='sfcWind'
   filename=variable+'_NAM22_'
elif variable=='t':
   var = 'tas'  
   filename=var+'_NAM22_'
elif variable=='pr':
   filename=variable+'_NAM22_'
   var = variable  
   filename=variable+'_NAM22_'
# -------- HELPER FUNCTIONS --------
def load_data(path, varname):
    nc = Dataset(path, 'r')
    data = nc.variables[varname][:]
    times = nc.variables['time'][:]
    time_units = nc.variables['time'].units
    try:
        time_calendar = nc.variables['time'].calendar
    except AttributeError:
        time_calendar = "standard"
    
    # decode times properly
    datetimes = num2date(times, units=time_units, calendar=time_calendar)
    
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    nc.close()
    return data, datetimes, lat, lon

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
#    print('mean list', mean_list)
    return np.stack([x[1] for x in mean_list])

def seasonal_change_and_ttest(hist_means, fut_means, var_type='t'):
    #print('shapes hist_measn', np.shape(hist_means))
    hist_vals = separate_years(hist_means)
    fut_vals = separate_years(fut_means)
    if var_type == 't':
        delta = fut_vals.mean(axis=0) - hist_vals.mean(axis=0)
    else:
        delta = ((fut_vals.mean(axis=0) - hist_vals.mean(axis=0)) / hist_vals.mean(axis=0)) * 100
 #   print('max and mins', np.amax(hist_vals), np.amin(hist_vals))
  #  print('max and minsfut', np.amax(fut_vals), np.amin(fut_vals))
   # print(np.nanstd(hist_vals), np.nanstd(fut_vals))
    #print("hist_vals shape:", hist_vals.shape)
    #print("fut_vals shape:", fut_vals.shape)
    #print("axis for t-test: 0")
    print("n years hist:", hist_vals.shape[0], "n years fut:", fut_vals.shape[0])

    ttest = scipy.stats.ttest_ind(hist_vals, fut_vals, axis=0)
    return delta, ttest.pvalue

# -------- LOAD DATA --------
hist_path = os.path.join(basepath, f"{filename}hist.nc")
fut_path = os.path.join(basepath, f"{filename}{period}.nc")

data_hist, time_hist,lat,lon = load_data(hist_path, var)
data_fut, time_fut,lat,lon = load_data(fut_path, var)

#print('shape of loaded data', np.shape(data_hist), np.shape(time_hist))
# -------- COMPUTE MEANS --------
seasonal_hist = compute_seasonal_means(data_hist, time_hist)
seasonal_fut = compute_seasonal_means(data_fut, time_fut)
#print(seasonal_hist)
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

print('pval', pval)
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

