import os
import datetime
import xarray as xr
import sys
# --------------------
# User settings
# --------------------
domain = "d03"

basepath = '/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
percentiles_path = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/percentiles/'
out_path = f'/gpfs/fs7/dfo/hpcmc/pfm/amh001/DATA/WRF/CompoundCounts/'

os.makedirs(out_path, exist_ok=True)

file_map = {'tas': f't_{domain}_daily', 'pr': f'pr_{domain}_daily'}

# --------------------
# Helper functions
# --------------------
def decode_time(ds, period):
    if "time" in ds:
        try:
            ds = xr.decode_cf(ds)
        except Exception:
            base_date = datetime.datetime(1986, 1, 1) if period == "hist" else datetime.datetime(2046, 1, 1)
            ds["time"] = [base_date + datetime.timedelta(hours=int(h)) for h in ds["time"].values]
    return ds

def get_season_months(season):
    return {
        "mam": [3,4,5],
        "jja": [6,7,8],
        "son": [9,10,11],
        "djf": [12,1,2]
    }[season]

def load_thresholds(var, base):
    """Load percentile thresholds for tas or pr, based on base period."""
    if var == "tas":
       perc= 90 if t_type=="warm" else 10
       fname = f"{var}{perc}p_{base}_roll5.nc" 
        #fname = f"WRF_tas90p_{base}.nc" if t_type == "warm" else f"WRF_tas10p_{base}.nc"
    elif var == "pr":
        perc = 75 if pr_type =="wet" else 25
        fname = f"{var}{perc}p_{base}_roll29.nc" 
    variablename = f"{var}_{perc}p"
    print('the filenames for perc ', fname)
    return xr.open_dataset(percentiles_path + fname)[variablename]

def count_events(tas, pr, t_thresh, pr_thresh, months, season):

    tas = tas.sel(time=tas["time.month"].isin(months))
    pr  = pr.sel(time=pr["time.month"].isin(months))

    dayofyear = tas["time.dayofyear"]

    t_daily_thresh  = t_thresh.isel(dayofyear=dayofyear - 1)
    pr_daily_thresh = pr_thresh.isel(dayofyear=dayofyear - 1)

    if t_type == "warm" and pr_type == "wet":
        events = (tas > t_daily_thresh) & (pr > pr_daily_thresh)
    elif t_type == "cold" and pr_type == "wet":
        events = (tas < t_daily_thresh) & (pr > pr_daily_thresh)
    elif t_type == "warm" and pr_type == "dry":
        events = (tas > t_daily_thresh) & (pr <= pr_daily_thresh)
    elif t_type == "cold" and pr_type == "dry":
        events = (tas < t_daily_thresh) & (pr <= pr_daily_thresh)

    if season == "djf":
        season_year = xr.where(
            events["time.month"] == 12,
            events["time.year"] + 1,
            events["time.year"]
        )
    else:
        season_year = events["time.year"]

    return (
        events
        .assign_coords(season_year=("time", season_year.values))
        .groupby("season_year")
        .sum("time")
    )
# --------------------
# Main
# --------------------
def process_season(season):
    months = get_season_months(season)

    # Paths
    tas_path = os.path.join(basepath, f"{'historical' if period=='hist' else period}/variables_complete/{file_map['tas']}.nc")
    pr_path  = os.path.join(basepath, f"{'historical' if period=='hist' else period}/variables_complete/{file_map['pr']}.nc")

    # Load data
    tas = xr.open_dataset(tas_path)["T2"]
    pr  = xr.open_dataset(pr_path)["pr"]
    #print(pr)

    pr = decode_time(pr, period)
    tas = decode_time(tas, period)
    # Always use historical thresholds
    
    t_thresh  = load_thresholds("tas", "hist")
    pr_thresh = load_thresholds("pr",  "hist")

    counts = count_events(tas, pr, t_thresh, pr_thresh, months,season)

    out_file = f"CanESM2-WRF-{t_type}_{pr_type}_{period}_{season}_yearly_lte.nc"
    counts.name = "count"

    counts.to_netcdf(os.path.join(out_path, out_file))
    return f"{season.upper()} done."

if __name__ == "__main__":
    print ("The script has the name %s" % (sys.argv[0]))
    #t_type = "warm"       # "warm" or "cold"
    #pr_type = "wet"       # "wet" or "dry"
    #period = "hist"      # "hist", "rcp45", "rcp85"
    period = (sys.argv[1])
    t_type = (sys.argv[2]) # )
    pr_type=(sys.argv[3])
    print(period, t_type, pr_type)
#variable = 't' # 't', 'pr', 'wind'
    for season in ["mam", "jja", "son", "djf"]:
        print(process_season(season))
    print("All seasons processed (counts only).")

