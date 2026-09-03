#!/home/amh001/space_fs7/software_2022/python/py_2024/bin/python
"""
Convert a WRF d03 daily file to CF-compliant, zlib-compressed form.

Usage:
  cf_convert.py --var pr|wspd|t2 --scenario historical|rcp45|rcp85 [--force]

Source (read-only):
  historical : WRF_FILES/original/{pr,t,wind}_d03_daily*.nc
  rcp45/rcp85: WRF_FILES/original/RCP{45,85}/{pr,t,wind}_d03_daily*.nc
Output:
  WRF_FILES/cf_compliant/data/{scenario}/CanESM2-WRF_{var}_d03_{scenario}_daily_{start}_{end}_300x300.nc

Guarantees:
  - metadata only: data values are bit-identical to the source (verified)
  - output is written to a .tmp file, verified, then atomically renamed
  - refuses to overwrite an existing output unless --force is given
  - date range in the output name is read from the file's time axis

The three source files have no units attribute for the data variable;
units below are declared by this script (mm/day, m s-1, K).
"""
import argparse
import os
import sys
from datetime import datetime, timezone

import netCDF4
import numpy as np

ROOT = "/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/WRF_FILES"
ORIG = os.path.join(ROOT, "original")
DATA = os.path.join(ROOT, "cf_compliant", "data")

# var: (data variable, standard_name, long_name, units)
# Source files are discovered by content (which variable they hold), not by
# name — the originals use inconsistent naming across scenarios.
VARS = {
    "pr":   ("pr",   "precipitation_amount", "daily total precipitation",     "mm/day"),
    "t2":   ("T2",   "air_temperature",      "daily mean 2-m air temperature", "K"),
    "wspd": ("wspd", "air_speed",            "daily mean wind speed",         "m s-1"),
}

def find_source(scenario, data_var):
    """Locate the source file in the scenario dir that contains data_var."""
    src_dir = ORIG if scenario == "historical" else os.path.join(ORIG, "RCP" + scenario[3:].upper())
    matches = []
    for fn in sorted(os.listdir(src_dir)):
        if not fn.endswith(".nc"):
            continue
        try:
            with netCDF4.Dataset(os.path.join(src_dir, fn), "r") as ds:
                if data_var in ds.variables:
                    matches.append(fn)
        except Exception:
            continue
    if len(matches) != 1:
        sys.exit(f"expected exactly one source with variable {data_var!r} in {src_dir}, found: {matches}")
    return os.path.join(src_dir, matches[0])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var", required=True, choices=sorted(VARS))
    ap.add_argument("--scenario", required=True, choices=("historical", "rcp45", "rcp85"))
    ap.add_argument("--force", action="store_true", help="overwrite existing output")
    args = ap.parse_args()

    data_var, std_name, long_name, units = VARS[args.var]
    scenario = args.scenario
    src = find_source(scenario, data_var)
    if not os.path.exists(src):
        sys.exit(f"source not found: {src}")

    # --- read source (header + full data) ---
    with netCDF4.Dataset(src, "r") as ds_in:
        nt, ny, nx = (ds_in.dimensions[d].size for d in ("time", "lat", "lon"))
        assert (nt, ny, nx) == (7305, 300, 300), f"unexpected dims {(nt, ny, nx)}"
        assert "time_bnds" in ds_in.variables, "expected time_bnds in source"
        assert data_var in ds_in.variables, f"expected {data_var} in source"

        time_data = ds_in["time"][:]
        lon_data  = ds_in["lon"][:]
        lat_data  = ds_in["lat"][:]
        data      = ds_in[data_var][:]
        time_units   = ds_in["time"].getncattr("units")
        time_calendar = ds_in["time"].getncattr("calendar")
        cell_methods = ds_in[data_var].getncattr("cell_methods") if "cell_methods" in ds_in[data_var].ncattrs() else None
        src_history  = ds_in.getncattr("history") if "history" in ds_in.ncattrs() else ""

    t0 = netCDF4.num2date(time_data[0], time_units, time_calendar)
    t1 = netCDF4.num2date(time_data[-1], time_units, time_calendar)
    y0, y1 = t0.year, t1.year

    out_dir = os.path.join(DATA, scenario)
    os.makedirs(out_dir, exist_ok=True)
    dst = os.path.join(out_dir, f"CanESM2-WRF_{args.var}_d03_{scenario}_daily_{y0}_{y1}_300x300.nc")
    if os.path.exists(dst) and not args.force:
        sys.exit(f"output exists: {dst} (use --force to overwrite)")

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %z")
    history = (f"{src_history}\n"
               f"{now}: CF-compliant copy of {os.path.basename(src)}: "
               f"renamed dims lon->y, lat->x; removed time_bnds, CDO, CDI, _CoordinateAxisType; "
               f"added standard_name, long_name, units ({units}) to {data_var}"
               + (f"; kept cell_methods ({cell_methods})" if cell_methods else "")
               + "; Conventions=CF-1.11; rebuilt with zlib compression (complevel=4).")

    print(f"src : {src}")
    print(f"dst : {dst}")
    print(f"time: {t0} -> {t1}  ({nt} days)  units: {time_units!r}")
    print(f"var : {data_var}  range: {data.min():.3f} .. {data.max():.3f}")

    # --- write ---
    tmp = dst + ".tmp"
    comp = dict(zlib=True, complevel=4)
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds_out:
        ds_out.createDimension("time", nt)
        ds_out.createDimension("y", ny)
        ds_out.createDimension("x", nx)
        ds_out.setncattr("Conventions", "CF-1.11")
        ds_out.setncattr("history", history)

        v = ds_out.createVariable("time", "f8", ("time",), **comp)
        v.standard_name = "time"; v.long_name = "time"
        v.units = time_units; v.calendar = time_calendar; v.axis = "T"
        v[:] = time_data

        v = ds_out.createVariable("lat", "f8", ("y", "x"), **comp)
        v.standard_name = "latitude"; v.long_name = "latitude"; v.units = "degrees_north"
        v[:] = lat_data

        v = ds_out.createVariable("lon", "f8", ("y", "x"), **comp)
        v.standard_name = "longitude"; v.long_name = "longitude"; v.units = "degrees_east"
        v[:] = lon_data

        v = ds_out.createVariable(data_var, "f8", ("time", "y", "x"), **comp)
        v.standard_name = std_name
        v.long_name = long_name
        v.units = units
        if cell_methods:
            v.cell_methods = cell_methods
        v.coordinates = "lat lon"
        v[:] = data

    # --- verify bit-identity, then replace ---
    with netCDF4.Dataset(src, "r") as a, netCDF4.Dataset(tmp, "r") as b:
        for name in ("time", "lon", "lat", data_var):
            assert np.array_equal(a[name][:], b[name][:]), f"MISMATCH in {name}!"
        print("verified: all data bit-identical to source")

    os.replace(tmp, dst)
    print(f"done: {dst}  ({os.path.getsize(dst)/1e9:.3f} GB)")

if __name__ == "__main__":
    main()
