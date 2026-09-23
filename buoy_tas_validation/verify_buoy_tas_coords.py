#!/usr/bin/env python
"""Verify buoy model-file naming: embedded (lon, lat) vs truth coords.

Reads truth coords from station_cdo/buoys/buoy_descs.csv and checks every
per-buoy model file's embedded 1-point coordinate:

  * data/wrf_stations/t_d01_{ECCC,NOAA}_buoy_<ID>.nc      (var T2)
  * data/wrf_stations/t_d02_st<ID>.nc                     (var T2)
  * data/wrf_stations/t_d03_st<ID>.nc                     (var T2)
  * data/cdo_extractions/canesm2_raw_buoys/tas_buoy_<ID>.nc (var tas)
  * data/cdo_extractions/canrcm4_buoys/tas_buoy_<ID>.nc   (var tas)

Usage:
    python verify_buoy_tas_coords.py [--root NOTEBOOKS] [--tol 0.05]

Exit code 0 = all files match their name; 1 = at least one mismatch.
Run this after any re-extraction so a rotated/misnamed file can never
silently enter the SAT bias analysis.
"""
import argparse, math, os, re, sys
import netCDF4 as nc

PATTERNS = [
    (r"t_d01_(ECCC|NOAA)_buoy_(\w+)\.nc$", "data/wrf_stations", 0.05),
    (r"t_d02_st(\w+)\.nc$",                "data/wrf_stations", 0.05),
    (r"t_d03_st(\w+)\.nc$",                "data/wrf_stations", 0.20),
    (r"tas_buoy_(\w+)\.nc$", "data/cdo_extractions/canesm2_raw_buoys", 0.05),
    (r"tas_buoy_(\w+)\.nc$", "data/cdo_extractions/canrcm4_buoys", 0.05),
]

def dmatch(l1, la1, l2, la2):
    return math.hypot(l1 - l2, (la1 - la2) * math.cos(math.radians((la1 + la2) / 2)))

def load_truth(path):
    truth = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split(",")
        if p[1] == "id":
            continue
        truth[p[1]] = (float(p[2]), float(p[3]))  # (lon, lat)
    return truth

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.getcwd())
    ap.add_argument("--tol", type=float, default=0.05,
                    help="max allowed offset for d01/d02/GCM files, deg "
                         "(d03 uses a fixed 0.20 deg tolerance for the 3 km cell)")
    args = ap.parse_args()
    root = args.root
    truth = load_truth(os.path.join(root, "station_cdo/buoys/buoy_descs.csv"))
    bad = 0
    checked = 0
    for pat, d, tol in PATTERNS:
        d = os.path.join(root, d)
        if not os.path.isdir(d):
            print("  (dir missing, skip) %s" % d)
            continue
        for f in sorted(os.listdir(d)):
            m = re.match(pat, f)
            if not m:
                continue
            bid = m.group(len(m.groups()))
            if bid not in truth:
                print("  UNLISTED  %s/%s (id not in buoy_descs.csv)" % (d, f))
                bad += 1
                continue
            ds = nc.Dataset(os.path.join(d, f))
            lon = float(ds.variables["lon"][0, 0])
            lat = float(ds.variables["lat"][0, 0])
            ds.close()
            tlon, tlat = truth[bid]
            ddeg = dmatch(lon, lat, tlon, tlat)
            checked += 1
            status = "OK" if ddeg <= tol else "MISMATCH"
            if ddeg > tol:
                bad += 1
            print("  %-72s (%9.4f,%10.4f) truth=(%9.4f,%10.4f) d=%.4f  %s"
                  % (os.path.join(os.path.basename(d), f), lat, lon, tlat, tlon, ddeg, status))
    print("%d files checked, %d mismatch(es)" % (checked, bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
