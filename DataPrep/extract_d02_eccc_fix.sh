#!/usr/bin/env bash
# Re-extract the 6 rotated ECCC D02 buoys with correct coords.
set -uo pipefail
cd "$(dirname "$0")"
CDO="${CDO:-/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo}"
SRC="/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/historical/variables_complete/t_d02_hourly.nc"
LOGDIR="station_cdo/logs"; OUTDIR="data/wrf_stations"
mkdir -p "$LOGDIR" "$OUTDIR"
for id in 46131 46132 46134 46146 46204 46206; do
  w="weight_files_gen/CanESM2-WRF-D02/weights_st_${id}.nc"
  out="$OUTDIR/t_d02_st${id}.nc"
  log="$LOGDIR/buoyt_d02_fix_${id}.log"
  if [ -s "$out" ]; then echo "SKIP existing $out" | tee "$log"; continue; fi
  # pick a desc with correct coords (ECCC or NOAA; both correct for these per buoy_descs.csv)
  desc=""
  for a in ECCC NOAA; do
    [ -f "station_cdo/buoys/descs/${a}_buoy_${id}.txt" ] && desc="station_cdo/buoys/descs/${a}_buoy_${id}.txt" && break
  done
  echo "CMD: $CDO -b F64 -f nc remap,$desc,$w $SRC $out" | tee "$log"
  "$CDO" -b F64 -f nc "remap,$desc,$w" "$SRC" "$out" 2>&1 >> "$log"
  echo "exit=$? for $id" | tee -a "$log"
done
echo "=== d02 eccc fix done ==="
