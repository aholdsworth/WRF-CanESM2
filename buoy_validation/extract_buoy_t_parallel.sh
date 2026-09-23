#!/usr/bin/env bash
# Parallel D02 (remaining 5) + D03 (all 14) T2 buoy extraction.
# Usage: extract_buoy_t_parallel.sh [NJOBS]
set -uo pipefail
cd "$(dirname "$0")"
CDO="${CDO:-/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo}"
NJOBS="${1:-8}"
SRC2="/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/historical/variables_complete/t_d02_hourly.nc"
SRC3="/gpfs/fs7/dfo/hpcmc/pfm/evg000/historical_verification/variables_complete/t_d03.nc"
LOGDIR="station_cdo/logs"; OUTDIR="data/wrf_stations"
mkdir -p "$LOGDIR" "$OUTDIR"

remap_one() {
  local dom="$1" id="$2" src="$3"
  local w="weight_files_gen/CanESM2-WRF-${dom^^}/weights_st_${id}.nc"
  local out="$OUTDIR/t_${dom}_st${id}.nc"
  local log="$LOGDIR/buoyt_par_${dom}_${id}.log"
  if [ -s "$out" ]; then echo "SKIP existing $(basename "$out")" | tee "$log"; return 0; fi
  local desc=""
  for a in NOAA ECCC; do
    [ -f "station_cdo/buoys/descs/${a}_buoy_${id}.txt" ] && desc="station_cdo/buoys/descs/${a}_buoy_${id}.txt" && break
  done
  if [ -z "$desc" ] || [ ! -f "$w" ] || [ ! -f "$src" ]; then
    echo "MISSING dep for $dom/$id (desc=$desc w=$w src=$src)" | tee "$log"; return 2
  fi
  {
    echo "CMD: $CDO -b F64 -f nc remap,$desc,$w $src $out"
    "$CDO" -b F64 -f nc "remap,$desc,$w" "$src" "$out" 2>&1
    echo "exit=$?"
  } > "$log" 2>&1
  echo "$id ($dom): $(tail -1 "$log")"
}
export -f remap_one
export CDO SRC2 SRC3 LOGDIR OUTDIR

# Build the work list: dom id
{
  for id in 46132 46134 46146 46204 46206; do echo "d02 $id"; done
  for id in 46005 46029 46041 46050 46087 46088 46089 46131 46132 46134 46146 46204 46206 46207; do echo "d03 $id"; done
} | xargs -P "$NJOBS" -n 2 bash -c 'remap_one "$0" "$1" "$( [ "$0" = d02 ] && echo "$SRC2" || echo "$SRC3" )" '
echo "=== parallel done ==="
