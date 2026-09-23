#!/usr/bin/env bash
# Fresh WRF d01 wspd extraction for the 7 ECCC buoy locations, from the PFM
# source wind_d01_hourly.nc + locally-generated 1-point nearest-cell weights
# (weight_files_gen/CanESM2-WRF-D01/weights_st_<ID>.nc).
#
# Why: the original wind_d01_ECCC_buoy_*.nc files had ROTATED filenames (6 of
# 7 held a different buoy's location), and a rename attempt left the set in an
# inconsistent state.  Re-extracting all 7 from source guarantees each file's
# name matches its embedded coords.  46131 (49.91,-124.99) had NO original
# file at all, so this also fills that gap.
#
# Method is consistent with d02/d03 (extract_buoy_wind_d02d03.sh): CDO remap
# with a 1-point desc + hand-written SCRIP weights, -b F64 -f nc.
#
# Output: data/wrf_stations/wind_d01_ECCC_buoy_<ID>.nc
# Resumable: skips existing non-empty outputs (use --force to redo).
set -uo pipefail
cd "$(dirname "$0")"

CDO="${CDO:-/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo}"
PFM="${PFM:-/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF}"
SRC="${WIND_D01_HOURLY:-$PFM/historical/variables_complete/wind_d01_hourly.nc}"
W="weight_files_gen/CanESM2-WRF-D01"
DESCDIR="station_cdo/buoys/descs"
LOGDIR="station_cdo/logs"
OUTDIR="data/wrf_stations"
mkdir -p "$LOGDIR" "$OUTDIR"

FORCE=0
[ "${1:-}" = "--force" ] && FORCE=1

IDS="46131 46132 46134 46146 46204 46206 46207"

remap_one() {
    local id="$1"
    local desc="$DESCDIR/ECCC_buoy_${id}.txt"
    local w="$W/weights_st_${id}.nc"
    local out="$OUTDIR/wind_d01_ECCC_buoy_${id}.nc"
    local log="$LOGDIR/buoywind_d01_${id}.log"
    if [ "$FORCE" -eq 0 ] && [ -s "$out" ]; then
        echo "SKIP existing $(basename "$out")" | tee "$log"
        echo "exit=0" >> "$log"
        return 0
    fi
    if [ ! -f "$desc" ]; then echo "NO DESC $desc" | tee "$log"; echo "exit=2" >> "$log"; return 2; fi
    if [ ! -f "$w" ]; then echo "MISSING WEIGHTS: $w" | tee "$log"; echo "exit=2" >> "$log"; return 2; fi
    {
        echo "CMD: $CDO -b F64 -f nc remap,$desc,$w $SRC $out"
        "$CDO" -b F64 -f nc "remap,$desc,$w" "$SRC" "$out" 2>&1
        echo "exit=$?"
    } > "$log" 2>&1
    echo "$id: $(tail -1 "$log")"
}

for id in $IDS; do
    remap_one "$id"
done

echo "=== done ==="
