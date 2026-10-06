#!/usr/bin/env bash
# Fresh WRF d01 T2 (near-surface air temperature) extraction for all 14 buoys,
# from the PFM source t_d01_hourly.nc (99x99, 175320 h, 1986-2005) + the
# locally-generated 1-point nearest-cell weights
# (weight_files_gen/CanESM2-WRF-D01/weights_st_<ID>.nc, rev2, verified).
#
# Why: the original t_d01_{ECCC,NOAA}_buoy_*.nc files pulled from PFM had
# ROTATED ECCC filenames (6 of 7 held a different buoy's location; 46131 had
# no file at all; 46207 duplicated).  Re-extracting all 14 from source
# guarantees each file's name matches its embedded coords, consistent with the
# wind re-extraction (extract_buoy_wind_d01.sh).
#
# Method: CDO remap with a 1-point desc + hand-written SCRIP weights,
# -b F64 -f nc.  Variable stays T2 (K).
#
# Output: data/wrf_stations/t_d01_{ECCC,NOAA}_buoy_<ID>.nc
# Resumable: skips existing non-empty outputs (use --force to redo).
set -uo pipefail
cd "$(dirname "$0")"

CDO="${CDO:-/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo}"
PFM="${PFM:-/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF}"
SRC="${T_D01_HOURLY:-$PFM/historical/variables_complete/t_d01_hourly.nc}"
W="weight_files_gen/CanESM2-WRF-D01"
DESCDIR="station_cdo/buoys/descs"
LOGDIR="station_cdo/logs"
OUTDIR="data/wrf_stations"
mkdir -p "$LOGDIR" "$OUTDIR"

FORCE=0
[ "${1:-}" = "--force" ] && FORCE=1

# 7 ECCC + 46005/46029/46041 (legacy NOAA, in packaged set) +
# 46050/46087/46088/46089 (new NOAA, model files for completeness)
remap_one() {
    local agency="$1" id="$2"
    local desc="$DESCDIR/${agency}_buoy_${id}.txt"
    local w="$W/weights_st_${id}.nc"
    local out="$OUTDIR/t_d01_${agency}_buoy_${id}.nc"
    local log="$LOGDIR/buoyt_d01_${id}.log"
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
    echo "$id ($agency): $(tail -1 "$log")"
}

for id in 46131 46132 46134 46146 46204 46206 46207; do
    remap_one "ECCC" "$id"
done
for id in 46005 46029 46041 46050 46087 46088 46089; do
    remap_one "NOAA" "$id"
done

echo "=== done ==="
