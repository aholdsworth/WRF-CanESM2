#!/usr/bin/env bash
# WRF d01 wspd extraction for the 5 USW WA/OR land stations (Olympia, Portland,
# SeaTac, Yakima, Astoria) from PFM wind_d01_hourly.nc + hand-written SCRIP
# bilinear 1-point weights (weight_files_gen/CanESM2-WRF-D01/weights_st_<ID>.nc,
# generated 2026-09-29 from files/NOAA_d03_stations_wind.csv truth coords).
# Same method family as extract_buoy_wind_d01.sh (2026-09-23 rev2).
# Output: data/wrf_stations/wind_d01_NOAA_<ID>.nc
set -uo pipefail
cd "$(dirname "$0")"
CDO=/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo
SRC=/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/historical/variables_complete/wind_d01_hourly.nc
W=weight_files_gen/CanESM2-WRF-D01
DESCDIR=station_cdo/land/descs
LOGDIR=station_cdo/logs
OUTDIR=data/wrf_stations
mkdir -p "$LOGDIR" "$OUTDIR"
FORCE=0
[ "${1:-}" = "--force" ] && FORCE=1
IDS="USW00024227 USW00024229 USW00024233 USW00024243 USW00094224"
for id in $IDS; do
    out="$OUTDIR/wind_d01_NOAA_${id}.nc"
    log="$LOGDIR/wind_d01_usw_${id}.log"
    if [ "$FORCE" -eq 0 ] && [ -s "$out" ]; then
        echo "SKIP $out"
        continue
    fi
    desc="$DESCDIR/NOAA_${id}.txt"
    wgt="$W/weights_st_${id}.nc"
    echo "CMD: $CDO -b F64 -f nc remap,$desc,$wgt $SRC $out" > "$log"
    "$CDO" -b F64 -f nc "remap,$desc,$wgt" "$SRC" "$out" >> "$log" 2>&1
    echo "exit=$?" >> "$log"
    echo "$id: $(tail -1 "$log")"
done
echo "=== done ==="
