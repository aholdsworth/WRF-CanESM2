#!/usr/bin/env bash
# Fresh WRF d02 + d03 T2 (near-surface air temperature) extraction for all
# 14 buoys, from the PFM sources:
#   d02: t_d02_hourly.nc (200x200 rotated curvilinear, 175320 h)
#        $PFM/historical/variables_complete/
#   d03: t_d03.nc (300x300, 175320 h)
#        $EVG/historical_verification/variables_complete/
# with the locally-generated 1-point nearest-cell weights
# (weight_files_gen/CanESM2-WRF-D0{2,3}/weights_st_<ID>.nc) built for the
# wind work against these exact source grids (verified identical).
#
# Method consistent with the wind extraction (extract_buoy_wind_d02d03.sh):
# CDO remap, 1-point desc + hand-written SCRIP weights, -b F64 -f nc.
# d03 weights use bilinear inside the domain and nearest-cell at the edge
# (46005); its d03 value is excluded from the analysis anyway (off-domain
# nearest-cell artifact, documented in the wind results).
#
# Output: data/wrf_stations/t_d0{2,3}_st<ID>.nc   (var T2, K)
# Resumable: skips existing non-empty outputs (use --force to redo).
set -uo pipefail
cd "$(dirname "$0")"

CDO="${CDO:-/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo}"
PFM="${PFM:-/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF}"
EVG="${EVG:-/gpfs/fs7/dfo/hpcmc/pfm/evg000}"
SRC2="$PFM/historical/variables_complete/t_d02_hourly.nc"
SRC3="$EVG/historical_verification/variables_complete/t_d03.nc"
DESCDIR="station_cdo/buoys/descs"
LOGDIR="station_cdo/logs"
OUTDIR="data/wrf_stations"
mkdir -p "$LOGDIR" "$OUTDIR"

FORCE=0
[ "${1:-}" = "--force" ] && FORCE=1

IDS="46005 46029 46041 46050 46087 46088 46089 46131 46132 46134 46146 46204 46206 46207"

remap_one() {
    local dom="$1" id="$2" src="$3"
    local w="weight_files_gen/CanESM2-WRF-${dom^^}/weights_st_${id}.nc"
    local out="$OUTDIR/t_${dom}_st${id}.nc"
    local log="$LOGDIR/buoyt_${dom}_${id}.log"
    if [ "$FORCE" -eq 0 ] && [ -s "$out" ]; then
        echo "SKIP existing $(basename "$out")" | tee "$log"
        echo "exit=0" >> "$log"
        return 0
    fi
    # pick a desc with the right coords (ECCC or NOAA; coords identical either way)
    local desc=""
    for a in ECCC NOAA; do
        [ -f "$DESCDIR/${a}_buoy_${id}.txt" ] && desc="$DESCDIR/${a}_buoy_${id}.txt" && break
    done
    if [ -z "$desc" ]; then echo "NO DESC for $id" | tee "$log"; echo "exit=2" >> "$log"; return 2; fi
    if [ ! -f "$w" ]; then echo "MISSING WEIGHTS: $w" | tee "$log"; echo "exit=2" >> "$log"; return 2; fi
    if [ ! -f "$src" ]; then echo "MISSING SOURCE: $src" | tee "$log"; echo "exit=2" >> "$log"; return 2; fi
    {
        echo "CMD: $CDO -b F64 -f nc remap,$desc,$w $src $out"
        "$CDO" -b F64 -f nc "remap,$desc,$w" "$src" "$out" 2>&1
        echo "exit=$?"
    } > "$log" 2>&1
    echo "$id ($dom): $(tail -1 "$log")"
}

for id in $IDS; do remap_one d02 "$id" "$SRC2"; done
for id in $IDS; do remap_one d03 "$id" "$SRC3"; done

echo "=== done ==="
