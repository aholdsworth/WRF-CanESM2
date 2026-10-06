#!/usr/bin/env bash
# Extract WRF d02 + d03 hourly wspd for the 14 buoy locations.
#
# d03: local 1-point CDO descs (station_cdo/buoys/descs, .mat truth coords) +
#      hand-written SCRIP weights in weight_files_gen/CanESM2-WRF-D03/
#      weights_st<ID>.nc (gen_buoy_wind_d03_weights.py; CDO 1.9.9 genbic
#      rejects geo_em).  Bilinear inside the domain; nearest-cell for points
#      at/beyond the domain edge (46005).
# d02: same local descs + hand-written SCRIP weights in
#      weight_files_gen/CanESM2-WRF-D02/weights_st<ID>.nc, built against the
#      200x200 curvilinear src grid of wind_d02_hourly.nc.  Eva's PFM d02
#      weights are NOT used: their ECCC dst points are the same rotated
#      assignment documented in station_cdo/buoys/BUOY_STATUS.md.
# Sources: variables_complete/{wind_d02_hourly,wspd_d03_hourly}.nc (175320 h).
# Output: data/wrf_stations/wind_d{02,03}_st<ID>.nc (f8, 1x1, mirrors
#         t_d{02,03}_st<ID>.nc naming).
# Resumable: skips existing non-empty outputs.  Plain loop + wait (no xargs).
#
# Env: NP (default 8)
set -uo pipefail
cd "$(dirname "$0")"

NP="${NP:-8}"
CDO=/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo
PFM=/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF
SRC_D02="$PFM/historical/variables_complete/wind_d02_hourly.nc"
SRC_D03="$PFM/historical/variables_complete/wspd_d03_hourly.nc"
W_D02="weight_files_gen/CanESM2-WRF-D02"
W_D03="weight_files_gen/CanESM2-WRF-D03"
LOGDIR="station_cdo/logs"
DESCDIR="station_cdo/buoys/descs"
mkdir -p "$LOGDIR" data/wrf_stations

# --- 1-point descs for every buoy in the table (regenerate, idempotent) ---
while IFS=, read -r agency id lon lat; do
    case "$agency" in \#*|agency) continue ;; esac
    [ -z "$agency" ] && continue
    [ "$id" = "desw1" ] && continue
    [ "$id" = "wpow1" ] && continue
    out="${DESCDIR}/${agency}_buoy_${id}.txt"
    {
        echo "# CDO grid file"
        echo "# Station ${id} (buoy) — ${agency}, regenerated from buoy_descs.csv"
        echo "#"
        echo "gridtype = curvilinear"
        echo "gridsize = 1"
        echo "xsize = 1"
        echo "ysize = 1"
        echo "#"
        echo "# Longitudes"
        echo "xvals = ${lon}"
        echo "#"
        echo "# Latitude"
        echo "yvals = ${lat}"
    } > "$out"
done < "station_cdo/buoys/buoy_descs.csv"

IDS="$(for f in "$DESCDIR"/*_buoy_46*.txt; do basename "$f" .txt | sed 's/.*_buoy_//'; done)"

remap_one() {
    local d="$1" id="$2"
    local desc="" w src out log
    for f in "$DESCDIR"/*_buoy_"$id".txt; do [ -f "$f" ] && desc="$f" && break; done
    if [ -z "${desc:-}" ]; then echo "NO DESC $id" >&2; return 2; fi
    case "$d" in
        d02) w="$W_D02/weights_st_${id}.nc"; src="$SRC_D02" ;;
        d03) w="$W_D03/weights_st_${id}.nc"; src="$SRC_D03" ;;
    esac
    out="data/wrf_stations/wind_${d}_st${id}.nc"
    log="$LOGDIR/buoywind_${d}_${id}.log"
    if [ -s "$out" ]; then echo "SKIP existing $(basename "$out")" > "$log"; echo "exit=0" >> "$log"; return 0; fi
    if [ ! -f "$w" ]; then echo "MISSING WEIGHTS: $w" > "$log"; echo "exit=2" >> "$log"; return 1; fi
    {
        echo "CMD: $CDO -b F64 -f nc remap,$desc,$w $src $out"
        "$CDO" -b F64 -f nc "remap,$desc,$w" "$src" "$out" 2>&1
        echo "exit=$?"
    } > "$log" 2>&1
}

run_jobs() {
    # Launch all 28 jobs as background processes, letting the OS scheduler
    # + GPFS handle concurrency.  CDO remap is I/O-bound on a 56 GB file;
    # 14-way concurrency is fine on a 128-core node.
    for id in $IDS; do
        for d in d02 d03; do
            remap_one "$d" "$id" &
        done
    done
    wait
}
run_jobs

echo
echo "=== SUMMARY ==="
ok=0; skip=0; fail=0
for id in $IDS; do
    for d in d02 d03; do
        log="$LOGDIR/buoywind_${d}_${id}.log"
        if [ -f "$log" ]; then
            if grep -q '^SKIP' "$log"; then skip=$((skip+1))
            elif grep -q 'exit=0' "$log"; then ok=$((ok+1))
            else fail=$((fail+1)); echo "FAILED: $d $id"; fi
        else
            fail=$((fail+1)); echo "NO LOG: $d $id"
        fi
    done
done
echo "buoys=$(echo "$IDS" | grep -c .) x 2: done=${ok} skipped=${skip} failed=${fail}"
