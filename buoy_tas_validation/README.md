# Buoy near-surface air-temperature (`tas`) validation

In progress — D01 re-extraction (this commit) is the first step.

## Status
- **D01 (75 km) re-extracted 2026-09-23** from PFM `t_d01_hourly.nc`
  (99x99, 175320 h, 1986-2005) with the verified D01 1-point nearest-cell
  weights (`weight_files_gen/CanESM2-WRF-D01/weights_st_<ID>.nc`, rev2).
  The original PFM `t_d01_{ECCC,NOAA}_buoy_*.nc` pull had **rotated ECCC
  filenames** (6 of 7 files held a different buoy's location; 46131 had no
  file; 46207 was duplicated). All 6 misnamed ECCC files are preserved as
  `t_d01_ECCC_buoy_<ID>.nc.rotated_20260923` in `data/wrf_stations/`.
  New files: 14 buoys (7 ECCC + 46005/46029/46041/46050/46087/46088/46089),
  all verified: embedded coords = truth (d = 0.0000 deg), 175320 steps,
  0 NaN, var `T2` (K).
- **CanESM2 (~100 km):** `data/cdo_extractions/canesm2_raw_buoys/tas_buoy_<ID>.nc`,
  daily, 1850-2005 (slice to 1986-2005 = days 49640:56940), 16 buoys, coords verified.
- **CanRCM4 (~25 km):** `data/cdo_extractions/canrcm4_buoys/tas_buoy_<ID>.nc`,
  monthly, 1986-2005, 16 buoys, coords verified.
- **D02 (15 km) / D03 (3 km):** no buoy `t` extractions yet — pending.
- **Obs:** buoy `.mat` files carry a SAT column (see
  `buoy_wind_bias_handoff_20260922.md`: "carry only wind + SAT; SAT/SST
  intentionally skipped" in the wind run). Needs the BuoyData mount to
  pin down column index/units.

## Files
| file | purpose |
|---|---|
| `extract_buoy_t_d01.sh` | Re-extract D01 T2 for all 14 buoys from PFM source. Resumable (`--force` to redo). Run from the NOTEBOOKS root (paths relative to it, like the wind scripts). |
| `verify_buoy_tas_coords.py` | Check embedded coords vs `station_cdo/buoys/buoy_descs.csv` for all 48 buoy tas files (16 D01 + 16 CanESM2 + 16 CanRCM4). Exit 0 = all match. **Run after any re-extraction.** |

## Coordinate convention
Truth coords = `station_cdo/buoys/buoy_descs.csv` (from the `.mat` headers /
canonical station inventory). A file is "correctly named" iff its embedded
1-point (lon, lat) matches its ID's truth coords within 0.05 deg. Never trust
PFM ECCC-buoy filenames — match by embedded coordinates.
