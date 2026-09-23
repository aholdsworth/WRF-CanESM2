# Buoy near-surface air-temperature (`tas`) validation

Buoy validation of near-surface air temperature for the CanESM2-WRF
nested-downscaling study, mirroring the wind validation
(`../buoy_wind_validation/`).

## Deliverables
- **Report text + tables:** `results/buoy_tas_results_for_report.md`
  (description paragraph, pooled & per-buoy tables; raw bias in K +
  bias % + MAE/RMSE).
- **Mathematical methods (LaTeX):** `results/buoy_tas_methods_math.md`.
- **Data:** `results/buoy_tas_bias.csv`, `results/buoy_tas_bias_pooled.csv`.

## Files
| file | purpose |
|---|---|
| `buoy_tas_bias.py` | Main analysis. Parses the buoy `.mat` "Temperature:Air" channel (deg C → K, ~15-min → hourly → 90% cascade), loads the five per-buoy model series (WRF D01/D02/D03 hourly T2, CanESM2/CanRCM4 daily tas), writes per-buoy and pooled skill metrics. Self-contained (inlines the wind pipeline's cascade; no external project imports). |
| `extract_buoy_t_d01.sh` | Re-extract D01 T2 for all 14 buoys from PFM `t_d01_hourly.nc` (the original PFM pull had rotated ECCC filenames). Resumable (`--force`). Run from the NOTEBOOKS root. |
| `extract_buoy_t_parallel.sh` | Re-extract D02 + D03 T2 for the 14 buoys in parallel (`xargs -P 8`) from PFM `t_d02_hourly.nc` / `t_d03.nc`, with coordinate-verified descriptors + SCRIP weights (the rotated Sep-15 files were re-extracted). Resumable (skips non-empty outputs). Run from the NOTEBOOKS root. |
| `verify_buoy_tas_coords.py` | Check embedded coords vs `station_cdo/buoys/buoy_descs.csv` for all buoy tas files (14 D01 + 14 D02 + 14 D03 + 16 CanESM2 + 16 CanRCM4). Exit 0 = all match. **Run after any re-extraction.** |

## Model columns (2026-09-23 state)
| column | source | cadence | status |
|---|---|---|---|
| WRF D01 (75 km) | `data/wrf_stations/t_d01_{ECCC,NOAA}_buoy_<ID>.nc` (T2, K) | hourly | ✅ re-extracted 2026-09-23, coords verified |
| WRF D02 (15 km) | `data/wrf_stations/t_d02_st<ID>.nc` (T2, K) | hourly | ✅ re-extracted 2026-09-23, coords verified |
| WRF D03 (3 km) | `data/wrf_stations/t_d03_st<ID>.nc` (T2, K) | hourly | ✅ re-extracted 2026-09-23, coords verified (46005 excluded: off-domain nearest cell) |
| CanESM2 (~100 km) | `data/cdo_extractions/canesm2_raw_buoys/tas_buoy_<ID>.nc` (K, 1850–2005 → 1986–2005) | daily | ✅ coords verified |
| CanRCM4 (~25 km) | `data/cdo_extractions/canrcm4_buoys/tas_buoy_<ID>.nc` (K, 1986–2005) | daily | ✅ coords verified |

## Observation
`MB_<ID>_HM.mat` (NOAA/MSC-OS met buoys): column 4 = "Temperature:Air",
deg C, ~15-min cadence, NaN = missing; column 0 = MATLAB datenum
(classic epoch 719529); header coords = truth. Sensor is in the buoy
met-package above the sea surface (package wind stresses referenced to
4.2 m); no height scaling is applied (log-profile correction is a momentum
correction). Restricted to 1986-01-01..2005-12-31.

## Conventions (locked, same as wind)
- 90 % missing-data cascade (day/month/season), year = mean of passing
  seasons, kept if ≥ `--year-gate` (default 3) of 4 seasons; `--tag
  _strict` = 4/4.
- Apples-to-apples: model annuals computed on the obs surviving years.
- Metrics on shared-year annuals: period means, raw bias (K), bias %,
  MAE, RMSE (n ≥ 2).
- Pooled headline = mean over the core set: buoys with ≥ `--min-years`
  (default 8) shared years. For temperature this is **9 buoys**: 46005,
  46041, 46050, 46131, 46132, 46146, 46204, 46206, 46207. The D03 pooled
  row covers 8 (46005 excluded — off-domain nearest cell).
- Resolution labels: D01 outer nest 75 km, D02 middle 15 km, D03 inner 3 km,
  CanESM2 parent GCM ~100 km (coarsest), CanRCM4 ~25 km.
- PFM ECCC-buoy filenames are NOT trusted — WRF files are matched by
  embedded coords (the script does this automatically).

## Usage
```
# from the NOTEBOOKS root (py_2024 python; obs root overridable):
python buoy_tas_bias.py --year-gate 3
python verify_buoy_tas_coords.py          # after any re-extraction
./extract_buoy_t_d01.sh --force           # redo D01 (needs PFM mount)
./extract_buoy_t_parallel.sh 8            # redo D02 + D03 (needs PFM mount)
```
Environment: py_2024 only. No network. Read-only obs + PFM.
