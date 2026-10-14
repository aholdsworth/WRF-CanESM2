# WRF-CanESM2

Analysis codes and notebooks for the manuscript
*"Projecting future changes over the coastal Pacific Northwest: climate and
extremes from a convection-permitting model."*

The analysis is based on WRF d03 (300 x 300 curvilinear grid over the
north-western North America coast, lon -133.8 to -116.6, lat 43.4 to 54.7)
driven by CanESM2, with a 1986-2005 historical period and 2046-2065 future
periods (RCP4.5 and RCP8.5). Results are compared against CanRCM4, raw
CanESM2, and NA-CORDEX/CMIP5 ensemble fields.

## Folders

| folder | contents |
|--------|----------|
| `Evaluation/` | All station/model-comparison skill code (1986–2005, 5 models): `climate_taylor/` (3-panel climate Taylor), `bias_figure/` (D03 bias figure), `bias_comparisons/` (multi-model bias maps), `elevation_bias/` (elevation-bias grid + shared skill-calc core: `seasonal_pooling.py`, station/skill/R² table producers), `buoy_validation/` (D03 buoy validation), `climate_land/` + `climate_buoy/` (land & buoy climatology tables), `extremes_validation/` (p75/p95 T & P validation). See its README. |
| `DataPrep/` | CF conversion (`cf_convert.py`) + station/buoy data extraction (`make_wrf_canon.py`, `extract_canrcm4_stations.py`, per-buoy `.sh`, climatology). See its README. |
| `mean_changes/` | Seasonal/annual mean-change scripts (deltas + t-tests for WRF d03, CanRCM4, CanESM2) and the `plot_means_*` notebooks (incl. CORDEX/CMIP5 comparison panels), annual averages, and elevation-binned change. |
| `extremes/` | Block-bootstrap scripts for extreme-percentile significance (per-gridcell and pooled land/ocean regional null test) and the `plot_extremes*` notebooks. (Station-level p75/p95 validation lives in `Evaluation/extremes_validation/`. |
| `CompoundEvents/` | Two-stage compound-extremes pipeline: per-DOY rolling percentile thresholds, then counts of warm/cold x wet/dry compound event days; drizzle-day counts; plotting notebooks. |
| `compression/` | Compression/PPC verification notebooks and scripts for the WRF d03 run. |
| `land_and_ocean/` | Regional-average (land/ocean) analysis notebooks. |
| `WRFDomainLib.py` | Shared helpers for the WRF d03 domain (coordinates, masks, plotting). |

## Where to read the results

Result tables and methods are committed in each subfolder's `results/`
directory (not in this README):

- Report-grade Table 1 (mean annual bias, RMSE, pooled R², 5 models ×
  T/P/wind) and Table 2 (buoys): built from the canonical pickles, see
  `Evaluation/elevation_bias/results/` and `Evaluation/buoy_validation/results/`.
- Elevation-bias diagnostic grid: `Evaluation/elevation_bias/results/`
  (`elev_bias_results_for_report.md`, `skill_summary.csv`,
  `skill_table_bias_mae_rmse.tex`) + methods in `elev_bias_methods_math.md`.
- Taylor diagrams, bias figures: `Evaluation/climate_taylor/` and
  `Evaluation/bias_figure/` (PNGs are regenerated locally, not committed).

These numbers supersede any earlier hand-derived values (the verification
pipeline had unit/gate bugs, all fixed 2026-09-25…29; see the per-folder
READMEs).

## How to run

All `Evaluation/` scripts run from the **analysis working tree** (this
repo is a *subfolder* of it) using the `py_2024` environment, with

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # OpenBLAS thread-spin gotcha

Dependency order (each step caches to `pickles/`):

1. `DataPrep/` — CF conversion (`cf_convert.py`) and per-station/per-buoy
   extractions (`make_wrf_canon.py`, `extract_canrcm4_stations.py`, buoy
   `.sh` scripts) → `data/wrf_canon/`, `data/cdo_extractions/`.
2. `make_station_tables.py` → `pickles/station_tables.pkl` — the **single
   canonical land table** (173 stations × 5 models); everything downstream
   takes station sets from it.
3. `make_skill_tables.py` (seconds; projection of the table),
   `make_bias_figure.py` + `bias_figure.ipynb`, `taylor_data.py` +
   `taylor_figure.py`, buoy `buoy_{tas,wind}_bias.py`,
   `extremes_validation/` — all read the table or the validated buoy pickles.

## Data

Data files are **not** in this repository. The canonical working-tree
locations (`data/`, `pickles/`, `files/`, `plotting_files/`, `report/`),
formats, sizes and checksums are documented in [`DATA.md`](DATA.md);
observations (NOAA GHCN-d, ECCC, BC BCH, buoys) and model output (WRF
d01–d03 driven by CanESM2, CanRCM4, raw CanESM2, NA-CORDEX/CMIP5) are
archived out-of-band.

## Notes

- Scripts in `DataPrep/`, `mean_changes/`, `extremes/`, and `CompoundEvents/`
  contain **absolute paths** to the original data on the HPC system where the
  analysis was run; edit the path constants at the top of each script before
  re-running.
- Each folder's README documents its expected inputs/outputs and the
  notebook-to-script relationships.
