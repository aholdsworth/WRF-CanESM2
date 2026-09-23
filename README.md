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
| `DataPrep/` | `cf_convert.py` — converts the raw WRF d03 daily files (pr, T2, wspd) to CF-compliant, zlib-compressed netCDF4 used for publication (bit-identical data, metadata only). See its README. |
| `mean_changes/` | Seasonal/annual mean-change scripts (deltas + t-tests for WRF d03, CanRCM4, CanESM2) and the `plot_means_*` notebooks (incl. CORDEX/CMIP5 comparison panels), annual averages, and elevation-binned change. |
| `extremes/` | Block-bootstrap scripts for extreme-percentile significance (per-gridcell and pooled land/ocean regional null test) and the `plot_extremes*` notebooks. |
| `elevation_bias/` | Elevation-bias diagnostic grid: the `elev_bias_grid.ipynb` walkthrough (5 rows × 4 cols: map, T-bias, PR-bias, wind-bias scatter vs elevation bias; PR also in absolute mm/day) plus the table producers (`make_station_tables.py`, `make_skill_tables.py`, `make_elev_bias_grid.py`, `plot_elev_bias_revised.py`, `seasonal_pooling.py`) and the result/methods tables in `results/`. |
| `CompoundEvents/` | Two-stage compound-extremes pipeline: per-DOY rolling percentile thresholds, then counts of warm/cold x wet/dry compound event days; drizzle-day counts; plotting notebooks. |
| `compression/` | Compression/PPC verification notebooks and scripts for the WRF d03 run. |
| `land_and_ocean/` | Regional-average (land/ocean) analysis notebooks. |
| `WRFDomainLib.py` | Shared helpers for the WRF d03 domain (coordinates, masks, plotting). |

## Elevation-bias diagnostic grid — headline tables

From `elevation_bias/` (1986–2005 mean, 175 stations; sign = model − obs).
Full tables, methods and the absolute-PR (mm/day) re-expression are in
[`elevation_bias/results/elev_bias_results_for_report.md`](elevation_bias/results/elev_bias_results_for_report.md).

**Elevation bias (m)** — model orography − observed station elevation
(extraction-consistent / CDO bicubic-weight orography):

| Row | n | median (m) | std (m) | mean (m) |
|-----|---|-----------:|--------:|---------:|
| D3 (WRF d03) | 164 | 41.6 | 165.5 | 89.8 |
| D15 (WRF d02) | 164 | 170.3 | 300.6 | 247.2 |
| D75 (WRF d01) | 164 | 304.6 | 318.1 | 331.7 |
| CanRCM4 | 164 | 242.3 | 348.8 | 294.3 |
| CanESM2 | 164 | 391.2 | 368.5 | 375.4 |

**Pooled seasonal skill** (bias / MAE / RMSE, season-pooled per-station;
15 rows = 5 datasets × {t, pr, wind}; see `elevation_bias/results/skill_summary.csv`):

| dataset | T bias (°C) | T RMSE (°C) | PR bias (%) | PR RMSE (%) | wind bias (%) |
|---------|------------:|------------:|------------:|------------:|--------------:|
| D03 (3 km) | −1.03 | 1.64 | +380 | 918 | +22.5 |
| D02 (15 km) | −2.04 | 2.37 | +831 | 2071 | +25.3 |
| D01 (~75 km) | −3.07 | 3.33 | +997 | 2243 | −2.8 |
| CanESM2 | −1.94 | 2.59 | +715 | 980 | +47.6 |
| CanRCM4 | −3.60 | 3.80 | +508 | 1272 | +48.4 |

All five rows are positively biased in elevation (coarser grid → higher
orography). Cold temperature bias grows with coarseness. The coarse Can* rows
show a large dry precipitation bias (median ≈ −90 %, ≈ −28 mm/day).

## Notes

- Scripts in `DataPrep/`, `mean_changes/`, `extremes/`, and `CompoundEvents/`
  contain **absolute paths** to the original data on the HPC system where the
  analysis was run; edit the path constants at the top of each script before
  re-running.
- Data files are not stored in this repository; each folder's README
  documents the expected inputs/outputs and the notebook-to-script
  relationships.
