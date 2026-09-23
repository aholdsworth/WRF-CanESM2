# Elevation-bias diagnostic grid

A 5 × 4 diagnostic figure and its tables for the manuscript
*"Projecting future changes over the coastal Pacific Northwest: climate and
extremes from a convection-permitting model."*

Each **row** is a model/domain and each **column** is a panel: (1) the row's
own model topography with station points coloured by that row's elevation
bias, then scatter panels of (2) temperature bias, (3) precipitation bias,
(4) wind-speed bias, all against **elevation bias (m)** on the y-axis.

* **Rows (5):** `D3` = WRF d03 (3 km, blue), `D15` = WRF d02 (15 km, orange),
  `D75` = WRF d01 (~75 km, green, emphasised), `CanRCM4` (purple, ~20–28 km
  rotated pole), `CanESM2` (red, ~100 km raw GCM).
* **Reference period:** 1986–2005 mean, 175 stations (164 land + 11 buoys).
* **Sign convention:** bias = model − observation (positive = model high).

The mathematical description (elevation/orography, the percent and absolute
precip biases, wind bias, the observations-side season gate, and skill
metrics) is in [`results/elev_bias_methods_math.md`](results/elev_bias_methods_math.md).
The report value tables are in
[`results/elev_bias_results_for_report.md`](results/elev_bias_results_for_report.md).

## The notebook

**[`elev_bias_grid.ipynb`](elev_bias_grid.ipynb)** is the walkthrough:
setup → load the station-table pickle → inspect the data → build the 5 × 4
grid (§5) → re-express the PR column in absolute units, mm/day (§5b and §5c,
which need no pickle rebuild) → optional per-dataset elevation-bias maps.
It is self-contained: it imports the companion scripts below for the
domain-file readers (`read_geo_em`, `read_orog`) and the WRF-domain helper
(`WRFDomainLib`, at the repo root) rather than re-implementing them.

Outputs (written to `<BASE>/figures/elevation_bias/`, PNGs are git-ignored):

| figure | cell |
|--------|------|
| `elevation_bias_grid.png` | §5 — the 5 × 4 percent-bias grid |
| `elevation_bias_grid_pr_abs.png` | §5b — PR bias in mm/day, 5 × 1 |
| `elevation_bias_grid_pr_abs_full.png` | §5c — full 5 × 4 with PR in mm/day |
| `obs_station_elevations.png` + 5 `*_elev_bias.png` | §6 — per-dataset maps |

## Scripts

| script | purpose |
|--------|---------|
| `make_station_tables.py` | Builds `pickles/station_tables.pkl` (dict of 5 DataFrames, 175 rows each) — per-station elevation + T/PR/wind bias for all 5 datasets. Reads the WRF per-station CDO extractions, the Can* station files, and the raw ECCC/NOAA/BCH observations. |
| `make_skill_tables.py` | Builds `pickles/skill_seasonal.pkl` / `skill_summary.csv` (.pkl) — the season-pooled per-station bias/MAE/RMSE and the 15-row summary. Uses the same obs-side gate as `seasonal_pooling.py`. |
| `make_elev_bias_grid.py` | Standalone script version of the §5 figure (same output the notebook writes). |
| `plot_elev_bias_revised.py` | Per-dataset elevation-bias **maps** (6 PNGs) and the shared data readers (`read_geo_em`, `read_orog`, the ECCC/BCH/NOAA/buoy loaders, `weighted_model_elevs`). Importing it renders no figures (plotting is under `if __name__ == '__main__'` / `_make_figures()`). |
| `seasonal_pooling.py` | The single source of truth for the observations-side season→year→period gate (90 % season completeness, ≥15 surviving years) and the `skill()` aggregation. Shared by `make_skill_tables.py` (and `mae_vs_eva.py` in the working tree). |

`WRFDomainLib.py` (repo root) supplies the WRF d01–d03 domain geometry
(`calc_wps_domain_info`, `reproject_corners`) used by the §6 per-dataset maps.

### Run order (to regenerate from raw data)

```
make_station_tables.py     # -> pickles/station_tables.pkl
make_skill_tables.py       # -> pickles/skill_seasonal.pkl, skill_summary.csv
make_elev_bias_grid.py     # -> figures/elevation_bias/elevation_bias_grid.png
# or:  jupyter nbconvert --execute elev_bias_grid.ipynb
```

The notebook only needs `pickles/station_tables.pkl` (and
`pickles/obs_pr_mm_annual.pkl` for the mm/day cells; it recomputes that cache
from raw obs if missing) plus the topography files in `plotting_files/` — so
a quick re-plot needs no raw-data rebuild.

## Inputs (not stored in this repo)

All paths are rooted at `BASE`, the analysis working tree. Set **`ELEV_BIAS_BASE`**
to relocate the tree (the scripts and the notebook all honour it); the default
is the HPC working tree used for the analysis:
`/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/NOTEBOOKS`.

| path (under `BASE`) | what |
|---------------------|------|
| `files/` | station lists: `ECCC_d03_stations.csv`, `BCH_d03_stations.csv`, `NOAA_d03_stations.csv`, `ECCC_buoys.csv`, `NOAA_buoys.csv` (coords + elevation) |
| `plotting_files/` | topography/geometry: `geo_em.d01/d02/d03.nc`, `orog_CanESM2.nc`, `orog_CanRCM4.nc`, `namelist.wps.txt` |
| `data/wrf_stations/` | WRF per-station CDO extractions (`t/pr/wind_d0X_*.nc`) |
| `data/eccc_obs/`, `data/noaa_obs/`, `data/bch_obs/` | raw observations (ECCC daily, NOAA GHCN-d, BCH p24) |
| `data/eccc_hourly_obs/` | ECCC hourly wind (for the wind bias) |
| `data/cdo_extractions/`, `data/canesm2_raw/`, `data/canrcm4/` | CanESM2 / CanRCM4 per-station extractions |
| `weight_files_gen/` | locally regenerated CDO bicubic weight files (Can* ECCC/BCH) |
| `station_cdo/eccc_bch_desc_map.json` | table-id → pipeline-sid map for ECCC/BCH weight lookup |

## Result tables

The committed tables live in `results/`:

* **`skill_summary.csv`** — the 15-row pooled seasonal skill summary
  (5 datasets × {t, pr, wind}: bias, MAE, RMSE, n_stations, n_years range),
  produced by `make_skill_tables.py`.
* **`elev_bias_results_for_report.md`** — the per-dataset bias/MAE/RMSE and
  elevation-bias summary tables (per-station medians/stds from
  `station_tables.pkl`, plus the absolute-PR mm/day table) for the report.
* **`elev_bias_methods_math.md`** — the equations (plain text + LaTeX).

A condensed version of the headline tables is also in the top-level
`README.md`.

## Caveats

* **WRF/Can* PR means** are inflated by a few extreme D75 ECCC sites
  (+2,000–7,000 %); the medians (≈ −13 / +73 / +145 %) are representative.
* **Coarse-grid over-elevation** is real physics, now shown at full magnitude
  by the extraction-consistent (bicubic) orography; the older nearest-GP
  "vertical line" stacking artefact is gone.
* **Wind** has the smallest station count (ECCC hourly-wind stations only).
* Full bug-fix log, verification gates and the NOAA-inches unit note are in
  the working-tree `elevation_bias_handoff.md` (not in this repo).
