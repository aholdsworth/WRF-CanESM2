# Evaluation — station/model-comparison skill code

All station- and model-vs-observation comparisons for the manuscript
*"Projecting future changes over the coastal Pacific Northwest: climate and
extremes from a convection-permitting model."* (1986–2005, WRF d03 driven by
CanESM2 vs D02/D01, CanRCM4, CanESM2). Each comparison lives in its own
subfolder; the shared skill-calculation modules are single copies (below).

## Subfolders

| subfolder | contents |
|-----------|----------|
| `climate_taylor/` | 3-panel climate Taylor diagram (T/PR/wind, 20-yr climatologies, 5 models) — `climate_all_taylor.py` + notebook, `taylor_figure.py`. |
| `bias_figure/` | D03 station-bias figure (climatology maps + gated bias dots): `bias_figure.ipynb` (plot) + `make_bias_figure.py` (compute) + `taylor_data.py` (station sets/annual series); the climatology builder `make_climo_d03.py` lives in `DataPrep/`. |
| `bias_comparisons/` | Multi-model station-bias maps (D02/D01/CanRCM4/CanESM2 rows × T/PR/wind cols): `bias_comparisons.ipynb` (fiddle) + `make_bias_comparisons.py` (twin). |
| `elevation_bias/` | Elevation-bias diagnostic grid + the **shared skill-calc core**: `seasonal_pooling.py` (Rule 1/2 pooling gate, `skill()`), `make_station_tables.py` (canonical 173-land-row table), `make_skill_tables.py` (skill-table projection), `make_r2_table.py` (report R² column), grid scripts + notebook, result tables in `results/`. |
| `buoy_validation/` | D03 buoy validation (8 tas, 7 wind; relaxed 75 %/80 % gates): `buoy_{tas,wind}_bias.py`, extraction/weight scripts, validated-table results in `results/`. |
| `climate_land/` | Land climatology table (pooled obs vs model, t/pr/wind × 5 models) → `report/climate_land_table.tex`. |
| `climate_buoy/` | Buoy climatology table (pooled obs vs model, t/wind × 5 models) → `report/climate_buoy_table.tex`. |
| `extremes_validation/` | 75th/95th-percentile (T & P) validation at nearly-complete stations (110): p95 skill table + Taylor; `extreme_p95.ipynb` + `make_extreme_p95_nb.py` (regenerator) + `test_p95_taylor.py`. |

## Shared skill-calc modules (single copies — do not duplicate)

| module | home | used by |
|--------|------|---------|
| `seasonal_pooling.py` | `elevation_bias/` | bias_figure, bias_comparisons, extremes_validation, make_skill_tables, make_r2_table, buoy scripts (via the workspace copy) |
| `taylor_data.py` | `bias_figure/` | make_bias_figure.py (station sets, WRF/Can* annual series) |
| `taylor_figure.py` | `climate_taylor/` | climate_all_taylor, p95 Taylor |

Cross-folder consumers add the sibling folder to `sys.path`
(`os.path.join(BASE, '..', 'elevation_bias')` etc.) — they never copy the
module.

## Conventions (all subfolders)

- **Rule 1 (annuals):** direct day→year mean; year valid only if ≥ 90 % of
  expected days present; period value only if ≥ 15 of 20 years survive.
- **Rule 2 (seasons):** meteorological DJF/MAM/JJA/SON; valid only if all 3
  months present (no day-fraction gate).
- **Like-for-like:** model restricted to the obs surviving years.
- **Buoy-only relaxed gate:** day ≥ 75 % of hours, year ≥ 80 % of days AND
  data in ≥ 10 of 12 months (buoy scripts only; land is strict 90 %).
- Bias = model − obs (°C for T; % for PR/wind in skill tables; m/s for buoy
  wind in the bias figures).
- Canonical land set: the 173-row `station_tables.pkl` (100 USC + 9 USW +
  ECCC + BCH; no buoys); per-variable finite sets via the D03 finite-bias
  columns (t 146 / pr 159 / wind 15).

## Data (not in this repo)

Everything runs from the analysis working tree (the repo's parent):
`data/` (wrf_canon, cdo_extractions, obs), `pickles/` (canonical caches),
`files/` (station lists), `plotting_files/` (geo_em, bathymetry, orog),
`report/` (tex outputs). Data extraction lives in `DataPrep/`.
Run with the `py_2024` environment and
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` (OpenBLAS thread-spin gotcha).
