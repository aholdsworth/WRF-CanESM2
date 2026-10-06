# bias_figure — station bias figure (WRF D03 vs observations, 1986–2005)

Publication figure: 3 columns (T / PR / wind) × 2 rows for the WRF d03
domain — top row the 1986–2005 observation climatology maps, bottom row
model-vs-obs bias (annual + seasonal points, gated like-for-like annuals
per the Rule-1 pooling conventions). Stations: ECCC + BCH + NOAA land
(100 USC + 9 USW, 173 canonical land rows) + 7 D03 wind buoys + 8 D03
temperature buoys (buoys drawn as diamonds). Bias = model − obs (degC for
T; % for PR and wind).

## Files

| file | role |
|------|------|
| `bias_figure.ipynb`        | PLOT step: loads the precomputed pickles and draws the panels (run top-to-bottom, `py_2024` env). |
| `make_bias_figure.py`      | COMPUTE step (15–20 min): reads all obs + WRF d03 per-station files, computes everything, saves `pickles/bias_figure.pkl` + `pickles/bias_filter_comparison.csv`. |
| `seasonal_pooling.py`      | (shared, at `elevation_bias/seasonal_pooling.py`) — Rule 1/2 pooling gate + `skill()` aggregation; the scripts here import it from that one canonical location. |
| `taylor_data.py`           | Station-set + WRF/Can* per-station file readers and annual series builders used by `make_bias_figure.py` (canonical 173-land-row sets from `station_tables.pkl`). |
| `make_climo_d03.py`        | Builds `pickles/climo_d03_hist.pkl` — the gridded 1986–2005 d03 climatology (t degC, pr mm/day, wind m/s) the notebook shows in the top row. |

The notebook also imports `WRFDomainLib` (at the **repo root** — `sys.path.insert(0,
repo_root)`, i.e. the parent of this folder's sibling layout in the main
workspace) for `calc_wps_domain_info` / `reproject_corners` (Lambert D03 box).

## Inputs (not in this repo)

Run from the main verification workspace (`NOTEBOOKS/`), which provides:

- `pickles/station_tables.pkl` — canonical 173-land-row table (Rule-1 root
  producer; D03 finite-bias columns define the per-variable land station
  sets: t 146 / pr 159 / wind 15).
- `pickles/buoy_wind_tables.pkl`, `pickles/buoy_tas_tables.pkl`,
  `pickles/buoy_tas_d03_per_buoy.pkl` — validated D03 buoy obs+model
  (the bias figure's buoy points are exactly the published
  buoy-validation set; no buoy value is parsed from `data/wrf_canon`).
- `pickles/skill_seasonal.pkl` — published skill biases (Eva's numbers).
- `data/wrf_canon/` — WRF per-station extractions, `{var}_{dom}_{AGENCY}_{ID}.nc`.
- `data/` obs trees (NOAA GHCN-d CSVs, ECCC, BCH), `files/` station lists.
- `pickles/climo_d03_hist.pkl` — from `make_climo_d03.py` (reads the
  CF-compliant d03 daily files under
  `…/AI_WORKDIR/WRF_FILES/cf_compliant/data/historical/`).
- `plotting_files/` — `geo_em.d02.nc`, `geo_em.d03.nc` (terrain base +
  grid), `new_bath.nc` (bathymetry), `namelist.wps.txt` (Lambert params).

## Run order

```
make_climo_d03.py      # -> pickles/climo_d03_hist.pkl   (once)
make_bias_figure.py    # -> pickles/bias_figure.pkl       (delete the pkl to force)
# then run bias_figure.ipynb top-to-bottom (figures/bias/*)
```

Outputs (git-ignored): `figures/bias/*.png|pdf`.
