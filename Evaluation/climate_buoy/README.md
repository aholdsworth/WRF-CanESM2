# climate_buoy — buoy climatology table (pooled obs vs model, 1986–2005)

Publication table `report/climate_buoy_table.tex` (+ `figures/climate/climate_buoy_map.png`):
per-variable (t / wind) × per-model (D3, D15, D75, CanESM2, CanRCM4) pooled
obs-vs-model climatology — obs / bias (model − obs) / RMSE / R² — for the
validated D03 buoy set (8 temperature buoys `TAS_BUOYS_D03`, 7 wind buoys
`WIND_BUOYS_D03`; no pr). Obs and model values come from the validated buoy
table pickles (never re-parsed from `data/wrf_canon` 46xxx files).

## Files

| file | role |
|------|------|
| `climate_buoy.ipynb` | Self-contained notebook: reads the cached annual series, pools (buoy × year), writes the table + map. Run top-to-bottom. |

## Inputs (not in this repo)

Run from the main verification workspace (`NOTEBOOKS/`), which provides:

- `pickles/taylor_annuals.pkl` — per-station × per-year (1986–2005) obs +
  model annuals for all 5 models (cache version 10; buoy entries from the
  validated buoy tables).
- `pickles/buoy_wind_tables.pkl` — validated buoy truth coordinates + model
  period values for the wind set.

## Outputs (git-ignored)

- `report/climate_buoy_table.tex` — the report table body
- `figures/climate/climate_buoy_map.png`
