# climate_land — land climatology table (pooled obs vs model, 1986–2005)

Publication table `report/climate_land_table.tex` (+ `figures/climate/climate_land_map.png`):
per-variable (t / pr / wind) × per-model (D3, D15, D75, CanESM2, CanRCM4) pooled
obs-vs-model climatology — obs / bias (model − obs) / RMSE / R² — for the
canonical 173-station land set (t 146 / pr 159 / wind 15 finite; buoys excluded).

## Files

| file | role |
|------|------|
| `climate_land.ipynb` | Self-contained notebook: reads the cached annual series, pools (station × year), writes the table + map. Run top-to-bottom. |

## Inputs (not in this repo)

Run from the main verification workspace (`NOTEBOOKS/`), which provides:

- `pickles/taylor_annuals.pkl` — per-station × per-year (1986–2005) obs +
  model annuals for all 5 models (cache version 10).
- `pickles/station_tables.pkl` — canonical 173-land-row table (station
  coordinates + finite set per variable).

## Outputs (git-ignored)

- `report/climate_land_table.tex` — the report table body
- `figures/climate/climate_land_map.png`
