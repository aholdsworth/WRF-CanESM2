# bias_comparisons — multi-model station bias maps (D02/D01/CanRCM4/CanESM2 vs obs)

Publication figure `figures/bias_comparisons.png`: 4 rows (models, top →
bottom: D02, D01, CanRCM4, CanESM2) × 3 columns (temperature bias °C,
precipitation bias %, wind bias %). Each panel = the row model's own
topography (terrain colormap, 0–3000 m) with per-station bias dots
(model − obs, gated Rule-1 reference) on the column's diverging scale
(`RdBu_r` / `BrBG` / `PuOr`, one shared colorbar per column).

Canonical observation set (2026-09-29, identical to `bias_figure` /
Taylor / elevation-bias grid): land = all ECCC d03 + all BCH + NOAA-land
(USC*) + the 5 USW WA/OR wind stations (173-row canonical land table);
buoys pinned to the validated D03 set per variable — t 8
(`TAS_BUOYS_D03`), wind 7 (`WIND_BUOYS_D03`), pr none.

## Files

| file | role |
|------|------|
| `bias_comparisons.ipynb`        | FIDDLE notebook: reproduces `bias_comparisons.png` from `pickles/station_tables.pkl` (no re-extraction). KNOBS cell (map window, zoom mode, topo scale, ranges, dot size/edge, mask-ocean) is where you play; run top-to-bottom. |
| `make_bias_comparisons.py`      | Standalone script twin (older PlateCarree layout): reads `bias_figure.pkl` (gated Fig-1 station set + obs reference), `data/wrf_canon` (D02/D01 rows), `seasonal_pooling.cdo_station_value` (Can* rows), and the validated buoy tables; writes `figures/bias_comparisons.png`. |
| `seasonal_pooling.py`           | Shared pooling source of truth; used by `make_bias_comparisons.py` via `sp.cdo_station_value` for the CanESM2/CanRCM4 per-station period values (M9/M10 CDO files). |

The notebook also imports `WRFDomainLib` from the **repo root**
(`calc_wps_domain_info`, for the D03 Lambert box).

## Inputs (not in this repo)

Run from the main verification workspace (`NOTEBOOKS/`), which provides:

- `pickles/station_tables.pkl` — canonical 173-land-row table (all
  station data the notebook reads: per-dataset period means + elevations).
- `pickles/bias_figure.pkl` + `pickles/buoy_{tas,wind}_tables.pkl` —
  used by `make_bias_comparisons.py` (gated Fig-1 set, buoy obs/model).
- `data/wrf_canon/` — WRF per-station extractions, `{var}_{dom}_{AGENCY}_{ID}.nc`
  (D02/D01 rows); `data/` M9/M10 CDO per-station files (Can* rows).
- `plotting_files/` — `geo_em.d01.nc`–`geo_em.d03.nc` (WRF HGT_M
  topography) and `orog_CanESM2.nc` / `orog_CanRCM4.nc` (Can* orography).

## Run

```
# fiddle & re-render (fast, reads only station_tables.pkl):
jupyter bias_comparisons.ipynb        # or run top-to-bottom

# standalone re-render:
/home/amh001/space_fs7/software_2022/python/py_2024/bin/python make_bias_comparisons.py
```

Output (git-ignored): `figures/bias_comparisons.png`.
