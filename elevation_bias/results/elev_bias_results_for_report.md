# Elevation-bias diagnostic grid — result tables

Summary tables for the 5-dataset elevation-bias diagnostic
(`../elev_bias_grid.ipynb`, §5/§5b/§5c). 1986–2005 mean, 175 stations
(164 land + 11 buoys; buoys contribute elevation only via sea-level).

**Row ↔ dataset** (per the manuscript convention): D3 = WRF d03 (3 km, blue),
D15 = WRF d02 (15 km, orange), D75 = WRF d01 (75 km, green, emphasised),
CanRCM4 (purple), CanESM2 (red).

Statistics are per-station over the non-NaN rows of
`pickles/station_tables.pkl` (the legacy catalogue the notebook loads).
The pooled, year/season-gated skill summary is in `skill_summary.csv`
(from `make_skill_tables.py`).

## Elevation bias (m) — model orography − observed station elevation

| Row     | n   | median (m) | std (m) | mean (m) |
|---------|-----|-----------:|--------:|---------:|
| D3      | 164 |        41.6 |   165.5 |     89.8 |
| D15     | 164 |       170.3 |   300.6 |    247.2 |
| D75     | 164 |       304.6 |   318.1 |    331.7 |
| CanRCM4 | 164 |       242.3 |   348.8 |    294.3 |
| CanESM2 | 164 |       391.2 |   368.5 |    375.4 |

All rows are positively biased (coarser model grid → higher orography than the
observed valley station). The WRF rows use the **extraction-consistent**
(CDO bicubic-weight) orography; CanESM2/CanRCM4 likewise. `elev_method`
column records `bicubic`/`corner` (WRF) or `bicubic`/`nearest` (Can*) per row.

## Temperature bias (°C) — 1986–2005 mean

| Row     | n   | median (°C) | std (°C) | mean (°C) |
|---------|-----|------------:|---------:|----------:|
| D3      | 147 |        −0.96 |     1.16 |    −1.11 |
| D15     | 147 |        −1.80 |     1.64 |    −2.06 |
| D75     | 146 |        −2.86 |     1.93 |    −3.07 |
| CanRCM4 | 147 |        −2.95 |     2.19 |    −3.52 |
| CanESM2 | 147 |        −1.91 |     2.05 |    −1.85 |

Cold bias grows with model coarseness (D75/CanRCM4 coldest).

## Precipitation bias (%) — 100·(model−obs)/obs, 1986–2005 mean

| Row     | n   | median (%) | std (%) | mean (%) |
|---------|-----|-----------:|--------:|---------:|
| D3      | 164 |       −13.1 |  3116.2 |  1181.3 |
| D15     | 164 |        72.5 |  5433.5 |  2377.2 |
| D75     | 163 |       145.2 |  5692.1 |  2742.5 |
| CanRCM4 | 164 |       −92.6 |   136.2 |   −25.6 |
| CanESM2 | 164 |       −90.3 |   322.7 |    20.8 |

The WRF means/stds are inflated by a few extreme D75 ECCC sites (order
+2,000–7,000 %) — the medians (−13 / +73 / +145 %) are the representative
values. The coarse Can* rows show a large **dry** bias (median ≈ −90 %), which
the absolute-PR table below confirms spans ≈ −28 mm/day.

## Precipitation bias, absolute (mm/day) — model − obs, 1986–2005 mean

Re-expressed from the percent table via the exact identity
`abs = (pct/100)·obs` with `obs` the 1986–2005 annual-mean precip in mm/day
(annual mm/yr ÷ 365). This is what §5b/§5c of the notebook plot.

| Row     | n   | median (mm/day) | ≈ mm/yr |
|---------|-----|----------------:|--------:|
| D3      | 164 |            −4.0 |   −1,465 |
| D15     | 164 |            44.0 |  +16,056 |
| D75     | 163 |            63.3 |  +23,107 |
| CanRCM4 | 164 |           −27.8 |  −10,143 |
| CanESM2 | 164 |           −26.8 |   −9,774 |

## Wind-speed bias (%) — 100·(model−obs)/obs, hourly-mean m/s

| Row     | n   | median (%) | std (%) | mean (%) |
|---------|-----|-----------:|--------:|---------:|
| D3      |  18 |       30.03 |   50.30 |   39.00 |
| D15     |  17 |       33.23 |   47.13 |   40.57 |
| D75     |  17 |        6.43 |   48.24 |   22.90 |
| CanRCM4 |  18 |       61.45 |   59.38 |   68.01 |
| CanESM2 |  18 |       65.03 |   53.44 |   70.75 |

Wind has the smallest n (ECCC hourly-wind stations only).

## Pooled seasonal skill summary (from `make_skill_tables.py`)

`skill_summary.csv` (15 rows = 5 datasets × {t, pr, wind}) with the pooled
median-of-stations bias / MAE / RMSE and surviving-year counts:

| dataset | var  | bias | MAE | RMSE | n_st | n_years |
|---------|------|-----:|----:|-----:|-----:|--------:|
| D03     | t    | −1.03 °C | 1.46 | 1.64 | 121 | 15–20 |
| D03     | pr   | +380 % | 891 % | 918 % | 129 | 15–20 |
| D03     | wind | +22.5 % | 0.62 | 0.68 | 10 | 16–20 |
| D02     | t    | −2.04 °C | 2.21 | 2.37 | 121 | 15–20 |
| D02     | pr   | +831 % | 2035 % | 2071 % | 129 | 15–20 |
| D02     | wind | +25.3 % | 0.79 | 0.85 | 10 | 16–20 |
| D01     | t    | −3.07 °C | 3.20 | 3.33 | 120 | 15–20 |
| D01     | pr   | +997 % | 2207 % | 2243 % | 128 | 15–20 |
| D01     | wind | −2.8 % | 0.56 | 0.62 | 10 | 16–20 |
| CanESM2 | t    | −1.94 °C | 2.39 | 2.59 | 121 | 15–20 |
| CanESM2 | pr   | +715 % | 956 % | 980 % | 129 | 15–20 |
| CanESM2 | wind | +47.6 % | 1.16 | 1.21 | 10 | 16–20 |
| CanRCM4 | t    | −3.60 °C | 3.64 | 3.80 | 121 | 15–20 |
| CanRCM4 | pr   | +508 % | 1247 % | 1272 % | 129 | 15–20 |
| CanRCM4 | wind | +48.4 % | 1.33 | 1.37 | 10 | 16–20 |

(The pr % here is the *season-pooled* per-station bias from
`make_skill_tables.py`, i.e. the same gated convention as the skill table; it
differs from the legacy `station_tables` medians above by the year/season
gate and the pooled-vs-median-of-stations step.)

---

*Generated from `pickles/station_tables.pkl` and `pickles/skill_summary.csv`
(2026-09-23 build, post NOAA-inches fix and bicubic-elevation method).
See `../README.md` for inputs and `elev_bias_methods_math.md` for the
equations.*
