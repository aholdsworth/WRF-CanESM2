# Buoy near-surface air-temperature (tas) validation — results for report (2026-09-23)

**Mathematical description of methods:** see
`buoy_tas_methods_math.md` (obs decoding, the 90 % missing-data cascade,
apples-to-apples year pairing, skill metrics, and pooling — each given in
plain text and LaTeX).

Source: `buoy_tas_bias.py --year-gate 3` (relaxed) →
`files/buoy_tas_bias{,_pooled}.csv` (regenerated copies in `results/`).
Three model columns: WRF D01 (outer nest, ~75 km, hourly T2 re-extracted
2026-09-23 with verified coordinates), CanESM2 (~100 km, daily `tas`, the
parent GCM and coarsest column) and CanRCM4 (~25 km, daily `tas`).
WRF D02 (15 km) / D03 (3 km) are pending — no per-buoy `t` extractions
exist yet.

## 1. Description (methods + results)

Fourteen North Pacific / BC-coast buoys with near-surface air-temperature
records (1986–2005) were used to validate near-surface air temperature from
three model columns: WRF D01 (outer nest, ~75 km), CanESM2 at its native
~100 km resolution (the parent GCM, the coarsest column), and CanRCM4
(~25 km, a comparable downscaling of CanESM2). The observation is the buoy
met-package air-temperature channel ("Temperature:Air", degrees C) from the
NOAA/MSC-OS `MB_<ID>_HM.mat` files, sampled at ~15-min cadence; it is
converted to kelvin (+273.15) and reduced to hourly (mean of the hour's
samples), daily and annual values with the same 90 % missing-data cascade
used for the wind analysis (a day is kept if ≥90 % of its hours are present
and no intra-day gap exceeds 3 h; a month if ≥90 % of its days survive; a
season if ≥2 of its 3 complete months survive; a year if ≥3 of its 4
seasons pass), and model annuals are computed on exactly the same surviving
years as the observations (apples-to-apples). The buoy air-temperature
sensor sits in the met-package above the sea surface (the package wind
stresses are referenced to 4.2 m), while the model fields are 2 m (WRF T2)
or near-surface (CanESM2/CanRCM4 `tas`); no height scaling is applied
(log-profile scaling is a momentum correction and does not transfer to a
conserved scalar). Per buoy and model we report the period means, raw bias
in K, bias % = 100·(model−obs)/obs, MAE and RMSE on the shared-year annual
values. The pooled headline averages the per-buoy statistics over the core
set of buoys with ≥8 shared surviving years: 46005, 46041, 46050, 46131,
46132, 46146, 46204, 46206 and 46207 (nine buoys — three more than the wind
core, because air temperature has far fewer gaps than the wind vector).

Pooled results (9-buoy core, K): WRF D01 bias −1.60 K (−0.56 %; MAE 1.67,
RMSE 1.81 K), CanESM2 −1.31 K (−0.46 %; 1.48, 1.65) and CanRCM4 −0.68 K
(−0.24 %; 0.94, 1.12). All three columns are slightly cold against the
buoys, and the cold bias weakens monotonically from the coarsest column
(D01) to the ~25 km CanRCM4, which is closest to neutral and has the
smallest error magnitude. The per-buoy pattern is not uniform: the two
buoys in the Queen Charlotte / Hecate Strait region (46131, 46146) show a
large −3.4 to −3.6 K (≈ −1.2 %) bias for both D01 and CanESM2 but only
−0.9 to −1.6 K for CanRCM4, whereas the west-coast buoys (46005, 46041,
46204, 46207) show the smallest biases. Note on D01 placement: the 75 km
D01 cells are ~0.7° across, so buoy positions fall up to ~0.35° from the
nearest cell centre; the D01 column therefore carries an extra sub-grid
placement uncertainty not present at finer resolution. A cross-check
against the independent land-station pipeline (120–121 stations) shows the
same sign and ordering: mean temperature bias −3.07 K (D01), −2.04 K (D02),
−1.03 K (D03), −1.94 K (CanESM2), −3.60 K (CanRCM4) — i.e. the coarse outer
column is systematically cold over land as well.

## 2. Table — pooled skill summary (9-buoy core, ≥8 shared years)

| Model | n buoys | Obs mean (K) | Model mean (K) | Bias (K) | Bias (%) | MAE (K) | RMSE (K) |
|---|---|---|---|---|---|---|---|
| WRF D01 (~75 km)  | 9 | 283.89 | 282.30 | −1.60 | −0.56 | 1.67 | 1.81 |
| CanESM2 (~100 km) | 9 | 283.89 | 282.58 | −1.31 | −0.46 | 1.48 | 1.65 |
| CanRCM4 (~25 km)  | 9 | 283.89 | 283.22 | −0.68 | −0.24 | 0.94 | 1.12 |

## 3. Table — per-buoy detail (core set; relaxed year gate: ≥3 of 4 seasons/yr)

| Buoy | Model | yrs | Obs (K) | Model (K) | Bias (K) | Bias (%) | MAE (K) | RMSE (K) |
|---|---|---|---|---|---|---|---|---|
| 46005 | WRF D01  | 13 | 284.53 | 283.83 | −0.71 | −0.25 | 1.02 | 1.14 |
| 46005 | CanESM2  | 13 | 284.53 | 284.26 | −0.28 | −0.10 | 0.92 | 1.13 |
| 46005 | CanRCM4  | 13 | 284.53 | 284.12 | −0.41 | −0.15 | 0.98 | 1.19 |
| 46041 | WRF D01  | 12 | 283.97 | 281.87 | −2.10 | −0.74 | 2.10 | 2.12 |
| 46041 | CanESM2  | 12 | 283.97 | 283.48 | −0.49 | −0.17 | 0.55 | 0.67 |
| 46041 | CanRCM4  | 12 | 283.97 | 283.72 | −0.25 | −0.09 | 0.42 | 0.51 |
| 46050 | WRF D01  |  8 | 284.54 | 283.60 | −0.94 | −0.33 | 0.94 | 1.12 |
| 46050 | CanESM2  |  8 | 284.54 | 284.12 | −0.43 | −0.15 | 0.69 | 0.85 |
| 46050 | CanRCM4  |  8 | 284.54 | 284.71 | +0.16 | +0.06 | 0.61 | 0.69 |
| 46131 | WRF D01  | 12 | 283.85 | 280.40 | −3.45 | −1.22 | 3.45 | 3.52 |
| 46131 | CanESM2  | 12 | 283.85 | 280.43 | −3.42 | −1.20 | 3.42 | 3.51 |
| 46131 | CanRCM4  | 12 | 283.85 | 282.99 | −0.86 | −0.30 | 0.97 | 1.13 |
| 46132 | WRF D01  |  8 | 283.74 | 282.86 | −0.88 | −0.31 | 0.92 | 1.21 |
| 46132 | CanESM2  |  8 | 283.74 | 282.98 | −0.76 | −0.27 | 0.96 | 1.23 |
| 46132 | CanRCM4  |  8 | 283.74 | 283.10 | −0.64 | −0.23 | 0.87 | 1.10 |
| 46146 | WRF D01  | 10 | 284.27 | 280.73 | −3.54 | −1.24 | 3.54 | 3.72 |
| 46146 | CanESM2  | 10 | 284.27 | 280.66 | −3.60 | −1.27 | 3.60 | 3.87 |
| 46146 | CanRCM4  | 10 | 284.27 | 282.68 | −1.59 | −0.56 | 1.78 | 2.11 |
| 46204 | WRF D01  | 11 | 282.87 | 282.29 | −0.59 | −0.21 | 0.89 | 1.03 |
| 46204 | CanESM2  | 11 | 282.87 | 281.71 | −1.17 | −0.41 | 1.18 | 1.34 |
| 46204 | CanRCM4  | 11 | 282.87 | 282.25 | −0.63 | −0.22 | 0.79 | 0.96 |
| 46206 | WRF D01  | 10 | 283.84 | 282.78 | −1.06 | −0.37 | 1.09 | 1.22 |
| 46206 | CanESM2  | 10 | 283.84 | 282.74 | −1.10 | −0.39 | 1.20 | 1.34 |
| 46206 | CanRCM4  | 10 | 283.84 | 282.97 | −0.86 | −0.30 | 0.98 | 1.15 |
| 46207 | WRF D01  |  8 | 283.41 | 282.31 | −1.10 | −0.39 | 1.10 | 1.20 |
| 46207 | CanESM2  |  8 | 283.41 | 282.83 | −0.58 | −0.21 | 0.80 | 0.95 |
| 46207 | CanRCM4  |  8 | 283.41 | 282.39 | −1.02 | −0.36 | 1.07 | 1.26 |

*Non-core buoys (in the per-buoy CSV, not in the pool): 46029 (7 obs years),
46134 (3), 46087/46088/46089 (1 each — deployed 2004–2005).*
