# Buoy wind-speed validation — results for report (2026-09-23)

**Mathematical description of methods:** see
`files/buoy_wind_methods_math.md` (height scaling, 90 % missing-data
cascade, bilinear extraction, skill metrics, and pooling — each given in
plain text and LaTeX).

Source: `buoy_wind_bias.py --year-gate 3` (relaxed) → `files/buoy_wind_bias{,_pooled}.csv`.
All five model columns now use consistently bilinear-extracted 10 m wind speed;
the WRF D01 column was re-extracted on 2026-09-23 (see DATA_INVENTORY.md M12)
after a transposed-grid bug in the first D01 weight generator was found and
fixed (D01 pooled bias moved from a spurious +23.5 % to −10.8 %).

## 1. Description (methods + results)

Twelve North Pacific / BC-coast buoys with 10-min wind-vector records
(1986–2005) were used to validate 10 m wind speed from five model columns:
WRF downscalings of CanESM2 on nested domains D01 (outer nest, ~75 km), D02
(middle nest, ~15 km) and D03 (inner nest, ~3 km), CanESM2 itself at its
native ~100 km resolution (the parent GCM, the coarsest column), and
CanRCM4 (~25 km, a comparable downscaling of CanESM2). Observed wind speed is the magnitude
of the two-component anemometer vector at 5 m, scaled up to 10 m with the
logarithmic wind profile, z0 = 2 mm (factor 1.0684, applied to the
observations only); model values are used natively at 10 m. Both records are
reduced to a common set of years with a 90 % missing-data cascade (a day is
kept if ≤10 % of its hours are present and no intra-day gap exceeds 3 h; a
month if ≤10 % of its days survive; a year if ≥3 of its 4 seasons pass), and
model annuals are computed on exactly the same surviving years as the
observations (apples-to-apples). Per buoy and model we report the period
means, bias % = 100·(model−obs)/obs, MAE and RMSE. The pooled headline
averages the per-buoy statistics over the core set of buoys with ≥8 shared
surviving years: 46005, 46041, 46131, 46146, 46204 and 46206.

Pooled results (6-buoy core): WRF D01 −10.8 % bias (MAE 1.35, RMSE 1.45 m s⁻¹),
WRF D02 +5.0 % (0.69, 0.82), WRF D03 +9.0 % (0.66, 0.78), CanESM2 −1.2 %
(0.96, 1.05) and CanRCM4 +16.2 % (1.11, 1.22). The mid-resolution WRF D02
shows the smallest absolute error of the five columns, and the 3 km WRF D03
the lowest MAE/RMSE; the coarsest column (CanESM2, ~100 km) is the most
accurate in mean bias (−1.2 %) but with larger error magnitude than the WRF
downscalings; CanRCM4 carries a large +16 % wind-speed bias. The D03 pooled
obs mean is 6.00 m s⁻¹ rather than 6.25 because buoy 46005 sits off the D03
nest domain and is dropped from that column only. Note on D01 placement: the
75 km D01 cells are ~0.7° across, so buoy positions fall up to ~0.35° from
the nearest cell centre; the D01 column therefore carries an extra sub-grid
placement uncertainty not present in D02/D03 (where buoy positions fall
well inside a single cell).

## 2. Table — pooled skill summary (6-buoy core, ≥8 shared years)

| Model | n buoys | Obs mean (m/s) | Model mean (m/s) | Bias (%) | MAE (m/s) | RMSE (m/s) |
|---|---|---|---|---|---|---|
| WRF D01 (~75 km)  | 6 | 6.25 | 5.58 | −10.8 | 1.35 | 1.45 |
| WRF D02 (~15 km)  | 6 | 6.25 | 6.56 |  +5.0 | 0.69 | 0.82 |
| WRF D03 (~3 km)   | 5 | 6.00 | 6.54 |  +9.0 | 0.66 | 0.78 |
| CanESM2 (~100 km) | 6 | 6.25 | 6.18 |  −1.2 | 0.96 | 1.05 |
| CanRCM4 (~25 km)  | 6 | 6.25 | 7.26 | +16.2 | 1.11 | 1.22 |

*WRF D03 core is 5 buoys (46005 excluded — off the D03 nest domain); its obs
mean is therefore over the 5 in-domain buoys.*

## 3. Table — per-buoy detail (relaxed year gate: ≥3 of 4 seasons/yr)

| Buoy | Model | yrs | Obs (m/s) | Model (m/s) | Bias (%) | MAE (m/s) | RMSE (m/s) |
|---|---|---|---|---|---|---|---|
| 46005 | WRF D01 | 10 | 7.48 | 7.87 |  +5.2 | 0.53 | 0.66 |
| 46005 | WRF D02 | 10 | 7.48 | 7.96 |  +6.5 | 0.58 | 0.72 |
| 46005 | WRF D03 |  — | 7.48 |  —   |  —    |  —   |  —   |
| 46005 | CanESM2 | 10 | 7.48 | 9.19 | +22.8 | 1.71 | 1.79 |
| 46005 | CanRCM4 | 10 | 7.48 | 8.55 | +14.4 | 1.08 | 1.20 |
| 46029 | WRF D01 |  6 | 5.84 | 5.64 |  −3.4 | 0.48 | 0.56 |
| 46029 | WRF D02 |  6 | 5.84 | 6.61 | +13.2 | 0.77 | 0.96 |
| 46029 | WRF D03 |  6 | 5.84 | 6.54 | +11.9 | 0.70 | 0.89 |
| 46029 | CanESM2 |  6 | 5.84 | 5.81 |  −0.5 | 0.50 | 0.60 |
| 46029 | CanRCM4 |  6 | 5.84 | 7.88 | +34.9 | 2.04 | 2.18 |
| 46041 | WRF D01 |  8 | 5.82 | 4.83 | −17.0 | 1.00 | 1.09 |
| 46041 | WRF D02 |  8 | 5.82 | 6.32 |  +8.5 | 0.49 | 0.70 |
| 46041 | WRF D03 |  8 | 5.82 | 6.36 |  +9.2 | 0.53 | 0.73 |
| 46041 | CanESM2 |  8 | 5.82 | 5.58 |  −4.1 | 0.40 | 0.45 |
| 46041 | CanRCM4 |  8 | 5.82 | 7.39 | +26.9 | 1.56 | 1.64 |
| 46050 | WRF D01 |  — | 6.15 |  —   |  —    |  —   |  —   |  (no d01 file yet — see note)
| 46050 | WRF D02 |  3 | 6.15 | 6.57 |  +6.8 | 0.42 | 0.49 |
| 46050 | WRF D03 |  3 | 6.15 | 6.53 |  +6.2 | 0.38 | 0.47 |
| 46050 | CanESM2 |  3 | 6.15 | 5.56 |  −9.5 | 0.59 | 0.70 |
| 46050 | CanRCM4 |  3 | 6.15 | 7.75 | +26.1 | 1.60 | 1.62 |
| 46087 | all     |  1 | 5.08 |  —   |  —    |  —   |  —   |
| 46088 | all     |  1 | 4.82 |  —   |  —    |  —   |  —   |
| 46089 | all     |  0 |  —   |  —   |  —    |  —   |  —   |
| 46131 | WRF D01 | 10 | 5.16 | 3.07 | −40.4 | 2.08 | 2.11 |
| 46131 | WRF D02 | 10 | 5.16 | 4.29 | −16.7 | 0.86 | 0.96 |
| 46131 | WRF D03 | 10 | 5.16 | 5.52 |  +7.0 | 0.48 | 0.58 |
| 46131 | CanESM2 | 10 | 5.16 | 4.55 | −11.8 | 0.62 | 0.71 |
| 46131 | CanRCM4 | 10 | 5.16 | 6.17 | +19.7 | 1.01 | 1.11 |
| 46132 | WRF D01 |  3 | 7.69 | 7.86 |  +2.2 | 0.42 | 0.54 |
| 46132 | WRF D02 |  3 | 7.69 | 8.23 |  +7.0 | 0.54 | 0.74 |
| 46132 | WRF D03 |  3 | 7.69 | 8.29 |  +7.8 | 0.60 | 0.78 |
| 46132 | CanESM2 |  3 | 7.69 | 7.75 |  +0.7 | 0.53 | 0.60 |
| 46132 | CanRCM4 |  3 | 7.69 | 8.81 | +14.6 | 1.12 | 1.32 |
| 46134 | all     |  1 | 3.00 |  —   |  —    |  —   |  —   |
| 46146 | WRF D01 | 10 | 5.22 | 2.47 | −52.7 | 2.75 | 2.77 |
| 46146 | WRF D02 | 10 | 5.22 | 5.88 | +12.6 | 0.66 | 0.75 |
| 46146 | WRF D03 | 10 | 5.22 | 6.05 | +16.0 | 0.84 | 0.91 |
| 46146 | CanESM2 | 10 | 5.22 | 3.85 | −26.2 | 1.36 | 1.42 |
| 46146 | CanRCM4 | 10 | 5.22 | 5.09 |  −2.5 | 0.36 | 0.45 |
| 46204 | WRF D01 |  9 | 7.69 | 7.86 |  +2.2 | 0.52 | 0.67 |
| 46204 | WRF D02 |  9 | 7.69 | 7.85 |  +2.0 | 0.52 | 0.67 |
| 46204 | WRF D03 |  9 | 7.69 | 7.84 |  +1.9 | 0.52 | 0.67 |
| 46204 | CanESM2 |  9 | 7.69 | 7.00 |  −9.0 | 0.74 | 0.87 |
| 46204 | CanRCM4 |  9 | 7.69 | 8.47 | +10.1 | 0.88 | 1.01 |
| 46206 | WRF D01 |  8 | 6.13 | 7.36 | +20.0 | 1.22 | 1.37 |
| 46206 | WRF D02 |  8 | 6.13 | 7.07 | +15.3 | 1.00 | 1.12 |
| 46206 | WRF D03 |  8 | 6.13 | 6.96 | +13.4 | 0.92 | 1.03 |
| 46206 | CanESM2 |  8 | 6.13 | 6.89 | +12.3 | 0.90 | 1.04 |
| 46206 | CanRCM4 |  8 | 6.13 | 7.89 | +28.7 | 1.76 | 1.92 |
| 46207 | WRF D01 |  6 | 7.71 | 7.93 |  +2.9 | 0.54 | 0.62 |
| 46207 | WRF D02 |  6 | 7.71 | 8.02 |  +4.1 | 0.54 | 0.66 |
| 46207 | WRF D03 |  6 | 7.71 | 8.04 |  +4.3 | 0.54 | 0.67 |
| 46207 | CanESM2 |  6 | 7.71 | 8.38 |  +8.7 | 0.69 | 0.99 |
| 46207 | CanRCM4 |  6 | 7.71 | 8.64 | +12.0 | 0.93 | 1.12 |

*46005/WRF D03: off the D03 nest domain → excluded from that column only.*

### Notes
- 46050, 46087, 46088, 46089 (the four "new NOAA" buoys): the d01
  re-extraction covered the 7 ECCC + the 3 legacy NOAA (46005/46029/46041)
  buoy files. Their d02/d03 files exist and are correct; only their d01
  columns are pending a follow-up extraction from `wind_d01_hourly.nc` (all
  four locations are inside the D01 grid). They do not enter the 6-buoy
  pooled core, so the headline table is unaffected.
- Metrics require ≥2 shared years; 46087/46088/46134 (1 yr each) are shown
  without metrics; 46089 has no surviving years.
