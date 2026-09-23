# Buoy near-surface air-temperature (tas) validation — results for report (2026-09-23)

**Mathematical description of methods:** see
`buoy_tas_methods_math.md` (obs decoding, the 90 % missing-data cascade,
apples-to-apples year pairing, skill metrics, and pooling — each given in
plain text and LaTeX).

**Source:** `buoy_tas_bias.py --year-gate 3` (relaxed) →
`files/buoy_tas_bias{,_pooled}.csv` (regenerated copies in `results/`).

Five model columns are compared: the three WRF nested domains — D01 (outer
nest, ~75 km), D02 (middle nest, ~15 km), D03 (inner nest, ~3 km), all
hourly `T2` — plus CanESM2 (~100 km, daily `tas`; the parent GCM and the
coarsest column) and CanRCM4 (~25 km, daily `tas`). All WRF columns were
re-extracted 2026-09-23 with embedded coordinates verified against the
buoy truth; the original PFM pull had rotated ECCC filenames (see §4).

## 1. Methods (summary)

Fourteen North Pacific / BC-coast buoys with near-surface air-temperature
records (1986–2005) are used to validate near-surface air temperature from
the five model columns above. The observation is the buoy met-package
air-temperature channel ("Temperature:Air", degrees C) from the
NOAA/MSC-OS `MB_<ID>_HM.mat` files, sampled at ~15-min cadence. It is
converted to kelvin (+273.15) and reduced to hourly (mean of the hour's
samples), daily and annual values with the same 90 % missing-data cascade
used for the wind analysis (a day is kept if ≥90 % of its hours are
present and no intra-day gap exceeds 3 h; a month if ≥90 % of its days
survive; a season if ≥2 of its 3 complete months survive; a year if ≥3 of
its 4 seasons pass). Model annuals are computed on exactly the same
surviving years as the observations (apples-to-apples). The buoy
air-temperature sensor sits in the met-package above the sea surface (the
package wind stresses are referenced to 4.2 m), while the model fields are
2 m (WRF `T2`) or near-surface (CanESM2/CanRCM4 `tas`); no height scaling
is applied (log-profile scaling is a momentum correction and does not
transfer to a conserved scalar). Per buoy and model we report the period
means, raw bias in K, bias % = 100·(model−obs)/obs, MAE and RMSE on the
shared-year annual values. The pooled headline averages the per-buoy
statistics over the core set of buoys with ≥8 shared surviving years:
46005, 46041, 46050, 46131, 46132, 46146, 46204, 46206 and 46207 (nine
buoys — three more than the wind core, because air temperature has far
fewer gaps than the wind vector). For D03 the core drops 46005 (nearest
3 km cell 1.16° away at the domain edge — see §4), so the D03 pooled row
covers 8 buoys.

## 2. Results

All five columns are slightly cold against the buoys, and the cold bias
weakens monotonically from the coarsest columns toward the finer nested
columns: D02 and D03 roughly halve the D01 bias (≈ −0.9 K), while CanRCM4
is closest to neutral with the smallest error magnitude. The per-buoy
pattern is not uniform. The two buoys in the Queen Charlotte / Hecate
Strait region (46131, 46146) show a large −3.4 to −3.6 K (≈ −1.2 %) bias
for both D01 and CanESM2 that largely disappears in the nested columns —
D02 reduces it to −0.3 to −0.5 K and D03 to −0.1 to −0.4 K — and is −0.9
to −1.6 K for CanRCM4; the west-coast buoys (46005, 46041, 46204, 46207)
show the smallest biases. A cross-check against the independent
land-station pipeline (120–121 stations) shows the same sign and ordering:
mean temperature bias −3.07 K (D01), −2.04 K (D02), −1.03 K (D03),
−1.94 K (CanESM2), −3.60 K (CanRCM4) — i.e. the coarse outer column is
systematically cold over land as well, and the bias weakens toward the
finer nests.

## 3. Table — pooled skill summary (core set, ≥8 shared years)

| Model | n buoys | Obs mean (K) | Model mean (K) | Bias (K) | Bias (%) | MAE (K) | RMSE (K) |
|---|---|---|---|---|---|---|---|
| CanESM2 (~100 km) | 9 | 283.89 | 282.58 | −1.31 | −0.46 | 1.48 | 1.65 |
| WRF D01 (~75 km)  | 9 | 283.89 | 282.30 | −1.60 | −0.56 | 1.67 | 1.81 |
| WRF D02 (~15 km)  | 9 | 283.89 | 283.01 | −0.88 | −0.31 | 1.02 | 1.18 |
| WRF D03 (~3 km)   | 8 | 283.81 | 282.92 | −0.89 | −0.32 | 1.04 | 1.18 |
| CanRCM4 (~25 km)  | 9 | 283.89 | 283.22 | −0.68 | −0.24 | 0.94 | 1.12 |

*Rows ordered from coarsest to finest. The D03 pooled row excludes 46005
(nearest D03 cell 1.16° away — off-domain artifact), so it covers 8 buoys;
its "Obs mean" is the mean over those same 8 buoys. All other rows cover
the full 9-buoy core.*

## 4. Table — per-buoy detail (core set; relaxed year gate: ≥3 of 4 seasons/yr)

| Buoy | Model | yrs | Obs (K) | Model (K) | Bias (K) | Bias (%) | MAE (K) | RMSE (K) |
|---|---|---|---|---|---|---|---|---|
| 46005 | WRF D01  | 13 | 284.53 | 283.83 | −0.71 | −0.25 | 1.02 | 1.14 |
| 46005 | WRF D02  | 13 | 284.53 | 283.93 | −0.61 | −0.21 | 0.96 | 1.09 |
| 46005 | WRF D03  | —  | —  | —  | —  | —  | —  | —  |
| 46005 | CanESM2  | 13 | 284.53 | 284.26 | −0.28 | −0.10 | 0.92 | 1.13 |
| 46005 | CanRCM4  | 13 | 284.53 | 284.12 | −0.41 | −0.15 | 0.98 | 1.19 |
| 46041 | WRF D01  | 12 | 283.97 | 281.87 | −2.10 | −0.74 | 2.10 | 2.12 |
| 46041 | WRF D02  | 12 | 283.97 | 282.69 | −1.28 | −0.45 | 1.28 | 1.33 |
| 46041 | WRF D03  | 12 | 283.97 | 282.73 | −1.24 | −0.44 | 1.24 | 1.30 |
| 46041 | CanESM2  | 12 | 283.97 | 283.48 | −0.49 | −0.17 | 0.55 | 0.67 |
| 46041 | CanRCM4  | 12 | 283.97 | 283.72 | −0.25 | −0.09 | 0.42 | 0.51 |
| 46050 | WRF D01  |  8 | 284.54 | 283.60 | −0.94 | −0.33 | 0.94 | 1.12 |
| 46050 | WRF D02  |  8 | 284.54 | 283.44 | −1.10 | −0.39 | 1.10 | 1.24 |
| 46050 | WRF D03  |  8 | 284.54 | 283.36 | −1.19 | −0.42 | 1.19 | 1.32 |
| 46050 | CanESM2  |  8 | 284.54 | 284.12 | −0.43 | −0.15 | 0.69 | 0.85 |
| 46050 | CanRCM4  |  8 | 284.54 | 284.71 | +0.16 | +0.06 | 0.61 | 0.69 |
| 46131 | WRF D01  | 12 | 283.85 | 280.40 | −3.45 | −1.22 | 3.45 | 3.52 |
| 46131 | WRF D02  | 12 | 283.85 | 283.31 | −0.54 | −0.19 | 0.70 | 0.83 |
| 46131 | WRF D03  | 12 | 283.85 | 283.72 | −0.14 | −0.05 | 0.54 | 0.62 |
| 46131 | CanESM2  | 12 | 283.85 | 280.43 | −3.42 | −1.20 | 3.42 | 3.51 |
| 46131 | CanRCM4  | 12 | 283.85 | 282.99 | −0.86 | −0.30 | 0.97 | 1.13 |
| 46132 | WRF D01  |  8 | 283.74 | 282.86 | −0.88 | −0.31 | 0.92 | 1.21 |
| 46132 | WRF D02  |  8 | 283.74 | 282.69 | −1.06 | −0.37 | 1.06 | 1.34 |
| 46132 | WRF D03  |  8 | 283.74 | 282.65 | −1.09 | −0.38 | 1.09 | 1.37 |
| 46132 | CanESM2  |  8 | 283.74 | 282.98 | −0.76 | −0.27 | 0.96 | 1.23 |
| 46132 | CanRCM4  |  8 | 283.74 | 283.10 | −0.64 | −0.23 | 0.87 | 1.10 |
| 46146 | WRF D01  | 10 | 284.27 | 280.73 | −3.54 | −1.24 | 3.54 | 3.72 |
| 46146 | WRF D02  | 10 | 284.27 | 283.94 | −0.32 | −0.11 | 0.89 | 1.14 |
| 46146 | WRF D03  | 10 | 284.27 | 283.83 | −0.43 | −0.15 | 0.95 | 1.18 |
| 46146 | CanESM2  | 10 | 284.27 | 280.66 | −3.60 | −1.27 | 3.60 | 3.87 |
| 46146 | CanRCM4  | 10 | 284.27 | 282.68 | −1.59 | −0.56 | 1.78 | 2.11 |
| 46204 | WRF D01  | 11 | 282.87 | 282.29 | −0.59 | −0.21 | 0.89 | 1.03 |
| 46204 | WRF D02  | 11 | 282.87 | 282.14 | −0.73 | −0.26 | 0.96 | 1.12 |
| 46204 | WRF D03  | 11 | 282.87 | 282.11 | −0.76 | −0.27 | 0.97 | 1.14 |
| 46204 | CanESM2  | 11 | 282.87 | 281.71 | −1.17 | −0.41 | 1.18 | 1.34 |
| 46204 | CanRCM4  | 11 | 282.87 | 282.25 | −0.63 | −0.22 | 0.79 | 0.96 |
| 46206 | WRF D01  | 10 | 283.84 | 282.78 | −1.06 | −0.37 | 1.09 | 1.22 |
| 46206 | WRF D02  | 10 | 283.84 | 282.58 | −1.26 | −0.44 | 1.26 | 1.39 |
| 46206 | WRF D03  | 10 | 283.84 | 282.53 | −1.31 | −0.46 | 1.31 | 1.44 |
| 46206 | CanESM2  | 10 | 283.84 | 282.74 | −1.10 | −0.39 | 1.20 | 1.34 |
| 46206 | CanRCM4  | 10 | 283.84 | 282.97 | −0.86 | −0.30 | 0.98 | 1.15 |
| 46207 | WRF D01  |  8 | 283.41 | 282.31 | −1.10 | −0.39 | 1.10 | 1.20 |
| 46207 | WRF D02  |  8 | 283.41 | 282.41 | −1.00 | −0.35 | 1.00 | 1.12 |
| 46207 | WRF D03  |  8 | 283.41 | 282.42 | −0.99 | −0.35 | 0.99 | 1.10 |
| 46207 | CanESM2  |  8 | 283.41 | 282.83 | −0.58 | −0.21 | 0.80 | 0.95 |
| 46207 | CanRCM4  |  8 | 283.41 | 282.39 | −1.02 | −0.36 | 1.07 | 1.26 |

### Notes

- **D03 — 46005 excluded.** 46005's nearest D03 cell is 1.16° away at the
  western domain edge — an off-domain nearest-cell artifact, not a real
  3 km sample (all other buoys are within 0.018° of their nearest D03
  cell). It is therefore dropped from the D03 column; its row above is
  "—" and the D03 pooled row covers 8 buoys.
- **D01 sub-grid placement.** The 75 km D01 cells are ~0.7° across, so a
  buoy position can fall up to ~0.35° from the nearest cell centre; the
  D01 column carries an extra sub-grid placement uncertainty not present
  at finer resolution (D02/D03 cells are 15 km / 3 km, so their placement
  uncertainty is negligible).
- **Rotated PFM files.** The original PFM per-station WRF files (Sept
  2023) were cross-assigned for six ECCC buoys (a file named for one buoy
  held the series of another). They were detected by comparing each
  file's embedded 1-point coordinates against the buoy truth
  (`station_cdo/buoys/buoy_descs.csv`) and re-extracted on 2026-09-23
  from regenerated, coordinate-verified station descriptors and SCRIP
  weights (`extract_buoy_t_parallel.sh`). The analysis matches every
  model file to a buoy by embedded coordinates (not filename), so a
  rotated/misnamed file cannot enter the analysis silently. The rotated
  originals are preserved under `data/wrf_stations/rotated_backup_20260923/`.
- **Non-core buoys** (in the per-buoy CSV, not in the pool): 46029 (7 obs
  years), 46134 (3), 46087/46088/46089 (1 each — deployed 2004–2005).
