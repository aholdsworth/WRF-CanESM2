# Extremes validation — 75th/95th percentile (T & P) at nearly-complete stations

Obs-vs-model comparison of the **75th and 95th percentile of daily temperature
(T)** and **daily precipitation (P)** at land stations with **nearly-complete
records**, so that obs and model are like-for-like.  The 95th percentile is
the statistic of the repo `extremes/` bootstrap pipeline (daily values pooled
over the record, P "as-is" with dry days included); the 75th is a diagnostic
for the upper half of the distribution — it separates models whose p95 looks
good because their *whole distribution* (mean-like) fit is better from models
that specifically reproduce the upper tail (see Figure 5, p75 bias vs p95
bias per station).

## Station set

110 land stations (ECCC + BCH + NOAA) from the canonical
`<BASE>/pickles/station_tables.pkl` (`rule1-canonical-land-v1`) with Rule-1
surviving years **>= 19 of 20 for both T and P**.

## Method

- **Like-for-like:** the surviving-year set is defined by the obs series
  (Rule-1: a year is valid only if >= 90 % of its expected days are present);
  each model series (complete) is restricted to the same years.
- **Percentiles:** 75th and 95th percentile of the daily values pooled over
  all days of the shared surviving years.  P includes dry days (the "pr
  as-is" convention of the `extremes/` bootstrap pipeline).
- **Bias:** model percentile minus obs percentile (positive = model high).
- **p95 skill (Sections 11-12):** RMSE = root-mean-square of the
  per-station (model - obs) p95; R² = pooled squared Pearson correlation of
  the *annual* p95 series (per-station p95 of each calendar year, pooled over
  station x year over the shared surviving years) — the same "pooled R²"
  convention as the repo's Taylor / skill-table pipeline.  The p95 Taylor
  diagram (Figure 7) is drawn in the `taylor_3panel` style (correlation arc,
  RMS contours, auto-scaled std axis, single-row legend): each point is a
  dataset's pooled (station x year) p95 series (std, r vs pooled obs p95),
  black star = obs reference (refstd over the union of shared station-years).

Model datasets: `D3` = WRF d03, `D15` = WRF d02, `D75` = WRF d01,
`CanRCM4`, `CanESM2` (CDO per-station files under
`data/cdo_extractions/{canrcm4_noaa113,canesm2_raw_noaa113}`).

## Notebooks

- `extreme_p95.ipynb` — self-contained walkthrough: station selection, daily
  series build (cached), per-station p75/p95, the p95 skill table
  (bias / RMSE / R²) and a p95 Taylor diagram, plus five figures.  Imports
  the shared obs loaders / Rule-1 gate from `<BASE>/seasonal_pooling.py` and
  inlines the model-file readers (it does not import
  `make_station_tables.py`).

Run (single-threaded — see the OpenBLAS thread-spin note in AGENTS.md):

    cd <BASE>
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        /home/amh001/space_fs7/software_2022/python/py_2024/bin/python \
        _run_extreme_p95_nb.py

The notebook is regenerated from `<BASE>/_mk_extreme_p95_nb.py` (edit the
builder, re-run it, then run the notebook).  The daily series are cached in
`<BASE>/pickles/extreme_p95_daily.pkl`; bump `CACHE_VERSION` in the notebook
(rebuild via the builder) if the readers or station set change.  A cache-hit
re-run reproduces the result CSVs exactly.

## Outputs

| output | location |
|--------|----------|
| `extreme_p95_table.csv` (per-station p75/p95 + bias) | `extremes_validation/results/` (tracked) |
| `extreme_p95_summary.csv` (per-dataset mean percentile + bias) | `extremes_validation/results/` (tracked) |
| `p75_obs_vs_model.png`, `p95_obs_vs_model.png` (2 x 5 scatter) | `<BASE>/figures/extremes_validation/` |
| `p75_bias_hist.png`, `p95_bias_hist.png` (2 x 5 bias histograms) | `<BASE>/figures/extremes_validation/` |
| `p75_vs_p95_bias.png` (p75 vs p95 bias per station) | `<BASE>/figures/extremes_validation/` |
| `extreme_p95_skill.csv` (p95 bias / RMSE / R² per dataset) | `extremes_validation/results/` (tracked) |
| `p95_taylor.png` (2-panel p95 Taylor diagram) | `<BASE>/figures/extremes_validation/` |
| `extreme_p95.pkl`, `extreme_p95_daily.pkl` (caches) | `<BASE>/pickles/` |

## Headline results (2026-10-01)

Mean of the per-station percentile across the 110 stations (bias = model -
obs): obs p75 means T 14.63 C / P 3.85 mm/day; obs p95 means T 19.67 C /
P 18.94 mm/day.

| dataset | T p75 bias (C) | T p95 bias (C) | P p75 bias (mm/d) | P p95 bias (mm/d) |
|---------|----------------|----------------|-------------------|-------------------|
| D3 (WRF d03)  | -0.52 | **+0.20** | -1.21 | **-1.76** |
| D15 (WRF d02) | -1.45 | -0.73 | -0.24 | +18.73 |
| D75 (WRF d01) | -2.49 | -1.69 | +0.41 | +24.69 |
| CanRCM4       | -2.95 | -2.02 | +0.97 | +0.79 |
| CanESM2       | -0.89 | +1.58 | +0.76 | -1.60 |

Reading the pair of percentiles together:

* **T:** the p75 bias tracks each model's whole-distribution (mean-like)
  bias.  CanESM2's p95 bias (+1.58 C) is *worse* than its p75 bias
  (-0.89 C) — it over-heats the warm tail specifically; its p95 was never
  actually "good", just competitive in absolute value.  WRF d03 is the only
  dataset whose p95 bias (+0.2 C) is smaller than its p75 bias (-0.5 C): it
  reproduces the warm tail.
* **P:** D3 and CanRCM4 stay within ~2 mm/day at *both* percentiles (good
  across the upper half of the distribution).  CanESM2 flips sign between
  p75 (+0.8) and p95 (-1.6): its heavy-rain tail is too weak relative to its
  median rain.  The coarser WRF domains are fine at p75 (|bias| < 1 mm/day)
  but grossly over-state the p95 (+19 / +25 mm/day): their extreme daily
  rain events are far too intense.  The D3 P p95 bias (median per-station
  -8.9 %) is consistent with the canonical Rule-1 annual-mean P bias on the
  same subset (-7.8 %).

## p95 skill (bias / RMSE / R²) — Section 11, Figure 6

RMSE over the 110 stations; R² = pooled Pearson R² of the annual p95 series
(station x year, shared surviving years).  Full table in
`results/extreme_p95_skill.csv`.

| var | dataset | bias | RMSE | R² |
|-----|---------|------|------|------|
| T (C) | D3 | +0.20 | 2.23 | 0.69 |
| T (C) | D15 | -0.73 | 2.19 | 0.60 |
| T (C) | D75 | -1.69 | 2.57 | 0.50 |
| T (C) | CanRCM4 | -2.02 | 3.02 | 0.48 |
| T (C) | CanESM2 | +1.58 | 2.75 | 0.33 |
| P (mm/d) | D3 | -1.76 | 4.59 | 0.73 |
| P (mm/d) | D15 | +18.73 | 24.87 | 0.68 |
| P (mm/d) | D75 | +24.69 | 31.04 | 0.56 |
| P (mm/d) | CanRCM4 | +0.79 | 7.23 | 0.54 |
| P (mm/d) | CanESM2 | -1.60 | 10.38 | 0.21 |

D3 is closest to the reference in both panels (best R², std near obs);
D15/D75 are over-dispersed in P (std ~2.2x obs) with a large wet bias but
still reasonable R² — the coarse domains track *which* stations/years are
wettest but over-state the intensity of the extremes; CanESM2 is the most
under-dispersed in P (std ~0.37x obs, R² 0.21) and the most mis-correlated
in T (R² 0.33), i.e. its warm T tail and weak P tail are genuinely different
extreme-day distributions, not small shifts.
