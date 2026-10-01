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

Model datasets: `D3` = WRF d03, `D15` = WRF d02, `D75` = WRF d01,
`CanRCM4`, `CanESM2` (CDO per-station files under
`data/cdo_extractions/{canrcm4_noaa113,canesm2_raw_noaa113}`).

## Notebooks

- `extreme_p95.ipynb` — self-contained walkthrough: station selection, daily
  series build (cached), per-station p75/p95, and five figures.  Imports the
  shared obs loaders / Rule-1 gate from `<BASE>/seasonal_pooling.py` and
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
