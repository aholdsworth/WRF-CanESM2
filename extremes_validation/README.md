# Extremes validation — 95th percentile (T & P) at nearly-complete stations

Obs-vs-model comparison of the **95th percentile of daily temperature (T)**
and **daily precipitation (P)** at land stations with **nearly-complete
records**, so that obs and model are like-for-like.

## Station set

110 land stations (ECCC + BCH + NOAA) from the canonical
`<BASE>/pickles/station_tables.pkl` (`rule1-canonical-land-v1`) with Rule-1
surviving years **>= 19 of 20 for both T and P**.

## Method

- **Like-for-like:** the surviving-year set is defined by the obs series
  (Rule-1: a year is valid only if >= 90 % of its expected days are present);
  each model series (complete) is restricted to the same years.
- **Percentile:** 95th percentile of the daily values pooled over all days of
  the shared surviving years.  P includes dry days (the "pr as-is"
  convention of the `extremes/` bootstrap pipeline).
- **Bias:** model p95 minus obs p95 (positive = model high).

Model datasets: `D3` = WRF d03, `D15` = WRF d02, `D75` = WRF d01,
`CanRCM4`, `CanESM2` (CDO per-station files under
`data/cdo_extractions/{canrcm4_noaa113,canesm2_raw_noaa113}`).

## Notebooks

- `extreme_p95.ipynb` — self-contained walkthrough: station selection, daily
  series build (cached), per-station p95, and the two figures.  Imports the
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
| `extreme_p95_table.csv` (per-station p95 + bias) | `extremes_validation/results/` (tracked) |
| `extreme_p95_summary.csv` (per-dataset mean p95 + bias) | `extremes_validation/results/` (tracked) |
| `p95_obs_vs_model.png` (2 x 5 scatter, obs vs model p95) | `<BASE>/figures/extremes_validation/` |
| `p95_bias_hist.png` (2 x 5 histogram of per-station p95 bias) | `<BASE>/figures/extremes_validation/` |
| `extreme_p95.pkl`, `extreme_p95_daily.pkl` (caches) | `<BASE>/pickles/` |

## Headline results (2026-09-30)

Mean of the per-station p95 across the 110 stations (bias = model - obs):
obs p95 means T 19.67 C, P 18.94 mm/day.

| dataset | T p95 bias (C) | P p95 bias (mm/day) |
|---------|----------------|---------------------|
| D3 (WRF d03)  | +0.20 | -1.76 |
| D15 (WRF d02) | -0.73 | +18.73 |
| D75 (WRF d01) | -1.69 | +24.69 |
| CanRCM4       | -2.02 | +0.79 |
| CanESM2       | +1.58 | -1.60 |

Only WRF d03 and CanESM2 reproduce the extreme daily rain within a few
mm/day; the coarser WRF domains grossly over-state it.  For T, WRF d03 is
nearly exact (median per-station bias -0.13 C).  The D3 P p95 bias
(median per-station -8.9 %) is consistent with the canonical Rule-1
annual-mean P bias on the same subset (-7.8 %).
