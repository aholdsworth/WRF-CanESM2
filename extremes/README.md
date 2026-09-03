# Extremes (bootstrap significance)

Bootstrap scripts for the extreme-percentile analysis of WRF d03
(CanESM2-WRF) against the historical period, and the plotting notebooks.

## Scripts

Both scripts work in **seasonal-year blocks** (DJF uses Dec of the previous
calendar year; blocks are the 3 months of each season-year) and run
1000 bootstrap iterations in parallel (joblib, 8 processes). Configuration
(variable, scenario) is set near the bottom of each script.

### bootstrap_gridcell_extreme_significance.py

Per-gridcell block bootstrap: resamples the historical season-year blocks,
computes the chosen percentile of each bootstrap sample, and derives the
distribution of the future-minus-historical percentile change.

Usage:

    python bootstrap_gridcell_extreme_significance.py {tmin|t|pr|wind} {yes|no} {yes|no}

- arg 1: variable (`tmin` = 5th percentile of T2; `t` = 95th percentile of
  T2; `pr`, `wind` as-is)
- arg 2: `yes`/`no` — express the change relative to the historical value
- arg 3: `yes`/`no` — subtract the historical median first

Scenario is `rcp85`, domain `d03`. Output (one pickle per season):

    DATA/WRF/Extremes/{var}_{season}_rcp85_boots_years_1000_perc_{perc}{rel}{med}.pkl

`rel` = `_rel` if relative; `med` = `_med50`/`_med75` if median-subtracted.

### bootstrap_region_mean_null_test.py

Pooled land/ocean regional null test: same block bootstrap, but pools grid
cells into land and ocean regions (mask from `domain/geo_em.d03.nc`
`LANDMASK`), with top-quartile filtering for precipitation. Used to test
whether the regional mean extreme change is significant.

Usage:

    python bootstrap_region_mean_null_test.py {tmin|t|pr|wind} {yes|no} {yes|no}

Output (one pickle per season):

    DATA/WRF/Extremes/{var}_{season}_{period}_region_null_boots_1000_perc_{perc}{rel}{med}.pkl

## Notebooks

- `plot_extremes.ipynb` (rcp45) and `plot_extremes-rcp85.ipynb` — load the
  per-season pickles from the gridcell script and map the percentile
  changes / significance.
- `plot_minusmed_pr_wind.ipynb` — plots of the median-subtracted variants.
