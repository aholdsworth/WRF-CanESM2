# Compound extremes (temperature x precipitation)

Two-stage pipeline for compound-extreme analysis of WRF d03
(CanESM2-WRF): build day-of-year percentile thresholds, then count
compound events (warm/cold x wet/dry) against the historical thresholds.

## Stage 1 — thresholds: `compute_percentiles_rolling.py`

Computes per-day-of-year percentiles with a small **rolling window**
(5 days for T2, 29 days for pr) to smooth the threshold curve:

- T2: 10th and 90th percentiles (cold / warm thresholds)
- pr: 25th and 75th percentiles (dry / wet thresholds)

Run for each period (currently only `rcp85` is enabled in the script's
`periods` list; set it to `["hist", "rcp45", "rcp85"]` to regenerate all).
Outputs:

    DATA/WRF/percentiles/{tas,pr}{perc}p_{period}_roll{5,29}.nc

Only the **historical** thresholds are used for counting (the script
hard-codes `base = "hist"`), so at minimum the `hist` files are required.

## Stage 2 — counting: `count_compound_events_lte.py`

For each season (MAM, JJA, SON, DJF), counts days where both conditions
hold, against the historical daily-varying thresholds:

- `warm`/`cold`: T2 above 90th / below 10th historical threshold
- `wet`/`dry`: pr above 75th / at-or-below 25th historical threshold

DJF uses the meteorological season-year (Dec counted in the following
year). Outputs (one file per type x period x season):

    DATA/WRF/CompoundCounts/CanESM2-WRF-{warm|cold}_{wet|dry}_{period}_{season}_yearly_lte.nc

Usage:

    python count_compound_events_lte.py {hist|rcp45|rcp85} {warm|cold} {wet|dry}

Note: the script requires `import datetime` (used by the fallback time
decoding in `decode_time`); this has been added.

## `count_drizzle_days.py`

Counts drizzle days (pr < 0.2 mm/day) per season for historical and
rcp45, with dask chunking:

    DATA/WRF/percentiles/drizzle_days_{hist|rcp45}_0.2.nc

## Notebooks

- `plotCompoundExtremes.ipynb` — loads the Stage-2 `CanESM2-WRF-*_yearly_lte.nc`
  counts (all four event types, all periods/seasons) and maps them.
- `percentiles.ipynb` — exploration of the Stage-1 percentile files.
