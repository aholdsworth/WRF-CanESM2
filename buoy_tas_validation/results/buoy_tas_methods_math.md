# Buoy `tas` validation — mathematical description of methods

Companion to `buoy_tas_results_for_report.md`. Mirrors the wind analysis
(`../buoy_wind_validation/results/buoy_wind_methods_math.md`) so the two
comparisons are directly comparable; differences are noted.

## 1. Observed near-surface air temperature

Each buoy `MB_<ID>_HM.mat` file stores, per ~15-min sample, a time column
(0), wind components (1–2), pressure (3), **air temperature (4, deg C)**
and further channels. Column 4 ("Temperature:Air") is used. Let
$t_i$ be the sample time and $x_i$ its value, with MATLAB datenum times
decoded against the classic epoch (datenum of 1970-01-01 = 719529):

$$ t_i = \mathrm{datenum}_i - 719529 \quad (\text{days since 1970-01-01}), $$

and samples are snapped to the hour they represent (the 15-min datenums
land ~3.3 µs off the exact hour). Missing values are NaN (the wind channel
uses exact-zero as missing; the temperature channel uses NaN).

**Units.** The observation is recorded in degrees Celsius; the model fields
(`T2`, `tas`) are in kelvin:

$$ x_i^{(K)} = x_i^{(C)} + 273.15 . $$

**Height.** The buoy air-temperature probe is in the met-package above the
sea surface (package wind stresses are referenced to 4.2 m); model values
are 2 m (WRF `T2`) or near-surface (CanESM2/CanRCM4 `tas`). No scaling is
applied: the logarithmic wind-profile correction used in the wind analysis
is a momentum correction and is not applicable to a conserved scalar.

**Aggregation.** Hourly value = mean of that hour's raw samples (typically
a pair at :00 and :59.999997):

$$ \bar{x}_h = \frac{1}{n_h}\sum_{i \in h} x_i^{(K)}, $$

daily value = mean of the 24 hourly values, and so on (Section 3).

## 2. Model series

Per-buoy 1×1-point series at the buoy's true coordinates (verified to
0.0000° against the `.mat`-header truth, `station_cdo/buoys/buoy_descs.csv`):

- **WRF D01 (75 km):** hourly `T2` (K), re-extracted 2026-09-23 from
  `t_d01_hourly.nc` (99×99) with 1-point nearest-cell SCRIP weights; the
  original PFM pull had rotated ECCC filenames and was replaced
  (`extract_buoy_t_d01.sh`).
- **CanESM2 (~100 km):** daily `tas` (K), 1850–2005, sliced to 1986–2005.
- **CanRCM4 (~25 km):** daily `tas` (K), 1986–2005.

D01 placement caveat: at ~75 km the buoy can lie up to ~0.35° from the
nearest cell centre, so the D01 series carries extra sub-grid placement
uncertainty.

## 3. Missing-data cascade (90 % gate at every level)

Identical to the wind analysis. Given a series at cadence $c$
(hours/day = 24 for obs and D01, 1 for the daily model columns):

- **Day.** Let $n_d$ be the number of present values in day $d$ out of
  $c_d$ expected. The day survives iff

$$ \frac{n_d}{c_d} \ge 0.9 \quad\text{and}\quad \max\ \text{intra-day gap} \le 3\ \text{h} \quad (c_d \ge 2), $$

  where the daily value is $\bar{x}_d = \frac{1}{n_d}\sum x$.

- **Month.** Month $m$ survives iff

$$ \frac{n_m}{c_m} \ge 0.9, $$

  with $n_m$ surviving days, $c_m$ calendar days; monthly value = mean of
  the surviving daily values.

- **Season.** Season $s$ survives iff at least 2 of its 3 complete months
  survive (a 1-of-3 season is 33 % complete < 90 %); seasonal value =
  mean of the surviving monthly values.

- **Year.** Year $y$ with passing seasons $S_y \subset \{\mathrm{DJF},
  \mathrm{MAM}, \mathrm{JJA}, \mathrm{SON}\}$:

$$ \bar{x}_y = \frac{1}{|S_y|}\sum_{s \in S_y} \bar{x}_s,
\qquad \text{kept iff } |S_y| \ge g, $$

  with $g = 3$ (relaxed default; $g = 4$ = locked spec, `--tag _strict`).

## 4. Apples-to-apples pairing and metrics

For buoy $b$ and model $m$, let $Y_{bm}$ be the intersection of the
surviving-year sets of obs and model (the pairing is obs-driven: the same
years enter every model column for a given buoy). With $n_{bm} = |Y_{bm}|$
shared annual values $\{(\bar{x}^{\,o}_{y}, \bar{x}^{\,m}_{y})\}_{y \in
Y_{bm}}$:

$$ \text{period means:}\quad \bar{x}^{o} = \frac{1}{n}\sum_{y\in Y}\bar{x}^{\,o}_{y},
\qquad \bar{x}^{m} = \frac{1}{n}\sum_{y\in Y}\bar{x}^{\,m}_{y}, $$

$$ \text{raw bias:}\quad \Delta = \bar{x}^{m} - \bar{x}^{o} \quad [\text{K}], $$

$$ \text{bias \%:}\quad 100 \cdot \frac{\Delta}{\bar{x}^{o}}, $$

$$ \text{MAE} = \frac{1}{n}\sum_{y\in Y}\left|\bar{x}^{\,m}_{y} - \bar{x}^{\,o}_{y}\right|,
\qquad
\text{RMSE} = \left(\frac{1}{n}\sum_{y\in Y}\left(\bar{x}^{\,m}_{y} - \bar{x}^{\,o}_{y}\right)^2\right)^{1/2}, $$

computed for $n \ge 2$.

## 5. Pooled headline

The core set is the buoys with $n \ge 8$ shared surviving years (for
temperature this is 9 buoys: 46005, 46041, 46050, 46131, 46132, 46146,
46204, 46206, 46207 — temperature gaps are much rarer than wind-vector
gaps, so three more buoys qualify than for wind). The pooled row for model
$m$ is the mean over core buoys $b$:

$$ \bar{x}^{o}_{\text{pool}} = \frac{1}{B}\sum_{b}\bar{x}^{o}_{b},\qquad
\bar{x}^{m}_{\text{pool}} = \frac{1}{B}\sum_{b}\bar{x}^{m}_{b}, $$

with pooled bias $\bar{x}^{m}_{\text{pool}} - \bar{x}^{o}_{\text{pool}}$,
pooled bias % $100(\bar{x}^{m}_{\text{pool}}-\bar{x}^{o}_{\text{pool}})/\bar{x}^{o}_{\text{pool}}$,
and pooled MAE/RMSE = mean of the per-buoy MAE/RMSE.

## 6. Differences from the wind analysis

1. Quantity: near-surface air temperature (K) instead of 10 m wind speed.
2. No height scaling (momentum-only correction).
3. Missingness: NaN (no exact-zero convention).
4. Model columns: D01 + CanESM2 + CanRCM4 only (D02/D03 buoy `t`
   extractions pending).
5. Core set: 9 buoys (temperature) vs 6 (wind).
