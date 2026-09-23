# Buoy wind-speed validation — mathematical description of methods

Companion to `buoy_wind_results_for_report.md`. Implemented in
`buoy_wind_bias.py` (bias/metrics/cascade) and
`gen_buoy_wind_d0{1,2,3}_weights.py` + CDO `remap` (spatial extraction).
Each subsection gives the plain-text (Unicode) form used here and an
equivalent LaTeX form for the paper.

**Notation.** Buoys indexed by `j` (14 stations); models by `m` ∈
{WRF D01, WRF D02, WRF D03, CanESM2, CanRCM4}. Surviving calendar years
indexed by `i` over the set `Y_j` (the **shared-year intersection** of the
observation and model record for buoy `j`). Quantities:

- `s^{obs}_{j,t}` — observed 10 m wind speed (m s⁻¹), hourly
- `s^{mod}_{j,t,m}` — model 10 m wind speed (m s⁻¹) at buoy `j`, column `m`
- `d_{j,d}`, `s^{m}_{j,y}` — daily / annual (year `y`) means
- `o_{j,i}`, `x_{j,i,m}` — annual means of obs / model `m` in year `i`

---

## 1. Observed wind speed: vector magnitude, missingness, height scaling

The anemometer gives the northward (`V`) and eastward (`U`) components.
Instantaneous speed is the vector magnitude; a component that is **exactly
zero** is treated as missing (instrument convention), as are NaN samples.

```
w_t = sqrt( U_t^2 + V_t^2 )        (U_t or V_t == 0 or NaN  =>  missing)
```
LaTeX: `$w_t = \sqrt{U_t^{2}+V_t^{2}}$, with $w_t$ undefined when
$U_t$ or $V_t$ is exactly zero or missing.`

The instruments log at 5 m; the model fields are 10 m. The observation is
raised to 10 m with the logarithmic wind profile (constant-stress, neutral
boundary layer, reference roughness `z0`):

```
s_t^{obs} = w_t · ln( h_m / z0 ) / ln( h_o / z0 ) ,
             h_o = 5 m,  h_m = 10 m,  z0 = 0.002 m
      factor  = ln(5)/ln(2.5) = 1.0684
```
LaTeX:
`$s_t^{\mathrm{obs}} = w_t\,\dfrac{\ln(h_m/z_0)}{\ln(h_o/z_0)}$,
$h_o=5\,\mathrm{m}$, $h_m=10\,\mathrm{m}$, $z_0=0.002\,\mathrm{m}$;
the factor $\ln(5)/\ln(2.5)=1.0684$ is applied to the observations only
(model values are already at 10 m).`

Hourly value = mean of the raw samples falling in that clock hour
(observation samples are ≈15-min cadence).

---

## 2. Missing-data cascade (90 % gate, observation-driven)

Both observation and model series are reduced to annual means by the same
four-level cascade. A level **survives** only if at least a fraction
`g = 0.90` of its parent level is complete (plus, for days, a gap rule).

**Hourly → daily.** For day `d`, with `N_d` present hours and a maximum
consecutive in-day missing run `R_d` (hours):

```
d_{j,d} = (1 / N_d) · Σ_{t ∈ d, present} s_{j,t}
day survives ⇔  N_d / 24 ≥ g   AND   R_d ≤ g_max = 3 h
```
LaTeX: `$d_{j,d} = \tfrac{1}{N_d}\sum_{t\in d} s_{j,t}$ over present hours;
day kept if $N_d/24 \ge g$ and $R_d \le 3\ \mathrm{h}$. (The intra-day gap
rule prevents a single burst of samples from masquerading as a daily mean.)`

**Daily → monthly.** Month `m` survives if it has ≥ 90 % of its calendar
days present:

```
month survives ⇔  (days present in m) / (calendar days in m) ≥ g
```

**Monthly → seasonal.** Season = mean of its 3 constituent monthly means;
it survives if ≥ 2 of its 3 months survived:

```
s^{sea}_{j,y,s} = mean( monthly means of the 3 months in season s of year y )
season survives ⇔  ≥ 2 of 3 months present
```

**Seasonal → annual.** The annual mean is the mean of the surviving seasons
of that year; the year survives if at least `k` of the 4 seasons survived
(`k = 3` in the headline/"relaxed" run; `k = 4` in the "strict" run):

```
a_{j,y} = mean( s^{sea}_{j,y,s} over surviving seasons s )
year y ∈ Y_j  ⇔  (surviving seasons in y) ≥ k ,   k = 3 (relaxed) / 4 (strict)
```

Because the gate is applied to the observation first and the surviving-year
set `Y_j` is then used for the model, **model annuals are computed on exactly
the years the observation survived** (apples-to-apples; the gate changes
which years enter, never the obs/model pairing).

---

## 3. Spatial extraction: bilinear remap to the buoy point

Each model field is interpolated from its native grid to the buoy location
`(φ_b, λ_b)` using hand-written SCRIP weights and CDO `remap`. On a
curvilinear grid the four surrounding **cell centres** `c1…c4`, with
coordinates `(Φ_k, Λ_k)` and values `w_k`, define a local quadrilateral.
The point's barycentric parameters `(u, v)` are solved in that local tangent
plane and clamped to `[0, 1]` (clamping gives a nearest-node fallback at the
domain edge / open ocean). The extracted value is the bilinear combination:

```
w_b = (1-u)(1-v) w_1 + u(1-v) w_2 + (1-u) v w_3 + u v w_4 ,
      w_1,w_2,w_3,w_4 = values at the 4 surrounding cell centres
```
LaTeX:
`$w_b = (1-u)(1-v)\,w_1 + u(1-v)\,w_2 + (1-u)v\,w_3 + uv\,w_4$,
with $(u,v)$ the solution of the point's coordinates in the local
cell-centre quadrilateral, $u,v \in [0,1]$.`

The weights sum to one, are non-negative, and depend only on the buoy
position and the source grid (time-invariant). The same bilinear convention
is used for all three WRF domains and for the ECCC/NOAA set, so the columns
are directly comparable. **Resolution caveat (D01):** the ~75 km D01 cells
are ≈0.7° across, so a buoy may lie up to ≈0.35° from its nearest D01 cell
centre and the value is interpolated between four cell centres; D02/D03
buoy positions fall well inside a single cell. This sub-grid placement
uncertainty is specific to the coarsest **WRF** column (CanESM2 at ~150 km
is the coarsest model overall).

---

## 4. Skill metrics

For buoy `j` and model `m`, over the `n = |Y_j|` shared surviving years,
with `o_{j,i}` and `x_{j,i,m}` the annual means:

```
bias %_{j,m} = 100 · ( (1/n) Σ_i x_{j,i,m}  −  (1/n) Σ_i o_{j,i} ) / ( (1/n) Σ_i o_{j,i} )
MAE_{j,m}    = (1/n) Σ_i | x_{j,i,m} − o_{j,i} |
RMSE_{j,m}   = sqrt( (1/n) Σ_i ( x_{j,i,m} − o_{j,i} )^2 )
```
LaTeX:
`$\mathrm{bias\%}_{j,m} = 100\,\dfrac{\bar x_{j,m}-\bar o_j}{\bar o_j}$,
$\mathrm{MAE}_{j,m} = \tfrac{1}{n}\sum_{i\in Y_j}\lvert x_{j,i,m}-o_{j,i}\rvert$,
$\mathrm{RMSE}_{j,m} = \sqrt{\tfrac{1}{n}\sum_{i\in Y_j}(x_{j,i,m}-o_{j,i})^{2}}$,
where $\bar x_{j,m}=\tfrac{1}{n}\sum_i x_{j,i,m}$ and
$\bar o_j=\tfrac{1}{n}\sum_i o_{j,i}$.`

Metrics are reported only when `n ≥ 2` shared years exist.

---

## 5. Pooling (headline row)

The **core set** is the subset of buoys with `n ≥ n_min = 8` shared surviving
years: `{46005, 46041, 46131, 46146, 46204, 46206}` (D03 core is 5 buoys,
since 46005 lies off the D03 nest domain). The pooled value for a column is
the **unweighted mean of the per-buoy statistics** over the core set `C_m`
for that column:

```
Bias%_m = (1/|C_m|) Σ_{j ∈ C_m} bias%_{j,m}   (equivalently 100·(X_m − O_m)/O_m
              with O_m, X_m the core means of obs / model)
MAE_m   = (1/|C_m|) Σ_{j ∈ C_m} MAE_{j,m}
RMSE_m  = (1/|C_m|) Σ_{j ∈ C_m} RMSE_{j,m}
```
LaTeX:
`$\widehat{\mathrm{Bias\%}}_m = \tfrac{1}{|C_m|}\sum_{j\in C_m}
\mathrm{bias\%}_{j,m}$,
$\widehat{\mathrm{MAE}}_m = \tfrac{1}{|C_m|}\sum_{j\in C_m}\mathrm{MAE}_{j,m}$,
$\widehat{\mathrm{RMSE}}_m = \tfrac{1}{|C_m|}\sum_{j\in C_m}\mathrm{RMSE}_{j,m}$.`
