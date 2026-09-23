# Elevation-bias diagnostic grid — methods & equations

Plain-text + LaTeX description of the methods behind
`../elev_bias_grid.ipynb` and the table/figure producers. Sign convention
throughout: **bias = model − observation** (positive = model high).

## 1. Datasets (rows)

| row label | dataset | model grid | orography source |
|-----------|---------|-----------|------------------|
| D3        | WRF d03 | 3 km, 300×300 curvilinear | `geo_em.d03.nc` `HGT_M` |
| D15       | WRF d02 | 15 km | `geo_em.d02.nc` `HGT_M` |
| D75       | WRF d01 | ~75 km, 99×99 | `geo_em.d01.nc` `HGT_M` |
| CanRCM4   | CanRCM4 | ~20–28 km rotated pole (NAM-22) | `orog_CanRCM4.nc` |
| CanESM2   | CanESM2 (raw) | ~100 km lat/lon | `orog_CanESM2.nc` |

## 2. Elevation bias (y-axis of every scatter)

$$\text{elev\_bias} = \text{model\_elev} - \text{obs\_elev}$$

`obs_elev` is the station elevation from `files/*.csv`. **Model orography is
made extraction-consistent**: each model series was extracted by `cdo remap`
with a specific bicubic weight file (recorded in the T file's `history`
attribute), so the model elevation is the exact same CDO remap applied to
orography:

$$\text{model\_elev} = \sum_k w_k \cdot \text{orog}[i_k, j_k],
\qquad (i_k, j_k) = \mathrm{divmod}(\text{src\_address}_k, n_x)$$

- $w_k$ = the **station column** (`remap_matrix[:, 0]`) of the weight file;
  weights sum to 1.0 for on-land stations. Negative bicubic weights are
  legitimate (CDO `genbic`, `normalization = none`) — no renormalisation.
- `orog` is the *same file the model data was remapped from* (`HGT_M` for WRF,
  `orog_CanESM2/CanRCM4.nc` for the coarse rows). Ocean/masked source cells
  drop out of the sum.

**Fallback** (logged per station, recorded in `elev_method`): no T file / no
`remap` in `history` / weight file missing → nearest **cell-centre** lookup
(WRF, `elev_method = corner`) or nearest-grid-point (Can*, `elev_method =
nearest`). Coverage: all 175 WRF stations bicubic except a handful of buoys /
`BCK`/`BLN`; Can* 164/175 land bicubic, 11 nearest-GP (9 ECCC + 2 BCH).

## 3. Temperature bias (°C)

1986–2005 annual means of T2 (2 m) then difference:

$$T\text{-bias} = \overline{T}_{model} - \overline{T}_{obs}$$

## 4. Precipitation bias (%)

$$\text{PR-bias}(\%) = 100 \cdot \frac{\overline{P}_{model} - \overline{P}_{obs}}{\overline{P}_{obs}}$$

on annual means (mm/day). **Unit note:** GHCN-d (NOAA) `PRCP` is in **inches**
and is multiplied by 25.4 → mm (a 2026-09-23 bug had wrongly removed this; see
handoff). ECCC obs are native mm; BCH p24 tenths-of-mm for {BLN, DLU, MIS,
NTY, WOL} are ÷10.

**Absolute PR bias (mm/day)** (plotted in §5b/§5c): recovered exactly from the
percent table by

$$P\text{-bias}_{abs} = \frac{\text{PR-bias}(\%)}{100} \cdot \overline{P}_{obs\,[mm/day]}$$

with $\overline{P}_{obs}$ the 1986–2005 annual-mean precip in mm/day
(annual mm/yr ÷ 365). No pickle rebuild is needed for this re-expression.

## 5. Wind-speed bias (%)

$$\text{Wind-bias}(\%) = 100 \cdot \frac{\overline{U}_{model} - \overline{U}_{obs}}{\overline{U}_{obs}}$$

on annual means of hourly wind speed (m/s). ECCC obs = hourly `Wind Spd
(km/h)` ×0.277778 → m/s, mean over all hourly values 1986–2005.

## 6. Observations-side gate (seasonal pooling)

Season-first pooling, equal weight at each level, day → season → year → period
(Jones & McAvoy 2010; WMO No. 100):

- a **season** value is NaN if < 90 % of its expected days/hours are present
  (expected computed from the actual calendar, leap years included);
- a **year** is NaN if any of its 4 seasons is NaN;
- the **period** value is the mean of surviving annuals, reported only if
  ≥ 15 of 20 years survive.

Like-for-like: the surviving-year set comes from the OBS series (the only one
with gaps); the (complete) model series is restricted to those years before
averaging. This gate lives in `seasonal_pooling.py` (shared by
`make_skill_tables.py` and `mae_vs_eva.py`).

## 7. Skill metrics (per station, over shared surviving years)

$$\text{bias} = \overline{M} - \overline{O}, \qquad
\text{MAE} = \overline{|M_y - O_y|}, \qquad
\text{RMSE} = \sqrt{\overline{(M_y - O_y)^2}}$$

`M_y`, `O_y` are the yearly model/obs means. NaN if n_years < 15.

## 8. The grid figure (notebook §5)

5 rows × 4 columns. Column 1 = each row's **own** model topography (shared
terrain colourbar 0–3000 m, same map extent) with station points coloured by
that row's elevation bias (PRGn, −1500…1500 m, `TwoSlopeNorm` at 0).
Columns 2–4 = scatter of T / PR / wind bias (x) vs elevation bias (y,
−1500…1500 m). Each scatter: OLS fit (`scipy.stats.linregress`), R² of the
same variable pair, origin crosshairs, panel tags (a)–(t). The map column uses
manual figure-level axes (double gap before the T column). No cartopy
coastlines (offline — gridlines only).

## 9. Known artefacts & caveats (documented, not bugs)

- **Can* coarse-grid over-elevation:** the coarse rotated-pole / GCM grid
  over-estimates valley-station orography → the high-elevation band in the
  scatter is real physics, now shown at full magnitude by the bicubic
  orography (the older "vertical line" stacking artefact from pure
  nearest-GP is gone).
- **D75 PR outliers:** a few dry high-elevation ECCC sites (+2,000–7,000 %)
  dominate D01's PR mean/std and stretch the shared PR axis; the median
  (≈ +145 %) is representative.
- **Can* PR −100 % floor:** raw GCM-scale precip vs observed at wet NOAA
  stations floors the ratio near −100 % (a real large dry bias), confirmed by
  the absolute-PR table (≈ −28 mm/day).

Full provenance, bug-fix log and verification gates: the working-tree
`elevation_bias_handoff.md` (not in this repo).
