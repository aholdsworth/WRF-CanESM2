# Buoy validation (10 m wind speed + near-surface air temperature)

Buoy validation of two near-surface quantities — **10 m wind speed** and
**near-surface air temperature (`tas`)** — for the CanESM2-WRF
nested-downscaling study (1986–2005), across five model columns: WRF D01
(outer nest, 75 km), D02 (middle nest, 15 km) and D03 (inner nest, 3 km)
downscalings of CanESM2, raw CanESM2 (the parent GCM, ~100 km — the coarsest
column), and CanRCM4 (~25 km, a comparable downscaling of CanESM2). This is
the buoy analogue of the land-station analysis (see `DATA_INVENTORY.md` M12
in the working tree). The wind and temperature pipelines share the same
conventions (90 % missing-data cascade, apples-to-apples shared-year
pairing, metrics on shared-year annuals) so the two comparisons are directly
comparable.

**Combined publication table (LaTeX):**
[`results/buoy_summary_table.tex`](results/buoy_summary_table.tex).

## Deliverables

| quantity | report text + tables | mathematical methods | data |
|---|---|---|---|
| wind | [`results/buoy_wind_results_for_report.md`](results/buoy_wind_results_for_report.md) | [`results/buoy_wind_methods_math.md`](results/buoy_wind_methods_math.md) | `results/buoy_wind_bias.csv`, `results/buoy_wind_bias_pooled.csv` |
| tas | [`results/buoy_tas_results_for_report.md`](results/buoy_tas_results_for_report.md) | [`results/buoy_tas_methods_math.md`](results/buoy_tas_methods_math.md) | `results/buoy_tas_bias.csv`, `results/buoy_tas_bias_pooled.csv` |

## Scripts

| script | purpose |
|--------|---------|
| `buoy_wind_bias.py` | Wind main analysis. Parses the NOAA `.mat` observations (vector magnitude, 5 m→10 m log-profile scaling, 90 % day/month/season/year cascade), loads the five per-buoy model series, writes per-buoy and pooled skill metrics. Self-contained. |
| `gen_buoy_wind_d01_weights.py` | Hand-written SCRIP 1-point **bilinear** remap weights for the D01 99×99 curvilinear grid. **Rev2 — fixes a transposed-grid bug** (see below). |
| `gen_buoy_wind_d03_weights.py` | Hand-written SCRIP 1-point bilinear weights for the D03 300×300 cell grid, from `geo_em.d03.nc`. |
| `extract_buoy_wind_d01.sh` | CDO `remap` extraction of the 7 ECCC D01 buoy wind files from the PFM `wind_d01_hourly.nc`. |
| `buoy_tas_bias.py` | Tas main analysis. Parses the buoy `.mat` "Temperature:Air" channel (deg C → K, ~15-min → hourly → 90 % cascade), loads the five per-buoy model series (WRF D01/D02/D03 hourly T2, CanESM2/CanRCM4 daily tas), writes per-buoy and pooled skill metrics. Self-contained (inlines the wind pipeline's cascade). |
| `extract_buoy_t_d01.sh` | Re-extract D01 T2 for all 14 buoys from PFM `t_d01_hourly.nc` (the original PFM pull had rotated ECCC filenames). Resumable (`--force`). Run from the NOTEBOOKS root. |
| `extract_buoy_t_parallel.sh` | Re-extract D02 + D03 T2 for the 14 buoys in parallel (`xargs -P 8`) from PFM `t_d02_hourly.nc` / `t_d03.nc`, with coordinate-verified descriptors + SCRIP weights (the rotated Sep-15 files were re-extracted). Resumable. Run from the NOTEBOOKS root. |
| `verify_buoy_tas_coords.py` | Check embedded coords vs `station_cdo/buoys/buoy_descs.csv` for all buoy tas files (14 D01 + 14 D02 + 14 D03 + 16 CanESM2 + 16 CanRCM4). Exit 0 = all match. **Run after any re-extraction.** |

Run order (wind): `gen_buoy_wind_d0{1,3}_weights.py` → `extract_buoy_wind_d01.sh`
(D01) plus the d02/d03/CanESM2/CanRCM4 per-buoy extraction → `buoy_wind_bias.py`.
Run order (tas): `extract_buoy_t_d01.sh` + `extract_buoy_t_parallel.sh` →
`verify_buoy_tas_coords.py` → `buoy_tas_bias.py`.

```
python buoy_wind_bias.py                 # relaxed year gate (>=3/4 seasons)
python buoy_wind_bias.py --year-gate 4 --tag _strict   # strict (>=4/4 seasons)
python buoy_tas_bias.py --year-gate 3
python verify_buoy_tas_coords.py         # after any re-extraction
```

### D01 transposed-grid bug (fixed)

The first D01 wind weight generator transposed the curvilinear grid and wrote
the declared source centres as `lat.T.ravel()`. CDO `remap` reads `wspd` in
**file C-order** (dim `lon` major), so every declared coordinate was the
*swapped neighbour's* physical location. Because the D01 grid is asymmetric
(max |lat − latᵀ| ≈ 5°), CDO interpolated at the wrong place (e.g. buoy 46134
landed near 48.9°, −138.2°, mean 8.45 m s⁻¹, instead of ≈3.4 m s⁻¹ near
48.7°, −124.4°). The correct layout — verified field-for-field against the
working D02/D03 weight files — is: `src_grid_center_lat/lon == the source
file's lat/lon in file C-order`, `src_grid_dims = [nlon, nlat]`. Rev2 works in
file order throughout. A second, independent bug (a `dst_grid_frac` created on
the `src_grid_size` dimension) caused a CDO segfault and is also fixed.
D02/D03 grids are near-symmetric, so the same class of bug never showed there.

## Model columns (2026-09-23 state)

| column | source | cadence | status |
|---|---|---|---|
| WRF D01 (75 km) | wind: `data/wrf_stations/wind_d01_ECCC_buoy_<ID>.nc`; tas: `data/wrf_stations/t_d01_{ECCC,NOAA}_buoy_<ID>.nc` (T2, K) | hourly | ✅ tas re-extracted 2026-09-23, coords verified |
| WRF D02 (15 km) | `data/wrf_stations/wind_d02_st<ID>.nc` / `t_d02_st<ID>.nc` | hourly | ✅ tas re-extracted 2026-09-23, coords verified |
| WRF D03 (3 km) | `data/wrf_stations/wind_d03_st<ID>.nc` / `t_d03_st<ID>.nc` | hourly | ✅ tas re-extracted 2026-09-23, coords verified (46005 excluded: off-domain nearest cell) |
| CanESM2 (~100 km) | `data/cdo_extractions/canesm2_raw_buoys/{sfcWind,tas}_buoy_<ID>.nc` (1850–2005 → 1986–2005) | daily | ✅ coords verified |
| CanRCM4 (~25 km) | `data/cdo_extractions/canrcm4_buoys/{sfcWind,tas}_buoy_<ID>.nc` (1986–2005) | daily | ✅ coords verified |

## Observation

`MB_<ID>_HM.mat` (NOAA/MSC-OS met buoys): wind = vector components
(magnitude, 5 m → 10 m log-profile scaling, z0 = 2 mm); temperature =
column 4 "Temperature:Air" (deg C → K, ~15-min cadence, NaN = missing);
column 0 = MATLAB datenum (classic epoch 719529); header coords = truth. The
sensor is in the buoy met-package above the sea surface (package wind
stresses referenced to 4.2 m); no height scaling is applied to temperature
(the log-profile correction is a momentum correction). Restricted to
1986-01-01..2005-12-31.

## Conventions (locked, identical for wind and tas)

- 90 % missing-data cascade (day/month/season); year = mean of passing
  seasons, kept if ≥ `--year-gate` (default 3) of 4 seasons; `--tag _strict`
  = 4/4.
- Apples-to-apples: model annuals computed on the obs surviving years.
- Metrics on shared-year annuals: period means, raw bias (m s⁻¹ / K), bias
  %, MAE, RMSE (n ≥ 2).
- Pooled headline = mean over the core set: buoys with ≥ `--min-years`
  (default 8) shared years.
  - **Wind: 6-buoy core** (46005, 46041, 46131, 46146, 46204, 46206); the D03
    pooled row covers 5 (46005 excluded — off the D03 nest domain).
  - **Tas: 9-buoy core** (46005, 46041, 46050, 46131, 46132, 46146, 46204,
    46206, 46207); the D03 pooled row covers 8 (46005 excluded — off-domain
    nearest cell).
- Resolution labels: D01 outer nest 75 km, D02 middle 15 km, D03 inner 3 km,
  CanESM2 parent GCM ~100 km (coarsest), CanRCM4 ~25 km.
- PFM ECCC-buoy filenames are NOT trusted — WRF files are matched by
  embedded coords (the scripts do this automatically).

## Inputs (not stored in this repo)

Set as environment variables (defaults are the HPC paths used for the
analysis):

| variable | default | what |
|----------|---------|------|
| `BUOY_OBS_ROOT` | `/gpfs/fs7/.../amh001/DATA/BuoyData` | `MB_<ID>_HM.mat` NOAA observations |
| `WIND_D01_HOURLY` | `/gpfs/fs7/.../CanESM2-WRF/historical/variables_complete/wind_d01_hourly.nc` | D01 wind source grid (99×99, dims `time, lon, lat`) |
| `GEO_EM_D03` | `../geo_em.d03.nc` (repo root) | D03 geometry |
| `CDO`, `PFM` | HPC paths | used by the extraction scripts |

`station_cdo/buoys/buoy_descs.csv` + `station_cdo/buoys/descs/*.txt` hold the
truth buoy coordinates and CDO single-point descriptors (committed here).

## Outputs

The bias scripts write (relative to this folder):

- `files/buoy_{wind,tas}_bias[_<tag>].csv` — per-buoy × model metrics
  (includes both raw `bias` and `bias_pct`)
- `files/buoy_{wind,tas}_bias_pooled[_<tag>].csv` — pooled core rows
- `pickles/buoy_wind_tables[_<tag>].pkl` — annual/seasonal tables (wind)

The canonical (relaxed) results are copied to [`results/`](results/).

## Buoys

The **10-buoy packaged set** is the 7 ECCC (46131, 46132, 46134, 46146,
46204, 46206, 46207) + 3 legacy NOAA (46005, 46029, 46041). The four "new
NOAA" buoys (46050, 46087, 46088, 46089) are **excluded from the wind
package**: their D01 wind files were not re-extracted (only D02/D03 files
exist). 46050 does enter the tas package (its D01/D02/CanESM2/CanRCM4 tas
files exist) but is off the D03 nest domain, so it is excluded from the D03
rows.

## Environment

Python: the `py_2024` env (3.9, `netCDF4`, `numpy`, `pandas`, `scipy`).
CDO 1.9.9 for the extraction scripts. No network access required.
