# Buoy wind-speed validation

Validation of 10 m wind speed against North Pacific / BC-coast buoy
observations (1986–2005), across five model columns: WRF d01 (outer nest,
75 km), d02 (middle nest, 15 km) and d03 (inner nest, 3 km) downscalings of
CanESM2, raw CanESM2 (the parent GCM, ~100 km — the coarsest column), and
CanRCM4 (~25 km, a comparable downscaling of CanESM2). This is the buoy
analogue of the
land-station analysis (see `DATA_INVENTORY.md` M12 in the working tree).

The **mathematical description of the methods** (height scaling, the 90 %
missing-data cascade, bilinear extraction, skill metrics, pooling) is in
[`results/buoy_wind_methods_math.md`](results/buoy_wind_methods_math.md) —
each formula in plain text and LaTeX. The report paragraph and value tables
are in [`results/buoy_wind_results_for_report.md`](results/buoy_wind_results_for_report.md).

## Scripts

| script | purpose |
|--------|---------|
| `buoy_wind_bias.py` | Main analysis. Parses the NOAA `.mat` observations (vector magnitude, 5 m→10 m log-profile scaling, 90 % day/month/season/year cascade), loads the five per-buoy model series, and writes per-buoy and pooled skill metrics. Self-contained (no external project imports). |
| `gen_buoy_wind_d01_weights.py` | Hand-written SCRIP 1-point **bilinear** remap weights for the D01 99×99 curvilinear grid. **Rev2 — fixes a transposed-grid bug** in the first version (see below). |
| `gen_buoy_wind_d03_weights.py` | Hand-written SCRIP 1-point bilinear weights for the D03 300×300 cell grid, from `geo_em.d03.nc`. |
| `extract_buoy_wind_d01.sh` | CDO `remap` extraction of the 7 ECCC D01 buoy wind files from the PFM `wind_d01_hourly.nc`. |

Run order: `gen_buoy_wind_d0{1,3}_weights.py` → `extract_buoy_wind_d01.sh`
(D01) plus the d02/d03/CESM2/RCM4 per-buoy extraction → `buoy_wind_bias.py`.

```
python buoy_wind_bias.py                 # relaxed year gate (>=3/4 seasons)
python buoy_wind_bias.py --year-gate 4 --tag _strict   # strict (>=4/4 seasons)
```

### D01 transposed-grid bug (fixed)

The first D01 weight generator transposed the curvilinear grid and wrote the
declared source centres as `lat.T.ravel()`. CDO `remap` reads `wspd` in **file
C-order** (dim `lon` major), so every declared coordinate was the *swapped
neighbour's* physical location. Because the D01 grid is asymmetric
(max |lat − latᵀ| ≈ 5°), CDO interpolated at the wrong place (e.g. buoy 46134
landed near 48.9°, −138.2°, mean 8.45 m s⁻¹, instead of ≈3.4 m s⁻¹ near
48.7°, −124.4°). The correct layout — verified field-for-field against the
working D02/D03 weight files — is: `src_grid_center_lat/lon == the source
file's lat/lon in file C-order`, `src_grid_dims = [nlon, nlat]`. Rev2 works in
file order throughout. A second, independent bug (a `dst_grid_frac` created on
the `src_grid_size` dimension) caused a CDO segfault and is also fixed.
D02/D03 grids are near-symmetric, so the same class of bug never showed there.

## Inputs (not stored in this repo)

Set these as environment variables (the defaults are the HPC paths used for
the analysis):

| variable | default | what |
|----------|---------|------|
| `BUOY_OBS_ROOT` | `/gpfs/fs7/.../amh001/DATA/BuoyData` | `MB_<ID>_HM.mat` NOAA observations |
| `WIND_D01_HOURLY` | `/gpfs/fs7/.../CanESM2-WRF/historical/variables_complete/wind_d01_hourly.nc` | D01 source grid (99×99, dims `time, lon, lat`) |
| `GEO_EM_D03` | `../geo_em.d03.nc` (repo root) | D03 geometry |
| `CDO`, `PFM` | HPC paths | used by `extract_buoy_wind_d01.sh` |

Per-buoy model series read from (relative to this folder):

- `data/wrf_stations/wind_d01_ECCC_buoy_<ID>.nc` (D01, matched by **embedded
  coordinates**, not filename), `wind_d02_st<ID>.nc`, `wind_d03_st<ID>.nc`
- `data/cdo_extractions/canesm2_raw_buoys/sfcWind_buoy_<ID>.nc` (daily)
- `data/cdo_extractions/canrcm4_buoys/sfcWind_buoy_<ID>.nc` (daily)

`station_cdo/buoys/buoy_descs.csv` + `station_cdo/buoys/descs/*.txt` hold the
truth buoy coordinates and CDO single-point descriptors (committed here).

## Outputs

`buoy_wind_bias.py` writes (relative to this folder):

- `files/buoy_wind_bias[_<tag>].csv` — per-buoy × model metrics
- `files/buoy_wind_bias_pooled[_<tag>].csv` — pooled core row
- `pickles/buoy_wind_tables[_<tag>].pkl` — annual/seasonal tables

The canonical (relaxed) results are copied to [`results/`](results/):
`buoy_wind_bias.csv`, `buoy_wind_bias_pooled.csv`, plus the report and methods
markdown.

## Buoys

The **10-buoy packaged set** is the 7 ECCC (46131, 46132, 46134, 46146,
46204, 46206, 46207) + 3 legacy NOAA (46005, 46029, 46041). The four "new
NOAA" buoys (46050, 46087, 46088, 46089) are **excluded from this
package**: their D01 files were not re-extracted (only their D02/D03 files
exist). They do not enter the 6-buoy pooled core, so the headline table is
unaffected.

## Environment

Python: the `py_2024` env (3.9, `netCDF4`, `numpy`, `pandas`, `scipy`).
CDO 1.9.9 for the D01 extraction. No network access required.
