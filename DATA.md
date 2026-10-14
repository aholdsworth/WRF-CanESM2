# Data archive manifest

The analysis working tree (the repo's **parent** directory,
`AI_WORKDIR/NOTEBOOKS/`) holds all data; nothing data-like is committed to
this repo. This file documents the archive contents so a checkout plus the
archive reproduces every result. The authoritative on-disk inventory
(regenerated 2026-09-21, regenerate with `make_inventory.py` in the working
tree) is `DATA_INVENTORY.md` in the working tree root.

## Sizes (2026-10-14)

| path (working tree) | size | contents |
|---|---:|---|
| `data/` | 22 G | all raw model + obs data (see below) |
| `pickles/` | 774 M | derived caches / canonical tables (regenerable from `data/` via the scripts in this repo) |
| `plotting_files/` | 232 M | static geo inputs: `geo_em.d0{1,2,3}.nc`, bathymetry, orography, `NEP36-original.nc` |
| `files/` | 3.7 M | station/buoy metadata CSVs (station lists per agency + variable) |
| `report/` | 1.2 M | `verification_report.tex` + table `.tex` outputs |
| `Eva_data/` | 797 M | Eva's original comparison inputs (reference only) |
| `weight_files_gen/` | 366 M | locally regenerated CDO bicubic remap weights (Can* ECCC/BCH; used by `station_geo._can_weight_file` as the non-PFM source) |
| `station_cdo/` | 31 M | CDO extraction scripts + buoy descs + `eccc_bch_desc_map.json` + extraction logs (the desc map and buoy descs are also tracked in the repo under `Evaluation/`) |
| `CanESM2_raw-20260922T164025Z-1-001.zip` | 241 M | original source archive (see below) |

## Model data (1986–2005, historical)

| dir | size | what |
|---|---:|---|
| `data/wrf_stations/` | 3.1 G | WRF d01/d02/d03 per-station hourly extractions (CDO remap+genbic), `{t,pr,wind}_{dom}_*.nc` |
| `data/wrf_canon/` | 3.2 G | **canonical** WRF per-station extractions, `{var}_{dom}_{AGENCY}_{ID}.nc` — the single WRF model source for all station comparisons (built by `DataPrep/make_wrf_canon.py`). **Note:** the 46xxx buoy files here are corrupt and must not be re-parsed (see `AGENTS.md` in the working tree). |
| `data/cdo_extractions/` | 2.2 G | CDO per-station/buoy M10 extractions: `canesm2_raw_noaa113/`, `canrcm4_noaa113/` (land), `canesm2_raw_buoys/`, `canrcm4_buoys/`, `wrf_d02d03_cdo/`, `*_weights/` |
| `data/canesm2_raw/` | 5.3 G | CanESM2 raw IBCS (CMIP5 `CanESM2.r1i1p1.historical`) gridded daily, T63 64×128: `tas_hist.nc` [K], `pr_hist.nc` [kg m⁻² s⁻¹], `wind_hist.nc` [m/s] |
| `data/canrcm4/` | 6.6 G | CanRCM4 gridded daily, **rotated-pole** NAM22 grid: `tas_NAM22_hist.nc`, `pr_NAM22_hist.nc`, `wind_NAM22_hist.nc` (2-D lat/lon — do not nearest-point in lat/lon) |
| `data/wrf_pfm_cache/` | 20 M | PFM WRF output cache fragments used by extraction scripts |

The "CanESM2" token refers to three distinct datasets — raw IBCS gridded,
WRF-downscaled d01–d03, and the CF-converted d03 — do not conflate them
(`DATA_INVENTORY.md` §0).

## Observations (1986–2005)

| dir | size | what |
|---|---:|---|
| `data/noaa_obs/` | 122 M | NOAA GHCN-d daily CSVs (5 5-year files). Units: PRCP **mm** (not inches), TAVG tenths °C, SNOW water equiv, AWND km/h. |
| `data/eccc_obs/` | 116 M | ECCC daily obs (`en_climate_daily` CSVs, 44 station dirs × 20 yrs): 'Mean Temp (C)', 'Total Rain (mm)'; no wind |
| `data/eccc_hourly_obs/` | 664 M | ECCC hourly obs (hourly wind + tas for the 17 ECCC land-wind stations) |
| `data/bch_obs/` | 4.6 M | BC Climate (BCH) 24-h precip (`BCH_p24Id.csv`, precip in **mm**) + daily T |
| `data/buoy_obs/` | 173 M | buoy obs: tas daily (ECCC CSV) + wind 1-h (NOAA `.mat` + ECCC CSV), 15 files |

## Source archives

- `CanESM2_raw-20260922T164025Z-1-001.zip` (241 M) — the original raw-data
  transfer zip. Checksums (verified 2026-10-14):
  - MD5: `d07bfccf71042013e8d01d05598b5772`
  - SHA-256: `d07c3dcfca6c304400ac6d9467b9f10ad29a9abd86471bbb3d4d912c6c15f450`

## Derived caches (`pickles/`, all regenerable)

Dependency order (see `README.md` "How to run"): `station_tables.pkl`
(canonical 173-row land table) → `skill_seasonal.pkl` / `skill_summary.pkl`
→ `bias_figure.pkl`, `taylor_annuals.pkl`, buoy table pickles → figures.
Older superseded pickles are preserved under `pickles/OLDWAYS/` and `*.bak*`
files (not part of the archive).

## Archiving / transfer

The 22 G `data/` tree is transferred out-of-band (no network egress in the
development sandbox; GitHub LFS is not configured for it). Before transfer:

```bash
cd AI_WORKDIR/NOTEBOOKS
find data pickles files plotting_files report weight_files_gen station_cdo -type f -print0 \
  | xargs -0 sha256sum > data_archive.sha256   # regenerate before each transfer
tar -cJf wrf-canesm2-data-20261014.tar.xz data pickles files plotting_files report weight_files_gen station_cdo
```

Verify on the receiving side with `sha256sum -c data_archive.sha256`
(checksum file itself excluded from the tarball by generating it first).

## Licensing / provenance

- WRF run + CanESM2/CanRCM4 gridded fields: produced by the authors' group
  on the HPC system (PFM/evg000 mounts); raw sources documented in
  `data/MANIFEST.txt`.
- NOAA GHCN-d: public domain (U.S. federal government data).
- ECCC obs: Environment and Climate Change Canada (open government data).
- BCH (BC Climate): Government of British Columbia.
- Buoy obs: NOAA / ECCC as per-station in `files/NOAA_buoys.csv`,
  `files/ECCC_buoys.csv`.
