# Data preparation

Scripts for turning the raw CanESM2-WRF / CanRCM4 / CanESM2 output into the
CF-compliant gridded files and the per-station / per-buoy CDO extractions
used by the `Evaluation/` pipelines.

## Station & buoy extraction (feeds `Evaluation/`)

| script | produces |
|--------|----------|
| `make_wrf_canon.py` | `data/wrf_canon/` — canonical WRF per-station extractions, one `{var}_{dom}_{AGENCY}_{ID}.nc` per station/domain (t/pr/wind, d01–d03); the single WRF model source for all station comparisons. |
| `extract_canrcm4_stations.py` | CanRCM4 per-station CDO files (M10) for the land + USW station set. |
| `extract_wrf_wind_951.py` | WRF wind re-extraction for station 951 (corrupted-cache recovery). |
| `make_climo_d03.py` | `pickles/climo_d03_hist.pkl` — gridded 1986–2005 d03 climatology (t/pr/wind) from the CF-compliant d03 daily files. |
| `make_buoy_wind_obs.py` | `pickles/buoy_wind_obs_annuals.pkl` — buoy wind obs annuals for the Taylor wind panel. |
| `extract_buoy_t_d02d03.sh`, `extract_buoy_wind_d02d03.sh` | per-buoy CDO extractions (tas/wspd) for the d02/d03 columns. |
| `extract_d02_eccc_fix.sh` | ECCC d02 re-extraction fix (missing-station recovery). |
| `extract_wrf_wind_d01_usw.sh` | WRF d01 wind extraction for the 5 USW WA/OR wind stations. |

The buoy d01 extraction scripts + weight generators live in
`Evaluation/buoy_validation/` (part of that pipeline).  All extraction
scripts read raw model output and write under the analysis working tree
(`data/`, `station_cdo/`, `pickles/`) — none of that data is in this repo.

## cf_convert.py

Converts WRF d03 daily files to CF-compliant, zlib-compressed netCDF4.
One script handles all variables and scenarios.

### Usage

    /home/amh001/space_fs7/software_2022/python/py_2024/bin/python \
        cf_convert.py --var {pr,t2,wspd} --scenario {historical,rcp45,rcp85} [--force]

- `--force` overwrites an existing output (default: refuse).
- All paths in the script are absolute (hardcoded to the original data
  location on the HPC system); edit the `ROOT` constant to point elsewhere.

### What it does (metadata only — data values are never changed)

1. Reads the source, asserts the expected structure (7305 x 300 x 300,
   `time_bnds` present).
2. Writes a temporary output:
   - dims renamed `lon -> y`, `lat -> x`; `time_bnds`/`bnds` dropped
   - CDO, CDI, `_CoordinateAxisType` attributes removed
   - data variable gets `standard_name`, `long_name`, `units`;
     `cell_methods` kept if present
   - `time` keeps its source `units`/`calendar`, gains `standard_name`/`axis`
   - `lat`/`lon` get `standard_name` + degree units
   - global `Conventions = CF-1.11`; source `history` preserved and appended
   - zlib compression, complevel 4
3. Verifies the output is **bit-identical** to the source (all variables,
   full arrays) before atomically replacing the output file.

### Source discovery

The source file is found by **content** (which data variable it holds) in:

| scenario   | source dir                          |
|------------|-------------------------------------|
| historical | `WRF_FILES/original/`               |
| rcp45      | `WRF_FILES/original/RCP45/`         |
| rcp85      | `WRF_FILES/original/RCP85/`         |

This is because the originals use inconsistent file names (e.g. RCP85 files
carry a `_rcp85` suffix, RCP45 files do not; the RCP85 wind file has no
`wspd` in its name).

### Variable metadata

| --var  | source var | standard_name          | long_name                        | units  |
|--------|------------|------------------------|----------------------------------|--------|
| pr     | pr         | precipitation_amount   | daily total precipitation        | mm/day |
| t2     | T2         | air_temperature        | daily mean 2-m air temperature   | K      |
| wspd   | wspd       | air_speed              | daily mean wind speed            | m s-1  |

Note: the source files carry **no units attribute** on the data variable;
the units above are declared by this script (confirmed with the data owner:
mm/day, K, m s-1).

### Output naming

    CanESM2-WRF_<var>_d03_<scenario>_daily_<start>_<end>_300x300.nc

The year range is read from the file's own time axis.
(historical: 1986-2005; rcp45/rcp85: 2046-2065; 7305 days each)
