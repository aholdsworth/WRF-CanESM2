# Data preparation

Scripts for turning the raw CanESM2-WRF d03 daily output into the
CF-compliant files used in the publication.

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
