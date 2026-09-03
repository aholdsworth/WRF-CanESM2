# WRF-CanESM2

Analysis codes and notebooks for the manuscript
*"Projecting future changes over the coastal Pacific Northwest: climate and
extremes from a convection-permitting model."*

The analysis is based on WRF d03 (300 x 300 curvilinear grid over the
north-western North America coast, lon -133.8 to -116.6, lat 43.4 to 54.7)
driven by CanESM2, with a 1986-2005 historical period and 2046-2065 future
periods (RCP4.5 and RCP8.5). Results are compared against CanRCM4, raw
CanESM2, and NA-CORDEX/CMIP5 ensemble fields.

## Folders

| folder | contents |
|--------|----------|
| `DataPrep/` | `cf_convert.py` — converts the raw WRF d03 daily files (pr, T2, wspd) to CF-compliant, zlib-compressed netCDF4 used for publication (bit-identical data, metadata only). See its README. |
| `mean_changes/` | Seasonal/annual mean-change scripts (deltas + t-tests for WRF d03, CanRCM4, CanESM2) and the `plot_means_*` notebooks (incl. CORDEX/CMIP5 comparison panels), annual averages, and elevation-binned change. |
| `extremes/` | Block-bootstrap scripts for extreme-percentile significance (per-gridcell and pooled land/ocean regional null test) and the `plot_extremes*` notebooks. |
| `CompoundEvents/` | Two-stage compound-extremes pipeline: per-DOY rolling percentile thresholds, then counts of warm/cold x wet/dry compound event days; drizzle-day counts; plotting notebooks. |
| `compression/` | Compression/PPC verification notebooks and scripts for the WRF d03 run. |
| `land_and_ocean/` | Regional-average (land/ocean) analysis notebooks. |
| `WRFDomainLib.py` | Shared helpers for the WRF d03 domain (coordinates, masks, plotting). |

## Notes

- Scripts in `DataPrep/`, `mean_changes/`, `extremes/`, and `CompoundEvents/`
  contain **absolute paths** to the original data on the HPC system where the
  analysis was run; edit the path constants at the top of each script before
  re-running.
- Data files are not stored in this repository; each folder's README
  documents the expected inputs/outputs and the notebook-to-script
  relationships.
