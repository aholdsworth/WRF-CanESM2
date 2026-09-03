# Mean changes (seasonal and annual)

Scripts and notebooks for the seasonal/annual mean-change analysis of WRF
d03 (CanESM2-WRF) relative to the historical period, with comparison fields
from CanRCM4, raw CanESM2, and ensemble data.

## Scripts

Each script computes seasonal (DJF/MAM/JJA/SON) and annual means of a
variable over the historical (1986-2005) and future periods, then saves the
deltas and t-test p-values as a compressed `.npz` file. The variable,
period, and domain are set in the CONFIGURATION block at the top of each
script (edit before running).

| script | data | output `.npz` |
|--------|------|---------------|
| `compute_seasonal_changes.py` | WRF d03 daily (T2, pr, wspd) | `seasonal_deltas_and_pvals_{t,pr,wind}_{period}_d03.npz` |
| `compute_seasonal_changes_CanRCM4.py` | CanRCM4 gridded daily (tas, pr, sfcWind) | `seasonal_deltas_and_pvals_{t,pr,wind}_{period}_CanRCM4.npz` |
| `compute_seasonal_changes_CanESM2.py` | raw CanESM2 gridded daily | `seasonal_deltas_and_pvals_{t,pr,wind}_{period}_CanESM2.npz` |

### Seasonal-cycle (Chen et al. 2019) scripts

| script | what it computes | output |
|--------|------------------|--------|
| `compute_seasonal_changes_cycles.py` | seasonal amplitude (peak-to-peak) change %, and circular day-of-max/min shifts vs historical, per grid cell | `DATA/WRF/MEANS/{var}_seasonalAmp_change_{period}.nc` — consumed by `../land_and_ocean/Regional_Averages_sc_chen.ipynb` |
| `compute_seas_amplitudes_chen.py` | same diagnostics, but area-weighted land/ocean averages (mask from `domain/geo_em.d03.nc`) | `DATA/WRF/MEANS/{var}_seasonalAmp_landocean_avg_{period}.csv` — consumed by `../land_and_ocean/Regional_Averages_seasonal.ipynb` |

Outputs go to the `out_file` path in each script (absolute GPFS path).
`t` = absolute change (K); `pr`/`wind` = relative change (%).

## Notebooks

- `plot_means_{temperature,precip,wind}(-RCP85).ipynb` — the main
  comparison figures. Load the `.npz` files above (WRF d03, CanRCM4,
  CanESM2) plus pre-computed **CORDEX** and **CMIP5** seasonal-change files
  from `.../ensemble_means/NA-CORDEX_ensemble/` and
  `.../ensemble_means/CMIP5_ensemble/`.
- `AnnualAverages.ipynb` — annual (rather than seasonal) mean changes.
- `Elevation_vs_change.ipynb` — bins monthly mean changes by elevation
  using WRF d03 data directly (self-contained; no external scripts).
- The Chen seasonal-cycle outputs are plotted in `../land_and_ocean/`
  (`Regional_Averages_sc_chen.ipynb`, `Regional_Averages_seasonal.ipynb`).

## CORDEX / CMIP5 panels (provenance note)

The CORDEX panel of the `plot_means_*` notebooks reads pre-computed files of
the form

    CORDEX_{tas,pr,wind}_change_{rcp45,rcp85}_{season}.nc
    CORDEX_{tas,pr,wind}_change_rcp45_{season}_robustness.nc
    CORDEX_{tas,pr,wind}_change_rcp45_{season}_consensus.nc

together with analogous CMIP5 ensemble files. These were **not** produced by
the scripts in this folder: they were generated in **R** using the
[loadeR](https://doi.org/10.1016/j.envsoft.2018.09.009) / `loadeR.2nc`
package (file history attribute shows loadeR v1.7.0). The R script used to
make them is not included in this repository; the resulting netCDF files are
referenced directly by absolute path from the notebooks.
