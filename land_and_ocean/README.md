# Land and ocean regional averages

Notebooks that aggregate the WRF d03 analysis results over land and ocean
(masks from `domain/geo_em.d03.nc` `LANDMASK`) and plot the regional time
series / changes.

## Notebooks

| notebook | consumes |
|----------|----------|
| `Regional_Averages.ipynb` | WRF d03 daily data directly (self-contained land/ocean averages) |
| `Regional_Averages_seasonal.ipynb` | `DATA/WRF/MEANS/{t,pr,wind}_seasonalAmp_landocean_avg_{rcp45,rcp85}.csv` — produced by `../mean_changes/compute_seas_amplitudes_chen.py` |
| `Regional_Averages_sc_chen.ipynb` | `DATA/WRF/MEANS/{var}_seasonalAmp_change_{period}.nc` — produced by `../mean_changes/compute_seasonal_changes_cycles.py` |
| `Regional_averages_extremes.ipynb` | the bootstrap pickles from `../extremes/` (`*_boots_years_1000_perc_*.pkl` and `*_region_null_boots_1000_perc_*.pkl`) |

The producing scripts live in `../mean_changes/` and `../extremes/`; see
those folders' READMEs for usage.
