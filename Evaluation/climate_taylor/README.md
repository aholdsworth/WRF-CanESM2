# climate_taylor — "climate" 3-panel Taylor diagram (T / PR / wind)

Publication figure: Taylor diagram of the **20-yr (1986-2005) station
climatology** for 5 models, in 3 panels (2 m temperature, annual
precipitation, 10 m wind speed).

This is the "climate-scale" version of the Taylor diagram: each station's
20 annual values (1986-2005) are first reduced to their 20-year mean
climatology, and the per-station climatologies are pooled across stations.
It complements the anomaly-free pooled (station x year) 3-panel figure in
the main verification workspace.

## Files

| file | role |
|------|------|
| `climate_all_taylor.py`  | standalone script that builds the figure (executed cells of the notebook) |
| `climate_all_taylor.ipynb` | same content as a Jupyter notebook |
| `taylor_figure.py`       | shared Taylor-diagram library (axisartist polar axes, auto-scaled RMS contours, model style constants) |

`climate_all_taylor.py` reuses from `taylor_figure.py`: `MODELS`, `LABELS`,
`COLORS`, `MARKERS`, `VARS`, `VAR_UNITS`, `FIG_W_MM` / `FIG_H`, `FONT` /
`FONT_S`, `TaylorDiagram`, `_round_axis_ceiling`, `_auto_levels`.

## Input (not in this repo)

The script reads one pickle of per-station x per-year annual arrays:

    pickles/taylor_annuals.pkl        (cache version 10, ~0.5 MB)

- `ta['t']`  : (obs, models) — obs 2 m annual MEAN, models same
- `ta['pr']` : (obs, models) — annual SUM, mm
- `ta['wind']`: (obs, models) — 10 m annual MEAN, m/s
- `models` keys: `D03, D02, D01, CanESM2, CanRCM4`

Station sets (v10, all 5 model columns share the same obs set —
land stations from the canonical 173-row land table, plus the filtered
D03 ocean buoys): **t 154** (146 land + 8 buoys) / **pr 159** (land only) /
**wind 23** (16 land + 7 buoys).

The pickle is built by `taylor_data.py` in the main verification workspace
(`NOTEBOOKS/`, not in this repo) from the raw WRF post-process
(`data/wrf_canon`), CanESM2/CanRCM4 CDO extractions
(`data/cdo_extractions/…`), NOAA GHCN-d / ECCC / BCH obs and the validated
buoy table pickles — i.e. it requires the full data tree, which is why the
cache is not committed here.

Locate the pickle by either of:
1. the default path, i.e. `NOTEBOOKS/pickles/taylor_annuals.pkl` two
   directories above this folder, or
2. the environment variable `TAYLOR_ANNUALS_PKL=/path/to/taylor_annuals.pkl`.

## Regenerating the figure

    export TAYLOR_ANNUALS_PKL=/…/NOTEBOOKS/pickles/taylor_annuals.pkl   # if not in the default place
    MPLBACKEND=Agg /home/amh001/space_fs7/software_2022/python/py_2024/bin/python climate_all_taylor.py

or run the notebook (any kernel with numpy / scipy / matplotlib):

    jupyter nbconvert --to notebook --execute --inplace climate_all_taylor.ipynb

Outputs (written next to the script, overwritten each run):
`figures/climate_all_taylor.png` (600 dpi) and `.pdf` (vector).

### Expected stats (v10 pickle, 2026-10-02)

refstd (obs climatologies): t 2.130 K, pr 882.3 mm/yr, wind 2.093 m/s

| var | model   | n   | std      | corr  |
|-----|---------|-----|----------|-------|
| t   | D03     | 154 | 2.442    | 0.885 |
| t   | D02     | 154 | 2.606    | 0.783 |
| t   | D01     | 154 | 2.476    | 0.684 |
| t   | CanESM2 | 154 | 2.420    | 0.669 |
| t   | CanRCM4 | 154 | 3.297    | 0.757 |
| pr  | D03     | 159 | 710.9    | 0.944 |
| pr  | D02     | 159 | 1246.4   | 0.925 |
| pr  | D01     | 159 | 1015.5   | 0.830 |
| pr  | CanESM2 | 159 | 315.6    | 0.333 |
| pr  | CanRCM4 | 159 | 820.4    | 0.830 |
| wind| D03     | 23  | 1.894    | 0.963 |
| wind| D02     | 23  | 1.907    | 0.940 |
| wind| D01     | 23  | 1.970    | 0.846 |
| wind| CanESM2 | 23  | 1.521    | 0.759 |
| wind| CanRCM4 | 23  | 1.992    | 0.894 |

## Method notes

- Climatologies: `mean` over the finite annual values per station
  (t/wind: annual mean; pr: annual sum in mm/yr).
- refstd = std of obs climatologies over the union of all models'
  shared stations; each model's std/corr are computed on its shared
  station subset.
- Diagram geometry (axis range, RMS contour levels) auto-scales per
  panel via `taylor_figure._round_axis_ceiling` / `_auto_levels`.
- The wind panel's contour labels fall back to the mid-radius ring
  labelling because the shared corr=0.55 diagonal misses every contour
  there (see the comment in the plotting cell).

## Provenance

- Scripts copied from `NOTEBOOKS/` (CanESM2-WRF verification workspace)
  on 2026-10-02: `taylor_figure.py` verbatim; the `.py`/`.ipynb` pair is
  the executed code of `climate_all_taylor.ipynb` (last edited
  2026-10-01, "windlabels" revision) with only the pickle/output paths
  adjusted.
- Data pipeline (main workspace): `taylor_data.py` CACHE_VERSION 10,
  canonical land table `station_tables.pkl` (173 rows, USW-fixed),
  buoy model/obs from the validated buoy table pickles; see
  `taylor_3panel_handoff.md` in the main workspace for the full
  changelog (v10 added the CanESM2/CanRCM4 filtered D03 buoys).
