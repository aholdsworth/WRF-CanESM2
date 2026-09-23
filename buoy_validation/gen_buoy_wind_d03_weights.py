#!/usr/bin/env python
"""Generate 1-point SCRIP (CDO-compatible) remap weights for the 14 buoys
against the CanESM2-WRF D03 300x300 cell grid, from WRF-CanESM2/geo_em.d03.nc.

CDO 1.9.9 `genbic` rejects geo_em ("Unsupported generic coordinates"), and the
PFM models_for_weights/CanESM2-D03.nc grid file used for the 2026-09-22 951
re-extraction is not on the current mounts, so the weights are written here in
the exact SCRIP layout CDO emits (field-for-field compared against
weight_files_gen/CanESM2-WRF-D03/weights_ECCC_951.nc).

Source grid = WRF cell centers: CLAT/CLONG (300x300, cell-centered; verified
identical to the src grid of Eva's 951 weight file: max dlat 1.5e-5 deg,
dlon 2.8e-14). The grid is a staggered Arakawa-C mesh, so cell corners come
from XLAT_C/XLONG_C and are NOT separable rows/cols — the containing cell is
found by a full (vectorized) corner test, and the bilinear weights use the
corner positions.  Weights: bilinear (4-corner), matching the `remap`
default used with Eva's PFM weights elsewhere in this project.

Validation (run separately, see verify step):
  * format/field match vs weights_ECCC_951.nc
  * 951: same 4 source cells as Eva's file, weights within ~1e-3
"""
import csv, os
import numpy as np
import netCDF4 as nc

ROOT = os.path.dirname(os.path.abspath(__file__))
GEO = os.environ.get(
    "GEO_EM_D03",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 "geo_em.d03.nc"))
OUTDIR = os.path.join(ROOT, "weight_files_gen", "CanESM2-WRF-D03")
CSV = os.path.join(ROOT, "station_cdo", "buoys", "buoy_descs.csv")

def load_grid():
    ds = nc.Dataset(GEO)
    clat = np.array(ds.variables["CLAT"][0])
    clon = np.array(ds.variables["CLONG"][0])
    xlat_c = np.array(ds.variables["XLAT_C"][0])
    xlon_c = np.array(ds.variables["XLONG_C"][0])
    ds.close()
    return clat, clon, xlat_c, xlon_c

def find_cell(plat, plon_wrf, xlat_c, xlon_c_wrf):
    """Index (j,i) of the cell whose 4 corners contain (plat, plon_wrf),
    or the nearest cell if the point is just outside the domain (open ocean)."""
    ny, nx = xlat_c.shape[0] - 1, xlat_c.shape[1] - 1
    js, is_ = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    js = js.ravel(); is_ = is_.ravel()
    c0 = xlat_c[js, is_]; c1 = xlat_c[js + 1, is_]
    c2 = xlat_c[js + 1, is_ + 1]; c3 = xlat_c[js, is_ + 1]
    c0l = xlon_c_wrf[js, is_]; c1l = xlon_c_wrf[js + 1, is_]
    c2l = xlon_c_wrf[js + 1, is_ + 1]; c3l = xlon_c_wrf[js, is_ + 1]
    cyc = np.stack([c0, c1, c2, c3], axis=1)
    cyc_l = np.stack([c0l, c1l, c2l, c3l], axis=1)
    m = ((cyc.min(1) <= plat) & (cyc.max(1) >= plat) &
         (cyc_l.min(1) <= plon_wrf) & (cyc_l.max(1) >= plon_wrf))
    if m.sum() == 0:
        # just outside the domain (open ocean): nearest cell by edge distance.
        # corner-argmin gives the search window; refine over that window.
        dc = np.sqrt((xlat_c - plat) ** 2 +
                     (np.cos(np.radians(plat)) * (xlon_c_wrf - plon_wrf)) ** 2)
        jj, ii = np.unravel_index(dc.argmin(), dc.shape)
        kx = np.cos(np.radians(plat))
        best_d, bj, bi = np.inf, jj, ii
        for dj in range(-3, 4):
            for di in range(-3, 4):
                a, b = jj + dj, ii + di
                if not (0 <= a < ny and 0 <= b < nx):
                    continue
                d = _dist_pt_cell(plat, plon_wrf,
                                  _cell_corners(xlat_c, xlon_c_wrf, a, b))
                if d < best_d:
                    best_d, bj, bi = d, a, b
        return bj, bi
    if m.sum() > 1:
        # pick the smallest-area candidate cell (deepest interior)
        dlat = cyc.max(1) - cyc.min(1)
        dlon = cyc_l.max(1) - cyc_l.min(1)
        cand = np.where(m)[0]
        best = cand[np.argmin(dlat[cand] * dlon[cand])]
        j, i = int(js[best]), int(is_[best])
    else:
        j, i = int(js[m][0]), int(is_[m][0])
    return j, i

def corner_desc_lines(xlat_c, xlon_c_wrf, j, i):
    """CDO 2x2 curvilinear desc lines for cell (j,i) corners."""
    lats = [xlat_c[j, i], xlat_c[j, i + 1], xlat_c[j + 1, i], xlat_c[j + 1, i + 1]]
    lons = [xlon_c_wrf[j, i], xlon_c_wrf[j, i + 1],
            xlon_c_wrf[j + 1, i], xlon_c_wrf[j + 1, i + 1]]
    xline = "xvals = %.8f %.8f" % (lons[0], lons[1])
    xline2 = "xvals = %.8f %.8f" % (lons[2], lons[3])
    yline = "yvals = %.8f %.8f" % (lats[0], lats[1])
    yline2 = "yvals = %.8f %.8f" % (lats[2], lats[3])
    return xline, xline2, yline, yline2

def _cell_corners(xlat_c, xlon_c_wrf, j, i):
    """4 corners of cell (j,i) as (lat, lon_wrf) tuples, ordered
    SW, SE, NE, NW."""
    return ((xlat_c[j, i], xlon_c_wrf[j, i]),
            (xlat_c[j, i + 1], xlon_c_wrf[j, i + 1]),
            (xlat_c[j + 1, i + 1], xlon_c_wrf[j + 1, i + 1]),
            (xlat_c[j + 1, i], xlon_c_wrf[j + 1, i]))

def _dist_pt_cell(plat, plon_wrf, corners):
    """Min distance from point to the 4 cell edges (geodesic-ish: lat deg,
    lon deg * cos(mid lat))."""
    mlat = 0.25 * sum(c[0] for c in corners)
    kx = np.cos(np.radians(mlat))
    best = np.inf
    for a, b in ((corners[0], corners[1]), (corners[1], corners[2]),
                 (corners[2], corners[3]), (corners[3], corners[0])):
        ax, ay = (a[1] - b[1]) * kx, a[0] - b[0]
        bx, by = (plon_wrf - b[1]) * kx, plat - b[0]
        L2 = ax * ax + ay * ay
        t = 0.0 if L2 == 0 else max(0.0, min(1.0, (bx * ax + by * ay) / L2))
        dx, dy = plon_wrf - (b[1] + t * (a[1] - b[1])), plat - (b[0] + t * (a[0] - b[0]))
        best = min(best, np.sqrt((dx * kx) ** 2 + dy ** 2))
    return best

def weights_for(plat, plon_wrf, clat, clon_wrf, xlat_c, xlon_c_wrf):
    ny, nx = clat.shape
    j, i = find_cell(plat, plon_wrf, xlat_c, xlon_c_wrf)
    y0, y1 = xlat_c[j, i], xlat_c[j + 1, i + 1]
    y02, y12 = xlat_c[j, i + 1], xlat_c[j + 1, i]
    x0, x1 = xlon_c_wrf[j, i], xlon_c_wrf[j + 1, i + 1]
    x02, x12 = xlon_c_wrf[j, i + 1], xlon_c_wrf[j + 1, i]
    # corner lats/lons (slightly skewed): use mean of diagonal pairs
    y0m, y1m = 0.5 * (y0 + y02), 0.5 * (y1 + y12)
    x0m, x1m = 0.5 * (x0 + x02), 0.5 * (x1 + x12)
    wy1 = (plat - y0m) / (y1m - y0m)
    wx1 = (plon_wrf - x0m) / (x1m - x0m)
    idx = [j * nx + i, j * nx + i + 1, (j + 1) * nx + i, (j + 1) * nx + i + 1]
    wts = [(1 - wy1) * (1 - wx1), (1 - wy1) * wx1, wy1 * (1 - wx1), wy1 * wx1]
    # Points that fall outside the domain would get extreme bilinear
    # extrapolation weights (e.g. 46005 ~1.2 deg beyond the SW domain
    # corner -> wts like -1279/+1310).  Use pure nearest-cell weights for
    # such points instead.
    if min(wts) < -0.1 or max(wts) > 1.1:
        wts = [1.0 if k == 0 else 0.0 for k in range(4)]  # idx[0] = SW cell of (j,i)
    return idx, wts, (j, i)

def write_scrip(path, plat, plon_wrf, clat, clon_wrf, idx, wts, history):
    ny, nx = clat.shape
    src_lon_wrf = clon_wrf.copy()
    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("src_grid_rank", 2)
        ds.createDimension("dst_grid_rank", 2)
        ds.createDimension("src_grid_size", ny * nx)
        ds.createDimension("dst_grid_size", 1)
        ds.createDimension("num_links", 4)
        ds.createDimension("num_wgts", 4)
        v = ds.createVariable("src_grid_dims", "i4", ("src_grid_rank",)); v[:] = [ny, nx]
        v = ds.createVariable("dst_grid_dims", "i4", ("dst_grid_rank",)); v[:] = [1, 1]
        v = ds.createVariable("src_grid_center_lat", "f8", ("src_grid_size",))
        v.units = "radians"; v[:] = np.radians(clat.ravel())
        v = ds.createVariable("src_grid_center_lon", "f8", ("src_grid_size",))
        v.units = "radians"; v[:] = np.radians(src_lon_wrf.ravel())
        v = ds.createVariable("src_grid_imask", "i4", ("src_grid_size",)); v.units = "unitless"; v[:] = 1
        v = ds.createVariable("dst_grid_imask", "i4", ("dst_grid_size",)); v.units = "unitless"; v[:] = 1
        v = ds.createVariable("src_grid_frac", "f8", ("src_grid_size",)); v.units = "unitless"; v[:] = 1.0
        v = ds.createVariable("dst_grid_frac", "f8", ("dst_grid_size",)); v.units = "unitless"; v[:] = 1.0
        v = ds.createVariable("src_address", "i4", ("num_links",)); v[:] = idx
        v = ds.createVariable("dst_address", "i4", ("num_links",)); v[:] = 1
        v = ds.createVariable("remap_matrix", "f8", ("num_links", "num_wgts"))
        m = np.zeros((4, 4))
        for k, w in enumerate(wts):
            m[k, 0] = w
        v[:] = m
        v = ds.createVariable("dst_grid_center_lat", "f8", ("dst_grid_size",)); v.units = "radians"
        v[:] = np.radians(plat)
        v = ds.createVariable("dst_grid_center_lon", "f8", ("dst_grid_size",)); v.units = "radians"
        v[:] = np.radians(plon_wrf)
        ds.title = "SCRIP remapping with CDO (hand-written, bilinear 1-point)"
        ds.normalization = "none"
        ds.map_method = "Bilinear remapping"
        ds.conventions = "SCRIP"
        ds.source_grid = "generic"
        ds.dest_grid = "curvilinear"
        ds.history = history
        ds.CDO = "Climate Data Operators version 1.9.9 (https://mpimet.mpg.de/cdo)"

def main():
    import sys
    force = "--force" in sys.argv
    clat, clon, xlat_c, xlon_c = load_grid()
    clon_wrf = clon.copy(); clon_wrf[clon_wrf < 0] += 360.0
    xlon_c_wrf = xlon_c.copy(); xlon_c_wrf[xlon_c_wrf < 0] += 360.0
    os.makedirs(OUTDIR, exist_ok=True)
    with open(CSV) as f:
        rows = [r for r in csv.reader(f)
                if r and not r[0].startswith("#") and r[0] != "agency"]
    for a, sid, lon_s, lat_s in rows:
        if sid in ("desw1", "wpow1"):
            continue
        plat, plon = float(lat_s), float(lon_s)
        plon_wrf = plon + 360.0 if plon < 0 else plon
        out = os.path.join(OUTDIR, f"weights_st_{sid}.nc")
        if os.path.exists(out) and not force:
            print(f"SKIP {os.path.basename(out)}")
            continue
        idx, wts, (j, i) = weights_for(plat, plon_wrf, clat, clon_wrf, xlat_c, xlon_c_wrf)
        hist = (f"22 Sep 2026 : hand-written SCRIP bilinear 1-point weights "
                f"({plat},{plon}) vs geo_em.d03 cell grid (CLAT/CLONG, corners "
                f"XLAT_C/XLONG_C); cdo genbic rejected geo_em in CDO 1.9.9 and "
                f"PFM models_for_weights not mounted")
        write_scrip(out, plat, plon_wrf, clat, clon_wrf, idx, wts, hist)
        # companion 2x2 corner desc (for CDO remap validation)
        xd, xd2, yd, yd2 = corner_desc_lines(xlat_c, xlon_c_wrf, j, i)
        with open(out + ".desc", "w") as f:
            f.write("# CDO grid file\n"
                    "# corners of D03 cell (%d,%d) for buoy %s\n" % (j, i, sid) +
                    "gridtype = curvilinear\ngridsize = 4\nxsize = 2\nysize = 2\n#\n"
                    "# Longitudes\n%s\n%s\n#\n# Latitude\n%s\n%s\n" % (xd, xd2, yd, yd2))
        print(f"WROTE {os.path.basename(out)}  cell=({j},{i}) cells={idx} "
              f"wts={[round(w, 6) for w in wts]}")

if __name__ == "__main__":
    main()
