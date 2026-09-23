#!/usr/bin/env python
"""Generate 1-point SCRIP (CDO-compatible) BILINEAR remap weights for the buoys
against the CanESM2-WRF D01 99x99 curvilinear grid, read directly from the
source file wind_d01_hourly.nc (lat/lon are 2-D curvilinear, dims (lon,lat)).

CRITICAL (2026-09-23, fixed): the first version of this script transposed the
grid and wrote the declared source grid as lat.T.ravel().  CDO remap reads
wspd in FILE C-ORDER (dim 'lon' major), so every declared coordinate was the
SWAPPED neighbour's (the D01 grid is asymmetric, max |lat-lat.T| ~ 5 deg), and
CDO interpolated at the wrong physical location (e.g. buoy 46134 landed on a
cell near 48.9, -138.2 -> mean 8.45 instead of ~2.9 near 48.7, -124.4).

Correct layout (verified identical to the working D02/D03 weight files):
  * src_grid_dims = [nlon, nlat] as in the source file
  * src_grid_center_lat/lon = the source file's lat/lon fields in FILE C-ORDER
    (flat index k -> wspd[:, k//nlat, k%nlat], i.e. cell (a=k//nlat, b=k%nlat))
  * flat index of grid node (a, b) = a*nlat + b

Consistent with gen_buoy_wind_d03_weights.py / the D02 weights: bilinear over
the 4 surrounding source cells when the point is inside the domain, falling
back to a single nearest cell when the point is at/beyond the domain edge
(open ocean).  D01 has no corner field (lat/lon are cell centers), so the
point's cell is the nearest cell (by geodesic-ish distance to cell center,
refined over a local window), and bilinear interpolation is taken between the
4 surrounding cell-center values.  The residual sub-grid placement bias
(buoy position vs ~0.7 deg d01 cell spacing) is documented as a resolution
caveat in the report, per the locked "Option A" decision.

Usage:  python gen_buoy_wind_d01_weights.py [--force]
Output: weight_files_gen/CanESM2-WRF-D01/weights_st_<ID>.nc
"""
import csv, os, sys
import numpy as np
import netCDF4 as nc

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get(
    "WIND_D01_HOURLY",
    "/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/"
    "historical/variables_complete/wind_d01_hourly.nc")
OUTDIR = os.path.join(ROOT, "weight_files_gen", "CanESM2-WRF-D01")
CSV = os.path.join(ROOT, "station_cdo", "buoys", "buoy_descs.csv")


def load_grid():
    ds = nc.Dataset(SRC)
    la = np.array(ds.variables["lat"])   # (lon, lat) FILE ORDER
    lo = np.array(ds.variables["lon"])   # (lon, lat) FILE ORDER
    ds.close()
    return la, lo


def _dist_cell(a, b, plat, plon, lat, lon):
    """Geodesic-ish distance (deg) from point to cell center (a,b)."""
    kx = np.cos(np.radians(0.5 * (lat[a, b] + plat)))
    return np.sqrt((lat[a, b] - plat) ** 2 +
                   (kx * (lon[a, b] - plon)) ** 2)


def _find_cell(plat, plon, lat, lon):
    """Nearest cell (a,b) in file order; refined in a local 5x5 window from
    the coarse argmin."""
    nlon, nlat = lat.shape
    kx = np.cos(np.radians(plat))
    d2 = (lat - plat) ** 2 + (kx * (lon - plon)) ** 2
    a0, b0 = np.unravel_index(d2.argmin(), d2.shape)
    best_d, a, b = np.inf, a0, b0
    for dj in range(-2, 3):
        for di in range(-2, 3):
            a2, b2 = a0 + di, b0 + dj
            if not (0 <= a2 < nlon and 0 <= b2 < nlat):
                continue
            d = _dist_cell(a2, b2, plat, plon, lat, lon)
            if d < best_d:
                best_d, a, b = d, a2, b2
    return a, b


def _find_block(plat, plon, lat, lon):
    """Find the 2x2 cell block (top-left corner (a,b)) whose cell-center quad
    contains the point.  Vectorized over all blocks, like find_cell() in
    gen_buoy_wind_d03_weights.py.  Returns (a, b) or None if no quad contains
    the point (domain edge / open ocean)."""
    nlon, nlat = lat.shape
    # block (a,b) spans cells (a,b),(a+1,b),(a,b+1),(a+1,b+1)
    ny, nx = nlat - 1, nlon - 1
    js, is_ = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    js = js.ravel(); is_ = is_.ravel()
    c0 = lat[js, is_];      c1 = lat[js + 1, is_]
    c2 = lat[js + 1, is_ + 1]; c3 = lat[js, is_ + 1]
    c0l = lon[js, is_];      c1l = lon[js + 1, is_]
    c2l = lon[js + 1, is_ + 1]; c3l = lon[js, is_ + 1]
    yc = np.stack([c0, c1, c2, c3], axis=1)
    xc = np.stack([c0l, c1l, c2l, c3l], axis=1)
    m = ((yc.min(1) <= plat) & (yc.max(1) >= plat) &
         (xc.min(1) <= plon) & (xc.max(1) >= plon))
    if m.sum() == 0:
        return None
    if m.sum() > 1:
        # pick the smallest-area candidate (deepest interior)
        dlat = yc.max(1) - yc.min(1)
        dlon = xc.max(1) - xc.min(1)
        cand = np.where(m)[0]
        best = cand[np.argmin(dlat[cand] * dlon[cand])]
        a, b = int(js[best]), int(is_[best])
    else:
        a, b = int(js[m][0]), int(is_[m][0])
    return a, b


def bilinear_weights(plat, plon, lat, lon):
    """Return (src_flat_indices[4], weights[4]) bilinear around (plat, plon)
    on the D01 grid IN FILE ORDER (flat k = a*nlat + b).

    The point's 2x2 cell block is found by a full vectorized corner-containment
    test (as in gen_buoy_wind_d03_weights.py), then (u,v) are solved in the
    block's cell-center quad and clamped to [0,1].  When no quad contains the
    point (domain edge / open ocean) the block is taken as the nearest cell's
    top-left, and clamping keeps the weights positive and summing to 1."""
    nlon, nlat = lat.shape
    blk = _find_block(plat, plon, lat, lon)
    if blk is None:
        a, b = _find_cell(plat, plon, lat, lon)
    else:
        a, b = blk
    # block cells: (a,b), (a+1,b), (a,b+1), (a+1,b+1); clamp at domain edge
    a1 = min(a + 1, nlon - 1)
    b1 = min(b + 1, nlat - 1)
    # 4 corner points of the local quad (cell centers)
    p00 = (lat[a, b], lon[a, b])
    p10 = (lat[a1, b], lon[a1, b])
    p01 = (lat[a, b1], lon[a, b1])
    p11 = (lat[a1, b1], lon[a1, b1])
    # local tangent plane: u along a (00->10), v along b (00->01)
    eu_y, eu_x = p10[0] - p00[0], p10[1] - p00[1]
    ev_y, ev_x = p01[0] - p00[0], p01[1] - p00[1]
    det = eu_y * ev_x - eu_x * ev_y
    if abs(det) < 1e-12:
        u = v = 0.0
    else:
        u = ((plat - p00[0]) * ev_x - (plon - p00[1]) * ev_y) / det
        v = ((plon - p00[1]) * eu_y - (plat - p00[0]) * eu_x) / det
    u = min(max(u, 0.0), 1.0)
    v = min(max(v, 0.0), 1.0)
    w00 = (1 - u) * (1 - v)   # cell (a,   b)
    w10 = u * (1 - v)         # cell (a+1, b)
    w01 = (1 - u) * v         # cell (a,   b+1)
    w11 = u * v               # cell (a+1, b+1)
    idx = [a * nlat + b, a1 * nlat + b, a * nlat + b1, a1 * nlat + b1]
    wts = [w00, w10, w01, w11]
    return idx, wts


def write_scrip(path, plat, plon, lat, lon, idx, wts):
    nlon, nlat = lat.shape
    src_lat = np.ascontiguousarray(lat).ravel()
    src_lon = np.ascontiguousarray(lon).ravel()
    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("src_grid_rank", 2)
        ds.createDimension("dst_grid_rank", 2)
        ds.createDimension("src_grid_size", nlon * nlat)
        ds.createDimension("dst_grid_size", 1)
        ds.createDimension("num_links", 4)
        ds.createDimension("num_wgts", 4)
        v = ds.createVariable("src_grid_dims", "i4", ("src_grid_rank",)); v[:] = [nlon, nlat]
        v = ds.createVariable("dst_grid_dims", "i4", ("dst_grid_rank",)); v[:] = [1, 1]
        v = ds.createVariable("src_grid_center_lat", "f8", ("src_grid_size",))
        v.units = "radians"; v[:] = np.radians(src_lat)
        v = ds.createVariable("src_grid_center_lon", "f8", ("src_grid_size",))
        v.units = "radians"; v[:] = np.radians(src_lon)
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
        v[:] = np.radians(plon)
        ds.title = "SCRIP remapping with CDO (hand-written, bilinear 1-point)"
        ds.normalization = "none"
        ds.map_method = "Bilinear remapping"
        ds.conventions = "SCRIP"
        ds.source_grid = "generic"
        ds.dest_grid = "curvilinear"
        ds.history = ("23 Sep 2026 (rev2, fixed transposed-grid bug) : hand-written "
                      "SCRIP bilinear 1-point weights (%g,%g) vs wind_d01_hourly.nc "
                      "99x99 curvilinear grid, FILE C-ORDER centers; consistent with "
                      "D02/D03 bilinear pipeline" % (plat, plon))
        ds.CDO = "Climate Data Operators version 1.9.9 (https://mpimet.mpg.de/cdo)"


def main():
    force = "--force" in sys.argv
    lat, lon = load_grid()
    nlon, nlat = lat.shape
    os.makedirs(OUTDIR, exist_ok=True)
    with open(CSV) as f:
        rows = [r for r in csv.reader(f)
                if r and not r[0].startswith("#") and r[0] != "agency"]
    for a, sid, lon_s, lat_s in rows:
        if sid in ("desw1", "wpow1"):
            continue
        out = os.path.join(OUTDIR, "weights_st_%s.nc" % sid)
        if os.path.exists(out) and not force:
            print("SKIP %s" % os.path.basename(out))
            continue
        plat, plon = float(lat_s), float(lon_s)
        idx, wts = bilinear_weights(plat, plon, lat, lon)
        write_scrip(out, plat, plon, lat, lon, idx, wts)
        cells = []
        for k, w in zip(idx, wts):
            if abs(w) > 1e-12:
                aa, bb = divmod(k, nlat)
                cells.append("(%d,%d @ %.2f,%.2f w=%.3f)" % (aa, bb, lat[aa, bb], lon[aa, bb], w))
        print("WROTE %s  wsum=%.4f cells: %s"
              % (os.path.basename(out), sum(wts), "; ".join(cells)))


if __name__ == "__main__":
    main()
