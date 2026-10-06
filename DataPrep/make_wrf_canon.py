#!/usr/bin/env python
"""
make_wrf_canon.py -- build data/wrf_canon/: ONE canonical naming convention for
all WRF per-station hourly files (t / pr / wind x d01 / d02 / d03).

Canonical name:  {var}_{dom}_{AGENCY}_{ID}.nc
    AGENCY = ECCC | NOAA | BCH
    ID     = ClimateID (ECCC land), USW/USC id (NOAA land), 3-letter (BCH),
             46xxx (buoys)

Sources (all copied byte-identical EXCEPT the 7 ECCC d01-buoy wind files,
which additionally get their corrupted embedded lat/lon corrected -- see
"Known issues" below):

  d01        data/wrf_pfm_cache/{var}_d01_{AGENCY}_{ID}.nc
             (also untagged {var}_d01_{ID}.nc for 3-letter BCH)
  d02 / d03  data/wrf_pfm_cache/{var}_{dom}_st{ID}.nc  -> {AGENCY}_{ID}
             data/wrf_pfm_cache/{var}_{dom}_BCH_{ID}.nc (unchanged)

Known issue fixed here (documented 2024-2025, wind_dir_buoy_audit.log,
buoy_wind_bias_handoff_20260922.md): the ECCC d01 BUOY wind files carry
ROTATED embedded coordinates (6 of 7 hold a different buoy's location).
We re-extract wspd from PFM source wind_d01_hourly.nc with the verified
1-point CDO weights (weight_files_gen/CanESM2-WRF-D01/weights_st_<ID>.nc,
truth coords from station_cdo/buoys/buoy_descs.csv) into
data/wrf_stations_tmp/wind_d01_ECCC_<ID>.nc, verify the embedded coords, and
then copy.  NOAA buoy d01 files (46029/46041) are copied as-is (coords OK).

Idempotent: re-runs verify and skip.  Run:
  /home/amh001/space_fs7/software_2022/python/py_2024/bin/python make_wrf_canon.py
"""
import filecmp, glob, os, re, shutil, subprocess, sys
import numpy as np
from netCDF4 import Dataset

BASE = os.path.dirname(os.path.abspath(__file__))
PFM = os.path.join(BASE, 'data', 'wrf_pfm_cache')
# Read-only PFM per-station extraction dir (d01 ECCC_202 untagged, and the
# d01 USW / d02+d03 st1001 land wind files that are NOT in wrf_pfm_cache).
PFM_STATIONS = ('/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
                'historical/variables_stations')
CANON = os.path.join(BASE, 'data', 'wrf_canon')
WRF_ST = os.path.join(BASE, 'data', 'wrf_stations')          # 2026-09-22 truth-matched d01 buoy files
TMP = os.path.join(BASE, 'data', 'wrf_stations_tmp')         # fresh truth re-extractions
TRUTH = {  # (lat, lon) -- from station_cdo/buoys/buoy_descs.csv (2026-09-22 verified)
    '46131': (49.91, -124.99), '46132': (49.74, -127.93),
    '46134': (48.65, -123.50), '46146': (49.34, -123.73),
    '46204': (51.38, -128.77), '46206': (48.83, -126.00),
    '46207': (50.88, -129.91), '46029': (46.163, -124.487),
    '46041': (47.352, -124.739)}
ECCC_BUOY = ['46131', '46132', '46134', '46146', '46204', '46206', '46207']
NOAA_BUOY = ['46029', '46041']
CDO = '/home/amh001/space_fs7/software_2022/python/cdo_env/bin/cdo'
PFM_SRC_WIND_D01 = ('/gpfs/fs7/dfo/hpcmc/pfm/spfm000/CanESM2-WRF/'
                    'historical/variables_complete/wind_d01_hourly.nc')
DESC = os.path.join(BASE, 'station_cdo', 'buoys', 'descs', 'ECCC_buoy_%s.txt')
WGT = os.path.join(BASE, 'weight_files_gen', 'CanESM2-WRF-D01', 'weights_st_%s.nc')


NOAA_BUOY_IDS = ('46005', '46029', '46041', '46050', '46087', '46088',
                 '46089')


def agency_of(sid):
    if sid.startswith('46'):
        return 'NOAA' if sid in NOAA_BUOY_IDS else 'ECCC'
    if sid in ('desw1', 'wpow1'):
        # legacy: these 3-letter ids are NOAA buoys
        return 'NOAA'
    if sid.startswith('USC') or sid.startswith('USW'):
        return 'NOAA'
    if sid.isdigit():
        return 'ECCC'
    return 'BCH'


def agency_of_buoy(sid):
    return 'NOAA' if sid in NOAA_BUOY_IDS else 'ECCC'


def canonical_name(var, dom, sid):
    if sid.startswith('46') and sid not in NOAA_BUOY_IDS:
        # ECCC buoys
        return '%s_%s_ECCC_%s.nc' % (var, dom, sid)
    return '%s_%s_%s_%s.nc' % (var, dom, agency_of(sid), sid)


def embedded_coords(path):
    nc = Dataset(path)
    la = float(np.asarray(nc.variables['lat'][:]).ravel()[0])
    lo = float(np.asarray(nc.variables['lon'][:]).ravel()[0])
    nc.close()
    return la, lo


def same(a, b, atol=0.02):
    return abs(a[0] - b[0]) <= atol and abs(a[1] - b[1]) <= atol


# ---------------------------------------------------------------------------
# Step 1: fresh D01 ECCC-buoy wspd extractions at truth coords (if source avail.)
# ---------------------------------------------------------------------------
def refresh_d01_buoy_wind():
    """Re-extract the 7 ECCC buoy d01 wind files at truth coords.
    Returns dict id -> path, using:
      * fresh extraction at TMP (if PFM source readable + CDO works), else
      * the 2026-09-22 verified WRF_ST files (same provenance, coord-checked).
    """
    out = {}
    src_ok = os.path.exists(PFM_SRC_WIND_D01)
    os.makedirs(TMP, exist_ok=True)
    print('PFM source wind_d01_hourly.nc readable: %s' % src_ok)
    for bid in ECCC_BUOY:
        tmp = os.path.join(TMP, 'wind_d01_ECCC_%s.nc' % bid)
        if src_ok and not (os.path.exists(tmp) and os.path.getsize(tmp) > 0):
            r = subprocess.run(
                [CDO, '-b', 'F64', '-f', 'nc',
                 'remap,%s,%s' % (DESC % bid, WGT % bid),
                 PFM_SRC_WIND_D01, tmp],
                capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(tmp):
                print('  %s: FRESH EXTRACTION FAILED (%s) -> fallback to 2026-09-22 file'
                      % (bid, r.stderr.strip().splitlines()[-1] if r.stderr else '?'))
        if src_ok and os.path.exists(tmp) and os.path.getsize(tmp) > 0:
            la, lo = embedded_coords(tmp)
            assert same((la, lo), TRUTH[bid]), \
                'fresh %s coords (%.3f,%.3f) != truth' % (bid, la, lo)
            out[bid] = tmp
            print('  %s: fresh extraction @ (%.3f, %.3f) OK' % (bid, la, lo))
        else:
            p = os.path.join(WRF_ST, 'wind_d01_ECCC_buoy_%s.nc' % bid)
            assert os.path.exists(p), 'no fallback file for %s' % bid
            la, lo = embedded_coords(p)
            assert same((la, lo), TRUTH[bid]), \
                '2026-09-22 file %s coords (%.3f,%.3f) != truth' % (bid, la, lo)
            out[bid] = p
            print('  %s: using 2026-09-22 verified file @ (%.3f, %.3f)' % (bid, la, lo))
    return out


# ---------------------------------------------------------------------------
# Supplemental sources for land files missing from wrf_pfm_cache.
#
# wrf_pfm_cache is incomplete for a handful of LAND stations: it has no d01
# ECCC_82 / ECCC_202 / 3-letter BCH t-pr, no d01 USW wind, and no d02/d03
# ECCC-1001 (st1001) wind.  Those live in data/wrf_stations (d01 t/pr land) and
# the read-only PFM variables_stations dir (d01 202 + all land wind).  Without
# this step a rebuild would regress them out of the canon, so we top the canon
# up from the supplemental sources (validated: embedded coords == station,
# 20-yr span, byte-identical).  Idempotent: only fills names not already made.
# ---------------------------------------------------------------------------
SUPP_BCH_D01 = ('ALU','BCK','BLN','CLO','CMU','CMX','COQ','CQM','CRU','DAI','DLU',
                'ECL','GOC','MIS','NTY','STA','WAH','WOL')


def fill_supplemental_land():
    """Copy/link land files absent from wrf_pfm_cache into the canon from the
    supplemental sources.  Returns (n_added, n_missing)."""
    added, missing = 0, []
    def place(canonical, src):
        nonlocal added
        dst = os.path.join(CANON, canonical)
        if os.path.exists(dst):
            return  # already made by the PFM scan
        if not os.path.exists(src):
            missing.append(canonical)
            return
        if os.path.exists(dst):
            os.remove(dst)
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)
        added += 1

    # d01 t/pr land: 18 BCH + ECCC_82 from wrf_stations (BCH/82 are tagged)
    for var in ('t', 'pr'):
        for code in SUPP_BCH_D01:
            f = '%s_d01_BCH_%s.nc' % (var, code)
            place(f, os.path.join(WRF_ST, f))
        f82 = '%s_d01_ECCC_82.nc' % var
        place(f82, os.path.join(WRF_ST, f82))
        # d01 ECCC_202: untagged in PFM (t_d01_202.nc / pr_d01_202.nc)
        place('%s_d01_ECCC_202.nc' % var,
              os.path.join(PFM_STATIONS, '%s_d01_202.nc' % var))

    # land wind: d02/d03 ECCC-1001 (st1001) from PFM
    for dom in ('d02', 'd03'):
        place('wind_%s_ECCC_1001.nc' % dom,
              os.path.join(PFM_STATIONS, 'wind_%s_st1001.nc' % dom))
    # d01 USW wind (5) from PFM (rename USW* -> NOAA_USW*)
    for sid in ('00024227', '00024229', '00024233', '00024243', '00094224'):
        place('wind_d01_NOAA_USW%s.nc' % sid,
              os.path.join(PFM_STATIONS, 'wind_d01_USW%s.nc' % sid))
    # d01 ECCC-1001 land wind from wrf_stations (in wrf_stations, not PFM)
    place('wind_d01_ECCC_1001.nc',
          os.path.join(WRF_ST, 'wind_d01_ECCC_1001.nc'))
    return added, missing


def main():
    os.makedirs(CANON, exist_ok=True)
    n_copy, n_fix, n_missing = 0, 0, []

    d01_buoy_src = refresh_d01_buoy_wind()

    # ---- scan PFM for all per-station files
    files = {}
    for f in os.listdir(PFM):
        m = re.match(r'^(t|pr|wind)_(d0[123])_(.+)\.nc$', f)
        if not m:
            continue
        var, dom, rest = m.groups()
        if dom == 'd01':
            # d01 files carry an explicit AGENCY_ prefix (ECCC/NOAA/BCH)
            # for everything EXCEPT the untagged 3-letter BCH land stations
            for ag in ('ECCC_', 'NOAA_', 'BCH_'):
                if rest.startswith(ag):
                    files.setdefault((var, dom), {})[rest[len(ag):]] = f
                    break
            else:
                # untagged d01 (3-letter BCH land stations)
                files.setdefault((var, dom), {})[rest] = f
        else:
            # d02/d03 source files: st{AGENCY}_{ID} (e.g. stECCC_82,
            # stNOAA_USC00351643, stBCH_ECCC_46131), {AGENCY}_{ID} (ECCC/NOAA
            # d02/d03), or untagged 3-letter BCH land stations
            r = rest
            if r.startswith('st'):
                r = r[2:]
            if r.startswith('ECCC_') or r.startswith('NOAA_') or r.startswith('BCH_'):
                r = r.split('_', 1)[1]
            files.setdefault((var, dom), {})[r] = f

    # ---- wipe the canon dir and rebuild from scratch (idempotent)
    for f in os.listdir(CANON):
        os.remove(os.path.join(CANON, f))

    for (var, dom), stmap in sorted(files.items()):
        for sid, fname in sorted(stmap.items()):
            src = os.path.join(PFM, fname)
            if dom == 'd01' and var == 'wind' and sid in ECCC_BUOY:
                # D01 ECCC buoy wind: use verified truth-coord extraction
                # (the wrf_pfm_cache copies have rotated embedded coords)
                src = d01_buoy_src[sid]
                n_fix += 1
            elif dom == 'd01' and var == 'wind' and sid in NOAA_BUOY:
                la, lo = embedded_coords(src)
                assert same((la, lo), TRUTH[sid]), \
                    'NOAA buoy %s d01 coords wrong: (%.3f,%.3f)' % (sid, la, lo)
            # 46xxx files are ECCC buoys (d02/d03 wind copies are d01-era
            # artifacts, 1 yr of data); canonical name always ECCC
            dst = os.path.join(CANON, canonical_name(var, dom, sid))
            # hardlink when possible (same fs -> zero extra disk, instant);
            # fallback to byte copy otherwise.  src is never modified in
            # place, so the link is safe.
            if os.path.exists(dst):
                os.remove(dst)
            try:
                os.link(src, dst)
            except OSError:
                shutil.copy2(src, dst)
            n_copy += 1

    print('\ncanonical dir: %s' % CANON)
    print('files copied: %d  (supplemental land: %d; d01 ECCC-buoy wind re-extracted: %d)'
          % (n_copy, n_supp, n_fix))
    # supplemental land files (d01 BCH/82/202 t-pr, d01 USW + d02/03 1001 wind)
    n_supp, supp_missing = fill_supplemental_land()
    if supp_missing:
        print('  WARN supplemental land sources missing:', supp_missing)
    n_copy += n_supp

    # summary per (var, dom)
    for var in ('t', 'pr', 'wind'):
        for dom in ('d01', 'd02', 'd03'):
            n = len(glob.glob(os.path.join(CANON, '%s_%s_*' % (var, dom))))
            print('  %-16s %3d' % ('%s_%s' % (var, dom), n))


if __name__ == '__main__':
    main()
