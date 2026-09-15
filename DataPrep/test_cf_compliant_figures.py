#!/home/amh001/space_fs7/software_2022/python/py_2024/bin/python
"""
Minimal test: does the CF-compliant data reproduce the paper's regional-mean
figures?

Loads the CF-compliant daily files (pr, t2, wspd x historical/rcp45/rcp85) plus
the land mask from geo_em.d03.nc, and replicates the *exact* pipeline in
`land_and_ocean/Regional_Averages.ipynb` (get_monthly_climatologies):

    daily -> resample(time="1MS").mean()          [monthly-mean time series]
          -> cos(lat) area-weighted land / ocean mean time series
          -> monthly climatology (mean over the years of that month)
          -> delta = fut - hist (absolute)  |  100*(fut-hist)/hist (relative)
          -> SEM = sqrt(hist_sem^2 + fut_sem^2)
          -> seasonal amplitude = mean(JJA) - mean(DJF)

The heavy lifting is done in numpy (fast) rather than xarray so that reading
9 x 7305x300x300 float64 files stays practical. The operations are the same
as the notebook: a daily->monthly mean is exactly the mean of all days in the
month (the extra groupby-month at the climatology step is a no-op for a
monthly-mean series), and the weighted regional mean is the cos(lat)-weighted
sum of masked cells divided by the sum of the weights.

Run:
    python DataPrep/test_cf_compliant_figures.py

Outputs (in Figures/):
    Regional_averages_relative_CF.png   paper-style figure to eyeball
    cf_regional_means.npz               every array, for a numeric diff later
"""
import os
import glob
import numpy as np
import netCDF4

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CC = "/gpfs/fs7/dfo/hpcmc/comda/fs2_comda/amh001/AI_WORKDIR/WRF_FILES/cf_compliant/data"
GEO_EM = os.path.join(REPO, "geo_em.d03.nc")
PLOT_DIR = os.path.join(REPO, "Figures")
os.makedirs(PLOT_DIR, exist_ok=True)

# logical var -> (CF data-variable name, CF filename token, use-relative-delta)
SPEC = {
    "t":    ("T2",   "t2",   False),
    "pr":   ("pr",   "pr",   True),
    "wind": ("wspd", "wspd", True),
}
SCEN = ("historical", "rcp45", "rcp85")


def cf_path(var, scenario):
    dvar, token, _ = SPEC[var]
    hits = glob.glob(f"{CC}/{scenario}/CanESM2-WRF_{token}_d03_{scenario}_daily_*.nc")
    assert len(hits) == 1, f"expected 1 CF file {scenario}/{var}, got {hits}"
    return hits[0]


def load_mask():
    """Land/ocean masks + cos-lat weights, on the (y, x) grid."""
    ds = netCDF4.Dataset(GEO_EM, "r")
    lm = ds["LANDMASK"][:].squeeze()          # (300,300)
    xlat = ds["XLAT_M"][:].squeeze()          # (300,300)
    ds.close()
    land = lm == 1
    weights = np.cos(np.deg2rad(xlat))
    return land, ~land, weights


def monthly_climatology_dailymean(path, dvar, land, ocean, weights, NMONTH=12):
    """Return (delta vs hist already handled by caller); here: monthly climatology
    land/ocean time series from one CF file, computed directly from daily data
    without materializing the monthly-mean 3D array twice.

    We read the full (time, y, x) array once, then:
      - daily -> monthly mean per (year, month)
      - region mean (cos-lat weighted, masked)
      - groupby month .mean / .std / .count for climatology + SEM
    """
    ds = netCDF4.Dataset(path, "r")
    data = ds[dvar][:]                      # (T, y, x) float64
    dt = ds["time"][:].astype("float64")
    tu = ds["time"].getncattr("units")
    cal = ds["time"].getncattr("calendar") if "calendar" in ds["time"].ncattrs() else "standard"
    ds.close()

    # decode to (month, day) to get day-of-month and calendar month
    dates = netCDF4.num2date(dt, tu, cal)
    month = np.array([d.month for d in dates], dtype="int64")       # 1..12
    day = np.array([d.day for d in dates], dtype="int64")
    T = data.shape[0]

    # (y,x) region scalars: land-weighted and ocean-weighted means over the grid
    wl = weights * land
    wo = weights * ocean
    wl_sum = wl.sum()
    wo_sum = wo.sum()
    # weighted sum over grid for each day
    day_land = (data * wl).sum(axis=(1, 2)) / wl_sum      # (T,)
    day_ocean = (data * wo).sum(axis=(1, 2)) / wo_sum     # (T,)

    # monthly mean time series (mean of days in each year-month)
    nmonth = NMONTH
    # use month*100+day? no: group by calendar month only (climatology), but SEM
    # needs the per-year monthly-mean std. The notebook does:
    #   hist_ts = resample(time='1MS').mean()   # (years, months) series
    #   .groupby('time.month').mean()           # mean over years for that month
    #   .groupby('time.month').std(ddof=1)/sqrt(count)
    # Build per-year-month means: key = year*12 + (month-1)
    years = np.array([d.year for d in dates], dtype="int64")
    ym = (years - years.min()) * 12 + (month - 1)          # 0..nyears*12-1
    nym = ym.max() + 1

    def per_ym_series(day_vals):
        s = np.bincount(ym, weights=day_vals, minlength=nym)
        c = np.bincount(ym, minlength=nym)
        return s / np.where(c > 0, c, 1), c                # (nym,), (nym,)

    land_ym, cnt = per_ym_series(day_land)
    ocean_ym, _ = per_ym_series(day_ocean)

    # group by calendar month (month index 0..11 across all years)
    # ym encodes year*12 + (month-1); month index of a given ym slot:
    ym_month = (np.arange(nym) % 12)                       # 0..11

    def clim_std(ym_series, ym_cnt):
        # mean over years for each calendar month
        s = np.bincount(ym_month, weights=ym_series * (ym_cnt > 0), minlength=nmonth)
        c = np.bincount(ym_month, weights=(ym_cnt > 0), minlength=nmonth)
        m = s / np.where(c > 0, c, 1)
        # std(ddof=1) over years for each calendar month
        s2 = np.bincount(ym_month, weights=(ym_series ** 2) * (ym_cnt > 0), minlength=nmonth)
        var = (s2 / c - m ** 2) * (c / np.where(c - 1 > 0, c - 1, 1))
        var = np.where(c - 1 > 0, var, 0.0)
        std = np.sqrt(np.maximum(var, 0.0))
        return m, std, c

    lm_, ls_, lc = clim_std(land_ym, cnt)
    om_, os_, oc = clim_std(ocean_ym, cnt)

    return dict(
        land_clim=lm_, ocean_clim=om_,
        land_std=ls_, ocean_std=os_,
        land_cnt=lc, ocean_cnt=oc,
        months=month,
    )


def get_monthly_climatologies(var):
    dvar, token, relative = SPEC[var]
    land, ocean, weights = load_mask()

    hist = monthly_climatology_dailymean(cf_path(var, "historical"), dvar, land, ocean, weights)
    r45 = monthly_climatology_dailymean(cf_path(var, "rcp45"), dvar, land, ocean, weights)
    r85 = monthly_climatology_dailymean(cf_path(var, "rcp85"), dvar, land, ocean, weights)

    def delta_sem(h, f):
        hc, fc = h["land_clim"], f["land_clim"]
        ho, fo = h["ocean_clim"], f["ocean_clim"]
        if relative:
            dl = 100.0 * (fc - hc) / hc
            do = 100.0 * (fo - ho) / ho
            sem_l = 100.0 * np.sqrt(h["land_std"] ** 2 + f["land_std"] ** 2) / np.abs(hc)
            sem_o = 100.0 * np.sqrt(h["ocean_std"] ** 2 + f["ocean_std"] ** 2) / np.abs(ho)
        else:
            dl = fc - hc
            do = fo - ho
            sem_l = np.sqrt(h["land_std"] ** 2 + f["land_std"] ** 2)
            sem_o = np.sqrt(h["ocean_std"] ** 2 + f["ocean_std"] ** 2)
        return dl, do, sem_l, sem_o

    d45l, d45o, s45l, s45o = delta_sem(hist, r45)
    d85l, d85o, s85l, s85o = delta_sem(hist, r85)

    return {
        "historical": {"land": hist["land_clim"], "ocean": hist["ocean_clim"]},
        "rcp45": {"land": r45["land_clim"], "ocean": r45["ocean_clim"]},
        "rcp85": {"land": r85["land_clim"], "ocean": r85["ocean_clim"]},
        "delta45": {"land": d45l, "ocean": d45o, "land_sem": s45l, "ocean_sem": s45o},
        "delta85": {"land": d85l, "ocean": d85o, "land_sem": s85l, "ocean_sem": s85o},
    }


def seasonal_amplitude(clim):
    summer = [5, 6, 7]     # Jun, Jul, Aug
    winter = [11, 0, 1]    # Dec, Jan, Feb
    return np.mean(clim[summer]) - np.mean(clim[winter])


def get_amps(results):
    hl = seasonal_amplitude(results["historical"]["land"])
    ho = seasonal_amplitude(results["historical"]["ocean"])
    r45l = seasonal_amplitude(results["rcp45"]["land"])
    r45o = seasonal_amplitude(results["rcp45"]["ocean"])
    r85l = seasonal_amplitude(results["rcp85"]["land"])
    r85o = seasonal_amplitude(results["rcp85"]["ocean"])
    return (100.0 * (r45l - hl) / hl, 100.0 * (r45o - ho) / ho,
            100.0 * (r85l - hl) / hl, 100.0 * (r85o - ho) / ho)


def plot_error(ax, res, months):
    caps, capw, capt, lw = 2.5, 0.75, 0.75, 2
    color3, color2 = "#364B9A", "#8CB6D0"
    color4, color1 = "#D55E00", "#E69F00"
    ax.errorbar(months, res["delta45"]["ocean"], yerr=res["delta45"]["ocean_sem"],
                label="Ocean RCP4.5", marker="o", color=color2, linewidth=lw,
                capsize=caps, capthick=capt, elinewidth=capw)
    ax.errorbar(months, res["delta45"]["land"], yerr=res["delta45"]["land_sem"],
                label="Land RCP4.5", marker="d", linestyle="-", color=color1,
                markerfacecolor="None", markeredgecolor=color1, linewidth=lw,
                capsize=caps, capthick=capt, elinewidth=capw)
    ax.errorbar(months, res["delta85"]["ocean"], yerr=res["delta85"]["ocean_sem"],
                label="Ocean RCP8.5", marker="o", color=color3, linewidth=lw,
                capsize=caps, capthick=capt, elinewidth=capw)
    ax.errorbar(months, res["delta85"]["land"], yerr=res["delta85"]["land_sem"],
                label="Land RCP8.5", marker="d", linestyle="-", color=color4,
                markerfacecolor="None", markeredgecolor=color4, linewidth=lw,
                capsize=caps, capthick=capt, elinewidth=capw)


def amp_bar(cax, amps):
    color3, color2 = "#364B9A", "#8CB6D0"
    color4, color1 = "#D55E00", "#E69F00"
    d45_l, d45_o, d85_l, d85_o = amps
    x = np.arange(2); width = 0.35
    cax.bar(x - width/2, [d45_l, d45_o], width, label="RCP 4.5", color=[color1, color2])
    cax.bar(x + width/2, [d85_l, d85_o], width, label="RCP8.5", color=[color4, color3])
    cax.axhline(0, linestyle="--", linewidth=1, color="black")
    cax.set_xticks(x); cax.set_xticklabels(["Land", "Ocean"])
    cax.spines["top"].set_visible(False); cax.spines["right"].set_visible(False)


def make_figure(res_t, res_pr, res_wind):
    months = np.arange(1, 13)
    month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    subpanel = ["(a)", "(c)", "(e)"]

    plt.figure(figsize=(8, 7))

    ax = plt.subplot(3, 1, 1)
    plot_error(ax, res_t, months)
    ax.set_ylabel(r"$\Delta$ T ($^{\circ}$C)")
    ax.grid(True); ax.set_xticks(months); ax.set_xticklabels([])
    ax.legend(fontsize=8, framealpha=0.5)
    ax.set_ylim([0, 5])
    ax.text(0.01, 0.95, subpanel[0], transform=ax.transAxes, va="top", zorder=1000)
    cax = make_axes_locatable(ax).append_axes("right", size="30%", pad=0.85)
    amp_bar(cax, get_amps(res_t))
    cax.set_ylabel(r"$\Delta_r$ $A_{\mathrm{T}}$ (JJA $-$ DJF) (%)")
    cax.text(0.05, 0.1, "(b)", transform=cax.transAxes, va="top")

    ax = plt.subplot(3, 1, 2)
    plot_error(ax, res_pr, months)
    ax.grid(True); ax.set_xticks(months); ax.set_xticklabels([])
    ax.axhline(0, linestyle="--", linewidth=1.5, color="black")
    plt.ylabel(r"$\Delta_r$ Pr (%)")
    plt.ylim([-55, 55])
    ax.text(0.01, 0.95, subpanel[1], transform=ax.transAxes, va="top", zorder=1000)
    cax = make_axes_locatable(ax).append_axes("right", size="30%", pad=0.75)
    amp_bar(cax, get_amps(res_pr))
    cax.set_ylabel(r"$\Delta_r$ $A_{\mathrm{pr}}$ (JJA $-$ DJF) (%)")
    cax.text(0.05, 0.1, "(d)", transform=cax.transAxes, va="top")

    ax = plt.subplot(3, 1, 3)
    plot_error(ax, res_wind, months)
    ax.axhline(0, linestyle="--", linewidth=1.5, color="black")
    ax.grid(True); ax.set_xticks(months); ax.set_xticklabels(month_labels)
    plt.xticks(months)
    plt.ylabel(r"$\Delta_r$ Wspd (%)")
    plt.xlabel("Time (Months)")
    plt.ylim([-15, 15])
    ax.text(0.01, 0.95, subpanel[2], transform=ax.transAxes, va="top", zorder=1000)
    cax = make_axes_locatable(ax).append_axes("right", size="30%", pad=0.75)
    amp_bar(cax, get_amps(res_wind))
    cax.set_ylabel(r"$\Delta_r$ $A_{\mathrm{Wspd}}$ (JJA $-$ DJF) (%)")
    cax.text(0.05, 0.1, "(f)", transform=cax.transAxes, va="top")

    fig_path = os.path.join(PLOT_DIR, "Regional_averages_relative_CF.png")
    plt.savefig(fig_path, format="png", bbox_inches="tight", transparent=False, dpi=200)
    plt.close()
    return fig_path


def main():
    res_t = get_monthly_climatologies("t")
    res_pr = get_monthly_climatologies("pr")
    res_wind = get_monthly_climatologies("wind")

    fig_path = make_figure(res_t, res_pr, res_wind)
    print(f"saved figure: {fig_path}")

    npz = {}
    for name, res in (("t", res_t), ("pr", res_pr), ("wind", res_wind)):
        for block in ("historical", "rcp45", "rcp85", "delta45", "delta85"):
            for key, val in res[block].items():
                npz[f"{name}_{block}_{key}"] = np.asarray(val)
    npz_path = os.path.join(PLOT_DIR, "cf_regional_means.npz")
    np.savez(npz_path, **npz)
    print(f"saved numbers: {npz_path}")

    month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    for name, res in (("T (K)", res_t), ("Pr (%)", res_pr), ("Wspd (%)", res_wind)):
        print(f"\n=== {name}  delta85 (future - historical) ===")
        print("month : " + "  ".join(month_labels))
        print("land  : " + np.array2string(res["delta85"]["land"], precision=2))
        print("ocean : " + np.array2string(res["delta85"]["ocean"], precision=2))


if __name__ == "__main__":
    main()
