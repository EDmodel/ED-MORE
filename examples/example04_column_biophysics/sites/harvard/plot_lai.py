#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_lai.py -- the Harvard site's leaf area by PFT against what was measured.

* The model: each PFT's daily leaf area index (mean_lai_pft_py), the broadleaves' and the conifers'.
* MODIS MCD15A3H (4-day, 500 m; the ORNL DAAC subset around Harvard Forest): the median of the 3 x 3
  pixels around the tower, good retrievals only (raw value at most 100, i.e. LAI 0-10, and QC bit 0
  clear: the main radiative-transfer algorithm, not the back-up).
* LAI-2000 surveys (Harvard Forest HF069, the 33 plots around the EMS tower, 1 km south-east of NEON's
  tower): the plant area index of the whole canopy and its deciduous part (the leaf-off reading
  subtracted). The instrument sees stems and branches too and does not correct for clumping, so its
  level is a lower bound on the leaves and its timing is the reference.

The figure is the mean seasonal cycle over the scored years. The statistics also hold each year's
leaf-out and leaf-fall dates: where the broadleaf (model), deciduous (LAI-2000, between the surveys) or
total (MODIS, a running median of five retrievals) leaf area crosses halfway between its leaf-off and
its summer level.

Usage: python plot_lai.py --data <harvard data> [--model output/cal] [--out lai.png]
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
import evaluation as ev  # noqa: E402  (the example's colours)

BROAD, CONIFER = "#1baf7a", "#eda100"
MODIS_GRID = 17                       # the ORNL subset is 17 x 17 pixels, the tower in the middle one
LEAF_OFF, SUMMER = (1, 2, 3), (7, 8)  # the months a year's leaf-off and summer levels are taken over


def read_model(prefix):
    """The model's daily leaf area index by PFT: columns broad (PFT 1) and conifer (PFT 2)."""
    with Dataset(prefix + "-D.nc") as ds:
        t = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:]))
        lai = np.asarray(ds["mean_lai_pft_py"][:], float)[:, 0, :]
    return pd.DataFrame({"broad": lai[:, 0], "conifer": lai[:, 1]}, index=t)


def read_modis(data):
    """MCD15A3H LAI, the median of the good retrievals among the 3 x 3 pixels around the tower."""
    lai = pd.read_csv(os.path.join(data, "modis", "mcd15a3h_Lai_500m.csv"))
    qc = pd.read_csv(os.path.join(data, "modis", "mcd15a3h_FparLai_QC.csv"))
    mid = MODIS_GRID // 2
    cols = [f"p{r * MODIS_GRID + c}" for r in (mid - 1, mid, mid + 1) for c in (mid - 1, mid, mid + 1)]
    raw, flag = lai[cols].to_numpy(float), qc[cols].to_numpy(float).astype(int)
    good = np.where((raw <= 100) & (flag & 1 == 0), raw * 0.1, np.nan)
    n = np.isfinite(good).sum(axis=1)
    med = np.full(len(good), np.nan)
    med[n > 0] = np.nanmedian(good[n > 0], axis=1)
    return pd.Series(med, index=pd.to_datetime(lai["calendar_date"])).dropna()


def read_lai2000(data):
    d = pd.read_csv(os.path.join(data, "hf_archive", "hf069", "hf069-02-lai-site.csv"), parse_dates=["date"])
    return d[d["site"] == "ems"].set_index("date")[["all_mean", "decid_mean"]]


def crossing(s, level, rising):
    """The first time the series passes `level` (upward if rising), interpolated between its records."""
    v, t = s.to_numpy(float), s.index.to_numpy("datetime64[ns]").astype("int64").astype(float)
    hit = np.flatnonzero((v[:-1] <= level) & (v[1:] > level) if rising else (v[:-1] >= level) & (v[1:] < level))
    if not len(hit):
        return None
    i = hit[0]
    return pd.Timestamp(int(t[i] + (level - v[i]) / (v[i + 1] - v[i]) * (t[i + 1] - t[i]))).normalize()


def halfway_dates(s, years):
    """Each year's leaf-out and leaf-fall: where the spring rise and the autumn fall pass halfway
    between the year's leaf-off level (January-March, else its lowest record before June) and its
    summer (July-August) mean."""
    out = {}
    for y in years:
        g = s[str(y)].dropna()
        off = g[g.index.month.isin(LEAF_OFF)]
        lo = off.mean() if len(off) else g[g.index.month < 6].min()
        hi = g[g.index.month.isin(SUMMER)].mean()
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi - lo < 1.0:
            continue
        mid = 0.5 * (lo + hi)
        up, down = crossing(g[g.index.month <= 7], mid, True), crossing(g[g.index.month >= 8], mid, False)
        out[int(y)] = {"leaf_out": str(up.date()) if up is not None else None,
                       "leaf_fall": str(down.date()) if down is not None else None}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", default=os.environ.get("HARVARD_DATA"), required=os.environ.get("HARVARD_DATA") is None)
    ap.add_argument("--model", default=os.path.join(HERE, "output", "cal"), help="the run's output prefix")
    ap.add_argument("--out", default=os.path.join(HERE, "lai.png"))
    ap.add_argument("--stats", default=os.path.join(HERE, "output", "lai_stats.json"))
    args = ap.parse_args(argv)

    model = read_model(args.model)
    span = (model.index[0], model.index[-1])
    modis = read_modis(args.data)[span[0]:span[1]]
    survey = read_lai2000(args.data)[span[0]:span[1]]
    total = model["broad"] + model["conifer"]
    years = sorted(set(model.index.year))
    later = [y for y in years if y > years[0]]          # the first year starts from the census, in leaf

    #----- the numbers -----------------------------------------------------------------------------#
    def by_months(s, months, yrs):
        return float(s[s.index.month.isin(months) & s.index.year.isin(yrs)].mean())
    stats = {"years_compared": later,
             "summer_jul_aug": {"model_total": by_months(total, SUMMER, later),
                                "model_broad": by_months(model["broad"], SUMMER, later),
                                "model_conifer": by_months(model["conifer"], SUMMER, later),
                                "modis": by_months(modis, SUMMER, later),
                                "lai2000_total": by_months(survey["all_mean"], SUMMER, later),
                                "lai2000_deciduous": by_months(survey["decid_mean"], SUMMER, later)},
             "leaf_off_jan_mar": {"model_total": by_months(total, LEAF_OFF, later),
                                  "modis": by_months(modis, LEAF_OFF, later),
                                  "lai2000_total": by_months(survey["all_mean"], (3, 4), later)},
             "halfway_model_broad": halfway_dates(model["broad"], later),
             "halfway_lai2000_deciduous": halfway_dates(survey["decid_mean"], later),
             "halfway_modis": halfway_dates(modis.rolling(5, center=True).median(), later)}
    for key in ("leaf_out", "leaf_fall"):              # the model's dates less the surveys' [days]
        d = [(pd.Timestamp(stats["halfway_model_broad"][y][key]) - pd.Timestamp(v[key])).days
             for y, v in stats["halfway_lai2000_deciduous"].items()
             if v.get(key) and stats["halfway_model_broad"].get(y, {}).get(key)]
        stats[f"{key}_model_minus_lai2000_days"] = d
    doy = lambda s: s.groupby(s.index.dayofyear).mean()  # noqa: E731
    clim_m = model[model.index.year.isin(later)]
    clim_b, clim_c = doy(clim_m["broad"]), doy(clim_m["conifer"])
    clim_modis = modis[modis.index.year.isin(later)]
    clim_survey = survey[survey.index.year.isin(later)]

    #----- the figure ------------------------------------------------------------------------------#
    plt.rcParams.update({"axes.edgecolor": ev.MUTED, "axes.labelcolor": ev.INK, "xtick.color": ev.MUTED,
                         "ytick.color": ev.MUTED, "axes.grid": True, "grid.color": ev.GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.stackplot(clim_c.index, clim_c.values, clim_b.values, colors=[CONIFER, BROAD], alpha=0.55, lw=0,
                 labels=["ED-MORE, conifers (PFT 2)", "ED-MORE, broadleaves (PFT 1)"])
    ax.plot(clim_b.index, clim_b.values + clim_c.values, color=ev.MODEL, lw=1.6, label="ED-MORE, total")
    ax.plot(clim_modis.index.dayofyear, clim_modis.values, ls="none", marker="o", ms=2.6, color=ev.MUTED,
            label="MODIS MCD15A3H (3 × 3 median)")
    ax.plot(clim_survey.index.dayofyear, clim_survey["all_mean"], ls="none", marker="s", ms=4, color=ev.INK,
            label="LAI-2000, whole canopy (HF069)")
    ax.plot(clim_survey.index.dayofyear, clim_survey["decid_mean"], ls="none", marker="s", ms=4, mfc="white",
            mec=ev.INK, label="LAI-2000, deciduous part")
    ticks = pd.date_range("2021-01-01", periods=12, freq="MS")
    ax.set(xlim=(1, 366), xticks=ticks.dayofyear, xticklabels=list(ev.MONTHS), ylim=(0, None),
           xlabel="month", ylabel="leaf area index [m² m⁻²]")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_title(f"Harvard Forest: leaf area by PFT, the mean seasonal cycle {later[0]}-{later[-1]}",
                 color=ev.INK, fontsize=11, loc="left")
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    os.makedirs(os.path.dirname(os.path.abspath(args.stats)), exist_ok=True)
    with open(args.stats, "w") as fh:
        json.dump(stats, fh, indent=2)
        fh.write("\n")
    print(f"figure: {args.out}")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
