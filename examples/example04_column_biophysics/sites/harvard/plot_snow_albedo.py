#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_snow_albedo.py -- the Harvard site's snow pack and albedo against what was measured, with the
default parameters (no calibration).

* Snow: the model's snow water equivalent (daily) against the Harvard Forest snow pillow (HF155: a
  15 m canopy opening 1.2 km west of the tower, cut so that branches do not intercept the snow; daily
  SWE, 2009-). ED-MORE holds no snow in the canopy, so all of its snow lies on the ground, as the
  pillow's does.
* Albedo: the midday (10-14 local, incoming shortwave above 200 W m-2) ratio of reflected to incoming
  shortwave, the model's at its canopy top against the tower's four-component radiometer at 39 m, both
  on the hours the radiometer measured, by calendar month and against the snow on the ground; and MODIS MCD43A3 (the 3 x 3 500 m pixels around
  the towers, shortwave white- and black-sky, good-quality retrievals) where it has been fetched.

Usage: python plot_snow_albedo.py --data <harvard data> [--model output/eval] [--out snow_albedo.png]
"""
import argparse
import glob
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
import evaluation as ev  # noqa: E402  (the example's colours and paths)
import tower_inputs as ti  # noqa: E402  (on the path evaluation set)

UTC_OFFSET_H = -5
MIDDAY = (10, 14)              # [h local] the hours an albedo is taken over
SW_MIN = 200.0                 # [W m-2] incoming shortwave an hour needs
SNOW_SWE = 5.0                 # [mm] a day with at least this much water in the pack is a snow day


def stamps(ds):
    return pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:], hour=ds["hour"][:]))


def read_model(prefix):
    """The model's daily SWE [mm] and snow depth [m], and its hourly incoming and reflected shortwave,
    on local time."""
    with Dataset(prefix + "-D.nc") as ds:
        t = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:]))
        daily = pd.DataFrame({"swe": np.asarray(ds["mean_water_layer_mass_py"][:], float).squeeze(),
                              "depth": np.asarray(ds["mean_water_layer_depth_py"][:], float).squeeze()}, index=t)
    frames = []
    for f in sorted(glob.glob(prefix + "-F-*.nc")):
        with Dataset(f) as ds:
            frames.append(pd.DataFrame({"sw_in": np.asarray(ds["mean_sw_in_py"][:], float).squeeze(),
                                        "sw_up": np.asarray(ds["mean_sw_up_py"][:], float).squeeze()},
                                       index=stamps(ds) + pd.Timedelta(hours=UTC_OFFSET_H)))
    return daily, pd.concat(frames).sort_index()


def read_tower_albedo(site_toml):
    """The tower's hourly incoming and reflected shortwave where both were measured, on local time."""
    table = ti.read_standard(ti.read_site(site_toml))
    sw = pd.DataFrame({"sw_in": table.values["SWdown"].where(table.measured["SWdown"]),
                       "sw_up": table.values["SW_out"].where(table.measured["SW_out"])})
    hourly = sw.resample("1h").mean().where(sw.resample("1h").count() == 2)
    return hourly.set_axis(hourly.index + pd.Timedelta(hours=UTC_OFFSET_H))


def midday_albedo(hourly):
    """Each day's midday albedo: reflected over incoming shortwave, summed over the bright midday hours."""
    h = hourly[(hourly.index.hour >= MIDDAY[0]) & (hourly.index.hour < MIDDAY[1]) & (hourly["sw_in"] > SW_MIN)]
    h = h.dropna()
    s = h.groupby(h.index.normalize()).agg(sw_in=("sw_in", "sum"), sw_up=("sw_up", "sum"), n=("sw_in", "size"))
    return (s["sw_up"] / s["sw_in"]).where(s["n"] >= 2)


def read_pillow(data):
    d = pd.read_csv(os.path.join(data, "hf_archive", "hf155", "hf155-01-daily.csv"), parse_dates=["date"])
    return d.set_index("date")["swe"]


def read_modis(data):
    """MODIS MCD43A3 shortwave albedo (the mean of white and black sky, good retrievals), the mean of the
    pixels per date; None until the AppEEARS request has been fetched."""
    files = sorted(glob.glob(os.path.join(data, "modis", "mcd43a3_appeears*", "*MCD43A3*results*.csv")))
    if not files:
        return None
    d = pd.read_csv(files[0], parse_dates=["Date"])
    bsa = [c for c in d.columns if c.endswith("Albedo_BSA_shortwave")][0]
    wsa = [c for c in d.columns if c.endswith("Albedo_WSA_shortwave")][0]
    qc = [c for c in d.columns if c.endswith("BRDF_Albedo_Band_Mandatory_Quality_shortwave")][0]
    ok = (d[qc] == 0) & d[bsa].between(0, 1) & d[wsa].between(0, 1)
    d = d[ok]
    return (0.5 * (d[bsa] + d[wsa])).groupby(d["Date"]).mean()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", default=os.environ.get("HARVARD_DATA"), required=os.environ.get("HARVARD_DATA") is None)
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval"), help="the run's output prefix")
    ap.add_argument("--site", default=os.path.join(HERE, "harvard_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "snow_albedo.png"))
    ap.add_argument("--stats", default=os.path.join(HERE, "output", "snow_albedo_stats.json"))
    args = ap.parse_args(argv)

    daily, hourly = read_model(args.model)
    pillow = read_pillow(args.data).reindex(daily.index)
    tower_sw = read_tower_albedo(args.site).reindex(hourly.index)
    measured = tower_sw.notna().all(axis=1)                       # the model on the radiometer's hours
    a_model = midday_albedo(hourly[measured])
    a_tower = midday_albedo(tower_sw[measured]).reindex(a_model.index)
    modis = read_modis(args.data)
    days = pd.DataFrame({"swe_model": daily["swe"], "swe_pillow": pillow, "a_model": a_model, "a_tower": a_tower})

    #----- the winters: November-April, labelled by the year they end in ------------------------------#
    winter = days[days.index.month.isin([11, 12, 1, 2, 3, 4])].copy()
    winter["winter"] = np.where(winter.index.month >= 11, winter.index.year + 1, winter.index.year)
    rows = {}
    for w, g in winter.groupby("winter"):
        if g["swe_pillow"].notna().sum() < 100:
            continue
        both = g.dropna(subset=["swe_model", "swe_pillow"])
        rows[int(w)] = {"peak_pillow": float(both["swe_pillow"].max()), "peak_model": float(both["swe_model"].max()),
                        "snow_days_pillow": int((both["swe_pillow"] >= SNOW_SWE).sum()),
                        "snow_days_model": int((both["swe_model"] >= SNOW_SWE).sum()),
                        "mean_pillow": float(both["swe_pillow"].mean()), "mean_model": float(both["swe_model"].mean())}
    #----- the albedo, by month and by snow -----------------------------------------------------------#
    both = days.dropna(subset=["a_model", "a_tower"])
    monthly = both.groupby(both.index.month)[["a_tower", "a_model"]].mean()
    snowy = both[both["swe_pillow"] >= SNOW_SWE]
    bare_winter = both[(both["swe_pillow"] < SNOW_SWE) & both.index.month.isin([12, 1, 2, 3])]
    summer = both[both.index.month.isin([6, 7, 8])]
    stats = {"winters": rows,
             "albedo_monthly": {int(k): {"tower": float(v["a_tower"]), "model": float(v["a_model"])} for k, v in monthly.iterrows()},
             "albedo_snow_on_ground": {"tower": float(snowy["a_tower"].mean()), "model": float(snowy["a_model"].mean()), "days": int(len(snowy))},
             "albedo_dec_mar_no_snow": {"tower": float(bare_winter["a_tower"].mean()), "model": float(bare_winter["a_model"].mean()),
                                        "days": int(len(bare_winter))},
             "albedo_jja": {"tower": float(summer["a_tower"].mean()), "model": float(summer["a_model"].mean()), "days": int(len(summer))}}
    if modis is not None:
        mm = modis.groupby(modis.index.month).mean()
        stats["albedo_monthly_modis"] = {int(k): float(v) for k, v in mm.items()}

    #----- the figure ---------------------------------------------------------------------------------#
    plt.rcParams.update({"axes.edgecolor": ev.MUTED, "axes.labelcolor": ev.INK, "xtick.color": ev.MUTED,
                         "ytick.color": ev.MUTED, "axes.grid": True, "grid.color": ev.GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig = plt.figure(figsize=(14, 8.2))
    ax_swe = fig.add_subplot(2, 1, 1)
    ax_swe.plot(days.index, days["swe_pillow"], color=ev.TOWER, lw=1.4, label="snow pillow (HF155, canopy opening)")
    ax_swe.plot(days.index, days["swe_model"], color=ev.MODEL, lw=1.4, label="ED-MORE (ground under the canopy)")
    ax_swe.set(ylabel="snow water equivalent [mm]", title="snow on the ground")
    ax_swe.legend(frameon=False, fontsize=9, loc="upper right")
    ax_m = fig.add_subplot(2, 2, 3)
    ax_m.plot(monthly.index, monthly["a_tower"], color=ev.TOWER, lw=2, marker="o", ms=3.5, label="tower (39 m)")
    ax_m.plot(monthly.index, monthly["a_model"], color=ev.MODEL, lw=1.8, marker="o", ms=3.5, label="ED-MORE")
    if modis is not None:
        ax_m.plot(mm.index, mm.values, color=ev.MUTED, lw=1.4, ls="--", marker="s", ms=3, label="MODIS MCD43A3")
    ax_m.set(xlim=(0.5, 12.5), xticks=range(1, 13), xticklabels=list(ev.MONTHS), xlabel="month (local)",
             ylabel="midday shortwave albedo", title="albedo by month")
    ax_m.legend(frameon=False, fontsize=9)
    ax_s = fig.add_subplot(2, 2, 4)
    w = both[both.index.month.isin([12, 1, 2, 3])]
    ax_s.scatter(w["swe_pillow"], w["a_tower"], s=8, color=ev.TOWER, alpha=0.6, label="tower")
    ax_s.scatter(w["swe_model"], w["a_model"], s=8, color=ev.MODEL, alpha=0.6, label="ED-MORE")
    ax_s.set(xlabel="snow water equivalent [mm] (pillow for the tower, the model's own)",
             ylabel="midday shortwave albedo", title="December-March: albedo against the snow pack")
    ax_s.legend(frameon=False, fontsize=9)
    fig.suptitle("Harvard Forest, default parameters: snow and albedo against the snow pillow, NEON's radiometer"
                 + (" and MODIS" if modis is not None else ""), color=ev.INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    os.makedirs(os.path.dirname(os.path.abspath(args.stats)), exist_ok=True)
    with open(args.stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    print(f"figure: {args.out}")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
