# SPDX-License-Identifier: Apache-2.0
"""Example 07's figure: the modelled PAR at NEON's sensors against what they measured, through the
year, down the profile, and under clear and overcast skies.

    python plot_example.py [--out output]     (reads output/model_vs_sensors.csv; writes canopy_radiation.png)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402

HERE = Path(__file__).resolve().parent
MODEL, OBS, MUTED = "#2a78d6", "#2b2b2b", "#8a8a85"
FLOOR = ["SP1", "SP3", "SP5"]
MIDDAY_UTC = (15, 19)            # 10:00-14:00 EST


def ratio(df, num, mask):
    """Sum of a sensor over the masked half hours, over the sum of the top's."""
    ok = mask & df[num].notna() & df["par_top"].notna()
    return df.loc[ok, num].sum() / df.loc[ok, "par_top"].sum() if ok.sum() >= 10 else np.nan


def with_floor_mean(df):
    """The floor sensors averaged over those reporting in each half hour, the model over the same
    sensors. December-March are left out: snow covers the line sensors."""
    df = df.copy()
    obs = df[[f"obs_{s}" for s in FLOOR]].to_numpy()
    have = np.isfinite(obs)
    for col, src in (("obs_floor", obs), ("mod_floor", df[[f"mod_{s}" for s in FLOOR]].to_numpy()),
                     ("mod_floor_p10", df[[f"mod_{s}_p10" for s in FLOOR]].to_numpy()),
                     ("mod_floor_p90", df[[f"mod_{s}_p90" for s in FLOOR]].to_numpy())):
        picked = np.where(have, src, np.nan)
        n = np.isfinite(picked).sum(axis=1)
        df[col] = np.where(n > 0, np.nansum(picked, axis=1) / np.maximum(n, 1), np.nan)
    snow = df.index.month.isin([12, 1, 2, 3])
    df.loc[snow, ["obs_floor", "mod_floor", "mod_floor_p10", "mod_floor_p90"]] = np.nan
    return df


def daily_midday(df, sensor):
    hour = df.index.hour
    mid = ((hour >= MIDDAY_UTC[0]) & (hour < MIDDAY_UTC[1]) & df[f"obs_{sensor}"].notna()
           & df[f"mod_{sensor}"].notna() & df["par_top"].notna())
    g = df[mid].groupby(df.index[mid].normalize())
    top = g["par_top"].sum()
    out = pd.DataFrame({k: g[f"{k}_{sensor}"].sum() / top for k in ("obs", "mod")})
    for q in ("p10", "p90"):
        out[q] = g[f"mod_{sensor}_{q}"].sum() / top
    return out[g.size() >= 6]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=HERE / "output")
    a = ap.parse_args()
    df = pd.read_csv(a.out / "model_vs_sensors.csv", index_col=0, parse_dates=True)
    df = with_floor_mean(df)
    heights = {"L5": 29.37, "L4": 23.63, "L3": 16.13, "L2": 5.13, "L1": 0.36, "floor": 0.05}

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6})
    fig = plt.figure(figsize=(12, 8.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05], hspace=0.38, wspace=0.32)

    #----- Top row: midday transmittance day by day, at three heights and the floor. --------------#
    for k, (s, title) in enumerate([("L3", "16 m (in the crowns)"), ("L2", "5 m (below the crowns)"),
                                    ("L1", "0.4 m"), ("floor", "forest floor (3 line sensors)")]):
        ax = fig.add_subplot(gs[0, k])
        d = daily_midday(df, s)
        ax.fill_between(d.index, d["p10"], d["p90"], color=MODEL, alpha=0.15, lw=0, label="model, patch 10-90 %")
        ax.plot(d.index, d["mod"], color=MODEL, lw=1.6, label="model, horizontal mean")
        ax.plot(d.index, d["obs"], "o", ms=2.6, color=OBS, label="observed")
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 1.5)
        ax.set_title(title, loc="left", fontsize=9)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonth=[1, 4, 7, 10]))
        if k == 0:
            ax.set_ylabel("midday PAR / PAR above the canopy")
            ax.legend(loc="lower left", frameon=False, fontsize=8)

    #----- Bottom left: the profile, leaf-off and leaf-on. ----------------------------------------#
    month = df.index.month
    sun = df["cosz"] > 0.35
    seasons = {"leaf-off, April": month == 4, "leaf-on, July-August": (month == 7) | (month == 8)}
    for k, (name, m) in enumerate(seasons.items()):
        ax = fig.add_subplot(gs[1, k])
        z = [heights[s] for s in heights]
        ax.plot([ratio(df, f"mod_{s}", m & sun & df[f"obs_{s}"].notna()) for s in heights], z, "-",
                color=MODEL, lw=1.6, label="model")
        ax.plot([ratio(df, f"obs_{s}", m & sun) for s in heights], z, "o", ms=6, color=OBS, label="observed")
        ax.set_xscale("log")
        ax.set_xlim(5e-3, 1.5)
        ax.axhspan(18, 26, color=MUTED, alpha=0.08, lw=0)
        ax.text(6e-3, 22, "crowns", color=MUTED, fontsize=8, va="center")
        ax.set_title(name, loc="left", fontsize=9)
        ax.set_xlabel("PAR / PAR above the canopy")
        if k == 0:
            ax.set_ylabel("height [m]")
            ax.legend(loc="upper left", frameon=False, fontsize=8)

    #----- Bottom right: clear against overcast, the floor and 16 m, both seasons. ----------------#
    ax = fig.add_subplot(gs[1, 2:])
    skies = {"clear": df["diffuse_frac"] < 0.35, "overcast": df["diffuse_frac"] > 0.8}
    labels, ys, y = [], [], 0.0
    for sname, sm in seasons.items():
        for s, sl in (("L3", "16 m"), ("floor", "floor")):
            for kname, km in skies.items():
                mask = sm & km & sun & df[f"mod_{s}"].notna()
                o, mo = ratio(df, f"obs_{s}", mask), ratio(df, f"mod_{s}", mask & df[f"obs_{s}"].notna())
                ax.plot([o, mo], [y, y], "-", color="#d8d7d1", lw=1.2, zorder=1)
                ax.plot(o, y, "o", ms=6, color=OBS, zorder=2, label="observed" if not ys else None)
                ax.plot(mo, y, "D", ms=5.5, color=MODEL, zorder=3, label="model" if not ys else None)
                labels.append(f"{sname.split(',')[0]}, {sl}, {kname}")
                ys.append(y)
                y += 1.0
            y += 0.6
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("PAR / PAR above the canopy")
    ax.set_title("clear (diffuse < 0.35) against overcast (diffuse > 0.8) skies", loc="left", fontsize=9)
    ax.legend(loc="lower right", frameon=False, fontsize=8)

    fig.suptitle("ED-MORE's two-stream against NEON's PAR sensors, Harvard Forest, 2019 (prototype)",
                 x=0.06, ha="left", fontsize=11)
    fig.savefig(a.out / "canopy_radiation.png", dpi=150, bbox_inches="tight")
    print("wrote", a.out / "canopy_radiation.png")


if __name__ == "__main__":
    main()
