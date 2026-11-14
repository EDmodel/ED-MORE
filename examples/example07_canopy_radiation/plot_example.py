# SPDX-License-Identifier: Apache-2.0
"""Example 07's figure: the modelled PAR at NEON's sensors against what they measured: down the profile,
under clear and overcast skies, with the sun's height, leaf-off and leaf-on.

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
MODEL, DEFAULT, OBS, GRID = "#2a78d6", "#eb6834", "#2b2b2b", "#d8d7d1"
SEASON_COLOR = {"leaf-off": "#eb6834", "leaf-on": "#2a78d6"}
LABEL = {"mod": "model, foliage at random in each layer", "def": "model, default clumping (leaf 0.8, wood 0.5)"}
FLOOR = ["SP1", "SP3", "SP5"]
LEVELS = ["L5", "L4", "L3", "L2", "L1", "floor"]


def with_floor_mean(df):
    """The floor sensors averaged over those reporting in each half hour, the model over the same ones."""
    df = df.copy()
    have = np.isfinite(df[[f"obs_{s}" for s in FLOOR]].to_numpy())
    pairs = [("obs_floor", "obs_{}")] + [(f"{o}_floor{q}", f"{o}_{{}}{q}") for o in ("mod", "def")
                                         for q in ("", "_p10", "_p90")]
    for out_col, column in pairs:
        src = np.where(have, df[[column.format(s) for s in FLOOR]].to_numpy(), np.nan)
        n = np.isfinite(src).sum(axis=1)
        df[out_col] = np.where(n > 0, np.nansum(src, axis=1) / np.maximum(n, 1), np.nan)
    return df


def band(df, sensor, mask, q, which="mod"):
    """A model's patch percentile q (p10, p90) as a ratio of sums, over the same half hours as ratio()."""
    ok = mask & df[f"obs_{sensor}"].notna() & df[f"{which}_{sensor}_{q}"].notna() & df["par_top"].notna()
    return df.loc[ok, f"{which}_{sensor}_{q}"].sum() / df.loc[ok, "par_top"].sum() if ok.sum() >= 8 else np.nan


def ratio(df, sensor, mask, which):
    """Sum of the sensor (observed or modelled) over the masked half hours, over the sum of the top's,
    where both observed and modelled values exist."""
    ok = mask & df[f"obs_{sensor}"].notna() & df[f"mod_{sensor}"].notna() & df["par_top"].notna()
    if ok.sum() < 8:
        return np.nan
    return df.loc[ok, f"{which}_{sensor}"].sum() / df.loc[ok, "par_top"].sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=HERE / "output")
    a = ap.parse_args()
    df = with_floor_mean(pd.read_csv(a.out / "model_vs_sensors.csv", index_col=0, parse_dates=True))
    z = pd.read_csv(a.out / "sensor_heights.csv", index_col=0)["height_m"].to_dict()
    z["floor"] = 0.05
    sun = df["cosz"] > 0.35
    clear, overcast = df["diffuse_frac"] < 0.35, df["diffuse_frac"] > 0.8

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6})
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 8.4), gridspec_kw=dict(hspace=0.42, wspace=0.34))

    #----- (a, b) The profile, leaf-off and leaf-on: the model's patch mean and its 10-90 % range. -----#
    for ax, season in zip(axs[0, :2], ("leaf-off", "leaf-on")):
        m = (df["window"] == season) & sun
        h = [z[s] for s in LEVELS]
        lo = [band(df, s, m, "p10") for s in LEVELS]
        hi = [band(df, s, m, "p90") for s in LEVELS]
        ax.fill_betweenx(h, lo, hi, color=MODEL, alpha=0.15, lw=0, label="its patches' 10-90 % range")
        ax.plot([ratio(df, s, m, "mod") for s in LEVELS], h, "-", color=MODEL, lw=1.8, label=LABEL["mod"])
        ax.plot([ratio(df, s, m, "def") for s in LEVELS], h, "--", color=DEFAULT, lw=1.5, label=LABEL["def"])
        ax.plot([ratio(df, s, m, "obs") for s in LEVELS], h, "o", ms=6.5, color=OBS, label="observed")
        ax.set_xscale("log")
        ax.set_xlim(3e-3, 1.5)
        ax.set_ylim(-1.5, 31)
        ax.axhspan(17, 27.5, color="#8a8a85", alpha=0.08, lw=0)
        ax.text(4e-3, 22, "crowns", color="#8a8a85", fontsize=8, va="center")
        ax.set_title(f"({'ab'[season == 'leaf-on']}) {season}: "
                     + ("1-25 April 2019, conifers only" if season == "leaf-off" else "July 2019, full leaves"),
                     loc="left", fontsize=9)
        ax.set_xlabel("PAR / PAR above the canopy")
        if season == "leaf-off":
            ax.set_ylabel("height [m]")
            ax.legend(loc="upper left", frameon=False, fontsize=7)

    #----- (c) Clear against overcast. ----------------------------------------------------------------#
    ax = axs[0, 2]
    ys, labels, y = [], [], 0.0
    for season in ("leaf-off", "leaf-on"):
        for s, sl in (("L3", "16 m"), ("L2", "5 m"), ("floor", "floor")):
            for kname, km in (("clear", clear), ("overcast", overcast)):
                m = (df["window"] == season) & sun & km
                o, mo, de = ratio(df, s, m, "obs"), ratio(df, s, m, "mod"), ratio(df, s, m, "def")
                ax.plot([min(o, mo, de), max(o, mo, de)], [y, y], "-", color=GRID, lw=1.2, zorder=1)
                ax.plot(de, y, "s", ms=5, color=DEFAULT, zorder=2, label="default clumping" if not ys else None)
                ax.plot(mo, y, "D", ms=5.5, color=MODEL, zorder=3, label="foliage at random" if not ys else None)
                ax.plot(o, y, "o", ms=6, color=OBS, zorder=4, label="observed" if not ys else None)
                ys.append(y)
                labels.append(f"{season.split('-')[1]}, {sl}, {kname}")
                y += 1.0
            y += 0.5
        y += 0.5
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xticks([0.01, 0.03, 0.1, 0.3, 1.0])
    ax.set_xticklabels(["0.01", "0.03", "0.1", "0.3", "1"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("PAR / PAR above the canopy")
    ax.set_title("(c) clear (diffuse < 0.35) against overcast (> 0.8)", loc="left", fontsize=9)
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)

    #----- (d, e) Clear-sky transmittance against the sun's height. ------------------------------------#
    bins = np.arange(0.3, 1.01, 0.1)
    for ax, season in zip(axs[1, :2], ("leaf-off", "leaf-on")):
        for s, sl, marker in (("L3", "16 m", "^"), ("L2", "5 m", "s"), ("floor", "floor", "o")):
            o, mo, c = [], [], []
            for b0, b1 in zip(bins[:-1], bins[1:]):
                m = (df["window"] == season) & clear & (df["cosz"] >= b0) & (df["cosz"] < b1)
                o.append(ratio(df, s, m, "obs"))
                mo.append(ratio(df, s, m, "mod"))
                c.append(0.5 * (b0 + b1))
            ax.plot(c, mo, "-", color=MODEL, lw=1.6)
            ax.plot(c, o, marker, ms=6, color=OBS, ls="none", label=f"observed, {sl}")
        ax.plot([], [], "-", color=MODEL, lw=1.6, label="model (foliage at random)")
        ax.set_yscale("log")
        ax.set_xlabel("cosine of the sun's zenith angle (half-hour mean)")
        ax.set_title(f"({'de'[season == 'leaf-on']}) {season}, clear sky: transmittance against the sun",
                     loc="left", fontsize=9)
        if season == "leaf-off":
            ax.set_ylabel("PAR / PAR above the canopy")
        else:
            ax.set_ylim(top=0.8)                                            # room for the legend
            ax.legend(loc="upper left", frameon=False, fontsize=7.5, ncol=2)

    #----- (f) Every half hour at 5 m. --------------------------------------------------------------------#
    ax = axs[1, 2]
    for season, color in SEASON_COLOR.items():
        m = (df["window"] == season) & sun & df["obs_L2"].notna() & df["mod_L2"].notna()
        ax.plot(df.loc[m, "obs_L2"] / df.loc[m, "par_top"], df.loc[m, "mod_L2"] / df.loc[m, "par_top"], "o",
                ms=3, alpha=0.5, color=color, label=season)
    ax.plot([1e-3, 2], [1e-3, 2], "-", color="#8a8a85", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(5e-3, 1.5)
    ax.set_ylim(5e-3, 1.5)
    ax.set_xlabel("observed PAR / PAR above the canopy, 5 m")
    ax.set_ylabel("modelled")
    ax.set_title("(f) every half hour at 5 m", loc="left", fontsize=9)
    ax.legend(loc="upper left", frameon=False, fontsize=8)

    fig.suptitle("ED-MORE's two-stream against NEON's PAR sensors, Harvard Forest, 2019 (prototype)",
                 x=0.06, ha="left", fontsize=11)
    fig.savefig(a.out / "canopy_radiation.png", dpi=150, bbox_inches="tight")
    print("wrote", a.out / "canopy_radiation.png")


if __name__ == "__main__":
    main()
