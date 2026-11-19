#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_lai.py -- the Tonzi site's leaf area by PFT against what was measured, the mean seasonal cycle.

* The model: each PFT's daily leaf area index (mean_lai_pft_py), the oaks' and the grass's.
* MODIS MCD15A3H (4-day, 500 m; the ORNL DAAC subset around the tower): the median of the 3 x 3 pixels
  around the tower, good retrievals only (raw value at most 100, i.e. LAI 0-10, and QC bit 0 clear: the
  main radiative-transfer algorithm, not the back-up). A 500 m pixel holds oaks, the grass between them
  and the grass beneath them.
* The clip LAI of the AmeriFlux BADM (direct harvests) of the Tonzi understory (LAI_U), the grass beneath
  and between the oaks, read with example 02's Tonzi reader.

The figure is the mean seasonal cycle over the scored years (2011 on: 2010 starts from the census, in
leaf); the observations are each record of those years, on its day of the year.

Usage: python plot_lai.py --modis <folder> --badm <folder> [--model output/eval] [--out lai.png]
  --modis: the folder of us_california_tonzi_ranch_mcd15a3h_{Lai_500m,FparLai_QC}.csv
  --badm:  the folder of AMF_US-Ton_BASE-BADM_27-5.zip
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
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "example02_canopy_phenology", "tonzi_stand"))
import evaluation as ev  # noqa: E402  (the example's colours)

OAK, GRASS = "#5a7d2a", "#e3b23c"
MODIS_GRID = 17                       # the ORNL subset is 17 x 17 pixels, the tower in the middle one
MODIS_SITE = "us_california_tonzi_ranch_mcd15a3h"
SEASONS = {"grass at its peak (March-April)": (3, 4), "oaks alone (July-August)": (7, 8),
           "grass filling (December-January)": (12, 1)}


def read_model(prefix):
    """The model's daily leaf area index by PFT: columns oak (PFT 1) and grass (PFT 2)."""
    with Dataset(prefix + "-D.nc") as ds:
        t = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:]))
        lai = np.asarray(ds["mean_lai_pft_py"][:], float)[:, 0, :]
    return pd.DataFrame({"oak": lai[:, 0], "grass": lai[:, 1]}, index=t)


def read_modis(folder):
    """MCD15A3H LAI, the median of the good retrievals among the 3 x 3 pixels around the tower."""
    lai = pd.read_csv(os.path.join(folder, f"{MODIS_SITE}_Lai_500m.csv"))
    qc = pd.read_csv(os.path.join(folder, f"{MODIS_SITE}_FparLai_QC.csv"))
    mid = MODIS_GRID // 2
    cols = [f"p{r * MODIS_GRID + c}" for r in (mid - 1, mid, mid + 1) for c in (mid - 1, mid, mid + 1)]
    raw, flag = lai[cols].to_numpy(float), qc[cols].to_numpy(float).astype(int)
    good = np.where((raw <= 100) & (flag & 1 == 0), raw * 0.1, np.nan)
    n = np.isfinite(good).sum(axis=1)
    med = np.full(len(good), np.nan)
    med[n > 0] = np.nanmedian(good[n > 0], axis=1)
    return pd.Series(med, index=pd.to_datetime(lai["calendar_date"])).dropna()


def read_clips(folder):
    """The BADM's clip LAI: {label: series}."""
    import prepare_observations as po                   # example 02's BADM reader
    out = {}
    for site, var, label in (("US-Ton", "LAI_U", "Tonzi understory"),):
        groups = po.badm(os.path.join(folder, f"AMF_{site}_BASE-BADM_27-5.zip"))
        rows = [(pd.Timestamp(v["LAI_DATE"][:8]), float(v[var])) for kind, v in groups.values()
                if kind == "GRP_LAI" and var in v and len(v.get("LAI_DATE", "")) >= 8]
        out[label] = pd.Series(dict(rows)).sort_index()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--modis", required=True, help="the folder of the MCD15A3H subset around the tower")
    ap.add_argument("--badm", required=True, help="the folder of the BASE-BADM zips")
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval"), help="the run's output prefix")
    ap.add_argument("--out", default=os.path.join(HERE, "lai.png"))
    ap.add_argument("--stats", default=os.path.join(HERE, "output", "lai_stats.json"))
    args = ap.parse_args(argv)

    model = read_model(args.model)
    later = sorted(y for y in set(model.index.year) if y > model.index[0].year)
    model = model[model.index.year.isin(later)]
    modis = read_modis(args.modis)
    modis = modis[modis.index.year.isin(later)]
    clips = {k: v[v.index.year.isin(later)] for k, v in read_clips(args.badm).items()}
    total = model["oak"] + model["grass"]

    def by_months(s, months):                          # None (JSON null) where nothing was measured
        mean = float(s[s.index.month.isin(months)].mean())
        return round(mean, 3) if np.isfinite(mean) else None
    stats = {"years": later}
    for name, months in SEASONS.items():
        stats[name] = {"model_oak": by_months(model["oak"], months), "model_grass": by_months(model["grass"], months),
                       "model_total": by_months(total, months), "modis": by_months(modis, months),
                       **{f"clips, {k}": by_months(v, months) for k, v in clips.items()}}
    doy = lambda s: s.groupby(s.index.dayofyear).mean()  # noqa: E731
    clim_o, clim_g = doy(model["oak"]), doy(model["grass"])

    plt.rcParams.update({"axes.edgecolor": ev.MUTED, "axes.labelcolor": ev.INK, "xtick.color": ev.MUTED,
                         "ytick.color": ev.MUTED, "axes.grid": True, "grid.color": ev.GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.stackplot(clim_g.index, clim_g.values, clim_o.values, colors=[GRASS, OAK], alpha=0.55, lw=0,
                 labels=["ED-MORE, grass (PFT 2)", "ED-MORE, oaks (PFT 1)"])
    ax.plot(clim_o.index, clim_o.values + clim_g.values, color=ev.MODEL, lw=1.6, label="ED-MORE, total")
    ax.plot(modis.index.dayofyear, modis.values, ls="none", marker="o", ms=2.6, color=ev.MUTED,
            label="MODIS MCD15A3H (3 × 3 median)")
    for label, s in clips.items():
        ax.plot(s.index.dayofyear, s.values, ls="none", marker="s", ms=4, mec=ev.INK, mfc=ev.INK,
                label=f"clips, {label} (BADM)")
    ticks = pd.date_range("2021-01-01", periods=12, freq="MS")
    ax.set(xlim=(1, 366), xticks=ticks.dayofyear, xticklabels=list(ev.MONTHS), ylim=(0, None),
           xlabel="month", ylabel="leaf area index [m² m⁻²]")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.set_title(f"Tonzi Ranch: leaf area by PFT, the mean seasonal cycle {later[0]}-{later[-1]}",
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
