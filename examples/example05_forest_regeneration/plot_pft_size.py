#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Plot the per-PFT and per-DBH-class axes of an ED-MORE diagnostic file.

ED-MORE aggregates every cohort quantity onto five axes; two of them -- PFT and DBH size class --
describe *who* the stand is made of rather than how much of it there is, and they are what a
demographic model exists to produce. This is their reference consumer.

It draws six panels:

    1. aboveground biomass by PFT, over time        -- successional composition
    2. leaf area index by PFT, over time            -- who holds the canopy
    3. stem density by DBH class, first vs last     -- the size distribution (the inverse-J)
    4. aboveground biomass by DBH class, first/last -- where the carbon actually is
    5. biomass by DBH class over time, stacked      -- how the structure moves
    6. growth and mortality by DBH class            -- the size dependence of the vital rates

and prints the two closure identities the aggregation contract promises
(docs/science/diagnostics.md eq. (3)):

    sum over PFTs of mean_agb_pft_py  ==  mean_agb_py
    sum over classes of mean_agb_size_py  ==  mean_agb_py

Those are worth printing rather than assuming. They are the statement that the two axes partition
the stand, and a half-open/closed class-edge error or a wrong reduction weight breaks them while
every individual number still looks plausible.

EMPTY SETS ARE NOT ZEROS. A *mean* over an empty PFT or size class (mean_mort_rate_size_py, mean_dbh_pft_py)
is the fill value, while a *sum* over one (mean_agb_size_py, mean_nplant_size_py) is a true zero -- otherwise the
identities above would break the moment a class emptied. netCDF4 returns masked arrays, so the fill
is masked out here rather than plotted as 1e37; a reader that strips the mask with np.array() gets
the raw fill and a ruined axis.

Usage:
    python plot_pft_size.py DIAG.nc [-o OUT.png]

DIAG.nc is an aggregation file from any tier that carries the cohort axes -- monthly or annual.
(The annual stream forbids the raw cohort and patch axes, because a window longer than a month
straddles the disturbance restructuring, but the PFT and size axes are reductions and are fine.)

Requires numpy, matplotlib and netCDF4.
"""
import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")            # headless: write a PNG, no display needed
import matplotlib.pyplot as plt
from netCDF4 import Dataset, MFDataset

#----- Classic ED / Moorcroft et al. 2001 PFT colours, shared with plot_site_timeseries.py. -----#
PFT_COLORS = ["green", "blue", "magenta"]


def pft_color(i):
    """Colour for PFT index i (0-based): the Moorcroft scheme, cycling for any extra PFTs."""
    return PFT_COLORS[i] if i < len(PFT_COLORS) else f"C{i}"


def open_records(path):
    """Open one diagnostic file, or a glob over a tier's per-chunk files, as one record series.

    The monthly and daily tiers write one file per chunk, so a run's series is normally several
    files; the annual tier writes one. MFDataset needs the classic (netCDF-3) format, so fall back
    to concatenating by hand when it refuses -- the reader should not care which format the run
    happened to write.
    """
    paths = sorted(glob.glob(path)) if any(c in path for c in "*?[") else [path]
    if len(paths) == 1:
        return Dataset(paths[0]), paths
    try:
        return MFDataset(paths), paths
    except Exception:
        return _Concat(paths), paths


class _Concat:
    """Minimal stand-in for MFDataset: concatenates each variable along its record axis."""

    def __init__(self, paths):
        self._ds = [Dataset(p) for p in paths]
        self.variables = {}
        self.dimensions = self._ds[0].dimensions
        for name in self._ds[0].variables:
            if all(name in d.variables for d in self._ds):
                self.variables[name] = _ConcatVar(name, self._ds)

    def close(self):
        for d in self._ds:
            d.close()

    def __getitem__(self, k):
        return self.variables[k]


class _ConcatVar:
    def __init__(self, name, ds):
        self._name, self._ds = name, ds
        v0 = ds[0].variables[name]
        self.dimensions, self.units = v0.dimensions, getattr(v0, "units", "")
        self.long_name = getattr(v0, "long_name", name)

    def __getitem__(self, idx):
        parts = [d.variables[self._name][:] for d in self._ds]
        #----- A variable with no record axis (the class edges, the PFT list) is the same in every  !
        #      file; concatenating it would fabricate a longer axis. -------------------------------#
        if "time" not in self.dimensions:
            return parts[0][idx]
        return np.ma.concatenate(parts, axis=0)[idx]


def get(ds, name):
    """Variable by name, on the run's one polygon, or None when the run's config did not enable it."""
    if name not in ds.variables:
        return None
    var = ds.variables[name]
    x = var[:]
    return x.take(0, axis=var.dimensions.index("polygon")) if "polygon" in var.dimensions else x


def class_labels(lower, upper):
    """Readable DBH-class labels from the file's own edge coordinates, never hard-coded."""
    return [f"{lo:g}–{hi:g}" for lo, hi in zip(lower, upper)]


def report_closure(ds, per_axis, site_name, axis_name, label):
    """Print sum-over-axis vs the site total, which the aggregation contract says must agree.

    Returns the worst relative discrepancy, or None when the run did not write both sides.
    """
    site = get(ds, site_name)
    if per_axis is None or site is None:
        print(f"   {label:<34} (not in this file -- skipped)")
        return None
    tot = per_axis.sum(axis=1)
    denom = np.maximum(np.abs(site), 1e-30)
    rel = np.abs(tot - site) / denom
    worst = float(np.max(rel))
    mark = "OK " if worst < 1e-10 else "!! "
    print(f"   {mark}{label:<32} worst relative gap = {worst:.3e}   "
          f"({axis_name} sum {tot[-1]:.6g} vs site {site[-1]:.6g})")
    return worst


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", help="diagnostic netCDF file, or a glob over a tier's files")
    ap.add_argument("-o", "--out", default="pft_size.png", help="output PNG")
    args = ap.parse_args()

    ds, paths = open_records(args.infile)
    print(f"reading {len(paths)} file(s): {os.path.basename(paths[0])}"
          + (f" .. {os.path.basename(paths[-1])}" if len(paths) > 1 else ""))

    #----- A decimal calendar year for the plot axis: `time` counts days since the run's start. -----#
    t = get(ds, "year") + (get(ds, "month") - 1 + (get(ds, "day") - 1) / 31.0) / 12.0
    pft_ids = get(ds, "pft")
    lower, upper = get(ds, "dbh_lower"), get(ds, "dbh_upper")
    if lower is None or upper is None:
        raise SystemExit("this file has no dbh_class axis -- enable the size axis in [output]")
    labels = class_labels(lower, upper)
    npft = len(pft_ids) if pft_ids is not None else 0

    mean_agb_pft_py = get(ds, "mean_agb_pft_py")
    mean_lai_pft_py = get(ds, "mean_lai_pft_py")
    mean_agb_size_py = get(ds, "mean_agb_size_py")
    npl_size = get(ds, "mean_nplant_size_py")
    grow_size = get(ds, "mean_agb_growth_size_py")
    mort_size = get(ds, "mean_mort_rate_size_py")

    #----- The identities first: a figure drawn from axes that do not partition the stand is a     !
    #      picture of a bug, so say whether they do before drawing anything. ----------------------#
    print(f"\nclosure of the reduction axes ({len(t)} records, {npft} PFT(s), {len(labels)} DBH classes):")
    report_closure(ds, mean_agb_pft_py, "mean_agb_py", "PFT", "sum_pft mean_agb_pft_py == mean_agb_py")
    report_closure(ds, mean_agb_size_py, "mean_agb_py", "class", "sum_class mean_agb_size_py == mean_agb_py")
    report_closure(ds, get(ds, "mean_nplant_size_py"), "mean_nplant_py", "class",
                   "sum_class mean_nplant_size_py == mean_nplant_py")
    report_closure(ds, get(ds, "mean_nplant_pft_py"), "mean_nplant_py", "PFT",
                   "sum_pft mean_nplant_pft_py == mean_nplant_py")
    if npft == 1:
        print("\n   NOTE: this run has ONE active PFT, so panels 1-2 are a single flat trace.\n"
              "         They are drawn anyway -- the composition story is what they are for, and a\n"
              "         one-PFT run is the degenerate case of it, not a different figure.")

    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    fig.suptitle("ED-MORE stand composition and size structure", fontsize=14)

    #----- 1-2. Composition by PFT. ------------------------------------------------------------#
    for ax, dat, name in ((axes[0, 0], mean_agb_pft_py, "aboveground biomass"),
                          (axes[0, 1], mean_lai_pft_py, "leaf area index")):
        if dat is None:
            ax.set_axis_off()
            continue
        for p in range(dat.shape[1]):
            ax.plot(t, dat[:, p], color=pft_color(p), lw=2, label=f"PFT {int(pft_ids[p])}")
        ax.set_xlabel("year")
        ax.set_ylabel(f"{name} [{ds.variables[('mean_agb_pft_py' if dat is mean_agb_pft_py else 'mean_lai_pft_py')].units}]")
        ax.set_title(f"{name} by PFT")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    #----- 3-4. The size distribution, first record against last. The stem-density panel is on a   !
    #      log y-axis because a healthy stand's size distribution spans orders of magnitude -- the  !
    #      inverse-J is the signal, and a linear axis shows only the smallest class. ---------------#
    x = np.arange(len(labels))
    for ax, dat, name, logy in ((axes[0, 2], npl_size, "stem density", True),
                                (axes[1, 0], mean_agb_size_py, "aboveground biomass", False)):
        if dat is None:
            ax.set_axis_off()
            continue
        ax.bar(x - 0.2, dat[0], width=0.4, color="0.7", label=f"{t[0]:.0f}")
        ax.bar(x + 0.2, dat[-1], width=0.4, color="darkgreen", label=f"{t[-1]:.0f}")
        if logy:
            ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("DBH class [cm]")
        ax.set_ylabel(f"{name} [{ds.variables['mean_nplant_size_py' if logy else 'mean_agb_size_py'].units}]")
        ax.set_title(f"{name} by size class")
        ax.legend(fontsize=8, title="year")
        ax.grid(alpha=0.3, axis="y")

    #----- 5. Where the biomass sits, over the run. Stacked, so the panel answers "is the stand     !
    #      moving carbon into bigger trees" rather than six separate line questions. ---------------#
    ax = axes[1, 1]
    if mean_agb_size_py is not None:
        ax.stackplot(t, np.ma.filled(mean_agb_size_py, 0.0).T,
                     labels=labels, colors=plt.cm.viridis(np.linspace(0.15, 0.95, len(labels))))
        ax.set_xlabel("year")
        ax.set_ylabel(f"aboveground biomass [{ds.variables['mean_agb_size_py'].units}]")
        ax.set_title("biomass by size class over time")
        ax.legend(fontsize=7, title="DBH [cm]", loc="upper left", ncol=2)
        ax.grid(alpha=0.3)

    #----- 6. The vital rates against size -- the demographic content of the size axis. Growth is   !
    #      a per-ground-area SUM and mortality an nplant-weighted MEAN, so an empty class gives 0   !
    #      for one and the fill value for the other. Both are plotted from the masked array, which  !
    #      simply leaves the empty class blank on the mortality axis. ------------------------------#
    ax = axes[1, 2]
    n_empty = 0
    if grow_size is not None:
        ax.bar(x, grow_size[-1], color="seagreen", alpha=0.8, label="AGB growth")
        ax.set_ylabel(f"AGB growth [{ds.variables['mean_agb_growth_size_py'].units}]", color="seagreen")
        ax.tick_params(axis="y", labelcolor="seagreen")
    if mort_size is not None:
        ax2 = ax.twinx()
        ax2.plot(x, mort_size[-1], "o-", color="firebrick", lw=2, label="mortality rate")
        ax2.set_ylabel(f"mortality rate [{ds.variables['mean_mort_rate_size_py'].units}]", color="firebrick")
        ax2.tick_params(axis="y", labelcolor="firebrick")
        n_empty = int(np.ma.getmaskarray(mort_size[-1]).sum())
    #----- The empty-class note goes under the axis, not inside it: the tallest growth bar is in    !
    #      the largest occupied class, which is exactly where an in-axes annotation would land. ----#
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    note = "DBH class [cm]"
    if mort_size is not None and n_empty:
        note += f"      ({n_empty} empty class: a MEAN there is fill, a SUM is a true 0)"
    ax.set_xlabel(note, fontsize=8)
    ax.set_title(f"vital rates by size class ({t[-1]:.0f})")
    ax.grid(alpha=0.3, axis="y")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out, dpi=130)
    print(f"\nwrote {args.out}")
    ds.close()


if __name__ == "__main__":
    main()
