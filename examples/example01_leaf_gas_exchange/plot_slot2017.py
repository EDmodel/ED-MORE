# SPDX-License-Identifier: Apache-2.0
"""Plot the Slot & Winter (2017) reproduction as one figure, from the CSVs that
reproduce_slot2017.py writes:

  <prefix>_aci.csv        ci, ac, aj, anet   (net rates for F. insipida; anet = min(ac, aj))
  <prefix>_<species>.csv  tleaf_c, vcmax, jmax, rlight, gs, anet

Left: the A-Ci demand curve. Right: five stacked leaf-temperature panels, one colour per species.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

COLORS = ["tab:green", "tab:blue", "tab:orange", "tab:red"]
UMOL = r"[$\mu$mol m$^{-2}$s$^{-1}$]"
TEMP_PANELS = [("vcmax",  r"V$_{cmax}$",  UMOL,                     (0, 400)),
               ("jmax",   r"J$_{max}$",   UMOL,                     (0, 260)),
               ("gs",     r"g$_s$",       r"[mol m$^{-2}$s$^{-1}$]", (0, 0.6)),
               ("anet",   r"A$_{net}$",   UMOL,                     (0, 26)),
               ("rlight", r"R$_{light}$", UMOL,                     (0, 1.6))]


def first_crossing(x, y):
    """First x where y changes sign (linear interpolation), or None."""
    k = np.flatnonzero(np.diff(np.sign(y)))
    if k.size == 0:
        return None
    k = k[0]
    return x[k] - y[k] * (x[k + 1] - x[k]) / (y[k + 1] - y[k])


def plot_aci(ax, d):
    ci = d["ci"]
    ax.plot(ci, d["ac"],   color="tab:red",  lw=1.7, label=r"RuBP carboxylation (A$_c$)")
    ax.plot(ci, d["aj"],   color="tab:blue", lw=1.7, label=r"RuBP regeneration (A$_j$)")
    ax.plot(ci, d["anet"], color="black",    lw=2.6, label=r"limiting rate (A$_{net}$)")
    ax.axhline(0.0, color="0.6", lw=0.8)
    #----- Compensation point, then the Ac/Aj crossover above it (below Gamma* both rates are
    #      negative and cross once more, which is not the transition the figure marks).
    gamma = first_crossing(ci, d["anet"])
    if gamma is not None:
        ax.plot([gamma], [0.0], "o", color="black", ms=4)
        ax.annotate(rf"$\Gamma$ = {gamma:.0f}", (gamma, 0.0), textcoords="offset points",
                    xytext=(8, -14), fontsize=9)
        above = ci > gamma
        xo = first_crossing(ci[above], (d["ac"] - d["aj"])[above])
        if xo is not None:
            ax.plot([xo], [np.interp(xo, ci, d["anet"])], marker="*", color="purple", ms=16,
                    zorder=5)
    ax.set_xlabel(r"intercellular CO$_2$   C$_i$ [$\mu$mol mol$^{-1}$]", fontsize=11)
    ax.set_ylabel(r"net photosynthesis  A [$\mu$mol m$^{-2}$ s$^{-1}$]", fontsize=11)
    ax.set_title("(a)  A–C$_i$ demand curve, F. insipida (paper Fig. 1b)", fontsize=12, loc="left")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9.5, loc="lower right")
    ax.text(0.035, 0.96, "V$_{cmax}$ = 161\nJ$_{max}$ = 238\n(corrected, 25 °C)",
            transform=ax.transAxes, va="top", fontsize=9.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))


def plot_slot2017(prefix, species, out):
    """Draw the A-Ci panel and the temperature panels for `species` (CSV name stems) into `out`."""
    fig = plt.figure(figsize=(14, 8))
    grid = fig.add_gridspec(5, 2, width_ratios=[1.55, 1.0], hspace=0.16, wspace=0.24,
                            left=0.06, right=0.985, top=0.90, bottom=0.09)
    plot_aci(fig.add_subplot(grid[:, 0]),
             np.genfromtxt(f"{prefix}_aci.csv", delimiter=",", names=True))

    data = {sp: np.genfromtxt(f"{prefix}_{sp}.csv", delimiter=",", names=True) for sp in species}
    for i, (col, name, unit, ylim) in enumerate(TEMP_PANELS):
        ax = fig.add_subplot(grid[i, 1])
        for sp, color in zip(species, COLORS):
            ax.plot(data[sp]["tleaf_c"], data[sp][col], color=color, lw=1.8,
                    label=sp.replace("_", ". "))
        ax.set(xlim=(25, 42), ylim=ylim)
        ax.set_ylabel(f"{name}\n{unit}", fontsize=8.5)
        ax.text(0.015, 0.9, f"({'bcdef'[i]})", transform=ax.transAxes, va="top", fontsize=10)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=9, ncol=4, loc="lower left", bbox_to_anchor=(0.0, 1.02),
                      frameon=False, handlelength=1.5, columnspacing=1.0)
        if i < len(TEMP_PANELS) - 1:
            ax.tick_params(labelbottom=False)
    ax.set_xlabel("leaf temperature [°C]   (b–f: paper Fig. 2)", fontsize=11)

    fig.suptitle("Slot & Winter (2017) reproduced with the MEDS leaf model", fontsize=13)
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
