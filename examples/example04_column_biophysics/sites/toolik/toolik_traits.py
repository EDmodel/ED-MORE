#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Toolik site's measured soil, for ed_config_eval.toml: what this prints is what its [soil] lists
and [init].soil_temp hold, with the source beside each value.

The soil column (soil_column) is the moist acidic tussock tundra's at NEON's TOOL megapit
(DP1.00096.001), below the pit's Oi (0-6 cm, recognisable dead moss and litter), which is the moss layer
on the soil ([soil].moss_depth). The column's depths are the pit's less those 6 cm.

  * The organic horizons, Oe (6-17 cm in the pit, hemic) and Oa (17-21 cm, sapric), have no bulk density
    or retention in the pit. Each takes the Imnavait Creek pit IMN-1's sample nearest it in decomposition
    (ORNL DAAC 2149, Bakian-Dogaheh et al. 2023; moist acidic tussock tundra 7 km south-east of the
    tower): the hemic Oe IMN-1-2 (7.5-15 cm, fibre 0.39), the sapric Oa IMN-1-4 (22-29 cm, fibre 0.20,
    organic matter 60 % against the Oa's 53 %); its bulk density and porosity, and the van Genuchten curve
    fitted to its paired water content and matric potential with the porosity as theta_sat. Ks is peat's by its decomposition
    (Letts et al. 2000: hemic 2e-6, sapric 1e-7 m/s). The thermal properties are CLM5's blend of the
    mineral soil below (the Bg's texture) with peat, by the organic fraction from the pit's carbon and the
    samples' bulk density (Lawrence & Slater 2008).
  * The mineral horizons, Bg (21-45 cm) and Cg (45-73 cm), take the model's own texture formulas
    (scripts/prepare_input/soil).
  * The permafrost below the thaw (Cgf), from the bottom of the layer holding NEON's deepest sensors that
    thaw (56 cm) to 18 m: the Cg's texture and chemistry at the dry density of the region's loamy permafrost,
    the median of the loamy layers below the thaw depth at Jorgenson's (2016) Toolik sites (Arctic Data
    Center, doi:10.18739/A2CP5C). Their ice lenses push the grains apart: texture and carbon would give
    them about twice that density, and two thirds of their volume is ice. The model has no excess ice, so
    the lower density stands in for it, as pore space the ice fills: a layer that thawed would keep that
    water as liquid. Nothing deeper than the pit was measured at the tower.

The start (soil_start): each layer's temperature on 1-2 July 2019, midsummer, when the spin-up starts, from
NEON's five soil plots (DP1.00041.001, 2-56 cm, each sensor at its recorded depth that day); below the
deepest sensor, a straight line to the deep temperature at its depth ([soil].deep_temp, deep_depth).

Usage:
  python toolik_traits.py --data <toolik data>    # the folder holding neon/, ornl_above/ and organic_layer/
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts", "prepare_input", "soil"))
import soil_horizons as sh              # noqa: E402

KPA_TO_M = 0.10197                      # [m of water head per kPa]
PIT_SOIL_TOP = 0.06                     # [m] the pit's Oi, the moss layer: the soil column starts below it
# [m] the column's layer bottoms below its top: 5-9 cm layers through the horizons and the active layer
# (CALM's end-of-season thaw is 45-53 cm below the moss surface), then growing about 1.4-fold to 18 m.
LAYER_BOTTOMS = [0.05, 0.11, 0.15, 0.23, 0.31, 0.39, 0.48, 0.57, 0.67, 0.80, 1.00, 1.30, 1.70, 2.30, 3.20,
                 4.50, 6.50, 9.50, 13.5, 18.0]
# The organic horizons: the IMN-1 sample nearest each in decomposition, the hemic Oe IMN-1-2 (fibre 0.39),
# the sapric Oa IMN-1-4 (fibre 0.20; organic matter 60 %, the Oa's 53 %), and peat's Ks by decomposition
# (Letts et al. 2000). IMN-1-3 (15-22 cm, fibre 0.33) fits a curve too steep for the frozen and the drying
# peat (n 3.7: its pores near saturation almost to residual).
ORGANIC = {"Oe": dict(samples=[2], ksat=2.0e-6),
           "Oa": dict(samples=[4], ksat=1.0e-7)}
MINERAL = ["Bg", "Cg"]
# [m] the permafrost table in the column: the bottom of the layer holding NEON's 56 cm sensors, the deepest
# that thaw at the tower (to about 1 C in August, 2019-2024); below it the Cg is the ice-rich permafrost
PERMAFROST_TOP = 0.57
# Jorgenson (2016)'s permafrost soil database: the Toolik sites' soil layers, their texture class (L: loamy
# mineral), dry density [g/cm3] and depths [cm], and each site's thaw depth [cm] (999: not reached)
JORGENSON = os.path.join("organic_layer", "jorgenson_thermokarst_soil_db_2009-2011",
                         "20092011_MTJ_Thermokarst_Permafrost_Soil.xlsx")
# The deep temperature: the nearest borehole of the northern foothills, Galbraith Lake (25 km south),
# -4.55 C at 20 m in 2014 after -5.0 in 2000 (Romanovsky, Cable & Kholodov 2015, GeoQuebec 2015, paper
# 479, Fig. 1b, read off the figure), carried to 2021 at its 2000-2014 trend, 0.3 C per decade.
DEEP_TEMP = 273.15 - 4.4                # [K]
DEEP_DEPTH = 20.0 - PIT_SOIL_TOP        # [m] the borehole's 20 m below the ground surface
START = ("2019-07-01", "2019-07-03")    # the days whose mean starts the soil: the spin-up starts in midsummer


def vg_theta_head(h, theta_sat, theta_res, alpha, n):
    """Van Genuchten water content at a suction head h [m, > 0]."""
    return theta_res + (theta_sat - theta_res) / (1.0 + (alpha * h) ** n) ** (1.0 - 1.0 / n)


def fit_imnavait(folder, samples):
    """One van Genuchten curve for IMN-1's `samples`, pooled: theta_sat their mean porosity, theta_res,
    alpha [1/m] and n fitted to their paired TEROS 21 matric potential and water content (-9 kPa and
    drier: the wet end is not measured)."""
    props = pd.read_csv(os.path.join(folder, "Soil_Physical_Properties.csv"))
    props = props[(props["soil_pit_number"] == "IMN-1") & props["pit_sample_id"].isin(samples)]
    psi = pd.read_csv(os.path.join(folder, "Hydrology_Matric_Potential.csv"))
    theta = pd.read_csv(os.path.join(folder, "Hydrology_Soil_Moisture.csv"))
    heads, waters = [], []
    for s in samples:
        col = f"IMN-1-{s}"
        h = np.abs(pd.to_numeric(psi[col], errors="coerce").to_numpy()) * KPA_TO_M
        t = pd.to_numeric(theta[col], errors="coerce").to_numpy()
        ok = np.isfinite(h) & np.isfinite(t) & (h > 0) & (h < 9e4 * KPA_TO_M) & (t > 0)
        heads.append(h[ok])
        waters.append(t[ok])
    h, t = np.concatenate(heads), np.concatenate(waters)
    theta_sat = float(props["porosity"].mean())
    fit = least_squares(lambda p: vg_theta_head(h, theta_sat, p[0], 10 ** p[1], p[2]) - t, [0.05, 0.0, 1.5],
                        bounds=([0.0, -3.0, 1.05], [0.9 * theta_sat, 3.0, 4.0]))
    return dict(theta_sat=theta_sat, theta_res=float(fit.x[0]), alpha=float(10 ** fit.x[1]), n=float(fit.x[2]),
                bd=float(props["field_bulk_density"].mean()), n_obs=len(t),
                rmse=float(np.sqrt(np.mean(fit.fun ** 2))))


def permafrost_dry_density(data):
    """The median dry density [g/cm3] of the loamy layers below each site's thaw depth at Jorgenson's
    Toolik sites, and how many layers it is over."""
    book = pd.read_excel(os.path.join(data, JORGENSON), sheet_name=["tbl_Site", "tbl_Soil_Stratigraphy"])
    layers = book["tbl_Soil_Stratigraphy"].merge(book["tbl_Site"][["SiteID", "ThawDepth"]], on="SiteID")
    for col in ("DepthTop", "Dens-Strat-gcc", "ThawDepth"):
        layers[col] = pd.to_numeric(layers[col], errors="coerce")
    frozen = layers[layers["SiteID"].astype(str).str.startswith("TOOL") & (layers["TextAggr1"] == "L")
                    & (layers["ThawDepth"] < 900) & (layers["DepthTop"] >= layers["ThawDepth"])]
    density = frozen["Dens-Strat-gcc"].dropna()
    return float(density.median()), len(density)


def horizons(data):
    """The pit's horizons below its Oi, in the column's depths, with every [soil] property, and the
    permafrost below the thaw (Cgf)."""
    pit = sh.read_megapit(glob.glob(os.path.join(data, "neon", "DP1.00096.001", "*"))[0])
    imnavait = os.path.join(data, "ornl_above", "Arctic_Soil_Properties")
    bg = pit.loc["Bg"]
    bg_theta_sat = sh.horizon_hydraulics(bg["sand"], bg["silt"], bg["clay"], bg["oc"], bg["cec"], bg["ph"],
                                         bg["bd"])["theta_sat"]          # the mineral's porosity, fragments aside
    rows = []
    for name in list(ORGANIC) + MINERAL:
        h = pit.loc[name]
        r = dict(horizon=name, top=h["top"] - PIT_SOIL_TOP, bottom=h["bottom"] - PIT_SOIL_TOP, oc=h["oc"])
        if name in ORGANIC:
            curve = fit_imnavait(imnavait, ORGANIC[name]["samples"])
            f_om = sh.organic_fraction(curve["bd"], h["oc"])
            r.update(curve, ksat=ORGANIC[name]["ksat"], f_om=f_om,
                     **sh.organic_thermal(f_om, bg["sand"], bg["clay"], bg_theta_sat, curve["theta_sat"]))
        else:
            hyd = sh.horizon_hydraulics(h["sand"], h["silt"], h["clay"], h["oc"], h["cec"], h["ph"], h["bd"],
                                        h["rock_volume"])
            r.update(hyd, bd=h["bd"], f_om=0.0, **sh.texture_thermal(h["sand"], h["clay"], hyd["theta_sat"]))
        rows.append(r)
    rows[-1]["bottom"] = PERMAFROST_TOP                                   # the Cg's thawing part
    cg = pit.loc["Cg"]
    bd_frozen, n_layers = permafrost_dry_density(data)
    hyd = sh.horizon_hydraulics(cg["sand"], cg["silt"], cg["clay"], cg["oc"], cg["cec"], cg["ph"], bd_frozen,
                                cg["rock_volume"])
    rows.append(dict(horizon="Cgf", top=PERMAFROST_TOP, bottom=LAYER_BOTTOMS[-1], oc=cg["oc"], bd=bd_frozen, f_om=0.0,
                     **hyd, **sh.texture_thermal(cg["sand"], cg["clay"], hyd["theta_sat"])))
    print(f"the permafrost below {PERMAFROST_TOP} m: dry density {bd_frozen:.2f} g/cm3, the median of {n_layers} "
          "loamy layers below the thaw depth at Jorgenson's (2016) Toolik sites")
    return pd.DataFrame(rows).set_index("horizon")


def soil_column(data):
    t = horizons(data)
    print("soil horizons below the pit's Oi (depths [m] in the column), NEON megapit (DP1.00096.001);")
    print("Oe and Oa: Imnavait IMN-1's curves (ORNL DAAC 2149), Letts et al.'s peat Ks, CLM5's thermal blend:")
    print(t[["top", "bottom", "oc", "bd", "f_om", "theta_sat", "theta_res", "alpha", "n", "ksat",
             "solid_conductivity", "dry_conductivity", "dry_heat_capacity"]].to_string(float_format=lambda v: f"{v:.4g}"))
    for name in ORGANIC:
        print(f"  {name}'s curve: {t.loc[name, 'n_obs']:.0f} paired readings, rmse {t.loc[name, 'rmse']:.3f} m3/m3")
    bottoms = np.array(LAYER_BOTTOMS)
    tops = np.concatenate([[0.0], bottoms[:-1]])
    layer_horizon = []
    for top, bottom in zip(tops, bottoms):
        inside = t[(t["top"] <= top + 1e-9) & (t["bottom"] >= bottom - 1e-9)]
        if len(inside) != 1:
            raise SystemExit(f"ERROR: the layer {top}-{bottom} m straddles a horizon boundary")
        layer_horizon.append(inside.index[0])
    print(f"the column: {len(bottoms)} layers to {bottoms[-1]} m, by horizon {layer_horizon}")
    keys = {"theta_sat": ("theta_sat", "{:.3f}"), "theta_res": ("theta_res", "{:.3f}"),
            "curve_par_a": ("alpha", "{:.3g}"), "curve_par_n": ("n", "{:.3g}"), "ksat": ("ksat", "{:.3g}"),
            "solid_conductivity": ("solid_conductivity", "{:.3g}"), "dry_conductivity": ("dry_conductivity", "{:.3g}"),
            "dry_heat_capacity": ("dry_heat_capacity", "{:.4g}")}
    print("[soil]")
    print(f"layer_bottoms = [{', '.join(f'{b:g}' for b in bottoms)}]")
    for key, (col, fmt) in keys.items():
        print(f"{key} = [{', '.join(fmt.format(t.loc[h, col]) for h in layer_horizon)}]")
    return tops, bottoms


def read_soil_temperature(data, start, end):
    """NEON's soil temperature between `start` and `end` (UTC): each sensor's mean [degC] and its depth
    [m] from the sensor positions on `start`; a sensor above the ground (zOffset >= 0) is left out."""
    folder = os.path.join(data, "neon", "DP1.00041.001", start[:7])
    positions = pd.read_csv(glob.glob(os.path.join(folder, "*sensor_positions*"))[0], dtype={"HOR.VER": str})
    day = pd.Timestamp(start, tz="UTC")
    rows = []
    for f in sorted(glob.glob(os.path.join(folder, "*ST_30_minute*.csv"))):
        m = re.search(r"\.(00\d)\.(50\d)\.030\.", os.path.basename(f))
        x = pd.read_csv(f, usecols=["startDateTime", "soilTempMean", "finalQF"])
        x = x[(x["finalQF"] == 0) & (x["startDateTime"] >= start) & (x["startDateTime"] < end)]
        p = positions[positions["HOR.VER"] == f"{m.group(1)}.{m.group(2)}"]
        p = p[(pd.to_datetime(p["start"]) <= day) & (pd.to_datetime(p["end"]).fillna(pd.Timestamp.max.tz_localize("UTC")) > day)]
        if x.empty or p.empty or p["zOffset"].iloc[0] >= 0:
            continue
        rows.append(dict(plot=int(m.group(1)), depth=-float(p["zOffset"].iloc[0]), temp=float(x["soilTempMean"].mean())))
    return pd.DataFrame(rows)


def soil_start(data, tops, bottoms):
    s = read_soil_temperature(data, *START)
    print(f"the soil on {START[0]} to {START[1]} (NEON DP1.00041.001), each plot's sensors [m, degC]:")
    for plot, g in s.groupby("plot"):
        print(f"  plot {plot}: " + ", ".join(f"{d:.2f} m {v:.2f}" for d, v in zip(g["depth"], g["temp"])))
    plots = [g.sort_values("depth") for _, g in s.groupby("plot")]
    deepest = s["depth"].max()

    def measured(z):
        """The plots' mean [K] at depth z, over the plots whose sensors reach it."""
        return np.mean([np.interp(z, g["depth"], g["temp"]) for g in plots if g["depth"].max() >= z]) + 273.15

    t_deepest = measured(deepest)
    nodes = 0.5 * (tops + bottoms)
    temp = np.array([measured(z) if z <= deepest else
                     t_deepest + (DEEP_TEMP - t_deepest) * (z - deepest) / (DEEP_DEPTH - deepest) for z in nodes])
    print(f"  at the layers' middles above {deepest:.2f} m, the mean of the plots reaching each; below, a line "
          f"from {t_deepest:.2f} K to {DEEP_TEMP:.2f} K at {DEEP_DEPTH} m")
    print("[init]")
    print(f"soil_temp = [{', '.join(f'{v:.2f}' for v in temp)}]")
    print("[soil] (the deep temperature)")
    print(f"deep_temp = {DEEP_TEMP:.2f}\ndeep_depth = {DEEP_DEPTH:.2f}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", required=True, help="the Toolik data folder (neon/, ornl_above/)")
    args = ap.parse_args(argv)
    tops, bottoms = soil_column(args.data)
    soil_start(args.data, tops, bottoms)


if __name__ == "__main__":
    main()
