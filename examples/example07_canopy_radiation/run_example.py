# SPDX-License-Identifier: Apache-2.0
"""Example 07 (prototype): ED-MORE's canopy radiation against NEON's vertical PAR sensors at Harvard Forest.

The two-stream (edmore.canopy) runs every daytime half hour of two windows of 2019: leaf-off (1-25
April: the snow gone, the buds not yet broken; only the conifers have leaves) and leaf-on (July). It is
driven by the PAR measured above the canopy (NEON's top level, 39 m), split into beam and diffuse by
the EMS tower's direct and diffuse PAR sensors. The stand around each sensor comes from the 2019 census
of the ForestGEO 35-ha plot: one patch per 20 m quadrat, one cohort per stem, with Harvard's height and
leaf-area allometry (HARVARD_PFT below). A sensor at height z reads the downwelling PAR at the top of
the tallest cohort layer at or below z (the ground's when none is), averaged over the patches by area:
the horizontal mean a 1-D layered canopy predicts.

    python run_example.py --data <data>        (<data>: the folder the fetch scripts fill)

Writes output/model_vs_sensors.csv (one row per half hour: drivers, each sensor observed and modelled,
and the modelled 10th/90th patch percentiles) and output/transmittance_summary.csv, and prints the table.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts" / "prepare_input" / "flux_tower"))
from edmore.canopy import Canopy               # noqa: E402
from edmore.config import RunConfig            # noqa: E402
import tower_conversions as conv               # noqa: E402  (the model's sun, in Python)

import harvard_data as hd                       # noqa: E402

PAR_PER_WATT = 4.57          # [umol J-1] PAR photons per joule; only the absolute fluxes depend on it
WINDOWS = {"leaf-off": ("2019-04-01", "2019-04-25"), "leaf-on": ("2019-07-01", "2019-07-31")}
#----- Harvard's allometry on the repository's PFT file: PFT 1 deciduous, PFT 2 evergreen (the conifers).
#      Heights: ED-MORE's log-quadratic curve fitted to NEON's measured trees at HARV (3,400 trees, each
#      DBH class weighted equally), capped at the canopy top the lidar shows. Deciduous leaf area: scaled
#      to a July LAI of 3.95 around the towers, where litterfall x leaf mass per area, MODIS and LAI-2000
#      agree; fit_allometry.py reproduces both. ------------------------------------------------------#
HARVARD_PFT = {
    "allometry.height_allometry": ["chave2014", "chave2014"],
    "allometry.chave_a":   [0.47289, 0.43288],
    "allometry.chave_b":   [1.0693, 0.71177],
    "allometry.chave_c":   [-0.097823, 0.00026528],
    "allometry.chave_rse": [0.25066, 0.25667],
    "allometry.hgt_max":   [28.0, 28.0],
    "allometry.lai_b1":    [0.1711, 0.23384770],
}
ALLOM_KEYS = ("height_allometry", "chave_a", "chave_b", "chave_c", "chave_rse", "chave_e", "b1Ht", "b2Ht",
              "hgt_max", "lai_b1", "lai_b2", "wai_b1", "wai_b2")
#----- Two optics, the run's two models: "mod", the foliage at random within each cohort layer (leaf and
#      wood clumping 1.0), since the 20 m patches already carry the gaps between crowns; and "def", the
#      repository's default clumping (leaf 0.8, wood 0.5). Every other optical key is the default. ----#
OPTICS = {"mod": ["pft.leaf_clumping=1.0", "pft.wood_clumping=1.0"], "def": []}


def harvard_config(out: Path, extra: list[str]) -> RunConfig:
    """The repository's main and PFT files with Harvard's allometry for PFTs 1 and 2, and any
    `key=value` of the PFT file set for every PFT (e.g. pft.leaf_clumping=1.0), written into `out`."""
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    for key, values in HARVARD_PFT.items():                 # PFTs past the second (unused) copy the first
        cfg.set(key, list(values) + [values[0]] * (cfg.n_pft - len(values)), file="pft")
    cfg.set("allometry.chave_e", 0.0, file="pft")            # the site's stress index: one value, all PFTs
    for item in extra:
        key, value = item.split("=")
        cfg.set(key, [float(value)] * cfg.n_pft, file="pft")
    cfg.write(out)
    return cfg


def summary(df: pd.DataFrame, sensors: list[str]) -> pd.DataFrame:
    """Sensor PAR over top PAR (ratio of sums) by window and sky, observed and modelled, over the half
    hours where the sensor, the top, the split and the model are all valid and the sun is above ~20 deg."""
    rows = []
    skies = {"all": np.ones(len(df), bool), "clear (diffuse < 0.35)": df["diffuse_frac"] < 0.35,
             "overcast (diffuse > 0.8)": df["diffuse_frac"] > 0.8}
    for wname in WINDOWS:
        for kname, kmask in skies.items():
            for s in sensors:
                ok = ((df["window"] == wname) & kmask & (df["cosz"] > 0.35) & df[f"obs_{s}"].notna()
                      & df["par_top"].notna() & df[f"mod_{s}"].notna())
                if ok.sum() < 20:
                    continue
                top = df.loc[ok, "par_top"].sum()
                rows.append(dict(window=wname, sky=kname, sensor=s, n=int(ok.sum()),
                                 observed=df.loc[ok, f"obs_{s}"].sum() / top,
                                 model=df.loc[ok, f"mod_{s}"].sum() / top,
                                 model_default_clumping=df.loc[ok, f"def_{s}"].sum() / top))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=HERE / "output")
    ap.add_argument("--cell", type=float, default=20.0, help="patch size [m]; 0 = one uniform patch")
    ap.add_argument("--radius", type=float, default=40.0, help="the stand around each sensor [m]")
    ap.add_argument("--pft-set", nargs="*", default=[], metavar="KEY=VALUE",
                    help="PFT-file keys set for every PFT, e.g. pft.leaf_clumping=1.0")
    a = ap.parse_args()
    t0 = time.time()

    canopies = {}
    for key, extra in OPTICS.items():
        cfg = harvard_config(a.out / f"config_{key}", extra + a.pft_set)
        canopies[key] = Canopy(a.out / f"config_{key}" / "main.toml")
    allom = {k: [cfg.get(f"allometry.{k}", file="pft", pft=p) for p in (1, 2)] for k in ALLOM_KEYS}

    #----- Sensors and drivers, in the two windows. -----------------------------------------------------#
    tower, heights = hd.neon_tower_par(a.data)
    floor, floor_xy = hd.neon_floor_par(a.data)
    obs = tower.join(floor, how="outer")
    window = pd.Series(np.nan, index=obs.index, dtype=object)
    for name, (t_start, t_end) in WINDOWS.items():
        window[(obs.index >= t_start) & (obs.index < pd.Timestamp(t_end) + pd.Timedelta(days=1))] = name
    obs = obs[window.notna()]
    idx = obs.index
    drv = pd.DataFrame(index=idx)
    drv["window"] = window[window.notna()]
    drv["par_top"] = obs["L6"]
    drv["diffuse_frac"] = hd.ems_diffuse_fraction(a.data, idx)
    drv["cosz"] = conv.window_mean_cosz(idx.values.astype("datetime64[s]"), 1800, hd.TOWER_LAT, hd.TOWER_LON)
    run = (drv["par_top"].notna() & drv["diffuse_frac"].notna() & (drv["cosz"] > 0.02)
           & (drv["par_top"] > 5.0)).to_numpy()
    vis = drv["par_top"].where(run, 0.0).to_numpy() / PAR_PER_WATT
    vis_beam = vis * (1.0 - drv["diffuse_frac"].fillna(1.0).to_numpy())
    vis_diff = vis - vis_beam
    cosz = drv["cosz"].to_numpy()
    leaf_on = (drv["window"] == "leaf-on").to_numpy()

    #----- Where each sensor stands, and at what height. -----------------------------------------------#
    x_tower, y_tower = hd.plot_xy(hd.TOWER_LAT, hd.TOWER_LON)
    sites = {"tower": ((x_tower, y_tower), {lv: heights[lv] for lv in ("L1", "L2", "L3", "L4", "L5")})}
    for sp, xy in floor_xy.items():
        sites[sp] = (xy, {sp: 0.0})
    stems = hd.live_stems(a.data)

    out = drv.copy()
    for s in list(heights) + list(floor_xy):
        out[f"obs_{s}"] = obs[s]
    for site, ((x0, y0), sensors) in sites.items():
        patches = hd.stand_around(stems, x0, y0, a.radius, allom, cell=a.cell)
        weights = np.array([p["area"] for p in patches])
        weights = weights / weights.sum()
        per_patch = {(o, s): np.full((len(patches), len(idx)), np.nan) for o in OPTICS for s in sensors}
        for k, p in enumerate(patches):
            for season_on in (False, True):
                hours = np.flatnonzero(run & (leaf_on == season_on))
                for key in per_patch:
                    per_patch[key][k, hours] = vis[hours]        # no stems: the open sky's light
                if p["height"].size == 0 or hours.size == 0:
                    continue
                lai = np.where(p["pft"] == 1, p["lai"] * float(season_on), p["lai"]) + 1e-9
                for o, canopy in canopies.items():
                    r = canopy.radiation(cosz[hours], vis_beam[hours], vis_diff[hours], vis_beam[hours],
                                         vis_diff[hours], p["pft"], p["height"], lai, p["wai"])
                    for s, z in sensors.items():
                        below = np.flatnonzero(p["height"] <= z)
                        if below.size == 0:
                            per_patch[o, s][k, hours] = r["ground_vis"]
                        else:
                            tallest = below[np.argmax(p["height"][below])]
                            per_patch[o, s][k, hours] = r["cohort_incid_vis"][:, tallest]
        for (o, s), m in per_patch.items():
            out[f"{o}_{s}"] = np.where(run, (weights[:, None] * m).sum(axis=0), np.nan) * PAR_PER_WATT
            out[f"{o}_{s}_p10"] = np.where(run, np.percentile(m, 10, axis=0), np.nan) * PAR_PER_WATT
            out[f"{o}_{s}_p90"] = np.where(run, np.percentile(m, 90, axis=0), np.nan) * PAR_PER_WATT
        lai_dec = sum(w * p["lai"][p["pft"] == 1].sum() for w, p in zip(weights, patches))
        lai_ev = sum(w * p["lai"][p["pft"] == 2].sum() for w, p in zip(weights, patches))
        wai = sum(w * p["wai"].sum() for w, p in zip(weights, patches))
        top = np.mean([p["height"].max() for p in patches if p["height"].size])
        print(f"{site}: {len(patches)} patches, {sum(p['height'].size for p in patches)} stems; LAI deciduous "
              f"{lai_dec:.2f}, evergreen {lai_ev:.2f}; WAI {wai:.2f}; tallest tree per patch {top:.1f} m", flush=True)

    a.out.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out / "model_vs_sensors.csv", index_label="start_utc")
    pd.Series(heights).to_csv(a.out / "sensor_heights.csv", header=["height_m"], index_label="sensor")
    sensors = ["L5", "L4", "L3", "L2", "L1"] + list(floor_xy)
    tab = summary(out, sensors)
    tab.to_csv(a.out / "transmittance_summary.csv", index=False)
    with pd.option_context("display.width", 160):
        tab["obs / model / default clumping"] = (tab["observed"].map("{:.3f}".format) + " / "
                                                 + tab["model"].map("{:.3f}".format) + " / "
                                                 + tab["model_default_clumping"].map("{:.3f}".format))
        print(tab.pivot_table(index="sensor", columns=["window", "sky"], values="obs / model / default clumping",
                              aggfunc="first").reindex(sensors).to_string())
    print("heights [m]: " + ", ".join(f"{k} {v:.2f}" for k, v in heights.items()))
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
