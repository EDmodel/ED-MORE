# SPDX-License-Identifier: Apache-2.0
"""Example 07 (prototype): ED-MORE's canopy radiation against NEON's vertical PAR sensors at Harvard Forest.

The two-stream (edmore.canopy) runs every daytime half hour of one year, driven by the PAR measured
above the canopy (NEON's top level, 39 m) and split into beam and diffuse by the EMS tower's direct
and diffuse PAR sensors. The stand around each sensor comes from the 2019 census of the ForestGEO
35-ha plot (one patch per 20 m cell, the plot's quadrat, one cohort per stem, ED-MORE's default
allometry, optionally with a site height cap); each deciduous species' leaf area follows the ground
observations of tagged trees (HF003) through the year. A sensor at height z reads the downwelling
PAR at the top of the tallest cohort layer at or below z (the ground's when none is), averaged over
the patches by area: the horizontal mean a 1-D layered canopy predicts.

    python run_example.py --data <data> --hgt-max 28      (<data>: the folder the fetch scripts fill)

Writes output/model_vs_sensors.csv (one row per half hour: drivers, each sensor observed and
modelled, and the modelled 10th/90th patch percentiles) and prints the transmittance summary.
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
RADIUS = 40.0                # [m] the stand around each sensor
ALLOM_KEYS = ("b1Ht", "b2Ht", "hgt_max", "lai_b1", "lai_b2", "wai_b1", "wai_b2")


def summary(df: pd.DataFrame, sensors: list[str]) -> pd.DataFrame:
    """Sensor PAR over top PAR (ratio of sums) by season and sky, observed and modelled, over the
    half hours where the sensor, the top and the split are all valid and the sun is above ~20 deg."""
    rows = []
    month = df.index.month
    seasons = {"leaf-off (Apr)": month == 4, "leaf-on (Jul-Aug)": (month == 7) | (month == 8),
               "autumn (Oct)": month == 10}
    skies = {"all": np.ones(len(df), bool), "clear (diffuse < 0.35)": df["diffuse_frac"] < 0.35,
             "overcast (diffuse > 0.8)": df["diffuse_frac"] > 0.8}
    for sname, smask in seasons.items():
        for kname, kmask in skies.items():
            for s in sensors:
                ok = (smask & kmask & (df["cosz"] > 0.35) & df[f"obs_{s}"].notna() & df["par_top"].notna()
                      & df[f"mod_{s}"].notna())
                if ok.sum() < 20:
                    continue
                top = df.loc[ok, "par_top"].sum()
                rows.append(dict(season=sname, sky=kname, sensor=s, n=int(ok.sum()),
                                 observed=df.loc[ok, f"obs_{s}"].sum() / top,
                                 model=df.loc[ok, f"mod_{s}"].sum() / top))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--year", type=int, default=2019)
    ap.add_argument("--out", type=Path, default=HERE / "output")
    ap.add_argument("--cell", type=float, default=20.0, help="patch size [m]; 0 = one uniform patch")
    ap.add_argument("--phenology", choices=("hf003", "modis"), default="hf003",
                    help="deciduous leaf area through the year: ground observations by species (HF003) "
                         "or MODIS LAI scaled 0-1")
    ap.add_argument("--radius", type=float, default=RADIUS, help="the stand around each sensor [m]")
    ap.add_argument("--hgt-max", type=float, default=None,
                    help="cap on tree height [m]; default the PFT file's hgt_max")
    a = ap.parse_args()
    t0 = time.time()

    #----- The configuration: the repository's defaults, with the PFT file's allometry and optics. ---#
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    canopy = Canopy(cfg.write(a.out / "config"))
    allom = {k: [cfg.get(f"allometry.{k}", file="pft", pft=p) for p in (1, 2)] for k in ALLOM_KEYS}
    if a.hgt_max is not None:
        allom["hgt_max"] = [a.hgt_max, a.hgt_max]

    #----- Sensors and drivers. ---------------------------------------------------------------------#
    tower, heights = hd.neon_tower_par(a.data)
    floor, floor_xy = hd.neon_floor_par(a.data)
    obs = tower.join(floor, how="outer")
    obs = obs[obs.index.year == a.year]
    idx = obs.index
    drv = pd.DataFrame(index=idx)
    drv["par_top"] = obs["L6"]
    drv["diffuse_frac"] = hd.ems_diffuse_fraction(a.data, idx)
    drv["cosz"] = conv.window_mean_cosz(idx.values.astype("datetime64[s]"), 1800, hd.TOWER_LAT, hd.TOWER_LON)
    run = drv["par_top"].notna() & drv["diffuse_frac"].notna() & (drv["cosz"] > 0.02) & (drv["par_top"] > 5.0)
    vis = drv["par_top"].where(run, 0.0).to_numpy() / PAR_PER_WATT
    vis_beam = vis * (1.0 - drv["diffuse_frac"].fillna(1.0).to_numpy())
    vis_diff = vis - vis_beam
    cosz = drv["cosz"].to_numpy()
    if a.phenology == "hf003":
        leaf_frac = hd.species_leaf_fraction(a.data, a.year)                  # days x species
    else:
        leaf_frac = hd.deciduous_fraction(a.data, a.year).to_frame("mean")    # days x 1
    frac_days = leaf_frac.index.normalize()
    frac_arr = leaf_frac.to_numpy()

    def frac_columns(species):
        """Each cohort's column in the leaf-fraction table: its species', else the mean."""
        cols = list(leaf_frac.columns)
        return np.array([cols.index(hd.hf003_code(s)) if hd.hf003_code(s) in cols else cols.index("mean")
                         for s in species], dtype=int)

    #----- Where each sensor stands, and at what height. -----------------------------------------------#
    x_tower, y_tower = hd.plot_xy(hd.TOWER_LAT, hd.TOWER_LON)
    sites = {"tower": ((x_tower, y_tower), {lv: heights[lv] for lv in ("L1", "L2", "L3", "L4", "L5")})}
    for sp, xy in floor_xy.items():
        sites[sp] = (xy, {sp: 0.0})
    stems = hd.live_stems(a.data)

    out = drv.copy()
    for s in list(heights) + list(floor_xy):
        if s in obs:
            out[f"obs_{s}"] = obs[s]
    day_of = idx.normalize()
    day_hours = [h for h in (np.flatnonzero((day_of == d) & run.to_numpy()) for d in np.unique(day_of))
                 if h.size]
    for site, ((x0, y0), sensors) in sites.items():
        patches = hd.stand_around(stems, x0, y0, a.radius, allom, cell=a.cell)
        weights = np.array([p["area"] for p in patches])
        weights = weights / weights.sum()
        per_patch = {s: np.full((len(patches), len(idx)), np.nan) for s in sensors}
        for k, p in enumerate(patches):
            fcol = frac_columns(p["sp"])
            for hours in day_hours:
                for s in sensors:
                    per_patch[s][k, hours] = vis[hours]          # no stems: the open sky's light
                if p["height"].size == 0:
                    continue
                f = frac_arr[frac_days.get_loc(day_of[hours[0]]), fcol]
                lai = np.where(p["pft"] == 1, p["lai"] * f, p["lai"]) + 1e-9
                r = canopy.radiation(cosz[hours], vis_beam[hours], vis_diff[hours], vis_beam[hours],
                                     vis_diff[hours], p["pft"], p["height"], lai, p["wai"])
                for s, z in sensors.items():
                    below = np.flatnonzero(p["height"] <= z)
                    if below.size == 0:
                        per_patch[s][k, hours] = r["ground_vis"]
                    else:
                        tallest = below[np.argmax(p["height"][below])]
                        per_patch[s][k, hours] = r["cohort_incid_vis"][:, tallest]
        for s in sensors:
            m = per_patch[s]
            out[f"mod_{s}"] = np.where(run, (weights[:, None] * m).sum(axis=0), np.nan) * PAR_PER_WATT
            out[f"mod_{s}_p10"] = np.where(run, np.nanpercentile(m, 10, axis=0), np.nan) * PAR_PER_WATT
            out[f"mod_{s}_p90"] = np.where(run, np.nanpercentile(m, 90, axis=0), np.nan) * PAR_PER_WATT
        lai_on = sum(w * (p["lai"].sum()) for w, p in zip(weights, patches))
        wai = sum(w * (p["wai"].sum()) for w, p in zip(weights, patches))
        lai_ev = sum(w * (p["lai"][p["pft"] == 2].sum()) for w, p in zip(weights, patches))
        print(f"{site}: {len(patches)} patches, {sum(p['height'].size for p in patches)} stems; "
              f"full-leaf LAI {lai_on:.2f} (evergreen {lai_ev:.2f}), WAI {wai:.2f}", flush=True)

    a.out.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out / "model_vs_sensors.csv", index_label="start_utc")
    sensors = ["L5", "L4", "L3", "L2", "L1"] + list(floor_xy)
    tab = summary(out, sensors)
    tab.to_csv(a.out / "transmittance_summary.csv", index=False)
    with pd.option_context("display.width", 140, "display.max_rows", 200, "display.float_format", "{:.3f}".format):
        print(tab.pivot_table(index=["season", "sky"], columns="sensor", values=["observed", "model"])
              .reindex(columns=sensors, level=1))
    print(f"heights [m]: " + ", ".join(f"{k} {v:.2f}" for k, v in heights.items()))
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
