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

Two optics run side by side: "mod", the clumping and leaf angle calibrate_optics.py fitted to the
sensors below the crowns (calibrated_optics.json), and "def", the repository's defaults.

    python run_example.py --data <data>        (<data>: the folder the fetch scripts fill)

Writes output/model_vs_sensors.csv (one row per half hour: drivers, each sensor observed and modelled,
and the modelled 10th/90th patch percentiles) and output/transmittance_summary.csv, and prints the table.
"""
from __future__ import annotations

import argparse
import json
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
CALIBRATED = HERE / "calibrated_optics.json"


def harvard_config(out: Path | None, extra: dict) -> RunConfig:
    """The repository's main and PFT files with Harvard's allometry for PFTs 1 and 2, and the PFT-file
    keys in `extra` ({"pft.leaf_clumping": 0.9, ...}) set for every PFT, written into `out` if given."""
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    for key, values in HARVARD_PFT.items():                 # PFTs past the second (unused) copy the first
        cfg.set(key, list(values) + [values[0]] * (cfg.n_pft - len(values)), file="pft")
    cfg.set("allometry.chave_e", 0.0, file="pft")            # the site's stress index: one value, all PFTs
    for key, value in extra.items():
        cfg.set(key, [float(value)] * cfg.n_pft, file="pft")
    if out is not None:
        cfg.write(out)
    return cfg


def prepare(data: Path, cell: float = 20.0, radius: float = 40.0) -> dict:
    """The drivers, the sensors and the stand around each, for the two windows."""
    cfg = harvard_config(None, {})
    allom = {k: [cfg.get(f"allometry.{k}", file="pft", pft=p) for p in (1, 2)] for k in ALLOM_KEYS}
    tower, heights = hd.neon_tower_par(data)
    floor, floor_xy = hd.neon_floor_par(data)
    obs = tower.join(floor, how="outer")
    window = pd.Series(np.nan, index=obs.index, dtype=object)
    for name, (t_start, t_end) in WINDOWS.items():
        window[(obs.index >= t_start) & (obs.index < pd.Timestamp(t_end) + pd.Timedelta(days=1))] = name
    obs = obs[window.notna()]
    idx = obs.index
    drv = pd.DataFrame(index=idx)
    drv["window"] = window[window.notna()]
    drv["par_top"] = obs["L6"]
    drv["diffuse_frac"] = hd.ems_diffuse_fraction(data, idx)
    drv["cosz"] = conv.window_mean_cosz(idx.values.astype("datetime64[s]"), 1800, hd.TOWER_LAT, hd.TOWER_LON)
    run = (drv["par_top"].notna() & drv["diffuse_frac"].notna() & (drv["cosz"] > 0.02)
           & (drv["par_top"] > 5.0)).to_numpy()
    vis = drv["par_top"].where(run, 0.0).to_numpy() / PAR_PER_WATT
    vis_beam = vis * (1.0 - drv["diffuse_frac"].fillna(1.0).to_numpy())
    x_tower, y_tower = hd.plot_xy(hd.TOWER_LAT, hd.TOWER_LON)
    sites = {"tower": ((x_tower, y_tower), {lv: heights[lv] for lv in ("L1", "L2", "L3", "L4", "L5")})}
    for sp, xy in floor_xy.items():
        sites[sp] = (xy, {sp: 0.0})
    stems = hd.live_stems(data)
    stands = {site: hd.stand_around(stems, x0, y0, radius, allom, cell=cell) for site, ((x0, y0), _) in sites.items()}
    for s in list(heights) + list(floor_xy):
        drv[f"obs_{s}"] = obs[s]
    return dict(drv=drv, run=run, vis=vis, vis_beam=vis_beam, vis_diff=vis - vis_beam,
                cosz=drv["cosz"].to_numpy(), leaf_on=(drv["window"] == "leaf-on").to_numpy(),
                sites=sites, stands=stands, heights=heights, floor=list(floor_xy))


def simulate(prep: dict, canopy: Canopy, sensors: list[str] | None = None, spread: bool = True) -> dict:
    """The modelled PAR [umol m-2 s-1] at each sensor (or only those in `sensors`): its patches' area
    mean, and with `spread` their 10th and 90th percentiles, one value per half hour."""
    run, vis, vb, vd, cosz = prep["run"], prep["vis"], prep["vis_beam"], prep["vis_diff"], prep["cosz"]
    n = len(run)
    out = {}
    for site, (_, site_sensors) in prep["sites"].items():
        site_sensors = {s: z for s, z in site_sensors.items() if sensors is None or s in sensors}
        if not site_sensors:
            continue
        patches = prep["stands"][site]
        weights = np.array([p["area"] for p in patches])
        weights = weights / weights.sum()
        per_patch = {s: np.full((len(patches), n), np.nan) for s in site_sensors}
        for k, p in enumerate(patches):
            for season_on in (False, True):
                hours = np.flatnonzero(run & (prep["leaf_on"] == season_on))
                for s in site_sensors:
                    per_patch[s][k, hours] = vis[hours]          # no stems: the open sky's light
                if p["height"].size == 0 or hours.size == 0:
                    continue
                lai = np.where(p["pft"] == 1, p["lai"] * float(season_on), p["lai"]) + 1e-9
                r = canopy.radiation(cosz[hours], vb[hours], vd[hours], vb[hours], vd[hours],
                                     p["pft"], p["height"], lai, p["wai"])
                for s, z in site_sensors.items():
                    below = np.flatnonzero(p["height"] <= z)
                    if below.size == 0:
                        per_patch[s][k, hours] = r["ground_vis"]
                    else:
                        tallest = below[np.argmax(p["height"][below])]
                        per_patch[s][k, hours] = r["cohort_incid_vis"][:, tallest]
        for s, m in per_patch.items():
            out[s] = np.where(run, (weights[:, None] * m).sum(axis=0), np.nan) * PAR_PER_WATT
            if spread:
                out[f"{s}_p10"] = np.where(run, np.percentile(m, 10, axis=0), np.nan) * PAR_PER_WATT
                out[f"{s}_p90"] = np.where(run, np.percentile(m, 90, axis=0), np.nan) * PAR_PER_WATT
    return out


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
                                 model_default=df.loc[ok, f"def_{s}"].sum() / top))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=HERE / "output")
    ap.add_argument("--cell", type=float, default=20.0, help="patch size [m]; 0 = one uniform patch")
    ap.add_argument("--radius", type=float, default=40.0, help="the stand around each sensor [m]")
    a = ap.parse_args()
    t0 = time.time()

    calibrated = json.loads(CALIBRATED.read_text())["keys"] if CALIBRATED.exists() else {}
    if not calibrated:
        print(f"no {CALIBRATED.name}: 'mod' runs the defaults too (run calibrate_optics.py)")
    prep = prepare(a.data, a.cell, a.radius)
    out = prep["drv"].copy()
    for name, keys in (("mod", calibrated), ("def", {})):
        harvard_config(a.out / f"config_{name}", keys)
        for s, v in simulate(prep, Canopy(a.out / f"config_{name}" / "main.toml")).items():
            out[f"{name}_{s}"] = v
    for site, patches in prep["stands"].items():
        w = np.array([p["area"] for p in patches]); w = w / w.sum()
        lai_dec = sum(x * p["lai"][p["pft"] == 1].sum() for x, p in zip(w, patches))
        lai_ev = sum(x * p["lai"][p["pft"] == 2].sum() for x, p in zip(w, patches))
        wai = sum(x * p["wai"].sum() for x, p in zip(w, patches))
        top = np.mean([p["height"].max() for p in patches if p["height"].size])
        print(f"{site}: {len(patches)} patches, {sum(p['height'].size for p in patches)} stems; LAI deciduous "
              f"{lai_dec:.2f}, evergreen {lai_ev:.2f}; WAI {wai:.2f}; tallest tree per patch {top:.1f} m")

    a.out.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out / "model_vs_sensors.csv", index_label="start_utc")
    pd.Series(prep["heights"]).to_csv(a.out / "sensor_heights.csv", header=["height_m"], index_label="sensor")
    sensors = ["L5", "L4", "L3", "L2", "L1"] + prep["floor"]
    tab = summary(out, sensors)
    tab.to_csv(a.out / "transmittance_summary.csv", index=False)
    with pd.option_context("display.width", 170):
        tab["obs / calibrated / default"] = (tab["observed"].map("{:.3f}".format) + " / "
                                             + tab["model"].map("{:.3f}".format) + " / "
                                             + tab["model_default"].map("{:.3f}".format))
        print(tab.pivot_table(index="sensor", columns=["window", "sky"], values="obs / calibrated / default",
                              aggfunc="first").reindex(sensors).to_string())
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
