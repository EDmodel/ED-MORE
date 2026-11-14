# SPDX-License-Identifier: Apache-2.0
"""Calibrate the foliage's clumping and leaf angle against NEON's PAR sensors below the crowns.

Three keys of every PFT, with the wood angle at its default (75 deg, from TLS cylinder models):
leaf clumping, the mean leaf angle and wood clumping. The targets are the light at the sensors the 1-D
layers can be compared with, 5 m, 0.4 m and the three floor sensors (the 16 and 24 m sensors sit
inside the crowns, where the layers' crown-depth assumption decides), as ratios of sums to the light
above the canopy, in both windows: under clear skies in three bands of the sun's height (which the
leaf angle shapes) and under overcast skies. The cost is the mean squared log of model over observed.
A grid, then Nelder-Mead from its best point; calibrated_optics.json keeps the result.

    python calibrate_optics.py --data <data>
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

import run_example as rx
from edmore.canopy import Canopy

SENSORS = ["L2", "L1", "SP1", "SP3", "SP5"]
KEYS = ("pft.leaf_clumping", "pft.leaf_angle_mean", "pft.wood_clumping")
BOUNDS = ((0.4, 1.0), (20.0, 75.0), (0.3, 1.0))
GRID = (np.arange(0.5, 1.01, 0.1), np.arange(25.0, 71.0, 10.0), np.array([0.4, 0.55, 0.7, 0.85, 1.0]))
COSZ_BANDS = ((0.35, 0.55), (0.55, 0.75), (0.75, 1.01))


def target_masks(drv) -> dict:
    """The half-hour sets each target sums over, per window: clear sky by sun band, and overcast."""
    masks = {}
    for w in rx.WINDOWS:
        base = (drv["window"] == w).to_numpy()
        for lo, hi in COSZ_BANDS:
            masks[f"{w}, clear, cos z {lo:.2f}-{hi:.2f}"] = base & (drv["diffuse_frac"] < 0.35).to_numpy() \
                & (drv["cosz"] >= lo).to_numpy() & (drv["cosz"] < hi).to_numpy()
        masks[f"{w}, overcast"] = base & (drv["diffuse_frac"] > 0.8).to_numpy() & (drv["cosz"] > 0.35).to_numpy()
    return masks


def residuals(prep, masks, modelled) -> dict:
    """ln(model / observed) of each (target, sensor) with at least 10 half hours of both."""
    drv, out = prep["drv"], {}
    top = drv["par_top"].to_numpy()
    for name, m in masks.items():
        for s in SENSORS:
            obs = drv[f"obs_{s}"].to_numpy()
            ok = m & np.isfinite(obs) & np.isfinite(modelled[s]) & np.isfinite(top)
            if ok.sum() >= 10:
                out[(name, s)] = np.log((modelled[s][ok].sum() / top[ok].sum()) / (obs[ok].sum() / top[ok].sum()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=rx.HERE / "output" / "calibration")
    a = ap.parse_args()
    t0 = time.time()
    prep = rx.prepare(a.data)
    masks = target_masks(prep["drv"])
    #----- Every evaluation is kept on disk, so an interrupted calibration resumes where it stopped. ---#
    a.work.mkdir(parents=True, exist_ok=True)
    store = a.work / "evaluations.json"
    saved = json.loads(store.read_text()) if store.exists() else []
    cache = {tuple(e["x"]): (e["cost"], {tuple(k.split("|")): v for k, v in e["res"].items()}) for e in saved}

    def cost(x):
        x = tuple(round(float(np.clip(v, lo, hi)), 6) for v, (lo, hi) in zip(x, BOUNDS))
        if x not in cache:
            rx.harvard_config(a.work, dict(zip(KEYS, x)))
            res = residuals(prep, masks, rx.simulate(prep, Canopy(a.work / "main.toml"), SENSORS, spread=False))
            cache[x] = (float(np.mean(np.square(list(res.values())))), res)
            saved.append(dict(x=list(x), cost=cache[x][0], res={"|".join(k): v for k, v in res.items()}))
            store.write_text(json.dumps(saved))
        return cache[x][0]

    default = cost((0.8, 45.0, 0.5))
    grid = [tuple(round(float(v), 6) for v in (c, t, w)) for c in GRID[0] for t in GRID[1] for w in GRID[2]]
    for i, x in enumerate(grid):
        cost(x)
        if i % 25 == 0:
            print(f"grid {i}/{len(grid)} ({time.time() - t0:.0f} s)", flush=True)
    start = np.array(min(grid, key=cost))
    steps = np.diag([-0.05, 5.0, -0.075])                    # inward from the upper bounds of the clumpings
    simplex = np.vstack([start, start + steps])
    fit = minimize(cost, start, method="Nelder-Mead",
                   options=dict(xatol=0.01, fatol=1e-4, maxfev=60, initial_simplex=simplex))
    best = tuple(round(float(np.clip(v, lo, hi)), 6) for v, (lo, hi) in zip(fit.x, BOUNDS))
    cost(best)

    #----- What the data decide: the cost's floor along each key, the others free (from the grid). ------#
    print(f"\ncost (mean squared log model/observed): default (0.8, 45, 0.5) {default:.4f}; best {cache[best][0]:.4f}")
    for k, key in enumerate(KEYS):
        prof = {v: min(c for x, (c, _) in cache.items() if abs(x[k] - v) < 1e-9) for v in GRID[k]}
        print(f"  {key}: profile " + ", ".join(f"{v:g}: {c:.4f}" for v, c in prof.items()))
    print("\nresiduals ln(model/observed), best / default:")
    for (name, s), r in cache[best][1].items():
        print(f"  {name:32s} {s:4s} {r:+.2f} / {cache[(0.8, 45.0, 0.5)][1].get((name, s), np.nan):+.2f}")
    result = dict(keys=dict(zip(KEYS, [round(v, 3) for v in best])), cost=cache[best][0], cost_default=default,
                  wood_angle_mean="the PFT file's default (75)", sensors=SENSORS, targets=list(masks),
                  evaluations=len(cache), seconds=round(time.time() - t0))
    rx.CALIBRATED.write_text(json.dumps(result, indent=1) + "\n")
    print(f"\n{rx.CALIBRATED.name}: {result['keys']} ({len(cache)} runs, {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
