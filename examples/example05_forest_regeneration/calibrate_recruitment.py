#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Calibrate each PFT's recruit carbon efficiency so the model's ingrowth past 1 cm is BCI's.

pft.recruit_carbon_efficiency is the share of a tree's reproduction carbon that becomes recruits at
min_cohort_height: seed set, germination and seedling survival together. Its target is the census:
the new stems of at least 1 cm per hectare and year, by PFT, over the five intervals of 1985-2010
(example 03's recruits.csv.gz).

A trial starts the coupled model of ed_config_regeneration.toml from the 1985 census, runs 20
years (1985-2004), and counts every cohort that crosses 1 cm dbh, weighted by its stems per m2 of
site, day by day; the last ten years are scored, when the seedlings born after the start have had
time to reach 1 cm. Each crossing is kept with its PFT and the LAI of its patch, so the model's
ingrowth can be set against the census's by canopy LAI too.

Recruits shade little, so a PFT's ingrowth is close to linear in its efficiency: one trial with no
reproduction recruits (efficiency 0: the seed rain alone) and one at efficiency E0 give each PFT's
efficiency for the census's ingrowth, and a third trial checks it.

Usage:
  python calibrate_recruitment.py trial --efficiency 0 0 0 --out output/recruitment/seed.npz
  python calibrate_recruitment.py trial --efficiency 1e-3 1e-3 1e-3 --out output/recruitment/e0.npz
  python calibrate_recruitment.py solve output/recruitment/seed.npz output/recruitment/e0.npz
  python calibrate_recruitment.py score output/recruitment/check.npz
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import calibrate_growth as cg             # noqa: E402  (the paths, the census start)
from edmore.config import RunConfig       # noqa: E402

RECRUITS = os.path.join(cg.EX03, "data", "recruits.csv.gz")
WORK = os.path.join(HERE, "output", "recruitment")
START, YEARS, SCORED = 1985, 20, 10       # model years 1985-2004; the last 10 scored
LAI_EDGES = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 99.0])
DBH_INGROWTH = 1.0                        # [cm] the census's smallest stem


def census_ingrowth():
    """New stems >= 1 cm [1/ha/yr] by PFT over 1985-2010, and by LAI within 20 m (LAI_EDGES)."""
    import pandas as pd
    r = pd.read_csv(RECRUITS).dropna(subset=["dt"])
    total = np.array([r.n_new[r.pft == p].sum() / r.exposure[r.pft == p].sum() * 1e4 for p in (1, 2, 3)])
    li = np.digitize(r.lai, LAI_EDGES) - 1
    by_lai = np.array([[r.n_new[(r.pft == p) & (li == j)].sum()
                        / max(r.exposure[(r.pft == p) & (li == j)].sum(), 1e-9) * 1e4
                        for j in range(len(LAI_EDGES) - 1)] for p in (1, 2, 3)])
    return total, by_lai


def trial_config(efficiency, threads, fitted=False, seed_rain=None):
    """The example's model from the 1985 census for YEARS years, with these recruit efficiencies;
    with `fitted`, the growth calibration's best parameters in place of the PFT file's."""
    cfg = RunConfig.load(cg.CONFIG)
    if fitted:
        cg.apply(cfg, cg.best_parameters())
    for key, value in {"run.start_time": f"{START}-01-01", "run.end_time": f"{START + YEARS}-01-01",
                       "run.n_threads": threads, "init.init_mode": 1, "init.census_file": cg.CENSUS,
                       "output.enabled": False, "output.state_write": False}.items():
        cfg.set(key, value)
    for p, e in enumerate(efficiency, start=1):
        cfg.set("pft.recruit_carbon_efficiency", float(e), file="pft", pft=p)
        if seed_rain is not None:
            cfg.set("pft.seed_rain_recruits", float(seed_rain[p - 1]), file="pft", pft=p)
    return cfg


def run_trial(cfg, out):
    """Every crossing of 1 cm: (year, PFT, stems per m2 of site, patch LAI) -> `out` (.npz); and the
    site's area-weighted area by patch-LAI class each year, the exposure the rates are per."""
    from edmore.model import Run
    folder = os.path.splitext(out)[0]
    main = cfg.write(folder)
    os.chdir(folder)
    events, exposure = [], np.zeros((YEARS, len(LAI_EDGES) - 1))
    with Run(main, verbose=False) as run:
        prev = {}
        for step in run:
            c = run.cohorts("global_id", "pft", "nplant", "dbh", "leaf_area", "owner_patch")
            area = run.patches("share")["share"]
            ip = c["owner_patch"] - 1
            lai = np.zeros(len(area))
            np.add.at(lai, ip, c["nplant"] * c["leaf_area"])
            year = step.date.year - START - (1 if (step.date.month, step.date.day) == (1, 1) else 0)
            if 0 <= year < YEARS:
                exposure[year] += np.bincount(np.digitize(lai, LAI_EDGES) - 1, weights=area,
                                              minlength=len(LAI_EDGES) - 1)[:len(LAI_EDGES) - 1] / 365.25
            for i, g in enumerate(c["global_id"]):
                if c["dbh"][i] >= DBH_INGROWTH and prev.get(int(g), DBH_INGROWTH) < DBH_INGROWTH:
                    events.append((year, c["pft"][i], area[ip[i]] * c["nplant"][i], lai[ip[i]]))
            prev = dict(zip(c["global_id"].tolist(), c["dbh"].tolist()))
        run.finalize()
    np.savez(out, events=np.array(events).reshape(-1, 4), exposure=exposure,
             efficiency=np.array(cfg.get("pft.recruit_carbon_efficiency", file="pft")),
             seed_rain=np.array(cfg.get("pft.seed_rain_recruits", file="pft")))


def ingrowth(path):
    """Model ingrowth past 1 cm [1/ha/yr] by PFT over the scored years, and by patch LAI."""
    z = np.load(path)
    e, expo = z["events"], z["exposure"][-SCORED:].sum(0)        # [yr] of site area per LAI class
    keep = e[:, 0] >= YEARS - SCORED
    total = np.array([e[keep & (e[:, 1] == p), 2].sum() / SCORED * 1e4 for p in (1, 2, 3)])
    li = np.digitize(e[:, 3], LAI_EDGES) - 1
    by_lai = np.array([[e[keep & (e[:, 1] == p) & (li == j), 2].sum() / max(expo[j], 1e-9) * 1e4
                        for j in range(len(LAI_EDGES) - 1)] for p in (1, 2, 3)])
    return total, by_lai, expo / SCORED


def report(path):
    target, target_lai = census_ingrowth()
    total, by_lai, share = ingrowth(path)
    print(f"{os.path.basename(path)}: efficiency {np.load(path)['efficiency']}")
    print("  ingrowth past 1 cm [1/ha/yr]  model " + " ".join(f"{v:7.1f}" for v in total)
          + "   census " + " ".join(f"{v:6.1f}" for v in target))
    print("  by LAI class " + " ".join(f"[{a:g},{b:g})" for a, b in zip(LAI_EDGES[:-1], LAI_EDGES[1:]))
          + f"   (model site area in each: {' '.join(f'{v:.2f}' for v in share)})")
    for p in (1, 2, 3):
        print(f"   PFT {p} model  " + " ".join(f"{v:7.1f}" for v in by_lai[p - 1])
              + "   census " + " ".join(f"{v:6.1f}" for v in target_lai[p - 1]))
    return total


def solve(seed_path, e0_path):
    """Each PFT's efficiency for the census's ingrowth, from the seed-rain-only and E0 trials."""
    target, _ = census_ingrowth()
    seed, _, _ = ingrowth(seed_path)
    e0_total, _, _ = ingrowth(e0_path)
    e0 = np.load(e0_path)["efficiency"]
    per_unit = (e0_total - seed) / e0
    efficiency = np.where(per_unit > 0, (target - seed) / per_unit, np.nan)
    for p in (1, 2, 3):
        print(f"PFT {p}: census {target[p-1]:6.1f}; seed rain alone {seed[p-1]:6.1f}; at {e0[p-1]:.1e} "
              f"{e0_total[p-1]:6.1f} -> efficiency {efficiency[p-1]:.3e}"
              + ("   (the seed rain alone exceeds the census)" if target[p-1] < seed[p-1] else ""))
    os.makedirs(WORK, exist_ok=True)
    with open(os.path.join(WORK, "solved.json"), "w") as fh:
        json.dump({"efficiency": [float(v) for v in efficiency], "census": target.tolist(),
                   "seed_rain_only": seed.tolist(), "at_e0": e0_total.tolist(), "e0": e0.tolist()}, fh, indent=1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("trial")
    t.add_argument("--efficiency", type=float, nargs=3, required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--threads", type=int, default=8)
    t.add_argument("--fitted", action="store_true", help="take the growth calibration's best parameters")
    t.add_argument("--seed-rain", type=float, nargs=3, help="[plant/m2/yr] in place of the PFT file's")
    s = sub.add_parser("solve")
    s.add_argument("seed")
    s.add_argument("e0")
    r = sub.add_parser("score")
    r.add_argument("trials", nargs="+")
    args = ap.parse_args(argv)
    if args.cmd == "trial":
        run_trial(trial_config(args.efficiency, args.threads, args.fitted, args.seed_rain),
                  os.path.abspath(args.out))
    elif args.cmd == "solve":
        solve(args.seed, args.e0)
    else:
        for path in args.trials:
            report(path)


if __name__ == "__main__":
    main()
