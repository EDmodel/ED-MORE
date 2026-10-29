#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Check, and calibrate, the PFTs' growth against BCI's census.

What pft_parameters.toml sets from data, and why, is in its comments: wood, size and mortality are
example 03's census fits; the leaves' top-of-canopy turnover rate (1 / their lifespan) and mass per
area are Panama's sun leaves', with Ma et al. (2025)'s light plasticity; fine-root carbon equals leaf
carbon; stem respiration is charged on sapwood volume at a rate set from stand totals in the
literature; the fine roots' own respiration, the share of growth put into reproduction and the sink
limit on diameter growth are set from field data too.

`check` runs that file as it is against the census and writes calibration.json. The search fits three
parameters, shared by the PFTs: the height at which trees start to reproduce (only trees above it
pay for reproduction, so it shapes growth by size), fine-root turnover as a multiple of the leaves'
(leaf_turnover_rate_toc), and one scale on the three PFTs' Vcmax25 (the growth level; the plasticity
slopes follow it). The shipped onset height, 5.27 m, is its fit; root turnover (1.5 x the leaves') and
Vcmax25 (example 04's, a scale of 1) are set from data instead.

The target is BCI's own growth over the census interval 2005-2010: each tree's dbh growth [cm/yr],
by PFT, size D and overtopping LAI L at the start of the interval, the trees that died within it
included with their last measured growth (example 03's growth table), averaged per class as
measured -- negative increments included, since they are measurement error that averages out. A
trial starts the coupled model of ed_config_regeneration.toml from the 2005 census (example 03's
census file), runs the same five years on ERA5-Land's weather, and measures each cohort's growth the
census way: its dbh change from the start to the end of the run (or to the last month it was seen),
over that time, credited to its D, L and stems per m2 of site at the start. A PFT's score is the mean
squared log-ratio of the model's class mean growth to the census's, over size classes and over
size-and-L classes (the two averaged). L is binned coarsely (below 3, 3-6, above 6): the census's L,
the leaf area of taller trees within 20 m, is a noisy measure of a tree's light, which flattens the
census's response to it. The loss is the PFTs' sum plus a penalty ((x - mean) / sd)^2 for each stand
total in TARGETS, from the trial's years 2-5.

The search: stage 1 is a Latin hypercube of 60 trials over the box, stage 2 one of 40 in a box an
eighth as wide around stage 1's best, and the final trial the best of both. Each trial takes about
12 minutes on 8 threads.

Usage:
  python calibrate_growth.py design --stage 1         # writes output/calibration/stage1/tasks.json
  python calibrate_growth.py run --stage 1 --workers 5 --threads 8   # its trials, here
  python calibrate_growth.py trial --stage 1 --task 7                # one trial (a cluster array task)
  python calibrate_growth.py fit --stage 1            # the stage's best -> best.json
  (the same for --stage 2, which designs around stage 1's best)
  python calibrate_growth.py final                    # the best trial, once -> calibration.json
  python calibrate_growth.py check                    # pft_parameters.toml as it is -> calibration.json

Needs example 03's census files: python ../example03_demography/prepare_census.py CENSUS_DIR
"""
import argparse
import calendar
import json
import os
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
EX03 = os.path.join(ROOT, "examples", "example03_demography")
sys.path.insert(0, os.path.join(ROOT, "python"))
from edmore.config import RunConfig         # noqa: E402

CONFIG = os.path.join(HERE, "ed_config_regeneration.toml")
CENSUS = os.path.join(EX03, "data", "bci_1985_census.csv")          # calibrate_recruitment.py's start
GROWTH_CENSUS = os.path.join(EX03, "data", "bci_2005_census.csv")   # the growth check's start
GROWTH_TABLE = os.path.join(EX03, "data", "growth.csv.gz")          # the census's growth, tree by tree
WORK = os.path.join(HERE, "output", "calibration")
RESULT = os.path.join(HERE, "calibration.json")

#: The fitted parameters and their boxes.
FIT = {"min_repro_height": (5.0, 30.0),             # [m] the height a cohort reproduces from
       "fineroot_turnover": (1.5, 5.0),          # [--] fine-root turnover / leaf_turnover_rate_toc
       "vcmax_scale": (0.7, 1.2)}                # [--] on the PFT file's Vcmax25
#: Stand totals of the trial's years 2-5, (mean, sd) each, added to the loss as ((x - mean) / sd)^2.
TARGETS = {"fineroot": (0.27, 0.11),    # fine roots' share of NPP: Malhi, Doughty & Galbraith (2011), 35 sites
           "gpp": (3.15, 0.30),         # [kgC/m2/yr] BCI's GPP (GOSIF; the tower's 2.8 is likely low)
           "root_resp": (0.56, 0.10),   # [kgC/m2/yr] roots + rhizosphere, Gigante and Manaus (pft_parameters.toml)
           "repro_c": (0.05, 0.025)}    # [kgC/m2/yr] flowers and fruit, BCI's litter traps (pft_parameters.toml)
N_TRIALS = {1: 60, 2: 40}

START, YEARS = 2005, 5               # the census interval 2005-2010 on ERA5-Land's own weather
D_EDGES = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 400.0])   # [cm] size classes (example 03's)
L_EDGES = np.array([0.0, 3.0, 6.0, 99.0])                              # overtopping LAI classes
EPS = 0.02                           # [cm/yr] keeps the log-ratio finite where growth is ~0
MIN_N = 5                            # cohorts a class needs to count
MIN_TREES = 20                       # census trees a class needs to count

with open(os.path.join(HERE, "pft_parameters.toml"), "rb") as _fh:
    _PFT_FILE = tomllib.load(_fh)
AGF = _PFT_FILE["allometry"]["aboveground_frac"][0]   # the coarse roots' share of wood is 1 - AGF
GRF = _PFT_FILE["pft"]["growth_resp_factor"][0]       # construction cost, charged on growth


def plasticity_slopes(vcmax25, rd_vcmax_ratio):
    """Ma et al. (2025)'s light-plasticity slopes of Vcmax25 and Rd25 [per unit overtopping LAI]."""
    return ([-(0.00242 * v + 0.06212) for v in vcmax25],
            [-(0.18974 * r * v + 0.11744) for v, r in zip(vcmax25, rd_vcmax_ratio)])


def apply(cfg, par):
    """Set the fitted parameters in a RunConfig that holds the PFT file's values."""
    vc = [v * par["vcmax_scale"] for v in cfg.get("pft.vcmax25", file="pft")]
    k_vc, k_rd = plasticity_slopes(vc, cfg.get("pft.rd_vcmax_ratio", file="pft"))
    cfg.set("pft.vcmax25", vc, file="pft")
    cfg.set("pft.kplastic_vcmax", k_vc, file="pft")
    cfg.set("pft.kplastic_rd", k_rd, file="pft")
    cfg.set("pft.fineroot_turnover_rate",
            [par["fineroot_turnover"] * r for r in cfg.get("pft.leaf_turnover_rate_toc", file="pft")], file="pft")
    cfg.set("pft.min_repro_height", [par["min_repro_height"]] * len(cfg.get("pft.wood_density", file="pft")),
            file="pft")
    return cfg


def trial_config(par, threads):
    """The example's model from the 2005 census, five years, yearly site totals only; with `par`
    (None: the PFT file as it is), the fitted parameters set."""
    cfg = RunConfig.load(CONFIG)
    for key, value in {"run.start_time": f"{START}-01-01", "run.end_time": f"{START + YEARS}-01-01",
                       "run.n_threads": threads, "init.init_mode": 1, "init.census_file": GROWTH_CENSUS,
                       "output.enabled": True, "output.prefix": "trial", "output.monthly.enabled": False,
                       "output.state_write": False}.items():
        cfg.set(key, value)
    return apply(cfg, par) if par else cfg


def run_trial(cfg, out):
    """Run `cfg` and measure every starting cohort's growth the census way -> `out` (.npz): rows
    (0, dbh, growth, overtopping LAI, stems per m2 of site, PFT), dbh, LAI and stems at the start and
    growth [cm/yr] = its dbh change to the last month it was seen, over that time (at least a year);
    and the carbon budget of years 2-5. A cohort's record stops before it fuses with another (its
    density rises): the fused cohort's dbh is the two cohorts' mean, a change that is not growth."""
    from edmore.model import Run
    folder = os.path.splitext(out)[0]
    main = cfg.write(folder)
    os.chdir(folder)
    cfg.remove_diagnostics(folder)                     # an earlier attempt's files: edmore never overwrites output
    start, last, density, fused = None, {}, {}, set()
    with Run(main, verbose=False) as run:
        for step in run:
            if start is None or step.date.day == 1:
                c = run.cohorts("global_id", "dbh", "overtopping_lai", "nplant", "pft", "owner_patch")
                t = step.date.toordinal() / 365.25
                area = run.patches("share")["share"][c["owner_patch"] - 1]
                if start is None:
                    start = {int(i): (d, l, a * n, p, t) for i, d, l, a, n, p in
                             zip(c["global_id"], c["dbh"], c["overtopping_lai"], area, c["nplant"], c["pft"])}
                for i, d, a, n in zip(c["global_id"], c["dbh"], area, c["nplant"]):
                    i = int(i)
                    if i in fused:
                        continue
                    if a * n > density.get(i, np.inf) * (1 + 1e-9):     # its stems rose: it fused
                        fused.add(i)
                        continue
                    density[i], last[i] = a * n, (d, t)
        run.finalize()
    rows = [(0.0, d0, (last[i][0] - d0) / (last[i][1] - t0), l0, w0, p)
            for i, (d0, l0, w0, p, t0) in start.items() if last[i][1] - t0 >= 1.0]
    np.savez(out, samples=np.array(rows), **carbon_budget(os.path.join(folder, "output", "trial-Y.nc")))


#: the stand totals a trial keeps, the mean of its years 2-5 [kgC/m2/yr]
BUDGET = ("mean_gpp_py", "mean_leaf_resp_py", "mean_stem_resp_py", "mean_root_resp_py", "mean_growth_resp_py",
          "mean_npp_leaf_py", "mean_npp_fineroot_py", "mean_npp_wood_py", "mean_npp_repro_py", "mean_npp_storage_py",
          "mean_root_exudate_py")


#: the fast loop's carbon fluxes, which the files give as mean rates [umol CO2 m-2 s-1]
UMOL_RATES = ("mean_gpp_py", "mean_leaf_resp_py", "mean_stem_resp_py", "mean_root_resp_py")
UMOL_TO_KGC = 1.20107e-8                              # [kgC/umol CO2], ed_constants' umol_2_kgc


def carbon_budget(path):
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        ds.set_auto_mask(False)
        yr = np.array([366.0 if calendar.isleap(int(y)) else 365.0 for y in ds["year"][:]]) * 86400.0

        def per_year(k):                              # [kgC/m2/yr]: a mean rate times its year's length
            x = np.asarray(ds[k][:], dtype=float)
            return x * UMOL_TO_KGC * yr.reshape((-1,) + (1,) * (x.ndim - 1)) if k in UMOL_RATES else x
        return {k: float(np.mean(per_year(k)[1:])) for k in BUDGET}


def shares(b):
    """NPP net of all plant respiration (as Malhi et al. measure it), its use efficiency, and the
    shares of reproduction and fine roots in it; GPP; and the roots' autotrophic respiration
    [kgC/m2/yr]: the fine roots' maintenance, the coarse roots' share of stem respiration, and the
    construction cost of fine- and coarse-root growth."""
    npp = sum(float(b[k]) for k in ("mean_npp_leaf_py", "mean_npp_fineroot_py", "mean_npp_wood_py",
                                    "mean_npp_repro_py", "mean_npp_storage_py"))
    g = float(b["mean_gpp_py"])
    root_resp = (float(b["mean_root_resp_py"]) + (1 - AGF) * float(b["mean_stem_resp_py"])
                 + GRF * (float(b["mean_npp_fineroot_py"]) + (1 - AGF) * float(b["mean_npp_wood_py"])))
    return {"cue": npp / g, "repro": float(b["mean_npp_repro_py"]) / npp, "repro_c": float(b["mean_npp_repro_py"]),
            "fineroot": float(b["mean_npp_fineroot_py"]) / npp, "gpp": g, "root_resp": root_resp,
            "stem_resp_above": AGF * float(b["mean_stem_resp_py"]), "leaf_npp": float(b["mean_npp_leaf_py"]),
            "wood_npp_above": AGF * float(b["mean_npp_wood_py"]), "exudate": float(b["mean_root_exudate_py"])}


# ----- scoring ----------------------------------------------------------------------------------------
def census_growth():
    """The census's growth over the interval the trial runs: dbh, overtopping LAI, PFT and growth
    [cm/yr] of every tree, as measured."""
    import pandas as pd
    g = pd.read_csv(GROWTH_TABLE)
    return g[g.interval == START]


def binned(dbh, lai, g, w, census, l_edges):
    """Rows (D class, L class, cohorts, model mean, census mean) of the classes with at least MIN_N
    cohorts and MIN_TREES census trees."""
    di, li = np.digitize(dbh, D_EDGES) - 1, np.digitize(lai, l_edges) - 1
    ci = np.digitize(census.dbh, D_EDGES) - 1
    cj = np.digitize(census.lai_over, l_edges) - 1
    rows = []
    for i in range(len(D_EDGES) - 1):
        for j in range(len(l_edges) - 1):
            b, c = (di == i) & (li == j), (ci == i) & (cj == j)
            if b.sum() >= MIN_N and c.sum() >= MIN_TREES:
                rows.append((i, j, b.sum(), np.average(g[b], weights=w[b]), census.g[c].mean()))
    return np.array(rows)


def score(samples, census=None):
    """Per PFT: the classes by size ("size") and by size and L ("light"), their mean squared
    log-ratios, the mean log-ratio by size ("bias"), and the loss (the two scores averaged)."""
    census = census_growth() if census is None else census
    q, dbh, g, lai, w, pft = samples.T
    out = {}
    for p in (1, 2, 3):
        m = (dbh >= D_EDGES[0]) & (pft == p)
        cp = census[census.pft == p]
        size = binned(dbh[m], lai[m], g[m], w[m], cp, np.array([0.0, 99.0]))
        light = binned(dbh[m], lai[m], g[m], w[m], cp, L_EDGES)
        #----- a starving cohort can lose wood, so a class mean can fall below zero: read it as none
        r_size = np.log((np.maximum(size[:, 3], 0) + EPS) / (np.maximum(size[:, 4], 0) + EPS))
        r_light = np.log((np.maximum(light[:, 3], 0) + EPS) / (np.maximum(light[:, 4], 0) + EPS))
        out[p] = {"size": size, "light": light, "bias": float(np.mean(r_size)),
                  "loss_size": float(np.mean(r_size ** 2)), "loss_light": float(np.mean(r_light ** 2))}
        out[p]["loss"] = 0.5 * (out[p]["loss_size"] + out[p]["loss_light"])
    return out


# ----- the search -------------------------------------------------------------------------------------
def stage_dir(stage):
    return os.path.join(WORK, f"stage{stage}")


def design(stage):
    """Write the stage's trials: a Latin hypercube over the box (stage 1) or around stage 1's best."""
    from scipy.stats import qmc
    lo = np.array([b[0] for b in FIT.values()])
    hi = np.array([b[1] for b in FIT.values()])
    if stage > 1:
        with open(os.path.join(stage_dir(stage - 1), "best.json")) as fh:
            best = json.load(fh)["par"]
        mid = np.array([best[k] for k in FIT])
        half = (hi - lo) / 16
        lo, hi = np.maximum(mid - half, lo), np.minimum(mid + half, hi)
    points = qmc.scale(qmc.LatinHypercube(d=len(FIT), seed=stage).random(N_TRIALS[stage]), lo, hi)
    tasks = [dict(zip(FIT, map(float, pt))) for pt in points]
    os.makedirs(stage_dir(stage), exist_ok=True)
    with open(os.path.join(stage_dir(stage), "tasks.json"), "w") as fh:
        json.dump(tasks, fh, indent=1)
    print(f"stage {stage}: {len(tasks)} trials -> {stage_dir(stage)}/tasks.json")


def trial(stage, task, threads):
    with open(os.path.join(stage_dir(stage), "tasks.json")) as fh:
        par = json.load(fh)[task - 1]
    run_trial(trial_config(par, threads), os.path.join(stage_dir(stage), f"t{task:03d}.npz"))


def run_stage(stage, workers, threads):
    with open(os.path.join(stage_dir(stage), "tasks.json")) as fh:
        n = len(json.load(fh))
    todo = [i for i in range(1, n + 1) if not os.path.exists(os.path.join(stage_dir(stage), f"t{i:03d}.npz"))]
    cmd = lambda i: [sys.executable, __file__, "trial", "--stage", str(stage), "--task", str(i),
                     "--threads", str(threads)]
    with ThreadPoolExecutor(workers) as pool:
        failed = [i for i, s in zip(todo, pool.map(lambda i: subprocess.call(cmd(i)), todo)) if s]
    if failed:
        raise SystemExit(f"ERROR: stage {stage} trials failed: {failed}")


def total_loss(z, s):
    """The PFTs' growth losses plus the stand totals' penalties."""
    b = shares(z)
    return float(sum(s[p]["loss"] for p in (1, 2, 3))
                 + sum(((b[k] - m) / sd) ** 2 for k, (m, sd) in TARGETS.items()))


def fit(stage):
    """The stage's trials ranked by loss; the best -> best.json."""
    with open(os.path.join(stage_dir(stage), "tasks.json")) as fh:
        tasks = json.load(fh)
    census = census_growth()
    results = []
    for i, par in enumerate(tasks, start=1):
        f = os.path.join(stage_dir(stage), f"t{i:03d}.npz")
        if os.path.exists(f):
            z = np.load(f)
            s = score(z["samples"], census)
            results.append({"task": i, "par": par, "loss": total_loss(z, s),
                            "size": [s[p]["loss_size"] for p in (1, 2, 3)],
                            "light": [s[p]["loss_light"] for p in (1, 2, 3)],
                            "bias": [s[p]["bias"] for p in (1, 2, 3)], "budget": shares(z)})
    print(f"stage {stage}: {len(results)} of {len(tasks)} trials scored")
    results.sort(key=lambda r: r["loss"])
    for r in results[:8]:
        print(f"  t{r['task']:03d} loss {r['loss']:.3f}  size " + " ".join(f"{v:.2f}" for v in r["size"])
              + "  light " + " ".join(f"{v:.2f}" for v in r["light"])
              + "  bias " + " ".join(f"{v:+.2f}" for v in r["bias"])
              + "  | CUE {cue:.2f} repro {repro_c:.3f} fineroot {fineroot:.0%} GPP {gpp:.2f} root resp {root_resp:.2f} leaf NPP {leaf_npp:.2f} wood above {wood_npp_above:.2f} |  ".format(**r["budget"])
              + " ".join(f"{k} {v:.3g}" for k, v in r["par"].items()))
    with open(os.path.join(stage_dir(stage), "best.json"), "w") as fh:
        json.dump(results[0], fh, indent=1)
    return results[0]


def best_parameters():
    """The best fitted parameters over the stages run so far."""
    bests = []
    for d in sorted(os.listdir(WORK)):
        f = os.path.join(WORK, d, "best.json")
        if d.startswith("stage") and os.path.exists(f):
            with open(f) as fh:
                bests.append(json.load(fh))
    return min(bests, key=lambda b: b["loss"])["par"]


def final(threads, par=None):
    """The best trial of all stages (or, with `par` {}, the PFT file as it is), run once and scored
    -> calibration.json."""
    par = best_parameters() if par is None else par
    out = os.path.join(WORK, "final.npz" if par else "check.npz")
    run_trial(trial_config(par, threads), out)
    z = np.load(out)
    s = score(z["samples"])
    cfg = trial_config(par, 1)
    about = ("carbon costs fitted to" if par else "pft_parameters.toml as it is, against")
    result = {"about": f"calibrate_growth.py: {about} BCI's census growth, 2005-2010",
              "fitted": par, "targets": TARGETS, "loss": total_loss(z, s),
              "carbon_budget": {k: round(v, 4) for k, v in shares(z).items()},
              "pft": {k: [round(v, 4) for v in cfg.get(f"pft.{k}", file="pft")]
                      for k in ("vcmax25", "fineroot_turnover_rate")},
              "score": {str(p): {k: round(s[p][k], 4) for k in ("loss", "loss_size", "loss_light", "bias")}
                        for p in (1, 2, 3)},
              # the final trial's mean growth [cm/yr] by size class: dbh from, dbh to, model, census
              "size_classes": {str(p): [[D_EDGES[int(r[0])], D_EDGES[int(r[0]) + 1], round(r[3], 4), round(r[4], 4)]
                                        for r in s[p]["size"]] for p in (1, 2, 3)},
              # and by size and overtopping LAI class: dbh from, dbh to, LAI from, LAI to, model, census
              "light_classes": {str(p): [[D_EDGES[int(r[0])], D_EDGES[int(r[0]) + 1], L_EDGES[int(r[1])],
                                          L_EDGES[int(r[1]) + 1], round(r[3], 4), round(r[4], 4)]
                                         for r in s[p]["light"]] for p in (1, 2, 3)}}
    with open(RESULT, "w") as fh:
        json.dump(result, fh, indent=1)
    print(json.dumps(result, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("design", "run", "trial", "fit", "final", "check"):
        sp = sub.add_parser(name)
        if name not in ("final", "check"):
            sp.add_argument("--stage", type=int, default=1)
        if name in ("run", "trial", "final", "check"):
            sp.add_argument("--threads", type=int, default=8, help="threads per trial")
        if name == "run":
            sp.add_argument("--workers", type=int, default=1, help="trials at once")
        if name == "trial":
            sp.add_argument("--task", type=int, required=True, help="1-based, as a Slurm array index")
    args = ap.parse_args(argv)
    for f in (GROWTH_CENSUS, GROWTH_TABLE):
        if not os.path.exists(f):
            raise SystemExit(f"ERROR: no {f}; run example 03's prepare_census.py first")
    if args.cmd == "design":
        design(args.stage)
    elif args.cmd == "run":
        run_stage(args.stage, args.workers, args.threads)
    elif args.cmd == "trial":
        trial(args.stage, args.task, args.threads)
    elif args.cmd == "fit":
        fit(args.stage)
    elif args.cmd == "final":
        final(args.threads)
    else:
        final(args.threads, par={})


if __name__ == "__main__":
    main()
