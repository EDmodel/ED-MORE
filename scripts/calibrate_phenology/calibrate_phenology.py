#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""calibrate_phenology -- fit a stand's leaf phenology to observations, with each cohort's own cues.

The cues a cohort reads in the coupled model (the canopy-air temperature around it, the hours of light
at its top) depend on the canopy, which depends on the phenology: the fit alternates. A coupled run of
the frozen stand gives every cohort's daily cues; the offline stand (the model's own kernel and leaf
rule, every cohort a column) is fitted on them; a new coupled run with the fitted values gives new
cues; until the values settle. README.md has the method and the worked example.

    calibrate_phenology.py iterate --config phenology.toml --work runs --workers 40
    calibrate_phenology.py run     --config phenology.toml --work runs            # one coupled run
    calibrate_phenology.py fit     --config phenology.toml --work runs --workers 40   # one fit on it
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
#----- the edmore package of this source tree: here, and in the coupled runs this starts
sys.path.insert(0, str(ROOT / "python"))
os.environ["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT / "python"), os.environ.get("PYTHONPATH")]))

import pandas as pd                    # noqa: E402

import coupled                         # noqa: E402
import fit as fitting                  # noqa: E402
import offline                         # noqa: E402
import pft_file                        # noqa: E402
import settings as settings_mod        # noqa: E402
from edmore.config import RunConfig    # noqa: E402
from targets import Target             # noqa: E402


def base_config(s) -> RunConfig:
    cfg = RunConfig.load(s["run"]["config"])
    coupled.check_base(cfg)
    return cfg


def with_values(cfg: RunConfig, keys, values: dict) -> RunConfig:
    out = cfg.copy()
    for k in keys:
        out.set(settings_mod.pft_file_key(k.name), float(values[k.label]), file="pft", pft=k.pft)
    return out


def drivers_for(s, cfg: RunConfig, values: dict, rundir: Path) -> Path:
    """The drivers of a coupled run with these values: made, or found when a finished run in
    `rundir` had the same configuration (its main and PFT files, the fitted values among them: a
    change anywhere else, such as the stand's leaf area, makes a new run)."""
    main = coupled.build(cfg, s["keys"], values, rundir, s["run"]["threads"])
    digest = coupled.config_digest(main)
    path, stamp = rundir / "drivers.npz", rundir / "config.sha256"
    if path.exists() and stamp.exists() and stamp.read_text().strip() == digest:
        print(f"  {rundir.name}: reusing the coupled run")
        return path
    print(f"  {rundir.name}: the coupled run ...", flush=True)
    coupled.run(main, s["run"]["runner"])
    coupled.read(rundir / "out", values, coupled.restart_file(cfg)).save(path)
    stamp.write_text(digest + "\n")
    return path


def stand_params(cfg: RunConfig, values: dict, keys, drivers: coupled.Drivers) -> dict:
    run_cfg = with_values(cfg, keys, values)
    return {p: offline.pft_params(run_cfg, p) for p in sorted(set(drivers.pft.tolist()))}


def check_closure(s, stand: offline.Stand) -> dict:
    """Stop when the offline stand does not reproduce the coupled run it came from."""
    worst = offline.closure(stand)
    for p, (diff, day) in worst.items():
        print(f"  closure, PFT {p}: the offline leaf cover is within {diff:.4f} of the coupled run's (worst {day})")
    bad = {p: w for p, w in worst.items() if w[0] > s["fit"]["closure"]}
    if bad:
        raise SystemExit(f"closure fails (> {s['fit']['closure']}) for PFT(s) {sorted(bad)}: the offline stand is "
                         "not the coupled run. A carbon-limited flush (too little storage) or a cohort that "
                         "changed would do this; see README.md, The closure check.")
    return {str(p): w[0] for p, w in worst.items()}


def one_fit(s, cfg, values, rundir, workers) -> dict:
    """The coupled run with `values`, its closure check, and the fit on its drivers."""
    path = drivers_for(s, cfg, values, rundir)
    d = coupled.Drivers.load(path)
    base = stand_params(cfg, values, s["keys"], d)
    closure = check_closure(s, offline.Stand(d, base))
    spec = fitting.Spec(str(path), base, s["keys"], s["targets"], s["run"]["score_from"])
    print(f"  fitting {len(s['keys'])} keys on {len(d.pft)} cohorts x {len(d.dates)} days ...", flush=True)
    res = fitting.fit(spec, values, s, workers)
    res["closure"] = closure
    report(s, res)
    (rundir / "fit.json").write_text(json.dumps(res, indent=1, default=float))
    write_table(spec, res["values"], rundir / "targets.csv")
    return res


def report(s, res) -> None:
    print(f"  {res['message']} ({res['nit']} generations, {res['nfev']} trials)")
    print(f"  {'target':28s} {'start':>9s} {'fitted':>9s}")
    for name, v in res["scores"].items():
        print(f"  {name:28s} {res['scores_start'][name]:9.3f} {v:9.3f}")
    print(f"  {'key':28s} {'start':>11s} {'fitted':>11s}")
    for k in s["keys"]:
        print(f"  {k.label:28s} {res['start'][k.label]:11.5g} {res['values'][k.label]:11.5g}")


def write_table(spec, values, path: Path) -> None:
    """The observations beside the model with these values."""
    p = fitting.Problem(spec)
    cover, loss = p.run(values)
    pd.concat([t.table(cover, loss) for t in p.targets]).to_csv(path, index=False)


def moved(keys, a: dict, b: dict) -> float:
    return max(abs(a[k.label] - b[k.label]) / (k.hi - k.lo) for k in keys)


def best_coupled(trials: list) -> tuple:
    """Of (values, coupled score, run) triples -- each the values a coupled run was made with and that
    run's score of them -- the one whose own run scored best."""
    return min(trials, key=lambda t: t[1])


def cmd_iterate(a, s) -> None:
    cfg = base_config(s)
    values = coupled.pft_values(cfg, s["keys"])
    work = Path(a.work).resolve()
    history, trials, drivers_run = [], [], None
    for i in range(int(s["fit"]["iterations"])):
        print(f"iteration {i + 1}")
        res = one_fit(s, cfg, values, work / f"run{i}", a.workers)
        coupled_score = float(sum(res["scores_start"].values()))     # the run's own values, on its own cues
        trials.append((values, coupled_score, f"run{i}"))
        step = moved(s["keys"], res["values"], values)
        history.append({"run": f"run{i}", "coupled_score": coupled_score, "score": res["score"], "moved": step,
                        "closure": res["closure"]})
        print(f"  the largest move: {step:.3f} of a key's range")
        settled = i > 0 and step < s["fit"]["settle"]
        values = res["values"]
        if settled:
            print("settled")
            break
    else:
        #----- A cue the canopy feeds back on (WATER: the leaves dry the soil whose potential is the cue) can
        #      keep the fits from settling: each overshoots on the last run's cues. The values kept are then
        #      those whose own coupled run scored best, not the last fit's, which no run has tried.
        values, best, drivers_run = best_coupled(trials)
        print(f"not settled after {s['fit']['iterations']} fits: the values written are those {drivers_run} was "
              f"made with, whose own run scored best ({best:.3f})")
    #----- drivers_run: the coupled run made with the values written, when one was (an unsettled fit)
    (work / "iterations.json").write_text(json.dumps({"history": history, "values": values,
                                                      "drivers_run": drivers_run}, indent=1))
    out = pft_file.write(Path(cfg.pft_path), work / "pft_parameters_fitted.toml", s["keys"], values,
                         f"From {s['path'].name}: {len(history)} fit(s); see {work.name}/iterations.json.")
    print(f"wrote {out}")


def cmd_run(a, s) -> None:
    cfg = base_config(s)
    values = json.loads(Path(a.values).read_text())["values"] if a.values else coupled.pft_values(cfg, s["keys"])
    rundir = Path(a.work).resolve() / "run"
    path = drivers_for(s, cfg, values, rundir)
    d = coupled.Drivers.load(path)
    check_closure(s, offline.Stand(d, stand_params(cfg, values, s["keys"], d)))
    for t in s["targets"]:
        Target(t, d, s["run"]["score_from"])           # each target finds its cohorts and observations


def cmd_fit(a, s) -> None:
    cfg = base_config(s)
    values = json.loads(Path(a.values).read_text())["values"] if a.values else coupled.pft_values(cfg, s["keys"])
    one_fit(s, cfg, values, Path(a.work).resolve() / "run", a.workers)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=("iterate", "run", "fit"))
    ap.add_argument("--config", required=True, help="the site's phenology calibration file")
    ap.add_argument("--work", required=True, help="the directory for the runs and the results")
    ap.add_argument("--workers", type=int, default=1, help="processes for the fit's trials")
    ap.add_argument("--values", help="run, fit: a fit.json whose values to use instead of the PFT file's")
    a = ap.parse_args(argv)
    try:
        s = settings_mod.load(a.config)
    except settings_mod.SettingsError as e:
        raise SystemExit(f"{a.config}: {e}")
    {"iterate": cmd_iterate, "run": cmd_run, "fit": cmd_fit}[a.command](a, s)


if __name__ == "__main__":
    main()
