#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""calibrate_fast.py -- fit ED-MORE's fast (sub-daily) parameters to a flux tower, with the stand frozen
at its initial structure (the protocol README.md links).

One joint fit (§6.3): Levenberg-Marquardt on every fitted key at once, over ten-day calibration
windows (every target) and seasonal soil-drying runs (LE), each a frozen run restarted from its own
state chain; one refresh of the chains, weights and sigma after the first convergence; then the final
gradient matrix for a rough uncertainty (the Laplace covariance at the MAP), the declared
alternatives, and the validation. The data rules (data_rules.py), the targets' observation models
(observation_models.py) and the priors from the site's climate (priors.py) come first. Every setting
of a calibration is documented, with its default, in calibration_reference.toml.

Commands (each reads a calibration's settings, e.g. examples/example04_column_biophysics/calibration.toml):
  report    the data report only: every target's records through every filter, the u* and closure
            diagnostics, sigma, and the windows and seasonal runs (printed as calibration.toml entries)
  check     the tool and the model's restart on one window: a chain, a trial with its parameter
            record checked, one gradient column that must move the output, a repeated trial that must
            match byte for byte, and a stand unchanged at the trial's end (CTest runs it)
  fit       the data report, the keys' coverage and screening, the joint fit, the uncertainty, the
            validation, report.md and the calibrated configs
  variants  the structural variants' fits side by side: each key's MAP and their spread
  worker    a node's worker for --queue (started inside the Slurm allocation)

Trials run through ED-MORE's Python API (`python -m edmore.model`, which needs libedmore.so) unless
--runner names a edmore executable; the two give the same output, bit for bit. The EEO vcmax25
prior uses ED-MORE's own leaf through libedmore (edmore.plant.leaf). A trial that fails stops the command
with its log.

Usage:
  calibrate_fast.py fit --config calibration.toml --variant interception_off --work runs/off --workers 40
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
#----- the edmore package of this source tree: edmore.config here, edmore.model in the trials this starts
sys.path.insert(0, str(ROOT / "python"))
os.environ["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT / "python"), os.environ.get("PYTHONPATH")]))

import calibrated_files                # noqa: E402
import chains                          # noqa: E402
import data_rules                      # noqa: E402
import fit                             # noqa: E402
import report                          # noqa: E402
import targets                         # noqa: E402
import trials                          # noqa: E402
import uncertainty                     # noqa: E402
from calibration import Calibration    # noqa: E402
from edmore.config import read_record  # noqa: E402
from workers import Task, make_workers, queue_worker   # noqa: E402



def repo_relative(path):
    """The path from the repository's root when it lies inside it, so a fit kept in the repository
    names no one's folders; otherwise the path as it is."""
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(path)

def log_to(path):
    fh = open(path, "a")
    lock = threading.Lock()

    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        with lock:
            print(line, flush=True)
            fh.write(line + "\n")
            fh.flush()
    return log


def open_work(args, log_name):
    """The calibration, its working directory and log, and the workers."""
    cal = Calibration(args.config, args.variant)
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    log = log_to(work / log_name)
    workers = make_workers("queue" if getattr(args, "queue", False) else "local",
                           getattr(args, "workers", 1), work / "queue")
    return cal, work, log, workers


def base_record(cal, work, workers, runner, log):
    """Run one day at the base configs to read the model's parameter record: the defaults of keys the
    base files leave to the model (cal.record)."""
    rdir = work / "record"
    csv = rdir / "out" / f"{trials.PREFIX}_parameters.csv"
    if not csv.exists():
        cfg = trials.with_keys(cal.base, [], [], cal.overrides)
        w = cal.windows[0]
        for k, v in {"run.start_time": trials.stamp(w.start), "run.end_time": trials.stamp(w.start + dt.timedelta(days=1)),
                     "run.slow_on": False, "run.n_threads": 1, "output.state_write": False, "output.enabled": True,
                     "output.dir": str(rdir / "out"), "output.prefix": trials.PREFIX, "output.daily.enabled": False,
                     "output.fast.enabled": False, "output.monthly.enabled": False,
                     "output.annual.enabled": False}.items():
            cfg.set(k, v)
        (rdir / "out").mkdir(parents=True, exist_ok=True)
        main = cfg.write(rdir)
        log("reading the parameter record from a one-day base run")
        status = workers.run([Task("record", trials.command(runner, main), str(rdir), str(rdir / "run.log"),
                                   trials.timeout_for(1, cal.fit_settings["timeout_per_day"]))])["record"][0]
        if status != "ok" or not csv.exists():
            raise trials.TrialError(f"the one-day base run failed: {status}\n{trials.log_tail(rdir / 'run.log')}")
    cal.record = read_record(csv, {rdir / "main.toml": "main", rdir / "pft.toml": "pft"})
    return cal.record


def windows_as_toml(cal) -> str:
    """The windows and seasonal runs as calibration.toml entries (to freeze the rules' choice)."""
    out = [f'[[windows.list]]\nname = "{w.name}"\nstart = "{w.start:%Y-%m-%d}"\nrole = "{w.role}"\n' for w in cal.windows]
    out += [f'[[seasonal_runs.list]]\nname = "{w.name}"\nstart = "{w.start:%Y-%m-%d}"\ndays = {w.days}\n'
            for w in cal.seasonal_runs]
    return "\n".join(out)


# ------------------------------------------------------------------------------------------------
# report, check
# ------------------------------------------------------------------------------------------------
def cmd_report(args):
    cal, work, log, workers = open_work(args, "report.log")
    cal.load_data()
    rep = report.data_report(cal, cal.rows_for(cal.fit_windows()), log)
    report.save_json(work / "data_report.json", rep)
    print(windows_as_toml(cal))
    workers.close()
    return 0


def cmd_check(args):
    """The tool and the model's restart on the first calibration window (module docstring)."""
    cal, work, log, workers = open_work(args, "check.log")
    cal.load_data()
    w = next(x for x in cal.windows if x.role == "cal")
    base_record(cal, work, workers, args.runner, log)
    keys = cal.keys()
    start = np.array([p.centre for p in keys])
    per_day = float(cal.fit_settings["timeout_per_day"])
    max_threads = min(int(cal.fit_settings["max_threads"]), workers.capacity)
    states = chains.run_chains([w], cal.chain_lead_days, cal.base, keys, start, work / "chains", workers,
                               args.runner, cal.overrides, per_day, log, max_threads)
    runner = trials.TrialRunner(keys, [w], cal.rows_for([w]), states, cal.base, cal.overrides, args.runner, workers,
                                work / "trials", cal.step, per_day, cal.fixed_observation_keys(keys), log=log,
                                max_threads=max_threads)
    name = args.key or keys[0].name
    fk = fit.FreeKeys(runner, [[p.name for p in keys].index(name)], start.copy())
    u0 = fk.u_prior
    r0 = fk.residuals([u0])[0]
    J, _ = fit.gradient_central(fk, u0, r0)
    moved = bool(np.any(J[:, 0] != 0.0))
    log(f"the gradient column of {name} moves the output: {moved} (|J| = {np.linalg.norm(J[:, 0]):.4g})")

    def run_again(subdir, write_state=False):
        threads = trials.threads_for(w.days, max_threads)
        tdir = trials.build_trial(cal.base, keys, start, w, states[w.name], work / subdir, cal.overrides, write_state,
                                  threads)
        status = workers.run([Task(subdir, trials.command(args.runner, tdir / "main.toml"), str(tdir),
                                   str(tdir / "run.log"), trials.timeout_for(w.days, per_day), threads)])[subdir][0]
        if status != "ok":
            raise trials.TrialError(f"{tdir.name}: {status}\n{trials.log_tail(tdir / 'run.log')}")
        return tdir, trials.finish(tdir, keys, start, cal.step)
    _, again = run_again("repeat")
    first = trials.load_series(trials.build_trial(cal.base, keys, start, w, states[w.name], work / "trials", cal.overrides))
    same = all(np.array_equal(first[v].values, again[v].values) for v in trials.TRIAL_VARIABLES)
    log(f"a repeated trial is byte-identical: {same}")
    sdir, _ = run_again("stand", write_state=True)
    differ = trials.stand_unchanged(states[w.name], sorted((sdir / "out").glob(f"{trials.PREFIX}-S-*.nc"))[-1])
    log(f"the stand is unchanged at the trial's end: {not differ}" + (f" (differ: {differ})" if differ else ""))
    for t, sc in targets.target_scores([runner.rows[w.name]], r0).items():
        log(f"  {t:10s} nrmse {sc['nrmse']:.3f}  n = {sc['n']}")
    workers.close()
    return 0 if (moved and same and not differ) else 1


# ------------------------------------------------------------------------------------------------
# fit
# ------------------------------------------------------------------------------------------------
def refresh_weights(runner, values, rows_list, log, scale_sigma_cap=None):
    """At these values: each target's effective-sample weight and, with scale_sigma_cap, each target's
    sigma scaled by the model's misfit (targets.scale_sigma). Returns (weights, scales)."""
    for wr in rows_list:
        for t in wr.targets:
            t.weight = 1.0
    r = runner.residuals([values])[0]
    scales = None
    if scale_sigma_cap is not None:                  # every weight is 1 here: the scale reads plain chi^2
        scales = targets.scale_sigma(rows_list, r, scale_sigma_cap)
        log(f"sigma scaled by the model's misfit (max(1, sqrt(chi2 per row)), at most {scale_sigma_cap}): "
            + ", ".join(f"{k} {v:.2f}" for k, v in scales.items()))
    weights = targets.set_weights(rows_list, r)
    log("effective-sample weights n_eff/n: " + ", ".join(
        f"{k} {np.mean([v for kk, v in weights.items() if kk.endswith('/' + k)]):.2f}"
        for k in sorted({kk.split('/')[1] for kk in weights})))
    return weights, scales


def cmd_fit(args):
    """The joint fit (module docstring)."""
    cal, work, log, workers = open_work(args, "fit.log")
    fs = cal.fit_settings
    t_start = time.time()
    cal.load_data()
    fit_windows, val_windows = cal.fit_windows(), cal.validation_windows()
    all_windows = cal.windows + cal.seasonal_runs
    log(f"{f'variant {cal.variant}: ' if cal.variant else ''}{len(fit_windows) - len(cal.seasonal_runs)} "
        f"calibration, {len(val_windows)} validation windows and {len(cal.seasonal_runs)} seasonal runs")
    base_record(cal, work, workers, args.runner, log)
    keys = cal.keys()
    start = np.array([p.centre for p in keys])
    default = np.array([p.default for p in keys])
    per_day = float(fs["timeout_per_day"])
    max_threads = min(int(fs["max_threads"]), workers.capacity)

    def chains_at(values, windows):
        return chains.run_chains(windows, cal.chain_lead_days, cal.base, keys, values, work / "chains", workers,
                                 args.runner, cal.overrides, per_day, log, max_threads)
    rows = cal.rows_for(all_windows)
    runner = trials.TrialRunner(keys, fit_windows, rows, chains_at(start, all_windows), cal.base, cal.overrides,
                                args.runner, workers, work / "trials", cal.step, per_day,
                                cal.fixed_observation_keys(keys), args.keep_netcdf, log, max_threads)
    states_start = dict(runner.states)
    fit_rows = [rows[w.name] for w in fit_windows]
    out = {"variant": cal.variant, "config": repo_relative(cal.path),
           "windows": {w.name: [str(w.start), w.role, w.days] for w in all_windows},
           "registry": {p.name: {"state": p.state, "kind": p.kind, "scope": p.scope, "reason": p.reason}
                        for p in cal.registry()},
           "keys": {p.name: {"kind": p.kind, "scope": p.scope, "range": [p.lo, p.hi], "default": p.default,
                             "prior_centre": p.centre, "prior_sd_u": p.sd_u,
                             "prior_source": p.prior.get("source", p.source)} for p in keys},
           "climate_priors": cal.data.climate_priors,
           "kattge_knorr": [[f, k, i, v] for (f, k, i), v in cal.kattge_knorr.items()]}
    cp = cal.data.climate_priors or {}
    if isinstance(cp.get("kattge_knorr"), dict):
        log(f"Kattge & Knorr at the growth temperature {cp['climate']['t_growth_c']:.1f} C: {cp['kattge_knorr']}")
    for k, e in (cp.get("eeo") or {}).items():
        log(f"EEO prior centre of {k}: {e['eeo']:.4g}" + (f" (the {fs['plant_type']} meta-analysis: "
            f"{e['meta']['centre']:.4g}, {e['meta_apart_sd']:.1f} sd apart{' -- FLAGGED' if e['meta_apart_sd'] > 2 else ''})"
            if e.get("meta") else ""))
    out["data"] = report.data_report(cal, rows, log)

    # ----- keys whose process the kept data never sample are fixed (best-practice plan §2.5)
    coverage = data_rules.process_coverage(rows, cal.data.obs, fit_windows, cal.daytime_sw)
    fixed = data_rules.fix_by_coverage(keys, coverage, int(fs["min_process_records"]))
    out["process_coverage"], out["fixed_by_coverage"] = coverage, fixed
    log(f"kept records sampling each process: {coverage}")
    for k, why in fixed.items():
        log(f"fixed by coverage: {k}: {why}")

    # ----- the weights, then the screening: the first (central) gradient matrix at the prior centres
    out["weights_start"], _ = refresh_weights(runner, start, fit_rows, log)
    candidates = [i for i, p in enumerate(keys) if p.name not in fixed]
    fk = fit.FreeKeys(runner, candidates, start.copy())
    u0 = fk.u_prior
    r0 = fk.residuals([u0])[0]
    log(f"start: Phi = {fk.cost(u0, r0):.6g} over {len(r0)} rows; screening {len(candidates)} keys "
        f"({2 * len(candidates) * len(fit_windows)} trials)")
    J0, smooth0 = fit.gradient_central(fk, u0, r0)
    keep, out["screening"] = fit.triage(fk, J0, smooth0, fit_rows, float(fs["uninformed_sd_ratio"]),
                                        float(fs["max_correlation"]))
    for k, why in out["screening"]["fixed"].items():
        log(f"screening fixes {k}: {why}")
    if not keep:
        raise SystemExit("the screening left no key to fit")
    free = [candidates[k] for k in keep]
    log("fitting " + ", ".join(keys[i].name for i in free))
    out["cost_start"] = fk.cost(u0, r0)
    out["scores_cal_start"] = targets.target_scores(fit_rows, r0)

    # ----- Levenberg-Marquardt; the refresh of the chains, weights and sigma; Levenberg-Marquardt again
    lm_settings = dict(max_iter=int(fs["max_iter"]), min_cost_drop=float(fs["min_cost_drop"]),
                       recompute_every=int(fs["recompute_every"]), log=log)
    fk = fit.FreeKeys(runner, free, start.copy())
    res = fit.lm(fk, fk.u_of(start), label="fit", J=J0[:, keep], r=r0, **lm_settings)
    values = fk.values(res["u"])
    out["lm"] = {"first": {"cost": res["cost"], "iterations": res["iterations"], "history": res["history"]}}
    log("refresh: the state chains at the current values, then the weights and sigma there")
    out["values_refresh"] = {p.name: float(values[i]) for i, p in enumerate(keys)}
    runner.states = chains_at(values, all_windows)
    out["weights_refresh"], out["sigma_scale_refresh"] = refresh_weights(
        runner, values, fit_rows, log, scale_sigma_cap=float(fs["sigma_scale_max"]))
    fk = fit.FreeKeys(runner, free, values.copy())
    res = fit.lm(fk, fk.u_of(values), label="fit after the refresh", **lm_settings)
    values_map = fk.values(res["u"])
    out["lm"]["after_refresh"] = {"cost": res["cost"], "iterations": res["iterations"], "history": res["history"]}

    after_the_fit(cal, fk, values_map, fit_rows, out, log)
    out["map"] = {p.name: float(values_map[i]) for i, p in enumerate(keys)}
    out["default"] = {p.name: float(default[i]) for i, p in enumerate(keys)}
    out["start"] = {p.name: float(start[i]) for i, p in enumerate(keys)}
    if val_windows:
        validate(cal, runner, default, start, values_map, val_windows, states_start, chains_at, out, log)
    out["trials"] = {"run": runner.n_trials, "median_s": float(np.median(runner.seconds)) if runner.seconds else None,
                     "wall_s": time.time() - t_start}
    calibrated_files.write(cal, keys, values_map, work)
    report.write_report(out, work / "report.md")
    report.save_json(work / "fit.json", out)
    v = out.get("validation")
    log(f"done: validation {'passed' if v and v['pass'] else 'failed' if v else 'not run'}"
        + (f", cost {out['cost_val']['default']:.6g} (default) -> {out['cost_val']['map']:.6g} (calibrated)" if v else "")
        + f"; {runner.n_trials} trials in {(time.time() - t_start) / 60:.1f} min")
    workers.close()
    return 0


def after_the_fit(cal, fk, values_map, fit_rows, out, log):
    """The final central gradient matrix at the MAP and what follows from it: the covariance and its
    intervals, the Gauss-Newton step left, the linearity check, the key table (prior z, bounds), the
    declared alternatives, the model/tower ratios and kappa."""
    runner = fk.runner
    u_map = fk.u_of(values_map)
    r_map = fk.residuals([u_map])[0]
    log(f"the final gradient matrix (central differences, {2 * len(fk.free) * len(runner.windows)} trials)")
    J_map, smooth_map = fit.gradient_central(fk, u_map, r_map)
    out["chi2_per_row"] = targets.chi2_per_row(fit_rows, r_map)
    cov_u = fit.posterior(J_map, fk.prior_sd)
    out["fitted"] = [p.name for p in fk.keys]
    out["smoothness_map"] = {p.name: float(v) for p, v in zip(fk.keys, smooth_map)}
    out["covariance"] = fit.value_covariance(fk, u_map, cov_u)
    out["correlation"] = fit.correlation(cov_u)
    out["intervals"] = fit.intervals(fk, u_map, cov_u)
    out["sigma_ratio"] = {k: v["sigma_ratio"] for k, v in out["intervals"].items()}
    for k, v in out["intervals"].items():
        log(f"  {k:22s} {v['map']:10.4g}  68 % [{v['i68'][0]:.4g}, {v['i68'][1]:.4g}]  95 % "
            f"[{v['i95'][0]:.4g}, {v['i95'][1]:.4g}]  sd ratio {v['sigma_ratio']:.2f}  ({v['kind']})")
    #----- how far the fit stopped from its own optimum: the Gauss-Newton step left, in posterior sd
    rest = uncertainty.gauss_newton_step(fk, u_map, J_map, r_map) / np.sqrt(np.diag(cov_u))
    out["unconverged_sd"] = {p.name: float(v) for p, v in zip(fk.keys, rest)}
    log("the Gauss-Newton step left at the MAP [posterior sd] (0 at a converged fit): " + ", ".join(
        f"{n} {v:+.2f}" for n, v in out["unconverged_sd"].items()))
    out["linearity"] = uncertainty.linearity(fk, u_map, cov_u, J_map, int(cal.settings["uncertainty"]["linearity_dirs"]), log)
    out["local_only"] = any(d["local_only"] for d in out["linearity"])
    out["key_table"] = uncertainty.key_table(fk, u_map, J_map, r_map, fit_rows)
    for k, v in out["key_table"].items():
        log(f"  prior z of {k:22s} {v['z']:+.2f} ({v['kind']}, pushed by {v['pushed_by']})"
            + (f"; near its {v['near_bound']} bound, pushed by {v['bound_pushed_by']}" if v["near_bound"] else ""))
    out["cost"] = {"start": out["cost_start"], "map": fk.cost(u_map, r_map)}
    out["scores_cal"] = {"start": out["scores_cal_start"], "map": targets.target_scores(fit_rows, r_map)}
    out["alternatives"] = uncertainty.alternatives(cal, fk, u_map, runner.windows, cov_u, J_map, r_map,
                                                   out["sigma_scale_refresh"], log)
    out["ratios"] = report.ratio_tables(cal, runner, values_map, runner.windows)
    out["kappa"] = report.kappa_report(cal.data.obs, runner.keys, values_map, cal.daytime_sw)
    if out["kappa"]:
        kr = out["kappa"]
        log(f"kappa {kr['kappa']:.3f}: the tower's respiration {kr['reco_tower']:.2f} -> {kr['reco_implied']:.2f}, "
            f"GPP {kr['gpp_tower']:.2f} -> {kr['gpp_implied']:.2f} umol m-2 s-1 over its measured records")


def validate(cal, runner, default, start, values_map, val_windows, states_start, chains_at, out, log):
    """The default and the calibrated values scored on the validation windows the fit never saw, each
    from chains run with its own values: the default's from chains at the base configuration (the
    start's chains ran at the prior centres, which differ where a prior is centred off the base value)."""
    states_map = dict(runner.states)
    if np.allclose(default, start):
        runner.states = states_start
    else:
        log("running the state chains with the base configuration for the validation's default")
        runner.states = chains_at(default, val_windows)
    r_default = runner.residuals([default], val_windows)[0]
    runner.states = states_map
    r_map = runner.residuals([values_map], val_windows)[0]
    val_rows = [runner.rows[w.name] for w in val_windows]
    out["scores_val"] = {"default": targets.target_scores(val_rows, r_default), "map": targets.target_scores(val_rows, r_map)}
    out["cost_val"] = {"default": float(r_default @ r_default), "map": float(r_map @ r_map)}
    out["validation"] = report.validation_verdict(out["scores_val"], out["cost_val"])


# ------------------------------------------------------------------------------------------------
# variants
# ------------------------------------------------------------------------------------------------
def cmd_variants(args):
    """The structural variants side by side (best-practice plan §5.3, §7.1): each key's MAP in every
    finished fit given, and the spread of the MAPs in the largest posterior sd."""
    fits = {}
    for path in args.fits:
        rep = json.loads(Path(path).read_text())
        fits[rep.get("variant") or Path(path).parent.name] = rep
    keys = sorted(set().union(*(set(r["intervals"]) for r in fits.values())))
    table = {}
    for k in keys:
        row = {v: r["intervals"][k] for v, r in fits.items() if k in r["intervals"]}
        maps = [x["map"] for x in row.values()]
        #----- the spread over the largest posterior sd (the 68 % interval's half width), so keys of any range compare
        widest = max((x["i68"][1] - x["i68"][0]) / 2.0 for x in row.values())
        table[k] = {"map": {v: x["map"] for v, x in row.items()},
                    "spread_sd": float((max(maps) - min(maps)) / widest) if widest > 0 else None}
    report.save_json(Path(args.out), {"variants": list(fits), "keys": table,
                                      "cost_val": {v: r.get("cost_val") for v, r in fits.items()},
                                      "validation": {v: (r.get("validation") or {}).get("pass") for v, r in fits.items()}})
    print(f"{'key':22s} " + " ".join(f"{v:>16s}" for v in fits) + "   spread [sd]")
    for k, e in table.items():
        print(f"{k:22s} " + " ".join(f"{e['map'].get(v, float('nan')):16.4g}" for v in fits)
              + f"   {e['spread_sd'] if e['spread_sd'] is not None else float('nan'):.2f}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, runs=True):
        p.add_argument("--config", required=True, help="the calibration's settings (calibration.toml)")
        p.add_argument("--variant", default=None, help="a structural variant declared in [variants]")
        p.add_argument("--work", required=True, help="the working directory (trials, chains, results)")
        if runs:
            p.add_argument("--runner", default=trials.PYTHON_RUNNER,
                           help='how a trial runs: "python" (the Python API, python -m edmore.model; the '
                                'default) or the path of a edmore executable')
            p.add_argument("--workers", type=int, default=os.cpu_count(), help="trials at once on this machine")
            p.add_argument("--queue", action="store_true",
                           help="hand the trials to `worker`s on other nodes through a queue under --work")
            p.add_argument("--keep-netcdf", action="store_true", help="keep the trials' netCDF output")
    p = sub.add_parser("report")
    common(p, runs=False)
    p.set_defaults(func=cmd_report)
    p = sub.add_parser("check")
    common(p)
    p.add_argument("--key", default=None, help="the key whose gradient column must move the output (default: the first)")
    p.set_defaults(func=cmd_check)
    p = sub.add_parser("fit")
    common(p)
    p.set_defaults(func=cmd_fit)
    p = sub.add_parser("variants")
    p.add_argument("fits", nargs="+", help="the fit.json of each variant's fit")
    p.add_argument("--out", required=True, help="the comparison's JSON")
    p.set_defaults(func=cmd_variants)
    p = sub.add_parser("worker")
    p.add_argument("--queue", required=True, help="the queue directory (<work>/queue)")
    p.add_argument("--slots", type=int, default=os.cpu_count())
    p.set_defaults(func=lambda a: queue_worker(a.queue, a.slots) or 0)
    args = ap.parse_args(argv)
    #----- a trial runs in its own directory, so an executable's path must not be relative
    if getattr(args, "runner", trials.PYTHON_RUNNER) != trials.PYTHON_RUNNER:
        args.runner = str(Path(args.runner).resolve())
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
