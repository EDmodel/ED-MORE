# SPDX-License-Identifier: Apache-2.0
"""What the tool reports (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2, §7.2): the data report before
the fit, and after it the model/tower ratios, kappa's implications, the validation's verdict and
report.md, the results for a reader. fit.json holds everything.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import data_rules
import targets
import trials


def save_json(path, obj):
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, Path):
            return str(o)
        if isinstance(o, tuple):
            return list(o)
        raise TypeError(type(o).__name__)
    Path(path).write_text(json.dumps(obj, indent=1, default=conv))


# ------------------------------------------------------------------------------------------------
# before the fit
# ------------------------------------------------------------------------------------------------
def data_report(cal, rows: dict, log) -> dict:
    """The data report (best-practice plan §2, §6.2): every target's records over the whole record
    through each step, the u* diagnostics and the rule each target got, the windows and seasonal
    runs, what the fit's windows span of the record's conditions, the closure, sigma, and the tower
    reader's checks. `rows` holds the fit's windows (calibration.Calibration.rows_for)."""
    d = cal.data
    filters = targets.record_report(d.obs, d.forcing_observed, cal.target_settings, cal.daytime_sw,
                                    d.sun_elevation, cal.utc_offset_h)
    for t, r in filters.items():
        log(f"data {t:9s} " + " -> ".join(f"{label} {n}" for label, n in r["counts"]))
    for t, u in d.ustar.items():
        g = u["diagnostic"]
        if g is not None:
            b = g.get("bootstrap", {})
            log(f"data {t} u* diagnostic ({g['driver']} classes): {g['outcome']}"
                + (f" at {g['threshold']}" if g.get("threshold") is not None else "")
                + (f", bootstrap 5-50-95 % {b['threshold_p05_p50_p95']}, outcomes {b['outcome_share']}"
                   if b.get("threshold_p05_p50_p95") else "") + (f" ({g['note']})" if g.get("note") else ""))
        log(f"data {t} u* filter: {u['ustar_min']} ({u['reason']})")
    for w in cal.windows + cal.seasonal_runs:
        log(f"data window {w.name:14s} {w.role:8s} {w.start:%Y-%m-%d} + {w.days} d")
    if d.seasonal.get("years"):
        log("data water deficit by year (deepest inside the leaf-on months): " + ", ".join(
            f"{y['year']} {y['deficit_mm']:.0f} mm on {y['end']}" for y in d.seasonal["years"]))
    if d.seasonal.get("note"):
        log(f"data seasonal runs: {d.seasonal['note']}")
    #----- what the fit's kept records span of the record's daytime conditions
    kept = {k: [] for k in ("sw_in", "tair", "vpd", "deficit_mm")}
    for w in cal.fit_windows():
        wr = rows[w.name]
        o = d.obs.reindex(wr.index).iloc[sorted({int(i) for t in wr.targets for i in t.rows})]
        for k in ("sw_in", "tair", "vpd"):
            kept[k] += o[k].tolist()
        kept["deficit_mm"] += d.deficit.reindex(pd.date_range(w.start, w.end, freq="D")).tolist()
    day = d.obs["sw_in"] > cal.daytime_sw
    record = {k: d.obs.loc[day, k].to_numpy() for k in ("sw_in", "tair", "vpd")}
    record["deficit_mm"] = d.deficit.to_numpy()
    coverage = data_rules.range_coverage(record, kept)
    log("data the fit's records span of the record's daytime range: " + ", ".join(
        f"{k} {v['share']:.0%}" for k, v in coverage.items() if v))
    cl = d.closure
    if cl.get("f_median") is not None:
        log(f"data closure: daily (H + LE) / (Rnet - G) median {cl['daily_closure_median']:.2f} over {cl['valid_days']} "
            f"days, f median {cl['f_median']:.2f} (10-90 % {cl['f_p10_p90'][0]:.2f}-{cl['f_p10_p90'][1]:.2f})")
    if cl["attribution"].get("rises_by_vpd_class"):
        log("data closure attribution (rise from the calmest to the most turbulent third, by VPD class): "
            + "; ".join(f"{k.upper()} " + " ".join(f"{x:+.0%}" for x in v)
                        for k, v in cl["attribution"]["rises_by_vpd_class"].items()))
    log(f"data closure shares: H {cl['shares']['h']}, LE {cl['shares']['le']} ({cl['reason']})")
    for t, e in d.sigma.items():
        log(f"data {t} sigma = {e['sigma_abs']:.3g} + {e['sigma_rel']:.3g} |{e['flux']}| ({e['source']}), "
            f"evaluated at the smoothed {e['flux']}")
    return {"filters": filters, "ustar": d.ustar, "windows": d.windows, "seasonal": d.seasonal,
            "range_coverage": coverage, "closure": cl, "sigma": d.sigma, "tower_checks": d.reader}


# ------------------------------------------------------------------------------------------------
# after the fit
# ------------------------------------------------------------------------------------------------
def ratio_tables(cal, runner, values, windows) -> dict:
    """Per target: the model/tower ratio of the means overall, by local hour and by incoming-shortwave
    quartile, on the fit's rows at these values (the cached trials)."""
    dirs = runner.run([values], windows)[0]
    ok = runner.observation_keys(values)
    sw_all = cal.data.obs["sw_in"]
    acc = {}
    for w, td in zip(windows, dirs):
        df = trials.load_series(td).reindex(runner.rows[w.name].index)
        hr = targets.local_hours(df.index, cal.utc_offset_h)
        sw = sw_all.reindex(df.index).to_numpy()
        for t in runner.rows[w.name].targets:
            a = acc.setdefault(t.name, {"m": [], "o": [], "h": [], "sw": []})
            a["m"].append(targets.model_values(t, df))
            a["o"].append(targets.observed(t, ok))
            a["h"].append(hr[t.rows])
            a["sw"].append(sw[t.rows])
    out = {}
    for name, a in acc.items():
        m, o, h, sw = (np.concatenate(a[k]) for k in ("m", "o", "h", "sw"))
        by_hour = {int(x): float(m[h == x].mean() / o[h == x].mean()) for x in np.unique(h)
                   if (h == x).sum() >= 10 and abs(o[h == x].mean()) > 1e-9}
        q = np.nanquantile(sw, [0.25, 0.5, 0.75])
        cls = np.digitize(sw, q)
        by_light = {f"SW quartile {int(c) + 1}": float(m[cls == c].mean() / o[cls == c].mean())
                    for c in np.unique(cls) if (cls == c).sum() >= 10 and abs(o[cls == c].mean()) > 1e-9}
        out[name] = {"overall": float(m.mean() / o.mean()) if abs(o.mean()) > 1e-9 else None,
                     "by_local_hour": by_hour, "by_light": by_light}
    return out


def kappa_report(obs: pd.DataFrame, keys, values, daytime_sw: float) -> dict | None:
    """kappa, and the respiration and GPP it implies over the tower's measured records. The
    respiration is scaled at every hour; GPP gains the scaled respiration only by day, as in GPP's
    observation model (there is no GPP at night)."""
    names = [p.name for p in keys]
    if "kappa" not in names:
        return None
    k = float(values[names.index("kappa")])
    m = obs["gpp"].notna() & obs["reco"].notna()
    day = obs.loc[m, "sw_in"] > daytime_sw
    reco, gpp = obs.loc[m, "reco"].mean(), obs.loc[m, "gpp"].mean()
    to_kgc = 12.011e-9 * 365.25 * 86400.0
    gpp_implied = (obs.loc[m, "gpp"] + (1.0 / k - 1.0) * obs.loc[m, "reco"].where(day, 0.0)).mean()
    return {"kappa": k, "records": int(m.sum()), "reco_tower": float(reco), "reco_implied": float(reco / k),
            "gpp_tower": float(gpp), "gpp_implied": float(gpp_implied),
            "gpp_tower_kgc_yr": float(gpp * to_kgc), "gpp_implied_kgc_yr": float(gpp_implied * to_kgc),
            "note": "means over the records with GPP and its respiration measured, day and night; GPP gains "
                    "the scaled respiration by day only"}


def validation_verdict(scores: dict, cost: dict, max_worse: float = 0.10) -> dict:
    """The validation's verdict: the calibrated values beat the default on the validation windows'
    cost, and no target's normalized RMSE is more than max_worse worse than the default's."""
    worse = {k: scores["map"][k]["nrmse"] / scores["default"][k]["nrmse"] - 1.0 for k in scores["map"]}
    return {"pass": bool(cost["map"] < cost["default"] and all(v <= max_worse for v in worse.values())),
            "nrmse_change": worse,
            "rule": f"the calibrated cost below the default's, and no target's RMSE more than {max_worse:.0%} worse"}


def write_report(report: dict, path: Path):
    """report.md: the fit's results for a reader."""
    variant = f", variant {report['variant']}" if report.get("variant") else ""
    L = [f"# Fast calibration: {Path(report['config']).parent.name}{variant}", ""]
    v = report.get("validation")
    L += ["## Validation", ""]
    if v:
        L += [f"**{'Passed' if v['pass'] else 'Failed'}**: {v['rule']}.", "",
              "| target | default | calibrated | change |", "|---|---|---|---|"]
        sv = report["scores_val"]
        L += [f"| {t} | {sv['default'][t]['nrmse']:.2f} | {sv['map'][t]['nrmse']:.2f} | {v['nrmse_change'][t]:+.0%} |"
              for t in sv["map"]]
        L += ["", f"RMSE in units of each target's sigma. Cost: default {report['cost_val']['default']:.6g}, "
                  f"calibrated {report['cost_val']['map']:.6g}."]
    else:
        L += ["No validation windows."]
    L += ["", "The full record with the slow tier on is the next check: run the calibrated configs (budgets "
              "closed; dry-season GPP and LE no worse than the default's)."]
    L += ["", "## Keys", "", "| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for k, iv in report.get("intervals", {}).items():
        kt = report["key_table"][k]
        note = []
        if kt["kind"] == "trait" and abs(kt["z"]) > 2.0:
            note.append(f"beyond 2 prior sd, pushed by {kt['pushed_by']}: diagnose it (plan §5.1) or relabel it effective")
        if kt["near_bound"]:
            note.append(f"near its {kt['near_bound']} bound, pushed by {kt['bound_pushed_by']}")
        L.append(f"| {k} | {kt['kind']} | {kt['scope']} | {iv['map']:.4g} | {iv['i68'][0]:.4g}-{iv['i68'][1]:.4g} | "
                 f"{iv['prior_centre']:.4g} | {kt['z']:+.2f} | {iv['sigma_ratio']:.2f} | {'; '.join(note)} | "
                 f"{str(kt['prior_source'])[:80]} |")
    fixed = {**report.get("fixed_by_coverage", {}), **report.get("screening", {}).get("fixed", {})}
    if fixed:
        L += ["", "Fixed at their priors: " + "; ".join(f"{k} ({why})" for k, why in fixed.items())]
    L += ["", "## Targets", "", "| target | chi2/n | sigma scale | model/tower | n (cal) |", "|---|---|---|---|---|"]
    rt = report.get("ratios", {})
    for t, c in report.get("chi2_per_row", {}).items():
        sc = report["sigma_scale_refresh"].get(t, 1.0)
        n = report["scores_cal"]["map"][t]["n"]
        L.append(f"| {t} | {c:.2f} | {sc:.2f} | {rt.get(t, {}).get('overall') or float('nan'):.3f} | {n} |")
    for t, r in rt.items():
        L.append(f"\n{t}, model/tower by local hour: " + ", ".join(f"{h} h {x:.2f}" for h, x in r["by_local_hour"].items()))
        L.append(f"{t}, by light: " + ", ".join(f"{c} {x:.2f}" for c, x in r["by_light"].items()))
    kr = report.get("kappa")
    if kr:
        L += ["", "## kappa", "", f"kappa = {kr['kappa']:.3f}: the tower's respiration {kr['reco_tower']:.2f} -> "
              f"{kr['reco_implied']:.2f} umol m-2 s-1; GPP {kr['gpp_tower']:.2f} -> {kr['gpp_implied']:.2f} "
              f"({kr['gpp_tower_kgc_yr']:.2f} -> {kr['gpp_implied_kgc_yr']:.2f} kgC m-2 yr-1), {kr['note']}."]
    L += ["", "## Uncertainty", "", "Laplace, from the final gradient matrix"
          + (" -- LOCAL ONLY (the linearity check failed)" if report.get("local_only") else "") + "."]
    for name, a in report.get("alternatives", {}).items():
        if "keys" not in a:
            L.append(f"- alternative {name}: {a.get('note')}")
            continue
        L.append(f"- alternative {name} ({a['what']}): largest linear shift {a['worst_shift_sd']:.2f} sd"
                 + (" -- refitted: " + ", ".join(f"{k} {x['map']:.4g}->{x['refit']:.4g}" for k, x in a["keys"].items())
                    if "refit" in a else ""))
    path.write_text("\n".join(L) + "\n")
