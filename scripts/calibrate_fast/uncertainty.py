# SPDX-License-Identifier: Apache-2.0
"""The uncertainty after the fit (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §7.1), from the final
gradient matrix J at the MAP (fit.py's notation):

  linearity      along the posterior's leading directions, +-1 sd: the actual change in the cost
                 against the quadratic prediction. Outside 0.5-2x, the covariance is "local only".
  unconverged    the Gauss-Newton step left at the MAP, per key in posterior sd (0 at a converged fit).
  the key table  each key's prior z, whether it sits near a bound, and the target that pushes it there.
  alternatives   the declared alternatives (GPP's u* filter, the closure shares, the partitioning):
                 each one's shift of the MAP, first the linear estimate from J (the cached trials),
                 and the fit rerun from the MAP under it when a key moves more than refit_beyond_sd.
"""
from __future__ import annotations

import math

import numpy as np

import fit
import observation_models
from targets import target_slices

ALTERNATIVES = ("gpp_ustar", "closure", "partitioning")


def linearity(fk, u, cov_u, J, n_dir=3, log=print) -> list:
    """Along the leading principal directions of the posterior, +-1 sd: the actual change in the
    cost against the quadratic prediction d'(J'J + P)d. The mean of the two sides over the prediction
    is the actual curvature (outside [0.5, 2] the covariance is local only); half their
    difference is the slope left where the fit stopped short of its optimum (0 at a converged fit)."""
    vals, vecs = np.linalg.eigh(cov_u)
    order = np.argsort(vals)[::-1][:n_dir]
    r0 = fk.residuals([u])[0]
    base = fk.cost(u, r0)
    H = J.T @ J + np.diag(1.0 / fk.prior_sd ** 2)
    pts, pred = [], []
    for k in order:
        d = math.sqrt(vals[k]) * vecs[:, k]
        pts += [u + d, u - d]
        pred.append(float(d @ H @ d))
    R = fk.residuals(pts)
    out = []
    for n, k in enumerate(order):
        acts = [fk.cost(pts[2 * n + s], R[2 * n + s]) - base for s in (0, 1)]
        ratios = [a / pred[n] for a in acts]
        curv = 0.5 * (ratios[0] + ratios[1])
        out.append({"direction": int(k), "sd_u": float(math.sqrt(vals[k])), "predicted": pred[n],
                    "delta_cost_plus": acts[0], "delta_cost_minus": acts[1], "curvature_ratio": curv,
                    "slope_ratio": 0.5 * (ratios[0] - ratios[1]), "local_only": not (0.5 <= curv <= 2.0),
                    "loadings": {p.name: float(vecs[j, k]) for j, p in enumerate(fk.keys)}})
        log(f"linearity along direction {k}: curvature {curv:.2f} x the quadratic's, slope left "
            f"{out[-1]['slope_ratio']:+.2f}" + (" (local only)" if out[-1]["local_only"] else ""))
    return out


def gauss_newton_step(fk, u, J, r) -> np.ndarray:
    """One Gauss-Newton step of the cost from u: du = -(J'J + P)^-1 (J'r + P (u - u_prior))."""
    P = np.diag(1.0 / fk.prior_sd ** 2)
    return -np.linalg.solve(J.T @ J + P, J.T @ r + P @ (np.asarray(u) - fk.u_prior))


def key_table(fk, u, J, r, rows_list, near=0.05) -> dict:
    """Per key at u: its prior z ((u - u_prior) / sd_prior), whether it lies within `near` of a bound,
    and the target whose gradient pushes it furthest from its prior and toward that bound."""
    slices = target_slices(rows_list)
    out = {}
    for j, p in enumerate(fk.keys):
        grad = {}
        for name, s in slices:
            grad[name] = grad.get(name, 0.0) + float(J[s, j] @ r[s])

        def pushing(direction):
            """The target whose descent moves u_j furthest in `direction` (descending moves u by -grad)."""
            return max(grad, key=lambda k: -grad[k] * direction) if grad and direction else None
        z = float((u[j] - p.u0) / p.sd_u)
        pos = p.position(p.to_value(u[j]))
        bound = "lower" if pos < near else "upper" if pos > 1.0 - near else None
        out[p.name] = {"z": z, "kind": p.kind, "scope": p.scope, "prior_centre": p.centre,
                       "prior_source": p.prior.get("source", p.source), "pushed_by": pushing(np.sign(z)),
                       "near_bound": bound,
                       "bound_pushed_by": pushing(-1.0 if bound == "lower" else 1.0) if bound else None,
                       "range": [p.lo, p.hi]}
    return out


def alternative_inputs(cal, name):
    """The targets and observations of one declared alternative (best-practice plan §7.1), or a note
    saying why it does not apply:
      gpp_ustar     GPP filtered at the daytime plateau of its own u* diagnostic, where the fit used
                    the provider's threshold
      closure       Bowen (H and LE scaled together), where the fit used the attribution's shares
      partitioning  the provider's daytime partitioning (GPP_DT, RECO_DT), where the site TOML declares it"""
    data, ts = cal.data, cal.target_settings
    if name == "gpp_ustar":
        u = data.ustar.get("gpp", {})
        diag = u.get("diagnostic") or {}
        if u.get("rule") != "provider" or diag.get("outcome") != "plateau":
            return None, None, "the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau"
        alt = {k: dict(v) for k, v in ts.items()}
        alt["gpp"]["ustar_min"] = diag["threshold"]
        return alt, data.obs, f"GPP u* >= {diag['threshold']} (the fit used {ts['gpp']['ustar_min']})"
    if name == "closure":
        if data.closure["shares"]["h"] is None:
            return None, None, "the fit already used Bowen"
        obs = data.obs.copy()
        obs["h_c"], obs["le_c"] = observation_models.corrected(obs, None, None)
        return ts, obs, "Bowen (H and LE scaled together) against the fit's attribution shares"
    if name == "partitioning":
        if data.obs["gpp_dt"].notna().sum() == 0 or data.obs["reco_dt"].notna().sum() == 0:
            return None, None, "the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)"
        obs = data.obs.copy()
        obs["gpp"], obs["reco"] = obs["gpp_dt"], obs["reco_dt"]
        return ts, obs, "the provider's daytime partitioning (GPP_DT, RECO_DT)"
    raise SystemExit(f"[uncertainty].alternatives: unknown {name!r}; known: {ALTERNATIVES}")


def alternatives(cal, fk, u_map, windows, cov_u, J_map, r_map, sigma_scale, log=print) -> dict:
    """The declared alternatives. Each one's shift of the MAP is first the linear estimate from the
    final gradient matrix -- the Gauss-Newton step on the alternative's rows minus that on the fit's
    own, from the cached trials. When a key moves by more than [uncertainty].refit_beyond_sd posterior
    sd, the fit is rerun from the MAP under the alternative and both MAPs are reported."""
    us, fs = cal.settings["uncertainty"], cal.fit_settings
    sd = np.sqrt(np.diag(cov_u))
    runner = fk.runner
    out = {}
    for name in us["alternatives"]:
        target_settings, obs, what = alternative_inputs(cal, name)
        if target_settings is None:
            out[name] = {"note": what}
            log(f"alternative {name}: skipped ({what})")
            continue
        #----- the fit's windows rebuilt under the alternative, with the fit's own weights and sigma scales
        alt = cal.rows_for(windows, target_settings, obs)
        for w in windows:
            weight = {t.name: t.weight for t in runner.rows[w.name].targets}
            for t in alt[w.name].targets:
                t.sigma = t.sigma * sigma_scale.get(t.name, 1.0)
                t.weight = weight.get(t.name, 1.0)
        saved = runner.rows
        runner.rows = {**saved, **alt}
        try:
            r_alt = fk.residuals([u_map], windows)[0]
            J_alt, _ = fit.gradient_central(fk, u_map, r_alt)
            du = gauss_newton_step(fk, u_map, J_alt, r_alt) - gauss_newton_step(fk, u_map, J_map, r_map)
            entry = {"what": what, "keys": {}}
            for k, p in enumerate(fk.keys):
                entry["keys"][p.name] = {"map": float(p.to_value(u_map[k])), "linear": float(p.to_value(u_map[k] + du[k])),
                                         "shift_sd": float(du[k] / sd[k])}
            worst = max(abs(v["shift_sd"]) for v in entry["keys"].values())
            log(f"alternative {name} ({what}): linear shift " + ", ".join(
                f"{n} {v['map']:.4g}->{v['linear']:.4g} ({v['shift_sd']:+.2f} sd)" for n, v in entry["keys"].items()))
            if worst > float(us["refit_beyond_sd"]):
                log(f"alternative {name}: a key moves {worst:.2f} sd > {us['refit_beyond_sd']}: refitting from the MAP")
                res = fit.lm(fk, u_map, max_iter=int(us["refit_max_iter"]), min_cost_drop=float(fs["min_cost_drop"]),
                             log=log, label=f"refit under {name}", J=J_alt, r=r_alt,
                             recompute_every=int(fs["recompute_every"]))
                for k, p in enumerate(fk.keys):
                    v = entry["keys"][p.name]
                    v["refit"] = float(p.to_value(res["u"][k]))
                    v["refit_shift_sd"] = float((res["u"][k] - u_map[k]) / sd[k])
                entry["refit"] = {"cost": res["cost"], "iterations": res["iterations"]}
            entry["worst_shift_sd"] = worst
            out[name] = entry
        finally:
            runner.rows = saved
    return out

