# SPDX-License-Identifier: Apache-2.0
"""The estimation (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §6.3): one joint Levenberg-Marquardt fit
of every free key on the stacked residuals of every target in every window, with a Gaussian prior
per key in the transformed space, the triage of the first gradient matrix, and the Laplace
covariance at the MAP.

In this module `u` is the keys' transformed values (parameters.py: unbounded, so the fit needs no
bounds), `J` the gradient matrix d(residuals)/du, and the objective (the "cost") is

    Phi(u) = || r(u) ||^2 + || (u - u_prior) / sd_prior ||^2

The residuals come from the trial runner (trials.TrialRunner), which runs one trial per (set of
values, window), every trial of a batch at once. The fit's iterations use one-sided differences,
with J carried between full recomputations by Broyden's rank-one update; the screening at the start
and the uncertainty at the end use central differences, which also measure each key's smoothness.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from parameters import GRADIENT_STEP_U, interval
from targets import target_slices


@dataclass
class FreeKeys:
    """The keys the fit moves (the rest stay at their values), their prior, and u -> values."""
    runner: object                   # trials.TrialRunner
    free: list                       # indices into runner.keys
    fixed_values: np.ndarray         # every key's value; the free ones are overwritten

    @property
    def keys(self):
        return [self.runner.keys[i] for i in self.free]

    @property
    def u_prior(self) -> np.ndarray:
        return np.array([p.u0 for p in self.keys])

    @property
    def prior_sd(self) -> np.ndarray:
        """The prior's sd in u, per free key."""
        return np.array([p.sd_u for p in self.keys])

    def values(self, u) -> np.ndarray:
        """Every key's value at u."""
        out = self.fixed_values.copy()
        for k, i in enumerate(self.free):
            out[i] = float(self.runner.keys[i].to_value(u[k]))
        return out

    def u_of(self, values) -> np.ndarray:
        return np.array([float(self.runner.keys[i].to_u(values[i])) for i in self.free])

    def prior_residuals(self, u) -> np.ndarray:
        return (np.asarray(u) - self.u_prior) / self.prior_sd

    def cost(self, u, r) -> float:
        return float(r @ r + self.prior_residuals(u) @ self.prior_residuals(u))

    def residuals(self, us: list, windows=None) -> list:
        return self.runner.residuals([self.values(u) for u in us], windows)


def gradient_central(fk: FreeKeys, u, r0, h=GRADIENT_STEP_U):
    """J by central differences in u, every column's two trials in one batch. Returns J and each
    column's smoothness: the agreement of its + and - one-sided slopes, 1 for a smooth key."""
    k = len(u)
    pts = []
    for j in range(k):
        e = np.zeros(k)
        e[j] = h
        pts += [u + e, u - e]
    R = fk.residuals(pts)
    J = np.zeros((len(r0), k))
    smooth = np.full(k, np.nan)
    for j in range(k):
        rp, rm = R[2 * j], R[2 * j + 1]
        J[:, j] = (rp - rm) / (2.0 * h)
        dp, dm = rp - r0, r0 - rm
        den = float(dm @ dm)
        smooth[j] = float(dp @ dm) / den if den > 0 else (1.0 if float(dp @ dp) == 0 else np.inf)
    return J, smooth


def gradient_forward(fk: FreeKeys, u, r0, h=GRADIENT_STEP_U) -> np.ndarray:
    """J by forward differences in u, every column's trial in one batch."""
    k = len(u)
    pts = []
    for j in range(k):
        e = np.zeros(k)
        e[j] = h
        pts.append(u + e)
    R = fk.residuals(pts)
    J = np.zeros((len(r0), k))
    for j in range(k):
        J[:, j] = (R[j] - r0) / h
    return J


def lm(fk: FreeKeys, u_start, max_iter=10, min_cost_drop=1e-3, lam0=1e-2, max_step=2.0, log=print, label="",
       J=None, r=None, recompute_every=3):
    """Levenberg-Marquardt from u_start. J is the given one (or forward differences at the start),
    carried by Broyden's rank-one update after each accepted step, and recomputed in full every
    `recompute_every` accepted steps and whenever an updated J finds no descent. Each iteration tries
    three damping values at once. It stops when the cost falls by less than min_cost_drop (a
    fraction), or after max_iter iterations."""
    u = np.array(u_start, dtype=float)
    if r is None:
        r = fk.residuals([u])[0]
    cost = fk.cost(u, r)
    if J is None:
        J = gradient_forward(fk, u, r)
    fresh, lam, since = True, lam0, 0
    hist = [{"iter": 0, "cost": cost, "lambda": lam, "u": u.tolist()}]
    log(f"{label} start: Phi = {cost:.6g}")
    it = 0
    while it < max_iter:
        it += 1
        Jf = np.vstack([J, np.diag(1.0 / fk.prior_sd)])
        rf = np.concatenate([r, fk.prior_residuals(u)])
        A, g = Jf.T @ Jf, Jf.T @ rf
        D = np.diag(np.maximum(np.diag(A), 1e-12))
        accepted = False
        for _ in range(4):
            lams = [lam / 10.0, lam, lam * 10.0]
            cands = []
            for lm_ in lams:
                d = np.linalg.solve(A + lm_ * D, -g)
                s = np.max(np.abs(d))
                cands.append(u + (d * max_step / s if s > max_step else d))
            R = fk.residuals(cands)
            costs = [fk.cost(c, rc) for c, rc in zip(cands, R)]
            b = int(np.argmin(costs))
            if costs[b] < cost:
                accepted = True
                break
            lam *= 100.0
        if not accepted:
            hist.append({"iter": it, "cost": cost, "lambda": lam, "u": u.tolist(), "accepted": False})
            if fresh:
                log(f"{label} iter {it}: no descent at any damping; stopping")
                break
            log(f"{label} iter {it}: no descent with the updated gradient matrix; recomputing it")
            J, fresh, since, lam = gradient_forward(fk, u, r), True, 0, lam0
            continue
        rel = (cost - costs[b]) / cost
        du, dr = cands[b] - u, R[b] - r
        J = J + np.outer(dr - J @ du, du) / float(du @ du)            # Broyden's rank-one update
        u, r, cost, lam = cands[b], R[b], costs[b], lams[b]
        fresh, since = False, since + 1
        hist.append({"iter": it, "cost": cost, "lambda": lam, "u": u.tolist(), "accepted": True})
        log(f"{label} iter {it}: Phi = {cost:.6g} (drop {100 * rel:.2f} %, lambda {lam:.3g})")
        if rel < min_cost_drop:
            break
        if since >= recompute_every and it < max_iter:
            J, fresh, since = gradient_forward(fk, u, r), True, 0
    return {"u": u, "r": r, "cost": cost, "J": J, "history": hist, "iterations": it}


def posterior(J, prior_sd) -> np.ndarray:
    """The Laplace covariance in u: (J'J + Sigma_prior^-1)^-1."""
    return np.linalg.inv(J.T @ J + np.diag(1.0 / prior_sd ** 2))


def value_covariance(fk: FreeKeys, u, cov_u):
    """The covariance in the keys' own units."""
    d = np.array([p.dvalue_du(uk) for p, uk in zip(fk.keys, u)])
    return cov_u * np.outer(d, d)


def correlation(cov):
    s = np.sqrt(np.diag(cov))
    return cov / np.outer(s, s)


def intervals(fk: FreeKeys, u, cov_u) -> dict:
    """Per key: the MAP, its sd in u, the 68 % and 95 % intervals in its own units (asymmetric,
    inside the range), and the posterior-to-prior sd ratio."""
    out = {}
    for k, p in enumerate(fk.keys):
        sd = math.sqrt(cov_u[k, k])
        out[p.name] = {"map": float(p.to_value(u[k])), "sd_u": sd, "i68": interval(p, u[k], sd, 1.0),
                       "i95": interval(p, u[k], sd, 1.96), "sigma_ratio": sd / p.sd_u,
                       "prior_centre": p.centre, "range": [p.lo, p.hi], "kind": p.kind}
    return out


def triage(fk: FreeKeys, J, smooth, rows_list, uninformed_sd_ratio=0.9, max_correlation=0.95):
    """The screening of the first gradient matrix (best-practice plan §4.1, §6.3): per key, its
    sensitivity per target, its posterior-to-prior sd ratio and its smoothness. A key is fixed at its
    prior's centre when its column is zero (dead), its response is rough (the two one-sided slopes
    disagree by more than 3x), the data barely inform it (sd ratio >= uninformed_sd_ratio), or it is
    the less informed of a pair correlated beyond max_correlation. Returns (indices kept, report)."""
    keys = fk.keys
    cov_u = posterior(J, fk.prior_sd)
    ratio = np.sqrt(np.diag(cov_u)) / fk.prior_sd
    corr = correlation(cov_u)
    sens = {}
    for j, p in enumerate(keys):
        per = {}
        for name, s in target_slices(rows_list):
            col = J[s, j] * fk.prior_sd[j]
            acc = per.setdefault(name, [0.0, 0])
            acc[0] += float(col @ col)
            acc[1] += len(col)
        sens[p.name] = {k: math.sqrt(v[0] / v[1]) for k, v in per.items() if v[1]}
    why = {}
    for j, p in enumerate(keys):
        if not np.any(J[:, j]):
            why[j] = "dead: its gradient column is zero (a key the windows never use)"
        elif np.isfinite(smooth[j]) and not 1 / 3 <= smooth[j] <= 3:
            why[j] = f"rough: the two one-sided slopes disagree (smoothness {smooth[j]:.2f}); make the model continuous"
        elif ratio[j] >= uninformed_sd_ratio:
            why[j] = f"uninformed: posterior/prior sd ratio {ratio[j]:.2f} >= {uninformed_sd_ratio}"
    keep = [j for j in range(len(keys)) if j not in why]
    #----- correlated pairs among the rest: fix the less informed one, until none is left
    while len(keep) > 1:
        c = correlation(posterior(J[:, keep], fk.prior_sd[keep]))
        np.fill_diagonal(c, 0.0)
        a, b = np.unravel_index(np.argmax(np.abs(c)), c.shape)
        if abs(c[a, b]) < max_correlation:
            break
        ja, jb = keep[a], keep[b]
        drop, other = (ja, jb) if ratio[ja] >= ratio[jb] else (jb, ja)
        why[drop] = f"collinear with {keys[other].name} (correlation {c[a, b]:+.3f}): the less informed of the two"
        keep.remove(drop)
    report = {"keys": [p.name for p in keys],
              "sigma_ratio": {p.name: float(ratio[j]) for j, p in enumerate(keys)},
              "smoothness": {p.name: float(smooth[j]) for j, p in enumerate(keys)},
              "sensitivity": sens,
              "collinear": [(keys[a].name, keys[b].name, float(corr[a, b]))
                            for a in range(len(keys)) for b in range(a + 1, len(keys)) if abs(corr[a, b]) > 0.9],
              "singular_values": (np.linalg.svd(J * fk.prior_sd, compute_uv=False) if J.size else np.zeros(0)).tolist(),
              "fitted": [keys[j].name for j in keep],
              "fixed": {keys[j].name: why[j] for j in sorted(why)}}
    return keep, report
