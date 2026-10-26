#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Fit the three vital-rate laws to the census tables; writes vital_rates.json for the ED-MORE driver.

    python fit_vital_rates.py          # reads data/ (prepare_census.py); a few minutes

Each law is fitted for each PFT, and each predicts a MEAN rate -- what a cohort carries:

  growth       [cm/yr]    g = [g_min + (g_max - g_min) / (1 + exp(-k (lnD - lnD0)))] exp(-b L)
  mortality    [1/yr]     m = gamma + alpha exp(-beta g)                       (Camac et al. 2018)
  recruitment  [1/m2/yr]  ln R = c0 + c1 LAI                                   (Poisson GLM)

D is dbh [cm]; L the overtopping LAI, the leaf area of taller trees within 20 m [m2/m2]; LAI the leaf
area index within 20 m of a quadrat's centre, all PFTs' together: shade, not seed source -- a PFT's own
leaf area nearby does not raise its recruitment in the census.

Growth rises with size from g_min to g_max, half-way at D0, and shade shrinks it by exp(-b L); it is
fitted to each tree's mean by Gamma quasi-likelihood, which needs growth above zero, so increments of
zero or less are set to a tenth of the smallest positive one. It is the growth of every tree alive at an
interval's start, as a cohort's is: a survivor's own, and a tree that died within the interval its last
measured growth (prepare_census.py) -- so from the 1990-1995 interval on, the first with a measurement
before it. Mortality takes g from the growth law
rather than the tree's measured growth, because a cohort's growth is the law's too; it is fitted by
maximum likelihood over each census interval, P(die) = 1 - exp(-m dt), so m is a rate per year.
Recruitment is fitted to each quadrat's rate weighted by its area x interval, which is the Poisson
model of the count. Skill is cross-validated on 1-ha blocks.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.linear_model import PoissonRegressor
from sklearn.metrics import r2_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PFTS = (1, 2, 3)
N_FOLD = 5


def growth_law(theta, dbh, lai_over):
    """[cm/yr] for theta = (g_min, g_max, D0, k, b)."""
    g_min, g_max, d0, k, b = theta
    return (g_min + (g_max - g_min) / (1.0 + np.exp(-k * (np.log(dbh) - np.log(d0))))) * np.exp(-b * lai_over)


def fit_growth(dbh, lai_over, g):
    """g_min, g_max, D0, k, b by Gamma quasi-likelihood (all kept positive), from three starting D0."""
    def deviance(log_theta):
        mu = growth_law(np.exp(log_theta), dbh, lai_over)
        return np.sum(g / mu + np.log(mu))
    fits = [minimize(deviance, np.log([0.05, 0.5, d0, 2.0, 0.01]), method="Nelder-Mead",
                     options={"maxiter": 8000, "xatol": 1e-7, "fatol": 1e-7}) for d0 in (7.0, 20.0, 50.0)]
    return [float(v) for v in np.exp(min(fits, key=lambda f: f.fun).x)]


def recruitment_terms(lai):
    return np.column_stack([lai])


def glm(model, X, y, weight=None):
    model.fit(X, y, sample_weight=weight)
    return [float(model.intercept_)] + [float(c) for c in model.coef_]


def log_link(coef, X):
    return np.exp(coef[0] + X @ np.asarray(coef[1:]))


def camac(theta, g):
    gamma, alpha, beta = theta
    return gamma + alpha * np.exp(-beta * g)


def fit_camac(g, dead, dt, start=(0.01, 0.05, 10.0)):
    """gamma, alpha, beta [1/yr, 1/yr, yr/cm] by maximum likelihood, all kept positive."""
    def negative_log_likelihood(log_theta):
        h = camac(np.exp(log_theta), g) * dt                  # expected deaths per tree in its interval
        return -(dead * np.log(-np.expm1(-h)) - (1 - dead) * h).sum()
    fit = minimize(negative_log_likelihood, np.log(start), method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 4000})
    return [float(v) for v in np.exp(fit.x)]


def fit_laws(g, s, r):
    """The three laws for each PFT, from the given rows."""
    laws = {"growth": [], "mortality": [], "recruitment": []}
    for p in PFTS:
        gp, sp, rp = g[g.pft == p], s[s.pft == p], r[r.pft == p]
        a = fit_growth(gp.dbh.to_numpy(), gp.lai_over.to_numpy(), gp.g_fit.to_numpy())
        laws["growth"].append(a)
        laws["mortality"].append(fit_camac(growth_law(a, sp.dbh.to_numpy(), sp.lai_over.to_numpy()),
                                           sp.dead.to_numpy(), sp.dt.to_numpy()))
        laws["recruitment"].append(glm(PoissonRegressor(alpha=0.0, solver="newton-cholesky", max_iter=1000),
                                       recruitment_terms(rp.lai), rp.rate, rp.exposure))
    return laws


def predict(laws, g, s, r):
    """Each row's growth [cm/yr], chance of dying within its interval, and recruitment [1/m2/yr]."""
    gp, ps, rp = np.zeros(len(g)), np.zeros(len(s)), np.zeros(len(r))
    for k, p in enumerate(PFTS):
        a, theta, c = laws["growth"][k], laws["mortality"][k], laws["recruitment"][k]
        m = (g.pft == p).to_numpy()
        gp[m] = growth_law(a, g.dbh[m], g.lai_over[m])
        m = (s.pft == p).to_numpy()
        ps[m] = -np.expm1(-camac(theta, growth_law(a, s.dbh[m], s.lai_over[m])) * s.dt[m])
        m = (r.pft == p).to_numpy()
        rp[m] = log_link(c, recruitment_terms(r.lai[m]))
    return gp, ps, rp


def class_mean_r2(obs, pred, df):
    """R2 of the means by PFT, size class and overtopping-LAI class -- what a cohort carries."""
    keys = [df.pft.values, np.digitize(df.dbh.values, [2, 5, 10, 20, 40, 80]),
            np.digitize(df.lai_over.values, [0.5, 1, 2, 3, 4, 5, 6])]
    t = pd.DataFrame({"o": obs, "p": pred}).groupby(keys).agg(o=("o", "mean"), p=("p", "mean"),
                                                              n=("o", "size"))
    t = t[t.n >= 50]
    return r2_score(t.o, t.p, sample_weight=t.n)


g = pd.read_csv(os.path.join(DATA, "growth.csv.gz"))
g = g[g.interval >= 1990].reset_index(drop=True)              # the dying trees' last growth exists from here
print(f"growth rows: {len(g)}, {g.died.sum()} of them trees that died within their interval "
      f"({g.died.mean() * 100:.1f} %)")
s = pd.read_csv(os.path.join(DATA, "survival.csv.gz"))
r = pd.read_csv(os.path.join(DATA, "recruits.csv.gz")).dropna(subset=["dt"])
epsilon = 0.1 * g.g[g.g > 0].min()
g["g_fit"] = np.maximum(g.g, epsilon)

# skill: each 1-ha block predicted by laws fitted without it
blocks = np.random.default_rng(0).permutation(np.unique(np.r_[g.block, s.block, r.block]))
fold = {b: i % N_FOLD for i, b in enumerate(blocks)}
g_cv, p_cv, r_cv = np.zeros(len(g)), np.zeros(len(s)), np.zeros(len(r))
for k in range(N_FOLD):
    held = [t.block.map(fold).to_numpy() == k for t in (g, s, r)]
    laws = fit_laws(g[~held[0]], s[~held[1]], r[~held[2]])
    pred = predict(laws, g[held[0]], s[held[1]], r[held[2]])
    g_cv[held[0]], p_cv[held[1]], r_cv[held[2]] = pred
w = r.exposure.to_numpy()
print(f"growth       R2 {r2_score(g.g, g_cv):.2f} for trees, {class_mean_r2(g.g, g_cv, g):.3f} for class "
      f"means; mean {g.g.mean():.4f} measured, {g.g_fit.mean():.4f} with <= 0 set to {epsilon:.5f}, "
      f"{g_cv.mean():.4f} predicted [cm/yr]")
print(f"mortality    AUC {roc_auc_score(s.dead, p_cv):.2f} for trees, R2 {class_mean_r2(s.dead, p_cv, s):.3f} "
      f"for class means; {s.dead.mean() * 100:.1f} % die, {p_cv.mean() * 100:.1f} % predicted")
print(f"recruitment  R2 {r2_score(r.rate, r_cv, sample_weight=w):.2f} for quadrats; plot total "
      f"{np.average(r.rate, weights=w) * 3e4:.0f} observed, {np.average(r_cv, weights=w) * 3e4:.0f} "
      f"predicted [stems/ha/yr]")

laws = fit_laws(g, s, r)
out = {
    "about": "Census-trained vital rates for examples/example03_demography (fit_vital_rates.py); "
             "one value per PFT. Recruitment's leaf area indices beyond 'range' take the range's edge.",
    "growth": {"law": "[g_min + (g_max - g_min) / (1 + exp(-k (lnD - lnD0)))] exp(-b L) [cm/yr], L the overtopping LAI",
               **{name: [t[i] for t in laws["growth"]] for i, name in enumerate(("g_min", "g_max", "D0", "k", "b"))},
               "epsilon": float(epsilon)},
    "mortality": {"law": "gamma + alpha exp(-beta growth) [1/yr]",
                  "gamma": [t[0] for t in laws["mortality"]], "alpha": [t[1] for t in laws["mortality"]],
                  "beta": [t[2] for t in laws["mortality"]]},
    "recruitment": {"law": "exp(c0 + c1 LAI) [1/m2/yr]", "coefficients": laws["recruitment"]},
    "range": {"lai": [0.0, float(r.lai.max())]},
}
with open(os.path.join(HERE, "vital_rates.json"), "w") as fh:
    json.dump(out, fh, indent=1)
for k, p in enumerate(PFTS):
    print(f"PFT {p}: growth {np.round(laws['growth'][k], 4)}; gamma, alpha, beta "
          f"{np.round(laws['mortality'][k], 4)}; recruitment {np.round(laws['recruitment'][k], 4)}")
print("wrote vital_rates.json")
