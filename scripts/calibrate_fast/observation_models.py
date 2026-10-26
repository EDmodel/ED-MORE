# SPDX-License-Identifier: Apache-2.0
"""The targets' observation models (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §3.1, §3.3): how the
tower's number relates to the true flux, each known bias as its own term, and the random error.

Closure (H and LE). The tower misses the same fraction of the turbulent flux at every hour:

    f_d       = median over +-window_days of the daily (Rnet - G) / (H + LE)
    H_true,h  = H_obs,h  + s_H  (f_d - 1) (H_obs,h + LE_obs,h)
    LE_true,h = LE_obs,h + s_LE (f_d - 1) (H_obs,h + LE_obs,h)        s_H + s_LE = 1

Over a whole day the ground and canopy storage roughly net out, so no estimate of them is needed
(G = 0 where the site TOML declares none). Only days with at least min_measured of their records
measured count, and the provider's gap-filled values fill those days' sums. The shares come from
the attribution test: within VPD classes, which of H/Rnet and LE/Rnet rises from the calmest to the
most turbulent third of the records (by more than min_rise, median over the VPD classes). One
rising takes the whole share, both rising means Bowen (s_H = H / (H + LE), both scaled by f), and
neither means as measured. [closure].shares = "none" keeps both as measured.

Respiration (GPP). The provider made GPP from its own respiration, GPP_tower = R_tower - NEE_obs, and
that respiration carries one multiplicative error, R_tower = kappa R_true. So the comparison is

    r = [GPP_model - GPP_tower - (1/kappa - 1) R_tower] / sigma_NEE

with kappa an observation key (never written to a MEDS config) and sigma NEE's random error. This
is the same as fitting daytime NEE closed with the tower's own respiration. Without kappa in the
fit, kappa is 1 and the tower's GPP is taken as given.

Random error, sigma = sigma_abs + sigma_rel |x|, from (in order) the provider's per-record random
uncertainty, the paired-day estimate (Hollinger & Richardson 2005: the same half hour on two days in
a row with similar light, temperature, VPD and wind), or the target's defaults. It is evaluated at
a smoothed observation -- the mean of the measured values at the same time of day within
+-smooth_days -- not at the record's own: a sigma that follows each record's own value gives the
randomly low ones more weight and biases the fit low.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

#: the flux whose random error sets each target's sigma (GPP's is NEE's: the respiration is a term)
ERROR_FLUX = {"le": "le", "h": "h", "gpp": "nee"}


# ------------------------------------------------------------------------------------------------
# closure
# ------------------------------------------------------------------------------------------------
def closure_factor(values: pd.DataFrame, measured: pd.DataFrame, utc_offset_h: float, settings: dict):
    """f_d per record (NaN where the window holds too few valid days), and the daily table's
    summary. `values` holds Rnet, H, LE (and G) as the provider gives them, gap-filled included;
    `measured` the masks."""
    need = [c for c in ("Rnet", "H", "LE") if c not in values]
    if need:
        raise SystemExit(f"the closure model needs {need} in the site TOML's [fluxes] "
                         "(or [closure].shares = \"none\": H and LE as measured)")
    day = (values.index + pd.Timedelta(hours=utc_offset_h)).floor("D")
    g = values["G"] if "G" in values else pd.Series(0.0, index=values.index)
    present = values[["Rnet", "H", "LE"]].notna().all(axis=1) & g.notna()
    meas = measured["H"] & measured["LE"]
    by = pd.DataFrame({"avail": (values["Rnet"] - g).where(present), "turb": (values["H"] + values["LE"]).where(present),
                       "present": present.astype(float), "meas": meas.astype(float)}).groupby(day)
    daily = pd.DataFrame({"avail": by["avail"].sum(min_count=1), "turb": by["turb"].sum(min_count=1),
                          "present": by["present"].mean(), "meas": by["meas"].mean()})
    ok = (daily["present"] >= 0.95) & (daily["meas"] >= float(settings["min_measured"])) & (daily["turb"] > 0)
    daily["ratio"] = (daily["avail"] / daily["turb"]).where(ok)
    w = int(settings["window_days"])
    f = daily["ratio"].asfreq("D").rolling(2 * w + 1, center=True, min_periods=int(settings["min_days"])).median()
    per_record = pd.Series(f.reindex(day).to_numpy(), index=values.index)
    rep = {"valid_days": int(ok.sum()), "days": int(len(daily)),
           "daily_closure_median": float((1.0 / daily["ratio"]).median()) if ok.any() else None,
           "f_median": float(f.median()) if f.notna().any() else None,
           "f_p10_p90": [float(f.quantile(0.1)), float(f.quantile(0.9))] if f.notna().any() else None}
    return per_record, rep


def attribution(obs: pd.DataFrame, settings: dict, min_rnet: float, daytime_sw: float) -> dict:
    """The closure attribution test (module docstring): within equal-count VPD classes, the rise
    of H/Rnet and LE/Rnet from the calmest third of the records (by u*) to the most turbulent."""
    m = (obs["sw_in"] > daytime_sw) & obs[["h", "le", "rnet", "ustar", "vpd"]].notna().all(axis=1)
    m &= obs["rnet"] >= min_rnet
    o = obs[m]
    n_vpd = int(settings["vpd_classes"])
    if len(o) < 30 * n_vpd:
        return {"note": f"{len(o)} daytime records with H, LE, Rnet, u* and VPD: too few for the test",
                "rises": {"h": None, "le": None}}
    vclass = pd.qcut(o["vpd"].rank(method="first"), n_vpd, labels=False)
    rises = {"h": [], "le": []}
    for v in range(n_vpd):
        g = o[vclass == v]
        u = pd.qcut(g["ustar"].rank(method="first"), 3, labels=False)
        lo, hi = g[u == 0], g[u == 2]
        for k in ("h", "le"):
            rises[k].append(float((hi[k].sum() / hi["rnet"].sum()) / (lo[k].sum() / lo["rnet"].sum()) - 1.0))
    med = {k: float(np.median(v)) for k, v in rises.items()}
    return {"rises_by_vpd_class": rises, "rises": med, "records": int(len(o))}


def closure_shares(rule: str, test: dict, min_rise: float):
    """(s_H, s_LE, reason), s_H = None meaning Bowen (each record's own H / (H + LE))."""
    if rule == "none":
        return 0.0, 0.0, "as measured (the calibration's choice)"
    if rule == "bowen":
        return None, None, "Bowen: H and LE scaled together (the calibration's choice; the FLUXNET convention)"
    if rule != "attribution":
        raise ValueError(f'[closure].shares must be "attribution", "bowen" or "none", not {rule!r}')
    r = test["rises"]
    if r["h"] is None:
        return 0.0, 0.0, "too few records for the attribution test: as measured"
    h_up, le_up = r["h"] > min_rise, r["le"] > min_rise
    if h_up and le_up:
        return None, None, f"both rise with u* (H {r['h']:+.0%}, LE {r['le']:+.0%}): Bowen"
    if h_up:
        return 1.0, 0.0, f"H/Rnet rises with u* ({r['h']:+.0%}), LE/Rnet does not ({r['le']:+.0%}): the gap is H's"
    if le_up:
        return 0.0, 1.0, f"LE/Rnet rises with u* ({r['le']:+.0%}), H/Rnet does not ({r['h']:+.0%}): the gap is LE's"
    return 0.0, 0.0, f"neither rises with u* (H {r['h']:+.0%}, LE {r['le']:+.0%}): as measured"


def corrected(obs: pd.DataFrame, s_h, s_le) -> tuple[pd.Series, pd.Series]:
    """H and LE corrected for the closure gap with shares s_H and s_LE (both 0: as measured; both
    None: Bowen, H and LE scaled together); where f_d is unknown the record has no corrected value."""
    if s_h == 0.0 and s_le == 0.0:
        return obs["h"], obs["le"]
    f = obs["closure_f"]
    if s_h is None:
        return obs["h"] * f, obs["le"] * f
    gap = (f - 1.0) * (obs["h"] + obs["le"])
    return obs["h"] + s_h * gap, obs["le"] + s_le * gap


# ------------------------------------------------------------------------------------------------
# random error
# ------------------------------------------------------------------------------------------------
def smoothed(x: pd.Series, utc_offset_h: float, days: int) -> pd.Series:
    """For each record, the mean of the measured values at the same time of day within +-days."""
    local = x.index + pd.Timedelta(hours=utc_offset_h)
    tod = (local - local.floor("D")).total_seconds().astype(int)
    table = pd.DataFrame({"v": x.to_numpy(), "day": local.floor("D"), "tod": tod}).pivot_table(
        index="day", columns="tod", values="v", aggfunc="mean", dropna=False).asfreq("D")
    sm = table.rolling(2 * days + 1, center=True, min_periods=1).mean()
    stacked = sm.stack(future_stack=True)
    return pd.Series(stacked.reindex(list(zip(local.floor("D"), tod))).to_numpy(), index=x.index)


def fit_sigma(level: np.ndarray, spread, n_bins: int):
    """sigma = a + b |x| fitted to a spread estimate per equal-count bin of |x|: (a, b, bins), or
    None when the data show no spread (a sigma of 0 would make a row infinitely certain). a is at
    least the smallest bin's spread, b at least 0."""
    order = np.argsort(level)
    bins = np.array_split(order, n_bins)
    xs = np.array([level[b].mean() for b in bins])
    ys = np.array([spread(b) for b in bins])
    if not np.all(np.isfinite(ys)) or ys.min() <= 0.0:
        return None
    ws = np.array([len(b) for b in bins], float)
    A = np.vstack([np.ones_like(xs), xs]).T * np.sqrt(ws)[:, None]
    a, b = np.linalg.lstsq(A, ys * np.sqrt(ws), rcond=None)[0]
    return float(max(a, ys.min())), float(max(b, 0.0)), xs.tolist(), ys.tolist()


def paired_sigma(obs: pd.DataFrame, col: str, settings: dict, step: float):
    """The paired-day estimate of col's random error: (a, b, report) or None with too few pairs.
    A pair is the same record on two days in a row, both measured, with similar light (PAR, or
    incoming shortwave), air temperature, VPD and wind; epsilon = (x1 - x2) / sqrt(2) and sigma per
    bin of the pair's mean |x| is epsilon's sd."""
    lag = int(round(86400.0 / step))
    x = obs[col].to_numpy()
    x2 = np.roll(x, -lag)
    ok = np.isfinite(x) & np.isfinite(x2)
    ok[-lag:] = False
    light, tol = ("par", float(settings["paired_dpar"])) if obs["par"].notna().any() else ("sw_in", float(settings["paired_dsw"]))
    for c, t in ((light, tol), ("tair", float(settings["paired_dt"])), ("vpd", float(settings["paired_dvpd"])),
                 ("wind", float(settings["paired_dwind"]))):
        if obs[c].notna().any():
            v = obs[c].to_numpy()
            ok &= np.abs(np.roll(v, -lag) - v) < t
    n = int(ok.sum())
    if n < int(settings["min_pairs"]):
        return None
    eps = (x[ok] - x2[ok]) / math.sqrt(2.0)
    level = 0.5 * np.abs(x[ok] + x2[ok])
    fitted = fit_sigma(level, lambda idx: float(np.std(eps[idx])), int(settings["bins"]))
    if fitted is None:
        return None
    a, b, xs, ys = fitted
    return a, b, {"pairs": n, "bin_level": xs, "bin_sigma": ys}


def provider_sigma(x: pd.Series, unc: pd.Series, settings: dict):
    """sigma = a + b |x| fitted to the provider's per-record random uncertainty (median per bin)."""
    ok = (x.notna() & unc.notna()).to_numpy()
    if ok.sum() < int(settings["min_pairs"]):
        return None
    level, u = np.abs(x.to_numpy()[ok]), unc.to_numpy()[ok]
    fitted = fit_sigma(level, lambda idx: float(np.median(u[idx])), int(settings["bins"]))
    if fitted is None:
        return None
    a, b, xs, ys = fitted
    return a, b, {"records": int(ok.sum()), "bin_level": xs, "bin_sigma": ys}


def set_sigmas(targets: dict, obs: pd.DataFrame, settings: dict, utc_offset_h: float, step: float) -> dict:
    """Each turbulent target's sigma_abs and sigma_rel from the first source that has enough records
    (the provider's uncertainty, the paired days, else the target's defaults), and the column its
    sigma is evaluated at (the smoothed observation of its error flux, obs["<flux>_smooth"]). Changes
    `targets` in place; returns the report."""
    report = {}
    for t, c in targets.items():
        if not c["on"] or t not in ERROR_FLUX:
            continue
        flux = ERROR_FLUX[t]
        est, used = None, "default"
        unc = obs.get(f"{flux}_randunc")
        if unc is not None and unc.notna().any():
            est, used = provider_sigma(obs[flux], unc, settings), "provider"
        if est is None:
            est, used = paired_sigma(obs, flux, settings, step), "paired days"
        if est is None:
            used = "default"
        else:
            c["sigma_abs"], c["sigma_rel"] = est[0], est[1]
        obs[f"{flux}_smooth"] = smoothed(obs[flux], utc_offset_h, int(settings["smooth_days"]))
        c["sigma_at"] = f"{flux}_smooth"
        report[t] = {"source": used, "flux": flux, "sigma_abs": c["sigma_abs"], "sigma_rel": c.get("sigma_rel", 0.0),
                     **({"estimate": est[2]} if est else {})}
    return report


def kappa_sd(obs: pd.DataFrame) -> float | None:
    """The kappa prior's sd from the gap between the provider's two partitionings (night-time RECO
    against daytime RECO_DT), on records with both: None when the site TOML declares only one."""
    ok = obs["reco"].notna() & obs["reco_dt"].notna()
    if ok.sum() == 0 or obs.loc[ok, "reco"].sum() <= 0:
        return None
    return float(max(0.05, abs(obs.loc[ok, "reco_dt"].sum() / obs.loc[ok, "reco"].sum() - 1.0)))
