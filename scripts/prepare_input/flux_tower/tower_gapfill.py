# SPDX-License-Identifier: Apache-2.0
"""tower_gapfill.py -- explicit, flagged gap filling for a tower forcing file
(scripts/prepare_input/flux_tower/README.md, "Gap filling"). ED-MORE never gap-fills, so every value the file carries is
either observed or filled here, and the file says which with a <Var>_qc flag:

  0  observed
  1  short gap (at most gapfill.short_gap_max records): linear in time; wind in the energy form;
     shortwave through the clearness index
  3  the model's own longwave synthesis regressed onto the tower, or the mean diurnal variation
     (2 is unused: a fill from another source, such as ERA5-Land, is the user's to make first)
  4  filled by the data provider (FLUXNET *_QC > 0)
  5  relative humidity recovered from the provider's VPD through its own saturation curve

Every function works on the tower's own records: numpy arrays aligned with the tower's UTC stamps,
each record the mean over its interval.
"""
import numpy as np
import pandas as pd

import tower_inputs as ti

conv = ti.conv
QC_OBSERVED, QC_SHORT, QC_SYNTH_OR_MDV, QC_PROVIDER, QC_FROM_VPD = 0, 1, 3, 4, 5
MIN_FIT_POINTS = 48            # a regression group smaller than this falls back to the pooled fit


def missing_runs(missing):
    """(first, last) index pairs of each run of True in `missing`."""
    m = np.concatenate([[False], np.asarray(missing, dtype=bool), [False]])
    edges = np.flatnonzero(np.diff(m.astype(int)))
    return list(zip(edges[::2], edges[1::2] - 1))


def fill_short(y, qc, max_len, kind="linear", mean_cosz=None):
    """Fill runs of at most max_len missing records bounded on both sides by values. `kind` is
    'linear', 'energy' (wind: linear in u^2) or 'shortwave' (linear in the clearness index
    SW / (S0 <cos z>), zero at night; needs mean_cosz)."""
    y = y.copy()
    qc = qc.copy()
    n = len(y)
    for a, b in missing_runs(~np.isfinite(y)):
        if b - a + 1 > max_len or a == 0 or b == n - 1:
            continue
        i0, i1 = a - 1, b + 1
        w = (np.arange(a, b + 1) - i0) / (i1 - i0)
        if kind == "energy":
            y[a:b + 1] = np.sqrt((1 - w) * y[i0] ** 2 + w * y[i1] ** 2)
        elif kind == "shortwave":
            toa = conv.SOLAR_CONSTANT * mean_cosz
            day = mean_cosz > conv.COSZ_BAR_MIN
            k0 = y[i0] / toa[i0] if day[i0] else np.nan
            k1 = y[i1] / toa[i1] if day[i1] else np.nan
            k0 = k1 if not np.isfinite(k0) else k0
            k1 = k0 if not np.isfinite(k1) else k1
            idx = np.arange(a, b + 1)
            if not np.isfinite(k0) and day[idx].any():
                continue                                 # daylight between two nights: leave it
            kt = np.clip((1 - w) * np.nan_to_num(k0) + w * np.nan_to_num(k1), 0.0, 1.0)
            y[a:b + 1] = np.where(day[idx], kt * toa[idx], 0.0)
        else:
            y[a:b + 1] = (1 - w) * y[i0] + w * y[i1]
        qc[a:b + 1] = QC_SHORT
    return y, qc


def fill_mean_diurnal(y, qc, records_per_day, windows_days=(7, 14, 30)):
    """Mean diurnal variation (Falge et al. 2001): each missing record takes the mean of the
    OBSERVED records at the same time of day within +-7 days, widening to 14 and 30 when there are
    none."""
    y = y.copy()
    qc = qc.copy()
    observed = np.isfinite(y) & (qc == QC_OBSERVED)
    n = len(y)
    for i in np.flatnonzero(~np.isfinite(y)):
        for days in windows_days:
            j = i + records_per_day * np.arange(-days, days + 1)
            j = j[(j >= 0) & (j < n)]
            j = j[observed[j]]
            if j.size:
                y[i] = float(np.mean(y[j]))
                qc[i] = QC_SYNTH_OR_MDV
                break
    return y, qc


def regression_groups(stamps_utc, mean_cosz, by_day_night):
    """Calendar month, and day or night when asked: the groups a regression is fitted within."""
    month = pd.DatetimeIndex(np.asarray(stamps_utc, dtype="datetime64[s]")).month.to_numpy()
    if not by_day_night:
        return month
    return month * 2 + (mean_cosz > conv.COSZ_BAR_MIN).astype(int)


def _design(x):
    """[1, x1, x2, ...] from one predictor (1-D) or several (2-D, one column each)."""
    x = np.asarray(x, dtype=float)
    x = x[:, None] if x.ndim == 1 else x
    return np.column_stack([np.ones(len(x)), x])


def fit_linear(y, x, fit, groups, fallback=None):
    """Least-squares coefficients per group for y = c0 + c1 x1 + ..., fitted where `fit`; a group
    with fewer than MIN_FIT_POINTS takes the pooled fit, and with fewer than that pooled, `fallback`
    (default: the identity on the first predictor). Returns the coefficients, the pooled fit, and
    whether the fallback was used."""
    design = _design(x)

    def ols(m):
        if m.sum() < MIN_FIT_POINTS:
            return None
        c, *_ = np.linalg.lstsq(design[m], y[m], rcond=None)
        return c
    if fallback is None:
        fallback = np.zeros(design.shape[1])
        fallback[1] = 1.0
    pooled = ols(fit)
    fell_back = pooled is None
    pooled = np.asarray(fallback, dtype=float) if fell_back else pooled
    coefficients = {}
    for g in np.unique(groups):
        c = ols(fit & (groups == g))
        coefficients[int(g)] = pooled if c is None else c
    return coefficients, pooled, fell_back


def apply_linear(x, groups, coefficients):
    design = _design(x)
    c = np.array([coefficients[int(g)] for g in groups])
    return np.sum(design * c, axis=1)


def fill_by_regression(y, qc, x, groups, code, lower=None, upper=None, fallback=None):
    """Fill the missing values of y from predictors x (the parts of the model's synthesis) by a
    linear regression fitted on the observed overlap within each group; with too little overlap to
    fit, the coefficients `fallback`. Returns the filled series, its qc, and the
    fit (for the report)."""
    y = y.copy()
    qc = qc.copy()
    usable = np.all(np.isfinite(_design(x)), axis=1)
    fit = np.isfinite(y) & usable & (qc == QC_OBSERVED)
    coefficients, pooled, fell_back = fit_linear(y, x, fit, groups, fallback)
    predicted = np.where(usable, apply_linear(np.nan_to_num(x), groups, coefficients), np.nan)
    if lower is not None or upper is not None:
        predicted = np.clip(predicted, lower, upper)
    fill = ~np.isfinite(y) & np.isfinite(predicted)
    y[fill] = predicted[fill]
    qc[fill] = code
    residual = y[fit] - predicted[fit]
    report = dict(fit_points=int(fit.sum()), filled=int(fill.sum()), fell_back_unfitted=bool(fell_back),
                  pooled_coefficients=[round(float(c), 6) for c in pooled],
                  overlap_bias=float(np.mean(residual)) if fit.any() else None,
                  overlap_rmse=float(np.sqrt(np.mean(residual ** 2))) if fit.any() else None)
    return y, qc, report


# ---------------------------------------------------------------------------------------------
# The longwave predictors: the model's synthesis.
# ---------------------------------------------------------------------------------------------
def clearness_held_through_night(sw, mean_cosz):
    """The clearness index the model's longwave synthesis sees: by day the interval's SW over its
    mean top-of-atmosphere flux, and after dark the last daytime value (dusk's cloudiness), seeded
    clear (1) before the first sunrise."""
    kt = np.where(mean_cosz > conv.COSZ_BAR_MIN,
                  np.clip(np.maximum(sw, 0.0) / (conv.SOLAR_CONSTANT * np.maximum(mean_cosz, 1e-30)), 0.0, 1.0),
                  np.nan)
    kt = pd.Series(kt).ffill().fillna(1.0).to_numpy()
    return kt


def synthesized_longwave(tair_k, rh, psurf_ground_pa, sw, mean_cosz, cloud_a=conv.LW_CLOUD_A):
    """ED-MORE's lwdown_source = "synthesize" from the tower's own temperature, humidity, pressure and
    shortwave (docs/science/forcing.md eq. (17))."""
    q = conv.rh_to_specific_humidity(rh, tair_k, psurf_ground_pa)
    kt = clearness_held_through_night(sw, mean_cosz)
    return conv.synthesize_lwdown(tair_k, q, psurf_ground_pa, kt, cloud_a)


def synthesis_predictors(tair_k, rh, psurf_ground_pa, sw, mean_cosz):
    """The two parts of the model's synthesis, eps_clear sigma T^4 and eps_clear sigma T^4 (1 - kt),
    as separate predictors. Regressing on both fits the cloud term's coefficient at the site
    instead of taking the model's 0.22, and on held-out records the two parts predict better than
    one regression on the whole synthesis."""
    q = conv.rh_to_specific_humidity(rh, tair_k, psurf_ground_pa)
    kt = clearness_held_through_night(sw, mean_cosz)
    clear = conv.synthesize_lwdown(tair_k, q, psurf_ground_pa, np.ones_like(q), 0.0)
    return np.column_stack([clear, clear * (1.0 - kt)])


# With too few observations to fit, the synthesis predictors fall back to the model's own
# synthesis: clear + 0.22 clear (1 - kt), which is what lwdown_source = "synthesize" would give.
SYNTHESIS_FALLBACK = (0.0, 1.0, conv.LW_CLOUD_A)


# ---------------------------------------------------------------------------------------------
# Re-centring the states: ED-MORE reads a state as a value at its stamp (docs/science/forcing.md §3).
# ---------------------------------------------------------------------------------------------
def recentre(y, qc, stamp, energy=False):
    """Interval means -> values AT the stamps, which is how ED-MORE reads a state. On a begin-stamped
    file the stamp t_k starts interval k, so x(t_k) = (x_{k-1} + x_k)/2; on an end-stamped one it ends
    it, so x(t_k) = (x_k + x_{k+1})/2. The first (begin) or last (end) record has one neighbour and
    keeps its mean. Wind averages its square, as the reader interpolates it. A value built from a
    filled interval carries that interval's flag."""
    y = np.asarray(y, dtype=float)
    z = y ** 2 if energy else y.copy()
    out = z.copy()
    flags = qc.copy()
    if stamp == "begin":
        out[1:] = 0.5 * (z[:-1] + z[1:])
        flags[1:] = np.maximum(qc[:-1], qc[1:])
    else:
        out[:-1] = 0.5 * (z[:-1] + z[1:])
        flags[:-1] = np.maximum(qc[:-1], qc[1:])
    return (np.sqrt(out) if energy else out), flags
