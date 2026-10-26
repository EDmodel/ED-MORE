# SPDX-License-Identifier: Apache-2.0
"""The data rules (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2.3, §2.5, §6.2). Each works at any
tower from its own data and forcing; the calibration's settings can override every outcome.

  u* per target      within classes of the target's driver (PAR for GPP, Rnet for LE and H), the
                     flux-to-driver ratio across u* classes. The plateau test (Papale et al. 2006)
                     finds the lowest class within `plateau_fraction` of the mean of the classes
                     above; a bootstrap over days gives its spread (Barr et al. 2013). The outcome:
                       plateau   filter at the threshold
                       flat      the lowest class already passes: no filter
                       rising    still rising at the top classes: no filter (the target's
                                 observation model has to handle it)
                       too_few   too few records for the test: no filter
  processes          for each process a key acts through (wet canopy, night, snow, drought), the
                     kept records that sample it. A key whose process has almost none is fixed.
  windows            `per_year` ten-day calibration windows, one per slot of the year inside the
                     leaf-on months, each the slot's best covered (measured turbulent fluxes,
                     observed forcing); validation windows in the same slots of other years.
  seasonal runs      up to two 120-day runs, each ending at a year's deepest cumulative water
                     deficit (rain minus Priestley-Taylor evaporation, from the forcing).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DRIVERS = {"gpp": "par", "le": "rnet", "h": "rnet"}
#: kept records with rain in them or the record before sample the wet canopy
WET_LAG_RECORDS = 1
SIGMA_SB = 5.670374419e-8


# ------------------------------------------------------------------------------------------------
# u*: a diagnostic per target
# ------------------------------------------------------------------------------------------------
def _plateau(ratios: np.ndarray, valid: np.ndarray, fraction: float):
    """The plateau test on one driver class's ratios by u* class: (outcome, index of the class).
    Classes with too few records are skipped."""
    idx = np.flatnonzero(valid)
    if len(idx) < 3:
        return "too_few", None
    r = ratios[idx]
    for k in range(len(idx) - 1):
        if r[k] >= fraction * r[k + 1:].mean():
            if k == 0:
                return "flat", idx[0]
            if k >= len(idx) - 2:                 # only the top pair passes: still rising
                return "rising", None
            return "plateau", idx[k]
    return "rising", None


def _classify(sums: np.ndarray, classes, min_n: int, fraction: float):
    """Outcome and threshold from per-(driver class, u* class) sums [flux, driver, n]."""
    outcomes, thresholds = [], []
    for d in range(sums.shape[0]):
        flux, drv, n = sums[d, :, 0], sums[d, :, 1], sums[d, :, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(drv > 0, flux / drv, np.nan)
        out, k = _plateau(ratios, (n >= min_n) & np.isfinite(ratios), fraction)
        outcomes.append(out)
        thresholds.append(None if k is None else (0.0 if out == "flat" else float(classes[k])))
    usable = [o for o in outcomes if o != "too_few"]
    if not usable:
        return "too_few", None, outcomes
    if sum(o == "rising" for o in usable) > len(usable) / 2:
        return "rising", None, outcomes
    thr = float(np.median([t for o, t in zip(outcomes, thresholds) if o in ("flat", "plateau")]))
    return ("flat" if thr == 0.0 else "plateau"), thr, outcomes


def ustar_diagnostic(obs: pd.DataFrame, target: str, settings: dict, daytime_sw: float,
                     utc_offset_h: float) -> dict:
    """The u* diagnostic of one target on the whole record (measured daytime records): its outcome,
    threshold, the bootstrap's 5-50-95 % thresholds and outcome shares, and the ratio table."""
    driver = DRIVERS[target]
    classes = list(settings["classes"])
    if obs[driver].notna().sum() == 0:
        return {"outcome": "too_few", "driver": driver, "note": f"no {driver} (the site TOML declares none)"}
    flux, drv, u = (obs[c].to_numpy() for c in (target, driver, "ustar"))
    ok = (obs["sw_in"].to_numpy() > daytime_sw) & np.isfinite(flux) & np.isfinite(drv) & np.isfinite(u)
    ok &= drv >= float(settings["min_driver"][driver])
    if ok.sum() < int(settings["min_records"]):
        return {"outcome": "too_few", "driver": driver, "records": int(ok.sum()),
                "note": f"fewer than {settings['min_records']} measured daytime records"}
    day = (obs.index + pd.Timedelta(hours=utc_offset_h)).floor("D").to_numpy()[ok]
    flux, drv, u = flux[ok], drv[ok], u[ok]
    n_dc = int(settings["driver_classes"])
    dclass = np.minimum((pd.Series(drv).rank(method="first").to_numpy() - 1) * n_dc // len(drv), n_dc - 1).astype(int)
    uclass = np.clip(np.searchsorted(classes, u, side="right") - 1, 0, len(classes) - 2)
    days, day_id = np.unique(day, return_inverse=True)
    n_uc = len(classes) - 1
    #----- per (day, driver class, u* class): the sums the test needs, so a bootstrap is a re-weighting
    cell = (day_id * n_dc + dclass) * n_uc + uclass
    size = len(days) * n_dc * n_uc
    per_day = np.stack([np.bincount(cell, weights=flux, minlength=size),
                        np.bincount(cell, weights=drv, minlength=size),
                        np.bincount(cell, minlength=size).astype(float)], axis=-1)
    per_day = per_day.reshape(len(days), n_dc, n_uc, 3)
    min_n, fraction = int(settings["min_class_records"]), float(settings["plateau_fraction"])
    outcome, thr, by_driver = _classify(per_day.sum(axis=0), classes, min_n, fraction)
    rng = np.random.default_rng(int(settings["seed"]))
    boot_thr, boot_out = [], []
    for _ in range(int(settings["n_boot"])):
        w = np.bincount(rng.integers(0, len(days), len(days)), minlength=len(days)).astype(float)
        o, t, _ = _classify(np.einsum("d,dcuk->cuk", w, per_day), classes, min_n, fraction)
        boot_out.append(o)
        boot_thr.append(np.nan if t is None else t)
    boot_thr = np.array(boot_thr)
    total = per_day.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        table = [[None if total[d, k, 2] < min_n else float(total[d, k, 0] / total[d, k, 1])
                  for k in range(n_uc)] for d in range(n_dc)]
    finite = boot_thr[np.isfinite(boot_thr)]
    return {"outcome": outcome, "threshold": thr, "records": int(len(flux)), "driver": driver,
            "by_driver_class": by_driver, "ratio": table, "classes": classes,
            "driver_class_edges": [float(np.min(drv[dclass == d])) for d in range(n_dc)] + [float(np.max(drv))],
            "bootstrap": {"n": int(settings["n_boot"]),
                          "threshold_p05_p50_p95": ([float(np.percentile(finite, q)) for q in (5, 50, 95)]
                                                    if len(finite) else None),
                          "outcome_share": {o: boot_out.count(o) / len(boot_out) for o in sorted(set(boot_out))}}}


def provider_threshold(provider: dict):
    """The provider's CO2 u* threshold (the site TOML's [provider]): a number, or {year: number};
    None when it gives none."""
    v = provider.get("ustar_threshold")
    if v is None:
        return None
    if isinstance(v, dict):
        return {int(k): float(x) for k, x in v.items()}
    return float(v)


def resolve_ustar(targets: dict, obs: pd.DataFrame, provider, settings: dict, daytime_sw: float,
                  utc_offset_h: float):
    """The targets' settings with every u* rule made a number (or a table of years, or None), and
    the report: each target's diagnostic, rule and reason. A target's ustar_min is a number,
    "provider" (the provider's CO2 threshold), "diagnostic" (its own diagnostic's outcome), or
    absent (no filter)."""
    out = {k: dict(v) for k, v in targets.items()}
    report = {}
    for t, c in out.items():
        if not c["on"]:
            continue
        rule = c.get("ustar_min")
        diag = ustar_diagnostic(obs, t, settings, daytime_sw, utc_offset_h) if t in DRIVERS else None
        entry = {"rule": rule, "diagnostic": diag}
        if rule == "provider":
            if provider is None:
                raise ValueError(f"[targets.{t}].ustar_min = \"provider\" needs the provider's u* threshold "
                                 "(the site TOML's [provider].ustar_threshold)")
            c["ustar_min"] = provider
            entry["reason"] = "the provider's CO2 threshold"
        elif rule == "diagnostic":
            if diag is None:
                raise ValueError(f"[targets.{t}].ustar_min = \"diagnostic\": the diagnostic runs on gpp, le and h only")
            c["ustar_min"] = diag["threshold"] if diag["outcome"] == "plateau" else None
            entry["reason"] = {"plateau": "a plateau: filtered at it", "flat": "the ratio is flat: no filter",
                               "rising": "still rising at the top u* classes: no filter (the observation model's job)",
                               "too_few": "too few records for the diagnostic: no filter"}[diag["outcome"]]
        elif rule is not None and not isinstance(rule, (int, float, dict)):
            raise ValueError(f"[targets.{t}].ustar_min must be a number, \"provider\" or \"diagnostic\", not {rule!r}")
        else:
            entry["reason"] = "set in the calibration's settings" if rule is not None else "no u* filter"
        entry["ustar_min"] = c.get("ustar_min")
        report[t] = entry
    return out, report


def threshold_at(index: pd.DatetimeIndex, thr) -> np.ndarray:
    """A u* threshold per record: a number, or the provider's table of years (a year it does not
    list takes the median of those it does)."""
    if isinstance(thr, dict):
        fill = float(np.median(list(thr.values())))
        return np.array([thr.get(int(y), fill) for y in index.year], dtype=float)
    return np.full(len(index), float(thr))


# ------------------------------------------------------------------------------------------------
# which processes the kept data sample
# ------------------------------------------------------------------------------------------------
PROCESSES = ("wet_canopy", "night", "snow", "drought")
TURBULENT = ("le", "h", "gpp", "ustar")


def process_coverage(rows: dict, obs: pd.DataFrame, windows, daytime_sw: float) -> dict:
    """For each process, the kept records that sample it, over the given windows:
      wet_canopy  turbulent-target records with rain in them or the record before
      night       turbulent-target records at night
      snow        albedo records within a week of freezing air
      drought     LE records of the seasonal runs (the soil-drying runs)"""
    counts = {p: 0 for p in PROCESSES}
    wet = obs["rain"].fillna(0.0) > 0.0
    wet = wet | wet.shift(WET_LAG_RECORDS, fill_value=False)
    frost = obs["days_since_frost"] <= 7
    for w in windows:
        wr = rows[w.name]
        kept = set()
        for t in wr.targets:
            if t.name in TURBULENT:
                kept.update(int(i) for i in t.rows)
            if t.name == "albedo":
                counts["snow"] += int(frost.reindex(wr.index, fill_value=False).to_numpy()[t.rows].sum())
            if t.name == "le" and w.role == "seasonal":
                counts["drought"] += len(t.rows)
        kept = np.array(sorted(kept), dtype=int)
        if len(kept):
            counts["wet_canopy"] += int(wet.reindex(wr.index, fill_value=False).to_numpy()[kept].sum())
            sw = obs["sw_in"].reindex(wr.index).to_numpy()[kept]
            counts["night"] += int((sw <= daytime_sw).sum())
    return counts


def days_since_frost(tair: pd.Series, utc_offset_h: float = 0.0) -> pd.Series:
    """For each record, the whole days since the last day whose air fell to freezing (0 degC);
    a large number where the record shows none."""
    local = tair.index + pd.Timedelta(hours=utc_offset_h)
    daily_min = tair.groupby(local.floor("D")).min().asfreq("D")
    frosty = daily_min <= 273.15
    last = pd.Series(frosty.index.where(frosty), index=frosty.index).ffill()
    since = (frosty.index - pd.DatetimeIndex(last)).days.to_numpy(dtype=float, na_value=np.nan)
    since = pd.Series(np.where(np.isnan(since), 1.0e6, since), index=frosty.index)
    return pd.Series(since.reindex(local.floor("D")).to_numpy(), index=tair.index)


def fix_by_coverage(keys, coverage: dict, min_records: int) -> dict:
    """Fix every key whose process the kept data sample in fewer than min_records records; returns
    {key: reason}."""
    fixed = {}
    for p in keys:
        if p.process and coverage[p.process] < min_records:
            fixed[p.name] = (f"the kept data sample its process ({p.process}) in {coverage[p.process]} "
                             f"records, fewer than {min_records}")
    return fixed


# ------------------------------------------------------------------------------------------------
# the windows
# ------------------------------------------------------------------------------------------------
def day_scores(obs: pd.DataFrame, forcing_observed: pd.Series, daytime_sw: float, utc_offset_h: float) -> pd.Series:
    """Each local day's coverage: the share of its daytime records with H, LE and GPP measured,
    times the share of all its records with observed forcing."""
    day = (obs.index + pd.Timedelta(hours=utc_offset_h)).floor("D")
    daytime = obs["sw_in"] > daytime_sw
    turb = daytime & obs[["h", "le", "gpp"]].notna().all(axis=1)
    frac_turb = turb.groupby(day).sum() / daytime.groupby(day).sum().replace(0, np.nan)
    frac_forc = forcing_observed.reindex(obs.index, fill_value=False).astype(float).groupby(day).mean()
    s = (frac_turb * frac_forc).fillna(0.0)
    return s.asfreq("D", fill_value=0.0)


def slot_of(doy: int, per_year: int) -> int:
    return min(int((doy - 1) * per_year / 365.25), per_year - 1)


def candidate_windows(scores: pd.Series, days: int, per_year: int, leaf_on: list, first_start) -> pd.DataFrame:
    """Every span of `days` inside one slot and inside the leaf-on months, starting on or after
    first_start: its start, year, slot and mean daily coverage."""
    roll = scores.rolling(days).mean().shift(-(days - 1))
    rows = []
    for start, sc in roll.dropna().items():
        if start < first_start:
            continue
        end = start + pd.Timedelta(days=days - 1)
        if end.year != start.year:
            continue
        s0, s1 = slot_of(start.dayofyear, per_year), slot_of(end.dayofyear, per_year)
        span = pd.date_range(start, end, freq="D")
        if s0 != s1 or not set(span.month) <= set(leaf_on):
            continue
        rows.append((start, start.year, s0, float(sc)))
    return pd.DataFrame(rows, columns=["start", "year", "slot", "coverage"])


def select_windows(scores: pd.Series, settings: dict, leaf_on: list, first_start):
    """The calibration and validation windows by the rule (module docstring), from each day's
    coverage (day_scores), and the report. Returns ([(name, start, role)], report)."""
    days, per_year, need = int(settings["days"]), int(settings["per_year"]), float(settings["min_coverage"])
    cand = candidate_windows(scores, days, per_year, leaf_on, first_start)
    cand = cand[cand["coverage"] >= need]
    chosen, report = [], {"slots": {}}
    for s in range(per_year):
        c = cand[cand["slot"] == s].sort_values(["coverage", "start"], ascending=[False, True])
        if c.empty:
            report["slots"][s + 1] = {"note": f"no window covers {need} or more"}
            continue
        best = c.iloc[0]
        chosen.append((f"cal{s + 1}_{best['year']}", best["start"].to_pydatetime(), "cal"))
        entry = {"cal": [str(best["start"].date()), round(best["coverage"], 3)]}
        v = c[c["year"] != best["year"]]
        if not v.empty:
            bv = v.iloc[0]
            chosen.append((f"val{s + 1}_{bv['year']}", bv["start"].to_pydatetime(), "val"))
            entry["val"] = [str(bv["start"].date()), round(bv["coverage"], 3)]
        report["slots"][s + 1] = entry
    return chosen, report


# ------------------------------------------------------------------------------------------------
# the water-deficit index and the seasonal runs
# ------------------------------------------------------------------------------------------------
def priestley_taylor_mm(tair_k, sw, lw, psurf, albedo=0.15, emissivity=0.97, alpha=1.26) -> np.ndarray:
    """Daily Priestley-Taylor potential evaporation [mm d-1] from daily mean air temperature [K],
    shortwave and longwave down [W m-2] and pressure [Pa]; net radiation from a fixed surface
    albedo and emissivity, ground heat taken as 0 over a day, and 0 when the net radiation is."""
    tc = np.asarray(tair_k) - 273.15
    es = 0.6108 * np.exp(17.27 * tc / (tc + 237.3))                # [kPa]
    delta = 4098.0 * es / (tc + 237.3) ** 2                        # [kPa K-1]
    gamma = 0.000665 * np.asarray(psurf) / 1000.0                   # [kPa K-1]
    rn = (1.0 - albedo) * np.asarray(sw) + np.asarray(lw) - emissivity * SIGMA_SB * np.asarray(tair_k) ** 4
    lam = 2.45e6                                                   # [J kg-1]
    return np.maximum(alpha * delta / (delta + gamma) * rn * 86400.0 / lam, 0.0)


def water_deficit(daily: pd.DataFrame) -> pd.Series:
    """The cumulative water deficit [mm, <= 0]: a bucket that loses the potential evaporation and
    gains the rain each day, and never holds more than it started with (D = min(0, D + P - PET))."""
    d, out = 0.0, []
    for p, e in zip(daily["rain_mm"].to_numpy(), daily["pet_mm"].to_numpy()):
        d = min(0.0, d + p - e)
        out.append(d)
    return pd.Series(out, index=daily.index, name="deficit_mm")


def forcing_daily(path, grid: int) -> pd.DataFrame:
    """Daily means of the forcing file (UTC days): air temperature, radiation, pressure, and the
    day's rain [mm] and Priestley-Taylor evaporation [mm]."""
    from netCDF4 import Dataset
    from tower import forcing_times
    with Dataset(path) as ds:
        when = forcing_times(ds)
        col = {v: np.asarray(ds[v][:, grid - 1], dtype=float) for v in ("Tair", "SWdown", "LWdown", "PSurf", "Rainf")}
    df = pd.DataFrame(col, index=when).resample("1D").mean().dropna()
    df["rain_mm"] = df["Rainf"] * 86400.0
    df["pet_mm"] = priestley_taylor_mm(df["Tair"], df["SWdown"], df["LWdown"], df["PSurf"])
    return df


def seasonal_runs(deficit: pd.Series, settings: dict, leaf_on: list, first_start, scores: pd.Series,
                  min_coverage: float):
    """Up to max_runs runs of `days` days, each ending at a year's deepest deficit inside the
    leaf-on months, deepest years first. A year is skipped when its deepest deficit is under
    min_deficit_mm, or when the run's mean daily coverage (day_scores) is under min_coverage -- a
    run the data barely see scores nothing. Returns ([(name, start, days)], report)."""
    days, n_max, need = int(settings["days"]), int(settings["max_runs"]), float(settings["min_deficit_mm"])
    d = deficit[deficit.index.month.isin(leaf_on)]
    found = []
    for year, g in d.groupby(d.index.year):
        end = g.idxmin()
        start = end - pd.Timedelta(days=days)
        span = pd.date_range(start, end, freq="D")
        cover = float(scores.reindex(span, fill_value=0.0).mean())
        found.append({"year": int(year), "end": str(end.date()), "deficit_mm": float(-g.min()),
                      "coverage": round(cover, 3),
                      "usable": bool(-g.min() >= need and start >= first_start and set(span.month) <= set(leaf_on)
                                     and cover >= min_coverage),
                      "start": start})
    usable = sorted([f for f in found if f["usable"]], key=lambda f: -f["deficit_mm"])[:n_max]
    runs = [(f"seasonal_{f['year']}", f["start"].to_pydatetime(), days)
            for f in sorted(usable, key=lambda f: f["start"])]
    report = {"years": [{k: v for k, v in f.items() if k != "start"} for f in found],
              "runs": [[n, str(s.date()), dd] for n, s, dd in runs]}
    if not runs:
        report["note"] = (f"no year's deficit reaches {need:g} mm inside the leaf-on months, after the chains' lead "
                          f"and with coverage {min_coverage:g}: no seasonal runs, the drought keys are fixed")
    return runs, report


# ------------------------------------------------------------------------------------------------
# what the windows cover
# ------------------------------------------------------------------------------------------------
def range_coverage(record: dict, kept: dict) -> dict:
    """For each variable, the share of the record's 1-99 % range that the kept records span."""
    out = {}
    for k, x in record.items():
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        y = np.asarray(kept.get(k, []), dtype=float)
        y = y[np.isfinite(y)]
        if len(x) < 10 or len(y) == 0:
            out[k] = None
            continue
        lo, hi = np.percentile(x, [1, 99])
        a, b = max(lo, y.min()), min(hi, y.max())
        out[k] = {"record": [float(lo), float(hi)], "kept": [float(y.min()), float(y.max())],
                  "share": float(max(0.0, b - a) / (hi - lo)) if hi > lo else 1.0}
    return out
