# SPDX-License-Identifier: Apache-2.0
"""The targets: each compares a series of observations with the same quantity of the offline stand,
on the cohorts it sees, each counted by its full-canopy leaf area (phenology_reference.toml,
[target]). A target's residuals are (model - observation) / sigma on the scored observations.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


class TargetError(ValueError):
    pass


def year_minmax(x: pd.Series) -> pd.Series:
    """Each year's values scaled between that year's least and most (NaN where they are equal)."""
    g = x.groupby(x.index.year)
    lo, hi = g.transform("min"), g.transform("max")
    return ((x - lo) / (hi - lo).where(hi > lo)).astype(float)


def crossing_doy(x: pd.Series, level: float, direction: str, doy_min: int, doy_max: int) -> float:
    """The first day of the year (interpolated) within [doy_min, doy_max] the series crosses
    `level` going `direction`; doy_max when it does not."""
    g = x[(x.index.day_of_year >= doy_min) & (x.index.day_of_year <= doy_max)].dropna()
    v, d = g.to_numpy(), g.index.day_of_year.to_numpy()
    if direction == "up":
        k = np.flatnonzero((v[:-1] < level) & (v[1:] >= level))
    else:
        k = np.flatnonzero((v[:-1] > level) & (v[1:] <= level))
    if not len(k):
        return float(doy_max)
    i = k[0]
    return float(d[i] + (level - v[i]) * (d[i + 1] - d[i]) / (v[i + 1] - v[i]))


class Target:
    """One [[target]] on one set of drivers: the cohorts it sees, its observations, and its
    model side."""

    def __init__(self, t: dict, drivers, score_from):
        self.t, self.name, self.kind = t, t["name"], t["kind"]
        d = drivers
        sel = np.isin(d.pft, t["pfts"]) & (d.height >= t["min_height"]) & (d.height < t["max_height"])
        self.cols = np.flatnonzero(sel)
        self.w = d.area[self.cols] * d.lai_full[self.cols]
        if not len(self.cols) or self.w.sum() <= 0.0:
            raise TargetError(f"target {self.name}: no cohort with leaves among PFTs {t['pfts']} "
                              f"from {t['min_height']} to {t['max_height']} m")
        self.full = float(self.w.sum())
        self.index = pd.DatetimeIndex(d.dates)
        start, end = pd.Timestamp(score_from), self.index[-1]
        if self.kind == "annual_loss":
            years = [y for y in sorted(set(self.index.year)) if pd.Timestamp(f"{y}-01-01") >= start
                     and pd.Timestamp(f"{y}-12-31") <= end]
            if not years:
                raise TargetError(f"target {self.name}: no whole year scored")
            self.years = years
            self.obs = np.full(len(years), float(t["value"]))
            return
        obs = pd.read_csv(t["file"])
        if self.kind == "date":
            if not {"year", "doy"} <= set(obs.columns):
                raise TargetError(f"{t['file']}: a date target's file has columns year, doy")
            obs = obs[(obs.year >= start.year) & (obs.year <= end.year)]
            obs = obs[np.array([pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=t["doy_min"] - 1) >= start
                                for y in obs.year], dtype=bool)]
            self.years, self.obs = obs.year.to_numpy(int), obs.doy.to_numpy(float)
        else:
            if not {"date", "value"} <= set(obs.columns):
                raise TargetError(f"{t['file']}: a {self.kind} target's file has columns date, value")
            obs = obs.assign(date=pd.to_datetime(obs.date)).dropna(subset=["value"])
            last = end - pd.Timedelta(days=t.get("window_days", 1) - 1)
            obs = obs[(obs.date >= start) & (obs.date <= last)]
            self.dates, self.obs = pd.DatetimeIndex(obs.date), obs.value.to_numpy(float)
            self.pos = self.index.get_indexer(self.dates)
            if (self.pos < 0).any():
                raise TargetError(f"target {self.name}: {self.dates[self.pos < 0][0].date()} is not a day of the run")
            if self.kind == "litter_fraction":
                self.leaf_years = self._leaf_years()
        if not len(self.obs):
            raise TargetError(f"target {self.name}: no observation from {start.date()} to {end.date()}")

    def _leaf_years(self):
        """Per observed leaf year: the record of the day before it starts and of its last observation."""
        out = []
        for y in sorted(set(self.dates.year)):
            begin = pd.Timestamp(f"{y}-{self.t['year_start']}")
            if begin < self.index[0] + pd.Timedelta(days=1):
                raise TargetError(f"target {self.name}: the leaf year {y} starts before the run")
            last = self.dates[self.dates.year == y].max()
            out.append((y, self.index.get_loc(begin - pd.Timedelta(days=1)), self.index.get_loc(last)))
        return out

    def model(self, cover: np.ndarray, loss: np.ndarray) -> np.ndarray:
        """The model's value at each observation."""
        lai = pd.Series(self.w @ cover[self.cols], index=self.index)
        t = self.t
        if self.kind == "fraction":
            x = lai / self.full if t["scale"] == "none" else year_minmax(lai)
            n = int(t["window_days"])
            x = x.rolling(n).mean().shift(-(n - 1)) if n > 1 else x
            return x.to_numpy()[self.pos]
        if self.kind == "fallen":
            peak = lai.groupby(lai.index.year).transform("max")
            fallen = (1.0 - lai / peak.where(peak > 0)).where(lai.index.day_of_year >= t["from_doy"], 0.0)
            return fallen.to_numpy()[self.pos]
        if self.kind == "date":
            x = year_minmax(lai)
            return np.array([crossing_doy(x[str(y)], t["level"], t["direction"], t["doy_min"], t["doy_max"])
                             for y in self.years])
        lost = self.w @ loss[self.cols]
        if self.kind == "annual_loss":
            per_year = pd.Series(lost, index=self.index).groupby(self.index.year).sum()
            return per_year.reindex(self.years).to_numpy() / self.full
        cum = np.cumsum(lost)
        out = np.empty(len(self.obs))
        for y, i0, i1 in self.leaf_years:
            m = self.dates.year == y
            out[m] = (cum[self.pos[m]] - cum[i0]) / max(cum[i1] - cum[i0], 1e-12)
        return out

    def residuals(self, cover: np.ndarray, loss: np.ndarray) -> np.ndarray:
        r = (self.model(cover, loss) - self.obs) / self.t["sigma"]
        return np.where(np.isfinite(r), r, 10.0)       # an undefined model value (an empty year) counts as far off

    def score(self, cover: np.ndarray, loss: np.ndarray) -> float:
        return float(self.t["weight"] * np.mean(self.residuals(cover, loss) ** 2))

    def table(self, cover: np.ndarray, loss: np.ndarray) -> pd.DataFrame:
        """The observations beside the model, for the reports and figures."""
        when = (self.dates if hasattr(self, "dates") else pd.Index(self.years, name="year"))
        return pd.DataFrame({"target": self.name, "when": when, "observed": self.obs,
                             "model": self.model(cover, loss)})
