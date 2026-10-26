#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 02: one ED-MORE phenology kernel, four forests, four leaf habits.

  Harvard Forest  deciduous broadleaf: warmth and day length
  Hyytiala        Scots pine (evergreen): warmth and the bright hours; needle fall stops at a floor
  BCI             a leaf exchanger: senescence on the bright hours of the dry season, canopy kept full
  Palo Verde      a drought-deciduous dry forest: predawn leaf water potential and day length

The kernel and the leaf rule are the compiled Fortran of the coupled model, reached through
edmore.plant.pheno, with carbon never limiting the flush; the cues, switches and leaf habits are
docs/science/plant_phenology.md §2-§6. Light at every site is ERA5-Land through ED-MORE's forcing
reader (drivers/era5_par_hours_*.csv.gz).

Usage (after fetch_phenology_data.py; build libedmore.so first, see the README):
  python run_phenology.py                  # run fitted_parameters.json: print scores, write the figure
  python run_phenology.py --fit            # refit every site first
  python run_phenology.py --fit bci        # refit one site
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import edmore  # noqa: F401  -- an installed package (pip install python/)
except ImportError:
    sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))  # or the source tree
import edmore.plant.pheno as pheno                                   # noqa: E402

DATA = os.path.join(HERE, "data")
DRIVERS = os.path.join(HERE, "drivers")
FIT_FILE = os.path.join(HERE, "fitted_parameters.json")
FIGURE = os.path.join(HERE, "canopy_phenology.png")
Cue = pheno.Cue


def data(name):
    return os.path.join(DATA, name)


def driver(name):
    return os.path.join(DRIVERS, name)


def crossing_day(series, first_doy, last_doy):
    """Per year, the (interpolated) day of year the series first rises through 0.5."""
    out = {}
    for year, g in series.groupby(series.index.year):
        g = g[(g.index.dayofyear >= first_doy) & (g.index.dayofyear <= last_doy)]
        k = np.flatnonzero(np.diff(np.sign(g.values - 0.5)) > 0)
        if len(k):
            i = k[0]
            x0, x1 = g.index[i].dayofyear, g.index[i + 1].dayofyear
            out[year] = x0 + (0.5 - g.values[i]) * (x1 - x0) / (g.values[i + 1] - g.values[i])
    return pd.Series(out)


def timing_skill(obs, mod):
    j = obs.index.intersection(mod.index)
    return float(np.corrcoef(obs[j], mod[j])[0, 1]), float(np.sqrt(np.mean((obs[j] - mod[j]) ** 2)))


def share_by_month(daily, years):
    """Each month's share of its calendar year's leaf fall, averaged over the years."""
    daily = daily[daily.index.year.isin(list(years))]
    m = daily.groupby([daily.index.year, daily.index.month]).sum()
    share = m / m.groupby(level=0).transform("sum")
    return share.groupby(level=1).mean().reindex(range(1, 13), fill_value=0.0)


class LightHours:
    """The hours a day of PAR above any par_min. Each column of drivers/era5_par_hours_<site>.csv.gz
    counts the 900 s steps of a UTC day above one value of a grid; between grid values the hours are
    interpolated in log par_min."""

    def __init__(self, key, dates):
        t = pd.read_csv(driver(f"era5_par_hours_{key}.csv.gz"), comment="#", parse_dates=["date"])
        t = t.set_index("date").reindex(dates)
        if t.isna().any().any():
            raise ValueError(f"{key}: the hours of light do not cover {dates[0]:%Y-%m-%d} - {dates[-1]:%Y-%m-%d}")
        self.log_grid = np.log([float(c.rsplit("_", 1)[1]) for c in t.columns])
        self.hours = t.values * 0.25

    def __call__(self, par_min):
        x = np.clip(np.log(max(par_min, 1e-9)), self.log_grid[0], self.log_grid[-1])
        i = min(int(np.searchsorted(self.log_grid, x, side="right")) - 1, len(self.log_grid) - 2)
        w = (x - self.log_grid[i]) / (self.log_grid[i + 1] - self.log_grid[i])
        return (1.0 - w) * self.hours[:, i] + w * self.hours[:, i + 1]


def collection_amounts(litter, rounds):
    """Model litter [canopy fraction] over each collection (start inclusive, end exclusive)."""
    cum = np.r_[0.0, np.cumsum(litter.values)]
    return cum[litter.index.get_indexer(rounds.end)] - cum[litter.index.get_indexer(rounds.start)]


def cumulative_errors(mod, obs, years):
    """Each calendar year's cumulative fraction of its fall at the collection ends: model minus traps.
    A single big collection cannot dominate it, and the traps' unit cancels."""
    out = []
    for y in sorted(set(years)):
        k = years == y
        cm, co = np.cumsum(mod[k]), np.cumsum(obs[k])
        if cm[-1] <= 0.0:                         # a trial with no litter that year
            return np.ones(1)
        out.append(cm / cm[-1] - co / co[-1])
    return np.concatenate(out)


class Site:
    """A site's daily drivers and observations; simulate() runs the kernel over them."""
    key = name = obs_label = ""
    lat = cover0 = 0.0
    window = ("", "")                             # the years the figure shows
    BASE = {}                                     # the cue masks
    FIXED = {}                                    # parameters held, not fitted
    dates = tair = psi = light = lai_obs = None

    def params(self, values):
        return pheno.Params(**{**self.BASE, **self.FIXED, **values})

    def simulate(self, values):
        p = self.params(values)
        ph = pheno.Phenology(p, leaf_cover=self.cover0)
        n = len(self.dates)
        uses_light = (int(p.flush_cue_mask) | int(p.shed_cue_mask)) & int(Cue.LIGHT)
        hours = self.light(p.par_min) if uses_light else np.zeros(n)
        tair = self.tair if self.tair is not None else np.full(n, 298.15)
        psi = self.psi if self.psi is not None else np.zeros(n)
        days = [ph.step(temp_day=t, par_hours=h, predawn_leaf_psi=w, doy=d)
                for t, h, w, d in zip(tair, hours, psi, self.dates.dayofyear)]
        return pd.DataFrame({"cover": [x.leaf_cover for x in days],
                             "litter": [x.senescence + x.background for x in days],
                             "flush": [x.leaf_flush_tendency for x in days],
                             "shed": [x.leaf_shed_tendency for x in days]}, index=self.dates)


class LitterSite(Site):
    """A site scored on litter traps: each calendar year's cumulative fall, and the canopies a year."""
    rounds = None                                 # DataFrame: start, end, days, mass (any unit)
    SCORED_YEARS = ()
    ANNUAL_LEAF_FALL = 1.0

    def annual(self, sim):
        return float(sim.litter.groupby(sim.index.year).sum().reindex(list(self.SCORED_YEARS)).mean())

    def litter_errors(self, sim):
        r = self.rounds[self.rounds.end.dt.year.isin(self.SCORED_YEARS)]
        return cumulative_errors(collection_amounts(sim.litter, r), r.mass.values, r.end.dt.year.values)

    def loss(self, sim):
        return float(np.mean(self.litter_errors(sim) ** 2)
                     + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2)

    def leaf_fall_by_month(self, sim):
        """Each month's share of the year's leaf fall: traps (each collection spread evenly over its
        days) and model."""
        r = self.rounds
        daily = pd.concat([pd.Series(m / max(d, 1), index=pd.date_range(s, e - pd.Timedelta(days=1)))
                           for s, e, d, m in zip(r.start, r.end, r.days, r.mass)])
        return share_by_month(daily, self.SCORED_YEARS), share_by_month(sim.litter, self.SCORED_YEARS)

    def scores(self, sim):
        mod = collection_amounts(sim.litter, self.rounds) / self.rounds.days.values
        obs = (self.rounds.mass / self.rounds.days).values
        return {"RMSE cumulative litter": float(np.sqrt(np.mean(self.litter_errors(sim) ** 2))),
                "litter rate r, collection by collection": float(np.corrcoef(mod, obs)[0, 1]),
                "leaf litter [canopies/yr]": self.annual(sim),
                "lowest leaf cover": float(sim.cover.min())}


def temperate_free():
    """The temperature and rate parameters fitted at Harvard Forest."""
    return {"flush_degree_days": (5.0, 1000.0),       # [K day] warmth centre
            "shed_base_temp": (278.0, 300.0),         # [K]
            "shed_degree_days": (5.0, 600.0),         # [K day] cold centre
            "flush_rate_max": (1 / 60, 1 / 3),        # [1/day]
            "shed_rate_max": (0.0005, 1 / 3),         # [1/day]
            "leaf_turnover_rate": (0.0, 0.6)}         # [1/yr]


class HarvardForest(Site):
    key, name = "harvard_forest", "Harvard Forest, deciduous broadleaf"
    obs_label = "HF003 tagged trees"
    lat, cover0 = 42.538, 0.0
    window = ("2016", "2018")
    first, last = 2003, 2023                          # 2003 grows the first canopy (spin-up)
    BASE = dict(flush_cue_mask=Cue.TEMP | Cue.LIGHT, shed_cue_mask=Cue.TEMP | Cue.LIGHT)
    #----- Leaf litter per year as a share of the canopy: a deciduous canopy is built once a year
    #      and every leaf falls. Without it the fit can keep flushing in October, when the warmth
    #      sum is still high, and drop nearly two canopies a year while matching every timing.
    ANNUAL_LEAF_FALL = 1.0
    CONIFERS = {"hemlock", "w.pine", "r.pine", "spruce"}
    NOT_LEAVES = ("twig", "bark", "fruit", "flower", "bud", "acorn", "cone", "seed", "non.leaf")

    def __init__(self):
        y0, y1 = self.first, self.last
        met = pd.read_csv(data("hf001-06-daily-m.csv"), parse_dates=["date"], usecols=["date", "airt"])
        met = met.set_index("date").loc[f"{y0}":f"{y1}"].asfreq("D").interpolate()
        self.dates, self.tair = met.index, met.airt.values + 273.15
        self.light = LightHours(self.key, self.dates)
        #----- MODIS LAI as canopy fullness: a 3-point running median, scaled each year between
        #      its winter level (the evergreen understorey and hemlock) and its summer top.
        lai = pd.read_csv(data("harvard_modis_lai.csv"), parse_dates=["date"]).set_index("date").lai
        lai = lai.dropna().rolling(3, center=True).median().dropna().loc[f"{y0 + 1}":f"{y1}"]
        rel = []
        for _, g in lai.groupby(lai.index.year):
            doy = g.index.dayofyear
            winter = g[(doy <= 90) | (doy >= 330)].median()
            summer = g[(doy >= 170) & (doy <= 240)].quantile(0.9)
            rel.append(((g - winter) / (summer - winter)).clip(-0.1, 1.2))
        self.lai_obs = pd.concat(rel)
        #----- HF003: the mean fraction of leaves fallen on the tagged trees, by date. ---------#
        fall = pd.read_csv(data("hf003-04-fall.csv"), parse_dates=["date"]).dropna(subset=["lfall"])
        fall = fall[(fall.date.dt.year > y0) & (fall.date.dt.year <= y1)]
        self.fall = fall.groupby("date").lfall.mean().clip(0, 100) / 100.0
        #----- HF069: broadleaf litter in the EMS baskets, cumulative within each leaf year. ----#
        lit = pd.read_csv(data("hf069-05-litter.csv"), parse_dates=["date"])
        lit = lit[(lit.site == "ems") & (lit.year > y0) & (lit.year < y1)]
        spp = lit.spp.str.lower()
        lit = lit[~spp.isin(self.CONIFERS) & ~spp.str.contains("|".join(self.NOT_LEAVES))]
        b = lit.groupby(["year", "date"]).wt_g.sum().reset_index()
        b["frac"] = b.groupby("year").wt_g.cumsum() / b.groupby("year").wt_g.transform("sum")
        self.basket = b

    def free(self):
        """Fitted parameters and bounds. par_min stays low (the hours of light ~ the day length)."""
        longest = pheno.daylength(self.lat, 172)
        return {**temperate_free(),
                "log10_par_min": (0.0, 2.0),                     # 1-100 umol/m2/s
                "flush_light_hours": (8.0, longest - 1.0),       # [h/day]
                "flush_light_sharpness": (0.5, 8.0),             # [1/h]
                "shed_light_hours": (6.0, longest),              # [h/day]
                "light_window": (1.0, 15.0)}                     # [day]

    def fallen(self, sim):
        """Fraction of the year's peak canopy lost, from August on (as HF003 counts it)."""
        peak = sim.cover.groupby(sim.index.year).transform("max")
        return (1.0 - sim.cover / peak.clip(lower=1e-6)).where(sim.index.dayofyear >= 213, 0.0)

    def basket_fraction(self, sim):
        """Cumulative fraction of each leaf year's model litter (from 1 April) at the basket dates."""
        cum, out = sim.litter.cumsum(), []
        for year, g in self.basket.groupby("year"):
            start, end = cum.asof(pd.Timestamp(f"{year}-04-01")), cum.asof(g.date.max())
            out.append((cum.reindex(g.date, method="nearest").values - start) / max(end - start, 1e-9))
        return np.concatenate(out)

    def errors(self, sim):
        return {"MODIS LAI": sim.cover.reindex(self.lai_obs.index).values - self.lai_obs.values,
                "leaf fall": self.fallen(sim).reindex(self.fall.index).values - self.fall.values,
                "baskets": self.basket_fraction(sim) - self.basket.frac.values}

    def annual(self, sim):
        return float(sim.litter.loc[f"{self.first + 1}":].sum() / (self.last - self.first))

    def loss(self, sim):
        return (sum(np.nanmean(e ** 2) for e in self.errors(sim).values())
                + 4.0 * (self.annual(sim) - self.ANNUAL_LEAF_FALL) ** 2)

    def leaf_fall_by_month(self, sim):
        """Each month's share of the year's leaf fall: the tagged trees' fraction fallen,
        interpolated daily from 0 on 1 August, and the model's litter."""
        parts = []
        for y, f in self.fall.groupby(self.fall.index.year):
            days = pd.date_range(f"{y}-08-01", f"{y}-12-31")
            f = pd.concat([pd.Series([0.0], index=[days[0]]), f[f.index > days[0]]])
            frac = f.reindex(f.index.union(days)).interpolate("time").reindex(days).ffill()
            parts.append(frac.diff().fillna(0.0).clip(lower=0.0))
        years = range(self.first + 1, self.last + 1)
        return share_by_month(pd.concat(parts), years), share_by_month(sim.litter, years)

    def scores(self, sim):
        out = {f"RMSE {k}": float(np.sqrt(np.nanmean(e ** 2))) for k, e in self.errors(sim).items()}
        spring = timing_skill(crossing_day(self.lai_obs, 60, 200), crossing_day(sim.cover, 60, 200))
        autumn = timing_skill(crossing_day(self.fall, 213, 366), crossing_day(self.fallen(sim), 213, 366))
        out["spring half-green r"], out["spring half-green RMSE [d]"] = spring
        out["autumn half-fallen r"], out["autumn half-fallen RMSE [d]"] = autumn
        out["leaf litter [canopies/yr]"] = self.annual(sim)
        return out


class Hyytiala(LitterSite):
    key, name = "hyytiala", "Hyytiälä, Scots pine"
    obs_label = "needle traps (ICOS)"
    lat, cover0 = 61.8475, 1.0
    window = ("2020", "2022")
    start, end = "2018-01-01", "2024-07-31"
    BASE = dict(flush_cue_mask=Cue.TEMP | Cue.LIGHT, shed_cue_mask=Cue.TEMP | Cue.LIGHT)
    #----- Held, not fitted: the needle-fall timing does not constrain the flush. ---------------#
    FIXED = dict(flush_degree_days=88.0, flush_rate_max=0.2866)
    SCORED_YEARS = (2019, 2020, 2021, 2022, 2023)
    #----- Needle fall per year as a share of the canopy: southern-Finnish Scots pine keeps
    #      3.4-4.2 needle cohorts and needles live about three years (Pensa & Jalkanen 1999,
    #      Silva Fennica 33:654), so about 0.3 of it falls each year.
    ANNUAL_LEAF_FALL = 0.30

    def __init__(self):
        met = pd.read_csv(data("hyytiala_fluxnet_dd.csv"), usecols=["TIMESTAMP", "TA_F"], na_values=[-9999])
        met["date"] = pd.to_datetime(met.TIMESTAMP.astype(str), format="%Y%m%d")
        tair = met.set_index("date").loc[self.start:self.end].TA_F.interpolate()
        self.dates, self.tair = tair.index, tair.values + 273.15
        self.light = LightHours(self.key, self.dates)
        lit = self.needle_litter(data("hyytiala_ancillary.csv"))
        first = self.dates[0] + pd.Timedelta(days=200)              # after a 200-day spin-up
        self.rounds = lit[(lit.start >= first) & (lit.end <= self.dates[-1])].reset_index(drop=True)

    def free(self):
        """Bounded so that no single cloudy spell can trigger needle fall: a light window of at
        least a week and gentle light switches."""
        return {"shed_base_temp": (278.0, 300.0), "shed_degree_days": (5.0, 600.0),
                "shed_rate_max": (0.0005, 1 / 3), "leaf_turnover_rate": (0.0, 0.6),
                "min_leaf_cover": (0.4, 0.95),
                "log10_par_min": (1.7, 2.78),                    # 50-600 umol/m2/s
                "flush_light_hours": (0.0, 16.0), "flush_light_sharpness": (0.5, 2.0),
                "shed_light_hours": (0.0, 16.0), "shed_light_sharpness": (-2.0, -0.25),
                "light_window": (7.0, 60.0)}

    @staticmethod
    def needle_litter(path):
        """ICOS ancillary litter: the plot-mean Scots pine foliage litter, per collection interval."""
        a = pd.read_csv(path, dtype=str)
        a = a[a.VARIABLE_GROUP == "GRP_LITTER"]
        g = a.pivot_table(index="GROUP_ID", columns="VARIABLE", values="DATAVALUE", aggfunc="first")
        g = g[(g.LITTER_DEBRIS_STATISTIC == "Mean") & (g.LITTER_DEBRIS_ORGAN == "Foliage")]
        g = g[~g.LITTER_DEBRIS_APPROACH.fillna("").str.contains("Yearly")]
        g = g[g.LITTER_DEBRIS_SPP.fillna("").str.startswith("Pinus sylvestris")]
        end = pd.to_datetime(g.LITTER_DEBRIS_DATE.fillna(g.LITTER_DEBRIS_DATE_END), format="%Y%m%d")
        start = pd.to_datetime(g.LITTER_DEBRIS_DATE_START, format="%Y%m%d", errors="coerce")
        out = pd.DataFrame({"end": end.values, "start": start.values,
                            "mass": g.LITTER_DEBRIS.astype(float).values * 1000.0})   # [g DM m-2]
        out = out.sort_values("end").drop_duplicates("end").reset_index(drop=True)
        out["start"] = out.start.fillna(out.end.shift(1))
        out["days"] = (out.end - out.start).dt.days
        return out.dropna(subset=["start"])


class BCI(LitterSite):
    """A light-driven leaf exchanger: it flushes while its water threshold is never reached and
    senesces on many bright hours, so its leaves turn over in the dry season and the canopy stays
    full."""
    key, name = "bci", "Barro Colorado Island, leaf exchanger"
    obs_label = "litter traps (GLiMP)"
    lat, cover0 = 9.15, 1.0
    window = ("2014", "2016")
    start, end = "2012-07-03", "2017-08-31"
    BASE = dict(flush_cue_mask=Cue.WATER, shed_cue_mask=Cue.LIGHT | Cue.WATER)
    #----- A water threshold far below the surrogate's range (-0.05 to -0.9 MPa): the exchanger
    #      never sheds for drought, and its water switch keeps it flushing. --------------------#
    FIXED = dict(leaf_psi_tlp=-1.5, flush_water_sum=3.0, flush_water_sharpness=6.667,
                 shed_water_sum=1.0, shed_water_sharpness=20.0)
    SCORED_YEARS = (2013, 2014, 2015, 2016)
    #----- Leaf litter per year as a share of the canopy: the leaves falling at BCI each year have
    #      about the canopy's area, 7.3 m2 per m2 of ground (Leigh 1999, ORNL NPP data set BRR).
    ANNUAL_LEAF_FALL = 1.0

    def __init__(self):
        tower = pd.read_csv(data("BCI_v5.1.csv"), parse_dates=["date"], na_values=["NaN"]).set_index("date")
        daily = tower[["tair", "SWC"]].resample("D").mean().interpolate()
        self.soil_psi_fit(daily.SWC)
        daily = daily.loc[self.start:self.end]
        self.dates = daily.index
        self.tair = daily.tair.values + 273.15
        self.psi = self.soil_psi(daily.SWC.values)
        self.light = LightHours(self.key, self.dates)
        #----- GLiMP control plots: mean fine litter per trap and collection. -------------------#
        lit = pd.read_csv(data("glimp_litterfall.csv"), encoding="utf-8-sig")
        lit = lit[lit.treat == "CT"]
        lit["date"] = pd.to_datetime(lit.date, dayfirst=True)
        c = lit.groupby("date").agg(mass=("mass.g.trap", "mean"), days=("days.accum", "median"))
        c["end"] = c.index
        c["start"] = c.end - pd.to_timedelta(c.days, unit="D")
        first = self.dates[0] + pd.Timedelta(days=180)              # after a half-year spin-up
        self.rounds = c[(c.start >= first) & (c.end <= self.dates[-1])].reset_index(drop=True)

    def soil_psi_fit(self, swc):
        """A Campbell retention curve, psi = -exp(a) theta_g^b, fitted to the soil water content
        and potential Kupers et al. (2019) sampled together; the tower's volumetric water content
        converts to gravimetric by its ratio to the plot's on their four sampling periods."""
        pairs = pd.concat([pd.read_csv(data(f), sep="\t") for f in
                           ("kupers_soil_moisture_mapping.txt", "kupers_soil_moisture_small_scale.txt")])
        pairs = pairs[(pairs.swp < -0.05) & (pairs.swc > 5.0)]      # the WP4C reads ~0 when wet
        self.campbell_b, self.campbell_a = np.polyfit(np.log(pairs.swc / 100.0), np.log(-pairs.swp), 1)
        plot = pd.read_csv(data("kupers_soil_moisture_mapping.txt"), sep="\t")
        plot["date"] = pd.to_datetime(plot.date, format="%m/%d/%Y")
        periods = plot.groupby("period").agg(date=("date", "median"), swc=("swc", "median"))
        tower = [swc.loc[t - pd.Timedelta(days=3): t + pd.Timedelta(days=3)].mean() for t in periods.date]
        self.bulk_density = float(np.median(np.array(tower) / (periods.swc.values / 100.0)))

    def soil_psi(self, swc):
        theta_g = swc / self.bulk_density
        return -np.exp(self.campbell_a + self.campbell_b * np.log(theta_g))

    def free(self):
        """Rates, the evergreen floor and the light trigger. par_min is held to bright light:
        free, it drifts to the photoperiod, a calendar with no year-to-year signal."""
        return {"flush_rate_max": (1 / 60, 1 / 3), "shed_rate_max": (0.0005, 1 / 3),
                "leaf_turnover_rate": (0.0, 3.0), "min_leaf_cover": (0.8, 0.95),
                "log10_par_min": (2.0, 3.08),                    # 100-1200 umol/m2/s
                "shed_light_hours": (0.0, 13.0),                 # [h/day]
                "shed_light_sharpness": (0.1, 4.0),              # > 0: many bright hours trigger senescence
                "light_window": (2.0, 60.0)}                     # [day]


class PaloVerde(LitterSite):
    """A drought-deciduous forest: flushing and senescence on the predawn leaf water potential and
    the hours of light. The water potential is hypothetical: an ED-MORE run's canopy predawn value for
    an evergreen stand at the site, so the fit asks what cues map that driver onto the observed
    canopy, not how this forest's own water status moves."""
    key, name = "palo_verde", "Palo Verde, drought-deciduous"
    obs_label = "litter traps (Xu et al. 2016)"
    lat, cover0 = 10.35, 1.0
    window = ("2010", "2012")
    start, end = "2008-01-01", "2013-12-31"
    BASE = dict(flush_cue_mask=Cue.WATER | Cue.LIGHT, shed_cue_mask=Cue.WATER | Cue.LIGHT)
    SCORED_YEARS = (2009, 2010, 2011, 2012, 2013)
    #----- About one canopy a year: 0.56 kg m-2 yr-1 of leaf litter over a canopy of LAI 5.5 at
    #      80-100 g m-2 of leaf is 1.0-1.3 canopies. Without the term the fit churns leaves. -----#
    ANNUAL_LEAF_FALL = 1.0

    def __init__(self):
        d = pd.read_csv(driver("palo_verde_predawn_psi.csv"), comment="#", parse_dates=["date"]).set_index("date")
        d = d.loc[self.start:self.end]
        self.dates, self.psi = d.index, d.predawn_leaf_psi.values
        self.light = LightHours(self.key, self.dates)
        #----- MODIS LAI by month, over its 95th percentile: canopy fullness. --------------------#
        lai = pd.read_csv(driver("palo_verde_modis_lai.csv"), comment="#", parse_dates=["month"]).set_index("month").lai
        self.lai_obs = (lai / lai.quantile(0.95)).clip(upper=1.1)
        #----- Leaf litter by month: each month is one collection. ------------------------------#
        r = pd.read_csv(driver("palo_verde_leaf_litter.csv"), comment="#", parse_dates=["month"])
        r["start"], r["end"] = r.month, r.month + pd.offsets.MonthBegin(1)
        r["days"] = (r.end - r.start).dt.days
        r["mass"] = r.leaf_litter_g_m2_day * r.days
        first = self.dates[0] + pd.Timedelta(days=150)              # after a five-month spin-up
        self.rounds = r[(r.start >= first) & (r.end <= self.dates[-1])].reset_index(drop=True)

    def free(self):
        """Water and light on both sides, with both light signs free."""
        return {"leaf_psi_tlp": (-3.0, -0.3),                    # [MPa]
                "flush_water_sum": (0.5, 60.0), "flush_water_sharpness": (0.5, 20.0),
                "shed_water_sum": (0.5, 60.0), "shed_water_sharpness": (0.5, 20.0),
                "flush_rate_max": (1 / 60, 1 / 3), "shed_rate_max": (0.0005, 1 / 3),
                "leaf_turnover_rate": (0.0, 1.0), "min_leaf_cover": (0.0, 0.5),
                "log10_par_min": (0.0, 3.08),                    # 1-1200 umol/m2/s
                "flush_light_hours": (0.0, 13.0), "flush_light_sharpness": (-8.0, 8.0),
                "shed_light_hours": (0.0, 13.0), "shed_light_sharpness": (-8.0, 8.0),
                "light_window": (1.0, 60.0)}

    def lai_errors(self, sim):
        cm = sim.cover.resample("MS").mean().reindex(self.lai_obs.index)
        return cm.values - self.lai_obs.values

    def loss(self, sim):
        return float(np.nanmean(self.lai_errors(sim) ** 2) + super().loss(sim))

    def scores(self, sim):
        return {"RMSE relative LAI": float(np.sqrt(np.nanmean(self.lai_errors(sim) ** 2))),
                **super().scores(sim)}


SITES = {cls.key: cls for cls in (HarvardForest, Hyytiala, BCI, PaloVerde)}
_LOADED = {}


def site(key):
    """Each site's data, loaded once per process (also in each fitting worker)."""
    if key not in _LOADED:
        _LOADED[key] = SITES[key]()
    return _LOADED[key]


def as_values(names, x):
    """Fitted vector -> parameter values; par_min is fitted as its log10."""
    v = dict(zip(names, map(float, x)))
    if "log10_par_min" in v:
        v["par_min"] = 10.0 ** v.pop("log10_par_min")
    return v


def _loss(x, key, names):
    s = site(key)
    try:
        value = s.loss(s.simulate(as_values(names, x)))
    except (ValueError, FloatingPointError):
        return 1e3
    return float(value) if np.isfinite(value) else 1e3


def fit(key, workers, seeds=4):
    """Differential evolution from several seeds, each polished by a bounded Nelder-Mead; the best
    wins, and the spread of the seeds is reported (a wide spread means the optimum is not unique)."""
    from scipy.optimize import differential_evolution, minimize
    free = site(key).free()
    names, bounds, runs = list(free), list(free.values()), []
    for seed in range(1, seeds + 1):
        r = differential_evolution(_loss, bounds, args=(key, names), popsize=25, maxiter=250, tol=1e-8,
                                   seed=seed, polish=False, workers=workers, updating="deferred")
        q = minimize(_loss, r.x, args=(key, names), method="Nelder-Mead", bounds=bounds,
                     options=dict(maxiter=2000, xatol=1e-5, fatol=1e-9))
        x, v = (q.x, q.fun) if q.fun < r.fun else (r.x, r.fun)
        runs.append((float(v), x))
        print(f"  {key} seed {seed}: {r.fun:.5f} -> polished {v:.5f}", flush=True)
    runs.sort(key=lambda t: t[0])
    print(f"  {key}: best {runs[0][0]:.5f}; seeds {', '.join(f'{v:.5f}' for v, _ in runs)}", flush=True)
    return {n: float(f"{v:.6g}") for n, v in as_values(names, runs[0][1]).items()}


def plot(runs, path):
    """Top row: each site's canopy over three years. Bottom row: when the leaves fall."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, len(runs), figsize=(4.2 * len(runs), 7.2), constrained_layout=True)
    letters = "abcdefgh"
    for j, (s, sim) in enumerate(runs):
        ax = axs[0, j]
        w = sim.loc[s.window[0]:s.window[1]]
        ax.plot(w.index, w.flush, color="tab:blue", lw=0.8, label="flush tendency")
        ax.plot(w.index, w.shed, color="tab:red", lw=0.8, ls="--", label="senescence tendency")
        ax.plot(w.index, w.cover, color="tab:green", lw=2, label="model leaf cover")
        if s.lai_obs is not None:
            lai = s.lai_obs.loc[s.window[0]:s.window[1]]
            ax.plot(lai.index, lai.values, "k.", ms=4, label="MODIS LAI, relative")
        ax.set(ylim=(-0.08, 1.15), title=f"({letters[j]}) {s.name}")
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax = axs[1, j]
        obs, mod = s.leaf_fall_by_month(sim)
        ax.bar(obs.index - 0.2, obs.values, 0.4, color="k", label="leaf fall, observed")
        ax.bar(mod.index + 0.2, mod.values, 0.4, color="tab:green", label="leaf fall, model")
        ax.set(xticks=range(1, 13), xticklabels="JFMAMJJASOND", ylim=(0, None),
               title=f"({letters[j + len(runs)]}) leaf fall: {s.obs_label}")
    axs[0, 0].set_ylabel("fraction of the full canopy")
    axs[1, 0].set_ylabel("share of the year's leaf fall")
    entries = {lab: h for ax in axs.flat for h, lab in zip(*ax.get_legend_handles_labels())}
    fig.legend(entries.values(), entries.keys(), loc="outside upper center", ncol=len(entries), fontsize=9)
    fig.savefig(path, dpi=130)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fit", nargs="*", metavar="SITE",
                    help="refit these sites (all if none given) and rewrite " + os.path.basename(FIT_FILE))
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="processes for the fit")
    args = ap.parse_args(argv)
    with open(FIT_FILE) as fh:
        values = json.load(fh)
    if args.fit is not None:
        for key in args.fit or list(SITES):
            values[key] = fit(key, args.workers)
            with open(FIT_FILE) as fh:                       # other sites may have been refitted meanwhile
                current = json.load(fh)
            current[key] = values[key]
            with open(FIT_FILE, "w") as fh:
                json.dump(current, fh, indent=2)
                fh.write("\n")
            values = current
    runs = []
    for key in SITES:
        s = site(key)
        sim = s.simulate(values[key])
        print(f"\n{s.name}")
        for name, value in {**s.FIXED, **values[key]}.items():
            print(f"  {name:26s} {value:g}")
        for name, value in s.scores(sim).items():
            print(f"  {name:46s} {value:.3f}")
        runs.append((s, sim))
    plot(runs, FIGURE)
    print(f"\nwrote {os.path.basename(FIGURE)}")


if __name__ == "__main__":
    main()
