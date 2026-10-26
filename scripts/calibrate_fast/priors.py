# SPDX-License-Identifier: Apache-2.0
"""The priors from the site's climate (protocol §3.5). Everything here
comes from the forcing, never from the tower's fluxes.

  growth climate   over the leaf-on months: the mean air temperature (Kattge & Knorr's growth
                   temperature), and the daytime means (incoming shortwave above 10 W m-2) of air
                   temperature, VPD, pressure and PAR; CO2 from the run's own CO2 setting
  stomatal_g1      the least-cost value: Medlyn's g1 equals xi = sqrt(beta (K + Gamma*) / (1.6 eta*)),
                   beta = 146 (Prentice et al. 2014; Stocker et al. 2020), at the daytime climate with
                   ED-MORE's own Rubisco kinetics; eta* is water's viscosity relative to 25 degC
  vcmax25          the coordination condition (Wang et al. 2017; Smith et al. 2019), solved with
                   ED-MORE's own leaf equations (edmore.plant.leaf): the vcmax25 at which the Rubisco- and
                   light-limited rates are equal at the daytime climate, with ci from the least-cost
                   ratio, the configured phi_psii, theta_J and leaf absorptance, and ED-MORE's
                   temperature response at the Kattge & Knorr values below
  Kattge & Knorr   (2007) for the growth temperature T [degC]: Jmax/Vcmax = 2.59 - 0.035 T,
                   ds_vcmax = 668.39 - 1.07 T, ds_jmax = 659.70 - 0.75 T

Each EEO value is a prior's centre with a wide range (log-sd 0.5): one rule for every plant type.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

R_GAS = 8.314462618          # [J mol-1 K-1]
T25 = 298.15                 # [K]
BETA_LEAST_COST = 146.0      # the unit-cost ratio of carboxylation to transpiration (Stocker et al. 2020)
PAR_PER_SW = 2.1             # [umol J-1] PAR photons per joule of total shortwave (ED-MORE's constant-path blend)
#: the O2 mole fraction and pressure gstar25 is stated at, as in ED-MORE's leaf kernel: Gamma* as a mole
#: fraction does not depend on pressure (docs/science/leaf_gas_exchange.md eq. (3))
O2_REF_GSTAR, P_REF_GSTAR = 0.209, 1.0e5          # [mol mol-1], [Pa]
EEO_KEYS = {"pft.stomatal_g1": "g1", "pft.vcmax25": "vcmax25"}
KATTGE_KNORR = {"pft.jmax_vcmax_ratio": (2.59, 0.035), "leaf_physiology.ds_vcmax": (668.39, 1.07),
                "leaf_physiology.ds_jmax": (659.70, 0.75)}
#: ED-MORE's leaf-physiology settings the priors read, and the model's own defaults
LEAF_DEFAULTS = {"kc25": 40.49, "ko25": 27840.0, "gstar25": 4.275, "ea_kc": 79430.0, "ea_ko": 36380.0,
                 "ea_gstar": 37830.0, "ea_vcmax": 65330.0, "ea_jmax": 43540.0, "hd_vcmax": 200000.0,
                 "hd_jmax": 200000.0, "ds_vcmax": 650.0, "ds_jmax": 640.0, "o2_mol_frac": 0.209,
                 "leaf_absorptance": 0.85, "phi_psii": 0.74}


def kattge_knorr(t_growth_c: float) -> dict:
    """{config key: value} of the shape keys at Kattge & Knorr's acclimated values."""
    return {k: a - b * t_growth_c for k, (a, b) in KATTGE_KNORR.items()}


def sat_vapor_pressure(t_k):
    """[Pa] over liquid water: Bolton (1980), as the model's forcing kernels."""
    tc = np.asarray(t_k, dtype=float) - 273.15
    return 611.2 * np.exp(17.67 * tc / (tc + 243.5))


def co2_ppm(base, years) -> float:
    """The run's CO2 [ppm] over these years: forcing.co2_const, or the mean of the CO2 file's
    yearly values."""
    if base.get("forcing.co2_source") != "file":
        return float(base.get("forcing.co2_const"))
    path = Path(base.get("forcing.co2_file"))
    if not path.is_absolute():
        path = (Path(base.path).parent / path).resolve()
    vals = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        m = re.match(r"^(\d{3,4})\s+([0-9.eE+-]+)$", line)
        if m:
            vals[int(m.group(1))] = float(m.group(2))
    got = [vals[y] for y in years if y in vals] or [vals[max(vals)]]
    return float(np.mean(got))


def growth_climate(forcing: str, grid: int, leaf_on: list, daytime_sw: float = 10.0) -> dict:
    """The growing season's climate from the forcing file (module docstring)."""
    from netCDF4 import Dataset, num2date
    with Dataset(forcing) as ds:
        t = ds["time"]
        when = pd.to_datetime([d.strftime("%Y-%m-%d %H:%M:%S")
                               for d in num2date(t[:], t.units, only_use_cftime_datetimes=False)])
        col = {v: np.asarray(ds[v][:, grid - 1], dtype=float) for v in ("Tair", "RHair", "PSurf", "SWdown")}
    df = pd.DataFrame(col, index=when)
    df = df[df.index.month.isin(leaf_on)]
    day = df[df["SWdown"] > daytime_sw]
    vpd = (1.0 - day["RHair"]) * sat_vapor_pressure(day["Tair"])
    return {"t_growth_c": float(df["Tair"].mean() - 273.15), "t_day_k": float(day["Tair"].mean()),
            "vpd_day_pa": float(vpd.mean()), "p_day_pa": float(day["PSurf"].mean()),
            "par_day": float(PAR_PER_SW * day["SWdown"].mean()), "years": sorted(set(df.index.year)),
            "leaf_on_months": list(leaf_on)}


def _arr(k25, ea, t_k):
    return k25 * math.exp(ea / R_GAS * (1.0 / T25 - 1.0 / t_k))


def viscosity_ratio(t_k: float) -> float:
    """Water's viscosity at t_k relative to 25 degC (Vogel's equation)."""
    eta = lambda t: 2.414e-5 * 10.0 ** (247.8 / (t - 140.0))             # noqa: E731  [Pa s]
    return eta(t_k) / eta(T25)


def kinetics(lp: dict, t_k: float, p_pa: float) -> dict:
    """ED-MORE's Rubisco kinetics as partial pressures [Pa] at t_k and the site's pressure."""
    gstar = _arr(lp["gstar25"], lp["ea_gstar"], t_k) * lp["o2_mol_frac"] / O2_REF_GSTAR * p_pa / P_REF_GSTAR
    kc, ko = _arr(lp["kc25"], lp["ea_kc"], t_k), _arr(lp["ko25"], lp["ea_ko"], t_k)
    o2 = lp["o2_mol_frac"] * p_pa
    return {"gstar": gstar, "kc": kc, "ko": ko, "o2": o2, "k": kc * (1.0 + o2 / ko)}


def least_cost_xi(lp: dict, clim: dict) -> float:
    """xi [Pa^0.5] of the least-cost hypothesis at the daytime climate."""
    k = kinetics(lp, clim["t_day_k"], clim["p_day_pa"])
    return math.sqrt(BETA_LEAST_COST * (k["k"] + k["gstar"]) / (1.6 * viscosity_ratio(clim["t_day_k"])))


def eeo_g1(lp: dict, clim: dict) -> float:
    """Medlyn's g1 [kPa^0.5] at the least-cost optimum: xi in Pa^0.5 / sqrt(1000)."""
    return least_cost_xi(lp, clim) / math.sqrt(1000.0)


def eeo_vcmax25(lp: dict, pft: dict, clim: dict, ca_ppm: float, temp_response: str = "peaked") -> float:
    """The coordination vcmax25 with ED-MORE's own leaf equations (module docstring). `pft` holds the
    PFT's theta_j and jmax_vcmax_ratio; `lp` the leaf-physiology settings (with Kattge & Knorr's ds)."""
    from scipy.optimize import brentq
    from edmore.plant import leaf
    t, p = clim["t_day_k"], clim["p_day_pa"]
    k = kinetics(lp, t, p)
    ca = ca_ppm * 1.0e-6 * p
    xi = least_cost_xi(lp, clim)
    chi = k["gstar"] / ca + (1.0 - k["gstar"] / ca) * xi / (xi + math.sqrt(clim["vpd_day_pa"]))
    to_mf = 1.0e6 / p                                       # partial pressure [Pa] -> mole fraction [umol mol-1]
    ci = chi * ca * to_mf

    def at_t(k25, which):
        if temp_response == "peaked":
            return leaf.peaked(k25, lp[f"ea_{which}"], lp[f"hd_{which}"], lp[f"ds_{which}"], t)
        return leaf.arrhenius(k25, lp[f"ea_{which}"], t)

    def gap(v25):
        vcmax, jmax = at_t(v25, "vcmax"), at_t(pft["jmax_vcmax_ratio"] * v25, "jmax")
        j = leaf.electron_transport_j(clim["par_day"], jmax, absorptance=lp["leaf_absorptance"],
                                      phi_psii=lp["phi_psii"], theta=pft["theta_j"])
        r = leaf.assimilation_demand_c3(ci, vcmax, j, gstar=k["gstar"] * to_mf, kc=k["kc"] * to_mf,
                                        ko=k["ko"] * to_mf, o2=k["o2"] * to_mf,
                                        colimitation=leaf.Colimitation.MINIMUM)
        return r.Ac - r.Aj
    return float(brentq(gap, 1.0, 500.0, xtol=1e-6))


def leaf_settings(value) -> dict:
    """ED-MORE's leaf-physiology settings through `value(key, default)` (the base config, else the
    model's parameter record, else the model's default)."""
    return {k: float(value(f"leaf_physiology.{k}", v)) for k, v in LEAF_DEFAULTS.items()}
