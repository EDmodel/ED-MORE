#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""A soil horizon's properties for a site's [soil] lists, by the model's own formulas, so that a site's
script prints what ED-MORE derives from the same texture:

  * the fine earth's van Genuchten curve and Ks from its texture (Hodnett & Tomasella 2002; Ottoni et
    al. 2019), with the coarse fragments after Saxton & Rawls (2006): ed_soil_hydr_lib.f90,
    soil_texture_to_hydraulics;
  * the thermal properties from sand, clay and porosity: ed_therm_lib.f90, soil_texture_to_thermal;
  * an organic horizon's thermal properties: the mineral's blended with peat's by the horizon's organic
    fraction, as CLM5 does (Lawrence & Slater 2008; CLM5 Technical Note, "Soil and snow thermal
    properties");
  * NEON's soil megapit (DP1.00096.001): its horizons' depths, texture, chemistry and bulk density.

The sites' scripts import it (sites/harvard/harvard_traits.py, sites/toolik/toolik_traits.py).
"""
import glob
import math
import os

import pandas as pd

GRAV_HEAD = 9.80665e-3                  # [MPa per m of water head]
PARTICLE_DENSITY = 2.65                 # [g/cm3] mineral grains and fragments (Saxton & Rawls 2006)
# The solids' conductivity by their quartz content (Johansen 1975; Peters-Lidard et al. 1998) and their heat
# capacity per sand and clay (CLM5), as ed_therm_lib.f90.
K_QUARTZ = 7.7                          # [W/m/K]  quartz
K_OTHER_MINERALS = 2.0                  # [W/m/K]  the other minerals
CV_SAND = 2.128e6                       # [J/m3/K] sand grains
CV_CLAY = 2.385e6                       # [J/m3/K] clay grains
CLM_GRAIN_DENSITY = 2700.0              # [kg/m3]
# Peat's end of CLM5's blend (Lawrence & Slater 2008).
OM_DENSITY_MAX = 130.0                  # [kg/m3] organic matter per layer volume at which a layer is all peat
OM_PER_C = 1.0 / 0.58                   # [g OM / g C] organic matter from organic carbon
K_SOLID_OM = 0.25                       # [W/m/K]  peat's solids
K_DRY_OM = 0.05                         # [W/m/K]  dry peat
CV_OM = 2.5e6                           # [J/m3/K] peat's solids


def vg_theta(head, theta_sat, theta_res, alpha, n):
    """Van Genuchten water content at a matric head [m, < 0]."""
    return theta_res + (theta_sat - theta_res) * (1.0 + (alpha * abs(head)) ** n) ** (1.0 / n - 1.0)


def fragment_volume(fragment_mass, bulk_density):
    """The coarse fragments' share of the volume from their share of the mass [g/g] and the fine
    earth's bulk density [g/cm3]: the two densities."""
    w = fragment_mass
    return (w / PARTICLE_DENSITY) / (w / PARTICLE_DENSITY + (1.0 - w) / bulk_density)


def horizon_hydraulics(sand, silt, clay, oc, cec, ph, bd, rock_volume=0.0):
    """A horizon's curve and Ks, as soil_texture_to_hydraulics: sand, silt, clay and organic carbon `oc`
    in % of the fine earth, `cec` [cmol/kg], `ph`, `bd` the fine earth's bulk density [g/cm3] and
    `rock_volume` the fragments' share of the volume."""
    alpha = math.exp((-2.294 - 3.526 * silt + 2.440 * oc - 0.076 * cec - 11.331 * ph + 0.019 * silt ** 2) / 100.0) \
        * GRAV_HEAD * 1000.0                                                  # [1/kPa] -> [1/m]
    n = math.exp((62.986 - 0.833 * clay - 0.529 * oc + 0.593 * ph + 0.007 * clay ** 2 - 0.014 * sand * silt) / 100.0)
    theta_sat = (81.799 + 0.099 * clay - 31.42 * bd + 0.018 * cec + 0.451 * ph - 0.0005 * sand * clay) / 100.0
    theta_res = max((22.733 - 0.164 * sand + 0.235 * cec - 0.831 * ph + 0.0018 * clay ** 2
                     + 0.0026 * sand * clay) / 100.0, 0.0)
    drainable = max(theta_sat - vg_theta(-0.033 / GRAV_HEAD, theta_sat, theta_res, alpha, n), 1e-30)
    ksat = 19.31 / 86400.0 * drainable ** 1.948
    ratio = bd / PARTICLE_DENSITY
    rock_mass = rock_volume / (ratio * (1.0 - rock_volume) + rock_volume)
    return dict(theta_sat=theta_sat * (1 - rock_volume), theta_res=theta_res * (1 - rock_volume), alpha=alpha, n=n,
                ksat=ksat * (1 - rock_mass) / (1 - rock_mass * (1 - 1.5 * ratio)), rock_volume=rock_volume)


def texture_thermal(sand, clay, theta_sat):
    """The solids' conductivity, the dry layer's and the dry layer's heat capacity, as
    soil_texture_to_thermal: sand and clay in % of the fine earth, the layer's porosity. The solids conduct
    by their quartz, the sand: silt and clay count with the other minerals."""
    quartz = sand / 100.0
    k_solid = K_QUARTZ ** quartz * K_OTHER_MINERALS ** (1.0 - quartz)
    bulk_density = CLM_GRAIN_DENSITY * (1.0 - theta_sat)
    k_dry = (0.135 * bulk_density + 64.7) / (CLM_GRAIN_DENSITY - 0.947 * bulk_density)
    dry_cvol = (CV_SAND * sand + CV_CLAY * clay) / (sand + clay) * (1.0 - theta_sat)
    return dict(solid_conductivity=k_solid, dry_conductivity=k_dry, dry_heat_capacity=dry_cvol)


def organic_fraction(bulk_density, carbon):
    """CLM5's organic fraction: the organic matter per layer volume over peat's, at most 1, from the
    bulk density [g/cm3] and the organic carbon [% of the dry mass]."""
    return min(bulk_density * 1000.0 * carbon / 100.0 * OM_PER_C / OM_DENSITY_MAX, 1.0)


def organic_thermal(f_om, sand, clay, mineral_theta_sat, theta_sat):
    """An organic horizon's thermal properties as CLM5 blends them: the mineral's (from the texture and
    the mineral porosity `mineral_theta_sat`) and peat's, weighted by the organic fraction; the solids'
    heat capacity per layer volume at the horizon's own porosity."""
    mineral = texture_thermal(sand, clay, mineral_theta_sat)
    cv_mineral = (CV_SAND * sand + CV_CLAY * clay) / (sand + clay)
    return dict(solid_conductivity=(1 - f_om) * mineral["solid_conductivity"] + f_om * K_SOLID_OM,
                dry_conductivity=(1 - f_om) * mineral["dry_conductivity"] + f_om * K_DRY_OM,
                dry_heat_capacity=((1 - f_om) * cv_mineral + f_om * CV_OM) * (1.0 - theta_sat))


def read_megapit(folder):
    """NEON's megapit (DP1.00096.001, one month's folder): one row per horizon, by name, with its top
    and bottom [m], sand, silt and clay [%], organic carbon `oc` [%], CEC, pH, the fine earth's bulk
    density `bd` [g/cm3] (NaN where the pit has none) and the fragments' volume share."""
    bio = pd.read_csv(glob.glob(os.path.join(folder, "*mgp_perbiogeosample*"))[0]).set_index("horizonName")
    bulk = pd.read_csv(glob.glob(os.path.join(folder, "*mgp_perbulksample*"))[0]).set_index("horizonName")
    rows = []
    for name, h in bio.iterrows():
        bd = float(bulk["bulkDensExclCoarseFrag"].get(name, float("nan")))
        fragments = (h["coarseFrag2To5"] + h["coarseFrag5To20"]) / 1000.0     # [g/g] of the < 20 mm fraction
        rows.append(dict(horizon=name, top=h["biogeoTopDepth"] / 100.0, bottom=h["biogeoBottomDepth"] / 100.0,
                         sand=h["sandTotal"], silt=h["siltTotal"], clay=h["clayTotal"], oc=h["carbonTot"] / 10.0,
                         cec=h["cecdNh4"], ph=h["phH2o"], bd=bd,
                         rock_volume=fragment_volume(fragments, bd) if bd == bd else float("nan")))
    return pd.DataFrame(rows).set_index("horizon").sort_values("top")
