# SPDX-License-Identifier: Apache-2.0
"""Census-trained vital rates for the demography engine.

fit_vital_rates.py fits three laws per PFT to the BCI censuses and writes their coefficients to
vital_rates.json. Each step this module reads the stand from the engine and returns the three arrays
``Site.apply_rates`` takes:

  growth       [cm/yr]        [g_min + (g_max - g_min) / (1 + exp(-k (lnD - lnD0)))] exp(-b L)
  mortality    [1/yr]         gamma + alpha exp(-beta growth)          (Camac et al. 2018)
  recruitment  [plants/m2/yr] exp(c0 + c1 LAI), per PFT and patch

D is the cohort's dbh and L its overtopping LAI, both the engine's own: the census described every tree
by the same allometry and the same index (the leaf area of taller trees). LAI is the patch's leaf
area index, all PFTs' together; beyond the census's range it takes its edge.

The census counts every death, including the canopy trees that fall and open a gap. The engine kills
those through treefall disturbance (canopy cohorts on the disturbed area die), so that rate is taken
off the mortality of cohorts tall enough to die in a gap.
"""
import json
import os

try:
    import tomllib                      # py3.11+
except ModuleNotFoundError:             # pragma: no cover -- py3.10 and older
    import tomli as tomllib

import numpy as np

COEFFICIENTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vital_rates.json")
FIELDS = ("dbh", "height", "nplant", "leaf_area", "overtopping_lai", "pft", "owner_patch")


class CensusLaws:
    def __init__(self, main_config):
        with open(main_config, "rb") as fh:
            demography = tomllib.load(fh)["demography"]
        self.treefall_rate = demography["patch_disturbance_rate"]
        self.treefall_height = demography["disturbance_survive_height"]
        with open(COEFFICIENTS) as fh:
            laws = json.load(fh)
        self.g = np.array([laws["growth"][k] for k in ("g_min", "g_max", "D0", "k", "b")]).T   # [pft, 5]
        self.camac = np.array([laws["mortality"][k] for k in ("gamma", "alpha", "beta")]).T
        self.c = np.array(laws["recruitment"]["coefficients"])                 # [pft, 2]
        self.lai_range = laws["range"]["lai"]
        self.n_pft = len(self.g)

    def growth(self, pft, dbh, lai_over):
        """[cm/yr] for 0-based ``pft``; the arguments broadcast against each other."""
        t = self.g[pft]
        size = 1.0 / (1.0 + np.exp(-t[..., 3] * (np.log(dbh) - np.log(t[..., 2]))))
        return (t[..., 0] + (t[..., 1] - t[..., 0]) * size) * np.exp(-t[..., 4] * lai_over)

    def mortality(self, pft, growth):
        """[1/yr] at the given growth [cm/yr]."""
        t = self.camac[pft]
        return t[..., 0] + t[..., 1] * np.exp(-t[..., 2] * growth)

    def recruitment(self, pft, lai):
        """[plants/m2/yr] in a patch of leaf area index ``lai``."""
        c = self.c[pft]
        return np.exp(c[..., 0] + c[..., 1] * np.clip(lai, *self.lai_range))

    def rates(self, site):
        """(growth[n], mortality[n], recruitment[n_pft, n_patch]) for the stand in ``site``."""
        s = {k: site.get(k) for k in FIELDS}
        pft = s["pft"].astype(int) - 1
        patch = s["owner_patch"].astype(int) - 1

        lai = np.zeros(max(site.n_patch, 1))
        np.add.at(lai, patch, s["nplant"] * s["leaf_area"])
        recruitment = self.recruitment(np.arange(self.n_pft)[:, None], lai[None, :])

        growth = self.growth(pft, s["dbh"], s["overtopping_lai"])
        mortality = self.mortality(pft, growth)
        tall = s["height"] >= self.treefall_height
        mortality[tall] = np.maximum(mortality[tall] - self.treefall_rate, 0.0)
        return growth, mortality, recruitment
