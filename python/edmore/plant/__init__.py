# SPDX-License-Identifier: Apache-2.0
"""edmore.plant — a Pythonic front end to the ED-MORE plant-ecophysiology model (libedmore.so).

The plant module goes beyond a single leaf process, so it is split into submodules (one per process
family), all backed by the one shared library, libedmore.so:

    edmore.plant.leaf — leaf gas exchange: FvCB C3 / Collatz C4 photosynthesis + Leuning / Medlyn /
                        Katul stomatal conductance, Arrhenius / peaked temperature response, and the
                        coupled A-gs-Ci solver.
    edmore.plant.pheno  — leaf phenology: temperature, hours-of-light and water cues + traits ->
                        flush and senescence tendencies -> leaf cover, by the coupled model's rule.

Both are exposed with dataclasses + enums so callers never touch ctypes:

    import edmore.plant.leaf as leaf
    flux = leaf.gas_exchange(par=1500.0, leaf_temp=298.15, vpd=1000.0, ca=400.0, params=leaf.c3_params())

    import edmore.plant.pheno as pheno
    ph = pheno.Phenology(pheno.temperate_deciduous())
    day = ph.step(temp_day=290.0, par_hours=13.0, doy=150)     # day.leaf_cover, day.senescence

Importing `edmore.plant` is cheap; the compiled library (see edmore._libedmore) loads lazily on first use.
"""
from __future__ import annotations

from . import leaf, pheno

__all__ = ["leaf", "pheno"]
