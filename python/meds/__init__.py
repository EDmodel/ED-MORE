# SPDX-License-Identifier: Apache-2.0
"""MEDS — Modular Ecosystem Dynamics Simulator (Python interface).

Four sub-modules backed by ONE compiled library (``libmeds.so``, see ``meds._libmeds``), and one
(``meds.config``) that needs no library:

    meds.plant       — plant ecophysiology kernels: leaf gas exchange (``meds.plant.leaf``) and
                       the leaf-phenology signal kernel (``meds.plant.pheno``). Stateless.
    meds.demography  — the demographic carbon slow loop alone, behind an opaque site handle:
                       no biophysics, no forcing, no output streams.
    meds.model       — the FULL coupled model, exactly what the ``meds_main`` executable runs:
                       fast biophysics + slow demography, live met forcing, netCDF output and
                       both conservation ledgers, with the time loop handed to Python.
    meds.config      — a run's main TOML and the PFT (plant trait) TOML it names: read them,
                       change keys, write them for a run, and read back the parameter record.
    meds.canopy      — the canopy's kernels over a frozen stand (the leaf solve over many leaves,
                       the two-stream), parameterized by a run's own configuration.

``meds.model`` is the whole model; ``meds.plant``, ``meds.demography`` and ``meds.canopy`` are pieces
of it exposed for their own sake.

Importing `meds` is cheap and does NOT load any compiled library; each sub-module loads the
shared library lazily on first import, so `import meds` works even without it built.
"""
__version__ = "0.4.0"
__all__ = ["plant", "demography", "model", "config", "canopy"]
