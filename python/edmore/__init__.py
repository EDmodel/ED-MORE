# SPDX-License-Identifier: Apache-2.0
"""ED-MORE's Python interface.

Four sub-modules backed by ONE compiled library (``libedmore.so``, see ``edmore._libedmore``), and
two (``edmore.config``, ``edmore.forcing``) that need no library:

    edmore.plant       — plant ecophysiology kernels: leaf gas exchange (``edmore.plant.leaf``) and
                         the leaf-phenology signal kernel (``edmore.plant.pheno``). Stateless.
    edmore.demography  — the demographic carbon slow loop alone, behind an opaque site handle:
                         no biophysics, no forcing, no output streams.
    edmore.model       — the FULL coupled model, exactly what the ``edmore`` executable runs:
                         fast biophysics + slow demography, live met forcing, netCDF output and
                         both conservation ledgers, with the time loop handed to Python.
    edmore.config      — a run's main TOML and the PFT (plant trait) TOML it names: read them,
                         change keys, write them for a run, and read back the parameter record.
    edmore.canopy      — the canopy's kernels over a frozen stand (the leaf solve over many leaves,
                         the two-stream), parameterized by a run's own configuration.
    edmore.forcing     — the writer of a forcing file ([forcing].format = "default"), which every
                         tool that makes one uses.

``edmore.model`` is the whole model; ``edmore.plant``, ``edmore.demography`` and ``edmore.canopy`` are pieces
of it exposed for their own sake.

Importing `edmore` is cheap and does NOT load any compiled library; each sub-module loads the
shared library lazily on first import, so `import edmore` works even without it built.
"""
__version__ = "0.4.0"
__all__ = ["plant", "demography", "model", "config", "canopy", "forcing"]
