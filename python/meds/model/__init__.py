# SPDX-License-Identifier: Apache-2.0
"""meds.model — drive the FULL coupled MEDS model from Python.

The other sub-packages expose pieces: ``meds.plant.leaf`` and ``meds.plant.pheno`` are
stateless kernels, ``meds.canopy`` runs the canopy's kernels over a frozen stand, and
``meds.demography`` drives the slow loop alone. This one drives
what the ``meds_main`` executable drives — coupled fast biophysics + slow demography,
live met forcing, netCDF output, both conservation ledgers — with the time loop in
Python:

    from meds.model import Run

    with Run("meds_config_main.toml") as run:
        for step in run:
            if step.is_new_year:
                print(step.date, run.total_agb, run.soil_carbon)

``Run.step`` takes one slow step of the same run (``meds_region``) the executable steps a
month at a time, so a Python-driven run and a ``meds_main`` run of the same config produce
byte-identical output. ``run(config)``
runs a config to its end, and ``python -m meds.model CONFIG`` does that from the shell in
place of ``meds_main CONFIG``. ``meds.config`` reads and changes the config files. Importing
this sub-module loads the shared library; ``import meds`` alone does not.
"""
from ._run import COMPLETED, Run, StepInfo, run

__all__ = ["COMPLETED", "Run", "StepInfo", "run"]
