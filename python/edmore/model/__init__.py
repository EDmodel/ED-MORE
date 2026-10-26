# SPDX-License-Identifier: Apache-2.0
"""edmore.model — drive the FULL coupled ED-MORE model from Python.

The other sub-packages expose pieces: ``edmore.plant.leaf`` and ``edmore.plant.pheno`` are
stateless kernels, ``edmore.canopy`` runs the canopy's kernels over a frozen stand, and
``edmore.demography`` drives the slow loop alone. This one drives
what the ``edmore`` executable drives — coupled fast biophysics + slow demography,
live met forcing, netCDF output, both conservation ledgers — with the time loop in
Python:

    from edmore.model import Run

    with Run("ed_config_main.toml") as run:
        for step in run:
            if step.is_new_year:
                print(step.date, run.total_agb, run.soil_carbon)

``Run.step`` takes one slow step of the same run (``ed_region``) the executable steps a
month at a time, so a Python-driven run and a ``edmore`` run of the same config produce
byte-identical output. ``run(config)``
runs a config to its end, and ``python -m edmore.model CONFIG`` does that from the shell in
place of ``edmore CONFIG``. ``edmore.config`` reads and changes the config files. Importing
this sub-module loads the shared library; ``import edmore`` alone does not.
"""
from ._run import COMPLETED, Run, StepInfo, run

__all__ = ["COMPLETED", "Run", "StepInfo", "run"]
