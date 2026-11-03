# SPDX-License-Identifier: Apache-2.0
"""The whole-column water and energy budgets close from the output alone (docs/science/diagnostics.md
§4, The budgets): on the daily and the sub-daily tier, every record's terms add up to its residual,

    water:  precip + shed_water - et - runoff - veg_surface_overflow - water_store_gain = -resid_water
    energy: rnet + heat_bottom + precip_enth - h - le - outflow_enth - heat_gain_*    = -resid_energy

and the ecosystem's carbon is its vegetation's plus its soil's. A term the files leave out, or one
booked on a different basis from the ledger's, breaks a budget by far more than round-off.

The run starts from the shipped census on a synthetic hourly forcing with an afternoon shower every
other day (test_restart_exact.synthetic_forcing), the canopy's surface water and the slow step on.
EDMORE_MAIN names the executable; CTest sets it.
"""
import subprocess

import numpy as np
from netCDF4 import Dataset

from test_restart_exact import EDMORE_MAIN, forced_config, synthetic_forcing

WATER_IN = ["mean_precip_py", "mean_shed_water_py"]
WATER_OUT = ["mean_et_py", "mean_runoff_py", "mean_veg_surface_overflow_py", "mean_water_store_gain_py"]
ENERGY_IN = ["mean_rnet_py", "mean_heat_bottom_py", "mean_precip_enth_py"]
ENERGY_OUT = ["mean_h_py", "mean_le_py", "mean_outflow_enth_py", "mean_heat_gain_soil_py", "mean_heat_gain_canopy_py",
              "mean_heat_gain_surface_py"]
BUDGETS = [("water", WATER_IN, WATER_OUT, "mean_resid_water_py"),
           ("energy", ENERGY_IN, ENERGY_OUT, "mean_resid_energy_py")]


def run(run_dir):
    forcing = run_dir / "forcing.nc"
    synthetic_forcing(forcing)
    io_config = run_dir / "io_budgets.toml"
    names = [n for _, i, o, r in BUDGETS for n in (*i, *o, r)]
    io_config.write_text("[variables]\n" + "".join(f'{n} = "F D"\n' for n in names))
    cfg = forced_config(run_dir, forcing, "2001-06-01", "2001-06-07")
    for key, value in {"fast.canopy_water_on": True, "output.state_write": False, "output.enabled": True,
                       "output.fast.enabled": True, "output.io_config": str(io_config),
                       "output.dir": str(run_dir / "out"), "output.prefix": "b"}.items():
        cfg.set(key, value)
    res = subprocess.run([EDMORE_MAIN, str(cfg.write(run_dir))], cwd=run_dir, capture_output=True,
                         text=True, timeout=600)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    return run_dir / "out"


def read(paths, name):
    return np.concatenate([np.ma.filled(Dataset(p)[name][:], np.nan).ravel() for p in paths])


def test_the_budgets_close_from_the_files(tmp_path):
    out = run(tmp_path)
    for tier in ("D", "F"):
        paths = sorted(out.glob(f"b-{tier}-*.nc"))
        assert paths, f"no {tier} file in {out}"
        assert read(paths, "mean_precip_py").max() > 0.0, f"{tier}: the forcing's showers"
        for budget, ins, outs, resid in BUDGETS:
            terms = {n: read(paths, n) for n in (*ins, *outs, resid)}
            net = sum(terms[n] for n in ins) - sum(terms[n] for n in outs)
            scale = sum(np.abs(terms[n]) for n in (*ins, *outs))
            gap = np.abs(net + terms[resid])
            assert np.all(gap <= 1.0e-9 * scale + 1.0e-14), \
                f"{tier}: the {budget} budget is off by up to {gap.max():.3e} (scale {scale.max():.3e})"
    paths = sorted(out.glob("b-D-*.nc"))
    eco, veg, soil = (read(paths, n) for n in ("mean_eco_carbon_py", "mean_veg_carbon_py", "mean_soilc_total_py"))
    assert np.allclose(eco, veg + soil, rtol=1e-12, atol=0), "ecosystem carbon = vegetation + soil"
