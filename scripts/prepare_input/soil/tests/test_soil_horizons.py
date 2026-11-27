# SPDX-License-Identifier: Apache-2.0
"""scripts/prepare_input/soil: a horizon's properties by the model's formulas, and CLM5's organic blend."""
import importlib.util
import math
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location("soil_horizons", Path(__file__).resolve().parents[1] / "soil_horizons.py")
sh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sh)


def test_thermal_is_the_model_formula():
    t = sh.texture_thermal(50.0, 20.0, 0.4)
    assert t["solid_conductivity"] == pytest.approx(7.7 ** 0.5 * 2.0 ** 0.5)          # quartz = the sand
    rho = 2700.0 * 0.6
    assert t["dry_conductivity"] == pytest.approx((0.135 * rho + 64.7) / (2700.0 - 0.947 * rho))
    assert t["dry_heat_capacity"] == pytest.approx((2.128e6 * 50 + 2.385e6 * 20) / 70 * 0.6)


def test_organic_blend_ends():
    mineral = sh.texture_thermal(45.0, 18.0, 0.45)
    assert sh.organic_thermal(0.0, 45.0, 18.0, 0.45, 0.45) == pytest.approx(mineral)
    peat = sh.organic_thermal(1.0, 45.0, 18.0, 0.45, 0.9)
    assert peat == pytest.approx(dict(solid_conductivity=0.25, dry_conductivity=0.05,
                                      dry_heat_capacity=2.5e6 * 0.1))


def test_organic_fraction_is_organic_matter_over_peats_and_capped():
    assert sh.organic_fraction(0.13, 58.0) == pytest.approx(1.0)      # 130 kg/m3 of organic matter
    assert sh.organic_fraction(0.065, 58.0) == pytest.approx(0.5)
    assert sh.organic_fraction(0.5, 40.0) == 1.0


def test_hydraulics_is_hodnett_tomasella_with_ottoni_ks():
    r = sh.horizon_hydraulics(40.0, 40.0, 20.0, 1.0, 10.0, 5.5, 1.3)
    assert r["theta_sat"] == pytest.approx((81.799 + 0.099 * 20 - 31.42 * 1.3 + 0.018 * 10 + 0.451 * 5.5
                                            - 0.0005 * 40 * 20) / 100)
    assert r["n"] == pytest.approx(math.exp((62.986 - 0.833 * 20 - 0.529 * 1 + 0.593 * 5.5 + 0.007 * 400
                                             - 0.014 * 1600) / 100))
    drained = r["theta_sat"] - sh.vg_theta(-0.033 / sh.GRAV_HEAD, r["theta_sat"], r["theta_res"], r["alpha"], r["n"])
    assert r["ksat"] == pytest.approx(19.31 / 86400 * drained ** 1.948)


def test_fragments_take_volume_and_slow_the_flow():
    clean = sh.horizon_hydraulics(40.0, 40.0, 20.0, 1.0, 10.0, 5.5, 1.3)
    stony = sh.horizon_hydraulics(40.0, 40.0, 20.0, 1.0, 10.0, 5.5, 1.3, rock_volume=0.2)
    assert stony["theta_sat"] == pytest.approx(0.8 * clean["theta_sat"])
    assert stony["ksat"] < clean["ksat"]
    assert sh.fragment_volume(0.0, 1.3) == 0.0
    assert 0.0 < sh.fragment_volume(0.3, 1.3) < 0.3                     # fragments are denser than the fine earth
