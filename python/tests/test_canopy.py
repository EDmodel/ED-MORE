# SPDX-License-Identifier: Apache-2.0
"""edmore.canopy and edmore.plant.leaf.gas_exchange_batch (needs libedmore.so; skipped without it): the
canopy's fast pieces on their own, the same ones the model's fast loop uses."""
from pathlib import Path

import numpy as np
import pytest

import edmore.plant.leaf as leaf
from edmore.config import RunConfig

ROOT = Path(__file__).resolve().parents[2]


def _lib_or_skip():
    try:
        leaf.self_test()
    except FileNotFoundError as exc:
        pytest.skip(f"libedmore.so not built: {exc}")


def test_batch_is_the_single_solve_leaf_by_leaf():
    _lib_or_skip()
    p = leaf.c3_params(vcmax25=50.0, jmax25=90.0)
    par = np.array([0.0, 50.0, 400.0, 1500.0])
    tl = np.array([295.0, 298.15, 301.0, 305.0])
    out = leaf.gas_exchange_batch(par, tl, 1200.0, 400.0, p, pressure=99000.0)
    for i in range(len(par)):
        one = leaf.gas_exchange(par=par[i], leaf_temp=tl[i], vpd=1200.0, ca=400.0, params=p, pressure=99000.0)
        assert out["A_gross"][i] == one.A_gross and out["gs"][i] == one.gs and out["ci"][i] == one.ci
        assert out["limitation"][i] == int(one.limitation) and out["converged"][i] == one.converged


def _config(tmp_path, **keys):
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    for k, v in keys.items():
        cfg.set(k.replace("__", "."), v)
    return cfg.write(tmp_path)


def test_plastic_traits_follow_the_switch(tmp_path):
    _lib_or_skip()
    from edmore.canopy import Canopy
    off = Canopy(_config(tmp_path / "off", trait_dynamics__trait_plasticity_on=False))
    on = Canopy(_config(tmp_path / "on", trait_dynamics__trait_plasticity_on=True))
    pft = np.array([1, 1, 1], dtype=np.int32)
    lai_above = np.array([0.0, 2.0, 5.0])
    v_off, _ = off.plastic_traits(pft, lai_above)
    v_on, _ = on.plastic_traits(pft, lai_above)                    # re-opens its own configuration
    assert np.all(v_off == v_off[0]) and v_on[0] == pytest.approx(v_off[0])
    assert v_on[0] > v_on[1] > v_on[2] > 0.0                        # Vcmax falls with the leaf area above
    assert off.n_pft >= 1


def test_a_config_rewritten_at_the_same_path_is_read_again(tmp_path):
    """A calibration loop rewrites its trial's config at one path: each new Canopy reads the file
    as it is now, not the configuration the library last opened from that path."""
    _lib_or_skip()
    from edmore.canopy import Canopy
    pft = np.array([1, 1], dtype=np.int32)
    lai_above = np.array([0.0, 5.0])
    path = _config(tmp_path, trait_dynamics__trait_plasticity_on=False)
    v_off, _ = Canopy(path).plastic_traits(pft, lai_above)
    assert path == _config(tmp_path, trait_dynamics__trait_plasticity_on=True)
    v_on, _ = Canopy(path).plastic_traits(pft, lai_above)
    assert v_off[0] == v_off[1] and v_on[0] > v_on[1]


def test_the_two_stream_conserves_the_beam(tmp_path):
    _lib_or_skip()
    from edmore.canopy import Canopy
    c = Canopy(_config(tmp_path))
    cosz = np.array([0.0, 0.3, 0.9])
    vb, vd = np.array([0.0, 100.0, 300.0]), np.array([0.0, 80.0, 100.0])
    nb, nd = vb * 1.1, vd * 1.1
    r = c.radiation(cosz, vb, vd, nb, nd, pft=[1, 1], height=[25.0, 5.0], lai=[3.0, 1.5], wai=[0.4, 0.1])
    assert r["cohort_leaf_vis"].shape == (3, 2)
    assert np.allclose(r["cohort_leaf_vis"].sum(axis=1), r["leaf_vis"])
    inc = vb + vd
    alb = 0.15                                                       # ground_albedo_vis in ed_config_main.toml
    #----- every VIS photon is reflected to the sky, absorbed by leaves or wood, or absorbed by the ground
    absorbed_ground = (1.0 - alb) * r["ground_vis"]
    assert np.allclose(r["up_vis"] + r["leaf_vis"] + r["wood_vis"] + absorbed_ground, inc, atol=1e-9)
    assert r["leaf_vis"][2] > r["leaf_vis"][1] > 0.0 and r["leaf_vis"][0] == 0.0
    #----- the taller cohort, on top, absorbs more per unit of leaf area than the one below
    per_leaf = r["cohort_leaf_vis"][2] / np.array([3.0, 1.5])
    assert per_leaf[0] > per_leaf[1]


def test_thermal_acclimation_is_refused(tmp_path):
    _lib_or_skip()
    from edmore.canopy import Canopy
    with pytest.raises(ValueError, match="acclimation"):
        Canopy(_config(tmp_path, leaf_physiology__thermal_acclimation=True))
