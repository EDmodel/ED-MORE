# SPDX-License-Identifier: Apache-2.0
"""scripts/prepare_input/pft_traits: the fits a PFT file's fitted traits come from, and the shipped values."""
import importlib.util
import math
from pathlib import Path

try:
    import tomllib
except ImportError:                     # Python < 3.11
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[4]
_spec = importlib.util.spec_from_file_location("pft_traits", Path(__file__).resolve().parents[1] / "pft_traits.py")
pft_traits = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pft_traits)


def test_wood_kmax_is_xu_and_held_at_the_fitted_range():
    assert pft_traits.wood_kmax(0.6) == math.exp(2.348 - 2.455 * 0.6 + 0.5 * 0.6186)
    assert pft_traits.wood_kmax(0.8) < pft_traits.wood_kmax(0.4)
    assert pft_traits.wood_kmax(0.2) == pft_traits.wood_kmax(0.35)
    assert pft_traits.wood_kmax(1.1) == pft_traits.wood_kmax(0.95)


def test_camac_hazard_at_the_reference_density_is_its_scale():
    gamma, alpha, beta = pft_traits.mortality(0.6)
    assert (gamma, alpha, beta) == (0.0094, 0.05, 18.72)


def test_slow_turnover_keeps_one_rate_through_the_canopy():
    k_vcmax, k_rd, k_sla, k_turnover = pft_traits.plasticity(40.0, 0.33)
    assert k_vcmax < 0 and k_rd == k_vcmax and k_sla > 0
    assert k_turnover == 0.0


def test_the_reference_pft_file_carries_the_fits():
    """The shipped values came from the model's own arithmetic, which may round exp and log
    differently in the last place."""
    with open(ROOT / "ed_config_pft.toml", "rb") as fh:
        cfg = tomllib.load(fh)
    pft = cfg["pft"]
    rho, vcmax, turnover = pft["wood_density"], pft["vcmax25"], pft["leaf_turnover_rate_toc"]
    close = lambda a, b: math.isclose(a, b, rel_tol=1e-13, abs_tol=1e-15)
    for i, r in enumerate(rho):
        assert close(cfg["hydraulics"]["wood_kmax"][i], pft_traits.wood_kmax(r))
        for name, value in zip(("mort_gamma", "mort_alpha", "mort_beta"), pft_traits.mortality(r)):
            assert close(pft[name][i], value)
        slopes = pft_traits.plasticity(vcmax[i], turnover[i])
        for name, value in zip(("kplastic_vcmax", "kplastic_rd", "kplastic_sla", "kplastic_leaf_turnover"), slopes):
            assert close(pft[name][i], value)
