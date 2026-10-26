# SPDX-License-Identifier: Apache-2.0
"""Smoke tests for meds.plant.pheno (needs libmeds.so built; see python/README.md).

Skips itself cleanly if the shared library hasn't been built, so `pytest` never hard-fails on a
machine that only has the Python sources.
"""
import pytest

import meds.plant.pheno as pheno


def _lib_or_skip():
    try:
        pheno.self_test()
    except FileNotFoundError as exc:
        pytest.skip(f"libmeds.so not built: {exc}")


def _drive(params, days, doy0=1, **env):
    """Run one habit for `days` steps under a constant environment; return the last Day."""
    ph = pheno.Phenology(params)
    day = None
    for i in range(days):
        day = ph.step(doy=(doy0 - 1 + i) % 365 + 1, **env)
    return day


def test_self_test_passes():
    _lib_or_skip()


def test_no_cues_hold_the_canopy_full():
    _lib_or_skip()
    day = _drive(pheno.Params(), 60, temp_day=260.0, par_hours=6.0, predawn_leaf_psi=-5.0)
    assert day.leaf_flush_tendency == 1.0 and day.leaf_shed_tendency == 0.0
    assert day.leaf_cover == 1.0 and day.senescence == 0.0


def test_deciduous_flushes_in_summer_and_senesces_in_autumn():
    _lib_or_skip()
    p = pheno.temperate_deciduous()
    summer = _drive(p, 60, doy0=150, temp_day=295.0, par_hours=15.0)
    assert summer.leaf_flush_tendency > 0.9 and summer.leaf_shed_tendency < 0.05
    autumn = _drive(p, 80, doy0=250, temp_day=275.0, par_hours=8.5)
    assert autumn.leaf_shed_tendency > 0.5 and autumn.leaf_flush_tendency < 0.2
    assert autumn.leaf_cover < 0.01                           # no floor: the canopy goes ~bare


def test_evergreen_is_a_leaf_cover_floor_not_a_flag():
    _lib_or_skip()
    # The deciduous cues with a floor: the same autumn senescence stops at min_leaf_cover
    # (plus the day's small residual flush).
    p = pheno.temperate_deciduous(min_leaf_cover=0.7)
    day = _drive(p, 80, doy0=250, temp_day=275.0, par_hours=8.5)
    assert 0.7 <= day.leaf_cover < 0.701


def test_drought_deciduous_sheds_when_dry_and_reflushes():
    _lib_or_skip()
    # Palo Verde's days: about 12.4 h of light in the wet season, 11.5 h in the dry season.
    p = pheno.drought_deciduous()                             # turgor-loss point -1.355 MPa
    ph = pheno.Phenology(p)
    for _ in range(60):
        wet = ph.step(predawn_leaf_psi=-0.3, par_hours=12.4)
    for _ in range(120):                                      # a four-month dry season
        dry = ph.step(predawn_leaf_psi=-3.5, par_hours=11.5)
    assert wet.leaf_shed_tendency < 0.05 and wet.leaf_cover > 0.99
    assert dry.leaf_shed_tendency > 0.9 and abs(dry.leaf_cover - p.min_leaf_cover) < 1e-3
    for _ in range(60):
        rewet = ph.step(predawn_leaf_psi=-0.3, par_hours=12.4)
    assert rewet.leaf_flush_tendency > 0.2 and rewet.leaf_cover > 0.99


def test_light_exchanging_turns_over_while_staying_full():
    _lib_or_skip()
    p = pheno.light_exchanging()
    dim = _drive(p, 60, par_hours=2.0)                        # few hours above par_min
    bright = _drive(p, 60, par_hours=10.0)
    assert dim.leaf_shed_tendency < 0.05 < bright.leaf_shed_tendency
    assert bright.leaf_flush_tendency > 0.99
    assert bright.senescence > 0.0 and bright.leaf_cover > 0.95   # exchanging, not thinning


def test_leaf_step_matches_the_carbon_rule():
    _lib_or_skip()
    p = pheno.Params(flush_rate_max=0.06, shed_rate_max=0.1, leaf_turnover_rate=0.5)
    cover, sen, bg = pheno.leaf_step(0.5, 0.5, 1.0, p)
    assert abs(sen - 0.05) < 1e-15                            # shed_rate_max * tendency * cover
    assert bg == 0.0                                          # no background turnover while leaves drop
    assert abs(cover - (0.5 - sen + 0.03)) < 1e-15            # + flush_rate_max * tendency
    # With no drop, the leaves age: turnover/yr * flush tendency * cover.
    cover, sen, bg = pheno.leaf_step(0.5, 0.5, 0.0, p)
    assert sen == 0.0 and abs(bg - 0.5 / 365.2425 * 0.5 * 0.5) < 1e-15
    # A bare canopy loses nothing however high the shed tendency.
    assert pheno.leaf_step(0.0, 0.0, 1.0, p) == (0.0, 0.0, 0.0)


def test_par_hours_counts_the_steps_above_par_min():
    import numpy as np
    par = np.array([[0.0, 3.0, 10.0, 400.0, 900.0, 4.0], [0.0] * 6])
    assert list(pheno.par_hours(par, 4.0, 5.0)) == [12.0, 0.0]     # three 4-h steps above 5
    assert list(pheno.par_hours(par, 4.0, 500.0)) == [4.0, 0.0]


def test_a_new_cohort_starts_its_light_memory_from_its_first_day():
    _lib_or_skip()
    ph = pheno.Phenology(pheno.temperate_deciduous())
    assert ph.state.light_hours_mean < 0.0
    ph.step(temp_day=290.0, par_hours=13.0, doy=150)
    assert ph.state.light_hours_mean == 13.0


def test_daylength_is_the_model_formula():
    _lib_or_skip()
    assert abs(pheno.daylength(0.0, 172) - 12.0) < 0.5
    assert pheno.daylength(61.85, 172) > 18.0                 # Hyytiala midsummer
    assert pheno.daylength(80.0, 355) == 0.0                  # polar night


def test_unknown_preset_override_is_rejected():
    with pytest.raises(TypeError):
        pheno.temperate_deciduous(k_flush_max=0.1)            # not a Params field
