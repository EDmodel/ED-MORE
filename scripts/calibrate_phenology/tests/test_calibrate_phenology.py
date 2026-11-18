# SPDX-License-Identifier: Apache-2.0
"""Tests for scripts/calibrate_phenology: the settings against phenology_reference.toml, the fitted
PFT file, each target kind on series made by hand, the coupled run's daily file read into drivers,
and an offline fit on synthetic cohorts that finds the values the observations were made with. No
ED-MORE run; the offline stand needs libedmore.so and skips without it."""
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import calibrate_phenology     # noqa: E402,F401  (puts the source tree's edmore package on the path)
import coupled                 # noqa: E402
import pft_file                # noqa: E402
import settings                # noqa: E402
from targets import Target, crossing_doy  # noqa: E402

YEARS = (2019, 2020, 2021)


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def site_file(tmp_path, fit='keys = ["shed_degree_days_pft1"]', target=None, extra="") -> Path:
    target = target or ('[[target]]\nname = "loss"\nkind = "annual_loss"\npfts = [1]\nsigma = 0.1\n'
                        'value = 1.0\n')
    return write(tmp_path / "phenology.toml",
                 f'[run]\nconfig = "main.toml"\nscore_from = "2020-01-01"\n\n[fit]\n{fit}\n{extra}\n{target}')


#----- settings -------------------------------------------------------------------------------------#
def test_settings_fill_the_defaults_and_resolve_paths(tmp_path):
    s = settings.load(site_file(tmp_path))
    assert s["run"]["config"] == (tmp_path / "main.toml").resolve()
    assert s["run"]["score_from"] == dt.date(2020, 1, 1)
    assert s["fit"]["iterations"] == 4 and s["run"]["runner"] == "python"
    (k,) = s["keys"]
    assert (k.name, k.pft, k.lo, k.hi) == ("shed_degree_days", 1, 5.0, 600.0)
    assert s["targets"][0]["weight"] == 1.0 and s["targets"][0]["max_height"] == 1.0e9


@pytest.mark.parametrize("fit, extra, message", [
    ('keys = ["par_min_pft1"]', "", "hour thresholds absorb it"),
    ('keys = ["flush_cue_mask_pft1"]', "", "leaf habit"),
    ('keys = ["leaf_psi_tlp_pft1"]', "", "not a fittable"),
    ('keys = ["shed_degree_days"]', "", "key_pftN"),
    ('keys = ["shed_degree_days_pft1"]\nwarm_start = true', "", "unknown setting"),
    ('keys = ["shed_degree_days_pft1"]', "[bounds]\nshed_degree_days_pft1 = [50.0, 10.0]", "lo < hi"),
])
def test_settings_stop_on_what_they_cannot_do(tmp_path, fit, extra, message):
    with pytest.raises(settings.SettingsError, match=message):
        settings.load(site_file(tmp_path, fit=fit, extra=extra))


def test_the_water_cue_is_fitted(tmp_path):
    s = settings.load(site_file(tmp_path, fit='keys = ["shed_water_sum_pft2", "flush_water_sharpness_pft2"]'))
    assert [(k.name, k.pft, k.lo, k.hi) for k in s["keys"]] == [("shed_water_sum", 2, 0.1, 300.0),
                                                                  ("flush_water_sharpness", 2, 0.02, 20.0)]


def test_a_water_cue_reads_the_leafs_turgor_loss_point():
    """The turgor-loss point the WATER cue sums against is the model's, from the leaf's pressure-volume
    curve: pi0 eps / (pi0 + eps), one value for every PFT or one per PFT."""
    pheno = lib_or_skip()
    import offline
    from edmore.config import RunConfig
    cfg = RunConfig({"run": {}}, {"pft": {"wood_density": [0.66, 0.5]},
                                  "hydraulics": {"leaf_pi0": [-2.0, -1.2], "leaf_elastic_mod": 12.0},
                                  "phenology": {"flush_cue_mask": [1, 4], "shed_cue_mask": [5, 4]}})
    grass = offline.pft_params(cfg, 2)
    assert grass.flush_cue_mask == pheno.Cue.WATER
    assert grass.leaf_psi_tlp == pytest.approx(-1.2 * 12.0 / (-1.2 + 12.0))
    assert offline.pft_params(cfg, 1).leaf_psi_tlp == pytest.approx(-2.0 * 12.0 / 10.0)


def test_a_bound_for_one_pft_overrides_the_default(tmp_path):
    s = settings.load(site_file(tmp_path, fit='keys = ["shed_degree_days_pft1", "shed_degree_days_pft2"]',
                                extra="[bounds]\nshed_degree_days_pft2 = [20.0, 80.0]"))
    assert [(k.lo, k.hi) for k in s["keys"]] == [(5.0, 600.0), (20.0, 80.0)]


@pytest.mark.parametrize("target, message", [
    ('[[target]]\nname = "x"\nkind = "budburst"\npfts = [1]\nsigma = 1.0\n', "kind must be"),
    ('[[target]]\nname = "x"\nkind = "annual_loss"\npfts = [1]\nsigma = 0.0\nvalue = 1.0\n', "sigma"),
    ('[[target]]\nname = "x"\nkind = "fallen"\npfts = [1]\nsigma = 0.1\n', "'file' is required"),
    ('[[target]]\nname = "x"\nkind = "fallen"\npfts = [1]\nsigma = 0.1\nfile = "a.csv"\nlevel = 0.5\n',
     "unknown setting"),
])
def test_targets_are_checked_by_kind(tmp_path, target, message):
    with pytest.raises(settings.SettingsError, match=message):
        settings.load(site_file(tmp_path, target=target))


#----- the fitted PFT file -------------------------------------------------------------------------#
def test_the_pft_file_keeps_its_comments(tmp_path):
    src = write(tmp_path / "pft.toml",
                "# a PFT file\n[pft]\nleaf_turnover_rate_toc = [0.0001, 0.33]   # [1/yr] two PFTs\n\n"
                "[phenology]\nshed_degree_days = [47.889, 50.0]   # [K day]\nflush_degree_days = [90.0, 100.0]\n")
    keys = settings.load(site_file(tmp_path, fit='keys = ["shed_degree_days_pft2", "leaf_turnover_rate_pft2"]'))["keys"]
    out = pft_file.write(src, tmp_path / "fitted.toml", keys,
                         {"shed_degree_days_pft2": 61.5, "leaf_turnover_rate_pft2": 0.25}, "a note")
    text = out.read_text()
    assert "shed_degree_days = [47.889, 61.5]   # [K day]" in text
    assert "leaf_turnover_rate_toc = [0.0001, 0.25]   # [1/yr] two PFTs" in text
    assert "flush_degree_days = [90.0, 100.0]" in text and "# a note" in text
    assert pft_file.set_element("[phenology]\nx = [1.5, 2.5]\n", "phenology.x", 0, 4.0) == "[phenology]\nx = [4.0, 2.5]\n"
    with pytest.raises(pft_file.PftFileError, match="no one-line array"):
        pft_file.set_element(src.read_text(), "phenology.light_window", 0, 1.0)


#----- the targets on series made by hand ----------------------------------------------------------#
def drivers(cover, loss=None, pft=None, height=None, area=None, lai_full=None):
    n_col, n_day = cover.shape
    dates = np.arange(np.datetime64(f"{YEARS[0]}-01-01"), np.datetime64(f"{YEARS[0]}-01-01") + np.timedelta64(n_day, "D"))
    one = np.ones(n_col)
    z = np.zeros((n_col, n_day))
    return coupled.Drivers(dates=dates, pft=np.asarray(pft if pft is not None else one, int),
                           height=np.asarray(height if height is not None else 20 * one, float),
                           area=np.asarray(area if area is not None else one, float),
                           lai_full=np.asarray(lai_full if lai_full is not None else one, float),
                           temp=z + 280.0, hours=z + 12.0, psi=z, turnover=z + 0.1, cover=cover)


def one_target(tmp_path, kind, d, rows=None, **kw):
    t = {"name": kind, "kind": kind, "pfts": [1], "min_height": 0.0, "max_height": 1e9, "sigma": 1.0,
         "weight": 1.0, **kw}
    if rows is not None:
        t["file"] = tmp_path / f"{kind}.csv"
        pd.DataFrame(rows).to_csv(t["file"], index=False)
    return Target(t, d, dt.date(2020, 1, 1))


def ramp_year():
    """A deciduous year: bare to 1 May, full from 1 June to 30 September, bare from 31 October."""
    days = pd.date_range(f"{YEARS[0]}-01-01", f"{YEARS[-1]}-12-31")
    doy = days.day_of_year.to_numpy()
    up = np.clip((doy - 121) / 31.0, 0.0, 1.0)
    down = np.clip((304 - doy) / 31.0, 0.0, 1.0)
    return days, np.minimum(up, down)


def test_fraction_weights_the_cohorts_by_full_leaf_area(tmp_path):
    days, c = ramp_year()
    cover = np.vstack([c, np.ones_like(c)])
    d = drivers(cover, area=[0.5, 0.5], lai_full=[3.0, 1.0])
    t = one_target(tmp_path, "fraction", d, {"date": ["2020-05-16", "2020-07-01"], "value": [0.0, 0.0]},
                   scale="none", window_days=1)
    i = days.get_loc(pd.Timestamp("2020-05-16"))
    assert np.allclose(t.model(cover, cover * 0), [(3 * c[i] + 1) / 4, 1.0])
    # A different height range sees only the other cohort.
    t2 = one_target(tmp_path, "fraction", drivers(cover, height=[25.0, 5.0]),
                    {"date": ["2020-05-16"], "value": [0.0]}, scale="none", window_days=1, max_height=10.0)
    assert np.allclose(t2.model(cover, cover * 0), [1.0])


def test_fraction_on_minmax_and_an_8_day_window(tmp_path):
    days, c = ramp_year()
    cover = (0.6 + 0.4 * c)[None, :]          # an evergreen-like floor: minmax scales it back to 0-1
    t = one_target(tmp_path, "fraction", drivers(cover), {"date": ["2020-05-01"], "value": [0.0]},
                   scale="year_minmax", window_days=8)
    i = days.get_loc(pd.Timestamp("2020-05-01"))
    assert np.isclose(t.model(cover, cover * 0)[0], c[i:i + 8].mean())


def test_fallen_counts_from_its_day_against_the_years_peak(tmp_path):
    days, c = ramp_year()
    cover = (0.5 * c)[None, :]                 # the peak is half a canopy: fallen is relative to it
    t = one_target(tmp_path, "fallen", drivers(cover), {"date": ["2020-07-01", "2020-10-15", "2020-11-15"],
                                                        "value": [0, 0, 0]}, from_doy=213)
    i = days.get_loc(pd.Timestamp("2020-10-15"))
    assert np.allclose(t.model(cover, cover * 0), [0.0, 1.0 - c[i], 1.0])


def test_date_finds_the_crossing_and_a_missing_one_is_late(tmp_path):
    days, c = ramp_year()
    cover = c[None, :]
    t = one_target(tmp_path, "date", drivers(cover), {"year": [2020, 2021], "doy": [140.0, 140.0]},
                   level=0.5, direction="up", doy_min=60, doy_max=200)
    assert np.allclose(t.model(cover, cover * 0), [121 + 15.5, 121 + 15.5])      # the ramp is by day of year
    down = one_target(tmp_path, "date", drivers(cover), {"year": [2020], "doy": [290.0]},
                      level=0.5, direction="down", doy_min=213, doy_max=366)
    assert np.isclose(down.model(cover, cover * 0)[0], 304 - 15.5)
    assert crossing_doy(pd.Series(1.0, index=days[days.year == 2020]), 0.5, "up", 60, 200) == 200.0


def test_annual_loss_and_litter_fraction(tmp_path):
    days, c = ramp_year()
    cover = c[None, :]
    loss = np.zeros_like(cover)
    autumn = (days.month == 10)
    loss[0, autumn] = 1.0 / 31.0                # one canopy a year, all in October
    d = drivers(cover, lai_full=[2.0])
    t = one_target(tmp_path, "annual_loss", d, value=1.0)
    assert t.years == [2020, 2021] and np.allclose(t.model(cover, loss), [1.0, 1.0])
    lit = one_target(tmp_path, "litter_fraction", d,
                     {"date": ["2020-10-16", "2020-11-30", "2021-10-31"], "value": [0, 0, 0]}, year_start="04-01")
    assert np.allclose(lit.model(cover, loss), [16 / 31, 1.0, 1.0])


def test_no_cohort_or_no_observation_stops(tmp_path):
    days, c = ramp_year()
    d = drivers(c[None, :])
    with pytest.raises(ValueError, match="no cohort"):
        one_target(tmp_path, "annual_loss", d, value=1.0, pfts=[2])
    with pytest.raises(ValueError, match="no observation"):
        one_target(tmp_path, "fallen", d, {"date": ["2019-10-01"], "value": [0.5]}, from_doy=213)


#----- the coupled run's daily file ---------------------------------------------------------------#
def daily_file(path, n_day=4, changed=False):
    from netCDF4 import Dataset
    count = np.array([2, 1])
    with Dataset(path, "w") as ds:
        for name, n in (("time", n_day), ("polygon", 1), ("cohort", 4), ("patch", 2)):
            ds.createDimension(name, n)
        t = pd.date_range("2020-12-30", periods=n_day)
        for name, v in (("year", t.year), ("month", t.month), ("day", t.day)):
            ds.createVariable(name, "i4", ("time",))[:] = np.asarray(v)
        c = ds.createVariable("cohort_count", "i4", ("time", "polygon", "patch"))
        c[:] = np.broadcast_to(count, (n_day, 1, 2))
        if changed:
            c[-1, 0, :] = [1, 1]
        ds.createVariable("cohort_offset", "i4", ("time", "polygon", "patch"))[:] = np.broadcast_to([1, 3], (n_day, 1, 2))
        ds.createVariable("area_pa", "f8", ("time", "polygon", "patch"))[:] = np.broadcast_to([0.7, 0.3], (n_day, 1, 2))
        values = {"pft_co": [1, 2, 1], "mean_height_co": [20.0, 8.0, 3.0], "mean_pheno_temp_co": [281.0, 280.0, 279.0],
                  "mean_pheno_light_hours_co": [12.0, 9.0, 6.0], "mean_dmax_psi_leaf_co": [-0.1, -0.2, -0.3],
                  "mean_leaf_turnover_rate_co": [0.1, 0.2, 0.3]}
        for name, v in values.items():
            x = ds.createVariable(name, "f8", ("time", "polygon", "cohort"), fill_value=-9999.0)
            x[:] = -9999.0
            x[:, 0, :3] = np.broadcast_to(v, (n_day, 3))
        lai = ds.createVariable("mean_lai_co", "f8", ("time", "polygon", "cohort"), fill_value=-9999.0)
        lai[:] = -9999.0
        lai[:, 0, :3] = np.outer(np.linspace(1.0, 0.5, n_day), [2.0, 1.0, 0.5])


def test_the_daily_file_gives_one_column_a_cohort(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    daily_file(out / "pheno-D.nc")
    d = coupled.read(out, {"x_pft1": 1.0})
    assert d.pft.tolist() == [1, 2, 1] and np.allclose(d.area, [0.7, 0.7, 0.3])
    assert np.allclose(d.lai_full, [2.0, 1.0, 0.5]) and np.allclose(d.cover[:, -1], 0.5)
    assert np.allclose(d.temp[:, 0], [281.0, 280.0, 279.0]) and d.doy.tolist() == [365, 366, 1, 2]
    d.save(tmp_path / "d.npz")
    back = coupled.Drivers.load(tmp_path / "d.npz")
    assert back.values == {"x_pft1": 1.0} and np.array_equal(back.cover, d.cover)


def test_a_changed_stand_stops(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    daily_file(out / "pheno-D.nc", changed=True)
    with pytest.raises(coupled.RunError, match="cohorts changed"):
        coupled.read(out, {})


#----- the offline stand and the fit (needs libedmore.so) ------------------------------------------#
def lib_or_skip():
    import edmore.plant.pheno as pheno
    try:
        pheno.self_test()
    except (FileNotFoundError, OSError) as exc:
        pytest.skip(f"libedmore.so not built: {exc}")
    return pheno


def synthetic_stand(pheno):
    """Twelve deciduous cohorts in two patches, cooler and darker lower down, three years."""
    days = pd.date_range(f"{YEARS[0]}-01-01", f"{YEARS[-1]}-12-31")
    doy = days.day_of_year.to_numpy()
    length = np.array([pheno.daylength(42.5, int(x)) for x in doy])
    height = np.tile(np.linspace(2.0, 28.0, 6), 2)
    temp = 283.15 + 12.0 * np.cos(2 * np.pi * (doy - 201) / 365.0)[None, :] - 0.05 * (28.0 - height)[:, None]
    hours = np.maximum(length[None, :] - 0.03 * (28.0 - height)[:, None], 0.0)
    n = len(height)
    d = coupled.Drivers(dates=days.to_numpy().astype("datetime64[D]"), pft=np.ones(n, int), height=height,
                        area=np.repeat([0.6, 0.4], 6), lai_full=0.2 + height / 10.0, temp=temp, hours=hours,
                        psi=np.zeros_like(temp), turnover=np.full_like(temp, 0.01), cover=np.ones_like(temp))
    return d


def test_the_fit_finds_the_values_the_observations_were_made_with(tmp_path):
    pheno = lib_or_skip()
    import fit
    import offline
    d = synthetic_stand(pheno)
    truth = pheno.temperate_deciduous(shed_degree_days=60.0, flush_degree_days=120.0, leaf_turnover_rate=0.01)
    start = {"shed_degree_days_pft1": 150.0, "flush_degree_days_pft1": 60.0}
    stand = offline.Stand(d, {1: truth})
    cover, loss = stand.run({1: truth})
    d.cover = cover                                       # the "coupled run" made with the truth
    assert max(w for w, _ in offline.closure(offline.Stand(d, {1: truth})).values()) < 1e-12

    days = pd.DatetimeIndex(d.dates)
    w = d.area * d.lai_full
    whole = pd.Series(w @ cover / w.sum(), index=days)
    fall_days = days[(days.year >= 2020) & (days.month >= 9) & (days.day % 5 == 0)]
    peak = whole.groupby(whole.index.year).transform("max")
    fallen = (1 - whole / peak)[fall_days]
    pd.DataFrame({"date": fall_days.strftime("%Y-%m-%d"), "value": fallen.values}).to_csv(tmp_path / "f.csv", index=False)
    spring = days[(days.year >= 2020) & (days.month.isin([4, 5, 6]))]
    pd.DataFrame({"date": spring.strftime("%Y-%m-%d"), "value": whole[spring].values}).to_csv(tmp_path / "s.csv",
                                                                                              index=False)
    s = settings.load(site_file(tmp_path, fit='keys = ["shed_degree_days_pft1", "flush_degree_days_pft1"]\n'
                                              'popsize = 8\nmaxiter = 60',
                                target='[[target]]\nname = "fallen"\nkind = "fallen"\npfts = [1]\nsigma = 0.05\n'
                                       'file = "f.csv"\n\n[[target]]\nname = "spring"\nkind = "fraction"\n'
                                       'pfts = [1]\nsigma = 0.05\nfile = "s.csv"\n'))
    path = tmp_path / "drivers.npz"
    d.save(path)
    base = {1: pheno.temperate_deciduous(shed_degree_days=150.0, flush_degree_days=60.0, leaf_turnover_rate=0.01)}
    spec = fit.Spec(str(path), base, s["keys"], s["targets"], s["run"]["score_from"])
    res = fit.fit(spec, start, s, workers=2)               # two workers: each starts fresh and loads the drivers
    assert res["score"] < 0.01 * sum(res["scores_start"].values())
    assert abs(res["values"]["shed_degree_days_pft1"] - 60.0) < 3.0
    assert abs(res["values"]["flush_degree_days_pft1"] - 120.0) < 6.0


def test_an_unsettled_fit_keeps_the_values_its_own_run_scored_best():
    trials = [({"k": 1.0}, 30.0, "run0"), ({"k": 2.0}, 12.0, "run1"), ({"k": 3.0}, 25.0, "run2")]
    assert calibrate_phenology.best_coupled(trials) == ({"k": 2.0}, 12.0, "run1")


def test_a_run_is_written_with_absolute_paths(tmp_path, monkeypatch):
    """The run starts in its own directory, so a relative work directory must not reach its files."""
    from edmore.config import RunConfig
    monkeypatch.chdir(tmp_path)
    base = RunConfig({"run": {}, "output": {}}, {"pft": {"wood_density": [0.5]}, "phenology": {"shed_degree_days": [50.0]}})
    keys = settings.load(site_file(tmp_path))["keys"]
    main = coupled.build(base, keys, {"shed_degree_days_pft1": 60.0}, Path("runs") / "run0", threads=1)
    cfg = RunConfig.load(main)
    for key in ("output.dir", "output.io_config"):
        assert Path(cfg.get(key)).is_absolute() and Path(cfg.get(key)).parent.exists()
    assert cfg.get("phenology.shed_degree_days", file="pft") == [60.0]


def test_a_run_is_reused_only_with_the_same_configuration(tmp_path):
    """The digest a finished run is found by covers the whole configuration, not only the fitted values."""
    from edmore.config import RunConfig
    keys = settings.load(site_file(tmp_path))["keys"]
    base = RunConfig({"run": {}, "output": {}}, {"pft": {"wood_density": [0.5], "lai_b1": [0.17]},
                                                 "phenology": {"shed_degree_days": [50.0]}})
    values = {"shed_degree_days_pft1": 60.0}
    d0 = coupled.config_digest(coupled.build(base, keys, values, tmp_path / "a", threads=1))
    assert coupled.config_digest(coupled.build(base, keys, values, tmp_path / "a", threads=1)) == d0
    leafier = base.copy()
    leafier.set("pft.lai_b1", 0.12, file="pft", pft=1)
    assert coupled.config_digest(coupled.build(leafier, keys, values, tmp_path / "a", threads=1)) != d0
