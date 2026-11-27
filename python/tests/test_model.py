# SPDX-License-Identifier: Apache-2.0
"""edmore.model runs a config as edmore does.

- A run through the Python API (`python -m edmore.model`) writes the same output files as edmore,
  bit for bit, and the same parameter record.
- Two runs in one Python process each write their own parameter record.
- Patch areas weight the cohorts to the site totals the model reports.
- A census or state file that cannot start the run stops it: opening the run raises.

Each run is two days from the pseudo census in data/census_example.csv on the fast loop's constant
reference climate. EDMORE_MAIN names the executable for the comparison (CTest sets it); without it
that test is skipped.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from edmore.config import RunConfig, read_record
from edmore.model import COMPLETED, run

ROOT = Path(__file__).resolve().parents[2]
EDMORE_MAIN = os.environ.get("EDMORE_MAIN")


def two_days(run_dir: Path, monthly: bool = False, overrides: dict | None = None) -> Path:
    """A two-day census run with hourly and daily output (and monthly, which keeps the slow cohort
    diagnostics), written into run_dir; `overrides` sets further keys."""
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    for key, value in {"run.start_time": "2001-06-01 00:00:00", "run.end_time": "2001-06-03 00:00:00",
                       "run.n_threads": 1, "fast.fast_biophysics_on": True, "forcing.forcing_on": False,
                       "init.init_mode": 1,
                       "init.census_file": str(ROOT / "data/census_example.csv"),
                       "output.state_write": False, "output.enabled": True, "output.dir": str(run_dir / "out"),
                       "output.prefix": "m", "output.fast.enabled": True, "output.daily.enabled": True,
                       "output.monthly.enabled": monthly, "output.annual.enabled": False,
                       **(overrides or {})}.items():
        cfg.set(key, value)
    (run_dir / "out").mkdir(parents=True, exist_ok=True)
    return cfg.write(run_dir)


def outputs(run_dir: Path) -> dict:
    """{file name: {variable: values}} for every netCDF file the run wrote."""
    out = {}
    for path in sorted((run_dir / "out").glob("*.nc")):
        with Dataset(path) as ds:
            out[path.name] = {v: np.array(ds[v][:]) for v in ds.variables}
    return out


def difference_report(a: dict, b: dict, logs: dict | None = None, limit: int = 30) -> str:
    """Where two runs' outputs part: for each differing variable, the first record (along its first
    axis, time) whose values differ and its largest difference, earliest first; then the log lines
    the two runs do not share. A failure on another machine then says where to look."""
    rows = []
    for f in sorted(a):
        for v in a[f]:
            x, y = a[f][v], b[f].get(v)
            if y is None or x.shape != y.shape:
                rows.append((f, -1, v, f"shapes {x.shape} and {None if y is None else y.shape}"))
                continue
            try:
                x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
            except (TypeError, ValueError):
                if not np.array_equal(x, y):
                    rows.append((f, 0, v, "non-numeric values differ"))
                continue
            neq = ~((x == y) | (np.isnan(x) & np.isnan(y)))
            if not neq.any():
                continue
            first = int(np.argwhere(neq)[0][0]) if x.ndim else 0
            d = np.abs(x - y)[neq]
            scale = np.maximum(np.abs(x[neq]), np.abs(y[neq]))
            rel = np.max(d / np.where(scale > 0, scale, 1.0))
            rows.append((f, first, v, f"max |diff| {d.max():.3e}, max relative {rel:.3e}"))
    rows.sort(key=lambda r: (r[1], r[0], r[2]))
    lines = [f"{len(rows)} variables differ; earliest first (file, first differing record, variable):"]
    lines += [f"  {f} [{first}] {v}: {what}" for f, first, v, what in rows[:limit]]
    if logs:
        (na, la), (nb, lb) = logs.items()
        only_a = [s for s in la if s not in set(lb)]
        only_b = [s for s in lb if s not in set(la)]
        lines.append(f"log lines only in {na}: {only_a[:15]}")
        lines.append(f"log lines only in {nb}: {only_b[:15]}")
    return "\n".join(lines)


def record(run_dir: Path) -> dict:
    """The run's parameter record, with the run's own directory written as <run>."""
    rec = read_record(run_dir / "out" / "m_parameters.csv",
                      {run_dir / "main.toml": "main", run_dir / "pft.toml": "pft"})
    return {k: (p, v.replace(str(run_dir), "<run>") if isinstance(v, str) else v)
            for k, (p, v) in rec.items()}


@pytest.mark.skipif(not EDMORE_MAIN or not Path(EDMORE_MAIN).exists(),
                    reason="EDMORE_MAIN is not set to an edmore executable")
def test_python_run_matches_edmore(tmp_path):
    logs = {}
    for name, argv in (("exe", [EDMORE_MAIN]), ("py", [sys.executable, "-m", "edmore.model"])):
        main = two_days(tmp_path / name)
        res = subprocess.run(argv + [str(main)], cwd=tmp_path / name, capture_output=True, text=True,
                             timeout=600)
        assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
        logs[name] = res.stdout
    for log in logs.values():
        assert log.rstrip().splitlines()[-1] == COMPLETED
    a, b = outputs(tmp_path / "exe"), outputs(tmp_path / "py")
    assert a.keys() == b.keys() and a
    differ = [(f, v) for f in a for v in a[f]
              if a[f][v].shape != b[f][v].shape or not np.array_equal(a[f][v], b[f][v], equal_nan=True)]
    runs = {n: [s.replace(str(tmp_path / n), "<run>") for s in log.splitlines()] for n, log in logs.items()}
    assert not differ, "edmore and edmore.model differ:\n" + difference_report(a, b, runs)
    assert record(tmp_path / "exe") == record(tmp_path / "py")


def test_the_difference_report_names_the_first_differing_record():
    same = np.arange(6.0).reshape(3, 2)
    late = same.copy()
    late[2, 1] += 1e-12
    early = same.copy()
    early[0, 0] += 1.0
    a = {"m-D.nc": {"same": same, "late": same, "early": same}}
    b = {"m-D.nc": {"same": same, "late": late, "early": early}}
    report = difference_report(a, b, {"exe": ["x", "only exe"], "py": ["x"]})
    lines = report.splitlines()
    assert lines[0].startswith("2 variables differ")
    assert "[0] early" in lines[1] and "[2] late" in lines[2]       # earliest first
    assert "only exe" in report and " same" not in report


def test_each_run_in_a_process_keeps_its_own_record(tmp_path):
    for name in ("first", "second"):
        run(two_days(tmp_path / name), verbose=False)
    sources = {k[0] for k in read_record(tmp_path / "second" / "out" / "m_parameters.csv")}
    assert not [s for s in sources if str(tmp_path / "first") in s], sources
    assert record(tmp_path / "first").keys() == record(tmp_path / "second").keys()


def test_patch_areas_weight_cohorts_to_the_site_totals(tmp_path):
    from edmore.model import Run
    with Run(two_days(tmp_path), verbose=False) as r:
        r.step()
        p = r.patches("area", "site", "share")
        c = r.cohorts("agb", "nplant", "owner_patch")
        assert len(p["area"]) == r.n_patch
        assert p["area"].sum() == pytest.approx(1.0, abs=1e-12)
        site_agb = (p["share"][c["owner_patch"] - 1] * c["nplant"] * c["agb"]).sum()
        assert site_agb == pytest.approx(r.total_agb, rel=1e-12)
        #----- one site: every patch in it, its area the whole polygon, each share the patch's area
        assert r.n_site == 1 and (p["site"] == 1).all()
        assert r.sites("area")["area"].tolist() == [1.0]
        assert (p["share"] == p["area"]).all()


def test_cohort_gpp_is_the_last_steps(tmp_path):
    """A cohort's gpp is its gross GPP per plant over the last slow step: never negative, and a
    day of the root config's sunshine gives the stand some."""
    from edmore.model import Run
    with Run(two_days(tmp_path), verbose=False) as r:
        r.step()
        area = r.patches("area")["area"]
        c = r.cohorts("gpp", "nplant", "owner_patch")
        assert (c["gpp"] >= 0.0).all()
        assert (area[c["owner_patch"] - 1] * c["nplant"] * c["gpp"]).sum() > 0.0
        npp = r.cohorts("npp", "npp_wood")
        assert np.isnan(npp["npp"]).all() or (np.isfinite(npp["npp"]).all() and np.isfinite(npp["npp_wood"]).all())


def test_cohort_carbon_closes(tmp_path):
    """A cohort's last-step GPP less its maintenance respiration is, per year, its NPP, growth
    respiration and root exudate: the allocator's identity, read through the Python API."""
    from edmore.model import Run
    with Run(two_days(tmp_path, monthly=True), verbose=False) as r:
        r.step()
        c = r.cohorts("gpp", "leaf_resp", "stem_resp", "root_resp", "npp", "npp_leaf", "npp_fineroot",
                      "npp_wood", "npp_storage", "npp_repro", "growth_resp", "exudate", "leaf_carbon")
    dt_yr = 1.0 / 365.2425                                     # a day, the slow step
    net = (c["gpp"] - c["leaf_resp"] - c["stem_resp"] - c["root_resp"]) / dt_yr
    parts = c["npp_leaf"] + c["npp_fineroot"] + c["npp_wood"] + c["npp_storage"] + c["npp_repro"]
    assert np.allclose(parts, c["npp"], rtol=1e-12, atol=1e-15)
    assert np.allclose(net, c["npp"] + c["growth_resp"] + c["exudate"], rtol=1e-9, atol=1e-12)
    assert (c["leaf_carbon"] > 0.0).all()


@pytest.mark.parametrize("mode, key, name", [(1, "init.census_file", "none.csv"),
                                             (2, "init.restart_file", "none.nc")])
def test_a_missing_initial_file_stops_the_run(tmp_path, mode, key, name):
    """A run asked to start from a census or a state file that is not there raises: it never starts
    from bare ground instead."""
    from edmore.model import Run
    main = two_days(tmp_path, overrides={"init.init_mode": mode, key: str(tmp_path / name)})
    with pytest.raises(RuntimeError, match="could not open"):
        Run(main, verbose=False)
