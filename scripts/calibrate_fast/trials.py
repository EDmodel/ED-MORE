# SPDX-License-Identifier: Apache-2.0
"""Trials: frozen ED-MORE runs over one window, each restarted from that window's state with one set of
key values (protocol §6).

A trial is a directory holding its own main and PFT TOML, made from the base configuration with
edmore.config by setting keys. Its name is a hash of the two, so a repeated set of values reuses the
finished run. A trial runs through the Python API (`python -m edmore.model`) or the edmore
executable; both write the same files.

A trial passes only when:
  - the run ends with "OK: simulation completed";
  - no whole-site budget breached its tolerance;
  - the model's parameter record lists every key the trial set, read from the trial's file, with
    the value written (a misspelt key or one the model does not read fails here, rather than
    silently running the default);
  - its fast output is on the tower's interval, with no missing or non-finite value.
A trial that fails stops the fit with its log: a model that cannot run a set of values inside the
keys' ranges has a bug to fix.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

import targets
from edmore.config import RunConfig, read_record
from workers import Task

#: the fast variables a trial writes: every target's model side (targets.py)
TRIAL_VARIABLES = ("mean_sw_in_py", "mean_sw_up_py", "mean_lw_up_py", "mean_rnet_py", "mean_le_py",
                   "mean_h_py", "mean_gpp_py", "mean_nee_py", "mean_ustar_py")
PREFIX = "t"
#: the runner that runs a config through the Python API instead of an executable
PYTHON_RUNNER = "python"
#: the line edmore and edmore.model print when a run ends well (edmore.model.COMPLETED)
COMPLETED = "OK: simulation completed"
#: the variables that describe the stand in a state file (the stand must not change in a trial)
STAND_VARIABLES = ("pft", "nplant", "dbh", "height", "leaf_area", "leaf_carbon", "fineroot_carbon",
                   "wood_carbon", "overtopping_lai", "patch_area")


class TrialError(RuntimeError):
    pass


@dataclass
class Window:
    name: str
    start: dt.datetime          # UTC, the run's start_time
    days: int
    role: str                   # "cal" (calibration), "val" (validation) or "seasonal" (a seasonal run)

    @property
    def end(self) -> dt.datetime:
        return self.start + dt.timedelta(days=self.days)


def stamp(t: dt.datetime) -> str:
    return t.strftime("%Y-%m-%d %H:%M:%S")


def digest(*objs) -> str:
    """A content hash that ignores key order: the same settings give the same digest however the
    tables were assembled."""
    text = json.dumps(objs, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha1(text.encode()).hexdigest()


#: a run takes one thread for every this many simulated days, up to [fit].max_threads: an iteration
#: waits for its longest runs (the seasonal runs, the chains), while a ten-day window gains little
DAYS_PER_THREAD = 15


def threads_for(days: float, max_threads: int) -> int:
    """The threads of a run of `days` simulated days. ED-MORE's output is the same at any count."""
    return max(1, min(int(max_threads), math.ceil(float(days) / DAYS_PER_THREAD)))


def timeout_for(days: float, per_day: float) -> float:
    """A run's timeout [s]: per_day for each simulated day, at least ten days' worth so a short run's
    start-up (reading the census) fits."""
    return float(per_day) * max(float(days), 10.0)


def with_keys(base: RunConfig, keys, values, overrides: dict | None = None) -> RunConfig:
    """A copy of the base configuration with the keys set to these values, then the calibration's
    overrides. An observation key (kappa) enters the residuals, never a run."""
    cfg = base.copy()
    for p, v in zip(keys, values):
        if p.file != "obs":
            cfg.set(p.key, float(v), file=p.file, pft=p.pft if p.file == "pft" else None)
    for k, v in (overrides or {}).items():
        cfg.set(k, v)
    return cfg


def build_trial(base: RunConfig, keys, values, window: Window, state_file, root: Path,
                overrides: dict | None = None, write_state: bool = False, threads: int = 1) -> Path:
    """Write a trial directory (or find the finished one) and return its path. The thread count is
    set after the name is made, so a finished trial is found whatever the count it ran with."""
    cfg = with_keys(base, keys, values, overrides)
    run = {"run.start_time": stamp(window.start), "run.end_time": stamp(window.end),
           "run.slow_on": False, "run.n_threads": 1,         # in the name as 1; the real count is set below
           "init.init_mode": 2, "init.restart_file": str(state_file), "init.reacclimate_traits": True,
           "output.state_write": write_state, "output.state_file_prefix": PREFIX,
           "output.state_interval": "1000 years", "output.enabled": True, "output.prefix": PREFIX,
           "output.fast.enabled": True, "output.fast.file_chunk": "year",
           "output.daily.enabled": False, "output.monthly.enabled": False,
           "output.annual.enabled": False}
    for k, v in run.items():
        cfg.set(k, v)
    tag = digest(cfg.main, cfg.pft)
    tdir = root / f"{window.name}-{tag[:16]}"
    if (tdir / "series.npz").exists():
        return tdir
    (tdir / "out").mkdir(parents=True, exist_ok=True)
    (tdir / "output_variables.toml").write_text(
        "[variables]\n" + "".join(f'{v} = "F"\n' for v in TRIAL_VARIABLES))
    for k, v in {"output.dir": str(tdir / "out"), "output.state_file_dir": str(tdir / "out"),
                 "output.io_config": str(tdir / "output_variables.toml"), "run.n_threads": int(threads)}.items():
        cfg.set(k, v)
    cfg.remove_diagnostics(tdir)        # an unfinished attempt's files: edmore never overwrites output
    cfg.write(tdir)
    return tdir


def command(runner: str, config: Path) -> list[str]:
    """The command that runs a config: through the Python API when the runner is "python", else
    the edmore executable the runner names."""
    if runner == PYTHON_RUNNER:
        return [sys.executable, "-m", "edmore.model", str(config)]
    return [runner, str(config)]


def log_tail(path: Path, n: int = 15) -> str:
    text = path.read_text(errors="replace") if path.exists() else ""
    return "\n".join(text.splitlines()[-n:])


# ----- the parameter record ---------------------------------------------------------------------
def check_record(tdir: Path, keys, values) -> None:
    rec = read_record(tdir / "out" / f"{PREFIX}_parameters.csv",
                      {tdir / "main.toml": "main", tdir / "pft.toml": "pft"})
    bad = []
    for p, v in zip(keys, values):
        if p.file == "obs":
            continue
        hit = rec.get((p.file, p.key, p.pft if p.file == "pft" else 0))
        if hit is None:
            bad.append(f"{p.name} ({p.key}): not read by the model")
        elif not hit[0]:
            bad.append(f"{p.name} ({p.key}): read as a default, not from the trial's file")
        elif not math.isclose(hit[1], float(v), rel_tol=1e-14, abs_tol=0.0):
            bad.append(f"{p.name} ({p.key}): the model read {hit[1]!r}, the trial wrote {float(v)!r}")
    if bad:
        raise TrialError(f"{tdir.name}: the parameter record does not match the trial:\n  " + "\n  ".join(bad))


# ----- the output --------------------------------------------------------------------------------
def read_series(out_dir: Path, step: float) -> pd.DataFrame:
    """The trial's fast records on the UTC start of each record; they must be `step` seconds apart,
    the tower's interval, or the pairing with the observations would be wrong."""
    frames = []
    for path in sorted(out_dir.glob(f"{PREFIX}-F-*.nc")):
        with Dataset(path) as ds:
            cols = {v: np.asarray(ds[v][:], dtype=float).squeeze() for v in TRIAL_VARIABLES}
            when = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:],
                                       hour=ds["hour"][:], minute=ds["minute"][:]))
        frames.append(pd.DataFrame(cols, index=when))
    if not frames:
        raise TrialError(f"no fast output in {out_dir}")
    df = pd.concat(frames).sort_index()
    gaps = np.unique(np.diff(df.index.values).astype("timedelta64[s]").astype(float))
    if len(df) > 1 and not np.allclose(gaps, step):
        raise TrialError(f"{out_dir}: the fast output is {gaps} s apart, not the tower's {step:g} s")
    return df.where(df.abs() < 1e30)


def finish(tdir: Path, keys, values, step: float, keep_netcdf: bool = False) -> pd.DataFrame:
    """Check a completed trial and keep its fast series (series.npz); raise TrialError if it
    failed. `step` is the tower's interval."""
    log = (tdir / "run.log").read_text(errors="replace") if (tdir / "run.log").exists() else ""
    if COMPLETED not in log:
        raise TrialError(f"{tdir.name}: the run did not complete\n{log_tail(tdir / 'run.log')}")
    #----- a set of values that breaks conservation is not a good run, however well it fits
    for which, fails in re.findall(r"budget\[(whole_\w+)\].*fails = (\d+)/", log):
        if int(fails) > 0:
            raise TrialError(f"{tdir.name}: the {which} budget breached its tolerance {fails} times")
    check_record(tdir, keys, values)
    df = read_series(tdir / "out", step)
    if df.isna().any().any():
        raise TrialError(f"{tdir.name}: missing or non-finite values in the fast output")
    np.savez(tdir / "series.npz", index=df.index.values.astype("datetime64[s]").astype(np.int64),
             **{v: df[v].values for v in TRIAL_VARIABLES})
    if not keep_netcdf:
        for f in (tdir / "out").glob(f"{PREFIX}-F-*.nc"):
            f.unlink()
    return df


def load_series(tdir: Path) -> pd.DataFrame:
    z = np.load(tdir / "series.npz")
    idx = pd.to_datetime(z["index"], unit="s")
    return pd.DataFrame({v: z[v] for v in TRIAL_VARIABLES}, index=idx)


def stand_unchanged(state_a, state_b) -> list:
    """The stand variables that differ between two state files (none: the trial left the stand as
    it found it). Each is read on the site's one polygon, so a chain state written before state files
    had the polygon axis compares with a trial's."""
    def site(var):
        x = np.asarray(var[:])
        return x[0] if var.dimensions[:1] == ("polygon",) else x
    with Dataset(state_a) as a, Dataset(state_b) as b:
        return [v for v in STAND_VARIABLES if v in a.variables
                and not np.array_equal(site(a[v]), site(b[v]))]


# ----- the runner -------------------------------------------------------------------------------
@dataclass
class TrialRunner:
    """Runs sets of key values over windows and returns their residuals."""
    keys: list                       # every key the fit may move (parameters.Param)
    windows: list                    # the windows residuals() runs by default
    rows: dict                       # window name -> targets.WindowRows
    states: dict                     # window name -> its state file
    base: RunConfig                  # the base main and PFT files
    overrides: dict
    runner: str                      # "python" (the Python API) or a edmore executable
    workers: object
    root: Path
    step: float                      # the tower's interval [s]: the trials' output must match it
    timeout_per_day: float
    fixed_observation_keys: dict = field(default_factory=dict)   # observation keys not fitted
    keep_netcdf: bool = False
    log: object = print
    max_threads: int = 1             # a run's threads: threads_for(its days, this)
    seconds: list = field(default_factory=list)
    n_trials: int = 0

    def observation_keys(self, values) -> dict:
        """The observation keys (kappa) of a set of values: they enter the residuals, never a run."""
        out = dict(self.fixed_observation_keys)
        out.update({p.key: float(v) for p, v in zip(self.keys, values) if p.file == "obs"})
        return out

    def run(self, value_sets: list, windows) -> list:
        """Run (or find finished) every (values, window) trial; returns, per set of values, its
        trial directories. A failed trial raises TrialError."""
        threads = {w.name: threads_for(w.days, self.max_threads) for w in windows}
        dirs = [[build_trial(self.base, self.keys, vals, w, self.states[w.name], self.root, self.overrides,
                             threads=threads[w.name]) for w in windows] for vals in value_sets]
        todo, seen = [], set()
        for vals, tds in zip(value_sets, dirs):
            for w, td in zip(windows, tds):
                if td not in seen and not (td / "series.npz").exists():
                    seen.add(td)
                    todo.append((td, vals, Task(td.name, command(self.runner, td / "main.toml"), str(td),
                                                str(td / "run.log"), timeout_for(w.days, self.timeout_per_day),
                                                threads[w.name])))
        status = self.workers.run([t for _, _, t in todo]) if todo else {}
        self.n_trials += len(todo)
        for td, vals, task in todo:
            st, secs = status[task.id]
            self.seconds.append(secs)
            if st != "ok":
                raise TrialError(f"{td.name}: {st}\n{self.describe(vals)}\n{log_tail(td / 'run.log')}")
            try:
                finish(td, self.keys, vals, self.step, self.keep_netcdf)
            except TrialError as e:
                raise TrialError(f"{e}\n{self.describe(vals)}") from None
        return dirs

    def residuals(self, value_sets: list, windows=None) -> list:
        """One stacked residual vector per set of values (over self.keys)."""
        windows = self.windows if windows is None else windows
        out = []
        for vals, tds in zip(value_sets, self.run(value_sets, windows)):
            ok = self.observation_keys(vals)
            parts = []
            for w, td in zip(windows, tds):
                try:
                    parts.append(targets.residual(self.rows[w.name], load_series(td), ok))
                except ValueError as e:
                    raise TrialError(f"{td.name}: {e}\n{self.describe(vals)}") from None
            out.append(np.concatenate(parts) if parts else np.zeros(0))
        return out

    def describe(self, values) -> str:
        return "the values: " + ", ".join(f"{p.name} {float(v):.6g}" for p, v in zip(self.keys, values))
