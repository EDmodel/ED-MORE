# SPDX-License-Identifier: Apache-2.0
"""The state chains (protocol §6.1).

Every window has its own chain: a frozen run from the initial stand (the base config's census or
state) that starts `chain_lead_days` before the window and writes a state at its start. Every trial
of that window restarts from it, so the soil water and temperature a window starts with are the
model's own, not the census's initialization. A chain runs with one set of key values (the start's,
and again the fit's at the refresh); a chain already run with the same values is reused.
"""
from __future__ import annotations

import datetime as dt
import threading
from pathlib import Path

from trials import TrialError, command, digest, log_tail, stamp, threads_for, timeout_for, with_keys
from workers import Task

PREFIX = "c"


def run_chain(window, lead_days: int, base, keys, values, root: Path, workers, runner, overrides,
              timeout_per_day: float, max_threads: int = 1) -> Path:
    """Run one window's chain; return the state file at the window's start."""
    start = window.start - dt.timedelta(days=lead_days)
    tag = digest(window.name, str(start), str(window.start),
                 [(p.name, float(v)) for p, v in zip(keys, values) if p.file != "obs"], overrides or {},
                 base.main, base.pft)
    cdir = root / f"chain-{window.name}-{tag[:12]}"
    state = cdir / "out" / f"{PREFIX}-S-{window.start.strftime('%Y%m%d%H%M%S')}.nc"
    if state.exists():
        return state
    cfg = with_keys(base, keys, values, overrides)
    threads = threads_for(lead_days, max_threads)
    for k, v in {"run.start_time": stamp(start), "run.end_time": stamp(window.start), "run.slow_on": False,
                 "run.n_threads": threads, "init.init_mode": 1,
                 "output.state_write": True, "output.state_file_dir": str(cdir / "out"),
                 "output.state_file_prefix": PREFIX, "output.state_interval": "1000 years",
                 "output.enabled": False}.items():
        cfg.set(k, v)
    (cdir / "out").mkdir(parents=True, exist_ok=True)
    main = cfg.write(cdir)
    res = workers.run([Task(f"chain-{window.name}", command(runner, main), str(cdir), str(cdir / "run.log"),
                            timeout_for(lead_days, timeout_per_day), threads)])
    status = next(iter(res.values()))[0]
    if status != "ok" or not state.exists():
        raise TrialError(f"the chain of {window.name} ({start} -> {window.start}) failed: {status}\n"
                         f"{log_tail(cdir / 'run.log')}")
    return state


def run_chains(windows, lead_days: int, base, keys, values, root: Path, workers, runner, overrides,
               timeout_per_day: float, log=print, max_threads: int = 1) -> dict:
    """Every window's state, all chains at once: {window name: state file}."""
    states, errors = {}, []

    def one(w):
        try:
            states[w.name] = run_chain(w, lead_days, base, keys, values, root, workers, runner, overrides,
                                       timeout_per_day, max_threads)
        except Exception as e:           # noqa: BLE001 -- raised below
            errors.append(e)
    log(f"state chains: {len(windows)} windows, each from the initial stand {lead_days} days before it")
    threads = [threading.Thread(target=one, args=(w,)) for w in windows]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    if errors:
        raise errors[0]
    return states
