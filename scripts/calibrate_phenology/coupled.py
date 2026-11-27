# SPDX-License-Identifier: Apache-2.0
"""The coupled run and the drivers it gives the offline fit.

A run is the site's own configuration with the phenology values of one iteration, its stand frozen
([demography].demography_on = false): every cohort keeps its size, its place and its row in the
output, so each one is a column of the offline stand. The daily cohort output gives each cohort the
cues it read (the canopy-air temperature around it, the hours of light at its top, its predawn leaf
water potential), its leaf turnover rate, and its leaf area.

A run from the census starts every cohort flushing in full leaf, as the offline stand starts by default.
A run from a state file ([init].init_mode = 2) starts each cohort from that file's phenology state, its
tendencies and its cue memory, which the offline stand then starts its column from, with the leaf cover
of the run's first day (a start in the dormant season changes no leaf on that day).
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from edmore.config import RunConfig
from settings import pft_file_key

PREFIX = "pheno"
#----- The daily cohort output a run writes; the layout (cohort_offset, cohort_count) comes with it.
VARIABLES = ("area_pa", "pft_co", "mean_height_co", "mean_lai_co", "mean_pheno_temp_co",
             "mean_pheno_light_hours_co", "mean_dmax_psi_leaf_co", "mean_leaf_turnover_rate_co")
#----- A state file's phenology memory per cohort, in the order of edmore.plant.pheno.State's fields.
STATE_VARIABLES = ("leaf_flush_tendency_co", "leaf_shed_tendency_co", "growing_degree_days_co",
                   "cold_degree_days_co", "wet_psi_sum_co", "dry_psi_sum_co", "light_hours_avg_co")


class RunError(RuntimeError):
    pass


def check_base(cfg: RunConfig) -> None:
    """The site's configuration must freeze the stand and run the slow tier, or the offline
    stand's columns would not be the coupled run's cohorts."""
    if cfg.get("demography.demography_on", default=True) is not False:
        raise RunError(f"{cfg.path}: [demography].demography_on must be false: the cohorts are the fit's "
                       "columns, and a growing stand would change them")
    if cfg.get("run.slow_on", default=True) is not True:
        raise RunError(f"{cfg.path}: [run].slow_on must be true: phenology runs on the slow tier")


def pft_values(cfg: RunConfig, keys) -> dict:
    """The values the configuration gives the keys, by label."""
    return {k.label: float(cfg.get(pft_file_key(k.name), file="pft", pft=k.pft)) for k in keys}


def build(base: RunConfig, keys, values: dict, rundir: Path, threads: int) -> Path:
    """Write the run's files into `rundir` and return its main file. The paths it writes are absolute:
    the run starts in its own directory."""
    rundir = Path(rundir).resolve()
    cfg = base.copy()
    for k in keys:
        cfg.set(pft_file_key(k.name), float(values[k.label]), file="pft", pft=k.pft)
    io = rundir / "output_variables.toml"
    rundir.mkdir(parents=True, exist_ok=True)
    io.write_text('[variables]\n"*" = false\n' + "".join(f'{v} = "D"\n' for v in VARIABLES))
    for k, v in {"output.enabled": True, "output.dir": str(rundir / "out"), "output.prefix": PREFIX,
                 "output.io_config": str(io), "output.overwrite": True, "output.state_write": False,
                 "output.daily.enabled": True, "output.daily.file_chunk": "all",
                 "output.fast.enabled": False, "output.monthly.enabled": False,
                 "output.annual.enabled": False, "run.n_threads": int(threads)}.items():
        cfg.set(k, v)
    return cfg.write(rundir)


def config_digest(main: Path) -> str:
    """A digest of a run's configuration: its main file and the PFT file written beside it."""
    return hashlib.sha256(main.read_bytes() + (main.parent / "pft.toml").read_bytes()).hexdigest()


def restart_file(cfg: RunConfig):
    """The state file a run starts from, or None for a run from the census."""
    if int(cfg.get("init.init_mode", default=1)) != 2:
        return None
    return Path(cfg.get("init.restart_file"))


def run(main: Path, runner: str) -> None:
    """Run a configuration; its log goes beside it."""
    cmd = [sys.executable, "-m", "edmore.model", str(main)] if runner == "python" else [runner, str(main)]
    log = main.parent / "run.log"
    with open(log, "w") as f:
        done = subprocess.run(cmd, cwd=main.parent, stdout=f, stderr=subprocess.STDOUT, env=os.environ.copy())
    if done.returncode != 0:
        tail = "".join(log.read_text(errors="replace").splitlines(keepends=True)[-15:])
        raise RunError(f"the coupled run failed ({' '.join(cmd)}):\n{tail}")


@dataclass
class Drivers:
    """A coupled run's cohorts as the offline stand's columns: (n_col,) traits, (n_col, n_day) series."""
    dates: np.ndarray            # datetime64[D] (n_day,)
    pft: np.ndarray              # (n_col,) 1-based
    height: np.ndarray           # [m]
    area: np.ndarray             # [-] the share of the site the cohort's patch covers
    lai_full: np.ndarray         # [m2/m2 patch] the cohort's full-canopy leaf area
    temp: np.ndarray             # [K]
    hours: np.ndarray            # [h/day]
    psi: np.ndarray              # [MPa]
    turnover: np.ndarray         # [1/yr]
    cover: np.ndarray            # [-] the coupled run's leaf / full leaf
    values: dict = field(default_factory=dict)   # the values the run was made with, by label
    state: np.ndarray = None     # (n_col, len(STATE_VARIABLES)) the phenology state a column starts from,
                                 #   a state file's; None: born flushing, as from the census
    cover0: np.ndarray = None    # (n_col,) [-] the leaf cover it starts from with `state`

    @property
    def doy(self) -> np.ndarray:
        d = self.dates.astype("datetime64[D]")
        return (d - d.astype("datetime64[Y]")).astype(int) + 1

    def save(self, path: Path) -> None:
        np.savez(path, **{k: v for k, v in self.__dict__.items() if k != "values" and v is not None},
                 value_labels=np.array(list(self.values)), value_numbers=np.array(list(self.values.values())))

    @classmethod
    def load(cls, path: Path) -> "Drivers":
        with np.load(path) as z:
            d = {k: z[k] for k in z.files if not k.startswith("value_")}
            values = dict(zip(z["value_labels"].tolist(), z["value_numbers"].tolist()))
        return cls(**d, values=values)


def read_state(path: Path, pft: np.ndarray) -> np.ndarray:
    """Each cohort's phenology state in the state file a run started from, (n_col, len(STATE_VARIABLES)).
    Stops unless the file holds the run's cohorts, in its order."""
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        if ds.dimensions["polygon"].size != 1:
            raise RunError(f"{path}: the phenology calibration runs one site (one polygon)")
        file_pft = np.ma.filled(ds["pft_co"][:], 0)[0].astype(int)
        if file_pft.tolist() != pft.tolist():
            raise RunError(f"{path}: its cohorts' PFTs {file_pft.tolist()} are not the run's {pft.tolist()}")
        return np.column_stack([np.ma.filled(ds[v][:], np.nan)[0].astype(float) for v in STATE_VARIABLES])


def read(out_dir: Path, values: dict, restart: Path = None) -> Drivers:
    """The drivers of a finished run, which started from the state file `restart` (None: from the
    census). Stops if the run had more than one site or its cohorts changed: the offline stand needs
    one fixed set of columns."""
    from netCDF4 import Dataset
    files = sorted(out_dir.glob(f"{PREFIX}-D*.nc"))
    if len(files) != 1:
        raise RunError(f"{out_dir}: expected one daily file, found {[f.name for f in files]}")
    with Dataset(files[0]) as ds:
        if ds.dimensions["polygon"].size != 1:
            raise RunError("the phenology calibration runs one site (one polygon)")
        dates = np.array([np.datetime64(f"{int(y):04d}-{int(m):02d}-{int(d):02d}")
                          for y, m, d in zip(ds["year"][:], ds["month"][:], ds["day"][:])])

        def get(name):
            return np.ma.filled(ds[name][:].astype(float), np.nan)[:, 0, :]

        count = np.ma.filled(ds["cohort_count"][:], 0)[:, 0, :].astype(int)
        offset = np.ma.filled(ds["cohort_offset"][:], 0)[:, 0, :].astype(int)
        if (count != count[0]).any() or (offset != offset[0]).any():
            raise RunError("the cohorts changed during the run: is [demography].demography_on false?")
        n = int(count[0].sum())
        patch = np.full(n, -1)
        for p, (o, c) in enumerate(zip(offset[0], count[0])):
            patch[o - 1:o - 1 + c] = p
        area = get("area_pa")[0][patch]
        pft = get("pft_co")[:, :n]
        if (pft != pft[0]).any():
            raise RunError("a cohort's PFT changed during the run")
        series = {v: get(v)[:, :n].T for v in ("mean_height_co", "mean_lai_co", "mean_pheno_temp_co",
                                                 "mean_pheno_light_hours_co", "mean_dmax_psi_leaf_co",
                                                 "mean_leaf_turnover_rate_co")}
    for v, x in series.items():
        if not np.isfinite(x).all():
            raise RunError(f"{files[0].name}: {v} has missing values")
    lai = series["mean_lai_co"]
    #----- The stand is frozen and each cohort reaches its full canopy in the run (the census starts in
    #      full leaf; a run from a state file passes a summer), so its largest leaf area is that canopy.
    lai_full = lai.max(axis=1)
    cover = np.divide(lai, lai_full[:, None], out=np.ones_like(lai), where=lai_full[:, None] > 0)
    state = read_state(restart, pft[0].astype(int)) if restart is not None else None
    return Drivers(dates=dates, pft=pft[0].astype(int), height=series["mean_height_co"][:, 0], area=area,
                   lai_full=lai_full, temp=series["mean_pheno_temp_co"],
                   hours=series["mean_pheno_light_hours_co"], psi=series["mean_dmax_psi_leaf_co"],
                   turnover=series["mean_leaf_turnover_rate_co"], cover=cover, values=dict(values),
                   state=state, cover0=cover[:, 0].copy() if state is not None else None)
