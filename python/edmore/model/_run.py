# SPDX-License-Identifier: Apache-2.0
"""High-level wrapper over the ED-MORE full-coupled-model C-API.

``Run`` is an opaque handle to a live simulation -- the same thing the ``edmore``
executable holds -- with the time loop handed back to Python:

    with Run("ed_config_main.toml") as run:
        for step in run:
            if step.is_new_year:
                print(step.date, run.total_agb, run.soil_carbon)

Each step is one slow step of the same run (``ed_region``) the executable steps a month
at a time, so a Python-driven run and a ``edmore`` run of the same config produce
byte-identical netCDF. Nothing is
re-implemented here; this is the executable's loop, in Python.

Array getters COPY out of the Fortran SoA into fresh numpy arrays (they never alias),
and the SoA is reordered by cohort/patch fusion and fission at the monthly cadence, so
a positional snapshot cached across a step is stale -- ``global_id`` is the only key
that survives a step.
"""
import ctypes
from datetime import date

import numpy as np

from ._ffi import lib

#----- Field ids match ed_c_api_demography's, so the two sub-packages read the same way. ----  #
_REAL = {"dbh": 0, "height": 1, "nplant": 2, "agb": 3, "leaf_area": 4,
         "overtopping_lai": 5, "growth_avg": 6, "wood_carbon": 7,
         "gpp": 8,   # [kgC/plant] gross GPP over the last slow step (a run's own; not a demography site's)
         "npp": 9, "npp_wood": 10,   # [kgC/plant/yr] the last slow step's NPP and its wood part; NaN
                                     # unless the output lists a site NPP total (which keeps them)
         #----- [kgC/plant] maintenance respiration over the last slow step, against that GPP (a run's own)
         "leaf_resp": 11, "stem_resp": 12, "root_resp": 13,
         #----- [kgC/plant] the live pools, and [1/yr] the baseline leaf turnover rate
         "leaf_carbon": 14, "fineroot_carbon": 15, "storage_carbon": 16, "leaf_turnover_rate": 17,
         #----- [kgC/plant/yr] the rest of the last slow step's allocation; NaN as npp is
         "npp_leaf": 18, "npp_fineroot": 19, "npp_storage": 20, "npp_repro": 21,
         "growth_resp": 22, "exudate": 23}
_INT = {"pft": 0, "owner_patch": 1, "global_id": 2}
_PATCH_REAL = {"area": 0, "age": 1, "share": 2}
_PATCH_INT = {"site": 0}
_SITE_REAL = {"area": 0}

#----- A step's status codes (the RUN_* parameters of ed_polygon). -----------------------  #
OK, FINISHED, ERR_NAN, ERR_AREA, ERR_BOUNDS = 0, 1, 2, 3, 4


class StepInfo:
    """What one slow step advanced to. Yielded by iterating a Run."""

    __slots__ = ("date", "istep", "iyear", "is_new_month", "is_new_year")

    def __init__(self, d, istep, iyear, is_new_month, is_new_year):
        self.date, self.istep, self.iyear = d, istep, iyear
        self.is_new_month, self.is_new_year = is_new_month, is_new_year

    def __repr__(self):
        return f"<StepInfo {self.date} step={self.istep}>"


class Run:
    """A live ED-MORE simulation driven from Python.

    Parameters
    ----------
    config : str or Path
        The run TOML -- the same file ``edmore`` takes on its command line.
    verbose : bool
        Keep the Fortran progress lines (init/forcing/summary/ledger) on stdout.
        They go to the process's stdout, not through Python's ``sys.stdout``.
    """

    def __init__(self, config, verbose=True):
        b = str(config).encode("utf-8")
        self.handle = lib.ed_run_open(b, len(b), 1 if verbose else 0)
        if self.handle == -1:
            raise RuntimeError("ED-MORE run registry full (max 4 concurrent runs)")
        if self.handle < 0:
            raise RuntimeError(f"ED-MORE could not open a run from {config!r}")
        self._finalized = False

    #----- Lifecycle ---------------------------------------------------------------------#
    def step(self):
        """Advance ONE slow step. Returns a StepInfo, or None once the run is over."""
        if self.is_done:
            return None
        prev = self.date
        status = lib.ed_run_step(self.handle)
        if status == ERR_NAN:
            #----- The Fortran driver RETURNS these rather than `error stop`ping, which is
            #      the only reason a bad state is an exception here instead of a dead
            #      interpreter. See the ed_region module header.
            raise RuntimeError(f"ED-MORE: NaN in the state or the energy/water ledger at {self.date}")
        if status == ERR_BOUNDS:
            #----- A cohort's density, carbon or water went negative or infinite, or a CENTURY
            #      pool went negative or past the divergence ceiling. The Fortran side has
            #      already printed which field and cohort (or pool and patch); for soil carbon
            #      the slow ledger's report names the PHASE that did it.
            raise RuntimeError(
                f"ED-MORE: a cohort or soil-carbon state out of its physical bounds at {self.date} "
                "(see the printed field and cohort, or pool and patch)")
        if status == FINISHED:
            return None
        if status < 0:
            raise RuntimeError(f"ED-MORE: invalid run handle {self.handle}")
        now = self.date
        return StepInfo(now, self.istep, self.iyear,
                        (now.year, now.month) != (prev.year, prev.month),
                        now.year != prev.year)

    def run_to_end(self, callback=None):
        """Step until end_time. ``callback(run, step)`` is called after each step."""
        n = 0
        for info in self:
            n += 1
            if callback is not None:
                callback(self, info)
        return n

    def finalize(self):
        """Terminal checkpoint, conservation reports, close the netCDF streams.

        The state stays readable afterwards -- freeing it is ``close()``.
        """
        if self._finalized:
            return 0
        status = lib.ed_run_finalize(self.handle)
        self._finalized = True
        if status == ERR_AREA:
            raise RuntimeError("ED-MORE: site area not conserved")
        return status

    def close(self):
        """Finalize if needed, then release the site's Fortran allocations."""
        if self.handle is None:
            return
        self.finalize()
        lib.ed_run_free(self.handle)
        self.handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def __iter__(self):
        while True:
            info = self.step()
            if info is None:
                return
            yield info

    #----- Where the calendar is ---------------------------------------------------------#
    @property
    def is_done(self):
        return bool(lib.ed_run_is_done(self.handle))

    @property
    def date(self):
        return date(lib.ed_run_year(self.handle),
                    lib.ed_run_month(self.handle),
                    lib.ed_run_day(self.handle))

    @property
    def istep(self):
        return lib.ed_run_istep(self.handle)

    @property
    def iyear(self):
        return lib.ed_run_iyear(self.handle)

    #----- Site aggregates: the same reducers the executable's summary line prints. -------#
    @property
    def n_patch(self):
        return lib.ed_run_n_patch(self.handle)

    @property
    def n_cohort(self):
        return lib.ed_run_n_cohort(self.handle)

    @property
    def n_site(self):
        return lib.ed_run_n_site(self.handle)

    @property
    def total_agb(self):
        """Above-ground biomass [kgC/m2], area-weighted over patches."""
        return lib.ed_run_total_agb(self.handle)

    @property
    def total_lai(self):
        """Leaf area index [m2/m2]."""
        return lib.ed_run_total_lai(self.handle)

    @property
    def total_nplant(self):
        """Stem density [plants/m2]."""
        return lib.ed_run_total_nplant(self.handle)

    @property
    def total_basal_area(self):
        """Basal area [m2/m2]."""
        return lib.ed_run_total_basal_area(self.handle)

    @property
    def soil_carbon(self):
        """All seven CENTURY pools [kgC/m2], area-weighted.

        Identically zero unless ``[soil_carbon].soil_carbon_on`` -- with the feature
        off, litter is DISCARDED rather than stored, so this stays flat at zero while
        the stand grows. That is worth watching from the loop.
        """
        return lib.ed_run_soil_carbon(self.handle)

    #----- Per-cohort state, copied out ---------------------------------------------------#
    def cohorts(self, *fields):
        """Copy per-cohort fields into numpy arrays: ``run.cohorts("dbh", "pft")``.

        With no arguments, returns every field. Valid names are the keys of
        ``Run.real_fields`` and ``Run.int_fields``.
        """
        if not fields:
            fields = tuple(_REAL) + tuple(_INT)
        n = self.n_cohort
        out = {}
        for f in fields:
            if f in _REAL:
                buf = np.zeros(max(n, 1), dtype=np.float64)
                lib.ed_run_get_real(self.handle, _REAL[f],
                                      buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            elif f in _INT:
                buf = np.zeros(max(n, 1), dtype=np.int32)
                lib.ed_run_get_int(self.handle, _INT[f],
                                     buf.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
            else:
                raise KeyError(f"unknown cohort field {f!r}; "
                               f"have {sorted(set(_REAL) | set(_INT))}")
            out[f] = buf[:n]
        return out

    #----- Per-patch state, copied out ----------------------------------------------------
    def patches(self, *fields):
        """Copy per-patch fields into numpy arrays: ``run.patches("share")``.

        ``area`` is the patch's fraction of its site, ``site`` that site (1-based), ``share`` the
        patch's fraction of the polygon (its area times its site's) and ``age`` its years since
        the disturbance that made it. A cohort's ``owner_patch`` (1-based) indexes these, so a
        polygon total of a cohort quantity x is ``sum(share[owner_patch - 1] * nplant * x)``.
        """
        if not fields:
            fields = tuple(_PATCH_REAL) + tuple(_PATCH_INT)
        n = self.n_patch
        out = {}
        for f in fields:
            if f in _PATCH_REAL:
                buf = np.zeros(max(n, 1), dtype=np.float64)
                lib.ed_run_get_patch_real(self.handle, _PATCH_REAL[f],
                                            buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            elif f in _PATCH_INT:
                buf = np.zeros(max(n, 1), dtype=np.int32)
                lib.ed_run_get_patch_int(self.handle, _PATCH_INT[f],
                                           buf.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
            else:
                raise KeyError(f"unknown patch field {f!r}; "
                               f"have {sorted(set(_PATCH_REAL) | set(_PATCH_INT))}")
            out[f] = buf[:n]
        return out

    #----- Per-site state, copied out -----------------------------------------------------
    def sites(self, *fields):
        """Copy per-site fields into numpy arrays: ``run.sites("area")``, each site's fraction of
        the polygon. A patch's ``site`` (1-based) indexes them."""
        if not fields:
            fields = tuple(_SITE_REAL)
        n = self.n_site
        out = {}
        for f in fields:
            if f not in _SITE_REAL:
                raise KeyError(f"unknown site field {f!r}; have {sorted(_SITE_REAL)}")
            buf = np.zeros(max(n, 1), dtype=np.float64)
            lib.ed_run_get_site_real(self.handle, _SITE_REAL[f],
                                       buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            out[f] = buf[:n]
        return out

    real_fields = tuple(_REAL)
    int_fields = tuple(_INT)
    patch_fields = tuple(_PATCH_REAL) + tuple(_PATCH_INT)
    site_fields = tuple(_SITE_REAL)

    def __repr__(self):
        if self.handle is None:
            return "<Run closed>"
        return (f"<Run {self.date} step={self.istep} "
                f"cohorts={self.n_cohort} patches={self.n_patch}>")


#: the line edmore prints when a run ends well; a tool that checks a run's log looks for it
COMPLETED = " OK: simulation completed, area conserved, no NaNs."


def run(config, verbose=True):
    """Run a config from its start to its ``end_time``, as the ``edmore`` executable does: the
    same steps, the same output files and the same closing line on stdout. Returns the number of
    slow steps. A NaN in the state, a cohort or soil-carbon state out of bounds, or a site whose
    area is not conserved raises RuntimeError, where the executable would stop."""
    with Run(config, verbose=verbose) as r:
        n = r.run_to_end()
        r.finalize()
    if verbose:
        print(COMPLETED, flush=True)
    return n
