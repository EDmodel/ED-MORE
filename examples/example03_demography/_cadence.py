# SPDX-License-Identifier: Apache-2.0
"""Shared slow-loop calendar cadence for the example drivers.

Yields the ``(is_new_month, is_new_year)`` flags MEDS's stepper folds into the
monthly/annual fuse-fission triggers, from a simple day counter — the same logic
the Fortran ``meds_stepper`` derives from the real calendar, sufficient for a
fixed-length spin-up example.
"""


def slow_steps(dt_years, n_years, nday=365):
    """Yield (istep, is_new_month, is_new_year) for an n_years spin-up."""
    step_days = max(1, round(dt_years * nday))
    nsteps = (n_years * nday) // step_days
    yday = prev_month = 0
    for istep in range(1, nsteps + 1):
        yday += step_days
        new_year = yday > nday
        if new_year:
            yday -= nday
        month = min(12, (yday - 1) * 12 // nday + 1)
        new_month = (month != prev_month) or new_year or istep == 1
        if istep == 1:
            new_year = True
        prev_month = month
        yield istep, new_month, new_year
