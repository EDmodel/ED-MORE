! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_budget_check -- the conservation ledgers (docs/science/numerical_scheme.md §5b):
! the closure test, a per-store accumulator, the merge of patches into a site and of steps
! into a run, and the end-of-run report. A residual is judged against the gross flux that
! crossed the store's boundary, never against the store, and the signed sum of the
! residuals is kept, so a small one-signed leak shows as a drift. Stores are per unit
! ground area: energy [J/m2], water [kg/m2], CO2 [umol/m2]. It sits in shared so every
! process module and the driver test conservation the same way.
!---~---
module ed_budget_check
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: budget_t, closure_ok, budget_imbalance, budget_accumulate, budget_check_stop
   public :: budget_check, budget_merge, budget_report, last_check_closed
   public :: BUDGET_RTOL_FLUX, BUDGET_ENERGY_RATE_FLOOR, BUDGET_WATER_RATE_FLOOR, BUDGET_CO2_RATE_FLOOR

   !--- Closure tolerances of every ledger (docs/science/numerical_scheme.md §5b).
   real(wp), parameter :: BUDGET_RTOL_FLUX         = 1.0e-6_wp  !--- [-] relative to the gross boundary flux
   real(wp), parameter :: BUDGET_ENERGY_RATE_FLOOR = 1.0e-3_wp  !--- [W/m2] absolute floor, as a rate
   real(wp), parameter :: BUDGET_WATER_RATE_FLOOR  = 1.0e-9_wp  !--- [kg/m2/s] absolute floor (~1e-4 mm/day)
   real(wp), parameter :: BUDGET_CO2_RATE_FLOOR    = 1.0e-6_wp  !--- [umol/m2/s] absolute floor

   !--- One store's ledger over a series of checks (fast steps, days, a run).
   type :: budget_t
      real(wp)    :: influx  = 0.0_wp    !--- inflow of the last check (an amount, from budget_check)
      real(wp)    :: outflux = 0.0_wp    !--- outflow of the last check
      real(wp)    :: resid   = 0.0_wp    !--- [store-unit/m2] residual of the last check
      real(wp)    :: worst   = 0.0_wp    !--- [store-unit/m2] largest |resid| of any check
      real(wp)    :: resid_sum  = 0.0_wp !--- [store-unit/m2] signed sum of resid (> 0: the store gains from nowhere)
      real(wp)    :: abs_sum    = 0.0_wp !--- [store-unit/m2] sum of |resid|
      real(wp)    :: flux_gross = 0.0_wp !--- [store-unit/m2] sum of |influx| + |outflux|
      real(wp)    :: elapsed    = 0.0_wp !--- [s] simulated time the checks span (when callers pass it)
      integer(ik) :: n_check = 0_ik      !--- checks made
      integer(ik) :: n_fail  = 0_ik      !--- checks that breached the tolerance
   end type budget_t

contains

   !---~---
   !   The closure test: |resid| <= rtol |scale| + atol, with scale the gross boundary
   ! flux of the window.
   !---~---
   elemental pure function closure_ok(resid, scale, rtol, atol) result(ok)
      real(wp), intent(in) :: resid, scale, rtol, atol
      logical              :: ok
      ok = abs(resid) <= rtol * abs(scale) + atol
   end function closure_ok

   !---~---
   !   One store's residual: (store1 - store0) - dt (influx - outflux); dt = 1 when the
   ! fluxes are already amounts.
   !---~---
   elemental pure function budget_imbalance(store0, store1, influx, outflux, dt) result(resid)
      real(wp), intent(in) :: store0, store1, influx, outflux, dt
      real(wp)             :: resid
      resid = (store1 - store0) - dt * (influx - outflux)
   end function budget_imbalance

   !---~---
   !   Fold one check into b: the window's fluxes and residual, the worst and summed
   ! residuals, the gross flux and the counters. elapsed [s] lets
   ! budget_report turn the summed residual into a mean leak rate.
   !---~---
   pure subroutine budget_accumulate(b, store0, store1, influx, outflux, dt, scale, rtol, atol, elapsed)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: store0, store1, influx, outflux, dt, scale, rtol, atol
      real(wp), optional, intent(in) :: elapsed
      b%influx  = influx
      b%outflux = outflux
      b%resid   = budget_imbalance(store0, store1, influx, outflux, dt)
      b%worst   = max(b%worst, abs(b%resid))
      b%resid_sum  = b%resid_sum  + b%resid
      b%abs_sum    = b%abs_sum    + abs(b%resid)
      b%flux_gross = b%flux_gross + dt * (abs(influx) + abs(outflux))
      if (present(elapsed)) b%elapsed = b%elapsed + elapsed
      b%n_check = b%n_check + 1_ik
      if (.not. closure_ok(b%resid, scale, rtol, atol)) b%n_fail = b%n_fail + 1_ik
   end subroutine budget_accumulate

   !---~---
   !   With debug set, print the residual and stop the run when it breaches the
   ! tolerance; without, do nothing.
   !---~---
   subroutine budget_check_stop(resid, scale, rtol, atol, label, debug)
      real(wp),         intent(in) :: resid, scale, rtol, atol
      character(len=*), intent(in) :: label
      logical,          intent(in) :: debug
      if (debug .and. .not. closure_ok(resid, scale, rtol, atol)) then
         write(*, '(a,a,a,es13.5,a,es13.5)') 'ed_budget_check: budget did NOT close [', &
            trim(label), '] resid = ', resid, '  tol = ', rtol * abs(scale) + atol
         error stop 'ed_budget_check: budget did not close'
      end if
   end subroutine budget_check_stop

   !---~---
   !   One ledger check, the form every fast-loop ledger uses: accumulate, then stop on
   ! a breach if halt is set. influx and outflux are amounts over dt_window.
   !---~---
   subroutine budget_check(b, store0, store1, influx, outflux, dt_window, rate_floor, label, halt)
      type(budget_t),   intent(inout) :: b
      real(wp),         intent(in)    :: store0, store1, influx, outflux
      real(wp),         intent(in)    :: dt_window    !--- [s] time the amounts span
      real(wp),         intent(in)    :: rate_floor   !--- [store-unit/m2/s] absolute tolerance floor, as a rate
      character(len=*), intent(in)    :: label
      logical,          intent(in)    :: halt         !--- stop the run on a breach (Debug)
      real(wp) :: scale, atol
      scale = abs(influx) + abs(outflux)
      atol  = rate_floor * dt_window
      call budget_accumulate(b, store0, store1, influx, outflux, 1.0_wp, scale, BUDGET_RTOL_FLUX, atol, &
                             elapsed=dt_window)
      call budget_check_stop(b%resid, scale, BUDGET_RTOL_FLUX, atol, label, halt)
   end subroutine budget_check

   !---~---
   !   Whether the last budget_check on b closed, by the tolerance it applied. A caller
   ! can check with halt = .false. and decide later: RK45 does, because a step it
   ! hands to the ESDIRK2 rescue is discarded and must not stop the run.
   !---~---
   pure logical function last_check_closed(b, dt_window, rate_floor) result(ok)
      type(budget_t), intent(in) :: b
      real(wp),       intent(in) :: dt_window, rate_floor
      ok = closure_ok(b%resid, abs(b%influx) + abs(b%outflux), BUDGET_RTOL_FLUX, rate_floor * dt_window)
   end function last_check_closed

   !---~---
   !   Fold c into b with weight w: a patch's area fraction when the patches are
   ! summed, 1 when a step's sum joins the run. Sums and elapsed scale with w, so the
   ! result stays per unit ground area and the patches report their common window
   ! once; counts add, and worst is the larger.
   !---~---
   pure subroutine budget_merge(b, c, w)
      type(budget_t), intent(inout) :: b
      type(budget_t), intent(in)    :: c
      real(wp),       intent(in)    :: w
      b%resid_sum  = b%resid_sum  + w * c%resid_sum
      b%abs_sum    = b%abs_sum    + w * c%abs_sum
      b%flux_gross = b%flux_gross + w * c%flux_gross
      b%elapsed    = b%elapsed    + w * c%elapsed
      b%worst      = max(b%worst, c%worst)
      b%n_check    = b%n_check + c%n_check
      b%n_fail     = b%n_fail  + c%n_fail
   end subroutine budget_merge

   !---~---
   !   One report line per store: the signed summed residual, the mean leak rate it
   ! implies over elapsed, the worst single check and the failed checks.
   !---~---
   subroutine budget_report(b, label, unit_store, unit_rate)
      type(budget_t),   intent(in) :: b
      character(len=*), intent(in) :: label, unit_store, unit_rate
      real(wp) :: rate
      rate = 0.0_wp
      if (b%elapsed > 0.0_wp) rate = b%resid_sum / b%elapsed
      write(*, '(a,a,a,es11.3,a,a,a,es11.3,a,a,a,es10.2,a,a,a,i0,a,i0)')                         &
         ' budget[', trim(label), ']  cumulative resid = ', b%resid_sum, ' ', trim(unit_store),   &
         '  mean leak = ', rate, ' ', trim(unit_rate), '  worst/check = ', b%worst, ' ',          &
         trim(unit_store), '  fails = ', b%n_fail, '/', b%n_check
   end subroutine budget_report

end module ed_budget_check
