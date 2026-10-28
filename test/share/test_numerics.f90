! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_numerics -- unit tests for the shared numerical primitives + the conservation checker.
!   1. thomas_solve: recover a known solution of a diagonally-dominant tridiagonal system, with
!      the ACTUAL arrays longer than the active count n (the explicit-shape interface).
!   2. quadratic_smaller_root: -> min(a,b) as theta -> 1; strictly below min for theta < 1.
!   3. adaptive_step_update: grow-capped (fmax), shrink-floored (fmin), identity-ish at err = 1.
!   4. ed_budget_check: closure predicate, imbalance algebra, accumulator fail-count, and the
!      no-op behaviour of the Debug hard-stop when debug = .false.
!   5. ascending_order: the reverse on a non-increasing input, ties included; the sort otherwise,
!      with the same tie rule (the larger index first).
!   6. conserved: the tolerance scales with the terms' absolute sum, not the net total.
!   7. pair_mean: the two-member weighted mean fusion blends with; weights of any size count.
!   8. equal_height_layers: runs of exactly equal heights are one layer, either way up.
!   9. The calendar (ed_time): leap years, day and month stepping, day of year, comparisons,
!      the decimal year, parsing and formatting, and the one solar declination that solar_cosz
!      and daylength both use.
!  10. The sun on a slope (ed_time::slope_incidence_cosz; docs/science/forcing.md §8): level
!      ground sees cos z, an equator-facing slope at a winter noon sees cos(z - slope), a
!      pole-facing one less, and an east-facing one more of the morning sun.
!  11. linear_interpolate: states are linear between the bracketing records.
!  12. cbrt: the real cube root, with the sign of its argument.
!---~---
program test_numerics
   use ed_test_assert, only : check, check_true, test_report
   use ed_kinds,          only : wp, ik
   use ed_numerics,       only : quadratic_smaller_root, adaptive_step_update, ascending_order, conserved, &
                                 pair_mean, equal_height_layers, linear_interpolate, cbrt
   use ed_linear_algebra,   only : thomas_solve
   use ed_time,           only : ed_time_t, is_leap_year, days_in_month, days_in_year, day_of_year,           &
                                 time_advance_days, time_advance_months, time_lt, time_eq, string_to_time,    &
                                 time_to_string, time_to_stamp, time_to_decimal_year, solar_cosz, daylength,  &
                                 solar_declination, solar_cosz_utc, seconds_into_day, slope_incidence_cosz
   use ed_budget_check,   only : budget_t, closure_ok, budget_imbalance, budget_accumulate,    &
                                 budget_check_stop, budget_check, budget_merge, budget_report, &
                                 BUDGET_RTOL_FLUX, last_check_closed
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none

   call test_thomas()
   call test_quadratic()
   call test_cbrt()
   call test_step_update()
   call test_budget()
   call test_ascending_order()
   call test_conserved()
   call test_pair_mean()
   call test_equal_height_layers()
   call test_calendar()
   call test_slope_sun()
   call test_interpolation()

   call test_report('test_numerics')

contains



   !--- 1. Tridiagonal solve of a 4-node system embedded in size-8 arrays.
   subroutine test_thomas()
      integer(ik), parameter :: NMAX = 8_ik, N = 4_ik
      real(wp) :: a(NMAX), b(NMAX), c(NMAX), d(NMAX), x(NMAX), xtrue(N)
      integer(ik) :: k
      a = 0.0_wp
      b = 0.0_wp
      c = 0.0_wp
      d = 0.0_wp
      x = 0.0_wp
      xtrue = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      !--- Symmetric diagonally-dominant matrix: b=4, a=c=-1. d = A*xtrue.
      do k = 1_ik, N
         a(k) = -1.0_wp
         b(k) = 4.0_wp
         c(k) = -1.0_wp
      end do
      d(1) = b(1)*xtrue(1) + c(1)*xtrue(2)
      do k = 2_ik, N-1_ik
         d(k) = a(k)*xtrue(k-1) + b(k)*xtrue(k) + c(k)*xtrue(k+1)
      end do
      d(N) = a(N)*xtrue(N-1) + b(N)*xtrue(N)
      call thomas_solve(a, b, c, d, x, N)          ! actual arrays are longer than n
      do k = 1_ik, N
         call check('thomas x', x(k), xtrue(k), 1.0e-12_wp)
      end do
   end subroutine test_thomas

   !--- 12. The real cube root: exact on cubes, odd in x, zero at zero.
   subroutine test_cbrt()
      call check('cbrt of 27', cbrt(27.0_wp), 3.0_wp, 1.0e-12_wp)
      call check('cbrt of -8 keeps the sign', cbrt(-8.0_wp), -2.0_wp, 1.0e-12_wp)
      call check('cbrt of 0', cbrt(0.0_wp), 0.0_wp, 0.0_wp)
      call check('cbrt inverts the cube', cbrt(0.37_wp ** 3), 0.37_wp, 1.0e-12_wp)
   end subroutine test_cbrt

   !--- 2. Smoothed-min quadratic.
   subroutine test_quadratic()
      real(wp) :: x1, x9
      x1 = quadratic_smaller_root(1.0_wp, 3.0_wp, 5.0_wp)     ! theta=1 -> exact min
      call check('quad theta=1 -> min', x1, 3.0_wp, 1.0e-12_wp)
      x9 = quadratic_smaller_root(0.9_wp, 3.0_wp, 5.0_wp)     ! co-limited: strictly below min, positive
      call check_true('quad theta<1 below min', x9 < 3.0_wp .and. x9 > 0.0_wp, x9)
   end subroutine test_quadratic

   !--- 3. Step-doubling factor clamps.
   subroutine test_step_update()
      real(wp) :: f_ok, f_grow, f_shrink
      f_ok     = adaptive_step_update(1.0_wp,    0.9_wp, 0.25_wp, 4.0_wp)   ! err=1 -> safety
      call check('step err=1 -> safety', f_ok, 0.9_wp, 1.0e-12_wp)
      f_grow   = adaptive_step_update(1.0e-4_wp, 0.9_wp, 0.25_wp, 4.0_wp)   ! tiny err -> capped at fmax
      call check('step grow capped', f_grow, 4.0_wp, 1.0e-12_wp)
      f_shrink = adaptive_step_update(100.0_wp,  0.9_wp, 0.25_wp, 4.0_wp)   ! big err -> floored at fmin
      call check('step shrink floored', f_shrink, 0.25_wp, 1.0e-12_wp)
   end subroutine test_step_update

   !--- 4. Conservation checker.
   subroutine test_budget()
      type(budget_t) :: b
      logical  :: ok_closed, ok_open
      real(wp) :: resid
      !--- Closure predicate (mixed rtol/atol).
      ok_closed = closure_ok(1.0e-10_wp, 100.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      ok_open   = closure_ok(1.0_wp,     100.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      call check_true('closure_ok closed', ok_closed, 1.0_wp)
      call check_true('closure_ok open flagged', .not. ok_open, 0.0_wp)
      !--- Imbalance algebra: a perfectly closed store has zero residual.
      resid = budget_imbalance(10.0_wp, 12.0_wp, 3.0_wp, 1.0_wp, 1.0_wp)     ! (12-10) - (3-1) = 0
      call check('imbalance closed = 0', resid, 0.0_wp, 1.0e-14_wp)
      !--- Accumulator: one closed + one open store -> n_fail = 1, worst > 0.
      call budget_accumulate(b, 10.0_wp, 12.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 12.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      call budget_accumulate(b, 10.0_wp, 15.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 12.0_wp, 1.0e-6_wp, 1.0e-8_wp)
      call check_true('accumulator n_check=2', b%n_check == 2_ik, real(b%n_check, wp))
      call check_true('accumulator n_fail=1', b%n_fail == 1_ik, real(b%n_fail, wp))
      call check('accumulator worst', b%worst, 3.0_wp, 1.0e-14_wp)   ! second store off by 3
      !--- Debug hard-stop is a NO-OP when debug = .false. (must return, not halt).
      call budget_check_stop(1.0e6_wp, 1.0_wp, 1.0e-6_wp, 1.0e-8_wp, 'no-op path', .false.)
      call check_true('budget_check_stop no-op returns', .true., 1.0_wp)
      !---~---
      !   SIGNED accumulation: the two residuals above were 0 and +3, so the cumulative signed
      ! sum is +3 and the absolute sum is 3 -- a running max alone could not tell +3 from -3.
      !---~---
      call check('accumulator resid_sum (signed)', b%resid_sum, 3.0_wp, 1.0e-14_wp)
      call check('accumulator abs_sum', b%abs_sum, 3.0_wp, 1.0e-14_wp)
      call check('accumulator flux_gross', b%flux_gross, 8.0_wp, 1.0e-14_wp)   ! 2 checks x (3+1)
      !---~---
      !   budget_check: FLUX-scaled tolerance. A residual of 1e-5 on a gross flux of 4 breaches
      ! rtol=1e-6 (tol 4e-6 + floor*dt); the same residual on a gross flux of 1e4 does not. A
      ! store-scaled rule would judge both against store1 = 12.
      !---~---
      block
         type(budget_t) :: c, d, r
         call budget_check(c, 10.0_wp, 12.0_wp + 1.0e-5_wp, 3.0_wp, 1.0_wp, 900.0_wp, 1.0e-12_wp, 'small flux', .false.)
         call check_true('budget_check: 1e-5 leak on a gross flux of 4 is flagged', c%n_fail == 1_ik, real(c%n_fail, wp))
         call budget_check(d, 10.0_wp, 12.0_wp + 1.0e-5_wp, 5.0e3_wp + 2.0_wp, 5.0e3_wp, 900.0_wp, 1.0e-12_wp, &
               'large flux', .false.)
         call check_true('budget_check: same leak on a gross flux of 1e4 passes', d%n_fail == 0_ik, real(d%n_fail, wp))
         call check('budget_check records elapsed', c%elapsed, 900.0_wp, 1.0e-12_wp)
         !---~---
         !   last_check_closed re-applies the same tolerance to the recorded amounts, so a caller
         ! can check with halt = .false. and decide later whether the breach counts (RK45).
         !---~---
         call check_true('last_check_closed: the small-flux breach is open',                          &
                         .not. last_check_closed(c, 900.0_wp, 1.0e-12_wp), 0.0_wp)
         call check_true('last_check_closed: the large-flux check is closed',                         &
                         last_check_closed(d, 900.0_wp, 1.0e-12_wp), 1.0_wp)
         !---~---
         !   merge: two "patches" of area 0.25 and 0.75 with residuals +1e-5 each -> site +1e-5;
         ! elapsed is the common window (weights sum to 1), counts add.
         !---~---
         r = budget_t()
         call budget_merge(r, c, 0.25_wp)
         call budget_merge(r, d, 0.75_wp)
         call check('budget_merge: area-weighted resid_sum', r%resid_sum, 1.0e-5_wp, 1.0e-12_wp)
         call check('budget_merge: elapsed is the common window', r%elapsed, 900.0_wp, 1.0e-12_wp)
         call check_true('budget_merge: counts add', r%n_check == 2_ik .and. r%n_fail == 1_ik, real(r%n_check, wp))
         call check_true('BUDGET_RTOL_FLUX is the documented 1e-6', BUDGET_RTOL_FLUX == 1.0e-6_wp, BUDGET_RTOL_FLUX)
         call budget_report(r, 'unit-test', 'J/m2', 'W/m2')
      end block
   end subroutine test_budget

   !---~---
   !   5. ascending_order. A tallest-first canopy with a tie (the model's case), a shuffled one
   ! with the same heights, and the empty and one-cohort canopies.
   !---~---
   subroutine test_ascending_order()
      real(wp)    :: tall_first(5), shuffled(5)
      integer(ik) :: ord(5), ord0(0), ord1(1)
      tall_first = [30.0_wp, 20.0_wp, 20.0_wp, 10.0_wp, 5.0_wp]
      call ascending_order(tall_first, 5_ik, ord)
      call check_true('ascending_order: tallest-first input gives the reverse, ties included',     &
                      all(ord == [5_ik, 4_ik, 3_ik, 2_ik, 1_ik]), real(ord(1), wp))
      shuffled = [20.0_wp, 5.0_wp, 30.0_wp, 20.0_wp, 10.0_wp]
      call ascending_order(shuffled, 5_ik, ord)
      call check_true('ascending_order: shuffled input is sorted, the larger index first on a tie', &
                      all(ord == [2_ik, 5_ik, 4_ik, 1_ik, 3_ik]), real(ord(3), wp))
      call ascending_order(shuffled, 0_ik, ord0)
      call check_true('ascending_order: an empty canopy', size(ord0) == 0, 0.0_wp)
      call ascending_order(shuffled, 1_ik, ord1)
      call check_true('ascending_order: one cohort', ord1(1) == 1_ik, real(ord1(1), wp))
   end subroutine test_ascending_order

   !---~---
   !   8. Cohorts of exactly equal height share a layer, numbered from the first cohort, in
   ! ascending (radiation, wind) or descending (overtopping LAI) order.
   !---~---
   subroutine test_equal_height_layers()
      integer(ik) :: layer(5), layer0(0), n_layer
      call equal_height_layers([2.0_wp, 2.0_wp, 7.5_wp, 20.0_wp, 20.0_wp], 5_ik, layer, n_layer)
      call check_true('equal_height_layers: ascending, two ties', all(layer == [1_ik, 1_ik, 2_ik, 3_ik, 3_ik])  &
                      .and. n_layer == 3_ik, real(n_layer, wp))
      call equal_height_layers([30.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 5.0_wp], 5_ik, layer, n_layer)
      call check_true('equal_height_layers: descending, three tied', all(layer == [1_ik, 2_ik, 2_ik, 2_ik, 3_ik])&
                      .and. n_layer == 3_ik, real(n_layer, wp))
      call equal_height_layers([20.0_wp, 20.0_wp * (1.0_wp + epsilon(1.0_wp))], 2_ik, layer(1:2), n_layer)
      call check_true('equal_height_layers: heights one ulp apart are two layers', n_layer == 2_ik,           &
                      real(n_layer, wp))
      call equal_height_layers([real(wp) ::], 0_ik, layer0, n_layer)
      call check_true('equal_height_layers: an empty canopy has no layer', n_layer == 0_ik, real(n_layer, wp))
   end subroutine test_equal_height_layers

   !--- The two-member weighted mean fusion blends with: weights of any size count.
   subroutine test_pair_mean()
      call check('pair_mean: equal weights', pair_mean(12.0_wp, 1.0_wp, 24.0_wp, 1.0_wp), 18.0_wp, 1.0e-15_wp)
      !---~---
      !   Two leafless cohorts' leaf carbon, ~1e-193, as in a BCI run. A threshold at 1e-30
      ! would keep the first SLA; the mean must blend as it does at any other size.
      !---~---
      call check('pair_mean: weights of 1e-193 still blend',                                         &
                 pair_mean(12.0_wp, 1.0e-193_wp, 24.0_wp, 3.0e-193_wp), 21.0_wp, 1.0e-14_wp)
      call check('pair_mean: no weight at all keeps the first value',                                &
                 pair_mean(12.0_wp, 0.0_wp, 24.0_wp, 0.0_wp), 12.0_wp, 0.0_wp)
      call check('pair_mean: all the weight on the second gives the second exactly',                 &
                 pair_mean(12.0_wp, 0.0_wp, 24.0_wp, 5.0_wp), 24.0_wp, 0.0_wp)
   end subroutine test_pair_mean

   !--- conserved: the tolerance scales with the terms' absolute sum, not the net total.
   subroutine test_conserved()
      real(wp), parameter :: ROUNDOFF_RTOL = 1.0e-12_wp   ! the round-off of summing a stored amount
      real(wp) :: nan_v
      !---~---
      !   A dry canopy whose film total is a tiny NEGATIVE number, reordered. Relative to the
      ! net total floored at 1e-30 this round-off (3.7e-40) would fail.
      !---~---
      call check_true('conserved: round-off of a near-zero mixed-sign total passes',                 &
                      conserved(-2.8665915636071925e-24_wp, -2.8665915636071928e-24_wp,             &
                                2.927102e-24_wp, ROUNDOFF_RTOL), 0.0_wp)
      call check_true('conserved: a real change of 1e-6 of the magnitude fails',                     &
                      .not. conserved(1.0_wp, 1.0_wp + 1.0e-6_wp, 1.0_wp, ROUNDOFF_RTOL), 0.0_wp)
      call check_true('conserved: nothing to conserve and nothing changed passes',                   &
                      conserved(0.0_wp, 0.0_wp, 0.0_wp, ROUNDOFF_RTOL), 0.0_wp)
      call check_true('conserved: something from nothing fails',                                     &
                      .not. conserved(0.0_wp, 1.0e-30_wp, 0.0_wp, ROUNDOFF_RTOL), 0.0_wp)
      !---~---
      !   A BCI fusion: leaf area at the bottom of the double range flushes to zero. The format
      ! resolves nothing below tiny/epsilon, so that is not a conservation failure.
      !---~---
      call check_true('conserved: an amount below what the format resolves may flush to zero',        &
                      conserved(1.4029406990646508e-307_wp, 0.0_wp, 1.4029406990646508e-307_wp,       &
                                ROUNDOFF_RTOL), 0.0_wp)
      nan_v = ieee_value(nan_v, ieee_quiet_nan)
      call check_true('conserved: a NaN fails', .not. conserved(1.0_wp, nan_v, 1.0_wp, ROUNDOFF_RTOL), 0.0_wp)
   end subroutine test_conserved


   !---~---
   !   9. The calendar: leap years (with the century rule), stepping by days and months, day of
   ! year, comparisons, the decimal year, string round trips, and one solar declination.
   !---~---
   subroutine test_calendar()
      type(ed_time_t)   :: a, b
      logical           :: ok
      real(wp)          :: cz_noon, cz_mid, cz_summer, cz_winter
      print '(a)', 'test_calendar:'
      call check_true('2000 is a leap year (div by 400)', is_leap_year(2000_ik))
      call check_true('1900 is not a leap year (div by 100, not 400)', .not. is_leap_year(1900_ik))
      call check_true('2004 is a leap year', is_leap_year(2004_ik))
      call check_true('2001 is not a leap year', .not. is_leap_year(2001_ik))
      call check_true('Feb 2000 has 29 days', days_in_month(2000_ik, 2_ik) == 29_ik)
      call check_true('Feb 1900 has 28 days', days_in_month(1900_ik, 2_ik) == 28_ik)
      call check_true('April has 30 days', days_in_month(2001_ik, 4_ik) == 30_ik)
      call check_true('January has 31 days', days_in_month(2001_ik, 1_ik) == 31_ik)
      call check_true('2000 has 366 days', days_in_year(2000_ik) == 366_ik)
      call check_true('2001 has 365 days', days_in_year(2001_ik) == 365_ik)

      !--- Day stepping across leap and non-leap February and the year boundary.
      a = time_advance_days(ed_time_t(2000_ik, 2_ik, 28_ik), 1_ik)
      call check_true('2000-02-28 + 1 day -> 02-29 (leap)', a%month == 2_ik .and. a%day == 29_ik)
      a = time_advance_days(ed_time_t(1900_ik, 2_ik, 28_ik), 1_ik)
      call check_true('1900-02-28 + 1 day -> 03-01 (non-leap)', a%month == 3_ik .and. a%day == 1_ik)
      a = time_advance_days(ed_time_t(2000_ik, 12_ik, 31_ik), 1_ik)
      call check_true('year roll-over on +1 day', a%year == 2001_ik .and. a%month == 1_ik .and. a%day == 1_ik)
      a = time_advance_days(ed_time_t(2000_ik, 1_ik, 1_ik), days_in_year(2000_ik))
      call check_true('+366 days from a leap Jan 1', a%year == 2001_ik .and. a%month == 1_ik .and. a%day == 1_ik)

      !--- Month stepping clamps the day into the shorter target month.
      a = time_advance_months(ed_time_t(2000_ik, 1_ik, 31_ik), 1_ik)
      call check_true('Jan 31 + 1 month -> Feb 29 (leap clamp)', a%month == 2_ik .and. a%day == 29_ik)
      a = time_advance_months(ed_time_t(2001_ik, 1_ik, 31_ik), 1_ik)
      call check_true('Jan 31 + 1 month -> Feb 28 (non-leap clamp)', a%month == 2_ik .and. a%day == 28_ik)
      a = time_advance_months(ed_time_t(2000_ik, 12_ik, 15_ik), 1_ik)
      call check_true('Dec + 1 month rolls the year', a%year == 2001_ik .and. a%month == 1_ik .and. a%day == 15_ik)
      a = time_advance_months(ed_time_t(2000_ik, 6_ik, 10_ik), 12_ik)
      call check_true('+12 months -> same day next year', a%year == 2001_ik .and. a%month == 6_ik .and. a%day == 10_ik)

      call check_true('Jan 1 is DOY 1', day_of_year(ed_time_t(2000_ik, 1_ik, 1_ik)) == 1_ik)
      call check_true('Mar 1 (leap) is DOY 61', day_of_year(ed_time_t(2000_ik, 3_ik, 1_ik)) == 61_ik)
      call check_true('Dec 31 (leap) is DOY 366', day_of_year(ed_time_t(2000_ik, 12_ik, 31_ik)) == 366_ik)
      call check_true('Dec 31 (non-leap) is DOY 365', day_of_year(ed_time_t(2001_ik, 12_ik, 31_ik)) == 365_ik)

      !--- The solar zenith cosine at the equinox and the solstices.
      cz_noon   = solar_cosz(ed_time_t(2001_ik, 3_ik, 21_ik),  43200.0_wp,    0.0_wp)   ! equinox, equator, noon
      cz_mid    = solar_cosz(ed_time_t(2001_ik, 3_ik, 21_ik),      0.0_wp,    0.0_wp)   ! equinox, equator, midnight
      cz_summer = solar_cosz(ed_time_t(2001_ik, 6_ik, 21_ik),    43200.0_wp, 45.0_wp)   ! N-summer solstice noon, 45N
      cz_winter = solar_cosz(ed_time_t(2001_ik,   12_ik, 21_ik), 43200.0_wp, 45.0_wp)   ! N-winter solstice noon, 45N
      call check_true('equinox equator noon: cosz ~ 1', abs(cz_noon - 1.0_wp) < 0.02_wp, cz_noon)
      call check_true('midnight: cosz floored to 0', cz_mid == 0.0_wp, cz_mid)
      call check_true('45N summer noon sun higher than winter', cz_summer > cz_winter .and. cz_winter > 0.0_wp)

      call check_true('2000 spans 366 days',                                                       &
                      time_eq(time_advance_days(ed_time_t(2000_ik, 1_ik, 1_ik), 366_ik), ed_time_t(2001_ik, 1_ik, 1_ik)))
      call check_true('2001 spans 365 days',                                                       &
                      time_eq(time_advance_days(ed_time_t(2001_ik, 1_ik, 1_ik), 365_ik), ed_time_t(2002_ik, 1_ik, 1_ik)))

      a = ed_time_t(2000_ik, 1_ik, 1_ik)
      b = ed_time_t(2000_ik, 1_ik, 2_ik)
      call check_true('Jan 1 < Jan 2', time_lt(a, b))
      call check_true('Jan 2 not < Jan 1', .not. time_lt(b, a))
      call check_true('a date is not before itself', .not. time_lt(a, a))
      call check_true('equal dates compare equal', time_eq(a, ed_time_t(2000_ik, 1_ik, 1_ik)))
      a = ed_time_t(2000_ik, 1_ik, 1_ik, 10_ik, 0_ik, 0_ik)
      b = ed_time_t(2000_ik, 1_ik, 1_ik, 11_ik, 0_ik, 0_ik)
      call check_true('same day, earlier hour is less', time_lt(a, b))

      !--- The decimal calendar year: mid-leap-year is exactly x.5.
      call check('Jan 1 -> integer year', time_to_decimal_year(ed_time_t(2000_ik, 1_ik, 1_ik)), 2000.0_wp, 1.0e-9_wp)
      call check('DOY 184 of 366 -> year + 0.5', time_to_decimal_year(ed_time_t(2000_ik, 7_ik, 2_ik)), 2000.5_wp,    &
                 1.0e-9_wp)

      !--- String round trips and validation.
      call string_to_time('2000-01-01', a, ok)
      call check_true('parse YYYY-MM-DD', ok .and. a%year == 2000_ik .and. a%month == 1_ik .and. a%day == 1_ik)
      call string_to_time('1999-12-31 23:59:59', a, ok)
      call check_true('parse date + time', ok .and. a%hour == 23_ik .and. a%minute == 59_ik .and. a%second == 59_ik)
      call string_to_time('2000/06/15', a, ok)
      call check_true('parse with / separators', ok .and. a%month == 6_ik .and. a%day == 15_ik)
      call string_to_time('2000-02-29', a, ok)
      call check_true('2000-02-29 is a valid (leap) date', ok)
      call string_to_time('2001-02-29', a, ok)
      call check_true('2001-02-29 is rejected (non-leap)', .not. ok)
      call string_to_time('2000-13-01', a, ok)
      call check_true('month 13 is rejected', .not. ok)
      call string_to_time('not-a-date', a, ok)
      call check_true('garbage is rejected', .not. ok)
      call check_true('14-char stamp', time_to_stamp(ed_time_t(2050_ik, 1_ik, 1_ik)) == '20500101000000')
      call check_true('stamp carries the time of day',                                             &
                      time_to_stamp(ed_time_t(1999_ik, 12_ik, 31_ik, 23_ik, 59_ik, 59_ik)) == '19991231235959')
      call check_true('human string', time_to_string(ed_time_t(2000_ik, 1_ik, 1_ik)) == '2000-01-01 00:00:00')

      !---~---
      !   ONE declination behind both solar consumers. solar_cosz and daylength both read
      ! solar_declination, Cooper's form (docs/science/forcing.md eq. (4)); two standard
      ! approximations, Cooper and White (1997), differ by a 1.25-day phase offset. Each
      ! consumer is INVERTED back to a declination rather than the formula re-implemented, so a
      ! divergence is caught even if the formula itself changes.
      !---~---
      block
         real(wp)    :: decl, cz, dl, lat, decl_from_cosz, decl_from_dl, arg
         integer(ik) :: d, doy
         real(wp),    parameter :: PI_ = 3.14159265358979323846_wp
         integer(ik), parameter :: DOYS(4) = [ 46_ik, 105_ik, 196_ik, 288_ik ]
         lat = 42.44_wp
         do d = 1_ik, size(DOYS, kind=ik)
            doy  = DOYS(d)
            decl = solar_declination(doy)
            !--- At solar NOON the hour angle is zero, so cosz = cos(lat - decl) exactly.
            cz = solar_cosz(date_of_doy(2025_ik, doy), 43200.0_wp, lat)
            decl_from_cosz = lat * PI_ / 180.0_wp - acos(min(1.0_wp, max(-1.0_wp, cz)))
            call check('solar_cosz uses solar_declination', decl_from_cosz, decl, 1.0e-12_wp)
            !--- daylength inverts to the same declination through its own acos.
            dl  = daylength(lat, doy)
            arg = cos(dl * PI_ / 24.0_wp)                    ! = -tan(lat) tan(decl)
            decl_from_dl = atan(-arg / tan(lat * PI_ / 180.0_wp))
            call check('daylength uses solar_declination', decl_from_dl, decl, 1.0e-12_wp)
         end do
         !---~---
         !   The shared function is the Cooper form: zero declination at its ascending crossing
         ! (doy 81), positive in NH summer, negative in NH winter.
         !---~---
         call check('declination ~ 0 at doy 81', solar_declination(81_ik), 0.0_wp, 1.0e-3_wp)
         call check_true('declination near +23.45 deg at the June solstice', solar_declination(172_ik) > 0.40_wp)
         call check_true('declination near -23.45 deg at the December solstice', solar_declination(355_ik) < -0.40_wp)
      end block
      !---~---
      !   The day length's polar branches: a polar day at 80 N in June, a polar night in
      ! December, and 12 h at the equator.
      !---~---
      call check_true('daylength: a polar day at 80 N in June is ~ 24 h', daylength(80.0_wp, 172_ik) > 23.5_wp)
      call check_true('daylength: a polar night at 80 N in December is ~ 0 h', daylength(80.0_wp, 355_ik) < 0.5_wp)
      call check_true('daylength: the equator has ~ 12 h', abs(daylength(0.0_wp, 172_ik) - 12.0_wp) < 0.5_wp)
   end subroutine test_calendar

   !--- 10. The sun on a slope: 45 N on the winter solstice, at 0 E (solar noon ~12 UTC).
   subroutine test_slope_sun()
      real(wp), parameter :: LAT = 45.0_wp, LON = 0.0_wp, SLOPE = 30.0_wp
      real(wp), parameter :: DEG = 3.14159265358979323846_wp / 180.0_wp
      type(ed_time_t)   :: noon, morning, afternoon
      real(wp)          :: cosz, flat, south, north, east_am, east_pm
      print '(a)', 'test_slope_sun:'
      noon      = ed_time_t(year=2021, month=12, day=21, hour=12)
      morning   = ed_time_t(year=2021, month=12, day=21, hour=10)
      afternoon = ed_time_t(year=2021, month=12, day=21, hour=14)
      cosz  = solar_cosz_utc(noon, seconds_into_day(noon), LAT, LON)
      flat  = slope_incidence_cosz(noon, seconds_into_day(noon), LAT, LON, 0.0_wp, 0.0_wp)
      south = slope_incidence_cosz(noon, seconds_into_day(noon), LAT, LON, SLOPE, 180.0_wp)
      north = slope_incidence_cosz(noon, seconds_into_day(noon), LAT, LON, SLOPE, 0.0_wp)
      call check('level ground sees the sun at cos z', flat, cosz, 1.0e-12_wp)
      !---~---
      !   At noon the sun is due south, so an equator-facing slope sees it at z - slope; the 12 UTC
      ! instant is a few minutes off solar noon at 0 E, hence the tolerance.
      !---~---
      call check('a 30 deg equator-facing slope sees the noon sun at cos(z - 30 deg)', south,          &
                 cos(acos(cosz) - SLOPE * DEG), 2.0e-3_wp)
      call check_true('a slope facing the pole sees it lower than level ground', north < cosz)
      east_am = slope_incidence_cosz(morning, seconds_into_day(morning), LAT, LON, SLOPE, 90.0_wp)
      east_pm = slope_incidence_cosz(afternoon, seconds_into_day(afternoon), LAT, LON, SLOPE, 90.0_wp)
      call check_true('a slope facing east sees more of the morning sun than of the afternoon''s', east_am > east_pm)
   end subroutine test_slope_sun

   !--- 11. Linear interpolation between two records.
   subroutine test_interpolation()
      print '(a)', 'test_interpolation:'
      call check('linear midpoint', linear_interpolate(10.0_wp, 20.0_wp, 0.5_wp), 15.0_wp, 1.0e-12_wp)
      call check('linear at the first record', linear_interpolate(10.0_wp, 20.0_wp, 0.0_wp), 10.0_wp, 1.0e-12_wp)
      call check('linear at the second record', linear_interpolate(10.0_wp, 20.0_wp, 1.0_wp), 20.0_wp, 1.0e-12_wp)
   end subroutine test_interpolation

   !--- Calendar date of a day-of-year, so the cosz inversion above can be driven by doy.
   pure type(ed_time_t) function date_of_doy(year, doy) result(t)
      integer(ik), intent(in) :: year, doy
      t = time_advance_days(ed_time_t(year = year, month = 1_ik, day = 1_ik), doy - 1_ik)
   end function date_of_doy
end program test_numerics
