! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_fast_integ_schemes -- the fast loop's two schemes and the march that drives them
! (ed_fast_integ_esdirk2, ed_fast_integ_rk45, ed_fast_integ_stepper::adaptive_march;
! docs/science/numerical_scheme.md §3, eqs. (1)-(4)).
!   1. ESDIRK2 (esdirk2_step) on the kernel record: second order on the decoupled CO2
!      twin; at least first order on the coupled soil-top temperature; plant water
!      physical and the embedded estimate bounded at the shipped 900 s; a day of the
!      adaptive march stays bounded.
!   2. ESDIRK2's ADAPTIVE MARCH (adaptive_march): a tighter rtol takes more steps, both
!      agree with a fine fixed march, and the adaptive march is the cheaper.
!   3. RK45 (rk45_step) on the kernel record: fifth order on the CO2 twin and on the
!      coupled soil-top temperature.
!   4. RK45's ADAPTIVE MARCH (adaptive_march's RK45 branch) on the fixture's frozen record,
!      as fast_step drives it: the error follows the tolerance, a tighter one takes more
!      steps, and only accepted steps reach the boundary-flux amounts, which close the
!      canopy air's three ledgers to round-off across a rejected trial.
!---~---
program test_fast_integ_schemes
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : tiny_num
   use ed_fast_params,          only : INTEG_ESDIRK2, INTEG_RK45, RK45_RTOL
   use ed_therm_lib,            only : cas_enthalpy_to_temp, cas_temp_to_enthalpy
   use ed_soil_heat,            only : soil_layer_temp
   use ed_fast_integ_types,     only : march_state_t, frozen_t, error_control_t, surface_tend_t, boundary_flux_t, &
                                       state_init, state_sub
   use ed_fast_integ_backward_euler,   only : backward_euler_stage, advance_plant_water
   use ed_fast_integ_esdirk2,   only : esdirk2_step
   use ed_fast_integ_rk45,      only : rk45_step, RK45_P_ORDER
   use ed_fast_integ_stepper,   only : adaptive_march
   use ed_fast_integ_control,   only : default_error_control, build_error_control, state_wrms_grouped
   use ed_test_fast_fixture,    only : fast_column_t, make_fast_column, reset_column, set_noon_forcing,           &
                                       frozen_record, kernel_record
   use ed_test_assert,          only : check_true, test_report
   implicit none

   integer(ik), parameter :: N = 2_ik            !--- the kernel record's cohorts
   type(fast_column_t)    :: col                 !--- the fixture column: its site's soil, and section 4's record
   type(error_control_t)  :: ec_default          !--- the type's defaults: no budget stop
   integer(ik)            :: nsl

   call make_fast_column(col)
   nsl = col%site%soil_profile%n_active

   call test_esdirk2()
   call test_esdirk2_adaptive()
   call test_rk45_order()
   call test_rk45_adaptive()

   call test_report('test_fast_integ_schemes')

contains

   !---~---
   !   1. ESDIRK2 (ARS(2,2,2), gamma = 1 - 1/sqrt(2); docs/science/numerical_scheme.md eqs.
   !    (1)-(2)). Its tableau is second order. On the canopy air's CO2 twin of the kernel
   !    record, whose biotic source and atmospheric conductance are both held, nothing lags,
   !    and the second order shows. The soil-top temperature is coupled to the canopy air
   !    through the backward-Euler stage, which holds the ground temperature, and the
   !    canopy-air conductances, at their values at the stage's start while its surface
   !    Newton solves the canopy air (ed_fast_integ_backward_euler): the surface sees the
   !    soil a stage late, and the observed order there is about one.
   !---~---
   subroutine test_esdirk2()
      type(march_state_t)   :: y, yref, y1, y2, y4, ynew, yerr, ylo, ytmp
      type(frozen_t)        :: frozen
      type(error_control_t) :: ec
      real(wp)    :: e1, e2, e4, p_lo, p_hi, tref, errnorm, tcas
      integer(ik) :: step, ns, nr
      logical     :: physical
      print '(a)', 'test_esdirk2:'
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)

      !---~---
      !   Self-convergence on the CO2 twin over 1800 s: the reference at 1800/64 s, the solutions at
      ! 450, 225 and 112.5 s.
      !---~---
      call march_esdirk2(y, frozen, 1800.0_wp / 64.0_wp, 64_ik, yref)
      tref = yref%cas_co2
      call march_esdirk2(y, frozen, 450.0_wp, 4_ik, y1)
      call march_esdirk2(y, frozen, 225.0_wp, 8_ik, y2)
      call march_esdirk2(y, frozen, 112.5_wp, 16_ik, y4)
      e1 = abs(y1%cas_co2 - tref)
      e2 = abs(y2%cas_co2 - tref)
      e4 = abs(y4%cas_co2 - tref)
      p_lo = log(e1 / max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2 / max(e4, tiny_num)) / log(2.0_wp)
      print '(a,3es10.3,a,2f6.2)', '   CO2 errors:', e1, e2, e4, ' ; observed order', p_lo, p_hi
      call check_true('esdirk2_step: second order on the decoupled CO2 twin (p >= 1.9)',                          &
                      p_lo >= 1.9_wp .and. p_hi >= 1.9_wp, p_hi)

      !---~---
      !   The same on the soil-top temperature. The first-order backward-Euler step's error at 225 s
      ! is printed beside it; its constant is the record's and can nearly cancel (it does with the
      ! longwave about the tissue's emission temperature), so it bounds nothing.
      !---~---
      call march_esdirk2(y, frozen, 1800.0_wp / 64.0_wp, 64_ik, yref)
      tref = soil_top_temp(yref)
      call march_esdirk2(y, frozen, 450.0_wp, 4_ik, y1)
      call march_esdirk2(y, frozen, 225.0_wp, 8_ik, y2)
      call march_esdirk2(y, frozen, 112.5_wp, 16_ik, y4)
      e1 = abs(soil_top_temp(y1) - tref)
      e2 = abs(soil_top_temp(y2) - tref)
      e4 = abs(soil_top_temp(y4) - tref)
      p_lo = log(e1 / max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2 / max(e4, tiny_num)) / log(2.0_wp)
      call march_be_euler(y, frozen, 225.0_wp, 8_ik, ytmp)
      print '(a,3es10.3,a,2f6.2,a,es10.3)', '   soil-top T errors:', e1, e2, e4, ' ; observed order', p_lo, p_hi, &
            ' ; first-order step at 225 s', abs(soil_top_temp(ytmp) - tref)
      call check_true('esdirk2_step: at least first order on the coupled soil-top temperature (p >= 0.9)',        &
                      p_lo >= 0.9_wp .and. p_hi >= 0.9_wp, min(p_lo, p_hi))

      !--- Plant water stays finite and positive at the shipped 900 s over 6 h.
      call state_init(y, N, nsl, y1)
      physical = .true.
      do step = 1_ik, 24_ik
         call esdirk2_step(y1, frozen, col%site, ec_default, N, nsl, 900.0_wp, ytmp, yerr, col%ws%stages)
         call state_init(ytmp, N, nsl, y1)
         physical = physical .and. all(y1%leaf_water_mass(1:N) == y1%leaf_water_mass(1:N))                         &
                    .and. all(y1%wood_water_mass(1:N) == y1%wood_water_mass(1:N))                                  &
                    .and. minval(y1%leaf_water_mass(1:N)) > 0.0_wp .and. minval(y1%wood_water_mass(1:N)) > 0.0_wp
      end do
      call check_true('esdirk2_step: plant water stays physical at dt = 900 s', physical, minval(y1%leaf_water_mass(1:N)))

      !---~---
      !   The embedded estimate at 900 s, in the production norm the march accepts a step by, is
      ! bounded (not detonating): measured against the first-order embedded solution, relative to
      ! the step's start (adaptive_march).
      !---~---
      ec = default_error_control(1.0e-3_wp)
      call esdirk2_step(y, frozen, col%site, ec_default, N, nsl, 900.0_wp, ynew, yerr, col%ws%stages)
      call state_sub(ynew, yerr, N, nsl, ylo)
      errnorm = state_wrms_grouped(ynew, ylo, y, N, nsl, ec%tols)
      print '(a,es10.3)', '   embedded estimate at 900 s, in tolerance units:', errnorm
      call check_true('esdirk2_step: the embedded estimate is bounded at dt = 900 s (< 50 tolerance units)',       &
                      errnorm < 50.0_wp, errnorm)

      !--- Stiffness: a day of the adaptive march stays physical and bounded.
      ec = default_error_control(1.0e-3_wp)
      call adaptive_march(INTEG_ESDIRK2, y, frozen, col%site, N, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, col%ws, ns, nr)
      tcas = cas_enthalpy_to_temp(ytmp%cas_enthalpy, ytmp%cas_shv)
      print '(a,i0,a,i0)', '   adaptive ESDIRK2, 24 h: steps = ', ns, ', rejects = ', nr
      call check_true('adaptive_march (ESDIRK2): a day stays bounded (280 < T_cas < 320 K)',                       &
                      tcas > 280.0_wp .and. tcas < 320.0_wp, tcas)
   end subroutine test_esdirk2

   !---~---
   !   2. ESDIRK2's ADAPTIVE MARCH under the production controller (adaptive_march), so the
   !    test fails for the reasons a run would: over 30 min from cool, dry canopy air, a
   !    tighter rtol takes more steps, both agree with a fixed march at 2 s, and the adaptive
   !    march takes fewer steps than that fixed one. The canopy is ventilated four times as
   !    strongly as the kernel record's.
   !---~---
   subroutine test_esdirk2_adaptive()
      type(march_state_t)   :: y, y1, y2, yr
      type(frozen_t)        :: frozen
      type(error_control_t) :: ec
      real(wp)    :: tc1, tc2, tcr
      integer(ik) :: ns1, ns2, nr1, nr2
      print '(a)', 'test_esdirk2_adaptive:'
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      frozen%cas%box%atm_conductance_enthalpy = 4.0_wp * frozen%cas%box%atm_conductance_enthalpy
      frozen%cas%box%atm_conductance_vapor    = 4.0_wp * frozen%cas%box%atm_conductance_vapor
      frozen%cas%box%atm_conductance_co2      = 4.0_wp * frozen%cas%box%atm_conductance_co2
      call set_cas(y, 290.0_wp, 0.009_wp)
      ec = default_error_control(1.0e-3_wp)
      call adaptive_march(INTEG_ESDIRK2, y, frozen, col%site, N, nsl, 1800.0_wp, ec, 50.0_wp, y1, col%ws, ns1, nr1)
      ec = default_error_control(1.0e-6_wp)
      call adaptive_march(INTEG_ESDIRK2, y, frozen, col%site, N, nsl, 1800.0_wp, ec, 50.0_wp, y2, col%ws, ns2, nr2)
      call march_esdirk2(y, frozen, 2.0_wp, 900_ik, yr)
      tc1 = cas_enthalpy_to_temp(y1%cas_enthalpy, y1%cas_shv)
      tc2 = cas_enthalpy_to_temp(y2%cas_enthalpy, y2%cas_shv)
      tcr = cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv)
      print '(a,i0,a,i0,a)', '   steps: rtol = 1e-3 -> ', ns1, ', rtol = 1e-6 -> ', ns2, '; fixed 2 s -> 900'
      call check_true('adaptive_march (ESDIRK2): rtol = 1e-3 agrees with the fine fixed march (< 0.5 K)',          &
                      abs(tc1 - tcr) < 0.5_wp, tc1 - tcr)
      call check_true('adaptive_march (ESDIRK2): rtol = 1e-6 agrees with the fine fixed march (< 0.5 K)',          &
                      abs(tc2 - tcr) < 0.5_wp, tc2 - tcr)
      call check_true('adaptive_march (ESDIRK2): a tighter rtol takes more steps', ns2 > ns1, real(ns2 - ns1, wp))
      call check_true('adaptive_march (ESDIRK2): cheaper than the fine fixed march (steps < 900)', ns1 < 900_ik,     &
                      real(ns1, wp))
      call check_true('adaptive_march (ESDIRK2): the march made progress (steps > 0)', ns1 > 0_ik, real(ns1, wp))
   end subroutine test_esdirk2_adaptive

   !---~---
   !   3. RK45 (Cash-Karp 5(4)): observed order. RK45 evaluates the surface and the soil heat
   !    together in every stage (docs/science/numerical_scheme.md §3), so the fifth order
   !    shows on the coupled soil-top temperature too, not only on the decoupled CO2 twin.
   !    A fifth-order error ~ C h^5 falls to round-off at a much coarser h than a second-
   !    order one, so h stays large enough for truncation to dominate at the finest point:
   !    over 2400 s, steps of 300, 150 and 75 s (300 s is under the ~360 s explicit limit the
   !    canopy air's ~130 s time scale sets; this compares accuracy, not stability), against
   !    a reference at 2400/1024 s, 32 times finer than the finest point. p = log2(e(h) /
   !    e(h/2)) approaches 5 as h shrinks; p >= 4.5 allows for pre-asymptotic effects.
   !---~---
   subroutine test_rk45_order()
      type(march_state_t) :: y, yref, y1, y2, y4, y8
      type(frozen_t)      :: frozen
      real(wp) :: e1, e2, e4, e8, p_lo, p_hi, p_fine, tref, noise
      print '(a)', 'test_rk45_order:'
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)

      !--- The CO2 twin: the tableau itself is fifth order on this right-hand side.
      call march_rk45(y, frozen, 2400.0_wp / 1024.0_wp, 1024_ik, yref)
      tref = yref%cas_co2
      call march_rk45(y, frozen, 300.0_wp, 8_ik, y1)
      call march_rk45(y, frozen, 150.0_wp, 16_ik, y2)
      call march_rk45(y, frozen, 75.0_wp, 32_ik, y4)
      e1 = abs(y1%cas_co2 - tref)
      e2 = abs(y2%cas_co2 - tref)
      e4 = abs(y4%cas_co2 - tref)
      p_lo = log(e1 / max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2 / max(e4, tiny_num)) / log(2.0_wp)
      print '(a,3es10.3,a,2f6.2)', '   CO2 errors:', e1, e2, e4, ' ; observed order', p_lo, p_hi
      call check_true('rk45_step: fifth order on the decoupled CO2 twin (p >= 4.5)',                              &
                      p_lo >= 4.5_wp .and. p_hi >= 4.5_wp, p_hi)

      !---~---
      !   The coupled soil-top temperature, ESDIRK2's first-order variable. Four points, and the
      ! order asserted on the two finest pairs: the most strongly coupled variable carries a
      ! visible h^6 term at the coarsest h, so p across the coarsest pair is not asymptotic and
      ! reads below or above 5 by whichever way that term points. The soil water moves at the
      ! network's steady rate, so the soil-top temperature is smooth and its error reaches the
      ! reference's own round-off (~1e3 eps T) by 75 s; the order is asserted above that floor.
      !---~---
      call march_rk45(y, frozen, 2400.0_wp / 1024.0_wp, 1024_ik, yref)
      tref = soil_top_temp(yref)
      call march_rk45(y, frozen, 300.0_wp, 8_ik, y1)
      call march_rk45(y, frozen, 150.0_wp, 16_ik, y2)
      call march_rk45(y, frozen, 75.0_wp, 32_ik, y4)
      call march_rk45(y, frozen, 37.5_wp, 64_ik, y8)
      e1 = abs(soil_top_temp(y1) - tref)
      e2 = abs(soil_top_temp(y2) - tref)
      e4 = abs(soil_top_temp(y4) - tref)
      e8 = abs(soil_top_temp(y8) - tref)
      p_lo   = log(e1 / max(e2, tiny_num)) / log(2.0_wp)
      p_hi   = log(e2 / max(e4, tiny_num)) / log(2.0_wp)
      p_fine = log(e4 / max(e8, tiny_num)) / log(2.0_wp)
      print '(a,4es10.3,a,3f6.2)', '   soil-top T errors:', e1, e2, e4, e8, ' ; observed order', p_lo, p_hi, p_fine
      noise = 1.0e3_wp * epsilon(1.0_wp) * tref
      call check_true('rk45_step: fifth order on the coupled soil-top temperature (p >= 4.5 above round-off)',     &
                      (p_hi >= 4.5_wp .or. e4 < noise) .and. (p_fine >= 4.5_wp .or. e8 < noise), p_fine)
   end subroutine test_rk45_order

   !---~---
   !   4. RK45's ADAPTIVE MARCH, as fast_step drives it: the fixture's frozen record for a
   !    900 s step at noon (live canopy-air conductances, the tissue heat store, the water
   !    network's flows), the run's error control at RK45's estimate order, from a first trial
   !    over the whole step. Its accuracy follows the tolerance: rtol 1e-4 (RK45's shipped
   !    one) and 1e-7, with atol_scale tightened alike so the absolute floor does not mask
   !    rtol, against 900 fixed 1 s steps. And the boundary-flux amounts sum only the accepted
   !    steps: then the canopy air's energy, water and CO2 change by exactly what crossed its
   !    boundary, across the first, rejected trial; summing a rejected trial breaks that by
   !    the trial's whole flux.
   !---~---
   subroutine test_rk45_adaptive()
      real(wp), parameter :: H = 900.0_wp
      type(march_state_t)   :: ys, y1, y2, yr
      type(frozen_t)        :: fz
      type(error_control_t) :: ec
      type(boundary_flux_t) :: acc
      real(wp)    :: d1, d2, cap_m, cap_c, r_enth, r_vap, r_co2
      integer(ik) :: ns1, ns2, nr1, nr2, n
      print '(a)', 'test_rk45_adaptive:'
      call reset_column(col)
      call set_noon_forcing(col)
      call frozen_record(col, fz, ys, H)
      n = col%static%n
      col%cfg%rtol = RK45_RTOL
      col%cfg%atol_scale = 1.0_wp
      ec = build_error_control(col%cfg)
      ec%p_order = RK45_P_ORDER
      call adaptive_march(INTEG_RK45, ys, fz, col%site, n, nsl, H, ec, H, y1, col%ws, ns1, nr1, acc=acc)
      col%cfg%rtol = 1.0e-3_wp * RK45_RTOL
      col%cfg%atol_scale = 1.0e-3_wp
      ec = build_error_control(col%cfg)
      ec%p_order = RK45_P_ORDER
      call adaptive_march(INTEG_RK45, ys, fz, col%site, n, nsl, H, ec, H, y2, col%ws, ns2, nr2)
      col%cfg%rtol = RK45_RTOL
      col%cfg%atol_scale = 1.0_wp
      call march_rk45_on(ys, fz, n, 1.0_wp, 900_ik, yr)
      d1 = abs(cas_enthalpy_to_temp(y1%cas_enthalpy, y1%cas_shv) - cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv))
      d2 = abs(cas_enthalpy_to_temp(y2%cas_enthalpy, y2%cas_shv) - cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv))
      print '(a,i0,a,i0,a,i0,a,i0,a,2es10.3)', '   RK45 steps (rejects): rtol 1e-4 -> ', ns1, ' (', nr1,             &
            '), 1e-7 -> ', ns2, ' (', nr2, '); canopy-air T error', d1, d2
      call check_true('adaptive_march (RK45): within 0.01 K of the fine march at its shipped rtol', d1 < 1.0e-2_wp, d1)
      call check_true('adaptive_march (RK45): a tolerance 1000 times tighter is at least 10 times closer',          &
                      d2 < 0.1_wp * d1 .or. d2 < 1.0e-9_wp, d2)
      call check_true('adaptive_march (RK45): a tighter tolerance takes more steps', ns2 > ns1, real(ns2 - ns1, wp))
      call check_true('adaptive_march (RK45): the first trial over the whole step was rejected (the next is live)',  &
                      nr1 > 0_ik, real(nr1, wp))
      cap_m = fz%cas%box%air_mass_capacity
      cap_c = fz%cas%box%air_molar_capacity
      r_enth = cap_m * (y1%cas_enthalpy - ys%cas_enthalpy) - (acc%cas_enth_in - acc%cas_enth_out)
      r_vap  = cap_m * (y1%cas_shv - ys%cas_shv) - (acc%cas_vap_in - acc%cas_vap_out)
      r_co2  = cap_c * (y1%cas_co2 - ys%cas_co2) - (acc%cas_co2_in - acc%cas_co2_out)
      call check_true('adaptive_march (RK45): only accepted steps are summed (canopy-air energy closes)',           &
                      abs(r_enth) < 1.0e-9_wp * max(abs(acc%cas_enth_in), 1.0_wp), r_enth)
      call check_true('adaptive_march (RK45): only accepted steps are summed (canopy-air water closes)',            &
                      abs(r_vap) < 1.0e-9_wp * max(abs(acc%cas_vap_in), 1.0e-6_wp), r_vap)
      call check_true('adaptive_march (RK45): only accepted steps are summed (canopy-air CO2 closes)',              &
                      abs(r_co2) < 1.0e-9_wp * max(abs(acc%cas_co2_in), 1.0_wp), r_co2)
   end subroutine test_rk45_adaptive

   !--- The soil-top temperature of a march state.
   function soil_top_temp(y) result(t)
      type(march_state_t), intent(in) :: y
      real(wp) :: t, fl
      call soil_layer_temp(y%soil_energy(1), y%theta(1), col%site%soil_profile%dry_heat_capacity(1),             &
                           col%site%soil_profile%layer_curve(1), t, fl)
   end function soil_top_temp

   !--- Canopy air at a temperature and humidity.
   subroutine set_cas(y, t, q)
      type(march_state_t), intent(inout) :: y
      real(wp),            intent(in)    :: t, q
      y%cas_shv = q
      y%cas_enthalpy = cas_temp_to_enthalpy(t, q)
   end subroutine set_cas

   !--- nstep fixed ESDIRK2 steps of dt on the kernel record, the embedded estimate discarded.
   subroutine march_esdirk2(y0, frozen, dt, nstep, y_out)
      type(march_state_t), intent(in)    :: y0
      type(frozen_t),      intent(in)    :: frozen
      real(wp),            intent(in)    :: dt
      integer(ik),         intent(in)    :: nstep
      type(march_state_t), intent(inout) :: y_out
      type(march_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call state_init(y0, N, nsl, y)
      do s = 1_ik, nstep
         call esdirk2_step(y, frozen, col%site, ec_default, N, nsl, dt, ytmp, yerr, col%ws%stages)
         call state_init(ytmp, N, nsl, y)
      end do
      call state_init(y, N, nsl, y_out)
   end subroutine march_esdirk2

   !---~---
   !   nstep first-order backward-Euler steps of dt: the gamma = 1 degraded ESDIRK2, one stage and
   ! the plant-water split on its own transpiration (test_fast_integ_kernels' probe).
   !---~---
   subroutine march_be_euler(y0, frozen, dt, nstep, y_out)
      type(march_state_t), intent(in)    :: y0
      type(frozen_t),      intent(in)    :: frozen
      real(wp),            intent(in)    :: dt
      integer(ik),         intent(in)    :: nstep
      type(march_state_t), intent(inout) :: y_out
      type(march_state_t)  :: y, ytmp
      type(surface_tend_t) :: surf_tend
      integer(ik) :: s
      call state_init(y0, N, nsl, y)
      do s = 1_ik, nstep
         call backward_euler_stage(y, frozen, col%site, ec_default, N, nsl, dt, ytmp, surf_tend)
         call advance_plant_water(y, frozen, N, dt, surf_tend%transp_c(1:N), ytmp)
         call state_init(ytmp, N, nsl, y)
      end do
      call state_init(y, N, nsl, y_out)
   end subroutine march_be_euler

   !--- nstep fixed Cash-Karp steps of dt on the kernel record, the fifth-order solution kept.
   subroutine march_rk45(y0, frozen, dt, nstep, y_out)
      type(march_state_t), intent(in)    :: y0
      type(frozen_t),      intent(in)    :: frozen
      real(wp),            intent(in)    :: dt
      integer(ik),         intent(in)    :: nstep
      type(march_state_t), intent(inout) :: y_out
      call march_rk45_on(y0, frozen, N, dt, nstep, y_out)
   end subroutine march_rk45

   subroutine march_rk45_on(y0, frozen, n, dt, nstep, y_out)
      type(march_state_t), intent(in)    :: y0
      type(frozen_t),      intent(in)    :: frozen
      integer(ik),         intent(in)    :: n, nstep
      real(wp),            intent(in)    :: dt
      type(march_state_t), intent(inout) :: y_out
      type(march_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call state_init(y0, n, nsl, y)
      do s = 1_ik, nstep
         call rk45_step(y, frozen, col%site, n, nsl, dt, ytmp, yerr)
         call state_init(ytmp, n, nsl, y)
      end do
      call state_init(y, n, nsl, y_out)
   end subroutine march_rk45_on

end program test_fast_integ_schemes
