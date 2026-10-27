! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_canopy_air_space -- the canopy-air box (ed_canopy_air_space;
! docs/science/canopy_air_space_biophysics.md eqs. (1)-(2)), stepped by
! cas_box_step_implicit.
!   1. ENTHALPY (eq. (2)): each step solves the box's backward-Euler enthalpy balance, and
!      the box warms toward a warmer atmosphere and its sources.
!   2. CO2 (eq. (1)): at the atmosphere's CO2 it holds still; displaced, it relaxes toward
!      the atmosphere without overshoot, more slowly over a weaker conductance; under a
!      constant source it settles at co2_atm + source/conductance; uptake draws it down
!      and respiration pushes it up.
!---~---
program test_canopy_air_space
   use ed_test_assert,        only : check, check_true, test_report
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : mmdry
   use ed_therm_lib,          only : cas_temp_to_enthalpy, cas_enthalpy_to_temp
   use ed_canopy_air_space,   only : cas_env_t, cas_box_step_implicit
   implicit none

   call test_cas()
   call test_cas_co2()

   call test_report('test_canopy_air_space')

contains

   !---~---
   !   1. The enthalpy twin (eq. (2)), the vapour and CO2 twins held inert (no sources and no
   ! vapour conductance; a nonzero molar capacity avoids a 0/0 in the unused CO2 branch). Each
   ! step's backward-Euler balance, M (h_new - h) = dt (S + g (h_atm - h_new)), is assembled
   ! here from the step's result.
   !---~---
   subroutine test_cas()
      type(cas_env_t) :: cas_env
      real(wp) :: enth, shv, temp, resid, enthalpy_atm, worst, enth_new, shv_new, co2_new
      real(wp), parameter :: CAS_MASS_CAPACITY = 1.2_wp * 20.0_wp     ! rho_air * can_depth  [kg/m2]
      real(wp), parameter :: GATM = 1.2_wp * 0.3_wp * 1.0_wp  ! rho_air * ustar * heat_transfer_factor(c3)  [kg/m2/s]
      integer(ik) :: step
      print '(a)', 'test_cas:'
      shv  = 0.012_wp
      enth = cas_temp_to_enthalpy(295.0_wp, shv)                  ! CAS starts at 295 K
      temp = 295.0_wp
      enthalpy_atm = cas_temp_to_enthalpy(300.0_wp, shv)              ! warmer atmosphere
      cas_env%air_mass_capacity        = CAS_MASS_CAPACITY
      cas_env%air_molar_capacity       = 1.0_wp                    ! inert CO2 twin (avoid 0/0)
      cas_env%atm_conductance_enthalpy = GATM
      cas_env%atm_conductance_vapor    = 0.0_wp
      cas_env%atm_conductance_co2      = 0.0_wp
      cas_env%atm_enthalpy             = enthalpy_atm
      cas_env%atm_specific_humidity    = 0.0_wp
      cas_env%atm_co2                  = 0.0_wp
      cas_env%surface_enthalpy_source  = 50.0_wp + 20.0_wp + 10.0_wp + 5.0_wp  ! cohort + ground sensible + vapour-enthalpy
      cas_env%surface_vapor_source     = 0.0_wp
      cas_env%biotic_co2_source        = 0.0_wp
      worst = 0.0_wp
      do step = 1_ik, 30_ik
         call cas_box_step_implicit(enth, shv, 0.0_wp, cas_env, 60.0_wp, enth_new, shv_new, co2_new)
         resid = cas_env%air_mass_capacity * (enth_new - enth)                                   &
                 - 60.0_wp * (cas_env%surface_enthalpy_source                                    &
                              + cas_env%atm_conductance_enthalpy * (cas_env%atm_enthalpy - enth_new))
         worst = max(worst, abs(resid))
         enth = enth_new
         shv = shv_new
         temp = cas_enthalpy_to_temp(enth, shv)
      end do
      call check_true('CAS energy residual ~ 0 (30 steps)', worst < 1.0e-6_wp, worst)
      call check_true('CAS warms toward atmosphere + surfaces', temp > 295.0_wp, temp)
   end subroutine test_cas

   !---~---
   !   2. The CO2 twin of the same box (docs/science/canopy_air_space_biophysics.md eq. (1)). At
   ! the atmosphere's CO2 with no source it holds still; displaced, it relaxes toward the
   ! atmosphere without overshoot at any step, and more slowly over a weaker conductance;
   ! under a constant source it settles at co2_atm + source/conductance; uptake draws it
   ! down and respiration pushes it up.
   !---~---
   subroutine test_cas_co2()
      real(wp), parameter :: CO2_ATM = 400.0_wp, RHO_AIR = 1.2_wp, USTAR = 0.3_wp, SOURCE = 5.0_wp
      real(wp) :: g_co2, c, c_weak
      integer(ik) :: step
      print '(a)', 'test_cas_co2:'
      g_co2 = RHO_AIR / mmdry * USTAR                    ! [mol/m2/s] dry air, transfer factor 1
      c = cas_co2_after(CO2_ATM, 0.0_wp, g_co2, 1800.0_wp)
      call check('no source at the atmosphere''s CO2: unchanged', c, CO2_ATM, 1.0e-9_wp)
      c = cas_co2_after(500.0_wp, 0.0_wp, g_co2, 60.0_wp)
      call check_true('displaced: relaxes toward the atmosphere', c > CO2_ATM .and. c < 500.0_wp, c)
      c_weak = cas_co2_after(500.0_wp, 0.0_wp, 0.5_wp * g_co2, 60.0_wp)
      call check_true('a weaker conductance relaxes more slowly', c_weak > c, c_weak - c)
      c = cas_co2_after(500.0_wp, 0.0_wp, g_co2, 1.0e9_wp)
      call check_true('a huge step does not overshoot the atmosphere (L-stable)', c >= CO2_ATM, c)
      call check('a huge step lands on the atmosphere', c, CO2_ATM, 1.0e-3_wp)
      c = CO2_ATM
      do step = 1_ik, 4000_ik
         c = cas_co2_after(c, SOURCE, g_co2, 60.0_wp)
      end do
      call check('a constant source settles at co2_atm + source/conductance', c, CO2_ATM + SOURCE / g_co2, 1.0e-4_wp)
      c = cas_co2_after(CO2_ATM, -15.0_wp, g_co2, 300.0_wp)
      call check_true('uptake draws the canopy air below the atmosphere', c < CO2_ATM, c)
      c = cas_co2_after(CO2_ATM, 10.0_wp, g_co2, 300.0_wp)
      call check_true('respiration pushes the canopy air above the atmosphere', c > CO2_ATM, c)
   end subroutine test_cas_co2

   !---~---
   !   One step of the CAS box's CO2, on a 20 m column of dry air at 1.2 kg/m3; the enthalpy and
   ! vapour twins are held inert (a nonzero mass capacity avoids their 0/0).
   !---~---
   function cas_co2_after(co2, source, g_co2, dt) result(co2_new)
      real(wp), intent(in) :: co2, source, g_co2, dt
      real(wp) :: co2_new, enth_new, shv_new
      type(cas_env_t) :: env
      env%air_mass_capacity   = 1.0_wp
      env%air_molar_capacity  = 1.2_wp / mmdry * 20.0_wp
      env%atm_conductance_co2 = g_co2
      env%atm_co2             = 400.0_wp
      env%biotic_co2_source   = source
      call cas_box_step_implicit(0.0_wp, 0.0_wp, co2, env, dt, enth_new, shv_new, co2_new)
   end function cas_co2_after

end program test_canopy_air_space
