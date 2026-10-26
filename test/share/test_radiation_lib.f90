! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_radiation_lib -- ed_radiation_lib, the radiation above the canopy                    !
! (docs/science/forcing.md):                                                                 !
!   1. the shortwave partition, Erbs (clearness index) and Weiss-Norman (band-specific):      !
!      four streams that sum to the total, mostly beam under a clear sky and mostly diffuse   !
!      under an overcast one, none at night, and the Weiss-Norman beam's pressure sensitivity; !
!   2. the longwave synthesis: Brutsaert and Idso-Jackson emissivities against hand values,   !
!      a cloud term that only adds, and the clearness index's night sentinel;                  !
!   3. the beam on a slope: level ground keeps its beam, no beam stays none, and the beam      !
!      normal to a low sun is held at the solar constant.                                      !
!   4. the longwave at another height (forcing.md §8): the file's longwave scales by the       !
!      clear-sky emissivity times T^4, so a warmer, moister site gets more, and dz = 0 none.  !
!==========================================================================================!
program test_radiation_lib
   use ed_test_assert,     only : check, check_true, test_report
   use ed_kinds,           only : wp
   use ed_radiation_lib,   only : LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO, SWPART_CLEARIDX, SWPART_WEISS_NORMAN,  &
                                  partition_shortwave, clearness_index, clear_sky_emissivity,            &
                                  synthesize_lwdown, slope_beam_factor, lapse_longwave
   implicit none

   call test_sw_partition()
   call test_sw_partition_weissnorman()
   call test_lwdown_synthesis()
   call test_slope_beam()
   call test_lapse_longwave()

   call test_report('test_radiation_lib')

contains

   !----- 4. shortwave partition (Erbs clearness-index). -----------------------------------!
   subroutine test_sw_partition()
      real(wp) :: pb, pd, nb, nd, cosz
      real(wp), parameter :: PS = 98000.0_wp
      print '(a)', '-- test 4: SW partition (Erbs) --'
      cosz = 0.8_wp
      call partition_shortwave(800.0_wp, cosz, PS, SWPART_CLEARIDX, pb, pd, nb, nd)
      call check('4 streams sum to total SW', pb+pd+nb+nd, 800.0_wp, 1.0e-9_wp)
      call check_true('all streams >= 0', pb>=0 .and. pd>=0 .and. nb>=0 .and. nd>=0, min(pb,pd,nb,nd))
      call check_true('clear sky -> mostly beam', (pb+nb) > (pd+nd), (pb+nb)-(pd+nd))
      ! overcast (low kt) -> mostly diffuse
      call partition_shortwave(150.0_wp, cosz, PS, SWPART_CLEARIDX, pb, pd, nb, nd)
      call check_true('overcast -> mostly diffuse', (pd+nd) > (pb+nb), (pd+nd)-(pb+nb))
      ! night -> all zero
      call partition_shortwave(0.0_wp, 0.0_wp, PS, SWPART_CLEARIDX, pb, pd, nb, nd)
      call check('night -> zero SW', pb+pd+nb+nd, 0.0_wp, 1.0e-30_wp)
   end subroutine test_sw_partition

   !----- Weiss-Norman band-specific partition: exact energy conservation, clear/overcast beam    !
   !      trend, dawn all-diffuse, and pressure sensitivity of the beam fraction.                   !
   subroutine test_sw_partition_weissnorman()
      real(wp) :: pb, pd, nb, nd, cosz, beam_hi, beam_lo
      print '(a)', '-- test 4b: SW partition (Weiss-Norman) --'
      cosz = 0.8_wp
      call partition_shortwave(800.0_wp, cosz, 98000.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      call check('WN: 4 streams sum to total SW', pb+pd+nb+nd, 800.0_wp, 1.0e-8_wp)
      call check_true('WN: all streams >= 0', pb>=0 .and. pd>=0 .and. nb>=0 .and. nd>=0, min(pb,pd,nb,nd))
      call check_true('WN: clear sky -> mostly beam', (pb+nb) > (pd+nd), (pb+nb)-(pd+nd))
      call check_true('WN: PAR beam > NIR-diffuse in clear sky', pb > nd, pb-nd)
      beam_hi = pb + nb
      ! overcast (low observed vs potential) -> beam collapses, diffuse dominates
      call partition_shortwave(120.0_wp, cosz, 98000.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      call check('WN overcast: streams sum to total', pb+pd+nb+nd, 120.0_wp, 1.0e-8_wp)
      call check_true('WN overcast -> mostly diffuse', (pd+nd) > (pb+nb), (pd+nd)-(pb+nb))
      ! grazing sun (cosz below WN_COSZ_MIN) -> all diffuse, still conserving
      call partition_shortwave(50.0_wp, 0.01_wp, 98000.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      call check('WN dawn: streams sum to total', pb+pd+nb+nd, 50.0_wp, 1.0e-9_wp)
      call check('WN dawn: no beam', pb+nb, 0.0_wp, 1.0e-30_wp)
      ! pressure sensitivity: psurf_pa enters the WN optical depth, so the split MUST change with it
      call partition_shortwave(800.0_wp, cosz, 101325.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      beam_lo = pb + nb
      call check_true('WN: surface pressure changes the beam split', abs(beam_hi-beam_lo) > 1.0e-6_wp, &
                      abs(beam_hi-beam_lo))
   end subroutine test_sw_partition_weissnorman

   !----- LONGWAVE SYNTHESIS (docs/science/forcing.md §11, eq. (17)). The clear-sky emissivities   !
   !      are closed-form, so these are known-answer checks against Brutsaert (1975) and Idso &    !
   !      Jackson (1969) evaluated by hand -- not against the code's own output. What they pin     !
   !      beyond arithmetic: the cloud term only ever ADDS, the night sentinel is distinguishable  !
   !      from an overcast sky, and both emissivity forms stay inside physical bounds.             !
   subroutine test_lwdown_synthesis()
      real(wp) :: eps, lw_clear, lw_cloudy, kt
      real(wp), parameter :: SIG = 5.670374419e-8_wp
      print '(a)', '-- 12. LWdown synthesis (Brutsaert / Idso + cloud term) --'

      !----- Brutsaert at 20 C, q = 0.010, 1013.25 hPa: eps = 1.24 (e/T)^(1/7) = 0.81985. -------!
      eps = clear_sky_emissivity(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp)
      call check('Brutsaert eps at 20 C, q = 0.010', eps, 0.81985_wp, 1.0e-4_wp)
      !----- Idso-Jackson is temperature only: 1 - 0.261 exp(-7.77e-4 dT^2) = 0.80866. ----------!
      eps = clear_sky_emissivity(LW_CLEAR_IDSO, 293.15_wp, 0.010_wp, 101325.0_wp)
      call check('Idso eps at 20 C (humidity-independent)', eps, 0.80866_wp, 1.0e-4_wp)
      call check('Idso ignores humidity entirely',                                          &
                       clear_sky_emissivity(LW_CLEAR_IDSO, 293.15_wp, 0.001_wp, 101325.0_wp),     &
                       clear_sky_emissivity(LW_CLEAR_IDSO, 293.15_wp, 0.030_wp, 101325.0_wp), 1.0e-12_wp)
      !----- Brutsaert does NOT: a moister sky emits more. --------------------------------------!
      call check_true('Brutsaert: a moister sky is more emissive',                                &
             clear_sky_emissivity(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.030_wp, 101325.0_wp) >         &
             clear_sky_emissivity(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.001_wp, 101325.0_wp))

      !----- The full synthesis, clear sky (kt = 1): LW = eps sigma T^4 = 343.33 W/m2. ----------!
      lw_clear = synthesize_lwdown(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp,          &
                                   1.0_wp, 0.22_wp)
      call check('clear-sky LWdown at 20 C', lw_clear, 0.81985_wp * SIG * 293.15_wp**4, 1.0e-2_wp)
      !----- Overcast (kt = 0.3) with a = 0.22: x (1 + 0.22*0.7) = 396.20 W/m2. -----------------!
      lw_cloudy = synthesize_lwdown(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp,         &
                                    0.3_wp, 0.22_wp)
      call check('cloudy-sky LWdown (kt = 0.3, a = 0.22)', lw_cloudy, 396.199_wp, 1.0e-2_wp)
      call check_true('the cloud term only ADDS to the clear-sky flux', lw_cloudy > lw_clear)
      !----- a = 0 recovers the clear-sky value for ANY kt: the term is genuinely optional. -----!
      call check('cloud coefficient 0 => clear sky at any kt',                              &
                       synthesize_lwdown(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp,    &
                                         0.0_wp, 0.0_wp), lw_clear, 1.0e-9_wp)

      !----- The clearness index: a real value by day, a NEGATIVE sentinel at night. A zero here  !
      !      would read as fully overcast and give every night the maximum cloud correction. -----!
      kt = clearness_index(500.0_wp, 0.5_wp)
      call check('clearness index by day = SW/(S0 cosz)', kt, 500.0_wp / (1361.0_wp*0.5_wp), 1.0e-9_wp)
      call check_true('clearness index is a NEGATIVE sentinel at night, not 0',                    &
                      clearness_index(0.0_wp, 0.0_wp) < 0.0_wp, clearness_index(0.0_wp, 0.0_wp))
      call check_true('clearness index is capped at 1', clearness_index(5000.0_wp, 0.5_wp) <= 1.0_wp)
   end subroutine test_lwdown_synthesis

   !----- The beam on a slope (docs/science/forcing.md §8). -----------------------------------------!
   subroutine test_slope_beam()
      real(wp), parameter :: S0 = 1361.0_wp
      print '(a)', '-- the beam on a slope --'
      call check_true('level ground keeps its beam', slope_beam_factor(300.0_wp, 0.5_wp, 0.5_wp) == 1.0_wp)
      call check_true('no beam has nothing to scale', slope_beam_factor(0.0_wp, 0.5_wp, 0.7_wp) == 1.0_wp)
      call check_true('a slope facing the sun takes cos i / cos z',                                  &
                      abs(slope_beam_factor(300.0_wp, 0.5_wp, 0.7_wp) - 1.4_wp) < 1.0e-12_wp)
      call check_true('the beam normal to a low sun is held at the solar constant',                  &
                      abs(slope_beam_factor(100.0_wp, 0.01_wp, 0.5_wp) - S0 * 0.5_wp / 100.0_wp) < 1.0e-9_wp)

   end subroutine test_slope_beam

   !----- 4. The longwave at a site 500 m below its cell, warmer, moister and at a higher pressure. --!
   subroutine test_lapse_longwave()
      real(wp), parameter :: T0 = 285.0_wp, P0 = 95000.0_wp, Q0 = 0.007_wp
      real(wp), parameter :: T1 = 288.25_wp, P1 = 100600.0_wp, Q1 = 0.0078_wp
      real(wp) :: lw1
      print '(a)', '-- 4. the longwave at another height --'
      lw1 = lapse_longwave(300.0_wp, LW_CLEAR_BRUTSAERT, T0, Q0, P0, T1, Q1, P1)
      call check('file longwave scales by eps*T^4', lw1,                                           &
                 300.0_wp * clear_sky_emissivity(LW_CLEAR_BRUTSAERT, T1, Q1, P1) * T1**4           &
                          / (clear_sky_emissivity(LW_CLEAR_BRUTSAERT, T0, Q0, P0) * T0**4), 1.0e-9_wp)
      call check_true('the warmer, moister site gets more longwave', lw1 > 300.0_wp, lw1 - 300.0_wp)
      call check('dz = 0 leaves longwave unchanged',                                               &
                 lapse_longwave(300.0_wp, LW_CLEAR_IDSO, T0, Q0, P0, T0, Q0, P0), 300.0_wp, 1.0e-9_wp)
   end subroutine test_lapse_longwave

end program test_radiation_lib
