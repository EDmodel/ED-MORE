! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_demography -- unit tests for the demographic laws (meds_demography;                 !
! docs/science/demography.md §1-§4). Each law is checked by a property it must have, not   !
! by its own formula typed again:                                                          !
!   1. GROWTH        : the prospective rate inverts the wood allometry; no wood growth, no   !
!                      diameter growth; a loss beyond the wood stops at no wood.             !
!   2. MORTALITY     : the hazard is gamma + alpha at zero growth and falls to gamma as      !
!                      growth rises; it reads the growth average once there is one.         !
!   3. TO THE FLOOR  : the decline takes the density exactly to the floor over its window,  !
!                      and is zero at or below the floor.                                   !
!   4. RECRUITMENT   : the efficiency is the open value at no LAI and halves at ln2/k; a     !
!                      cohort's recruitment is linear in its density and zero without seed.  !
!   5. DECISIONS     : a pool founds a cohort at the threshold and not below; the disturbed  !
!                      share composes over two half-intervals; a cohort at the survival     !
!                      height dies in the gap.                                              !
!==========================================================================================!
program test_demography
   use meds_test_assert,   only : check_close, check_true, test_report
   use meds_kinds,         only : wp
   use meds_allometry_lib, only : allometry_t, size_to_wood_carbon, dbh_to_height
   use meds_demography,    only : npp_to_growth, camac_mortality, decline_to_floor_mortality,          &
                                  recruit_efficiency, npp_to_recruitment, founds_cohort,               &
                                  disturbed_fraction, survives_treefall
   implicit none

   real(wp), parameter :: RHO = 0.6_wp                       ! [g/cm3] wood density
   type(allometry_t), parameter :: ALLOM = allometry_t(hgt_max = 35.0_wp, aboveground_frac = 0.7_wp)

   call test_growth()
   call test_mortality()
   call test_decline_to_floor()
   call test_recruitment()
   call test_decisions()

   call test_report('test_demography')

contains

   !----- 1. The prospective growth rate: wood carbon for 20 cm plus the wood of 21 cm less 20 cm !
   !      gives 1 cm over the step, through the allometry's own inverse. -----------------------!
   subroutine test_growth()
      real(wp) :: w20, w21, w0
      w20 = size_to_wood_carbon(ALLOM, 20.0_wp, dbh_to_height(ALLOM, 20.0_wp), RHO)
      w21 = size_to_wood_carbon(ALLOM, 21.0_wp, dbh_to_height(ALLOM, 21.0_wp), RHO)
      call check_close('growth: the wood of 1 cm grows 1 cm/yr', npp_to_growth(ALLOM, w20, w21 - w20, 20.0_wp,   &
                       1.0_wp, RHO), 1.0_wp, 1.0e-6_wp)
      call check_close('growth: half the step, twice the rate', npp_to_growth(ALLOM, w20, w21 - w20, 20.0_wp,   &
                       0.5_wp, RHO), 2.0_wp, 2.0e-6_wp)
      call check_close('growth: no wood growth, no diameter growth', npp_to_growth(ALLOM, w20, 0.0_wp, 20.0_wp,  &
                       1.0_wp, RHO), 0.0_wp, 1.0e-6_wp)
      !----- A loss beyond the wood stops at no wood: the rate is that of losing all of it. ---!
      w0 = npp_to_growth(ALLOM, w20, -w20, 20.0_wp, 1.0_wp, RHO)
      call check_close('growth: a loss beyond the wood stops at no wood',                                    &
                       npp_to_growth(ALLOM, w20, -2.0_wp * w20, 20.0_wp, 1.0_wp, RHO), w0, 1.0e-12_wp)
      call check_true ('growth: losing wood shrinks the diameter', w0 < 0.0_wp)
   end subroutine test_growth

   !----- 2. The Camac hazard: its two limits, its monotone fall, and which growth it reads. --!
   subroutine test_mortality()
      real(wp), parameter :: GAMMA = 0.01_wp, ALPHA = 0.2_wp, BETA = 1.5_wp
      real(wp) :: m0, m1, m2
      m0 = camac_mortality(0.0_wp, 9.0_wp, .true., GAMMA, ALPHA, BETA)
      call check_close('mortality: gamma + alpha at zero growth', m0, GAMMA + ALPHA, 1.0e-14_wp)
      call check_close('mortality: falls to gamma as growth rises',                                         &
                       camac_mortality(50.0_wp, 0.0_wp, .true., GAMMA, ALPHA, BETA), GAMMA, 1.0e-14_wp)
      m1 = camac_mortality(0.5_wp, 0.0_wp, .true., GAMMA, ALPHA, BETA)
      m2 = camac_mortality(1.0_wp, 0.0_wp, .true., GAMMA, ALPHA, BETA)
      call check_true ('mortality: faster growth, lower hazard', m0 > m1 .and. m1 > m2)
      call check_close('mortality: with no average yet it reads the prospective rate',                      &
                       camac_mortality(9.0_wp, 0.0_wp, .false., GAMMA, ALPHA, BETA), m0, 1.0e-14_wp)
   end subroutine test_mortality

   !----- 3. The decline to the floor: n exp(-m window) is the floor; nothing at the floor. --!
   subroutine test_decline_to_floor()
      real(wp), parameter :: FLOOR = 1.0e-8_wp, WINDOW = 0.25_wp
      real(wp) :: m
      m = decline_to_floor_mortality(0.3_wp, FLOOR, WINDOW)
      call check_close('floor: the decline lands on the floor over its window', 0.3_wp * exp(-m * WINDOW),   &
                       FLOOR, 1.0e-20_wp)
      call check_close('floor: no decline at the floor', decline_to_floor_mortality(FLOOR, FLOOR, WINDOW),    &
                       0.0_wp, 0.0_wp)
      call check_close('floor: no decline below it', decline_to_floor_mortality(0.5_wp * FLOOR, FLOOR, WINDOW), &
                       0.0_wp, 0.0_wp)
   end subroutine test_decline_to_floor

   !----- 4. The recruit efficiency and a cohort's recruitment. -----------------------------!
   subroutine test_recruitment()
      real(wp), parameter :: E_OPEN = 0.3_wp, K = 0.5_wp
      real(wp) :: r1, r2
      call check_close('recruit: the open efficiency at no LAI', recruit_efficiency(E_OPEN, K, 0.0_wp),       &
                       E_OPEN, 1.0e-15_wp)
      call check_close('recruit: half of it at LAI = ln2/k', recruit_efficiency(E_OPEN, K, log(2.0_wp) / K),  &
                       0.5_wp * E_OPEN, 1.0e-15_wp)
      r1 = npp_to_recruitment(0.1_wp, 0.02_wp, 1.0_wp / 365.0_wp, 0.15_wp, 0.01_wp)
      r2 = npp_to_recruitment(0.2_wp, 0.02_wp, 1.0_wp / 365.0_wp, 0.15_wp, 0.01_wp)
      call check_close('recruit: twice the plants, twice the recruits', r2, 2.0_wp * r1, 1.0e-12_wp * r2)
      call check_close('recruit: no seed, no recruits', npp_to_recruitment(0.1_wp, 0.0_wp, 1.0_wp / 365.0_wp,  &
                       0.15_wp, 0.01_wp), 0.0_wp, 0.0_wp)
      !----- The recruits carry the carbon the efficiency passes on: rec * dt * c_min. -------!
      call check_close('recruit: the recruits carry the carbon passed on', r1 * (1.0_wp / 365.0_wp) * 0.01_wp,  &
                       0.1_wp * 0.02_wp * 0.15_wp, 1.0e-15_wp)
   end subroutine test_recruitment

   !----- 5. The founding test, the disturbed share and treefall survival. --------------------!
   subroutine test_decisions()
      real(wp) :: f1, fh
      call check_true ('founding: a pool at the threshold founds a cohort', founds_cohort(1.0e-2_wp, 1.0e-2_wp))
      call check_true ('founding: a smaller pool carries over', .not. founds_cohort(0.99e-2_wp, 1.0e-2_wp))
      f1 = disturbed_fraction(0.014_wp, 1.0_wp)
      fh = disturbed_fraction(0.014_wp, 0.5_wp)
      call check_close('disturbance: two half-years disturb what one year does', 1.0_wp - (1.0_wp - fh)**2, f1, &
                       1.0e-15_wp)
      call check_close('disturbance: no rate, no gaps', disturbed_fraction(0.0_wp, 1.0_wp), 0.0_wp, 0.0_wp)
      call check_true ('disturbance: never the whole stand', disturbed_fraction(100.0_wp, 1.0_wp) <= 1.0_wp)
      call check_true ('treefall: a shorter cohort survives', survives_treefall(9.99_wp, 10.0_wp))
      call check_true ('treefall: one at the height dies', .not. survives_treefall(10.0_wp, 10.0_wp))
   end subroutine test_decisions

end program test_demography
