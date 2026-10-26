! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_canopy_water -- the canopy's interception film (ed_canopy_water;                    !
! docs/science/vegetation_energy_dynamics.md, "The canopy's interception film").           !
!   1. INTERCEPTION (intercept_canopy_layer): a dry canopy fills below its capacity with   !
!      no drip and an exact water balance; a full one drips, holds its capacity and is     !
!      wholly wet; a film above its capacity, as leaves fall from under it, drips the      !
!      excess, and throughfall plus the change in storage is still the rain.               !
!==========================================================================================!
program test_canopy_water
   use ed_test_assert,    only : check, check_true, test_report
   use ed_kinds,          only : wp
   use ed_canopy_water,   only : intercept_canopy_layer
   implicit none

   call test_interception()

   call test_report('test_canopy_water')

contains

   !----- 1. One cohort of LAI 2 and SAI 0.5 under 0.1 mm/s of rain for 10 minutes, its film holding    !
   !      0.1 kg/m2 per unit of leaf and of wood area (0.25 kg/m2 in all), the interception fraction's  !
   !      extinction 0.5 and scale 1. ----------------------------------------------------------------!
   subroutine test_interception()
      real(wp) :: lw, tf, dr, sw, bal, rain, dt
      print '(a)', 'test_interception:'
      dt   = 600.0_wp
      rain = 1.0e-4_wp
      !----- Dry canopy fills below capacity: no drip, exact balance. -----!
      lw = 0.0_wp
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, dt, 0.1_wp, 0.1_wp, 0.5_wp, 1.0_wp,        &
                                  tf, dr, sw)
      bal = tf + lw / dt                                ! throughfall + storage rate
      call check('interception water balance', bal, rain, 1.0e-12_wp)
      call check_true('no drip below capacity', dr < 1.0e-30_wp, dr)
      call check_true('leaf_water within capacity', lw <= 0.1_wp * 2.5_wp + 1.0e-12_wp, lw)
      !----- Saturated canopy overflows: drip appears, storage capped. -----!
      lw = 0.25_wp                                       ! = 0.1*LAI + 0.1*SAI (full)
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, dt, 0.1_wp, 0.1_wp, 0.5_wp, 1.0_wp,        &
                                  tf, dr, sw)
      call check_true('drip when full', dr > 0.0_wp, dr)
      call check_true('capacity respected', lw <= 0.25_wp + 1.0e-12_wp, lw)
      call check_true('sigma_w saturates to 1', abs(sw - 1.0_wp) < 1.0e-9_wp, sw)
      !----- A film already ABOVE capacity (leaf area lost under a full film) drips the excess:       !
      !      storage capped, and throughfall + storage change = rain -- nothing discarded. ---------!
      lw = 0.30_wp                                       ! capacity is 0.25: 0.05 too much
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, dt, 0.1_wp, 0.1_wp, 0.5_wp, 1.0_wp,        &
                                  tf, dr, sw)
      call check_true('over capacity: storage capped', lw <= 0.25_wp + 1.0e-12_wp, lw)
      bal = tf + (lw - 0.30_wp) / dt
      call check('over capacity: the excess drips, the water balance holds', bal, rain, 1.0e-12_wp)
   end subroutine test_interception

end program test_canopy_water
