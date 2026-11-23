! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_ground_surface_lib -- the laws of the ground's surface layers (ed_ground_surface_lib).
!   1. POROUS PHASE (porous_layer_temp, porous_energy_state, porous_layer_energy,
!      porous_energy_at_temp): the all-liquid closed form; energy -> temperature -> energy
!      round-trips on the liquid line, on the freezing curve and on the plateau of water above
!      saturation at t_3ple; the apparent heat capacity against a central difference.
!   2. POROUS LAYER (porous_layer_depth, porous_layer_capacity, porous_layer_fill,
!      porous_layer_resistance, porous_layer_matter): depth and capacity from the matter; the
!      resistance's dry and saturated limits; a layer only from min_depth.
!   3. WATER LAYER, MASS (water_layer_free_water, water_layer_depth, matrix_depth_update,
!      combined_density): a pond is all free water, a layer under its hold has none, one over it
!      the liquid less h/(1-h) of its ice; rain fills the matrix's pores and adds no depth, a
!      pond is as deep as its water; snowfall adds its volume, ice that forms fills the pores
!      first, melt shrinks the matrix at its density; volumes add.
!   4. WATER LAYER, SURFACE (water_layer_cover, water_layer_albedo, water_layer_resistance): the
!      cover of 1 and 10 cm at 250 kg/m3, and of 1 and 5 cm of ponded water; the albedo's fresh
!      and aged ends; the resistance of the matrix and of the water standing above it.
!---~---
program test_ground_surface_lib
   use ed_test_assert,        only : check_true, check_close, test_report
   use ed_kinds,              only : wp
   use ed_constants,          only : rho_h2o, rho_ice, t_3ple, cp_liq, tsupercool_liq, k_water
   use ed_soil_hydr_lib,      only : soil_layer_curve_t, build_soil_layer_curve, SOIL_RETENTION_VG,            &
                                     soil_liquid_capacity
   use ed_ground_surface_lib, only : porous_layer_temp, porous_energy_state, porous_layer_energy,               &
                                     porous_energy_at_temp, porous_layer_depth, porous_layer_capacity,           &
                                     porous_layer_fill, porous_layer_resistance, porous_layer_matter,            &
                                     water_layer_free_water, water_layer_depth, matrix_depth_update,             &
                                     combined_density, water_layer_cover, water_layer_albedo,                    &
                                     water_layer_resistance
   implicit none

   call test_porous_phase()
   call test_porous_layer()
   call test_water_layer_mass()
   call test_water_layer_surface()

   call test_report('test_ground_surface_lib')

contains

   !--- 1. The porous medium's phase.
   subroutine test_porous_phase()
      type(soil_layer_curve_t) :: c
      real(wp), parameter :: CDRY = 1.2e6_wp, THETA = 0.35_wp
      real(wp) :: temps(4), e, t, fl, de, e_lo, e_hi, fl2, liq, dliq, wmass, h
      integer  :: i
      call build_soil_layer_curve(SOIL_RETENTION_VG, 0.4314_wp, 0.0981_wp, 0.7607_wp, 1.3554_wp, 4.54e-6_wp,   &
                                  0.05_wp, c)
      wmass = THETA * rho_h2o
      !--- The liquid line in closed form.
      e = porous_layer_energy(290.0_wp, THETA, CDRY, c)
      call porous_layer_temp(e, THETA, CDRY, c, t, fl)
      call check_close('the all-liquid line in closed form', t,                                             &
                       (e + wmass * cp_liq * tsupercool_liq) / (CDRY + wmass * cp_liq), 1.0e-9_wp)
      call check_close('all liquid above t_3ple', fl, 1.0_wp, 0.0_wp)
      !--- Round trips on both sides of the freezing point, the freezing curve's liquid below it.
      temps = [295.0_wp, t_3ple, 272.0_wp, 260.0_wp]
      do i = 1, size(temps)
         e = porous_layer_energy(temps(i), THETA, CDRY, c)
         call porous_energy_state(e, THETA, CDRY, c, t, fl, de)
         call check_close('energy -> temperature round-trips', t, temps(i), 1.0e-7_wp)
         if (temps(i) < t_3ple) then
            call soil_liquid_capacity(c, temps(i), liq, dliq)
            call check_close('below freezing the curve holds the liquid', fl * THETA, min(liq, THETA), 1.0e-9_wp)
         end if
      end do
      !--- Water above saturation freezes on a plateau at t_3ple.
      call porous_energy_at_temp(t_3ple, 0.5_wp, CDRY, c, e_lo, de, fl2)
      e_hi = porous_layer_energy(t_3ple + 1.0e-6_wp, 0.5_wp, CDRY, c)
      call porous_layer_temp(0.5_wp * (e_lo + e_hi), 0.5_wp, CDRY, c, t, fl)
      call check_close('water beyond saturation freezes at t_3ple', t, t_3ple, 1.0e-9_wp)
      call check_true('and is part liquid there', fl > fl2 .and. fl < 1.0_wp, fl)
      !--- The apparent heat capacity against a central difference on the freezing curve.
      h = 1.0e-4_wp
      call porous_energy_at_temp(265.0_wp, THETA, CDRY, c, e, de, fl)
      call check_close('the apparent heat capacity is dE/dT', de,                                           &
                       (porous_layer_energy(265.0_wp + h, THETA, CDRY, c)                                    &
                        - porous_layer_energy(265.0_wp - h, THETA, CDRY, c)) / (2.0_wp * h), 1.0e-4_wp * de)
   end subroutine test_porous_phase

   !--- 2. A porous layer: litter at 45 kg/m3, theta_sat 0.98, k 0.10 dry and 0.1294 saturated.
   subroutine test_porous_layer()
      real(wp), parameter :: RHO_B = 45.0_wp, THETA_SAT = 0.98_wp, K_DRY = 0.10_wp, K_WET = 0.1294_wp
      real(wp), parameter :: MIN_DEPTH = 1.0e-4_wp
      real(wp) :: matter, depth, cap
      matter = 0.9_wp                                              ! 2 cm of litter
      depth  = porous_layer_depth(matter, RHO_B)
      cap    = porous_layer_capacity(matter, RHO_B, THETA_SAT)
      call check_close('depth is matter over bulk density', depth, 0.02_wp, 1.0e-15_wp)
      call check_close('capacity is theta_sat over the depth', cap, THETA_SAT * rho_h2o * 0.02_wp, 1.0e-12_wp)
      call check_close('fill is the share of the capacity', porous_layer_fill(0.5_wp * cap, cap), 0.5_wp, 1.0e-15_wp)
      call check_close('fill stops at saturation', porous_layer_fill(2.0_wp * cap, cap), 1.0_wp, 0.0_wp)
      call check_close('a dry layer conducts at k_dry', porous_layer_resistance(depth, 0.0_wp, K_DRY, K_WET),     &
                       depth / K_DRY, 1.0e-15_wp)
      call check_close('a saturated one at k_wet', porous_layer_resistance(depth, 1.0_wp, K_DRY, K_WET),          &
                       depth / K_WET, 1.0e-15_wp)
      call check_close('below min_depth no layer', porous_layer_matter(0.5_wp * MIN_DEPTH * RHO_B, RHO_B,         &
                       MIN_DEPTH), 0.0_wp, 0.0_wp)
      call check_close('from min_depth all of it', porous_layer_matter(MIN_DEPTH * RHO_B, RHO_B, MIN_DEPTH),      &
                       MIN_DEPTH * RHO_B, 0.0_wp)
   end subroutine test_porous_layer

   !--- 3. The water layer's mass: holding at LEAF-3's 1:9, the matrix at fresh snow's 250 kg/m3.
   subroutine test_water_layer_mass()
      real(wp), parameter :: H = 0.1_wp, RHO_FRESH = 250.0_wp
      real(wp) :: w, ice, liq
      w = 2.0_wp
      call check_close('a pond is all free water', water_layer_free_water(w, 1.0_wp, H), w, 1.0e-15_wp)
      call check_close('a layer under its hold has none', water_layer_free_water(w, 0.05_wp, H), 0.0_wp, 0.0_wp)
      w   = 10.0_wp
      liq = 3.0_wp
      ice = w - liq
      call check_close('over the hold: the liquid less h/(1-h) of the ice', water_layer_free_water(w, liq / w, H), &
                       liq - H / (1.0_wp - H) * ice, 1.0e-12_wp)
      !--- Depth: 100 kg/m2 of snow is a 0.4 m matrix.
      call check_close('rain fills the pores and adds no depth', water_layer_depth(0.4_wp, 100.0_wp, 10.0_wp),     &
                       0.4_wp, 0.0_wp)
      call check_close('a pond is as deep as its water', water_layer_depth(0.0_wp, 0.0_wp, 50.0_wp), 0.05_wp,      &
                       1.0e-15_wp)
      call check_close('snowfall adds its volume', matrix_depth_update(0.0_wp, 0.0_wp, 25.0_wp, 25.0_wp,          &
                       RHO_FRESH), 0.1_wp, 1.0e-15_wp)
      call check_close('ice that forms fills the pores first', matrix_depth_update(0.4_wp, 100.0_wp, 110.0_wp,     &
                       0.0_wp, RHO_FRESH), 0.4_wp, 0.0_wp)
      call check_close('a pond that freezes is solid ice', matrix_depth_update(0.0_wp, 0.0_wp, 50.0_wp, 0.0_wp,   &
                       RHO_FRESH), 50.0_wp / rho_ice, 1.0e-15_wp)
      call check_close('melt shrinks the matrix at its density', matrix_depth_update(0.4_wp, 100.0_wp, 50.0_wp,   &
                       0.0_wp, RHO_FRESH), 0.2_wp, 1.0e-15_wp)
      call check_close('volumes add', combined_density(10.0_wp, RHO_FRESH, 10.0_wp, rho_ice),                     &
                       20.0_wp / (10.0_wp / RHO_FRESH + 10.0_wp / rho_ice), 1.0e-12_wp)
      call check_close('no mass takes the other''s density', combined_density(0.0_wp, rho_ice, 5.0_wp, RHO_FRESH),&
                       RHO_FRESH, 0.0_wp)
   end subroutine test_water_layer_mass

   !---~---
   !   4. The water layer's surface. The cover at the shipped scale factor 2.5, exponent 1 and
   ! the ground's roughness 0.01 m: tanh(0.16) and tanh(1.6) for 1 and 10 cm of fresh snow, and
   ! tanh(0.04) and tanh(0.2) for 1 and 5 cm of ponded water.
   !---~---
   subroutine test_water_layer_surface()
      real(wp), parameter :: A = 2.5_wp, M = 1.0_wp, Z0 = 0.01_wp
      call check_close('1 cm of fresh snow covers 0.16', water_layer_cover(0.01_wp, 250.0_wp, A, M, Z0),           &
                       tanh(0.16_wp), 1.0e-15_wp)
      call check_close('10 cm of fresh snow covers 0.92', water_layer_cover(0.10_wp, 250.0_wp, A, M, Z0),          &
                       tanh(1.6_wp), 1.0e-15_wp)
      call check_close('1 cm of ponded water covers 0.04', water_layer_cover(0.01_wp, 1000.0_wp, A, M, Z0),        &
                       tanh(0.04_wp), 1.0e-15_wp)
      call check_close('5 cm of ponded water covers 0.20', water_layer_cover(0.05_wp, 1000.0_wp, A, M, Z0),        &
                       tanh(0.2_wp), 1.0e-15_wp)
      call check_close('dry snow is fresh', water_layer_albedo(0.0_wp, 0.85_wp, 0.518_wp), 0.85_wp, 0.0_wp)
      call check_close('wet snow is aged', water_layer_albedo(1.0_wp, 0.85_wp, 0.518_wp), 0.518_wp, 0.0_wp)
      call check_close('the matrix alone', water_layer_resistance(0.4_wp, 0.4_wp, 0.15_wp), 0.4_wp / 0.15_wp,     &
                       0.0_wp)
      call check_close('and the water standing above it', water_layer_resistance(0.15_wp, 0.1_wp, 0.15_wp),       &
                       0.1_wp / 0.15_wp + 0.05_wp / k_water, 1.0e-12_wp)
   end subroutine test_water_layer_surface

end program test_ground_surface_lib
