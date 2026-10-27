! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_canopy_aerodynamics -- the canopy's turbulent exchange (ed_canopy_aerodynamics;
! docs/science/canopy_aerodynamics.md §2-§4, eqs. (3), (7)-(11)).
!   1. NEUTRAL (mo_surface_layer): with the air and the canopy air alike, u* and the
!      scalar factor are the log law's, kappa u / ln(zldis/z0m) and kappa / ln(zldis/z0m),
!      and the same a hair either side of neutral.
!   2. STABILITY (mo_surface_layer): unstable gives a larger u* than neutral and stable a
!      smaller one; a stable night over a 30 m forest keeps zeta within zeta_max_stable
!      and its coupling rises with the wind; under a weak wind and hot canopy air the
!      sweeps settle on the root a 200-sweep solve finds.
!   3. CANOPY-TOP WIND (canopy_aerodynamics, eq. (8)): neutral, it is the log law's share
!      of the wind at zref, ln((h - d)/z0) / ln((zref - d)/z0).
!   4. BOUNDARY LAYER (canopy_aerodynamics, eqs. (10)-(11)): a leaf's conductance is above
!      its floor and rises with the wind; an isothermal 5 cm leaf at 20 C has Monteith &
!      Unsworth's laminar flat-plate value with the air's tabulated viscosity and
!      diffusivity.
!   5. TWO COHORTS (canopy_aerodynamics, eq. (9)): u* and the ground conductance are
!      positive, a closed canopy couples the ground less than bare soil, and the wind
!      falls top to bottom by each closed crown's extinction.
!   6. OPEN CANOPY (canopy_aerodynamics, eq. (7)): the ground conductance is the bare
!      value at W = 1, the dense-canopy value at W = 0, and rises monotonically between;
!      W = exp(-(LAI + WAI)); a bare patch couples as CLM5 couples bare soil.
!   7. TIED COHORTS (canopy_aerodynamics, eq. (9)): a crown split into two cohorts of half
!      the leaf area meets the wind the whole crown does, and so does the cohort below.
!   8. CALLER ORDER (canopy_aerodynamics): cohorts handed over tallest first, with the
!      order the gather forms, the tall one meets the top's wind and the larger boundary
!      layer; the wood's boundary layer follows the wood's temperature, not the leaves'.
!   9. THE REFERENCE LEVEL (canopy_aerodynamics, §2): canopy air at the temperature and
!      humidity of the air at the canopy-air top is neutral, with no temperature scale,
!      whatever the reference height.
!  10. SNOW BURIAL (canopy_burial_fraction, canopy_roughness, eq. (1)): no snow buries
!      nothing and leaves the canopy's own roughness and displacement; the buried share
!      rises with the pack's depth, continuously, to ED2's cap for a pack as deep as the
!      canopy or deeper; as burial goes to 1 the roughness and displacement go to the
!      snow surface's, continuously; a canopy below the snow's roughness takes the snow's.
!  11. CANOPY-AIR CONDUCTANCES (canopy_air_conductances, cas_atm_conductances, eq. (6)):
!      neutral, the log law's rho kappa^2 u / ln^2((zref - d)/z0) for enthalpy and vapour,
!      and the dry air's moles for CO2; more wind or a rougher canopy conducts more, a
!      buried canopy less, stable canopy air less and unstable more; at the patch kernel's
!      own state the re-solve gives what cas_atm_conductances makes of the kernel's u*.
!---~---
program test_canopy_aerodynamics
   use ed_test_assert,           only : check, check_true, test_report
   use ed_kinds,                 only : wp, ik
   use ed_constants,             only : von_karman, mmdry
   use ed_biophysics_params,     only : aero_params_t
   use ed_therm_lib,             only : cas_molar_density
   use ed_canopy_aerodynamics,   only : aero_env_t, aero_args_t, aero_out_t, ensure_aero_out_capacity,          &
                                        canopy_aerodynamics, mo_surface_layer, open_canopy_fraction,             &
                                        set_aero_env_atm, set_aero_env_canopy, canopy_burial_fraction,          &
                                        canopy_roughness, canopy_air_conductances, cas_atm_conductances
   use ed_numerics,              only : ascending_order
   implicit none

   call test_neutral()
   call test_stability_order()
   call test_canopy_top_wind()
   call test_boundary_layer()
   call test_two_cohorts()
   call test_open_canopy()
   call test_tied_layer()
   call test_caller_order()
   call test_reference_level()
   call test_snow_burial()
   call test_cas_conductances()

   call test_report('test_canopy_aerodynamics')

contains

   !--- 1. Neutral log-law recovery.
   subroutine test_neutral()
      type(aero_params_t) :: cfg
      real(wp) :: ustar, heat_transfer_factor, zeta, rib, obu, uref, z0m, zldis, ln, us_side
      uref  = 5.0_wp
      z0m   = 0.1_wp
      zldis = 10.0_wp                                       ! displace = 0 => zldis = zref
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp, 0.010_wp, ustar, heat_transfer_factor, zeta, rib, obu)
      ln = log(zldis / z0m)
      call check('neutral ustar', ustar, von_karman * uref / ln, 1.0e-5_wp)
      call check('neutral heat_transfer_factor', heat_transfer_factor, von_karman / ln, 1.0e-6_wp)
      call check('neutral zeta ~ 0', zeta, 0.0_wp, 1.0e-6_wp)
      !---~---
      !   A hair either side of neutral gives the neutral answer: the convective velocity
      ! fades with the buoyancy flux, so the wind has no jump there.
      !---~---
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp + 1.0e-9_wp, 0.010_wp, us_side, heat_transfer_factor, zeta, rib, obu)
      call check('ustar a hair unstable of neutral', us_side, ustar, 1.0e-6_wp)
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp - 1.0e-9_wp, 0.010_wp, us_side, heat_transfer_factor, zeta, rib, obu)
      call check('ustar a hair stable of neutral', us_side, ustar, 1.0e-6_wp)
   end subroutine test_neutral

   !--- 2. Stability ordering: unstable > neutral > stable ustar.
   subroutine test_stability_order()
      type(aero_params_t) :: cfg
      real(wp) :: us_u, us_n, us_s, t1, z, r, o, g_1
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 298.15_wp, 0.010_wp,          &
                            298.15_wp, 0.010_wp, us_n, t1, z, r, o)                 ! neutral
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 296.15_wp, 0.010_wp,          &
                            300.15_wp, 0.010_wp, us_u, t1, z, r, o)                 ! CAS warmer -> unstable
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 300.15_wp, 0.010_wp,          &
                            296.15_wp, 0.010_wp, us_s, t1, z, r, o)                 ! air warmer -> stable
      call check_true('unstable ustar > neutral', us_u > us_n, us_u - us_n)
      call check_true('stable ustar < neutral',   us_s < us_n, us_n - us_s)
      !--- Canopy air 3 K colder than the air at zref over a 30 m forest (d 18.9 m, z0 3.9 m).
      call mo_surface_layer(cfg, 1.0_wp, 40.0_wp, 18.9_wp, 3.9_wp, 290.15_wp, 0.008_wp,         &
                            287.15_wp, 0.008_wp, us_s, t1, z, r, o)
      g_1 = us_s * t1
      call check_true('stable sweeps within zeta_max_stable', z <= cfg%zeta_max_stable, z)
      call mo_surface_layer(cfg, 2.0_wp, 40.0_wp, 18.9_wp, 3.9_wp, 290.15_wp, 0.008_wp,         &
                            287.15_wp, 0.008_wp, us_s, t1, z, r, o)
      call check_true('stable coupling rises with wind', us_s * t1 > 1.5_wp * g_1, us_s * t1 - g_1)
      !---~---
      !   Free convection: 0.5 m/s over a 30 m forest, the canopy air 5 K warmer, where four fixed
      ! sweeps missed the coupling by 19 %.
      !---~---
      call mo_surface_layer(cfg, 0.5_wp, 35.0_wp, 18.9_wp, 3.9_wp, 298.15_wp, 0.010_wp,         &
                            303.15_wp, 0.010_wp, us_s, t1, z, r, o)
      g_1 = us_s * t1
      cfg%n_iter_mo = 200_ik
      call mo_surface_layer(cfg, 0.5_wp, 35.0_wp, 18.9_wp, 3.9_wp, 298.15_wp, 0.010_wp,         &
                            303.15_wp, 0.010_wp, us_s, t1, z, r, o)
      call check('free convection: the sweeps settle on the root', g_1, us_s * t1, 1.0e-3_wp * us_s * t1)
   end subroutine test_stability_order

   !---~---
   !   A neutral patch at 20 C: the air and the canopy air alike, and the ground at their
   ! temperature, under a wind u_ref at zref = 30 m.
   !---~---
   function neutral_env(u_ref) result(env)
      real(wp), intent(in) :: u_ref   !--- [m/s]
      type(aero_env_t)     :: env
      env%u_ref     = u_ref
      env%zref      = 30.0_wp
      env%theta_atm = 293.15_wp
      env%can_temp  = 293.15_wp
      env%t_ground  = 293.15_wp
      env%shv_atm   = 0.010_wp
      env%can_shv   = 0.010_wp
   end function neutral_env

   !---~---
   !   One cohort of 5 cm leaves at the top of a 20 m canopy, at the canopy air's temperature,
   ! under the neutral wind u_ref.
   !---~---
   subroutine one_leaf_cohort(u_ref, out)
      real(wp),         intent(in)    :: u_ref   !--- [m/s]
      type(aero_out_t), intent(inout) :: out
      type(aero_params_t) :: cfg
      type(aero_args_t)   :: geom
      integer(ik), parameter :: ORDER(1) = [1_ik]
      geom%veg_height   = 20.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call ensure_aero_out_capacity(out, 1_ik)
      call canopy_aerodynamics(cfg, neutral_env(u_ref), geom, 1_ik, ORDER, [20.0_wp], [0.5_wp], [293.15_wp],  &
                               [293.15_wp], [0.05_wp], [0.02_wp], out)
   end subroutine one_leaf_cohort

   !---~---
   !   3. The canopy-top wind (eq. (8)): neutral, u* and the profile are the log law's, so the
   ! wind at the canopy top is the wind at zref scaled by ln((h - d)/z0) / ln((zref - d)/z0).
   !---~---
   subroutine test_canopy_top_wind()
      type(aero_params_t) :: cfg
      type(aero_out_t)    :: out
      real(wp), parameter :: U_REF = 3.0_wp, H = 20.0_wp, ZREF = 30.0_wp
      real(wp) :: z0, d, uh_log
      call one_leaf_cohort(U_REF, out)
      z0 = cfg%z0m_ratio * H
      d  = cfg%d_ratio * H
      uh_log = U_REF * log((H - d) / z0) / log((ZREF - d) / z0)
      call check('neutral canopy-top wind is the log law''s', out%uh, uh_log, 1.0e-12_wp * uh_log)
      call check_true('the canopy-top wind is positive and below the wind at zref',                 &
                      out%uh > 0.0_wp .and. out%uh < U_REF, out%uh)
   end subroutine test_canopy_top_wind

   !---~---
   !   4. The leaf boundary layer (eqs. (10)-(11)), through the kernel. An isothermal leaf has no
   ! free convection, so at the low Reynolds number of a 5 cm leaf in ~1.5 m/s its conductance is
   ! the laminar flat plate's, Nu = 0.60 Re^0.5 (Monteith & Unsworth 2013, Principles of
   ! Environmental Physics, ch. 10), g = Nu kappa_T / l, with the air's viscosity 1.51e-5 and
   ! thermal diffusivity 2.15e-5 m2/s at 20 C (their Appendix A.3). Within 1 %: the linear
   ! fits of the two air properties sit within 0.4 % of the table.
   !---~---
   subroutine test_boundary_layer()
      type(aero_params_t) :: cfg
      type(aero_out_t)    :: lo, hi, mid
      real(wp), parameter :: NU_20 = 1.51e-5_wp, KAPPA_20 = 2.15e-5_wp, WIDTH = 0.05_wp   ! [m2/s], [m2/s], [m]
      real(wp) :: g_plate
      call one_leaf_cohort(1.0_wp, lo)
      call one_leaf_cohort(6.0_wp, hi)
      call check_true('gbh positive above floor', lo%leaf_gbh(1) > cfg%gbhmos_min, lo%leaf_gbh(1))
      call check_true('gbh increases with wind', hi%leaf_gbh(1) > lo%leaf_gbh(1), hi%leaf_gbh(1) - lo%leaf_gbh(1))
      call one_leaf_cohort(3.0_wp, mid)
      g_plate = 0.60_wp * sqrt(mid%wind(1) * WIDTH / NU_20) * KAPPA_20 / WIDTH
      call check('an isothermal 5 cm leaf: the laminar flat plate at 20 C', mid%leaf_gbh(1), g_plate,         &
                 1.0e-2_wp * g_plate)
   end subroutine test_boundary_layer

   !--- 5. Two cohorts under the default forcing.
   subroutine test_two_cohorts()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik), parameter :: N = 2_ik
      integer(ik), parameter :: ORDER(N) = [1_ik, 2_ik]   ! the arrays are bottom to top already
      real(wp), parameter :: HEIGHT(N) = [8.0_wp, 18.0_wp], LAI(N) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: TL(N) = [299.0_wp, 300.0_wp], TW(N) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: LW(N) = [0.04_wp, 0.04_wp], BD(N) = [0.02_wp, 0.02_wp]
      geom%veg_height   = 18.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp

      call ensure_aero_out_capacity(out, N)
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)

      call check_true('ustar > 0', out%ustar > 0.0_wp, out%ustar)
      call check_true('ggnet > 0', out%ggnet > 0.0_wp, out%ggnet)
      call check_true('ggnet <= ggbare (closed canopy)', out%ggnet <= out%ggbare + 1.0e-12_wp, out%ggbare - out%ggnet)
      call check_true('wind decreases top->bottom', out%wind(2) > out%wind(1), out%wind(2) - out%wind(1))
      !---~---
      !   Each closed crown takes the wind down by exp(-LAI/4) to its mid-depth and exp(-LAI/2)
      ! to its base (docs/science/canopy_aerodynamics.md eq. (9)).
      !---~---
      call check('top cohort: mid-crown wind = uh exp(-LAI/4)', out%wind(2),                        &
                 max(cfg%ugbmin, out%uh * exp(-0.25_wp * 2.0_wp)), 1.0e-12_wp)
      call check('bottom cohort: under the top crown, then to its own mid-depth', out%wind(1),       &
                 max(cfg%ugbmin, out%uh * exp(-0.50_wp * 2.0_wp) * exp(-0.25_wp * 1.5_wp)), 1.0e-12_wp)
   end subroutine test_two_cohorts

   !--- 6. The open-canopy blend of the ground conductance (eq. (7)).
   subroutine test_open_canopy()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik), parameter :: N = 2_ik, NW = 11_ik
      integer(ik), parameter :: ORDER(N) = [1_ik, 2_ik]   ! the arrays are bottom to top already
      real(wp), parameter :: HEIGHT(N) = [8.0_wp, 18.0_wp], LAI(N) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: TL(N) = [299.0_wp, 300.0_wp], TW(N) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: LW(N) = [0.04_wp, 0.04_wp], BD(N) = [0.02_wp, 0.02_wp]
      real(wp) :: gg(NW), r_edmore, r_clm, nu
      integer(ik) :: k
      geom%veg_height = 18.0_wp
      geom%snowfac    = 0.0_wp
      env%t_ground = env%can_temp + 2.0_wp           ! a warm ground: the dense term's stability is 0
      call ensure_aero_out_capacity(out, N)

      !--- The open fraction: 1 on bare ground, exp(-PAI) under plants.
      call check('a patch with no leaf or wood area is fully open', open_canopy_fraction(0.0_wp, 0.0_wp), &
                 1.0_wp, 0.0_wp)
      call check('W = exp(-(LAI + WAI))', open_canopy_fraction(5.0_wp, 1.0_wp), exp(-6.0_wp), 1.0e-15_wp)

      !--- The two limits, exactly, and a monotone rise between them.
      do k = 1_ik, NW
         geom%opencan_frac = real(k - 1_ik, wp) / real(NW - 1_ik, wp)
         call canopy_aerodynamics(cfg, env, geom, N, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)
         gg(k) = out%ggnet
         if (k == 1_ik) call check('W = 0: the dense-canopy conductance', out%ggnet, out%ggveg, 0.0_wp)
         if (k == NW)   call check('W = 1: the bare-ground conductance', out%ggnet, out%ggbare, 0.0_wp)
      end do
      call check_true('the bare-ground conductance exceeds the dense one', out%ggbare > out%ggveg,  &
                      out%ggbare - out%ggveg)
      call check_true('the ground conductance rises monotonically with W', all(gg(2:) > gg(:NW - 1)), &
                      minval(gg(2:) - gg(:NW - 1)))
      !---~---
      !   At LAI 6 (W = 0.0025) the blend is the linear weighting of the two terms, and the
      ! ground stays far closer to the dense value than to the bare one.
      !---~---
      geom%opencan_frac = open_canopy_fraction(6.0_wp, 0.0_wp)
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)
      call check('LAI 6: W ggbare + (1 - W) ggveg', out%ggnet,                                       &
                 geom%opencan_frac * out%ggbare + (1.0_wp - geom%opencan_frac) * out%ggveg, 1.0e-15_wp)
      call check_true('LAI 6: the ground couples nearly as under the dense form alone',                &
                      out%ggnet - out%ggveg < 0.01_wp * (out%ggbare - out%ggveg), out%ggnet / out%ggveg)
      !---~---
      !   A bare patch at the soil's roughness, neutral: 1/(u* temp1) + 1/ggbare is CLM5's
      ! (ln(zldis/z0m) + 0.13 (z0m u*/nu)^0.45) / (kappa u*), the resistance from the soil to zref.
      !---~---
      geom%veg_height   = cfg%soil_rough / cfg%z0m_ratio
      geom%opencan_frac = 1.0_wp
      env%can_temp  = 298.15_wp
      env%theta_atm = 298.15_wp
      env%t_ground  = 298.15_wp
      env%can_shv   = env%shv_atm
      call canopy_aerodynamics(cfg, env, geom, 0_ik, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)
      r_edmore = 1.0_wp / (out%ustar * out%heat_transfer_factor) + 1.0_wp / out%ggnet
      nu     = cfg%kin_visc0 * (1.0_wp + cfg%dkin_visc * (298.15_wp - cfg%t_ref_air))
      r_clm  = (log((env%zref - out%displace) / out%rough)                                        &
                + 0.13_wp * (out%rough * out%ustar / nu) ** 0.45_wp) / (von_karman * out%ustar)
      call check('a bare patch is coupled as CLM5 couples bare soil', r_edmore, r_clm, 1.0e-10_wp * r_clm)
   end subroutine test_open_canopy

   !---~---
   !   7. Tied cohorts share one wind layer (docs/science/canopy_aerodynamics.md eq. (9)), as they
   ! share one radiation layer: test_two_cohorts' stand with its top crown split into two cohorts
   ! of half the leaf area each must give every cohort the wind the unsplit stand gives it.
   !---~---
   subroutine test_tied_layer()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: one, split
      real(wp), parameter :: H2(2) = [8.0_wp, 18.0_wp], LAI2(2) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: TL2(2) = [299.0_wp, 300.0_wp], TW2(2) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: LW2(2) = [0.04_wp, 0.04_wp], BD2(2) = [0.02_wp, 0.02_wp]
      real(wp), parameter :: H3(3) = [8.0_wp, 18.0_wp, 18.0_wp], LAI3(3) = [1.5_wp, 1.0_wp, 1.0_wp]
      real(wp), parameter :: TL3(3) = [299.0_wp, 300.0_wp, 300.0_wp], TW3(3) = [298.5_wp, 299.0_wp, 299.0_wp]
      real(wp), parameter :: LW3(3) = [0.04_wp, 0.04_wp, 0.04_wp], BD3(3) = [0.02_wp, 0.02_wp, 0.02_wp]
      integer(ik), parameter :: ORD2(2) = [1_ik, 2_ik], ORD3(3) = [1_ik, 2_ik, 3_ik]   ! bottom to top
      geom%veg_height   = 18.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call ensure_aero_out_capacity(one, 2_ik)
      call ensure_aero_out_capacity(split, 3_ik)
      call canopy_aerodynamics(cfg, env, geom, 2_ik, ORD2, H2, LAI2, TL2, TW2, LW2, BD2, one)
      call canopy_aerodynamics(cfg, env, geom, 3_ik, ORD3, H3, LAI3, TL3, TW3, LW3, BD3, split)
      call check('tied: both halves of the top crown meet its wind', split%wind(2), one%wind(2), 1.0e-14_wp)
      call check('tied: the two halves meet the same wind', split%wind(3), split%wind(2), 0.0_wp)
      call check('tied: the cohort below meets the wind it met under the whole crown', split%wind(1),     &
                 one%wind(1), 1.0e-14_wp)
      call check('tied: the same leaf boundary layer', split%leaf_gbh(3), one%leaf_gbh(2), 1.0e-14_wp)
   end subroutine test_tied_layer

   !---~---
   !   8. The caller's cohort order. canopy_aerodynamics takes the cohorts bottom to top through
   ! the order it is handed, from arrays the gather leaves tallest first, so the tall cohort
   ! (index 1) meets the canopy top's wind and the larger boundary layer and the short one the
   ! attenuated wind; handing the gather order straight to the kernel inverts this. The wood's
   ! boundary layer follows the wood's own temperature (free convection is driven by the
   ! surface-air difference), which catches one temperature array passed for both, letting a
   ! warm leaf raise the wood's conductance.
   !---~---
   subroutine test_caller_order()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out, ref, warm_leaf, warm_wood
      integer(ik), parameter :: N = 2_ik
      real(wp), parameter :: HEIGHT(N) = [18.0_wp, 6.0_wp], LAI(N) = [2.0_wp, 2.0_wp]   ! tallest first
      real(wp), parameter :: LW(N) = [0.04_wp, 0.04_wp], BD(N) = [0.02_wp, 0.02_wp]
      integer(ik) :: order(N)
      real(wp)    :: t_air(N), t_warm(N)
      call ascending_order(HEIGHT, N, order)                 ! as the gather forms it
      env%u_ref = 3.0_wp
      geom%veg_height   = 18.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call ensure_aero_out_capacity(out, N)
      call canopy_aerodynamics(cfg, env, geom, N, order, HEIGHT, LAI, [300.0_wp, 300.0_wp], [300.0_wp, 300.0_wp], &
                               LW, BD, out)
      call check_true('caller order: the tall cohort (index 1) meets more wind', out%wind(1) > out%wind(2),        &
                      out%wind(1) - out%wind(2))
      call check_true('caller order: the tall cohort has the larger leaf boundary layer',                          &
                      out%leaf_gbw(1) > out%leaf_gbw(2), out%leaf_gbw(1) - out%leaf_gbw(2))
      t_air  = env%can_temp
      t_warm = env%can_temp + 8.0_wp
      call ensure_aero_out_capacity(ref, N)
      call ensure_aero_out_capacity(warm_leaf, N)
      call ensure_aero_out_capacity(warm_wood, N)
      call canopy_aerodynamics(cfg, env, geom, N, order, HEIGHT, LAI, t_air, t_air, LW, BD, ref)
      call canopy_aerodynamics(cfg, env, geom, N, order, HEIGHT, LAI, t_warm, t_air, LW, BD, warm_leaf)
      call canopy_aerodynamics(cfg, env, geom, N, order, HEIGHT, LAI, t_air, t_warm, LW, BD, warm_wood)
      call check_true('the wood''s boundary layer ignores the leaf temperature',                                   &
                      all(warm_leaf%wood_gbh(1:N) == ref%wood_gbh(1:N)), maxval(abs(warm_leaf%wood_gbh(1:N) - ref%wood_gbh(1:N))))
      call check_true('warm leaves raise the leaf boundary layer (free convection)',                               &
                      all(warm_leaf%leaf_gbh(1:N) > ref%leaf_gbh(1:N)), minval(warm_leaf%leaf_gbh(1:N) - ref%leaf_gbh(1:N)))
      call check_true('warm wood raises the wood boundary layer (free convection)',                                &
                      all(warm_wood%wood_gbh(1:N) > ref%wood_gbh(1:N)), minval(warm_wood%wood_gbh(1:N) - ref%wood_gbh(1:N)))
   end subroutine test_caller_order

   !---~---
   !   9. The reference level (docs/science/canopy_aerodynamics.md §2, "The reference level"): the
   ! surface layer compares the canopy air with the air at the canopy-air top, both referenced
   ! to that level, so canopy air at the temperature and humidity of the air there is neutral,
   ! with no temperature scale, whatever zref is. Referencing only the air to the ground would
   ! solve it as stable, tstar = heat_transfer_factor*(g/cp)*zref.
   !---~---
   subroutine test_reference_level()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik) :: order(1)
      call ascending_order([38.0_wp], 1_ik, order)
      geom%veg_height   = 38.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call set_aero_env_atm(env, 3.0_wp, 43.0_wp, 300.0_wp, 0.015_wp)
      call set_aero_env_canopy(env, 300.0_wp, 0.015_wp, 300.0_wp)
      call ensure_aero_out_capacity(out, 1_ik)
      call canopy_aerodynamics(cfg, env, geom, 1_ik, order, [38.0_wp], [5.0_wp], [300.0_wp], [300.0_wp], [0.04_wp], &
                               [0.02_wp], out)
      call check('reference level: canopy air like the air above has no temperature scale', out%tstar, 0.0_wp, 0.0_wp)
      call check('reference level: canopy air like the air above is neutral', out%zeta, 0.0_wp, 0.0_wp)
   end subroutine test_reference_level

   !---~---
   !   10. Snow burial (docs/science/canopy_aerodynamics.md eq. (1)): a pack buries the share of
   ! the canopy's height its depth reaches, capped at ED2's 0.99 so a buried canopy keeps a
   ! trace of its own roughness, and the roughness and displacement blend from the canopy's
   ! own (z0 = 0.13 h, d = 0.63 h) toward the snow surface's (snow_rough, no displacement). A
   ! 2 m shrub under a pack swept from bare ground to 1.5 times its height.
   !---~---
   subroutine test_snow_burial()
      type(aero_params_t) :: cfg
      type(aero_args_t)   :: geom
      integer(ik), parameter :: NS = 61_ik
      real(wp), parameter :: H = 2.0_wp, CAP = 0.99_wp      ! [m], [-] ED2's cap on the burial fraction
      real(wp) :: depth(NS), f(NS), z0(NS), d(NS), z0_own, d_own, z0_snow, d_snow, step
      integer(ik) :: k
      geom%veg_height   = H
      geom%opencan_frac = 0.0_wp
      step = 1.5_wp * H / real(NS - 1_ik, wp)
      do k = 1_ik, NS
         depth(k) = step * real(k - 1_ik, wp)
         f(k) = canopy_burial_fraction(depth(k), H)
         geom%snowfac = f(k)
         call canopy_roughness(cfg, geom, z0(k), d(k))
      end do
      z0_own = cfg%z0m_ratio * H
      d_own  = cfg%d_ratio * H
      !--- No snow: nothing buried, and the canopy's own roughness and displacement.
      call check('no snow buries nothing', f(1), 0.0_wp, 0.0_wp)
      call check('no snow: the canopy''s own roughness, 0.13 h', z0(1), z0_own, 1.0e-15_wp * z0_own)
      call check('no snow: the canopy''s own displacement, 0.63 h', d(1), d_own, 1.0e-15_wp * d_own)
      !--- The buried share rises with the pack, as the share of the height it reaches, to the cap.
      call check_true('burial rises monotonically with the pack''s depth', all(f(2:) >= f(:NS - 1_ik)),       &
                      minval(f(2:) - f(:NS - 1_ik)))
      call check('a pack half the canopy''s height buries half of it', canopy_burial_fraction(0.5_wp * H, H),  &
                 0.5_wp, 1.0e-15_wp)
      call check_true('a pack as deep as the canopy, or deeper, buries it to the cap',                         &
                      all(abs(pack(f, depth >= H) - CAP) <= 1.0e-15_wp), maxval(pack(f, depth >= H)))
      call check_true('the burial is continuous: a step in depth moves it by at most the step over h',      &
                      all(abs(f(2:) - f(:NS - 1_ik)) <= step / H * (1.0_wp + 1.0e-12_wp)),                     &
                      maxval(abs(f(2:) - f(:NS - 1_ik))))
      !--- The roughness and displacement fall monotonically and continuously with the burial.
      call check_true('a deeper pack never roughens the canopy', all(z0(2:) <= z0(:NS - 1_ik)),               &
                      maxval(z0(2:) - z0(:NS - 1_ik)))
      call check_true('a deeper pack never raises the displacement', all(d(2:) <= d(:NS - 1_ik)),             &
                      maxval(d(2:) - d(:NS - 1_ik)))
      call check_true('the roughness is continuous in the pack''s depth',                                      &
                      all(abs(z0(2:) - z0(:NS - 1_ik)) <= (z0_own - cfg%snow_rough) * step / H                 &
                          * (1.0_wp + 1.0e-12_wp)), maxval(abs(z0(2:) - z0(:NS - 1_ik))))
      call check_true('under the deepest pack the canopy keeps only a trace of its own roughness',            &
                      z0(NS) > cfg%snow_rough .and. z0(NS) - cfg%snow_rough <= (1.0_wp - CAP) * z0_own        &
                      * (1.0_wp + 1.0e-12_wp) .and. d(NS) <= (1.0_wp - CAP) * d_own * (1.0_wp + 1.0e-12_wp),   &
                      z0(NS) - cfg%snow_rough)
      !--- Burial 1 is the snow surface itself, and burial approaching 1 approaches it.
      geom%snowfac = 1.0_wp
      call canopy_roughness(cfg, geom, z0_snow, d_snow)
      call check('fully buried: the snow surface''s roughness', z0_snow, cfg%snow_rough, 1.0e-15_wp)
      call check('fully buried: no displacement', d_snow, 0.0_wp, 1.0e-15_wp)
      geom%snowfac = 1.0_wp - 1.0e-6_wp
      call canopy_roughness(cfg, geom, z0(1), d(1))
      call check_true('nearly buried: within a millionth of the canopy''s own of the snow surface',            &
                      abs(z0(1) - z0_snow) <= 1.0e-6_wp * z0_own * (1.0_wp + 1.0e-9_wp)                        &
                      .and. abs(d(1) - d_snow) <= 1.0e-6_wp * d_own * (1.0_wp + 1.0e-9_wp), z0(1) - z0_snow)
      !--- A canopy below the snow surface's roughness (5 mm of grass) takes the snow's.
      geom%veg_height = 0.005_wp
      geom%snowfac    = 0.0_wp
      call canopy_roughness(cfg, geom, z0(1), d(1))
      call check('a canopy below the snow''s roughness takes the snow''s', z0(1), cfg%snow_rough, 0.0_wp)
   end subroutine test_snow_burial

   !---~---
   !   11. The canopy-air <-> atmosphere conductances (docs/science/canopy_aerodynamics.md eq. (6)),
   ! re-solved at a canopy-air state (canopy_air_conductances) from the surface layer and the
   ! roughness, under a wind u at zref = 30 m (neutral_env). Neutral, u* = kappa u / ln and the
   ! profile factor is kappa / ln, ln = ln((zref - d)/z0), so g = rho kappa^2 u / ln^2; vapour
   ! takes the heat factor (z0q = z0h) and CO2 counts the dry air's moles, rho (1 - q) / M_dry.
   !---~---
   subroutine test_cas_conductances()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      real(wp), parameter :: RHO = 1.2_wp, H = 20.0_wp, ZREF = 30.0_wp, U = 3.0_wp
      real(wp), parameter :: PACK_DEPTH(4) = [0.0_wp, 0.5_wp, 1.0_wp, 2.0_wp]   ! [m]
      real(wp) :: gh, gv, gc, g_log, ln, g(4), gk_h, gk_v, gk_c
      integer(ik) :: k, order(1)
      env = neutral_env(U)
      geom%veg_height   = H
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call canopy_air_conductances(cfg, env, geom, RHO, env%theta_atm, env%shv_atm, gh, gv, gc)
      ln    = log((ZREF - cfg%d_ratio * H) / (cfg%z0m_ratio * H))
      g_log = RHO * von_karman ** 2 * U / ln ** 2
      call check('neutral: g = rho kappa^2 u / ln^2((zref - d)/z0)', gh, g_log, 1.0e-12_wp * g_log)
      call check('vapour takes the enthalpy''s conductance (z0q = z0h)', gv, gh, 0.0_wp)
      call check('CO2: the same exchange counted in the dry air''s moles, (1 - q)/M_dry per kg', gc,           &
                 gh * (1.0_wp - env%shv_atm) / mmdry, 1.0e-12_wp * gc)
      !--- More wind conducts more: neutral, in proportion to it.
      do k = 1_ik, 4_ik
         env = neutral_env(U * 2.0_wp ** (k - 1_ik))
         call canopy_air_conductances(cfg, env, geom, RHO, env%theta_atm, env%shv_atm, g(k), gv, gc)
      end do
      call check_true('a stronger wind conducts more', all(g(2:) > g(:3)), minval(g(2:) - g(:3)))
      call check('neutral, the conductance is linear in the wind', g(4), 8.0_wp * g(1), 1.0e-12_wp * g(4))
      !--- A rougher canopy (taller, under the same zref) conducts more.
      env = neutral_env(U)
      do k = 1_ik, 4_ik
         geom%veg_height = 2.5_wp * 2.0_wp ** (k - 1_ik)
         call canopy_air_conductances(cfg, env, geom, RHO, env%theta_atm, env%shv_atm, g(k), gv, gc)
      end do
      call check_true('a rougher canopy conducts more (h 2.5, 5, 10, 20 m)', all(g(2:) > g(:3)),               &
                      minval(g(2:) - g(:3)))
      !--- A 2 m shrub under a deepening pack (0, 0.5, 1, 2 m) conducts less.
      geom%veg_height = 2.0_wp
      do k = 1_ik, 4_ik
         geom%snowfac = canopy_burial_fraction(PACK_DEPTH(k), geom%veg_height)
         call canopy_air_conductances(cfg, env, geom, RHO, env%theta_atm, env%shv_atm, g(k), gv, gc)
      end do
      call check_true('a buried canopy conducts less, the more so the deeper the pack', all(g(2:) < g(:3)),   &
                      maxval(g(2:) - g(:3)))
      !---~---
      !   Canopy air 3 K colder than the air above (stable) conducts less than neutral, 3 K warmer
      ! (unstable) more.
      !---~---
      geom%veg_height = H
      geom%snowfac    = 0.0_wp
      call canopy_air_conductances(cfg, env, geom, RHO, env%theta_atm - 3.0_wp, env%shv_atm, g(1), gv, gc)
      call canopy_air_conductances(cfg, env, geom, RHO, env%theta_atm + 3.0_wp, env%shv_atm, g(2), gv, gc)
      call check_true('stable canopy air conducts less than neutral', g(1) < gh, gh - g(1))
      call check_true('unstable canopy air conducts more than neutral', g(2) > gh, g(2) - gh)
      !---~---
      !   The patch kernel's u* and profile factors at its own, unstable canopy air give, through
      ! cas_atm_conductances, the conductances the re-solve finds at that state: the frozen
      ! record's and the march's are one exchange.
      !---~---
      env%can_temp = env%theta_atm + 3.0_wp
      env%can_shv  = 0.012_wp
      call ascending_order([H], 1_ik, order)
      call ensure_aero_out_capacity(out, 1_ik)
      call canopy_aerodynamics(cfg, env, geom, 1_ik, order, [H], [4.0_wp], [env%can_temp], [env%can_temp],       &
                               [0.05_wp], [0.02_wp], out)
      call cas_atm_conductances(RHO, cas_molar_density(RHO, env%can_shv), out%ustar, out%heat_transfer_factor,  &
                                out%vapour_transfer_factor, gk_h, gk_v, gk_c)
      call canopy_air_conductances(cfg, env, geom, RHO, env%can_temp, env%can_shv, gh, gv, gc)
      call check('the kernel''s u* through cas_atm_conductances: the re-solved enthalpy conductance', gk_h, gh,   &
                 1.0e-12_wp * gh)
      call check('... and the re-solved vapour and CO2 conductances', abs(gk_v - gv) / gv + abs(gk_c - gc) / gc,  &
                 0.0_wp, 1.0e-12_wp)
   end subroutine test_cas_conductances

end program test_canopy_aerodynamics
