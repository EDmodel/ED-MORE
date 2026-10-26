! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_optics_lib -- the canopy's optical kernels (meds_optics_lib;                        !
! docs/science/canopy_radiation_transfer.md §1-§2, eqs. (1)-(3)). A PFT's leaf-angle       !
! moments reach the Beta LIDF through meds_pft_derive's moments_to_beta_params, as at load.!
!   1. LEAF-ANGLE LIMITS (beta_lidf, leaf_bf): a uniform inclination gives bf = 0.5, a     !
!      planophile canopy bf near 1 and an erectophile one bf near 0.                       !
!   2. NORMALIZATION (beta_lidf): narrow distributions, whose density underflows, still    !
!      sum to one, and a zero spread at 45 deg puts the leaf area in its class.            !
!   3. DELTA LIDF (moments_to_beta_params): a zero spread gives a finite Beta (p, q).      !
!   4. G(MU) (gfun_direct, eq. (3)): the Ross G-function against its closed forms: leaves  !
!      at one inclination theta_L project mu cos(theta_L) while theta_L + theta_s <= 90    !
!      deg (G = mu for horizontal leaves); near-vertical leaves project (2/pi) sin theta_s;!
!      a spherical LIDF projects G = 1/2 toward every sun.                                 !
!==========================================================================================!
program test_optics_lib
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : pi
   use meds_optics_lib,  only : beta_lidf, leaf_bf, gfun_direct, N_LEAF_CLASS
   use meds_pft_derive,  only : moments_to_beta_params
   implicit none

   call test_leaf_angle_limits()
   call test_lidf_normalization()
   call test_lidf_delta()
   call test_g_function()

   call test_report('test_optics_lib')

contains

   !----- 1. bf = <cos^2 theta_leaf> (eq. (2)) at its limits. ---------------------------------!
   subroutine test_leaf_angle_limits()
      real(wp) :: p, q, bf
      real(wp) :: lidf(N_LEAF_CLASS)
      print '(a)', 'test_leaf_angle_limits:'
      !----- Uniform inclination (Beta(1,1)): <cos^2> = 0.5 exactly. ------------------------!
      lidf = beta_lidf(1.0_wp, 1.0_wp)
      bf   = leaf_bf(lidf)
      call check('bf(uniform LAD)', bf, 0.5_wp, 2.0e-2_wp)
      !----- Strongly planophile (small mean angle): bf -> ~1 (near-horizontal). ------------!
      call moments_to_beta_params(8.0_wp, 6.0_wp, p, q)
      lidf = beta_lidf(p, q)
      bf   = leaf_bf(lidf)
      call check_true('bf(planophile) > 0.90', bf > 0.90_wp, bf)
      !----- Strongly erectophile (near-vertical): bf -> ~0. --------------------------------!
      call moments_to_beta_params(82.0_wp, 6.0_wp, p, q)
      lidf = beta_lidf(p, q)
      bf   = leaf_bf(lidf)
      call check_true('bf(erectophile) < 0.10', bf < 0.10_wp, bf)
   end subroutine test_leaf_angle_limits

   !----- 2. Narrow distributions, whose unnormalized density underflows, still sum to one; a     !
   !      spread of 0 at 45 deg puts the leaf area in the 40-50 deg class. -----------------------!
   subroutine test_lidf_normalization()
      real(wp) :: p, q
      real(wp) :: lidf(N_LEAF_CLASS)
      print '(a)', 'test_lidf_normalization:'
      call moments_to_beta_params(64.0_wp, 3.0_wp, p, q)
      lidf = beta_lidf(p, q)
      call check('LIDF sums to one (64 +- 3 deg)', sum(lidf), 1.0_wp, 1.0e-12_wp)
      call moments_to_beta_params(45.0_wp, 0.0_wp, p, q)
      lidf = beta_lidf(p, q)
      call check('LIDF sums to one (45 +- 0 deg)', sum(lidf), 1.0_wp, 1.0e-12_wp)
      call check_true('45 +- 0 deg in its class', lidf(5) > 0.99_wp, lidf(5))
   end subroutine test_lidf_normalization

   !----- 3. A delta leaf-angle distribution (std_deg = 0) must give FINITE Beta (p, q): the       !
   !      variance floor keeps kappa finite, where mt(1-mt)/0 would give a NaN LIDF. --------------!
   subroutine test_lidf_delta()
      real(wp) :: p, q
      print '(a)', 'test_lidf_delta:'
      call moments_to_beta_params(45.0_wp, 0.0_wp, p, q)     ! mean=45deg -> mt=0.5; std=0 -> vt floored
      call check('delta LIDF p finite (variance floored)', p, 124999.5_wp, 1.0e-3_wp)
      call check('delta LIDF q finite (variance floored)', q, 124999.5_wp, 1.0e-3_wp)
   end subroutine test_lidf_delta

   !----- 4. The Ross G-function (eq. (3)) against its closed forms (Ross 1981, The Radiation      !
   !      Regime and Architecture of Plant Stands, ch. 3). With all the leaf area at one           !
   !      inclination theta_L and theta_L + theta_s <= 90 deg, no leaf is seen from its back and   !
   !      G = mu cos(theta_L) exactly: on the grid's flattest class (5 deg) that is the horizontal !
   !      limit G = mu, short of it by cos(5 deg). As theta_L -> 90 deg, G -> (2/pi) sin theta_s; !
   !      the steepest class (89 deg) lies within 1e-3 of it, an error of order cos^2(theta_L).   !
   !      A spherical LIDF, each class holding cos(lo) - cos(hi) of the leaf area, projects G = 1/2 !
   !      toward every sun; the classes' midpoint rule holds that within 2.5e-3. ----------------!
   subroutine test_g_function()
      !----- [deg] SCOPE's 13 inclination classes: upper bounds, 10 deg wide to 80, then 2 deg. ---!
      real(wp), parameter :: CLASS_TOP_DEG(N_LEAF_CLASS) =                                          &
         [ 10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp, 50.0_wp, 60.0_wp, 70.0_wp, 80.0_wp,                     &
           82.0_wp, 84.0_wp, 86.0_wp, 88.0_wp, 90.0_wp ]
      real(wp), parameter :: DEG = pi / 180.0_wp
      real(wp) :: lidf(N_LEAF_CLASS), mu, sin_s, g_vertical, worst_flat, worst_vert, worst_sph, lo
      integer(ik) :: i, k
      print '(a)', 'test_g_function:'

      !----- Horizontal leaves: the 5 deg class alone, suns from 84 deg to the zenith. -------------!
      lidf    = 0.0_wp
      lidf(1) = 1.0_wp
      worst_flat = 0.0_wp
      do i = 1_ik, 10_ik
         mu = 0.1_wp * real(i, wp)
         worst_flat = max(worst_flat, abs(gfun_direct(lidf, mu) - mu * cos(5.0_wp * DEG)))
      end do
      call check('leaves at 5 deg: G = mu cos(5 deg), the horizontal limit', worst_flat, 0.0_wp, 1.0e-14_wp)

      !----- Vertical leaves: the 89 deg class alone, suns from 84 deg to 26 deg. ------------------!
      lidf     = 0.0_wp
      lidf(N_LEAF_CLASS) = 1.0_wp
      worst_vert = 0.0_wp
      do i = 1_ik, 9_ik
         mu    = 0.1_wp * real(i, wp)
         sin_s = sqrt(1.0_wp - mu * mu)
         g_vertical = (2.0_wp / pi) * sin_s
         worst_vert = max(worst_vert, abs(gfun_direct(lidf, mu) - g_vertical) / g_vertical)
      end do
      call check('leaves at 89 deg: G = (2/pi) sin(theta_s), the vertical limit', worst_vert, 0.0_wp, 1.0e-3_wp)

      !----- A spherical LIDF: the same projection toward every sun. -------------------------------!
      lo = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         lidf(k) = cos(lo * DEG) - cos(CLASS_TOP_DEG(k) * DEG)
         lo = CLASS_TOP_DEG(k)
      end do
      call check('the spherical LIDF holds all the leaf area', sum(lidf), 1.0_wp, 1.0e-14_wp)
      worst_sph = 0.0_wp
      do i = 1_ik, 20_ik
         mu = 0.05_wp * real(i, wp)
         worst_sph = max(worst_sph, abs(gfun_direct(lidf, mu) - 0.5_wp))
      end do
      call check('a spherical LIDF: G = 1/2 toward every sun', worst_sph, 0.0_wp, 2.5e-3_wp)
   end subroutine test_g_function

end program test_optics_lib
