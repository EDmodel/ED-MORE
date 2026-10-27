! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_optics_lib -- the optical kernels of the canopy radiation: the leaf-inclination
! distribution (LIDF), a Beta distribution on the SCOPE 13-class grid; its moments bf and
! G(mu); and the single-scatter pair of a leaf or wood surface (Goel & Strebel 1984;
! Verhoef 1984; van der Tol et al. 2009). The equations are
! docs/science/canopy_radiation_transfer.md §1-§2, eqs. (1)-(4); ed_pft_derive builds
! the per-PFT optics from these kernels and ed_canopy_radiation the per-cohort.
!---~---
module ed_optics_lib
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : pi, pio2, tiny_num
   implicit none
   private

   public :: N_LEAF_CLASS
   public :: beta_lidf, leaf_bf, gfun_direct
   public :: scatter_pair

   !--- SCOPE's 13 inclination classes: 10 degrees wide to 80, then 2 degrees to 90.
   integer(ik), parameter :: N_LEAF_CLASS = 13_ik
   !--- Class midpoints [deg].
   real(wp), parameter :: LEAF_CLASS_DEG(N_LEAF_CLASS) =                                      &
      [  5.0_wp, 15.0_wp, 25.0_wp, 35.0_wp, 45.0_wp, 55.0_wp, 65.0_wp, 75.0_wp,               &
        81.0_wp, 83.0_wp, 85.0_wp, 87.0_wp, 89.0_wp ]
   !--- Class upper bounds [deg]; class 1 starts at 0.
   real(wp), parameter :: LEAF_BND_DEG(N_LEAF_CLASS) =                                        &
      [ 10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp, 50.0_wp, 60.0_wp, 70.0_wp, 80.0_wp,               &
        82.0_wp, 84.0_wp, 86.0_wp, 88.0_wp, 90.0_wp ]

contains

   !---~---
   !   The Beta LIDF: the leaf-area share of each inclination class, for shape
   ! parameters (p, q) (docs/science/canopy_radiation_transfer.md eq. (1)), by a
   ! midpoint rule: no special function, and the singular endpoints are never hit.
   ! The density is summed relative to its largest value on the grid, so a narrow
   ! distribution, whose unnormalized density underflows everywhere, still sums to
   ! one; the Beta normalizing constant cancels.
   !---~---
   pure function beta_lidf(p, q) result(lidf)
      real(wp), intent(in) :: p, q
      real(wp)             :: lidf(N_LEAF_CLASS)
      integer(ik), parameter :: NSUB = 64_ik           ! midpoint panels per class
      real(wp)    :: tlo, dt, acc, log_max
      integer(ik) :: k, j
      log_max = -huge(1.0_wp)
      tlo     = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         dt = (LEAF_BND_DEG(k) / 90.0_wp - tlo) / real(NSUB, wp)
         do j = 1_ik, NSUB
            log_max = max(log_max, beta_log_kernel(tlo + (real(j, wp) - 0.5_wp) * dt, p, q))
         end do
         tlo = LEAF_BND_DEG(k) / 90.0_wp
      end do
      tlo = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         dt  = (LEAF_BND_DEG(k) / 90.0_wp - tlo) / real(NSUB, wp)
         acc = 0.0_wp
         do j = 1_ik, NSUB
            acc = acc + exp(beta_log_kernel(tlo + (real(j, wp) - 0.5_wp) * dt, p, q) - log_max) * dt
         end do
         lidf(k) = acc
         tlo     = LEAF_BND_DEG(k) / 90.0_wp
      end do
      lidf = lidf / sum(lidf)                           ! the largest panel adds its whole dt
   end function beta_lidf

   !---~---
   !   The log of the unnormalized Beta density, (p-1) ln t + (q-1) ln(1-t), with t
   ! kept off 0 and 1 so an integrable endpoint singularity stays finite.
   !---~---
   pure function beta_log_kernel(t, p, q) result(f)
      real(wp), intent(in) :: t, p, q
      real(wp)             :: f, tt
      tt = min(max(t, 1.0e-9_wp), 1.0_wp - 1.0e-9_wp)   ! clamp-ok: the beta density's domain (0, 1)
      f  = (p - 1.0_wp) * log(tt) + (q - 1.0_wp) * log(1.0_wp - tt)
   end function beta_log_kernel

   !--- The midpoint inclination of class k [rad].
   pure function leaf_class_angle(k) result(theta)
      integer(ik), intent(in) :: k
      real(wp)                :: theta
      theta = LEAF_CLASS_DEG(k) * pi / 180.0_wp
   end function leaf_class_angle

   !---~---
   !   bf = <cos^2 theta_leaf> over the LIDF (docs/science/canopy_radiation_transfer.md
   ! eq. (2)).
   !---~---
   pure function leaf_bf(lidf) result(bf)
      real(wp), intent(in) :: lidf(:)
      real(wp)             :: bf
      integer(ik)          :: k
      bf = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         bf = bf + lidf(k) * cos(leaf_class_angle(k)) ** 2
      end do
   end function leaf_bf

   !---~---
   !   The Ross G-function G(mu), the mean projection of leaf area toward the sun, by
   ! SCOPE's volscatt projection (docs/science/canopy_radiation_transfer.md eq. (3)).
   ! cosz is mu, the cosine of the solar zenith angle.
   !---~---
   pure function gfun_direct(lidf, cosz) result(gee)
      real(wp), intent(in) :: lidf(:)
      real(wp), intent(in) :: cosz
      real(wp)             :: gee
      real(wp)             :: sinz, cs, ss, as, bts, chi_s
      integer(ik)          :: k
      sinz = sqrt(max(0.0_wp, 1.0_wp - cosz * cosz))
      gee  = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         cs = cos(leaf_class_angle(k)) * cosz
         ss = sin(leaf_class_angle(k)) * sinz
         as = max(ss, cs)
         if (as < tiny_num) then
            chi_s = 0.0_wp
         else
            bts   = acos(min(max(-cs / as, -1.0_wp), 1.0_wp))
            chi_s = (2.0_wp / pi) * ((bts - pio2) * cs + sin(bts) * ss)
         end if
         gee = gee + lidf(k) * chi_s
      end do
   end function gfun_direct

   !---~---
   !   Single-scatter albedo omega = rho + tau and asymmetry g = bf (rho - tau)/omega,
   ! with g = 0 for a perfect absorber (docs/science/canopy_radiation_transfer.md eq.
   ! (4)).
   !---~---
   pure subroutine scatter_pair(rho, tau, bf, omega, g)
      real(wp), intent(in)  :: rho, tau, bf
      real(wp), intent(out) :: omega, g
      omega = rho + tau
      if (omega > tiny_num) then
         g = bf * (rho - tau) / omega
      else
         g = 0.0_wp
      end if
   end subroutine scatter_pair

end module ed_optics_lib
