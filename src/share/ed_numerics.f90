! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_numerics -- small numerical helpers every layer shares: clamps, a smooth step, the
! smooth minimum of two limits, a rounded max(x, 0), the step-size factor, the Newton
! solvers' step safeguards, the Illinois steps of a bracketed root, linear interpolation,
! weighted means, an ordering and the canopy layers it gives, and the test that an
! operator kept what it only moves. The linear solvers are in ed_linear_algebra.
!---~---
module ed_numerics
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : safe_exp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private

   public :: quadratic_smaller_root, adaptive_step_update, smooth_max0
   public :: sufficient_decrease, fraction_to_boundary
   public :: root_bracket_t, regula_falsi_point, illinois_narrow
   public :: logistic, clamp01, clamp, linear_interpolate, weighted_mean
   public :: ascending_order, equal_height_layers
   public :: conserved, pair_mean

   !---~---
   !   A line-search step is accepted when the merit falls by at least this share of
   ! the step taken (Armijo's sufficient-decrease condition).
   !---~---
   real(wp), parameter :: ARMIJO_C = 1.0e-4_wp

   !---~---
   !   A bracket [a, b] on a root of a continuous function, with residuals fa and fb of
   ! opposite signs, for the Illinois method: regula falsi that halves the kept end's
   ! residual when the same end is kept twice, so both ends close in. The caller
   ! evaluates the function: regula_falsi_point gives the next point to try, and
   ! illinois_narrow narrows the bracket with that point's residual.
   !---~---
   type :: root_bracket_t
      real(wp)    :: a  = 0.0_wp      ! the left end
      real(wp)    :: fa = 0.0_wp      ! the residual at a (the Illinois step may have halved it)
      real(wp)    :: b  = 0.0_wp      ! the right end
      real(wp)    :: fb = 0.0_wp      ! the residual at b (likewise)
      integer(ik) :: kept = 0_ik      ! the end kept last time: -1 left, +1 right, 0 none
   end type root_bracket_t

contains

   !---~---
   !   ord(k) is the index of the k-th smallest of x(1:n), the larger index first among
   ! equals (the cohorts from the canopy bottom up). The cohorts are kept tallest
   ! first, so the model's input is reversed in O(n); any other input takes a
   ! selection sort, which gives the same order on a non-increasing x.
   !---~---
   pure subroutine ascending_order(x, n, ord)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: x(n)
      integer(ik), intent(out) :: ord(n)
      logical     :: used(n)
      real(wp)    :: xmin
      integer(ik) :: j, k, imin
      if (.not. any(x(1:n-1) < x(2:n))) then
         ord = [(n - k + 1_ik, k = 1_ik, n)]
         return
      end if
      used = .false.
      do k = 1_ik, n
         imin = 0_ik
         xmin = huge(1.0_wp)
         do j = 1_ik, n
            if (.not. used(j) .and. x(j) <= xmin) then
               xmin = x(j)
               imin = j
            end if
         end do
         ord(k) = imin
         used(imin) = .true.
      end do
   end subroutine ascending_order

   !---~---
   !   The canopy layers of cohorts in height order, either way up. Cohorts of exactly
   ! equal height are co-dominant and share one layer: layer(k) numbers the run of
   ! equal heights that height(k) is in, from 1 at height(1), and n_layer counts the
   ! runs. The radiation, the wind and the overtopping LAI all group cohorts by it.
   ! Exact equality is the right test: every PFT recruits at one shared size, so
   ! recruits tie bit for bit.
   !---~---
   pure subroutine equal_height_layers(height, n, layer, n_layer)
      integer(ik), intent(in)  :: n            ! number of cohorts
      real(wp),    intent(in)  :: height(n)    ! [m] ascending or descending
      integer(ik), intent(out) :: layer(n)     ! each cohort's layer, 1 at height(1)
      integer(ik), intent(out) :: n_layer      ! number of layers, 0 when n = 0
      integer(ik) :: k
      n_layer = 0_ik
      do k = 1_ik, n
         if (k == 1_ik) then
            n_layer = 1_ik
         else if (height(k) /= height(k-1_ik)) then
            n_layer = n_layer + 1_ik
         end if
         layer(k) = n_layer
      end do
   end subroutine equal_height_layers

   !---~---
   !   The smaller root of theta x^2 - (a + b) x + a b = 0: a smooth min(a, b), equal
   ! to it at theta = 1 (docs/science/leaf_gas_exchange.md §2).
   !---~---
   elemental pure function quadratic_smaller_root(theta, a, b) result(x)
      real(wp), intent(in) :: theta, a, b
      real(wp)             :: x, s, disc
      s    = a + b
      disc = max(s * s - 4.0_wp * theta * a * b, 0.0_wp)
      x    = (s - sqrt(disc)) / (2.0_wp * max(theta, tiny(theta)))    ! theta in (0,1]; floor guards /0
   end function quadratic_smaller_root

   !---~---
   !   The factor a step size is multiplied by: safety * err^(-1/(p_order + 1)),
   ! clamped to [fmin, fmax] (docs/science/numerical_scheme.md eq. (4)).
   !---~---
   elemental pure function adaptive_step_update(err, safety, fmin, fmax, p_order) result(factor)
      real(wp),    intent(in)           :: err, safety, fmin, fmax
      integer(ik), intent(in), optional :: p_order   ! order of the embedded error estimate (default 1)
      real(wp)             :: factor, expo
      integer(ik) :: p
      p = 1_ik
      if (present(p_order)) p = p_order
      expo = -1.0_wp / real(p + 1_ik, wp)
      factor = min(fmax, max(fmin, safety * max(err, tiny(err)) ** expo))   ! err = 0 gives fmax
   end function adaptive_step_update

   !---~---
   !   max(x, 0) with its corner rounded over [0, width] by the cubic 2x^2/w - x^3/w^2,
   ! which matches both sides in value and slope, so Newton never meets a corner.
   ! Every one-way water flow uses it: root uptake, the leaves' supply limit and
   ! ground evaporation (docs/science/soil_plant_water.md eq. (1)).
   !---~---
   elemental subroutine smooth_max0(x, width, f, df)
      real(wp), intent(in)  :: x, width
      real(wp), intent(out) :: f, df    ! the rounded max and its slope
      real(wp) :: t
      if (x <= 0.0_wp) then
         f = 0.0_wp
         df = 0.0_wp
      else if (x < width) then
         t  = x / width
         f  = x * t * (2.0_wp - t)
         df = t * (4.0_wp - 3.0_wp * t)
      else
         f = x
         df = 1.0_wp
      end if
   end subroutine smooth_max0

   !---~---
   !   Regula falsi's point on a bracket: where the line through its two ends crosses
   ! zero, or the midpoint when rounding puts that outside (a, b).
   !---~---
   elemental pure real(wp) function regula_falsi_point(br) result(c)
      type(root_bracket_t), intent(in) :: br
      c = (br%a * br%fb - br%b * br%fa) / (br%fb - br%fa)
      if (.not. (c > br%a .and. c < br%b)) c = 0.5_wp * (br%a + br%b)
   end function regula_falsi_point

   !---~---
   !   The bracket narrowed to the side of c, with fc the residual at c, by the Illinois
   ! rule: the end kept a second time in a row has its residual halved, which pulls the
   ! next point towards it so that end moves too.
   !---~---
   elemental pure subroutine illinois_narrow(br, c, fc)
      type(root_bracket_t), intent(inout) :: br
      real(wp),             intent(in)    :: c, fc
      if (br%fa * fc < 0.0_wp) then           ! the root is in (a, c): the left end is kept
         br%b  = c
         br%fb = fc
         if (br%kept == -1_ik) br%fa = 0.5_wp * br%fa
         br%kept = -1_ik
      else                                    ! the root is in (c, b): the right end is kept
         br%a  = c
         br%fa = fc
         if (br%kept == 1_ik) br%fb = 0.5_wp * br%fb
         br%kept = 1_ik
      end if
   end subroutine illinois_narrow

   !---~---
   !   Armijo's test for a backtracking line search: true when the merit at a trial
   ! step of `step_share` of the Newton step has fallen to (1 - ARMIJO_C step_share)
   ! of the merit at the start.
   !---~---
   elemental pure logical function sufficient_decrease(merit_trial, merit, step_share) result(ok)
      real(wp), intent(in) :: merit_trial, merit, step_share
      ok = merit_trial <= (1.0_wp - ARMIJO_C * step_share) * merit
   end function sufficient_decrease

   !---~---
   !   The largest share of a Newton step dx that keeps x, a quantity that is zero at
   ! its boundary, from covering more than `share` of its way there (the
   ! fraction-to-the-boundary rule); a step away from the boundary is not limited.
   !---~---
   elemental pure real(wp) function fraction_to_boundary(x, dx, share) result(lam)
      real(wp), intent(in) :: x, dx, share
      if (dx < 0.0_wp) then
         lam = share * x / (-dx)
      else
         lam = huge(x)
      end if
   end function fraction_to_boundary

   !---~---
   !   The logistic 1 / (1 + exp(-z)), a smooth step from 0 to 1
   ! (docs/science/plant_phenology.md eq. (2)).
   !---~---
   elemental pure function logistic(z) result(f)
      real(wp), intent(in) :: z
      real(wp)             :: f
      f = 1.0_wp / (1.0_wp + safe_exp(-z))
   end function logistic

   !--- x clamped to [0, 1].
   elemental pure function clamp01(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = min(1.0_wp, max(0.0_wp, x))
   end function clamp01

   !--- x clamped to [lo, hi], with lo <= hi.
   elemental pure function clamp(x, lo, hi) result(y)
      real(wp), intent(in) :: x, lo, hi
      real(wp)             :: y
      y = min(hi, max(lo, x))
   end function clamp

   !--- The value a fraction w of the way from a to b: (1 - w) a + w b.
   elemental pure real(wp) function linear_interpolate(a, b, w) result(x)
      real(wp), intent(in) :: a, b, w
      x = (1.0_wp - w) * a + w * b
   end function linear_interpolate

   !---~---
   !   The sum of x(k) w(k) over k = 1..n, in that order: a weighted mean when the
   ! weights sum to one.
   !---~---
   pure function weighted_mean(x, w, n) result(xbar)
      real(wp),    intent(in) :: x(:), w(:)
      integer(ik), intent(in) :: n
      real(wp) :: xbar
      integer(ik) :: k
      xbar = 0.0_wp
      do k = 1_ik, n
         xbar = xbar + x(k) * w(k)
      end do
   end function weighted_mean

   !---~---
   !   True when an amount an operator only moves is unchanged to round-off: |after
   ! before| <= rtol * magnitude, plus the smallest difference the number format
   ! resolves. A NaN fails (docs/science/numerical_scheme.md §9).
   !---~---
   elemental pure logical function conserved(before, after, magnitude, rtol) result(ok)
      real(wp), intent(in) :: before, after   ! the amount before and after the operator
      real(wp), intent(in) :: magnitude       ! sum of the absolute values of its terms
      real(wp), intent(in) :: rtol            ! relative tolerance
      ok = .not. (ieee_is_nan(before) .or. ieee_is_nan(after) .or. ieee_is_nan(magnitude))
      if (ok) ok = abs(after - before) <= rtol * magnitude + tiny(magnitude) / epsilon(magnitude)
   end function conserved

   !---~---
   !   The weighted mean of two values, written so it stays between them for any
   ! weights however small; with no weight at all it keeps a. There is no size
   ! threshold: a leafless cohort's tiny leaf carbon still weighs its traits in a
   ! fusion.
   !---~---
   elemental pure real(wp) function pair_mean(a, wa, b, wb) result(m)
      real(wp), intent(in) :: a, wa   ! the first value and its weight (>= 0)
      real(wp), intent(in) :: b, wb   ! the second value and its weight (>= 0)
      if (wa + wb > 0.0_wp) then
         m = a + (b - a) * (wb / (wa + wb))
      else
         m = a
      end if
   end function pair_mean

end module ed_numerics
