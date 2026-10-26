! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_linear_algebra -- the linear solvers and the matrix exponential every layer shares   !
! (docs/science/numerical_scheme.md §9):                                                    !
!   thomas_solve          a tridiagonal system (the soil heat column)                       !
!   lu_factor, lu_solve   a small dense system (the soil-carbon steady state)               !
!   arrowhead_*           a block arrowhead system (the soil-plant water network)           !
!   matrix_exp            exp(A) (the soil-carbon step)                                     !
!                                                                                          !
! Sizes are passed, not taken from the arrays, so a caller can hand over fixed-size storage  !
! longer than the system it solves.                                                         !
!==========================================================================================!
module meds_linear_algebra
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: thomas_solve
   public :: lu_factor, lu_solve
   public :: arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   public :: matrix_exp

   !----- A block arrowhead system, filled by the caller, and its factorisation. --------------!
   type :: arrowhead_system_t
      integer(ik) :: nc = 0_ik                !< core unknowns
      integer(ik) :: nb = 0_ik                !< blocks
      integer(ik) :: mb = 0_ik                !< unknowns per block
      real(wp),    allocatable :: a(:,:)      !< (nc, nc)     core block
      real(wp),    allocatable :: b(:,:)      !< (nc, nb)     core rows, each block's first column
      real(wp),    allocatable :: c(:,:)      !< (nc, nb)     each block's first row, over the core columns
      real(wp),    allocatable :: d(:,:,:)    !< (mb, mb, nb) the diagonal blocks
      real(wp),    allocatable :: lu(:,:)     !< (nc, nc)     the LU of the core after the blocks are folded in
      real(wp),    allocatable :: dinv(:,:,:) !< (mb, mb, nb) the inverses of the diagonal blocks
      integer(ik), allocatable :: piv(:)      !< (nc)         the core LU's row interchanges
   end type arrowhead_system_t

contains

   !----- Solve a(k) x(k-1) + b(k) x(k) + c(k) x(k+1) = d(k), k = 1..n, by the Thomas sweep   !
   !      (a(1) and c(n) are not used). No pivoting: the systems it solves are diagonally      !
   !      dominant. A pivot is kept away from zero, so a degenerate input gives a finite x. ---!
   pure subroutine thomas_solve(a, b, c, d, x, n)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: a(n), b(n), c(n), d(n)
      real(wp),    intent(out) :: x(n)
      real(wp)    :: cp(n), dp(n), denom
      integer(ik) :: k
      x = 0.0_wp
      denom = sign(max(abs(b(1)), tiny_num), b(1))
      cp(1) = c(1) / denom
      dp(1) = d(1) / denom
      do k = 2_ik, n
         denom = b(k) - a(k) * cp(k-1)
         denom = sign(max(abs(denom), tiny_num), denom)
         cp(k) = c(k) / denom
         dp(k) = (d(k) - a(k) * dp(k-1)) / denom
      end do
      x(n) = dp(n)
      do k = n - 1_ik, 1_ik, -1_ik
         x(k) = dp(k) - cp(k) * x(k+1)
      end do
   end subroutine thomas_solve

   !----- LU factorisation with partial pivoting, in place: L unit lower triangular below the   !
   !      diagonal, U upper triangular on and above it. ok is false on a zero or NaN pivot. ---!
   pure subroutine lu_factor(a, piv, n, ok)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: a(n, n)
      integer(ik), intent(out)   :: piv(n)    !< the row interchanged with row k at step k
      logical,     intent(out)   :: ok
      integer(ik) :: i, j, k, p
      real(wp)    :: row(n), f
      ok = .true.
      do k = 1_ik, n
         p = k
         do i = k + 1_ik, n
            if (abs(a(i, k)) > abs(a(p, k))) p = i
         end do
         piv(k) = p
         if (a(p, k) == 0.0_wp .or. a(p, k) /= a(p, k)) then
            ok = .false.
            return
         end if
         if (p /= k) then
            row = a(k, :)
            a(k, :) = a(p, :)
            a(p, :) = row
         end if
         do i = k + 1_ik, n
            f = a(i, k) / a(k, k)
            a(i, k) = f
            if (f == 0.0_wp) cycle
            do j = k + 1_ik, n
               a(i, j) = a(i, j) - f * a(k, j)
            end do
         end do
      end do
   end subroutine lu_factor

   !----- Solve with lu_factor's result: x is the right-hand side on entry, the solution on   !
   !      exit. ----------------------------------------------------------------------------!
   pure subroutine lu_solve(a, piv, n, x)
      integer(ik), intent(in)    :: n
      real(wp),    intent(in)    :: a(n, n)
      integer(ik), intent(in)    :: piv(n)
      real(wp),    intent(inout) :: x(n)
      integer(ik) :: i, j
      real(wp)    :: t
      do i = 1_ik, n
         if (piv(i) /= i) then
            t = x(i)
            x(i) = x(piv(i))
            x(piv(i)) = t
         end if
      end do
      do i = 2_ik, n
         do j = 1_ik, i - 1_ik
            x(i) = x(i) - a(i, j) * x(j)
         end do
      end do
      do i = n, 1_ik, -1_ik
         do j = i + 1_ik, n
            x(i) = x(i) - a(i, j) * x(j)
         end do
         x(i) = x(i) / a(i, i)
      end do
   end subroutine lu_solve

   !----- Size a system for nc core unknowns and nb blocks of mb, keeping the storage when the  !
   !      shape is unchanged, and zero its blocks. -----------------------------------------!
   pure subroutine arrowhead_allocate(s, nc, nb, mb)
      type(arrowhead_system_t), intent(inout) :: s
      integer(ik),              intent(in)    :: nc, nb, mb
      if (s%nc /= nc .or. s%nb /= nb .or. s%mb /= mb .or. .not. allocated(s%a)) then
         if (allocated(s%a))   deallocate(s%a, s%b, s%c, s%d, s%lu, s%dinv, s%piv)
         allocate(s%a(nc, nc), s%b(nc, nb), s%c(nc, nb), s%d(mb, mb, nb), s%lu(nc, nc), s%dinv(mb, mb, nb),  &
                  s%piv(nc))
         s%nc = nc
         s%nb = nb
         s%mb = mb
      end if
      s%a = 0.0_wp
      s%b = 0.0_wp
      s%c = 0.0_wp
      s%d = 0.0_wp
   end subroutine arrowhead_allocate

   !----- Factor the system with add_core and add_block added to the diagonals: invert each     !
   !      block, fold it into the core, and LU-factor the core. The blocks are left as they     !
   !      were, so they can be factored again with other additions. ok is false if a block or  !
   !      the core is singular. -----------------------------------------------------------!
   pure subroutine arrowhead_factor(s, add_core, add_block, ok)
      type(arrowhead_system_t), intent(inout) :: s
      real(wp),                 intent(in)    :: add_core(:), add_block(:, :)   !< (nc), (mb, nb)
      logical,                  intent(out)   :: ok
      integer(ik) :: ib, i, j
      real(wp)    :: e11, cj
      ok = .true.
      s%lu = s%a
      do i = 1_ik, s%nc
         s%lu(i, i) = s%lu(i, i) + add_core(i)
      end do
      do ib = 1_ik, s%nb
         s%dinv(:, :, ib) = s%d(:, :, ib)
         do i = 1_ik, s%mb
            s%dinv(i, i, ib) = s%dinv(i, i, ib) + add_block(i, ib)
         end do
         call invert_small(s%dinv(:, :, ib), s%mb, ok)
         if (.not. ok) return
         e11 = s%dinv(1, 1, ib)
         do j = 1_ik, s%nc
            cj = e11 * s%c(j, ib)
            if (cj == 0.0_wp) cycle
            do i = 1_ik, s%nc
               s%lu(i, j) = s%lu(i, j) - s%b(i, ib) * cj
            end do
         end do
      end do
      call lu_factor(s%lu, s%piv, s%nc, ok)
   end subroutine arrowhead_factor

   !----- Solve the factored system for one right-hand side (r_core, r_block). ---------------!
   pure subroutine arrowhead_solve(s, r_core, r_block, x_core, x_block)
      type(arrowhead_system_t), intent(in)  :: s
      real(wp),                 intent(in)  :: r_core(:), r_block(:, :)   !< (nc), (mb, nb)
      real(wp),                 intent(out) :: x_core(:), x_block(:, :)
      real(wp)    :: rc(s%nc), cx
      integer(ik) :: ib, i, l
      rc = r_core(1:s%nc)
      !----- Each block's own solution, folded onto the core through its first unknown. -------!
      do ib = 1_ik, s%nb
         do i = 1_ik, s%mb
            x_block(i, ib) = 0.0_wp
            do l = 1_ik, s%mb
               x_block(i, ib) = x_block(i, ib) + s%dinv(i, l, ib) * r_block(l, ib)
            end do
         end do
         if (x_block(1, ib) /= 0.0_wp) rc = rc - s%b(:, ib) * x_block(1, ib)
      end do
      call lu_solve(s%lu, s%piv, s%nc, rc)
      x_core(1:s%nc) = rc
      !----- Back to the blocks, through the core's solution. ---------------------------------!
      do ib = 1_ik, s%nb
         cx = sum(s%c(:, ib) * rc)
         x_block(1:s%mb, ib) = x_block(1:s%mb, ib) - s%dinv(1:s%mb, 1, ib) * cx
      end do
   end subroutine arrowhead_solve

   !----- Invert a small block in place: closed form for 1x1 and 2x2, Gauss-Jordan otherwise.  !
   !      ok is false if it is singular. ---------------------------------------------------!
   pure subroutine invert_small(m, n, ok)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: m(n, n)
      logical,     intent(out)   :: ok
      real(wp)    :: det, m11
      ok = .true.
      if (n == 1_ik) then
         if (m(1, 1) == 0.0_wp .or. m(1, 1) /= m(1, 1)) then
            ok = .false.
            return
         end if
         m(1, 1) = 1.0_wp / m(1, 1)
      else if (n == 2_ik) then
         det = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
         if (det == 0.0_wp .or. det /= det) then
            ok = .false.
            return
         end if
         m11 = m(1, 1)
         m(1, 1) =  m(2, 2) / det
         m(2, 2) = m11 / det
         m(1, 2) = -m(1, 2) / det
         m(2, 1) = -m(2, 1) / det
      else
         call invert_gauss_jordan(m, n, ok)
      end if
   end subroutine invert_small

   pure subroutine invert_gauss_jordan(m, n, ok)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: m(n, n)
      logical,     intent(out)   :: ok
      real(wp)    :: aug(n, 2*n), piv_val, f, row(2*n)
      integer(ik) :: i, j, p
      aug = 0.0_wp
      aug(:, 1:n) = m
      do i = 1_ik, n
         aug(i, n + i) = 1.0_wp
      end do
      ok = .true.
      do j = 1_ik, n
         p = j
         do i = j + 1_ik, n
            if (abs(aug(i, j)) > abs(aug(p, j))) p = i
         end do
         piv_val = aug(p, j)
         if (piv_val == 0.0_wp .or. piv_val /= piv_val) then
            ok = .false.
            return
         end if
         if (p /= j) then
            row = aug(j, :)
            aug(j, :) = aug(p, :)
            aug(p, :) = row
         end if
         aug(j, :) = aug(j, :) / piv_val
         do i = 1_ik, n
            if (i == j) cycle
            f = aug(i, j)
            if (f /= 0.0_wp) aug(i, :) = aug(i, :) - f * aug(j, :)
         end do
      end do
      m = aug(:, n+1:2*n)
   end subroutine invert_gauss_jordan

   !----- exp(A) by scaling and squaring: halve A s times until its largest row sum is at most  !
   !      1/2, then exp_scaled_squared. --------------------------------------------------!
   pure subroutine matrix_exp(a_in, e, m)
      integer(ik), intent(in)  :: m
      real(wp),    intent(in)  :: a_in(m, m)
      real(wp),    intent(out) :: e(m, m)
      real(wp)    :: nrm, scale
      integer(ik) :: i, j, s
      nrm = 0.0_wp
      do i = 1_ik, m
         scale = 0.0_wp
         do j = 1_ik, m
            scale = scale + abs(a_in(i, j))
         end do
         nrm = max(nrm, scale)
      end do
      s = 0_ik
      scale = 1.0_wp
      do while (nrm * scale > 0.5_wp)
         s = s + 1_ik
         scale = scale * 0.5_wp
      end do
      call exp_scaled_squared(a_in, e, m, s)
   end subroutine matrix_exp

   !----- exp(A) as the 12-term Taylor series of A/2^s, squared s times. -------------------!
   pure subroutine exp_scaled_squared(a_in, e, m, s)
      integer(ik), intent(in)  :: m, s
      real(wp),    intent(in)  :: a_in(m, m)
      real(wp),    intent(out) :: e(m, m)
      integer(ik), parameter :: N_TAYLOR = 12_ik
      real(wp)    :: a(m, m), term(m, m), tmp(m, m)
      integer(ik) :: i, k, sq
      a = a_in * (2.0_wp ** (-s))
      e = 0.0_wp
      term = 0.0_wp
      do i = 1_ik, m
         e(i, i)    = 1.0_wp
         term(i, i) = 1.0_wp
      end do
      do k = 1_ik, N_TAYLOR
         call matmul_square(term, a, tmp, m)
         term = tmp / real(k, wp)
         e = e + term
      end do
      do sq = 1_ik, s
         call matmul_square(e, e, tmp, m)
         e = tmp
      end do
   end subroutine exp_scaled_squared

   !----- c = a b for square matrices, written out so no array temporary is made              !
   !      (docs/building.md). ------------------------------------------------------------!
   pure subroutine matmul_square(a, b, c, m)
      integer(ik), intent(in)  :: m
      real(wp),    intent(in)  :: a(m, m), b(m, m)
      real(wp),    intent(out) :: c(m, m)
      integer(ik) :: i, j, k
      real(wp)    :: acc
      do j = 1_ik, m
         do i = 1_ik, m
            acc = 0.0_wp
            do k = 1_ik, m
               acc = acc + a(i, k) * b(k, j)
            end do
            c(i, j) = acc
         end do
      end do
   end subroutine matmul_square

end module meds_linear_algebra
