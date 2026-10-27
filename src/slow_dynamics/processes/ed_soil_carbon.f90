! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_soil_carbon -- the daily step of a patch's soil carbon pools, which obey the carbon
! matrix ODE dX/dt = B*I(t) + A*xi(t)*K*X(t) [kgC/m2/day] (docs/science/soil_carbon.md
! eq. (1)); the chemistry (xi, A, K, Rh) is ed_respiration's.
!   * necromass_to_litter, add_necromass -- dead plant carbon to litter (§5).
!   * build_litter_input             -- the patch's litter onto the source term u (§5).
!   * soil_carbon_step               -- one day, by forward Euler or the exact exponential,
!                                       on the fast loop's accumulated decay (§6); the new
!                                       pools and the checks come back in
!                                       soil_carbon_out_t (§7, §9).
!   * solve_soil_carbon_steady_state -- the semi-analytic spin-up (SASU, §8).
!   * soil_carbon_bad_pool           -- plausibility, as distinct from conservation (§9).
! The pools (soil_carbon_state_t) and the litter (litter_input_t) are patch state, in
! ed_patch_state_vars.
!---~---
module ed_soil_carbon
   use ed_kinds,               only : wp, ik
   use ed_patch_state_vars,    only : litter_input_t, soil_carbon_state_t
   use ed_biogeochem_params,   only : N_SOIL_POOL, soil_carbon_params_t, DECOMP_STEP_EXPM,               &
                                      IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL,        &
                                      IP_MICR, IP_SLOW, IP_PASSIVE
   use ed_constants,           only : tiny_num
   use ed_linear_algebra,      only : matrix_exp, lu_factor, lu_solve
   use ed_respiration,         only : decomposition_transfer_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private

   public :: soil_carbon_out_t
   public :: necromass_to_litter, add_necromass, build_litter_input, soil_carbon_step
   public :: solve_soil_carbon_steady_state, soil_carbon_steady_pools
   public :: soil_carbon_bad_pool, SOILC_MAX_PLAUSIBLE, SOILC_NEG_TOL

   !---~---
   !   One daily step's result: the new pools, which the caller commits, and its
   ! carbon-mass checks (soil_carbon.md §7, §9).
   !---~---
   type :: soil_carbon_out_t
      type(soil_carbon_state_t) :: pools     !--- the pools and lignin after the step
      !---~---
      !   rh_seam_gap: zero by construction on every step, since the patch structure
      ! changes before the fast window opens; do not widen a tolerance for it (§7).
      !---~---
      real(wp) :: rh_seam_gap  = 0.0_wp   !--- [kgC/m2/day] the step's Rh less the fast loop's (set by the caller)
      real(wp) :: lignin_resid = 0.0_wp   !--- [kgC/m2/day] max_s |dL_s - (lignin_in_s - (1-surv_s)*L_s)|; tracer check
      !---~---
      !   lambda, the largest fraction of any pool the step withdraws, says whether
      ! freezing the pool is sound; past 1 forward Euler goes negative (§6).
      !---~---
      real(wp)    :: lambda_max  = 0.0_wp   !--- [-] max_j dvec_j = max_j (xi_int_j * K_j)
      integer(ik) :: lambda_pool = 0_ik     !--- which pool attained it (IP_* index)
   end type soil_carbon_out_t

   !---~---
   !   What "physically impossible" means for a pool (soil_carbon.md §9): a carbon pool
   ! is never negative beyond round-off, and 1e4 kgC/m2 (~500x the richest real soil)
   ! catches divergence, not implausible science.
   !---~---
   real(wp), parameter :: SOILC_NEG_TOL       = 1.0e-9_wp   !--- [kgC/m2] round-off allowance below zero
   real(wp), parameter :: SOILC_MAX_PLAUSIBLE = 1.0e4_wp    !--- [kgC/m2] per-pool divergence ceiling

contains

   !---~---
   !   Split per-area dead carbon [kgC/m2] by tissue: canopy (leaves and reproduction) all
   ! above ground and fine roots all below, both labile by f_labile_leaf; wood and storage
   ! above/below by aboveground_frac, wood labile by f_labile_stem and storage by
   ! f_labile_leaf. Only the structural streams carry lignin. Pass 0 for a tissue that
   ! does not apply (the turnover call has no wood or storage).
   !---~---
   pure subroutine necromass_to_litter(canopy_c, root_c, wood_c, storage_c,                          &
                                       f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac,   &
                                       labile_grnd, labile_soil, struct_grnd, struct_soil,             &
                                       lignin_grnd, lignin_soil)
      real(wp), intent(in)  :: canopy_c, root_c, wood_c, storage_c
      real(wp), intent(in)  :: f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac
      real(wp), intent(out) :: labile_grnd, labile_soil, struct_grnd, struct_soil
      real(wp), intent(out) :: lignin_grnd, lignin_soil
      real(wp) :: canopy_lab, canopy_str, root_lab, root_str, woody_lab, woody_str

      canopy_lab = canopy_c * f_labile_leaf
      canopy_str = canopy_c * (1.0_wp - f_labile_leaf)
      root_lab   = root_c * f_labile_leaf
      root_str   = root_c * (1.0_wp - f_labile_leaf)
      woody_lab  = wood_c * f_labile_stem + storage_c * f_labile_leaf
      woody_str  = wood_c * (1.0_wp - f_labile_stem) + storage_c * (1.0_wp - f_labile_leaf)

      labile_grnd = canopy_lab + woody_lab * aboveground_frac
      labile_soil = root_lab   + woody_lab * (1.0_wp - aboveground_frac)
      struct_grnd = canopy_str + woody_str * aboveground_frac
      struct_soil = root_str   + woody_str * (1.0_wp - aboveground_frac)
      lignin_grnd = lignin_frac * struct_grnd
      lignin_soil = lignin_frac * struct_soil
   end subroutine necromass_to_litter

   !---~---
   !   Add one group of dead tissue [kgC/m2] to a patch's litter accumulator: the pool
   ! streams (reproduction going with the canopy) and the same carbon by tissue.
   !---~---
   pure subroutine add_necromass(lit, leaf_c, root_c, wood_c, storage_c, repro_c,                   &
                                 f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac)
      type(litter_input_t), intent(inout) :: lit
      real(wp),             intent(in)    :: leaf_c, root_c, wood_c, storage_c, repro_c
      real(wp),             intent(in)    :: f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac
      real(wp) :: lab_g, lab_s, str_g, str_s, lig_g, lig_s

      call necromass_to_litter(leaf_c + repro_c, root_c, wood_c, storage_c,                         &
                               f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac,         &
                               lab_g, lab_s, str_g, str_s, lig_g, lig_s)
      lit%labile_grnd = lit%labile_grnd + lab_g
      lit%labile_soil = lit%labile_soil + lab_s
      lit%struct_grnd = lit%struct_grnd + str_g
      lit%struct_soil = lit%struct_soil + str_s
      lit%lignin_grnd = lit%lignin_grnd + lig_g
      lit%lignin_soil = lit%lignin_soil + lig_s
      lit%leaf_c      = lit%leaf_c      + leaf_c
      lit%fineroot_c  = lit%fineroot_c  + root_c
      lit%wood_c      = lit%wood_c      + wood_c
      lit%storage_c   = lit%storage_c   + storage_c
      lit%repro_c     = lit%repro_c     + repro_c
   end subroutine add_necromass

   !---~---
   !   Place the already-split litter streams (four carbon streams and two lignin fluxes)
   ! into the source term u. The per-cohort splits are made before summing to the patch,
   ! so this does no chemistry of its own (docs/science/soil_carbon.md §5).
   !---~---
   pure subroutine build_litter_input(litter, u, lignin_in)
      type(litter_input_t), intent(in)  :: litter
      real(wp),             intent(out) :: u(N_SOIL_POOL)   !--- [kgC/m2/day] input to each pool (1-4 nonzero)
      real(wp),             intent(out) :: lignin_in(2)     !--- [kgC/m2/day] lignin flux to struct_grnd/struct_soil
      u = 0.0_wp
      u(IP_FAST_GRND)   = litter%labile_grnd
      u(IP_FAST_SOIL)   = litter%labile_soil
      u(IP_STRUCT_GRND) = litter%struct_grnd
      u(IP_STRUCT_SOIL) = litter%struct_soil
      lignin_in(1) = litter%lignin_grnd
      lignin_in(2) = litter%lignin_soil
   end subroutine build_litter_input

   !---~---
   !   Advance the pools one day, X_new = X + u + A*diag(xi_int*K)*X, with xi_int the fast
   ! loop's accumulated decay integral of the day rather than xi from daily-mean drivers,
   ! so the pool loss equals the day's respired carbon by construction
   ! (docs/science/soil_carbon.md §6, eqs. (5)-(6)). DECOMP_STEP_EULER is forward Euler;
   ! DECOMP_STEP_EXPM the exact exponential of the frozen operator, for large accelerated
   ! steps. rh_today = litter in - pool change for either solver; it is an audit value
   ! and is not fed back to the canopy air. u is one day's input [kgC/m2/day]. The new
   ! pools come back in out%pools; the caller commits them.
   !---~---
   subroutine soil_carbon_step(pools, u, lignin_in, xi_int, params, rh_today, out)
      type(soil_carbon_state_t),  intent(in)    :: pools
      real(wp),                   intent(in)    :: u(N_SOIL_POOL), lignin_in(2)
      real(wp),                   intent(in)    :: xi_int(N_SOIL_POOL)   !--- [day] INT_day xi_j dt (fast accumulator)
      type(soil_carbon_params_t), intent(in)    :: params
      real(wp),                   intent(out)   :: rh_today              !--- [kgC/m2/day] actual pool loss = fast loop's Rh
      type(soil_carbon_out_t),    intent(out)   :: out                   !--- the new pools and the checks
      real(wp) :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: x0(N_SOIL_POOL), x1(N_SOIL_POOL), dvec(N_SOIL_POOL), loss(N_SOIL_POOL)
      real(wp) :: lig0(2), lig1(2), d_s(2), surv(2), litter_in, dc_pool, lig_resid
      integer(ik) :: i, j

      call decomposition_transfer_matrix(pools, params, a_mat, k_diag, er)
      x0 = pools%carbon
      lig0(1) = pools%lignin(1)
      lig0(2) = pools%lignin(2)

      !--- Per-pool day-integrated decay fraction d_j = xi_int_j * K_j [-].
      do j = 1_ik, N_SOIL_POOL
         dvec(j) = xi_int(j) * k_diag(j)
      end do

      !--- Each donor decays at its own rate, and only then is its loss split, A*(xi*K*X) (§4).
      select case (params%step_solver)
      case (DECOMP_STEP_EXPM)
         call step_expm(a_mat, dvec, u, x0, x1)
      case default   ! DECOMP_STEP_EULER
         do j = 1_ik, N_SOIL_POOL
            loss(j) = dvec(j) * x0(j)                 ! [kgC/m2] day-integrated loss of donor j
         end do
         do i = 1_ik, N_SOIL_POOL
            x1(i) = x0(i) + u(i)
            do j = 1_ik, N_SOIL_POOL
               x1(i) = x1(i) + a_mat(i, j) * loss(j)  ! a_ii=-1 removes donor i's own loss; a_ij adds transfers
            end do
         end do
      end select

      !---~---
      !   Lignin decays in lockstep with the structural carbon that holds it, by the
      ! solver's own survival fraction: 1 - d_s under Euler, exp(-d_s) under EXPM;
      ! the linear form under EXPM would drive lignin negative at large d_s (§6).
      !---~---
      d_s(1) = dvec(IP_STRUCT_GRND)
      d_s(2) = dvec(IP_STRUCT_SOIL)
      select case (params%step_solver)
      case (DECOMP_STEP_EXPM)
         surv(1) = exp(-d_s(1))
         surv(2) = exp(-d_s(2))
      case default
         surv(1) = 1.0_wp - d_s(1)
         surv(2) = 1.0_wp - d_s(2)
      end select
      lig1(1) = lig0(1) * surv(1) + lignin_in(1)
      lig1(2) = lig0(2) * surv(2) + lignin_in(2)

      !--- The new state, which the caller commits.
      out%pools = pools
      out%pools%carbon = x1
      out%pools%lignin(1) = lig1(1)
      out%pools%lignin(2) = lig1(2)

      !--- carbon-mass audit (both solvers): rh_today = litter_in - dC_pool (transfers cancel).
      litter_in = sum(u)
      dc_pool   = sum(x1) - sum(x0)
      rh_today  = litter_in - dc_pool
      !---~---
      !   No pool-mass residual: rh_today is defined from the pool change, so it is zero by
      ! definition; the independent check is test_soil_carbon's (soil_carbon.md §9).
      !---~---
      !--- lambda: the largest fraction of any pool this step withdraws (soil_carbon_out_t).
      out%lambda_max  = 0.0_wp
      out%lambda_pool = 0_ik
      do i = 1_ik, N_SOIL_POOL
         if (dvec(i) > out%lambda_max) then
            out%lambda_max  = dvec(i)
            out%lambda_pool = i
         end if
      end do
      out%rh_seam_gap = 0.0_wp                                       ! the caller sets it against the fast loop's Rh
      !--- Lignin balance dL_s == lignin_in_s - (1-surv_s)*L_s, zero for either solver.
      lig_resid = 0.0_wp
      do i = 1_ik, 2_ik
         lig_resid = max(lig_resid, abs( (lig1(i) - lig0(i)) - (lignin_in(i) - (1.0_wp - surv(i)) * lig0(i)) ))
      end do
      out%lignin_resid = lig_resid
   end subroutine soil_carbon_step

   !---~---
   !   Semi-analytic spin-up (SASU; docs/science/soil_carbon.md eq. (7)): the steady state
   ! under climatological-mean drivers, X_ss = -(A*xi_bar*K)^{-1} * u_bar, solved on the
   ! active sub-block {j: K_j > 0} only, since inert pools (K = 0) make the full system
   ! singular. Inert pools are held at zero.
   !---~---
   pure subroutine solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, pools_ss)
      real(wp),                  intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)
      real(wp),                  intent(in)  :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL)
      type(soil_carbon_state_t), intent(out) :: pools_ss
      real(wp)    :: x_ss(N_SOIL_POOL)
      call steady_state_vector(a_mat, k_diag, xi_bar, u_bar, x_ss)
      pools_ss = soil_carbon_state_t()                       ! zero all (inert pools + lignin stay 0)
      pools_ss%carbon = x_ss
   end subroutine solve_soil_carbon_steady_state

   !---~---
   !   The spin-up's steady pools: SASU for the configured climatological xi and litter input
   ! ([soil_carbon].spinup_*), with the transfer matrix of empty pools (lignin fraction 0).
   !---~---
   function soil_carbon_steady_pools(params) result(pools)
      type(soil_carbon_params_t), intent(in) :: params
      type(soil_carbon_state_t) :: pools
      type(litter_input_t) :: lit
      real(wp) :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL), lignin_bar(2)
      call decomposition_transfer_matrix(soil_carbon_state_t(), params, a_mat, k_diag, er)
      xi_bar = params%spinup_xi
      lit%labile_grnd = params%spinup_labile_grnd
      lit%labile_soil = params%spinup_labile_soil
      lit%struct_grnd = params%spinup_struct_grnd
      lit%struct_soil = params%spinup_struct_soil
      call build_litter_input(lit, u_bar, lignin_bar)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, pools)
   end function soil_carbon_steady_pools

   !--- Solve M[active]*x = -u[active] on the active sub-block by LU with partial pivoting.
   pure subroutine steady_state_vector(a_mat, k_diag, xi_bar, u_bar, x_ss)
      real(wp), intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)
      real(wp), intent(in)  :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL)
      real(wp), intent(out) :: x_ss(N_SOIL_POOL)
      real(wp)    :: xik(N_SOIL_POOL), m_red(N_SOIL_POOL, N_SOIL_POOL), b_red(N_SOIL_POOL)
      real(wp)    :: x_red(N_SOIL_POOL)
      integer(ik) :: idx(N_SOIL_POOL), piv(N_SOIL_POOL), na, i, j
      logical     :: ok

      xik   = xi_bar * k_diag                             ! diag(xi*K)
      m_red = 0.0_wp
      b_red = 0.0_wp
      idx   = 0_ik
      !--- active index set {j: K_j > 0}.
      na = 0_ik
      do j = 1_ik, N_SOIL_POOL
         if (k_diag(j) > tiny_num) then
            na = na + 1_ik
            idx(na) = j
         end if
      end do
      !--- reduced block M = (A*diag(xi*K))[active,active], RHS b = -(u)[active].
      do j = 1_ik, na
         b_red(j) = -u_bar(idx(j))
         do i = 1_ik, na
            m_red(i, j) = a_mat(idx(i), idx(j)) * xik(idx(j))
         end do
      end do
      x_red = b_red
      call lu_factor(m_red(1:na, 1:na), piv(1:na), na, ok)
      if (.not. ok) error stop 'ed_soil_carbon: the steady-state soil carbon system is singular'
      call lu_solve(m_red(1:na, 1:na), piv(1:na), na, x_red(1:na))
      !--- scatter back (inert pools stay 0).
      x_ss = 0.0_wp
      do j = 1_ik, na
         x_ss(idx(j)) = x_red(j)
      end do
   end subroutine steady_state_vector

   !---~---
   !   The exact daily step: with M = A*diag(xi_int*K) and the day's input u held constant,
   ! [X;1] evolves as exp(Aug)*[X0;1], Aug = [[M, u],[0, 0]], which carries the input
   ! without inverting M, singular when a pool is inert (docs/science/soil_carbon.md eq. (6)).
   !---~---
   pure subroutine step_expm(a_mat, dvec, u, x0, x1)
      real(wp), intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL)
      real(wp), intent(in)  :: dvec(N_SOIL_POOL)          !--- day-integrated decay fraction xi_int*K [-]
      real(wp), intent(in)  :: u(N_SOIL_POOL), x0(N_SOIL_POOL)
      real(wp), intent(out) :: x1(N_SOIL_POOL)
      integer(ik), parameter :: M = N_SOIL_POOL + 1_ik
      real(wp)    :: aug(M, M), e(M, M), y0(M), y1(M)
      integer(ik) :: i, j

      !--- assemble Aug = [[A*diag(dvec), u],[0,0]].
      aug = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         do i = 1_ik, N_SOIL_POOL
            aug(i, j) = a_mat(i, j) * dvec(j)
         end do
      end do
      do i = 1_ik, N_SOIL_POOL
         aug(i, M) = u(i)
      end do
      call matrix_exp(aug, e, M)
      !--- y1 = e * [x0; 1].
      y0(1:N_SOIL_POOL) = x0
      y0(M)             = 1.0_wp
      do i = 1_ik, M
         y1(i) = 0.0_wp
         do j = 1_ik, M
            y1(i) = y1(i) + e(i, j) * y0(j)
         end do
      end do
      x1 = y1(1:N_SOIL_POOL)
   end subroutine step_expm

   !---~---
   !   0 when every pool is physically possible, else the index of the first that is not:
   ! not-a-number, negative beyond round-off, or above the divergence ceiling. An index, so
   ! the caller's message can name the pool, which points to the operator that broke it.
   ! Cheap enough to call at every operator boundary.
   !---~---
   pure function soil_carbon_bad_pool(pools) result(k)
      type(soil_carbon_state_t), intent(in) :: pools
      integer(ik)                     :: k
      real(wp)    :: x(N_SOIL_POOL)
      integer(ik) :: j
      x = pools%carbon
      k = 0_ik
      do j = 1_ik, N_SOIL_POOL
         !---~---
         !   NaN first, and alone, by classification: a comparison with a NaN can raise IEEE
         ! invalid (ordered ones always, and `x /= x` under ifx at -O0), which a trapping
         ! build stops on, and `.or.` does not short-circuit.
         !---~---
         if (ieee_is_nan(x(j))) then
            k = j
            return
         end if
         if (x(j) < -SOILC_NEG_TOL .or. x(j) > SOILC_MAX_PLAUSIBLE) then
            k = j
            return
         end if
      end do
   end function soil_carbon_bad_pool

end module ed_soil_carbon
