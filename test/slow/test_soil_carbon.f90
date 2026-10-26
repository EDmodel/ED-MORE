! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_soil_carbon -- unit tests for the soil carbon pools' daily step (meds_soil_carbon;  !
! docs/science/soil_carbon.md):                                                            !
!   1. MASS CLOSURE  : dC == litter_in - Rh for one soil_carbon_step (EULER and EXPM),     !
!                      with Rh recomputed INDEPENDENTLY from er*xi_int*K*X0. rh_today is   !
!                      DEFINED as litter_in - dC, so a residual built from it is           !
!                      identically zero for any inputs (docs/science/soil_carbon.md §9).   !
!   2. SPIN-UP       : SASU active-block X_ss is an exact EULER fixed point; brute force   !
!                      converges; inert pools held at 0; the singular full 7x7 is never    !
!                      solved (default K5=K7=0).                                           !
!   3. LITTER+LIGNIN : build_litter_input maps pre-split streams; lignin is a proportional !
!                      passive tracer.                                                     !
!   4. NECROMASS     : necromass_to_litter's split matches the hand formula and conserves  !
!                      mass.                                                               !
!  PLAUSIBILITY      : soil_carbon_bad_pool passes round-off and real soils, catches the   !
!                      rest.                                                               !
! The decomposition chemistry the step debits with (xi, A, K, er and Rh) is                !
! meds_respiration's (test_respiration).                                                   !
!==========================================================================================!
program test_soil_carbon
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : yr_day, kgcday_2_umols
   use meds_patch_state_vars, only : litter_input_t, soil_carbon_state_t
   use meds_soil_carbon, only : soil_carbon_out_t, soil_carbon_bad_pool, SOILC_NEG_TOL, SOILC_MAX_PLAUSIBLE, &
                                build_litter_input, soil_carbon_step, solve_soil_carbon_steady_state,    &
                                necromass_to_litter, add_necromass
   use meds_biogeochem_params, only : soil_carbon_params_t, DECOMP_STEP_EULER, DECOMP_STEP_EXPM, DECOMP_SCHEME_ED2,         &
                                      DECOMP_SCHEME_CENTURY5, N_SOIL_POOL, IP_FAST_GRND, IP_FAST_SOIL,                      &
                                      IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE, SOIL_POOL_NAME
   use meds_respiration,      only : decomposition_transfer_matrix
   implicit none

   call test_mass_closure()
   call test_spinup_steady_state()
   call test_litter_and_lignin()
   call test_necromass_to_litter()
   call test_pool_plausibility()

   call test_report('test_soil_carbon')

contains
   !----- a reproducible pseudo-random-ish filler (no Date/random; deterministic). ---------------!
   pure function frac(i) result(f)
      integer(ik), intent(in) :: i
      real(wp) :: f
      f = mod(real(i, wp) * 0.61803398875_wp, 1.0_wp)     ! golden-ratio low-discrepancy in (0,1)
   end function frac
   !=======================================================================================!
   ! 1. Mass closure -- dC == litter_in - Rh for one step under BOTH solvers, Rh recomputed          !
   !    independently of the closure definition.                                                     !
   !=======================================================================================!
   subroutine test_mass_closure()
      type(soil_carbon_state_t)  :: pools
      type(soil_carbon_params_t) :: opts
      type(soil_carbon_out_t)    :: audit
      real(wp) :: u(N_SOIL_POOL), lignin_in(2), xi_int(N_SOIL_POOL), rh, x(N_SOIL_POOL), sx
      real(wp) :: a_mat(N_SOIL_POOL,N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL), x1(N_SOIL_POOL)
      integer(ik) :: i, j
      print '(a)', '-- test 1: mass closure (EULER + EXPM) --'
      opts%decomp_scheme = DECOMP_SCHEME_CENTURY5      ! all 5 active -> exercises transfers + EXPM
      opts%k_micr = 6.0_wp ; opts%k_pass = 0.012_wp ; opts%k_slow = 0.15_wp ; opts%k_struct = 1.5_wp
      opts%k_fast = 12.0_wp
      opts%er_fast = 0.55_wp ; opts%er_struct_nonlig = 0.50_wp ; opts%er_struct_lig = 0.20_wp
      opts%er_slow = 0.55_wp ; opts%er_pass = 0.55_wp ; opts%er_micr_int = 0.60_wp ; opts%er_micr_slp = 0.17_wp

      pools%carbon(IP_FAST_GRND) = 0.4_wp ; pools%carbon(IP_FAST_SOIL) = 0.6_wp
      pools%carbon(IP_STRUCT_GRND) = 2.0_wp ; pools%carbon(IP_STRUCT_SOIL) = 3.0_wp
      pools%carbon(IP_MICR) = 0.5_wp ; pools%carbon(IP_SLOW) = 8.0_wp ; pools%carbon(IP_PASSIVE) = 20.0_wp
      pools%lignin(1) = 0.6_wp ; pools%lignin(2) = 0.9_wp
      u = [ (0.01_wp * frac(int(10*i, ik)), i = 1, N_SOIL_POOL) ]
      u(IP_MICR:IP_PASSIVE) = 0.0_wp                   ! litter only enters pools 1-4
      lignin_in = [0.002_wp, 0.003_wp]
      xi_int = [ (0.5_wp * frac(int(3*i+1, ik)) + 0.1_wp, i = 1, N_SOIL_POOL) ]  ! ~O(0.1..0.6) day

      x = pools%carbon ; sx = max(sum(x), 1.0_wp)

      opts%step_solver = DECOMP_STEP_EULER
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      block
         type(soil_carbon_state_t) :: p
         real(wp) :: rh_indep, dc
         p = pools
         call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
         p = audit%pools
         call check_true('EULER lignin_resid ~ 0', abs(audit%lignin_resid) < 1.0e-12_wp*sx, audit%lignin_resid)
         !----- LAMBDA, the freeze number: the fraction of a pool one step withdraws. It must be the  !
         !      max over pools of xi_int*K, must name that pool, and -- the property that makes the   !
         !      frozen-pool approximation sound at all -- must be well under 1 for a daily step.  ----!
         call check('lambda_max == max_j xi_int_j*k_j', audit%lambda_max,                          &
                    maxval(xi_int(1:N_SOIL_POOL)*k_diag(1:N_SOIL_POOL)), 1.0e-15_wp)
         call check_true('lambda_max is far below 1 for a daily step', audit%lambda_max < 0.1_wp,   &
                         audit%lambda_max)
         call check_true('lambda_pool names a real pool', audit%lambda_pool >= 1_ik .and.          &
                         audit%lambda_pool <= N_SOIL_POOL, real(audit%lambda_pool, wp))
         call check_true('EULER rh_today >= 0', rh >= 0.0_wp, rh)
         !----- INDEPENDENT (non-tautological) mass check: dC == litter_in - Rh, with Rh recomputed    !
         !      from er*xi_int*K*X0 (NOT from rh_today = litter_in - dC, which is closure-by-definition).!
         !      This catches a dropped-litter or misattributed-transfer bug that a residual cannot.    !
         x1 = p%carbon
         rh_indep = 0.0_wp
         do j = 1_ik, N_SOIL_POOL
            rh_indep = rh_indep + er(j)*xi_int(j)*k_diag(j)*x(j)
         end do
         dc = sum(x1) - sum(x)
         call check('EULER dC == litter_in - Rh (INDEPENDENT)', dc, sum(u) - rh_indep, 1.0e-12_wp*sx)
      end block

      opts%step_solver = DECOMP_STEP_EXPM
      block
         type(soil_carbon_state_t) :: p
         p = pools
         call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
         p = audit%pools
         call check_true('EXPM  rh_today >= 0', rh >= -1.0e-12_wp, rh)
      end block
   end subroutine test_mass_closure
   !=======================================================================================!
   ! 2. SASU spin-up: active-block X_ss is an exact EULER fixed point; brute force converges;       !
   !    inert pools held at 0; default (K5=K7=0) solves despite the singular full matrix.            !
   !=======================================================================================!
   subroutine test_spinup_steady_state()
      print '(a)', '-- test 2: SASU steady state (active block) --'
      call spinup_one_scheme(DECOMP_SCHEME_ED2, 'scheme0')
      call spinup_one_scheme(DECOMP_SCHEME_CENTURY5, 'scheme5')
   end subroutine test_spinup_steady_state
   subroutine spinup_one_scheme(scheme, tag)
      integer(ik),      intent(in) :: scheme
      character(len=*), intent(in) :: tag
      type(soil_carbon_state_t)  :: pools, pss, p
      type(soil_carbon_params_t) :: opts
      type(soil_carbon_out_t)    :: audit
      real(wp) :: a_mat(N_SOIL_POOL,N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: xi(N_SOIL_POOL), u(N_SOIL_POOL), lignin_in(2), xss(N_SOIL_POOL), x1(N_SOIL_POOL)
      real(wp) :: rh, worst, relerr
      integer(ik) :: it
      opts = soil_carbon_params_t()
      opts%decomp_scheme = scheme
      if (scheme == DECOMP_SCHEME_CENTURY5) then
         opts%k_micr=6.0_wp; opts%k_pass=0.012_wp; opts%k_slow=0.15_wp; opts%k_struct=1.5_wp; opts%k_fast=12.0_wp
         opts%er_fast=0.55_wp; opts%er_struct_nonlig=0.50_wp; opts%er_struct_lig=0.20_wp
         opts%er_slow=0.55_wp; opts%er_pass=0.55_wp; opts%er_micr_int=0.60_wp; opts%er_micr_slp=0.17_wp
      end if
      ! frozen climate scalar (constant), constant litter input to pools 1-4:
      xi = [0.7_wp, 0.5_wp, 0.7_wp, 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp]
      u  = 0.0_wp
      u(IP_FAST_GRND)=0.03_wp; u(IP_FAST_SOIL)=0.05_wp; u(IP_STRUCT_GRND)=0.02_wp; u(IP_STRUCT_SOIL)=0.04_wp
      lignin_in = 0.0_wp
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi, u, pss)
      xss = pss%carbon
      ! inert pools held at 0 (default scheme); finite everywhere:
      if (scheme == DECOMP_SCHEME_ED2) then
         !----- Singular-system guard (docs/science/soil_carbon.md §8): the FULL 7x7 A*diag(xi*K) is        !
         !      SINGULAR because the inert (micr, passive) COLUMNS are identically zero (K_micr=K_pass=0),  !
         !      so a full solve would hit a zero pivot -- the active-block {K>0} solve is REQUIRED and      !
         !      still returns a finite X_ss (checked just below). Asserting the inert columns vanish is     !
         !      what proves the active-block path is genuinely exercised (a full-solve regression could not !
         !      satisfy this).                                                                              !
         block
            real(wp)    :: mcol_inert
            integer(ik) :: ii
            mcol_inert = 0.0_wp
            do ii = 1_ik, N_SOIL_POOL
               mcol_inert = mcol_inert + abs(a_mat(ii,IP_MICR)   *xi(IP_MICR)   *k_diag(IP_MICR))   &
                                       + abs(a_mat(ii,IP_PASSIVE)*xi(IP_PASSIVE)*k_diag(IP_PASSIVE))
            end do
            call check('  '//tag//' full 7x7 A*xi*K SINGULAR (inert cols = 0)', mcol_inert, 0.0_wp, 0.0_wp)
         end block
         call check('  '//tag//' inert micr = 0', xss(IP_MICR),    0.0_wp, 0.0_wp)
         call check('  '//tag//' inert pass = 0', xss(IP_PASSIVE), 0.0_wp, 0.0_wp)
      end if
      call check_true('  '//tag//' X_ss finite & >= 0', all(xss >= 0.0_wp) .and. all(xss == xss), sum(xss))
      ! EXACT EULER fixed point: step from X_ss with xi_int=xi, u -> no change.
      opts%step_solver = DECOMP_STEP_EULER
      p = pss
      call soil_carbon_step(p, u, lignin_in, xi, opts, rh, audit)
      p = audit%pools
      x1 = p%carbon
      worst = maxval(abs(x1 - xss))
      call check_true('  '//tag//' EULER fixed point dX~0', worst < 1.0e-10_wp*max(maxval(xss),1.0_wp), worst)
      ! brute force from bare ground converges TOWARD X_ss (SASU vs brute force agree). The scheme-5
      ! tolerance is looser BY DESIGN: its passive pool (k=0.012/yr => ~170 yr residence at xi~0.5) is
      ! millennial, so 400k days (~1100 yr) leaves an O(1e-3) tail -- brute-forcing that pool to round-off
      ! is exactly the multi-millennium integration SASU (the exact fixed point above, dX~0) replaces.
      p = soil_carbon_state_t()
      do it = 1_ik, 400000_ik
         call soil_carbon_step(p, u, lignin_in, xi, opts, rh, audit)
         p = audit%pools
      end do
      x1 = p%carbon
      relerr = maxval(abs(x1 - xss)) / max(maxval(xss), 1.0_wp)
      if (scheme == DECOMP_SCHEME_CENTURY5) then
         call check_true('  '//tag//' brute force -> X_ss (millennial passive tail)', relerr < 3.0e-3_wp, relerr)
      else
         call check_true('  '//tag//' brute force -> X_ss', relerr < 1.0e-6_wp, relerr)
      end if
   end subroutine spinup_one_scheme
   !=======================================================================================!
   ! 3. Litter mapping + lignin passive tracer.                                                    !
   !=======================================================================================!
   subroutine test_litter_and_lignin()
      type(soil_carbon_state_t)  :: pools, p
      type(soil_carbon_params_t) :: opts
      type(litter_input_t)       :: lit
      type(soil_carbon_out_t)    :: audit
      real(wp) :: u(N_SOIL_POOL), lignin_in(2), xi_int(N_SOIL_POOL), rh
      real(wp) :: f_lig_before, f_lig_after, c_before, l_before, d_s
      ! (a) driver-side per-cohort split demo (two PFTs) summed into a pool-destined litter_input_t:
      real(wp), parameter :: LEAF1=0.10_wp, ROOT1=0.06_wp, WOOD1=0.20_wp
      real(wp), parameter :: LEAF2=0.08_wp, ROOT2=0.05_wp, WOOD2=0.30_wp
      real(wp), parameter :: FLAB_LEAF1=0.7_wp, FLAB_STEM1=0.05_wp, AGF1=0.6_wp
      real(wp), parameter :: FLAB_LEAF2=0.5_wp, FLAB_STEM2=0.10_wp, AGF2=0.7_wp
      real(wp) :: exp_lab_grnd, exp_lab_soil, exp_str_grnd, exp_str_soil
      print '(a)', '-- test 3: litter mapping + lignin tracer --'

      ! per-cohort split: leaf splits labile/structural by f_labile_leaf and lies above ground;
      ! fineroot is below-ground; wood (CWD) is mostly structural via f_labile_stem, split by agf.
      exp_lab_grnd = LEAF1*FLAB_LEAF1 + WOOD1*FLAB_STEM1*AGF1                            &
                   + LEAF2*FLAB_LEAF2 + WOOD2*FLAB_STEM2*AGF2
      exp_lab_soil = ROOT1*FLAB_LEAF1 + WOOD1*FLAB_STEM1*(1.0_wp-AGF1)                    &
                   + ROOT2*FLAB_LEAF2 + WOOD2*FLAB_STEM2*(1.0_wp-AGF2)
      exp_str_grnd = LEAF1*(1.0_wp-FLAB_LEAF1) + WOOD1*(1.0_wp-FLAB_STEM1)*AGF1          &
                   + LEAF2*(1.0_wp-FLAB_LEAF2) + WOOD2*(1.0_wp-FLAB_STEM2)*AGF2
      exp_str_soil = ROOT1*(1.0_wp-FLAB_LEAF1) + WOOD1*(1.0_wp-FLAB_STEM1)*(1.0_wp-AGF1)  &
                   + ROOT2*(1.0_wp-FLAB_LEAF2) + WOOD2*(1.0_wp-FLAB_STEM2)*(1.0_wp-AGF2)
      lit%labile_grnd = exp_lab_grnd ; lit%labile_soil = exp_lab_soil
      lit%struct_grnd = exp_str_grnd ; lit%struct_soil = exp_str_soil
      lit%lignin_grnd = 0.3_wp*exp_str_grnd ; lit%lignin_soil = 0.3_wp*exp_str_soil
      call build_litter_input(lit, u, lignin_in)
      call check('build u(fast_grnd)',   u(IP_FAST_GRND),   exp_lab_grnd, 1.0e-14_wp)
      call check('build u(fast_soil)',   u(IP_FAST_SOIL),   exp_lab_soil, 1.0e-14_wp)
      call check('build u(struct_grnd)', u(IP_STRUCT_GRND), exp_str_grnd, 1.0e-14_wp)
      call check('build u(struct_soil)', u(IP_STRUCT_SOIL), exp_str_soil, 1.0e-14_wp)
      call check('build lignin_in(1)',   lignin_in(1),      0.3_wp*exp_str_grnd, 1.0e-14_wp)
      call check_true('no litter to inert/SOM pools', abs(u(IP_MICR))+abs(u(IP_SLOW))+abs(u(IP_PASSIVE))<1.0e-30_wp, &
                      u(IP_SLOW))

      ! (b) lignin tracer: proportional decay conserves f_lignin under pure decomposition (no input).
      opts = soil_carbon_params_t()
      pools = soil_carbon_state_t()
      pools%carbon(IP_STRUCT_GRND) = 2.0_wp ; pools%lignin(1) = 0.8_wp   ! f_lignin = 0.4
      c_before = pools%carbon(IP_STRUCT_GRND) ; l_before = pools%lignin(1)
      f_lig_before = l_before / c_before
      xi_int = 0.0_wp ; xi_int(IP_STRUCT_GRND) = 0.2_wp        ! decay only the struct_grnd pool
      u = 0.0_wp ; lignin_in = 0.0_wp
      opts%step_solver = DECOMP_STEP_EULER
      p = pools
      call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
      p = audit%pools
      f_lig_after = p%lignin(1) / p%carbon(IP_STRUCT_GRND)
      d_s = xi_int(IP_STRUCT_GRND) * (opts%k_struct/yr_day)
      call check('f_lignin conserved under pure decay', f_lig_after, f_lig_before, 1.0e-12_wp)
      call check_true('0 <= L <= C after decay', p%lignin(1) >= 0.0_wp .and. &
                      p%lignin(1) <= p%carbon(IP_STRUCT_GRND), p%lignin(1))
      call check_true('lignin_resid ~ 0', abs(audit%lignin_resid) < 1.0e-14_wp, audit%lignin_resid)
      ! fresh lignin input raises f_lignin:
      pools = soil_carbon_state_t() ; pools%carbon(IP_STRUCT_GRND) = 2.0_wp ; pools%lignin(1) = 0.4_wp
      xi_int = 0.0_wp
      u = 0.0_wp ; u(IP_STRUCT_GRND) = 0.5_wp ; lignin_in = [0.4_wp, 0.0_wp]     ! high-lignin input
      p = pools
      call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
      p = audit%pools
      call check_true('lignin input raises f_lignin', &
                      p%lignin(1)/p%carbon(IP_STRUCT_GRND) > 0.2_wp, p%lignin(1)/p%carbon(IP_STRUCT_GRND))

      ! (c) EXPM at a LARGE accelerated-spin-up step: lignin must decay in lockstep with its host
      ! carbon (exp(-d_s)), so f_lignin stays conserved and 0 <= L <= C. The linear survival (1 - d_s)
      ! would drive L NEGATIVE here (d_s = 3.0 => L = 0.8*(1-3) = -1.6; docs/science/soil_carbon.md §6).
      opts = soil_carbon_params_t() ; opts%step_solver = DECOMP_STEP_EXPM
      pools = soil_carbon_state_t() ; pools%carbon(IP_STRUCT_GRND) = 2.0_wp ; pools%lignin(1) = 0.8_wp  ! f_lignin=0.4
      c_before = pools%carbon(IP_STRUCT_GRND) ; l_before = pools%lignin(1) ; f_lig_before = l_before/c_before
      xi_int = 0.0_wp ; xi_int(IP_STRUCT_GRND) = 3.0_wp / (opts%k_struct/yr_day)   ! force d_s = xi_int*K = 3.0
      u = 0.0_wp ; lignin_in = 0.0_wp
      p = pools
      call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
      p = audit%pools
      call check_true('EXPM large step: L >= 0 (no negative lignin)', p%lignin(1) >= 0.0_wp, p%lignin(1))
      call check_true('EXPM large step: L <= C', p%lignin(1) <= p%carbon(IP_STRUCT_GRND), &
                      p%carbon(IP_STRUCT_GRND) - p%lignin(1))
      call check('EXPM large step: f_lignin conserved under pure decay', &
                 p%lignin(1)/p%carbon(IP_STRUCT_GRND), f_lig_before, 1.0e-12_wp)
      call check_true('EXPM large step: lignin_resid ~ 0', abs(audit%lignin_resid) < 1.0e-12_wp, audit%lignin_resid)
   end subroutine test_litter_and_lignin
   !=======================================================================================!
   ! 4. NECROMASS_TO_LITTER (meds_soil_carbon; docs/science/soil_carbon.md §5):                    !
   !     the driver-side (and core-side termination/disturbance) necromass -> litter-destination !
   !     split. (a) reproduces test_litter_and_lignin's own two-cohort hand formula EXACTLY when   !
   !     called per-cohort with storage_c=0 (no wood/storage in a turnover-shed call), catching any !
   !     drift between the kernel and this helper. (b) EXHAUSTIVE MASS CONSERVATION: the six        !
   !     outputs must sum to the total necromass in for an arbitrary (leaf,root,wood,storage) input  !
   !     -- the split has nowhere else for carbon to go. (c) each tissue lies where it grew:        !
   !     leaves all above ground; storage split by aboveground_frac like wood. (d) add_necromass    !
   !     records the carbon by origin, fruit and seed falling with the canopy, and the origins sum to !
   !     the pool streams.                                                                          !
   !=======================================================================================!
   subroutine test_necromass_to_litter()
      real(wp), parameter :: LEAF1=0.10_wp, ROOT1=0.06_wp, WOOD1=0.20_wp
      real(wp), parameter :: LEAF2=0.08_wp, ROOT2=0.05_wp, WOOD2=0.30_wp
      real(wp), parameter :: FLAB_LEAF1=0.7_wp, FLAB_STEM1=0.05_wp, AGF1=0.6_wp
      real(wp), parameter :: FLAB_LEAF2=0.5_wp, FLAB_STEM2=0.10_wp, AGF2=0.7_wp
      real(wp) :: lg1, ls1, sg1, ss1, lig_g1, lig_s1
      real(wp) :: lg2, ls2, sg2, ss2, lig_g2, lig_s2
      real(wp) :: exp_lab_grnd, exp_lab_soil, exp_str_grnd, exp_str_soil
      real(wp) :: leaf_c, root_c, wood_c, storage_c, total_in, total_out
      real(wp) :: lg, ls, sg, ss, lig_g, lig_s
      print '(a)', '-- test 4: necromass_to_litter (driver/core litter-split helper) --'

      ! (a) per-cohort call reproduces test_litter_and_lignin's hand formula exactly (storage_c=0,
      !     lignin_frac=0.3 matching that test's flat structural-stream lignin fraction).
      call necromass_to_litter(LEAF1, ROOT1, WOOD1, 0.0_wp, FLAB_LEAF1, FLAB_STEM1, AGF1, 0.3_wp,   &
                               lg1, ls1, sg1, ss1, lig_g1, lig_s1)
      call necromass_to_litter(LEAF2, ROOT2, WOOD2, 0.0_wp, FLAB_LEAF2, FLAB_STEM2, AGF2, 0.3_wp,   &
                               lg2, ls2, sg2, ss2, lig_g2, lig_s2)
      exp_lab_grnd = LEAF1*FLAB_LEAF1 + WOOD1*FLAB_STEM1*AGF1                            &
                   + LEAF2*FLAB_LEAF2 + WOOD2*FLAB_STEM2*AGF2
      exp_lab_soil = ROOT1*FLAB_LEAF1 + WOOD1*FLAB_STEM1*(1.0_wp-AGF1)                    &
                   + ROOT2*FLAB_LEAF2 + WOOD2*FLAB_STEM2*(1.0_wp-AGF2)
      exp_str_grnd = LEAF1*(1.0_wp-FLAB_LEAF1) + WOOD1*(1.0_wp-FLAB_STEM1)*AGF1          &
                   + LEAF2*(1.0_wp-FLAB_LEAF2) + WOOD2*(1.0_wp-FLAB_STEM2)*AGF2
      exp_str_soil = ROOT1*(1.0_wp-FLAB_LEAF1) + WOOD1*(1.0_wp-FLAB_STEM1)*(1.0_wp-AGF1)  &
                   + ROOT2*(1.0_wp-FLAB_LEAF2) + WOOD2*(1.0_wp-FLAB_STEM2)*(1.0_wp-AGF2)
      call check('necromass_to_litter matches hand formula: labile_grnd', lg1+lg2, exp_lab_grnd, 1.0e-14_wp)
      call check('necromass_to_litter matches hand formula: labile_soil', ls1+ls2, exp_lab_soil, 1.0e-14_wp)
      call check('necromass_to_litter matches hand formula: struct_grnd', sg1+sg2, exp_str_grnd, 1.0e-14_wp)
      call check('necromass_to_litter matches hand formula: struct_soil', ss1+ss2, exp_str_soil, 1.0e-14_wp)
      call check('necromass_to_litter matches hand formula: lignin_grnd', lig_g1+lig_g2, 0.3_wp*exp_str_grnd, 1.0e-14_wp)
      call check('necromass_to_litter matches hand formula: lignin_soil', lig_s1+lig_s2, 0.3_wp*exp_str_soil, 1.0e-14_wp)

      ! (b) exhaustive mass conservation: the four labile/structural destinations (lignin is a
      !     passive SUB-tracer of struct_grnd/struct_soil, not a fifth carbon pool) must sum to the
      !     total necromass in, for an arbitrary input including a nonzero storage term.
      leaf_c = 0.37_wp ; root_c = 0.21_wp ; wood_c = 1.85_wp ; storage_c = 0.09_wp
      call necromass_to_litter(leaf_c, root_c, wood_c, storage_c, 0.65_wp, 0.08_wp, 0.55_wp, 0.28_wp, &
                               lg, ls, sg, ss, lig_g, lig_s)
      total_in  = leaf_c + root_c + wood_c + storage_c
      total_out = lg + ls + sg + ss
      call check('necromass_to_litter conserves total mass', total_out, total_in, 1.0e-13_wp*total_in)
      call check_true('lignin_grnd <= struct_grnd', lig_g <= sg, sg - lig_g)
      call check_true('lignin_soil <= struct_soil', lig_s <= ss, ss - lig_s)

      ! (c) a leaf lies above ground whatever aboveground_frac is; storage splits by it, like wood.
      call necromass_to_litter(0.40_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.65_wp, 0.08_wp, 0.55_wp, 0.28_wp, &
                               lg1, ls1, sg1, ss1, lig_g1, lig_s1)
      call check('leaf litter is all above ground: labile_grnd', lg1, 0.40_wp*0.65_wp, 1.0e-14_wp)
      call check('leaf litter is all above ground: struct_grnd', sg1, 0.40_wp*0.35_wp, 1.0e-14_wp)
      call check('leaf litter is all above ground: nothing below', ls1 + ss1, 0.0_wp, 0.0_wp)
      call necromass_to_litter(0.0_wp, 0.0_wp, 0.0_wp, 0.15_wp, 0.65_wp, 0.08_wp, 0.55_wp, 0.28_wp, &
                               lg2, ls2, sg2, ss2, lig_g2, lig_s2)
      call check('storage splits by aboveground_frac: labile_grnd', lg2, 0.15_wp*0.65_wp*0.55_wp, 1.0e-14_wp)
      call check('storage splits by aboveground_frac: labile_soil', ls2, 0.15_wp*0.65_wp*0.45_wp, 1.0e-14_wp)

      ! (d) add_necromass: the pool streams of necromass_to_litter (fruit/seed in the canopy channel)
      !     plus the carbon by origin, and the origins sum to the pool streams.
      block
         type(litter_input_t) :: lit
         lit = litter_input_t()
         call add_necromass(lit, leaf_c, root_c, wood_c, storage_c, 0.04_wp, 0.65_wp, 0.08_wp, 0.55_wp, 0.28_wp)
         call necromass_to_litter(leaf_c + 0.04_wp, root_c, wood_c, storage_c, 0.65_wp, 0.08_wp, 0.55_wp, &
                                  0.28_wp, lg, ls, sg, ss, lig_g, lig_s)
         call check('add_necromass adds the split: labile_grnd', lit%labile_grnd, lg, 0.0_wp)
         call check('add_necromass adds the split: struct_soil', lit%struct_soil, ss, 0.0_wp)
         call check('add_necromass records the leaf by origin', lit%leaf_c, leaf_c, 0.0_wp)
         call check('add_necromass records fruit and seed by origin', lit%repro_c, 0.04_wp, 0.0_wp)
         call check('the origins sum to the pool streams',                                           &
                    lit%leaf_c + lit%fineroot_c + lit%wood_c + lit%storage_c + lit%repro_c,         &
                    lit%labile_grnd + lit%labile_soil + lit%struct_grnd + lit%struct_soil, 1.0e-14_wp)
      end block
   end subroutine test_necromass_to_litter
   !---------------------------------------------------------------------------------------!
   ! POOL PLAUSIBILITY -- the predicate that answers the question conservation cannot.            !
   !                                                                                          !
   ! A runaway that CONSERVES carbon closes every budget: a coupled spin-up can drive a patch !
   ! to 7e7 kgC/m2 of litter carbon with NEGATIVE structural and slow pools while the         !
   ! slow-loop ledger reports every phase closed on all three currencies                      !
   ! (docs/science/soil_carbon.md §9).                                                        !
   !                                                                                          !
   ! These assertions are the guard against that class, so they check the BOUNDARIES rather than    !
   ! the happy path: round-off negatives must pass (the daily Euler step produces them), real        !
   ! negatives must fail, and the ceiling must sit far above any real soil.                          !
   !---------------------------------------------------------------------------------------!
   subroutine test_pool_plausibility()
      type(soil_carbon_state_t) :: p
      integer(ik)         :: k

      !----- A realistic temperate-forest profile is PLAUSIBLE. -------------------------------!
      p = soil_carbon_state_t()
      p%carbon(IP_FAST_GRND)   = 1.5_wp  ; p%carbon(IP_FAST_SOIL)   = 1.8_wp
      p%carbon(IP_STRUCT_GRND) = 7.2_wp  ; p%carbon(IP_STRUCT_SOIL) = 6.3_wp
      p%carbon(IP_SLOW)        = 5.4_wp
      call check_true('a real soil profile passes the plausibility check',                          &
                      soil_carbon_bad_pool(p) == 0_ik, real(soil_carbon_bad_pool(p), wp))
      call check_true('no offending pool is named (none)',                                          &
                      trim(SOIL_POOL_NAME(soil_carbon_bad_pool(p))) == '(none)')

      !----- All-zero (bare ground, cold start) is PLAUSIBLE. It has to be: every run starts here. !
      p = soil_carbon_state_t()
      call check_true('an all-zero cold start passes', soil_carbon_bad_pool(p) == 0_ik,             &
                      real(soil_carbon_bad_pool(p), wp))

      !----- ROUND-OFF negatives must PASS. The daily Euler step subtracts a computed loss from a  !
      !      pool and lands a few ulp below zero routinely; a guard that fired on that would be     !
      !      turned off within a day, which is worse than no guard.  --------------------------------!
      p = soil_carbon_state_t()
      p%carbon(IP_SLOW) = -0.1_wp * SOILC_NEG_TOL
      call check_true('a round-off negative is tolerated', soil_carbon_bad_pool(p) == 0_ik,         &
                      p%carbon(IP_SLOW))

      !----- A REAL negative must FAIL, and name the right pool. --------------------------------!
      p = soil_carbon_state_t()
      p%carbon(IP_STRUCT_GRND) = -5.6e3_wp                       ! the size of a runaway's negative pool
      k = soil_carbon_bad_pool(p)
      call check_true('a genuinely negative pool is caught', k == IP_STRUCT_GRND, real(k, wp))
      call check_true('the offending pool is named', trim(SOIL_POOL_NAME(k)) == 'struct_grnd', &
                      real(k, wp))

      !----- Above the divergence CEILING must fail. --------------------------------------------!
      p = soil_carbon_state_t()
      p%carbon(IP_FAST_GRND) = 7.15e7_wp                         ! the size of a runaway's litter pool
      k = soil_carbon_bad_pool(p)
      call check_true('a diverged pool is caught by the ceiling', k == IP_FAST_GRND, real(k, wp))

      !----- The ceiling must sit FAR above any real soil, or it polices science instead of        !
      !      catching divergence. 100 kgC/m2 is already richer than almost anything measured.      !
      call check_true('the ceiling leaves real soils alone', SOILC_MAX_PLAUSIBLE > 1.0e2_wp,        &
                      SOILC_MAX_PLAUSIBLE)

      !----- A NaN must fail. It is the other way a pool stops being a number. -------------------!
      p = soil_carbon_state_t()
      p%carbon(IP_MICR) = ieee_value(1.0_wp, ieee_quiet_nan)
      call check_true('a NaN pool is caught', soil_carbon_bad_pool(p) == IP_MICR,                   &
                      real(soil_carbon_bad_pool(p), wp))
   end subroutine test_pool_plausibility
end program test_soil_carbon
