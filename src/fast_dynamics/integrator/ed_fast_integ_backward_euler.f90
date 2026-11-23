! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_backward_euler -- the fast loop's backward-Euler solves, two systems
! side by side. One ESDIRK2 stage (docs/science/numerical_scheme.md §3 `esdirk2`, eq. (1))
! solves the canopy air against the leaf, wood and ground temperatures by a small Newton
! whose Jacobian is the box's own terms less the surface's source slopes, taken by forward
! differences and reused over a patch's day: backward_euler_stage is the stage,
! newton_surface_solve its inner solve, surface_source_slopes and surface_jacobian the
! Jacobian's two parts, cas_box_commit the box commit; advance_plant_water and
! advance_surf_water_full move plant and surface water once per step.
! The soil-plant water network's step is implicit Euler on adaptive sub-steps
! (docs/science/soil_plant_water.md §3): water_step, its Newton sub-step, the Newton
! matrix and the sub-step's commit. The two Newton solves share only the line search and
! the step bound of ed_numerics.
!---~---
module ed_fast_integ_backward_euler
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use ed_kinds, only : wp, ik
   use ed_constants,   only : tiny_num
   use ed_soil_heat,   only : soil_heat_env_t, soil_heat_out_t, soil_energy_step_implicit, set_soil_heat_env
   use ed_ground_surface_lib, only : porous_layer_temp
   use ed_patch_state_vars,   only : soil_energy_state_t
   use ed_fast_integ_derivs, only : surface_derivs, cas_conductances, water_derivs, node_balance, wood_potentials
   use ed_fast_integ_control, only : water_scales
   use ed_canopy_air_space,   only : cas_box_step_implicit, cas_env_t
   use ed_fast_integ_types,   only : march_state_t, frozen_t, surface_state_t, cas_boundary_t,                  &
                               surface_tend_t, stage_bflux_t, error_control_t, tol_set_t, fast_static_state_t,    &
                               water_drive_t, water_tend_t, water_state_t, water_flows_t, water_work_t,           &
                               GRP_THETA, POND, WOOD, LEAF, COHORT_STORES, WORK_WATER_NSUB, WORK_WATER_NEWTON,    &
                               WORK_WATER_REJECT, WORK_WATER_NFAIL, WORK_WATER_FORCED_NEWTON,                     &
                               WORK_WATER_FORCED_FAST, WORK_WATER_FORCED_CUT, WORK_WATER_FORCED_COHORT,           &
                               WORK_WATER_NODE_UNSOLVED, WORK_WATER_FAST_POND, WORK_WATER_FAST_SLOT,              &
                               WORK_WATER_FAST_LAYER, WORK_WATER_CUT_POND, WORK_WATER_CUT_SLOT,                   &
                               WORK_WATER_CUT_LAYER, state_init, surface_source_slopes_t, water_organic_t
   use ed_site_state_vars, only : site_surface_t
   use ed_therm_lib, only : cas_enthalpy_to_temp, internal_energy_liquid, liquid_energy_to_temp,                   &
                           temp_to_internal_energy
   use ed_numerics,     only : sufficient_decrease, fraction_to_boundary
   use ed_soil_hydr_lib,   only : pond_water_to_head, soil_layer_curve_t, curve_kinks
   use ed_soil_water,      only : ground_evaporation_setup
   use ed_linear_algebra, only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   use ed_soil_params, only : soil_profile_t, soil_params_t
   use ed_pft_params, only : hydraulics_table_t
   implicit none
   private

   public :: backward_euler_stage, advance_plant_water, advance_surf_water_full
   public :: water_step
   public :: water_jacobian_error   ! used only by the tests
   public :: newton_surface_solve, surface_source_slopes, surface_jacobian, surface_residual_scaled, &
             SURFACE_NEWTON_MAX                                                        ! used only by the tests

   !---~---
   !   The canopy-air Newton (newton_surface_solve). A stage's balance counts as solved when its
   ! enthalpy and humidity residuals, as changes of the canopy air over the stage, are within
   ! these absolute tolerances plus SURFACE_RTOL of the state: 50 J/kg is about 0.05 K of
   ! canopy-air temperature. A solve takes fresh source slopes at most SURFACE_NEWTON_MAX times,
   ! once more when its last step was on held ones; it halves a step at most SURFACE_LS_MAX times
   ! and stops at SURFACE_FEVAL_MAX evaluations of the surface, which also bounds its steps on held
   ! slopes.
   !---~---
   real(wp),    parameter :: SURFACE_RTOL      = 1.0e-7_wp
   real(wp),    parameter :: SURFACE_ATOL_ENTH = 5.0e1_wp    ! [J/kg]
   real(wp),    parameter :: SURFACE_ATOL_SHV  = 1.0e-6_wp   ! [kg/kg]
   integer(ik), parameter :: SURFACE_NEWTON_MAX = 4_ik
   integer(ik), parameter :: SURFACE_LS_MAX     = 6_ik
   integer(ik), parameter :: SURFACE_FEVAL_MAX  = 24_ik
   !---~---
   !   Held source slopes are taken afresh after a step on them that converges poorly: one whose
   ! line search had to shorten it, or one that cut the scaled residual by less than
   ! SLOPES_MIN_FALL. Slopes taken at the iterate cut it by far more than that near the solution,
   ! so a smaller fall says the held ones no longer describe the surface.
   !---~---
   real(wp),    parameter :: SLOPES_MIN_FALL = 4.0_wp

   !---~---
   !   The water network's Newton. The balance counts as solved when every store's imbalance
   ! over the sub-step is below NEWTON_TOL of the step's absolute water tolerance -- far inside
   ! the error the sub-step itself is allowed. A Newton step may take a soil layer or a tissue
   ! at most STEP_MAX_SHARE of the way to empty (fraction to the boundary), a soil layer's
   ! unknown no closer than SOIL_TAU_FLOOR, and move its saturation or relative water content by
   ! at most STEP_MAX_DRY. The floor lies far below any frozen layer's melt plateau
   ! (SE_MIN x MELT_DRY_SHARE) and far above the denormals whose powers the conductivity's slope
   ! cannot take.
   !---~---
   real(wp),    parameter :: NEWTON_TOL     = 1.0e-2_wp
   integer(ik), parameter :: NEWTON_MAX     = 25_ik
   integer(ik), parameter :: BACKTRACK_MAX  = 12_ik
   real(wp),    parameter :: STEP_MAX_SHARE = 0.9_wp
   real(wp),    parameter :: STEP_MAX_DRY   = 0.2_wp
   real(wp),    parameter :: SOIL_TAU_FLOOR = 1.0e-12_wp
   !---~---
   !   The node classes (docs/science/soil_plant_water.md §3, "Node classes"). A core node whose
   ! time to relax to its neighbours, its storage over its flows' slope, is under FAST_SHARE of the
   ! sub-step is eliminated: its storage is then below the Newton tolerance's share of its flows,
   ! so a Newton step treats it as algebraic anyway. Its own balance is solved to
   ! NODE_SOLVE_SHARE of the Newton tolerance, so the slack it leaves its neighbour is negligible,
   ! in at most NODE_SOLVE_MAX steps.
   !---~---
   real(wp),    parameter :: FAST_SHARE       = NEWTON_TOL
   real(wp),    parameter :: NODE_SOLVE_SHARE = 1.0e-3_wp
   integer(ik), parameter :: NODE_SOLVE_MAX   = 100_ik
   !--- [-] The largest scaled imbalance the merit squares: far beyond any step the line search keeps.
   real(wp),    parameter :: MERIT_X_MAX      = 1.0e100_wp
   !---~---
   !   Its sub-steps: the error controller of an order-one method (exponent -1/2), its growth and
   ! shrink limits, the factor a sub-step whose Newton failed is cut by, and the shortest
   ! sub-step, as a share of the fast step, below which a sub-step is taken as it stands.
   !---~---
   real(wp),    parameter :: SAFETY = 0.9_wp, GROW_MAX = 4.0_wp, SHRINK_MIN = 0.2_wp
   real(wp),    parameter :: NEWTON_FAIL_CUT = 0.25_wp
   real(wp),    parameter :: H_MIN_SHARE = 1.0e-5_wp
   integer(ik), parameter :: TRIALS_MAX = 2000_ik
   !--- [kg/m2] Below this a pond holds no water to carry a temperature.
   real(wp),    parameter :: POND_TINY = 1.0e-9_wp

contains

   !---~---
   !   backward_euler_stage -- one backward-Euler stage Y = base + dt*f(Y) of the stiff
   ! block (each ESDIRK2 stage, dt = gamma*h, is one call): the CAS twins, implicit in the
   ! atmospheric exchange through the coupled 2x2 Newton, and soil heat (implicit,
   ! Thomas). Soil water moves at the water network's steady rate; plant water is passed
   ! through and advanced once (advance_plant_water).
   !---~---
   subroutine backward_euler_stage(y, frozen, site, ec, n, nsl, dt, slopes, y_out, surf_tend, bf)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      type(error_control_t), intent(in)  :: ec     ! its budget stop ([fast].debug_error)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(surface_source_slopes_t), intent(inout) :: slopes   ! the Newton's source slopes, reused and refreshed
      !---~---
      !   The stage's state and its own surface tendencies, the ones its CAS source and bf used; the
      ! caller b-weights transp_c into the plant water, so the leaves' water debit and the canopy
      ! air's vapour credit agree (docs/science/numerical_scheme.md §5a). Both filled in place.
      !---~---
      type(march_state_t),  intent(inout) :: y_out
      type(surface_tend_t),  intent(inout) :: surf_tend
      type(stage_bflux_t), optional, intent(out) :: bf  ! per-stage boundary-flux rates for the ESDIRK2 ledger

      type(surface_state_t)      :: y_stage
      type(cas_boundary_t)       :: cas_stage   ! the CAS boundary with this stage's conductances
      type(soil_energy_state_t) :: se
      type(soil_heat_env_t)     :: eforc
      type(soil_heat_out_t)        :: eflux
      real(wp)    :: t_ground, fliq1, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp)    :: enth1, shv1, e_infil, e_drain, t_cas1, water_heat(nsl)
      real(wp)    :: co21, enth_newton, shv_newton
      integer(ik) :: newton_steps
      logical     :: converged

      call state_init(y, n, nsl, y_out)

      !--- diagnose the soil-top temperature so the ground skin sees the current state.
      call porous_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                &
                           site%soil_profile%layer_curve(1), t_ground, fliq1)

      cas_mass_capacity = frozen%cas%box%air_mass_capacity
      cas_molar_capacity = frozen%cas%box%air_molar_capacity
      cas_stage = frozen%cas              ! plain scalars: a cheap copy, overridden with this stage's conductances below

      !---~---
      !   Re-solve the surface layer at this stage's canopy-air state, so the stage pays for the
      ! ventilation its own temperature earns (docs/science/numerical_scheme.md §2); the box commit,
      ! the Newton and the ledger all read these values. cas_stage carries them with mo_live cleared,
      ! so surface_derivs does not re-solve on each Newton residual evaluation.
      !---~---
      call cas_conductances(frozen%cas, y%cas_enthalpy, y%cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      cas_stage%box%atm_conductance_enthalpy = g_atm_heat
      cas_stage%box%atm_conductance_vapor = g_atm_vapour
      cas_stage%box%atm_conductance_co2 = g_atm_co2
      cas_stage%mo_live = .false.

      !---~---
      !   The coupled 2x2 Newton for the leaves and the canopy air; its final surf_tend drives the box
      ! commit and the soil sinks, one flux per interface. Its solution and step count are the
      ! tests'; the box commit reads surf_tend.
      !---~---
      call newton_surface_solve(y, cas_stage, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, &
                                slopes, surf_tend, enth_newton, shv_newton, newton_steps, converged)
      !---~---
      !   One box commit for all three twins, in flux form, so the budgets close for any surf_tend:
      ! CO2 rides the same box as enthalpy and humidity.
      !---~---
      call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, cas_stage, dt, enth1, shv1, co21)
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      y_out%cas_co2      = co21

      !--- Soil heat: implicit backward Euler, Thomas solve.
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      !---~---
      !   The enthalpy of the water the network moved: across the faces at its time-mean flows, in
      ! from the pond at t_infil, out with root uptake and exfiltration at each layer's start-of-step
      ! temperature (what the wood and the pond were credited), and out with drainage
      ! (docs/science/soil_plant_water.md §5). The leaf pays its transpiration's full vapour
      ! enthalpy in surface_derivs.
      !---~---
      call set_soil_heat_env(nsl, surf_tend%g_top, y%theta, frozen%roots%root_heat, frozen%hydrology%exfil_enth,  &
                             frozen%hydrology%w_flux_frozen, frozen%hydrology%infiltration,                       &
                             frozen%hydrology%t_infil, frozen%hydrology%drainage, frozen%hydrology%t_bot, eforc,  &
                             water_heat, e_drain)
      e_infil = frozen%hydrology%infiltration * internal_energy_liquid(frozen%hydrology%t_infil)
      call soil_energy_step_implicit(se, eforc, site%soil_profile, site%soil, dt, eflux)
      !---~---
      !   [fast].debug_error stops the run when the soil heat's own budget misses by more than
      ! soil_energy_budget_atol [K] of the column's heat capacity.
      !---~---
      if (ec%debug_error .and. abs(eflux%energy_resid) > ec%soil_energy_budget_atol * eflux%heat_capacity)    &
         error stop 'soil heat step: energy budget did not close'
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !---~---
      !   Soil water is the network's, its end-of-step theta committed once by update_dynamic_state
      ! (docs/science/numerical_scheme.md §5b). The stages move theta at its steady rate, so each
      ! layer holds the water whose enthalpy the faces brought in; a steady rate is integrated exactly
      ! and adds nothing to the error estimate.
      !---~---
      y_out%theta(1:nsl) = y%theta(1:nsl) + dt * frozen%hydrology%theta_rate(1:nsl)

      !--- Plant water mass is passed through (advance_plant_water advances it).
      y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)

      !---~---
      !   This stage's boundary-flux rates for the ledger, all at the committed state, so the amounts
      ! esdirk2_step and adaptive_march b-weight and accumulate telescope to closure.
      !---~---
      if (present(bf)) then
         t_cas1 = cas_enthalpy_to_temp(enth1, shv1)         ! committed CAS temperature, for the reported H
            bf%cas_enth_in  = surf_tend%src_enth + g_atm_heat*frozen%cas%box%atm_enthalpy
            bf%cas_enth_out = g_atm_heat*enth1
            bf%cas_vap_in   = surf_tend%src_vap  + g_atm_vapour*frozen%cas%box%atm_specific_humidity
            bf%cas_vap_out  = g_atm_vapour*shv1
            bf%cas_co2_in   = frozen%cas%nee_biotic + g_atm_co2*frozen%cas%box%atm_co2
            bf%cas_co2_out  = g_atm_co2*y_out%cas_co2
            !---~---
            !   The bottom face the soil solve committed (zero on the adiabatic base, the conduction
            ! to deep_temp under the Dirichlet anchor) crosses both the soil's and the column's book.
            !---~---
            bf%soil_enth_in = surf_tend%g_top - eflux%bottom_heat + e_infil
            bf%soil_enth_out= sum(water_heat(1:nsl)) + e_drain
            !---~---
            !   Soil water's storage change and flows are the network's, booked once per step by
            ! close_ledgers: the per-stage soil-water amounts are zero, and the whole-column
            ! water here is only the CAS-vapour exchange.
            !---~---
            bf%soil_wat_in  = 0.0_wp
            bf%soil_wat_out = 0.0_wp
            !---~---
            !   ground_rad is the tiles' radiation, cover-weighted (rnet_soil when bare). Infiltration,
            ! exfiltration and uptake move enthalpy between tracked stores and cancel; rain and runoff
            ! are added once at the outer level (docs/science/numerical_scheme.md §5b).
            !---~---
            bf%whole_enth_in= surf_tend%coh_rnet + surf_tend%ground_rad - eflux%bottom_heat
            !---~---
            !   The condensate is not a boundary loss: the caller deposits it in soil layer 1 with the
            ! enthalpy surface_derivs debited from the CAS (whole_cond, whole_cond_enth).
            !---~---
            bf%whole_enth_out= g_atm_heat*(enth1 - frozen%cas%box%atm_enthalpy) + e_drain
            bf%whole_wat_in = 0.0_wp
            bf%whole_wat_out = g_atm_vapour*(shv1 - frozen%cas%box%atm_specific_humidity)
            bf%whole_cond   = surf_tend%cond                     ! deposited into soil layer 1, not lost
            bf%whole_cond_enth = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
            bf%stores       = surf_tend%stores
            bf%stack_newton   = surf_tend%stack_newton
            bf%stack_fallback = surf_tend%stack_fallback
            bf%veg_surface_evap    = sum(surf_tend%leaf_surface_evap(1:n)) + sum(surf_tend%wood_surface_evap(1:n))
            bf%h_canopy     = surf_tend%h_canopy
            bf%h_ground     = surf_tend%h_ground
            !---~---
            !   The reported H, against theta_atm (T at zref), the ledger's basis, at the canopy
            ! air's specific heat.
            !---~---
            bf%atm_heat_out = g_atm_heat*frozen%cas%cp*(t_cas1 - frozen%cas%env%theta_atm)
            bf%atm_vap_out  = g_atm_vapour*(shv1  - frozen%cas%box%atm_specific_humidity)
            bf%rad_net      = surf_tend%coh_rnet + surf_tend%ground_rad
      end if
   end subroutine backward_euler_stage

   !---~---
   !   newton_surface_solve -- a direct 2x2 Newton solve of the coupled backward-Euler canopy-air block
   ! {R_H, R_q} = 0 in (cas_enthalpy, cas_shv), whose surface sources depend on both through tcas,
   ! qcas and qsat; robust near saturation, where a leaf <-> CAS Picard iteration is not. The
   ! Jacobian (surface_jacobian) is the box's own terms, exact at every stage, less the surface's
   ! source slopes (surface_source_slopes), which capture the strong VPD self-limiting d src_vap/dq.
   ! The slopes change slowly over a patch's day, so a solve steps on the set it is handed and takes
   ! one afresh at the iterate only when none is held or a step on the held set converges poorly
   ! (docs/science/numerical_scheme.md §3, "The stage Newton"). A singular-Jacobian fallback, a
   ! line search and an evaluation cap; it never stops the run. Returns surf_tend at the solution,
   ! from which backward_euler_stage commits the box, the solution itself, the steps taken and
   ! whether it converged.
   !---~---
   subroutine newton_surface_solve(y, cas, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, &
                                   slopes, surf_tend, cas_enthalpy, cas_shv, steps, ok)
      type(march_state_t),   intent(in)    :: y
      type(cas_boundary_t),   intent(in)    :: cas      ! the stage's CAS boundary (its own conductances)
      type(frozen_t),  intent(in)    :: frozen   ! tissue / veg_surface / ground / snow for surface_derivs
      real(wp),               intent(in)    :: t_ground
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, cas_mass_capacity, g_atm_heat, g_atm_vapour
      type(surface_source_slopes_t), intent(inout) :: slopes   ! reused when held, taken afresh here when not
      type(surface_tend_t),   intent(inout) :: surf_tend   ! the evaluation at the solution, on return
      real(wp),               intent(out)   :: cas_enthalpy, cas_shv   ! [J/kg], [kg/kg] the solution
      integer(ik),            intent(out)   :: steps      ! the Newton steps taken
      logical,                intent(out)   :: ok         ! the balance is solved

      type(surface_state_t) :: y_stage
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt, scaled, scaled_new
      integer(ik) :: ls, nfeval, fresh
      logical     :: reused
      real(wp),    parameter :: DETEPS = 1.0e-30_wp

      H0 = y%cas_enthalpy
      q0 = y%cas_shv
      Hk = H0
      qk = q0
      steps = 0_ik
      fresh = 0_ik
      ok = .false.
      reused = .false.
      y_stage%cas_co2 = y%cas_co2
      y_stage%cas_enthalpy = Hk
      y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%veg_surface, frozen%ground,                            &
         t_ground, n, surf_tend)
      nfeval = 1_ik
      R_H = cas_mass_capacity*(Hk - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%box%atm_enthalpy - Hk)
      R_q = cas_mass_capacity*(qk - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%box%atm_specific_humidity  - qk)
      scaled = surface_residual_scaled(R_H, R_q, Hk, qk, dt, cas_mass_capacity)

      do
         if (scaled <= 1.0_wp) then
            ok = .true.
            exit
         end if
         if (nfeval >= SURFACE_FEVAL_MAX) exit
         if (.not. slopes%valid) then
            !--- Fresh slopes at most SURFACE_NEWTON_MAX times, and once more after a step on held ones.
            if (fresh > SURFACE_NEWTON_MAX .or. (fresh == SURFACE_NEWTON_MAX .and. .not. reused)) exit
            call surface_source_slopes(Hk, qk, y%cas_co2, cas, frozen, t_ground, surf_tend, n, slopes)
            nfeval = nfeval + 2_ik
            fresh = fresh + 1_ik
            reused = .false.
         else
            reused = .true.
         end if
         steps = steps + 1_ik
         call surface_jacobian(slopes, cas_mass_capacity/dt, g_atm_heat, g_atm_vapour, J11, J12, J21, J22)
         detJ = J11*J22 - J12*J21
         if (detJ <= DETEPS*abs(J11*J22) .or. detJ <= 0.0_wp) then       ! singular / sign-flipped guard
            delH = -R_H / max(J11, tiny_num)                             ! damped-diagonal (Picard-like) fallback
            delq = -R_q / max(J22, tiny_num)
         else
            delH = (-R_H*J22 + R_q*J12)/detJ
            delq = (-R_q*J11 + R_H*J21)/detJ
         end if
         lam = 1.0_wp
         rn0 = R_H*R_H + R_q*R_q
         Ht = Hk
         qt = qk
         RHt = R_H
         Rqt = R_q
         do ls = 1_ik, SURFACE_LS_MAX
            Ht = Hk + lam*delH
            qt = qk + lam*delq
            !---~---
            !   No supersaturation clamp on the state: a jump there inflates the embedded error near
            ! RH = 1. Transient supersaturation is relaxed by surface_derivs' condensation sink.
            !---~---
            y_stage%cas_enthalpy = Ht
            y_stage%cas_shv = qt
            call surface_derivs(y_stage, cas, frozen%tissue, frozen%veg_surface, frozen%ground,                      &
               t_ground, n, surf_tend)
            nfeval = nfeval + 1_ik
            RHt = cas_mass_capacity*(Ht - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%box%atm_enthalpy - Ht)
            Rqt = cas_mass_capacity*(qt - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%box%atm_specific_humidity  - qt)
            if (sufficient_decrease(RHt*RHt + Rqt*Rqt, rn0, lam)) exit
            lam = 0.5_wp*lam
         end do
         Hk = Ht
         qk = qt
         R_H = RHt
         R_q = Rqt
         scaled_new = surface_residual_scaled(R_H, R_q, Hk, qk, dt, cas_mass_capacity)
         !---~---
         !   A step on held slopes that converged poorly -- shortened by the line search, or a fall
         ! of the residual below SLOPES_MIN_FALL while still unsolved -- has them taken afresh.
         !---~---
         if (reused .and. scaled_new > 1.0_wp .and. (ls > 1_ik .or. SLOPES_MIN_FALL*scaled_new > scaled))          &
            slopes%valid = .false.
         scaled = scaled_new
      end do
      cas_enthalpy = Hk
      cas_shv = qk
      !---~---
      !   surf_tend holds the last evaluation, which is at (Hk, qk): the start when the first test
      ! passes, else the accepted line-search trial. The caller commits the box from it.
      !---~---
   end subroutine newton_surface_solve

   !---~---
   !   The canopy-air balance's residuals as changes of the canopy air over the stage, each over its
   ! tolerance; the larger of the two. At or below one the balance is solved.
   !---~---
   pure real(wp) function surface_residual_scaled(r_h, r_q, h, q, dt, cas_mass_capacity) result(scaled)
      real(wp), intent(in) :: r_h, r_q                 ! [W/m2], [kg/m2/s] the residuals
      real(wp), intent(in) :: h, q                     ! [J/kg], [kg/kg] the canopy air's state
      real(wp), intent(in) :: dt, cas_mass_capacity    ! [s], [kg/m2]
      scaled = max(abs(r_h)*dt/cas_mass_capacity / (SURFACE_ATOL_ENTH + SURFACE_RTOL*abs(h)),                   &
                   abs(r_q)*dt/cas_mass_capacity / (SURFACE_ATOL_SHV + SURFACE_RTOL*abs(q)))
   end function surface_residual_scaled

   !---~---
   !   The backward-Euler canopy-air box commit, through the shared box kernel
   ! cas_box_step_implicit (ed_canopy_air_space; docs/science/canopy_air_space_biophysics.md
   ! eq. (2)): one implementation of the box, not a copy of its formulas here.
   !---~---
   pure subroutine cas_box_commit(h0, q0, c0, surf_tend, cas, dt, h1, q1, c1)
      real(wp),               intent(in)  :: h0, q0, c0, dt
      type(surface_tend_t),   intent(in)  :: surf_tend
      type(cas_boundary_t),   intent(in)  :: cas        ! the stage's boundary: its box carries the stage's conductances
      real(wp),               intent(out) :: h1, q1, c1
      type(cas_env_t) :: env
      env = cas%box
      env%surface_enthalpy_source = surf_tend%src_enth
      env%surface_vapor_source    = surf_tend%src_vap
      env%biotic_co2_source       = cas%nee_biotic
      call cas_box_step_implicit(h0, q0, c0, env, dt, h1, q1, c1)
   end subroutine cas_box_commit

   !---~---
   !   The surface's source slopes at (Hk, qk), by forward differences of surface_derivs from its
   ! evaluation there (surf_tend): how src_enth and src_vap respond to the canopy air's enthalpy
   ! and humidity. Two evaluations of the surface; the set is marked held.
   !---~---
   subroutine surface_source_slopes(Hk, qk, co2, cas, frozen, t_ground, surf_tend, n, slopes)
      real(wp),               intent(in)    :: Hk, qk, co2
      type(cas_boundary_t),   intent(in)    :: cas
      type(frozen_t),         intent(in)    :: frozen
      real(wp),               intent(in)    :: t_ground
      type(surface_tend_t),   intent(in)    :: surf_tend         ! the evaluation at (Hk, qk)
      integer(ik),            intent(in)    :: n
      type(surface_source_slopes_t), intent(inout) :: slopes
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq
      !--- The difference steps: the square root of the machine epsilon of each variable's scale.
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      y_stage%cas_co2 = co2
      y_stage%cas_enthalpy = Hk + dH
      y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%veg_surface, frozen%ground,                            &
         t_ground, n, sfp)
      slopes%dsrc_enth_dh = (sfp%src_enth - surf_tend%src_enth)/dH
      slopes%dsrc_vap_dh  = (sfp%src_vap - surf_tend%src_vap)/dH
      y_stage%cas_enthalpy = Hk
      y_stage%cas_shv = qk + dq
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%veg_surface, frozen%ground,                            &
         t_ground, n, sfp)
      slopes%dsrc_enth_dq = (sfp%src_enth - surf_tend%src_enth)/dq
      slopes%dsrc_vap_dq  = (sfp%src_vap - surf_tend%src_vap)/dq
      slopes%valid = .true.
   end subroutine surface_source_slopes

   !---~---
   !   The 2x2 Jacobian of (R_H, R_q) in (H, q): the box's own terms, its capacity over the stage
   ! (capacity_rate, cas_mass_capacity/dt) and its conductances to the atmosphere, both exact at
   ! every stage, less the surface's source slopes.
   !---~---
   pure subroutine surface_jacobian(slopes, capacity_rate, g_atm_heat, g_atm_vapour, J11, J12, J21, J22)
      type(surface_source_slopes_t), intent(in) :: slopes
      real(wp),               intent(in)  :: capacity_rate, g_atm_heat, g_atm_vapour   ! [kg/m2/s]
      real(wp),               intent(out) :: J11, J12, J21, J22
      J11 = capacity_rate + g_atm_heat - slopes%dsrc_enth_dh
      J12 =                            - slopes%dsrc_enth_dq
      J21 =                            - slopes%dsrc_vap_dh
      J22 = capacity_rate + g_atm_vapour - slopes%dsrc_vap_dq
   end subroutine surface_jacobian

   !---~---
   !   advance_plant_water -- each cohort's leaf and wood water over an ESDIRK2 step at the network's
   ! time-mean rates: sapflow and uptake in; sapflow and the b-weighted transpiration the stages
   ! charged the canopy air out. The stages hold transpiration at the network's, so a dt_fast step
   ! ends where the network did. Nothing is floored: the network keeps every store inside its curve.
   !---~---
   pure subroutine advance_plant_water(y, frozen, n, dt, transp_c_bw, y_out)
      type(march_state_t),  intent(in)    :: y
      type(frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: dt
      real(wp),              intent(in)    :: transp_c_bw(n)  ! [kg/m2 ground/s] per cohort, b-weighted over the step
      type(march_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      do i = 1_ik, n
         y_out%leaf_water_mass(i) = y%leaf_water_mass(i)                                                  &
              + dt * (frozen%plant%sapflow_frozen(i) - transp_c_bw(i) / max(frozen%plant%nplant(i), tiny_num))
         y_out%wood_water_mass(i) = y%wood_water_mass(i)                                                  &
              + dt * (frozen%plant%uptake_frozen(i) - frozen%plant%sapflow_frozen(i))
      end do
   end subroutine advance_plant_water

   !---~---
   !   advance_surf_water_full -- the canopy's surface water over the full dt, the surface-water
   ! analogue of advance_plant_water: frozen interception in, b-weighted surface-water evaporation
   ! out (the number the CAS was credited). Unclamped: update_dynamic_state clamps it and books the
   ! overflow and deficit.
   !---~---
   pure subroutine advance_surf_water_full(y, frozen, n, dt, leaf_surface_evap_bw, wood_surface_evap_bw, y_out)
      type(march_state_t),  intent(in)    :: y
      type(frozen_t), intent(in)    :: frozen
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt
      real(wp),               intent(in)    :: leaf_surface_evap_bw(n), wood_surface_evap_bw(n)
      type(march_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      !---~---
      !   Not floored: the caller's deficit bookkeeping must see negative surface water to credit the
      ! shortfall to the whole-column outflow; a floor here would create that water.
      !---~---
      do i = 1_ik, n
         y_out%leaf_surf_water(i) = y%leaf_surf_water(i) + dt*(frozen%veg_surface%intercept_leaf(i) - leaf_surface_evap_bw(i))
         y_out%wood_surf_water(i) = y%wood_surf_water(i) + dt*(frozen%veg_surface%intercept_wood(i) - wood_surface_evap_bw(i))
      end do
   end subroutine advance_surf_water_full

   !---~---
   !   The soil-plant water network's step (docs/science/soil_plant_water.md §3): implicit
   ! Euler on adaptive sub-steps, each a Newton solve (eq. (2)) in every store's own
   ! unknown, with a block-arrowhead linear solve (ed_linear_algebra) and a filtered
   ! local-error estimate. Water moves by the flows themselves, so it is conserved to
   ! round-off whatever Newton's residual.
   !---~---

   !---~---
   !   water_step -- the network over one fast step, from the state set_water_state set, under
   ! the GRP_THETA tolerances. On return the state holds the end-of-step water, unknowns and
   ! potentials, the flows what each path moved, and work_count the network's work: accepted
   ! sub-steps, Newton iterations, rejected sub-steps, the sub-steps taken although they did
   ! not converge or meet the tolerance (counted, not hidden) with the reason and the worst
   ! store's class, the node balances left unsolved, and the core nodes by class. The drive is
   ! held over the step, so the ground evaporation's constant terms are set into it once, here.
   !---~---
   subroutine water_step(dt, static, profile, soil, hyd, drive, tols, water, work_count)
      real(wp),                  intent(in)    :: dt
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(inout) :: drive
      type(tol_set_t),           intent(in)    :: tols
      type(water_work_t),        intent(inout) :: water
      integer(ik),               intent(inout) :: work_count(:)   ! the step report's work counters
      real(wp)    :: t, h, h_min, err, rtol, atol
      real(wp)    :: scale_soil(water%nsoil), scale_coh(COHORT_STORES, water%ncoh)
      real(wp)    :: tau_soil(water%nsoil), tau_coh(COHORT_STORES, water%ncoh)
      integer(ik) :: trials, iters, i_sol, worst_at, n_unsolved
      logical     :: converged, last, spilled
      logical     :: fast(water%nsoil), cut(water%nsoil)   ! the sub-step's eliminated and detached core nodes

      rtol = tols%rtol(GRP_THETA)
      atol = tols%atol(GRP_THETA)
      call clear_flows(water%flows)
      work_count(WORK_WATER_NSUB:WORK_WATER_NFAIL) = 0_ik
      work_count(WORK_WATER_FORCED_NEWTON:WORK_WATER_CUT_LAYER) = 0_ik
      call arrowhead_allocate(water%sys, water%nsoil, water%ncoh, COHORT_STORES)
      call water_scales(static, profile, hyd, drive%org, scale_soil, scale_coh)
      call ground_evaporation_setup(profile, soil, drive%ground)
      water%start_ready = .false.

      t = 0.0_wp
      h = dt
      h_min = H_MIN_SHARE * dt
      trials = 0_ik
      do while (t < dt * (1.0_wp - 1.0e-12_wp))
         h = min(h, dt - t)
         last = (h <= h_min) .or. (trials >= TRIALS_MAX)
         trials = trials + 1_ik
         tau_soil = water%state%tau_soil
         tau_coh = water%state%tau_coh
         call water_newton_substep(static, profile, soil, hyd, drive, h, atol, scale_soil, scale_coh,           &
                                   tau_soil, tau_coh, water, converged, iters, i_sol, fast, cut, worst_at, n_unsolved)
         work_count(WORK_WATER_NEWTON) = work_count(WORK_WATER_NEWTON) + iters
         work_count(WORK_WATER_NODE_UNSOLVED) = work_count(WORK_WATER_NODE_UNSOLVED) + n_unsolved
         if (.not. converged .and. .not. last) then
            work_count(WORK_WATER_REJECT) = work_count(WORK_WATER_REJECT) + 1_ik
            h = max(h * NEWTON_FAIL_CUT, h_min)
            cycle
         end if
         !---~---
         !   A sub-step whose start already balanced took no iteration: it ends where it began, its
         ! inflows unchanged, so its error is zero, and it has no factorisation of its own to
         ! filter with.
         !---~---
         if (.not. converged) then
            err = huge(1.0_wp)
         else if (iters == 0_ik) then
            err = 0.0_wp
         else
            call water_substep_error(water%sys_start, water%sys, h, water%ev(water%i_start), water%ev(i_sol),   &
                                     scale_soil, scale_coh, rtol, atol, fast .or. cut, err, worst_at)
         end if
         if (err > 1.0_wp .and. .not. last) then
            work_count(WORK_WATER_REJECT) = work_count(WORK_WATER_REJECT) + 1_ik
            h = max(h * max(SHRINK_MIN, SAFETY * err ** (-0.5_wp)), h_min)
            cycle
         end if
         if (err > 1.0_wp) call count_forced(converged, worst_at)
         call count_classes(fast, WORK_WATER_FAST_POND, WORK_WATER_FAST_SLOT, WORK_WATER_FAST_LAYER)
         call count_classes(cut, WORK_WATER_CUT_POND, WORK_WATER_CUT_SLOT, WORK_WATER_CUT_LAYER)
         !--- A sub-step's balance that is not finite is never committed: the run stops.
         if (.not. finite_evaluation(water%ev(i_sol), tau_soil, tau_coh))                                     &
            error stop 'water_step: a sub-step of the soil-plant water network is not finite'
         call water_commit_substep(drive, h, tau_soil, tau_coh, water%ev(i_sol), water%state, water%flows,  &
                                   fast, cut, spilled)
         !---~---
         !   The next sub-step starts where this one ended: its start is this one's solution, whose
         ! Jacobian sys still holds -- unless the spill moved the pond, or a node was detached, whose
         ! cut faces the next sub-step's classes must not read.
         !---~---
         water%i_start = i_sol
         water%start_ready = .not. (spilled .or. any(cut))
         t = t + h
         work_count(WORK_WATER_NSUB) = work_count(WORK_WATER_NSUB) + 1_ik
         if (err <= 1.0_wp)                                                                                   &
            h = h * min(GROW_MAX, max(SHRINK_MIN, SAFETY * max(err, 1.0e-10_wp) ** (-0.5_wp))) ! clamp-ok: controller floor
      end do
      water%state%psi_soil = water%ev(water%i_start)%psi_soil
      water%state%psi_coh  = water%ev(water%i_start)%psi_coh
   contains
      !--- An accepted sub-step's nodes of a class, by kind: the pond, the organic slots, the soil layers.
      subroutine count_classes(mask, i_pond, i_slot, i_layer)
         logical,     intent(in) :: mask(:)
         integer(ik), intent(in) :: i_pond, i_slot, i_layer   ! their counters in work_count
         if (mask(POND)) work_count(i_pond) = work_count(i_pond) + 1_ik
         work_count(i_slot)  = work_count(i_slot)  + count(mask(2:1+drive%org%n), kind=ik)
         work_count(i_layer) = work_count(i_layer) + count(mask(2+drive%org%n:), kind=ik)
      end subroutine count_classes

      !---~---
      !   A sub-step taken although unconverged: why (its Newton had not converged, or it failed the
      ! error filter), and the class of its worst store (worst_at: a core node, or nsoil + the
      ! cohort stores' column-major index; 0 if unknown).
      !---~---
      subroutine count_forced(newton_ok, at)
         logical,     intent(in) :: newton_ok
         integer(ik), intent(in) :: at
         work_count(WORK_WATER_NFAIL) = work_count(WORK_WATER_NFAIL) + 1_ik
         if (.not. newton_ok) work_count(WORK_WATER_FORCED_NEWTON) = work_count(WORK_WATER_FORCED_NEWTON) + 1_ik
         if (at > water%nsoil) then
            work_count(WORK_WATER_FORCED_COHORT) = work_count(WORK_WATER_FORCED_COHORT) + 1_ik
         else if (at >= 1_ik) then
            if (fast(at)) work_count(WORK_WATER_FORCED_FAST) = work_count(WORK_WATER_FORCED_FAST) + 1_ik
            if (cut(at))  work_count(WORK_WATER_FORCED_CUT)  = work_count(WORK_WATER_FORCED_CUT) + 1_ik
         end if
      end subroutine count_forced
   end subroutine water_step

   !--- Whether a sub-step's evaluation and unknowns are all finite numbers.
   pure logical function finite_evaluation(ev, tau_soil, tau_coh) result(ok)
      type(water_tend_t), intent(in) :: ev
      real(wp),           intent(in) :: tau_soil(:), tau_coh(:,:)
      ok = all(ieee_is_finite(tau_soil)) .and. all(ieee_is_finite(tau_coh))                                   &
           .and. all(ieee_is_finite(ev%water_soil)) .and. all(ieee_is_finite(ev%inflow_soil))                 &
           .and. all(ieee_is_finite(ev%water_coh)) .and. all(ieee_is_finite(ev%inflow_coh))
   end function finite_evaluation

   !---~---
   !   One sub-step: Newton on the end-of-sub-step water balance, with the core nodes in their
   ! classes (docs/science/soil_plant_water.md §3, "Node classes"). At the sub-step's start, from
   ! its evaluation, a node whose time scale is far below the sub-step is fast (classify_fast):
   ! every evaluation first solves its own balance with its neighbours held (nonlinear
   ! elimination), and Newton keeps its row with a zero residual, so the factorised matrix gives
   ! the rest the Schur-eliminated step. At every iteration, a connected group of nodes the solve
   ! cannot see is detached (inert_groups): its faces carry nothing and each member solves its
   ! own balance with its own sources. A step that would carry a node across a kink of its curve
   ! stops there for the iteration (chop_at_kinks). The node balances solve_node leaves unsolved
   ! are counted (n_unsolved).
   !---~---
   subroutine water_newton_substep(static, profile, soil, hyd, drive, h, atol, scale_soil, scale_coh,             &
                                   tau_soil, tau_coh, water, converged, iters, i_sol, fast, cut, worst_at, n_unsolved)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(in)    :: drive
      real(wp),                  intent(in)    :: h, atol, scale_soil(:), scale_coh(:,:)
      real(wp),                  intent(inout) :: tau_soil(:), tau_coh(:,:)   ! start of sub-step in, solution out
      type(water_work_t),        intent(inout) :: water
      logical,                   intent(out)   :: converged
      integer(ik),               intent(out)   :: iters
      integer(ik),               intent(out)   :: i_sol   ! which of water%ev holds the solution
      logical,                   intent(out)   :: fast(:), cut(:)   ! (nsoil) the eliminated and detached nodes
      integer(ik),               intent(out)   :: worst_at     ! the store of the last worst imbalance (water_residual)
      integer(ik),               intent(out)   :: n_unsolved   ! node balances solve_node left unsolved
      real(wp) :: res_soil(water%nsoil), res_coh(COHORT_STORES, water%ncoh)
      real(wp) :: step_soil(water%nsoil), step_coh(COHORT_STORES, water%ncoh)
      real(wp) :: try_soil(water%nsoil), try_coh(COHORT_STORES, water%ncoh)
      real(wp) :: res2_soil(water%nsoil), res2_coh(COHORT_STORES, water%ncoh)
      real(wp) :: merit, merit_try, worst, worst_try, lam, lam_max, tol
      integer(ik) :: it, ib, i_try, worst_at_try
      real(wp) :: step_free(water%nsoil)   ! the step before chop_at_kinks
      logical  :: ok, cut_new(water%nsoil)
      logical  :: chopped                  ! chop_at_kinks shortened the step

      converged = .false.
      iters = 0_ik
      worst_at = 0_ik
      n_unsolved = 0_ik
      fast = .false.
      cut  = .false.
      tol  = NEWTON_TOL * atol / h
      !---~---
      !   The network at the start, evaluated unless it already is; the iterates are evaluated into
      ! the other slot, and the solution is the start until one is taken.
      !---~---
      i_sol = water%i_start
      i_try = 3_ik - water%i_start
      if (.not. water%start_ready) then
         call water_derivs(static, profile, soil, hyd, drive, tau_soil, tau_coh, water%ev(i_sol), water%sys, .true.)
         water%start_ready = .true.
      end if
      call classify_fast(water%sys, water%ev(i_sol), h, tol, scale_soil, fast)
      call inert_groups(water%ev(i_sol), h, tol, scale_soil, cut)
      fast = fast .and. .not. cut
      if (any(fast) .or. any(cut)) then
         water%start_ready = .false.
         call evaluate(tau_soil, tau_coh, i_sol)
      end if
      call water_residual(water%state, h, water%ev(i_sol), atol, scale_soil, scale_coh, res_soil, res_coh,      &
                          merit, worst, worst_at, cut)
      if (merit /= merit) return
      do it = 1_ik, NEWTON_MAX
         if (worst <= 1.0_wp) exit
         iters = it
         !--- A group the solve has lost sight of is detached; one whose links have returned rejoins.
         call inert_groups(water%ev(i_sol), h, tol, scale_soil, cut_new)
         if (any(cut_new .neqv. cut)) then
            cut  = cut_new
            fast = fast .and. .not. cut
            call evaluate(tau_soil, tau_coh, i_sol)
            call water_residual(water%state, h, water%ev(i_sol), atol, scale_soil, scale_coh, res_soil, res_coh,  &
                                merit, worst, worst_at, cut)
            if (merit /= merit) return
            if (worst <= 1.0_wp) exit
         end if
         call detach_rows(water%sys, cut, h, scale_soil)
         call water_newton_factor(water%sys, h, water%ev(i_sol), ok)
         if (.not. ok) return
         if (it == 1_ik) water%sys_start = water%sys   ! the start's factorisation, for the error filter
         call arrowhead_solve(water%sys, merge(0.0_wp, -res_soil, fast .or. cut), -res_coh, step_soil, step_coh)
         !--- An eliminated or detached node's own solve, not the linear step, gives its unknown.
         step_soil = merge(0.0_wp, step_soil, fast .or. cut)
         step_free = step_soil
         call chop_at_kinks(profile, drive, tau_soil, step_soil)
         chopped = any(step_soil /= step_free)
         lam_max = step_limit(profile, drive%org, tau_soil, tau_coh, step_soil, step_coh)
         lam = lam_max
         !--- A trial's evaluation replaces the Jacobian sys held, so the start's is gone.
         water%start_ready = .false.
         do ib = 1_ik, BACKTRACK_MAX
            try_soil = tau_soil + lam * step_soil
            try_coh = tau_coh + lam * step_coh
            call evaluate(try_soil, try_coh, i_try)
            call water_residual(water%state, h, water%ev(i_try), atol, scale_soil, scale_coh, res2_soil,          &
                                res2_coh, merit_try, worst_try, worst_at_try, cut)
            if (merit_try == merit_try .and. sufficient_decrease(merit_try, merit, lam)) exit
            !---~---
            !   A step chopped at a kink is not Newton's, and the decrease Armijo's test expects
            ! of Newton's it cannot give: landing a node on a kink a hair away moves the merit by
            ! a hair. It is taken if the merit does not grow, so the next step may cross.
            !---~---
            if (chopped .and. merit_try <= merit) exit
            lam = 0.5_wp * lam
         end do
         if (merit_try /= merit_try) return
         tau_soil = try_soil
         tau_coh = try_coh
         res_soil = res2_soil
         res_coh = res2_coh
         merit = merit_try
         worst = worst_try
         worst_at = worst_at_try
         i_sol = i_try
      end do
      !---~---
      !   The error estimate solves with the first iteration's factorisation, at the start, and
      ! the last one's, one Newton step behind the solution, not a fresh one at it: its filter
      ! needs the stiff exchanges' scale, not the exact matrix, and a factorisation per sub-step
      ! is saved.
      !---~---
      converged = worst <= 1.0_wp

   contains

      !---~---
      !   The network at an iterate: each eliminated or detached node's unknown from its own balance
      ! (solve_node), the others' held, then every store's flows and the Jacobian, the detached
      ! nodes' faces cut.
      !---~---
      subroutine evaluate(ts, tc, i_ev)
         real(wp),    intent(inout) :: ts(:), tc(:,:)
         integer(ik), intent(in)    :: i_ev
         real(wp)    :: psi_wood(water%ncoh)
         integer(ik) :: c
         logical     :: solved
         if (any(fast) .or. any(cut)) then
            call wood_potentials(static, hyd, tc, psi_wood)
            do c = 1_ik, water%nsoil
               if (.not. (fast(c) .or. cut(c))) cycle
               call solve_node(static, profile, soil, drive, c, ts, water%state%water_soil(c), h,                &
                               NODE_SOLVE_SHARE * NEWTON_TOL * atol * scale_soil(c) / h, cut, psi_wood, solved)
               if (.not. solved) n_unsolved = n_unsolved + 1_ik
            end do
         end if
         call water_derivs(static, profile, soil, hyd, drive, ts, tc, water%ev(i_ev), water%sys, .true., cut)
      end subroutine evaluate

   end subroutine water_newton_substep

   !---~---
   !   The fast core nodes at a sub-step's start (docs/science/soil_plant_water.md §3, "Node
   ! classes"): a node's storage S = |dW/dtau| over the flow part of its Jacobian column, L, the sum
   ! of how far its unknown moves every balance, is the time it takes to relax to its neighbours,
   ! whatever the unknown's scale. Below FAST_SHARE of the sub-step, with a link the solve can see
   ! (seen_links, as inert_groups reads them), its storage is under the Newton tolerance's share of
   ! its flows, and it is eliminated. Of two adjacent fast nodes only the faster is eliminated.
   !---~---
   pure subroutine classify_fast(sys, ev, h, tol, scale_soil, fast)
      type(arrowhead_system_t), intent(in)  :: sys
      type(water_tend_t),       intent(in)  :: ev
      real(wp),                 intent(in)  :: h, tol, scale_soil(:)
      logical,                  intent(out) :: fast(:)
      real(wp)    :: big_l, t_relax(sys%nc)
      integer(ik) :: c, i, ib
      logical     :: seen(sys%nc), linked
      call seen_links(ev, tol, scale_soil, seen)
      do c = 1_ik, sys%nc
         big_l = 0.0_wp
         do i = 1_ik, sys%nc
            big_l = big_l + abs(sys%a(i, c))
         end do
         do ib = 1_ik, sys%nb
            big_l = big_l + abs(sys%c(c, ib))
         end do
         linked = seen(c) .or. ev%root_link(c) >= tol * scale_soil(c)
         if (c > 1_ik) linked = linked .or. seen(c - 1_ik)
         t_relax(c) = abs(ev%dwater_soil(c)) / max(big_l, tiny_num)
         fast(c) = linked .and. t_relax(c) < FAST_SHARE * h
      end do
      do c = 1_ik, sys%nc - 1_ik
         if (fast(c) .and. fast(c + 1_ik)) then
            if (t_relax(c) <= t_relax(c + 1_ik)) then
               fast(c + 1_ik) = .false.
            else
               fast(c) = .false.
            end if
         end if
      end do
   end subroutine classify_fast

   !---~---
   !   The chain's faces the solve can see: face c, joining nodes c and c + 1, is seen when a whole
   ! unit of either end's unknown moves either end's balance by the convergence tolerance (each
   ! weighed as water_residual weighs that balance), from its slopes whether or not it is cut.
   !---~---
   pure subroutine seen_links(ev, tol, scale_soil, seen)
      type(water_tend_t), intent(in)  :: ev
      real(wp),           intent(in)  :: tol, scale_soil(:)
      logical,            intent(out) :: seen(:)
      integer(ik) :: c, nc
      nc = size(seen, kind=ik)
      do c = 1_ik, nc - 1_ik
         seen(c) = max(abs(ev%face_dqa(c)), abs(ev%face_dqb(c))) >= tol * min(scale_soil(c), scale_soil(c + 1_ik))
      end do
      seen(nc) = .false.
   end subroutine seen_links

   !---~---
   !   The inert core nodes (docs/science/soil_plant_water.md §3, "Node classes"): the graph of the
   ! links the solve can see, each weighed as water_residual weighs a balance (a link is seen
   ! when a whole unit of one end's unknown moves the other end's balance by the convergence
   ! tolerance), read from every face's and root's slopes whether or not they are cut. A connected
   ! group none of whose members holds storage the solve can see, with no link it can see to a
   ! node outside it or to a cohort, is detached as one unit. Its matrix would be singular, or so
   ! close to it that a step sends its unknowns off by kilometres: an empty pond on sealed frozen
   ! ground, a frozen layer whose pores are full, and those with frozen organic slots between.
   !---~---
   pure subroutine inert_groups(ev, h, tol, scale_soil, cut)
      type(water_tend_t), intent(in)  :: ev
      real(wp),           intent(in)  :: h, tol, scale_soil(:)
      logical,            intent(out) :: cut(:)
      integer(ik) :: nc, c, first
      logical     :: seen_link(size(cut)), active
      nc = size(cut, kind=ik)
      call seen_links(ev, tol, scale_soil, seen_link)
      cut = .false.
      !--- The chain's groups: runs of nodes joined by seen faces.
      first = 1_ik
      active = .false.
      do c = 1_ik, nc
         if (abs(ev%dwater_soil(c)) / h >= tol * scale_soil(c) .or. ev%root_link(c) >= tol * scale_soil(c))      &
            active = .true.
         if (.not. seen_link(c)) then
            if (.not. active) cut(first:c) = .true.
            first  = c + 1_ik
            active = .false.
         end if
      end do
   end subroutine inert_groups

   !--- A detached node's row and column cleared and its diagonal set, so the rest solve without it.
   pure subroutine detach_rows(sys, cut, h, scale_soil)
      type(arrowhead_system_t), intent(inout) :: sys
      logical,                  intent(in)    :: cut(:)
      real(wp),                 intent(in)    :: h, scale_soil(:)
      integer(ik) :: c
      do c = 1_ik, sys%nc
         if (.not. cut(c)) cycle
         sys%a(c, :) = 0.0_wp
         sys%a(:, c) = 0.0_wp
         sys%b(c, :) = 0.0_wp
         sys%c(c, :) = 0.0_wp
         sys%a(c, c) = scale_soil(c) / h
      end do
   end subroutine detach_rows

   !---~---
   !   solve_node -- core node c's unknown from its own balance (node_balance), every other store
   ! held: Newton inside a bracket, bisecting when a step leaves the bracket or falls short. The
   ! bracket is found by stepping out from the current unknown, so of several roots the nearest
   ! is taken; on the dry side of the curve the steps are geometric, the unknown being a
   ! saturation that may span many decades. Solved to tol [kg/m2/s], far inside Newton's; solved
   ! says whether it was.
   !---~---
   subroutine solve_node(static, profile, soil, drive, c, ts, w0, h, tol, cut, psi_wood, solved)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(water_drive_t),       intent(in)    :: drive
      integer(ik),               intent(in)    :: c
      real(wp),                  intent(inout) :: ts(:)
      real(wp),                  intent(in)    :: w0, h, tol
      logical,                   intent(in)    :: cut(:)        ! (nsoil) the detached nodes
      real(wp),                  intent(in)    :: psi_wood(:)
      logical,                   intent(out)   :: solved
      real(wp)    :: x, r, dr, a, b, ra, rb, dx, x_new
      integer(ik) :: it
      logical     :: have_a, have_b
      solved = .true.
      x = ts(c)
      call balance_at(x, r, dr)
      if (abs(r) <= tol) return
      have_a = r < 0.0_wp
      have_b = .not. have_a
      a = x
      b = x
      ra = r
      rb = r
      !--- Step out towards the root until the residual changes sign.
      dx = -r / merge(dr, 1.0_wp, dr > 0.0_wp)
      do it = 1_ik, NODE_SOLVE_MAX
         if (have_a .and. have_b) exit
         if (have_a) then
            x_new = x + max(abs(dx), 1.0e-3_wp * max(abs(x), 1.0_wp))
         else
            x_new = x - max(abs(dx), 1.0e-3_wp * max(abs(x), 1.0_wp))
            if (c /= POND) x_new = max(x_new, max(SOIL_TAU_FLOOR, 1.0e-3_wp * x))
         end if
         x = x_new
         call balance_at(x, r, dr)
         if (abs(r) <= tol) then
            ts(c) = x
            return
         end if
         if (r < 0.0_wp) then
            have_a = .true.
            a = x
            ra = r
         else
            have_b = .true.
            b = x
            rb = r
         end if
         dx = 4.0_wp * dx
         if (c /= POND .and. .not. have_a .and. x <= SOIL_TAU_FLOOR) exit
      end do
      if (.not. (have_a .and. have_b)) then
         ts(c) = x
         solved = .false.
         return
      end if
      !--- Newton inside the bracket.
      x = merge(a, b, abs(ra) < abs(rb))
      call balance_at(x, r, dr)
      do it = 1_ik, NODE_SOLVE_MAX
         if (abs(r) <= tol) exit
         x_new = x - r / merge(dr, tiny_num, dr > 0.0_wp)
         if (.not. (dr > 0.0_wp) .or. x_new <= min(a, b) .or. x_new >= max(a, b)) then
            if (c /= POND .and. min(a, b) > 0.0_wp .and. max(a, b) > 1.0e2_wp * min(a, b)) then
               x_new = sqrt(a * b)
            else
               x_new = 0.5_wp * (a + b)
            end if
         end if
         if (abs(x_new - x) <= 1.0e-15_wp * max(abs(x), 1.0e-12_wp)) exit   ! clamp-ok: the scale of a unit-free unknown near zero
         x = x_new
         call balance_at(x, r, dr)
         if (r < 0.0_wp) then
            a = x
         else
            b = x
         end if
      end do
      ts(c) = x
      solved = abs(r) <= tol
   contains
      subroutine balance_at(xx, rr, drr)
         real(wp), intent(in)  :: xx
         real(wp), intent(out) :: rr, drr
         ts(c) = xx
         call node_balance(static, profile, soil, drive, c, ts, w0, h, cut, psi_wood, rr, drr)
      end subroutine balance_at
   end subroutine solve_node

   !---~---
   !   A Newton step that would carry a node across a kink of its curve's storage (curve_kinks:
   ! Campbell's air entry) lands on the kink instead, that node's step alone (Appleyard's
   ! chopping; Wang & Tchelepi 2013), which the line search takes if the merit does not grow;
   ! from the kink, the next step may cross.
   !---~---
   pure subroutine chop_at_kinks(profile, drive, tau_soil, step_soil)
      type(soil_profile_t), intent(in)    :: profile
      type(water_drive_t),  intent(in)    :: drive
      real(wp),             intent(in)    :: tau_soil(:)
      real(wp),             intent(inout) :: step_soil(:)
      integer(ik) :: c, n_org
      n_org = drive%org%n
      do c = 2_ik, size(tau_soil, kind=ik)
         if (step_soil(c) == 0.0_wp) cycle
         if (c <= 1_ik + n_org) then
            call chop(drive%org%curve(c - 1_ik), tau_soil(c), step_soil(c))
         else
            call chop(profile%layer_curve(c - 1_ik - n_org), tau_soil(c), step_soil(c))
         end if
      end do
   contains
      pure subroutine chop(cv, tau, d)
         type(soil_layer_curve_t), intent(in)    :: cv
         real(wp),                 intent(in)    :: tau
         real(wp),                 intent(inout) :: d
         real(wp)    :: kinks(1)
         integer(ik) :: n, i
         call curve_kinks(cv, kinks, n)
         do i = 1_ik, n
            call land(kinks(i), tau, d)
         end do
      end subroutine chop
      pure subroutine land(kink, tau, d)
         real(wp), intent(in)    :: kink, tau
         real(wp), intent(inout) :: d
         !--- On the kink already (to round-off, its scale floored near zero): the step may cross it.
         if (abs(tau - kink) <= 1.0e-12_wp * max(abs(kink), 1.0e-12_wp)) return   ! clamp-ok: scale near zero
         if ((tau - kink) * (tau + d - kink) < 0.0_wp) d = kink - tau
      end subroutine land
   end subroutine chop_at_kinks

   !--- One store's water-balance residual [kg/m2/s]: (W(tau) - W_start)/h - inflow.
   elemental real(wp) function store_residual(water, water_start, h, inflow) result(res)
      real(wp), intent(in) :: water, water_start, h, inflow
      res = (water - water_start) / h - inflow
   end function store_residual

   !---~---
   !   Every store's residual and, in the same pass, two measures of it, each store's imbalance
   ! over the sub-step against NEWTON_TOL of its absolute water tolerance: worst, below one,
   ! means the balance is solved; merit, the sum of squares, is what the line search asks to
   ! fall -- a Newton step always points down it, which is not true of the worst store alone.
   ! The soil's squares and the cohorts' are summed apart, then added. worst_at is the worst
   ! store: a core node, or nsoil + the cohort stores' column-major index. An imbalance that is not
   ! a finite number makes the merit NaN, which every caller refuses, before any bound applies.
   !---~---
   pure subroutine water_residual(state, h, ev, atol, scale_soil, scale_coh, res_soil, res_coh, merit, worst,    &
                                  worst_at, skip)
      type(water_state_t), intent(in)  :: state
      real(wp),            intent(in)  :: h, atol, scale_soil(:), scale_coh(:,:)
      type(water_tend_t),  intent(in)  :: ev
      real(wp),            intent(out) :: res_soil(:), res_coh(:,:), merit, worst
      integer(ik),         intent(out) :: worst_at
      logical, optional,   intent(in)  :: skip(:)   ! detached nodes: their water is their sources', whatever their curve
      real(wp)    :: tol, x, merit_coh
      integer(ik) :: c, i, j, nsoil
      logical     :: finite
      tol = NEWTON_TOL * atol
      merit = 0.0_wp
      worst = 0.0_wp
      worst_at = 0_ik
      finite = .true.
      nsoil = size(res_soil, kind=ik)
      do c = 1_ik, nsoil
         res_soil(c) = store_residual(ev%water_soil(c), state%water_soil(c), h, ev%inflow_soil(c))
         if (present(skip)) then
            if (skip(c)) res_soil(c) = 0.0_wp
         end if
         x = res_soil(c) * h / (tol * scale_soil(c))
         finite = finite .and. ieee_is_finite(x)
         x = sign(min(abs(x), MERIT_X_MAX), x)   ! clamp-ok: a wild trial's square stays finite, and is refused
         merit = merit + x * x
         if (abs(x) > worst) worst_at = c
         worst = max(worst, abs(x))
      end do
      merit_coh = 0.0_wp
      do i = 1_ik, size(res_coh, 2, kind=ik)
         do j = 1_ik, COHORT_STORES
            res_coh(j, i) = store_residual(ev%water_coh(j, i), state%water_coh(j, i), h, ev%inflow_coh(j, i))
            x = res_coh(j, i) * h / (tol * scale_coh(j, i))
            finite = finite .and. ieee_is_finite(x)
            x = sign(min(abs(x), MERIT_X_MAX), x)   ! clamp-ok: a wild trial's square stays finite, and is refused
            merit_coh = merit_coh + x * x
            if (abs(x) > worst) worst_at = nsoil + COHORT_STORES * (i - 1_ik) + j
            worst = max(worst, abs(x))
         end do
      end do
      merit = merit + merit_coh
      if (.not. finite) merit = ieee_value(merit, ieee_quiet_nan)
   end subroutine water_residual

   !---~---
   !   The largest share of a Newton step every store can take: no organic slot, soil layer or
   ! tissue is driven more than STEP_MAX_SHARE of the way to empty, and none on the dry side
   ! moves its saturation or relative water content by more than STEP_MAX_DRY. A wet slot's or
   ! layer's unknown is its potential, so only the first applies to it: it keeps a step that
   ! would carry it past the branch point and below empty inside the curve. The pond's head is
   ! not bounded.
   !---~---
   pure real(wp) function step_limit(profile, org, tau_soil, tau_coh, step_soil, step_coh) result(lam)
      type(soil_profile_t),  intent(in) :: profile
      type(water_organic_t), intent(in) :: org
      real(wp),              intent(in) :: tau_soil(:), tau_coh(:,:), step_soil(:), step_coh(:,:)
      integer(ik) :: k, i, j, c
      lam = 1.0_wp
      do j = 1_ik, org%n
         c = 1_ik + j
         lam = min(lam, fraction_to_boundary(tau_soil(c) - SOIL_TAU_FLOOR, step_soil(c), STEP_MAX_SHARE))
         if (tau_soil(c) <= org%curve(j)%se_switch) call cap_dry_change(tau_soil(c), step_soil(c), lam)
      end do
      do k = 1_ik, profile%n_active
         c = 1_ik + org%n + k
         lam = min(lam, fraction_to_boundary(tau_soil(c) - SOIL_TAU_FLOOR, step_soil(c), STEP_MAX_SHARE))
         if (tau_soil(c) <= profile%layer_curve(k)%se_switch) call cap_dry_change(tau_soil(c), step_soil(c), lam)
      end do
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         do j = 1_ik, COHORT_STORES
            lam = min(lam, fraction_to_boundary(tau_coh(j, i), step_coh(j, i), STEP_MAX_SHARE))
            call cap_dry_change(tau_coh(j, i), step_coh(j, i), lam)
         end do
      end do
   contains
      !---~---
      !   The cap on the change applies below full saturation (an RWC or Se of one), where the
      ! potential steepens towards empty; above it the curve is gentle and Newton moves freely.
      !---~---
      pure subroutine cap_dry_change(tau, d, lam)
         real(wp), intent(in)    :: tau, d
         real(wp), intent(inout) :: lam
         if (tau < 1.0_wp .and. abs(d) > STEP_MAX_DRY) lam = min(lam, STEP_MAX_DRY / abs(d))
      end subroutine cap_dry_change
   end function step_limit

   !--- The storage term dW/dtau / h each store adds to the Newton matrix's diagonal.
   pure subroutine water_storage_diagonal(h, ev, add_soil, add_coh)
      real(wp),           intent(in)  :: h
      type(water_tend_t), intent(in)  :: ev
      real(wp),           intent(out) :: add_soil(:), add_coh(:,:)
      add_soil = ev%dwater_soil / h
      add_coh  = ev%dwater_coh / h
   end subroutine water_storage_diagonal

   !---~---
   !   Factor the Newton matrix: the flows' Jacobian of the latest evaluation (sys) plus the
   ! storage terms of ev.
   !---~---
   pure subroutine water_newton_factor(sys, h, ev, ok)
      type(arrowhead_system_t), intent(inout) :: sys
      real(wp),                 intent(in)    :: h
      type(water_tend_t),       intent(in)    :: ev
      logical,                  intent(out)   :: ok
      real(wp) :: add_soil(size(ev%dwater_soil)), add_coh(COHORT_STORES, size(ev%dwater_coh, 2))
      call water_storage_diagonal(h, ev, add_soil, add_coh)
      call arrowhead_factor(sys, add_soil, add_coh, ok)
   end subroutine water_newton_factor


   !---~---
   !   The sub-step's local error relative to its tolerance (<= 1 to accept): implicit Euler's
   ! h (F_end - F_start)/2, filtered twice through (I - h dF/dW)^-1 so a stiff exchange that
   ! has already settled does not read as error (docs/science/soil_plant_water.md §3). Each
   ! store is weighed against atol of its scale plus rtol of its water. One pass filters through
   ! the sub-step's first factorisation, at its start, the other through its last, near its end,
   ! each with its own state's storage slopes: an exchange stiff at either end is filtered, as
   ! one that shuts as it settles (a refilling tissue's one-way uptake at the foot of its ramp)
   ! is stiff only at the start. The last factorisation keeps the coupling rows of the solution
   ! the last evaluation wrote: close to the matrix at the end, which is all a filter needs.
   !---~---
   pure subroutine water_substep_error(sys0, sys1, h, ev0, ev1, scale_soil, scale_coh, rtol, atol, skip, err,     &
                                       worst_at)
      type(arrowhead_system_t), intent(in) :: sys0, sys1 ! the Newton matrix, factored at the first and last iteration
      real(wp),                 intent(in) :: h
      type(water_tend_t),       intent(in) :: ev0, ev1   ! the network at the sub-step's start and end
      real(wp),                 intent(in) :: scale_soil(:), scale_coh(:,:), rtol, atol
      logical,                  intent(in) :: skip(:)    ! eliminated and detached nodes: their error is others'
      real(wp),                 intent(out) :: err       ! [-] the error over its tolerance
      integer(ik),              intent(out) :: worst_at  ! the store of the largest error, as water_residual's
      real(wp) :: v_soil(size(ev1%water_soil)), v_coh(COHORT_STORES, size(ev1%water_coh, 2))
      real(wp) :: e_soil(size(ev1%water_soil)), e_coh(COHORT_STORES, size(ev1%water_coh, 2))
      integer(ik) :: loc(2)
      real(wp) :: y_soil(size(ev1%water_soil)), y_coh(COHORT_STORES, size(ev1%water_coh, 2))
      v_soil = merge(0.0_wp, h * (ev1%inflow_soil - ev0%inflow_soil), skip)
      v_coh  = h * (ev1%inflow_coh - ev0%inflow_coh)
      call arrowhead_solve(sys0, v_soil / h, v_coh / h, y_soil, y_coh)
      v_soil = merge(0.0_wp, ev0%dwater_soil * y_soil, skip)
      v_coh  = ev0%dwater_coh * y_coh
      call arrowhead_solve(sys1, v_soil / h, v_coh / h, y_soil, y_coh)
      v_soil = merge(0.0_wp, ev1%dwater_soil * y_soil, skip)
      v_coh  = ev1%dwater_coh * y_coh
      e_soil = 0.5_wp * abs(v_soil) / (atol * scale_soil + rtol * abs(ev1%water_soil))
      err = maxval(e_soil)
      worst_at = maxloc(e_soil, dim=1, kind=ik)
      if (size(v_coh, 2) > 0) then
         e_coh = 0.5_wp * abs(v_coh) / (atol * scale_coh + rtol * abs(ev1%water_coh))
         if (maxval(e_coh) > err) then
            loc = maxloc(e_coh, kind=ik)
            worst_at = size(e_soil, kind=ik) + COHORT_STORES * (loc(2) - 1_ik) + loc(1)
         end if
         err = max(err, maxval(e_coh))
      end if
      if (err /= err) err = huge(1.0_wp)
   end subroutine water_substep_error

   !---~---
   !   Accept a sub-step: advance each store's water by the flows themselves, spill the pond above
   ! its capacity, carry the pond's enthalpy and the organic slots', and add what each path moved
   ! to the step's amounts. A store whose own balance was solved apart -- an eliminated core node,
   ! a cohort's leaves on their one path -- takes its water from its own curve, and the face to
   ! its strongest neighbour (the wood, for the leaves) carries the flow that closes its balance;
   ! the neighbour absorbs the slack (docs/science/soil_plant_water.md §3, "Conservation").
   !---~---
   subroutine water_commit_substep(drive, h, tau_soil, tau_coh, ev, state, flows, fast, cut, spilled)
      type(water_drive_t), intent(in)    :: drive
      real(wp),            intent(in)    :: h, tau_soil(:), tau_coh(:,:)
      type(water_tend_t),  intent(in)    :: ev
      type(water_state_t), intent(inout) :: state
      type(water_flows_t), intent(inout) :: flows
      logical,             intent(in)    :: fast(:), cut(:)   ! the eliminated and the detached core nodes
      logical,             intent(out)   :: spilled   ! the pond spilled, so its unknown moved
      real(wp) :: w_mix, e_mix, t_mix, u_mix, amount, spill, t_pond, pond_flux, sap(size(tau_coh, 2)), u_up, q
      real(wp) :: in_m, in_h, out_m, u_out, share   ! the organic layer's water in and out over the sub-step
      real(wp) :: qf(size(tau_soil)), dinfl(size(tau_soil)), slack
      integer(ik) :: i, nsl, n_org, nsoil, soff, c, g
      logical     :: snap(size(tau_soil))   ! the eliminated nodes taken to their curve's water

      nsl   = size(flows%face, kind=ik)
      n_org = drive%org%n
      nsoil = size(tau_soil, kind=ik)
      soff  = 1_ik + n_org
      !--- Down the chain: face c joins node c to c + 1.
      qf = 0.0_wp
      qf(1) = ev%pond_flux
      qf(2:1+n_org) = ev%org_face(1:n_org)
      qf(soff+1:soff+nsl-1) = ev%face_flux(1:nsl-1)
      !---~---
      !   Each eliminated node at its curve's water: the face to its strongest neighbour carries the
      ! flow that closes its balance there -- a neighbour that is not detached and holds the water
      ! the slack would take from it. If neither neighbour can, the node keeps the water its flows
      ! gave it, its slack Newton's, as any store's.
      !---~---
      dinfl = 0.0_wp
      snap  = fast
      if (.not. fast(POND)) call drain_only_what_is_held()
      do c = 1_ik, nsoil
         if (.not. fast(c)) cycle
         slack = (ev%water_soil(c) - state%water_soil(c)) / h - (ev%inflow_soil(c) + dinfl(c))
         if (c == 1_ik) then
            g = 2_ik
         else if (c == nsoil) then
            g = nsoil - 1_ik
         else if (abs(ev%face_dqb(c - 1_ik)) >= abs(ev%face_dqa(c))) then
            g = c - 1_ik
         else
            g = c + 1_ik
         end if
         if (.not. can_give(g)) then
            g = 2_ik * c - g
            if (g < 1_ik .or. g > nsoil) then
               snap(c) = .false.
               cycle
            end if
            if (.not. can_give(g)) then
               snap(c) = .false.
               cycle
            end if
         end if
         dinfl(c) = dinfl(c) + slack
         dinfl(g) = dinfl(g) - slack
         if (g == c - 1_ik) then
            qf(g) = qf(g) + slack
         else
            qf(c) = qf(c) - slack
         end if
      end do
      if (fast(POND) .and. .not. snap(POND)) call drain_only_what_is_held()
      pond_flux = qf(1)
      !--- The flow that closes each leaf's balance at its curve's water.
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         sap(i) = (ev%water_coh(LEAF, i) - state%water_coh(LEAF, i)) / h + ev%transp(i)
      end do

      !---~---
      !   The pond's enthalpy, in the order its water moves: the water reaching it joins it first,
      ! so what infiltrates leaves at the mixture's temperature; water the node below pushes up arrives
      ! at that node's temperature, the organic layer's or the top layer's; then the spill leaves
      ! at the pond's own. The organic slots gain what the pond gives them and lose what they push
      ! up, and the flow across the top layer's top face carries the giver's temperature.
      !---~---
      w_mix = state%water_soil(POND) + h * drive%rain
      e_mix = state%e_pond + h * drive%rain * internal_energy_liquid(drive%t_rain)
      flows%rain_enth = flows%rain_enth + h * drive%rain * internal_energy_liquid(drive%t_rain)
      t_mix = drive%t_rain
      if (w_mix > POND_TINY) t_mix = liquid_energy_to_temp(e_mix / w_mix)
      u_mix = internal_energy_liquid(t_mix)
      if (n_org == 0_ik) then
         u_up = internal_energy_liquid(drive%t_soil(1))
         if (pond_flux >= 0.0_wp) then
            amount = h * pond_flux
            flows%infil = flows%infil + amount
            flows%infil_enth = flows%infil_enth + amount * u_mix
            state%e_pond = e_mix - amount * u_mix
         else
            amount = -h * pond_flux
            flows%exfil = flows%exfil + amount
            flows%exfil_enth = flows%exfil_enth + amount * u_up
            state%e_pond = e_mix + amount * u_up
         end if
      else
         !---~---
         !   The organic layer, one thermal node: water passing through it in a sub-step leaves at
         ! the temperature it came in at, as far as what came in covers what leaves, and at the
         ! layer's own for the rest -- so the water a thin litter passes on does not heat or cool
         ! it by its difference from the layer's start, a store of a few hundred joules a kelvin.
         ! The layer's own water leaves as its energy holds it, its liquid share at its temperature
         ! and the rest ice, as its spill does: on the melt plateau the freezing curve lets nearly
         ! frozen litter drain, and water valued as liquid would take latent heat the litter does
         ! not have. What receives the ice share refreezes it, the soil by its freezing curve, the
         ! pond in the snow stage.
         !---~---
         q = h * qf(soff)
         in_m  = h * max(pond_flux, 0.0_wp) + max(-q, 0.0_wp)
         in_h  = h * max(pond_flux, 0.0_wp) * u_mix + max(-q, 0.0_wp) * internal_energy_liquid(drive%t_soil(1))
         out_m = h * max(-pond_flux, 0.0_wp) + max(q, 0.0_wp)
         u_out = temp_to_internal_energy(0.0_wp, 1.0_wp, drive%org%t_org, drive%org%fliq)
         if (out_m > 0.0_wp .and. in_m > 0.0_wp) then
            share = min(1.0_wp, in_m / out_m)
            u_out = share * in_h / in_m + (1.0_wp - share) * u_out
         end if
         if (pond_flux >= 0.0_wp) then
            amount = h * pond_flux
            flows%infil = flows%infil + amount
            flows%infil_enth = flows%infil_enth + amount * u_mix
            state%e_pond = e_mix - amount * u_mix
            flows%org_enth = flows%org_enth + amount * u_mix
         else
            amount = -h * pond_flux
            flows%exfil = flows%exfil + amount
            flows%exfil_enth = flows%exfil_enth + amount * u_out
            state%e_pond = e_mix + amount * u_out
            flows%org_enth = flows%org_enth - amount * u_out
         end if
         if (q >= 0.0_wp) then
            flows%top_in = flows%top_in + q
            flows%top_in_enth = flows%top_in_enth + q * u_out
            flows%org_enth = flows%org_enth - q * u_out
         else
            flows%top_out = flows%top_out - q
            flows%top_out_enth = flows%top_out_enth - q * internal_energy_liquid(drive%t_soil(1))
            flows%org_enth = flows%org_enth - q * internal_energy_liquid(drive%t_soil(1))
         end if
      end if

      !---~---
      !   Water, by the flows; the eliminated nodes and the leaves at their curves, their neighbours
      ! and the woods across the faces that closed them.
      !---~---
      state%water_soil = state%water_soil + h * (ev%inflow_soil + dinfl)
      state%water_coh  = state%water_coh + h * ev%inflow_coh
      do c = 1_ik, nsoil
         if (snap(c)) state%water_soil(c) = ev%water_soil(c)
      end do
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         state%water_coh(LEAF, i) = ev%water_coh(LEAF, i)
         state%water_coh(WOOD, i) = state%water_coh(WOOD, i) - h * (sap(i) - ev%sapflow(i))
      end do
      state%tau_soil = tau_soil
      state%tau_coh = tau_coh

      !--- Spill: water above the pond's capacity runs off at the pond's temperature.
      spill = max(0.0_wp, state%water_soil(POND) - drive%pond_capacity)
      spilled = spill > 0.0_wp
      if (spilled) then
         t_pond = liquid_energy_to_temp(state%e_pond / state%water_soil(POND))
         flows%spill = flows%spill + spill
         flows%spill_enth = flows%spill_enth + spill * internal_energy_liquid(t_pond)
         state%e_pond = state%e_pond - spill * internal_energy_liquid(t_pond)
         state%water_soil(POND) = drive%pond_capacity
         state%tau_soil(POND) = pond_water_to_head(drive%pond_capacity)
      end if
      if (state%water_soil(POND) <= POND_TINY) then
         state%e_pond = 0.0_wp
      end if

      !--- What each path moved.
      flows%face(1:nsl-1) = flows%face(1:nsl-1) + h * qf(soff+1:soff+nsl-1)
      flows%org_face(1:n_org) = flows%org_face(1:n_org) + h * qf(2:1+n_org)
      flows%org_evap(1:n_org) = flows%org_evap(1:n_org) + h * ev%org_evap(1:n_org)
      flows%drain   = flows%drain   + h * ev%drain
      flows%evap    = flows%evap    + h * ev%evap
      flows%uptake  = flows%uptake  + h * ev%uptake
      flows%sapflow = flows%sapflow + h * sap
      flows%transp  = flows%transp  + h * ev%transp
   contains
      !---~---
      !   The pond gives no more than it holds. Its water is its flows', and Newton, converged to
      ! its tolerance, may leave it drained a hair past empty; the face it drained through then
      ! carries that much less, the node below taking that much less as its own flows give it,
      ! through its slack if it is eliminated, so the water and its enthalpy leave the pond as they
      ! arrive below. Left past empty, the pond's energy would be cut to nothing with the water it
      ! owed: at Toolik 2.9e-5 kg/m2 infiltrated beyond a pond of 0.022, 26 J/m2 out of nowhere.
      !---~---
      subroutine drain_only_what_is_held()
         real(wp) :: past, cut_flux
         past = -(state%water_soil(POND) + h * (ev%inflow_soil(POND) + dinfl(POND)))
         if (past <= 0.0_wp .or. qf(1) <= 0.0_wp) return
         cut_flux    = min(past / h, qf(1))
         qf(1)       = qf(1) - cut_flux
         dinfl(POND) = dinfl(POND) + cut_flux
         dinfl(2)    = dinfl(2) - cut_flux
      end subroutine drain_only_what_is_held

      !--- Whether neighbour gg, not detached, holds the water the slack would take from it.
      logical function can_give(gg)
         integer(ik), intent(in) :: gg
         can_give = .not. cut(gg)
         if (can_give) can_give = state%water_soil(gg) + h * (ev%inflow_soil(gg) + dinfl(gg) - slack) >= 0.0_wp
      end function can_give
   end subroutine water_commit_substep

   pure subroutine clear_flows(flows)
      type(water_flows_t), intent(inout) :: flows
      flows%face = 0.0_wp
      flows%uptake = 0.0_wp
      flows%sapflow = 0.0_wp
      flows%transp = 0.0_wp
      flows%infil = 0.0_wp
      flows%infil_enth = 0.0_wp
      flows%exfil = 0.0_wp
      flows%exfil_enth = 0.0_wp
      flows%drain = 0.0_wp
      flows%evap = 0.0_wp
      flows%spill = 0.0_wp
      flows%spill_enth = 0.0_wp
      flows%rain_enth = 0.0_wp
      flows%org_face = 0.0_wp
      flows%org_evap = 0.0_wp
      flows%org_enth = 0.0_wp
      flows%top_in = 0.0_wp
      flows%top_in_enth = 0.0_wp
      flows%top_out = 0.0_wp
      flows%top_out_enth = 0.0_wp
   end subroutine clear_flows

   !---~---
   !   For the tests: how far the assembled Jacobian is from the residual's own derivative, at
   ! the current state for a sub-step of h -- the largest difference between J v and a
   ! central difference of the residual along v, relative to the largest entry of either. The
   ! drive is set up as water_step sets it.
   !---~---
   real(wp) function water_jacobian_error(static, profile, soil, hyd, drive, water, h) result(err)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(inout) :: drive
      type(water_work_t),        intent(inout) :: water
      real(wp),                  intent(in)    :: h
      type(water_tend_t) :: ev
      real(wp) :: v_soil(water%nsoil), v_coh(COHORT_STORES, water%ncoh)
      real(wp) :: jv_soil(water%nsoil), jv_coh(COHORT_STORES, water%ncoh)
      real(wp) :: rp_soil(water%nsoil), rp_coh(COHORT_STORES, water%ncoh)
      real(wp) :: rm_soil(water%nsoil), rm_coh(COHORT_STORES, water%ncoh)
      real(wp) :: add_soil(water%nsoil), add_coh(COHORT_STORES, water%ncoh)
      real(wp) :: eps, big
      integer(ik) :: c, i, j, l
      ev = water%ev(1)                       ! scratch of the network's shape
      call arrowhead_allocate(water%sys, water%nsoil, water%ncoh, COHORT_STORES)
      call ground_evaporation_setup(profile, soil, drive%ground)
      associate (st => water%state)
         do c = 1_ik, water%nsoil
            v_soil(c) = 1.0e-2_wp * (1.0_wp + 0.1_wp * real(c, wp)) * max(abs(st%tau_soil(c)), 1.0e-3_wp) ! clamp-ok: test step
         end do
         do i = 1_ik, water%ncoh
            do j = 1_ik, COHORT_STORES
               v_coh(j, i) = -1.0e-2_wp * (1.0_wp + 0.07_wp * real(i + j, wp)) * st%tau_coh(j, i)
            end do
         end do
         call water_derivs(static, profile, soil, hyd, drive, st%tau_soil, st%tau_coh, ev, water%sys, .true.)
         call water_storage_diagonal(h, ev, add_soil, add_coh)
         jv_soil = matmul(water%sys%a, v_soil) + add_soil * v_soil
         jv_coh = 0.0_wp
         do i = 1_ik, water%ncoh
            jv_soil = jv_soil + water%sys%b(:, i) * v_coh(WOOD, i)
            do l = 1_ik, COHORT_STORES
               jv_coh(l, i) = sum(water%sys%d(l, :, i) * v_coh(:, i)) + add_coh(l, i) * v_coh(l, i)
            end do
            jv_coh(WOOD, i) = jv_coh(WOOD, i) + sum(water%sys%c(:, i) * v_soil)
         end do
         eps = 1.0e-5_wp
         call water_derivs(static, profile, soil, hyd, drive, st%tau_soil + eps * v_soil, st%tau_coh + eps * v_coh, &
                           ev, water%sys, .false.)
         rp_soil = store_residual(ev%water_soil, st%water_soil, h, ev%inflow_soil)
         rp_coh  = store_residual(ev%water_coh, st%water_coh, h, ev%inflow_coh)
         call water_derivs(static, profile, soil, hyd, drive, st%tau_soil - eps * v_soil, st%tau_coh - eps * v_coh, &
                           ev, water%sys, .false.)
         rm_soil = store_residual(ev%water_soil, st%water_soil, h, ev%inflow_soil)
         rm_coh  = store_residual(ev%water_coh, st%water_coh, h, ev%inflow_coh)
      end associate
      rp_soil = (rp_soil - rm_soil) / (2.0_wp * eps)
      rp_coh = (rp_coh - rm_coh) / (2.0_wp * eps)
      big = max(maxval(abs(jv_soil)), maxval(abs(rp_soil)), tiny_num)
      err = maxval(abs(jv_soil - rp_soil)) / big
      if (water%ncoh > 0_ik) then
         big = max(maxval(abs(jv_coh)), maxval(abs(rp_coh)), tiny_num)
         err = max(err, maxval(abs(jv_coh - rp_coh)) / big)
      end if
   end function water_jacobian_error

end module ed_fast_integ_backward_euler
