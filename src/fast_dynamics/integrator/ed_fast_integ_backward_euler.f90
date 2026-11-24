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
   use ed_constants,   only : tiny_num, cp_air, cp_vap, tsupercool_vap
   use ed_soil_heat,   only : soil_heat_env_t, soil_heat_out_t, soil_energy_step_implicit, set_soil_heat_env
   use ed_ground_surface_lib, only : porous_layer_temp
   use ed_patch_state_vars,   only : soil_energy_state_t
   use ed_fast_integ_derivs, only : surface_derivs, cas_conductances, water_derivs, node_balance, wood_potentials, &
                                    unknowns_from_water
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
                               WORK_WATER_CUT_LAYER, WORK_WATER_FORCED_HELD, WORK_WATER_LINE_SEARCH,              &
                               WORK_WATER_BLOCK_SWEEP, WORK_WATER_HELD_OUT, state_init, surface_source_slopes_t,  &
                               water_porous_t
   use ed_site_state_vars, only : site_surface_t
   use ed_therm_lib, only : cas_enthalpy_to_temp, cas_heat_capacity
   use ed_numerics,     only : sufficient_decrease, fraction_to_boundary, root_bracket_t, regula_falsi_point,       &
                               illinois_narrow
   use ed_soil_hydr_lib,   only : pond_water_to_head, soil_layer_curve_t, curve_kinks, on_curve_corner
   use ed_soil_water,      only : ground_evaporation_setup
   use ed_linear_algebra, only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve, thomas_solve
   use ed_soil_params, only : soil_profile_t, soil_params_t
   use ed_pft_params, only : hydraulics_table_t
   implicit none
   private

   public :: backward_euler_stage, advance_plant_water, advance_surf_water_full
   public :: water_step
   public :: water_jacobian_error, substep_in_bounds, classify_fast_leaves, inert_groups, solve_node   ! used only by the tests
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
   !   [-] The most of its way to zero a canopy-air Newton step may take the canopy air's humidity
   ! or its absolute temperature (fraction to the boundary, as STEP_MAX_SHARE in the water
   ! network's Newton), so every trial the line search evaluates is a state the canopy air can
   ! hold.
   !---~---
   real(wp),    parameter :: SURFACE_STEP_MAX_SHARE = 0.9_wp
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
   !---~---
   !   When none of its BACKTRACK_MAX trials meets Armijo's test or lowers the merit, the line search
   ! halves on BACKTRACK_DEEP more times looking for any fall, to 2^-39 of the step: a step past
   ! which only round-off moves a smooth merit.
   !---~---
   integer(ik), parameter :: BACKTRACK_DEEP = 28_ik
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
   integer(ik), parameter :: NODE_SOLVE_MAX   = 200_ik
   !---~---
   !   The classes' hysteresis: a node keeps the class it had at the last accepted sub-step until
   ! the measure that put it there has moved CLASS_HYSTERESIS times past the threshold, so a node
   ! sitting at a threshold does not change class at round-off. A numerical choice: the classes
   ! change how a balance is solved, not what solves it.
   !---~---
   real(wp),    parameter :: CLASS_HYSTERESIS = 2.0_wp
   !---~---
   !   solve_node's bracket search: its first step Newton's, but at least NODE_STEP_MIN of the
   ! unknown (or of one), so it moves the unknown however flat the balance; how far one step may
   ! go towards zero on the dry side, to NODE_DRY_SHARE of the unknown, three decades of
   ! saturation a step, so a dry node crosses many decades geometrically and never steps to zero;
   ! and the steps out it takes before it gives up, NODE_BRACKET_MAX, each NODE_STEP_GROWTH times
   ! the last, which reach 4^40 (about 1e24) times the first, beyond any range of the unknown, so a
   ! balance not bracketed then has no change of sign. In the bracket, a spread wider than
   ! NODE_SPREAD_MAX (two decades) is halved geometrically, where the arithmetic midpoint would
   ! take about three points a decade, and a bracket narrower than NODE_ROUNDOFF epsilons of the
   ! unknown has closed to its last bits. A block of eliminated nodes takes at most
   ! BLOCK_NEWTON_MAX Newton steps or sweeps: a smooth block converges in under ten. Numerical
   ! choices: none moves a converged answer.
   !---~---
   real(wp),    parameter :: NODE_STEP_MIN    = 1.0e-3_wp
   real(wp),    parameter :: NODE_DRY_SHARE   = 1.0e-3_wp
   integer(ik), parameter :: NODE_BRACKET_MAX = 40_ik
   real(wp),    parameter :: NODE_SPREAD_MAX  = 1.0e2_wp
   real(wp),    parameter :: NODE_ROUNDOFF    = 8.0_wp
   real(wp),    parameter :: NODE_STEP_GROWTH = 4.0_wp
   integer(ik), parameter :: BLOCK_NEWTON_MAX = 50_ik
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
      real(wp)    :: enth1, shv1, e_drain, t_cas1, water_heat(nsl)
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
      !   The enthalpy of the water the network moved: across the soil's top in g_top, the ground
      ! surface's stack having valued it at its giver's temperature; across the faces below at its
      ! time-mean flows; out with root uptake at each layer's start-of-step temperature (what the
      ! wood was credited), and out with drainage (docs/science/soil_plant_water.md §5). The leaf
      ! pays its transpiration's full vapour enthalpy in surface_derivs.
      !---~---
      call set_soil_heat_env(nsl, surf_tend%g_top, y%theta, frozen%roots%root_heat,                              &
                             frozen%hydrology%w_flux_frozen, frozen%hydrology%drainage, frozen%hydrology%t_bot,    &
                             eforc, water_heat, e_drain)
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
            bf%soil_enth_in = surf_tend%g_top - eflux%bottom_heat
            bf%soil_enth_out= sum(water_heat(1:nsl)) + e_drain
            !---~---
            !   Soil water's storage change and flows are the network's, booked once per step by
            ! close_ledgers: the per-stage soil-water amounts are zero, and the whole-column
            ! water here is only the CAS-vapour exchange.
            !---~---
            bf%soil_wat_in  = 0.0_wp
            bf%soil_wat_out = 0.0_wp
            !---~---
            !   ground_rad is the tiles' radiation, cover-weighted (rnet_soil when bare). The water the
            ! network moves and the uptake carry enthalpy between tracked stores and cancel; rain and
            ! runoff are added once at the outer level (docs/science/numerical_scheme.md §5b).
            !---~---
            bf%whole_enth_in= surf_tend%coh_rnet + surf_tend%ground_rad - eflux%bottom_heat
            !---~---
            !   The condensate is not a boundary loss: the update puts it in the water layer with the
            ! enthalpy surface_derivs debited from the CAS (whole_cond, whole_cond_enth).
            !---~---
            bf%whole_enth_out= g_atm_heat*(enth1 - frozen%cas%box%atm_enthalpy) + e_drain
            bf%whole_wat_in = 0.0_wp
            bf%whole_wat_out = g_atm_vapour*(shv1 - frozen%cas%box%atm_specific_humidity)
            bf%whole_cond   = surf_tend%cond                     ! put in the water layer, not lost
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
         lam = surface_domain_step(Hk, qk, delH, delq)
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
   !   [-] The share of a canopy-air Newton step (dh, dq from the state h, q) the line search starts
   ! from: the whole step, or less where it would take the canopy air's humidity or its absolute
   ! temperature more than SURFACE_STEP_MAX_SHARE of its way to zero. The temperature along the
   ! step, T = (h + q c_v T_v) / c_p(q) (cas_enthalpy_to_temp), is a ratio of two functions linear
   ! in the share, so the share at which it reaches its floor is in closed form.
   !---~---
   pure real(wp) function surface_domain_step(h, q, dh, dq) result(lam)
      real(wp), intent(in) :: h, q, dh, dq   ! [J/kg], [kg/kg] the state and the step
      real(wp) :: num0, dnum, den0, dden, t_floor
      lam = min(1.0_wp, fraction_to_boundary(q, dq, SURFACE_STEP_MAX_SHARE))
      num0 = h + q * cp_vap * tsupercool_vap
      dnum = dh + dq * cp_vap * tsupercool_vap
      den0 = cas_heat_capacity(q)
      dden = dq * (cp_vap - cp_air)
      t_floor = (1.0_wp - SURFACE_STEP_MAX_SHARE) * num0 / den0
      if (dnum - t_floor * dden < 0.0_wp) lam = min(lam, (t_floor * den0 - num0) / (dnum - t_floor * dden))
   end function surface_domain_step

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
   !   A sub-step forced at the shortest length commits a solution of its balance (Newton
   ! converged, only the error filter failed) as any other. A sub-step whose Newton fails at the
   ! shortest length is not taken: the core node whose balance was worst is held out of the
   ! network for the rest of the step, as an inert node is detached (inert_groups) -- its faces
   ! carry nothing and it takes only its own sources, the pond the water reaching it, a porous or
   ! soil layer its evaporation -- and the sub-step is retaken from its start at the rest of the
   ! step. The rest of the network, the plants' supply from the soil among it, goes on. Only when
   ! the worst store is a cohort's, or a node already held out, is an unconverged Newton iterate
   ! committed, if what it commits is a state of the network (substep_in_bounds): every store
   ! finite and at or above empty, its water within the error tolerance of its curve's at its
   ! unknown; otherwise the network holds still for the rest of the step, every store as it began
   ! the sub-step but the pond, which takes the water reaching it. Either way each unknown is then
   ! read again from its committed water, so the next sub-step starts on its curve.
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
      integer(ik) :: trials, iters, i_sol, worst_at, n_unsolved, n_search, n_sweep
      logical     :: converged, last, spilled, forced, kept
      logical     :: reread   ! the next sub-step's unknowns are read again from the water committed
      logical     :: fast(water%nsoil), cut(water%nsoil)   ! the sub-step's eliminated and detached core nodes
      logical     :: held_out(water%nsoil)                 ! the core nodes held out for the rest of the step
      logical     :: was_fast(water%nsoil), was_cut(water%nsoil)   ! the classes of the last accepted sub-step
      logical     :: fast_leaf(water%ncoh), was_fast_leaf(water%ncoh)   ! each cohort's leaves fast, now and then
      type(water_state_t) :: start_state                   ! a forced sub-step's start, to hold still at
      type(water_flows_t) :: start_flows

      rtol = tols%rtol(GRP_THETA)
      atol = tols%atol(GRP_THETA)
      call clear_flows(water%flows)
      work_count(WORK_WATER_NSUB:WORK_WATER_NFAIL) = 0_ik
      work_count(WORK_WATER_FORCED_NEWTON:WORK_WATER_CUT_LAYER) = 0_ik
      work_count(WORK_WATER_FORCED_HELD:WORK_WATER_HELD_OUT) = 0_ik
      call arrowhead_allocate(water%sys, water%nsoil, water%ncoh, COHORT_STORES)
      call water_scales(static, profile, hyd, drive%porous, scale_soil, scale_coh)
      call ground_evaporation_setup(profile, soil, drive%ground)
      water%start_ready = .false.
      forced = .false.
      reread = .false.
      held_out = .false.
      was_fast = .false.
      was_cut  = .false.
      was_fast_leaf = .false.

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
                                   tau_soil, tau_coh, held_out, was_fast, was_cut, was_fast_leaf, water, converged,     &
                                   iters, i_sol, fast, cut, fast_leaf, worst_at, n_unsolved, n_search, n_sweep)
         work_count(WORK_WATER_NEWTON) = work_count(WORK_WATER_NEWTON) + iters
         work_count(WORK_WATER_NODE_UNSOLVED) = work_count(WORK_WATER_NODE_UNSOLVED) + n_unsolved
         work_count(WORK_WATER_LINE_SEARCH) = work_count(WORK_WATER_LINE_SEARCH) + n_search
         work_count(WORK_WATER_BLOCK_SWEEP) = work_count(WORK_WATER_BLOCK_SWEEP) + n_sweep
         if (.not. converged .and. .not. last) then
            work_count(WORK_WATER_REJECT) = work_count(WORK_WATER_REJECT) + 1_ik
            h = max(h * NEWTON_FAIL_CUT, h_min)
            cycle
         end if
         !--- Failed at the shortest length: the worst core node held out, and the sub-step retaken at the rest of the step.
         if (.not. converged .and. worst_at >= 1_ik .and. worst_at <= water%nsoil) then
            if (.not. held_out(worst_at)) then
               held_out(worst_at) = .true.
               work_count(WORK_WATER_HELD_OUT) = work_count(WORK_WATER_HELD_OUT) + 1_ik
               work_count(WORK_WATER_REJECT) = work_count(WORK_WATER_REJECT) + 1_ik
               h = dt - t
               cycle
            end if
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
         forced = err > 1.0_wp
         if (forced) call count_forced(converged, worst_at)
         call count_classes(fast, WORK_WATER_FAST_POND, WORK_WATER_FAST_SLOT, WORK_WATER_FAST_LAYER)
         call count_classes(cut, WORK_WATER_CUT_POND, WORK_WATER_CUT_SLOT, WORK_WATER_CUT_LAYER)
         was_fast = fast
         was_cut  = cut .and. .not. held_out
         was_fast_leaf = fast_leaf
         if (converged) then
            !--- A sub-step's balance that is not finite is never committed: the run stops.
            if (.not. finite_evaluation(water%ev(i_sol), tau_soil, tau_coh))                                  &
               error stop 'water_step: a sub-step of the soil-plant water network is not finite'
            call water_commit_substep(drive, h, tau_soil, tau_coh, water%ev(i_sol), water%state, water%flows,  &
                                      fast, cut, fast_leaf, spilled)
         else
            start_state = water%state
            start_flows = water%flows
            kept = finite_evaluation(water%ev(i_sol), tau_soil, tau_coh)
            if (kept) then
               call water_commit_substep(drive, h, tau_soil, tau_coh, water%ev(i_sol), water%state, water%flows, &
                                         fast, cut, fast_leaf, spilled)
               kept = substep_in_bounds(water%state, water%ev(i_sol), scale_soil, scale_coh, rtol, atol)
            end if
            if (.not. kept) then
               !---~---
               !   Held still: every store stays as it began the sub-step, but the pond, which takes
               ! the water reaching it whatever the rest does (and spills it past its capacity). The
               ! start held would fail again as it just did, so the hold lasts the rest of the step.
               !---~---
               water%state = start_state
               water%flows = start_flows
               h = dt - t
               call hold_pond(h)
               work_count(WORK_WATER_FORCED_HELD) = work_count(WORK_WATER_FORCED_HELD) + 1_ik
            end if
         end if
         !---~---
         !   The next sub-step starts where this one ended: its start is this one's solution, whose
         ! Jacobian sys still holds -- unless the spill moved the pond, or a node was detached, whose
         ! cut faces the next sub-step's classes must not read, or the sub-step was taken although its
         ! Newton had not converged, whose unknowns are read again from the water it committed. A
         ! sub-step forced by the error filter alone is a converged solution, its unknowns its water's:
         ! read again, an eliminated node's unknown would lose the balance its own solve gave it, and
         ! the next sub-step's error would see that balance's flow return at every length.
         !---~---
         water%i_start = i_sol
         reread = forced .and. .not. converged
         water%start_ready = .not. (spilled .or. any(cut) .or. reread)
         if (reread) call unknowns_from_water(static, profile, hyd, drive%ice(1:water%nsl), drive%porous, water%state)
         t = t + h
         work_count(WORK_WATER_NSUB) = work_count(WORK_WATER_NSUB) + 1_ik
         if (err <= 1.0_wp)                                                                                   &
            h = h * min(GROW_MAX, max(SHRINK_MIN, SAFETY * max(err, 1.0e-10_wp) ** (-0.5_wp))) ! clamp-ok: controller floor
      end do
      !--- The potentials an unconverged last sub-step left are its unknowns' as read again from its water.
      if (reread)                                                                                             &
         call water_derivs(static, profile, soil, hyd, drive, water%state%tau_soil, water%state%tau_coh,         &
                           water%ev(water%i_start), water%sys, .false.)
      water%state%psi_soil = water%ev(water%i_start)%psi_soil
      water%state%psi_coh  = water%ev(water%i_start)%psi_coh
   contains
      !--- The pond over a held span hh: the water reaching it, and its spill past its capacity.
      subroutine hold_pond(hh)
         real(wp), intent(in) :: hh   ! [s]
         real(wp) :: spill
         water%state%water_soil(POND) = water%state%water_soil(POND) + hh * drive%rain
         spill = max(0.0_wp, water%state%water_soil(POND) - drive%pond_capacity)
         spilled = spill > 0.0_wp
         water%flows%spill = water%flows%spill + spill
         water%state%water_soil(POND) = water%state%water_soil(POND) - spill
      end subroutine hold_pond

      !--- An accepted sub-step's nodes of a class, by kind: the pond, the porous layers, the soil layers.
      subroutine count_classes(mask, i_pond, i_slot, i_layer)
         logical,     intent(in) :: mask(:)
         integer(ik), intent(in) :: i_pond, i_slot, i_layer   ! their counters in work_count
         if (mask(POND)) work_count(i_pond) = work_count(i_pond) + 1_ik
         work_count(i_slot)  = work_count(i_slot)  + count(mask(2:1+drive%porous%n), kind=ik)
         work_count(i_layer) = work_count(i_layer) + count(mask(2+drive%porous%n:), kind=ik)
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

   !---~---
   !   Whether what a forced sub-step's unconverged iterate committed is a state of the network: every
   ! store's water finite and at or above empty, and within the error tolerance (atol of its scale
   ! plus rtol of its water) of what its curve holds at the iterate's unknown, the amount by which
   ! the iterate's balance missed. A store below empty or past what its curve holds there is the
   ! sign of an iterate gone wild.
   !---~---
   pure logical function substep_in_bounds(state, ev, scale_soil, scale_coh, rtol, atol) result(ok)
      type(water_state_t), intent(in) :: state   ! the water the iterate committed
      type(water_tend_t),  intent(in) :: ev      ! the network at the iterate
      real(wp),            intent(in) :: scale_soil(:), scale_coh(:,:), rtol, atol
      ok = all(ieee_is_finite(state%water_soil)) .and. all(ieee_is_finite(state%water_coh))
      if (.not. ok) return
      ok = all(state%water_soil >= 0.0_wp) .and. all(state%water_coh >= 0.0_wp)                               &
           .and. all(abs(state%water_soil - ev%water_soil) <= atol * scale_soil + rtol * abs(state%water_soil))  &
           .and. all(abs(state%water_coh - ev%water_coh) <= atol * scale_coh + rtol * abs(state%water_coh))
   end function substep_in_bounds

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
   ! its evaluation, a node whose time scale is far below the sub-step is fast (classify_fast), and
   ! adjacent fast nodes that relax at once make one block (fast_blocks); every
   ! evaluation first solves each block's balances together with its neighbours held (solve_block,
   ! nonlinear elimination), and Newton keeps their rows with a zero residual, so the factorised
   ! matrix gives the rest the Schur-eliminated step. At every iteration, a connected group of
   ! nodes the solve cannot see is detached (inert_groups), as is every node water_step holds out
   ! (held_out): its faces carry nothing and each member solves its own balance with its own
   ! sources. A step that would carry a node across a kink of its curve stops there for the
   ! iteration (chop_at_kinks). The line search takes no step that raises the merit: when every
   ! trial would, Newton stops unconverged (n_search). The node balances left unsolved and the
   ! blocks solved by sweeps are counted (n_unsolved, n_sweep).
   !---~---
   subroutine water_newton_substep(static, profile, soil, hyd, drive, h, atol, scale_soil, scale_coh,             &
                                   tau_soil, tau_coh, held_out, was_fast, was_cut, was_fast_leaf, water, converged,     &
                                   iters, i_sol, fast, cut, fast_leaf, worst_at, n_unsolved, n_search, n_sweep)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(in)    :: drive
      real(wp),                  intent(in)    :: h, atol, scale_soil(:), scale_coh(:,:)
      real(wp),                  intent(inout) :: tau_soil(:), tau_coh(:,:)   ! start of sub-step in, solution out
      logical,                   intent(in)    :: held_out(:)   ! (nsoil) core nodes held out, detached whatever their links
      logical,                   intent(in)    :: was_fast(:), was_cut(:)   ! (nsoil) the last accepted sub-step's classes
      logical,                   intent(in)    :: was_fast_leaf(:)          ! (ncoh) its cohorts' leaves' class
      type(water_work_t),        intent(inout) :: water
      logical,                   intent(out)   :: converged
      integer(ik),               intent(out)   :: iters
      integer(ik),               intent(out)   :: i_sol   ! which of water%ev holds the solution
      logical,                   intent(out)   :: fast(:), cut(:)   ! (nsoil) the eliminated and detached nodes
      logical,                   intent(out)   :: fast_leaf(:)      ! (ncoh) the cohorts whose leaves are fast
      integer(ik),               intent(out)   :: worst_at     ! the store of the last worst imbalance (water_residual)
      integer(ik),               intent(out)   :: n_unsolved   ! node balances left unsolved
      integer(ik),               intent(out)   :: n_search     ! line searches that found no step lowering the merit
      integer(ik),               intent(out)   :: n_sweep      ! block solves that fell back to sweeps
      real(wp) :: res_soil(water%nsoil), res_coh(COHORT_STORES, water%ncoh)
      real(wp) :: step_soil(water%nsoil), step_coh(COHORT_STORES, water%ncoh)
      real(wp) :: try_soil(water%nsoil), try_coh(COHORT_STORES, water%ncoh)
      real(wp) :: res2_soil(water%nsoil), res2_coh(COHORT_STORES, water%ncoh)
      real(wp) :: merit, merit_try, merit_best, worst, worst_try, lam, lam_max, lam_best, tol
      integer(ik) :: it, ib, i_try, worst_at_try
      real(wp) :: step_free(water%nsoil)   ! the step before chop_at_kinks
      real(wp) :: t_relax(water%nsoil)     ! [s] each core node's time to relax to its neighbours
      logical  :: at_once(water%nsoil)     ! each core node that relaxes at once (classify_fast)
      logical  :: ok, cut_new(water%nsoil), taken
      logical  :: chopped                  ! chop_at_kinks shortened the step

      converged = .false.
      iters = 0_ik
      worst_at = 0_ik
      n_unsolved = 0_ik
      n_search = 0_ik
      n_sweep = 0_ik
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
      call classify_fast(drive, tau_soil, water%sys, water%ev(i_sol), h, tol, scale_soil, was_fast, fast, t_relax,   &
                         at_once)
      call classify_fast_leaves(water%sys, water%ev(i_sol), h, atol, scale_coh, was_fast_leaf, fast_leaf)
      call inert_groups(water%ev(i_sol), h, tol, scale_soil, was_cut, held_out, cut)
      cut  = cut .or. held_out
      fast = fast .and. .not. cut
      call fast_blocks(t_relax, at_once, fast)
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
         call inert_groups(water%ev(i_sol), h, tol, scale_soil, cut .and. .not. held_out, held_out, cut_new)
         cut_new = cut_new .or. held_out
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
         !--- A trial's evaluation replaces the Jacobian sys held, so the start's is gone.
         water%start_ready = .false.
         call arrowhead_solve(water%sys, merge(0.0_wp, -res_soil, fast .or. cut), -res_coh, step_soil, step_coh)
         !--- An eliminated or detached node's own solve, not the linear step, gives its unknown.
         step_soil = merge(0.0_wp, step_soil, fast .or. cut)
         step_free = step_soil
         call chop_at_kinks(profile, drive, tau_soil, step_soil)
         chopped = any(step_soil /= step_free)
         lam_max = step_limit(profile, drive%porous, tau_soil, tau_coh, step_soil, step_coh)
         lam = lam_max
         taken = .false.
         merit_best = huge(1.0_wp)
         lam_best = 0.0_wp
         do ib = 1_ik, BACKTRACK_MAX + BACKTRACK_DEEP
            call trial(lam)
            if (merit_try == merit_try) then
               if (ib <= BACKTRACK_MAX) then
                  if (sufficient_decrease(merit_try, merit, lam)) then
                     taken = .true.
                     exit
                  end if
                  !---~---
                  !   A step chopped at a kink is not Newton's, and the decrease Armijo's test
                  ! expects of Newton's it cannot give: landing a node on a kink a hair away moves
                  ! the merit by a hair. It is taken if the merit does not grow, so the next step
                  ! may cross.
                  !---~---
                  if (chopped .and. merit_try <= merit) then
                     taken = .true.
                     exit
                  end if
                  if (merit_try < merit_best) then
                     merit_best = merit_try
                     lam_best   = lam
                  end if
               else if (merit_try < merit) then
                  taken = .true.
                  exit
               end if
            end if
            if (ib == BACKTRACK_MAX .and. merit_best < merit) exit
            lam = 0.5_wp * lam
         end do
         !---~---
         !   No trial met Armijo's test: the one of them that lowered the merit most is taken, else
         ! the first shorter one that lowers it at all, the merit falling along Newton's direction
         ! for every short enough step where the balances are smooth; if none does, Newton stops
         ! here, unconverged, rather than take a step that raises it.
         !---~---
         if (.not. taken) then
            if (.not. (merit_best < merit)) then
               n_search = n_search + 1_ik
               if (i_sol == i_try) call trial(0.0_wp)   ! the iterate's evaluation back in its slot
               return
            end if
            call trial(lam_best)
         end if
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

      !--- The network at lam of the Newton step: its residuals and their measures.
      subroutine trial(share)
         real(wp), intent(in) :: share
         try_soil = tau_soil + share * step_soil
         try_coh = tau_coh + share * step_coh
         call evaluate(try_soil, try_coh, i_try)
         call water_residual(water%state, h, water%ev(i_try), atol, scale_soil, scale_coh, res2_soil,             &
                             res2_coh, merit_try, worst_try, worst_at_try, cut)
      end subroutine trial

      !---~---
      !   The network at an iterate: each detached node's unknown from its own balance (solve_node)
      ! and each block of eliminated nodes' from their balances together (solve_block), the others'
      ! held, then every store's flows and the Jacobian, the detached nodes' faces cut.
      !---~---
      subroutine evaluate(ts, tc, i_ev)
         real(wp),    intent(inout) :: ts(:), tc(:,:)
         integer(ik), intent(in)    :: i_ev
         real(wp)    :: psi_wood(water%ncoh), node_tol(water%nsoil)
         integer(ik) :: c, c2
         logical     :: solved, swept
         if (any(fast) .or. any(cut)) then
            call wood_potentials(static, hyd, tc, psi_wood)
            node_tol = NODE_SOLVE_SHARE * NEWTON_TOL * atol * scale_soil / h
            c = 1_ik
            do while (c <= water%nsoil)
               if (cut(c)) then
                  call solve_node(static, profile, soil, drive, c, ts, water%state%water_soil(c), h, node_tol(c),   &
                                  cut, psi_wood, solved)
                  if (.not. solved) n_unsolved = n_unsolved + 1_ik
                  c = c + 1_ik
               else if (fast(c)) then
                  c2 = c
                  do while (c2 < water%nsoil)
                     if (.not. fast(c2 + 1_ik)) exit
                     c2 = c2 + 1_ik
                  end do
                  call solve_block(static, profile, soil, drive, c, c2, ts, water%state%water_soil, h, node_tol,    &
                                   cut, psi_wood, solved, swept)
                  if (.not. solved) n_unsolved = n_unsolved + 1_ik
                  if (swept) n_sweep = n_sweep + 1_ik
                  c = c2 + 1_ik
               else
                  c = c + 1_ik
               end if
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
   ! its flows, and it is eliminated. So is a node whose net outflow at the start would empty it
   ! within the sub-step, its relaxation time taken as nil: its storage there, a trace pond's over
   ! its depression, says nothing of the empty store it becomes, which a Newton step linearised at
   ! the start would carry past the head where its face meets the node below. So, too, is a node
   ! that starts on a corner of its curve (on_curve_corner): a Newton step linearised there holds
   ! one side of the corner only, its own bracketed solve both. Of a run of adjacent fast nodes,
   ! those that relax at once are eliminated together, as blocks (fast_blocks). A node eliminated
   ! at the last accepted sub-step (was_fast) stays so until its time scale passes
   ! CLASS_HYSTERESIS times the threshold, and its link is seen at a CLASS_HYSTERESIS-th of the
   ! tolerance (seen_links).
   !---~---
   pure subroutine classify_fast(drive, tau_soil, sys, ev, h, tol, scale_soil, was_fast, fast, t_relax, at_once)
      type(water_drive_t),      intent(in)  :: drive
      real(wp),                 intent(in)  :: tau_soil(:)
      type(arrowhead_system_t), intent(in)  :: sys
      type(water_tend_t),       intent(in)  :: ev
      real(wp),                 intent(in)  :: h, tol, scale_soil(:)
      logical,                  intent(in)  :: was_fast(:)
      logical,                  intent(out) :: fast(:)
      real(wp),                 intent(out) :: t_relax(:)   ! [s] each node's time to relax to its neighbours
      logical,                  intent(out) :: at_once(:)   ! it holds no water where it starts, or starts on a corner
      real(wp)    :: big_l
      integer(ik) :: c, i, ib
      logical     :: seen(sys%nc), linked
      call seen_links(ev, tol, scale_soil, 1.0_wp / face_hysteresis(was_fast(1:sys%nc)), seen)
      do c = 1_ik, sys%nc
         big_l = 0.0_wp
         do i = 1_ik, sys%nc
            big_l = big_l + abs(sys%a(i, c))
         end do
         do ib = 1_ik, sys%nb
            big_l = big_l + abs(sys%c(c, ib))
         end do
         linked = seen(c) .or. ev%root_link(c) >= tol * scale_soil(c) / hysteresis(was_fast(c))
         if (c > 1_ik) linked = linked .or. seen(c - 1_ik)
         t_relax(c) = abs(ev%dwater_soil(c)) / max(big_l, tiny_num)
         at_once(c) = ev%dwater_soil(c) == 0.0_wp
         if (c > 1_ik .and. c <= 1_ik + drive%porous%n) then
            at_once(c) = at_once(c) .or. on_curve_corner(drive%porous%ice(c - 1_ik), tau_soil(c))
         else if (c > 1_ik + drive%porous%n) then
            at_once(c) = at_once(c) .or. on_curve_corner(drive%ice(c - 1_ik - drive%porous%n), tau_soil(c))
         end if
         !--- So does one that empties within the sub-step.
         if (at_once(c) .or. -h * ev%inflow_soil(c) > ev%water_soil(c)) t_relax(c) = 0.0_wp
         fast(c) = linked .and. t_relax(c) < FAST_SHARE * h * hysteresis(was_fast(c))
      end do
   end subroutine classify_fast

   !--- [-] The factor by which a class's threshold moves for a node that held the class (CLASS_HYSTERESIS), else 1.
   elemental real(wp) function hysteresis(held) result(f)
      logical, intent(in) :: held
      f = merge(CLASS_HYSTERESIS, 1.0_wp, held)
   end function hysteresis

   !--- [-] Each face's hysteresis (hysteresis): face c's when either of its nodes, c and c + 1, held the class.
   pure function face_hysteresis(held) result(f)
      logical, intent(in) :: held(:)
      real(wp) :: f(size(held))
      integer(ik) :: c
      do c = 1_ik, size(held, kind=ik)
         f(c) = hysteresis(held(c))
         if (c < size(held, kind=ik)) f(c) = hysteresis(held(c) .or. held(c + 1_ik))
      end do
   end function face_hysteresis

   !---~---
   !   Each cohort's leaves at a sub-step's start, classified as a core node is (classify_fast): their
   ! time to relax, their storage over the flow part of their Jacobian column (the sapflow's slope
   ! and the transpiration's), fast below FAST_SHARE of the sub-step, with the same hysteresis.
   ! Leaves relax at once that hold no water to store (leafless), that their outflow at the start
   ! would empty within the sub-step, or that hold less than their cohort's absolute tolerance: the
   ! cohort's balance is solved to a share of the water its whole plant holds (water_scales), so the
   ! water such leaves take from their flows is not resolved and can leave their curve, even past
   ! empty. A fast leaf takes its water from its curve at the commit, the sapflow closing its
   ! balance; a slower one keeps the water its flows give it, as any store does.
   !---~---
   pure subroutine classify_fast_leaves(sys, ev, h, atol, scale_coh, was_fast, fast)
      type(arrowhead_system_t), intent(in)  :: sys
      type(water_tend_t),       intent(in)  :: ev
      real(wp),                 intent(in)  :: h, atol
      real(wp),                 intent(in)  :: scale_coh(:,:)   ! (2, ncoh) [kg/m2] each cohort's water scale
      logical,                  intent(in)  :: was_fast(:)      ! (ncoh) the last accepted sub-step's
      logical,                  intent(out) :: fast(:)          ! (ncoh)
      real(wp)    :: big_l, t_relax
      integer(ik) :: i
      do i = 1_ik, sys%nb
         big_l   = abs(sys%d(LEAF, LEAF, i)) + abs(sys%d(WOOD, LEAF, i))
         t_relax = abs(ev%dwater_coh(LEAF, i)) / max(big_l, tiny_num)
         if (ev%dwater_coh(LEAF, i) == 0.0_wp .or. -h * ev%inflow_coh(LEAF, i) > ev%water_coh(LEAF, i)        &
             .or. ev%water_coh(LEAF, i) < atol * scale_coh(LEAF, i) * hysteresis(was_fast(i))) t_relax = 0.0_wp
         fast(i) = t_relax < FAST_SHARE * h * hysteresis(was_fast(i))
      end do
   end subroutine classify_fast_leaves

   !---~---
   !   The blocks of eliminated nodes (docs/science/soil_plant_water.md §3, "Node classes"). In a
   ! run of adjacent fast nodes, the members that relax at once (at_once, from classify_fast: their
   ! curve holds no water where they start, as an empty pond's or a frozen layer's whose pores ice
   ! and liquid fill, below zero head; or they start on a corner of their curve, on_curve_corner)
   ! stay eliminated, adjacent ones together as one block (solve_block): no Newton row of their own
   ! could hold them. The other members leave the run for the Newton, whose rows hold them; a run of
   ! such members alone keeps its fastest one eliminated (a node that empties within the sub-step
   ! first), so no eliminated node is solved against a neighbour held that moves as fast.
   !---~---
   pure subroutine fast_blocks(t_relax, at_once, fast)
      real(wp),                 intent(in)    :: t_relax(:)
      logical,                  intent(in)    :: at_once(:)
      logical,                  intent(inout) :: fast(:)
      integer(ik) :: c1, c2, nc
      nc = size(fast, kind=ik)
      c1 = 1_ik
      do while (c1 <= nc)
         if (.not. fast(c1)) then
            c1 = c1 + 1_ik
            cycle
         end if
         c2 = c1
         do while (c2 < nc)
            if (.not. fast(c2 + 1_ik)) exit
            c2 = c2 + 1_ik
         end do
         if (c2 > c1) then
            if (any(at_once(c1:c2))) then
               fast(c1:c2) = at_once(c1:c2)
            else
               fast(c1:c2) = .false.
               fast(c1 - 1_ik + minloc(t_relax(c1:c2), dim=1, kind=ik)) = .true.
            end if
         end if
         c1 = c2 + 1_ik
      end do
   end subroutine fast_blocks

   !---~---
   !   The chain's faces the solve can see: face c, joining nodes c and c + 1, is seen when a whole
   ! unit of either end's unknown moves either end's balance by the convergence tolerance (each
   ! weighed as water_residual weighs that balance), from its slopes whether or not it is cut, the
   ! tolerance moved by each face's factor (the classes' hysteresis).
   !---~---
   pure subroutine seen_links(ev, tol, scale_soil, factor, seen)
      type(water_tend_t), intent(in)  :: ev
      real(wp),           intent(in)  :: tol, scale_soil(:)
      real(wp),           intent(in)  :: factor(:)   ! (nc) [-] face c's tolerance over the convergence tolerance
      logical,            intent(out) :: seen(:)
      integer(ik) :: c, nc
      nc = size(seen, kind=ik)
      do c = 1_ik, nc - 1_ik
         seen(c) = max(abs(ev%face_dqa(c)), abs(ev%face_dqb(c)))                                               &
                   >= tol * min(scale_soil(c), scale_soil(c + 1_ik)) * factor(c)
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
   ! ground, a frozen layer whose pores are full, and those with frozen porous layers between. A
   ! node detached before (was_cut) rejoins only once its storage or a link out is
   ! CLASS_HYSTERESIS times what the solve can see. A node held out for the rest of the step
   ! (held_out) is outside the network: its faces carry nothing, so they join no group, and its
   ! storage makes none active; a node whose only link was to it, an empty pond over a layer held
   ! out, is a group of its own, and detached if it holds nothing.
   !---~---
   pure subroutine inert_groups(ev, h, tol, scale_soil, was_cut, held_out, cut)
      type(water_tend_t), intent(in)  :: ev
      real(wp),           intent(in)  :: h, tol, scale_soil(:)
      logical,            intent(in)  :: was_cut(:)
      logical,            intent(in)  :: held_out(:)   ! the nodes held out for the rest of the step
      logical,            intent(out) :: cut(:)
      integer(ik) :: nc, c, first
      logical     :: seen_link(size(cut)), active
      real(wp)    :: seen_at
      nc = size(cut, kind=ik)
      call seen_links(ev, tol, scale_soil, face_hysteresis(was_cut(1:nc)), seen_link)
      do c = 1_ik, nc - 1_ik
         if (held_out(c) .or. held_out(c + 1_ik)) seen_link(c) = .false.
      end do
      cut = .false.
      !--- The chain's groups: runs of nodes joined by seen faces.
      first = 1_ik
      active = .false.
      do c = 1_ik, nc
         seen_at = tol * scale_soil(c) * hysteresis(was_cut(c))
         if (.not. held_out(c) .and. (abs(ev%dwater_soil(c)) / h >= seen_at .or. ev%root_link(c) >= seen_at))      &
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
   ! held, by a bracketed root search. The balance rises with the unknown (more water held, more
   ! flowing out), so the root lies on the side its sign points to: the search steps out that way
   ! from the current unknown, each step four times the last, Newton's at first (geometric towards
   ! zero on the dry side, the unknown being a saturation that may span many decades), so of
   ! several roots the nearest is taken. In the bracket it takes Newton's step where that lands
   ! inside and the last one at least halved the residual, else the Illinois method's point
   ! (ed_numerics), which closes in from both ends; it halves a bracket that has not halved over
   ! two points (a balance flat on one side, as an empty pond's over a dry litter, holds the
   ! Illinois point at that side), geometrically where the bracket spans decades. It is solved
   ! when the residual is within tol [kg/m2/s], far inside Newton's, or when the bracket has
   ! closed to round-off of the unknown with the residual within Newton's tolerance, the best the
   ! arithmetic gives; solved says whether it was, the unknown left at the best point found
   ! either way.
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
      type(root_bracket_t) :: br
      real(wp)    :: x, r, dr, x_new, dx, x_best, r_best, dr_best, r_last, width
      integer(ik) :: it, stalled
      logical     :: bracketed, newton_ok
      x = ts(c)
      x_best = x
      r_best = huge(1.0_wp)
      dr_best = 0.0_wp
      call balance_at(x, r, dr)
      if (solved) return
      !--- Step out towards the root until the residual changes sign.
      br%a = x
      br%fa = r
      bracketed = .false.
      dx = -sign(NODE_STEP_MIN * max(abs(x), 1.0_wp), r)
      if (dr > 0.0_wp) then
         if (abs(r / dr) > abs(dx)) dx = -r / dr
      end if
      do it = 1_ik, NODE_BRACKET_MAX
         x_new = x + dx
         if (c /= POND) x_new = max(x_new, NODE_DRY_SHARE * x, SOIL_TAU_FLOOR)
         if (x_new == x) exit                           ! at the dry floor, the residual still of one sign
         x = x_new
         call balance_at(x, r, dr)
         if (solved) then
            ts(c) = x
            return
         end if
         if (r * br%fa < 0.0_wp) then
            bracketed = .true.
            exit
         end if
         br%a  = x
         br%fa = r
         dx = NODE_STEP_GROWTH * dx
      end do
      if (.not. bracketed) then
         ts(c) = x_best
         return
      end if
      !--- The bracket's ends in order, a < b.
      br%b  = x
      br%fb = r
      if (br%a > br%b) then
         x_new = br%a
         br%a  = br%b
         br%b  = x_new
         x_new = br%fa
         br%fa = br%fb
         br%fb = x_new
      end if
      br%kept = 0_ik
      x  = x_best
      r  = r_best
      dr = dr_best
      newton_ok = .true.
      width = br%b - br%a
      stalled = 0_ik
      do it = 1_ik, NODE_SOLVE_MAX
         if (c /= POND .and. br%a > 0.0_wp .and. br%b > NODE_SPREAD_MAX * br%a) then
            x_new = sqrt(br%a * br%b)
         else if (stalled >= 2_ik) then
            x_new = 0.5_wp * (br%a + br%b)
         else
            x_new = regula_falsi_point(br)
            if (newton_ok .and. dr > 0.0_wp) then
               if (x - r / dr > br%a .and. x - r / dr < br%b) x_new = x - r / dr
            end if
         end if
         r_last = abs(r)
         x = x_new
         call balance_at(x, r, dr)
         if (solved) exit
         call illinois_narrow(br, x, r)
         newton_ok = abs(r) <= 0.5_wp * r_last
         !--- A bracket that has not halved over two points is halved outright at the third.
         if (br%b - br%a > 0.5_wp * width) then
            stalled = stalled + 1_ik
         else
            stalled = 0_ik
            width   = br%b - br%a
         end if
         !---~---
         !   Closed to round-off: the best the arithmetic gives, solved if within Newton's own
         ! tolerance (tol is NODE_SOLVE_SHARE of it), else the balance jumps there.
         !---~---
         if (br%b - br%a <= NODE_ROUNDOFF * epsilon(1.0_wp) * max(abs(br%a), abs(br%b))) then
            solved = abs(r_best) * NODE_SOLVE_SHARE <= tol
            exit
         end if
      end do
      ts(c) = x_best
   contains
      !--- The balance at xx, its slope, and the best point so far; solved once it is within tol.
      subroutine balance_at(xx, rr, drr)
         real(wp), intent(in)  :: xx
         real(wp), intent(out) :: rr, drr
         ts(c) = xx
         call node_balance(static, profile, soil, drive, c, ts, w0, h, cut, psi_wood, rr, drr)
         if (abs(rr) < abs(r_best)) then
            x_best  = xx
            r_best  = rr
            dr_best = drr
         end if
         solved = abs(rr) <= tol
      end subroutine balance_at
   end subroutine solve_node

   !---~---
   !   solve_block -- the unknowns of a block of adjacent eliminated core nodes, c1 to c2, from their
   ! own balances solved together, the nodes outside held: Newton on the block's tridiagonal
   ! Jacobian (node_balance's slopes to the nodes above and below; thomas_solve), each step kept
   ! inside every node's curve as Newton's own steps are (node_step_limit) and backtracked until
   ! the block's merit, the sum of its squared imbalances each over its tol, falls by Armijo's
   ! test. A block topped by the pond, which then holds no water (fast_blocks), has the pond
   ! eliminated inside it as the block is from the network: every evaluation first solves the
   ! pond's own balance (solve_node) with the node below it held, so the block's Newton runs on
   ! the nodes below, and the pond passes down what reaches it, its slopes folded into the node
   ! below's (the Schur complement). The pond's face onto a thin node below is far stiffer than
   ! that node's storage, and the pond's potential is linear in its unknown where the node's is
   ! not: kept in the Newton, the two locked heads part by the curvature of the node's curve at
   ! any step that moves them, and the merit grows at every step the line search tries. Where
   ! Newton stalls, a sweep of each node's own bracketed solve (solve_node), down the block and
   ! back up, moves it on, and Newton resumes from there (swept). A single node is solve_node's
   ! alone. solved says whether every balance is within its tol.
   !---~---
   subroutine solve_block(static, profile, soil, drive, c1, c2, ts, w0, h, tol, cut, psi_wood, solved, swept)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(water_drive_t),       intent(in)    :: drive
      integer(ik),               intent(in)    :: c1, c2
      real(wp),                  intent(inout) :: ts(:)
      real(wp),                  intent(in)    :: w0(:), h, tol(:)   ! (nsoil) each node's start water; [s]; [kg/m2/s]
      logical,                   intent(in)    :: cut(:)
      real(wp),                  intent(in)    :: psi_wood(:)
      logical,                   intent(out)   :: solved, swept
      real(wp)    :: r(c2 - c1 + 1_ik), lo(c2 - c1 + 1_ik), di(c2 - c1 + 1_ik), up(c2 - c1 + 1_ik)
      real(wp)    :: dx(c2 - c1 + 1_ik), x(c2 - c1 + 1_ik), merit, merit_try, lam
      integer(ik) :: n, it, ib, j, c, k1, m
      logical     :: taken, one
      n = c2 - c1 + 1_ik
      swept = .false.
      if (n == 1_ik) then
         call solve_node(static, profile, soil, drive, c1, ts, w0(c1), h, tol(c1), cut, psi_wood, solved)
         return
      end if
      !--- The block's Newton runs on k1 .. c2, m nodes: below the pond where the pond tops it.
      k1 = c1
      if (c1 == POND) k1 = c1 + 1_ik
      m = c2 - k1 + 1_ik
      call block_balance(r, lo, di, up, merit, solved)
      do it = 1_ik, BLOCK_NEWTON_MAX
         if (solved) return
         taken = .false.
         if (merit == merit) then
            call thomas_solve(lo(1:m), di(1:m), up(1:m), -r(1:m), dx(1:m), m)
            if (all(ieee_is_finite(dx(1:m)))) then
               lam = 1.0_wp
               do j = 1_ik, m
                  call node_step_limit(profile, drive%porous, k1 + j - 1_ik, ts(k1 + j - 1_ik), dx(j), lam)
               end do
               x = ts(c1:c2)
               do ib = 1_ik, BACKTRACK_MAX
                  ts(k1:c2) = x(k1 - c1 + 1_ik:) + lam * dx(1:m)
                  call block_balance(r, lo, di, up, merit_try, solved)
                  if (merit_try == merit_try) then
                     if (sufficient_decrease(merit_try, merit, lam)) then
                        taken = .true.
                        exit
                     end if
                  end if
                  lam = 0.5_wp * lam
               end do
               if (.not. taken) ts(c1:c2) = x
            end if
         end if
         if (taken) then
            merit = merit_try
         else
            !---~---
            !   Newton stalled: one sweep of each node's own solve, down the block and back up, and
            ! Newton again from there.
            !---~---
            swept = .true.
            do j = 1_ik, 2_ik * n
               c = merge(c1 + j - 1_ik, c2 - (j - n) + 1_ik, j <= n)
               call solve_node(static, profile, soil, drive, c, ts, w0(c), h, tol(c), cut, psi_wood, one)
            end do
            call block_balance(r, lo, di, up, merit, solved)
         end if
      end do
   contains
      !---~---
      !   The Newton nodes' balances at ts and their slopes in their own and their block neighbours'
      ! unknowns (the Jacobian's sub-, main and super-diagonal, rows 1 .. m), the merit over them,
      ! and whether every balance of the block is within its tol. A pond on top is solved first,
      ! the node below held, and its slopes folded into that node's own: the node's balance moves
      ! with its unknown as the pond's solution follows it.
      !---~---
      subroutine block_balance(rr, a, b, cc, merit_out, ok)
         real(wp), intent(out) :: rr(:), a(:), b(:), cc(:), merit_out
         logical,  intent(out) :: ok
         real(wp)    :: r_pond, dr_pond, dr_pond_down, dr_pond_up
         integer(ik) :: k
         logical     :: pond_ok
         pond_ok = .true.
         if (k1 > c1) then
            call solve_node(static, profile, soil, drive, c1, ts, w0(c1), h, tol(c1), cut, psi_wood, pond_ok)
            call node_balance(static, profile, soil, drive, c1, ts, w0(c1), h, cut, psi_wood, r_pond, dr_pond,       &
                              dr_pond_up, dr_pond_down)
            pond_ok = abs(r_pond) <= tol(c1)
         end if
         do k = 1_ik, m
            call node_balance(static, profile, soil, drive, k1 + k - 1_ik, ts, w0(k1 + k - 1_ik), h, cut, psi_wood, &
                              rr(k), b(k), a(k), cc(k))
         end do
         if (k1 > c1 .and. dr_pond > 0.0_wp) b(1) = b(1) - a(1) * dr_pond_down / dr_pond
         a(1)  = 0.0_wp
         cc(m) = 0.0_wp
         merit_out = sum((rr(1:m) / tol(k1:c2)) ** 2)
         ok = pond_ok .and. all(abs(rr(1:m)) <= tol(k1:c2))
         if (.not. all(ieee_is_finite(rr(1:m)))) merit_out = ieee_value(merit_out, ieee_quiet_nan)
      end subroutine block_balance
   end subroutine solve_block

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
      integer(ik) :: c, n_porous
      n_porous = drive%porous%n
      do c = 2_ik, size(tau_soil, kind=ik)
         if (step_soil(c) == 0.0_wp) cycle
         if (c <= 1_ik + n_porous) then
            call chop(drive%porous%curve(c - 1_ik), tau_soil(c), step_soil(c))
         else
            call chop(profile%layer_curve(c - 1_ik - n_porous), tau_soil(c), step_soil(c))
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
   !   The largest share of a Newton step every store can take: no porous layer, soil layer or
   ! tissue is driven more than STEP_MAX_SHARE of the way to empty, and none on the dry side
   ! moves its saturation or relative water content by more than STEP_MAX_DRY. A wet porous or
   ! soil layer's unknown is its potential, so only the first applies to it: it keeps a step that
   ! would carry it past the branch point and below empty inside the curve. The pond's head is
   ! not bounded.
   !---~---
   pure real(wp) function step_limit(profile, porous, tau_soil, tau_coh, step_soil, step_coh) result(lam)
      type(soil_profile_t),  intent(in) :: profile
      type(water_porous_t),  intent(in) :: porous
      real(wp),              intent(in) :: tau_soil(:), tau_coh(:,:), step_soil(:), step_coh(:,:)
      integer(ik) :: i, j, c
      lam = 1.0_wp
      do c = 2_ik, 1_ik + porous%n + profile%n_active
         call node_step_limit(profile, porous, c, tau_soil(c), step_soil(c), lam)
      end do
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         do j = 1_ik, COHORT_STORES
            lam = min(lam, fraction_to_boundary(tau_coh(j, i), step_coh(j, i), STEP_MAX_SHARE))
            call cap_dry_change(tau_coh(j, i), step_coh(j, i), lam)
         end do
      end do
   end function step_limit

   !---~---
   !   step_limit's rule for core node c, a porous or soil layer (the pond's head is not bounded):
   ! lam, the share of its step d taken from its unknown tau, held to what keeps it inside its curve.
   !---~---
   pure subroutine node_step_limit(profile, porous, c, tau, d, lam)
      type(soil_profile_t), intent(in)    :: profile
      type(water_porous_t), intent(in)    :: porous
      integer(ik),          intent(in)    :: c
      real(wp),             intent(in)    :: tau, d
      real(wp),             intent(inout) :: lam
      real(wp) :: se_switch
      if (c == POND) return
      if (c <= 1_ik + porous%n) then
         se_switch = porous%curve(c - 1_ik)%se_switch
      else
         se_switch = profile%layer_curve(c - 1_ik - porous%n)%se_switch
      end if
      lam = min(lam, fraction_to_boundary(tau - SOIL_TAU_FLOOR, d, STEP_MAX_SHARE))
      if (tau <= se_switch) call cap_dry_change(tau, d, lam)
   end subroutine node_step_limit

   !---~---
   !   The cap on the change applies below full saturation (an RWC or Se of one), where the
   ! potential steepens towards empty; above it the curve is gentle and Newton moves freely.
   !---~---
   pure subroutine cap_dry_change(tau, d, lam)
      real(wp), intent(in)    :: tau, d
      real(wp), intent(inout) :: lam
      if (tau < 1.0_wp .and. abs(d) > STEP_MAX_DRY) lam = min(lam, STEP_MAX_DRY / abs(d))
   end subroutine cap_dry_change

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
   ! its capacity, and add what each path moved to the step's amounts. A store whose own balance
   ! was solved apart -- a block of eliminated core nodes, a cohort's fast leaves on their one path
   ! (classify_fast_leaves) -- takes its water from its own curve, and the faces from the block to
   ! its neighbour across its stronger link (the wood, for the leaves) carry the flows that close
   ! each member's balance; the neighbour absorbs the block's slack (docs/science/soil_plant_water.md
   ! §3, "Conservation"). Water only: the heat it carries is the energy solves'.
   !---~---
   subroutine water_commit_substep(drive, h, tau_soil, tau_coh, ev, state, flows, fast, cut, fast_leaf, spilled)
      type(water_drive_t), intent(in)    :: drive
      real(wp),            intent(in)    :: h, tau_soil(:), tau_coh(:,:)
      type(water_tend_t),  intent(in)    :: ev
      type(water_state_t), intent(inout) :: state
      type(water_flows_t), intent(inout) :: flows
      logical,             intent(in)    :: fast(:), cut(:)   ! the eliminated and the detached core nodes
      logical,             intent(in)    :: fast_leaf(:)      ! the cohorts whose leaves are fast
      logical,             intent(out)   :: spilled   ! the pond spilled, so its unknown moved
      real(wp) :: spill, sap(size(tau_coh, 2))
      real(wp) :: qf(size(tau_soil)), dinfl(size(tau_soil)), slack(size(tau_soil)), total, carried
      integer(ik) :: i, n_porous, nsoil, c, c1, c2, g
      logical     :: snap(size(tau_soil))   ! the eliminated nodes taken to their curve's water

      n_porous = drive%porous%n
      nsoil    = size(tau_soil, kind=ik)
      !--- Down the chain: face c joins node c to c + 1.
      qf = ev%face
      !---~---
      !   Each block of eliminated nodes at its curves' water: the faces from it to the neighbour
      ! across its stronger link carry the flows that close each member's balance there -- a
      ! neighbour that is not detached and holds the water the block's slack would take from it;
      ! else the other neighbour. If neither can, the block keeps the water its flows gave it, its
      ! slack Newton's, as any store's.
      !---~---
      dinfl = 0.0_wp
      snap  = fast
      if (.not. fast(POND)) call drain_only_what_is_held()
      c1 = 1_ik
      do while (c1 <= nsoil)
         if (.not. fast(c1)) then
            c1 = c1 + 1_ik
            cycle
         end if
         c2 = c1
         do while (c2 < nsoil)
            if (.not. fast(c2 + 1_ik)) exit
            c2 = c2 + 1_ik
         end do
         do c = c1, c2
            slack(c) = (ev%water_soil(c) - state%water_soil(c)) / h - (ev%inflow_soil(c) + dinfl(c))
         end do
         total = sum(slack(c1:c2))
         !--- The neighbour across the stronger link (0 for none), else the other one.
         g = 0_ik
         if (c1 > 1_ik .and. c2 < nsoil) then
            g = merge(c1 - 1_ik, c2 + 1_ik, abs(ev%face_dqb(c1 - 1_ik)) >= abs(ev%face_dqa(c2)))
         else if (c1 > 1_ik) then
            g = c1 - 1_ik
         else if (c2 < nsoil) then
            g = c2 + 1_ik
         end if
         if (g > 0_ik) then
            if (.not. can_give(g, total)) then
               g = merge(c2 + 1_ik, c1 - 1_ik, g == c1 - 1_ik)
               if (g < 1_ik .or. g > nsoil) then
                  g = 0_ik
               else if (.not. can_give(g, total)) then
                  g = 0_ik
               end if
            end if
         end if
         if (g == 0_ik) then
            snap(c1:c2) = .false.
            c1 = c2 + 1_ik
            cycle
         end if
         carried = 0.0_wp
         if (g == c2 + 1_ik) then                ! down the block: face c carries what the members above it lack
            do c = c1, c2
               carried  = carried + slack(c)
               qf(c)    = qf(c) - carried
               dinfl(c) = dinfl(c) + slack(c)
            end do
         else                                   ! up the block: face c - 1 carries what the members below it lack
            do c = c2, c1, -1_ik
               carried    = carried + slack(c)
               qf(c - 1_ik) = qf(c - 1_ik) + carried
               dinfl(c)   = dinfl(c) + slack(c)
            end do
         end if
         dinfl(g) = dinfl(g) - total
         c1 = c2 + 1_ik
      end do
      if (fast(POND) .and. .not. snap(POND)) call drain_only_what_is_held()
      !--- The flow that closes each fast leaf's balance at its curve's water; a slower leaf's is its sapflow's.
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         sap(i) = ev%sapflow(i)
         if (fast_leaf(i)) sap(i) = (ev%water_coh(LEAF, i) - state%water_coh(LEAF, i)) / h + ev%transp(i)
      end do

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
         if (.not. fast_leaf(i)) cycle
         state%water_coh(LEAF, i) = ev%water_coh(LEAF, i)
         state%water_coh(WOOD, i) = state%water_coh(WOOD, i) - h * (sap(i) - ev%sapflow(i))
      end do
      state%tau_soil = tau_soil
      state%tau_coh = tau_coh

      !--- Spill: water above the pond's capacity runs off.
      spill = max(0.0_wp, state%water_soil(POND) - drive%pond_capacity)
      spilled = spill > 0.0_wp
      if (spilled) then
         flows%spill = flows%spill + spill
         state%water_soil(POND) = drive%pond_capacity
         state%tau_soil(POND) = pond_water_to_head(drive%pond_capacity)
      end if

      !--- What each path moved.
      flows%face  = flows%face + h * qf
      flows%infil = flows%infil + h * max(qf(1), 0.0_wp)
      flows%exfil = flows%exfil + h * max(-qf(1), 0.0_wp)
      flows%porous_evap(1:n_porous) = flows%porous_evap(1:n_porous) + h * ev%porous_evap(1:n_porous)
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
      ! through its slack if it is eliminated.
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

      !--- Whether neighbour gg, not detached, holds the water a slack of amount would take from it.
      logical function can_give(gg, amount)
         integer(ik), intent(in) :: gg
         real(wp),    intent(in) :: amount   ! [kg/m2/s]
         can_give = .not. cut(gg)
         if (can_give) can_give = state%water_soil(gg) + h * (ev%inflow_soil(gg) + dinfl(gg) - amount) >= 0.0_wp
      end function can_give
   end subroutine water_commit_substep

   pure subroutine clear_flows(flows)
      type(water_flows_t), intent(inout) :: flows
      flows%face = 0.0_wp
      flows%uptake = 0.0_wp
      flows%sapflow = 0.0_wp
      flows%transp = 0.0_wp
      flows%infil = 0.0_wp
      flows%exfil = 0.0_wp
      flows%drain = 0.0_wp
      flows%evap = 0.0_wp
      flows%spill = 0.0_wp
      flows%porous_evap = 0.0_wp
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
