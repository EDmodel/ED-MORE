! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_integ_backward_euler -- the fast loop's backward-Euler solves, two systems     !
! side by side. One ESDIRK2 stage (docs/science/numerical_scheme.md §3 `esdirk2`, eq. (1))     !
! solves the canopy air against the leaf, wood and ground temperatures by a small Newton   !
! on a forward-difference Jacobian: backward_euler_stage is the stage,                     !
! newton_surface_solve and jac_surface its inner solve, cas_box_commit the box commit;     !
! advance_plant_water and advance_surf_water_full move plant and film water once per step. !
! The soil-plant water network's step is implicit Euler on adaptive sub-steps              !
! (docs/science/soil_plant_water.md §3): water_step, its Newton sub-step, the Newton       !
! matrix and the sub-step's commit. The two Newton solves share only the line search and   !
! the step bound of meds_numerics.                                                         !
!==========================================================================================!
module meds_fast_integ_backward_euler
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, cp_air
   use meds_soil_heat, only : soil_heat_env_t, soil_heat_out_t, soil_energy_step_implicit, set_soil_heat_env,     &
                              soil_layer_temp
   use meds_patch_state_vars,   only : soil_energy_state_t
   use meds_fast_integ_derivs, only : surface_derivs, cas_conductances, water_derivs
   use meds_fast_integ_control, only : water_scales
   use meds_canopy_air_space, only : cas_box_step_implicit, cas_env_t
   use meds_fast_integ_types, only : march_state_t, frozen_t, surface_state_t, cas_boundary_t,                  &
                               surface_tend_t, stage_bflux_t, error_control_t, tol_set_t, fast_static_state_t,    &
                               water_drive_t, water_tend_t, water_state_t, water_flows_t, water_work_t,           &
                               GRP_THETA, POND, WOOD, LEAF, COHORT_STORES, WORK_WATER_NSUB, WORK_WATER_NEWTON,    &
                               WORK_WATER_REJECT, WORK_WATER_NFAIL, state_init
   use meds_site_state_vars, only : site_surface_t
   use meds_therm_lib, only : cas_enthalpy_to_temp, internal_energy_liquid, liquid_energy_to_temp
   use meds_numerics,     only : sufficient_decrease, fraction_to_boundary
   use meds_soil_hydr_lib, only : pond_water_to_head
   use meds_soil_water,    only : ground_evaporation_setup
   use meds_linear_algebra, only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   use meds_soil_params, only : soil_profile_t, soil_params_t
   use meds_pft_params, only : hydraulics_table_t
   implicit none
   private

   public :: backward_euler_stage, advance_plant_water, advance_surf_water_full
   public :: water_step
   public :: water_jacobian_error   ! used only by the tests

   !----- The water network's Newton. The balance counts as solved when every store's imbalance    !
   !      over the sub-step is below NEWTON_TOL of the step's absolute water tolerance -- far inside !
   !      the error the sub-step itself is allowed. A Newton step may take a soil layer or a tissue  !
   !      at most STEP_MAX_SHARE of the way to empty (fraction to the boundary) and move its         !
   !      saturation or relative water content by at most STEP_MAX_DRY. ------------------------!
   real(wp),    parameter :: NEWTON_TOL     = 1.0e-2_wp
   integer(ik), parameter :: NEWTON_MAX     = 25_ik
   integer(ik), parameter :: BACKTRACK_MAX  = 12_ik
   real(wp),    parameter :: STEP_MAX_SHARE = 0.9_wp
   real(wp),    parameter :: STEP_MAX_DRY   = 0.2_wp
   !----- Its sub-steps: the error controller of an order-one method (exponent -1/2), its growth and !
   !      shrink limits, the factor a sub-step whose Newton failed is cut by, and the shortest      !
   !      sub-step, as a share of the fast step, below which a sub-step is taken as it stands. ----!
   real(wp),    parameter :: SAFETY = 0.9_wp, GROW_MAX = 4.0_wp, SHRINK_MIN = 0.2_wp
   real(wp),    parameter :: NEWTON_FAIL_CUT = 0.25_wp
   real(wp),    parameter :: H_MIN_SHARE = 1.0e-5_wp
   integer(ik), parameter :: TRIALS_MAX = 2000_ik
   !----- [kg/m2] Below this a pond holds no water to carry a temperature. ------------------------!
   real(wp),    parameter :: POND_TINY = 1.0e-9_wp

contains

   !---------------------------------------------------------------------------------------!
   ! backward_euler_stage -- one backward-Euler stage Y = base + dt*f(Y) of the stiff      !
   ! block (each ESDIRK2 stage, dt = gamma*h, is one call): the CAS twins, implicit in the !
   ! atmospheric exchange through the coupled 2x2 Newton, and soil heat (implicit,         !
   ! Thomas). Soil water moves at the water network's steady rate; plant water is passed   !
   ! through and advanced once (advance_plant_water).                                      !
   !---------------------------------------------------------------------------------------!
   subroutine backward_euler_stage(y, frozen, site, ec, n, nsl, dt, y_out, surf_tend, bf)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      type(error_control_t), intent(in)  :: ec     !< its budget stop ([fast].debug_error)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      !----- The stage's state and its own surface tendencies, the ones its CAS source and bf used; the   !
      !      caller b-weights transp_c into the plant water, so the leaves' water debit and the canopy     !
      !      air's vapour credit agree (docs/science/numerical_scheme.md §5a). Both filled in place. -----!
      type(march_state_t),  intent(inout) :: y_out
      type(surface_tend_t),  intent(inout) :: surf_tend
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux rates for the ESDIRK2 ledger

      type(surface_state_t)      :: y_stage
      type(cas_boundary_t)       :: cas_stage   !< the CAS boundary with this stage's conductances
      type(soil_energy_state_t) :: se
      type(soil_heat_env_t)     :: eforc
      type(soil_heat_out_t)        :: eflux
      real(wp)    :: t_ground, fliq1, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp)    :: enth1, shv1, e_infil, e_drain, t_cas1, water_heat(nsl)
      real(wp)    :: co21
      integer(ik) :: nfeval
      logical     :: ok

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      call soil_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                &
                           site%soil_profile%layer_curve(1), t_ground, fliq1)

      cas_mass_capacity = frozen%cas%box%air_mass_capacity
      cas_molar_capacity = frozen%cas%box%air_molar_capacity
      cas_stage = frozen%cas              ! plain scalars: a cheap copy, overridden with this stage's conductances below

      !----- Re-solve the surface layer at this stage's canopy-air state, so the stage pays for the       !
      !      ventilation its own temperature earns (docs/science/numerical_scheme.md §2); the box commit,  !
      !      the Newton and the ledger all read these values. cas_stage carries them with mo_live cleared,  !
      !      so surface_derivs does not re-solve on each Newton residual evaluation. ---------------------!
      call cas_conductances(frozen%cas, y%cas_enthalpy, y%cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      cas_stage%box%atm_conductance_enthalpy = g_atm_heat
      cas_stage%box%atm_conductance_vapor = g_atm_vapour
      cas_stage%box%atm_conductance_co2 = g_atm_co2
      cas_stage%mo_live = .false.

      !----- The coupled 2x2 Newton for the leaves and the canopy air; its final surf_tend drives the box !
      !      commit and the soil sinks, one flux per interface. --------------------------------------!
      call newton_surface_solve(y, cas_stage, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, &
                                surf_tend, nfeval, ok)
      !----- One box commit for all three twins, in flux form, so the budgets close for any surf_tend: !
      !      CO2 rides the same box as enthalpy and humidity. -----------------------------------------!
      call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, cas_stage, dt, enth1, shv1, co21)
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      y_out%cas_co2      = co21

      !----- Soil heat: implicit backward Euler, Thomas solve. -------------------------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      !----- The enthalpy of the water the network moved: across the faces at its time-mean flows, in     !
      !      from the pond at t_infil, out with root uptake and exfiltration at each layer's start-of-step  !
      !      temperature (what the wood and the pond were credited), and out with drainage                !
      !      (docs/science/soil_plant_water.md §5). The leaf pays its transpiration's full vapour         !
      !      enthalpy in surface_derivs. ---------------------------------------------------------------!
      call set_soil_heat_env(nsl, surf_tend%g_top, y%theta, frozen%roots%root_heat, frozen%hydrology%exfil_enth,  &
                             frozen%hydrology%w_flux_frozen, frozen%hydrology%infiltration,                       &
                             frozen%hydrology%t_infil, frozen%hydrology%drainage, frozen%hydrology%t_bot, eforc,  &
                             water_heat, e_drain)
      e_infil = frozen%hydrology%infiltration * internal_energy_liquid(frozen%hydrology%t_infil)
      call soil_energy_step_implicit(se, eforc, site%soil_profile, site%soil, dt, eflux)
      !----- [fast].debug_error stops the run when the soil heat's own budget misses by more than  !
      !      soil_energy_budget_atol [K] of the column's heat capacity. -----------------------------!
      if (ec%debug_error .and. abs(eflux%energy_resid) > ec%soil_energy_budget_atol * eflux%heat_capacity)    &
         error stop 'soil heat step: energy budget did not close'
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- Soil water is the network's, its end-of-step theta committed once by update_dynamic_state   !
      !      (docs/science/numerical_scheme.md §5b). The stages move theta at its steady rate, so each      !
      !      layer holds the water whose enthalpy the faces brought in; a steady rate is integrated exactly !
      !      and adds nothing to the error estimate. ---------------------------------------------------!
      y_out%theta(1:nsl) = y%theta(1:nsl) + dt * frozen%hydrology%theta_rate(1:nsl)

      !----- The pond is passed through: update_dynamic_state commits the network's end-of-step pond, !
      !      and the pond is outside the error norm. --------------------------------------------------!
      y_out%w_surface      = y%w_surface
      y_out%w_surface_energy = y%w_surface_energy

      !----- Plant water mass is passed through (advance_plant_water advances it). -----------!
      y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)

      !----- This stage's boundary-flux rates for the ledger, all at the committed state, so the amounts   !
      !      esdirk2_step and adaptive_march b-weight and accumulate telescope to closure. -----------!
      if (present(bf)) then
         t_cas1 = cas_enthalpy_to_temp(enth1, shv1)         ! committed CAS temperature, for the reported H
            bf%cas_enth_in  = surf_tend%src_enth + g_atm_heat*frozen%cas%box%atm_enthalpy
            bf%cas_enth_out = g_atm_heat*enth1
            bf%cas_vap_in   = surf_tend%src_vap  + g_atm_vapour*frozen%cas%box%atm_specific_humidity
            bf%cas_vap_out  = g_atm_vapour*shv1
            bf%cas_co2_in   = frozen%cas%nee_biotic + g_atm_co2*frozen%cas%box%atm_co2
            bf%cas_co2_out  = g_atm_co2*y_out%cas_co2
            !----- The bottom face the soil solve committed (zero on the adiabatic base, the conduction  !
            !      to deep_temp under the Dirichlet anchor) crosses both the soil's and the column's book. !
            bf%soil_enth_in = surf_tend%g_top - eflux%bottom_heat + e_infil
            bf%soil_enth_out= sum(water_heat(1:nsl)) + e_drain
            !----- Soil water's storage change and flows are the network's, booked once per step by     !
            !      close_ledgers: the per-stage soil-water amounts are zero, and the whole-column          !
            !      water here is only the CAS-vapour exchange. ------------------------------------------!
            bf%soil_wat_in  = 0.0_wp
            bf%soil_wat_out = 0.0_wp
            !----- ground_rad is the tiles' radiation, cover-weighted (rnet_soil when bare). Infiltration,  !
            !      exfiltration and uptake move enthalpy between tracked stores and cancel; rain and runoff  !
            !      are added once at the outer level (docs/science/numerical_scheme.md §5b). -------------!
            bf%whole_enth_in= surf_tend%coh_rnet + surf_tend%ground_rad - eflux%bottom_heat
            !----- The condensate is not a boundary loss: the caller deposits it in soil layer 1 with the  !
            !      enthalpy surface_derivs debited from the CAS (whole_cond, whole_cond_enth). -----------!
            bf%whole_enth_out= g_atm_heat*(enth1 - frozen%cas%box%atm_enthalpy) + e_drain
            bf%whole_wat_in = 0.0_wp
            bf%whole_wat_out = g_atm_vapour*(shv1 - frozen%cas%box%atm_specific_humidity)
            bf%whole_cond   = surf_tend%cond                     ! deposited into soil layer 1, not lost
            bf%whole_cond_enth = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
            bf%snow_enth_in = surf_tend%snow_enth_in
            bf%snow_vap_out = surf_tend%snow_vap_out
            !----- The reported H, against theta_atm (T at zref), the ledger's basis. ----------------------!
            bf%atm_heat_out = g_atm_heat*cp_air*(t_cas1 - frozen%cas%env%theta_atm)
            bf%atm_vap_out  = g_atm_vapour*(shv1  - frozen%cas%box%atm_specific_humidity)
            bf%rad_net      = surf_tend%coh_rnet + surf_tend%ground_rad
      end if
   end subroutine backward_euler_stage

   !---------------------------------------------------------------------------------------!
   ! newton_surface_solve -- a direct 2x2 Newton solve of the coupled backward-Euler canopy-air block   !
   ! {R_H, R_q} = 0 in (cas_enthalpy, cas_shv), whose surface sources depend on both through tcas,    !
   ! qcas and qsat; robust near saturation, where a leaf <-> CAS Picard iteration is not. The forward- !
   ! difference Jacobian (jac_surface) captures the strong VPD self-limiting d src_vap/dq. A singular- !
   ! Jacobian fallback, a line search and an evaluation cap; it never stops the run. Returns surf_tend !
   ! at the solution, from which backward_euler_stage commits the box.                                 !
   !---------------------------------------------------------------------------------------!
   subroutine newton_surface_solve(y, cas, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, &
                                   surf_tend, nfeval, ok)
      type(march_state_t),   intent(in)    :: y
      type(cas_boundary_t),   intent(in)    :: cas      !< the stage's CAS boundary (its own conductances)
      type(frozen_t),  intent(in)    :: frozen   !< tissue / film / ground / snow for surface_derivs
      real(wp),               intent(in)    :: t_ground
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, cas_mass_capacity, g_atm_heat, g_atm_vapour
      type(surface_tend_t),   intent(inout) :: surf_tend   !< the evaluation at the solution, on return
      integer(ik),            intent(out)   :: nfeval
      logical,                intent(out)   :: ok

      type(surface_state_t) :: y_stage
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt
      integer(ik) :: it, ls
      real(wp),    parameter :: RTOL_N = 1.0e-7_wp, ATOL_H = 5.0e1_wp, ATOL_Q = 1.0e-6_wp
      real(wp),    parameter :: DETEPS = 1.0e-30_wp
      integer(ik), parameter :: NEWT_MAX = 4_ik, LS_MAX = 6_ik, FEVAL_CAP = 24_ik

      H0 = y%cas_enthalpy
      q0 = y%cas_shv
      Hk = H0
      qk = q0
      nfeval = 0_ik
      ok = .false.
      y_stage%cas_co2 = y%cas_co2
      y_stage%cas_enthalpy = Hk
      y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground,                            &
         t_ground, n, surf_tend)
      nfeval = nfeval + 1_ik
      R_H = cas_mass_capacity*(Hk - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%box%atm_enthalpy - Hk)
      R_q = cas_mass_capacity*(qk - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%box%atm_specific_humidity  - qk)

      do it = 1_ik, NEWT_MAX
         if ( abs(R_H)*dt/cas_mass_capacity <= ATOL_H + RTOL_N*abs(Hk) .and.                                  &
              abs(R_q)*dt/cas_mass_capacity <= ATOL_Q + RTOL_N*abs(qk) ) then
            ok = .true.
            exit
         end if
         call jac_surface(Hk, qk, y%cas_co2, cas, frozen, t_ground, surf_tend, n, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt, &
                          J11, J12, J21, J22, nfeval)
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
         do ls = 1_ik, LS_MAX
            Ht = Hk + lam*delH
            qt = qk + lam*delq
            !----- No supersaturation clamp on the state: a jump there inflates the embedded error near  !
            !      RH = 1. Transient supersaturation is relaxed by surface_derivs' condensation sink. ----!
            y_stage%cas_enthalpy = Ht
            y_stage%cas_shv = qt
            call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground,                      &
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
         if (nfeval >= FEVAL_CAP) exit
      end do
      !----- surf_tend holds the last evaluation, which is at (Hk, qk): the start when the first test  !
      !      passes, else the accepted line-search trial. The caller commits the box from it. --------!
   end subroutine newton_surface_solve

   !---------------------------------------------------------------------------------------!
   ! The backward-Euler canopy-air box commit, through the shared box kernel                !
   ! cas_box_step_implicit (meds_canopy_air_space; docs/science/canopy_air_space_biophysics.md !
   ! eq. (2)): one implementation of the box, not a copy of its formulas here.              !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_box_commit(h0, q0, c0, surf_tend, cas, dt, h1, q1, c1)
      real(wp),               intent(in)  :: h0, q0, c0, dt
      type(surface_tend_t),   intent(in)  :: surf_tend
      type(cas_boundary_t),   intent(in)  :: cas        !< the stage's boundary: its box carries the stage's conductances
      real(wp),               intent(out) :: h1, q1, c1
      type(cas_env_t) :: env
      env = cas%box
      env%surface_enthalpy_source = surf_tend%src_enth
      env%surface_vapor_source    = surf_tend%src_vap
      env%biotic_co2_source       = cas%nee_biotic
      call cas_box_step_implicit(h0, q0, c0, env, dt, h1, q1, c1)
   end subroutine cas_box_commit

   !----- 2x2 numerical Jacobian of (R_H, R_q) w.r.t. (H, q) by forward-differencing surface_derivs. --!
   subroutine jac_surface(Hk, qk, co2, cas, frozen, t_ground, surf_tend, n, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt, J11, &
                          J12, J21, J22, &
                          nfeval)
      real(wp),               intent(in)    :: Hk, qk, co2, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt
      type(cas_boundary_t),   intent(in)    :: cas
      type(frozen_t),  intent(in)    :: frozen
      real(wp),               intent(in)    :: t_ground
      type(surface_tend_t),   intent(in)    :: surf_tend         ! base eval at (Hk,qk)
      integer(ik),            intent(in)    :: n
      real(wp),               intent(out)   :: J11, J12, J21, J22
      integer(ik),            intent(inout) :: nfeval
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq, dse_dH, dsv_dH, dse_dq, dsv_dq
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      y_stage%cas_co2 = co2
      y_stage%cas_enthalpy = Hk + dH
      y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground,                            &
         t_ground, n, sfp)
      nfeval = nfeval + 1_ik
      dse_dH = (sfp%src_enth - surf_tend%src_enth)/dH
      dsv_dH = (sfp%src_vap - surf_tend%src_vap)/dH
      y_stage%cas_enthalpy = Hk
      y_stage%cas_shv = qk + dq
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground,                            &
         t_ground, n, sfp)
      nfeval = nfeval + 1_ik
      dse_dq = (sfp%src_enth - surf_tend%src_enth)/dq
      dsv_dq = (sfp%src_vap - surf_tend%src_vap)/dq
      J11 = cas_mass_capacity/dt + g_atm_heat - dse_dH
      J12 =              - dse_dq
      J21 =              - dsv_dH
      J22 = cas_mass_capacity/dt + g_atm_vapour - dsv_dq
   end subroutine jac_surface

   !---------------------------------------------------------------------------------------!
   ! advance_plant_water -- each cohort's leaf and wood water over an ESDIRK2 step at the network's  !
   ! time-mean rates: sapflow and uptake in; sapflow and the b-weighted transpiration the stages     !
   ! charged the canopy air out. The stages hold transpiration at the network's, so a dt_fast step  !
   ! ends where the network did. Nothing is floored: the network keeps every store inside its curve. !
   !---------------------------------------------------------------------------------------!
   pure subroutine advance_plant_water(y, frozen, n, dt, transp_c_bw, y_out)
      type(march_state_t),  intent(in)    :: y
      type(frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: dt
      real(wp),              intent(in)    :: transp_c_bw(n)  !< [kg/m2 ground/s] per cohort, b-weighted over the step
      type(march_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      do i = 1_ik, n
         y_out%leaf_water_mass(i) = y%leaf_water_mass(i)                                                  &
              + dt * (frozen%plant%sapflow_frozen(i) - transp_c_bw(i) / max(frozen%plant%nplant(i), tiny_num))
         y_out%wood_water_mass(i) = y%wood_water_mass(i)                                                  &
              + dt * (frozen%plant%uptake_frozen(i) - frozen%plant%sapflow_frozen(i))
      end do
   end subroutine advance_plant_water

   !---------------------------------------------------------------------------------------!
   ! advance_surf_water_full -- the canopy film over the full dt, the film analogue of                 !
   ! advance_plant_water: frozen interception in, b-weighted film evaporation out (the number the CAS  !
   ! was credited). Unclamped: update_dynamic_state clamps it and books the overflow and deficit.     !
   !---------------------------------------------------------------------------------------!
   pure subroutine advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
      type(march_state_t),  intent(in)    :: y
      type(frozen_t), intent(in)    :: frozen
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt
      real(wp),               intent(in)    :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      type(march_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      !----- Not floored: the caller's deficit bookkeeping must see a negative film to credit the       !
      !      shortfall to the whole-column outflow; a floor here would create that water. ------------!
      do i = 1_ik, n
         y_out%leaf_surf_water(i) = y%leaf_surf_water(i) + dt*(frozen%film%intercept_leaf(i) - film_evap_leaf_bw(i))
         y_out%wood_surf_water(i) = y%wood_surf_water(i) + dt*(frozen%film%intercept_wood(i) - film_evap_wood_bw(i))
      end do
   end subroutine advance_surf_water_full

   !=======================================================================================!
   !  The soil-plant water network's step (docs/science/soil_plant_water.md §3): implicit    !
   !  Euler on adaptive sub-steps, each a Newton solve (eq. (2)) in every store's own        !
   !  unknown, with a block-arrowhead linear solve (meds_linear_algebra) and a filtered      !
   !  local-error estimate. Water moves by the flows themselves, so it is conserved to       !
   !  round-off whatever Newton's residual.                                                   !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! water_step -- the network over one fast step, from the state set_water_state set, under  !
   ! the GRP_THETA tolerances. On return the state holds the end-of-step water, unknowns and  !
   ! potentials, the flows what each path moved, and work_count the network's work: accepted  !
   ! sub-steps, Newton iterations, rejected sub-steps, and sub-steps taken although they did  !
   ! not converge or meet the tolerance (counted, not hidden). The drive is held over the     !
   ! step, so the ground evaporation's constant terms are set into it once, here.             !
   !---------------------------------------------------------------------------------------!
   subroutine water_step(dt, static, profile, soil, hyd, drive, tols, water, work_count)
      real(wp),                  intent(in)    :: dt
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(inout) :: drive
      type(tol_set_t),           intent(in)    :: tols
      type(water_work_t),        intent(inout) :: water
      integer(ik),               intent(inout) :: work_count(:)   !< the step report's work counters
      real(wp)    :: t, h, h_min, err, rtol, atol
      real(wp)    :: scale_soil(water%nsoil), scale_coh(COHORT_STORES, water%ncoh)
      real(wp)    :: tau_soil(water%nsoil), tau_coh(COHORT_STORES, water%ncoh)
      integer(ik) :: trials, iters, i_sol, nsub, nnewton, nreject, nfail
      logical     :: converged, last, spilled

      rtol = tols%rtol(GRP_THETA)
      atol = tols%atol(GRP_THETA)
      call clear_flows(water%flows)
      nsub    = 0_ik
      nnewton = 0_ik
      nreject = 0_ik
      nfail   = 0_ik
      call arrowhead_allocate(water%sys, water%nsoil, water%ncoh, COHORT_STORES)
      call water_scales(static, profile, hyd, scale_soil, scale_coh)
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
                                   tau_soil, tau_coh, water, converged, iters, i_sol)
         nnewton = nnewton + iters
         if (.not. converged .and. .not. last) then
            nreject = nreject + 1_ik
            h = max(h * NEWTON_FAIL_CUT, h_min)
            cycle
         end if
         !----- A sub-step whose start already balanced took no iteration: it ends where it began, its  !
         !      inflows unchanged, so its error is zero, and it has no factorisation of its own to      !
         !      filter with. ---------------------------------------------------------------------------!
         if (.not. converged) then
            err = huge(1.0_wp)
         else if (iters == 0_ik) then
            err = 0.0_wp
         else
            err = water_substep_error(water%sys, h, water%ev(water%i_start), water%ev(i_sol), scale_soil,       &
                                      scale_coh, rtol, atol)
         end if
         if (err > 1.0_wp .and. .not. last) then
            nreject = nreject + 1_ik
            h = max(h * max(SHRINK_MIN, SAFETY * err ** (-0.5_wp)), h_min)
            cycle
         end if
         if (err > 1.0_wp) nfail = nfail + 1_ik
         call water_commit_substep(soil, drive, h, tau_soil, tau_coh, water%ev(i_sol), water%state, water%flows,  &
                                   spilled)
         !----- The next sub-step starts where this one ended: its start is this one's solution, whose  !
         !      Jacobian sys still holds -- unless the spill moved the pond. ----------------------------!
         water%i_start = i_sol
         water%start_ready = .not. spilled
         t = t + h
         nsub = nsub + 1_ik
         if (err <= 1.0_wp)                                                                                   &
            h = h * min(GROW_MAX, max(SHRINK_MIN, SAFETY * max(err, 1.0e-10_wp) ** (-0.5_wp))) ! clamp-ok: controller floor
      end do
      water%state%psi_soil = water%ev(water%i_start)%psi_soil
      water%state%psi_coh  = water%ev(water%i_start)%psi_coh
      work_count(WORK_WATER_NSUB)   = nsub
      work_count(WORK_WATER_NEWTON) = nnewton
      work_count(WORK_WATER_REJECT) = nreject
      work_count(WORK_WATER_NFAIL)  = nfail
   end subroutine water_step

   !----- One sub-step: Newton on the end-of-sub-step water balance. ---------------------------!
   subroutine water_newton_substep(static, profile, soil, hyd, drive, h, atol, scale_soil, scale_coh,             &
                                   tau_soil, tau_coh, water, converged, iters, i_sol)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(in)    :: drive
      real(wp),                  intent(in)    :: h, atol, scale_soil(:), scale_coh(:,:)
      real(wp),                  intent(inout) :: tau_soil(:), tau_coh(:,:)   !< start of sub-step in, solution out
      type(water_work_t),        intent(inout) :: water
      logical,                   intent(out)   :: converged
      integer(ik),               intent(out)   :: iters
      integer(ik),               intent(out)   :: i_sol   !< which of water%ev holds the solution
      real(wp) :: res_soil(water%nsoil), res_coh(COHORT_STORES, water%ncoh)
      real(wp) :: step_soil(water%nsoil), step_coh(COHORT_STORES, water%ncoh)
      real(wp) :: try_soil(water%nsoil), try_coh(COHORT_STORES, water%ncoh)
      real(wp) :: res2_soil(water%nsoil), res2_coh(COHORT_STORES, water%ncoh)
      real(wp) :: merit, merit_try, worst, worst_try, lam, lam_max
      integer(ik) :: it, ib, i_try
      logical  :: ok

      converged = .false.
      iters = 0_ik
      !----- The network at the start, evaluated unless it already is; the iterates are evaluated into  !
      !      the other slot, and the solution is the start until one is taken. ------------------------!
      i_sol = water%i_start
      i_try = 3_ik - water%i_start
      if (.not. water%start_ready) then
         call water_derivs(static, profile, soil, hyd, drive, tau_soil, tau_coh, water%ev(i_sol), water%sys, .true.)
         water%start_ready = .true.
      end if
      call water_residual(water%state, h, water%ev(i_sol), atol, scale_soil, scale_coh, res_soil, res_coh,      &
                          merit, worst)
      if (merit /= merit) return
      do it = 1_ik, NEWTON_MAX
         if (worst <= 1.0_wp) exit
         iters = it
         call water_newton_factor(water%sys, h, water%ev(i_sol), ok)
         if (.not. ok) return
         call arrowhead_solve(water%sys, -res_soil, -res_coh, step_soil, step_coh)
         lam_max = step_limit(profile, tau_soil, tau_coh, step_soil, step_coh)
         lam = lam_max
         !----- A trial's evaluation replaces the Jacobian sys held, so the start's is gone. ---------!
         water%start_ready = .false.
         do ib = 1_ik, BACKTRACK_MAX
            try_soil = tau_soil + lam * step_soil
            try_coh = tau_coh + lam * step_coh
            call water_derivs(static, profile, soil, hyd, drive, try_soil, try_coh, water%ev(i_try), water%sys,   &
                              .true.)
            call water_residual(water%state, h, water%ev(i_try), atol, scale_soil, scale_coh, res2_soil,          &
                                res2_coh, merit_try, worst_try)
            if (merit_try == merit_try .and. sufficient_decrease(merit_try, merit, lam)) exit
            lam = 0.5_wp * lam
         end do
         if (merit_try /= merit_try) return
         tau_soil = try_soil
         tau_coh = try_coh
         res_soil = res2_soil
         res_coh = res2_coh
         merit = merit_try
         worst = worst_try
         i_sol = i_try
      end do
      !----- The error estimate solves with the last iteration's factorisation, one Newton step       !
      !      behind the solution, not a fresh one at it: its filter needs the stiff exchanges' scale, !
      !      not the exact matrix, and a factorisation per sub-step is saved. ------------------------!
      converged = worst <= 1.0_wp
   end subroutine water_newton_substep

   !----- One store's water-balance residual [kg/m2/s]: (W(tau) - W_start)/h - inflow. ------------!
   elemental real(wp) function store_residual(water, water_start, h, inflow) result(res)
      real(wp), intent(in) :: water, water_start, h, inflow
      res = (water - water_start) / h - inflow
   end function store_residual

   !----- Every store's residual and, in the same pass, two measures of it, each store's imbalance !
   !      over the sub-step against NEWTON_TOL of its absolute water tolerance: worst, below one,  !
   !      means the balance is solved; merit, the sum of squares, is what the line search asks to  !
   !      fall -- a Newton step always points down it, which is not true of the worst store alone. !
   !      The soil's squares and the cohorts' are summed apart, then added. -----------------------!
   pure subroutine water_residual(state, h, ev, atol, scale_soil, scale_coh, res_soil, res_coh, merit, worst)
      type(water_state_t), intent(in)  :: state
      real(wp),            intent(in)  :: h, atol, scale_soil(:), scale_coh(:,:)
      type(water_tend_t),  intent(in)  :: ev
      real(wp),            intent(out) :: res_soil(:), res_coh(:,:), merit, worst
      real(wp)    :: tol, x, merit_coh
      integer(ik) :: c, i, j
      tol = NEWTON_TOL * atol
      merit = 0.0_wp
      worst = 0.0_wp
      do c = 1_ik, size(res_soil, kind=ik)
         res_soil(c) = store_residual(ev%water_soil(c), state%water_soil(c), h, ev%inflow_soil(c))
         x = res_soil(c) * h / (tol * scale_soil(c))
         merit = merit + x * x
         worst = max(worst, abs(x))
      end do
      merit_coh = 0.0_wp
      do i = 1_ik, size(res_coh, 2, kind=ik)
         do j = 1_ik, COHORT_STORES
            res_coh(j, i) = store_residual(ev%water_coh(j, i), state%water_coh(j, i), h, ev%inflow_coh(j, i))
            x = res_coh(j, i) * h / (tol * scale_coh(j, i))
            merit_coh = merit_coh + x * x
            worst = max(worst, abs(x))
         end do
      end do
      merit = merit + merit_coh
   end subroutine water_residual

   !----- The largest share of a Newton step every store can take: no soil layer or tissue is      !
   !      driven more than STEP_MAX_SHARE of the way to empty, and none on the dry side moves its      !
   !      saturation or relative water content by more than STEP_MAX_DRY. A wet soil layer's unknown   !
   !      is its potential, so only the first applies to it: it keeps a step that would carry it past  !
   !      the branch point and below empty inside the curve. The pond's head is not bounded. -------!
   pure real(wp) function step_limit(profile, tau_soil, tau_coh, step_soil, step_coh) result(lam)
      type(soil_profile_t), intent(in) :: profile
      real(wp),             intent(in) :: tau_soil(:), tau_coh(:,:), step_soil(:), step_coh(:,:)
      integer(ik) :: k, i, j
      lam = 1.0_wp
      do k = 1_ik, profile%n_active
         lam = min(lam, fraction_to_boundary(tau_soil(1_ik + k), step_soil(1_ik + k), STEP_MAX_SHARE))
         if (tau_soil(1_ik + k) <= profile%layer_curve(k)%se_switch)                                           &
            call cap_dry_change(tau_soil(1_ik + k), step_soil(1_ik + k), lam)
      end do
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         do j = 1_ik, COHORT_STORES
            lam = min(lam, fraction_to_boundary(tau_coh(j, i), step_coh(j, i), STEP_MAX_SHARE))
            call cap_dry_change(tau_coh(j, i), step_coh(j, i), lam)
         end do
      end do
   contains
      !----- The cap on the change applies below full saturation (an RWC or Se of one), where the     !
      !      potential steepens towards empty; above it the curve is gentle and Newton moves freely. ---!
      pure subroutine cap_dry_change(tau, d, lam)
         real(wp), intent(in)    :: tau, d
         real(wp), intent(inout) :: lam
         if (tau < 1.0_wp .and. abs(d) > STEP_MAX_DRY) lam = min(lam, STEP_MAX_DRY / abs(d))
      end subroutine cap_dry_change
   end function step_limit

   !----- The storage term dW/dtau / h each store adds to the Newton matrix's diagonal. -----------!
   pure subroutine water_storage_diagonal(h, ev, add_soil, add_coh)
      real(wp),           intent(in)  :: h
      type(water_tend_t), intent(in)  :: ev
      real(wp),           intent(out) :: add_soil(:), add_coh(:,:)
      add_soil = ev%dwater_soil / h
      add_coh  = ev%dwater_coh / h
   end subroutine water_storage_diagonal

   !----- Factor the Newton matrix: the flows' Jacobian of the latest evaluation (sys) plus the     !
   !      storage terms of ev. ------------------------------------------------------------------!
   pure subroutine water_newton_factor(sys, h, ev, ok)
      type(arrowhead_system_t), intent(inout) :: sys
      real(wp),                 intent(in)    :: h
      type(water_tend_t),       intent(in)    :: ev
      logical,                  intent(out)   :: ok
      real(wp) :: add_soil(size(ev%dwater_soil)), add_coh(COHORT_STORES, size(ev%dwater_coh, 2))
      call water_storage_diagonal(h, ev, add_soil, add_coh)
      call arrowhead_factor(sys, add_soil, add_coh, ok)
   end subroutine water_newton_factor

   !----- The sub-step's local error relative to its tolerance (<= 1 to accept): implicit Euler's  !
   !      h (F_end - F_start)/2, filtered twice through (I - h dF/dW)^-1 so a stiff exchange that  !
   !      has already settled does not read as error (docs/science/soil_plant_water.md §3). Each   !
   !      store is weighed against atol of its scale plus rtol of its water. The Newton matrix is  !
   !      the last iteration's factorisation, with the coupling rows of the solution the last      !
   !      evaluation wrote: close to the matrix at the end, which is all a filter needs. ----------!
   pure real(wp) function water_substep_error(sys, h, ev0, ev1, scale_soil, scale_coh, rtol, atol) result(err)
      type(arrowhead_system_t), intent(in) :: sys        !< the Newton matrix, factored at the last iteration
      real(wp),                 intent(in) :: h
      type(water_tend_t),       intent(in) :: ev0, ev1   !< the network at the sub-step's start and end
      real(wp),                 intent(in) :: scale_soil(:), scale_coh(:,:), rtol, atol
      real(wp) :: v_soil(size(ev1%water_soil)), v_coh(COHORT_STORES, size(ev1%water_coh, 2))
      real(wp) :: y_soil(size(ev1%water_soil)), y_coh(COHORT_STORES, size(ev1%water_coh, 2))
      integer(ik) :: pass
      v_soil = h * (ev1%inflow_soil - ev0%inflow_soil)
      v_coh  = h * (ev1%inflow_coh - ev0%inflow_coh)
      do pass = 1, 2
         call arrowhead_solve(sys, v_soil / h, v_coh / h, y_soil, y_coh)
         v_soil = ev1%dwater_soil * y_soil
         v_coh  = ev1%dwater_coh * y_coh
      end do
      err = maxval(0.5_wp * abs(v_soil) / (atol * scale_soil + rtol * abs(ev1%water_soil)))
      if (size(v_coh, 2) > 0)                                                                                 &
         err = max(err, maxval(0.5_wp * abs(v_coh) / (atol * scale_coh + rtol * abs(ev1%water_coh))))
      if (err /= err) err = huge(1.0_wp)
   end function water_substep_error

   !----- Accept a sub-step: advance each store's water by the flows themselves, spill the pond above !
   !      its capacity, carry the pond's enthalpy, and add what each path moved to the step's amounts.  !
   !      A store on one path (the pond, a cohort's leaves) takes its water from its own curve, and its !
   !      path carries the flow that closes its balance; the store across the path absorbs Newton's    !
   !      slack (docs/science/soil_plant_water.md §3, "Conservation"). ----------------------------!
   subroutine water_commit_substep(soil, drive, h, tau_soil, tau_coh, ev, state, flows, spilled)
      type(soil_params_t), intent(in)    :: soil
      type(water_drive_t), intent(in)    :: drive
      real(wp),            intent(in)    :: h, tau_soil(:), tau_coh(:,:)
      type(water_tend_t),  intent(in)    :: ev
      type(water_state_t), intent(inout) :: state
      type(water_flows_t), intent(inout) :: flows
      logical,             intent(out)   :: spilled   !< the pond spilled, so its unknown moved
      real(wp) :: w_mix, e_mix, t_mix, u_mix, amount, spill, t_pond, pond_flux, sap(size(tau_coh, 2))
      integer(ik) :: i, nsl

      nsl = size(flows%face, kind=ik)
      !----- The flows that close the pond's and each leaf's balance at their curves' water. -------!
      pond_flux = (state%water_soil(POND) + h * drive%rain - ev%water_soil(POND)) / h
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         sap(i) = (ev%water_coh(LEAF, i) - state%water_coh(LEAF, i)) / h + ev%transp(i)
      end do

      !----- The pond's enthalpy, in the order its water moves: the rain joins it first, so what    !
      !      infiltrates leaves at the mixture's temperature; water the top layer pushes up arrives   !
      !      at that layer's temperature; then the spill leaves at the pond's own. --------------!
      w_mix = state%water_soil(POND) + h * drive%rain
      e_mix = state%e_pond + h * drive%rain * internal_energy_liquid(drive%t_rain)
      flows%rain_enth = flows%rain_enth + h * drive%rain * internal_energy_liquid(drive%t_rain)
      t_mix = drive%t_rain
      if (w_mix > POND_TINY) t_mix = liquid_energy_to_temp(e_mix / w_mix)
      u_mix = internal_energy_liquid(t_mix)
      if (pond_flux >= 0.0_wp) then
         amount = h * pond_flux
         flows%infil = flows%infil + amount
         flows%infil_enth = flows%infil_enth + amount * u_mix
         state%e_pond = e_mix - amount * u_mix
      else
         amount = -h * pond_flux
         flows%exfil = flows%exfil + amount
         flows%exfil_enth = flows%exfil_enth + amount * internal_energy_liquid(drive%t_soil(1))
         state%e_pond = e_mix + amount * internal_energy_liquid(drive%t_soil(1))
      end if

      !----- Water, by the flows; the pond and the leaves at their curves, the top layer and the woods !
      !      across the paths that closed them. -----------------------------------------------!
      state%water_soil = state%water_soil + h * ev%inflow_soil
      state%water_coh  = state%water_coh + h * ev%inflow_coh
      state%water_soil(POND) = ev%water_soil(POND)
      state%water_soil(2) = state%water_soil(2) + h * (pond_flux - ev%pond_flux)
      do i = 1_ik, size(tau_coh, 2, kind=ik)
         state%water_coh(LEAF, i) = ev%water_coh(LEAF, i)
         state%water_coh(WOOD, i) = state%water_coh(WOOD, i) - h * (sap(i) - ev%sapflow(i))
      end do
      state%tau_soil = tau_soil
      state%tau_coh = tau_coh

      !----- Spill: water above the pond's capacity runs off at the pond's temperature. --------!
      spill = max(0.0_wp, state%water_soil(POND) - soil%w_pond_max)
      spilled = spill > 0.0_wp
      if (spilled) then
         t_pond = liquid_energy_to_temp(state%e_pond / state%water_soil(POND))
         flows%spill = flows%spill + spill
         flows%spill_enth = flows%spill_enth + spill * internal_energy_liquid(t_pond)
         state%e_pond = state%e_pond - spill * internal_energy_liquid(t_pond)
         state%water_soil(POND) = soil%w_pond_max
         state%tau_soil(POND) = pond_water_to_head(soil%w_pond_max)
      end if
      if (state%water_soil(POND) <= POND_TINY) then
         state%e_pond = 0.0_wp
      end if

      !----- What each path moved. --------------------------------------------------------!
      flows%face(1:nsl-1) = flows%face(1:nsl-1) + h * ev%face_flux(1:nsl-1)
      flows%drain   = flows%drain   + h * ev%drain
      flows%evap    = flows%evap    + h * ev%evap
      flows%uptake  = flows%uptake  + h * ev%uptake
      flows%sapflow = flows%sapflow + h * sap
      flows%transp  = flows%transp  + h * ev%transp
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
   end subroutine clear_flows

   !----- For the tests: how far the assembled Jacobian is from the residual's own derivative, at  !
   !      the current state for a sub-step of h -- the largest difference between J v and a        !
   !      central difference of the residual along v, relative to the largest entry of either. The !
   !      drive is set up as water_step sets it. --------------------------------------------------!
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

end module meds_fast_integ_backward_euler
