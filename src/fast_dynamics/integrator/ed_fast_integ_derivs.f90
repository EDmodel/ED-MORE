! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_derivs -- the right-hand sides the fast loop's schemes evaluate, free of
! side effects: no flux hidden in a backward-Euler denominator and no store mutated, so
! any tableau can combine the evaluations (docs/science/numerical_scheme.md §3).
!   * surface_derivs -- leaf and wood energy balance, ground skin, the three canopy-air
!     twins;
!   * march_derivs -- surface_derivs plus the soil-heat column and the soil and plant
!     water, which move at the water network's frozen rates (§2);
!   * cas_conductances -- the canopy-air <-> atmosphere conductances at a live canopy-air
!     state;
!   * water_derivs -- the soil-plant water network's flows, inflows and Jacobian at one
!     state, and set_water_state its state from the model's
!     (docs/science/soil_plant_water.md §1-§2).
!---~---
module ed_fast_integ_derivs
   use ed_kinds, only : wp, ik
   use ed_constants,   only : tiny_num, rho_h2o, grav_head
   use ed_therm_lib,   only : cas_enthalpy_to_temp, sat_specific_humidity, sat_specific_humidity_temp_deriv
   use ed_soil_heat,   only : soil_heat_env_t, soil_energy_time_deriv, set_soil_heat_env, soil_layer_temp
   use ed_surface, only : surface_env_t, ground_surface_exchange
   use ed_patch_state_vars,   only : soil_energy_state_t
   use ed_canopy_air_space, only : cas_env_t, cas_box_time_deriv, cas_condensation
   use ed_canopy_aerodynamics, only : canopy_air_conductances
   use ed_plant_energy,   only : plant_energy_env_t, canopy_tissue_energy
   use ed_canopy_water,   only : canopy_water_out_t
   use ed_fast_integ_types, only : surface_state_t, surface_tend_t, cas_boundary_t,                                  &
                               march_state_t, frozen_t, march_tend_t, fast_static_state_t,                    &
                               water_drive_t, water_tend_t, water_work_t, ensure_water_capacity, POND, WOOD, LEAF, &
                              water_organic_t
   use ed_site_state_vars, only : site_surface_t
   use ed_soil_hydr_lib,   only : soil_water_relations, soil_theta_to_unknown, pond_water_relations,              &
                                  pond_water_to_head, soil_ice_t
   use ed_plant_hydr_lib,   only : tissue_water_relations, tissue_water_to_rwc, tissue_rwc_to_unknown,            &
                                   tissue_psi_to_rwc
   use ed_plant_hydraulics,   only : supply_limited_transpiration, root_water_uptake_profile, xylem_sapflow,        &
                                     root_water_uptake
   use ed_soil_water, only : soil_water_flux, ground_evaporation_rate, bottom_flux, evaporation_supply
   use ed_linear_algebra, only : arrowhead_system_t
   use ed_soil_params, only : soil_profile_t, soil_params_t
   use ed_pft_params, only : hydraulics_table_t
   implicit none
   private

   public :: surface_derivs, march_derivs, cas_conductances
   public :: set_water_state, water_derivs, node_balance, core_node_relations, wood_potentials

   !---~---
   !   [-] The drive (head gradient) over which a wetting face's conductivity blends from the
   ! downstream node's to the upstream one's, through their arithmetic mean at no flow: a
   ! centimetre of head per metre, far below the gradients of a wetting or drying front and well
   ! above round-off. A numerical choice.
   !---~---
   real(wp), parameter :: WET_BLEND_WIDTH = 1.0e-2_wp

contains

   !---~---
   !   refresh_cas_conductances -- re-solve only the bulk Monin-Obukhov surface layer at a live canopy-
   ! air state for the three CAS <-> atmosphere conductances (docs/science/canopy_aerodynamics.md
   ! eq. (6)). The per-cohort boundary layers stay frozen (docs/science/numerical_scheme.md §2), as do
   ! the geometry and the reference-level state: only (T_cas, q_cas) is live. Vapour and CO2 reuse the
   ! heat transfer factor (z0q = z0h), CO2 on the molar capacity at the live humidity.
   !---~---
   pure subroutine refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp) :: tcas
      tcas = cas_enthalpy_to_temp(cas_enthalpy, cas_shv)
      !--- The pre-pass solved at tcas too, so only the state differs from its solve.
      call canopy_air_conductances(cas%aero_cfg, cas%env, cas%geom, cas%rho, tcas, cas_shv,              &
                                   g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine refresh_cas_conductances

   !---~---
   !   cas_conductances -- the one place that decides whether to re-solve the surface layer: a fresh
   ! solve when the bundle carries live MO inputs, else the bundle's own values. One rule keeps a
   ! tendency and the boundary flux charged against it in agreement (explicit_stage_bflux, RK45).
   !---~---
   pure subroutine cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      if (cas%mo_live) then
         call refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      else
         g_atm_heat = cas%box%atm_conductance_enthalpy
         g_atm_vapour = cas%box%atm_conductance_vapor
         g_atm_co2 = cas%box%atm_conductance_co2
      end if
   end subroutine cas_conductances

   !---~---
   !   surface_derivs -- the CAS surface-block RHS (leaf/wood energy balance + ground skin + the three
   ! CAS twins), side-effect-free. Every tissue flux is linearized about the current tcas, which is
   ! also the longwave emission base (docs/science/vegetation_energy_dynamics.md eqs. (1)-(2)); the
   ! twins' tendencies are docs/science/canopy_air_space_biophysics.md eq. (1).
   !---~---
   pure subroutine surface_derivs(y, cas, tissue, veg_surface, ground, t_ground, n, f)
      type(surface_state_t),        intent(in)  :: y
      type(cas_boundary_t),         intent(in)  :: cas      ! CAS <-> atmosphere boundary (a stage may pass its own conductances)
      type(plant_energy_env_t),     intent(in)  :: tissue   ! per-cohort leaf/wood energy coefficients
      type(canopy_water_out_t),     intent(in)  :: veg_surface     ! the canopy's surface water
      type(surface_env_t),          intent(in)  :: ground   ! what the ground surface reads, the snow pack among it
      real(wp),               intent(in)  :: t_ground   ! [K] soil-top temperature at this evaluation (live, not frozen)
      integer(ik),            intent(in)  :: n
      !---~---
      !   Filled in place; the per-cohort arrays are allocated only when missing or resized, so a
      ! Newton solve that evaluates this many times on one record allocates them once.
      !---~---
      type(surface_tend_t),   intent(inout) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt
      real(wp)    :: coh_h, coh_qw, coh_transp, coh_rnet, coh_surface_evap   ! the canopy's, into the canopy air
      type(surface_env_t) :: senv      ! the ground surface's surroundings at this evaluation
      real(wp)    :: gah_l, gaw_l, gac_l      ! the conductances this evaluation uses
      type(cas_env_t) :: cas_env

      call ensure_surface_tend(f, n)

      tcas   = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, cas%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, cas%press)

      !---~---
      !   Every cohort's wood and leaf at this canopy air (docs/science/vegetation_energy_dynamics.md);
      ! with canopy_water_on off, f_wet_c = 0 and the surface water's pathway vanishes.
      !---~---
      call canopy_tissue_energy(tissue, veg_surface, tcas, qcas, qsat_c, dqdt, cas%rho, n, f%leaf_temp, f%wood_temp,     &
                                f%transp_c, f%leaf_surface_evap, f%wood_surface_evap, coh_h, coh_qw, coh_transp,        &
                                coh_surface_evap, coh_rnet)

      !---~---
      !   The ground surface at this evaluation's canopy air and soil top, as the tissues' is at
      ! theirs: the snowfac-weighted pack and the bare soil (ed_surface).
      !---~---
      senv            = ground
      senv%can_temp   = tcas
      senv%can_shv    = qcas
      senv%rho_air    = cas%rho
      senv%cp_cas     = cas%cp
      senv%press      = cas%press
      senv%t_soil_top = t_ground
      call ground_surface_exchange(senv, f%h_ground, f%le_ground, f%g_top, f%ground_rad, f%snow_enth_in,         &
                                   f%snow_vap_out, f%org_enth_in, f%org_vap_out)

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + coh_surface_evap + ground%soil_evap + f%snow_vap_out
      if (ground%has_org) f%src_vap = f%src_vap + f%org_vap_out
      f%coh_rnet   = coh_rnet
      f%h_canopy   = coh_h

      !---~---
      !   The condensation sink (dew, fog; ed_canopy_air_space), taken from the canopy air's vapour
      ! and, as the condensate's liquid enthalpy, from its enthalpy.
      !---~---
      call cas_condensation(cas%box%air_mass_capacity, y%cas_shv, qsat_c, tcas, f%cond, f%cond_enth)
      f%src_vap  = f%src_vap  - f%cond
      f%src_enth = f%src_enth - f%cond_enth

      !---~---
      !   CAS box tendencies (shared kernel; the condensation adjustment above is folded into
      ! src_enth/src_vap, so the box math is scheme-agnostic).
      !---~---
      cas_env = cas%box
      cas_env%surface_enthalpy_source = f%src_enth
      cas_env%surface_vapor_source    = f%src_vap
      cas_env%biotic_co2_source       = cas%nee_biotic
      !---~---
      !   The CAS <-> atmosphere conductances follow the live canopy-air state: RK45 re-solves them
      ! here. ESDIRK2 arrives with mo_live cleared on its stage copy, having solved once for the
      ! stage, so its Newton's residual evaluations reuse that value.
      !---~---
      call cas_conductances(cas, y%cas_enthalpy, y%cas_shv, gah_l, gaw_l, gac_l)
      cas_env%atm_conductance_enthalpy = gah_l
      cas_env%atm_conductance_vapor    = gaw_l
      cas_env%atm_conductance_co2      = gac_l
      call cas_box_time_deriv(y%cas_enthalpy, y%cas_shv, y%cas_co2, cas_env,                 &
                                 f%d_cas_enthalpy, f%d_cas_shv, f%d_cas_co2)
   end subroutine surface_derivs

   !--- Allocate a surface-tendency record's per-cohort arrays for n cohorts, unless they are.
   pure subroutine ensure_surface_tend(f, n)
      type(surface_tend_t), intent(inout) :: f
      integer(ik),          intent(in)    :: n
      if (allocated(f%leaf_temp) .and. allocated(f%wood_temp) .and. allocated(f%transp_c) .and.   &
          allocated(f%leaf_surface_evap) .and. allocated(f%wood_surface_evap)) then
         if (size(f%leaf_temp) == n .and. size(f%wood_temp) == n .and. size(f%transp_c) == n .and. &
             size(f%leaf_surface_evap) == n .and. size(f%wood_surface_evap) == n) return
      end if
      if (allocated(f%leaf_temp))      deallocate(f%leaf_temp)
      if (allocated(f%wood_temp))      deallocate(f%wood_temp)
      if (allocated(f%transp_c))       deallocate(f%transp_c)
      if (allocated(f%leaf_surface_evap)) deallocate(f%leaf_surface_evap)
      if (allocated(f%wood_surface_evap)) deallocate(f%wood_surface_evap)
      allocate(f%leaf_temp(n), f%wood_temp(n), f%transp_c(n), f%leaf_surface_evap(n), f%wood_surface_evap(n))
   end subroutine ensure_surface_tend

   !---~---
   !   march_derivs -- the whole-column RHS an RK45 stage evaluates; commits nothing. Diagnoses the
   ! soil-top temperature from the state, runs the surface block, and assembles the soil-heat, soil-
   ! water and plant-water-mass tendencies. Water moves at the soil-plant water network's time-mean
   ! rates and the leaves transpire what the network gave them, so every stage sees the same water
   ! (docs/science/soil_plant_water.md §4, docs/science/numerical_scheme.md §5b).
   !---~---
   pure subroutine march_derivs(y, frozen, site, n, nsl, f, sf_out)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      integer(ik),           intent(in)  :: n      ! number of cohorts
      integer(ik),           intent(in)  :: nsl    ! number of active soil layers
      type(march_tend_t),   intent(out) :: f
      !---~---
      !   This stage's own surface tendencies, for RK45's whole-column boundary-flux ledger (as
      ! backward_euler_stage's surf_tend serves ESDIRK2's).
      !---~---
      type(surface_tend_t),  intent(out), optional :: sf_out

      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_state_t)  :: soil_energy
      type(soil_heat_env_t)     :: eforc
      real(wp)                   :: t_ground, fliq1, transp_i, e_drain, water_heat(nsl)
      integer(ik)                :: i

      allocate(f%d_leaf_water_mass(n), f%d_wood_water_mass(n), f%leaf_temp(n))
      allocate(f%d_leaf_surf_water(n), f%d_wood_surf_water(n))

      !--- The soil-top temperature from the current state, for the ground skin.
      call soil_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                &
                           site%soil_profile%layer_curve(1), t_ground, fliq1)

      !--- 1. Surface block (leaf + ground + CAS twins).
      y_stage%cas_enthalpy = y%cas_enthalpy
      y_stage%cas_shv = y%cas_shv
      y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%veg_surface, frozen%ground,                 &
                          t_ground, n, surf_tend)
      f%d_cas_enthalpy = surf_tend%d_cas_enthalpy
      f%d_cas_shv = surf_tend%d_cas_shv
      f%d_cas_co2 = surf_tend%d_cas_co2
      f%g_top = surf_tend%g_top
      f%leaf_temp(1:n) = surf_tend%leaf_temp(1:n)

      !---~---
      !   2. Soil water at the network's steady rate, which any consistent tableau integrates
      ! exactly, so RK45's theta ends where the network's does; so do drainage and uptake.
      !---~---
      f%dtheta_dt(1:nsl) = frozen%hydrology%theta_rate(1:nsl)
      f%drainage_rate    = frozen%hydrology%drainage
      f%uptake_rate      = frozen%roots%uptake

      !---~---
      !   3. Soil heat: g_top from the surface, and the enthalpy of the network's flows, the numbers
      ! the ESDIRK2 stage uses (docs/science/numerical_scheme.md §5b, soil_plant_water.md §5). The
      ! leaf pays its transpiration's full vapour enthalpy in surface_derivs.
      !---~---
      soil_energy%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      call set_soil_heat_env(nsl, surf_tend%g_top, y%theta, frozen%roots%root_heat, frozen%hydrology%exfil_enth,  &
                             frozen%hydrology%w_flux_frozen, frozen%hydrology%infiltration,                       &
                             frozen%hydrology%t_infil, f%drainage_rate, frozen%hydrology%t_bot, eforc,            &
                             water_heat, e_drain)
      call soil_energy_time_deriv(soil_energy, eforc, site%soil_profile, site%soil,   &
                                  f%dedt, bottom_face=f%g_bottom)

      !---~---
      !   4. Plant water mass per plant: the network's sapflow and uptake in, and out the
      ! transpiration this stage charged the canopy air, the network's own when the stages hold it.
      !---~---
      do i = 1_ik, n
         transp_i = surf_tend%transp_c(i) / max(frozen%plant%nplant(i), tiny_num)
         f%d_leaf_water_mass(i) = frozen%plant%sapflow_frozen(i) - transp_i
         f%d_wood_water_mass(i) = frozen%plant%uptake_frozen(i)  - frozen%plant%sapflow_frozen(i)
      end do

      !---~---
      !   5. Canopy surface water: frozen interception in, per-stage surface-water evaporation
      ! out, both per m2 of ground. Zero with canopy_water_on off.
      !---~---
      f%d_leaf_surf_water(1:n) = frozen%veg_surface%intercept_leaf(1:n) - surf_tend%leaf_surface_evap(1:n)
      f%d_wood_surf_water(1:n) = frozen%veg_surface%intercept_wood(1:n) - surf_tend%wood_surface_evap(1:n)
      if (present(sf_out)) sf_out = surf_tend
   end subroutine march_derivs

   !---~---
   !   The soil-plant water network (docs/science/soil_plant_water.md §1-§2): the pond, the
   ! soil layers, and each cohort's wood and leaves. The stores' curves and the paths'
   ! flows are the process libraries'; these join them.
   !---~---

   !---~---
   !   set_water_state -- the network's state from the model's: each organic slot's water and
   ! each soil layer's water content and ice, the pond's water and enthalpy, and each cohort's
   ! leaf and wood water [kg per plant]; each store's unknown is read from its water. A
   ! leafless cohort's leaves hold none, so their unknown starts at the wood's potential less
   ! the head to the crown, where an empty leaf comes to rest.
   !---~---
   subroutine set_water_state(static, profile, hyd, theta, ice, w_pond, e_pond, leaf_water, wood_water, water, org)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(hydraulics_table_t),  intent(in)    :: hyd
      real(wp),                  intent(in)    :: theta(:)                       ! (nsl) [m3/m3] liquid and ice
      type(soil_ice_t),          intent(in)    :: ice(:)                         ! (nsl) each layer's ice
      real(wp),                  intent(in)    :: w_pond, e_pond                 ! [kg/m2], [J/m2]
      real(wp),                  intent(in)    :: leaf_water(:), wood_water(:)   ! (ncoh) [kg/plant]
      type(water_work_t),        intent(inout) :: water
      type(water_organic_t), optional, intent(in) :: org                         ! the organic slots, if any
      integer(ik) :: k, i, j, n_org, soff
      real(wp)    :: psi, dpsi, w, dw, kk, dk, head1, elev
      n_org = 0_ik
      if (present(org)) n_org = org%n
      soff = 1_ik + n_org
      call ensure_water_capacity(water, profile%n_active, n_org, static%n)
      associate (st => water%state)
         do j = 1_ik, n_org
            st%tau_soil(1_ik + j)   = soil_theta_to_unknown(org%curve(j), org%ice(j),                            &
                                                            (org%water(j) + org%deficit(j)) / (rho_h2o * org%depth(j)))
            st%water_soil(1_ik + j) = org%water(j)
         end do
         do k = 1_ik, profile%n_active
            st%tau_soil(soff + k)   = soil_theta_to_unknown(profile%layer_curve(k), ice(k), theta(k))
            st%water_soil(soff + k) = rho_h2o * profile%dz(k) * theta(k)
         end do
         !---~---
         !   The pond starts at its water's head, but no higher than the head at which it neither
         ! gives nor takes water from the node below it, the top slot or the top layer: over dry
         ! ground a trace pond empties within the step, and its own head would start Newton at a
         ! flow thousands of times what it holds.
         !---~---
         st%water_soil(POND) = w_pond
         st%e_pond = e_pond
         if (n_org > 0_ik) then
            call soil_water_relations(org%curve(1), org%ice(1), st%tau_soil(2), psi, dpsi, w, dw, kk, dk)
            elev = -0.5_wp * org%depth(1)
         else
            call soil_water_relations(profile%layer_curve(1), ice(1), st%tau_soil(2), psi, dpsi, w, dw, kk, dk)
            elev = profile%z_node(1)
         end if
         head1 = psi / grav_head
         st%tau_soil(POND) = min(pond_water_to_head(w_pond), head1 + elev)
         do i = 1_ik, static%n
            associate (hp => hyd%pft(static%pft(i)))
               st%tau_coh(WOOD, i)   = tissue_rwc_to_unknown(                                                  &
                                          tissue_water_to_rwc(hp%wood_curve, static%bwood_water(i), wood_water(i)))
               st%water_coh(WOOD, i) = wood_water(i) * static%nplant(i)
               if (hp%leaf_curve%water_sat * static%bleaf(i) > tiny_num) then
                  st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(                                               &
                                             tissue_water_to_rwc(hp%leaf_curve, static%bleaf(i), leaf_water(i)))
                  st%water_coh(LEAF, i) = leaf_water(i) * static%nplant(i)
               else
                  call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i),            &
                                              st%tau_coh(WOOD, i), psi, dpsi, w, dw)
                  st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(tissue_psi_to_rwc(psi - static%crown_head(i),   &
                                             hp%leaf_curve%pi0, hp%leaf_curve%elastic_mod))
                  st%water_coh(LEAF, i) = 0.0_wp
               end if
            end associate
         end do
      end associate
   end subroutine set_water_state

   !---~---
   !   water_derivs -- the network at one state: each store's potential and water from its
   ! curve; the flow along every path (pond to the top node, down the chain of organic slots and
   ! soil layers, the slots' and the top layer's evaporation, drainage, each cohort's uptake,
   ! sapflow and transpiration); the net inflow into each store; and, when asked, the flows'
   ! Jacobian in the stores' unknowns (sys). The core is a chain, face c joining node c to node
   ! c + 1. A node detached from the solve (cut) exchanges nothing with its neighbours: every
   ! face it touches and its roots carry nothing, though each face's slopes are still kept (in
   ! ev's face_dqa and face_dqb, and root_link), so the solve can see when its links return.
   !---~---
   subroutine water_derivs(static, profile, soil, hyd, drive, tau_soil, tau_coh, ev, sys, with_jacobian, cut)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(in)    :: drive
      real(wp),                  intent(in)    :: tau_soil(:), tau_coh(:,:)
      type(water_tend_t),        intent(inout) :: ev
      type(arrowhead_system_t),  intent(inout) :: sys   ! the flows' Jacobian, written when with_jacobian
      logical,                   intent(in)    :: with_jacobian
      logical,         optional, intent(in)    :: cut(:) ! (nsoil) nodes detached from the solve
      integer(ik) :: k, i, ca, cb, c, nsl, nsoil, n_org, soff
      real(wp) :: q, dq_pa, dq_pb, dq_ka, dq_kb, dqa, dqb, u, s, ds_w, ds_l
      real(wp) :: e, de_dtheta, de_dpsi, de, qbot, dqbot_dpsi, dqbot_dk, theta1, span, supply, dsupply
      real(wp) :: psi_w, dpsi_w, u_sum, dqb_sum, du_dpsi(profile%n_active), du_dk(profile%n_active)
      real(wp) :: de_own, w_b, dist
      logical  :: wet, is_cut(size(tau_soil))

      nsl   = profile%n_active
      n_org = drive%org%n
      soff  = 1_ik + n_org
      nsoil = soff + nsl
      is_cut = .false.
      if (present(cut)) is_cut = cut(1:nsoil)

      ev%inflow_soil = 0.0_wp
      ev%inflow_coh = 0.0_wp
      !---~---
      !   The Jacobian's entries are written afresh: the core's band, which the flows add into, is
      ! cleared; every other entry an evaluation touches is set outright below, and the rest stay
      ! zero from arrowhead_allocate.
      !---~---
      if (with_jacobian) then
         do c = 1_ik, nsoil
            sys%a(c, c) = 0.0_wp
            if (c < nsoil) then
               sys%a(c, c + 1_ik) = 0.0_wp
               sys%a(c + 1_ik, c) = 0.0_wp
            end if
         end do
      end if

      !--- Storage.
      do c = 1_ik, nsoil
         call core_node_relations(profile, drive, c, tau_soil(c), ev%psi_soil(c), ev%dpsi_soil(c), ev%water_soil(c), &
                                  ev%dwater_soil(c), ev%k_soil(c), ev%dk_soil(c))
      end do
      do i = 1_ik, static%n
         associate (hp => hyd%pft(static%pft(i)))
            call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i), tau_coh(WOOD, i), &
                                        ev%psi_coh(WOOD, i), ev%dpsi_coh(WOOD, i), ev%water_coh(WOOD, i),        &
                                        ev%dwater_coh(WOOD, i))
            call tissue_water_relations(hp%leaf_curve, static%bleaf(i), static%nplant(i), tau_coh(LEAF, i),       &
                                        ev%psi_coh(LEAF, i), ev%dpsi_coh(LEAF, i), ev%water_coh(LEAF, i),        &
                                        ev%dwater_coh(LEAF, i))
         end associate
      end do

      !--- Down the chain: the pond to the top node, through the organic slots, between soil layers.
      ev%org_face = 0.0_wp
      do c = 1_ik, nsoil - 1_ik
         ca = c
         cb = c + 1_ik
         call core_face_geometry(profile, drive, c, w_b, dist, wet)
         call core_face_flux(wet, ev%k_soil(ca), ev%k_soil(cb), w_b, ev%psi_soil(ca), ev%psi_soil(cb), dist,          &
                             q, dq_pa, dq_pb, dq_ka, dq_kb)
         dqa = dq_pa * ev%dpsi_soil(ca) + dq_ka * ev%dk_soil(ca)
         dqb = dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb)
         ev%face_dqa(c) = dqa
         ev%face_dqb(c) = dqb
         if (is_cut(ca) .or. is_cut(cb)) then
            q   = 0.0_wp
            dqa = 0.0_wp
            dqb = 0.0_wp
         end if
         if (c == 1_ik) then
            ev%pond_flux = q
         else if (c <= n_org + 1_ik) then
            ev%org_face(c - 1_ik) = q
         else
            ev%face_flux(c - soff) = q
         end if
         call add_soil_flow(ca, cb, q, dqa, dqb)
      end do

      !---~---
      !   Rain into the pond; the slots' and the top layer's evaporation; drainage from the bottom
      ! one. The dry surface layer reads the liquid and dsl_theta_init of the ice, so its onset is
      ! dsl_theta_init of the pores the ice leaves (CLM5's effective porosity).
      !---~---
      ev%inflow_soil(POND) = ev%inflow_soil(POND) + drive%rain
      !---~---
      !   The organic layer's own evaporation, from its top slot by the soil's law (Philip's pore
      ! humidity, implicit in the slot's potential, the dried top's path in series); the slots below
      ! feed it through their faces.
      !---~---
      ev%org_evap = 0.0_wp
      if (n_org > 0_ik) then
         c = 2_ik
         call slot_evaporation(drive, soil, ev%water_soil(c), ev%dwater_soil(c), ev%psi_soil(c), ev%dpsi_soil(c),  &
                               e, de_own)
         ev%org_evap(1) = e
         ev%inflow_soil(c) = ev%inflow_soil(c) - e
         if (with_jacobian) sys%a(c, c) = sys%a(c, c) + de_own
      end if
      c = soff + 1_ik
      theta1 = ev%water_soil(c) / (rho_h2o * profile%dz(1)) - (1.0_wp - soil%dsl_theta_init) * drive%ice(1)%theta_ice
      call ground_evaporation_rate(theta1, ev%psi_soil(c), drive%ground, soil, e, de_dtheta, de_dpsi)
      !--- The layer supplies it while it holds water above residual, ice included.
      span = profile%theta_sat(1) - profile%theta_res(1)
      call evaporation_supply((ev%water_soil(c) / (rho_h2o * profile%dz(1)) - profile%theta_res(1)) / span,     &
                              supply, dsupply)
      de_dtheta = supply * de_dtheta + e * dsupply / span
      de_dpsi   = supply * de_dpsi
      e         = supply * e
      ev%evap = e
      ev%inflow_soil(c) = ev%inflow_soil(c) - e
      if (with_jacobian) sys%a(c, c) = sys%a(c, c)                                                              &
                         + de_dtheta * ev%dwater_soil(c) / (rho_h2o * profile%dz(1))                               &
                         + de_dpsi * ev%dpsi_soil(c)
      c = soff + nsl
      call bottom_flux(profile, soil, nsl, ev%psi_soil(c), ev%k_soil(c), drive%ice(nsl)%impedance, qbot,           &
                       dqbot_dpsi, dqbot_dk)
      ev%drain = rho_h2o * qbot
      ev%inflow_soil(c) = ev%inflow_soil(c) - ev%drain
      if (with_jacobian) sys%a(c, c) = sys%a(c, c)                                                              &
                         + rho_h2o * (dqbot_dpsi * ev%dpsi_soil(c) + dqbot_dk * ev%dk_soil(c))

      !---~---
      !   Each cohort: uptake from its layers, sapflow to its leaves, and transpiration. The
      ! cohort's uptake comes from its whole root profile at once, each layer's conductance
      ! times its roots' contact with it (none once it freezes); its own sums are kept in
      ! scalars over its layers and written once. The roots take nothing from the organic slots,
      ! nor from a layer detached from the solve.
      !---~---
      ev%root_link = 0.0_wp
      do i = 1_ik, static%n
         psi_w = ev%psi_coh(WOOD, i)
         dpsi_w = ev%dpsi_coh(WOOD, i)
         call root_water_uptake_profile(nsl, static%g_rhizo(1:nsl, i) * drive%root_contact(1:nsl),                  &
                                        ev%k_soil(soff+1:nsoil), ev%psi_soil(soff+1:nsoil), profile%z_node(1:nsl),    &
                                        psi_w,                                                                        &
                                        ev%uptake(1:nsl, i), du_dpsi, du_dk)
         u_sum = 0.0_wp
         dqb_sum = 0.0_wp
         do k = 1_ik, nsl
            c = soff + k
            ev%root_link(c) = max(ev%root_link(c), abs(du_dpsi(k) * ev%dpsi_soil(c) + du_dk(k) * ev%dk_soil(c)),   &
                                  abs(du_dpsi(k) * dpsi_w))
            if (is_cut(c)) ev%uptake(k, i) = 0.0_wp
            if (static%g_rhizo(k, i) <= 0.0_wp .or. is_cut(c)) then
               if (with_jacobian) then
                  sys%b(c, i) = 0.0_wp
                  sys%c(c, i) = 0.0_wp
               end if
               cycle
            end if
            u = ev%uptake(k, i)
            ev%inflow_soil(c) = ev%inflow_soil(c) - u
            u_sum = u_sum + u
            if (with_jacobian) then
               dqa = du_dpsi(k) * ev%dpsi_soil(c) + du_dk(k) * ev%dk_soil(c)
               dqb = -du_dpsi(k) * dpsi_w
               sys%a(c, c) = sys%a(c, c) + dqa
               sys%b(c, i) = dqb
               sys%c(c, i) = -dqa
               dqb_sum = dqb_sum - dqb
            end if
         end do
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) + u_sum
         associate (hp => hyd%pft(static%pft(i)))
            call xylem_sapflow(static%k_plant(i), ev%psi_coh(WOOD, i), ev%psi_coh(LEAF, i), static%crown_head(i),  &
                               hp%wood_psi50, hp%wood_kexp, hp%vuln_table, s, ds_w, ds_l)
         end associate
         ev%sapflow(i) = s
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) - s
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) + s
         call supply_limited_transpiration(drive%e_demand(i), tau_coh(LEAF, i), hyd%pft(static%pft(i))%rwc_close,   &
                                           e, de)
         ev%transp(i) = e
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) - e
         if (with_jacobian) then
            dqa = ds_w * ev%dpsi_coh(WOOD, i)
            dqb = ds_l * ev%dpsi_coh(LEAF, i)
            sys%d(WOOD, WOOD, i) = dqb_sum + dqa
            sys%d(WOOD, LEAF, i) = dqb
            sys%d(LEAF, WOOD, i) = -dqa
            sys%d(LEAF, LEAF, i) = -dqb + de
         end if
      end do

   contains

      !---~---
      !   A flow q from soil node a to soil node b, with dq/dtau_a and dq/dtau_b: the inflows
      ! and, when wanted, the Jacobian of the residual (W - W0)/h - inflow.
      !---~---
      subroutine add_soil_flow(a, b, q, dqa, dqb)
         integer(ik), intent(in) :: a, b
         real(wp),    intent(in) :: q, dqa, dqb
         ev%inflow_soil(a) = ev%inflow_soil(a) - q
         ev%inflow_soil(b) = ev%inflow_soil(b) + q
         if (.not. with_jacobian) return
         sys%a(a, a) = sys%a(a, a) + dqa
         sys%a(a, b) = sys%a(a, b) + dqb
         sys%a(b, a) = sys%a(b, a) - dqa
         sys%a(b, b) = sys%a(b, b) - dqb
      end subroutine add_soil_flow

   end subroutine water_derivs

   !---~---
   !   Core node c's potential, water and conductivity at its unknown, with their slopes: the
   ! pond (its conductivity the node below's, saturated, at that node's ice impedance, since its
   ! path lies in it), an organic slot, or a soil layer.
   !---~---
   pure subroutine core_node_relations(profile, drive, c, tau, psi, dpsi, w, dw, kc, dkc)
      type(soil_profile_t), intent(in)  :: profile
      type(water_drive_t),  intent(in)  :: drive
      integer(ik),          intent(in)  :: c
      real(wp),             intent(in)  :: tau
      real(wp),             intent(out) :: psi, dpsi, w, dw, kc, dkc
      integer(ik) :: soff
      soff = 1_ik + drive%org%n
      if (c == POND) then
         call pond_water_relations(tau, psi, dpsi, w, dw)
         if (drive%org%n > 0_ik) then
            kc = drive%org%curve(1)%ksat * drive%org%ice(1)%impedance
         else
            kc = profile%ksat(1) * drive%ice(1)%impedance
         end if
         dkc = 0.0_wp
      else if (c <= soff) then
         call soil_water_relations(drive%org%curve(c - 1_ik), drive%org%ice(c - 1_ik), tau, psi, dpsi, w, dw, kc, dkc)
         w = w - drive%org%deficit(c - 1_ik)      ! a slot below its residual: its curve's water less the lack
      else
         call soil_water_relations(profile%layer_curve(c - soff), drive%ice(c - soff), tau, psi, dpsi, w, dw, kc, dkc)
      end if
   end subroutine core_node_relations

   !---~---
   !   Face c of the chain, joining node c to node c + 1: the weight of node c + 1's conductivity
   ! in the face's, the distance between the two nodes, which is also the drop in elevation, and
   ! whether it wets (core_face_flux): the pond's face to the top slot and the faces between
   ! slots. The bottom slot's face into the top soil layer takes the soil's thickness-weighted
   ! mean, as a soil face, so the soil's conductivity limits what enters it: upstream, a saturated
   ! moss would drain into a clay at the moss's conductivity, tens of kilograms a second. The pond
   ! counts as a saturated layer as thick as the node below it.
   !---~---
   pure subroutine core_face_geometry(profile, drive, c, w_b, dist, wet)
      type(soil_profile_t), intent(in)  :: profile
      type(water_drive_t),  intent(in)  :: drive
      integer(ik),          intent(in)  :: c
      real(wp),             intent(out) :: w_b, dist
      logical,              intent(out) :: wet
      integer(ik) :: n_org, j
      n_org = drive%org%n
      if (c == 1_ik) then
         w_b = 0.5_wp
         if (n_org > 0_ik) then
            dist = 0.5_wp * drive%org%depth(1)
         else
            dist = -profile%z_node(1)
         end if
         wet = n_org > 0_ik
      else if (c <= n_org + 1_ik) then
         j = c - 1_ik
         if (j < n_org) then
            w_b  = drive%org%depth(j) / (drive%org%depth(j) + drive%org%depth(j + 1_ik))
            dist = 0.5_wp * (drive%org%depth(j) + drive%org%depth(j + 1_ik))
         else
            w_b  = drive%org%depth(j) / (drive%org%depth(j) + profile%dz(1))
            dist = 0.5_wp * drive%org%depth(j) - profile%z_node(1)
         end if
         wet = j < n_org        ! the bottom slot's face into the soil: the soil's conductivity limits it
      else
         j = c - 1_ik - n_org
         w_b  = profile%face_weight(j)
         dist = profile%dz_node(j)
         wet  = .false.
      end if
   end subroutine core_face_geometry

   !---~---
   !   The flow down a face of the chain [kg/m2/s] and its slopes. A soil face takes the
   ! thickness-weighted geometric mean of its nodes' conductivities (soil_water_flux). A face
   ! touching an organic slot wets: it takes the upstream node's conductivity (Forsyth, Wu &
   ! Pruess 1995), so water reaches a dry slot, or leaves a wet one into dry ground, at the wet
   ! side's conductivity rather than at the vanishing one of the dry side; near no flow, over a
   ! drive of WET_BLEND_WIDTH either way, the two blend smoothly through their arithmetic mean,
   ! so Newton meets no switch.
   !---~---
   pure subroutine core_face_flux(wet, k_a, k_b, w_b, psi_a, psi_b, dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
      logical,  intent(in)  :: wet
      real(wp), intent(in)  :: k_a, k_b, w_b, psi_a, psi_b, dist
      real(wp), intent(out) :: q, dq_pa, dq_pb, dq_ka, dq_kb
      real(wp) :: drive, x, sb, dsb, k_face
      if (.not. wet) then
         call soil_water_flux(k_a, k_b, w_b, psi_a, psi_b, dist, dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
         return
      end if
      drive = ((psi_a - psi_b) / grav_head + dist) / dist
      x     = max(-1.0_wp, min(1.0_wp, drive / WET_BLEND_WIDTH))
      sb    = 0.25_wp * (x + 1.0_wp) ** 2 * (2.0_wp - x)              ! 0 upstream b, 1 upstream a
      dsb   = 0.75_wp * (1.0_wp - x * x) / WET_BLEND_WIDTH
      k_face = k_b + (k_a - k_b) * sb
      q      = rho_h2o * k_face * drive
      dq_pa  = rho_h2o * (k_face + (k_a - k_b) * dsb * drive) / (grav_head * dist)
      dq_pb  = -dq_pa
      dq_ka  = rho_h2o * drive * sb
      dq_kb  = rho_h2o * drive * (1.0_wp - sb)
   end subroutine core_face_flux

   !---~---
   !   The top organic slot's evaporation [kg/m2/s] and its slope in the slot's unknown, by the
   ! soil's own law (ground_evaporation_rate): Philip's pore humidity at the slot's potential,
   ! smoothly stopping where it meets the canopy air's (no dew), through the dried top's path in
   ! series with the ground's conductance; and, as the top soil layer's, only while it holds water
   ! above its residual, ice included (evaporation_supply). A drying slot comes to rest at its
   ! equilibrium with the canopy air's humidity.
   !---~---
   pure subroutine slot_evaporation(drive, soil, w, dw, psi, dpsi, e, de)
      type(water_drive_t), intent(in)  :: drive
      type(soil_params_t), intent(in)  :: soil
      real(wp),            intent(in)  :: w, dw, psi, dpsi   ! the top slot's water [kg/m2], potential [MPa], slopes
      real(wp),            intent(out) :: e, de
      real(wp) :: theta, de_dtheta, de_dpsi, span, supply, dsupply
      theta = w / (rho_h2o * drive%org%depth(1))
      call ground_evaporation_rate(theta, psi, drive%org%air, soil, e, de_dtheta, de_dpsi)
      span = drive%org%curve(1)%theta_sat - drive%org%curve(1)%theta_res
      call evaporation_supply((theta - drive%org%curve(1)%theta_res) / span, supply, dsupply)
      de_dtheta = supply * de_dtheta + e * dsupply / span
      de_dpsi   = supply * de_dpsi
      e         = supply * e
      de        = de_dtheta * dw / (rho_h2o * drive%org%depth(1)) + de_dpsi * dpsi
   end subroutine slot_evaporation

   !--- Each cohort's wood potential [MPa] at its unknown, which a soil node's uptake reads.
   subroutine wood_potentials(static, hyd, tau_coh, psi_wood)
      type(fast_static_state_t), intent(in)  :: static
      type(hydraulics_table_t),  intent(in)  :: hyd
      real(wp),                  intent(in)  :: tau_coh(:,:)
      real(wp),                  intent(out) :: psi_wood(:)
      real(wp)    :: dpsi, w, dw
      integer(ik) :: i
      do i = 1_ik, static%n
         call tissue_water_relations(hyd%pft(static%pft(i))%wood_curve, static%bwood_water(i), static%nplant(i),   &
                                     tau_coh(WOOD, i), psi_wood(i), dpsi, w, dw)
      end do
   end subroutine wood_potentials

   !---~---
   !   node_balance -- core node c's own backward-Euler balance over a sub-step of h, R = (W -
   ! W0)/h - inflow [kg/m2/s], and dR/dtau_c, at its unknown tau with every other store held at
   ! its state in tau_soil and its neighbours' relations given: the faces to the nodes above and
   ! below (none to a detached node, cut), and its own sources -- rain into the pond, the top
   ! slot's evaporation, the top layer's evaporation,
   ! the bottom's drainage and each cohort's uptake at its wood's potential. The nonlinear
   ! elimination of a fast node and the closed form of a detached one solve R = 0 on it.
   !---~---
   subroutine node_balance(static, profile, soil, drive, c, tau_soil, w0, h, cut, psi_wood, r, dr)
      type(fast_static_state_t), intent(in)  :: static
      type(soil_profile_t),      intent(in)  :: profile
      type(soil_params_t),       intent(in)  :: soil
      type(water_drive_t),       intent(in)  :: drive
      integer(ik),               intent(in)  :: c
      real(wp),                  intent(in)  :: tau_soil(:)  ! (nsoil) the node's trial unknown among the others'
      real(wp),                  intent(in)  :: w0, h        ! [kg/m2] its water at the sub-step's start; [s]
      logical,                   intent(in)  :: cut          ! detached: no faces, no roots
      real(wp),                  intent(in)  :: psi_wood(:)  ! (ncoh) [MPa]
      real(wp),                  intent(out) :: r, dr
      real(wp) :: psi, dpsi, w, dw, kc, dkc, q, dq_pa, dq_pb, dq_ka, dq_kb, w_b, dist, inflow, dinflow
      real(wp) :: e, de, de_dtheta, de_dpsi, theta1, span, supply, dsupply, qbot, dqbot_dpsi, dqbot_dk
      real(wp) :: u, du_ds, du_dw, du_dk, nbr(4), x(6)
      integer(ik) :: nsl, soff, nsoil, k, i
      logical :: wet
      nsl   = profile%n_active
      soff  = 1_ik + drive%org%n
      nsoil = soff + nsl
      call core_node_relations(profile, drive, c, tau_soil(c), psi, dpsi, w, dw, kc, dkc)
      !--- The neighbours' potentials and conductivities, held.
      nbr = 0.0_wp
      if (c > 1_ik) then
         call core_node_relations(profile, drive, c - 1_ik, tau_soil(c - 1_ik), x(1), x(2), x(3), x(4), x(5), x(6))
         nbr(1) = x(1)
         nbr(2) = x(5)
      end if
      if (c < nsoil) then
         call core_node_relations(profile, drive, c + 1_ik, tau_soil(c + 1_ik), x(1), x(2), x(3), x(4), x(5), x(6))
         nbr(3) = x(1)
         nbr(4) = x(5)
      end if
      inflow  = 0.0_wp
      dinflow = 0.0_wp
      if (.not. cut) then
         if (c > 1_ik) then
            call core_face_geometry(profile, drive, c - 1_ik, w_b, dist, wet)
            call core_face_flux(wet, nbr(2), kc, w_b, nbr(1), psi, dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
            inflow  = inflow + q
            dinflow = dinflow + dq_pb * dpsi + dq_kb * dkc
         end if
         if (c < nsoil) then
            call core_face_geometry(profile, drive, c, w_b, dist, wet)
            call core_face_flux(wet, kc, nbr(4), w_b, psi, nbr(3), dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
            inflow  = inflow - q
            dinflow = dinflow - (dq_pa * dpsi + dq_ka * dkc)
         end if
      end if
      if (c == POND) inflow = inflow + drive%rain
      if (c == 2_ik .and. soff > 1_ik) then
         call slot_evaporation(drive, soil, w, dw, psi, dpsi, e, de)
         inflow  = inflow - e
         dinflow = dinflow - de
      end if
      if (c == soff + 1_ik) then
         theta1 = w / (rho_h2o * profile%dz(1)) - (1.0_wp - soil%dsl_theta_init) * drive%ice(1)%theta_ice
         call ground_evaporation_rate(theta1, psi, drive%ground, soil, e, de_dtheta, de_dpsi)
         span = profile%theta_sat(1) - profile%theta_res(1)
         call evaporation_supply((w / (rho_h2o * profile%dz(1)) - profile%theta_res(1)) / span, supply, dsupply)
         de_dtheta = supply * de_dtheta + e * dsupply / span
         de_dpsi   = supply * de_dpsi
         inflow  = inflow - supply * e
         dinflow = dinflow - (de_dtheta * dw / (rho_h2o * profile%dz(1)) + de_dpsi * dpsi)
      end if
      if (c == nsoil) then
         call bottom_flux(profile, soil, nsl, psi, kc, drive%ice(nsl)%impedance, qbot, dqbot_dpsi, dqbot_dk)
         inflow  = inflow - rho_h2o * qbot
         dinflow = dinflow - rho_h2o * (dqbot_dpsi * dpsi + dqbot_dk * dkc)
      end if
      if (c > soff .and. .not. cut) then
         k = c - soff
         do i = 1_ik, static%n
            if (static%g_rhizo(k, i) <= 0.0_wp) cycle
            call root_water_uptake(static%g_rhizo(k, i) * drive%root_contact(k), kc, psi, profile%z_node(k),        &
                                   psi_wood(i), u, du_ds, du_dw, du_dk)
            inflow  = inflow - u
            dinflow = dinflow - (du_ds * dpsi + du_dk * dkc)
         end do
      end if
      r  = (w - w0) / h - inflow
      dr = dw / h - dinflow
   end subroutine node_balance

end module ed_fast_integ_derivs
