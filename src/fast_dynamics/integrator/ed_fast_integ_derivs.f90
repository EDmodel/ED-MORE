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
   use ed_soil_heat,   only : soil_heat_env_t, soil_energy_time_deriv, set_soil_heat_env
   use ed_ground_surface_lib, only : porous_layer_temp
   use ed_surface, only : surface_env_t, ground_surface_exchange, N_SURFACE_STORE
   use ed_patch_state_vars,   only : soil_energy_state_t
   use ed_canopy_air_space, only : cas_env_t, cas_box_time_deriv, cas_condensation
   use ed_canopy_aerodynamics, only : canopy_air_conductances
   use ed_plant_energy,   only : plant_energy_env_t, canopy_tissue_energy
   use ed_canopy_water,   only : canopy_water_out_t
   use ed_fast_integ_types, only : surface_state_t, surface_tend_t, cas_boundary_t,                                  &
                               march_state_t, frozen_t, march_tend_t, fast_static_state_t,                    &
                               water_drive_t, water_tend_t, water_work_t, ensure_water_capacity, POND, WOOD, LEAF, &
                               water_porous_t, water_state_t
   use ed_site_state_vars, only : site_surface_t
   use ed_soil_hydr_lib,   only : soil_water_relations, soil_theta_to_unknown, pond_water_relations,              &
                                  pond_water_to_head, soil_ice_t
   use ed_plant_hydr_lib,   only : tissue_water_relations, tissue_water_to_rwc, tissue_rwc_to_unknown,            &
                                   tissue_psi_to_rwc
   use ed_plant_hydraulics,   only : supply_limited_transpiration, root_water_uptake_profile, xylem_sapflow,        &
                                     root_water_uptake
   use ed_soil_water, only : soil_water_flux, bottom_flux, dry_layer_resistance, ground_vapour, store_vapour_path,  &
                             vapour_path_t
   use ed_linear_algebra, only : arrowhead_system_t
   use ed_numerics,       only : smooth_max0
   use ed_soil_params, only : soil_profile_t, soil_params_t
   use ed_pft_params, only : hydraulics_table_t
   use ed_ground_surface_params, only : N_POROUS
   implicit none
   private

   public :: surface_derivs, canopy_demand, march_derivs, cas_conductances
   public :: set_water_state, unknowns_from_water, water_derivs, node_balance, core_node_relations, wood_potentials

   !---~---
   !   [-] The drive (head gradient) over which the pond's face turns from the conductivity of the
   ! node below to the pond's as the flow turns down (core_face_flux): a centimetre of head per
   ! metre, far below the gradients of a wetting front and well above round-off. A numerical
   ! choice.
   !---~---
   real(wp), parameter :: WET_BLEND_WIDTH = 1.0e-2_wp
   !--- The most stores that evaporate: every porous layer and the top soil layer.
   integer(ik), parameter :: N_EVAP_MAX = N_POROUS + 1_ik

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

      real(wp)    :: tcas, qcas, qsat_c
      real(wp)    :: coh_h, coh_qw, coh_transp, coh_rnet, coh_surface_evap   ! the canopy's, into the canopy air
      type(surface_env_t) :: senv      ! the ground surface's surroundings at this evaluation
      real(wp)    :: gah_l, gaw_l, gac_l      ! the conductances this evaluation uses
      type(cas_env_t) :: cas_env
      integer(ik) :: k

      call canopy_tissue_part(y, cas, tissue, veg_surface, n, f, tcas, qcas, qsat_c, coh_h, coh_qw, coh_transp,       &
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
      call ground_surface_exchange(senv, f%h_ground, f%le_ground, f%g_top, f%ground_rad, f%stores, f%stack_newton,  &
                                   f%stack_fallback)

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + coh_surface_evap + ground%soil_evap
      do k = 1_ik, N_SURFACE_STORE
         f%src_vap = f%src_vap + f%stores%vapour_out(k)
      end do
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
   !---~---
   !   The canopy's demand at one canopy-air state: each cohort's wood and leaf energy balance, its
   ! transpiration and its surface water's evaporation, without the ground, whose exchange the
   ! tissues' does not depend on (surface_derivs is this and the ground and the canopy air).
   !---~---
   pure subroutine canopy_demand(y, cas, tissue, veg_surface, n, f)
      type(surface_state_t),    intent(in)    :: y
      type(cas_boundary_t),     intent(in)    :: cas
      type(plant_energy_env_t), intent(in)    :: tissue
      type(canopy_water_out_t), intent(in)    :: veg_surface
      integer(ik),              intent(in)    :: n
      type(surface_tend_t),     intent(inout) :: f
      real(wp) :: tcas, qcas, qsat_c, coh_h, coh_qw, coh_transp, coh_surface_evap, coh_rnet
      call canopy_tissue_part(y, cas, tissue, veg_surface, n, f, tcas, qcas, qsat_c, coh_h, coh_qw, coh_transp,       &
                              coh_surface_evap, coh_rnet)
   end subroutine canopy_demand

   !---~---
   !   Every cohort's wood and leaf at this canopy air (docs/science/vegetation_energy_dynamics.md);
   ! with canopy_water_on off, f_wet_c = 0 and the surface water's pathway vanishes. Returns the
   ! canopy air's temperature, humidity and saturation the tissues saw, and the canopy's sums.
   !---~---
   pure subroutine canopy_tissue_part(y, cas, tissue, veg_surface, n, f, tcas, qcas, qsat_c, coh_h, coh_qw,       &
                                      coh_transp, coh_surface_evap, coh_rnet)
      type(surface_state_t),    intent(in)    :: y
      type(cas_boundary_t),     intent(in)    :: cas
      type(plant_energy_env_t), intent(in)    :: tissue
      type(canopy_water_out_t), intent(in)    :: veg_surface
      integer(ik),              intent(in)    :: n
      type(surface_tend_t),     intent(inout) :: f
      real(wp),                 intent(out)   :: tcas, qcas, qsat_c                    ! [K], [kg/kg], [kg/kg]
      real(wp),                 intent(out)   :: coh_h, coh_qw, coh_transp, coh_surface_evap, coh_rnet
      real(wp) :: dqdt
      call ensure_surface_tend(f, n)
      tcas   = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, cas%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, cas%press)
      call canopy_tissue_energy(tissue, veg_surface, tcas, qcas, qsat_c, dqdt, cas%rho, n, f%leaf_temp, f%wood_temp,     &
                                f%transp_c, f%leaf_surface_evap, f%wood_surface_evap, coh_h, coh_qw, coh_transp,        &
                                coh_surface_evap, coh_rnet)
   end subroutine canopy_tissue_part

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
      call porous_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                &
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
      !   3. Soil heat: g_top from the surface, the heat of the water crossing the soil's top
      ! included, and the enthalpy of the network's flows below it, the numbers the ESDIRK2 stage
      ! uses (docs/science/numerical_scheme.md §5b, soil_plant_water.md §5). The leaf pays its
      ! transpiration's full vapour enthalpy in surface_derivs.
      !---~---
      soil_energy%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      call set_soil_heat_env(nsl, surf_tend%g_top, y%theta, frozen%roots%root_heat,                              &
                             frozen%hydrology%w_flux_frozen, f%drainage_rate, frozen%hydrology%t_bot, eforc,      &
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
   !   set_water_state -- the network's state from the model's: each porous layer's water and
   ! each soil layer's water content and ice, the pond's water, and each cohort's leaf and wood
   ! water [kg per plant]; each store's unknown is read from its water (unknowns_from_water).
   !---~---
   subroutine set_water_state(static, profile, hyd, theta, ice, w_pond, leaf_water, wood_water, water, porous)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(hydraulics_table_t),  intent(in)    :: hyd
      real(wp),                  intent(in)    :: theta(:)                       ! (nsl) [m3/m3] liquid and ice
      type(soil_ice_t),          intent(in)    :: ice(:)                         ! (nsl) each layer's ice
      real(wp),                  intent(in)    :: w_pond                         ! [kg/m2]
      real(wp),                  intent(in)    :: leaf_water(:), wood_water(:)   ! (ncoh) [kg/plant]
      type(water_work_t),        intent(inout) :: water
      type(water_porous_t), optional, intent(in) :: porous                       ! the porous layers, if any
      type(water_porous_t) :: none
      integer(ik) :: k, i, j, n_porous, soff
      n_porous = 0_ik
      if (present(porous)) n_porous = porous%n
      soff = 1_ik + n_porous
      call ensure_water_capacity(water, profile%n_active, n_porous, static%n)
      associate (st => water%state)
         st%water_soil(POND) = w_pond
         do j = 1_ik, n_porous
            st%water_soil(1_ik + j) = porous%water(j)
         end do
         do k = 1_ik, profile%n_active
            st%water_soil(soff + k) = rho_h2o * profile%dz(k) * theta(k)
         end do
         do i = 1_ik, static%n
            st%water_coh(WOOD, i) = wood_water(i) * static%nplant(i)
            st%water_coh(LEAF, i) = leaf_water(i) * static%nplant(i)
         end do
      end associate
      if (present(porous)) then
         call unknowns_from_water(static, profile, hyd, ice, porous, water%state)
      else
         call unknowns_from_water(static, profile, hyd, ice, none, water%state)
      end if
   end subroutine set_water_state

   !---~---
   !   unknowns_from_water -- each store's unknown read from the water it holds, the inverse of its
   ! curve: a porous or soil layer's with its ice held, a porous layer below its residual at its
   ! curve's dry end (its deficit carried as the offset), each cohort's tissues from their water per
   ! plant. A leafless cohort's leaves hold none, so their unknown is the wood's potential less the
   ! head to the crown, where an empty leaf comes to rest. The pond's is its water's head, but no
   ! higher than the head at which it neither gives nor takes water from the node below it, the top
   ! porous layer or the top soil layer: over dry ground a trace pond empties within the step, and
   ! its own head would start Newton at a flow thousands of times what it holds.
   !---~---
   subroutine unknowns_from_water(static, profile, hyd, ice, porous, st)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(soil_ice_t),          intent(in)    :: ice(:)     ! (nsl) each soil layer's ice
      type(water_porous_t),      intent(in)    :: porous     ! the porous layers in the network
      type(water_state_t),       intent(inout) :: st
      integer(ik) :: k, i, j, soff
      real(wp)    :: psi, dpsi, w, dw, kk, dk, head1, elev
      soff = 1_ik + porous%n
      do j = 1_ik, porous%n
         st%tau_soil(1_ik + j) = soil_theta_to_unknown(porous%curve(j), porous%ice(j),                               &
                                                       (st%water_soil(1_ik + j) + porous%deficit(j))                 &
                                                       / (rho_h2o * porous%depth(j)))
      end do
      do k = 1_ik, profile%n_active
         st%tau_soil(soff + k) = soil_theta_to_unknown(profile%layer_curve(k), ice(k),                             &
                                                       st%water_soil(soff + k) / (rho_h2o * profile%dz(k)))
      end do
      if (porous%n > 0_ik) then
         call soil_water_relations(porous%curve(1), porous%ice(1), st%tau_soil(2), psi, dpsi, w, dw, kk, dk)
         elev = -0.5_wp * porous%depth(1)
      else
         call soil_water_relations(profile%layer_curve(1), ice(1), st%tau_soil(2), psi, dpsi, w, dw, kk, dk)
         elev = profile%z_node(1)
      end if
      head1 = psi / grav_head
      st%tau_soil(POND) = min(pond_water_to_head(st%water_soil(POND)), head1 + elev)
      do i = 1_ik, static%n
         associate (hp => hyd%pft(static%pft(i)))
            st%tau_coh(WOOD, i) = tissue_rwc_to_unknown(tissue_water_to_rwc(hp%wood_curve, static%bwood_water(i),   &
                                                                            st%water_coh(WOOD, i)                    &
                                                                            / max(static%nplant(i), tiny_num)))
            if (hp%leaf_curve%water_sat * static%bleaf(i) > tiny_num) then
               st%tau_coh(LEAF, i) = tissue_rwc_to_unknown(tissue_water_to_rwc(hp%leaf_curve, static%bleaf(i),      &
                                                                               st%water_coh(LEAF, i)                 &
                                                                               / max(static%nplant(i), tiny_num)))
            else
               call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i),                &
                                           st%tau_coh(WOOD, i), psi, dpsi, w, dw)
               st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(tissue_psi_to_rwc(psi - static%crown_head(i),         &
                                                             hp%leaf_curve%pi0, hp%leaf_curve%elastic_mod))
               st%water_coh(LEAF, i) = 0.0_wp
            end if
         end associate
      end do
   end subroutine unknowns_from_water

   !---~---
   !   water_derivs -- the network at one state: each store's potential and water from its
   ! curve; the flow along every path (pond to the top node, down the chain of porous layers and
   ! soil layers, the top porous and soil layers' evaporation, drainage, each cohort's uptake,
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
      integer(ik) :: k, i, ca, cb, c, nsl, nsoil, n_porous, soff
      real(wp) :: q, dq_pa, dq_pb, dq_ka, dq_kb, dqa, dqb, u, s, ds_w, ds_l
      real(wp) :: e, de, qbot, dsrc, e_vap(N_EVAP_MAX), de_vap(N_EVAP_MAX, N_EVAP_MAX)
      real(wp) :: psi_w, dpsi_w, u_sum, dqb_sum, du_dpsi(profile%n_active), du_dk(profile%n_active)
      real(wp) :: w_b, dist
      integer(ik) :: m
      logical  :: wet, is_cut(size(tau_soil))

      nsl   = profile%n_active
      n_porous = drive%porous%n
      soff  = 1_ik + n_porous
      nsoil = soff + nsl
      is_cut = .false.
      if (present(cut)) is_cut = cut(1:nsoil)

      ev%inflow_soil = 0.0_wp
      ev%inflow_coh = 0.0_wp
      !---~---
      !   The Jacobian's entries are written afresh: the core's band and the evaporating stores'
      ! block (nodes 2 .. soff + 1), which the flows add into, are cleared; every other entry an
      ! evaluation touches is set outright below, and the rest stay zero from arrowhead_allocate.
      !---~---
      m = n_porous + 1_ik
      if (with_jacobian) then
         do c = 1_ik, nsoil
            sys%a(c, c) = 0.0_wp
            if (c < nsoil) then
               sys%a(c, c + 1_ik) = 0.0_wp
               sys%a(c + 1_ik, c) = 0.0_wp
            end if
         end do
         sys%a(2:soff + 1_ik, 2:soff + 1_ik) = 0.0_wp
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

      !--- Down the chain: the pond to the top node, through the porous layers, between soil layers.
      ev%face = 0.0_wp
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
         ev%face(c) = q
         call add_soil_flow(ca, cb, q, dqa, dqb)
      end do

      !--- Each node's own sources: the water reaching the pond, and the drainage.
      do c = 1_ik, nsoil
         dsrc = 0.0_wp
         call node_sources(profile, soil, drive, c, ev%psi_soil(c), ev%dpsi_soil(c), ev%k_soil(c), ev%dk_soil(c),     &
                           ev%inflow_soil(c), dsrc, qbot)
         if (with_jacobian) sys%a(c, c) = sys%a(c, c) + dsrc
         if (c == nsoil) ev%drain = qbot
      end do

      !--- The porous layers' and the top soil layer's vapour, through the one leg to the canopy air.
      call network_vapour(profile, soil, drive, ev%psi_soil(2:soff + 1_ik), ev%dpsi_soil(2:soff + 1_ik),             &
                          ev%water_soil(2:soff + 1_ik), ev%dwater_soil(2:soff + 1_ik), e_vap, de_vap)
      ev%porous_evap = 0.0_wp
      ev%porous_evap(1:n_porous) = e_vap(1:n_porous)
      ev%evap = e_vap(m)
      ev%inflow_soil(2:soff + 1_ik) = ev%inflow_soil(2:soff + 1_ik) - e_vap(1:m)
      if (with_jacobian) sys%a(2:soff + 1_ik, 2:soff + 1_ik) = sys%a(2:soff + 1_ik, 2:soff + 1_ik) + de_vap(1:m, 1:m)

      !---~---
      !   Each cohort: uptake from its layers, sapflow to its leaves, and transpiration. The
      ! cohort's uptake comes from its whole root profile at once, each layer's conductance
      ! times its roots' contact with it (none once it freezes); its own sums are kept in
      ! scalars over its layers and written once. The roots take nothing from the porous layers,
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
   ! path lies in it), a porous layer, or a soil layer.
   !---~---
   pure subroutine core_node_relations(profile, drive, c, tau, psi, dpsi, w, dw, kc, dkc)
      type(soil_profile_t), intent(in)  :: profile
      type(water_drive_t),  intent(in)  :: drive
      integer(ik),          intent(in)  :: c
      real(wp),             intent(in)  :: tau
      real(wp),             intent(out) :: psi, dpsi, w, dw, kc, dkc
      integer(ik) :: soff
      soff = 1_ik + drive%porous%n
      if (c == POND) then
         call pond_water_relations(tau, psi, dpsi, w, dw)
         if (drive%porous%n > 0_ik) then
            kc = drive%porous%curve(1)%ksat * drive%porous%ice(1)%impedance
         else
            kc = profile%ksat(1) * drive%ice(1)%impedance
         end if
         dkc = 0.0_wp
      else if (c <= soff) then
         call soil_water_relations(drive%porous%curve(c - 1_ik), drive%porous%ice(c - 1_ik), tau, psi, dpsi, w, dw,   &
                                   kc, dkc)
         w = w - drive%porous%deficit(c - 1_ik)   ! a layer below its residual: its curve's water less the lack
      else
         call soil_water_relations(profile%layer_curve(c - soff), drive%ice(c - soff), tau, psi, dpsi, w, dw, kc, dkc)
      end if
   end subroutine core_node_relations

   !---~---
   !   Face c of the chain, joining node c to node c + 1: the weight of node c + 1's conductivity
   ! in the face's, the distance between the two nodes, which is also the drop in elevation, and
   ! whether it wets (core_face_flux): the pond's face onto the top porous layer alone. Every face
   ! below the pond takes the thickness-weighted geometric mean of its nodes' conductivities, as a
   ! soil face does, one conductance whichever way the water flows: between the porous layers, a
   ! frozen litter over a thawed moss draws water up through its own ice-impeded conductivity, not
   ! through the moss's, and gives it back the same way; into the top soil layer, the soil's
   ! conductivity limits what enters it. The pond counts as a saturated layer as thick as the node
   ! below it.
   !---~---
   pure subroutine core_face_geometry(profile, drive, c, w_b, dist, wet)
      type(soil_profile_t), intent(in)  :: profile
      type(water_drive_t),  intent(in)  :: drive
      integer(ik),          intent(in)  :: c
      real(wp),             intent(out) :: w_b, dist
      logical,              intent(out) :: wet
      integer(ik) :: n_porous, j
      n_porous = drive%porous%n
      if (c == 1_ik) then
         w_b = 0.5_wp
         if (n_porous > 0_ik) then
            dist = 0.5_wp * drive%porous%depth(1)
         else
            dist = -profile%z_node(1)
         end if
         wet = n_porous > 0_ik
      else if (c <= n_porous + 1_ik) then
         j = c - 1_ik
         if (j < n_porous) then
            w_b  = drive%porous%depth(j) / (drive%porous%depth(j) + drive%porous%depth(j + 1_ik))
            dist = 0.5_wp * (drive%porous%depth(j) + drive%porous%depth(j + 1_ik))
         else
            w_b  = drive%porous%depth(j) / (drive%porous%depth(j) + profile%dz(1))
            dist = 0.5_wp * drive%porous%depth(j) - profile%z_node(1)
         end if
         wet = .false.
      else
         j = c - 1_ik - n_porous
         w_b  = profile%face_weight(j)
         dist = profile%dz_node(j)
         wet  = .false.
      end if
   end subroutine core_face_geometry

   !---~---
   !   The flow down a face of the chain [kg/m2/s] and its slopes. A face takes the
   ! thickness-weighted geometric mean of its nodes' conductivities (soil_water_flux), but the
   ! pond's face onto a porous layer wets: water flowing down it takes the pond's conductivity, the
   ! saturated layer's, the upstream node's (Forsyth, Wu & Pruess 1995), so the pond's water
   ! enters a dry litter at its saturated conductivity rather than at the vanishing one of its dry
   ! pores; water flowing up takes the layer's own. The flow is the layer's conductivity times the
   ! drive, plus the difference of the two conductivities times the drive's rounded positive part
   ! (smooth_max0 over WET_BLEND_WIDTH): it rises with the drive at every drive, its slope between
   ! the two conductivities', so Newton meets neither a switch nor a fold.
   !---~---
   pure subroutine core_face_flux(wet, k_a, k_b, w_b, psi_a, psi_b, dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
      logical,  intent(in)  :: wet
      real(wp), intent(in)  :: k_a, k_b, w_b, psi_a, psi_b, dist
      real(wp), intent(out) :: q, dq_pa, dq_pb, dq_ka, dq_kb
      real(wp) :: drive, down, ddown
      if (.not. wet) then
         call soil_water_flux(k_a, k_b, w_b, psi_a, psi_b, dist, dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
         return
      end if
      drive = ((psi_a - psi_b) / grav_head + dist) / dist
      call smooth_max0(drive, WET_BLEND_WIDTH, down, ddown)    ! the drive's rounded downward part
      q      = rho_h2o * (k_b * drive + (k_a - k_b) * down)
      dq_pa  = rho_h2o * (k_b + (k_a - k_b) * ddown) / (grav_head * dist)
      dq_pb  = -dq_pa
      dq_ka  = rho_h2o * down
      dq_kb  = rho_h2o * (drive - down)
   end subroutine core_face_flux

   !---~---
   !   Core node c's own sources, besides its faces, its roots and its vapour, added to its net
   ! inflow [kg/m2/s], and their slope in its unknown added to dout, the slope of its outflow: the
   ! water reaching the pond, and the drainage from the bottom layer, which drain returns.
   !---~---
   pure subroutine node_sources(profile, soil, drive, c, psi, dpsi, kc, dkc, inflow, dout, drain)
      type(soil_profile_t), intent(in)    :: profile
      type(soil_params_t),  intent(in)    :: soil
      type(water_drive_t),  intent(in)    :: drive
      integer(ik),          intent(in)    :: c
      real(wp),             intent(in)    :: psi, dpsi, kc, dkc               ! its potential, conductivity, slopes
      real(wp),             intent(inout) :: inflow, dout                     ! [kg/m2/s], its slope in the unknown
      real(wp),             intent(out)   :: drain                            ! [kg/m2/s]
      real(wp)    :: qbot, dqbot_dpsi, dqbot_dk
      integer(ik) :: nsl
      nsl   = profile%n_active
      drain = 0.0_wp
      if (c == POND) inflow = inflow + drive%rain
      if (c == 1_ik + drive%porous%n + nsl) then
         call bottom_flux(profile, soil, nsl, psi, kc, drive%ice(nsl)%impedance, qbot, dqbot_dpsi, dqbot_dk)
         drain  = rho_h2o * qbot
         inflow = inflow - drain
         dout   = dout + rho_h2o * (dqbot_dpsi * dpsi + dqbot_dk * dkc)
      end if
   end subroutine node_sources

   !---~---
   !   network_vapour -- the vapour of the network's evaporating stores at one state, the porous
   ! layers (nodes 2 .. soff) and the top soil layer (node soff + 1), top down, given each one's
   ! potential and water with their slopes in its unknown: each store's path to the ground's
   ! surface node, through its own dry surface layer (dry_layer_resistance), which its own water
   ! sets, and the open pores of the porous layers above it (store_vapour_path), and their
   ! evaporation through the one leg to the canopy air (ground_vapour), with its slopes in every
   ! evaporating store's unknown, de(i, k) store i's in store k's. A store's dry surface layer
   ! reads its liquid and dsl_theta_init of its ice (a porous layer's dsl_theta_init being one),
   ! so its onset is dsl_theta_init of the pores the ice leaves (CLM5's effective porosity).
   !---~---
   pure subroutine network_vapour(profile, soil, drive, psi, dpsi, w, dw, e, de)
      type(soil_profile_t), intent(in)  :: profile
      type(soil_params_t),  intent(in)  :: soil
      type(water_drive_t),  intent(in)  :: drive
      real(wp),             intent(in)  :: psi(:), dpsi(:), w(:), dw(:)   ! (porous%n + 1) [MPa], [kg/m2] and slopes
      real(wp),             intent(out) :: e(:), de(:,:)                  ! [kg/m2/s], its slopes
      type(vapour_path_t) :: path(N_EVAP_MAX)
      real(wp)    :: r_dry, dr_dry
      integer(ik) :: j, m
      m = drive%porous%n + 1_ik
      do j = 1_ik, drive%porous%n
         call dry_layer_resistance(w(j) / (rho_h2o * drive%porous%depth(j)), drive%porous%depth(j),                  &
                                   drive%porous%air(j), r_dry, dr_dry)
         path(j) = store_vapour_path(drive%porous%air(j), psi(j), dpsi(j), w(j), dw(j), drive%porous%depth(j),        &
                                     drive%porous%curve(j)%theta_res, drive%porous%curve(j)%theta_sat, r_dry, dr_dry)
      end do
      call dry_layer_resistance(w(m) / (rho_h2o * profile%dz(1)) - (1.0_wp - soil%dsl_theta_init) * drive%ice(1)%theta_ice, &
                                soil%dsl_dmax, drive%ground, r_dry, dr_dry)
      path(m) = store_vapour_path(drive%ground, psi(m), dpsi(m), w(m), dw(m), profile%dz(1), profile%theta_res(1),      &
                                  profile%theta_sat(1), r_dry, dr_dry)
      call ground_vapour(m, path(1:m), drive%ground, e(1:m), de(1:m, 1:m))
   end subroutine network_vapour

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
   ! below (none to or from a detached node, as water_derivs cuts them), and its own sources --
   ! the water reaching the pond, the evaporation of a porous layer or the top soil layer, the
   ! other evaporating stores held (network_vapour), the bottom's drainage and each cohort's uptake
   ! at its wood's potential; and, when
   ! asked, dR/dtau of the node above (dr_up) and below (dr_down), which a block of eliminated
   ! nodes solves with. The nonlinear elimination of a fast node or block and the closed form of a
   ! detached one solve R = 0 on it.
   !---~---
   subroutine node_balance(static, profile, soil, drive, c, tau_soil, w0, h, cut, psi_wood, r, dr, dr_up, dr_down)
      type(fast_static_state_t), intent(in)  :: static
      type(soil_profile_t),      intent(in)  :: profile
      type(soil_params_t),       intent(in)  :: soil
      type(water_drive_t),       intent(in)  :: drive
      integer(ik),               intent(in)  :: c
      real(wp),                  intent(in)  :: tau_soil(:)  ! (nsoil) the node's trial unknown among the others'
      real(wp),                  intent(in)  :: w0, h        ! [kg/m2] its water at the sub-step's start; [s]
      logical,                   intent(in)  :: cut(:)       ! (nsoil) the detached nodes: no faces, no roots
      real(wp),                  intent(in)  :: psi_wood(:)  ! (ncoh) [MPa]
      real(wp),                  intent(out) :: r, dr
      real(wp), optional,        intent(out) :: dr_up, dr_down
      real(wp) :: psi, dpsi, w, dw, kc, dkc, q, dq_pa, dq_pb, dq_ka, dq_kb, w_b, dist, inflow, dinflow
      real(wp) :: qbot, dsrc, d_up, d_down
      real(wp) :: u, du_ds, du_dw, du_dk, above(6), below(6)
      real(wp) :: vp(N_EVAP_MAX), vdp(N_EVAP_MAX), vw(N_EVAP_MAX), vdw(N_EVAP_MAX), kdum, dkdum
      real(wp) :: e_vap(N_EVAP_MAX), de_vap(N_EVAP_MAX, N_EVAP_MAX)
      integer(ik) :: nsl, soff, nsoil, k, i, m, j
      logical :: wet
      nsl   = profile%n_active
      soff  = 1_ik + drive%porous%n
      nsoil = soff + nsl
      call core_node_relations(profile, drive, c, tau_soil(c), psi, dpsi, w, dw, kc, dkc)
      inflow  = 0.0_wp
      dinflow = 0.0_wp
      d_up    = 0.0_wp
      d_down  = 0.0_wp
      !--- The faces to the neighbours, held at their unknowns: (potential, its slope, water, its slope, K, its slope).
      if (c > 1_ik) then
         if (.not. (cut(c) .or. cut(c - 1_ik))) then
            call core_node_relations(profile, drive, c - 1_ik, tau_soil(c - 1_ik), above(1), above(2), above(3),     &
                                     above(4), above(5), above(6))
            call core_face_geometry(profile, drive, c - 1_ik, w_b, dist, wet)
            call core_face_flux(wet, above(5), kc, w_b, above(1), psi, dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
            inflow  = inflow + q
            dinflow = dinflow + dq_pb * dpsi + dq_kb * dkc
            d_up    = -(dq_pa * above(2) + dq_ka * above(6))
         end if
      end if
      if (c < nsoil) then
         if (.not. (cut(c) .or. cut(c + 1_ik))) then
            call core_node_relations(profile, drive, c + 1_ik, tau_soil(c + 1_ik), below(1), below(2), below(3),     &
                                     below(4), below(5), below(6))
            call core_face_geometry(profile, drive, c, w_b, dist, wet)
            call core_face_flux(wet, kc, below(5), w_b, psi, below(1), dist, q, dq_pa, dq_pb, dq_ka, dq_kb)
            inflow  = inflow - q
            dinflow = dinflow - (dq_pa * dpsi + dq_ka * dkc)
            d_down  = dq_pb * below(2) + dq_kb * below(6)
         end if
      end if
      dsrc = 0.0_wp
      call node_sources(profile, soil, drive, c, psi, dpsi, kc, dkc, inflow, dsrc, qbot)
      !--- An evaporating store's vapour, the others held: its slopes in its neighbours' unknowns join theirs.
      if (c > 1_ik .and. c <= soff + 1_ik) then
         m = soff
         do j = 1_ik, m
            if (j + 1_ik == c) then
               vp(j)  = psi
               vdp(j) = dpsi
               vw(j)  = w
               vdw(j) = dw
            else
               call core_node_relations(profile, drive, j + 1_ik, tau_soil(j + 1_ik), vp(j), vdp(j), vw(j), vdw(j),    &
                                        kdum, dkdum)
            end if
         end do
         call network_vapour(profile, soil, drive, vp(1:m), vdp(1:m), vw(1:m), vdw(1:m), e_vap, de_vap)
         j = c - 1_ik
         inflow = inflow - e_vap(j)
         dsrc   = dsrc + de_vap(j, j)
         if (j > 1_ik) d_up   = d_up + de_vap(j, j - 1_ik)
         if (j < m)    d_down = d_down + de_vap(j, j + 1_ik)
      end if
      dinflow = dinflow - dsrc
      if (c > soff .and. .not. cut(c)) then
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
      if (present(dr_up))   dr_up   = d_up
      if (present(dr_down)) dr_down = d_down
   end subroutine node_balance

end module ed_fast_integ_derivs
