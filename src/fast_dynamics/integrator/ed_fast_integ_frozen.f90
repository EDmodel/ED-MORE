! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_frozen -- the frozen record one fast step marches against: what is
! computed once at the start of the dt_fast, from that moment's state and forcing, and
! held fixed across the march (docs/science/numerical_scheme.md §2, Category 0).
! build_frozen runs it as named steps:
!   * the ground's stage: the stores the step starts from, the canopy's interception of the
!     rain, the porous layers' spill and the water layer's stage; then the radiation through
!     the canopy onto that ground (ed_fast_integ_radiation);
!   * the processes' rates at the start state: canopy aerodynamics, the root zone, leaf
!     gas exchange, stem and fine-root maintenance respiration, heterotrophic respiration,
!     the biotic CO2 source, and the canopy-air capacities and atmosphere conductances;
!   * the tissue coefficients, the rest of the canopy-air and ground boundary, the
!     transpiration demand and the ground inflow;
!   * the soil-plant water network's start state, and the march's.
! The water network is stepped next (water_step), and record_water_flows turns what it
! moved into the record's rates (docs/science/soil_plant_water.md §4-§5). Both schemes
! march against the same record, so they see the same GPP.
!---~---
module ed_fast_integ_frozen
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : tiny_num, rho_h2o, pi, r_gas, cp_organic,                            &
                                     kgcday_2_umols
   use ed_diag_vars,          only : CO_ANET, CO_AGROSS, CO_GSW, CO_GBW, CO_CI, CO_CS, CO_RD, CO_BETA_STOM,     &
                                     CO_BETA_NONSTOM, CO_LEAF_TEMP, CO_WOOD_TEMP, CO_LEAF_VPD, CO_PSI_LEAF,      &
                                     CO_ABS_PAR, CO_ABS_SW, CO_ABS_LW, CO_WIND, CO_LEAF_WATER, CO_WOOD_WATER,    &
                                     CO_GPP_RATE, CO_LEAF_PAR, CO_GB_MOL, CO_PSI_PREDAWN,    &
                                     CO_PSI_WOOD, CO_PLC, CO_SAPFLOW, CO_ROOT_UPTAKE, CO_WET_FRACTION
   use ed_plant_hydr_lib,     only : conductance_retained_fraction, root_soil_contact
   use ed_plant_hydraulics,   only : leaf_node_psi
   use ed_soil_derive,        only : top_soil_psi
   use ed_soil_params,        only : soil_profile_t
   use ed_soil_hydr_lib,      only : soil_ice_state
   use ed_config,             only : ed_config_t
   use ed_site_state_vars,    only : site_surface_t
   use ed_canopy_aerodynamics, only : aero_env_t, aero_args_t, aero_out_t, set_aero_env_atm, set_aero_env_canopy
   use ed_pft_params,       only : leaf_photo_table_t, hydraulics_table_t
   use ed_fast_integ_types,       only : fast_dynamic_state_t, fast_static_state_t, fast_forcing_t, fast_output_t,     &
                                     march_state_t, frozen_t, surface_state_t,                            &
                                     surface_tend_t, cas_boundary_t,                                              &
                                     plant_water_t, step_exports_t,                                              &
                                     water_work_t, POND, WOOD, WORK_WATER_NFAIL
   use ed_fast_integ_derivs, only : canopy_demand, set_water_state
   use ed_state_defaults,   only : is_unset_negative
   use ed_biophysics_params, only : aero_params_t
   use ed_patch_state_vars, only : cas_state_t, soil_carbon_state_t, cas_capacities, cas_capacity_change,          &
                                   porous_layer_state_t, porous_layer_curve, porous_layers_diagnose,              &
                                   porous_layers_present, porous_layers_spill, ground_top_temp, ground_surface_store
   use ed_ground_surface_params, only : ground_surface_layer_params_t, porous_material_t, N_POROUS
   use ed_biogeochem_params, only : N_SOIL_POOL
   use ed_canopy_aerodynamics, only : canopy_aerodynamics, cas_atm_conductances
   use ed_plant_energy,   only : plant_energy_env_t, sensible_heat_coeff, leaf_transp_coeff, leaf_surface_evap_coeff, &
                                 tissue_hcap_per_dt, EFFAREA_HEAT, EFFAREA_TRANSP, EFFAREA_EVAP
   use ed_canopy_water,   only : canopy_water_out_t, canopy_interception, limit_surface_evaporation
   use ed_leaf_gas_exchange, only : leaf_gas_exchange_batch, cohort_leaf_par
   use ed_respiration,   only : stem_maintenance_respiration, fine_root_maintenance_respiration, root_zone_temp_scale, &
                                decomposition_env_scalar, decomposition_transfer_matrix, heterotrophic_respiration
   use ed_surface, only : surface_env_t, porous_node_env_t, water_layer_stage, water_layer_cover_of, water_layer_exposed
   use ed_soil_water, only : porous_cover_path, slot_evaporation_setup
   use ed_therm_lib,          only : internal_energy_liquid,                                                     &
                                     cas_enthalpy_to_temp, cas_heat_capacity, air_vpd, air_density,             &
                                     soil_thermal_conductivity, specific_humidity_to_vapor_pressure
   use ed_numerics,           only : weighted_mean
   use ed_ground_surface_lib, only : porous_layer_fill, porous_layer_depth, porous_layer_capacity,              &
                                     porous_layer_resistance
   use ed_fast_integ_radiation, only : canopy_ground_radiation
   implicit none
   private

   public :: build_frozen, record_water_flows

contains

   !---~---
   !   build_frozen -- the frozen record both schemes march against, as named steps, and the
   ! water network's and the march's start states. It only reads the dynamic state: the
   ! ground's stage advances a copy of its stores, which the update commits. With the sky
   ! live, the radiation through the canopy fills the forcing's absorbed radiation after the
   ! stage. The record is the thread's, reused step to step (prepare_frozen). out receives the
   ! aerodynamics, the per-cohort rates, the step's GPP, NEE and decomposition scalars, and,
   ! with out%capture, the per-cohort diagnostics (docs/science/diagnostics.md §1, stage [2]).
   !---~---
   subroutine build_frozen(dt_fast, cfg, site, static, forcing, dynamic, out, water, frozen, y)
      real(wp),                   intent(in)    :: dt_fast
      type(ed_config_t),          intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(inout) :: forcing
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(fast_output_t),        intent(inout) :: out
      type(water_work_t),         intent(inout) :: water    ! the soil-plant water network
      type(frozen_t),             intent(inout) :: frozen
      type(march_state_t),       intent(out)   :: y
      !---~---
      !   The canopy air's temperature and humidity at the start state, its pressure, density and
      ! specific heat, and the ground's temperature, which the later steps read.
      !---~---
      real(wp)    :: tcas, qcas, press, rho, cp_cas, t_ground
      real(wp)    :: can_dmol      ! [mol/m3] the canopy air's dry-air molar density
      !--- The processes' patch totals [umol/m2 ground/s] and the root zone's drivers.
      real(wp)    :: gpp, ra_leaf, ra_stem, ra_root, rh, soil_temp_root, theta_mean
      real(wp)    :: t_surface        ! [K] the ground's surface: the top porous layer's, else the top soil node's
      real(wp)    :: t_soil_top       ! [K] the top soil node's
      integer(ik) :: n, nsl

      n   = static%n
      nsl = site%soil_profile%n_active
      call prepare_frozen(n, nsl, frozen)

      !---~---
      !   The ground's stage, then the radiation through the canopy onto the ground it leaves, so
      ! the two-stream's ground and the surface's balance have one cover.
      !---~---
      call ground_stage(dt_fast, cfg, site, static, forcing, dynamic, out, frozen)
      if (forcing%sky_live) call canopy_ground_radiation(cfg, static, dynamic, frozen%ground%exposed,           &
                                                         frozen%ground%snowfac, forcing)

      !---~---
      !   The processes' rates at the start state. The ground's surface is the top porous layer
      ! where there is one (docs/science/surface_organic_layer.md): the aerodynamics and the
      ! litter's decomposition read its temperature, the soil's own evaporation the top soil node's.
      !---~---
      t_soil_top = dynamic%soil_energy%soil_temp(1)
      t_surface  = ground_top_temp(dynamic%ground_surface, t_soil_top)
      call refresh_canopy_aerodynamics(cfg%aero, frozen%cas%env, static, forcing, dynamic%cas,               &
                                       t_surface, dynamic%leaf_temp, dynamic%wood_temp,                       &
                                       out%aero, tcas, qcas, press, rho, cp_cas, t_ground)
      call root_zone_environment(site%soil_profile, static%root_frac_patch(1:nsl), dynamic%soil_energy%soil_temp, &
                                 dynamic%soil_water%theta, soil_temp_root, theta_mean)
      call canopy_leaf_gas_exchange(static%leaf_photo, cfg%pft%hydraulics, site%soil_profile, static, forcing,     &
                                    out%aero, dynamic, tcas, qcas, press, rho, cp_cas, gpp, ra_leaf,             &
                                    frozen%tissue%h_coeff_leaf, frozen%tissue%g_transp_leaf, out%gpp_coh,        &
                                    out%leaf_resp_coh, out%capture, out%cdiag)
      call canopy_maintenance_respiration(cfg%leaf_physiology%ea_rd, cfg%leaf_physiology%hd_rd,                   &
                                          cfg%leaf_physiology%ds_rd, site%soil_profile,                           &
                                          static%root_frac_patch(1:nsl), static, dynamic%wood_temp,               &
                                          dynamic%soil_energy%soil_temp, ra_stem, ra_root, out%stem_resp_coh,     &
                                          out%root_resp_coh)
      if (porous_layers_present(dynamic%ground_surface%porous)) then
         call patch_heterotrophic_respiration(cfg, site%soil_profile, static%soil_carbon, t_ground, soil_temp_root, &
                                              theta_mean, rh, out%report%exports,                                 &
                                              surf_rel=porous_layers_wetness(dynamic%ground_surface%porous,        &
                                                                             cfg%ground_surface))
      else
         call patch_heterotrophic_respiration(cfg, site%soil_profile, static%soil_carbon, t_ground, soil_temp_root, &
                                              theta_mean, rh, out%report%exports)
      end if
      !---~---
      !   The biotic CO2 source (docs/science/plant_respiration.md eq. (7)), positive to the canopy
      ! air. slow_co2_rate is the autotrophic part this loop cannot compute, the allocator's daily
      ! growth respiration on realized growth, handed down as a frozen rate with the
      ! opposite-signed starvation correction (§5-§6 there).
      !---~---
      frozen%cas%nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp + static%slow_co2_rate
      out%report%exports%gpp_last = gpp
      out%report%exports%nee_last = frozen%cas%nee_biotic
      !---~---
      !   The canopy air's capacities, and its exchange with the atmosphere at the start state
      ! (docs/science/canopy_aerodynamics.md eq. (6)): heat on heat_transfer_factor, vapour and
      ! CO2 on vapour_transfer_factor (equal, z0q = z0h).
      !---~---
      call cas_capacities(rho, qcas, dynamic%cas%can_depth, can_dmol, frozen%cas%box%air_mass_capacity,         &
                          frozen%cas%box%air_molar_capacity)
      call cas_atm_conductances(rho, can_dmol, out%aero%ustar, out%aero%heat_transfer_factor,                  &
                                out%aero%vapour_transfer_factor, frozen%cas%box%atm_conductance_enthalpy,     &
                                frozen%cas%box%atm_conductance_vapor, frozen%cas%box%atm_conductance_co2)
      !--- The content change from the capacities moving since the last step: a budget term.
      call cas_capacity_change(dynamic%cas, frozen%cas%box%air_mass_capacity, frozen%cas%box%air_molar_capacity, &
                               out%report%exports%cas_capacity_energy, out%report%exports%cas_capacity_water,    &
                               out%report%exports%cas_capacity_co2)

      !---~---
      !   The thermal resistance from the ground's surface down to the top soil node, its depth over
      ! the top layer's conductivity, which the stack's bottom node conducts through.
      !---~---
      associate (prof => site%soil_profile)
         frozen%ground%r_soil_top = max(-prof%z_node(1), tiny_num)                                                &
                                    / max(soil_thermal_conductivity(dynamic%soil_water%theta(1),                  &
                                                                    dynamic%soil_energy%soil_fliq(1),             &
                                                                    prof%theta_sat(1), prof%solid_conductivity(1), &
                                                                    prof%dry_conductivity(1)), tiny_num)
      end associate

      call tissue_coefficients(dt_fast, cfg, static, forcing, dynamic, out%aero, rho, cp_cas, frozen%tissue,          &
                               frozen%veg_surface, frozen%plant)
      call buried_tissue_store(frozen%tissue, static%n, dt_fast, frozen%ground%buried_cap, frozen%ground%buried_energy)
      call canopy_air_boundary(cfg%aero, static%geom, out%aero, forcing, rho, cp_cas, press, frozen%cas, frozen%ground)
      call ground_inflow(site%soil_profile, forcing, dynamic, out%aero%ggnet, tcas, qcas, rho, &
                         t_soil_top, nsl, cfg%ground_surface, frozen)
      frozen%water_drive%pond_capacity = cfg%ground_surface%water%pond_capacity
      frozen%ground%ground_emiss = cfg%soil%ground_emissivity
      call transpiration_demand(dt_fast, cfg%canopy_water_on, n, dynamic, frozen)
      call porous_network_drive(cfg%ground_surface, frozen)

      !--- The start states: the water network's, from the stores' state^n water and an empty pond, and the march's.
      call set_water_state(static, site%soil_profile, cfg%pft%hydraulics, dynamic%soil_water%theta(1:nsl),       &
                           frozen%water_drive%ice(1:nsl), 0.0_wp,                                                &
                           dynamic%leaf_water_mass(1:n), dynamic%wood_water_mass(1:n), water,                    &
                           frozen%water_drive%porous)
      call pack_march_start(dynamic, site%soil_profile%dz(1), frozen%stage_soil_heat, n, nsl, y)
   end subroutine build_frozen

   !---~---
   !   ground_stage -- the ground's stores before the step's radiation and its march: their water
   ! and energy as the step finds them, which the whole-column ledgers start from; the canopy's
   ! interception of the rain, which the ground does not change; the porous layers' spill, what a
   ! layer holds beyond its saturation joining the water layer; and the water layer's stage
   ! (water_layer_stage): the rain through the canopy, all the snowfall, and the water the
   ! tissues shed join it, its free water is set apart for the water network, a film takes the
   ! temperature of the node it lies on, the top porous layer or the top soil layer, and its
   ! cover is set. The water layer's exchange with the canopy air and the soil is the march's.
   !---~---
   subroutine ground_stage(dt_fast, cfg, site, static, forcing, dynamic, out, frozen)
      real(wp),                   intent(in)    :: dt_fast
      type(ed_config_t),          intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(fast_output_t),        intent(inout) :: out
      type(frozen_t),             intent(inout) :: frozen
      real(wp)    :: rain_through   ! [kg/m2/s] the rain through the canopy
      integer(ik) :: k
      frozen%surface0 = ground_surface_store(dynamic%ground_surface)
      call intercept_rain(dt_fast, cfg, static, forcing, dynamic, frozen%veg_surface, rain_through)
      frozen%ground%pack   = dynamic%ground_surface%water
      frozen%ground%params = cfg%ground_surface%water
      frozen%ground%dt     = dt_fast
      frozen%porous = dynamic%ground_surface%porous
      call porous_layers_spill(frozen%porous, cfg%ground_surface, frozen%ground%pack%ground_store_t)
      k = top_porous(frozen%porous)
      if (k > 0_ik) then
         associate (top => frozen%porous(k))
            call water_layer_stage(rain_through, forcing%snowfall, forcing%air_temp, static%shed_water_rate,        &
                                   static%shed_water_enth, dt_fast, cfg%ground_surface%water,                      &
                                   cfg%ground_surface%min_depth, top%temp, top%fliq, frozen%ground%pack,            &
                                   frozen%stage, below=top%ground_store_t, below_dry_cap=cp_organic * top%matter, &
                                   below_depth=porous_layer_depth(top%matter, cfg%ground_surface%porous(k)%bulk_density), &
                                   below_curve=porous_layer_curve(cfg%ground_surface%porous(k), top))
            !--- A film's heat comes from the top porous layer.
            top%energy = top%energy - frozen%stage%film_heat
         end associate
         call porous_layers_diagnose(frozen%porous, cfg%ground_surface)
         frozen%stage_soil_heat = 0.0_wp
      else
         call water_layer_stage(rain_through, forcing%snowfall, forcing%air_temp, static%shed_water_rate,           &
                                static%shed_water_enth, dt_fast, cfg%ground_surface%water,                         &
                                cfg%ground_surface%min_depth, dynamic%soil_energy%soil_temp(1),                     &
                                dynamic%soil_energy%soil_fliq(1), frozen%ground%pack, frozen%stage)
         !--- A film's heat comes from the top soil layer, at the march's start.
         frozen%stage_soil_heat = frozen%stage%film_heat
      end if
      !---~---
      !   The layer exposed, what it holds and its free water together, sets its cover of the
      ! ground and its optics (water_layer_exposed).
      !---~---
      frozen%ground%exposed = water_layer_exposed(frozen%ground%pack, frozen%stage%free)
      frozen%ground%snowfac = water_layer_cover_of(frozen%ground%exposed, frozen%stage%is_layer,                  &
                                                   cfg%ground_surface%water, cfg%aero%soil_rough)
      !---~---
      !   The free water reaches the water network's pond over the step at a steady rate, so the
      ! pond starts from empty at the head the ground below gives it. The network moves its water;
      ! its heat is the water layer's store's.
      !---~---
      frozen%water_drive%rain = frozen%stage%free%water / dt_fast
      frozen%ground%free_temp  = frozen%stage%free_temp
      frozen%ground%free_start = frozen%stage%free
      out%report%exports%water_layer_cover = frozen%ground%snowfac
   end subroutine ground_stage

   !---~---
   !   The frozen record for a step, reused from the last: every scalar back at its default, as
   ! a fresh record has it, and the per-cohort and per-layer arrays kept when their sizes still
   ! fit. The plant water's heat terms start at zero: the start-of-step surface solves read them
   ! (as q_extra, whose default is 0) before the water network's flows give them their values.
   !---~---
   subroutine prepare_frozen(n, nsl, frozen)
      integer(ik),    intent(in)    :: n, nsl
      type(frozen_t), intent(inout) :: frozen
      type(frozen_t) :: defaults   ! default-initialised and never written: the defaults
      type(frozen_t) :: arrays     ! holds the arrays while the scalars are reset
      call move_frozen_arrays(frozen, arrays)
      frozen = defaults
      call move_frozen_arrays(arrays, frozen)
      call fit(frozen%tissue%h_coeff_leaf, n)
      call fit(frozen%tissue%g_transp_leaf, n)
      call fit(frozen%tissue%abs_sw, n)
      call fit(frozen%tissue%abs_lw, n)
      call fit(frozen%tissue%lai, n)
      call fit(frozen%tissue%h_coeff_w, n)
      call fit(frozen%tissue%abs_sw_wood, n)
      call fit(frozen%tissue%abs_lw_wood, n)
      call fit(frozen%tissue%wai, n)
      call fit(frozen%tissue%wood_dry_hcap, n)
      call fit(frozen%tissue%wood_wmass, n)
      call fit(frozen%tissue%leaf_dry_hcap, n)
      call fit(frozen%tissue%leaf_wmass, n)
      call fit(frozen%tissue%buried, n)
      call fit(frozen%tissue%leaf_hcap_per_dt, n)
      call fit(frozen%tissue%wood_hcap_per_dt, n)
      call fit(frozen%tissue%t_leaf0, n)
      call fit(frozen%tissue%t_wood0, n)
      call fit(frozen%tissue%qwflux_wl, n)
      call fit(frozen%tissue%q_wood_net, n)
      call fit(frozen%tissue%leaf_water_enth, n)
      call fit(frozen%tissue%wood_water_enth, n)
      call fit(frozen%tissue%leaf_emiss, n)
      call fit(frozen%tissue%wood_emiss, n)
      call fit(frozen%tissue%transp_net, n)
      call fit(frozen%veg_surface%f_wet_c, n)
      call fit(frozen%veg_surface%g_surf_leaf, n)
      call fit(frozen%veg_surface%g_surf_wood, n)
      call fit(frozen%veg_surface%intercept_leaf, n)
      call fit(frozen%veg_surface%intercept_wood, n)
      call fit(frozen%veg_surface%transfer_lw, n)
      call fit(frozen%roots%root_share, nsl)
      call fit(frozen%roots%qloss_frozen, n)
      call fit(frozen%roots%root_heat, nsl)
      call fit(frozen%plant%nplant, n)
      call fit(frozen%plant%sapflow_frozen, n)
      call fit(frozen%plant%uptake_frozen, n)
      call fit(frozen%hydrology%theta1, nsl)
      call fit(frozen%water_drive%e_demand, n)
      frozen%tissue%qwflux_wl(1:n)  = 0.0_wp
      frozen%tissue%q_wood_net(1:n) = 0.0_wp
   end subroutine prepare_frozen

   !---~---
   !   Move every array of the frozen record from one record to another (the list prepare_frozen
   ! sizes). An array this list missed would be freed by the reset and allocated afresh.
   !---~---
   subroutine move_frozen_arrays(from, to)
      type(frozen_t), intent(inout) :: from, to
      call move_alloc(from%tissue%h_coeff_leaf, to%tissue%h_coeff_leaf)
      call move_alloc(from%tissue%g_transp_leaf, to%tissue%g_transp_leaf)
      call move_alloc(from%tissue%abs_sw, to%tissue%abs_sw)
      call move_alloc(from%tissue%abs_lw, to%tissue%abs_lw)
      call move_alloc(from%tissue%lai, to%tissue%lai)
      call move_alloc(from%tissue%h_coeff_w, to%tissue%h_coeff_w)
      call move_alloc(from%tissue%abs_sw_wood, to%tissue%abs_sw_wood)
      call move_alloc(from%tissue%abs_lw_wood, to%tissue%abs_lw_wood)
      call move_alloc(from%tissue%wai, to%tissue%wai)
      call move_alloc(from%tissue%wood_dry_hcap, to%tissue%wood_dry_hcap)
      call move_alloc(from%tissue%wood_wmass, to%tissue%wood_wmass)
      call move_alloc(from%tissue%leaf_dry_hcap, to%tissue%leaf_dry_hcap)
      call move_alloc(from%tissue%leaf_wmass, to%tissue%leaf_wmass)
      call move_alloc(from%tissue%buried, to%tissue%buried)
      call move_alloc(from%tissue%leaf_hcap_per_dt, to%tissue%leaf_hcap_per_dt)
      call move_alloc(from%tissue%wood_hcap_per_dt, to%tissue%wood_hcap_per_dt)
      call move_alloc(from%tissue%t_leaf0, to%tissue%t_leaf0)
      call move_alloc(from%tissue%t_wood0, to%tissue%t_wood0)
      call move_alloc(from%tissue%qwflux_wl, to%tissue%qwflux_wl)
      call move_alloc(from%tissue%q_wood_net, to%tissue%q_wood_net)
      call move_alloc(from%tissue%leaf_water_enth, to%tissue%leaf_water_enth)
      call move_alloc(from%tissue%wood_water_enth, to%tissue%wood_water_enth)
      call move_alloc(from%tissue%leaf_emiss, to%tissue%leaf_emiss)
      call move_alloc(from%tissue%wood_emiss, to%tissue%wood_emiss)
      call move_alloc(from%tissue%transp_net, to%tissue%transp_net)
      call move_alloc(from%veg_surface%f_wet_c, to%veg_surface%f_wet_c)
      call move_alloc(from%veg_surface%g_surf_leaf, to%veg_surface%g_surf_leaf)
      call move_alloc(from%veg_surface%g_surf_wood, to%veg_surface%g_surf_wood)
      call move_alloc(from%veg_surface%intercept_leaf, to%veg_surface%intercept_leaf)
      call move_alloc(from%veg_surface%intercept_wood, to%veg_surface%intercept_wood)
      call move_alloc(from%veg_surface%transfer_lw, to%veg_surface%transfer_lw)
      call move_alloc(from%roots%root_share, to%roots%root_share)
      call move_alloc(from%roots%qloss_frozen, to%roots%qloss_frozen)
      call move_alloc(from%roots%root_heat, to%roots%root_heat)
      call move_alloc(from%plant%nplant, to%plant%nplant)
      call move_alloc(from%plant%sapflow_frozen, to%plant%sapflow_frozen)
      call move_alloc(from%plant%uptake_frozen, to%plant%uptake_frozen)
      call move_alloc(from%hydrology%theta1, to%hydrology%theta1)
      call move_alloc(from%water_drive%e_demand, to%water_drive%e_demand)
   end subroutine move_frozen_arrays

   !--- Make a an array of m elements, keeping it when it already is one.
   pure subroutine fit(a, m)
      real(wp), allocatable, intent(inout) :: a(:)
      integer(ik),           intent(in)    :: m
      if (allocated(a)) then
         if (size(a, kind=ik) == m) return
         deallocate(a)
      end if
      allocate(a(m))
   end subroutine fit

   !---~---
   !   refresh_canopy_aerodynamics -- the aerodynamic environment: its atmosphere from the
   ! forcing at the canopy-air top, the reference level the surface layer is solved to, and its
   ! canopy air diagnosed from the CAS state; then the whole-canopy aerodynamics solve (friction
   ! velocity, scalar profile factors, ground conductance, per-cohort boundary layers;
   ! docs/science/canopy_aerodynamics.md §2-§4). The canopy air's density and specific heat
   ! are held for the step, at its start state.
   !---~---
   subroutine refresh_canopy_aerodynamics(aero_cfg, aenv, static, forcing, cas, soil_temp_top,    &
                                          leaf_temp, wood_temp, aero, tcas, qcas, press, rho, cp_cas, t_ground)
      type(aero_params_t),       intent(in)    :: aero_cfg
      type(aero_env_t),          intent(inout) :: aenv
      type(fast_static_state_t), intent(in)    :: static     ! the cohorts, the canopy geometry and its top
      type(fast_forcing_t),      intent(in)    :: forcing
      type(cas_state_t),         intent(in)    :: cas
      real(wp),                  intent(in)    :: soil_temp_top    ! [K] top soil-node temperature (ground skin)
      real(wp),                  intent(in)    :: leaf_temp(:)     ! [K] per-cohort leaf temperature
      real(wp),                  intent(in)    :: wood_temp(:)     ! [K] per-cohort wood temperature
      type(aero_out_t),          intent(inout) :: aero
      real(wp),                  intent(out)   :: tcas, qcas, press, rho, cp_cas, t_ground
      call set_aero_env_atm(aenv, forcing%wind, static%z_top, forcing%air_temp, forcing%shv_atm)
      tcas = cas_enthalpy_to_temp(cas%can_enthalpy, cas%can_shv)
      qcas = cas%can_shv
      press = forcing%press
      !---~---
      !   The canopy air's own density, its temperature and humidity at the pressure of its top,
      ! and its specific heat at that humidity, the c_p of every sensible flux it exchanges.
      !---~---
      rho = air_density(tcas, press, qcas)
      cp_cas = cas_heat_capacity(qcas)
      t_ground = soil_temp_top
      call set_aero_env_canopy(aenv, tcas, qcas, t_ground)
      call canopy_aerodynamics(aero_cfg, aenv, static%geom, static%n, static%order, static%height, static%lai, &
                               leaf_temp, wood_temp, static%leaf_width, static%branch_diam, aero)
   end subroutine refresh_canopy_aerodynamics

   !---~---
   !   root_zone_environment -- root-weighted soil temperature and the column-mean soil moisture
   ! that drive fine-root and heterotrophic respiration.
   !---~---
   pure subroutine root_zone_environment(soil, root_frac, soil_temp, theta, soil_temp_root, theta_mean)
      type(soil_profile_t), intent(in)  :: soil
      real(wp),             intent(in)  :: root_frac(:)      ! [-]     the patch's fine roots by layer
      real(wp),             intent(in)  :: soil_temp(:)      ! [K]     per-layer soil temperature
      real(wp),             intent(in)  :: theta(:)          ! [m3/m3] per-layer soil moisture
      real(wp),             intent(out) :: soil_temp_root    ! [K]     root-fraction-weighted mean
      real(wp),             intent(out) :: theta_mean        ! [m3/m3] depth-weighted column mean
      integer(ik) :: k, nsl
      nsl = soil%n_active
      soil_temp_root = weighted_mean(soil_temp(1:nsl), root_frac, nsl)
      theta_mean = 0.0_wp
      do k = 1_ik, nsl
         theta_mean     = theta_mean     + theta(k) * soil%dz(k)
      end do
      theta_mean = theta_mean / max(-soil%soil_layer_z(nsl+1_ik), tiny_num)
   end subroutine root_zone_environment

   !---~---
   !   canopy_leaf_gas_exchange -- per-cohort leaf gas exchange (GPP, gs, Rd) through the batch kernel,
   ! the patch GPP and leaf-respiration totals, and the frozen leaf-energy coefficients h_coeff_leaf
   ! and g_transp_leaf (docs/science/leaf_gas_exchange.md). psi_leaf is read once per dt_fast from the
   ! state^n leaf water. The stomatal limb (eq. (10), §4.2) takes the cohort's predawn leaf potential,
   ! dmax_psi_leaf, through the kernel's psi keyword; a cohort with no history yet (a recruit, a run's
   ! first step) starts from its patch's top-layer soil potential rather than at fully turgid.
   ! The molar conductances the leaf solve works in and the velocities the energy balance takes
   ! convert at one molar density, the canopy air's, P/(R T_cas): the air whose density
   ! multiplies them in the latent fluxes, so a leaf's transpiration is M_w g_mol de/P whatever
   ! its temperature (docs/science/vegetation_energy_dynamics.md eq. (2)).
   !---~---
   subroutine canopy_leaf_gas_exchange(leaf_photo, hyd_table, soil, static, forcing, aero, dynamic,      &
                                       tcas, qcas, press, rho, cp_cas, gpp, ra_leaf, h_coeff_leaf,       &
                                       g_transp_leaf, gpp_coh, leaf_resp_coh, capture, cdiag)
      type(leaf_photo_table_t),   intent(in)    :: leaf_photo   ! the site's per-PFT leaf parameters
      type(hydraulics_table_t),   intent(in)    :: hyd_table    ! per-PFT leaf PV curve
      type(soil_profile_t),       intent(in)    :: soil         ! surface-layer retention (dmax_psi seed)
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing      ! absorbed PAR per cohort
      type(aero_out_t),           intent(in)    :: aero         ! leaf boundary-layer conductances, wind
      type(fast_dynamic_state_t), intent(in)    :: dynamic      ! leaf/wood temperature + water, CAS CO2, theta(1)
      real(wp),                   intent(in)    :: tcas, qcas, press   ! [K], [kg/kg], [Pa] the canopy air's
      real(wp),                   intent(in)    :: rho, cp_cas         ! [kg/m3], [J/kg/K] its density and specific heat
      real(wp),                   intent(out)   :: gpp, ra_leaf          ! [umol/m2 ground/s] patch totals
      real(wp),                   intent(out)   :: h_coeff_leaf(:), g_transp_leaf(:)
      real(wp),                   intent(inout) :: gpp_coh(:), leaf_resp_coh(:)   ! [umol/plant/s]
      !---~---
      !   Per-cohort diagnostic capture; off, the extra fields are never requested from the batch
      ! kernel.
      !---~---
      logical,                    intent(in)    :: capture
      real(wp),                   intent(inout) :: cdiag(:,:)   ! (N_CO_FAST_ROWS, ncoh) instantaneous values

      !--- Bare-array batch I/O for the per-cohort physiology kernel.
      real(wp) :: par_arr(static%n), vpd_arr(static%n), gb_arr(static%n), psi_leaf_arr(static%n)
      real(wp) :: dmax_psi_arr(static%n), dmax_psi_seed
      real(wp) :: a_gross_arr(static%n), gs_arr(static%n), rd_arr(static%n)
      real(wp) :: a_net_arr(static%n), ci_arr(static%n), cs_arr(static%n)
      real(wp) :: bstom_arr(static%n), bnstom_arr(static%n)
      real(wp) :: e_air, rho_mol, gsw_ms
      integer(ik) :: i, n

      n = static%n
      gpp = 0.0_wp
      ra_leaf = 0.0_wp
      gpp_coh(1:n)       = 0.0_wp
      leaf_resp_coh(1:n) = 0.0_wp
      e_air   = specific_humidity_to_vapor_pressure(qcas, press)
      rho_mol = press / (r_gas * tcas)                             ! [mol/m3] the canopy air's molar density
      do i = 1_ik, n
         par_arr(i)      = cohort_leaf_par(forcing%abs_par(i), static%lai(i), forcing%par_per_w)
         vpd_arr(i)      = air_vpd(dynamic%leaf_temp(i), e_air)
         gb_arr(i)       = aero%leaf_gbw(i) * rho_mol
         psi_leaf_arr(i) = leaf_node_psi(dynamic%leaf_water_mass(i), static%bleaf(i),                 &
               dynamic%wood_water_mass(i), static%bwood_water(i), static%height(i),              &
               hyd_table%pft(static%pft(i)))
      end do
      dmax_psi_seed = top_soil_psi(soil, dynamic%soil_water%theta(1))
      do i = 1_ik, n
         if (is_unset_negative(static%dmax_psi_leaf(i))) then
            dmax_psi_arr(i) = dmax_psi_seed                    ! unset: seed from the soil
         else
            dmax_psi_arr(i) = static%dmax_psi_leaf(i)
         end if
      end do
      if (capture) then
         call leaf_gas_exchange_batch(n, par_arr, dynamic%leaf_temp(1:n), vpd_arr, dynamic%cas%can_co2, press, e_air, &
                                      psi_leaf_arr, gb_arr, leaf_photo, static%pft(1:n),                  &
                                      static%vcmax25(1:n), static%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n),                                            &
                                      a_net=a_net_arr, ci=ci_arr, cs=cs_arr,                            &
                                      beta_stom=bstom_arr, beta_nonstom=bnstom_arr)
         do i = 1_ik, n
            cdiag(CO_ANET,         i) = a_net_arr(i)
            cdiag(CO_AGROSS,       i) = a_gross_arr(i)
            cdiag(CO_GSW,          i) = gs_arr(i)
            cdiag(CO_GBW,          i) = aero%leaf_gbw(i)
            cdiag(CO_CI,           i) = ci_arr(i)
            cdiag(CO_CS,           i) = cs_arr(i)
            cdiag(CO_RD,           i) = rd_arr(i)
            cdiag(CO_BETA_STOM,    i) = bstom_arr(i)
            cdiag(CO_BETA_NONSTOM, i) = bnstom_arr(i)
            cdiag(CO_LEAF_TEMP,    i) = dynamic%leaf_temp(i)
            cdiag(CO_WOOD_TEMP,    i) = dynamic%wood_temp(i)
            cdiag(CO_LEAF_VPD,     i) = vpd_arr(i)
            cdiag(CO_PSI_LEAF,     i) = psi_leaf_arr(i)
            cdiag(CO_ABS_PAR,      i) = forcing%abs_par(i)
            cdiag(CO_ABS_SW,       i) = forcing%abs_sw(i)
            cdiag(CO_ABS_LW,       i) = forcing%abs_lw(i)
            cdiag(CO_WIND,         i) = aero%wind(i)
            cdiag(CO_LEAF_WATER,   i) = dynamic%leaf_water_mass(i)
            cdiag(CO_WOOD_WATER,   i) = dynamic%wood_water_mass(i)
            cdiag(CO_LEAF_PAR,     i) = par_arr(i)
            cdiag(CO_GB_MOL,       i) = gb_arr(i)
            cdiag(CO_PSI_PREDAWN,  i) = dmax_psi_arr(i)
         end do
      else
         call leaf_gas_exchange_batch(n, par_arr, dynamic%leaf_temp(1:n), vpd_arr, dynamic%cas%can_co2, press, e_air, &
                                      psi_leaf_arr, gb_arr, leaf_photo, static%pft(1:n),                  &
                                      static%vcmax25(1:n), static%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n))
      end if
      do i = 1_ik, n
         gsw_ms  = gs_arr(i) / rho_mol
         !---~---
         !   GPP is gross A, and Rd is charged below either way, so the canopy air receives the
         ! leaf's net exchange A_g - Rd whatever the water-stress controls did to A_g.
         !---~---
         gpp     = gpp     + a_gross_arr(i) * static%leaf_area(i) * static%exposed(i) * static%nplant(i)
         gpp_coh(i) = a_gross_arr(i) * static%leaf_area(i) * static%exposed(i)   ! the leaves above the snow
         if (capture) cdiag(CO_GPP_RATE, i) = gpp_coh(i)
         ra_leaf = ra_leaf + rd_arr(i)      * static%leaf_area(i) * static%nplant(i)
         leaf_resp_coh(i) = rd_arr(i) * static%leaf_area(i)
         h_coeff_leaf(i) = sensible_heat_coeff(EFFAREA_HEAT * static%lai(i), aero%leaf_gbh(i), rho, cp_cas)
         g_transp_leaf(i)    = leaf_transp_coeff(EFFAREA_TRANSP, static%lai(i), aero%leaf_gbw(i), gsw_ms)
      end do
   end subroutine canopy_leaf_gas_exchange

   !---~---
   !   canopy_maintenance_respiration -- stem and fine-root maintenance respiration, per cohort
   ! [umol/plant/s] and as patch totals [umol/m2 ground/s] (docs/science/plant_respiration.md
   ! eqs. (3)-(4)), both on the peaked response with leaf Rd's parameters. Elemental kernels: the
   ! array actuals drive the broadcast; the patch-uniform soil temperature scale broadcasts.
   !---~---
   subroutine canopy_maintenance_respiration(ea_rd, hd_rd, ds_rd, soil, root_frac, static, wood_temp, soil_temp, &
                                             ra_stem, ra_root, stem_resp_coh, root_resp_coh)
      real(wp),                  intent(in)    :: ea_rd, hd_rd, ds_rd ! leaf Rd's parameters, in the peaked form
      real(wp),                  intent(in)    :: root_frac(:)        ! [-] the patch's fine roots by layer
      type(soil_profile_t),      intent(in)    :: soil                ! active layer count
      type(fast_static_state_t), intent(in)    :: static
      real(wp),                  intent(in)    :: wood_temp(:)        ! [K] per-cohort wood temperature
      real(wp),                  intent(in)    :: soil_temp(:)        ! [K] PER-LAYER soil temperature
      real(wp),                  intent(out)   :: ra_stem, ra_root    ! [umol/m2 ground/s]
      real(wp),                  intent(inout) :: stem_resp_coh(:), root_resp_coh(:)   ! [umol/plant/s]
      real(wp) :: stem_resp_arr(static%n), root_resp_arr(static%n), tscale_root
      integer(ik) :: i, n
      n = static%n
      !---~---
      !   The temperature response, summed over layers (docs/science/plant_respiration.md eq. (4)),
      ! is patch-uniform: evaluated once and broadcast over the cohorts.
      !---~---
      tscale_root = root_zone_temp_scale(soil_temp, root_frac, soil%n_active, ea_rd, hd_rd, ds_rd)
      ra_stem = 0.0_wp
      ra_root = 0.0_wp
      call stem_maintenance_respiration(wood_temp(1:n), static%diameter(1:n), static%height(1:n),   &
                                   static%wai(1:n), static%nplant(1:n),                          &
                                   static%aboveground_frac(1:n), static%is_woody(1:n),          &
                                   static%stem_resp_factor25(1:n), static%sap_volume(1:n),      &
                                   static%stem_resp_sapwood25(1:n), ea_rd, hd_rd, ds_rd,          &
                                   stem_resp_arr(1:n))
      call fine_root_maintenance_respiration(tscale_root, static%broot(1:n),                        &
                                   static%root_resp_factor25(1:n), root_resp_arr(1:n))
      do i = 1_ik, n
         ra_stem = ra_stem + stem_resp_arr(i) * static%nplant(i)
         ra_root = ra_root + root_resp_arr(i) * static%nplant(i)
         stem_resp_coh(i) = stem_resp_arr(i)
         root_resp_coh(i) = root_resp_arr(i)
      end do
   end subroutine canopy_maintenance_respiration

   !---~---
   !   patch_heterotrophic_respiration -- Rh [umol/m2/s] from the CENTURY matrix on the patch's frozen
   ! pools (docs/science/soil_carbon.md §7, eq. (2)), so the day's fast Rh equals the daily step's
   ! pool debit by construction; xi goes to exports%xi_step for the caller to accumulate. Zero with
   ! soil carbon off.
   !---~---
   subroutine patch_heterotrophic_respiration(cfg, soil, soil_carbon, t_ground, soil_temp_root,          &
                                              theta_mean, rh, exports, surf_rel)
      type(ed_config_t),    intent(in)    :: cfg
      type(soil_profile_t), intent(in)    :: soil          ! the surface layer's retention limits
      type(soil_carbon_state_t),  intent(in)    :: soil_carbon
      real(wp),             intent(in)    :: t_ground, soil_temp_root, theta_mean
      real(wp),             intent(out)   :: rh              ! [umol/m2/s]
      type(step_exports_t), intent(inout) :: exports
      real(wp), optional,   intent(in)    :: surf_rel        ! [-] the porous layers' wetness (pools 1 and 3)
      real(wp) :: xi(N_SOIL_POOL), a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: rh_kgc                                     ! [kgC/m2/day] the matrix's own currency
      if (cfg%soil_carbon%soil_carbon_on) then
         call decomposition_env_scalar(t_ground, soil_temp_root, theta_mean, soil%theta_res(1),                    &
                                  soil%theta_sat(1), soil_carbon, cfg%soil_carbon, xi, surf_rel)
         call decomposition_transfer_matrix(soil_carbon, cfg%soil_carbon, a_mat, k_diag, er)
         !---~---
         !   The matrix returns kgC/m2/day, the pools' currency, which exports%rh_matrix_step keeps for
         ! the seam check; the canopy air gets umol/m2/s (docs/science/soil_carbon.md §7).
         !---~---
         rh_kgc = heterotrophic_respiration(a_mat, k_diag, xi, soil_carbon)
         exports%xi_step        = xi
         exports%rh_matrix_step = rh_kgc                           ! [kgC/m2/day] -- the seam's currency
         rh = rh_kgc * kgcday_2_umols                              ! [umol/m2/s]  -- what the CAS gets
      else
         !---~---
         !   No soil carbon, no soil respiration: the pools are empty and take no litter, so any Rh
         ! here would be carbon from nothing (docs/science/soil_carbon.md §Parameters).
         !---~---
         rh = 0.0_wp
      end if
   end subroutine patch_heterotrophic_respiration

   !---~---
   !   intercept_rain -- the step's canopy interception (canopy_interception) with the canopy
   ! water store on, of the rain only, whatever lies on the ground: the canopy's surface water
   ! is liquid, so the snowfall passes through the canopy to the ground (snow held on the
   ! branches is not represented). A crown the snow buries has its exposed share's leaf and
   ! wood area (crown_exposed_fraction).
   !---~---
   subroutine intercept_rain(dt_fast, cfg, static, forcing, dynamic, veg_surface, rain_through)
      real(wp),                     intent(in)    :: dt_fast
      type(ed_config_t),            intent(in)    :: cfg
      type(fast_static_state_t),    intent(in)    :: static
      type(fast_forcing_t),         intent(in)    :: forcing
      type(fast_dynamic_state_t),   intent(in)    :: dynamic
      type(canopy_water_out_t),     intent(inout) :: veg_surface
      real(wp),                     intent(out)   :: rain_through   ! [kg/m2/s] the rain that reaches the ground
      integer(ik) :: n
      n = static%n
      veg_surface%intercept_leaf = 0.0_wp
      veg_surface%intercept_wood = 0.0_wp
      veg_surface%transfer_lw = 0.0_wp
      veg_surface%f_wet_c = 0.0_wp
      veg_surface%g_surf_leaf = 0.0_wp
      veg_surface%g_surf_wood = 0.0_wp
      rain_through = forcing%rainfall
      if (cfg%canopy_water_on)                                                                                 &
         call canopy_interception(n, forcing%rainfall, dynamic%leaf_surf_water,                                 &
                                  dynamic%wood_surf_water, static%lai, static%wai,                              &
                                  cfg%pft%leaf_surf_water_max(static%pft(1:n)),                                 &
                                  cfg%pft%wood_surf_water_max(static%pft(1:n)), cfg%intercept_k,                &
                                  cfg%intercept_alpha, dt_fast, veg_surface%intercept_leaf, veg_surface%intercept_wood,       &
                                  veg_surface%transfer_lw, veg_surface%f_wet_c, rain_through)
   end subroutine intercept_rain

   !---~---
   !   tissue_coefficients -- each cohort's leaf and wood energy-balance inputs: geometry,
   ! radiation, heat capacities and the start-of-step temperatures
   ! (docs/science/vegetation_energy_dynamics.md), and its surface-water evaporation conductances.
   !---~---
   subroutine tissue_coefficients(dt_fast, cfg, static, forcing, dynamic, aero, rho, cp_cas, tissue, veg_surface, plant)
      real(wp),                     intent(in)    :: dt_fast
      type(ed_config_t),            intent(in)    :: cfg
      type(fast_static_state_t),    intent(in)    :: static
      type(fast_forcing_t),         intent(in)    :: forcing
      type(fast_dynamic_state_t),   intent(in)    :: dynamic
      type(aero_out_t),             intent(in)    :: aero
      real(wp),                     intent(in)    :: rho, cp_cas   ! [kg/m3], [J/kg/K] the canopy air's
      type(plant_energy_env_t),     intent(inout) :: tissue
      type(canopy_water_out_t),     intent(inout) :: veg_surface
      type(plant_water_t),          intent(inout) :: plant
      integer(ik) :: i
      do i = 1_ik, static%n
         tissue%lai(i)    = static%lai(i)
         tissue%abs_sw(i) = forcing%abs_sw(i)
         tissue%abs_lw(i) = forcing%abs_lw(i)
         !---~---
         !   Wood heat capacity over the step: the dry heat of all the wood plus the water network's
         ! whole wood store (docs/science/vegetation_energy_dynamics.md eq. (8)). The surface water
         ! the wood holds joins at the commit, which reads the temperature back with the end-of-step
         ! surface water, so surface water far heavier than its tissue that evaporates within the
         ! step cannot leave its heat to the tissue alone ("The canopy's surface water").
         !---~---
         tissue%wood_dry_hcap(i) = static%wood_dry_hcap(i)
         tissue%wood_wmass(i)    = max(dynamic%wood_water_mass(i), 0.0_wp) * static%nplant(i)
         !--- Leaf heat capacity, the same parts: its water read from dynamic.
         tissue%leaf_dry_hcap(i) = static%leaf_dry_hcap(i)
         tissue%leaf_wmass(i)    = max(dynamic%leaf_water_mass(i), 0.0_wp) * static%nplant(i)
         !---~---
         !   The wood's radiation and sensible coefficient: its whole energy budget in surface_derivs
         ! (docs/science/vegetation_energy_dynamics.md §Leaf and wood energy balance).
         !---~---
         tissue%wai(i)         = static%wai(i)
         tissue%h_coeff_w(i)   = sensible_heat_coeff(pi * static%wai(i), aero%wood_gbh(i), rho, cp_cas)
         tissue%abs_sw_wood(i) = forcing%abs_sw_wood(i)
         tissue%abs_lw_wood(i) = forcing%abs_lw_wood(i)
         !---~---
         !   Tissue heat store a = cap/dt_fast, relaxing from the start-of-step temperature
         ! (docs/science/vegetation_energy_dynamics.md eqs. (3)-(5)): an algebraic closure in every
         ! stage evaluation, not a tableau unknown. The share a snow pack buries is the water
         ! layer's node's (buried_tissue_store), so the balance runs on the exposed share's.
         !---~---
         tissue%buried(i)           = 1.0_wp - static%exposed(i)
         tissue%leaf_hcap_per_dt(i) = static%exposed(i)                                                        &
                                      * tissue_hcap_per_dt(tissue%leaf_dry_hcap(i), tissue%leaf_wmass(i), dt_fast)
         tissue%wood_hcap_per_dt(i) = static%exposed(i)                                                        &
                                      * tissue_hcap_per_dt(tissue%wood_dry_hcap(i), tissue%wood_wmass(i), dt_fast)
         tissue%t_leaf0(i) = dynamic%leaf_temp(i)
         tissue%t_wood0(i) = dynamic%wood_temp(i)
         tissue%leaf_emiss(i) = static%leaf_emissivity(i)
         tissue%wood_emiss(i) = static%wood_emissivity(i)
         plant%nplant(i)   = static%nplant(i)
         !---~---
         !   Surface-water evaporation conductances, from the boundary layers; none with the canopy
         ! water store off.
         !---~---
         if (cfg%canopy_water_on) then
            veg_surface%g_surf_leaf(i) = leaf_surface_evap_coeff(EFFAREA_EVAP, static%lai(i), aero%leaf_gbw(i))
            veg_surface%g_surf_wood(i) = leaf_surface_evap_coeff(EFFAREA_EVAP, static%wai(i), aero%wood_gbw(i))
         end if
      end do
   end subroutine tissue_coefficients

   !---~---
   !   canopy_air_boundary -- the rest of the canopy-air boundary: what a per-stage surface-
   ! layer re-solve needs besides the live canopy air and its environment (set by
   ! refresh_canopy_aerodynamics), the air above, and the ground's radiation and
   ! conductance. The re-solve takes its roughness from the same geometry, so at the state^n
   ! canopy air it reproduces the conductances exactly.
   !---~---
   subroutine canopy_air_boundary(aero_cfg, geom, aero, forcing, rho, cp_cas, press, cas, ground)
      type(aero_params_t),     intent(in)    :: aero_cfg
      type(aero_args_t),       intent(in)    :: geom
      type(aero_out_t),        intent(in)    :: aero
      type(fast_forcing_t),    intent(in)    :: forcing
      real(wp),                intent(in)    :: rho, cp_cas, press   ! [kg/m3], [J/kg/K], [Pa] the canopy air's
      type(cas_boundary_t),    intent(inout) :: cas
      type(surface_env_t),     intent(inout) :: ground
      cas%aero_cfg     = aero_cfg
      cas%geom         = geom
      !---~---
      !   Live inputs license the per-stage re-solve; a record not built here (a test fixture)
      ! leaves mo_live false, and its box conductances are used as given.
      !---~---
      cas%mo_live      = .true.
      cas%box%atm_enthalpy = forcing%enthalpy_atm
      cas%box%atm_specific_humidity = forcing%shv_atm
      cas%box%atm_co2 = forcing%co2_atm
      cas%rho   = rho
      cas%cp    = cp_cas
      cas%press = press
      ground%rnet_soil = forcing%rnet_soil
      ground%sw_snow   = forcing%sw_snow
      ground%sw_ground = forcing%sw_ground
      ground%lw_down   = forcing%lw_down_ground
      ground%ggnet     = aero%ggnet
   end subroutine canopy_air_boundary

   !---~---
   !   transpiration_demand -- the start-of-step surface solve, at the state^n canopy air: each
   ! cohort's transpiration demand on the water network (docs/science/soil_plant_water.md §4).
   ! It takes the canopy air's exchange the pre-pass solved at that same state, so the
   ! surface layer is not solved again, and the ground is not solved at all: the tissues'
   ! exchange does not depend on it (canopy_demand).
   ! The evaporation conductance of each tissue's surface water is first scaled so its
   ! potential evaporation over dt_fast cannot exceed what that water holds plus this step's
   ! interception; then the leaf is credited the liquid enthalpy of the sapflow that replaces
   ! its transpiration, as the stages will, estimated as if the sapflow matched it, and the
   ! surface is solved again.
   !---~---
   subroutine transpiration_demand(dt_fast, canopy_water_on, n, dynamic, frozen)
      real(wp),                   intent(in)    :: dt_fast
      logical,                    intent(in)    :: canopy_water_on
      integer(ik),                intent(in)    :: n
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(frozen_t),      intent(inout) :: frozen
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sf0
      type(cas_boundary_t)  :: cas_start   ! the boundary with the pre-pass's conductances
      real(wp)    :: t_up_wl
      integer(ik) :: i
      y_stage%cas_enthalpy = dynamic%cas%can_enthalpy
      y_stage%cas_shv = dynamic%cas%can_shv
      y_stage%cas_co2 = dynamic%cas%can_co2
      cas_start = frozen%cas
      cas_start%mo_live = .false.
      call canopy_demand(y_stage, cas_start, frozen%tissue, frozen%veg_surface, n, sf0)
      if (canopy_water_on)                                                                                    &
         call limit_surface_evaporation(n, dt_fast, dynamic%leaf_surf_water, dynamic%wood_surf_water,              &
                                     sf0%leaf_surface_evap, sf0%wood_surface_evap, frozen%veg_surface)
      do i = 1_ik, n
         t_up_wl = merge(dynamic%wood_temp(i), dynamic%leaf_temp(i), sf0%transp_c(i) >= 0.0_wp)
         frozen%tissue%qwflux_wl(i) = sf0%transp_c(i) * internal_energy_liquid(t_up_wl)
      end do
      call canopy_demand(y_stage, cas_start, frozen%tissue, frozen%veg_surface, n, sf0)
      frozen%water_drive%e_demand(1:n) = sf0%transp_c(1:n)
   end subroutine transpiration_demand

   !---~---
   !   The tissue a snow pack buries, as the water layer's node holds it over the step
   ! (docs/science/vegetation_energy_dynamics.md, "Buried tissue"): each cohort's buried share of
   ! its leaves' and wood's heat capacity, its tissue and its water (tissue_hcap_per_dt), and their
   ! heat at their start temperatures, from zero kelvin, as the node's phase law counts a dry
   ! heat capacity's.
   !---~---
   pure subroutine buried_tissue_store(tissue, n, dt_fast, cap, energy)
      type(plant_energy_env_t), intent(in)  :: tissue
      integer(ik),              intent(in)  :: n
      real(wp),                 intent(in)  :: dt_fast
      real(wp),                 intent(out) :: cap, energy   ! [J/m2/K], [J/m2]
      real(wp)    :: c_leaf, c_wood
      integer(ik) :: i
      cap    = 0.0_wp
      energy = 0.0_wp
      do i = 1_ik, n
         if (.not. (tissue%buried(i) > 0.0_wp)) cycle
         c_leaf = tissue%buried(i) * tissue_hcap_per_dt(tissue%leaf_dry_hcap(i), tissue%leaf_wmass(i), dt_fast) * dt_fast
         c_wood = tissue%buried(i) * tissue_hcap_per_dt(tissue%wood_dry_hcap(i), tissue%wood_wmass(i), dt_fast) * dt_fast
         cap    = cap + c_leaf + c_wood
         energy = energy + c_leaf * tissue%t_leaf0(i) + c_wood * tissue%t_wood0(i)
      end do
   end subroutine buried_tissue_store

   !---~---
   !   The porous layers as nodes of the water network (docs/science/surface_organic_layer.md,
   ! Water), top down, every layer that holds matter: its retention curve at its depth and its
   ! water; its ice what its energy holds on its freezing curve, the soil's law, held over the step
   ! as a frozen soil layer's, so the network's liquid is the energy's. The water it lacks of its
   ! residual, if it holds less (litter just arrived, dry), is carried as an offset, so its
   ! potential is its curve's dry end and its mass exact. Each layer evaporates by the soil's law
   ! to the canopy air at its own temperature (slot_evaporation_setup), through its own dried top,
   ! which its water sets over the step, and the open pores of the layers above it, held over the
   ! step (porous_cover_path), as the top soil layer's vapour crosses them all.
   !---~---
   subroutine porous_network_drive(params, frozen)
      type(ground_surface_layer_params_t), intent(in)    :: params
      type(frozen_t),                      intent(inout) :: frozen
      type(porous_material_t) :: m
      real(wp)    :: depth, theta, ice, path
      integer(ik) :: k
      associate (o => frozen%water_drive%porous, layer => frozen%porous)
         o%n = 0_ik
         path = 0.0_wp
         do k = 1_ik, N_POROUS
            if (layer(k)%matter <= 0.0_wp) cycle
            m     = params%porous(k)
            depth = porous_layer_depth(layer(k)%matter, m%bulk_density)
            o%n = o%n + 1_ik
            o%slot(o%n)    = k
            o%depth(o%n)   = depth
            o%water(o%n)   = layer(k)%water
            o%curve(o%n)   = porous_layer_curve(m, layer(k))
            theta = layer(k)%water / (rho_h2o * depth)
            o%deficit(o%n) = max(0.0_wp, m%theta_res - theta) * rho_h2o * depth
            o%air(o%n)     = slot_evaporation_setup(frozen%water_drive%ground, layer(k)%temp, path, m%theta_sat)
            path  = path + porous_cover_path(depth, theta, m%theta_sat)
            ice   = (1.0_wp - layer(k)%fliq) * theta
            theta = max(theta, m%theta_res)
            o%ice(o%n) = soil_ice_state(o%curve(o%n), ice, theta - ice)
         end do
      end associate
   end subroutine porous_network_drive

   !---~---
   !   ground_inflow -- what drives the water network besides the transpiration demand and the
   ! pond, which the water layer's free water already is (docs/science/soil_plant_water.md §5):
   ! each layer's temperature and ice, its roots' contact, and the ground evaporation's boundary.
   ! Also the liquid enthalpy of the rain the canopy's surface water exchanges with, which the
   ! start-of-step surface solves read, so it runs before them.
   !---~---
   subroutine ground_inflow(profile, forcing, dynamic, ggnet, tcas, qcas, rho, t_ground, nsl, params, frozen)
      type(soil_profile_t),       intent(in)    :: profile
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      real(wp),                   intent(in)    :: ggnet         ! [m/s] ground <-> canopy-air conductance
      real(wp),                   intent(in)    :: tcas, qcas, rho, t_ground
      integer(ik),                intent(in)    :: nsl
      type(ground_surface_layer_params_t), intent(in)    :: params    ! the porous layers' materials
      type(frozen_t),      intent(inout) :: frozen
      real(wp) :: t_rain         ! [K] the temperature of the rain the canopy's surface water exchanges with
      real(wp) :: depth          ! [m] a porous layer's depth
      type(porous_material_t) :: m
      integer(ik) :: k
      associate (drive => frozen%water_drive)
         !--- The soil's vapour diffuses up through the porous layers, through their open pores.
         drive%ground%cover_path = 0.0_wp
         do k = 1_ik, N_POROUS
            if (frozen%porous(k)%matter <= 0.0_wp) cycle
            m     = params%porous(k)
            depth = porous_layer_depth(frozen%porous(k)%matter, m%bulk_density)
            drive%ground%cover_path = drive%ground%cover_path                                                    &
                                      + porous_cover_path(depth, frozen%porous(k)%water / (rho_h2o * depth), m%theta_sat)
         end do
         drive%ground%snow_free_frac  = 1.0_wp - frozen%ground%snowfac
         drive%ground%t_ground        = t_ground
         drive%ground%q_air           = qcas
         drive%ground%rho_air         = rho
         drive%ground%press           = forcing%press
         !--- Bare-soil aerodynamic resistance; snow_free_frac above weights the flux, not this.
         drive%ground%r_aero = 1.0_wp / max(ggnet, tiny_num)
         drive%t_soil(1:nsl) = dynamic%soil_energy%soil_temp(1:nsl)
         !---~---
         !   Each layer's ice beside the liquid it starts with: the ice is held over the step, and
         ! melts for the water drawn below that liquid (docs/science/soil_biophysics.md, "Frozen soil").
         !---~---
         drive%ice(1:nsl) = soil_ice_state(profile%layer_curve(1:nsl),                                            &
                                           (1.0_wp - dynamic%soil_energy%soil_fliq(1:nsl))                        &
                                           * dynamic%soil_water%theta(1:nsl),                                     &
                                           dynamic%soil_energy%soil_fliq(1:nsl) * dynamic%soil_water%theta(1:nsl))
         !--- Each layer's roots lose their contact with it as it freezes (docs/science/plant_hydraulics.md eq. (14)).
         drive%root_contact(1:nsl) = root_soil_contact(dynamic%soil_energy%soil_temp(1:nsl))
         !---~---
         !   The rain the surface water captures and drips back is liquid at the air's temperature
         ! (precip_enthalpy's rain); the canopy air's when no rain falls.
         !---~---
         t_rain = tcas
         if (forcing%rainfall > tiny_num) t_rain = forcing%air_temp
         frozen%veg_surface%rain_liquid_enthalpy = internal_energy_liquid(t_rain)
      end associate
   end subroutine ground_inflow

   !---~---
   !   The march's start state, the state^n canopy air, soil energy and water, plant water and
   ! surface water; the top soil layer less the heat a film took from it at the stage.
   !---~---
   subroutine pack_march_start(dynamic, dz1, film_heat, n, nsl, y)
      type(fast_dynamic_state_t), intent(in)  :: dynamic
      real(wp),                   intent(in)  :: dz1         ! [m] the top soil layer's thickness
      real(wp),                   intent(in)  :: film_heat   ! [J/m2] what a film took from it
      integer(ik),                intent(in)  :: n, nsl
      type(march_state_t),       intent(out) :: y
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n), y%leaf_surf_water(n), y%wood_surf_water(n))
      y%cas_enthalpy = dynamic%cas%can_enthalpy
      y%cas_shv = dynamic%cas%can_shv
      y%cas_co2 = dynamic%cas%can_co2
      y%soil_energy(1:nsl) = dynamic%soil_energy%soil_energy(1:nsl)
      y%soil_energy(1)     = y%soil_energy(1) - film_heat / dz1
      y%theta(1:nsl)       = dynamic%soil_water%theta(1:nsl)
      y%leaf_water_mass(1:n) = dynamic%leaf_water_mass(1:n)
      y%wood_water_mass(1:n) = dynamic%wood_water_mass(1:n)
      y%leaf_surf_water(1:n) = dynamic%leaf_surf_water(1:n)
      y%wood_surf_water(1:n) = dynamic%wood_surf_water(1:n)
   end subroutine pack_march_start

   !---~---
   !   record_water_flows -- what the water network moved over the step, as the frozen record's
   ! rates (docs/science/soil_plant_water.md §4-§5): ground evaporation, the soil's boundary
   ! flows, its faces and end-of-step water, and the free water's; the per-layer face closure;
   ! root uptake by layer and by cohort with the enthalpy it carries; sapflow and transpiration;
   ! the hydraulic diagnostics; the plant water's advective enthalpy; and the porous layers as the
   ! ground surface's stack reads them (porous_nodes).
   !---~---
   subroutine record_water_flows(dt_fast, cfg, site, static, dynamic, water, out, frozen)
      real(wp),                   intent(in)    :: dt_fast
      type(ed_config_t),          intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(water_work_t),         intent(in)    :: water
      type(fast_output_t),        intent(inout) :: out
      type(frozen_t),      intent(inout) :: frozen
      real(wp)    :: share_tot, face_resid, f_in, f_out, u_liq_soil_k(site%soil_profile%n_active)
      real(wp)    :: t_up_wl, u_liq_up, u_liq_leaf, u_liq_wood, sapflow_gnd, uptake_gnd, transp_gnd
      real(wp)    :: top                 ! [kg/m2] net into the top soil layer across its top face
      integer(ik) :: i, k, n, nsl, soff
      n    = static%n
      nsl  = site%soil_profile%n_active
      soff = 1_ik + water%n_porous
      associate (fl => water%flows, st => water%state)
         !--- Ground evaporation, the soil's boundary flows and the free water, as the network moved them.
         top = fl%face(soff)
         frozen%ground%soil_evap         = fl%evap / dt_fast
         frozen%hydrology%pond_infil     = fl%infil / dt_fast
         frozen%hydrology%q_top          = (top - fl%evap) / (rho_h2o * dt_fast)
         frozen%hydrology%drainage       = fl%drain / dt_fast
         frozen%hydrology%runoff_surf    = fl%spill / dt_fast
         out%report%exports%runoff_surf_export = fl%spill
         out%report%exports%drainage_export    = fl%drain
         frozen%hydrology%runoff_enth    = fl%spill / dt_fast * internal_energy_liquid(frozen%stage%free_temp)
         frozen%hydrology%t_bot          = dynamic%soil_energy%soil_temp(nsl)
         frozen%hydrology%free_water1    = st%water_soil(POND)
         frozen%hydrology%w_flux_frozen(1:nsl)   = 0.0_wp
         frozen%hydrology%w_flux_frozen(1:nsl-1) = fl%face(soff+1:soff+nsl-1) / (rho_h2o * dt_fast)
         do k = 1_ik, nsl
            frozen%hydrology%theta1(k) = st%water_soil(soff + k) / (rho_h2o * site%soil_profile%dz(k))
         end do
         frozen%hydrology%theta_rate(1:nsl) = (frozen%hydrology%theta1(1:nsl) - dynamic%soil_water%theta(1:nsl))  &
                                              / dt_fast

         !---~---
         !   Per-layer face closure: each layer's water change against the flows the soil energy will
         ! advect. It reads round-off, and checks that the faces handed on are the ones that moved the
         ! water; a network sub-step taken unconverged counts as a failed check.
         !---~---
         face_resid = 0.0_wp
         do k = 1_ik, nsl
            f_in = 0.0_wp
            f_out = 0.0_wp
            if (k > 1_ik)   f_in  = fl%face(soff + k - 1_ik)
            if (k < nsl)    f_out = fl%face(soff + k)
            if (k == nsl)   f_out = fl%drain
            if (k == 1_ik)  f_in  = top - fl%evap
            face_resid = face_resid + abs((frozen%hydrology%theta1(k) - dynamic%soil_water%theta(k))           &
                                          * site%soil_profile%dz(k) * rho_h2o - (f_in - f_out - sum(fl%uptake(k, 1:n))))
         end do
         associate (ledger => out%report%ledgers%soil_face_mass)
            ledger%resid   = face_resid
            ledger%worst   = max(ledger%worst, face_resid)
            ledger%abs_sum = ledger%abs_sum + face_resid
            ledger%n_check = ledger%n_check + 1_ik
            if (out%report%work(WORK_WATER_NFAIL) > 0_ik) ledger%n_fail = ledger%n_fail + 1_ik
         end associate

         !---~---
         !   Root uptake in total, by layer and by cohort, and its enthalpy at each layer's start-of-step
         ! temperature (docs/science/soil_plant_water.md §5), the number the wood is credited below. A
         ! patch whose roots took nothing places its zero sink by its root profile.
         !---~---
         frozen%roots%uptake = sum(fl%uptake(1:nsl, 1:n)) / dt_fast
         do k = 1_ik, nsl
            u_liq_soil_k(k) = internal_energy_liquid(dynamic%soil_energy%soil_temp(k))
            frozen%roots%root_share(k) = sum(fl%uptake(k, 1:n))
            frozen%roots%root_heat(k)  = sum(fl%uptake(k, 1:n)) * u_liq_soil_k(k) / dt_fast
         end do
         share_tot = sum(frozen%roots%root_share(1:nsl))
         if (share_tot > tiny_num) then
            frozen%roots%root_share(1:nsl) = frozen%roots%root_share(1:nsl) / share_tot
         else
            frozen%roots%root_share(1:nsl) = static%root_frac_patch(1:nsl)
         end if
         do i = 1_ik, n
            frozen%roots%qloss_frozen(i)   = sum(fl%uptake(1:nsl, i) * u_liq_soil_k(1:nsl)) / dt_fast
            frozen%plant%sapflow_frozen(i) = fl%sapflow(i) / (max(static%nplant(i), tiny_num) * dt_fast)
            frozen%plant%uptake_frozen(i)  = sum(fl%uptake(1:nsl, i)) / (max(static%nplant(i), tiny_num) * dt_fast)
            frozen%tissue%transp_net(i)    = fl%transp(i) / dt_fast
         end do
         frozen%tissue%transp_held = .true.

         !---~---
         !   Per-cohort hydraulic diagnostics, from the same network the physics commits: sapflow and
         ! root uptake per plant, and the wood's end-of-step potential and loss of conductance.
         !---~---
         if (out%capture) then
            do i = 1_ik, n
               out%cdiag(CO_SAPFLOW,     i) = frozen%plant%sapflow_frozen(i)
               out%cdiag(CO_ROOT_UPTAKE, i) = frozen%plant%uptake_frozen(i)
               out%cdiag(CO_WET_FRACTION, i) = frozen%veg_surface%f_wet_c(i)
               out%cdiag(CO_PSI_WOOD,    i) = st%psi_coh(WOOD, i)
               associate (hp => cfg%pft%hydraulics%pft(static%pft(i)))
                  out%cdiag(CO_PLC, i) = 1.0_wp - conductance_retained_fraction(st%psi_coh(WOOD, i), hp%wood_psi50, &
                                                                                hp%wood_kexp)
               end associate
            end do
         end if

         !---~---
         !   Advective enthalpy of the plant water (docs/science/soil_plant_water.md eqs. (4)-(5)): the
         ! network's flows, with upwind temperatures frozen at state^n, since the stage temperatures are
         ! computed from these terms; each tissue's water change carries its start-of-step enthalpy.
         !---~---
         do i = 1_ik, n
            t_up_wl    = merge(dynamic%wood_temp(i), dynamic%leaf_temp(i), fl%sapflow(i) >= 0.0_wp)
            u_liq_up   = internal_energy_liquid(t_up_wl)
            u_liq_leaf = internal_energy_liquid(dynamic%leaf_temp(i))
            u_liq_wood = internal_energy_liquid(dynamic%wood_temp(i))
            sapflow_gnd = fl%sapflow(i) / dt_fast                                  ! [kg/m2 ground/s]
            uptake_gnd  = sum(fl%uptake(1:nsl, i)) / dt_fast
            transp_gnd  = fl%transp(i) / dt_fast
            frozen%tissue%qwflux_wl(i)  = sapflow_gnd * u_liq_up - (sapflow_gnd - transp_gnd) * u_liq_leaf
            frozen%tissue%q_wood_net(i) = frozen%roots%qloss_frozen(i) - uptake_gnd * u_liq_wood                &
                                          - sapflow_gnd * (u_liq_up - u_liq_wood)
            frozen%tissue%wood_water_enth(i) = (uptake_gnd - sapflow_gnd) * u_liq_wood
            frozen%tissue%leaf_water_enth(i) = (sapflow_gnd - transp_gnd) * u_liq_leaf
         end do
      end associate
      call porous_nodes(cfg%ground_surface, dt_fast, water, frozen)
   end subroutine record_water_flows

   !---~---
   !   The porous layers and the free water as the ground surface's stack reads them over the
   ! march (porous_node_env_t; surface_env_t's free water), from what the water network moved:
   ! each layer's record holds its end water with its own evaporation not yet gone, which the
   ! update takes with the march's energy, and its node the end water it keeps, its vapour
   ! taking its whole enthalpy from the layer's energy; its dry heat capacity, its energy and
   ! temperature at the step's start, its thermal resistance at its end water, its evaporation
   ! over the step, and the network's net flow out of its bottom; the free water's net flow into
   ! the ground and its runoff. The heat that water carries is the stack's
   ! (ground_surface_exchange).
   !---~---
   subroutine porous_nodes(params, dt_fast, water, frozen)
      type(ground_surface_layer_params_t), intent(in)    :: params
      real(wp),                            intent(in)    :: dt_fast
      type(water_work_t),                  intent(in)    :: water
      type(frozen_t),                      intent(inout) :: frozen
      type(porous_material_t) :: m
      integer(ik) :: j, k
      associate (fl => water%flows, env => frozen%ground)
         env%free_down = fl%face(POND) / dt_fast
         env%runoff    = fl%spill / dt_fast
         env%free_end  = water%state%water_soil(POND)
         env%porous    = porous_node_env_t()
         do j = 1_ik, water%n_porous
            k = frozen%water_drive%porous%slot(j)
            m = params%porous(k)
            associate (layer => frozen%porous(k), node => env%porous(k))
               layer%water     = water%state%water_soil(1_ik + j) + fl%porous_evap(j)
               node%dry_cap    = cp_organic * layer%matter
               node%depth      = frozen%water_drive%porous%depth(j)
               node%curve      = frozen%water_drive%porous%curve(j)
               node%water      = water%state%water_soil(1_ik + j)
               node%energy     = layer%energy
               node%temp       = layer%temp
               node%fliq       = layer%fliq
               node%resistance = porous_layer_resistance(porous_layer_depth(layer%matter, m%bulk_density),          &
                                                         porous_layer_fill(layer%water,                              &
                                                         porous_layer_capacity(layer%matter, m%bulk_density,         &
                                                                               m%theta_sat)), m%k_dry, m%k_wet)
               node%vapour     = fl%porous_evap(j) / dt_fast
               node%down       = fl%face(1_ik + j) / dt_fast
            end associate
         end do
      end associate
   end subroutine porous_nodes

   !--- The top porous layer, the first top down that holds dry matter; 0 with none.
   pure integer(ik) function top_porous(porous) result(k_top)
      type(porous_layer_state_t), intent(in) :: porous(:)
      integer(ik) :: k
      k_top = 0_ik
      do k = size(porous, kind=ik), 1_ik, -1_ik
         if (porous(k)%matter > 0.0_wp) k_top = k
      end do
   end function top_porous

   !---~---
   !   [-] The porous layers' wetness, all layers together, which the above-ground litter's
   ! decomposition reads: their saturation above residual, the water they hold beyond their
   ! residual over what they hold beyond it saturated, as the soil's wetness is read.
   !---~---
   pure real(wp) function porous_layers_wetness(porous, params) result(wet)
      type(porous_layer_state_t),          intent(in) :: porous(:)
      type(ground_surface_layer_params_t), intent(in) :: params
      real(wp)    :: w, w_sat, w_res
      integer(ik) :: k
      w     = 0.0_wp
      w_sat = 0.0_wp
      w_res = 0.0_wp
      do k = 1_ik, N_POROUS
         if (porous(k)%matter <= 0.0_wp) cycle
         w     = w + porous(k)%water
         w_sat = w_sat + porous_layer_capacity(porous(k)%matter, params%porous(k)%bulk_density,                  &
                                               params%porous(k)%theta_sat)
         w_res = w_res + porous_layer_capacity(porous(k)%matter, params%porous(k)%bulk_density,                  &
                                               params%porous(k)%theta_res)
      end do
      wet = (w - w_res) / max(w_sat - w_res, tiny_num)
   end function porous_layers_wetness

end module ed_fast_integ_frozen
