! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_fast_integ_frozen -- the frozen record one fast step marches against: what is         !
! computed once at the start of the dt_fast, from that moment's state and forcing, and     !
! held fixed across the march (docs/science/numerical_scheme.md §2, Category 0).           !
! build_frozen runs it as named steps:                                                     !
!   * the processes' rates at the start state: canopy aerodynamics, the root zone, leaf    !
!     gas exchange, stem and fine-root maintenance respiration, heterotrophic respiration, !
!     the biotic CO2 source, and the canopy-air capacities and atmosphere conductances;    !
!   * the snow stage, canopy interception, the tissue coefficients, the rest of the        !
!     canopy-air and ground boundary, the transpiration demand and the ground inflow;      !
!   * the soil-plant water network's start state, and the march's.                         !
! The water network is stepped next (water_step), and record_water_flows turns what it     !
! moved into the record's rates (docs/science/soil_plant_water.md §4-§5). Both schemes     !
! march against the same record, so they see the same GPP.                                 !
!==========================================================================================!
module ed_fast_integ_frozen
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : tiny_num, cp_air, rho_h2o, pi, r_gas,                      &
                                     kgcday_2_umols
   use ed_diag_vars,        only :   CO_ANET, CO_AGROSS, CO_GSW, CO_GBW, CO_CI, CO_CS, CO_RD, CO_BETA_STOM,     &
                                     CO_BETA_NONSTOM, CO_LEAF_TEMP, CO_WOOD_TEMP, CO_LEAF_VPD, CO_PSI_LEAF,      &
                                     CO_ABS_PAR, CO_ABS_SW, CO_ABS_LW, CO_WIND, CO_LEAF_WATER, CO_WOOD_WATER,    &
                                     CO_GPP_RATE, CO_LEAF_PAR, CO_GB_MOL, CO_PSI_PREDAWN, CO_CA, CO_PRESSURE,    &
                                     CO_PSI_WOOD, CO_PLC, CO_SAPFLOW, CO_ROOT_UPTAKE
   use ed_plant_hydr_lib,     only : conductance_retained_fraction
   use ed_plant_hydraulics,   only : leaf_node_psi
   use ed_surface, only : snow_out_t
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
   use ed_fast_integ_derivs, only : surface_derivs, set_water_state
   use ed_state_defaults,   only : is_unset_negative
   use ed_biophysics_params, only : aero_params_t
   use ed_patch_state_vars, only : cas_state_t, soil_carbon_state_t, cas_capacities
   use ed_biogeochem_params, only :   N_SOIL_POOL
   use ed_canopy_aerodynamics,   only : canopy_aerodynamics, cas_atm_conductances
   use ed_plant_energy, only :   plant_energy_env_t, sensible_heat_coeff, leaf_transp_coeff, leaf_film_coeff, &
                                 tissue_hcap_per_dt, EFFAREA_HEAT, EFFAREA_TRANSP, EFFAREA_EVAP
   use ed_canopy_water, only   : canopy_water_out_t, canopy_interception, limit_film_evaporation
   use ed_leaf_gas_exchange,   only : leaf_gas_exchange_batch, cohort_leaf_par
   use ed_respiration, only :   stem_maintenance_respiration, fine_root_maintenance_respiration, root_zone_temp_scale, &
                                decomposition_env_scalar, decomposition_transfer_matrix, heterotrophic_respiration
   use ed_surface, only : surface_env_t, prepare_snow_stage
   use ed_soil_water, only : pond_inflow
   use ed_therm_lib,        only :   internal_energy_liquid, liquid_energy_to_temp,                              &
                                     cas_enthalpy_to_temp, air_vpd,                                             &
                                     soil_thermal_conductivity, specific_humidity_to_vapor_pressure
   use ed_numerics,         only :   weighted_mean
   implicit none
   private

   public :: build_frozen, record_water_flows

contains

   !---------------------------------------------------------------------------------------!
   ! build_frozen -- the frozen record both schemes march against, as named steps, and the !
   ! water network's and the march's start states. It only reads the dynamic state: the    !
   ! snow stage advances a copy of the pack, which the update commits. The record is the   !
   ! thread's, reused step to step (prepare_frozen). out receives the aerodynamics, the    !
   ! per-cohort rates, the step's GPP, NEE and decomposition scalars, and, with            !
   ! out%capture, the per-cohort diagnostics (docs/science/diagnostics.md §1, stage [2]).  !
   !---------------------------------------------------------------------------------------!
   subroutine build_frozen(dt_fast, cfg, site, static, forcing, dynamic, out, water, frozen, y)
      real(wp),                   intent(in)    :: dt_fast
      type(ed_config_t),          intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(fast_output_t),        intent(inout) :: out
      type(water_work_t),         intent(inout) :: water    !< the soil-plant water network
      type(frozen_t),             intent(inout) :: frozen
      type(march_state_t),       intent(out)   :: y
      !----- The canopy air's temperature and humidity at the start state, the air's pressure and     !
      !      density, and the ground's temperature, which the later steps read. ---------------------!
      real(wp)    :: tcas, qcas, press, rho, t_ground
      real(wp)    :: can_dmol      !< [mol/m3] the canopy air's dry-air molar density
      !----- The processes' patch totals [umol/m2 ground/s] and the root zone's drivers. ----------!
      real(wp)    :: gpp, ra_leaf, ra_stem, ra_root, rh, soil_temp_root, theta_mean
      real(wp)    :: throughfall   !< [kg/m2/s] the rain and snow through the canopy
      integer(ik) :: n, nsl

      n   = static%n
      nsl = site%soil_profile%n_active
      call prepare_frozen(n, nsl, frozen)

      !----- The processes' rates at the start state. --------------------------------------------!
      call refresh_canopy_aerodynamics(cfg%aero, frozen%cas%env, static, forcing, dynamic%cas,               &
                                       dynamic%soil_energy%soil_temp(1), dynamic%leaf_temp, dynamic%wood_temp, &
                                       out%aero, tcas, qcas, press, rho, t_ground)
      call root_zone_environment(site%soil_profile, static%root_frac_patch(1:nsl), dynamic%soil_energy%soil_temp, &
                                 dynamic%soil_water%theta, soil_temp_root, theta_mean)
      call canopy_leaf_gas_exchange(static%leaf_photo, cfg%pft%hydraulics, site%soil_profile, static, forcing,     &
                                    out%aero, dynamic, qcas, press, rho, gpp, ra_leaf, frozen%tissue%h_coeff_leaf, &
                                    frozen%tissue%g_transp_leaf, out%gpp_coh, out%leaf_resp_coh, out%capture,     &
                                    out%cdiag)
      call canopy_maintenance_respiration(cfg%leaf_physiology%ea_rd, cfg%leaf_physiology%hd_rd,                   &
                                          cfg%leaf_physiology%ds_rd, site%soil_profile,                           &
                                          static%root_frac_patch(1:nsl), static, dynamic%wood_temp,               &
                                          dynamic%soil_energy%soil_temp, ra_stem, ra_root, out%stem_resp_coh,     &
                                          out%root_resp_coh)
      call patch_heterotrophic_respiration(cfg, site%soil_profile, static%soil_carbon, t_ground, soil_temp_root, &
                                           theta_mean, rh, out%report%exports)
      !----- The biotic CO2 source (docs/science/plant_respiration.md eq. (7)), positive to the canopy  !
      !      air. slow_co2_rate is the autotrophic part this loop cannot compute, the allocator's daily !
      !      growth respiration on realized growth, handed down as a frozen rate with the              !
      !      opposite-signed starvation correction (§5-§6 there). ----------------------------------!
      frozen%cas%nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp + static%slow_co2_rate
      out%report%exports%gpp_last = gpp
      out%report%exports%nee_last = frozen%cas%nee_biotic
      !----- The canopy air's capacities, and its exchange with the atmosphere at the start state    !
      !      (docs/science/canopy_aerodynamics.md eq. (6)): heat on heat_transfer_factor, vapour and  !
      !      CO2 on vapour_transfer_factor (equal, z0q = z0h). -----------------------------------!
      call cas_capacities(rho, qcas, dynamic%cas%can_depth, can_dmol, frozen%cas%box%air_mass_capacity,         &
                          frozen%cas%box%air_molar_capacity)
      call cas_atm_conductances(rho, can_dmol, out%aero%ustar, out%aero%heat_transfer_factor,                  &
                                out%aero%vapour_transfer_factor, frozen%cas%box%atm_conductance_enthalpy,     &
                                frozen%cas%box%atm_conductance_vapor, frozen%cas%box%atm_conductance_co2)

      !----- The snow stage (docs/science/snow_biophysics.md §The snow store): the step's snow and rain  !
      !      into the pack and the liquid it holds out to the pond, before the water network, which takes !
      !      that water this step. The pack's exchange with the canopy air and the soil is the march's.   !
      !      A no-op without a pack. -----------------------------------------------------------------!
      frozen%ground%pack   = dynamic%snow
      frozen%ground%params = cfg%snow
      frozen%ground%dt     = dt_fast
      call prepare_snow_stage(forcing%snowfall, forcing%rainfall, forcing%air_temp, cfg%aero%soil_rough,          &
                              frozen%ground, frozen%snow)
      !----- The top soil node the ground surface conducts to, under the snow's base. -----------------!
      associate (prof => site%soil_profile)
         frozen%ground%dz_soil_top = max(-prof%z_node(1), tiny_num)
         frozen%ground%k_soil_top  = soil_thermal_conductivity(dynamic%soil_water%theta(1),                       &
                                                               dynamic%soil_energy%soil_fliq(1), prof%theta_sat(1), &
                                                               prof%solid_conductivity(1), prof%dry_conductivity(1))
      end associate

      call intercept_rain(dt_fast, cfg, static, forcing, dynamic, frozen%snow%exists, frozen%film, throughfall)
      call tissue_coefficients(dt_fast, cfg, static, forcing, dynamic, out%aero, rho, frozen%tissue, frozen%film, &
                               frozen%plant)
      call canopy_air_boundary(cfg%aero, static%geom, out%aero, forcing, rho, press, frozen%cas, frozen%ground)
      call ground_inflow(static, site%soil_profile, forcing, dynamic, out%aero%ggnet, throughfall, tcas, qcas, rho, &
                         t_ground, nsl, frozen)
      call transpiration_demand(dt_fast, cfg%canopy_water_on, n, dynamic, t_ground, frozen)

      !----- The start states: the water network's, from the stores' state^n water, and the march's. --!
      call set_water_state(static, site%soil_profile, cfg%pft%hydraulics, dynamic%soil_water%theta(1:nsl),       &
                           frozen%water_drive%ice(1:nsl),                                                        &
                           dynamic%soil_water%w_surface, dynamic%soil_water%w_surface_energy,                    &
                           dynamic%leaf_water_mass(1:n), dynamic%wood_water_mass(1:n), water)
      call pack_march_start(dynamic, n, nsl, y)
   end subroutine build_frozen

   !----- The frozen record for a step, reused from the last: every scalar back at its default, as  !
   !      a fresh record has it, and the per-cohort and per-layer arrays kept when their sizes still  !
   !      fit. The plant water's heat terms start at zero: the start-of-step surface solves read them !
   !      (as q_extra, whose default is 0) before the water network's flows give them their values. -!
   subroutine prepare_frozen(n, nsl, frozen)
      integer(ik),    intent(in)    :: n, nsl
      type(frozen_t), intent(inout) :: frozen
      type(frozen_t) :: defaults   !< default-initialised and never written: the defaults
      type(frozen_t) :: arrays     !< holds the arrays while the scalars are reset
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
      call fit(frozen%film%f_wet_c, n)
      call fit(frozen%film%g_film_leaf, n)
      call fit(frozen%film%g_film_w, n)
      call fit(frozen%film%intercept_leaf, n)
      call fit(frozen%film%intercept_wood, n)
      call fit(frozen%film%transfer_lw, n)
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

   !----- Move every array of the frozen record from one record to another (the list prepare_frozen !
   !      sizes). An array this list missed would be freed by the reset and allocated afresh. -------!
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
      call move_alloc(from%film%f_wet_c, to%film%f_wet_c)
      call move_alloc(from%film%g_film_leaf, to%film%g_film_leaf)
      call move_alloc(from%film%g_film_w, to%film%g_film_w)
      call move_alloc(from%film%intercept_leaf, to%film%intercept_leaf)
      call move_alloc(from%film%intercept_wood, to%film%intercept_wood)
      call move_alloc(from%film%transfer_lw, to%film%transfer_lw)
      call move_alloc(from%roots%root_share, to%roots%root_share)
      call move_alloc(from%roots%qloss_frozen, to%roots%qloss_frozen)
      call move_alloc(from%roots%root_heat, to%roots%root_heat)
      call move_alloc(from%plant%nplant, to%plant%nplant)
      call move_alloc(from%plant%sapflow_frozen, to%plant%sapflow_frozen)
      call move_alloc(from%plant%uptake_frozen, to%plant%uptake_frozen)
      call move_alloc(from%hydrology%theta1, to%hydrology%theta1)
      call move_alloc(from%water_drive%e_demand, to%water_drive%e_demand)
   end subroutine move_frozen_arrays

   !----- Make a an array of m elements, keeping it when it already is one. ----------------------!
   pure subroutine fit(a, m)
      real(wp), allocatable, intent(inout) :: a(:)
      integer(ik),           intent(in)    :: m
      if (allocated(a)) then
         if (size(a, kind=ik) == m) return
         deallocate(a)
      end if
      allocate(a(m))
   end subroutine fit

   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_aerodynamics -- the aerodynamic environment: its atmosphere from the        !
   ! forcing at the canopy-air top, the reference level the surface layer is solved to, and its  !
   ! canopy air diagnosed from the CAS state; then the whole-canopy aerodynamics solve (friction !
   ! velocity, scalar profile factors, ground conductance, per-cohort boundary layers;           !
   ! docs/science/canopy_aerodynamics.md §2-§4).                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_aerodynamics(aero_cfg, aenv, static, forcing, cas, soil_temp_top,    &
                                          leaf_temp, wood_temp, aero, tcas, qcas, press, rho, t_ground)
      type(aero_params_t),       intent(in)    :: aero_cfg
      type(aero_env_t),          intent(inout) :: aenv
      type(fast_static_state_t), intent(in)    :: static     !< the cohorts, the canopy geometry and its top
      type(fast_forcing_t),      intent(in)    :: forcing
      type(cas_state_t),         intent(in)    :: cas
      real(wp),                  intent(in)    :: soil_temp_top    !< [K] top soil-node temperature (ground skin)
      real(wp),                  intent(in)    :: leaf_temp(:)     !< [K] per-cohort leaf temperature
      real(wp),                  intent(in)    :: wood_temp(:)     !< [K] per-cohort wood temperature
      type(aero_out_t),          intent(inout) :: aero
      real(wp),                  intent(out)   :: tcas, qcas, press, rho, t_ground
      call set_aero_env_atm(aenv, forcing%wind, static%z_top, forcing%air_temp, forcing%shv_atm)
      tcas = cas_enthalpy_to_temp(cas%can_enthalpy, cas%can_shv)
      qcas = cas%can_shv
      press = forcing%press
      rho = forcing%rho_air
      t_ground = soil_temp_top
      call set_aero_env_canopy(aenv, tcas, qcas, t_ground)
      call canopy_aerodynamics(aero_cfg, aenv, static%geom, static%n, static%order, static%height, static%lai, &
                               leaf_temp, wood_temp, static%leaf_width, static%branch_diam, aero)
   end subroutine refresh_canopy_aerodynamics

   !---------------------------------------------------------------------------------------!
   ! root_zone_environment -- root-weighted soil temperature and the column-mean soil moisture     !
   ! that drive fine-root and heterotrophic respiration.                                             !
   !---------------------------------------------------------------------------------------!
   pure subroutine root_zone_environment(soil, root_frac, soil_temp, theta, soil_temp_root, theta_mean)
      type(soil_profile_t), intent(in)  :: soil
      real(wp),             intent(in)  :: root_frac(:)      !< [-]     the patch's fine roots by layer
      real(wp),             intent(in)  :: soil_temp(:)      !< [K]     per-layer soil temperature
      real(wp),             intent(in)  :: theta(:)          !< [m3/m3] per-layer soil moisture
      real(wp),             intent(out) :: soil_temp_root    !< [K]     root-fraction-weighted mean
      real(wp),             intent(out) :: theta_mean        !< [m3/m3] depth-weighted column mean
      integer(ik) :: k, nsl
      nsl = soil%n_active
      soil_temp_root = weighted_mean(soil_temp(1:nsl), root_frac, nsl)
      theta_mean = 0.0_wp
      do k = 1_ik, nsl
         theta_mean     = theta_mean     + theta(k) * soil%dz(k)
      end do
      theta_mean = theta_mean / max(-soil%soil_layer_z(nsl+1_ik), tiny_num)
   end subroutine root_zone_environment

   !---------------------------------------------------------------------------------------!
   ! canopy_leaf_gas_exchange -- per-cohort leaf gas exchange (GPP, gs, Rd) through the batch kernel,  !
   ! the patch GPP and leaf-respiration totals, and the frozen leaf-energy coefficients h_coeff_leaf   !
   ! and g_transp_leaf (docs/science/leaf_gas_exchange.md). psi_leaf is read once per dt_fast from the !
   ! state^n leaf water. The stomatal limb (eq. (10), §4.2) takes the cohort's predawn leaf potential, !
   ! dmax_psi_leaf, through the kernel's psi keyword; a cohort with no history yet (a recruit, a run's !
   ! first step) starts from its patch's top-layer soil potential rather than at fully turgid.        !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_leaf_gas_exchange(leaf_photo, hyd_table, soil, static, forcing, aero, dynamic,      &
                                       qcas, press, rho, gpp, ra_leaf, h_coeff_leaf, g_transp_leaf,        &
                                       gpp_coh, leaf_resp_coh, capture, cdiag)
      type(leaf_photo_table_t),   intent(in)    :: leaf_photo   !< the site's per-PFT leaf parameters
      type(hydraulics_table_t),   intent(in)    :: hyd_table    !< per-PFT leaf PV curve
      type(soil_profile_t),       intent(in)    :: soil         !< surface-layer retention (dmax_psi seed)
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing      !< absorbed PAR per cohort
      type(aero_out_t),           intent(in)    :: aero         !< leaf boundary-layer conductances, wind
      type(fast_dynamic_state_t), intent(in)    :: dynamic      !< leaf/wood temperature + water, CAS CO2, theta(1)
      real(wp),                   intent(in)    :: qcas, press, rho
      real(wp),                   intent(out)   :: gpp, ra_leaf          !< [umol/m2 ground/s] patch totals
      real(wp),                   intent(out)   :: h_coeff_leaf(:), g_transp_leaf(:)
      real(wp),                   intent(inout) :: gpp_coh(:), leaf_resp_coh(:)   !< [umol/plant/s]
      !----- Per-cohort diagnostic capture; off, the extra fields are never requested from the batch  !
      !      kernel. ----------------------------------------------------------------------------!
      logical,                    intent(in)    :: capture
      real(wp),                   intent(inout) :: cdiag(:,:)   !< (N_CO_FAST_ROWS, ncoh) instantaneous values

      !----- Bare-array batch I/O for the per-cohort physiology kernel. ----------------------------!
      real(wp) :: par_arr(static%n), vpd_arr(static%n), gb_arr(static%n), rho_mol_arr(static%n), &
            psi_leaf_arr(static%n)
      real(wp) :: dmax_psi_arr(static%n), dmax_psi_seed
      real(wp) :: a_gross_arr(static%n), gs_arr(static%n), rd_arr(static%n)
      real(wp) :: a_net_arr(static%n), ci_arr(static%n), cs_arr(static%n)
      real(wp) :: bstom_arr(static%n), bnstom_arr(static%n)
      real(wp) :: e_air, gsw_ms
      integer(ik) :: i, n

      n = static%n
      gpp = 0.0_wp
      ra_leaf = 0.0_wp
      gpp_coh(1:n)       = 0.0_wp
      leaf_resp_coh(1:n) = 0.0_wp
      e_air = specific_humidity_to_vapor_pressure(qcas, press)     ! loop-invariant
      do i = 1_ik, n
         rho_mol_arr(i)  = press / (r_gas * dynamic%leaf_temp(i))
         par_arr(i)      = cohort_leaf_par(forcing%abs_par(i), static%lai(i), forcing%par_per_w)
         vpd_arr(i)      = air_vpd(dynamic%leaf_temp(i), e_air)
         gb_arr(i)       = aero%leaf_gbw(i) * rho_mol_arr(i)
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
         call leaf_gas_exchange_batch(n, par_arr, dynamic%leaf_temp(1:n), vpd_arr, dynamic%cas%can_co2, press, &
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
            cdiag(CO_CA,           i) = dynamic%cas%can_co2
            cdiag(CO_PRESSURE,     i) = press
         end do
      else
         call leaf_gas_exchange_batch(n, par_arr, dynamic%leaf_temp(1:n), vpd_arr, dynamic%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, leaf_photo, static%pft(1:n),                  &
                                      static%vcmax25(1:n), static%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n))
      end if
      do i = 1_ik, n
         gsw_ms  = gs_arr(i) / max(rho_mol_arr(i), tiny_num)
         !----- GPP is gross A, and Rd is charged below either way, so the canopy air receives the !
         !      leaf's net exchange A_g - Rd whatever the water-stress controls did to A_g. -------!
         gpp     = gpp     + a_gross_arr(i) * static%leaf_area(i) * static%nplant(i)
         gpp_coh(i) = a_gross_arr(i) * static%leaf_area(i)
         if (capture) cdiag(CO_GPP_RATE, i) = a_gross_arr(i) * static%leaf_area(i)
         ra_leaf = ra_leaf + rd_arr(i)      * static%leaf_area(i) * static%nplant(i)
         leaf_resp_coh(i) = rd_arr(i) * static%leaf_area(i)
         h_coeff_leaf(i) = sensible_heat_coeff(EFFAREA_HEAT * static%lai(i), aero%leaf_gbh(i), rho, cp_air)
         g_transp_leaf(i)    = leaf_transp_coeff(EFFAREA_TRANSP, static%lai(i), aero%leaf_gbw(i), gsw_ms)
      end do
   end subroutine canopy_leaf_gas_exchange

   !---------------------------------------------------------------------------------------!
   ! canopy_maintenance_respiration -- stem and fine-root maintenance respiration, per cohort       !
   ! [umol/plant/s] and as patch totals [umol/m2 ground/s] (docs/science/plant_respiration.md       !
   ! eqs. (3)-(4)), both on the peaked response with leaf Rd's parameters. Elemental kernels: the   !
   ! array actuals drive the broadcast; the patch-uniform soil temperature scale broadcasts.        !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_maintenance_respiration(ea_rd, hd_rd, ds_rd, soil, root_frac, static, wood_temp, soil_temp, &
                                             ra_stem, ra_root, stem_resp_coh, root_resp_coh)
      real(wp),                  intent(in)    :: ea_rd, hd_rd, ds_rd !< leaf Rd's parameters, in the peaked form
      real(wp),                  intent(in)    :: root_frac(:)        !< [-] the patch's fine roots by layer
      type(soil_profile_t),      intent(in)    :: soil                !< active layer count
      type(fast_static_state_t), intent(in)    :: static
      real(wp),                  intent(in)    :: wood_temp(:)        !< [K] per-cohort wood temperature
      real(wp),                  intent(in)    :: soil_temp(:)        !< [K] PER-LAYER soil temperature
      real(wp),                  intent(out)   :: ra_stem, ra_root    !< [umol/m2 ground/s]
      real(wp),                  intent(inout) :: stem_resp_coh(:), root_resp_coh(:)   !< [umol/plant/s]
      real(wp) :: stem_resp_arr(static%n), root_resp_arr(static%n), tscale_root
      integer(ik) :: i, n
      n = static%n
      !----- The temperature response, summed over layers (docs/science/plant_respiration.md eq. (4)), !
      !      is patch-uniform: evaluated once and broadcast over the cohorts. ---------------------!
      tscale_root = root_zone_temp_scale(soil_temp, root_frac, soil%n_active, ea_rd, hd_rd, ds_rd)
      ra_stem = 0.0_wp
      ra_root = 0.0_wp
      call stem_maintenance_respiration(wood_temp(1:n), static%dbh(1:n), static%height(1:n),   &
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

   !---------------------------------------------------------------------------------------!
   ! patch_heterotrophic_respiration -- Rh [umol/m2/s] from the CENTURY matrix on the patch's frozen  !
   ! pools (docs/science/soil_carbon.md §7, eq. (2)), so the day's fast Rh equals the daily step's    !
   ! pool debit by construction; xi goes to exports%xi_step for the caller to accumulate. Zero with   !
   ! soil carbon off.                                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine patch_heterotrophic_respiration(cfg, soil, soil_carbon, t_ground, soil_temp_root,          &
                                              theta_mean, rh, exports)
      type(ed_config_t),    intent(in)    :: cfg
      type(soil_profile_t), intent(in)    :: soil          !< the surface layer's retention limits
      type(soil_carbon_state_t),  intent(in)    :: soil_carbon
      real(wp),             intent(in)    :: t_ground, soil_temp_root, theta_mean
      real(wp),             intent(out)   :: rh              !< [umol/m2/s]
      type(step_exports_t), intent(inout) :: exports
      real(wp) :: xi(N_SOIL_POOL), a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: rh_kgc                                     !< [kgC/m2/day] the matrix's own currency
      if (cfg%soil_carbon%soil_carbon_on) then
         call decomposition_env_scalar(t_ground, soil_temp_root, theta_mean, soil%theta_res(1),                    &
                                  soil%theta_sat(1), soil_carbon, cfg%soil_carbon, xi)
         call decomposition_transfer_matrix(soil_carbon, cfg%soil_carbon, a_mat, k_diag, er)
         !----- The matrix returns kgC/m2/day, the pools' currency, which exports%rh_matrix_step keeps for !
         !      the seam check; the canopy air gets umol/m2/s (docs/science/soil_carbon.md §7). ----------!
         rh_kgc = heterotrophic_respiration(a_mat, k_diag, xi, soil_carbon)
         exports%xi_step        = xi
         exports%rh_matrix_step = rh_kgc                           ! [kgC/m2/day] -- the seam's currency
         rh = rh_kgc * kgcday_2_umols                              ! [umol/m2/s]  -- what the CAS gets
      else
         !----- No soil carbon, no soil respiration: the pools are empty and take no litter, so any Rh   !
         !      here would be carbon from nothing (docs/science/soil_carbon.md §Parameters). ----------!
         rh = 0.0_wp
      end if
   end subroutine patch_heterotrophic_respiration

   !---------------------------------------------------------------------------------------!
   ! intercept_rain -- the step's canopy interception (canopy_interception) with the canopy  !
   ! water store on. None under a pack, which takes the rain and the snow; throughfall        !
   ! carries sub-threshold snowfall, which forms no pack, as liquid.                          !
   !---------------------------------------------------------------------------------------!
   subroutine intercept_rain(dt_fast, cfg, static, forcing, dynamic, pack_exists, film, throughfall)
      real(wp),                     intent(in)    :: dt_fast
      type(ed_config_t),            intent(in)    :: cfg
      type(fast_static_state_t),    intent(in)    :: static
      type(fast_forcing_t),         intent(in)    :: forcing
      type(fast_dynamic_state_t),   intent(in)    :: dynamic
      logical,                      intent(in)    :: pack_exists
      type(canopy_water_out_t),     intent(inout) :: film
      real(wp),                     intent(out)   :: throughfall   !< [kg/m2/s] what reaches the ground
      integer(ik) :: n
      n = static%n
      film%intercept_leaf = 0.0_wp
      film%intercept_wood = 0.0_wp
      film%transfer_lw = 0.0_wp
      film%f_wet_c = 0.0_wp
      film%g_film_leaf = 0.0_wp
      film%g_film_w = 0.0_wp
      throughfall = forcing%rainfall + forcing%snowfall
      if (cfg%canopy_water_on .and. .not. pack_exists)                                                         &
         call canopy_interception(n, forcing%rainfall + forcing%snowfall, dynamic%leaf_surf_water,              &
                                  dynamic%wood_surf_water, static%lai, static%wai,                              &
                                  cfg%pft%leaf_surf_water_max(static%pft(1:n)),                                 &
                                  cfg%pft%wood_surf_water_max(static%pft(1:n)), cfg%intercept_k,                &
                                  cfg%intercept_alpha, dt_fast, film%intercept_leaf, film%intercept_wood,       &
                                  film%transfer_lw, film%f_wet_c, throughfall)
   end subroutine intercept_rain

   !---------------------------------------------------------------------------------------!
   ! tissue_coefficients -- each cohort's leaf and wood energy-balance inputs: geometry,      !
   ! radiation, heat capacities and the start-of-step temperatures                            !
   ! (docs/science/vegetation_energy_dynamics.md), and its film-evaporation conductances.     !
   !---------------------------------------------------------------------------------------!
   subroutine tissue_coefficients(dt_fast, cfg, static, forcing, dynamic, aero, rho, tissue, film, plant)
      real(wp),                     intent(in)    :: dt_fast
      type(ed_config_t),            intent(in)    :: cfg
      type(fast_static_state_t),    intent(in)    :: static
      type(fast_forcing_t),         intent(in)    :: forcing
      type(fast_dynamic_state_t),   intent(in)    :: dynamic
      type(aero_out_t),             intent(in)    :: aero
      real(wp),                     intent(in)    :: rho
      type(plant_energy_env_t),     intent(inout) :: tissue
      type(canopy_water_out_t),     intent(inout) :: film
      type(plant_water_t),          intent(inout) :: plant
      integer(ik) :: i
      do i = 1_ik, static%n
         tissue%lai(i)    = static%lai(i)
         tissue%abs_sw(i) = forcing%abs_sw(i)
         tissue%abs_lw(i) = forcing%abs_lw(i)
         !----- Wood heat capacity over the step: the dry heat of all the wood plus the water network's    !
         !      whole wood store (docs/science/vegetation_energy_dynamics.md eq. (8)). The film the wood   !
         !      holds joins at the commit, which reads the temperature back with the end-of-step film, so  !
         !      a film far heavier than its tissue that evaporates within the step cannot leave its heat   !
         !      to the tissue alone ("The canopy's interception film"). -----------------------------------!
         tissue%wood_dry_hcap(i) = static%wood_dry_hcap(i)
         tissue%wood_wmass(i)    = max(dynamic%wood_water_mass(i), 0.0_wp) * static%nplant(i)
         !----- Leaf heat capacity, the same parts: its water read from dynamic. ----------------------------!
         tissue%leaf_dry_hcap(i) = static%leaf_dry_hcap(i)
         tissue%leaf_wmass(i)    = max(dynamic%leaf_water_mass(i), 0.0_wp) * static%nplant(i)
         !----- The wood's radiation and sensible coefficient: its whole energy budget in surface_derivs   !
         !      (docs/science/vegetation_energy_dynamics.md §Leaf and wood energy balance). -------------!
         tissue%wai(i)         = static%wai(i)
         tissue%h_coeff_w(i)   = sensible_heat_coeff(pi * static%wai(i), aero%wood_gbh(i), rho, cp_air)
         tissue%abs_sw_wood(i) = forcing%abs_sw_wood(i)
         tissue%abs_lw_wood(i) = forcing%abs_lw_wood(i)
         !----- Tissue heat store a = cap/dt_fast, relaxing from the start-of-step temperature             !
         !      (docs/science/vegetation_energy_dynamics.md eqs. (3)-(5)): an algebraic closure in every  !
         !      stage evaluation, not a tableau unknown. -------------------------------------------------!
         tissue%leaf_hcap_per_dt(i) = tissue_hcap_per_dt(tissue%leaf_dry_hcap(i), tissue%leaf_wmass(i), dt_fast)
         tissue%wood_hcap_per_dt(i) = tissue_hcap_per_dt(tissue%wood_dry_hcap(i), tissue%wood_wmass(i), dt_fast)
         tissue%t_leaf0(i) = dynamic%leaf_temp(i)
         tissue%t_wood0(i) = dynamic%wood_temp(i)
         tissue%leaf_emiss(i) = static%leaf_emissivity(i)
         tissue%wood_emiss(i) = static%wood_emissivity(i)
         plant%nplant(i)   = static%nplant(i)
         !----- Film-evaporation conductances, from the boundary layers; none with the canopy water store  !
         !      off. -------------------------------------------------------------------------------------!
         if (cfg%canopy_water_on) then
            film%g_film_leaf(i) = leaf_film_coeff(EFFAREA_EVAP, static%lai(i), aero%leaf_gbw(i))
            film%g_film_w(i) = leaf_film_coeff(EFFAREA_EVAP, static%wai(i), aero%wood_gbw(i))
         end if
      end do
   end subroutine tissue_coefficients

   !---------------------------------------------------------------------------------------!
   ! canopy_air_boundary -- the rest of the canopy-air boundary: what a per-stage surface-     !
   ! layer re-solve needs besides the live canopy air and its environment (set by            !
   ! refresh_canopy_aerodynamics), the air above, and the ground's radiation and              !
   ! conductance. The re-solve takes its roughness from the same geometry, so at the state^n  !
   ! canopy air it reproduces the conductances exactly.                                       !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_air_boundary(aero_cfg, geom, aero, forcing, rho, press, cas, ground)
      type(aero_params_t),     intent(in)    :: aero_cfg
      type(aero_args_t),       intent(in)    :: geom
      type(aero_out_t),        intent(in)    :: aero
      type(fast_forcing_t),    intent(in)    :: forcing
      real(wp),                intent(in)    :: rho, press
      type(cas_boundary_t),    intent(inout) :: cas
      type(surface_env_t),     intent(inout) :: ground
      cas%aero_cfg     = aero_cfg
      cas%geom         = geom
      !----- Live inputs license the per-stage re-solve; a record not built here (a test fixture)    !
      !      leaves mo_live false, and its box conductances are used as given. ----------------------!
      cas%mo_live      = .true.
      cas%box%atm_enthalpy = forcing%enthalpy_atm
      cas%box%atm_specific_humidity = forcing%shv_atm
      cas%box%atm_co2 = forcing%co2_atm
      cas%rho   = rho
      cas%press = press
      ground%rnet_soil = forcing%rnet_soil
      ground%sw_snow   = forcing%sw_snow
      ground%lw_down   = forcing%lw_down_ground
      ground%ggnet     = aero%ggnet
   end subroutine canopy_air_boundary

   !---------------------------------------------------------------------------------------!
   ! transpiration_demand -- the start-of-step surface solve, at the state^n canopy air: each !
   ! cohort's transpiration demand on the water network (docs/science/soil_plant_water.md §4).!
   ! It takes the canopy air's exchange the pre-pass solved at that same state, so the        !
   ! surface layer is not solved again.                                                       !
   ! Each film's evaporation conductance is first scaled so its potential evaporation over    !
   ! dt_fast cannot exceed what the film holds plus this step's interception; then the leaf   !
   ! is credited the liquid enthalpy of the sapflow that replaces its transpiration, as the    !
   ! stages will, estimated as if the sapflow matched it, and the surface is solved again.    !
   !---------------------------------------------------------------------------------------!
   subroutine transpiration_demand(dt_fast, canopy_water_on, n, dynamic, t_ground, frozen)
      real(wp),                   intent(in)    :: dt_fast
      logical,                    intent(in)    :: canopy_water_on
      integer(ik),                intent(in)    :: n
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      real(wp),                   intent(in)    :: t_ground
      type(frozen_t),      intent(inout) :: frozen
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sf0
      type(cas_boundary_t)  :: cas_start   !< the boundary with the pre-pass's conductances
      real(wp)    :: t_up_wl
      integer(ik) :: i
      y_stage%cas_enthalpy = dynamic%cas%can_enthalpy
      y_stage%cas_shv = dynamic%cas%can_shv
      y_stage%cas_co2 = dynamic%cas%can_co2
      cas_start = frozen%cas
      cas_start%mo_live = .false.
      call surface_derivs(y_stage, cas_start, frozen%tissue, frozen%film, frozen%ground, t_ground, n, sf0)
      if (canopy_water_on)                                                                                    &
         call limit_film_evaporation(n, dt_fast, dynamic%leaf_surf_water, dynamic%wood_surf_water,              &
                                     sf0%film_evap_leaf, sf0%film_evap_wood, frozen%film)
      do i = 1_ik, n
         t_up_wl = merge(dynamic%wood_temp(i), dynamic%leaf_temp(i), sf0%transp_c(i) >= 0.0_wp)
         frozen%tissue%qwflux_wl(i) = sf0%transp_c(i) * internal_energy_liquid(t_up_wl)
      end do
      call surface_derivs(y_stage, cas_start, frozen%tissue, frozen%film, frozen%ground, t_ground, n, sf0)
      frozen%water_drive%e_demand(1:n) = sf0%transp_c(1:n)
   end subroutine transpiration_demand

   !---------------------------------------------------------------------------------------!
   ! ground_inflow -- what drives the water network besides the transpiration demand       !
   ! (docs/science/soil_plant_water.md §5, "The ground inflow"): the water reaching the    !
   ! pond, the temperature that values it and the enthalpy it brings across the column's   !
   ! boundary (pond_inflow), each layer's temperature and ice, and the ground evaporation's!
   ! boundary. Also the liquid enthalpy of the rain the canopy films exchange water with,  !
   ! which the start-of-step surface solves read, so it runs before them.                  !
   !---------------------------------------------------------------------------------------!
   subroutine ground_inflow(static, profile, forcing, dynamic, ggnet, throughfall, tcas, qcas, rho, t_ground, nsl, &
                            frozen)
      type(fast_static_state_t),  intent(in)    :: static
      type(soil_profile_t),       intent(in)    :: profile
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      real(wp),                   intent(in)    :: ggnet         !< [m/s] ground <-> canopy-air conductance
      real(wp),                   intent(in)    :: throughfall   !< [kg/m2/s]
      real(wp),                   intent(in)    :: tcas, qcas, rho, t_ground
      integer(ik),                intent(in)    :: nsl
      type(frozen_t),      intent(inout) :: frozen
      real(wp) :: t_water_base   !< [K] the meltwater's or the rain's, before the shed water joins it
      associate (snow => frozen%snow, drive => frozen%water_drive)
         call pond_inflow(snow%exists, snow%melt_rate, snow%t_melt, throughfall, forcing%rainfall, forcing%snowfall,  &
                          forcing%air_temp, tcas, static%shed_water_rate, static%shed_water_enth, drive%rain,         &
                          drive%t_rain, t_water_base, frozen%hydrology%ground_inflow_enth)
         drive%ground%snow_free_frac  = 1.0_wp - frozen%ground%snowfac
         drive%ground%t_ground        = t_ground
         drive%ground%q_air           = qcas
         drive%ground%rho_air         = rho
         drive%ground%press           = forcing%press
         !----- Bare-soil aerodynamic resistance; snow_free_frac above weights the flux, not this. -------!
         drive%ground%r_aero = 1.0_wp / max(ggnet, tiny_num)
         drive%t_soil(1:nsl) = dynamic%soil_energy%soil_temp(1:nsl)
         !----- Each layer's ice, held over the step (docs/science/soil_biophysics.md, "Frozen soil"). ----!
         drive%ice(1:nsl) = soil_ice_state(profile%layer_curve(1:nsl),                                            &
                                           (1.0_wp - dynamic%soil_energy%soil_fliq(1:nsl))                        &
                                           * dynamic%soil_water%theta(1:nsl))
         !----- The rain the canopy films capture and drip back is the throughfall's, at its temperature. -!
         frozen%film%rain_liquid_enthalpy = internal_energy_liquid(t_water_base)
      end associate
   end subroutine ground_inflow

   !----- The march's start state, the state^n canopy air, soil energy and water, plant water and  !
   !      films; the pond rides it (the stages pass it through; each scheme commits the network's). -!
   subroutine pack_march_start(dynamic, n, nsl, y)
      type(fast_dynamic_state_t), intent(in)  :: dynamic
      integer(ik),                intent(in)  :: n, nsl
      type(march_state_t),       intent(out) :: y
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n), y%leaf_surf_water(n), y%wood_surf_water(n))
      y%cas_enthalpy = dynamic%cas%can_enthalpy
      y%cas_shv = dynamic%cas%can_shv
      y%cas_co2 = dynamic%cas%can_co2
      y%soil_energy(1:nsl) = dynamic%soil_energy%soil_energy(1:nsl)
      y%theta(1:nsl)       = dynamic%soil_water%theta(1:nsl)
      y%w_surface          = dynamic%soil_water%w_surface
      y%w_surface_energy   = dynamic%soil_water%w_surface_energy
      y%leaf_water_mass(1:n) = dynamic%leaf_water_mass(1:n)
      y%wood_water_mass(1:n) = dynamic%wood_water_mass(1:n)
      y%leaf_surf_water(1:n) = dynamic%leaf_surf_water(1:n)
      y%wood_surf_water(1:n) = dynamic%wood_surf_water(1:n)
   end subroutine pack_march_start

   !---------------------------------------------------------------------------------------!
   ! record_water_flows -- what the water network moved over the step, as the frozen record's !
   ! rates (docs/science/soil_plant_water.md §4-§5): ground evaporation, the soil's boundary   !
   ! flows, its faces and end-of-step water, and the pond; the per-layer face closure; root    !
   ! uptake by layer and by cohort with the enthalpy it carries; sapflow and transpiration;    !
   ! the hydraulic diagnostics; and the plant water's advective enthalpy.                      !
   !---------------------------------------------------------------------------------------!
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
      integer(ik) :: i, k, n, nsl
      n   = static%n
      nsl = site%soil_profile%n_active
      associate (fl => water%flows, st => water%state)
         !----- Ground evaporation, the soil's boundary flows and the pond, as the network moved them. ----!
         frozen%ground%soil_evap         = fl%evap / dt_fast
         frozen%hydrology%infiltration   = fl%infil / dt_fast
         frozen%hydrology%exfiltration   = fl%exfil / dt_fast
         frozen%hydrology%exfil_enth     = fl%exfil_enth / dt_fast
         frozen%hydrology%q_top          = (fl%infil - fl%exfil - fl%evap) / (rho_h2o * dt_fast)
         frozen%hydrology%drainage       = fl%drain / dt_fast
         frozen%hydrology%runoff_surf    = fl%spill / dt_fast
         out%report%exports%runoff_surf_export = fl%spill
         out%report%exports%drainage_export    = fl%drain
         frozen%hydrology%runoff_enth    = fl%spill_enth / dt_fast
         !----- Infiltration leaves the pond, so the top face is valued at the temperature that carries the  !
         !      enthalpy the pond gave up (docs/science/numerical_scheme.md §5b, "The pond holds heat"). ----!
         frozen%hydrology%t_infil = frozen%water_drive%t_rain
         if (fl%infil > tiny_num) frozen%hydrology%t_infil = liquid_energy_to_temp(fl%infil_enth / fl%infil)
         frozen%hydrology%t_bot           = dynamic%soil_energy%soil_temp(nsl)
         frozen%hydrology%w_surface1      = st%water_soil(POND)
         frozen%hydrology%w_surface_energy1 = st%e_pond
         frozen%hydrology%w_flux_frozen(1:nsl)   = 0.0_wp
         frozen%hydrology%w_flux_frozen(1:nsl-1) = fl%face(1:nsl-1) / (rho_h2o * dt_fast)
         do k = 1_ik, nsl
            frozen%hydrology%theta1(k) = st%water_soil(1_ik + k) / (rho_h2o * site%soil_profile%dz(k))
         end do
         frozen%hydrology%theta_rate(1:nsl) = (frozen%hydrology%theta1(1:nsl) - dynamic%soil_water%theta(1:nsl))  &
                                              / dt_fast

         !----- Per-layer face closure: each layer's water change against the flows the soil energy will      !
         !      advect. It reads round-off, and checks that the faces handed on are the ones that moved the     !
         !      water; a network sub-step taken unconverged counts as a failed check. -------------------------!
         face_resid = 0.0_wp
         do k = 1_ik, nsl
            f_in = 0.0_wp
            f_out = 0.0_wp
            if (k > 1_ik)   f_in  = fl%face(k-1)
            if (k < nsl)    f_out = fl%face(k)
            if (k == nsl)   f_out = fl%drain
            if (k == 1_ik)  f_in  = fl%infil - fl%exfil - fl%evap
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

         !----- Root uptake in total, by layer and by cohort, and its enthalpy at each layer's start-of-step  !
         !      temperature (docs/science/soil_plant_water.md §5), the number the wood is credited below. A   !
         !      patch whose roots took nothing places its zero sink by its root profile. --------------------!
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

         !----- Per-cohort hydraulic diagnostics, from the same network the physics commits: sapflow and   !
         !      root uptake per plant, and the wood's end-of-step potential and loss of conductance. -------!
         if (out%capture) then
            do i = 1_ik, n
               out%cdiag(CO_SAPFLOW,     i) = frozen%plant%sapflow_frozen(i)
               out%cdiag(CO_ROOT_UPTAKE, i) = frozen%plant%uptake_frozen(i)
               out%cdiag(CO_PSI_WOOD,    i) = st%psi_coh(WOOD, i)
               associate (hp => cfg%pft%hydraulics%pft(static%pft(i)))
                  out%cdiag(CO_PLC, i) = 1.0_wp - conductance_retained_fraction(st%psi_coh(WOOD, i), hp%wood_psi50, &
                                                                                hp%wood_kexp)
               end associate
            end do
         end if

         !----- Advective enthalpy of the plant water (docs/science/soil_plant_water.md eqs. (4)-(5)): the      !
         !      network's flows, with upwind temperatures frozen at state^n, since the stage temperatures are   !
         !      computed from these terms; each tissue's water change carries its start-of-step enthalpy. -----!
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
   end subroutine record_water_flows

end module ed_fast_integ_frozen
