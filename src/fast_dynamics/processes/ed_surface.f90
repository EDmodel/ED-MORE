! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_surface -- the ground's surface: the bare-ground skin and the one-layer water store
! (snow, ice and the liquid they hold; the pond is its free water), blended by snowfac,
! each giving the soil-top boundary and the canopy-air fluxes
! (docs/science/snow_biophysics.md):
!   * ground_surface_exchange -- the blend at one canopy-air and soil-top state;
!   * ground_surface_fluxes -- the bare ground's sensible and latent fluxes;
!   * prepare_snow_stage -- the pack before a step: the pond it takes, its precipitation,
!     drainage and cover;
!   * snow_cover_fraction / snow_accumulate / snow_take_pond / snow_drain_meltwater -- the
!     layer's mass;
!   * snow_blend_ground_optics / ground_tile_radiation -- the ground's optics under a
!     pack, and the radiation each tile absorbs at its own;
!   * surface_node_balance -- one node's balance over the step, the pack's and the porous
!     layers' alike; snow_base_conductance / snow_take_exchange -- the pack's.
!---~---
module ed_surface
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan_boltzmann, &
                                       latent_heat_fusion, tsupercool_liq, rho_ice
   use ed_therm_lib,            only : internal_energy_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_liquid, precip_enthalpy, liquid_energy_to_temp
   use ed_patch_state_vars,     only : ground_store_t, water_layer_state_t, ground_store_present, store_move,  &
                                       store_take_exchange
   use ed_ground_surface_params, only : water_material_t
   use ed_ground_surface_lib,   only : water_layer_cover, water_layer_albedo, combined_density,             &
                                       water_layer_free_water, water_layer_resistance, water_layer_bulk_density
   implicit none
   private

   public :: ground_surface_exchange, ground_surface_fluxes, prepare_snow_stage
   public :: snow_cover_fraction, snow_accumulate, snow_take_pond, pond_freezes, snow_drain_meltwater
   public :: snow_blend_ground_optics
   public :: ground_tile_radiation
   public :: snow_take_exchange
   public :: surface_env_t, snow_out_t, surface_node_t, node_flux_t, store_exchange_t
   public :: STORE_SNOW, STORE_ORGANIC, N_SURFACE_STORE
   public :: surface_node_balance, conductance_in_series

   !---~---
   !   The ground surface's stores that exchange with the canopy air apart from the soil: the
   ! water layer (the snow pack) and the organic layer. Each one's exchange is one entry of a
   ! store_exchange_t, which every record of the march carries: the energy the store gains and
   ! the vapour it gives the canopy air, as rates in an evaluation and as amounts once the march
   ! has summed them. The update commits each store's amounts to it.
   !---~---
   integer(ik), parameter :: STORE_SNOW = 1_ik, STORE_ORGANIC = 2_ik, N_SURFACE_STORE = 2_ik
   type :: store_exchange_t
      real(wp) :: energy_in(N_SURFACE_STORE)  = 0.0_wp   ! [W/m2] or [J/m2]       energy into the store
      real(wp) :: vapour_out(N_SURFACE_STORE) = 0.0_wp   ! [kg/m2/s] or [kg/m2]   vapour from it to the canopy air
   end type store_exchange_t

   !---~---
   !   The ground surface's surroundings, read by the snow balance and the bare ground alike:
   ! the radiation each tile absorbs and its conductance to the canopy air, the canopy air's
   ! and the soil top's state, the thermal resistances below the surface stores, the bare
   ! ground's evaporation, and the pack it exchanges with, the pack's cover and parameters, and
   ! the step its balance relaxes over (prepare_snow_stage sets the pack and the cover).
   !---~---
   type :: surface_env_t
      real(wp) :: rnet_soil = 0.0_wp                    ! [W/m2] the bare soil's net all-wave radiation per unit of its area
      real(wp) :: sw_snow   = 0.0_wp                    ! [W/m2] the shortwave the snow absorbs per unit of its area
      real(wp) :: lw_down   = 0.0_wp                    ! [W/m2] the longwave reaching the ground
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp ! [K],[kg/kg] the canopy air's temperature and specific humidity
      real(wp) :: ggnet    = 0.0_wp                     ! [m/s] ground to canopy-air conductance, heat and vapour alike
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp ! [kg/m3],[Pa] the canopy air's density and pressure
      real(wp) :: cp_cas   = cp_air                     ! [J/kg/K] its specific heat (cas_heat_capacity); dry until set
      real(wp) :: t_soil_top = 273.16_wp                ! [K] top soil node's temperature
      !---~---
      !   The thermal resistances below the surface stores, set with the step's record
      ! (build_frozen): from the ground's surface down to the top soil node, its depth over the top
      ! layer's conductivity; and from the pack's base down to the node it conducts into, the
      ! organic layer's middle where there is one, else the top soil node.
      !---~---
      real(wp) :: r_soil_top = 0.0_wp                   ! [m2 K/W]
      real(wp) :: r_below    = 0.0_wp                   ! [m2 K/W]
      real(wp) :: soil_evap  = 0.0_wp                   ! [kg/m2/s] the bare ground's evaporation, already area-weighted
      !---~---
      !   The porous layers above the soil (docs/science/surface_organic_layer.md), one node: its
      ! dry heat capacity, water, energy and thermal resistance, and its emissivity, the ground's.
      ! With it, rnet_soil is its net radiation. Its evaporation is the water network's, held over
      ! the step as the soil's is, per unit ground area, and its vapour's resistance through its
      ! dried top is the network's input.
      !---~---
      logical  :: has_org        = .false.
      real(wp) :: org_dry_cap    = 0.0_wp               ! [J/m2/K]
      real(wp) :: org_water      = 0.0_wp               ! [kg/m2]
      real(wp) :: org_energy     = 0.0_wp               ! [J/m2]
      real(wp) :: org_dry_resist = 0.0_wp               ! [s/m] its vapour through its dried top (porous_layers_sums)
      real(wp) :: org_resistance = 0.0_wp               ! [m2 K/W] top to bottom
      real(wp) :: org_evap_held  = 0.0_wp               ! [kg/m2/s] its evaporation, per unit ground area (< 0: dew)
      real(wp) :: ground_emiss   = 0.95_wp              ! [-]
      type(water_layer_state_t) :: pack                 ! the pack after the step's precipitation and drainage
      real(wp) :: pack_density = 0.0_wp                 ! [kg/m3] its bulk density over the step (prepare_snow_stage)
      real(wp) :: snowfac    = 0.0_wp                   ! [-] its Niu-Yang cover fraction (0 without a pack)
      type(water_material_t) :: params                  ! the water layer's material
      real(wp) :: dt         = 0.0_wp                   ! [s] the step
   end type surface_env_t

   !---~---
   !   What the pack's stage returns, per unit ground area: the water it drained, the drainage as
   ! the pond's inflow over the step, and the precipitation's enthalpy, for the ledgers. All zero
   ! without a pack.
   !---~---
   type :: snow_out_t
      !---~---
      !   The water the pack drained, its free liquid at the plateau's energy, and all it held when
      ! it melted out, on its way to the pond: a store the pack's moves fill.
      !---~---
      type(ground_store_t) :: drained                   ! [kg/m2], [J/m2]
      !--- The stage's hand-off.
      logical  :: exists    = .false.                   ! a pack is present after the step's precipitation
      real(wp) :: melt_rate = 0.0_wp                    ! [kg/m2/s] meltwater and dump to the pond over the step
      !---~---
      !   The drained water's temperature, at which its liquid internal energy equals its
      ! enthalpy, so the pond receives exactly that enthalpy; t_3ple when none drains.
      !---~---
      real(wp) :: t_melt    = t_3ple                    ! [K]
      real(wp) :: acc_enth  = 0.0_wp                    ! [J/m2]  precipitation (rain and snow) enthalpy that entered the pack
   end type snow_out_t

   !---~---
   !   One node of the ground's surface, the pack or the porous layers, as its balance over a step
   ! reads it (surface_node_balance): its state at the step's start (its temperature and liquid
   ! fraction read off its energy), the share of the ground it exchanges with the canopy air
   ! over, the net radiation it absorbs per unit of that share at its temperature, its vapour,
   ! either the canopy air's (live, from its own saturation) or the water network's rate held
   ! over the step, per unit ground area; the heat the node above gives it, and its conductance
   ! to the node below, both per unit ground area.
   !---~---
   type :: surface_node_t
      real(wp) :: dry_cap     = 0.0_wp      ! [J/m2/K] heat capacity of its dry matter (none for water)
      real(wp) :: water       = 0.0_wp      ! [kg/m2]  liquid and ice
      real(wp) :: energy      = 0.0_wp      ! [J/m2]   internal energy
      real(wp) :: temp        = t_3ple      ! [K]
      real(wp) :: fliq        = 0.0_wp      ! [-]
      real(wp) :: area        = 0.0_wp      ! [-]      the share of the ground it exchanges over
      real(wp) :: rad         = 0.0_wp      ! [W/m2]   net all-wave radiation per unit of that share
      real(wp) :: emiss       = 1.0_wp      ! [-]      its emissivity
      logical  :: vapour_held = .false.     ! its vapour is the water network's
      real(wp) :: vapour_rate = 0.0_wp      ! [kg/m2/s] that rate (< 0: dew)
      real(wp) :: heat_above  = 0.0_wp      ! [W/m2]   from the node above, into it
      real(wp) :: g_below     = 0.0_wp      ! [W/m2/K] to the node below
      real(wp) :: t_below     = t_3ple      ! [K]      the node below's temperature
   end type surface_node_t

   !---~---
   !   What a node's balance returns, per unit ground area, positive away from the node: its net
   ! radiation, its sensible and latent heat and vapour to the canopy air, its conduction to the
   ! node below, and its net gain; and the temperature it ends at.
   !---~---
   type :: node_flux_t
      real(wp) :: t_surf    = t_3ple                    ! [K]
      real(wp) :: rad       = 0.0_wp                    ! [W/m2]
      real(wp) :: h         = 0.0_wp                    ! [W/m2]
      real(wp) :: le        = 0.0_wp                    ! [W/m2]
      real(wp) :: w_flux    = 0.0_wp                    ! [kg/m2/s]
      real(wp) :: g_below   = 0.0_wp                    ! [W/m2]
      real(wp) :: energy_in = 0.0_wp                    ! [W/m2]
   end type node_flux_t

contains

   !--- Ground / surface skin

   !---~---
   !   ground_surface_exchange -- the ground surface at one canopy-air and soil-top state: the
   ! pack over its cover, its balance relaxed over the step with its own heat store
   ! (docs/science/snow_biophysics.md, "The snow store"), conducting into the porous layers or
   ! the soil top; the porous layers, where there are any, over the rest, conducting into the
   ! soil top; and else the bare soil over the rest. Each node's balance is surface_node_balance,
   ! taken top down. The fluxes are zero for a node that is absent, so with no snow and no layer
   ! this is exactly bare ground. soil_evap already carries its snow-free area factor
   ! (ground_evap_env_t%snow_free_frac), so le_soil is added whole.
   !---~---
   pure subroutine ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores)
      type(surface_env_t),    intent(in)  :: env          !        the ground surface's surroundings
      real(wp),               intent(out) :: h_ground     ! [W/m2] sensible heat to the canopy air
      real(wp),               intent(out) :: le_ground    ! [W/m2] latent heat to the canopy air
      real(wp),               intent(out) :: g_top        ! [W/m2] net energy into the soil top
      real(wp),               intent(out) :: ground_rad   ! [W/m2] net radiation, snow's and bare soil's
      type(store_exchange_t), intent(out) :: stores       !        each surface store's exchange (0 without it)
      type(surface_node_t) :: org
      type(node_flux_t)    :: sfx, ofx
      real(wp) :: h_bare, le_soil, cover
      stores = store_exchange_t()
      cover  = max(0.0_wp, min(1.0_wp, env%snowfac))
      if (env%has_org) then
         !--- The porous layers are the ground's surface: the pack's base conducts into them, and they into the soil top.
         org%dry_cap = env%org_dry_cap
         org%water   = env%org_water
         org%energy  = env%org_energy
         call internal_energy_to_temp(org%energy, org%water, org%dry_cap, org%temp, org%fliq)
         call snow_node_balance(env, cover, org%temp, sfx)
         org%area        = 1.0_wp - cover
         org%rad         = env%rnet_soil
         org%emiss       = env%ground_emiss
         org%vapour_held = .true.
         org%vapour_rate = env%org_evap_held
         org%heat_above  = sfx%g_below
         org%g_below     = conductance_in_series(0.5_wp * env%org_resistance, env%r_soil_top)
         org%t_below     = env%t_soil_top
         call surface_node_balance(org, env, ofx)
         call ground_surface_fluxes(env%t_soil_top, env%can_temp, env%ggnet, env%rho_air, env%cp_cas,               &
                                    env%soil_evap, h_bare, le_soil)
         h_ground   = sfx%h  + ofx%h
         le_ground  = sfx%le + ofx%le + le_soil
         g_top      = ofx%g_below - le_soil
         ground_rad = sfx%rad + ofx%rad
         stores%energy_in(STORE_ORGANIC)  = ofx%energy_in
         stores%vapour_out(STORE_ORGANIC) = ofx%w_flux
      else
         call snow_node_balance(env, cover, env%t_soil_top, sfx)
         call ground_surface_fluxes(env%t_soil_top, env%can_temp, env%ggnet, env%rho_air, env%cp_cas,               &
                                    env%soil_evap, h_bare, le_soil)
         h_ground  = sfx%h  + (1.0_wp - env%snowfac) * h_bare
         le_ground = sfx%le + le_soil
         g_top     = sfx%g_below                                                                         &
                     + (1.0_wp - env%snowfac) * (env%rnet_soil - h_bare)                                 &
                     - le_soil
         ground_rad = sfx%rad + (1.0_wp - env%snowfac) * env%rnet_soil
      end if
      stores%energy_in(STORE_SNOW)  = sfx%energy_in
      stores%vapour_out(STORE_SNOW) = sfx%w_flux
   end subroutine ground_surface_exchange

   !---~---
   !   The pack as a node of the ground's surface, its balance over the step (surface_node_balance):
   ! over its cover, at its own temperature's net radiation, its vapour the canopy air's, its base
   ! conducting through the whole pack and r_below into the node below at t_below. With no pack
   ! or no cover it exchanges nothing, and its surface is the node below's.
   !---~---
   pure subroutine snow_node_balance(env, cover, t_below, flux)
      type(surface_env_t), intent(in)  :: env
      real(wp),            intent(in)  :: cover      ! [-] the pack's cover
      real(wp),            intent(in)  :: t_below    ! [K] the node its base conducts into
      type(node_flux_t),   intent(out) :: flux
      type(surface_node_t) :: pack
      flux = node_flux_t()
      flux%t_surf = t_below
      if (.not. ground_store_present(env%pack) .or. env%pack%water <= env%params%tiny_snow_mass                  &
          .or. cover <= 0.0_wp) return
      pack%water  = env%pack%water
      pack%energy = env%pack%energy
      call internal_energy_to_temp(pack%energy, pack%water, 0.0_wp, pack%temp, pack%fliq)
      pack%area    = cover
      pack%rad     = env%sw_snow + env%params%emissivity * (env%lw_down - stefan_boltzmann * pack%temp ** 4)
      pack%emiss   = env%params%emissivity
      pack%g_below = cover * snow_base_conductance(env%pack%matrix_depth, env, env%params)
      pack%t_below = t_below
      call surface_node_balance(pack, env, flux)
   end subroutine snow_node_balance

   !---~---
   !   Bare-ground fluxes to the canopy air, instantaneous: sensible through ggnet at the
   ! canopy air's density and specific heat, latent from the soil evaporation (already
   ! area-weighted) at the vapour enthalpy the canopy air receives. ground_surface_exchange
   ! assembles G_top and the snow blend.
   !---~---
   elemental pure subroutine ground_surface_fluxes(t_ground, t_cas, ggnet, rho, cp_cas, soil_evap, h_bare, le_soil)
      real(wp), intent(in)  :: t_ground, t_cas, ggnet       ! [K], [K], [m/s]
      real(wp), intent(in)  :: rho, cp_cas, soil_evap       ! [kg/m3], [J/kg/K], [kg/m2/s]
      real(wp), intent(out) :: h_bare, le_soil              ! [W/m2], [W/m2]
      h_bare  = sensible_heat_slope(ggnet, rho, cp_cas) * (t_ground - t_cas)
      le_soil = soil_evap * enthalpy_vapor(t_ground)
   end subroutine ground_surface_fluxes

   !---~---
   !   The sensible heat a ground surface store gives the canopy air per kelvin it is warmer
   ! [W/m2/K]: its conductance to the canopy air times the air's heat capacity per volume. The
   ! bare soil, the pack and the organic layer all exchange through it.
   !---~---
   elemental real(wp) function sensible_heat_slope(ggnet, rho, cp_cas) result(slope)
      real(wp), intent(in) :: ggnet, rho, cp_cas   ! [m/s], [kg/m3], [J/kg/K]
      slope = ggnet * rho * cp_cas
   end function sensible_heat_slope

   !--- [W/m2/K] How a surface's thermal emission grows with its temperature: d(eps sigma T^4)/dT.
   elemental real(wp) function emission_slope(emiss, temp) result(slope)
      real(wp), intent(in) :: emiss, temp   ! [-], [K]
      slope = 4.0_wp * emiss * stefan_boltzmann * temp ** 3
   end function emission_slope

   !--- [W/m2/K] The conductance of two thermal resistances [m2 K/W] in series.
   elemental real(wp) function conductance_in_series(r_upper, r_lower) result(g)
      real(wp), intent(in) :: r_upper, r_lower
      g = 1.0_wp / max(r_upper + r_lower, tiny_num)
   end function conductance_in_series

   !--- Snow -- mass side: cover fraction, accumulation, meltwater drainage.

   !---~---
   !   Niu-Yang (2007) snow-cover fraction [-], tanh(depth/scale): it ramps the ground optics
   ! and the soil-top blend and area-weights the pack's exchange
   ! (docs/science/snow_biophysics.md). The scale is set by the ground's roughness, the
   ! micro-relief a thin pack has to fill, and grows with the layer's density.
   !---~---
   pure function snow_cover_fraction(swe, snow_depth, density, params, ground_rough) result(snowfac)
      real(wp),            intent(in) :: swe, snow_depth   ! [kg/m2], [m]
      real(wp),            intent(in) :: density           ! [kg/m3] the layer's bulk density
      type(water_material_t), intent(in) :: params
      real(wp),            intent(in) :: ground_rough   ! [m] the ground's roughness ([aerodynamics].soil_rough)
      real(wp) :: snowfac
      if (swe <= params%tiny_snow_mass) then
         snowfac = 0.0_wp
      else
         snowfac = water_layer_cover(snow_depth, density, params%cover_scale_factor, params%cover_density_exponent, ground_rough)
      end if
   end function snow_cover_fraction

   !---~---
   !   The ground's optics under a pack, for the two-stream: its albedo and emissivity blended
   ! toward the pack's by its cover, the one the pack's balance exchanges over, so a partial pack
   ! is continuous, and the temperature at which the blended emissivity emits what the soil and
   ! the pack emit at their own (docs/science/snow_biophysics.md, "The ground's radiation").
   ! Unchanged without cover.
   !---~---
   pure subroutine snow_blend_ground_optics(snow, fc, params, albedo_vis, albedo_nir, emiss, temp)
      type(water_layer_state_t), intent(in) :: snow
      real(wp),            intent(in)    :: fc           ! [-] the pack's cover (prepare_snow_stage)
      type(water_material_t), intent(in)    :: params
      real(wp),            intent(inout) :: albedo_vis   ! [-] VIS albedo: the bare soil's on entry
      real(wp),            intent(inout) :: albedo_nir   ! [-] NIR albedo: the bare soil's on entry
      real(wp),            intent(inout) :: emiss        ! [-] thermal emissivity: the bare soil's on entry
      real(wp),            intent(inout) :: temp         ! [K] emission temperature: the top soil node's on entry
      real(wp) :: snow_vis, snow_nir, emit_soil, emit_snow
      if (.not. (fc > 0.0_wp)) return
      call snow_albedos(snow, params, snow_vis, snow_nir)
      emit_soil  = emiss * temp ** 4                                ! [K4] emission over sigma
      emit_snow  = params%emissivity * snow%temp ** 4
      albedo_vis = (1.0_wp - fc) * albedo_vis + fc * snow_vis
      albedo_nir = (1.0_wp - fc) * albedo_nir + fc * snow_nir
      emiss      = (1.0_wp - fc) * emiss + fc * params%emissivity
      temp       = (((1.0_wp - fc) * emit_soil + fc * emit_snow) / emiss) ** 0.25_wp
   end subroutine snow_blend_ground_optics

   !---~---
   !   The pack's own shortwave albedos, fresh to aged by its surface liquid fraction; a pack
   ! this step's snowfall has yet to make is fresh.
   !---~---
   pure subroutine snow_albedos(snow, params, albedo_vis, albedo_nir)
      type(water_layer_state_t), intent(in) :: snow
      type(water_material_t), intent(in)  :: params
      real(wp),            intent(out) :: albedo_vis, albedo_nir   ! [-]
      real(wp) :: fl
      fl = 0.0_wp
      if (ground_store_present(snow)) fl = snow%fliq
      albedo_vis = water_layer_albedo(fl, params%albedo_vis_fresh, params%albedo_vis_aged)
      albedo_nir = water_layer_albedo(fl, params%albedo_nir_fresh, params%albedo_nir_aged)
   end subroutine snow_albedos

   !---~---
   !   The radiation each tile of the ground absorbs per unit of its own area, from what
   ! reaches the ground in each band (docs/science/snow_biophysics.md, "The ground's
   ! radiation"): the bare soil's net all-wave radiation at its albedo and emissivity and the
   ! top soil node's temperature, and the shortwave the snow absorbs at its own albedo. The
   ! snow's longwave is formed in its balance, at its own temperature. Weighted by the cover,
   ! the two give the ground absorption of the two-stream, whose optics
   ! snow_blend_ground_optics blends.
   !---~---
   pure subroutine ground_tile_radiation(dn_vis, dn_nir, dn_lw, soil_albedo_vis, soil_albedo_nir, soil_emiss,  &
                                         soil_temp, snow, params, rnet_soil, sw_snow)
      real(wp),            intent(in)  :: dn_vis, dn_nir, dn_lw   ! [W/m2] each band reaching the ground
      real(wp),            intent(in)  :: soil_albedo_vis, soil_albedo_nir, soil_emiss   ! [-] the bare soil's optics
      real(wp),            intent(in)  :: soil_temp               ! [K] the top soil node's temperature
      type(water_layer_state_t), intent(in) :: snow
      type(water_material_t), intent(in)  :: params
      real(wp),            intent(out) :: rnet_soil, sw_snow      ! [W/m2]
      real(wp) :: snow_vis, snow_nir
      rnet_soil = dn_vis * (1.0_wp - soil_albedo_vis) + dn_nir * (1.0_wp - soil_albedo_nir)              &
                  + soil_emiss * (dn_lw - stefan_boltzmann * soil_temp ** 4)
      call snow_albedos(snow, params, snow_vis, snow_nir)
      sw_snow = dn_vis * (1.0_wp - snow_vis) + dn_nir * (1.0_wp - snow_nir)
   end subroutine ground_tile_radiation

   !---~---
   !   Snowfall joins the pack as ice and rain-on-snow as liquid, with their precip_enthalpy. A
   ! layer is created only from min_new_snow_mass of new snow; with no layer and less, the
   ! pack stays empty and the precipitation reaches the ground with the throughfall. The
   ! snowfall adds its depth at the fresh-snow density; rain fills the pores and adds none. The
   ! pack's density over the stage mixes by volume with the snowfall's.
   !---~---
   pure subroutine snow_accumulate(snow, density, snowfall, rainf, air_temp, dt, params)
      type(water_layer_state_t), intent(inout) :: snow               ! the pack
      real(wp),            intent(inout) :: density                  ! [kg/m3] its bulk density over the stage
      real(wp),            intent(in)    :: snowfall, rainf, air_temp, dt   ! [kg/m2/s], [kg/m2/s], [K], [s]
      type(water_material_t), intent(in)    :: params
      real(wp) :: dm_snow, dm_rain
      logical  :: has_layer
      dm_snow = max(0.0_wp, snowfall) * dt
      dm_rain = max(0.0_wp, rainf) * dt
      has_layer = ground_store_present(snow) .or. (dm_snow >= params%min_new_snow_mass)
      if (.not. has_layer) return                          ! bare ground and sub-threshold snow: no pack
      if (.not. ground_store_present(snow)) density = params%fresh_density
      if (density /= params%fresh_density)                                                        &
         density = combined_density(snow%water, density, dm_snow, params%fresh_density)
      snow%water        = snow%water + dm_snow + dm_rain
      snow%energy       = snow%energy + precip_enthalpy(dm_rain, dm_snow, air_temp)
      snow%matrix_depth = snow%matrix_depth + dm_snow / params%fresh_density   ! rain fills pores / refreezes -> no bulk depth
   end subroutine snow_accumulate

   !---~---
   !   Percolation: the free liquid above the holding capacity moves to the drained store, and so
   ! does all the pack holds when it melts out.
   !---~---
   pure subroutine snow_drain_meltwater(snow, density, params, melt)
      type(water_layer_state_t), intent(inout) :: snow     ! the pack
      real(wp),            intent(inout) :: density  ! [kg/m3] its bulk density over the stage
      type(water_material_t), intent(in)    :: params
      type(snow_out_t),   intent(out)   :: melt     ! the drained and dumped water, mass and enthalpy
      real(wp) :: wmass_free
      melt = snow_out_t()
      if (.not. ground_store_present(snow)) return

      !--- Vanished or trace layer: dump the residual to the pond and revert to bare ground.
      if (snow%water <= params%tiny_snow_mass) then
         call store_move(snow, melt%drained, snow%water, snow%energy)
         call clear_layer(snow, density)
         return
      end if

      call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)

      !--- Drain the free liquid above the holding capacity (LEAF-3's 1:9).
      wmass_free = water_layer_free_water(snow%water, snow%fliq, params%liquid_holding_frac)
      if (wmass_free > 0.0_wp) then
         call store_move(snow, melt%drained, wmass_free, wmass_free * internal_energy_liquid(t_3ple))  ! at the plateau
         snow%matrix_depth = snow%water / density
      end if

      !--- Below the melt-out threshold after draining: dump the remainder to the pond.
      if (snow%water <= params%tiny_snow_mass) then
         call store_move(snow, melt%drained, snow%water, snow%energy)
         call clear_layer(snow, density)
      else
         call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)
      end if
   end subroutine snow_drain_meltwater

   !--- Revert to the empty store: no water, no energy, no depth.
   pure subroutine clear_layer(snow, density)
      type(water_layer_state_t), intent(inout) :: snow
      real(wp),                  intent(out)   :: density   ! [kg/m3]
      type(water_layer_state_t) :: empty                     ! default-initialised and never written
      snow    = empty
      density = 0.0_wp
   end subroutine clear_layer

   !---~---
   !   Whether the pond is below its freezing point: its energy is less than its water's as
   ! liquid at the triple point, so it holds ice (a pond of sub-threshold snow, of freezing
   ! rain, or one that has lost heat).
   !---~---
   pure logical function pond_freezes(w_pond, e_pond, params)
      real(wp),            intent(in) :: w_pond, e_pond   ! [kg/m2], [J/m2]
      type(water_material_t), intent(in) :: params
      pond_freezes = w_pond > params%tiny_snow_mass .and. e_pond < w_pond * internal_energy_liquid(t_3ple)
   end function pond_freezes

   !---~---
   !   The pond joins the water layer, its water and energy; with no layer, it starts one, as
   ! ice. The pond is left empty; the layer's drainage gives back what it does not hold.
   !---~---
   pure subroutine snow_take_pond(snow, density, pond)
      type(water_layer_state_t), intent(inout) :: snow      ! the layer
      real(wp),             intent(inout) :: density        ! [kg/m3] its bulk density over the stage
      type(ground_store_t), intent(inout) :: pond           ! the pond
      if (.not. ground_store_present(snow)) then
         density           = rho_ice
         snow%matrix_depth = pond%water / rho_ice
      end if
      call store_move(pond, snow, pond%water, pond%energy)
   end subroutine snow_take_pond

   !---~---
   !   prepare_snow_stage -- the water layer before a step (docs/science/snow_biophysics.md, "The
   ! snow store"): the pond joins it, the step's snowfall and rain-on-snow join it, and the
   ! liquid it holds beyond its capacity drains, a paired (mass, enthalpy) transfer, with its
   ! cover fraction. A pond below its freezing point starts the layer, as ice, before the
   ! precipitation, which then joins it; otherwise the pond joins a pack the precipitation leaves.
   ! The pond joins by freezing: the layer keeps no more liquid than it held before, so what it
   ! keeps of the pond is what its cold freezes, and that adds its depth as ice. The drained water
   ! refills the pond first, up to what it held, and the rest reaches it over the step as meltwater. env's pack
   ! enters as the pack at the step's start and leaves advanced, with env's cover set; the pond
   ! enters as the step's start and leaves as the pond the step starts from; st returns the
   ! drainage and the stage's hand-off. The exchange is then with this pack
   ! (ground_surface_exchange), and the pack takes it afterwards (snow_take_exchange).
   !---~---
   pure subroutine prepare_snow_stage(snowfall, rainfall, t_air, ground_rough, env, st, pond)
      real(wp),            intent(in)    :: snowfall       ! [kg/m2/s] frozen precipitation
      real(wp),            intent(in)    :: rainfall       ! [kg/m2/s] liquid precipitation (rain-on-snow)
      real(wp),            intent(in)    :: t_air          ! [K]       reference-level air temperature
      real(wp),            intent(in)    :: ground_rough   ! [m]       the ground's roughness ([aerodynamics].soil_rough)
      type(surface_env_t), intent(inout) :: env            ! its pack advanced, its cover set; params and dt read
      type(snow_out_t),    intent(out)   :: st             !           the drainage and the stage's hand-off
      type(ground_store_t), intent(inout) :: pond          !           the pond: its water and energy
      type(ground_store_t) :: pond0                        !           the pond at the stage's start
      real(wp) :: e_before, acc_enth, pond_kept
      real(wp) :: back, back_enth, liq_before, t_now, fl_now, extra
      logical  :: exists, took_pond, joined

      env%pack_density = water_layer_bulk_density(env%pack%water, env%pack%matrix_depth)
      pond0     = pond
      took_pond = .false.
      joined    = .false.
      if (.not. ground_store_present(env%pack) .and. pond_freezes(pond%water, pond%energy, env%params)) then
         call snow_take_pond(env%pack, env%pack_density, pond)   ! the frozen pond starts the layer
         took_pond = .true.
      end if
      e_before = env%pack%energy
      call snow_accumulate(env%pack, env%pack_density, snowfall, rainfall, t_air, env%dt, env%params)
      acc_enth = env%pack%energy - e_before              ! precipitation (rain and snow) enthalpy into the pack
      exists   = ground_store_present(env%pack)          ! a pack took the snow and the rain
      liq_before = 0.0_wp                                ! the liquid the layer held before a pond joined it
      if (exists .and. .not. took_pond .and. pond%water > 0.0_wp) then
         call internal_energy_to_temp(env%pack%energy, env%pack%water, 0.0_wp, t_now, fl_now)
         liq_before = fl_now * env%pack%water
         call snow_take_pond(env%pack, env%pack_density, pond)   ! the pond joins the pack
         took_pond = .true.
         joined    = .true.
      end if
      env%snowfac = 0.0_wp

      if (exists) then
         !---~---
         !   Meltwater goes to the pond with its enthalpy, valued at t_melt, and nowhere else; pack
         ! and pond are both tracked stores, so the transfer cancels in the column's ledger.
         !---~---
         call snow_drain_meltwater(env%pack, env%pack_density, env%params, st)
         if (took_pond .and. ground_store_present(env%pack)) then
            !---~---
            !   The pond joins the layer by freezing: its liquid does not rise into the pores, so the
            ! layer keeps no more liquid than it held before; the rest, at the plateau, drains.
            !---~---
            call internal_energy_to_temp(env%pack%energy, env%pack%water, 0.0_wp, t_now, fl_now)
            extra = fl_now * env%pack%water - liq_before
            if (extra > 0.0_wp) then
               call store_move(env%pack, st%drained, extra, extra * internal_energy_liquid(t_3ple))
               env%pack%matrix_depth = env%pack%water / env%pack_density
               if (env%pack%water <= env%params%tiny_snow_mass) then
                  call store_move(env%pack, st%drained, env%pack%water, env%pack%energy)
                  call clear_layer(env%pack, env%pack_density)
               else
                  call internal_energy_to_temp(env%pack%energy, env%pack%water, 0.0_wp, env%pack%temp, env%pack%fliq)
               end if
            end if
         end if
         if (took_pond) then
            !--- The drained water, at one temperature, refills the pond up to what it held.
            back      = min(st%drained%water, pond0%water)
            back_enth = 0.0_wp
            if (back > 0.0_wp) back_enth = st%drained%energy * (back / st%drained%water)
            call store_move(st%drained, pond, back, back_enth)
            pond_kept = pond0%water - back
            !--- What a pack keeps of the pond that joined it adds its depth as ice.
            if (joined .and. ground_store_present(env%pack) .and. pond_kept > 0.0_wp) then
               env%pack_density      = combined_density(env%pack%water - pond_kept, env%pack_density,      &
                                                        pond_kept, rho_ice)
               env%pack%matrix_depth = env%pack%water / env%pack_density
            end if
         end if
         st%melt_rate = st%drained%water / env%dt
         st%t_melt    = t_3ple
         if (st%drained%water > 0.0_wp) st%t_melt = liquid_energy_to_temp(st%drained%energy / st%drained%water)
         !---~---
         !   snowfac is snow, 1 - snowfac bare soil; the pack's exchange scales by it, so a patchy
         ! pack barely exchanges and its fluxes come snowfac-weighted.
         !---~---
         if (env%pack%water > env%params%tiny_snow_mass)                                               &
            env%snowfac = snow_cover_fraction(env%pack%water, env%pack%matrix_depth, env%pack_density,      &
                                              env%params, ground_rough)
      end if
      st%exists   = exists
      st%acc_enth = acc_enth
   end subroutine prepare_snow_stage

   !--- Snow -- energy side: surface fluxes, base conduction, backward-Euler energy advance.

   !---~---
   !   The pack's conductance to the node below it [W/m2/K]: the whole snow depth and r_below in
   ! series. The pack's one temperature is its surface's, the temperature it emits and
   ! exchanges heat and vapour at, so the heat between it and the ground crosses the whole pack. A
   ! deepening pack decouples the soil (docs/science/snow_biophysics.md).
   !---~---
   pure function snow_base_conductance(snow_depth, env, params) result(gcond)
      real(wp),            intent(in) :: snow_depth   ! [m]
      type(surface_env_t), intent(in) :: env
      type(water_material_t), intent(in) :: params
      real(wp) :: gcond
      gcond = conductance_in_series(water_layer_resistance(snow_depth, snow_depth, params%conductivity_at_fresh), env%r_below)
   end function snow_base_conductance

   !---~---
   !   One node's exchange over dt, per unit ground area, from a linearized backward-Euler balance
   ! of its energy that is not committed (docs/science/snow_biophysics.md, "The snow store";
   ! docs/science/surface_organic_layer.md, Heat): the balance's one root T, in closed form on the
   ! melt plateau or else on the ice or liquid line, with its dry heat capacity on each, and the
   ! fluxes at T. Over its share of the ground it takes its net radiation and gives sensible heat
   ! and vapour, a kilogram of which carries enthalpy_vapor less the water's own energy; it takes
   ! the heat from the node above and conducts into the node below. Implicit in its own
   ! temperature, so a node of little heat capacity relaxes without overshooting. Its own vapour,
   ! once the linearization takes it past the water it holds, is held there and the balance
   ! solved again, so a node that loses all its water to the air ends with none of its energy;
   ! a held vapour is the water network's, per unit ground area. The node takes the exchange's
   ! time integral in its own commit.
   !---~---
   pure subroutine surface_node_balance(node, env, flux)
      type(surface_node_t), intent(in)  :: node
      type(surface_env_t),  intent(in)  :: env      ! the canopy air, the conductance to it and the step
      type(node_flux_t),    intent(out) :: flux
      real(wp) :: a, dt, w, cd, e_old, t_n, fliq_n, h_n, dh_dt, w_n, dw_dt, vh, demit_dt, u_w, l_vap, g_n
      real(wp) :: r0, drdt, u_freeze, u_melt, e_bal, t_new, fl_new, dtemp, w_s, w_cap, h_s, g_s, emit, u_w_new
      integer  :: pass
      flux  = node_flux_t()
      a     = node%area
      dt    = env%dt
      w     = node%water
      cd    = node%dry_cap
      e_old = node%energy
      t_n   = node%temp
      fliq_n = node%fliq

      !--- The fluxes at T^n and their slopes; a kilogram of vapour takes its enthalpy less the water's energy.
      dh_dt = sensible_heat_slope(env%ggnet, env%rho_air, env%cp_cas)
      h_n   = dh_dt * (t_n - env%can_temp)
      w_n   = 0.0_wp
      dw_dt = 0.0_wp
      vh    = 0.0_wp
      if (node%vapour_held) then
         vh = node%vapour_rate
      else
         w_n   = env%ggnet * env%rho_air * (sat_specific_humidity(t_n, env%press, fliq_n) - env%can_shv)
         dw_dt = env%ggnet * env%rho_air * sat_specific_humidity_temp_deriv(t_n, env%press, fliq_n)
      end if
      demit_dt = emission_slope(node%emiss, t_n)
      u_w   = water_specific_energy(t_n, fliq_n)                ! from the datum: a trace of water stays exact
      l_vap = enthalpy_vapor(t_n) - u_w
      g_n   = node%g_below * (t_n - node%t_below)
      w_cap = huge(1.0_wp)
      if (a > 0.0_wp .and. .not. node%vapour_held) w_cap = w / (dt * a)   ! the vapour that empties it over dt

      !--- Solve e(T) = E^n + dt (r0 + drdt (T - T^n)): on the plateau first, else on the ice or the liquid line.
      do pass = 1, 2
         r0    = a * (node%rad - h_n - l_vap * w_n) - l_vap * vh + node%heat_above - g_n
         drdt  = -(a * (demit_dt + dh_dt + l_vap * dw_dt) + node%g_below)   ! < 0
         u_freeze = (cd + w * cp_ice) * t_3ple
         u_melt   = u_freeze + w * latent_heat_fusion
         t_new    = t_3ple
         e_bal    = e_old + dt * (r0 + drdt * (t_3ple - t_n))
         fl_new   = min(1.0_wp, max(0.0_wp, (e_bal - u_freeze) / max(w * latent_heat_fusion, tiny_num)))
         if (e_bal < u_freeze) then
            t_new  = (e_old + dt * (r0 - drdt * t_n)) / (cd + w * cp_ice - dt * drdt)
            fl_new = 0.0_wp
         else if (e_bal > u_melt) then
            t_new  = (e_old + w * cp_liq * tsupercool_liq + dt * (r0 - drdt * t_n)) / (cd + w * cp_liq - dt * drdt)
            fl_new = 1.0_wp
         end if
         dtemp = t_new - t_n
         w_s   = w_n + dw_dt * dtemp
         if (w_s <= w_cap) exit
         w_n   = w_cap                                          ! held, and solved again
         dw_dt = 0.0_wp
      end do
      w_s = min(w_s, w_cap)

      !--- The fluxes at T, on the linearization.
      h_s     = h_n + dh_dt * dtemp
      g_s     = g_n + node%g_below * dtemp
      emit    = demit_dt * dtemp                                ! the emission added as the node moves from T^n to T
      u_w_new = water_specific_energy(t_new, fl_new)
      flux%t_surf    = t_new
      flux%rad       = a * (node%rad - emit)
      flux%h         = a * h_s
      flux%w_flux    = a * w_s + vh
      flux%le        = flux%w_flux * (l_vap + u_w_new)           ! the latent heat plus the leaving water's own energy
      flux%g_below   = g_s
      flux%energy_in = flux%rad - flux%h - flux%le + node%heat_above - g_s
   end subroutine surface_node_balance

   !--- [J/kg] Water's specific internal energy at a temperature and liquid fraction (the inverter's datum).
   elemental function water_specific_energy(temp, fliq) result(u)
      real(wp), intent(in) :: temp, fliq   ! [K], [-]
      real(wp) :: u
      u = fliq * internal_energy_liquid(temp) + (1.0_wp - fliq) * cp_ice * temp
   end function water_specific_energy

   !---~---
   !   The pack takes an exchange's amounts, energy in and vapour out, keeping e(T) per
   ! kilogram, so its temperature follows from the new energy. A pack whose mass is gone is
   ! left at t_3ple for snow_drain_meltwater to empty.
   !---~---
   pure subroutine snow_take_exchange(snow, density, energy_amount, vapour_amount, params)
      type(water_layer_state_t), intent(inout) :: snow                      ! the pack
      real(wp),            intent(in)    :: density                         ! [kg/m3] its bulk density over the stage
      real(wp),            intent(in)    :: energy_amount, vapour_amount    ! [J/m2] in, [kg/m2] out
      type(water_material_t), intent(in) :: params
      if (.not. ground_store_present(snow)) return
      call store_take_exchange(snow, energy_amount, vapour_amount)
      snow%matrix_depth = snow%water / density
      if (snow%water > params%tiny_snow_mass) then
         call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)
      else
         snow%temp = t_3ple                               ! vanished; snow_drain_meltwater dumps the rest
         snow%fliq = 0.0_wp
      end if
   end subroutine snow_take_exchange

end module ed_surface
