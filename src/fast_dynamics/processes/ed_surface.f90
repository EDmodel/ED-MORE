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
!   * snow_surface_fluxes / snow_base_conductance / snow_surface_balance /
!     snow_take_exchange -- its energy.
!---~---
module ed_surface
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan_boltzmann, &
                                       latent_heat_fusion, tsupercool_liq, rho_ice
   use ed_therm_lib,            only : internal_energy_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_liquid, precip_enthalpy, liquid_energy_to_temp
   use ed_patch_state_vars,     only : water_layer_state_t, ground_store_present
   use ed_ground_surface_params, only : water_material_t
   use ed_numerics,             only : within_roundoff
   use ed_ground_surface_lib,   only : water_layer_cover, water_layer_albedo, combined_density,             &
                                       water_layer_free_water, water_layer_resistance, water_layer_bulk_density
   implicit none
   private

   public :: ground_surface_exchange, ground_surface_fluxes, prepare_snow_stage
   public :: snow_cover_fraction, snow_accumulate, snow_take_pond, pond_freezes, snow_drain_meltwater
   public :: snow_blend_ground_optics
   public :: ground_tile_radiation
   public :: snow_take_exchange
   public :: surface_env_t, snow_out_t, organic_out_t, store_exchange_t
   public :: STORE_SNOW, STORE_ORGANIC, N_SURFACE_STORE
   public :: organic_surface_balance, conductance_in_series

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
      !   The organic layer above the soil (docs/science/surface_organic_layer.md), one node:
      ! its dry heat capacity, water, energy, thermal resistance and its vapour's resistance
      ! through its dried top at the step's start, and its emissivity, the ground's. With it,
      ! rnet_soil is its net radiation. Once the water network has taken the layer's own
      ! evaporation (org_held), the balance holds that rate, as the soil's is held.
      !---~---
      logical  :: has_org        = .false.
      real(wp) :: org_dry_cap    = 0.0_wp               ! [J/m2/K]
      real(wp) :: org_water      = 0.0_wp               ! [kg/m2]
      real(wp) :: org_energy     = 0.0_wp               ! [J/m2]
      real(wp) :: org_dry_resist = 0.0_wp               ! [s/m] its vapour through its dried top (organic_layer_sums)
      real(wp) :: org_resistance = 0.0_wp               ! [m2 K/W] top to bottom
      logical  :: org_held       = .false.              ! its evaporation is the water network's
      real(wp) :: org_evap_held  = 0.0_wp               ! [kg/m2/s] that rate, per unit ground area (< 0: dew)
      real(wp) :: ground_emiss   = 0.95_wp              ! [-]
      type(water_layer_state_t) :: pack                 ! the pack after the step's precipitation and drainage
      real(wp) :: pack_density = 0.0_wp                 ! [kg/m3] its bulk density over the step (prepare_snow_stage)
      real(wp) :: snowfac    = 0.0_wp                   ! [-] its Niu-Yang cover fraction (0 without a pack)
      type(water_material_t) :: params                  ! the water layer's material
      real(wp) :: dt         = 0.0_wp                   ! [s] the step
   end type surface_env_t

   !---~---
   !   What the snow kernels return, per unit ground area. snow_surface_balance fills the
   ! exchange (positive upward, away from the pack) and snow_drain_meltwater the meltwater
   ! handed to the pond; each resets the whole record, so the other's fields read zero.
   ! prepare_snow_stage returns the drainage with the stage's hand-off: the meltwater as the
   ! pond's inflow, and the precipitation's enthalpy and the pack before the stage, for the
   ! ledgers. All zero without a pack.
   !---~---
   type :: snow_out_t
      !--- The surface balance.
      real(wp) :: t_surf  = 273.16_wp                   ! [K] diagnosed snow-surface temperature
      real(wp) :: h_snow  = 0.0_wp, le_snow = 0.0_wp    ! [W/m2] sensible + latent (sublimation/evap) to the canopy air
      real(wp) :: w_flux  = 0.0_wp                      ! [kg/m2/s] vapour to the canopy air (>0 sublimation, <0 deposition)
      real(wp) :: g_base  = 0.0_wp                      ! [W/m2] conduction from the snow base into the soil top
      real(wp) :: rnet    = 0.0_wp                      ! [W/m2] net all-wave radiation absorbed by the snow surface
      real(wp) :: energy_in = 0.0_wp                    ! [W/m2] what the pack gains: rnet - h_snow - le_snow - g_base
      real(wp) :: snowfac = 0.0_wp                      ! [-] Niu-Yang07 snow-cover / burial fraction
      !--- The meltwater drainage, a store-to-store hand-off to the pond.
      real(wp) :: melt_mass = 0.0_wp                    ! [kg/m2] free liquid drained this step -> the pond
      real(wp) :: melt_enth = 0.0_wp                    ! [J/m2]  matching enthalpy (internal_energy_liquid(t_3ple))
      real(wp) :: dump_mass = 0.0_wp                    ! [kg/m2] full-melt-out residual water dumped to the pond
      real(wp) :: dump_enth = 0.0_wp                    ! [J/m2]  full-melt-out residual energy dumped to the pond
      !--- The stage's hand-off.
      logical  :: exists    = .false.                   ! a pack is present after the step's precipitation
      real(wp) :: melt_rate = 0.0_wp                    ! [kg/m2/s] meltwater and dump to the pond over the step
      !---~---
      !   The drained water's temperature, at which its liquid internal energy equals its
      ! enthalpy, so the pond receives exactly that enthalpy; t_3ple when none drains.
      !---~---
      real(wp) :: t_melt    = t_3ple                    ! [K]
      real(wp) :: acc_enth  = 0.0_wp                    ! [J/m2]  precipitation (rain and snow) enthalpy that entered the pack
      !---~---
      !   The pack before the stage and the pond water it took, for the ledgers: the pond the
      ! step starts from is what the pack left it, so the two together are the stores at the
      ! step's start.
      !---~---
      real(wp) :: swe0      = 0.0_wp                    ! [kg/m2]
      real(wp) :: enth0     = 0.0_wp                    ! [J/m2]  internal energy
      real(wp) :: pond_mass = 0.0_wp                    ! [kg/m2] pond water the pack took and kept
      real(wp) :: pond_enth = 0.0_wp                    ! [J/m2]  its energy
   end type snow_out_t

   !---~---
   !   What the organic layer's balance returns, per unit ground area, positive away from it: its
   ! net radiation, sensible and latent heat and vapour on the snow-free ground, the conduction
   ! into the soil top, and its net gain.
   !---~---
   type :: organic_out_t
      real(wp) :: t_surf    = 273.16_wp                 ! [K]
      real(wp) :: rad       = 0.0_wp                    ! [W/m2]
      real(wp) :: h         = 0.0_wp                    ! [W/m2]
      real(wp) :: le        = 0.0_wp                    ! [W/m2]
      real(wp) :: w_flux    = 0.0_wp                    ! [kg/m2/s]
      real(wp) :: g_soil    = 0.0_wp                    ! [W/m2]
      real(wp) :: energy_in = 0.0_wp                    ! [W/m2]
   end type organic_out_t

contains

   !--- Ground / surface skin

   !---~---
   !   ground_surface_exchange -- the ground surface at one canopy-air and soil-top state:
   ! the snowfac-weighted pack, its balance relaxed over the step with its own heat store
   ! (docs/science/snow_biophysics.md, "The snow store"), and the bare soil on the rest.
   ! The fluxes come snowfac-weighted and are zero with no pack or no cover, so with no
   ! snow this is exactly bare ground, with no branch. soil_evap already carries its
   ! snow-free area factor (ground_evap_env_t%snow_free_frac), so le_soil is added whole.
   !---~---
   pure subroutine ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores)
      type(surface_env_t),    intent(in)  :: env          !        the ground surface's surroundings
      real(wp),               intent(out) :: h_ground     ! [W/m2] sensible heat to the canopy air
      real(wp),               intent(out) :: le_ground    ! [W/m2] latent heat to the canopy air
      real(wp),               intent(out) :: g_top        ! [W/m2] net energy into the soil top
      real(wp),               intent(out) :: ground_rad   ! [W/m2] net radiation, snow's and bare soil's
      type(store_exchange_t), intent(out) :: stores       !        each surface store's exchange (0 without it)
      type(snow_out_t)    :: sfx
      type(organic_out_t) :: ofx
      real(wp) :: h_bare, le_soil, t_org, fl_org
      stores = store_exchange_t()
      if (env%has_org) then
         !---~---
         !   The organic layer is the ground's surface: the pack's base conducts into it, at its
         ! temperature at the step's start, and it into the soil top.
         !---~---
         call internal_energy_to_temp(env%org_energy, env%org_water, env%org_dry_cap, t_org, fl_org)
         call snow_surface_balance(env, t_org, sfx)
         call organic_surface_balance(env, sfx%g_base, ofx)
         call ground_surface_fluxes(env%t_soil_top, env%can_temp, env%ggnet, env%rho_air, env%cp_cas,               &
                                    env%soil_evap, h_bare, le_soil)
         h_ground   = sfx%h_snow  + ofx%h
         le_ground  = sfx%le_snow + ofx%le + le_soil
         g_top      = ofx%g_soil - le_soil
         ground_rad = sfx%rnet + ofx%rad
         stores%energy_in(STORE_ORGANIC)  = ofx%energy_in
         stores%vapour_out(STORE_ORGANIC) = ofx%w_flux
      else
         call snow_surface_balance(env, env%t_soil_top, sfx)
         call ground_surface_fluxes(env%t_soil_top, env%can_temp, env%ggnet, env%rho_air, env%cp_cas,               &
                                    env%soil_evap, h_bare, le_soil)
         h_ground  = sfx%h_snow  + (1.0_wp - env%snowfac) * h_bare
         le_ground = sfx%le_snow + le_soil
         g_top     = sfx%g_base                                                                         &
                     + (1.0_wp - env%snowfac) * (env%rnet_soil - h_bare)                                 &
                     - le_soil
         ground_rad = sfx%rnet + (1.0_wp - env%snowfac) * env%rnet_soil
      end if
      stores%energy_in(STORE_SNOW)  = sfx%energy_in
      stores%vapour_out(STORE_SNOW) = sfx%w_flux
   end subroutine ground_surface_exchange

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
   ! toward the pack's by the snow-cover fraction, so a partial pack is continuous, and the
   ! temperature at which the blended emissivity emits what the soil and the pack emit at
   ! their own (docs/science/snow_biophysics.md, "The ground's radiation"). Unchanged without
   ! a pack.
   !---~---
   pure subroutine snow_blend_ground_optics(snow, params, ground_rough, albedo_vis, albedo_nir, emiss, temp)
      type(water_layer_state_t), intent(in) :: snow
      type(water_material_t), intent(in)    :: params
      real(wp),            intent(in)    :: ground_rough   ! [m] the ground's roughness ([aerodynamics].soil_rough)
      real(wp),            intent(inout) :: albedo_vis   ! [-] VIS albedo: the bare soil's on entry
      real(wp),            intent(inout) :: albedo_nir   ! [-] NIR albedo: the bare soil's on entry
      real(wp),            intent(inout) :: emiss        ! [-] thermal emissivity: the bare soil's on entry
      real(wp),            intent(inout) :: temp         ! [K] emission temperature: the top soil node's on entry
      real(wp) :: fc, snow_vis, snow_nir, emit_soil, emit_snow
      if (.not. (snow%water > params%tiny_snow_mass)) return
      fc = snow_cover_fraction(snow%water, snow%matrix_depth, water_layer_bulk_density(snow%water, snow%matrix_depth), &
                               params, ground_rough)
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
   !   Percolation: drain the free liquid above the holding capacity as a paired (mass,
   ! enthalpy) hand-off to the pond (melt_*), and dump the residual at melt-out (dump_*).
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
         melt%dump_mass  = max(0.0_wp, snow%water)
         melt%dump_enth  = snow%energy
         call clear_layer(snow, density)
         return
      end if

      call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)

      !--- Drain the free liquid above the holding capacity (LEAF-3's 1:9).
      wmass_free = water_layer_free_water(snow%water, snow%fliq, params%liquid_holding_frac)
      if (wmass_free > 0.0_wp) then
         melt%melt_mass    = wmass_free
         melt%melt_enth    = wmass_free * internal_energy_liquid(t_3ple)   ! drains at the plateau
         snow%water        = snow%water - wmass_free
         snow%energy       = snow%energy - melt%melt_enth
         snow%matrix_depth = snow%water / density
      end if

      !--- Below the melt-out threshold after draining: dump the remainder to the pond.
      if (snow%water <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%water)
         melt%dump_enth  = snow%energy
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
   pure subroutine snow_take_pond(snow, density, w_pond, e_pond)
      type(water_layer_state_t), intent(inout) :: snow      ! the layer
      real(wp),           intent(inout) :: density          ! [kg/m3] its bulk density over the stage
      real(wp),           intent(inout) :: w_pond, e_pond   ! [kg/m2], [J/m2] the pond
      if (.not. ground_store_present(snow)) then
         density           = rho_ice
         snow%matrix_depth = w_pond / rho_ice
      end if
      snow%water  = snow%water + w_pond
      snow%energy = snow%energy + e_pond
      w_pond = 0.0_wp
      e_pond = 0.0_wp
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
   pure subroutine prepare_snow_stage(snowfall, rainfall, t_air, ground_rough, env, st, w_pond, e_pond)
      real(wp),            intent(in)    :: snowfall       ! [kg/m2/s] frozen precipitation
      real(wp),            intent(in)    :: rainfall       ! [kg/m2/s] liquid precipitation (rain-on-snow)
      real(wp),            intent(in)    :: t_air          ! [K]       reference-level air temperature
      real(wp),            intent(in)    :: ground_rough   ! [m]       the ground's roughness ([aerodynamics].soil_rough)
      type(surface_env_t), intent(inout) :: env            ! its pack advanced, its cover set; params and dt read
      type(snow_out_t),    intent(out)   :: st             !           the drainage and the stage's hand-off
      real(wp),            intent(inout) :: w_pond         ! [kg/m2]   the pond's water
      real(wp),            intent(inout) :: e_pond         ! [J/m2]    the pond's energy
      real(wp) :: swe0, enth0, e_before, acc_enth, melt_enth, w_pond0, e_pond0, drained, drained_enth
      real(wp) :: back, back_enth, liq_before, t_now, fl_now, extra
      logical  :: exists, took_pond, joined

      swe0      = env%pack%water
      enth0     = env%pack%energy
      env%pack_density = water_layer_bulk_density(env%pack%water, env%pack%matrix_depth)
      w_pond0   = w_pond
      e_pond0   = e_pond
      took_pond = .false.
      joined    = .false.
      if (.not. ground_store_present(env%pack) .and. pond_freezes(w_pond, e_pond, env%params)) then
         call snow_take_pond(env%pack, env%pack_density, w_pond, e_pond)   ! the frozen pond starts the layer
         took_pond = .true.
      end if
      e_before = env%pack%energy
      call snow_accumulate(env%pack, env%pack_density, snowfall, rainfall, t_air, env%dt, env%params)
      acc_enth = env%pack%energy - e_before              ! precipitation (rain and snow) enthalpy into the pack
      exists   = ground_store_present(env%pack)          ! a pack took the snow and the rain
      liq_before = 0.0_wp                                ! the liquid the layer held before a pond joined it
      if (exists .and. .not. took_pond .and. w_pond > 0.0_wp) then
         call internal_energy_to_temp(env%pack%energy, env%pack%water, 0.0_wp, t_now, fl_now)
         liq_before = fl_now * env%pack%water
         call snow_take_pond(env%pack, env%pack_density, w_pond, e_pond)   ! the pond joins the pack
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
               st%melt_mass = st%melt_mass + extra
               st%melt_enth = st%melt_enth + extra * internal_energy_liquid(t_3ple)
               env%pack%water        = env%pack%water - extra
               env%pack%energy       = env%pack%energy - extra * internal_energy_liquid(t_3ple)
               env%pack%matrix_depth = env%pack%water / env%pack_density
               if (env%pack%water <= env%params%tiny_snow_mass) then
                  st%dump_mass = st%dump_mass + max(0.0_wp, env%pack%water)
                  st%dump_enth = st%dump_enth + env%pack%energy
                  call clear_layer(env%pack, env%pack_density)
               else
                  call internal_energy_to_temp(env%pack%energy, env%pack%water, 0.0_wp, env%pack%temp, env%pack%fliq)
               end if
            end if
         end if
         if (took_pond) then
            !--- The drained water, at one temperature, refills the pond up to what it held.
            drained      = st%melt_mass + st%dump_mass
            drained_enth = st%melt_enth + st%dump_enth
            back      = min(drained, w_pond0)
            back_enth = 0.0_wp
            if (back > 0.0_wp) back_enth = drained_enth * (back / drained)
            w_pond = back
            e_pond = back_enth
            st%melt_mass = drained - back
            st%melt_enth = drained_enth - back_enth
            st%dump_mass = 0.0_wp
            st%dump_enth = 0.0_wp
            st%pond_mass = w_pond0 - back
            st%pond_enth = e_pond0 - back_enth
            !--- What a pack keeps of the pond that joined it adds its depth as ice.
            if (joined .and. ground_store_present(env%pack) .and. st%pond_mass > 0.0_wp) then
               env%pack_density      = combined_density(env%pack%water - st%pond_mass, env%pack_density,      &
                                                        st%pond_mass, rho_ice)
               env%pack%matrix_depth = env%pack%water / env%pack_density
            end if
         end if
         st%melt_rate = (st%melt_mass + st%dump_mass) / env%dt
         melt_enth    = st%melt_enth + st%dump_enth
         st%t_melt    = t_3ple
         if (st%melt_mass + st%dump_mass > 0.0_wp)                                                  &
            st%t_melt = liquid_energy_to_temp(melt_enth / (st%melt_mass + st%dump_mass))
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
      st%swe0     = swe0 + st%pond_mass
      st%enth0    = enth0 + st%pond_enth
   end subroutine prepare_snow_stage

   !--- Snow -- energy side: surface fluxes, base conduction, backward-Euler energy advance.

   !---~---
   !   Sensible and latent fluxes per unit snow area at a snow-surface temperature (positive
   ! up), the latent at the vapour enthalpy the canopy air receives, and saturated over ice
   ! by the pack's own liquid fraction (docs/science/snow_biophysics.md eq. (1)).
   !---~---
   pure subroutine snow_surface_fluxes(t_surf, fliq, env, h_flux, w_flux, le_flux)
      real(wp),         intent(in)  :: t_surf    ! [K] snow-surface temperature
      real(wp),         intent(in)  :: fliq      ! [-] liquid fraction of the snow surface
      type(surface_env_t), intent(in)  :: env
      real(wp),         intent(out) :: h_flux, w_flux, le_flux   ! [W/m2], [kg/m2/s], [W/m2]
      h_flux = sensible_heat_slope(env%ggnet, env%rho_air, env%cp_cas) * (t_surf - env%can_temp)
      w_flux = env%ggnet * env%rho_air * (sat_specific_humidity(t_surf, env%press, fliq) - env%can_shv)
      le_flux = w_flux * enthalpy_vapor(t_surf)
   end subroutine snow_surface_fluxes

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
   !   The pack's exchange with its surroundings over dt, area-weighted, from a linearized
   ! backward-Euler balance of its enthalpy that is not committed: the balance's one root T,
   ! in closed form on the melt plateau or else on the ice or liquid line, and the fluxes at
   ! T (docs/science/snow_biophysics.md, "The snow store"). The pack's base conducts through
   ! r_below to the node below it, at t_below. With no pack or no cover the fluxes are zero and
   ! t_surf is t_below. The pack takes the exchange's time integral in snow_take_exchange.
   !---~---
   pure subroutine snow_surface_balance(env, t_below, flux)
      type(surface_env_t), intent(in)  :: env       ! with the pack, its cover (which scales all boundary exchange), params, dt
      real(wp),            intent(in)  :: t_below   ! [K] the node the pack's base conducts into
      type(snow_out_t),    intent(out) :: flux      ! the exchange
      real(wp) :: t_n, fliq_n, e_old, swe0, gcond, af, dtaf
      real(wp) :: h_n, w_n, le_n, g_n, dh_dt, dw_dt, demit_dt, l_vap, r0, drdt
      real(wp) :: u_freeze, u_melt, e_bal, t_new, dtemp, h_s, w_s, le_s, g_s, emit, rad_n
      real(wp) :: dt, area_frac

      dt        = env%dt
      area_frac = env%snowfac
      associate (snow => env%pack, params => env%params)
         flux = snow_out_t()
         af = max(0.0_wp, min(1.0_wp, area_frac))
         if (.not. ground_store_present(snow) .or. snow%water <= params%tiny_snow_mass .or. af <= 0.0_wp) then
            flux%t_surf = t_below                                     ! no snow surface (or zero cover) to balance
            return
         end if

         swe0  = snow%water
         e_old = snow%energy
         dtaf  = dt * af
         flux%snowfac = af
         call internal_energy_to_temp(e_old, swe0, 0.0_wp, t_n, fliq_n)
         gcond = snow_base_conductance(snow%matrix_depth, env, params)

         !---~---
         !   The fluxes at T^n and their slopes. A kilogram of vapour carries enthalpy_vapor; less
         ! the pack's specific energy, that is the latent heat it takes from the pack
         ! (sublimation on ice, evaporation on water).
         !---~---
         call snow_surface_fluxes(t_n, fliq_n, env, h_n, w_n, le_n)
         g_n      = gcond * (t_n - t_below)
         dh_dt    = sensible_heat_slope(env%ggnet, env%rho_air, env%cp_cas)
         dw_dt    = env%ggnet * env%rho_air * sat_specific_humidity_temp_deriv(t_n, env%press, fliq_n)
         demit_dt = emission_slope(params%emissivity, t_n)
         !---~---
         !   The snow's net radiation at T^n: what it absorbs at its own albedo and emissivity,
         ! less its own emission.
         !---~---
         rad_n    = env%sw_snow + params%emissivity * (env%lw_down - stefan_boltzmann * t_n ** 4)
         l_vap    = enthalpy_vapor(t_n) - e_old / swe0
         r0       = rad_n - h_n - l_vap * w_n - g_n
         drdt     = -(demit_dt + dh_dt + l_vap * dw_dt + gcond)                ! < 0

         !---~---
         !   Solve swe e(T) = E^n + dt af (r0 + drdt (T - T^n)): on the plateau first, else on the
         ! ice or the liquid line.
         !---~---
         u_freeze = swe0 * cp_ice * t_3ple
         u_melt   = u_freeze + swe0 * latent_heat_fusion
         t_new    = t_3ple
         e_bal    = e_old + dtaf * (r0 + drdt * (t_3ple - t_n))
         if (e_bal < u_freeze) then                                             ! frozen: e = cp_ice T
            t_new = (e_old + dtaf * (r0 - drdt * t_n)) / (swe0 * cp_ice - dtaf * drdt)
            e_bal = swe0 * cp_ice * t_new
         else if (e_bal > u_melt) then                                          ! liquid: e = cp_liq (T - T_sc)
            t_new = (e_old + swe0 * cp_liq * tsupercool_liq + dtaf * (r0 - drdt * t_n))                   &
                    / (swe0 * cp_liq - dtaf * drdt)
            e_bal = swe0 * internal_energy_liquid(t_new)
         end if

         !--- The fluxes at T, on the linearization.
         dtemp = t_new - t_n
         h_s   = h_n + dh_dt * dtemp
         w_s   = w_n + dw_dt * dtemp
         g_s   = g_n + gcond * dtemp
         emit  = demit_dt * dtemp                     ! the emission added as the surface moves from T^n to T
         if (w_s * dtaf < swe0) then
            le_s = w_s * (l_vap + e_bal / swe0)       ! the latent heat plus the leaving mass's own energy
         else                                         ! the whole pack sublimates, with all the energy left
            w_s  = swe0 / dtaf
            le_s = e_old / dtaf + rad_n - emit - h_s - g_s
         end if

         !---~---
         !   The fluxes per unit ground area, positive away from the pack: h_snow + le_snow to the
         ! canopy air's enthalpy, w_flux to its vapour and g_base to the soil top. energy_in is
         ! the pack's net gain.
         !---~---
         flux%t_surf    = t_new
         flux%h_snow    = af * h_s
         flux%le_snow   = af * le_s
         flux%w_flux    = af * w_s
         flux%g_base    = af * g_s
         flux%rnet      = af * (rad_n - emit)                     ! snow's area-weighted net radiation into the store
         flux%energy_in = flux%rnet - flux%h_snow - flux%le_snow - flux%g_base
      end associate
   end subroutine snow_surface_balance

   !---~---
   !   The organic layer's exchange over dt, from a linearized backward-Euler balance of its
   ! energy, as the pack's (snow_surface_balance), with its dry heat capacity on each line
   ! (docs/science/surface_organic_layer.md): on the snow-free ground its net radiation, sensible
   ! heat and evaporation, its water's vapour crossing its dried top in series with the ground's
   ! conductance; under the pack the pack's base flux g_snow, held; and over the whole patch
   ! conduction into the soil top, its half-depth and the top soil node in series. Implicit in its
   ! own temperature, so a thin layer of little heat capacity relaxes without overshooting.
   !---~---
   pure subroutine organic_surface_balance(env, g_snow, flux)
      type(surface_env_t), intent(in)  :: env
      real(wp),            intent(in)  :: g_snow     ! [W/m2] the pack's base flux into the layer
      type(organic_out_t), intent(out) :: flux
      real(wp) :: af, dt, w, cd, e_old, t_n, fliq_n, g_os, h_n, dh_dt, w_n, dw_dt, demit_dt
      real(wp) :: rad_n, u_w, l_vap, g_n, r0, drdt, u_freeze, u_melt, e_bal, t_new, dtemp, h_s, w_s, g_s, emit
      real(wp) :: u_w_new, le_s, fl_new, w_cap, g_vap
      integer  :: pass
      flux = organic_out_t()
      af    = 1.0_wp - max(0.0_wp, min(1.0_wp, env%snowfac))     ! the snow-free ground
      dt    = env%dt
      w     = env%org_water
      cd    = env%org_dry_cap
      e_old = env%org_energy
      call internal_energy_to_temp(e_old, w, cd, t_n, fliq_n)
      g_os = conductance_in_series(0.5_wp * env%org_resistance, env%r_soil_top)   ! its lower half and the top soil node

      !--- The fluxes at T^n and their slopes; a kilogram of vapour takes its enthalpy less the water's energy.
      dh_dt = sensible_heat_slope(env%ggnet, env%rho_air, env%cp_cas)
      h_n   = dh_dt * (t_n - env%can_temp)
      !---~---
      !   Evaporation and dew alike, from the water below its dried top: a full layer exchanges as
      ! open water, a drying one through the dry top's resistance, an empty one not at all (w_cap).
      ! Held at the water network's rate once it has taken it from the slots.
      !---~---
      if (env%org_held) then
         w_n = 0.0_wp
         if (af > 0.0_wp) w_n = env%org_evap_held / af
         dw_dt = 0.0_wp
      else
         g_vap = 1.0_wp / (1.0_wp / max(env%ggnet, tiny_num) + env%org_dry_resist)
         w_n   = g_vap * env%rho_air * (sat_specific_humidity(t_n, env%press, fliq_n) - env%can_shv)
         dw_dt = g_vap * env%rho_air * sat_specific_humidity_temp_deriv(t_n, env%press, fliq_n)
      end if
      demit_dt = emission_slope(env%ground_emiss, t_n)
      rad_n = env%rnet_soil                                      ! at T^n: its own emission included
      u_w   = water_specific_energy(t_n, fliq_n)                ! from the datum: a trace of water stays exact
      l_vap = enthalpy_vapor(t_n) - u_w
      g_n   = g_os * (t_n - env%t_soil_top)
      w_cap = huge(1.0_wp)
      if (af > 0.0_wp .and. .not. env%org_held) w_cap = w / (dt * af)   ! the evaporation that empties it over dt

      !---~---
      !   Solve e(T) = E^n + dt (r0 + drdt (T - T^n)): on the plateau first, else on the ice or the
      ! liquid line. An evaporation the linearization takes past the water the layer holds is held
      ! there and the balance solved again, so the energy the layer takes is the root's: with a few
      ! joules per kelvin, any mismatch would be a large temperature error.
      !---~---
      do pass = 1, 2
         r0    = af * (rad_n - h_n - l_vap * w_n) + g_snow - g_n
         drdt  = -(af * (demit_dt + dh_dt + l_vap * dw_dt) + g_os)  ! < 0
         u_freeze = (cd + w * cp_ice) * t_3ple
         u_melt   = u_freeze + w * latent_heat_fusion
         t_new    = t_3ple
         e_bal    = e_old + dt * (r0 + drdt * (t_3ple - t_n))
         fl_new   = min(1.0_wp, max(0.0_wp, (e_bal - u_freeze) / max(w * latent_heat_fusion, tiny_num)))
         if (e_bal < u_freeze) then
            t_new  = (e_old + dt * (r0 - drdt * t_n)) / (cd + w * cp_ice - dt * drdt)
            e_bal  = (cd + w * cp_ice) * t_new
            fl_new = 0.0_wp
         else if (e_bal > u_melt) then
            t_new  = (e_old + w * cp_liq * tsupercool_liq + dt * (r0 - drdt * t_n)) / (cd + w * cp_liq - dt * drdt)
            e_bal  = cd * t_new + w * internal_energy_liquid(t_new)
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
      h_s   = h_n + dh_dt * dtemp
      g_s   = g_n + g_os * dtemp
      emit  = demit_dt * dtemp
      u_w_new = water_specific_energy(t_new, fl_new)
      le_s = w_s * (l_vap + u_w_new)                           ! the latent heat plus the leaving water's own energy

      flux%t_surf    = t_new
      flux%rad       = af * (rad_n - emit)
      flux%h         = af * h_s
      flux%le        = af * le_s
      flux%w_flux    = af * w_s
      flux%g_soil    = g_s
      flux%energy_in = flux%rad - flux%h - flux%le + g_snow - g_s
   end subroutine organic_surface_balance

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
      if (.not. (snow%water - vapour_amount >= 0.0_wp                                                      &
                 .or. within_roundoff(snow%water - vapour_amount, snow%water + abs(vapour_amount))))         &
         error stop 'snow_take_exchange: the pack gave the canopy air more vapour than it held'
      snow%energy       = snow%energy + energy_amount
      snow%water        = max(0.0_wp, snow%water - vapour_amount)   ! clamp-ok: round-off (asserted)
      snow%matrix_depth = snow%water / density
      if (snow%water > params%tiny_snow_mass) then
         call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)
      else
         snow%temp = t_3ple                               ! vanished; snow_drain_meltwater dumps the rest
         snow%fliq = 0.0_wp
      end if
   end subroutine snow_take_exchange

end module ed_surface
