! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_surface -- the ground's surface: the bare-ground skin and the one-layer snow          !
! store, blended by snowfac, each giving the soil-top boundary and the canopy-air fluxes   !
! (docs/science/snow_biophysics.md):                                                       !
!   * ground_surface_exchange -- the blend at one canopy-air and soil-top state;           !
!   * ground_surface_fluxes -- the bare ground's sensible and latent fluxes;               !
!   * prepare_snow_stage -- the pack before a step: its precipitation, drainage and cover; !
!   * snow_cover_fraction / snow_accumulate / snow_drain_meltwater -- the snow's mass;     !
!   * snow_blend_ground_optics / ground_tile_radiation -- the ground's optics under a      !
!     pack, and the radiation each tile absorbs at its own;                                !
!   * snow_surface_fluxes / snow_base_conductance / snow_surface_balance /                 !
!     snow_take_exchange -- its energy.                                                    !
!==========================================================================================!
module ed_surface
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, stefan_boltzmann, &
                                       latent_heat_fusion, tsupercool_liq
   use ed_therm_lib,          only :   internal_energy_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_liquid, precip_enthalpy, liquid_energy_to_temp
   use ed_patch_state_vars,     only : snow_state_t
   use ed_biophysics_params,    only : snow_params_t
   implicit none
   private

   !----- [kg/m3] The snow density Niu & Yang (2007)'s cover-fraction scale is referenced to. ----!
   real(wp), parameter :: NY07_RHO_REF = 100.0_wp

   public :: ground_surface_exchange, ground_surface_fluxes, prepare_snow_stage
   public :: snow_cover_fraction, snow_accumulate, snow_drain_meltwater, snow_blend_ground_optics
   public :: ground_tile_radiation
   public :: snow_take_exchange
   public :: surface_env_t, snow_out_t

   !----- The ground surface's surroundings, read by the snow balance and the bare ground alike:   !
   !      the radiation each tile absorbs and its conductance to the canopy air, the canopy air's  !
   !      and the soil top's state, the bare ground's evaporation, and the pack it exchanges with, !
   !      the pack's cover and parameters, and the step its balance relaxes over                   !
   !      (prepare_snow_stage sets the pack and the cover). ---------------------------------------!
   type :: surface_env_t
      real(wp) :: rnet_soil = 0.0_wp                    !< [W/m2] the bare soil's net all-wave radiation per unit of its area
      real(wp) :: sw_snow   = 0.0_wp                    !< [W/m2] the shortwave the snow absorbs per unit of its area
      real(wp) :: lw_down   = 0.0_wp                    !< [W/m2] the longwave reaching the ground
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp !< [K],[kg/kg] the canopy air's temperature and specific humidity
      real(wp) :: ggnet    = 0.0_wp                     !< [m/s] ground to canopy-air conductance, heat and vapour alike
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp !< [kg/m3],[Pa] the canopy air's density and pressure
      real(wp) :: t_soil_top = 273.16_wp                !< [K] top soil node's temperature, the snow base's conduction target
      real(wp) :: k_soil_top = 0.0_wp                   !< [W/m/K] top soil layer's thermal conductivity
      real(wp) :: dz_soil_top = 0.05_wp                 !< [m] top soil node's depth |z_node(1)| below the snow base
      real(wp) :: soil_evap  = 0.0_wp                   !< [kg/m2/s] the bare ground's evaporation, already area-weighted
      type(snow_state_t)  :: pack                       !< the pack after the step's precipitation and drainage
      real(wp) :: snowfac    = 0.0_wp                   !< [-] its Niu-Yang cover fraction (0 without a pack)
      type(snow_params_t) :: params                     !< the snow parameters
      real(wp) :: dt         = 0.0_wp                   !< [s] the step
   end type surface_env_t

   !----- What the snow kernels return, per unit ground area. snow_surface_balance fills the       !
   !      exchange (positive upward, away from the pack) and snow_drain_meltwater the meltwater    !
   !      handed to the pond; each resets the whole record, so the other's fields read zero.       !
   !      prepare_snow_stage returns the drainage with the stage's hand-off: the meltwater as the  !
   !      pond's inflow, and the precipitation's enthalpy and the pack before the stage, for the   !
   !      ledgers. All zero without a pack. -------------------------------------------------------!
   type :: snow_out_t
      !----- The surface balance. -----------------------------------------------------------------!
      real(wp) :: t_surf  = 273.16_wp                   !< [K] diagnosed snow-surface temperature
      real(wp) :: h_snow  = 0.0_wp, le_snow = 0.0_wp    !< [W/m2] sensible + latent (sublimation/evap) to the canopy air
      real(wp) :: w_flux  = 0.0_wp                      !< [kg/m2/s] vapour to the canopy air (>0 sublimation, <0 deposition)
      real(wp) :: g_base  = 0.0_wp                      !< [W/m2] conduction from the snow base into the soil top
      real(wp) :: rnet    = 0.0_wp                      !< [W/m2] net all-wave radiation absorbed by the snow surface
      real(wp) :: energy_in = 0.0_wp                    !< [W/m2] what the pack gains: rnet - h_snow - le_snow - g_base
      real(wp) :: snowfac = 0.0_wp                      !< [-] Niu-Yang07 snow-cover / burial fraction
      !----- The meltwater drainage, a store-to-store hand-off to the pond. -----------------------!
      real(wp) :: melt_mass = 0.0_wp                    !< [kg/m2] free liquid drained this step -> the pond
      real(wp) :: melt_enth = 0.0_wp                    !< [J/m2]  matching enthalpy (internal_energy_liquid(t_3ple))
      real(wp) :: dump_mass = 0.0_wp                    !< [kg/m2] full-melt-out residual water dumped to the pond
      real(wp) :: dump_enth = 0.0_wp                    !< [J/m2]  full-melt-out residual energy dumped to the pond
      !----- The stage's hand-off. ----------------------------------------------------------------!
      logical  :: exists    = .false.                   !< a pack is present after the step's precipitation
      real(wp) :: melt_rate = 0.0_wp                    !< [kg/m2/s] meltwater and dump to the pond over the step
      !----- The drained water's temperature, at which its liquid internal energy equals its       !
      !      enthalpy, so the pond receives exactly that enthalpy; t_3ple when none drains. -------!
      real(wp) :: t_melt    = t_3ple                    !< [K]
      real(wp) :: acc_enth  = 0.0_wp                    !< [J/m2]  precipitation (rain and snow) enthalpy that entered the pack
      real(wp) :: swe0      = 0.0_wp                    !< [kg/m2] pack mass before the stage, for the ledger
      real(wp) :: enth0     = 0.0_wp                    !< [J/m2]  pack internal energy before the stage, for the ledger
   end type snow_out_t

contains

   !=======================================================================================!
   !  Ground / surface skin                                                                !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! ground_surface_exchange -- the ground surface at one canopy-air and soil-top state:   !
   ! the snowfac-weighted pack, its balance relaxed over the step with its own heat store  !
   ! (docs/science/snow_biophysics.md, "The snow store"), and the bare soil on the rest.   !
   ! The fluxes come snowfac-weighted and are zero with no pack or no cover, so with no    !
   ! snow this is exactly bare ground, with no branch. soil_evap already carries its       !
   ! snow-free area factor (ground_evap_env_t%snow_free_frac), so le_soil is added whole.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, snow_energy_in,           &
                                           snow_vapour_out)
      type(surface_env_t),    intent(in)  :: env               !<           the ground surface's surroundings
      real(wp),               intent(out) :: h_ground          !< [W/m2]    sensible heat to the canopy air
      real(wp),               intent(out) :: le_ground         !< [W/m2]    latent heat to the canopy air
      real(wp),               intent(out) :: g_top             !< [W/m2]    net energy into the soil top
      real(wp),               intent(out) :: ground_rad        !< [W/m2]    net radiation, snow's and bare soil's
      real(wp),               intent(out) :: snow_energy_in    !< [W/m2]    energy into the pack
      real(wp),               intent(out) :: snow_vapour_out   !< [kg/m2/s] vapour from the pack to the canopy air
      type(snow_out_t) :: sfx
      real(wp) :: h_bare, le_soil
      call snow_surface_balance(env, sfx)
      call ground_surface_fluxes(env%t_soil_top, env%can_temp, env%ggnet, env%rho_air, env%soil_evap, h_bare, le_soil)
      h_ground  = sfx%h_snow  + (1.0_wp - env%snowfac) * h_bare
      le_ground = sfx%le_snow + le_soil
      g_top     = sfx%g_base                                                                         &
                  + (1.0_wp - env%snowfac) * (env%rnet_soil - h_bare)                                 &
                  - le_soil
      ground_rad      = sfx%rnet + (1.0_wp - env%snowfac) * env%rnet_soil
      snow_energy_in  = sfx%energy_in
      snow_vapour_out = sfx%w_flux
   end subroutine ground_surface_exchange

   !---------------------------------------------------------------------------------------!
   ! Bare-ground fluxes to the canopy air, instantaneous: sensible through ggnet, latent   !
   ! from the soil evaporation (already area-weighted) at the vapour enthalpy the canopy   !
   ! air receives. ground_surface_exchange assembles G_top and the snow blend.             !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine ground_surface_fluxes(t_ground, t_cas, ggnet, rho, soil_evap, h_bare, le_soil)
      real(wp), intent(in)  :: t_ground, t_cas, ggnet, rho, soil_evap   !< [K], [K], [m/s], [kg/m3], [kg/m2/s]
      real(wp), intent(out) :: h_bare, le_soil                          !< [W/m2], [W/m2]
      h_bare  = ggnet * rho * cp_air * (t_ground - t_cas)
      le_soil = soil_evap * enthalpy_vapor(t_ground)
   end subroutine ground_surface_fluxes

   !=======================================================================================!
   !  Snow -- mass side: cover fraction, accumulation, meltwater drainage.                 !
   !=======================================================================================!

   !----- Niu-Yang (2007) snow-cover fraction [-], tanh(depth/scale): it ramps the ground optics   !
   !      and the soil-top blend and area-weights the pack's exchange                              !
   !      (docs/science/snow_biophysics.md). The scale is set by the ground's roughness, the       !
   !      micro-relief a thin pack has to fill. ---------------------------------------------------!
   pure function snow_cover_fraction(swe, snow_depth, params, ground_rough) result(snowfac)
      real(wp),            intent(in) :: swe, snow_depth   !< [kg/m2], [m]
      type(snow_params_t), intent(in) :: params
      real(wp),            intent(in) :: ground_rough   !< [m] the ground's roughness ([aerodynamics].soil_rough)
      real(wp) :: snowfac, scale
      if (swe <= params%tiny_snow_mass) then
         snowfac = 0.0_wp
      else
         scale   = max(params%ny07_a * ground_rough * (params%rho_snow / NY07_RHO_REF) ** params%ny07_m, tiny_num)
         snowfac = tanh(snow_depth / scale)
      end if
   end function snow_cover_fraction

   !----- The ground's optics under a pack, for the two-stream: its albedo and emissivity blended  !
   !      toward the pack's by the snow-cover fraction, so a partial pack is continuous, and the   !
   !      temperature at which the blended emissivity emits what the soil and the pack emit at     !
   !      their own (docs/science/snow_biophysics.md, "The ground's radiation"). Unchanged without !
   !      a pack. ---------------------------------------------------------------------------------!
   pure subroutine snow_blend_ground_optics(snow, params, ground_rough, albedo_vis, albedo_nir, emiss, temp)
      type(snow_state_t),  intent(in)    :: snow
      type(snow_params_t), intent(in)    :: params
      real(wp),            intent(in)    :: ground_rough   !< [m] the ground's roughness ([aerodynamics].soil_rough)
      real(wp),            intent(inout) :: albedo_vis   !< [-] VIS albedo: the bare soil's on entry
      real(wp),            intent(inout) :: albedo_nir   !< [-] NIR albedo: the bare soil's on entry
      real(wp),            intent(inout) :: emiss        !< [-] thermal emissivity: the bare soil's on entry
      real(wp),            intent(inout) :: temp         !< [K] emission temperature: the top soil node's on entry
      real(wp) :: fc, snow_vis, snow_nir, emit_soil, emit_snow
      if (.not. (snow%nlayer >= 1_ik .and. snow%swe(1) > params%tiny_snow_mass)) return
      fc = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), params, ground_rough)
      call snow_albedos(snow, params, snow_vis, snow_nir)
      emit_soil  = emiss * temp ** 4                                ! [K4] emission over sigma
      emit_snow  = params%snow_emiss * snow%snow_temp(1) ** 4
      albedo_vis = (1.0_wp - fc) * albedo_vis + fc * snow_vis
      albedo_nir = (1.0_wp - fc) * albedo_nir + fc * snow_nir
      emiss      = (1.0_wp - fc) * emiss + fc * params%snow_emiss
      temp       = (((1.0_wp - fc) * emit_soil + fc * emit_snow) / emiss) ** 0.25_wp
   end subroutine snow_blend_ground_optics

   !----- The pack's own shortwave albedos, fresh to aged by its surface liquid fraction; a pack   !
   !      this step's snowfall has yet to make is fresh. ------------------------------------------!
   pure subroutine snow_albedos(snow, params, albedo_vis, albedo_nir)
      type(snow_state_t),  intent(in)  :: snow
      type(snow_params_t), intent(in)  :: params
      real(wp),            intent(out) :: albedo_vis, albedo_nir   !< [-]
      real(wp) :: fl
      fl = 0.0_wp
      if (snow%nlayer >= 1_ik) fl = snow%snow_fliq(1)
      albedo_vis = (1.0_wp - fl) * params%albedo_vis_fresh + fl * params%albedo_vis_aged
      albedo_nir = (1.0_wp - fl) * params%albedo_nir_fresh + fl * params%albedo_nir_aged
   end subroutine snow_albedos

   !----- The radiation each tile of the ground absorbs per unit of its own area, from what        !
   !      reaches the ground in each band (docs/science/snow_biophysics.md, "The ground's          !
   !      radiation"): the bare soil's net all-wave radiation at its albedo and emissivity and the !
   !      top soil node's temperature, and the shortwave the snow absorbs at its own albedo. The   !
   !      snow's longwave is formed in its balance, at its own temperature. Weighted by the cover, !
   !      the two give the ground absorption of the two-stream, whose optics                       !
   !      snow_blend_ground_optics blends. --------------------------------------------------------!
   pure subroutine ground_tile_radiation(dn_vis, dn_nir, dn_lw, soil_albedo_vis, soil_albedo_nir, soil_emiss,  &
                                         soil_temp, snow, params, rnet_soil, sw_snow)
      real(wp),            intent(in)  :: dn_vis, dn_nir, dn_lw   !< [W/m2] each band reaching the ground
      real(wp),            intent(in)  :: soil_albedo_vis, soil_albedo_nir, soil_emiss   !< [-] the bare soil's optics
      real(wp),            intent(in)  :: soil_temp               !< [K] the top soil node's temperature
      type(snow_state_t),  intent(in)  :: snow
      type(snow_params_t), intent(in)  :: params
      real(wp),            intent(out) :: rnet_soil, sw_snow      !< [W/m2]
      real(wp) :: snow_vis, snow_nir
      rnet_soil = dn_vis * (1.0_wp - soil_albedo_vis) + dn_nir * (1.0_wp - soil_albedo_nir)              &
                  + soil_emiss * (dn_lw - stefan_boltzmann * soil_temp ** 4)
      call snow_albedos(snow, params, snow_vis, snow_nir)
      sw_snow = dn_vis * (1.0_wp - snow_vis) + dn_nir * (1.0_wp - snow_nir)
   end subroutine ground_tile_radiation

   !----- Snowfall joins the pack as ice and rain-on-snow as liquid, with their precip_enthalpy. A !
   !      layer is created only from min_new_snow_mass of new snow; with no layer and less, the    !
   !      pack stays empty and the precipitation reaches the ground with the throughfall. ---------!
   pure subroutine snow_accumulate(snow, snowfall, rainf, air_temp, dt, params)
      type(snow_state_t), intent(inout) :: snow                      !< the pack
      real(wp),            intent(in)    :: snowfall, rainf, air_temp, dt   !< [kg/m2/s], [kg/m2/s], [K], [s]
      type(snow_params_t), intent(in)    :: params
      real(wp) :: dm_snow, dm_rain
      logical  :: has_layer
      dm_snow = max(0.0_wp, snowfall) * dt
      dm_rain = max(0.0_wp, rainf) * dt
      has_layer = (snow%nlayer >= 1_ik) .or. (dm_snow >= params%min_new_snow_mass)
      if (.not. has_layer) return                          ! bare ground and sub-threshold snow: no pack
      snow%swe(1)         = snow%swe(1) + dm_snow + dm_rain
      snow%snow_energy(1) = snow%snow_energy(1) + precip_enthalpy(dm_rain, dm_snow, air_temp)
      snow%snow_depth(1)  = snow%snow_depth(1) + dm_snow / params%rho_snow   ! rain fills pores / refreezes -> no bulk depth
      snow%nlayer         = 1_ik
   end subroutine snow_accumulate

   !----- Percolation: drain the free liquid above the holding capacity as a paired (mass,         !
   !      enthalpy) hand-off to the pond (melt_*), and dump the residual at melt-out (dump_*). ----!
   pure subroutine snow_drain_meltwater(snow, params, melt)
      type(snow_state_t), intent(inout) :: snow     !< the pack
      type(snow_params_t), intent(in)    :: params
      type(snow_out_t),   intent(out)   :: melt     !< the drained and dumped water, mass and enthalpy
      real(wp) :: wmass_free, hold
      melt = snow_out_t()
      if (snow%nlayer == 0_ik) return

      !----- Vanished or trace layer: dump the residual to the pond and revert to bare ground. ----!
      if (snow%swe(1) <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%swe(1))
         melt%dump_enth  = snow%snow_energy(1)
         call clear_layer(snow)
         return
      end if

      call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))

      !----- Drain the free liquid above the holding capacity (LEAF-3's 1:9). ---------------------!
      hold       = max(1.0_wp - params%liquid_holding_frac, tiny_num)
      wmass_free = max(0.0_wp, snow%swe(1) * (snow%snow_fliq(1) - params%liquid_holding_frac) / hold)
      wmass_free = min(wmass_free, snow%swe(1))
      if (wmass_free > 0.0_wp) then
         melt%melt_mass      = wmass_free
         melt%melt_enth      = wmass_free * internal_energy_liquid(t_3ple)   ! drains at the plateau
         snow%swe(1)         = snow%swe(1) - wmass_free
         snow%snow_energy(1) = snow%snow_energy(1) - melt%melt_enth
         snow%snow_depth(1)  = snow%swe(1) / params%rho_snow
      end if

      !----- Below the melt-out threshold after draining: dump the remainder to the pond. ---------!
      if (snow%swe(1) <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%swe(1))
         melt%dump_enth  = snow%snow_energy(1)
         call clear_layer(snow)
      else
         call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      end if
   end subroutine snow_drain_meltwater

   !----- Revert to the empty "no snow" state (nlayer = 0). ---------------------------------------!
   pure subroutine clear_layer(snow)
      type(snow_state_t), intent(inout) :: snow
      snow%swe(1)         = 0.0_wp
      snow%snow_energy(1) = 0.0_wp
      snow%snow_depth(1)  = 0.0_wp
      snow%snow_fliq(1)   = 0.0_wp
      snow%snow_temp(1)   = t_3ple
      snow%nlayer         = 0_ik
   end subroutine clear_layer

   !---------------------------------------------------------------------------------------!
   ! prepare_snow_stage -- the pack before a step: the step's snowfall and rain-on-snow,   !
   ! then the drainage of the liquid it holds (the last step's melt and this step's rain)  !
   ! to the pond as a paired (mass, enthalpy) transfer, and its cover fraction             !
   ! (docs/science/snow_biophysics.md, "The snow store"). env's pack enters as the pack at !
   ! the step's start and leaves advanced, with env's cover set; st returns the drainage   !
   ! and the stage's hand-off. The exchange is then with this pack                         !
   ! (ground_surface_exchange), and the pack takes it afterwards (snow_take_exchange).     !
   !---------------------------------------------------------------------------------------!
   pure subroutine prepare_snow_stage(snowfall, rainfall, t_air, ground_rough, env, st)
      real(wp),            intent(in)    :: snowfall       !< [kg/m2/s] frozen precipitation
      real(wp),            intent(in)    :: rainfall       !< [kg/m2/s] liquid precipitation (rain-on-snow)
      real(wp),            intent(in)    :: t_air          !< [K]       reference-level air temperature
      real(wp),            intent(in)    :: ground_rough   !< [m]       the ground's roughness ([aerodynamics].soil_rough)
      type(surface_env_t), intent(inout) :: env            !< its pack advanced, its cover set; params and dt read
      type(snow_out_t),    intent(out)   :: st             !<           the drainage and the stage's hand-off
      real(wp) :: swe0, enth0, acc_enth, melt_enth
      logical  :: exists

      swe0  = env%pack%swe(1)
      enth0 = env%pack%snow_energy(1)
      call snow_accumulate(env%pack, snowfall, rainfall, t_air, env%dt, env%params)
      acc_enth = env%pack%snow_energy(1) - enth0   ! precipitation (rain and snow) enthalpy into the pack
      exists   = env%pack%nlayer >= 1_ik           ! a pack took the snow and the rain
      env%snowfac = 0.0_wp

      if (exists) then
         !----- Meltwater goes to the pond with its enthalpy, valued at t_melt, and nowhere else; pack !
         !      and pond are both tracked stores, so the transfer cancels in the column's ledger. -----!
         call snow_drain_meltwater(env%pack, env%params, st)
         st%melt_rate = (st%melt_mass + st%dump_mass) / env%dt
         melt_enth    = st%melt_enth + st%dump_enth
         st%t_melt    = t_3ple
         if (st%melt_mass + st%dump_mass > 0.0_wp)                                                  &
            st%t_melt = liquid_energy_to_temp(melt_enth / (st%melt_mass + st%dump_mass))
         !----- snowfac is snow, 1 - snowfac bare soil; the pack's exchange scales by it, so a patchy  !
         !      pack barely exchanges and its fluxes come snowfac-weighted. ---------------------------!
         if (env%pack%nlayer >= 1_ik .and. env%pack%swe(1) > env%params%tiny_snow_mass)                 &
            env%snowfac = snow_cover_fraction(env%pack%swe(1), env%pack%snow_depth(1), env%params, ground_rough)
      end if
      st%exists   = exists
      st%acc_enth = acc_enth
      st%swe0     = swe0
      st%enth0    = enth0
   end subroutine prepare_snow_stage

   !=======================================================================================!
   !  Snow -- energy side: surface fluxes, base conduction, backward-Euler energy advance. !
   !=======================================================================================!

   !----- Sensible and latent fluxes per unit snow area at a snow-surface temperature (positive    !
   !      up), the latent at the vapour enthalpy the canopy air receives, and saturated over ice   !
   !      by the pack's own liquid fraction (docs/science/snow_biophysics.md eq. (1)). ------------!
   pure subroutine snow_surface_fluxes(t_surf, fliq, env, h_flux, w_flux, le_flux)
      real(wp),         intent(in)  :: t_surf    !< [K] snow-surface temperature
      real(wp),         intent(in)  :: fliq      !< [-] liquid fraction of the snow surface
      type(surface_env_t), intent(in)  :: env
      real(wp),         intent(out) :: h_flux, w_flux, le_flux   !< [W/m2], [kg/m2/s], [W/m2]
      h_flux = env%ggnet * env%rho_air * cp_air * (t_surf - env%can_temp)
      w_flux = env%ggnet * env%rho_air * (sat_specific_humidity(t_surf, env%press, fliq) - env%can_shv)
      le_flux = w_flux * enthalpy_vapor(t_surf)
   end subroutine snow_surface_fluxes

   !----- Snow-base to soil-top conductance [W/m2/K]: the half snow layer and the top soil node in !
   !      series, so a deepening pack decouples the soil (docs/science/snow_biophysics.md). -------!
   pure function snow_base_conductance(snow_depth, env, params) result(gcond)
      real(wp),            intent(in) :: snow_depth   !< [m]
      type(surface_env_t),    intent(in) :: env
      type(snow_params_t), intent(in) :: params
      real(wp) :: gcond, r_series
      r_series = 0.5_wp * max(snow_depth, tiny_num) / max(params%k_snow, tiny_num)                  &
               + max(env%dz_soil_top, tiny_num) / max(env%k_soil_top, tiny_num)
      gcond    = 1.0_wp / max(r_series, tiny_num)
   end function snow_base_conductance

   !----- The pack's exchange with its surroundings over dt, area-weighted, from a linearized      !
   !      backward-Euler balance of its enthalpy that is not committed: the balance's one root T,  !
   !      in closed form on the melt plateau or else on the ice or liquid line, and the fluxes at  !
   !      T (docs/science/snow_biophysics.md, "The snow store"). With no pack or no cover the      !
   !      fluxes are zero and t_surf is the soil top's temperature. The pack takes the exchange's  !
   !      time integral in snow_take_exchange. ----------------------------------------------------!
   pure subroutine snow_surface_balance(env, flux)
      type(surface_env_t), intent(in)  :: env    !< with the pack, its cover (which scales all boundary exchange), params, dt
      type(snow_out_t),    intent(out) :: flux   !< the exchange
      real(wp) :: t_n, fliq_n, e_old, swe0, gcond, af, dtaf
      real(wp) :: h_n, w_n, le_n, g_n, dh_dt, dw_dt, demit_dt, l_vap, r0, drdt
      real(wp) :: u_freeze, u_melt, e_bal, t_new, dtemp, h_s, w_s, le_s, g_s, emit, rad_n
      real(wp) :: dt, area_frac

      dt        = env%dt
      area_frac = env%snowfac
      associate (snow => env%pack, params => env%params)
         flux = snow_out_t()
         af = max(0.0_wp, min(1.0_wp, area_frac))
         if (snow%nlayer == 0_ik .or. snow%swe(1) <= params%tiny_snow_mass .or. af <= 0.0_wp) then
            flux%t_surf = env%t_soil_top                              ! no snow surface (or zero cover) to balance
            return
         end if

         swe0  = snow%swe(1)
         e_old = snow%snow_energy(1)
         dtaf  = dt * af
         flux%snowfac = af
         call internal_energy_to_temp(e_old, swe0, 0.0_wp, t_n, fliq_n)
         gcond = snow_base_conductance(snow%snow_depth(1), env, params)

         !----- The fluxes at T^n and their slopes. A kilogram of vapour carries enthalpy_vapor; less  !
         !      the pack's specific energy, that is the latent heat it takes from the pack             !
         !      (sublimation on ice, evaporation on water). -------------------------------------------!
         call snow_surface_fluxes(t_n, fliq_n, env, h_n, w_n, le_n)
         g_n      = gcond * (t_n - env%t_soil_top)
         dh_dt    = env%ggnet * env%rho_air * cp_air
         dw_dt    = env%ggnet * env%rho_air * sat_specific_humidity_temp_deriv(t_n, env%press, fliq_n)
         demit_dt = 4.0_wp * params%snow_emiss * stefan_boltzmann * t_n ** 3
         !----- The snow's net radiation at T^n: what it absorbs at its own albedo and emissivity, !
         !      less its own emission. ------------------------------------------------------------!
         rad_n    = env%sw_snow + params%snow_emiss * (env%lw_down - stefan_boltzmann * t_n ** 4)
         l_vap    = enthalpy_vapor(t_n) - e_old / swe0
         r0       = rad_n - h_n - l_vap * w_n - g_n
         drdt     = -(demit_dt + dh_dt + l_vap * dw_dt + gcond)                ! < 0

         !----- Solve swe e(T) = E^n + dt af (r0 + drdt (T - T^n)): on the plateau first, else on the  !
         !      ice or the liquid line. ---------------------------------------------------------------!
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

         !----- The fluxes at T, on the linearization. ------------------------------------------------!
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

         !----- The fluxes per unit ground area, positive away from the pack: h_snow + le_snow to the  !
         !      canopy air's enthalpy, w_flux to its vapour and g_base to the soil top. energy_in is   !
         !      the pack's net gain. ------------------------------------------------------------------!
         flux%t_surf    = t_new
         flux%h_snow    = af * h_s
         flux%le_snow   = af * le_s
         flux%w_flux    = af * w_s
         flux%g_base    = af * g_s
         flux%rnet      = af * (rad_n - emit)                     ! snow's area-weighted net radiation into the store
         flux%energy_in = flux%rnet - flux%h_snow - flux%le_snow - flux%g_base
      end associate
   end subroutine snow_surface_balance

   !----- The pack takes an exchange's amounts, energy in and vapour out, keeping e(T) per         !
   !      kilogram, so its temperature follows from the new energy. A pack whose mass is gone is   !
   !      left at t_3ple for snow_drain_meltwater to empty. ---------------------------------------!
   pure subroutine snow_take_exchange(snow, energy_amount, vapour_amount, params)
      type(snow_state_t),  intent(inout) :: snow                            !< the pack
      real(wp),            intent(in)    :: energy_amount, vapour_amount    !< [J/m2] in, [kg/m2] out
      type(snow_params_t), intent(in)    :: params
      if (snow%nlayer == 0_ik) return
      snow%snow_energy(1) = snow%snow_energy(1) + energy_amount
      snow%swe(1)         = max(0.0_wp, snow%swe(1) - vapour_amount)
      snow%snow_depth(1)  = snow%swe(1) / params%rho_snow
      if (snow%swe(1) > params%tiny_snow_mass) then
         call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      else
         snow%snow_temp(1) = t_3ple                               ! vanished; snow_drain_meltwater dumps the rest
         snow%snow_fliq(1) = 0.0_wp
      end if
   end subroutine snow_take_exchange

end module ed_surface
