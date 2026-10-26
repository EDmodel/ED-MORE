! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_fast_forcing -- turns the forcing source into each fast step's drivers. The window's   !
! polygon samples (sample_window), each site's sample after its terrain (site_forcing), and !
! per patch and sub-step the air at the patch's canopy-air top and the two-stream canopy   !
! radiation (patch_forcing; docs/science/forcing.md, canopy_radiation_transfer.md §4).     !
!==========================================================================================!
module ed_fast_forcing
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : tiny_num, par_w_2_umol
   use ed_config,               only : ed_config_t
   use ed_therm_lib,            only : cas_temp_to_enthalpy, precip_phase
   use ed_time,                 only : ed_time_t, time_advance_seconds, seconds_into_day, slope_incidence_cosz
   use ed_met_types,            only : met_source_t, met_cursor_t, met_forcing_t
   use ed_met,                  only : met_advance, met_instant, met_to_cas_top, met_at_site
   use ed_radiation_lib,        only : slope_beam_factor
   use ed_polygon_state_vars,   only : polygon_t
   use ed_pft_params,           only : N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use ed_canopy_radiation, only :   rad_env_t, rad_out_t, alloc_rad_env, ground_optics_env_t
   use ed_canopy_radiation,     only : canopy_radiation, ground_optics
   use ed_surface,              only : snow_blend_ground_optics, ground_tile_radiation
   use ed_fast_integ_types,         only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t
   implicit none
   private

   public :: sample_window, site_forcing, patch_forcing

contains

   !----- This slow step's polygon samples, one per sub-step: met_advance in increasing time, and     !
   !      met_instant, a pure function of (reader state, time). Made once per window outside the    !
   !      patch loop: they are the same for every patch, and the file reader stays out of the       !
   !      parallel region. ----------------------------------------------------------------------!
   subroutine sample_window(met_src, met_cur, poly, cfg, step_start, met_sample, t_sample)
      type(met_source_t),  intent(in)    :: met_src
      type(met_cursor_t),  intent(inout) :: met_cur
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      type(ed_time_t),     intent(in)    :: step_start      !< calendar time at the start of the slow step
      type(met_forcing_t), intent(out)   :: met_sample(:)   !< (sub-step)
      type(ed_time_t),     intent(out)   :: t_sample(:)     !< (sub-step) each sample's time
      integer(ik) :: isub
      met_cur%kt_last_day = poly%kt_last_day          ! the polygon's memory of the sky, carried over
      do isub = 1_ik, size(met_sample, kind=ik)
         t_sample(isub) = time_advance_seconds(step_start,                                          &
                    (real(isub, wp) - 1.0_wp + cfg%forcing_sample_frac) * cfg%dt_fast)
         call met_advance(met_src, met_cur, t_sample(isub))
         met_sample(isub) = met_instant(met_src, met_cur, t_sample(isub))
      end do
      poly%kt_last_day = met_cur%kt_last_day
   end subroutine sample_window

   !----- Each site's sample, moved to its elevation and its beam to its slope (met_at_site,          !
   !      docs/science/forcing.md §8); a flat site at the polygon's elevation takes it as is. ----------!
   subroutine site_forcing(poly, cfg, met_cur, met_sample, t_sample, met_site)
      type(polygon_t),     intent(in)  :: poly
      type(ed_config_t),   intent(in)  :: cfg
      type(met_cursor_t),  intent(in)  :: met_cur
      type(met_forcing_t), intent(in)  :: met_sample(:)   !< (sub-step) the polygon's samples
      type(ed_time_t),     intent(in)  :: t_sample(:)     !< (sub-step) their times
      type(met_forcing_t), intent(out) :: met_site(:,:)   !< (sub-step, site)
      integer(ik) :: is, isub
      do is = 1_ik, poly%site%n
         associate (surf => poly%site%surface(is))
            do isub = 1_ik, size(met_sample, kind=ik)
               met_site(isub, is) = met_sample(isub)
               if (surf%elevation_m == poly%elevation_m .and. surf%slope_deg == 0.0_wp) cycle
               associate (m => met_sample(isub))
                  met_site(isub, is) = met_at_site(m, cfg%forcing, t_sample(isub)%month,                &
                     surf%elevation_m - poly%elevation_m,                                             &
                     slope_beam_factor(m%par_beam + m%nir_beam, m%cosz,                                  &
                                       slope_incidence_cosz(t_sample(isub), seconds_into_day(t_sample(isub)), &
                                                            met_cur%latitude_deg, met_cur%longitude_deg,      &
                                                            surf%slope_deg, surf%aspect_deg)))
               end associate
            end do
         end associate
      end do
   end subroutine site_forcing

   !---------------------------------------------------------------------------------------!
   ! The patch's drivers for one sub-step from its site's sample `met`: the air, moved to the  !
   ! patch's canopy-air top along its log profile when to_top (met_to_cas_top,                !
   ! docs/science/forcing.md §8; without a forcing source the record already is the air       !
   ! there), and the two-stream radiation, which the height does not change.                  !
   !---------------------------------------------------------------------------------------!
   subroutine patch_forcing(cfg, static, dynamic, met, to_top, forcing)
      type(ed_config_t),          intent(in)    :: cfg
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(met_forcing_t),        intent(in)    :: met
      logical,                    intent(in)    :: to_top
      type(fast_forcing_t),       intent(inout) :: forcing   !< sized to the largest patch (ensure_forcing_capacity)
      type(met_forcing_t) :: met_top
      if (to_top) then
         met_top = met_to_cas_top(met, cfg%forcing, static%z_top, static%displace, static%rough)
         call fill_forcing(forcing, met_top)
      else
         call fill_forcing(forcing, met)
      end if
      call apply_rt_forcing(forcing, cfg, static, dynamic, met)
   end subroutine patch_forcing

   !----- The patch's air from the sub-step's record: its enthalpy, humidity and CO2, the wind,       !
   !      pressure and density the aerodynamics reads, and the precipitation; apply_rt_forcing then   !
   !      fills the radiation. ---------------------------------------------------------------------!
   subroutine fill_forcing(forcing, met)
      type(fast_forcing_t), intent(inout) :: forcing
      type(met_forcing_t),  intent(in)    :: met
      forcing%enthalpy_atm = cas_temp_to_enthalpy(met%air_temp, met%air_shv)
      forcing%shv_atm      = met%air_shv
      forcing%co2_atm      = met%co2
      forcing%wind         = met%wind
      forcing%press        = met%air_pressure
      forcing%rho_air      = met%air_density
      !----- The rain/snow split at the record's temperature: with forcing, the patch's canopy-air top, !
      !      where the precipitation's enthalpy is valued too (docs/science/forcing.md §7). ------------!
      call precip_phase(met%precip, met%air_temp, forcing%rainfall, forcing%snowfall)
      forcing%air_temp     = met%air_temp                ! rainfall enthalpy reference (snow/rain-on-snow)
   end subroutine fill_forcing

   !----- Run the patch's two-stream canopy radiation and fill each cohort's absorbed SW,            !
   !      incident-equivalent PAR and net LW, the ground's SW and LW, and the top-of-canopy fluxes  !
   !      (docs/science/canopy_radiation_transfer.md §4). Cohorts arrive tallest first and the      !
   !      solver takes them bottom to top, in static%order, so outputs go back through it. ---------!
   subroutine apply_rt_forcing(forcing, cfg, static, dynamic, met)
      type(fast_forcing_t),       intent(inout) :: forcing
      type(ed_config_t),          intent(in)    :: cfg
      type(fast_static_state_t),  intent(in)    :: static     !< the cohorts (gather order, top first)
      type(fast_dynamic_state_t), intent(in)    :: dynamic    !< lagged tissue temperatures, top soil node, snow
      type(met_forcing_t),        intent(in)    :: met
      integer(ik) :: j, i, ig, ncoh
      integer(ik) :: pft_bt(static%n)
      real(wp)    :: lai_bt(static%n), wai_bt(static%n), hgt_bt(static%n)
      real(wp)    :: tleaf_bt(static%n), twood_bt(static%n)
      real(wp)    :: soil_albedo(N_RAD_BAND_DEFAULT), soil_temp_top
      type(rad_env_t)   :: rf
      type(rad_out_t)      :: flux
      type(ground_optics_env_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)

      ncoh = static%n
      !----- Bare-soil albedo per band from [soil] (no longwave reflection) and the ground's emission  !
      !      temperature, the top soil node's. -------------------------------------------------------!
      soil_albedo   = [cfg%soil%ground_albedo_vis, cfg%soil%ground_albedo_nir, 0.0_wp]
      soil_temp_top = dynamic%soil_energy%soil_temp(1)
      !----- A bare patch needs no special case: canopy_radiation's empty canopy gives the ground's. -!

      !----- The cohorts bottom -> top. ----------------------------------------------------------!
      do j = 1_ik, ncoh
         i = static%order(j)
         pft_bt(j) = static%pft(i)
         lai_bt(j) = static%lai(i)
         hgt_bt(j) = static%height(i)
         wai_bt(j) = static%wai(i)
         !----- The tissues emit at their start-of-sub-step temperatures. ---------------------------!
         tleaf_bt(j) = dynamic%leaf_temp(i)
         twood_bt(j) = dynamic%wood_temp(i)
      end do

      !----- rad_env_t from met (all W/m2, direct assignment). --------------------------------!
      call alloc_rad_env(rf, N_RAD_BAND_DEFAULT)
      rf%cosz = met%cosz
      rf%incid_beam(RAD_VIS) = met%par_beam
      rf%incid_diff(RAD_VIS) = met%par_diffuse
      rf%incid_beam(RAD_NIR) = met%nir_beam
      rf%incid_diff(RAD_NIR) = met%nir_diffuse
      rf%incid_beam(RAD_LW)  = 0.0_wp
      rf%incid_diff(RAD_LW)  = met%lw_down

      !----- Ground optics: bare-soil albedo/emissivity from [soil], skin temperature from the top soil node. --!
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = soil_albedo
      surf%soil_emiss  = cfg%soil%ground_emissivity
      surf%soil_temp   = soil_temp_top
      !----- A pack blends them toward its own (ed_surface). -----------------------------------  !
      call snow_blend_ground_optics(dynamic%snow, cfg%snow, cfg%aero%soil_rough, surf%soil_albedo(RAD_VIS),   &
                                    surf%soil_albedo(RAD_NIR), surf%soil_emiss, surf%soil_temp)
      he = [.false., .false., .true.]
      call ground_optics(surf, N_RAD_BAND_DEFAULT, he, rf%grnd_refl, rf%grnd_emiss)

      call canopy_radiation(cfg%pft%optics, rf, ncoh, pft_bt, hgt_bt, lai_bt, wai_bt, tleaf_bt, twood_bt, flux)

      !----- inverse-scatter: RT index j (bottom->top) maps to gather index static%order(j). ------!
      do j = 1_ik, ncoh
         ig = static%order(j)
         forcing%abs_sw(ig)      = flux%abs_leaf(RAD_VIS, j) + flux%abs_leaf(RAD_NIR, j)   ! absorbed leaf SW
         forcing%abs_par(ig)     = flux%abs_leaf(RAD_VIS, j)                                          &
                                   / max(cfg%leaf_physiology%leaf_absorptance, tiny_num)   ! incident PAR
         forcing%abs_lw(ig)      = flux%abs_leaf(RAD_LW, j)                  ! net leaf LW, its own emission included
         forcing%abs_sw_wood(ig) = flux%abs_wood(RAD_VIS, j) + flux%abs_wood(RAD_NIR, j)   ! absorbed wood SW
         forcing%abs_lw_wood(ig) = flux%abs_wood(RAD_LW, j)                  ! net wood LW, its own emission included
         forcing%par_top(ig)     = flux%incid_top(RAD_VIS, j)                ! PAR at the cohort's top
      end do
      forcing%abs_sw_ground = (flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS))                      &
                            + (flux%dn_ground(RAD_NIR) - flux%up_ground(RAD_NIR))
      forcing%abs_lw_ground = flux%dn_ground(RAD_LW) - flux%up_ground(RAD_LW)   ! net ground LW, emission included
      !----- What each tile takes of the radiation reaching the ground, at its own optics. ----------!
      call ground_tile_radiation(flux%dn_ground(RAD_VIS), flux%dn_ground(RAD_NIR), flux%dn_ground(RAD_LW),     &
                                 soil_albedo(RAD_VIS), soil_albedo(RAD_NIR), cfg%soil%ground_emissivity,       &
                                 soil_temp_top, dynamic%snow, cfg%snow, forcing%rnet_soil, forcing%sw_snow)
      forcing%lw_down_ground = flux%dn_ground(RAD_LW)
      forcing%par_per_w     = par_w_2_umol                     ! true VIS absorbed -> photon flux
      !----- Top-of-canopy fluxes: each band's albedo times the incident it was formed from, kept as   !
      !      fluxes (docs/science/canopy_radiation_transfer.md, Surface radiative fluxes on the       !
      !      diagnostic path). ------------------------------------------------------------------!
      forcing%sw_in_vis = rf%incid_beam(RAD_VIS) + rf%incid_diff(RAD_VIS)
      forcing%sw_in_nir = rf%incid_beam(RAD_NIR) + rf%incid_diff(RAD_NIR)
      forcing%sw_up_vis = flux%albedo(RAD_VIS) * forcing%sw_in_vis
      forcing%sw_up_nir = flux%albedo(RAD_NIR) * forcing%sw_in_nir
      forcing%lw_up     = flux%albedo(RAD_LW)  * (rf%incid_beam(RAD_LW) + rf%incid_diff(RAD_LW))
   end subroutine apply_rt_forcing

end module ed_fast_forcing
