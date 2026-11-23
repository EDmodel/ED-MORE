! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_forcing -- turns the forcing source into each fast step's drivers. The window's
! polygon samples (sample_window), each site's sample after its terrain (site_forcing), and
! per patch and sub-step the air at the patch's canopy-air top and the sky above it
! (patch_forcing; docs/science/forcing.md). The step runs the sky through the canopy itself.
!---~---
module ed_fast_forcing
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_therm_lib,            only : cas_temp_to_enthalpy, precip_phase
   use ed_time,                 only : ed_time_t, time_advance_seconds, seconds_into_day, slope_incidence_cosz, day_of_year,   &
                                       earth_sun_distance_factor
   use ed_met_types,            only : met_source_t, met_cursor_t, met_forcing_t
   use ed_met,                  only : met_advance, met_instant, met_to_cas_top, met_at_site
   use ed_radiation_lib,        only : slope_beam_factor
   use ed_polygon_state_vars,   only : polygon_t
   use ed_fast_integ_types,         only : fast_static_state_t, fast_forcing_t
   implicit none
   private

   public :: sample_window, site_forcing, patch_forcing

contains

   !---~---
   !   This slow step's polygon samples, one per sub-step: met_advance in increasing time, and
   ! met_instant, a pure function of (reader state, time). Made once per window outside the
   ! patch loop: they are the same for every patch, and the file reader stays out of the
   ! parallel region.
   !---~---
   subroutine sample_window(met_src, met_cur, poly, cfg, step_start, met_sample, t_sample)
      type(met_source_t),  intent(in)    :: met_src
      type(met_cursor_t),  intent(inout) :: met_cur
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      type(ed_time_t),     intent(in)    :: step_start      ! calendar time at the start of the slow step
      type(met_forcing_t), intent(out)   :: met_sample(:)   ! (sub-step)
      type(ed_time_t),     intent(out)   :: t_sample(:)     ! (sub-step) each sample's time
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

   !---~---
   !   Each site's sample, moved to its elevation and its beam to its slope (met_at_site,
   ! docs/science/forcing.md §8); a flat site at the polygon's elevation takes it as is.
   !---~---
   subroutine site_forcing(poly, cfg, met_cur, met_sample, t_sample, met_site)
      type(polygon_t),     intent(in)  :: poly
      type(ed_config_t),   intent(in)  :: cfg
      type(met_cursor_t),  intent(in)  :: met_cur
      type(met_forcing_t), intent(in)  :: met_sample(:)   ! (sub-step) the polygon's samples
      type(ed_time_t),     intent(in)  :: t_sample(:)     ! (sub-step) their times
      type(met_forcing_t), intent(out) :: met_site(:,:)   ! (sub-step, site)
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
                                                            surf%slope_deg, surf%aspect_deg),       &
                                       earth_sun_distance_factor(day_of_year(t_sample(isub)))))
               end associate
            end do
         end associate
      end do
   end subroutine site_forcing

   !---~---
   !   The patch's drivers for one sub-step from its site's sample `met`: the air, moved to the
   ! patch's canopy-air top along its log profile when to_top (met_to_cas_top,
   ! docs/science/forcing.md §8; without a forcing source the record already is the air
   ! there), and the sky, which the height does not change.
   !---~---
   subroutine patch_forcing(cfg, static, met, to_top, forcing)
      type(ed_config_t),          intent(in)    :: cfg
      type(fast_static_state_t),  intent(in)    :: static
      type(met_forcing_t),        intent(in)    :: met
      logical,                    intent(in)    :: to_top
      type(fast_forcing_t),       intent(inout) :: forcing   ! sized to the largest patch (ensure_forcing_capacity)
      type(met_forcing_t) :: met_top
      if (to_top) then
         met_top = met_to_cas_top(met, cfg%forcing, static%z_top, static%displace, static%rough)
         call fill_forcing(forcing, met_top)
      else
         call fill_forcing(forcing, met)
      end if
      forcing%sky_live    = .true.
      forcing%cosz        = met%cosz
      forcing%par_beam    = met%par_beam
      forcing%par_diffuse = met%par_diffuse
      forcing%nir_beam    = met%nir_beam
      forcing%nir_diffuse = met%nir_diffuse
      forcing%lw_down     = met%lw_down
   end subroutine patch_forcing

   !---~---
   !   The patch's air from the sub-step's record: its enthalpy, humidity and CO2, the wind,
   ! pressure the patch reads, and the precipitation.
   !---~---
   subroutine fill_forcing(forcing, met)
      type(fast_forcing_t), intent(inout) :: forcing
      type(met_forcing_t),  intent(in)    :: met
      forcing%enthalpy_atm = cas_temp_to_enthalpy(met%air_temp, met%air_shv)
      forcing%shv_atm      = met%air_shv
      forcing%co2_atm      = met%co2
      forcing%wind         = met%wind
      forcing%press        = met%air_pressure
      !---~---
      !   The rain/snow split at the record's temperature: with forcing, the patch's canopy-air top,
      ! where the precipitation's enthalpy is valued too (docs/science/forcing.md §7).
      !---~---
      call precip_phase(met%precip, met%air_temp, forcing%rainfall, forcing%snowfall)
      forcing%air_temp     = met%air_temp                ! rainfall enthalpy reference (snow/rain-on-snow)
   end subroutine fill_forcing

end module ed_fast_forcing
