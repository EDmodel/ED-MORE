! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_met -- the meteorological-forcing reader (docs/science/forcing.md §1-§3, §8, §10-§12).     !
! met_open opens an ED-MORE forcing file or the ED_ERA5land archive into the                        !
! met_source_t a run shares; met_cursor_init binds a polygon's met_cursor_t to it; met_prefetch     !
! loads what a step reads; met_advance and met_instant step and sample a cursor with no file        !
! access; met_at_site moves a sample to a site and met_to_cas_top to a patch's canopy-air top.      !
! read_record is the one ingest for both file sources (ed_met_default,                  !
! ed_met_era5land): it checks, converts, lapses and partitions a record, and a missing        !
! value stops the run. CO2 comes from [forcing].co2_source on model time, never from the met file  !
! (met_co2); a run without forcing opens only the CO2 (met_open_co2).                              !
!==========================================================================================!
module ed_met
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use iso_c_binding, only : c_int
   use ed_kinds, only : wp, ik
   use ed_numerics, only : linear_interpolate
   use ed_therm_lib, only :   air_density, lapse_air_temperature, lapse_pressure, lapse_specific_humidity,        &
                              dry_adiabat_temperature, sat_specific_humidity, rh_to_specific_humidity
   use ed_time, only :   ed_time_t, time_advance_seconds, seconds_into_day,   time_to_string, solar_cosz_utc,     &
                         mean_solar_cosz
   use ed_forcing_params, only :   forcing_params_t, MET_BACKEND_DEFAULT, MET_BACKEND_ERA5LAND,                   &
                                   METAVG_END, METAVG_BEGIN, CLAMP_ERROR, LW_SYNTHESIZE,                          &
                                   CO2_SOURCE_FILE, HEIGHT_ABOVE_ZERO_PLANE, WIND_EXPOSURE_OPEN_TERRAIN
   use ed_radiation_lib,  only :   COSZ_MIN, SWPART_PASSTHROUGH, lapse_longwave, partition_shortwave, clearness_index,      &
                                   synthesize_lwdown
   use ed_met_types, only : met_forcing_t, met_record_t, met_source_t, met_cursor_t,        &
                                  met_cells_t, HUMIDITY_RHAIR, HUMIDITY_TDEW, ED_FIELD
   use ed_met_co2, only : co2_series_read, co2_series_at, co2_series_covers,                 &
                               co2_series_end, co2_series_free
   use ed_netcdf_c, only : nc_close, nc_check
   use ed_met_share, only :   MET_OK, MET_ERR_CO2_FILE, MET_ERR_CO2_NOT_COVERED, file_lookup_sec,    &
                              humidity_field, add_line
   use ed_met_default, only   : default_file_open, default_file_record
   use ed_met_era5land, only : era5land_open, era5land_prefetch,                     &
                                       era5land_record
   implicit none
   private

   public :: met_open, met_open_co2, met_co2, met_cursor_init, met_advance, met_instant, met_close, met_prefetch
   public :: met_at_site, met_to_cas_top, cas_top_wind_factor
   public :: cosz_reconstruct_factor, disaggregate_sw

   integer(ik), parameter :: N_COSZ_SUB   = 10_ik      !< sub-samples per forcing interval for <cosz>_win
   real(wp),    parameter :: COSZ_BAR_MIN = 1.0e-3_wp  !< [-] window-mean-cosz floor for a valid reconstruction

contains

   !=======================================================================================!
   !  OPEN: the CO2, then the source: an ED-MORE forcing file                                 !
   !  (default_file_open) or the archive (era5land_open). A rejection prints the reader's    !
   !  reason, then returns its code with `stat` present, so the tests can exercise it, and    !
   !  stops the run without.                                                                 !
   !  run_start and run_end bound what a non-recycling run reads and the CO2 must cover.      !
   !=======================================================================================!
   subroutine met_open(src, fcfg, stat, run_start, run_end, cells)
      type(met_source_t),     intent(inout)  :: src
      type(forcing_params_t), intent(in)    :: fcfg
      integer(ik), optional,  intent(out)   :: stat
      type(ed_time_t),   optional, intent(in) :: run_start, run_end
      type(met_cells_t), optional, intent(in) :: cells   !< a region's cells (archive only); absent: the site's
      integer(ik)       :: vstat
      character(len=:), allocatable :: why

      if (present(stat)) stat = MET_OK
      why = ''

      src%fcfg       = fcfg
      src%backend    = fcfg%format
      src%grid_index = fcfg%grid_index
      src%dt_forcing = fcfg%timestep
      src%has_wind_vector = .false.

      !----- CO2 first, so every backend takes it the same way. ---------------------------------!
      call open_co2(src, run_start, run_end, vstat, why)
      if (vstat == MET_OK) then
         select case (fcfg%format)
         case (MET_BACKEND_ERA5LAND)
            call era5land_open(src, vstat, why, run_start, run_end, cells)
         case default                                              ! MET_BACKEND_DEFAULT
            call default_file_open(src, vstat, why, run_start, run_end)
         end select
      end if
      if (vstat /= MET_OK) then
         !----- Close the file first: an open netCDF handle keeps an HDF5 lock on the path, and a  !
         !      later open of the same file fails with a bare "permission denied". ---------------!
         call met_close(src)
         write(*,'(a)') why
         if (present(stat)) then
            stat = vstat
            return
         end if
         !----- A fixed stop code after the reason: ifx garbles a stop code that is not a constant. -!
         error stop 'met_open: the forcing cannot drive this run (see the messages above)'
      end if
   end subroutine met_open

   !----- A run without forcing reads only the CO2 its reference climate takes (docs/science/   !
   !      forcing.md §12): the series with co2_source = "file", nothing with "const". A rejection    !
   !      prints the reason and returns its code with `stat` present, and stops the run without. ---!
   subroutine met_open_co2(src, fcfg, run_start, run_end, stat)
      type(met_source_t),     intent(inout) :: src
      type(forcing_params_t), intent(in)    :: fcfg
      type(ed_time_t),        intent(in)    :: run_start, run_end
      integer(ik), optional,  intent(out)   :: stat
      integer(ik) :: vstat
      character(len=:), allocatable :: why
      why = ''
      src%fcfg = fcfg
      call open_co2(src, run_start, run_end, vstat, why)
      if (present(stat)) stat = vstat
      if (vstat == MET_OK) return
      write(*,'(a)') why
      if (present(stat)) return
      error stop 'met_open_co2: the CO2 cannot drive this run (see the messages above)'
   end subroutine met_open_co2

   !----- The free atmosphere's CO2 [umol/mol] at `now`, one way for every run (docs/science/      !
   !      forcing.md eq. (18)): the series with co2_source = "file", else co2_const. --------------!
   pure real(wp) function met_co2(src, now) result(co2)
      type(met_source_t), intent(in) :: src
      type(ed_time_t),    intent(in) :: now
      if (src%fcfg%co2_source == CO2_SOURCE_FILE) then
         co2 = co2_series_at(src%co2, now)
      else
         co2 = src%fcfg%co2_const
      end if
   end function met_co2

   !----- CO2 (docs/science/forcing.md §12): with co2_source = "file", read the ED-MORE CO2     !
   !      file and check it covers the run. ED-MORE does not extrapolate CO2, so a run the file   !
   !      does not cover stops here rather than decades in. The run span is absent only in unit  !
   !      tests. --------------------------------------------------------------------------------!
   subroutine open_co2(src, run_start, run_end, stat, why)
      type(met_source_t),            intent(inout) :: src
      type(ed_time_t), optional,     intent(in)    :: run_start, run_end
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(inout) :: why
      character(len=512) :: msg
      logical :: ok
      stat = MET_OK
      if (src%fcfg%co2_source /= CO2_SOURCE_FILE) return
      call co2_series_read(src%fcfg%co2_file, src%co2, ok, msg)
      if (.not. ok) then
         call add_line(why, ' met_open: [forcing].co2_file = '//trim(src%fcfg%co2_file))
         call add_line(why, '   '//trim(msg))
         stat = MET_ERR_CO2_FILE
         return
      end if
      if (present(run_start) .and. present(run_end)) then
         if (.not. co2_series_covers(src%co2, run_start, run_end)) then
            call add_line(why, ' met_open: the CO2 file does not cover the run: '//trim(src%fcfg%co2_file))
            call add_line(why, '   run  '//time_to_string(run_start)//' .. '//time_to_string(run_end))
            call add_line(why, '   file '//time_to_string(src%co2%first_start)//' .. '//                  &
                               time_to_string(co2_series_end(src%co2)))
            call add_line(why, '   ED-MORE does not extrapolate CO2: extend the file, or shorten the run.')
            stat = MET_ERR_CO2_NOT_COVERED
            return
         end if
      end if
   end subroutine open_co2

   !=======================================================================================!
   !  CURSOR: bind one polygon to the open source -- its cell and location -- and load its first  !
   !  bracket. The grid elevation (the lapse origin) is the cell's orography for the archive and   !
   !  [site].grid_elevation for an ED-MORE forcing file.                                           !
   !=======================================================================================!
   subroutine met_cursor_init(src, cur, cell, latitude_deg, longitude_deg, elevation_m)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(out) :: cur
      integer(ik),        intent(in)  :: cell
      real(wp),           intent(in)  :: latitude_deg, longitude_deg, elevation_m
      cur%cell = cell
      cur%latitude_deg = latitude_deg
      cur%longitude_deg = longitude_deg
      cur%elevation_m  = elevation_m
      cur%grid_elevation_m = src%fcfg%grid_elevation
      if (src%backend == MET_BACKEND_ERA5LAND) cur%grid_elevation_m = src%cells%elevation(cell)
      call load_bracket(src, cur, src%rec_first)
   end subroutine met_cursor_init

   !=======================================================================================!
   !  PREFETCH: load what the step starting at `step_start` reads, before its compute phase.   !
   !  An ED-MORE forcing file was read at open; the archive loads a month at a time.           !
   !=======================================================================================!
   subroutine met_prefetch(src, step_start)
      type(met_source_t), intent(inout)  :: src
      type(ed_time_t),    intent(in)    :: step_start
      if (src%backend == MET_BACKEND_ERA5LAND) call era5land_prefetch(src, step_start)
   end subroutine met_prefetch

   !=======================================================================================!
   !  ADVANCE: slide this cursor's bracket so rec_prev%when <= now < rec_next%when: a start     !
   !  before the first record holds or stops (start_clamp), the recycle seam wraps, and a run   !
   !  past the last record clamps to the last interval (docs/science/forcing.md §2).            !
   !=======================================================================================!
   subroutine met_advance(src, cur, now)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      type(ed_time_t),    intent(in)    :: now
      real(wp)    :: now_sec
      integer(ik) :: i
      logical     :: at_seam, past_end

      !----- The same recycle-mapped seconds on the file axis that met_instant uses. ----------!
      now_sec = file_lookup_sec(src, now)

      !----- The recycle seam: the window's last interval wraps to its first record              !
      !      (docs/science/forcing.md §9). The seam is the window's edge, not the file's. --------!
      !      Nested, not one .and.: Fortran need not short-circuit, and irec_cycle_last is 0 unless recycling.
      at_seam = .false.
      if (src%n_cycle_years >= 1_ik .and. src%irec_cycle_last > src%irec_cycle_first)               &
         at_seam = now_sec >= src%time_sec(src%irec_cycle_last)
      !----- A non-recycling run past the end clamps to the last interval. ----------------------!
      past_end = .not. src%fcfg%recycle .and. src%nrec >= 2_ik .and. now_sec >= src%time_sec(src%nrec)

      if (now_sec < src%time_sec(1)) then
         !----- A start before the first record. Hold loads records 1-2 and met_instant clamps  !
         !      w_next to 0; rec_next is left alone, or a stale copy would stop the reload at t1. !
         if (src%fcfg%start_clamp == CLAMP_ERROR) then
            error stop 'met_advance: model start precedes the first forcing record (start_clamp=error)'
         end if
         if (cur%irec_prev /= 1_ik) call load_bracket(src, cur, 1_ik)
      else if (at_seam) then
         if (.not. cur%at_wrap_seam) call load_wrap_bracket(src, cur)
      else if (past_end) then
         if (cur%irec_prev /= src%nrec - 1_ik) call load_bracket(src, cur, src%nrec - 1_ik)
      else
         !----- Search from the current bracket (cheap for a forward march). ------------------!
         i = max(1_ik, cur%irec_prev)
         do while (i < src%nrec - 1_ik)
            if (src%time_sec(i + 1_ik) <= now_sec) then
               i = i + 1_ik
            else
               exit
            end if
         end do
         do while (i > 1_ik)
            if (src%time_sec(i) > now_sec) then
               i = i - 1_ik
            else
               exit
            end if
         end do
         if (i /= cur%irec_prev) call load_bracket(src, cur, i)
      end if

      !----- The daytime clearness the longwave synthesis holds after dark, on whichever bracket  !
      !      the instant takes: updated here, outside the patch loop, because met_instant is       !
      !      intent(in) (docs/science/forcing.md §11). ---------------------------------------------!
      if (src%fcfg%lwdown_source == LW_SYNTHESIZE) call remember_clearness(src, cur, now)
   end subroutine met_advance

   !----- Update cur%kt_last_day when the sun is up, from the interval's mean shortwave against   !
   !      the sun at `now` (an approximation, docs/science/forcing.md §11). -----------------------!
   subroutine remember_clearness(src, cur, now)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      type(ed_time_t),    intent(in)    :: now
      real(wp) :: cosz_now, sw_total, kt
      type(met_record_t) :: r
      r = interval_mean_record(src, cur)
      cosz_now = solar_cosz_utc(now, seconds_into_day(now), cur%latitude_deg, cur%longitude_deg)
      sw_total = r%par_beam + r%par_diffuse + r%nir_beam + r%nir_diffuse
      kt = clearness_index(sw_total, cosz_now)
      if (kt >= 0.0_wp) cur%kt_last_day = kt         ! negative = night: keep dusk's value
   end subroutine remember_clearness

   !----- The record carrying the means over the interval that contains the model instant: the  !
   !      later of the bracket on an avg_convention = "end" file, the earlier on a "begin" one     !
   !      (docs/science/forcing.md §3). Every interval mean is taken from here. -----------------!
   pure function interval_mean_record(src, cur) result(r)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(met_record_t) :: r
      select case (src%fcfg%avg_convention)
      case (METAVG_BEGIN)
         r = cur%rec_prev
      case default
         r = cur%rec_next            ! METAVG_END (ERA5-Land) and any other
      end select
   end function interval_mean_record

   !=======================================================================================!
   !  INSTANT: the bracket at the model instant `now` (docs/science/forcing.md §3-§5): states     !
   !  linear, wind in energy form, precipitation held, shortwave reconstructed from the interval   !
   !  means, cosz and rho_air derived. A pure function, so it is safe to sample ahead.             !
   !=======================================================================================!
   pure function met_instant(src, cur, now) result(met)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      type(ed_time_t),    intent(in) :: now
      type(met_forcing_t) :: met
      type(met_record_t)  :: mean_rec
      type(ed_time_t)     :: mws
      real(wp) :: now_sec, tprev, tnext, w_next, cosz_now, factor, win_start_sec
      associate (f => src%fcfg, p => cur%rec_prev, n => cur%rec_next)

      !----- The solar zenith at `now` (solar_cosz_utc turns UTC into apparent solar time). ---!
      cosz_now = solar_cosz_utc(now, seconds_into_day(now), cur%latitude_deg, cur%longitude_deg)
      met%cosz = cosz_now

      !----- CO2 on model time, one way for every backend (docs/science/forcing.md eq. (18)). ----!
      met%co2 = met_co2(src, now)

      !----- The weight within the bracket, on the same recycle-mapped seconds as met_advance. -!
      now_sec = file_lookup_sec(src, now)
      if (cur%at_wrap_seam) then                    ! cycle-boundary: window's last rec -> its first
         tprev = src%time_sec(src%irec_cycle_last)
         tnext = tprev + src%dt_forcing
      else
         tprev = src%time_sec(cur%irec_prev)
         tnext = src%time_sec(min(cur%irec_prev + 1_ik, src%nrec))
      end if
      if (tnext > tprev) then
         w_next = min(1.0_wp, max(0.0_wp, (now_sec - tprev) / (tnext - tprev)))
      else
         w_next = 0.0_wp                                          ! degenerate (start-clamp hold)
      end if

      !----- States linear; wind in energy form (docs/science/forcing.md eq. (1)). ----------!
      met%air_temp   = linear_interpolate(p%air_temp,   n%air_temp,   w_next)
      met%air_shv     = linear_interpolate(p%air_shv,     n%air_shv,     w_next)
      met%air_pressure = linear_interpolate(p%air_pressure, n%air_pressure, w_next)
      met%lw_down   = linear_interpolate(p%lw_down,   n%lw_down,   w_next)
      met%wind     = interpolate_wind_energy(p%wind, n%wind, w_next)   ! floored at aerodynamics.ubmin where used
      !----- The vector interpolates linearly, which keeps its direction; never floored. --------!
      if (src%has_wind_vector) then
         met%wind_u = linear_interpolate(p%wind_u, n%wind_u, w_next)
         met%wind_v = linear_interpolate(p%wind_v, n%wind_v, w_next)
         met%has_wind_vector = .true.
      end if

      !----- Rain and shortwave both come from the record whose interval contains `now`. -----!
      mean_rec = interval_mean_record(src, cur)

      !----- Precipitation: the interval's total rate, held across it; each patch splits it. -!
      met%precip = mean_rec%rain

      !----- Shortwave: each mean stream times cosz(now)/<cosz>_win (docs/science/forcing.md      !
      !      eq. (5)), the window anchored on the model calendar (mws), so both see the same sun. !
      mws           = time_advance_seconds(now, tprev - now_sec)   ! model instant at the window start
      win_start_sec = seconds_into_day(mws)
      factor = cosz_reconstruct_factor(mws, win_start_sec,                                       &
                                       (tnext - tprev) / real(N_COSZ_SUB, wp), tnext - tprev,   &
                                       cur%latitude_deg, cur%longitude_deg)
      met%par_beam    = disaggregate_sw(mean_rec%par_beam,    cosz_now, factor)
      met%par_diffuse = disaggregate_sw(mean_rec%par_diffuse, cosz_now, factor)
      met%nir_beam    = disaggregate_sw(mean_rec%nir_beam,    cosz_now, factor)
      met%nir_diffuse = disaggregate_sw(mean_rec%nir_diffuse, cosz_now, factor)

      !----- The longwave synthesis (docs/science/forcing.md §11), after the shortwave: this      !
      !      instant's clearness index by day, the remembered one after dark. ----------------------!
      if (f%lwdown_source == LW_SYNTHESIZE) then
         block
            real(wp) :: sw_now, kt_now
            sw_now = met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse
            kt_now = clearness_index(sw_now, cosz_now)
            if (kt_now < 0.0_wp) kt_now = cur%kt_last_day
            met%lw_down = synthesize_lwdown(f%lw_clear_form, met%air_temp, met%air_shv, met%air_pressure,   &
                                           kt_now, f%lw_cloud_a)
         end block
      end if

      met%air_density = air_density(met%air_temp, met%air_pressure, met%air_shv)
      end associate
   end function met_instant

   !----- Wind interpolates as its square, so its kinetic energy varies linearly              !
   !      (docs/science/forcing.md eq. (1)). No floor here: the aerodynamics floors the speed  !
   !      at aerodynamics.ubmin where it uses it. --------------------------------------------!
   elemental function interpolate_wind_energy(prev, next, w_next) result(u)
      real(wp), intent(in) :: prev, next, w_next
      real(wp) :: u
      u = sqrt(max((1.0_wp - w_next) * prev * prev + w_next * next * next, 0.0_wp))
   end function interpolate_wind_energy

   !=======================================================================================!
   !  SHORTWAVE DISAGGREGATION (docs/science/forcing.md eqs. (5)-(6)): the factor 1/<cosz>_win  !
   !  that turns an interval mean F_avg into F(now) = F_avg*cosz(now)*factor, whose mean over the !
   !  same window is F_avg. It is the reciprocal of the mean cosine (night sub-samples as 0), not  !
   !  the mean secant, which is biased high at low sun. A fully night window returns 0.           !
   !=======================================================================================!
   pure function cosz_reconstruct_factor(t, win_start_sec, dt_sub, dt_win,                      &
                                         latitude_deg, longitude_deg) result(factor)
      type(ed_time_t),   intent(in) :: t
      real(wp),          intent(in) :: win_start_sec, dt_sub, dt_win
      real(wp),          intent(in) :: latitude_deg, longitude_deg
      real(wp)    :: factor, cosz_bar
      cosz_bar = mean_solar_cosz(t, win_start_sec, dt_sub, dt_win, latitude_deg, longitude_deg)
      if (cosz_bar > COSZ_BAR_MIN) then
         factor = 1.0_wp / cosz_bar
      else
         factor = 0.0_wp                                  ! fully-night window; F_avg is 0 too
      end if
   end function cosz_reconstruct_factor

   !----- Instantaneous SW stream from its interval mean: F_now = F_avg*cosz(now)*factor. ------!
   !      Night (cosz(now) <= COSZ_MIN) -> 0 automatically; factor already 0 for a night window. !
   elemental function disaggregate_sw(f_avg, cosz_now, factor) result(f_now)
      real(wp), intent(in) :: f_avg, cosz_now, factor
      real(wp) :: f_now
      if (cosz_now > COSZ_MIN) then
         f_now = f_avg * cosz_now * factor
      else
         f_now = 0.0_wp
      end if
   end function disaggregate_sw

   !=======================================================================================!
   !  AT A SITE (docs/science/forcing.md §8): the polygon's sample taken by the terrain step to a  !
   !  site dz above the polygon, when the run lapses its forcing, and its direct beam scaled by     !
   !  beam_factor (slope_beam_factor) to the site's slope. The diffuse light, the wind and the      !
   !  precipitation are the polygon's.                                                             !
   !=======================================================================================!
   pure function met_at_site(met, f, month, dz, beam_factor) result(site)
      type(met_forcing_t),    intent(in) :: met
      type(forcing_params_t), intent(in) :: f
      integer(ik),            intent(in) :: month
      real(wp),               intent(in) :: dz            !< [m] site above polygon
      real(wp),               intent(in) :: beam_factor   !< [-] slope_beam_factor
      type(met_forcing_t) :: site
      real(wp) :: gamma
      site = met
      if (f%apply_elevation_lapse .and. dz /= 0.0_wp) then
         gamma         = f%lapse_rate_tair(month)
         site%air_pressure = lapse_pressure(met%air_pressure, met%air_temp, dz, gamma)
         site%air_temp   = lapse_air_temperature(met%air_temp, dz, gamma)
         site%air_shv     = lapse_specific_humidity(met%air_shv, met%air_temp, met%air_pressure, site%air_temp, site%air_pressure)
         site%lw_down   = lapse_longwave(met%lw_down, f%lw_clear_form, met%air_temp, met%air_shv, met%air_pressure,     &
                                        site%air_temp, site%air_shv, site%air_pressure)
         site%air_density  = air_density(site%air_temp, site%air_pressure, site%air_shv)
      end if
      site%par_beam = met%par_beam * beam_factor
      site%nir_beam = met%nir_beam * beam_factor
   end function met_at_site

   !=======================================================================================!
   !  AT A PATCH'S CANOPY-AIR TOP (docs/science/forcing.md §8, eq. (15)): from the forcing's own  !
   !  heights, per patch and sub-step.                                                             !
   !=======================================================================================!
   !----- The factor taking the forcing's wind to z_top over a patch of roughness z0 and          !
   !      displacement d: an open-terrain wind is first returned to its blending height z_b, then  !
   !      the patch's neutral log profile carries it to z_top, heights above d floored at 2 z0. ---!
   pure real(wp) function cas_top_wind_factor(f, z_top, displace, rough) result(factor)
      type(forcing_params_t), intent(in) :: f
      real(wp),               intent(in) :: z_top      !< [m] canopy-air top, above the ground
      real(wp),               intent(in) :: displace   !< [m] zero-plane displacement
      real(wp),               intent(in) :: rough      !< [m] roughness length (> 0)
      real(wp) :: z_b, h_b, h_top
      factor = 1.0_wp
      z_b = f%wind_height
      if (f%wind_exposure == WIND_EXPOSURE_OPEN_TERRAIN) then
         factor = log(f%wind_blending_height / f%wind_exposure_z0) / log(f%wind_height / f%wind_exposure_z0)
         z_b    = f%wind_blending_height
      end if
      h_b = z_b
      if (f%height_above /= HEIGHT_ABOVE_ZERO_PLANE) h_b = z_b - displace
      h_b    = max(h_b, 2.0_wp * rough)
      h_top  = max(z_top - displace, 2.0_wp * rough)
      factor = factor * log(h_top / rough) / log(h_b / rough)
   end function cas_top_wind_factor

   !----- A forcing record at a patch's canopy-air top: wind (and its vector) by the patch's factor, !
   !      temperature along the dry adiabat from the forcing's temperature height (a reanalysis's    !
   !      is above d), so the aerodynamics' theta is exactly the forcing's; air density re-derived.  !
   !      Humidity is conserved and the radiation, precipitation and pressure are unchanged; the     !
   !      rain/snow split is made from this record's temperature where the column's forcing is      !
   !      filled (fill_forcing).                                                                     !
   pure function met_to_cas_top(met, f, z_top, displace, rough) result(top)
      type(met_forcing_t),    intent(in) :: met
      type(forcing_params_t), intent(in) :: f
      real(wp),               intent(in) :: z_top      !< [m] canopy-air top, above the ground
      real(wp),               intent(in) :: displace   !< [m] zero-plane displacement
      real(wp),               intent(in) :: rough      !< [m] roughness length (> 0)
      type(met_forcing_t) :: top
      real(wp) :: factor, z_t
      top    = met
      factor = cas_top_wind_factor(f, z_top, displace, rough)
      top%wind = factor * met%wind
      if (met%has_wind_vector) then
         top%wind_u = factor * met%wind_u
         top%wind_v = factor * met%wind_v
      end if
      z_t = f%tq_height
      if (f%height_above == HEIGHT_ABOVE_ZERO_PLANE) z_t = displace + f%tq_height
      top%air_temp  = dry_adiabat_temperature(met%air_temp, z_top - z_t)
      top%air_density = air_density(top%air_temp, top%air_pressure, top%air_shv)
   end function met_to_cas_top

   !=======================================================================================!
   !  CLOSE.                                                                                     !
   !=======================================================================================!
   subroutine met_close(src)
      type(met_source_t), intent(inout)  :: src
      integer(c_int) :: st
      if (src%backend == MET_BACKEND_DEFAULT .and. src%ncid >= 0_ik) then
         st = nc_close(int(src%ncid, c_int))
         call nc_check(st, 'met_close: nc_close')
         src%ncid = -1_ik
      end if
      if (allocated(src%time_sec)) deallocate(src%time_sec)
      if (allocated(src%month_year))  deallocate(src%month_year, src%month_month, src%month_rec0)
      if (allocated(src%buffer%values)) deallocate(src%buffer%values)
      if (allocated(src%carry))  deallocate(src%carry)
      if (allocated(src%head))   deallocate(src%head)
      if (allocated(src%series)) deallocate(src%series)
      src%series_col = 0_ik
      call co2_series_free(src%co2)
      src%buffer%year = 0_ik
      src%carry_rec = 0_ik
      src%n_head = 0_ik
      src%n_loads = 0_ik
   end subroutine met_close

   !----- Load rec_prev = record(irec), rec_next = record(irec+1) (clamped at the end). A slide   !
   !      by one record reuses the old upper record as the new lower one. ------------------------!
   subroutine load_bracket(src, cur, irec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      integer(ik),        intent(in)    :: irec
      if (cur%irec_prev >= 1_ik .and. irec == cur%irec_prev + 1_ik .and. .not. cur%at_wrap_seam) then
         cur%rec_prev = cur%rec_next
      else
         call read_record(src, cur, irec, cur%rec_prev)
      end if
      call read_record(src, cur, min(irec + 1_ik, src%nrec), cur%rec_next)
      cur%irec_prev    = irec
      cur%at_wrap_seam = .false.
   end subroutine load_bracket

   !----- The seam bracket: rec_prev = the window's last record, rec_next = its first, so the   !
   !      last interval interpolates across the wrap instead of clamping. -----------------------!
   subroutine load_wrap_bracket(src, cur)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      call read_record(src, cur, src%irec_cycle_last,  cur%rec_prev)
      call read_record(src, cur, src%irec_cycle_first, cur%rec_next)
      cur%irec_prev    = src%irec_cycle_last
      cur%at_wrap_seam = .true.
   end subroutine load_wrap_bracket

   !----- The one ingest: record irec at this polygon's cell, checked, converted, lapsed and    !
   !      partitioned (docs/science/forcing.md §6-§8). ------------------------------------------!
   subroutine read_record(src, cur, irec, rec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(out)   :: rec
      real(wp)    :: sw_total, humidity_value
      rec%when   = time_advance_seconds(src%base_time, src%time_sec(irec))
      !----- The values as the source stores them; from here on, one ingest for both. -----------!
      if (src%backend == MET_BACKEND_ERA5LAND) then
         call era5land_record(src, cur, irec, rec, humidity_value, sw_total)
      else
         call default_file_record(src, irec, rec, humidity_value, sw_total)
      end if
      call assert_finite(rec%air_temp, 'Tair', irec, src%grid_index)
      call assert_finite(humidity_value, ED_FIELD(humidity_field(src%humidity)), irec, src%grid_index)
      call assert_finite(rec%air_pressure, 'PSurf', irec, src%grid_index)
      !----- Specific humidity at the forcing's own temperature and pressure: a dewpoint and a     !
      !      relative humidity both go through the model's liquid-water saturation curve           !
      !      (docs/science/forcing.md eq. (11)). ---------------------------------------------------!
      select case (src%humidity)
      case (HUMIDITY_RHAIR)
         rec%air_shv = rh_to_specific_humidity(humidity_value, rec%air_temp, rec%air_pressure)
      case (HUMIDITY_TDEW)
         rec%air_shv = sat_specific_humidity(humidity_value, rec%air_pressure)
      case default
         rec%air_shv = humidity_value
      end select
      call assert_finite(rec%wind, 'Wind', irec, src%grid_index)
      call assert_finite(rec%rain, 'Rainf', irec, src%grid_index)
      if (src%fcfg%lwdown_source /= LW_SYNTHESIZE)                                                  &
         call assert_finite(rec%lw_down, 'LWdown', irec, src%grid_index)

      !----- The terrain lapse, from the cell's elevation to the polygon's                        !
      !      (docs/science/forcing.md §8 Step 1, eqs. (13)-(14)). ----------------------------------!
      if (src%fcfg%apply_elevation_lapse) then
         block
            real(wp) :: dz, gamma, t_grid, p_grid, q_grid
            dz     = cur%elevation_m - cur%grid_elevation_m       ! + when the polygon is higher
            gamma  = src%fcfg%lapse_rate_tair(rec%when%month)
            t_grid = rec%air_temp
            p_grid = rec%air_pressure
            q_grid = rec%air_shv   ! before it is overwritten
            rec%air_pressure = lapse_pressure(p_grid, t_grid, dz, gamma)
            rec%air_temp   = lapse_air_temperature(t_grid, dz, gamma)
            !----- Relative humidity is held: a measured one is used as is, any other is formed at the cell. !
            if (src%humidity == HUMIDITY_RHAIR) then
               rec%air_shv  = rh_to_specific_humidity(humidity_value, rec%air_temp, rec%air_pressure)
            else
               rec%air_shv  = lapse_specific_humidity(q_grid, t_grid, p_grid, rec%air_temp, rec%air_pressure)
            end if
            !----- Only a file's longwave: a synthesized one is built later from the lapsed T and q. -!
            if (src%fcfg%lwdown_source /= LW_SYNTHESIZE)                                           &
               rec%lw_down = lapse_longwave(rec%lw_down, src%fcfg%lw_clear_form, t_grid, q_grid, p_grid, &
                                           rec%air_temp, rec%air_shv, rec%air_pressure)
         end block
      end if

      if (src%fcfg%sw_partition == SWPART_PASSTHROUGH) then
         call assert_finite(rec%par_beam,    'SWdown_par_beam',    irec, src%grid_index)   ! required source
         call assert_finite(rec%par_diffuse, 'SWdown_par_diffuse', irec, src%grid_index)   ! fields -> no gap-fill
         call assert_finite(rec%nir_beam,    'SWdown_nir_beam',    irec, src%grid_index)
         call assert_finite(rec%nir_diffuse, 'SWdown_nir_diffuse', irec, src%grid_index)
      else
         call assert_finite(sw_total, 'SWdown', irec, src%grid_index)
         call partition_shortwave(sw_total, record_mean_cosz(src, cur, rec%when), rec%air_pressure,    &
                                  src%fcfg%sw_partition,                                          &
                                  rec%par_beam, rec%par_diffuse, rec%nir_beam, rec%nir_diffuse)
      end if
   end subroutine read_record

   !----- The sun the shortwave partition uses for one record: the mean of max(cosz, 0) over the  !
   !      interval it averages, which ends at the stamp or begins there                         !
   !      (docs/science/forcing.md §6). --------------------------------------------------------!
   real(wp) function record_mean_cosz(src, cur, when) result(cosz)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(ed_time_t),    intent(in) :: when
      real(wp) :: stamp_sec, start_sec, factor
      stamp_sec = seconds_into_day(when)
      start_sec = stamp_sec                                    ! avg_convention = "begin"
      if (src%fcfg%avg_convention == METAVG_END) start_sec = stamp_sec - src%dt_forcing
      factor = cosz_reconstruct_factor(when, start_sec, src%dt_forcing / real(N_COSZ_SUB, wp),       &
                                       src%dt_forcing, cur%latitude_deg, cur%longitude_deg)
      cosz = 0.0_wp
      if (factor > 0.0_wp) cosz = 1.0_wp / factor              ! factor = 1 / <cosz>_interval
   end function record_mean_cosz

   !----- ED-MORE never gap-fills: a NaN in a required field halts the run                       !
   !      (docs/science/forcing.md §10). -----------------------------------------------------!
   subroutine assert_finite(x, name, irec, grid)
      real(wp),         intent(in) :: x
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: irec, grid
      if (ieee_is_nan(x)) then
         write(*,'(4a,i0,a,i0)') 'met_driver: missing/NaN value in required field ', trim(name), &
               ' at time record ', ' ', irec, ', grid ', grid
         error stop 'met_driver: forcing has a data gap (ED-MORE does not gap-fill; fix the source upstream)'
      end if
   end subroutine assert_finite

end module ed_met
