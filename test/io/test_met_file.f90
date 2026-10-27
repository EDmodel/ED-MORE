! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_met_file -- the "default" forcing file and what ed_met does with it
! (docs/science/forcing.md):
!   - its reader: a netCDF round trip, gaps marked by _FillValue or missing_value and a
!     packed variable, the interval that contains an instant, the sunrise
!     interval's shortwave, the file and the [forcing] block agreeing, start-hold and the
!     window edges, recycling a window over many years (Feb 29 included), and the nearest
!     grid point;
!   - the cosz reconstruction that conserves an interval's mean shortwave;
!   - the vertical moves: a site above its polygon (met_at_site) and the forcing height to a
!     canopy-air top (met_to_cas_top);
!   - the CO2 series: its file, its period middles and every refusal;
!   - the contract a flux-tower file relies on: the three humidity forms, UTC, the heights it
!     states, rain and shortwave read from the interval that contains the instant on an end-
!     and a begin-stamped file, and the model's relative humidity returning the file's RHair.
! The laws it reads with are the libraries' tests': test_therm_lib, test_radiation_lib.
!---~---
program test_met_file
   use ed_test_assert,        only : check, check_true, test_report
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : t_3ple, h2o_dry_ratio, grav, cp_air
   use ed_time,               only : ed_time_t, time_advance_seconds, solar_cosz_utc
   use ed_therm_lib,          only : sat_vapor_pressure, air_density, sat_specific_humidity, &
                                     rh_to_specific_humidity, precip_phase, specific_humidity_to_vpd
   use ed_radiation_lib,      only : SWPART_CLEARIDX, SWPART_PASSTHROUGH, clearness_index
   use ed_forcing_params,     only : forcing_params_t, MET_BACKEND_DEFAULT, METAVG_END, METAVG_BEGIN, CLAMP_HOLD, &
                                     CLAMP_ERROR, GRIDMATCH_NEAREST, CO2_SOURCE_FILE, HEIGHT_ABOVE_ZERO_PLANE, &
                                     HEIGHT_ABOVE_GROUND, WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL, &
                                     LW_SYNTHESIZE
   use ed_met_types,          only : met_forcing_t, met_source_t, met_cursor_t, FLD_TAIR, FLD_QAIR, FLD_PSURF
   use ed_met,                only : cosz_reconstruct_factor, disaggregate_sw, cas_top_wind_factor, met_to_cas_top, &
                                     met_open, met_cursor_init, met_advance, met_instant, met_close, met_at_site
   use ed_met_grid,           only : great_circle_distance, nearest_grid_index
   use ed_met_share,          only : MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS, MET_ERR_START_NOT_A_RECORD, &
                                     MET_ERR_WINDOW_NOT_COVERED, MET_ERR_DT_MISMATCH, MET_ERR_ATTR_MISMATCH, &
                                     MET_ERR_CO2_FILE, MET_ERR_CO2_NOT_COVERED, MET_ERR_CO2_IN_MET_FILE, &
                                     MET_ERR_NOT_UTC, MET_ERR_HUMIDITY
   use ed_netcdf_c
   use iso_c_binding,         only : c_int, c_size_t, c_double
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none

   character(len=*), parameter :: NCFILE = 'test_met_file_tmp.nc'
   !---~---
   !   The file of test 7c: Tair's record GAP_REC at its _FillValue, Qair's next record at its
   ! missing_value, and PSurf (98000 Pa) stored packed.
   !---~---
   integer,  parameter :: GAP_REC      = 20
   real(wp), parameter :: TAIR_FILL    = -9999.0_wp
   real(wp), parameter :: QAIR_MISSING = -1.0_wp
   real(wp), parameter :: PSURF_SCALE  = 10.0_wp, PSURF_OFFSET = 500.0_wp
   !--- The forcing file the CO2 tests open: 365 daily records of 2021 (write_yearfile).
   character(len=*), parameter :: CO2_MET_FILE = 'test_met_co2_year_tmp.nc'
   character(len=1024) :: co2_default_file = ''

   character(len=*), parameter :: TOWER_FILE = 'test_met_file_tower_tmp.nc'
   real(wp),    parameter :: DT           = 1800.0_wp    ! [s] half-hourly, as a tower reports
   integer,     parameter :: NT           = 96           ! two days of records
   integer,     parameter :: RAIN_REC     = 30           ! the one record with rain
   real(wp),    parameter :: RAIN_RATE    = 2.0e-3_wp    ! [kg/m2/s] its rate
   integer,     parameter :: SATURATED_REC = 10          ! the record whose RHair is exactly 1
   real(wp),    parameter :: TOWER_HEIGHT = 41.0_wp      ! [m] above the ground (BCI)
   real(wp),    parameter :: PSURF        = 98800.0_wp   ! [Pa]
   type(ed_time_t)        :: base


   call test_cosz_reconstruction()
   call test_netcdf_roundtrip()
   call test_fill_and_packing()
   call test_sunrise_interval()
   call test_nearest_grid()
   call test_cas_top()
   call test_multiyear_cycling()
   call test_recycle_anchor_phase()
   call test_file_config_agreement()
   if (command_argument_count() >= 1) call get_command_argument(1, co2_default_file)
   call test_co2_series(trim(co2_default_file))

   base = ed_time_t(year=2020_ik, month=7_ik, day=1_ik)
   call test_humidity_forms()
   call test_utc_required()
   call test_stated_heights()
   call test_flux_interval(METAVG_END)
   call test_flux_interval(METAVG_BEGIN)
   call test_tower_round_trip()
   call test_met_at_site()

   call test_report('test_met_file')

contains

   !---~---
   !   The site's cursor into an opened source: record-by-record access goes through a cursor
   ! (docs/science/forcing.md §2); a site run has one, at cell 1 and the [site] location.
   !---~---
   subroutine site_cursor(src, cur, fc)
      type(met_source_t),     intent(in)  :: src
      type(met_cursor_t),     intent(out) :: cur
      type(forcing_params_t), intent(in)  :: fc
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude,                   &
                           fc%elevation)
   end subroutine site_cursor







   !--- 5. cosz reconstruction conserves the interval mean (the load-bearing method).
   subroutine test_cosz_reconstruction()
      type(ed_time_t) :: t
      real(wp) :: factor, f_avg, win_start, dt_win, dt_sub, sec, cosz, fsum
      integer(ik) :: i, n
      real(wp), parameter :: LAT = 42.44_wp, LON = -76.50_wp
      print '(a)', '-- test 5: cosz reconstruction (interval-mean conserving) --'
      t = ed_time_t(year=2020_ik, month=7_ik, day=1_ik)
      !---~---
      !   A SUNRISE window (10-11 UTC ~ 05-06 local at Ithaca), where cosz rises steeply from ~0.07
      ! this is where ED2's <sec z> mistake is worst (docs/science/forcing.md §5). Solar noon is ~17 UTC here.
      !---~---
      f_avg = 500.0_wp
      win_start = 10.0_wp*3600.0_wp
      dt_win = 3600.0_wp
      n = 12_ik
      dt_sub = dt_win / real(n, wp)
      factor = cosz_reconstruct_factor(t, win_start, dt_sub, dt_win, LAT, LON)  ! = 1/<cosz>
      !--- (a) the CORRECT 1/<cosz> factor conserves the interval mean EXACTLY over the same subsamples:
      fsum = 0.0_wp
      do i = 1_ik, n
         sec  = win_start + (real(i,wp)-0.5_wp)*dt_sub
         cosz = solar_cosz_utc(t, sec, LAT, LON)
         fsum = fsum + disaggregate_sw(f_avg, cosz, factor)
      end do
      call check('1/<cosz> conserves interval mean (sunrise window)', fsum/real(n,wp), f_avg, 1.0e-8_wp)
      !--- (b) a fully-night window -> factor 0 -> all SW routes to 0
      factor = cosz_reconstruct_factor(t, 4.0_wp*3600.0_wp, dt_sub, dt_win, LAT, LON) ! 04-05 UTC = night
      call check('night window -> factor 0', factor, 0.0_wp, 1.0e-30_wp)
   end subroutine test_cosz_reconstruction

   !--- 7. NetCDF round-trip: write a (time=25, grid=2) file, read it back.
   subroutine test_netcdf_roundtrip()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met_noon, met_night, met_g2
      type(ed_time_t)        :: base, now
      real(wp) :: sw_g1
      print '(a)', '-- test 7: NetCDF round-trip (multi-grid) --'
      base = ed_time_t(year=2020_ik, month=7_ik, day=1_ik)
      call write_synthetic_forcing(NCFILE, base)

      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%grid_index = 1_ik
      fc%timestep = 3600.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = 42.44_wp
      fc%longitude = -76.50_wp
      fc%recycle = .false.
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      call check_true('opened file: nrec=25', src%nrec == 25_ik, real(src%nrec, wp))
      call check_true('opened file: ngrid=2', src%ngrid == 2_ik, real(src%ngrid, wp))
      call check_true('a file with Wind only supplies no wind vector', .not. src%has_wind_vector)

      !--- noon (17 UTC): SW positive, Tair follows the diurnal input
      now = time_advance_seconds(base, 17.0_wp*3600.0_wp)
      call met_advance(src, cur, now)
      met_noon = met_instant(src, cur, now)
      call check_true('noon SWdown > 0', met_noon%sw_down() > 100.0_wp, met_noon%sw_down())
      call check_true('noon cosz > 0.5', met_noon%cosz > 0.5_wp, met_noon%cosz)
      call check_true('Tair in [285,305]', met_noon%air_temp > 285.0_wp .and. met_noon%air_temp < 305.0_wp, met_noon%air_temp)
      sw_g1 = met_noon%sw_down()

      !--- night (04 UTC): SW exactly 0
      now = time_advance_seconds(base, 4.0_wp*3600.0_wp)
      call met_advance(src, cur, now)
      met_night = met_instant(src, cur, now)
      call check('night SWdown = 0', met_night%sw_down(), 0.0_wp, 1.0e-12_wp)
      call check_true('an instant from a Wind-only file has no vector', .not. met_night%has_wind_vector)
      call met_close(src)

      !---~---
      !   A file carrying u10 and v10 supplies the vector, and the speed comes from it
      ! (docs/science/forcing.md §1).
      !---~---
      call write_synthetic_forcing(NCFILE, base, with_components=.true.)
      fc%grid_index = 1_ik
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      call check_true('a file with u10 and v10 supplies the wind vector', src%has_wind_vector)
      now = time_advance_seconds(base, 6.0_wp*3600.0_wp)
      call met_advance(src, cur, now)
      met_night = met_instant(src, cur, now)
      call check('wind_u = u10', met_night%wind_u, 3.0_wp, 1.0e-12_wp)
      call check('wind_v = v10', met_night%wind_v, -4.0_wp, 1.0e-12_wp)
      call check('speed from the components, not Wind', met_night%wind, 5.0_wp, 1.0e-12_wp)
      call met_close(src)
      call write_synthetic_forcing(NCFILE, base)

      !--- multi-grid: grid_index=2 carries a distinct (scaled) SW series
      fc%grid_index = 2_ik
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      now = time_advance_seconds(base, 17.0_wp*3600.0_wp)
      call met_advance(src, cur, now)
      met_g2 = met_instant(src, cur, now)
      call check_true('grid 2 SW differs from grid 1', abs(met_g2%sw_down() - sw_g1) > 1.0_wp, met_g2%sw_down()-sw_g1)
      call met_close(src)

      call test_edge_paths(base)
   end subroutine test_netcdf_roundtrip

   !---~---
   !   7c. A file's gaps and packing, by the netCDF attribute conventions: a value at its
   ! variable's _FillValue or missing_value is a gap (NaN), which stops the step that reads it
   ! (docs/science/forcing.md §10), and a packed variable is unpacked.
   !---~---
   subroutine test_fill_and_packing()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      type(ed_time_t)        :: base, now
      integer(ik) :: jt, jq, jp
      print '(a)', '-- test 7c: gaps by _FillValue and missing_value, and a packed variable --'
      base = ed_time_t(year=2020_ik, month=7_ik, day=1_ik)
      call write_synthetic_forcing(NCFILE, base, conventions=.true.)
      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%grid_index = 1_ik
      fc%timestep = 3600.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = 42.44_wp
      fc%longitude = -76.50_wp
      fc%recycle = .false.
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      jt = src%series_col(FLD_TAIR)
      jq = src%series_col(FLD_QAIR)
      jp = src%series_col(FLD_PSURF)
      call check_true('Tair at its _FillValue reads as a gap', ieee_is_nan(src%series(GAP_REC, jt)))
      call check_true('Qair at its missing_value reads as a gap', ieee_is_nan(src%series(GAP_REC + 1, jq)))
      call check_true('every other record reads as data', count(ieee_is_nan(src%series(:, jt))) == 1 .and.   &
                      count(ieee_is_nan(src%series(:, jq))) == 1)
      call check('PSurf is unpacked: scale_factor, then add_offset', maxval(abs(src%series(:, jp) - 98000.0_wp)), &
                 0.0_wp, 1.0e-9_wp)
      now = time_advance_seconds(base, 3.5_wp*3600.0_wp)          ! a step away from the gaps
      call met_advance(src, cur, now)
      met = met_instant(src, cur, now)
      call check_true('a step away from the gaps reads Tair', met%air_temp > 285.0_wp .and. met%air_temp < 305.0_wp, &
                      met%air_temp)
      call met_close(src)
      call write_synthetic_forcing(NCFILE, base)
   end subroutine test_fill_and_packing

   !---~---
   !   7b. A SUNRISE HOUR KEEPS ITS LIGHT. The shortwave partition's sun is the interval's mean
   ! cosz, so an hour that is lit only after its midpoint still splits its shortwave into the
   ! four streams, and the instants inside it average back to the hour's mean
   ! (docs/science/forcing.md §5-§6). The site longitude is searched so that one hour of the
   ! file has a dark midpoint and a lit end.
   !---~---
   subroutine test_sunrise_interval()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      type(ed_time_t)        :: base, now
      real(wp), parameter    :: SW = 30.0_wp, LAT = 42.44_wp
      integer  :: h, k, il, hsun
      real(wp) :: lon, cz_mid, factor, mean_sw
      print '(a)', '-- test 7b: a sunrise hour keeps its shortwave --'
      base = ed_time_t(year=2020_ik, month=7_ik, day=1_ik)
      hsun = -1
      search: do il = 0, 60
         lon = -76.5_wp - 0.25_wp * real(il, wp)                     ! one minute of sun time per step
         do h = 2, 24                                                  ! the hour ending at stamp h (UTC)
            cz_mid = solar_cosz_utc(base, (real(h, wp) - 0.5_wp) * 3600.0_wp, LAT, lon)
            factor = cosz_reconstruct_factor(base, real(h - 1, wp) * 3600.0_wp, 360.0_wp, 3600.0_wp, LAT, lon)
            if (cz_mid <= 1.0e-3_wp .and. factor > 0.0_wp .and. &
                solar_cosz_utc(base, real(h, wp) * 3600.0_wp, LAT, lon) > 0.02_wp) then
               hsun = h
               exit search
            end if
         end do
      end do search
      call check_true('a sunrise hour with a dark midpoint exists', hsun > 0, real(hsun, wp))
      call write_synthetic_forcing(NCFILE, base, sw_flat=SW)
      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%grid_index = 1_ik
      fc%timestep = 3600.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = LAT
      fc%longitude = lon
      fc%recycle = .false.
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      mean_sw = 0.0_wp
      do k = 1, 60                                                     ! one-minute mid-points
         now = time_advance_seconds(base, (real(hsun - 1, wp) + (real(k, wp) - 0.5_wp) / 60.0_wp) * 3600.0_wp)
         call met_advance(src, cur, now)
         met = met_instant(src, cur, now)
         mean_sw = mean_sw + met%sw_down() / 60.0_wp
      end do
      call check('sunrise hour: its instants average to the hour''s SW', mean_sw, SW, 0.1_wp * SW)
      call met_close(src)
      call write_synthetic_forcing(NCFILE, base)
   end subroutine test_sunrise_interval

   !---~---
   !   13. THE FILE AND THE CONFIG MUST AGREE. The reader checks three things against the
   ! config: the record spacing (dt_forcing), and the two self-describing global attributes the
   ! prep script writes. A file that disagrees with its config would be silently mis-timed or
   ! mis-partitioned -- no later check can see either, because both produce a plausible run.
   !---~---
   subroutine test_file_config_agreement()
      type(met_source_t)     :: src
      type(forcing_params_t) :: fc
      type(ed_time_t)        :: base
      integer(ik)            :: st
      print '(a)', '-- test 13: the forcing file must agree with [forcing] --'
      base = ed_time_t(year=2020_ik, month=7_ik, day=1_ik)
      call write_synthetic_forcing(NCFILE, base)

      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%grid_index = 1_ik
      fc%timestep = 3600.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = 42.44_wp
      fc%longitude = -76.50_wp
      fc%recycle = .false.

      !--- the honest config opens cleanly. A validator that rejects the good case is useless.
      call met_open(src, fc, stat=st)
      call check_true('a matching config opens (MET_OK)', st == MET_OK, real(st, wp))
      call met_close(src)

      !---~---
      !   dt_forcing lying about the cadence. It matters everywhere: dt_forcing places
      ! interval midpoints, disaggregates shortwave and brackets the recycle seam.
      !---~---
      fc%timestep = 1800.0_wp
      call met_open(src, fc, stat=st)
      call check_true('a half-hourly config against an hourly file is rejected',                   &
                      st == MET_ERR_DT_MISMATCH, real(st, wp))
      !--- A rejected open leaves nothing open: a caller's met_close must not close the file twice.
      call check_true('a rejected open releases the file', src%ncid < 0_ik, real(src%ncid, wp))
      call met_close(src)
      fc%timestep = 3600.0_wp

      !--- the file says its flux means end at the stamp; the config says they begin there.
      fc%avg_convention = METAVG_BEGIN
      call met_open(src, fc, stat=st)
      call check_true('avg_convention contradicting the file attribute is rejected',               &
                      st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      fc%avg_convention = METAVG_END

      !--- the file carries TOTAL shortwave; passthrough would read component fields it lacks.
      fc%sw_partition = SWPART_PASSTHROUGH
      call met_open(src, fc, stat=st)
      call check_true('passthrough against a total-shortwave file is rejected',                    &
                      st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      fc%sw_partition = SWPART_CLEARIDX
   end subroutine test_file_config_agreement

   !--- 8. Edge paths: the CLAMP_HOLD start and the recycle-window checks.
   subroutine test_edge_paths(base)
      type(ed_time_t), intent(in) :: base
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      type(ed_time_t)        :: now
      integer(ik) :: st
      real(wp) :: tair0, tair1, expect_hold
      print '(a)', '-- test 8: edge paths (start-hold + recycle) --'
      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%grid_index = 1_ik
      fc%timestep = 3600.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = 42.44_wp
      fc%longitude = -76.50_wp
      !--- synthetic Tair(hh) = 288 + 8 sin(2 pi (hh-15)/24)  -- reproduce the file's values for the checks:
      tair0 = 288.0_wp + 8.0_wp*sin(2.0_wp*3.14159265_wp*(0.0_wp-15.0_wp)/24.0_wp)
      tair1 = 288.0_wp + 8.0_wp*sin(2.0_wp*3.14159265_wp*(1.0_wp-15.0_wp)/24.0_wp)

      !---~---
      !   (a) CLAMP_HOLD: start before the first record holds rec1; then MARCHING into [t0,t1) must
      !     interpolate, not hold rec1 for the whole first interval.
      !---~---
      fc%start_clamp = CLAMP_HOLD
      fc%recycle = .false.
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      now = time_advance_seconds(base, -12.0_wp*3600.0_wp)          ! 12 h before the first record
      call met_advance(src, cur, now)
      met = met_instant(src, cur, now)
      call check('hold before t0 -> record 1', met%air_temp, tair0, 1.0e-4_wp)
      now = time_advance_seconds(base, 0.5_wp*3600.0_wp)            ! 00:30 -> midpoint of [t0,t1)
      call met_advance(src, cur, now)
      met = met_instant(src, cur, now)
      expect_hold = 0.5_wp*(tair0 + tair1)
      call check('march into [t0,t1) interpolates (not stuck at rec1)', met%air_temp, expect_hold, 1.0e-4_wp)
      call met_close(src)

      !---~---
      !   (b) A 25-record (24 h span) file CANNOT be recycled: a sub-year window drifts both
      !     hour-of-day and day-of-year on every wrap: an absolute-seconds span-wrap would
      !     SILENTLY drive a multi-decade run with sub-daily shortwave ~10 h out of phase. It must
      !     be REJECTED (docs/science/forcing.md §9).
      !---~---
      fc%start_clamp = CLAMP_ERROR
      fc%recycle = .true.
      fc%recycle_start = base
      fc%recycle_end   = time_advance_seconds(base, 24.0_wp*3600.0_wp)   ! 1 day, not a whole year
      call met_open(src, fc, stat=st)
      call check_true('sub-year recycle window rejected (not silently span-wrapped)',            &
                      st == MET_ERR_WINDOW_NOT_WHOLE_YEARS, real(st, wp))

      !---~---
      !   (c) A whole-year window whose start does NOT match a record stamp is rejected too
      !     ED-MORE never guesses the start of the day (e.g. a config saying 00:00 against an
      !     end-of-interval ERA5-Land file whose records are stamped 01:00).
      !---~---
      fc%recycle_start = time_advance_seconds(base, 1800.0_wp)           ! 00:30, between records
      fc%recycle_end   = ed_time_t(2021_ik, 7_ik, 1_ik, 0_ik, 30_ik)
      call met_open(src, fc, stat=st)
      call check_true('recycle_start off the record grid rejected',                              &
                      st == MET_ERR_START_NOT_A_RECORD, real(st, wp))

      !--- (d) A whole-year window ON the record grid, but the file only holds 24 h of it.
      fc%recycle_start = base
      fc%recycle_end   = ed_time_t(2021_ik, 7_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('window not covered by the file rejected',                                 &
                      st == MET_ERR_WINDOW_NOT_COVERED, real(st, wp))
      call met_close(src)
   end subroutine test_edge_paths

   !---~---
   !   The recycle mapping must preserve hour-of-day and day-of-year EXACTLY, for an anchor
   ! anywhere in the calendar -- not only Jan-1. A classifier that accepts only Jan-1 00:00
   ! files and span-wraps everything else shifts hour-of-day whenever the span is not a whole
   ! number of days (docs/science/forcing.md §9, eq. (16)).
   !---~---
   subroutine test_recycle_anchor_phase()
      character(len=*), parameter :: YF = 'test_met_anchorfile_tmp.nc'
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: m_ref, m_map
      integer(ik) :: st
      print '(a)', '-- test: recycle anchor + sub-daily phase --'
      call write_yearfile(YF, 2021_ik)                       ! 365 daily records from 2021-01-01 00:00
      fc%format = MET_BACKEND_DEFAULT
      fc%path = YF
      fc%grid_index = 1_ik
      fc%timestep = 86400.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .true.

      !---~---
      !   A MID-YEAR window: 2021-03-01 .. 2022-03-01. Plain year substitution would map a
      ! January instant to 2021-01, which is OUTSIDE this window; the anchored form sends it
      ! to 2022-01 instead. The file only reaches 2021-12-31, so this window is (correctly)
      ! rejected as not covered -- which is itself the check that the coverage test bites.
      !---~---
      fc%recycle_start = ed_time_t(2021_ik,   3_ik, 1_ik)
      fc%recycle_end   = ed_time_t(2022_ik,   3_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('mid-year window beyond the file is rejected',                             &
                      st == MET_ERR_WINDOW_NOT_COVERED, real(st, wp))

      !---~---
      !   Whole-year window on the file: every model year must read the SAME calendar day, for
      ! model years both before and after the window (the modulo must not go negative).
      !---~---
      fc%recycle_start = ed_time_t(2021_ik,   1_ik, 1_ik)
      fc%recycle_end   = ed_time_t(2022_ik,   1_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('whole-year window accepted', st == MET_OK, real(st, wp))
      call site_cursor(src, cur, fc)
      call met_advance(src, cur, ed_time_t(2021_ik,9_ik,17_ik))
      m_ref = met_instant(src, cur, ed_time_t(2021_ik,9_ik,17_ik))
      call met_advance(src, cur, ed_time_t(2049_ik,9_ik,17_ik))
      m_map = met_instant(src, cur, ed_time_t(2049_ik,9_ik,17_ik))
      call check('29 yr later reads the same calendar day', m_map%air_temp, m_ref%air_temp, 1.0e-9_wp)
      call met_advance(src, cur, ed_time_t(2013_ik,9_ik,17_ik))
      m_map = met_instant(src, cur, ed_time_t(2013_ik,9_ik,17_ik))
      call check('a model year BEFORE the window maps in too', m_map%air_temp, m_ref%air_temp, 1.0e-9_wp)
      call met_close(src)
   end subroutine test_recycle_anchor_phase

   !--- Write a synthetic (time=25 hourly, grid=2) ED-MORE forcing NetCDF via ed_netcdf_c.
   subroutine write_synthetic_forcing(path, base, with_components, with_co2air, sw_flat, conventions)
      character(len=*),  intent(in) :: path
      type(ed_time_t),   intent(in) :: base
      logical, optional, intent(in) :: with_components   ! also write u10 = 3, v10 = -4 (speed 5, not Wind's 3)
      logical, optional, intent(in) :: with_co2air       ! also write CO2air, which the reader rejects
      real(wp), optional, intent(in) :: sw_flat          ! every record's SWdown [W/m2], night and day
      logical, optional, intent(in) :: conventions       ! the gaps and packing of test 7c
      integer, parameter :: NT = 25, NG = 2
      integer(c_int) :: st, ncid, td, gd, vt, vla, vlo, vv(7), vu, vw, vc
      integer(c_int) :: dims2(2), dims1(1)
      real(c_double) :: tsec(NT), lat(NG), lon(NG)
      real(c_double) :: dat(NG, NT)        ! (grid, time) column-major == C [time][grid]
      integer :: it, ig, k
      real(wp) :: hh, sw_mean
      character(len=8), parameter :: VNAMES(7) = ['Tair    ','Qair    ','PSurf   ','Wind    ', &
                                                  'Rainf   ','SWdown  ','LWdown  ']
      integer(c_size_t) :: start2(2), count2(2), start1(1), count1(1)
      logical :: co2air, packed

      co2air = .false.
      if (present(with_co2air)) co2air = with_co2air
      packed = .false.
      if (present(conventions)) packed = conventions
      st = nc_create_f(path, NC_NETCDF4, ncid)
      call nc_check(st, 'write: create')
      st = nc_def_dim_f(ncid, 'time', int(NT, c_size_t), td)
      call nc_check(st, 'write: time dim')
      st = nc_def_dim_f(ncid, 'grid', int(NG, c_size_t), gd)
      call nc_check(st, 'write: grid dim')
      dims1(1) = td
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt)
      call nc_check(st, 'write: time var')
      st = nc_put_att_text_f(ncid, vt, 'units', int(len_trim('seconds since 2020-07-01 00:00:00'), c_size_t), &
                             'seconds since 2020-07-01 00:00:00')
                             call nc_check(st, 'write: units')
      dims1(1) = gd
      st = nc_def_var_f(ncid, 'latitude',  NC_DOUBLE, 1, dims1, vla)
      call nc_check(st, 'write: lat var')
      st = nc_def_var_f(ncid, 'longitude', NC_DOUBLE, 1, dims1, vlo)
      call nc_check(st, 'write: lon var')
      dims2(1) = td
      dims2(2) = gd                        ! [time, grid] (slowest first)
      do k = 1, 7
         st = nc_def_var_f(ncid, trim(VNAMES(k)), NC_DOUBLE, 2, dims2, vv(k))
         call nc_check(st, 'write: var '//trim(VNAMES(k)))
      end do
      if (co2air) then
         st = nc_def_var_f(ncid, 'CO2air', NC_DOUBLE, 2, dims2, vc)
         call nc_check(st, 'write: var CO2air')
      end if
      if (packed) then
         st = nc_put_att_double_f(ncid, vv(1), '_FillValue', NC_DOUBLE, TAIR_FILL)
         call nc_check(st, 'write: fill')
         st = nc_put_att_double_f(ncid, vv(2), 'missing_value', NC_DOUBLE, QAIR_MISSING)
         call nc_check(st, 'write: missing_value')
         st = nc_put_att_double_f(ncid, vv(3), 'scale_factor', NC_DOUBLE, PSURF_SCALE)
         call nc_check(st, 'write: scale')
         st = nc_put_att_double_f(ncid, vv(3), 'add_offset', NC_DOUBLE, PSURF_OFFSET)
         call nc_check(st, 'write: offset')
      end if
      if (present(with_components)) then
         if (with_components) then
            st = nc_def_var_f(ncid, 'u10', NC_DOUBLE, 2, dims2, vu)
            call nc_check(st, 'write: var u10')
            st = nc_def_var_f(ncid, 'v10', NC_DOUBLE, 2, dims2, vw)
            call nc_check(st, 'write: var v10')
         end if
      end if
      !---~---
      !   The two SELF-DESCRIBING global attributes the prep script writes. The fixture has to
      ! carry them or it is not the file the reader validates against, and the file/config
      ! agreement checks would be skipped here while firing in production -- the
      ! fixture-must-mirror-the-driver trap. This file is hourly, ERA5-Land-style, with total
      ! shortwave.
      !---~---
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'avg_convention', int(len_trim('end'), c_size_t), 'end')
      call nc_check(st, 'write: avg_convention')
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'sw_input_kind', int(len_trim('total'), c_size_t), 'total')
      call nc_check(st, 'write: sw_input_kind')
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'time_zone', 3_c_size_t, 'UTC')
      call nc_check(st, 'write: time_zone')
      st = nc_enddef(ncid)
      call nc_check(st, 'write: enddef')

      do it = 1, NT
         tsec(it) = real(it - 1, c_double) * 3600.0_c_double
      end do
      lat = [42.5_c_double, 42.4_c_double]
      lon = [-76.6_c_double, -76.5_c_double]
      start1(1) = 0_c_size_t
      count1(1) = int(NT, c_size_t)
      st = nc_put_vara_double(ncid, vt, start1, count1, tsec)
      call nc_check(st, 'write: time vals')
      count1(1) = int(NG, c_size_t)
      st = nc_put_vara_double(ncid, vla, start1, count1, lat)
      call nc_check(st, 'write: lat vals')
      st = nc_put_vara_double(ncid, vlo, start1, count1, lon)
      call nc_check(st, 'write: lon vals')

      start2 = [0_c_size_t, 0_c_size_t]
      count2 = [int(NT, c_size_t), int(NG, c_size_t)]
      do k = 1, 7
         do it = 1, NT
            hh = real(it - 1, wp)                      ! UTC hour 0..24
            do ig = 1, NG
               select case (k)
               case (1)
                  dat(ig,it) = 288.0_wp + 8.0_wp*sin(2.0_wp*3.14159265_wp*(hh-15.0_wp)/24.0_wp)         ! Tair
               case (2)
                  dat(ig,it) = 0.008_wp                ! Qair
               case (3)
                  dat(ig,it) = 98000.0_wp              ! PSurf
               case (4)
                  dat(ig,it) = 3.0_wp                  ! Wind
               case (5)
                  dat(ig,it) = 0.0_wp                  ! Rainf
               case (6)                                 ! SWdown: daytime bump, scaled per grid cell
                  sw_mean = max(0.0_wp, 900.0_wp*sin(3.14159265_wp*(hh-11.0_wp)/12.0_wp))
                  if (hh < 11.0_wp .or. hh > 23.0_wp) sw_mean = 0.0_wp
                  dat(ig,it) = sw_mean * (1.0_wp + 0.3_wp*real(ig-1, wp))
                  if (present(sw_flat)) dat(ig,it) = sw_flat
               case (7)
                  dat(ig,it) = 320.0_wp                ! LWdown
               end select
            end do
         end do
         if (packed) then
            select case (k)
            case (1)
               dat(:, GAP_REC)     = TAIR_FILL
            case (2)
               dat(:, GAP_REC + 1) = QAIR_MISSING
            case (3)
               dat = (dat - PSURF_OFFSET) / PSURF_SCALE
            end select
         end if
         st = nc_put_vara_double(ncid, vv(k), start2, count2, dat)
         call nc_check(st, 'write: vals '//trim(VNAMES(k)))
      end do
      if (present(with_components)) then
         if (with_components) then
            dat = 3.0_wp
            st = nc_put_vara_double(ncid, vu, start2, count2, dat)
            call nc_check(st, 'write: vals u10')
            dat = -4.0_wp
            st = nc_put_vara_double(ncid, vw, start2, count2, dat)
            call nc_check(st, 'write: vals v10')
         end if
      end if
      if (co2air) then
         dat = 415.0_wp
         st = nc_put_vara_double(ncid, vc, start2, count2, dat)
         call nc_check(st, 'write: vals CO2air')
      end if
      st = nc_close(ncid)
      call nc_check(st, 'write: close')
   end subroutine write_synthetic_forcing

   !--- Nearest-grid match: pure kernel (argmin + great-circle) and the reader override.
   subroutine test_nearest_grid()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      real(wp)    :: lon(3), lat(3), d
      integer(ik) :: idx
      print '(a)', '-- test: nearest-grid match --'
      lon = [-80.0_wp, -76.5_wp, -70.0_wp]
      lat = [40.0_wp, 42.44_wp, 45.0_wp]
      call check('argmin picks the coincident cell', real(nearest_grid_index(-76.5_wp,42.44_wp,lon,lat),wp), 2.0_wp, 0.5_wp)
      call check('argmin picks the SW cell',         real(nearest_grid_index(-79.0_wp,40.5_wp,lon,lat),wp), 1.0_wp, 0.5_wp)
      call check('great-circle self-distance = 0', great_circle_distance(-76.5_wp,42.44_wp,-76.5_wp,42.44_wp), 0.0_wp, 1.0e-6_wp)
      call check('great-circle 1 deg lat ~ 111 km', great_circle_distance(0.0_wp,0.0_wp,0.0_wp,1.0_wp), 111195.0_wp, 500.0_wp)
      !--- reader override: the 2-grid file has lon/lat cells (-76.6,42.5) and (-76.5,42.4).
      call write_synthetic_forcing(NCFILE, ed_time_t(2020_ik,7_ik,1_ik))
      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%timestep = 3600.0_wp
      fc%grid_match = GRIDMATCH_NEAREST
      fc%grid_index = 1_ik                                          ! grid_index deliberately WRONG for cell 2
      fc%latitude = 42.41_wp
      fc%longitude = -76.49_wp                              ! nearest to cell 2
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      call check('reader resolves nearest -> grid_index 2', real(src%grid_index,wp), 2.0_wp, 0.5_wp)
      call met_close(src)
      fc%latitude = 42.51_wp
      fc%longitude = -76.61_wp                              ! nearest to cell 1
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      call check('reader resolves nearest -> grid_index 1', real(src%grid_index,wp), 1.0_wp, 0.5_wp)
      call met_close(src)
   end subroutine test_nearest_grid

   !---~---
   !   The move to a patch's canopy-air top (docs/science/forcing.md §8). A 25 m canopy:
   ! d = 15.75 m, z0 = 3.25 m, top z_c = 30 m.
   !---~---
   subroutine test_cas_top()
      real(wp), parameter :: ZC = 30.0_wp, D = 15.75_wp, Z0 = 3.25_wp, GC = grav / cp_air
      type(forcing_params_t) :: f
      type(met_forcing_t)    :: m, top
      real(wp) :: fac, rain_top, snow_top, rain_ref, snow_ref
      print '(a)', '-- test: forcing to the canopy-air top --'
      !--- ERA5-Land: open-terrain 10 m wind, heights above the zero plane.
      f%tq_height = 2.0_wp
      f%wind_height = 10.0_wp
      f%height_above = HEIGHT_ABOVE_ZERO_PLANE
      f%wind_exposure = WIND_EXPOSURE_OPEN_TERRAIN
      f%wind_exposure_z0 = 0.03_wp
      f%wind_blending_height = 40.0_wp
      fac = cas_top_wind_factor(f, ZC, D, Z0)
      call check('ERA5 wind: open-terrain step, then the patch profile', fac,                     &
                 log(40.0_wp/0.03_wp) / log(10.0_wp/0.03_wp) * log((ZC - D)/Z0) / log(40.0_wp/Z0), 1.0e-12_wp)
      call check('ERA5 wind: ~0.73 u10 at the top of a 25 m canopy', fac, 0.729_wp, 1.0e-3_wp)
      m%air_temp = 290.0_wp
      top = met_to_cas_top(m, f, ZC, D, Z0)
      call check('temperature conserves theta from d + 2 m', top%air_temp, 290.0_wp - GC * (ZC - (D + 2.0_wp)), 1.0e-12_wp)
      call check('theta at the top equals theta at the forcing height',                            &
                 top%air_temp + GC * ZC, 290.0_wp + GC * (D + 2.0_wp), 1.0e-9_wp)
      !--- A flux tower: local wind at 40 m above the ground; the 2 z0 floor on a tall stand.
      f%height_above = HEIGHT_ABOVE_GROUND
      f%wind_exposure = WIND_EXPOSURE_LOCAL
      f%tq_height = 40.0_wp
      f%wind_height = 40.0_wp
      call check('tower wind: the patch profile from 40 m above the ground', cas_top_wind_factor(f, ZC, D, Z0), &
                 log((ZC - D)/Z0) / log((40.0_wp - D)/Z0), 1.0e-12_wp)
      top = met_to_cas_top(m, f, ZC, D, Z0)
      call check('tower temperature: along the dry adiabat from 40 m', top%air_temp, 290.0_wp + GC * 10.0_wp, 1.0e-12_wp)
      call check('tower at the top: no change', cas_top_wind_factor(f, 40.0_wp, D, Z0), 1.0_wp, 1.0e-15_wp)
      call check('heights floored at 2 z0 above d', cas_top_wind_factor(f, D + 1.0_wp, D, Z0),       &
                 log(2.0_wp) / log((40.0_wp - D)/Z0), 1.0e-12_wp)
      !--- The whole record: wind and vector scaled, temperature moved, rho re-derived, the rest kept.
      f%height_above = HEIGHT_ABOVE_ZERO_PLANE
      f%wind_exposure = WIND_EXPOSURE_OPEN_TERRAIN
      f%tq_height = 2.0_wp
      f%wind_height = 10.0_wp
      m%wind = 5.0_wp
      m%wind_u = 3.0_wp
      m%wind_v = -4.0_wp
      m%has_wind_vector = .true.
      m%air_temp = 290.0_wp
      m%air_shv = 0.009_wp
      m%air_pressure = 98000.0_wp
      top = met_to_cas_top(m, f, ZC, D, Z0)
      fac = cas_top_wind_factor(f, ZC, D, Z0)
      call check('record: wind speed scaled', top%wind, 5.0_wp * fac, 1.0e-12_wp)
      call check('record: vector scaled (u)', top%wind_u, 3.0_wp * fac, 1.0e-12_wp)
      call check('record: vector scaled (v)', top%wind_v, -4.0_wp * fac, 1.0e-12_wp)
      call check('record: temperature moved', top%air_temp, 290.0_wp - GC * (ZC - (D + 2.0_wp)), 1.0e-12_wp)
      call check('record: humidity conserved', top%air_shv, 0.009_wp, 1.0e-15_wp)
      call check('record: pressure stays at the ground', top%air_pressure, 98000.0_wp, 1.0e-12_wp)
      call check('record: shortwave unchanged', top%sw_down(), m%sw_down(), 1.0e-12_wp)
      call check('record: air density re-derived', top%air_density, air_density(top%air_temp, 98000.0_wp, 0.009_wp), 1.0e-12_wp)
      !---~---
      !   The record carries total precipitation, unchanged by the move; the rain/snow split is made
      ! from the moved temperature where the column's forcing is filled. At freezing, the cooler top
      ! air turns some of the rain to snow.
      !---~---
      m%air_temp = t_3ple
      m%precip = 1.0e-3_wp
      top = met_to_cas_top(m, f, ZC, D, Z0)
      call check('record: total precipitation unchanged', top%precip, 1.0e-3_wp, 1.0e-15_wp)
      call precip_phase(m%precip, m%air_temp, rain_ref, snow_ref)
      call precip_phase(top%precip, top%air_temp, rain_top, snow_top)
      call check_true('record: at freezing, the cooler top air turns some rain to snow', snow_top > snow_ref, &
                      snow_top - snow_ref)
   end subroutine test_cas_top

   !--- Multi-year CALENDAR recycling + Feb-29 reconciliation (whole-year daily file).
   subroutine test_multiyear_cycling()
      character(len=*), parameter :: YF = 'test_met_yearfile_tmp.nc'
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: m_ref, m_map
      print '(a)', '-- test: multi-year cycling + Feb-29 --'
      call write_yearfile(YF, 2021_ik)                             ! 365 daily records, Jan-1 2021 aligned
      fc%format = MET_BACKEND_DEFAULT
      fc%path = YF
      fc%grid_index = 1_ik
      fc%timestep = 86400.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .true.
      !--- The cycle is DECLARED, never inferred from the file.
      fc%recycle_start = ed_time_t(2021_ik,   1_ik, 1_ik)
      fc%recycle_end   = ed_time_t(2022_ik,   1_ik, 1_ik)
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      call check('n_cycle_years = 1', real(src%n_cycle_years,wp), 1.0_wp, 0.5_wp)
      call check('cycle anchor year = 2021', real(src%cycle_anchor%year,wp), 2021.0_wp, 0.5_wp)
      call check('cycle first record = record 1', real(src%irec_cycle_first,wp), 1.0_wp, 0.5_wp)
      !--- recycle identity: model 2023-07-01 reads the 2021-07-01 record (day-of-year exact).
      call met_advance(src, cur, ed_time_t(2021_ik,7_ik,1_ik))
      m_ref = met_instant(src, cur, ed_time_t(2021_ik,7_ik,1_ik))
      call met_advance(src, cur, ed_time_t(2023_ik,7_ik,1_ik))
      m_map = met_instant(src, cur, ed_time_t(2023_ik,7_ik,1_ik))
      call check('recycle: 2023-07-01 reads the 2021-07-01 record', m_map%air_temp, m_ref%air_temp, 1.0e-9_wp)
      !--- Feb-29 reconciliation: model 2024-02-29 (leap) maps to file 2021-02-28.
      call met_advance(src, cur, ed_time_t(2021_ik,2_ik,28_ik))
      m_ref = met_instant(src, cur, ed_time_t(2021_ik,2_ik,28_ik))
      call met_advance(src, cur, ed_time_t(2024_ik,2_ik,29_ik))
      m_map = met_instant(src, cur, ed_time_t(2024_ik,2_ik,29_ik))
      call check('Feb-29 (leap model) maps to file Feb-28', m_map%air_temp, m_ref%air_temp, 1.0e-9_wp)
      call met_close(src)
   end subroutine test_multiyear_cycling

   !---~---
   !   15. Prescribed CO2 (docs/science/forcing.md §12, eq. (18)): co2_const on every backend,
   ! or an ED-MORE CO2 file read on MODEL time -- each value at its period's middle, linear
   ! between middles, the end values held over the outer half-periods, a run the file does not
   ! cover rejected at open.
   !---~---
   subroutine test_co2_series(default_file)
      character(len=*), intent(in) :: default_file      ! the shipped series ('' -> its check is skipped)
      character(len=*), parameter :: CF = 'test_met_co2_tmp.txt'
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met, m_ref, m_map
      integer(ik)            :: st
      print '(a)', '-- test 15: prescribed CO2 --'

      !--- "const": co2_const (not the type's 420 default).
      call write_yearfile(CO2_MET_FILE, 2021_ik)
      fc = co2_test_forcing()
      fc%co2_const = 284.3_wp
      call met_open(src, fc)
      call site_cursor(src, cur, fc)
      met = met_instant(src, cur, ed_time_t(2021_ik, 7_ik, 1_ik))
      call check('const: CO2 = co2_const', met%co2, 284.3_wp, 1.0e-12_wp)
      call met_close(src)

      !---~---
      !   Annual means, with comments and a blank line. 2000 is a leap year, so its middle is
      ! Jul 2 00:00; 2001's is Jul 2 12:00. On 2001-01-01 the weight is 183 d / 365.5 d.
      !---~---
      call write_text(CF, [character(len=48) :: '# a comment', 'timestep  1 year', 'units umol/mol', &
                           '2000  300.0   # a trailing comment', '2001  310.0', '', '2002  330.0'])
      call open_file_co2(CF, ed_time_t(2000_ik,1_ik,1_ik), ed_time_t(2003_ik,1_ik,1_ik), src, cur, st)
      call check_true('year: the file opens and covers 2000-01-01 .. 2003-01-01', st == MET_OK, real(st, wp))
      call check('year: at the middle of 2001, the 2001 value', co2_now(src, cur,                    &
                 ed_time_t(2001_ik,7_ik,2_ik,12_ik)), 310.0_wp, 1.0e-9_wp)
      call check('year: linear between middles', co2_now(src, cur, ed_time_t(2001_ik,1_ik,1_ik)), &
                 300.0_wp + 10.0_wp * 183.0_wp / 365.5_wp, 1.0e-9_wp)
      call check('year: held before the first middle', co2_now(src, cur,                           &
                 ed_time_t(2000_ik,3_ik,1_ik)), 300.0_wp, 1.0e-12_wp)
      call check('year: held after the last middle', co2_now(src, cur,                             &
                 ed_time_t(2002_ik,12_ik,31_ik,23_ik)), 330.0_wp, 1.0e-12_wp)
      call met_close(src)
      call open_file_co2(CF, ed_time_t(1999_ik,12_ik,31_ik), ed_time_t(2003_ik,1_ik,1_ik), src, cur, st)
      call check_true('year: a run starting before the file is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))
      call open_file_co2(CF, ed_time_t(2000_ik,1_ik,1_ik), ed_time_t(2003_ik,1_ik,1_ik,1_ik), src, cur, st)
      call check_true('year: a run ending after the file is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))

      !---~---
      !   Monthly, across a leap February: Jan's middle is the 16th 12:00 (15.5 d), Feb's the
      ! 15th 12:00 (14.5 d), so Feb 1 00:00 is 15.5/30 of the way.
      !---~---
      call write_text(CF, [character(len=48) :: 'timestep 1 month', 'units umol/mol', '2024-01 400', &
                           '2024-02 410', '2024-03 420'])
      call open_file_co2(CF, ed_time_t(2024_ik,1_ik,1_ik), ed_time_t(2024_ik,4_ik,1_ik), src, cur, st)
      call check_true('month: opens', st == MET_OK, real(st, wp))
      call check('month: at the middle of a leap February', co2_now(src, cur,                       &
                 ed_time_t(2024_ik,2_ik,15_ik,12_ik)), 410.0_wp, 1.0e-9_wp)
      call check('month: Feb 1, linear by calendar days', co2_now(src, cur,                         &
                 ed_time_t(2024_ik,2_ik,1_ik)), 400.0_wp + 10.0_wp * 15.5_wp / 30.0_wp, 1.0e-9_wp)
      call met_close(src)

      !--- 5-year means from the year 1000: both spans are 1826 d, so 1005-01-01 is halfway.
      call write_text(CF, [character(len=48) :: 'timestep 5 year', 'units umol/mol', '1000 280', '1005 290'])
      call open_file_co2(CF, ed_time_t(1000_ik,1_ik,1_ik), ed_time_t(1010_ik,1_ik,1_ik), src, cur, st)
      call check_true('5 year: opens and covers 1000 .. 1010', st == MET_OK, real(st, wp))
      call check('5 year: halfway between the two middles', co2_now(src, cur,                       &
                 ed_time_t(1005_ik,1_ik,1_ik)), 285.0_wp, 1.0e-9_wp)
      call met_close(src)

      !--- Daily, and 30-minute (a flux tower's cadence).
      call write_text(CF, [character(len=48) :: 'timestep 1 day', 'units umol/mol', '2024-02-28 400', &
                           '2024-02-29 402'])
      call open_file_co2(CF, ed_time_t(2024_ik,2_ik,28_ik), ed_time_t(2024_ik,3_ik,1_ik), src, cur, st)
      call check_true('day: opens (a leap day row)', st == MET_OK, real(st, wp))
      call check('day: midnight between two daily means', co2_now(src, cur,                         &
                 ed_time_t(2024_ik,2_ik,29_ik)), 401.0_wp, 1.0e-9_wp)
      call met_close(src)
      call write_text(CF, [character(len=48) :: 'timestep 30 minute', 'units umol/mol',             &
                           '2024-07-01T00:00 400', '2024-07-01T00:30 402', '2024-07-01T01:00 404'])
      call open_file_co2(CF, ed_time_t(2024_ik,7_ik,1_ik), ed_time_t(2024_ik,7_ik,1_ik,1_ik,30_ik), src, cur, st)
      call check_true('30 minute: opens and covers 00:00 .. 01:30', st == MET_OK, real(st, wp))
      call check('30 minute: between the 00:15 and 00:45 middles', co2_now(src, cur,                &
                 ed_time_t(2024_ik,7_ik,1_ik,0_ik,30_ik)), 401.0_wp, 1.0e-9_wp)
      call check('30 minute: held over the last half-period', co2_now(src, cur,                     &
                 ed_time_t(2024_ik,7_ik,1_ik,1_ik,25_ik)), 404.0_wp, 1.0e-12_wp)
      call met_close(src)
      call open_file_co2(CF, ed_time_t(2024_ik,7_ik,1_ik), ed_time_t(2024_ik,7_ik,1_ik,1_ik,31_ik), src, cur, st)
      call check_true('30 minute: a run past 01:30 is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))

      !--- Files that break format 1 are rejected, never guessed at.
      call expect_bad('no timestep line', [character(len=48) :: 'units umol/mol', '2000 300', '2001 310'])
      call expect_bad('a second timestep line', [character(len=48) :: 'timestep 1 year', 'timestep 1 year', &
                      'units umol/mol', '2000 300', '2001 310'])
      call expect_bad('an unknown unit', [character(len=48) :: 'timestep 1 week', 'units umol/mol', &
                      '2000 300', '2001 310'])
      call expect_bad('units other than umol/mol', [character(len=48) :: 'timestep 1 year', 'units ppm', &
                      '2000 300', '2001 310'])
      call expect_bad('a keyword after the data', [character(len=48) :: 'timestep 1 year', 'units umol/mol', &
                      '2000 300', '2001 310', 'units umol/mol'])
      call expect_bad('a start at the wrong precision', [character(len=48) :: 'timestep 1 year',       &
                      'units umol/mol', '2000-01 300', '2001-01 310'])
      call expect_bad('an empty date part', [character(len=48) :: 'timestep 1 month', 'units umol/mol', &
                      '2000-01 300', '2000- 310'])
      call expect_bad('a gap', [character(len=48) :: 'timestep 1 year', 'units umol/mol', '2000 300', '2002 310'])
      call expect_bad('one row', [character(len=48) :: 'timestep 1 year', 'units umol/mol', '2000 300'])
      call expect_bad('a value that is not a number', [character(len=48) :: 'timestep 1 year',         &
                      'units umol/mol', '2000 300', '2001 abc'])
      call expect_bad('a value that is not positive', [character(len=48) :: 'timestep 1 year',         &
                      'units umol/mol', '2000 300', '2001 -5'])
      call expect_bad('a third column', [character(len=48) :: 'timestep 1 year', 'units umol/mol',     &
                      '2000 300 1', '2001 310 1'])
      call open_file_co2('no_such_co2_file.txt', ed_time_t(2000_ik,1_ik,1_ik),                      &
                         ed_time_t(2001_ik,1_ik,1_ik), src, cur, st)
      call check_true('bad file: a missing file is rejected', st == MET_ERR_CO2_FILE, real(st, wp))

      !---~---
      !   The point of reading CO2 on model time: recycled met repeats, the CO2 does not. Model
      ! 2023-07-02 12:00 reads the 2021 met record but the 2023 CO2.
      !---~---
      call write_text(CF, [character(len=48) :: 'timestep 1 year', 'units umol/mol', '2021 400',    &
                           '2022 410', '2023 420', '2024 430'])
      fc = co2_test_forcing()
      fc%recycle = .true.
      fc%recycle_start = ed_time_t(2021_ik, 1_ik, 1_ik)
      fc%recycle_end = ed_time_t(2022_ik, 1_ik, 1_ik)
      fc%co2_source = CO2_SOURCE_FILE
      fc%co2_file = CF
      call met_open(src, fc, stat=st, run_start=ed_time_t(2021_ik,1_ik,1_ik),                      &
                    run_end=ed_time_t(2025_ik,1_ik,1_ik))
      call check_true('recycle: opens with a CO2 file', st == MET_OK, real(st, wp))
      if (st == MET_OK) then
         call site_cursor(src, cur, fc)
         call met_advance(src, cur, ed_time_t(2021_ik,7_ik,2_ik,12_ik))
         m_ref = met_instant(src, cur, ed_time_t(2021_ik,7_ik,2_ik,12_ik))
         call met_advance(src, cur, ed_time_t(2023_ik,7_ik,2_ik,12_ik))
         m_map = met_instant(src, cur, ed_time_t(2023_ik,7_ik,2_ik,12_ik))
         call check('recycle: the met repeats (2023 reads the 2021 record)', m_map%air_temp, m_ref%air_temp, 1.0e-9_wp)
         call check('recycle: 2021 CO2', m_ref%co2, 400.0_wp, 1.0e-9_wp)
         call check('recycle: the CO2 does not repeat (2023 CO2)', m_map%co2, 420.0_wp, 1.0e-9_wp)
         call met_close(src)
      end if

      !--- A met file that carries CO2air is rejected: nothing would read it.
      call write_synthetic_forcing(NCFILE, ed_time_t(year=2020_ik, month=7_ik, day=1_ik), with_co2air=.true.)
      fc = forcing_params_t()
      fc%format = MET_BACKEND_DEFAULT
      fc%path = NCFILE
      fc%grid_index = 1_ik
      fc%timestep = 3600.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .false.
      call met_open(src, fc, stat=st)
      call check_true('CO2air: a met file carrying it is rejected', st == MET_ERR_CO2_IN_MET_FILE, real(st, wp))
      call write_synthetic_forcing(NCFILE, ed_time_t(year=2020_ik, month=7_ik, day=1_ik))

      !---~---
      !   The shipped default series (data/co2/): it loads, covers 1000 .. 2023, and reads
      ! back the CMIP7 values at the years' middles (1700 and 1850 are not leap years).
      !---~---
      if (len_trim(default_file) == 0) then
         print '(a)', '   NOTE: no path to the shipped CO2 file was given; its check is skipped'
      else
         call open_file_co2(default_file, ed_time_t(1000_ik,1_ik,1_ik), ed_time_t(2023_ik,1_ik,1_ik), &
                            src, cur, st)
         call check_true('shipped: loads and covers 1000-01-01 .. 2023-01-01', st == MET_OK, real(st, wp))
         if (st == MET_OK) then
            call check('shipped: 1000 (held)', co2_now(src, cur, ed_time_t(1000_ik,3_ik,1_ik)), 282.437_wp, 1.0e-9_wp)
            call check('shipped: 1700', co2_now(src, cur, ed_time_t(1700_ik,7_ik,2_ik,12_ik)),   277.537_wp, 1.0e-9_wp)
            call check('shipped: 1850', co2_now(src, cur, ed_time_t(1850_ik,7_ik,2_ik,12_ik)),   284.297_wp, 1.0e-9_wp)
            call check('shipped: 2022 (held)', co2_now(src, cur, ed_time_t(2022_ik,12_ik,31_ik)), 417.320_wp, 1.0e-9_wp)
            call met_close(src)
         end if
         call open_file_co2(default_file, ed_time_t(1000_ik,1_ik,1_ik), ed_time_t(2023_ik,1_ik,2_ik), &
                            src, cur, st)
         call check_true('shipped: a run past 2022 is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))
      end if
   end subroutine test_co2_series

   !--- Open CO2_MET_FILE with co2_file = path over the run [t0, t1]; `st` is met_open's.
   subroutine open_file_co2(path, t0, t1, src, cur, st)
      character(len=*),   intent(in)    :: path
      type(ed_time_t),    intent(in)    :: t0, t1
      type(met_source_t), intent(inout) :: src
      type(met_cursor_t), intent(out)   :: cur
      integer(ik),        intent(out)   :: st
      type(forcing_params_t) :: fc
      fc = co2_test_forcing()
      fc%co2_source = CO2_SOURCE_FILE
      fc%co2_file = path
      call met_open(src, fc, stat=st, run_start=t0, run_end=t1)
      if (st == MET_OK) call site_cursor(src, cur, fc)
   end subroutine open_file_co2

   !---~---
   !   The CO2 tests' [forcing]: CO2_MET_FILE, read whole, not recycled. The CO2 is read on model
   ! time, so the file's year need not be the run's.
   !---~---
   function co2_test_forcing() result(fc)
      type(forcing_params_t) :: fc
      fc%format = MET_BACKEND_DEFAULT
      fc%path = CO2_MET_FILE
      fc%grid_index = 1_ik
      fc%timestep = 86400.0_wp
      fc%avg_convention = METAVG_END
      fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .false.
   end function co2_test_forcing

   real(wp) function co2_now(src, cur, t) result(co2)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(ed_time_t),    intent(in) :: t
      type(met_forcing_t) :: met
      !---~---
      !   A series that failed to open has no rows. -1 fails every check that follows, where
      ! met_instant would read its unallocated arrays and the segfault would swallow the
      ! buffered output that says why the open failed.
      !---~---
      if (src%co2%n < 2_ik) then
         co2 = -1.0_wp
         return
      end if
      met = met_instant(src, cur, t)
      co2 = met%co2
   end function co2_now

   !--- A CO2 file with `lines` must be rejected at open (MET_ERR_CO2_FILE).
   subroutine expect_bad(what, lines)
      character(len=*), intent(in) :: what, lines(:)
      character(len=*), parameter  :: BF = 'test_met_co2_bad_tmp.txt'
      type(met_source_t) :: src
      type(met_cursor_t) :: cur
      integer(ik)        :: st
      call write_text(BF, lines)
      call open_file_co2(BF, ed_time_t(2000_ik,1_ik,1_ik), ed_time_t(2001_ik,1_ik,1_ik), src, cur, st)
      call check_true('bad file: '//what, st == MET_ERR_CO2_FILE, real(st, wp))
      if (st == MET_OK) call met_close(src)
   end subroutine expect_bad

   subroutine write_text(path, lines)
      character(len=*), intent(in) :: path, lines(:)
      integer :: u, j
      open(newunit=u, file=path, status='replace', action='write')
      do j = 1, size(lines)
         write(u, '(a)') trim(lines(j))
      end do
      close(u)
   end subroutine write_text

   !--- Write a whole-year DAILY (365 records, grid=1) forcing file; Tair encodes the day index.
   subroutine write_yearfile(path, year)
      character(len=*), intent(in) :: path
      integer(ik),      intent(in) :: year
      integer, parameter :: NT = 365, NG = 1
      integer(c_int)    :: st, ncid, td, gd, vt, vla, vlo, vv(7), dims2(2), dims1(1)
      integer(c_size_t) :: start2(2), count2(2), start1(1), count1(1)
      real(c_double)    :: tsec(NT), la(NG), lo(NG), dat(NG, NT)
      integer :: it, k
      character(len=8), parameter :: VNAMES(7) = ['Tair    ','Qair    ','PSurf   ','Wind    ', &
                                                  'Rainf   ','SWdown  ','LWdown  ']
      character(len=40) :: units
      write(units,'(a,i0,a)') 'seconds since ', year, '-01-01 00:00:00'
      st = nc_create_f(path, NC_NETCDF4, ncid)
      call nc_check(st, 'yf: create')
      st = nc_def_dim_f(ncid, 'time', int(NT,c_size_t), td)
      call nc_check(st, 'yf: time dim')
      st = nc_def_dim_f(ncid, 'grid', int(NG,c_size_t), gd)
      call nc_check(st, 'yf: grid dim')
      dims1(1) = td
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt)
      call nc_check(st, 'yf: time var')
      st = nc_put_att_text_f(ncid, vt, 'units', int(len_trim(units),c_size_t), trim(units))
      call nc_check(st, 'yf: units')
      dims1(1) = gd
      st = nc_def_var_f(ncid, 'latitude',  NC_DOUBLE, 1, dims1, vla)
      call nc_check(st, 'yf: lat var')
      st = nc_def_var_f(ncid, 'longitude', NC_DOUBLE, 1, dims1, vlo)
      call nc_check(st, 'yf: lon var')
      dims2(1) = td
      dims2(2) = gd
      do k = 1, 7
         st = nc_def_var_f(ncid, trim(VNAMES(k)), NC_DOUBLE, 2, dims2, vv(k))
         call nc_check(st, 'yf: var')
      end do
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'time_zone', 3_c_size_t, 'UTC')
      call nc_check(st, 'yf: time_zone')
      st = nc_enddef(ncid)
      call nc_check(st, 'yf: enddef')
      do it = 1, NT
         tsec(it) = real(it-1, c_double) * 86400.0_c_double
      end do
      la = 42.44_c_double
      lo = -76.50_c_double
      start1(1) = 0_c_size_t
      count1(1) = int(NT,c_size_t)
      st = nc_put_vara_double(ncid, vt, start1, count1, tsec)
      call nc_check(st, 'yf: time vals')
      count1(1) = int(NG,c_size_t)
      st = nc_put_vara_double(ncid, vla, start1, count1, la)
      call nc_check(st, 'yf: lat vals')
      st = nc_put_vara_double(ncid, vlo, start1, count1, lo)
      call nc_check(st, 'yf: lon vals')
      start2 = [0_c_size_t, 0_c_size_t]
      count2 = [int(NT,c_size_t), int(NG,c_size_t)]
      do k = 1, 7
         do it = 1, NT
            select case (k)
            case (1)
               dat(1,it) = 280.0_wp + real(it-1, wp)           ! Tair encodes the day index (unique per day)
            case (2)
               dat(1,it) = 0.006_wp                            ! Qair
            case (3)
               dat(1,it) = 99000.0_wp                          ! PSurf
            case (4)
               dat(1,it) = 2.0_wp                              ! Wind
            case (5)
               dat(1,it) = 0.0_wp                              ! Rainf
            case (6)
               dat(1,it) = 200.0_wp                            ! SWdown (daily mean)
            case (7)
               dat(1,it) = 300.0_wp                            ! LWdown
            end select
         end do
         st = nc_put_vara_double(ncid, vv(k), start2, count2, dat)
         call nc_check(st, 'yf: vals')
      end do
      st = nc_close(ncid)
      call nc_check(st, 'yf: close')
   end subroutine write_yearfile


   !---~---
   !   The [forcing] block a tower file is read with: BCI's location, 41 m heights above the
   ! ground, a wind measured over this canopy.
   !---~---
   function tower_config(convention) result(fc)
      integer(ik), intent(in) :: convention
      type(forcing_params_t)  :: fc
      fc%format = MET_BACKEND_DEFAULT
      fc%path = TOWER_FILE
      fc%grid_index = 1_ik
      fc%timestep = DT
      fc%avg_convention = convention
      fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = 9.1568_wp
      fc%longitude = -79.8486_wp
      fc%elevation = 150.0_wp
      fc%recycle = .false.
      fc%tq_height = TOWER_HEIGHT
      fc%wind_height = TOWER_HEIGHT
      fc%height_above = HEIGHT_ABOVE_GROUND
      fc%wind_exposure = WIND_EXPOSURE_LOCAL
   end function tower_config

   !--- The fixture's values at record k (1-based), so a test can name what it expects.
   pure real(wp) function tair_at(k) result(v)
      integer, intent(in) :: k
      v = 299.0_wp + 3.0_wp * sin(2.0_wp * 3.14159265_wp * (0.5_wp * real(k - 1, wp) - 13.0_wp) / 24.0_wp)
   end function tair_at

   pure real(wp) function rh_at(k) result(v)
      integer, intent(in) :: k
      v = 0.80_wp - 0.15_wp * sin(2.0_wp * 3.14159265_wp * (0.5_wp * real(k - 1, wp) - 13.0_wp) / 24.0_wp)
      if (k == SATURATED_REC) v = 1.0_wp
   end function rh_at

   !---~---
   !   Distinct daytime means, rising record by record, so a read from the wrong interval is off
   ! by 10 W/m2 rather than by a rounding error. The daylit records at BCI are 11:00-23:30 UTC.
   !---~---
   pure real(wp) function sw_at(k) result(v)
      integer, intent(in) :: k
      real(wp) :: hh
      hh = modulo(0.5_wp * real(k - 1, wp), 24.0_wp)
      v = 0.0_wp
      if (hh >= 13.0_wp .and. hh <= 21.0_wp) v = 100.0_wp + 10.0_wp * real(k, wp)
   end function sw_at

   pure real(wp) function stamp_seconds(k) result(s)
      integer, intent(in) :: k
      s = real(k - 1, wp) * DT
   end function stamp_seconds

   !---~---
   !   A one-point ED-MORE forcing file, half-hourly from 2020-07-01 00:00 UTC. `humidity` names the
   ! variable(s) it carries ('RHair', 'Tdew', 'Qair', 'RHair+Qair' or 'none'); an empty
   ! `time_zone` or `height_above` leaves that attribute out, and so does stated_height <= 0.
   !---~---
   subroutine write_tower_file(convention, humidity, time_zone, rh_scale, stated_height, height_above)
      character(len=*), intent(in) :: convention, humidity, time_zone, height_above
      real(wp),         intent(in) :: rh_scale, stated_height
      character(len=8)  :: names(9)
      integer(c_int)    :: st, ncid, td, gd, vt, vla, vlo, vid(9), dims1(1), dims2(2)
      integer(c_size_t) :: s1(1), c1(1), s2(2), c2(2)
      real(c_double)    :: tsec(NT), dat(1, NT), one(1)
      integer :: nvar, j, k
      names(1:6) = [character(len=8) :: 'Tair', 'PSurf', 'Wind', 'Rainf', 'SWdown', 'LWdown']
      nvar = 6
      select case (humidity)
      case ('RHair')
         nvar = 7
         names(7) = 'RHair'
      case ('Tdew')
         nvar = 7
         names(7) = 'Tdew'
      case ('Qair')
         nvar = 7
         names(7) = 'Qair'
      case ('RHair+Qair')
         nvar = 8
         names(7) = 'RHair'
         names(8) = 'Qair'
      end select
      st = nc_create_f(TOWER_FILE, NC_NETCDF4, ncid)
      call nc_check(st, 'tower: create')
      st = nc_def_dim_f(ncid, 'time', int(NT, c_size_t), td)
      call nc_check(st, 'tower: time dim')
      st = nc_def_dim_f(ncid, 'grid', 1_c_size_t, gd)
      call nc_check(st, 'tower: grid dim')
      dims1(1) = td
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt)
      call nc_check(st, 'tower: time var')
      st = nc_put_att_text_f(ncid, vt, 'units', 33_c_size_t, 'seconds since 2020-07-01 00:00:00')
      call nc_check(st, 'tower: units')
      dims1(1) = gd
      st = nc_def_var_f(ncid, 'latitude',  NC_DOUBLE, 1, dims1, vla)
      call nc_check(st, 'tower: lat')
      st = nc_def_var_f(ncid, 'longitude', NC_DOUBLE, 1, dims1, vlo)
      call nc_check(st, 'tower: lon')
      dims2(1) = td
      dims2(2) = gd
      do j = 1, nvar
         st = nc_def_var_f(ncid, trim(names(j)), NC_DOUBLE, 2, dims2, vid(j))
         call nc_check(st, 'tower: var')
      end do
      !---~---
      !   The global attributes make_tower_forcing.py writes; the fixture must carry them or the
      ! checks under test would be skipped here and fire in production.
      !---~---
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'avg_convention', int(len(convention), c_size_t), convention)
      call nc_check(st, 'tower: avg_convention')
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'sw_input_kind', 5_c_size_t, 'total')
      call nc_check(st, 'tower: sw_input_kind')
      if (len(time_zone) > 0) then
         st = nc_put_att_text_f(ncid, NC_GLOBAL, 'time_zone', int(len(time_zone), c_size_t), time_zone)
         call nc_check(st, 'tower: time_zone')
      end if
      if (stated_height > 0.0_wp) then
         st = nc_put_att_double_f(ncid, NC_GLOBAL, 'tq_height_m', NC_DOUBLE, real(stated_height, c_double))
         call nc_check(st, 'tower: tq_height_m')
         st = nc_put_att_double_f(ncid, NC_GLOBAL, 'wind_height_m', NC_DOUBLE, real(stated_height, c_double))
         call nc_check(st, 'tower: wind_height_m')
      end if
      if (len(height_above) > 0) then
         st = nc_put_att_text_f(ncid, NC_GLOBAL, 'height_above', int(len(height_above), c_size_t), height_above)
         call nc_check(st, 'tower: height_above')
      end if
      st = nc_enddef(ncid)
      call nc_check(st, 'tower: enddef')

      do k = 1, NT
         tsec(k) = real(stamp_seconds(k), c_double)
      end do
      s1(1) = 0_c_size_t
      c1(1) = int(NT, c_size_t)
      st = nc_put_vara_double(ncid, vt, s1, c1, tsec)
      call nc_check(st, 'tower: time vals')
      c1(1) = 1_c_size_t
      one = 9.1568_c_double
      st = nc_put_vara_double(ncid, vla, s1, c1, one)
      call nc_check(st, 'tower: lat vals')
      one = -79.8486_c_double
      st = nc_put_vara_double(ncid, vlo, s1, c1, one)
      call nc_check(st, 'tower: lon vals')
      s2 = 0_c_size_t
      c2 = [int(NT, c_size_t), 1_c_size_t]
      do j = 1, nvar
         do k = 1, NT
            select case (trim(names(j)))
            case ('Tair')
               dat(1, k) = tair_at(k)
            case ('PSurf')
               dat(1, k) = PSURF
            case ('Wind')
               dat(1, k) = 2.5_wp
            case ('Rainf')
               dat(1, k) = merge(RAIN_RATE, 0.0_wp, k == RAIN_REC)
            case ('SWdown')
               dat(1, k) = sw_at(k)
            case ('LWdown')
               dat(1, k) = 450.0_wp
            case ('RHair')
               dat(1, k) = rh_scale * rh_at(k)
            case ('Tdew')
               dat(1, k) = tair_at(k) - 3.0_wp
            case ('Qair')
               dat(1, k) = 0.018_wp
            end select
         end do
         st = nc_put_vara_double(ncid, vid(j), s2, c2, dat)
         call nc_check(st, 'tower: vals')
      end do
      st = nc_close(ncid)
      call nc_check(st, 'tower: close')
   end subroutine write_tower_file

   !--- The forcing at an instant `sec` seconds after the file's first stamp.
   function sample(src, cur, sec) result(met)
      type(met_source_t), intent(in)    :: src
      type(met_cursor_t), intent(inout) :: cur
      real(wp),           intent(in)    :: sec
      type(met_forcing_t) :: met
      type(ed_time_t)     :: now
      now = time_advance_seconds(base, sec)
      call met_advance(src, cur, now)
      met = met_instant(src, cur, now)
   end function sample

   !---~---
   !   1. Each humidity form becomes q through the model's own conversion, at the record's own
   ! temperature and pressure; a file with none, with two, or with RHair in percent is refused.
   !---~---
   subroutine test_humidity_forms()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      integer(ik) :: st
      integer, parameter :: K = 40
      print '(a)', '-- test 1: humidity forms --'
      fc = tower_config(METAVG_BEGIN)

      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('an RHair file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude, fc%elevation)
      met = sample(src, cur, stamp_seconds(K))
      call check('RHair -> q by rh_to_specific_humidity at the record''s T and P', met%air_shv,           &
                 rh_to_specific_humidity(rh_at(K), tair_at(K), PSURF), 1.0e-15_wp)
      call met_close(src)

      call write_tower_file('begin', 'Tdew', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a Tdew file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude, fc%elevation)
      met = sample(src, cur, stamp_seconds(K))
      call check('Tdew -> q by sat_specific_humidity', met%air_shv,                              &
                 sat_specific_humidity(tair_at(K) - 3.0_wp, PSURF), 1.0e-15_wp)
      call met_close(src)

      call write_tower_file('begin', 'Qair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a Qair file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude, fc%elevation)
      met = sample(src, cur, stamp_seconds(K))
      call check('Qair is used as it is', met%air_shv, 0.018_wp, 1.0e-15_wp)
      call met_close(src)

      call write_tower_file('begin', 'RHair+Qair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a file with two humidity variables is refused', st == MET_ERR_HUMIDITY, real(st, wp))
      call met_close(src)

      call write_tower_file('begin', 'none', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a file with no humidity variable is refused', st == MET_ERR_HUMIDITY, real(st, wp))
      call met_close(src)

      call write_tower_file('begin', 'RHair', 'UTC', 100.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('RHair written in percent is refused', st == MET_ERR_HUMIDITY, real(st, wp))
      call met_close(src)
   end subroutine test_humidity_forms

   !--- 2. ED-MORE runs in UTC: a file that does not say time_zone = "UTC" is refused.
   subroutine test_utc_required()
      type(met_source_t)     :: src
      type(forcing_params_t) :: fc
      integer(ik) :: st
      print '(a)', '-- test 2: the UTC requirement --'
      fc = tower_config(METAVG_BEGIN)
      call write_tower_file('begin', 'RHair', '', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a file without time_zone is refused', st == MET_ERR_NOT_UTC, real(st, wp))
      call met_close(src)
      call write_tower_file('begin', 'RHair', 'UTC-05:00', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a local-time file is refused', st == MET_ERR_NOT_UTC, real(st, wp))
      call met_close(src)
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a UTC file opens', st == MET_OK, real(st, wp))
      call met_close(src)
   end subroutine test_utc_required

   !--- 3. The heights a file states must be the ones [forcing] moves the forcing from.
   subroutine test_stated_heights()
      type(met_source_t)     :: src
      type(forcing_params_t) :: fc
      integer(ik) :: st
      print '(a)', '-- test 3: the stated heights --'
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      fc = tower_config(METAVG_BEGIN)
      fc%tq_height = 2.0_wp
      call met_open(src, fc, stat=st)
      call check_true('tq_height against a file stating 41 m is refused', st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      call met_close(src)
      fc = tower_config(METAVG_BEGIN)
      fc%wind_height = 10.0_wp
      call met_open(src, fc, stat=st)
      call check_true('wind_height against a file stating 41 m is refused', st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      call met_close(src)
      fc = tower_config(METAVG_BEGIN)
      fc%height_above = HEIGHT_ABOVE_ZERO_PLANE
      call met_open(src, fc, stat=st)
      call check_true('height_above against a file stating "ground" is refused',                     &
                      st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      call met_close(src)
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, 0.0_wp, '')
      fc = tower_config(METAVG_BEGIN)
      fc%tq_height = 2.0_wp
      call met_open(src, fc, stat=st)
      call check_true('a file that states no heights is not checked', st == MET_OK, real(st, wp))
      call met_close(src)
   end subroutine test_stated_heights

   !---~---
   !   4. Rain and shortwave are interval means, read from the record whose interval CONTAINS
   ! the instant: record k covers [t_k, t_k+dt) on a begin-stamped file and (t_k-dt, t_k] on an
   ! end-stamped one. A read one record off puts every storm in the wrong half-hour; the
   ! shortwave check averages the instantaneous flux over the reconstruction's own ten
   ! midpoints, which returns the interval mean exactly when the right record is used.
   !---~---
   subroutine test_flux_interval(convention)
      integer(ik), intent(in) :: convention
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      real(wp) :: t0, sw_mean
      integer  :: i, k_sw
      integer(ik) :: st
      if (convention == METAVG_BEGIN) then
         print '(a)', '-- test 4b: rain and shortwave intervals, begin-stamped --'
         call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
         t0 = stamp_seconds(RAIN_REC)                        ! the rain interval starts at its stamp
      else
         print '(a)', '-- test 4a: rain and shortwave intervals, end-stamped --'
         call write_tower_file('end', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
         t0 = stamp_seconds(RAIN_REC) - DT                   ! ... ends at its stamp
      end if
      fc = tower_config(convention)
      call met_open(src, fc, stat=st)
      call check_true('the file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude, fc%elevation)
      met = sample(src, cur, t0 - 0.5_wp * DT)
      call check('no rain in the interval before', met%precip, 0.0_wp, 1.0e-15_wp)
      met = sample(src, cur, t0 + 0.25_wp * DT)
      call check('the rain falls in its own interval (early)', met%precip, RAIN_RATE, 1.0e-15_wp)
      met = sample(src, cur, t0 + 0.75_wp * DT)
      call check('the rain falls in its own interval (late)', met%precip, RAIN_RATE, 1.0e-15_wp)
      met = sample(src, cur, t0 + 1.5_wp * DT)
      call check('no rain in the interval after', met%precip, 0.0_wp, 1.0e-15_wp)

      !--- A midday interval, 17:00-17:30 UTC (about solar noon at BCI).
      t0 = 17.0_wp * 3600.0_wp
      k_sw = nint(t0 / DT) + 1                               ! begin: the record stamped 17:00
      if (convention == METAVG_END) k_sw = k_sw + 1          ! end:   the record stamped 17:30
      sw_mean = 0.0_wp
      do i = 1, 10
         met = sample(src, cur, t0 + (real(i, wp) - 0.5_wp) * DT / 10.0_wp)
         sw_mean = sw_mean + (met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse) / 10.0_wp
      end do
      call check('the shortwave of the interval is its own record''s mean', sw_mean, sw_at(k_sw), 1.0e-9_wp)
      call met_close(src)

      !---~---
      !   The clearness index the longwave synthesis holds through the night is the containing
      ! interval's too: on a begin file, the next record's would be the dark one after sunset.
      !---~---
      fc%lwdown_source = LW_SYNTHESIZE
      call met_open(src, fc, stat=st)
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude, fc%elevation)
      met = sample(src, cur, t0 + 0.4_wp * DT)
      call check('the remembered clearness is the containing interval''s', cur%kt_last_day,              &
                 clearness_index(sw_at(k_sw), solar_cosz_utc(time_advance_seconds(base, t0 + 0.4_wp * DT),  &
                                 t0 + 0.4_wp * DT, fc%latitude, fc%longitude)), 1.0e-12_wp)
      call met_close(src)
   end subroutine test_flux_interval

   !--- 5. The round trip of a tower file.
   subroutine test_tower_round_trip()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      real(wp) :: vpd, rh_model, worst
      integer  :: k
      integer(ik) :: st
      print '(a)', '-- test 5: the round trip of a tower file --'
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      fc = tower_config(METAVG_BEGIN)
      call met_open(src, fc, stat=st)
      call check_true('the tower file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude, fc%elevation)

      !---~---
      !   At every stamp the model's own inverse returns the tower's relative humidity. This is
      ! what storing RHair buys: the provider's saturation curve never enters.
      !---~---
      worst = 0.0_wp
      do k = 1, NT - 1
         if (k == SATURATED_REC) cycle
         met = sample(src, cur, stamp_seconds(k))
         vpd = specific_humidity_to_vpd(met%air_temp, met%air_shv, met%air_pressure)
         rh_model = 1.0_wp - vpd / sat_vapor_pressure(met%air_temp)
         worst = max(worst, abs(rh_model - rh_at(k)))
      end do
      call check('the model''s RH at the forcing temperature is RHair at every stamp', worst, 0.0_wp, 1.0e-12_wp)
      met = sample(src, cur, stamp_seconds(SATURATED_REC))
      vpd = specific_humidity_to_vpd(met%air_temp, met%air_shv, met%air_pressure)
      call check('a saturated tower record reads back as VPD = 0', vpd, 0.0_wp, 1.0e-9_wp)

      call met_close(src)
   end subroutine test_tower_round_trip

   !---~---
   !   A site's forcing within its polygon (ed_met::met_at_site; docs/science/forcing.md §8): a
   ! site above its polygon is colder by the lapse rate times the height, at a lower pressure
   ! and the same relative humidity, under less longwave; one at the polygon's elevation on
   ! level ground keeps the polygon's forcing exactly; without the lapse only the beam scales.
   !---~---
   subroutine test_met_at_site()
      type(forcing_params_t) :: f
      type(met_forcing_t)    :: met, site
      print '(a)', '-- a site within its polygon --'
      !--- 3. A site above its polygon.
      f%apply_elevation_lapse = .true.
      f%lapse_rate_tair = 0.0065_wp
      met%air_temp = 290.0_wp
      met%air_pressure = 95000.0_wp
      met%air_shv = 0.010_wp
      met%lw_down = 350.0_wp
      met%par_beam = 200.0_wp
      met%nir_beam = 220.0_wp
      met%par_diffuse = 60.0_wp
      met%air_density = air_density(met%air_temp, met%air_pressure, met%air_shv)
      site = met_at_site(met, f, 1, 500.0_wp, 1.0_wp)
      call check_true('500 m up is 500 m x the lapse rate colder', abs(site%air_temp - (290.0_wp - 3.25_wp)) < 1.0e-12_wp)
      call check_true('and at a lower pressure', site%air_pressure < met%air_pressure)
      call check_true('at the same relative humidity',                                               &
                      abs(relative_humidity(site) - relative_humidity(met)) < 1.0e-9_wp)
      call check_true('under less longwave', site%lw_down < met%lw_down)
      call check_true('with its air density from its own air',                                       &
                      site%air_density == air_density(site%air_temp, site%air_pressure, site%air_shv))
      call check_true('the light is the polygon''s on level ground',                                  &
                      site%par_beam == met%par_beam .and. site%par_diffuse == met%par_diffuse)
      site = met_at_site(met, f, 1, 0.0_wp, 1.0_wp)
      call check_true('a site at the polygon''s elevation on level ground keeps its forcing exactly',  &
                      site%air_temp == met%air_temp .and. site%air_pressure == met%air_pressure .and.                &
                      site%air_shv == met%air_shv .and. site%lw_down == met%lw_down .and.                        &
                      site%par_beam == met%par_beam .and. site%nir_beam == met%nir_beam .and.            &
                      site%air_density == met%air_density)
      f%apply_elevation_lapse = .false.
      site = met_at_site(met, f, 1, 500.0_wp, 1.5_wp)
      call check_true('without the lapse the air is the polygon''s',                                  &
                      site%air_temp == met%air_temp .and. site%air_pressure == met%air_pressure)
      call check_true('and the slope scales the beam alone',                                          &
                      site%par_beam == 1.5_wp * met%par_beam .and. site%nir_beam == 1.5_wp * met%nir_beam &
                      .and. site%par_diffuse == met%par_diffuse)

   end subroutine test_met_at_site

   !--- Relative humidity of a forcing sample: its vapour pressure over saturation.
   real(wp) function relative_humidity(m)
      type(met_forcing_t), intent(in) :: m
      relative_humidity = m%air_shv * m%air_pressure / (h2o_dry_ratio + (1.0_wp - h2o_dry_ratio) * m%air_shv)  &
                          / sat_vapor_pressure(m%air_temp)
   end function relative_humidity

end program test_met_file
