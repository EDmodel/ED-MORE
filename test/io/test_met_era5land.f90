! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_met_era5land -- the ED_ERA5land archive backend of the met reader                   !
! (docs/science/forcing.md §1, §2 and §9). A tiny synthetic archive is written through the !
! netCDF C API: a 4 x 18 grid (10 deg x 20 deg) spanning the antimeridian and two 16-cell  !
! chunk columns wide, every month of 2021 and January 2022, two no-data cells, and values  !
! that are analytic functions of cell and hour, so every check below has a known answer.   !
! It is read back through met_open / met_advance / met_instant and through the reader's    !
! own selection and month-load routines.                                                   !
!==========================================================================================!
program test_met_era5land
   use meds_test_assert,     only : check, check_true, test_report
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : h2o_dry_ratio
   use meds_time,            only : meds_time_t
   use meds_forcing_params,  only : forcing_params_t, MET_BACKEND_ERA5LAND, METAVG_END, CLAMP_ERROR
   use meds_radiation_lib,   only : SWPART_CLEARIDX, SWPART_PASSTHROUGH, clear_sky_emissivity
   use meds_met_types,   only : met_source_t, met_cursor_t, met_forcing_t, met_cells_t, met_month_t
   use meds_therm_lib,       only : sat_specific_humidity, lapse_pressure, sat_vapor_pressure
   use iso_c_binding,        only : c_double
   use meds_met,      only : met_open, met_cursor_init, met_advance, met_instant, met_close, met_prefetch
   use meds_met_share, only : MET_OK, MET_ERR_ARCHIVE, MET_ERR_ATTR_MISMATCH
   use meds_met_era5land,    only : era5land_path, era5land_default_template, era5land_load_month, ERA_OK,  &
                                    ERA_ERR_NAN, ERA_ERR_FILL, ERA_NVAR, ERA_TAIR, ERA_TDEW, ERA_PSURF,      &
                                    ERA_U10, ERA_V10, ERA_LWDOWN
   use meds_met_grid,        only : select_site_cell, select_box_cells, GRID_OK, GRID_ERR_NO_CELL
   use meds_test_archive, only : T0, field, hour_index, write_archive, write_month, archive_grid
   implicit none

   character(len=*), parameter :: DIR = 'era5land_test_tmp'

   call write_archive(DIR)

   call test_paths()
   call test_site_selection()
   call test_box_selection()
   call test_nan_rejected()
   call test_fill_rejected()
   call test_driver_months()
   call test_driver_recycle()
   call test_driver_rejections()
   call test_region_cursors()
   call test_driver_recycle_anchors()

   call test_report('test_met_era5land')

contains

   !----- What the driver does around a step: prefetch the day's forcing (the step's I/O), then   !
   !      advance and sample inside it with no file access. Returns the instantaneous forcing.  ----!
   function step_sample(src, cur, t) result(met)
      type(met_source_t), intent(inout) :: src
      type(met_cursor_t), intent(inout) :: cur
      type(meds_time_t),  intent(in)    :: t
      type(met_forcing_t) :: met
      integer(ik) :: loads
      call met_prefetch(src, meds_time_t(t%year, t%month, t%day))
      loads = src%n_loads
      call met_advance(src, cur, t) ; met = met_instant(src, cur, t)
      if (src%n_loads /= loads) error stop 'test_met_era5land: the reader loaded a month inside a step'
   end function step_sample

   !----- The site's cursor into an opened source: cell 1, at the config's location. ---------------!
   subroutine site_cursor(src, cur, fc)
      type(met_source_t),     intent(in)  :: src
      type(met_cursor_t),     intent(out) :: cur
      type(forcing_params_t), intent(in)  :: fc
      call met_cursor_init(src, cur, 1_ik, fc%latitude, fc%longitude,                   &
                           fc%elevation)
   end subroutine site_cursor

   !=======================================================================================!
   !  1. Templates.                                                                              !
   !=======================================================================================!
   subroutine test_paths()
      print '(a)', '-- era5land 1: path templates --'
      call check_true('default template',                                                        &
         trim(era5land_path(era5land_default_template(), '/data/x', 'Tair', 2021_ik, 3_ik))       &
         == '/data/x/ED_ERA5land_Tair_202103.nc')
      call check_true('custom template with a repeated token',                                   &
         trim(era5land_path('{data_path}/{yyyy}/{var}_{yyyy}{mm}.nc', '/d', 'u10', 2021_ik, 12_ik)) &
         == '/d/2021/u10_202112.nc')
   end subroutine test_paths

   !=======================================================================================!
   !  2. Site -> cell: its own cell when valid, else the nearest valid one within the limit.     !
   !=======================================================================================!
   subroutine test_site_selection()
      type(met_cells_t) :: dom
      real(wp)    :: d
      integer(ik) :: st
      print '(a)', '-- era5land 2: site to cell --'
      call select_site_cell(archive_grid(DIR//'/ED_ERA5land_static.nc'), 35.0_wp, -120.0_wp, 50.0_wp, dom, d, st)
      call check_true('a site on a valid cell takes that cell', st == GRID_OK .and. dom%ncell == 1_ik &
                      .and. dom%row(1) == 1_ik .and. dom%col(1) == 3_ik)
      call check('its distance is 0', d, 0.0_wp, 1.0e-6_wp)
      call check('its elevation comes from the static file', dom%elevation(1), 103.0_wp, 1.0e-6_wp)
      call check_true('a site reads a 1 x 1 box', dom%nchunk == 1_ik .and. dom%chunk_nrow(1) == 1_ik  &
                      .and. dom%chunk_ncol(1) == 1_ik .and. dom%chunk_row(1) == 1_ik .and. dom%chunk_col(1) == 3_ik)

      !----- (0,5) has no data. Its nearest valid neighbour is (1,5), 10 deg of latitude south   !
      !      (1112 km), not (0,4) or (0,6), 20 deg of longitude away at 45 N (about 1570 km). ------!
      call select_site_cell(archive_grid(DIR//'/ED_ERA5land_static.nc'), 45.0_wp, -80.0_wp, 2000.0_wp, dom, d, st)
      call check_true('a no-data site takes the nearest valid cell', st == GRID_OK                &
                      .and. dom%row(1) == 1_ik .and. dom%col(1) == 5_ik)
      call check('and reports its distance', d, 1111.95_wp, 0.1_wp)
      call select_site_cell(archive_grid(DIR//'/ED_ERA5land_static.nc'), 45.0_wp, -80.0_wp, 1000.0_wp, dom, d, st)
      call check_true('no valid cell within the limit is an error', st == GRID_ERR_NO_CELL)
   end subroutine test_site_selection

   !=======================================================================================!
   !  3. Box: valid cells in row-major order, across the antimeridian; and a month load of      !
   !     them, which reads one chunk column per touched chunk.                                   !
   !=======================================================================================!
   subroutine test_box_selection()
      type(met_cells_t) :: dom
      type(met_month_t)  :: buf
      character(len=512) :: msg
      real(wp)    :: across(4), inside(4)
      integer(ik) :: st, g
      print '(a)', '-- era5land 3: box selection and month load --'
      across = [40.0_wp, 150.0_wp, 20.0_wp, -150.0_wp] ; inside = [50.0_wp, -130.0_wp, 10.0_wp, -90.0_wp]
      !----- [N, W, S, E] = [40, 150, 20, -150]: rows 35 and 25 N; columns 160 E, then -180 and  !
      !      -160 across 180. (2,0) has no data, so 5 cells.                                        !
      call select_box_cells(archive_grid(DIR//'/ED_ERA5land_static.nc'), across, dom, st)
      call check_true('box across 180: 5 valid cells', st == GRID_OK .and. dom%ncell == 5_ik)
      call check_true('row-major, west to east across 180',                                      &
                      all(dom%row == [1_ik, 1_ik, 1_ik, 2_ik, 2_ik]) .and.                        &
                      all(dom%col == [17_ik, 0_ik, 1_ik, 17_ik, 1_ik]))
      call check_true('two chunk columns touched', dom%nchunk == 2_ik)
      !----- Each chunk reads only its cells' box: rows 1-2 x cols 0-1 west of 180, rows 1-2 x col 17 east. !
      call check_true('each chunk reads only the box of its cells',                                &
                      all(dom%chunk_row == [1_ik, 1_ik]) .and. all(dom%chunk_col == [0_ik, 17_ik]) .and. &
                      all(dom%chunk_nrow == [2_ik, 2_ik]) .and. all(dom%chunk_ncol == [2_ik, 1_ik]))
      call select_box_cells(archive_grid(DIR//'/ED_ERA5land_static.nc'), inside, dom, st)
      call check_true('box not across 180: 4 rows x 2 columns', st == GRID_OK .and. dom%ncell == 8_ik)
      !----- A region's selection rule: column 1 is 30% land, so land_fraction_min = 0.5 drops     !
      !      its two cells and keeps the order of the rest. -----------------------------------------!
      call select_box_cells(archive_grid(DIR//'/ED_ERA5land_static.nc'), across, dom, st, land_fraction_min=0.5_wp)
      call check_true('land_fraction_min drops the mostly-water cells',                            &
                      st == GRID_OK .and. dom%ncell == 3_ik .and.                                   &
                      all(dom%row == [1_ik, 1_ik, 2_ik]) .and. all(dom%col == [17_ik, 0_ik, 17_ik]))
      call check_true('the box keeps each cell''s land fraction', all(dom%land_fraction == 1.0_wp))

      call select_box_cells(archive_grid(DIR//'/ED_ERA5land_static.nc'), across, dom, st)
      call era5land_load_month(era5land_default_template(), DIR, dom, 2021_ik, 2_ik, buf, st, msg)
      call check_true('February loads', st == ERA_OK .and. buf%nt == 672_ik)
      g = hour_index(meds_time_t(2021_ik, 2_ik, 1_ik, 1_ik))
      call check('first hour, cell (2,17), Tair', real(buf%values(1, 4, ERA_TAIR), wp),          &
                 field(ERA_TAIR, 2_ik, 17_ik, g), 1.0e-9_wp)
      call check('last hour, cell (1,0), u10', real(buf%values(672, 2, ERA_U10), wp),            &
                 field(ERA_U10, 1_ik, 0_ik, g + 671_ik), 1.0e-9_wp)
      call check('mid-month, cell (1,1), PSurf', real(buf%values(300, 3, ERA_PSURF), wp),         &
                 field(ERA_PSURF, 1_ik, 1_ik, g + 299_ik), 1.0e-9_wp)
   end subroutine test_box_selection

   !=======================================================================================!
   !  4. MEDS never gap-fills (docs/science/forcing.md §10): a selected cell with no data  !
   !     is an error, never a fill. The second static file claims data at (0,5), where the !
   !     month files hold NaN.                                                             !
   !=======================================================================================!
   subroutine test_nan_rejected()
      type(met_cells_t) :: dom
      type(met_month_t)  :: buf
      character(len=512) :: msg
      real(wp)    :: d
      integer(ik) :: st
      print '(a)', '-- era5land 4: NaN in a selected cell --'
      call select_site_cell(archive_grid(DIR//'/static_allvalid.nc'), 45.0_wp, -80.0_wp, 50.0_wp, dom, d, st)
      call era5land_load_month(era5land_default_template(), DIR, dom, 2021_ik, 1_ik, buf, st, msg)
      call check_true('a NaN in a selected cell is rejected', st == ERA_ERR_NAN)
      call check_true('and the message names the variable', index(msg, 'Tair') > 0)
   end subroutine test_nan_rejected

   !----- 4b. A missing value is recognised as NaN only, so a file whose _FillValue is a number would  !
   !      pass that number off as data: it is refused. One month of every variable is                !
   !      written again, into a directory of its own, declaring _FillValue = -9999. -----------------!
   subroutine test_fill_rejected()
      character(len=*), parameter :: DIR2 = 'era5land_test_fill_tmp'
      type(met_cells_t) :: dom
      type(met_month_t)  :: buf
      character(len=512) :: msg
      real(wp)    :: d
      integer(ik) :: st, v
      print '(a)', '-- era5land 4b: a numeric _FillValue --'
      call execute_command_line('rm -rf '//DIR2//' && mkdir -p '//DIR2)
      do v = 1_ik, ERA_NVAR
         call write_month(DIR2, v, 2021_ik, 1_ik, fill=-9999.0_c_double)
      end do
      call select_site_cell(archive_grid(DIR//'/ED_ERA5land_static.nc'), 35.0_wp, -120.0_wp, 50.0_wp, dom, d, st)
      call era5land_load_month(era5land_default_template(), DIR2, dom, 2021_ik, 1_ik, buf, st, msg)
      call check_true('a month file with a numeric _FillValue is refused', st == ERA_ERR_FILL)
      call check_true('and the message says so', index(msg, '_FillValue') > 0)
      call era5land_load_month(era5land_default_template(), DIR, dom, 2021_ik, 1_ik, buf, st, msg)
      call check_true('the same month without the attribute loads', st == ERA_OK)
   end subroutine test_fill_rejected

   !----- A reader config on the synthetic archive at the valid cell (1,3), 35 N 120 W. ----------!
   function archive_config() result(fc)
      type(forcing_params_t) :: fc
      fc%format = MET_BACKEND_ERA5LAND ; fc%data_path = DIR ; fc%max_distance_km = 50.0_wp
      fc%timestep = 3600.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%latitude = 35.0_wp ; fc%longitude = -120.0_wp
      fc%recycle = .false. ; fc%start_clamp = CLAMP_ERROR
   end function archive_config

   !=======================================================================================!
   !  5. The driver over a month boundary: conversions, the seam bracket, the wind vector, and  !
   !     the static elevation.                                                                   !
   !=======================================================================================!
   subroutine test_driver_months()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, g
      real(wp)    :: u1, u2, v1, v2, t0, p0, q0, lw0, t1, p1, rh0, rh1
      integer     :: k
      print '(a)', '-- era5land 5: driver across a month boundary --'
      fc = archive_config()
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),               &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      call check_true('opens (MET_OK)', st == MET_OK)
      call site_cursor(src, cur, fc)
      call check_true('reads January and February: 1416 hourly records', src%nrec == 1416_ik)
      call check('grid_elevation from the static file', cur%grid_elevation_m, 103.0_wp, 1.0e-6_wp)
      call check_true('the archive supplies the wind vector', src%has_wind_vector)

      !----- At a stamp every weight is 0, so the record comes through untouched. ---------------!
      t = meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik) ; g = hour_index(t)
      met = step_sample(src, cur, t)
      call check('Tair at a stamp', met%air_temp, field(ERA_TAIR, 1_ik, 3_ik, g), 1.0e-9_wp)
      call check('Qair = the kernel of the stored Tdew and PSurf', met%air_shv,                      &
                 sat_specific_humidity(field(ERA_TDEW, 1_ik, 3_ik, g),                    &
                                               field(ERA_PSURF, 1_ik, 3_ik, g)), 1.0e-15_wp)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; v1 = field(ERA_V10, 1_ik, 3_ik, g)
      call check('wind_u = u10', met%wind_u, u1, 1.0e-12_wp)
      call check('wind_v = v10', met%wind_v, v1, 1.0e-12_wp)
      call check('speed from the components', met%wind, sqrt(u1**2 + v1**2), 1.0e-12_wp)
      call check_true('the instant carries the vector', met%has_wind_vector)

      !----- Inside January's file: 31 Jan 23:30 sits between 23:00 and 00:00 (both January). --!
      t = meds_time_t(2021_ik, 1_ik, 31_ik, 23_ik, 30_ik) ; g = hour_index(meds_time_t(2021_ik, 1_ik, 31_ik, 23_ik))
      met = step_sample(src, cur, t)
      call check('31 Jan 23:30 interpolates inside the January file', met%air_temp,                 &
                 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, g) + field(ERA_TAIR, 1_ik, 3_ik, g + 1_ik)), 1.0e-9_wp)

      !----- The seam: 1 Feb 00:00 is January's last record, 01:00 February's first. -----------!
      t = meds_time_t(2021_ik, 2_ik, 1_ik, 0_ik, 30_ik) ; g = hour_index(meds_time_t(2021_ik, 2_ik, 1_ik))
      met = step_sample(src, cur, t)
      call check('1 Feb 00:30 brackets January''s last and February''s first record', met%air_temp, &
                 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, g) + field(ERA_TAIR, 1_ik, 3_ik, g + 1_ik)), 1.0e-9_wp)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; u2 = field(ERA_U10, 1_ik, 3_ik, g + 1_ik)
      v1 = field(ERA_V10, 1_ik, 3_ik, g) ; v2 = field(ERA_V10, 1_ik, 3_ik, g + 1_ik)
      call check('the vector interpolates linearly (u)', met%wind_u, 0.5_wp * (u1 + u2), 1.0e-12_wp)
      call check('the vector interpolates linearly (v)', met%wind_v, 0.5_wp * (v1 + v2), 1.0e-12_wp)
      call check_true('the energy-form speed is at least the vector''s length',                    &
                      met%wind >= sqrt(met%wind_u**2 + met%wind_v**2) - 1.0e-12_wp)
      !----- January loaded at open, February at its first day's prefetch; the seam's 00:00          !
      !      record came from January's buffer, not from a second read. -----------------------------!
      call check_true('each month read once: January, then February', src%n_loads == 2_ik,          &
                      real(src%n_loads, wp))
      call met_close(src)

      !----- The TERRAIN lapse (docs/science/forcing.md §8, eqs. (13)-(14)): a site 500 m above its  !
      !      cell, in January, with twelve distinct monthly rates. T by January's rate, P               !
      !      hydrostatically, q at constant relative humidity, file LW by eps*T^4. The reader leaves    !
      !      the wind alone: the move to the canopy-air top is per patch, in the fast loop. ------------!
      fc%apply_elevation_lapse = .true.
      fc%lapse_rate_tair = [(0.0040_wp + 0.0002_wp * real(k, wp), k = 1, 12)]   ! January: 0.0042 K/m
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),               &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      fc%elevation = src%cells%elevation(1) + 500.0_wp
      call site_cursor(src, cur, fc)
      t = meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik) ; g = hour_index(t)
      met = step_sample(src, cur, t)
      t0  = field(ERA_TAIR, 1_ik, 3_ik, g) ; p0 = field(ERA_PSURF, 1_ik, 3_ik, g)
      q0  = sat_specific_humidity(field(ERA_TDEW, 1_ik, 3_ik, g), p0)
      lw0 = field(ERA_LWDOWN, 1_ik, 3_ik, g)
      t1  = t0 - 0.0042_wp * 500.0_wp ; p1 = lapse_pressure(p0, t0, 500.0_wp, 0.0042_wp)
      call check('lapse: T by January''s rate', met%air_temp, t1, 1.0e-9_wp)
      call check('lapse: P hydrostatic', met%air_pressure, p1, 1.0e-6_wp)
      rh0 = q0 * p0 / (h2o_dry_ratio + (1.0_wp - h2o_dry_ratio) * q0) / sat_vapor_pressure(t0)
      rh1 = met%air_shv * met%air_pressure / (h2o_dry_ratio + (1.0_wp - h2o_dry_ratio) * met%air_shv)            &
            / sat_vapor_pressure(met%air_temp)
      call check('lapse: relative humidity held', rh1, rh0, 1.0e-12_wp)
      call check_true('lapse: the colder site holds less vapour', met%air_shv < q0, met%air_shv - q0)
      call check('lapse: file longwave scaled by eps*T^4', met%lw_down,                              &
                 lw0 * clear_sky_emissivity(fc%lw_clear_form, met%air_temp, met%air_shv, met%air_pressure)      &
                     * met%air_temp**4 / (clear_sky_emissivity(fc%lw_clear_form, t0, q0, p0) * t0**4), 1.0e-9_wp)
      call check_true('lapse: the colder, drier site gets less longwave', met%lw_down < lw0, met%lw_down - lw0)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; v1 = field(ERA_V10, 1_ik, 3_ik, g)
      call check('the reader leaves the wind as delivered (u)', met%wind_u, u1, 1.0e-12_wp)
      call check('the reader leaves the wind as delivered (v)', met%wind_v, v1, 1.0e-12_wp)
      call met_close(src)
   end subroutine test_driver_months

   !=======================================================================================!
   !  6. Recycling across the months of a declared window, including the year-end seam, which   !
   !     brackets December's last record (2022-01-01 00:00) and January's first.                 !
   !=======================================================================================!
   subroutine test_driver_recycle()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, g
      print '(a)', '-- era5land 6: recycling across months --'
      fc = archive_config()
      fc%recycle = .true.
      fc%recycle_start = T0 ; fc%recycle_end = meds_time_t(2022_ik, 1_ik, 1_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('a whole-year window opens without a run period', st == MET_OK)
      call site_cursor(src, cur, fc)
      call check_true('the window is 8760 records', src%nrec == 8760_ik)
      t = meds_time_t(2035_ik, 6_ik, 15_ik, 12_ik) ; g = hour_index(meds_time_t(2021_ik, 6_ik, 15_ik, 12_ik))
      met = step_sample(src, cur, t)
      call check('2035-06-15 12:00 reads 2021-06-15 12:00', met%air_temp, field(ERA_TAIR, 1_ik, 3_ik, g), 1.0e-9_wp)
      t = meds_time_t(2036_ik, 1_ik, 1_ik, 0_ik, 30_ik)
      met = step_sample(src, cur, t)
      call check('the year-end seam brackets December''s last and January''s first record',       &
                 met%air_temp, 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, 8759_ik) + field(ERA_TAIR, 1_ik, 3_ik, 0_ik)), &
                 1.0e-9_wp)
      !----- Walk the model from 2036-12-30 across the wrap day by day: the wrap step needs December's  !
      !      last record and January, and the reader reads January once more, not December again. ---!
      block
         integer(ik) :: loads
         t = meds_time_t(2036_ik, 12_ik, 30_ik, 12_ik) ; met = step_sample(src, cur, t)
         t = meds_time_t(2036_ik, 12_ik, 31_ik, 23_ik, 30_ik) ; met = step_sample(src, cur, t)
         loads = src%n_loads
         t = meds_time_t(2037_ik, 1_ik, 1_ik, 0_ik, 30_ik) ; met = step_sample(src, cur, t)
         call check('across the wrap: December''s last and January''s first record', met%air_temp,   &
                    0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, 8759_ik) + field(ERA_TAIR, 1_ik, 3_ik, 0_ik)), 1.0e-9_wp)
         call check_true('the wrap reads January once, December''s record comes from the buffer',  &
                         src%n_loads == loads + 1_ik, real(src%n_loads - loads, wp))
      end block
      call met_close(src)
   end subroutine test_driver_recycle

   !=======================================================================================!
   !  7. Rejected at open: a month the archive lacks, total shortwave read as components, and  !
   !     a site with no valid cell in reach.                                                     !
   !=======================================================================================!
   subroutine test_driver_rejections()
      type(met_source_t)     :: src
      type(forcing_params_t) :: fc
      integer(ik) :: st
      print '(a)', '-- era5land 7: rejections at open --'
      fc = archive_config()
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 11_ik, 1_ik),               &
                    run_end=meds_time_t(2022_ik, 3_ik, 1_ik))
      call check_true('a run past the archive''s last month is rejected', st == MET_ERR_ARCHIVE)
      fc%sw_partition = SWPART_PASSTHROUGH
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 3_ik, 1_ik),                &
                    run_end=meds_time_t(2021_ik, 4_ik, 1_ik))
      call check_true('passthrough against total shortwave is rejected', st == MET_ERR_ATTR_MISMATCH)
      fc = archive_config()
      fc%latitude = 45.0_wp ; fc%longitude = -80.0_wp
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 3_ik, 1_ik),                &
                    run_end=meds_time_t(2021_ik, 4_ik, 1_ik))
      call check_true('a site with no valid cell within max_distance_km is rejected', st == MET_ERR_ARCHIVE)
   end subroutine test_driver_rejections

   !=======================================================================================!
   !  8. One source, many cursors (docs/science/forcing.md §2): a region opens the archive       !
   !     once for its cells, and each polygon walks it with its own cursor. Each cursor must read !
   !     what a site run at its cell reads, bit for bit, and one prefetch must serve them all.    !
   !=======================================================================================!
   subroutine test_region_cursors()
      integer(ik), parameter :: PICK(2) = [2_ik, 5_ik], PR(2) = [1_ik, 2_ik], PC(2) = [0_ik, 1_ik]
      type(met_source_t)     :: region, site(2)
      type(met_cursor_t)     :: cur(2), cur_site(2)
      type(met_cells_t)      :: cells
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met, met_site
      type(meds_time_t)      :: t, times(3)
      integer(ik) :: st, k, j, loads
      logical     :: same
      print '(a)', '-- era5land 8: one region source, one cursor per polygon --'
      !----- The box across 180 from test 3: (1,17), (1,0), (1,1), (2,17), (2,1). Polygons at the   !
      !      second and fifth cells, (1,0) and (2,1); a site run at each cell for comparison. ------!
      call select_box_cells(archive_grid(DIR//'/ED_ERA5land_static.nc'), [40.0_wp, 150.0_wp, 20.0_wp, -150.0_wp], &
                               cells, st)
      fc = archive_config()
      call met_open(region, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),            &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik), cells=cells)
      call check_true('a region opens on its cells (MET_OK)', st == MET_OK .and. region%cells%ncell == 5_ik)
      do k = 1_ik, 2_ik
         call met_cursor_init(region, cur(k), PICK(k), cells%lat(PICK(k)), cells%lon(PICK(k)),    &
                              fc%elevation)
         fc = archive_config()
         fc%latitude = cells%lat(PICK(k)) ; fc%longitude = cells%lon(PICK(k))
         call met_open(site(k), fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),        &
                       run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
         call site_cursor(site(k), cur_site(k), fc)
      end do
      call check('the first cursor takes its cell''s orography', cur(1)%grid_elevation_m, 100.0_wp, 1.0e-6_wp)
      call check('the second cursor takes its cell''s orography', cur(2)%grid_elevation_m, 201.0_wp, 1.0e-6_wp)

      !----- A stamp, a mid-month interpolation and the February seam. ----------------------------!
      times = [meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik), meds_time_t(2021_ik, 1_ik, 31_ik, 23_ik, 30_ik), &
               meds_time_t(2021_ik, 2_ik, 1_ik, 0_ik, 30_ik)]
      same = .true.
      do j = 1_ik, size(times)
         t = times(j)
         call met_prefetch(region, meds_time_t(t%year, t%month, t%day))
         loads = region%n_loads
         do k = 1_ik, 2_ik
            call met_advance(region, cur(k), t) ; met = met_instant(region, cur(k), t)
            met_site = step_sample(site(k), cur_site(k), t)
            same = same .and. met%air_temp == met_site%air_temp .and. met%air_shv == met_site%air_shv       &
                 .and. met%air_pressure == met_site%air_pressure .and. met%precip == met_site%precip    &
                 .and. met%wind == met_site%wind                                                &
                 .and. met%wind_u == met_site%wind_u .and. met%wind_v == met_site%wind_v        &
                 .and. met%lw_down == met_site%lw_down .and. met%par_beam == met_site%par_beam    &
                 .and. met%par_diffuse == met_site%par_diffuse                                  &
                 .and. met%nir_beam == met_site%nir_beam                                        &
                 .and. met%nir_diffuse == met_site%nir_diffuse .and. met%cosz == met_site%cosz
            if (j == 1_ik) call check('at a stamp the cursor reads its own cell''s Tair', met%air_temp, &
                                      field(ERA_TAIR, PR(k), PC(k), hour_index(t)), 1.0e-9_wp)
         end do
         call check_true('one prefetch serves every cursor (no read inside the step)',           &
                         region%n_loads == loads)
      end do
      call check_true('each polygon reads what a site run at its cell reads, bit for bit', same)
      call check_true('the region read each month once', region%n_loads == 2_ik)
      call met_close(region) ; call met_close(site(1))

      !----- A source closed and opened again starts clean (the C API reuses a run's source): on the  !
      !      second cell, a fresh cursor's first bracket at 1 Feb 00:30 needs January's last record,  !
      !      which would be the first cell's if met_close left the carried record behind.  -------------!
      fc = archive_config()
      fc%latitude = cells%lat(PICK(2)) ; fc%longitude = cells%lon(PICK(2))
      call met_open(site(1), fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),             &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      call site_cursor(site(1), cur_site(1), fc)
      t = times(3)
      met = step_sample(site(1), cur_site(1), t)
      call check('a reopened source reads its new cell''s carried record', met%air_temp,               &
                 0.5_wp * (field(ERA_TAIR, PR(2), PC(2), hour_index(meds_time_t(2021_ik, 2_ik, 1_ik)))  &
                 + field(ERA_TAIR, PR(2), PC(2), hour_index(meds_time_t(2021_ik, 2_ik, 1_ik)) + 1_ik)), 1.0e-9_wp)
      call met_close(site(1)) ; call met_close(site(2))
   end subroutine test_region_cursors

   !=======================================================================================!
   !  9. Recycle windows starting at other record stamps. Only a window starting at 01:00 puts  !
   !     the seam at the start of a daily step. A midnight anchor puts it at the end of one, and  !
   !     any other hour inside one, whose remaining hours read the window's first day while the  !
   !     buffer holds the window's last month. Neither may load inside the step. The third        !
   !     window starts at 00:00 on the 1st, as a region's may: its first record is carried into   !
   !     the next month without reading the month before it again.                               !
   !=======================================================================================!
   subroutine test_driver_recycle_anchors()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_params_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, loads, d
      print '(a)', '-- era5land 9: recycle windows from other record stamps --'
      fc = archive_config()
      fc%recycle = .true.

      !----- 2021-01-02 00:00 .. 2022-01-02 00:00: the seam, 2022-01-01 23:00 + 1 h, ends a step. ---!
      fc%recycle_start = meds_time_t(2021_ik, 1_ik, 2_ik) ; fc%recycle_end = meds_time_t(2022_ik, 1_ik, 2_ik)
      call met_open(src, fc, stat=st)
      call check_true('a midnight anchor opens (MET_OK)', st == MET_OK)
      call site_cursor(src, cur, fc)
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 1_ik, 12_ik))
      call check('the day before a midnight anchor reads the window''s last day', met%air_temp,          &
                 field(ERA_TAIR, 1_ik, 3_ik, hour_index(meds_time_t(2022_ik, 1_ik, 1_ik, 12_ik))), 1.0e-9_wp)
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 1_ik, 23_ik, 30_ik))
      call check('a midnight anchor''s seam: the window''s last and first records', met%air_temp,         &
                 seam_mean(meds_time_t(2022_ik, 1_ik, 1_ik, 23_ik), fc%recycle_start), 1.0e-9_wp)
      loads = src%n_loads
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 2_ik, 0_ik, 30_ik))
      call check('after it, the window''s first hour', met%air_temp, pair_mean(fc%recycle_start), 1.0e-9_wp)
      call check_true('the next day loads the window''s first month, once', src%n_loads == loads + 1_ik,  &
                      real(src%n_loads - loads, wp))
      call met_close(src)

      !----- 2021-01-01 06:00 .. 2022-01-01 06:00: the seam, 2022-01-01 05:00 + 1 h, falls inside a day  !
      !      whose last 18 hours are the window's first day. -------------------------------------------!
      fc%recycle_start = meds_time_t(2021_ik, 1_ik, 1_ik, 6_ik)
      fc%recycle_end   = meds_time_t(2022_ik, 1_ik, 1_ik, 6_ik)
      call met_open(src, fc, stat=st)
      call check_true('a 06:00 anchor opens (MET_OK)', st == MET_OK)
      call site_cursor(src, cur, fc)
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 1_ik, 2_ik, 30_ik))
      call check('before a 06:00 anchor''s seam: the window''s last month', met%air_temp,                &
                 pair_mean(meds_time_t(2022_ik, 1_ik, 1_ik, 2_ik)), 1.0e-9_wp)
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 1_ik, 5_ik, 30_ik))
      call check('a 06:00 anchor''s seam, inside the step', met%air_temp,                                 &
                 seam_mean(meds_time_t(2022_ik, 1_ik, 1_ik, 5_ik), fc%recycle_start), 1.0e-9_wp)
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 1_ik, 12_ik, 30_ik))
      call check('the same step after the seam reads the window''s first day', met%air_temp,           &
                 pair_mean(meds_time_t(2021_ik, 1_ik, 1_ik, 12_ik)), 1.0e-9_wp)
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 1_ik, 23_ik, 30_ik))
      call check('... through its last hour', met%air_temp, pair_mean(meds_time_t(2021_ik, 1_ik, 1_ik, 23_ik)), &
                 1.0e-9_wp)
      loads = src%n_loads
      met = step_sample(src, cur, meds_time_t(2036_ik, 1_ik, 2_ik, 0_ik, 30_ik))
      call check('the next day reads the window''s first month', met%air_temp,                           &
                 pair_mean(meds_time_t(2021_ik, 1_ik, 2_ik)), 1.0e-9_wp)
      call check_true('... loaded once', src%n_loads == loads + 1_ik, real(src%n_loads - loads, wp))
      call met_close(src)

      !----- 2021-02-01 00:00 .. 2022-02-01 00:00, walked as a region walks a month: one prefetch at  !
      !      the month's start, and no later day of the month may load. The seam ends the month's last  !
      !      day; the next month takes the window's first record from the kept day. -------------------!
      fc%recycle_start = meds_time_t(2021_ik, 2_ik, 1_ik) ; fc%recycle_end = meds_time_t(2022_ik, 2_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('an anchor at 00:00 on the 1st opens (MET_OK)', st == MET_OK)
      call site_cursor(src, cur, fc)
      call met_prefetch(src, meds_time_t(2036_ik, 1_ik, 1_ik))
      loads = src%n_loads
      do d = 1_ik, 31_ik
         call met_prefetch(src, meds_time_t(2036_ik, 1_ik, d))
      end do
      t = meds_time_t(2036_ik, 1_ik, 31_ik, 23_ik, 30_ik)
      call met_advance(src, cur, t) ; met = met_instant(src, cur, t)
      call check('the seam of a window from the 1st ends the month', met%air_temp,                         &
                 seam_mean(meds_time_t(2022_ik, 1_ik, 31_ik, 23_ik), fc%recycle_start), 1.0e-9_wp)
      call check_true('the month loaded nothing after its first day', src%n_loads == loads)
      met = step_sample(src, cur, meds_time_t(2036_ik, 2_ik, 1_ik, 0_ik, 30_ik))
      call check('the next month opens on the window''s first hour', met%air_temp,                        &
                 pair_mean(fc%recycle_start), 1.0e-9_wp)
      call check_true('... loading its own month only', src%n_loads == loads + 1_ik,                     &
                      real(src%n_loads - loads, wp))
      call met_close(src)
   end subroutine test_driver_recycle_anchors

   !----- The site cell's Tair half an hour after stamp t: the mean of its record and the next. ----!
   real(wp) function pair_mean(t) result(x)
      type(meds_time_t), intent(in) :: t
      x = 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, hour_index(t)) + field(ERA_TAIR, 1_ik, 3_ik, hour_index(t) + 1_ik))
   end function pair_mean

   !----- Its Tair half an hour into the seam: the mean of the window's last record and its first. --!
   real(wp) function seam_mean(t_last, t_first) result(x)
      type(meds_time_t), intent(in) :: t_last, t_first
      x = 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, hour_index(t_last)) + field(ERA_TAIR, 1_ik, 3_ik, hour_index(t_first)))
   end function seam_mean

end program test_met_era5land
