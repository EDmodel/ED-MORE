! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_met_era5land -- format = "ERA5land": the ERA5-Land archive, monthly files of          !
! hourly fields over its global grid (docs/science/forcing.md §1-§2; the files' layout is  !
! MEDS_FORCING_DESIGN.md §14), laid end to end as one hourly axis over the months a run    !
! needs. era5land_open opens it, era5land_prefetch loads a month before the step that      !
! reads it, and era5land_record hands the shared ingest one record's values as stored.     !
! era5land_grid is its grid (ed_met_grid). A month load checks every file, and ED-MORE     !
! never gap-fills: a NaN in a selected cell is an error, never a fill.                     !
!==========================================================================================!
module ed_met_era5land
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use iso_c_binding,       only : c_int, c_size_t, c_double, c_float
   use ed_kinds,            only : wp, ik
   use ed_constants,        only : day_sec
   use ed_time,             only : ed_time_t, time_advance_seconds, seconds_between, seconds_into_day,    &
                                   time_to_string, days_in_month, time_units_base
   use ed_forcing_params,   only : MET_PATH_LEN, ARCHIVE_DT_SEC, forcing_params_t
   use ed_radiation_lib,    only : SWPART_PASSTHROUGH
   use ed_met_types,        only : met_record_t, met_source_t, met_cursor_t, met_grid_t, met_cells_t,       &
                                   met_month_t, HUMIDITY_TDEW
   use ed_met_grid,         only : select_site_cell, GRID_OK, GRID_ERR_OPEN, GRID_ERR_NO_CELL
   use ed_netcdf_c,         only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f, nc_get_att_text_f,            &
                                   nc_get_att_double_f, nc_get_vara_double, nc_get_vara_float, nc_close,     &
                                   NC_NOERR, NC_NOWRITE
   use ed_met_share,      only :   MET_OK, MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, REC_MATCH_TOL,             &
                                   validate_recycle_window, file_lookup_sec, first_record_after,             &
                                   validate_file_against_config, add_line
   implicit none
   private

   public :: era5land_open, era5land_prefetch, era5land_record
   public :: era5land_grid, era5land_path, era5land_default_template, era5land_month_hours, era5land_load_month
   public :: ERA_NVAR, ERA_TAIR, ERA_TDEW, ERA_PSURF, ERA_U10, ERA_V10, ERA_RAINF, ERA_SWDOWN, ERA_LWDOWN
   public :: ERA_VAR_NAME, ERA_VAR_UNITS, ERA_CHUNK
   public :: ERA_OK, ERA_ERR_OPEN, ERA_ERR_GRID, ERA_ERR_TIME, ERA_ERR_UNITS, ERA_ERR_NAN, ERA_ERR_FILL

   !----- The archive's eight variables, their file names and the units the builder writes. -----!
   integer(ik), parameter :: ERA_NVAR = 8_ik
   integer(ik), parameter :: ERA_TAIR = 1_ik, ERA_TDEW = 2_ik, ERA_PSURF = 3_ik, ERA_U10 = 4_ik,     &
                             ERA_V10 = 5_ik, ERA_RAINF = 6_ik, ERA_SWDOWN = 7_ik, ERA_LWDOWN = 8_ik
   character(len=6),  parameter :: ERA_VAR_NAME(ERA_NVAR) =                                        &
      ['Tair  ', 'Tdew  ', 'PSurf ', 'u10   ', 'v10   ', 'Rainf ', 'SWdown', 'LWdown']
   character(len=10), parameter :: ERA_VAR_UNITS(ERA_NVAR) =                                       &
      ['K         ', 'K         ', 'Pa        ', 'm s-1     ', 'm s-1     ', 'kg m-2 s-1',           &
       'W m-2     ', 'W m-2     ']
   !----- The archive's time coordinate is seconds since this instant. ---------------------------!
   type(ed_time_t), parameter :: ERA_EPOCH = ed_time_t(1970_ik, 1_ik, 1_ik, 0_ik, 0_ik, 0_ik)
   integer(ik), parameter :: ERA_CHUNK = 16_ik     !< [cells] the files' chunk edge
   real(wp),    parameter :: STAMP_TOL = 0.5_wp    !< [s] slack for "this record is that stamp"

   !----- Status codes. --------------------------------------------------------------------------!
   integer(ik), parameter :: ERA_OK          = 0_ik
   integer(ik), parameter :: ERA_ERR_OPEN    = 1_ik   !< a file is missing or unreadable
   integer(ik), parameter :: ERA_ERR_GRID    = 2_ik   !< a file's grid differs from the static file's
   integer(ik), parameter :: ERA_ERR_TIME    = 3_ik   !< a month file's stamps are not that month's hours
   integer(ik), parameter :: ERA_ERR_UNITS   = 4_ik   !< a variable's units attribute is not the expected one
   integer(ik), parameter :: ERA_ERR_NAN     = 5_ik   !< a selected cell has a missing value
   integer(ik), parameter :: ERA_ERR_FILL    = 6_ik   !< a variable's _FillValue is a number, not NaN

contains

   !=======================================================================================!
   !  OPEN the ED_ERA5land archive: its cells and month axis (open_archive), the           !
   !  declared recycle window checked against that axis, and the axis's first month loaded,  !
   !  with the window's first day kept for the recycle seam. A rejection returns its status   !
   !  and the reason; the caller closes the source.                                        !
   !=======================================================================================!
   subroutine era5land_open(src, stat, why, run_start, run_end, cells)
      type(met_source_t),            intent(inout) :: src
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(out)   :: why    !< what stops the run when stat /= MET_OK
      type(ed_time_t),   optional,   intent(in)    :: run_start, run_end
      type(met_cells_t), optional,   intent(in)    :: cells  !< a region's cells; absent: the site's
      why = ''
      call open_archive(src, run_start, run_end, stat, why, cells)
      if (stat == MET_OK) call validate_recycle_window(src, stat, why)
      if (stat /= MET_OK) return
      call load_axis_month(src, 1_ik, stat, why)
      if (stat /= MET_OK) return
      call keep_window_head(src)
      src%rec_first = 1_ik
   end subroutine era5land_open

   !=======================================================================================!
   !  PREFETCH from the archive: the month the step starting at `step_start` reads, and the    !
   !  record before it (00:00 on the 1st, or the window's last record at the recycle wrap),     !
   !  which comes from the outgoing month, so a run loads each month once                  !
   !  (docs/science/forcing.md §2). It assumes the daily steps from midnight validate_config   !
   !  requires of the archive.                                                             !
   !=======================================================================================!
   subroutine era5land_prefetch(src, step_start)
      type(met_source_t), intent(inout)  :: src
      type(ed_time_t),    intent(in)    :: step_start
      real(wp)    :: s0
      integer(ik) :: r, k, kp, p
      logical     :: wrap
      s0 = file_lookup_sec(src, step_start)
      wrap = .false.
      if (src%n_cycle_years >= 1_ik .and. src%irec_cycle_last > src%irec_cycle_first) then
         wrap = s0 >= src%time_sec(src%irec_cycle_last) - REC_MATCH_TOL
      end if
      if (wrap) then                                    ! the window's first month, after its last record
         k = axis_month_of(src, src%irec_cycle_first)
         p = src%irec_cycle_last
      else
         r = first_record_after(src, s0)                ! clamped to the axis
         k = axis_month_of(src, r)
         p = src%month_rec0(k)
      end if
      if (p > 0_ik .and. src%carry_rec /= p) then       ! the record before the month
         if (.not. allocated(src%carry)) allocate(src%carry(src%cells%ncell, ERA_NVAR))
         kp = axis_month_of(src, p)
         if (.not. month_loaded(src, kp) .and. in_window_head(src, p)) then
            src%carry = src%head(p - src%irec_cycle_first + 1_ik, :, :)   ! a midnight anchor's record
         else
            if (.not. month_loaded(src, kp)) call load_month_or_stop(src, kp)
            src%carry = src%buffer%values(p - src%month_rec0(kp), :, :)
         end if
         src%carry_rec = p
      end if
      if (.not. month_loaded(src, k)) call load_month_or_stop(src, k)
   end subroutine era5land_prefetch

   !----- One record's values at the cursor's cell, as the archive stores them: the dewpoint    !
   !      (converted in the shared ingest), the wind components with their speed, and the total  !
   !      shortwave, never the four streams (open_archive refuses "passthrough"). ---------------!
   subroutine era5land_record(src, cur, irec, rec, humidity_value, sw_total)
      type(met_source_t), intent(in)    :: src
      type(met_cursor_t), intent(in)    :: cur
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(inout) :: rec
      real(wp),           intent(out)   :: humidity_value, sw_total
      integer(ik) :: h
      call locate_record(src, irec, h)
      rec%air_temp   = archive_value(src, cur, h, ERA_TAIR)
      rec%air_pressure = archive_value(src, cur, h, ERA_PSURF)
      humidity_value = archive_value(src, cur, h, ERA_TDEW)
      rec%wind_u   = archive_value(src, cur, h, ERA_U10)
      rec%wind_v   = archive_value(src, cur, h, ERA_V10)
      rec%wind     = sqrt(rec%wind_u**2 + rec%wind_v**2)
      rec%rain    = archive_value(src, cur, h, ERA_RAINF)            ! total rainfall rate [kg/m2/s]
      rec%lw_down   = archive_value(src, cur, h, ERA_LWDOWN)
      sw_total     = archive_value(src, cur, h, ERA_SWDOWN)
   end subroutine era5land_record

   !=======================================================================================!
   !  ED_ERA5land OPEN: resolve the paths, take a region's cells or bind the site to its cell     !
   !  (whose orography is the grid elevation the terrain lapse starts from), and lay the months   !
   !  the run needs end to end as one hourly axis. Every month file is checked for existence      !
   !  here, so a gap in the archive stops the run before it starts rather than decades into it.   !
   !=======================================================================================!
   subroutine open_archive(src, run_start, run_end, stat, why, cells)
      type(met_source_t),            intent(inout) :: src
      type(ed_time_t), optional,     intent(in)    :: run_start, run_end
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(inout) :: why
      type(met_cells_t), optional,   intent(in)    :: cells
      character(len=MET_PATH_LEN) :: path
      character(len=MET_PATH_LEN + 160) :: buf
      real(wp)       :: distance_km, first
      integer(ik)    :: est, y0, m0, y1, m1, nmonth, k, v, h, y, m, nmiss
      integer(c_int) :: st, ncid
      logical        :: exists

      stat = MET_OK
      associate (f => src%fcfg)
      src%file_template = f%file_template
      if (len_trim(src%file_template) == 0) src%file_template = era5land_default_template()
      src%grid = era5land_grid(f)

      if (f%sw_partition == SWPART_PASSTHROUGH) then
         call add_line(why, ' met_open: the ED_ERA5land archive carries TOTAL shortwave, but')
         call add_line(why, '   [forcing].sw_partition = "passthrough" expects the four component streams.')
         stat = MET_ERR_ATTR_MISMATCH
         return
      end if

      !----- The cells: a region's, selected by its caller, or the site's own. -------------------!
      est = GRID_OK
      if (present(cells)) then
         src%cells = cells
      else
         call select_site_cell(src%grid, f%latitude, f%longitude, f%max_distance_km, src%cells, distance_km, est)
      end if
      if (est /= GRID_OK) then
         select case (est)
         case (GRID_ERR_NO_CELL)
            write(buf,'(a,f0.1,a,f0.4,a,f0.4,a)') ' met_open: no valid ERA5-Land cell within ',        &
                  f%max_distance_km, ' km of the site (', f%latitude, ', ', f%longitude, ')'
            call add_line(why, trim(buf))
         case (GRID_ERR_OPEN)
            call add_line(why, ' met_open: cannot open the ED_ERA5land static file '//trim(src%grid%static_file))
         case default
            call add_line(why, ' met_open: the static file '//trim(src%grid%static_file)//                  &
                            ' does not hold a regular lat/lon grid with valid and elevation')
         end select
         stat = MET_ERR_ARCHIVE
         return
      end if
      src%grid_index = 1_ik
      src%ngrid = src%cells%ncell
      src%has_wind_vector = .true.
      src%humidity = HUMIDITY_TDEW                           ! the archive stores the 2 m dewpoint
      !----- A caller that chose the cells reports them; the site's own choice is reported here. ---!
      if (.not. present(cells))                                                                 &
         write(*,'(a,f0.3,a,f0.3,a,f7.2,a,f0.1,a)') ' force : ED_ERA5land cell (', src%cells%lat(1), ', ',  &
               src%cells%lon(1), '), ', distance_km, ' km from the site; orography ', src%cells%elevation(1), ' m'

      !----- The months: the declared recycle window, or the run period. ----------------------!
      if (f%recycle) then
         call archive_month_of(f%recycle_start, .false., y0, m0)
         call archive_month_of(time_advance_seconds(f%recycle_end, -src%dt_forcing), .false., y1, m1)
      else
         if (.not. (present(run_start) .and. present(run_end))) then
            call add_line(why, ' met_open: format = "ERA5land" without recycling needs the run period')
            stat = MET_ERR_ARCHIVE
            return
         end if
         call archive_month_of(run_start, .false., y0, m0)
         call archive_month_of(run_end,   .true.,  y1, m1)
      end if
      nmonth = (y1 - y0) * 12_ik + (m1 - m0) + 1_ik
      if (nmonth < 1_ik) then
         call add_line(why, ' met_open: the ED_ERA5land month range is empty')
         stat = MET_ERR_ARCHIVE
         return
      end if
      allocate(src%month_year(nmonth), src%month_month(nmonth), src%month_rec0(nmonth))
      y = y0
      m = m0
      src%nrec = 0_ik
      do k = 1_ik, nmonth
         src%month_year(k) = y
         src%month_month(k) = m
         src%month_rec0(k) = src%nrec
         src%nrec = src%nrec + era5land_month_hours(y, m)
         m = m + 1_ik
         if (m > 12_ik) then
            m = 1_ik
            y = y + 1_ik
         end if
      end do
      src%base_time = ERA_EPOCH
      if (allocated(src%time_sec)) deallocate(src%time_sec)
      allocate(src%time_sec(src%nrec))
      do k = 1_ik, nmonth
         first = seconds_between(ERA_EPOCH, ed_time_t(src%month_year(k), src%month_month(k), 1_ik, 1_ik))
         do h = 1_ik, era5land_month_hours(src%month_year(k), src%month_month(k))
            src%time_sec(src%month_rec0(k) + h) = first + ARCHIVE_DT_SEC * real(h - 1_ik, wp)
         end do
      end do

      nmiss = 0_ik
      do k = 1_ik, nmonth
         do v = 1_ik, ERA_NVAR
            path = era5land_path(src%file_template, f%data_path, ERA_VAR_NAME(v),                  &
                                 src%month_year(k), src%month_month(k))
            inquire(file=trim(path), exist=exists)
            if (.not. exists) then
               nmiss = nmiss + 1_ik
               if (nmiss == 1_ik) call add_line(why, ' met_open: ED_ERA5land file missing: '//trim(path))
            end if
         end do
      end do
      if (nmiss > 0_ik) then
         write(buf,'(a,i0,a,i0,a)') '   ', nmiss, ' of ', nmonth * ERA_NVAR, ' files the run needs are missing.'
         call add_line(why, trim(buf))
         stat = MET_ERR_ARCHIVE
         return
      end if

      !----- The first month file describes the archive; check it against the config. ----------!
      path = era5land_path(src%file_template, f%data_path, 'Tair', y0, m0)
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) then
         call add_line(why, ' met_open: cannot open '//trim(path))
         stat = MET_ERR_ARCHIVE
         return
      end if
      call validate_file_against_config(src, ncid, stat, why)
      st = nc_close(ncid)
      src%buffer%year = 0_ik
      src%carry_rec = 0_ik
      src%n_loads = 0_ik
      end associate
   end subroutine open_archive

   !----- The archive month whose file holds the record at t, taken to the hour below (or above,   !
   !      round_up) on the hourly grid. Files are end-stamped, so 00:00 on the 1st belongs to the   !
   !      month before. -----------------------------------------------------------------------------!
   subroutine archive_month_of(t, round_up, year, month)
      type(ed_time_t),   intent(in)  :: t
      logical,           intent(in)  :: round_up
      integer(ik),       intent(out) :: year, month
      type(ed_time_t) :: before
      real(wp) :: hours
      hours = seconds_between(ERA_EPOCH, t) / ARCHIVE_DT_SEC
      if (round_up) then
         hours = real(ceiling(hours), wp)
      else
         hours = real(floor(hours), wp)
      end if
      before = time_advance_seconds(ERA_EPOCH, hours * ARCHIVE_DT_SEC - 1.0_wp)
      year = before%year
      month = before%month
   end subroutine archive_month_of

   !----- The axis month holding record irec (binary search: an axis can span decades). ----------!
   pure integer(ik) function axis_month_of(src, irec) result(k)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: irec
      integer(ik) :: lo, hi, mid
      lo = 1_ik
      hi = size(src%month_rec0, kind=ik)          ! month_rec0(lo) < irec always holds
      do while (lo < hi)
         mid = (lo + hi + 1_ik) / 2_ik
         if (src%month_rec0(mid) < irec) then
            lo = mid
         else
            hi = mid - 1_ik
         end if
      end do
      k = lo
   end function axis_month_of

   pure logical function month_loaded(src, k) result(yes)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: k
      yes = src%buffer%year == src%month_year(k) .and. src%buffer%month == src%month_month(k)
   end function month_loaded

   pure logical function in_window_head(src, irec) result(yes)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: irec
      yes = irec >= src%irec_cycle_first .and. irec < src%irec_cycle_first + src%n_head
   end function in_window_head

   !----- Keep the recycle window's first day, from axis month 1 just loaded: its first record    !
   !      through the next midnight, all in that month's file. The seam falls inside a daily step  !
   !      unless the window starts at 01:00, and that step reads the window's last record, in the  !
   !      axis's last month, and then these. A midnight anchor needs its one record.               !
   subroutine keep_window_head(src)
      type(met_source_t), intent(inout) :: src
      real(wp)    :: time_of_day
      integer(ik) :: h1
      src%n_head = 0_ik
      if (allocated(src%head)) deallocate(src%head)
      if (src%n_cycle_years < 1_ik) return
      time_of_day = seconds_into_day(src%cycle_anchor)
      src%n_head = 1_ik
      if (time_of_day > REC_MATCH_TOL)                                                            &
         src%n_head = 1_ik + nint((day_sec - time_of_day) / src%dt_forcing, ik)
      h1 = src%irec_cycle_first - src%month_rec0(1)
      allocate(src%head(src%n_head, src%cells%ncell, ERA_NVAR))
      src%head = src%buffer%values(h1:h1 + src%n_head - 1_ik, :, :)
   end subroutine keep_window_head

   !----- Read axis month k into the buffer: the only place the archive is read after open. A      !
   !      month that fails its checks returns MET_ERR_ARCHIVE and the reason. -----------------------!
   subroutine load_axis_month(src, k, stat, why)
      type(met_source_t),            intent(inout) :: src
      integer(ik),                   intent(in)    :: k
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(inout) :: why
      character(len=MET_PATH_LEN + 128) :: message
      integer(ik) :: est
      stat = MET_OK
      call era5land_load_month(src%file_template, src%fcfg%data_path, src%cells, src%month_year(k), &
                               src%month_month(k), src%buffer, est, message)
      if (est /= ERA_OK) then
         call add_line(why, ' met_driver: ED_ERA5land month failed its checks: '//trim(message))
         stat = MET_ERR_ARCHIVE
         return
      end if
      src%n_loads = src%n_loads + 1_ik
   end subroutine load_axis_month

   !----- A month a step needs, in the I/O phase: one that fails its checks stops the run, since   !
   !      ED-MORE never gap-fills (docs/science/forcing.md §10). ------------------------------------!
   subroutine load_month_or_stop(src, k)
      type(met_source_t), intent(inout) :: src
      integer(ik),        intent(in)    :: k
      character(len=:), allocatable :: why
      integer(ik) :: stat
      why = ''
      call load_axis_month(src, k, stat, why)
      if (stat /= MET_OK) then
         write(*,'(a)') why
         error stop 'met_driver: an ED_ERA5land month cannot be read (ED-MORE does not gap-fill)'
      end if
   end subroutine load_month_or_stop

   !----- Where record irec's values are: its hour h > 0 in the loaded month, h = 0 for the carried !
   !      record, h < 0 for the -h-th record of the window's first day. Anything else was not       !
   !      prefetched, which is a programming error, not a data gap. -------------------------------!
   subroutine locate_record(src, irec, h)
      type(met_source_t), intent(in)   :: src
      integer(ik),        intent(in)  :: irec
      integer(ik),        intent(out) :: h
      integer(ik) :: k
      k = axis_month_of(src, irec)
      if (month_loaded(src, k)) then
         h = irec - src%month_rec0(k)
      else if (irec == src%carry_rec) then
         h = 0_ik
      else if (in_window_head(src, irec)) then
         h = src%irec_cycle_first - 1_ik - irec
      else
         write(*,'(a,i0,2a)') ' met_driver: forcing record ', irec, ' at ',                            &
               time_to_string(time_advance_seconds(src%base_time, src%time_sec(irec)))
         error stop 'met_driver: a forcing record was not prefetched before the step (internal error)'
      end if
   end subroutine locate_record

   pure real(wp) function archive_value(src, cur, h, var) result(val)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      integer(ik),        intent(in) :: h, var
      if (h == 0_ik) then
         val = real(src%carry(cur%cell, var), wp)
      else if (h < 0_ik) then
         val = real(src%head(-h, cur%cell, var), wp)
      else
         val = real(src%buffer%values(h, cur%cell, var), wp)
      end if
   end function archive_value

   !=======================================================================================!
   !  Paths. Tokens {data_path}, {var}, {yyyy} and {mm}; the defaults are the archive's own    !
   !  flat layout.                                                                         !
   !=======================================================================================!
   function era5land_path(template, data_path, var, year, month) result(path)
      character(len=*), intent(in) :: template, data_path, var
      integer(ik),      intent(in) :: year, month
      character(len=MET_PATH_LEN)  :: path
      character(len=4) :: yyyy
      character(len=2) :: mm
      character(len=:), allocatable :: s
      write(yyyy, '(i4.4)') year
      write(mm, '(i2.2)') month
      s = replace_all(trim(template), '{data_path}', trim(data_path))
      s = replace_all(s, '{var}', trim(var))
      s = replace_all(s, '{yyyy}', yyyy)
      s = replace_all(s, '{mm}', mm)
      path = s
   end function era5land_path

   pure function era5land_default_template() result(t)
      character(len=MET_PATH_LEN) :: t
      t = '{data_path}/ED_ERA5land_{var}_{yyyy}{mm}.nc'
   end function era5land_default_template

   pure function era5land_default_static() result(t)
      character(len=MET_PATH_LEN) :: t
      t = '{data_path}/ED_ERA5land_static.nc'
   end function era5land_default_static

   !----- The archive's grid: the static file a [forcing] block names (static_file, or the default, -!
   !      in the archive's folder) and the files' 16 x 16 chunks. ---------------------------------!
   function era5land_grid(f) result(grid)
      type(forcing_params_t), intent(in) :: f
      type(met_grid_t) :: grid
      grid%static_file = f%static_file
      if (len_trim(grid%static_file) == 0) grid%static_file = era5land_default_static()
      grid%static_file = era5land_path(grid%static_file, f%data_path, '', 0_ik, 0_ik)
      grid%chunk = ERA_CHUNK
   end function era5land_grid

   !----- Every occurrence of `token` in `s` replaced by `value`. ----------------------------------!
   pure function replace_all(s, token, value) result(out)
      character(len=*), intent(in)  :: s, token, value
      character(len=:), allocatable :: out
      integer :: k, pos
      out = ''
      pos = 1
      do
         k = index(s(pos:), token)
         if (k == 0) exit
         out = out // s(pos:pos + k - 2) // value
         pos = pos + k - 1 + len(token)
      end do
      out = out // s(pos:)
   end function replace_all

   !----- Records in an archive month file: 01:00 on the 1st .. 00:00 on the 1st of the next month. !
   elemental integer(ik) function era5land_month_hours(year, month) result(nt)
      integer(ik), intent(in) :: year, month
      nt = nint(day_sec / ARCHIVE_DT_SEC, ik) * days_in_month(year, month)
   end function era5land_month_hours

   !=======================================================================================!
   !  MONTH LOAD (docs/science/forcing.md §1): every variable of one archive month for every    !
   !  domain cell. Each file's grid, stamps and units are checked against what the reader  !
   !  assumes; each touched chunk column is one hyperslab read of (month hours, <=16, <=16).     !
   !  `message` names the offending file, cell or hour on a nonzero status.                !
   !=======================================================================================!
   subroutine era5land_load_month(template, data_path, dom, year, month, buf, stat, message)
      character(len=*),   intent(in)    :: template, data_path
      type(met_cells_t),  intent(in)    :: dom
      integer(ik),        intent(in)    :: year, month
      type(met_month_t),  intent(inout) :: buf
      integer(ik),        intent(out)   :: stat
      character(len=*),   intent(out)   :: message
      character(len=MET_PATH_LEN) :: path
      character(len=64)  :: units
      integer(c_int)     :: ncid, vid, st
      integer(c_size_t)  :: start3(3), count3(3)
      real(c_float), allocatable :: block(:,:,:)
      real(c_double)     :: fill
      integer(ik) :: nt, v, k, m, c, nr, nc
      integer(ik) :: bad(3)

      stat = ERA_OK
      message = ''
      nt = era5land_month_hours(year, month)
      buf%year = 0_ik
      if (allocated(buf%values)) then
         if (any(shape(buf%values) /= [nt, dom%ncell, ERA_NVAR])) deallocate(buf%values)
      end if
      if (.not. allocated(buf%values)) allocate(buf%values(nt, dom%ncell, ERA_NVAR))

      do v = 1_ik, ERA_NVAR
         path = era5land_path(template, data_path, ERA_VAR_NAME(v), year, month)
         st = nc_open_f(trim(path), NC_NOWRITE, ncid)
         if (st /= NC_NOERR) then
            stat = ERA_ERR_OPEN
            message = trim(path)
            return
         end if
         call check_month_file(ncid, dom, year, month, nt, stat)
         if (stat /= ERA_OK) then
            st = nc_close(ncid)
            message = trim(path)
            return
         end if
         st = nc_inq_varid_f(ncid, trim(ERA_VAR_NAME(v)), vid)
         if (st == NC_NOERR) st = nc_get_att_text_f(ncid, vid, 'units', units)
         if (st /= NC_NOERR .or. trim(units) /= trim(ERA_VAR_UNITS(v))) then
            st = nc_close(ncid)
            stat = ERA_ERR_UNITS
            message = trim(path)//': units "'//trim(units)//'", expected "'//trim(ERA_VAR_UNITS(v))//'"'
            return
         end if
         !----- The reader recognises a missing value as NaN and nothing else (below). A file whose   !
         !      _FillValue is a number would slip that number into the forcing as data, so it is     !
         !      refused. The archive stores NaN, and a file without the attribute has no fill.  -----!
         st = nc_get_att_double_f(ncid, vid, '_FillValue', fill)
         if (st == NC_NOERR) then
            if (.not. ieee_is_nan(fill)) then
               st = nc_close(ncid)
               stat = ERA_ERR_FILL
               write(message, '(2a,es12.4,a)') trim(path), ': _FillValue is ', fill, ', not NaN'
               return
            end if
         end if
         do k = 1_ik, dom%nchunk
            nr = dom%chunk_nrow(k)
            nc = dom%chunk_ncol(k)
            allocate(block(nc, nr, nt))                      ! C order [time][lat][lon]
            start3 = [0_c_size_t, int(dom%chunk_row(k), c_size_t), int(dom%chunk_col(k), c_size_t)]
            count3 = [int(nt, c_size_t), int(nr, c_size_t), int(nc, c_size_t)]
            st = nc_get_vara_float(ncid, vid, start3, count3, block)   ! the archive stores float32
            if (st /= NC_NOERR) then
               deallocate(block)
               st = nc_close(ncid)
               stat = ERA_ERR_OPEN
               message = trim(path)//': read failed'
               return
            end if
            do m = dom%chunk_first(k), dom%chunk_first(k + 1_ik) - 1_ik
               c = dom%by_chunk(m)
               buf%values(:, c, v) = block(dom%col(c) - dom%chunk_col(k) + 1_ik,                     &
                                           dom%row(c) - dom%chunk_row(k) + 1_ik, :)
            end do
            deallocate(block)
         end do
         st = nc_close(ncid)
      end do

      !----- ED-MORE never gap-fills (docs/science/forcing.md §10): the static mask promises data  !
      !      in every selected cell. ieee_is_nan, because a -fpe0 build can trap x /= x on a NaN,    !
      !      and findloc over an integer mask, not a logical one (a compiler trap, CLAUDE.md). ------!
      if (any(ieee_is_nan(buf%values))) then
         bad = int(findloc(merge(1, 0, ieee_is_nan(buf%values)), 1), ik)
         write(message, '(a,a,i4.4,a,i2.2,a,i0,a,f9.3,a,f9.3,a)') trim(ERA_VAR_NAME(bad(3))), ' ',         &
               year, '-', month, ' hour ', bad(1), ' at cell (', dom%lat(bad(2)), ', ', dom%lon(bad(2)), ')'
         stat = ERA_ERR_NAN
         return
      end if
      buf%year = year
      buf%month = month
      buf%nt = nt
   end subroutine era5land_load_month

   !----- A month file must have the static grid and exactly the month's end-stamped hours. -------!
   subroutine check_month_file(ncid, dom, year, month, nt, stat)
      integer(c_int),     intent(in)  :: ncid
      type(met_cells_t),  intent(in)  :: dom
      integer(ik),        intent(in)  :: year, month, nt
      integer(ik),        intent(out) :: stat
      integer(c_int)     :: st, vid
      integer(c_size_t)  :: n, start1(1), count1(1)
      character(len=256) :: units
      type(ed_time_t)    :: base
      real(c_double), allocatable :: tsec(:)
      real(wp)           :: first
      logical            :: ok
      integer(ik)        :: k

      stat = ERA_ERR_GRID
      st = nc_inq_dimlen_f(ncid, 'lat', n)
      if (st /= NC_NOERR .or. int(n, ik) /= dom%nlat) return
      st = nc_inq_dimlen_f(ncid, 'lon', n)
      if (st /= NC_NOERR .or. int(n, ik) /= dom%nlon) return
      stat = ERA_ERR_TIME
      st = nc_inq_dimlen_f(ncid, 'time', n)
      if (st /= NC_NOERR .or. int(n, ik) /= nt) return
      st = nc_inq_varid_f(ncid, 'time', vid)
      if (st /= NC_NOERR) return
      st = nc_get_att_text_f(ncid, vid, 'units', units)
      if (st /= NC_NOERR) return
      call time_units_base(units, base, ok)
      if (.not. ok) return
      allocate(tsec(nt))
      start1 = 0_c_size_t
      count1 = int(nt, c_size_t)
      st = nc_get_vara_double(ncid, vid, start1, count1, tsec)
      if (st /= NC_NOERR) return
      first = seconds_between(base, ed_time_t(year, month, 1_ik, 1_ik, 0_ik, 0_ik))
      do k = 1_ik, nt
         if (abs(tsec(k) - (first + ARCHIVE_DT_SEC * real(k - 1_ik, wp))) > STAMP_TOL) return
      end do
      stat = ERA_OK
   end subroutine check_month_file

end module ed_met_era5land
