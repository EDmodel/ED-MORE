! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met_share -- what the two forcing file sources share (docs/science/forcing.md   !
! §1, §2 and §9): the reader's status codes, the time axis with its declared recycle window     !
! (validate_recycle_window, recycle_model_to_file), and the check of a file against the         !
! [forcing] config that describes it (validate_file_against_config). A MEDS forcing file and    !
! the ED_ERA5land archive both lay their records on one axis of seconds since a base time, so    !
! meds_met brackets and recycles both with one code. A check that fails returns its status and  !
! its reason, line by line (add_line); met_open prints it.                                       !
!==========================================================================================!
module meds_met_share
   use iso_c_binding, only : c_int, c_double
   use meds_kinds, only : wp, ik
   use meds_time, only : meds_time_t, time_advance_seconds, seconds_between, is_leap_year,        &
                         time_to_string, whole_years_between
   use meds_forcing_params, only : MET_BACKEND_DEFAULT, METAVG_BEGIN, HEIGHT_ABOVE_GROUND, MAX_RECYCLE_YEARS
   use meds_radiation_lib,  only : SWPART_PASSTHROUGH
   use meds_met_types, only : met_source_t, HUMIDITY_RHAIR, HUMIDITY_TDEW, FLD_QAIR,          &
                                  FLD_RHAIR, FLD_TDEW
   use meds_netcdf_c, only : nc_inq_varid_f, nc_get_att_text_f, nc_get_att_double_f, NC_NOERR,    &
                             NC_GLOBAL
   implicit none
   private

   public :: MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS, MET_ERR_START_NOT_A_RECORD,                  &
             MET_ERR_WINDOW_NOT_COVERED, MET_ERR_DT_MISMATCH, MET_ERR_AXIS_NOT_UNIFORM,           &
             MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, MET_ERR_CO2_FILE, MET_ERR_CO2_NOT_COVERED,   &
             MET_ERR_CO2_IN_MET_FILE, MET_ERR_NOT_UTC, MET_ERR_HUMIDITY, MET_ERR_FILE
   public :: REC_MATCH_TOL, validate_recycle_window, file_lookup_sec, first_record_after,         &
             validate_file_against_config, humidity_field, add_line

   !----- met_open status codes (see the `stat` argument). ---------------------------------!
   integer(ik), parameter :: MET_OK                          = 0_ik
   integer(ik), parameter :: MET_ERR_WINDOW_NOT_WHOLE_YEARS  = 1_ik   !< recycle_end - recycle_start /= N years
   integer(ik), parameter :: MET_ERR_START_NOT_A_RECORD      = 2_ik   !< recycle_start is not a record stamp
   integer(ik), parameter :: MET_ERR_WINDOW_NOT_COVERED      = 3_ik   !< file stops short of recycle_end
   integer(ik), parameter :: MET_ERR_DT_MISMATCH             = 4_ik   !< [forcing].timestep /= the file's record spacing
   integer(ik), parameter :: MET_ERR_AXIS_NOT_UNIFORM        = 5_ik   !< the file's time axis is ragged
   integer(ik), parameter :: MET_ERR_ATTR_MISMATCH           = 6_ik   !< a file global attribute contradicts the config
   integer(ik), parameter :: MET_ERR_ARCHIVE                 = 7_ik   !< archive: static file, site cell or a month missing
   integer(ik), parameter :: MET_ERR_CO2_FILE                = 8_ik   !< the CO2 file is unreadable or breaks format 1
   integer(ik), parameter :: MET_ERR_CO2_NOT_COVERED         = 9_ik   !< the CO2 file does not cover the run
   integer(ik), parameter :: MET_ERR_CO2_IN_MET_FILE         = 10_ik  !< the met file carries CO2air
   integer(ik), parameter :: MET_ERR_NOT_UTC                 = 11_ik  !< the file does not say time_zone = "UTC"
   integer(ik), parameter :: MET_ERR_HUMIDITY                = 12_ik  !< no, two, or out-of-range humidity variables
   integer(ik), parameter :: MET_ERR_FILE                    = 13_ik  !< the file cannot be read as a MEDS forcing file

   !----- Tolerance for "this record stamp is that instant" [s]: the axis is float seconds, so an  !
   !      exact == would be brittle, and half a second is far below any forcing interval. ---------!
   real(wp), parameter :: REC_MATCH_TOL = 0.5_wp
   !----- A file's declared measurement height agrees with the config within this [m]. -----------!
   real(wp), parameter :: HEIGHT_MATCH_TOL = 0.01_wp

contains

   !=======================================================================================!
   !  Check the declared recycle window ([forcing].recycle_start/recycle_end) against the file  !
   !  just opened. MEDS never infers the window from the file: it is told where the cycle starts  !
   !  and how long it is, and checks that the file agrees (docs/science/forcing.md §9).          !
   !=======================================================================================!
   subroutine validate_recycle_window(src, stat, why)
      type(met_source_t),            intent(inout) :: src
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(inout) :: why     !< the reason, when stat /= MET_OK
      character(len=320) :: buf
      real(wp)    :: start_sec, end_sec
      integer(ik) :: i
      logical     :: covered

      stat = MET_OK
      src%cycle_anchor = src%fcfg%recycle_start
      src%n_cycle_years = 0_ik
      src%irec_cycle_first = 0_ik
      src%irec_cycle_last = 0_ik
      if (.not. src%fcfg%recycle) return

      !----- The whole-year span, checked here as well as in validate_config, so a driver built   !
      !      directly from a forcing_params_t (as the unit tests build one) cannot skip it. -------!
      src%n_cycle_years = whole_years_between(src%fcfg%recycle_start, src%fcfg%recycle_end,     &
                                              MAX_RECYCLE_YEARS)
      if (src%n_cycle_years < 1_ik) then
         write(buf,'(4a)') ' met_open: recycle window ', time_to_string(src%fcfg%recycle_start),   &
                           ' .. ', time_to_string(src%fcfg%recycle_end)
         call add_line(why, trim(buf))
         call add_line(why, '   is not an exact whole number of calendar years.')
         stat = MET_ERR_WINDOW_NOT_WHOLE_YEARS
         return
      end if

      start_sec = seconds_between(src%base_time, src%fcfg%recycle_start)
      end_sec   = seconds_between(src%base_time, src%fcfg%recycle_end)

      !----- recycle_start must be a record stamp: a config saying 00:00:00 against a file stamped !
      !      01:00:00 is reported with both values, not re-derived. ----------------------------------!
      do i = 1_ik, src%nrec
         if (abs(src%time_sec(i) - start_sec) <= REC_MATCH_TOL) then
            src%irec_cycle_first = i
            exit
         end if
      end do
      if (src%irec_cycle_first == 0_ik) then
         call add_line(why, ' met_open: forcing.recycle_start = '//time_to_string(src%fcfg%recycle_start))
         call add_line(why, '   does not match any record in the forcing file. The file record stamps run')
         write(buf,'(4a)') '   from ', time_to_string(time_advance_seconds(src%base_time, src%time_sec(1))), &
                           ' to ',    time_to_string(time_advance_seconds(src%base_time, src%time_sec(src%nrec)))
         call add_line(why, trim(buf))
         call add_line(why, '   Set recycle_start to an actual record stamp -- MEDS does not guess it. Note')
         call add_line(why, '   avg_convention="end" files (ERA5-Land) stamp the END of each interval, so the')
         call add_line(why, '   first record of a calendar year is 01:00:00 for hourly data, not 00:00:00.')
         stat = MET_ERR_START_NOT_A_RECORD
         return
      end if

      !----- The file must hold a record in the last interval below recycle_end; otherwise the   !
      !      seam bracket would interpolate across a hole. ----------------------------------------!
      do i = src%nrec, src%irec_cycle_first, -1_ik
         if (src%time_sec(i) < end_sec - REC_MATCH_TOL) then
            src%irec_cycle_last = i
            exit
         end if
      end do
      !----- Two steps, not one .and.: Fortran need not short-circuit, and irec_cycle_last may be 0. !
      covered = src%irec_cycle_last > src%irec_cycle_first
      if (covered) covered = src%time_sec(src%irec_cycle_last) >= end_sec - src%dt_forcing - REC_MATCH_TOL
      if (.not. covered) then
         write(buf,'(4a)') ' met_open: the forcing file does not cover the declared recycle window ', &
                           time_to_string(src%fcfg%recycle_start), ' .. ',                         &
                           time_to_string(src%fcfg%recycle_end)
         call add_line(why, trim(buf))
         call add_line(why, '   last record in the file: '//                                         &
                            time_to_string(time_advance_seconds(src%base_time, src%time_sec(src%nrec))))
         stat = MET_ERR_WINDOW_NOT_COVERED
         return
      end if

      write(*,'(a,i0,4a)') ' force : recycling ', src%n_cycle_years, ' calendar year(s), ',      &
            time_to_string(src%fcfg%recycle_start), ' .. ', time_to_string(src%fcfg%recycle_end)
   end subroutine validate_recycle_window

   !----- A model instant's seconds on the file axis, for both met_advance and met_instant so they !
   !      agree: mapped into the declared window when recycling (recycle_model_to_file), otherwise !
   !      the instant itself. ---------------------------------------------------------------------!
   pure function file_lookup_sec(src, now) result(s)
      type(met_source_t), intent(in)  :: src
      type(meds_time_t),  intent(in) :: now
      real(wp) :: s
      if (src%fcfg%recycle .and. src%n_cycle_years >= 1_ik) then
         s = seconds_between(src%base_time, recycle_model_to_file(src, now))
      else
         s = seconds_between(src%base_time, now)
      end if
   end function file_lookup_sec

   !----- Map a model instant into the declared recycle window by substituting the year, relative !
   !      to the window's anchor, with month, day and time of day kept exact (docs/science/forcing.md !
   !      eq. (16)). Feb 29 reads Feb 28 when the target file year is not a leap year. -------------!
   pure function recycle_model_to_file(src, now) result(now_file)
      type(met_source_t), intent(in)  :: src
      type(meds_time_t),  intent(in) :: now
      type(meds_time_t) :: now_file
      integer(ik) :: yf, off
      off = merge(0_ik, 1_ik, mdhms_at_or_after(now, src%cycle_anchor))
      yf  = src%cycle_anchor%year + off                                                          &
          + modulo(now%year - src%cycle_anchor%year - off, src%n_cycle_years)
      now_file = now
      now_file%year = yf
      if (now%month == 2_ik .and. now%day == 29_ik .and. .not. is_leap_year(yf)) now_file%day = 28_ik
   end function recycle_model_to_file

   !----- Compare two instants on (month, day, hour, minute, second) only: their place in the    !
   !      calendar year, whatever the year. --------------------------------------------------------!
   pure logical function mdhms_at_or_after(a, b) result(yes)
      type(meds_time_t), intent(in) :: a, b
      integer(ik) :: ka, kb
      ka = ((a%month*32_ik + a%day)*24_ik + a%hour)*3600_ik + a%minute*60_ik + a%second
      kb = ((b%month*32_ik + b%day)*24_ik + b%hour)*3600_ik + b%minute*60_ik + b%second
      yes = ka >= kb
   end function mdhms_at_or_after

   !----- The first record strictly after file time s (the last record if none). ------------------!
   pure integer(ik) function first_record_after(src, s) result(r)
      type(met_source_t), intent(in)  :: src
      real(wp),           intent(in) :: s
      integer(ik) :: lo, hi, mid
      lo = 1_ik
      hi = src%nrec
      do while (lo < hi)
         mid = (lo + hi) / 2_ik
         if (src%time_sec(mid) > s + REC_MATCH_TOL) then
            hi = mid
         else
            lo = mid + 1_ik
         end if
      end do
      r = lo
   end function first_record_after

   !----- The file variable that carries each humidity form. -----------------------------------!
   pure integer(ik) function humidity_field(form) result(fld)
      integer(ik), intent(in) :: form
      select case (form)
      case (HUMIDITY_RHAIR)
         fld = FLD_RHAIR
      case (HUMIDITY_TDEW)
         fld = FLD_TDEW
      case default
         fld = FLD_QAIR
      end select
   end function humidity_field

   !---------------------------------------------------------------------------------------!
   ! validate_file_against_config -- the forcing file describes itself, and so does the config;   !
   ! when they disagree the run is wrong in a way no later check can see, so stop here. The        !
   ! record spacing matters most; the attributes are checked when present (docs/science/forcing.md !
   ! §1).                                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine validate_file_against_config(src, ncid, vstat, why)
      type(met_source_t),            intent(in)    :: src
      integer(c_int),                intent(in)    :: ncid
      integer(ik),                   intent(out)   :: vstat
      character(len=:), allocatable, intent(inout) :: why     !< the reason, when vstat /= MET_OK
      character(len=320) :: buf
      character(len=64) :: attr
      real(wp)          :: dt_file
      integer(c_int)    :: st
      integer(ik)       :: i

      vstat = MET_OK

      !----- (a) The record spacing: the first interval against the config, then every other     !
      !      interval against the first (MEDS reads a uniform axis only). -------------------------!
      if (src%nrec >= 2_ik) then
         dt_file = src%time_sec(2) - src%time_sec(1)
         if (abs(dt_file - src%dt_forcing) > REC_MATCH_TOL) then
            write(buf,'(a,f12.3,a)') ' met_open: [forcing].timestep = ', src%dt_forcing, ' s'
            call add_line(why, trim(buf))
            write(buf,'(a,f12.3,a)') '   but the file record spacing is ', dt_file, ' s.'
            call add_line(why, trim(buf))
            call add_line(why, '   The timestep places interval midpoints, disaggregates shortwave and')
            call add_line(why, '   brackets the recycle seam. A wrong value mis-times all three.')
            vstat = MET_ERR_DT_MISMATCH
            return
         end if
         do i = 3_ik, src%nrec
            if (abs((src%time_sec(i) - src%time_sec(i-1)) - dt_file) > REC_MATCH_TOL) then
               write(buf,'(a,i0,a)') ' met_open: forcing record spacing changes at record ', i, '.'
               call add_line(why, trim(buf))
               call add_line(why, '   MEDS assumes a uniform time axis.')
               vstat = MET_ERR_AXIS_NOT_UNIFORM
               return
            end if
         end do
      end if

      !----- (b) avg_convention: which end of the interval the file's flux means belong to. ------!
      st = nc_get_att_text_f(ncid, NC_GLOBAL, 'avg_convention', attr)
      if (st == NC_NOERR .and. len_trim(attr) > 0) then
         if (trim(attr) /= trim(avg_convention_name(src%fcfg%avg_convention))) then
            call add_line(why, ' met_open: file avg_convention = "'//trim(attr)//                &
                               '", config = "'//trim(avg_convention_name(src%fcfg%avg_convention))//'"')
            vstat = MET_ERR_ATTR_MISMATCH
            return
         end if
      end if

      !----- (c) sw_input_kind: "total" needs a partition, "components" must not be partitioned.   !
      st = nc_get_att_text_f(ncid, NC_GLOBAL, 'sw_input_kind', attr)
      if (st == NC_NOERR .and. len_trim(attr) > 0) then
         if (trim(attr) == 'total' .and. src%fcfg%sw_partition == SWPART_PASSTHROUGH) then
            call add_line(why, ' met_open: the file carries TOTAL shortwave (sw_input_kind = "total")')
            call add_line(why, '   but [forcing].sw_partition = "passthrough", which expects the four')
            call add_line(why, '   component streams. The run would read fields the file does not have.')
            vstat = MET_ERR_ATTR_MISMATCH
            return
         end if
         if (trim(attr) == 'components' .and. src%fcfg%sw_partition /= SWPART_PASSTHROUGH) then
            call add_line(why, ' met_open: the file already carries the four shortwave components')
            call add_line(why, '   (sw_input_kind = "components") but [forcing].sw_partition asks for a')
            call add_line(why, '   split. Partitioning an already-split stream double-counts the beam.')
            vstat = MET_ERR_ATTR_MISMATCH
            return
         end if
      end if

      !----- (d) The clock: a MEDS forcing file must say time_zone = "UTC" (docs/science/forcing.md !
      !      §1); the archive is UTC by construction. --------------------------------------------!
      if (src%backend == MET_BACKEND_DEFAULT) then
         st = nc_get_att_text_f(ncid, NC_GLOBAL, 'time_zone', attr)
         if (st /= NC_NOERR) attr = '(absent)'
         if (trim(attr) /= 'UTC') then
            call add_line(why, ' met_open: the forcing file''s time_zone attribute is "'//trim(attr)//'".')
            call add_line(why, '   MEDS runs in UTC: build the file on a UTC clock and set time_zone = "UTC".')
            vstat = MET_ERR_NOT_UTC
            return
         end if
      end if

      !----- (e) The heights the file states, against those the config moves the forcing from    !
      !      (docs/science/forcing.md §8). ------------------------------------------------------!
      call check_height_attribute(ncid, 'tq_height_m',        src%fcfg%tq_height,   vstat, why)
      if (vstat /= MET_OK) return
      call check_height_attribute(ncid, 'wind_height_m',      src%fcfg%wind_height, vstat, why)
      if (vstat /= MET_OK) return
      call check_height_attribute(ncid, 'wind_meas_height_m', src%fcfg%wind_height, vstat, why)
      if (vstat /= MET_OK) return
      st = nc_get_att_text_f(ncid, NC_GLOBAL, 'height_above', attr)
      if (st == NC_NOERR .and. len_trim(attr) > 0) then
         if (trim(attr) /= trim(height_above_name(src%fcfg%height_above))) then
            call add_line(why, ' met_open: file height_above = "'//trim(attr)//'", config = "'//      &
                               trim(height_above_name(src%fcfg%height_above))//'"')
            vstat = MET_ERR_ATTR_MISMATCH
            return
         end if
      end if

      !----- (f) CO2air: CO2 comes from [forcing].co2_source, never the met file, so a file that  !
      !      carries it is rejected rather than ignored (docs/science/forcing.md §12). ------------!
      block
         integer(c_int) :: vid
         if (nc_inq_varid_f(ncid, 'CO2air', vid) == NC_NOERR) then
            call add_line(why, ' met_open: the forcing file carries CO2air, but CO2 comes from')
            call add_line(why, '   [forcing].co2_source ("const" or "file"), never from the met file.')
            call add_line(why, '   Remove it (ncks -x -v CO2air in.nc out.nc) or rebuild the file with')
            call add_line(why, '   make_forcing_file.py, which does not write it.')
            vstat = MET_ERR_CO2_IN_MET_FILE
            return
         end if
      end block
   end subroutine validate_file_against_config

   !----- A file's stated measurement height (a global attribute in m), when present, against the   !
   !      configured one.                                                                             !
   subroutine check_height_attribute(ncid, name, configured, vstat, why)
      integer(c_int),                intent(in)    :: ncid
      character(len=*),              intent(in)    :: name
      real(wp),                      intent(in)    :: configured
      integer(ik),                   intent(out)   :: vstat
      character(len=:), allocatable, intent(inout) :: why
      character(len=320) :: buf
      real(c_double) :: stated
      vstat = MET_OK
      if (nc_get_att_double_f(ncid, NC_GLOBAL, name, stated) /= NC_NOERR) return
      if (abs(real(stated, wp) - configured) > HEIGHT_MATCH_TOL) then
         write(buf,'(3a,f0.3,a,f0.3,a)') ' met_open: the forcing file states ', name, ' = ', stated,     &
                                         ' m, but [forcing] says ', configured, ' m.'
         call add_line(why, trim(buf))
         vstat = MET_ERR_ATTR_MISMATCH
      end if
   end subroutine check_height_attribute

   !----- One more line of a reader's reason for refusing a file; the opener prints it. ---------!
   pure subroutine add_line(why, line)
      character(len=:), allocatable, intent(inout) :: why
      character(len=*),              intent(in)    :: line
      if (.not. allocated(why)) why = ''
      if (len(why) > 0) why = why//new_line('a')
      why = why//line
   end subroutine add_line

   pure function height_above_name(code) result(nm)
      integer(ik), intent(in) :: code
      character(len=10) :: nm
      nm = 'zero_plane'
      if (code == HEIGHT_ABOVE_GROUND) nm = 'ground'
   end function height_above_name

   !----- The config code's own spelling, so the mismatch message quotes both sides in one vocabulary. !
   pure function avg_convention_name(code) result(nm)
      integer(ik), intent(in) :: code
      character(len=8) :: nm
      select case (code)
      case (METAVG_BEGIN)
         nm = 'begin'
      case default
         nm = 'end'
      end select
   end function avg_convention_name

end module meds_met_share
