! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_met_co2 -- the prescribed free-atmosphere CO2 (docs/science/forcing.md §12): read a       !
! ED-MORE CO2 file, format 1, into a co2_series_t (co2_series_read) and look it up at a model      !
! instant (co2_series_at). Each row is a period's mean, placed at the period's middle: the lookup  !
! is linear between middles and holds the end values over the outer half-periods (eq. (18)). A    !
! missing period is an error, never interpolated. met_open checks once that the file covers the    !
! run (co2_series_covers). The lookup is on model time, so the CO2 rises while recycled met repeats. !
!==========================================================================================!
module ed_met_co2
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ed_kinds,             only : wp, ik
   use ed_time,              only : ed_time_t, string_to_time, time_to_string, time_eq, seconds_between,    &
                                    time_advance_seconds, time_advance_days, time_advance_months,           &
                                    time_advance_years
   use ed_met_types,   only : co2_series_t
   use ed_numerics,        only : linear_interpolate
   implicit none
   private

   public :: co2_series_read, co2_series_at, co2_series_covers, co2_series_end, co2_series_free

   !----- The period units a `timestep` line may name, and the separators a period start carries  !
   !      at each unit's precision, in order (a year is digits alone).  ---------------------------!
   integer(ik), parameter :: UNIT_YEAR = 1_ik, UNIT_MONTH = 2_ik, UNIT_DAY = 3_ik,               &
                             UNIT_HOUR = 4_ik, UNIT_MINUTE = 5_ik
   character(len=6), parameter :: UNIT_NAME(5) = [character(len=6) ::                            &
                                  'year', 'month', 'day', 'hour', 'minute']
   character(len=4), parameter :: START_SEPS(5) = [character(len=4) :: '', '-', '--', '--T', '--T:']
   character(len=*), parameter :: START_FORM(5) = [character(len=16) ::                          &
                                  'YYYY', 'YYYY-MM', 'YYYY-MM-DD', 'YYYY-MM-DDThh', 'YYYY-MM-DDThh:mm']
   integer, parameter :: LINE_LEN = 1024

contains

   !=======================================================================================!
   !  READ: parse and check an ED-MORE CO2 file. On failure `ok` is false, `s` is empty and     !
   !  `msg` says why, naming the line; the caller decides whether that stops the run (met_open   !
   !  does).                                                                                     !
   !=======================================================================================!
   subroutine co2_series_read(path, s, ok, msg)
      character(len=*),   intent(in)  :: path
      type(co2_series_t), intent(out) :: s
      logical,            intent(out) :: ok
      character(len=*),   intent(out) :: msg
      character(len=LINE_LEN) :: raw, line, tok, rest, ntok, utok, mid, extra
      type(ed_time_t) :: t_start, t_next
      real(wp)    :: v
      integer     :: u, ios, iline
      integer(ik) :: nrow, i, step_n, step_unit
      logical     :: seen_step, seen_units

      ok = .false.
      msg = ''
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         msg = 'cannot open the file'
         return
      end if

      !----- Pass 1: count the data rows, so the series is allocated once. --------------------!
      nrow = 0_ik
      do
         read(u, '(a)', iostat=ios) raw
         if (ios /= 0) exit
         call clean_line(raw, line)
         if (len_trim(line) == 0) cycle
         call split_token(line, tok, rest)
         if (trim(tok) /= 'timestep' .and. trim(tok) /= 'units') nrow = nrow + 1_ik
      end do
      if (nrow < 2_ik) then
         close(u)
         msg = 'fewer than two data rows'
         return
      end if
      allocate(s%mid_sec(nrow), s%co2(nrow))
      s%n = nrow

      !----- Pass 2: the keyword lines, then the rows, each checked as it is read. --------------!
      rewind(u)
      seen_step = .false.
      seen_units = .false.
      step_n = 0_ik
      step_unit = 0_ik
      i = 0_ik
      iline = 0
      do
         read(u, '(a)', iostat=ios) raw
         if (ios /= 0) exit
         iline = iline + 1
         call clean_line(raw, line)
         if (len_trim(line) == 0) cycle
         call split_token(line, tok, rest)
         select case (trim(tok))
         case ('timestep')
            if (seen_step .or. i > 0_ik) then
               call line_msg(msg, iline, 'a timestep line must come once, before the first data row')
               exit
            end if
            call split_token(rest, ntok, mid)
            call split_token(mid, utok, extra)
            read(ntok, *, iostat=ios) step_n
            !----- `utok` untrimmed: findloc of a shorter value returns 0 (docs/building.md).           !
            step_unit = int(findloc(UNIT_NAME, utok, dim=1), ik)
            if (ios /= 0 .or. step_n < 1_ik .or. step_unit == 0_ik .or. len_trim(extra) > 0) then
               call line_msg(msg, iline, 'expected "timestep <n> <unit>", n a whole number >= 1 and '// &
                             'unit one of year, month, day, hour, minute')
               exit
            end if
            seen_step = .true.
         case ('units')
            if (seen_units .or. i > 0_ik) then
               call line_msg(msg, iline, 'a units line must come once, before the first data row')
               exit
            end if
            if (trim(rest) /= 'umol/mol') then
               call line_msg(msg, iline, 'expected "units umol/mol", the only unit accepted')
               exit
            end if
            seen_units = .true.
         case default
            if (.not. (seen_step .and. seen_units)) then
               call line_msg(msg, iline, 'a data row comes before the timestep and units lines')
               exit
            end if
            call parse_start(tok, step_unit, t_start, ok)
            if (.not. ok) then
               call line_msg(msg, iline, 'the period start "'//trim(tok)//'" is not a valid '//         &
                             trim(START_FORM(step_unit))//' (the precision of "timestep ... '//        &
                             trim(UNIT_NAME(step_unit))//'")')
               exit
            end if
            call split_token(rest, ntok, extra)
            read(ntok, *, iostat=ios) v
            if (ios /= 0 .or. len_trim(extra) > 0) then
               call line_msg(msg, iline, 'expected "<period start> <value>"')
               exit
            end if
            if (.not. ieee_is_finite(v) .or. v <= 0.0_wp) then
               call line_msg(msg, iline, 'the value must be finite and above zero')
               exit
            end if
            i = i + 1_ik
            if (i == 1_ik) then
               s%first_start = t_start
            else if (.not. time_eq(t_start, t_next)) then
               call line_msg(msg, iline, 'the period starts at '//time_to_string(t_start)//            &
                             ', but the one before it ends at '//time_to_string(t_next)//              &
                             '. Rows must be consecutive; ED-MORE never gap-fills')
               exit
            end if
            t_next = advance_period(t_start, step_n, step_unit)
            s%mid_sec(i) = 0.5_wp * (seconds_between(s%first_start, t_start)                       &
                                     + seconds_between(s%first_start, t_next))
            s%co2(i)     = v
         end select
      end do
      close(u)

      ok = len_trim(msg) == 0
      if (.not. ok) then
         call co2_series_free(s)
         return
      end if
      s%span_end_sec = seconds_between(s%first_start, t_next)
   end subroutine co2_series_read

   !=======================================================================================!
   !  LOOKUP: the CO2 at model instant `now` [umol/mol]. Linear between period middles; the end  !
   !  values are held beyond the first and last middles (a covered run only reaches those over   !
   !  the file's outer half-periods).                                                            !
   !=======================================================================================!
   pure real(wp) function co2_series_at(s, now) result(co2)
      type(co2_series_t), intent(in) :: s
      type(ed_time_t),    intent(in) :: now
      real(wp)    :: x
      integer(ik) :: lo, hi, k
      x = seconds_between(s%first_start, now)
      if (x <= s%mid_sec(1)) then
         co2 = s%co2(1)
         return
      end if
      if (x >= s%mid_sec(s%n)) then
         co2 = s%co2(s%n)
         return
      end if
      lo = 1_ik
      hi = s%n                                   ! mid_sec(lo) <= x < mid_sec(hi)
      do while (hi - lo > 1_ik)
         k = (lo + hi) / 2_ik
         if (s%mid_sec(k) <= x) then
            lo = k
         else
            hi = k
         end if
      end do
      co2 = linear_interpolate(s%co2(lo), s%co2(hi), (x - s%mid_sec(lo)) / (s%mid_sec(hi) - s%mid_sec(lo)))
   end function co2_series_at

   !----- Does the series cover [t0, t1] -- from the first period's start to the last one's end? --!
   pure logical function co2_series_covers(s, t0, t1) result(yes)
      type(co2_series_t), intent(in) :: s
      type(ed_time_t),    intent(in) :: t0, t1
      yes = seconds_between(s%first_start, t0) >= 0.0_wp .and.                                   &
            seconds_between(s%first_start, t1) <= s%span_end_sec
   end function co2_series_covers

   !----- The end of the last period (for messages). ----------------------------------------!
   pure function co2_series_end(s) result(t)
      type(co2_series_t), intent(in) :: s
      type(ed_time_t) :: t
      t = time_advance_seconds(s%first_start, s%span_end_sec)
   end function co2_series_end

   subroutine co2_series_free(s)
      type(co2_series_t), intent(inout) :: s
      if (allocated(s%mid_sec)) deallocate(s%mid_sec)
      if (allocated(s%co2))     deallocate(s%co2)
      s%n = 0_ik
      s%span_end_sec = 0.0_wp
   end subroutine co2_series_free

   !----- The start of the next period: n calendar units on. --------------------------------!
   pure function advance_period(t, n, unit) result(t2)
      type(ed_time_t),   intent(in) :: t
      integer(ik),       intent(in) :: n, unit
      type(ed_time_t) :: t2
      select case (unit)
      case (UNIT_YEAR)
         t2 = time_advance_years(t, n)
      case (UNIT_MONTH)
         t2 = time_advance_months(t, n)
      case (UNIT_DAY)
         t2 = time_advance_days(t, n)
      case (UNIT_HOUR)
         t2 = time_advance_seconds(t, 3600.0_wp * real(n, wp))
      case default
         t2 = time_advance_seconds(t, 60.0_wp * real(n, wp))
      end select
   end function advance_period

   !----- A period start at exactly the unit's precision: digit groups joined by the unit's       !
   !      separators in order (e.g. "--T" for an hour), none empty, then a valid calendar instant. -!
   subroutine parse_start(label, unit, t, ok)
      character(len=*),  intent(in)  :: label
      integer(ik),       intent(in)  :: unit
      type(ed_time_t),   intent(out) :: t
      logical,           intent(out) :: ok
      character(len=4) :: seps
      integer :: j, nsep
      logical :: prev_digit
      ok = .false.
      seps = ''
      nsep = 0
      prev_digit = .false.
      t = ed_time_t()
      do j = 1, len_trim(label)
         select case (label(j:j))
         case ('0':'9')
            prev_digit = .true.
         case ('-', 'T', ':')
            if (.not. prev_digit .or. nsep >= len(seps)) return    ! an empty group, or too many
            nsep = nsep + 1
            seps(nsep:nsep) = label(j:j)
            prev_digit = .false.
         case default
            return
         end select
      end do
      if (.not. prev_digit .or. seps /= START_SEPS(unit)) return
      call string_to_time(label, t, ok)
   end subroutine parse_start

   !----- Tabs to blanks, the '#' comment dropped, left-justified. --------------------------!
   pure subroutine clean_line(raw, line)
      character(len=*), intent(in)  :: raw
      character(len=*), intent(out) :: line
      integer :: j
      line = raw
      do j = 1, len_trim(line)
         if (line(j:j) == char(9)) line(j:j) = ' '
      end do
      j = index(line, '#')
      if (j > 0) line(j:) = ''
      line = adjustl(line)
   end subroutine clean_line

   !----- The first blank-delimited token of `line`, and the rest left-justified. ------------!
   pure subroutine split_token(line, tok, rest)
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: tok, rest
      character(len=len(line)) :: buf
      integer :: k
      buf = adjustl(line)
      k = index(buf, ' ')
      if (k == 0) then
         tok = buf
         rest = ''
      else
         tok = buf(1:k-1)
         rest = adjustl(buf(k:))
      end if
   end subroutine split_token

   subroutine line_msg(msg, iline, text)
      character(len=*), intent(out) :: msg
      integer,          intent(in)  :: iline
      character(len=*), intent(in)  :: text
      write(msg, '(a,i0,2a)') 'line ', iline, ': ', text
   end subroutine line_msg

end module ed_met_co2
