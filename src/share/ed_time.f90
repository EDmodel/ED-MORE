! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_time -- calendar time and the sun's position. A run goes from a start date to an      !
! end date on the proleptic Gregorian calendar, leap years included, and every output      !
! record carries its date. This module does the date arithmetic, through the Julian Day    !
! Number, and the solar geometry that follows from the date and the clock: the             !
! declination, the zenith cosine at local solar or UTC time, its mean over a window, the   !
! sun on a slope and the daylength (docs/science/forcing.md §4).                           !
!==========================================================================================!
module ed_time
   use ed_kinds,       only : wp, ik, lk
   use ed_constants,   only : pi, day_sec, deg_to_rad, tiny_num
   implicit none
   private

   public :: ed_time_t, calendar_step_t, calendar_step
   public :: is_leap_year, days_in_month, days_in_year, day_of_year
   public :: time_advance_days, time_advance_months, years_between
   public :: time_advance_years, whole_years_between
   public :: seconds_into_day, seconds_between, time_advance_seconds
   public :: time_lt, time_eq, time_valid
   public :: string_to_time, time_units_base, time_to_string, time_to_stamp, time_to_decimal_year
   public :: solar_cosz, daylength, doy_effective, solar_declination
   public :: solar_cosz_utc, utc_to_solar_seconds, mean_solar_cosz, slope_incidence_cosz

   !----- A calendar instant; the defaults make ed_time_t() a valid date. ---------------  !
   type :: ed_time_t
      integer(ik) :: year   = 2000_ik
      integer(ik) :: month  = 1_ik
      integer(ik) :: day    = 1_ik
      integer(ik) :: hour   = 0_ik
      integer(ik) :: minute = 0_ik
      integer(ik) :: second = 0_ik
   end type ed_time_t

   !----- One step of the calendar: its start and end, its length in days, and the boundaries it !
   !      crosses. A step that reaches a new month also reaches a new day, and one that reaches a !
   !      new year a new month. calendar_step builds it. -------------------------------------!
   type :: calendar_step_t
      type(ed_time_t)   :: start, end
      integer(ik)       :: days = 1_ik
      logical           :: new_day = .false., new_month = .false., new_year = .false.
   end type calendar_step_t

   integer(ik), parameter :: MONTH_LEN(12) = [31_ik, 28_ik, 31_ik, 30_ik, 31_ik, 30_ik,       &
                                               31_ik, 31_ik, 30_ik, 31_ik, 30_ik, 31_ik]

contains

   !----- The proleptic-Gregorian leap-year rule. -----------------------------------------!
   elemental logical function is_leap_year(year) result(leap)
      integer(ik), intent(in) :: year
      leap = (mod(year, 4_ik) == 0_ik .and. mod(year, 100_ik) /= 0_ik) .or.                   &
             mod(year, 400_ik) == 0_ik
   end function is_leap_year

   elemental integer(ik) function days_in_month(year, month) result(d)
      integer(ik), intent(in) :: year, month
      d = MONTH_LEN(month)
      if (month == 2_ik .and. is_leap_year(year)) d = 29_ik
   end function days_in_month

   elemental integer(ik) function days_in_year(year) result(d)
      integer(ik), intent(in) :: year
      d = 365_ik
      if (is_leap_year(year)) d = 366_ik
   end function days_in_year

   !----- The day of the year, 1 on 1 January. --------------------------------------------!
   pure integer(ik) function day_of_year(t) result(doy)
      type(ed_time_t), intent(in) :: t
      type(ed_time_t) :: jan1
      jan1 = ed_time_t(year = t%year)          ! 1 January, 00:00, of the same year
      doy  = int(to_jdn(t) - to_jdn(jan1), ik) + 1_ik
   end function day_of_year

   !----- The solar declination [rad], Cooper (1969): 23.45 deg sin(2 pi (284 + doy)/365)  !
   !      (docs/science/forcing.md eq. (4)). It is the model's only declination:           !
   !      solar_cosz and daylength both use it. -------------------------------------------!
   elemental real(wp) function solar_declination(doy) result(decl)
      integer(ik), intent(in) :: doy          !< day of year (1 = 1 January)
      decl = 23.45_wp * pi / 180.0_wp * sin(2.0_wp * pi * (284.0_wp + real(doy, wp)) / 365.0_wp)
   end function solar_declination

   !----- The cosine of the solar zenith angle, 0 at night (docs/science/forcing.md        !
   !      eq. (4)). t_sec is local apparent solar time, which the forcing makes from UTC,  !
   !      the longitude and the equation of time (eqs. (2)-(3)). --------------------------!
   pure real(wp) function solar_cosz(t, t_sec, latitude_deg) result(cosz)
      type(ed_time_t),   intent(in) :: t
      real(wp),          intent(in) :: t_sec          !< [s] local solar seconds into the day (0 = midnight)
      real(wp),          intent(in) :: latitude_deg   !< [deg] site latitude (+N)
      real(wp) :: deg2rad, lat, decl, hourangle, frac_day
      deg2rad   = pi / 180.0_wp
      lat       = latitude_deg * deg2rad
      decl      = solar_declination(day_of_year(t))
      frac_day  = t_sec / day_sec                  ! 0 at midnight, 0.5 at solar noon
      hourangle = 2.0_wp * pi * (frac_day - 0.5_wp)   ! -pi (midnight) .. 0 (noon) .. +pi
      cosz      = sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hourangle)
      cosz      = max(cosz, 0.0_wp)
   end function solar_cosz

   !----- The equation of time [s], Spencer (1971) (docs/science/forcing.md eq. (3)). -----!
   pure function equation_of_time(t) result(eot_sec)
      type(ed_time_t), intent(in) :: t
      real(wp) :: eot_sec, b
      b = 2.0_wp * pi * real(day_of_year(t) - 1_ik, wp) / 365.0_wp        ! Spencer 1971
      eot_sec = 229.18_wp * ( 0.000075_wp + 0.001868_wp*cos(b) - 0.032077_wp*sin(b)             &
                            - 0.014615_wp*cos(2.0_wp*b) - 0.040849_wp*sin(2.0_wp*b) ) * 60.0_wp
   end function equation_of_time

   !----- Local apparent solar seconds into the day at a UTC instant: the UTC seconds, plus    !
   !      240 s per degree of longitude, plus the equation of time (docs/science/forcing.md    !
   !      eqs. (2)-(3)). Every forcing clock is UTC. -------------------------------------------!
   pure function utc_to_solar_seconds(t, sec_utc, longitude_deg) result(s)
      type(ed_time_t),   intent(in) :: t              !< record/model time (for declination + EoT via doy)
      real(wp),          intent(in) :: sec_utc        !< [s] seconds into the UTC day
      real(wp),          intent(in) :: longitude_deg
      real(wp) :: s
      s = sec_utc + longitude_deg * 240.0_wp + equation_of_time(t)
   end function utc_to_solar_seconds

   !----- cos(solar zenith) at a UTC instant, through its apparent solar time. -------------!
   pure function solar_cosz_utc(t, sec_utc, latitude_deg, longitude_deg) result(cosz)
      type(ed_time_t),   intent(in) :: t
      real(wp),          intent(in) :: sec_utc, latitude_deg, longitude_deg
      real(wp) :: cosz, sec_solar
      sec_solar = utc_to_solar_seconds(t, sec_utc, longitude_deg)
      cosz = solar_cosz(t, sec_solar, latitude_deg)
   end function solar_cosz_utc

   !----- The mean of max(cos z, 0) over a window of UTC seconds [win_start_sec, win_start_sec +  !
   !      dt_win), by the midpoint rule on about dt_win / dt_sub sub-samples: night counts as 0    !
   !      (docs/science/forcing.md eq. (6)). ------------------------------------------------------!
   pure real(wp) function mean_solar_cosz(t, win_start_sec, dt_sub, dt_win, latitude_deg, longitude_deg)  &
                                          result(cosz_bar)
      type(ed_time_t),   intent(in) :: t
      real(wp),          intent(in) :: win_start_sec, dt_sub, dt_win
      real(wp),          intent(in) :: latitude_deg, longitude_deg
      real(wp)    :: cosz_sum, sec, ci
      integer(ik) :: nsub, i
      nsub = max(1_ik, nint(dt_win / max(dt_sub, tiny_num), ik))
      cosz_sum = 0.0_wp
      do i = 1_ik, nsub                                   ! midpoint rule over the FULL window
         sec = win_start_sec + (real(i, wp) - 0.5_wp) * (dt_win / real(nsub, wp))
         ci  = solar_cosz_utc(t, sec, latitude_deg, longitude_deg)
         cosz_sum = cosz_sum + max(ci, 0.0_wp)            ! night sub-samples contribute exactly 0
      end do
      cosz_bar = cosz_sum / real(nsub, wp)
   end function mean_solar_cosz

   !----- The sun on a slope (docs/science/forcing.md §8): the cosine of its incidence on a            !
   !      surface tilted `slope_deg` from horizontal and facing `aspect_deg` (clockwise from north),   !
   !      cos i = cos z cos s + sin z sin s cos(A - aspect), with A the sun's azimuth, also clockwise  !
   !      from north. 0 with the sun below the horizon or behind the slope. ---------------------------!
   pure real(wp) function slope_incidence_cosz(t, sec_utc, latitude_deg, longitude_deg, slope_deg,     &
                                               aspect_deg) result(cosi)
      type(ed_time_t),   intent(in) :: t
      real(wp),          intent(in) :: sec_utc, latitude_deg, longitude_deg, slope_deg, aspect_deg
      real(wp) :: lat, decl, hour_angle, cosz, sinz, azimuth, slope
      lat        = latitude_deg * deg_to_rad
      decl       = solar_declination(day_of_year(t))
      hour_angle = 2.0_wp * pi * (utc_to_solar_seconds(t, sec_utc, longitude_deg) / day_sec - 0.5_wp)
      cosz       = sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hour_angle)
      cosi       = 0.0_wp
      if (cosz <= 0.0_wp) return
      sinz    = sqrt(max(0.0_wp, 1.0_wp - cosz * cosz))
      azimuth = atan2(-cos(decl) * sin(hour_angle), sin(decl) * cos(lat) - cos(decl) * sin(lat) * cos(hour_angle))
      slope   = slope_deg * deg_to_rad
      cosi    = max(0.0_wp, cosz * cos(slope) + sinz * sin(slope) * cos(azimuth - aspect_deg * deg_to_rad))
   end function slope_incidence_cosz

   !----- Daylength [h], 24/pi acos(-tan(lat) tan(decl)): 0 in the polar night, 24 in the  !
   !      polar day. ----------------------------------------------------------------------!
   elemental real(wp) function daylength(lat_deg, doy) result(dl)
      real(wp),    intent(in) :: lat_deg
      integer(ik), intent(in) :: doy
      real(wp) :: latr, decl, arg
      latr = lat_deg * pi / 180.0_wp
      decl = solar_declination(doy)
      arg  = -tan(latr) * tan(decl)
      if (arg >= 1.0_wp) then
         dl = 0.0_wp                 ! polar night
      else if (arg <= -1.0_wp) then
         dl = 24.0_wp                ! polar day
      else
         dl = 24.0_wp / pi * acos(arg)
      end if
   end function daylength

   !----- The day of year, shifted by half a year in the southern hemisphere, so a         !
   !      seasonal cue written for the north serves both hemispheres. ---------------------!
   elemental integer(ik) function doy_effective(doy, hemis_north) result(de)
      integer(ik), intent(in) :: doy
      logical,     intent(in) :: hemis_north
      if (hemis_north) then
         de = doy
      else
         de = modulo(doy - 1_ik + 182_ik, 365_ik) + 1_ik
      end if
   end function doy_effective

   !=======================================================================================!
   ! Julian Day Number conversions (Fliegel & Van Flandern 1968). Day numbers are          !
   ! 64-bit, so a span in seconds stays exact over centuries: 32-bit seconds overflow      !
   ! after 68 years.                                                                       !
   !=======================================================================================!
   pure integer(lk) function to_jdn(t) result(jdn)
      type(ed_time_t), intent(in) ::   t
      integer(lk) :: y, m, d
      y = int(t%year, lk) ; m = int(t%month, lk) ; d = int(t%day, lk)
      jdn = (1461_lk * (y + 4800_lk + (m - 14_lk) / 12_lk)) / 4_lk                            &
          + (367_lk * (m - 2_lk - 12_lk * ((m - 14_lk) / 12_lk))) / 12_lk                     &
          - (3_lk * ((y + 4900_lk + (m - 14_lk) / 12_lk) / 100_lk)) / 4_lk                    &
          + d - 32075_lk
   end function to_jdn

   pure subroutine from_jdn(jdn, year, month, day)
      integer(lk), intent(in)  :: jdn
      integer(ik), intent(out) :: year, month, day
      integer(lk) :: l, n, i, j, k
      l = jdn + 68569_lk
      n = (4_lk * l) / 146097_lk
      l = l - (146097_lk * n + 3_lk) / 4_lk
      i = (4000_lk * (l + 1_lk)) / 1461001_lk
      l = l - (1461_lk * i) / 4_lk + 31_lk
      j = (80_lk * l) / 2447_lk
      k = l - (2447_lk * j) / 80_lk
      l = j / 11_lk
      j = j + 2_lk - 12_lk * l
      i = 100_lk * (n - 49_lk) + i + l
      year = int(i, ik) ; month = int(j, ik) ; day = int(k, ik)
   end subroutine from_jdn

   !=======================================================================================!
   ! Date stepping.                                                                        !
   !=======================================================================================!
   !----- Advance by ndays whole days; the time of day is kept. ---------------------------!
   pure function time_advance_days(t, ndays) result(t2)
      type(ed_time_t),   intent(in) :: t
      integer(ik),       intent(in) :: ndays
      type(ed_time_t) :: t2
      t2 = t
      call from_jdn(to_jdn(t) + int(ndays, lk), t2%year, t2%month, t2%day)
   end function time_advance_days

   !----- The calendar step of `days` days from `start`. -----------------------------------!
   pure function calendar_step(start, days) result(step)
      type(ed_time_t),   intent(in) :: start
      integer(ik),       intent(in) :: days
      type(calendar_step_t) :: step
      step%start     = start
      step%end       = time_advance_days(start, days)
      step%days      = days
      step%new_year  = step%end%year /= start%year
      step%new_month = step%new_year .or. step%end%month /= start%month
      step%new_day   = step%new_month .or. step%end%day /= start%day
   end function calendar_step

   !----- Advance by whole months; a day past the end of the target month becomes its last !
   !      day (31 January + 1 month is 28 or 29 February). --------------------------------!
   pure function time_advance_months(t, nmonths) result(t2)
      type(ed_time_t),   intent(in) :: t
      integer(ik),       intent(in) :: nmonths
      type(ed_time_t) :: t2
      integer(ik) :: total
      t2 = t
      total    = t%year * 12_ik + (t%month - 1_ik) + nmonths
      t2%year  = total / 12_ik
      t2%month = mod(total, 12_ik) + 1_ik
      t2%day   = min(t%day, days_in_month(t2%year, t2%month))
   end function time_advance_months

   !----- b - a in decimal calendar years, for progress messages. -------------------------!
   pure real(wp) function years_between(a, b) result(y)
      type(ed_time_t), intent(in)   :: a, b
      y = time_to_decimal_year(b) - time_to_decimal_year(a)
   end function years_between

   !=======================================================================================!
   ! Arithmetic to the second, across day boundaries (the forcing's interpolation window). !
   !=======================================================================================!
   pure integer(lk) function sec_of_day(t) result(s)
      type(ed_time_t), intent(in)   :: t
      s = int(t%hour, lk) * 3600_lk + int(t%minute, lk) * 60_lk + int(t%second, lk)
   end function sec_of_day

   !----- Seconds since midnight of the instant's own day (0 .. 86399). -------------------!
   pure real(wp) function seconds_into_day(t) result(s)
      type(ed_time_t), intent(in) :: t
      s = real(sec_of_day(t), wp)
   end function seconds_into_day

   !----- Seconds from a to b (positive if b is later). -----------------------------------!
   pure real(wp) function seconds_between(a, b) result(s)
      type(ed_time_t), intent(in) ::   a, b
      s = real((to_jdn(b) - to_jdn(a)) * 86400_lk + sec_of_day(b) - sec_of_day(a), wp)
   end function seconds_between

   !----- Advance by dsec seconds, rounded to the nearest second; dsec may be negative or  !
   !      span many days. -----------------------------------------------------------------!
   pure function time_advance_seconds(t, dsec) result(t2)
      type(ed_time_t),   intent(in) :: t
      real(wp),          intent(in) :: dsec
      type(ed_time_t) :: t2
      integer(lk) :: total, ndays, rem
      total = sec_of_day(t) + nint(dsec, lk)             ! seconds from midnight of t's day
      ndays = int(floor(real(total, wp) / day_sec), lk)
      rem   = total - ndays * 86400_lk                   ! 0 .. 86399, also when total < 0
      t2 = time_advance_days(t, int(ndays, ik))
      t2%hour   = int(rem / 3600_lk, ik)
      t2%minute = int(mod(rem, 3600_lk) / 60_lk, ik)
      t2%second = int(mod(rem, 60_lk), ik)
   end function time_advance_seconds

   !=======================================================================================!
   ! Comparisons, on the full instant.                                                     !
   !=======================================================================================!
   pure logical function time_lt(a, b) result(yes)
      type(ed_time_t),   intent(in) :: a, b
      integer(lk) :: ja, jb
      ja = to_jdn(a) ; jb = to_jdn(b)
      if (ja /= jb) then
         yes = ja < jb
      else
         yes = sec_of_day(a) < sec_of_day(b)
      end if
   end function time_lt

   pure logical function time_eq(a, b) result(yes)
      type(ed_time_t),   intent(in) :: a, b
      yes = to_jdn(a) == to_jdn(b) .and. sec_of_day(a) == sec_of_day(b)
   end function time_eq

   !----- Advance by n calendar years, keeping the month, day and time of day; 29 February !
   !      in a non-leap target year becomes 28 February, as in the forcing's recycling     !
   !      (docs/science/forcing.md §9). ---------------------------------------------------!
   pure function time_advance_years(t, n) result(out)
      type(ed_time_t),   intent(in) :: t
      integer(ik),       intent(in) :: n
      type(ed_time_t) :: out
      out = t
      out%year = t%year + n
      if (t%month == 2_ik .and. t%day == 29_ik .and. .not. is_leap_year(out%year)) out%day = 28_ik
   end function time_advance_years

   !----- n when b is exactly n calendar years after a (1 <= n <= nmax), else 0. The       !
   !      forcing's recycle window must span whole years, or every wrap shifts the hour of !
   !      day and the season (docs/science/forcing.md §9). --------------------------------!
   pure function whole_years_between(a, b, nmax) result(n)
      type(ed_time_t),   intent(in) :: a, b
      integer(ik),       intent(in) :: nmax
      integer(ik) :: n, k
      n = 0_ik
      do k = 1_ik, nmax
         if (time_eq(time_advance_years(a, k), b)) then ; n = k ; return ; end if
         if (time_lt(b, time_advance_years(a, k))) return       ! past b: not a whole number of years
      end do
   end function whole_years_between

   !----- True when the month, day, hour, minute and second are in range (a parsed config  !
   !      date); any year passes. ---------------------------------------------------------!
   pure logical function time_valid(t) result(ok)
      type(ed_time_t),   intent(in) :: t
      ok = t%month >= 1_ik  .and. t%month <= 12_ik                                           &
           .and. t%day  >= 1_ik  .and. t%day   <= days_in_month(t%year, max(min(t%month,12_ik),1_ik)) &
           .and. t%hour   >= 0_ik .and. t%hour   <= 23_ik                                    &
           .and. t%minute >= 0_ik .and. t%minute <= 59_ik                                    &
           .and. t%second >= 0_ik .and. t%second <= 59_ik
   end function time_valid

   !=======================================================================================!
   ! String <-> time.                                                                      !
   !=======================================================================================!
   !----- Parse "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" ('/', 'T' and ',' also separate). A  !
   !      missing trailing field defaults (month and day to 1, time to 0). ok is false on  !
   !      an unreadable or out-of-range date, so the caller can use its default. ----------!
   subroutine string_to_time(s, t, ok)
      character(len=*),  intent(in)  :: s
      type(ed_time_t),   intent(out) :: t
      logical,           intent(out) :: ok
      character(len=len(s)) :: buf
      integer(ik) :: i, ntok, ios
      integer(ik) :: y, mo, d, h, mi, se

      t = ed_time_t()
      y = 2000_ik ; mo = 1_ik ; d = 1_ik ; h = 0_ik ; mi = 0_ik ; se = 0_ik

      !----- Separators become blanks; a list-directed read then takes the numbers. -------!
      buf = s
      do i = 1_ik, int(len_trim(buf), ik)
         select case (buf(i:i))
         case ('-', '/', ':', 'T', 't', ',') ; buf(i:i) = ' '
         end select
      end do

      ntok = count_tokens(buf)
      ios  = 0_ik
      select case (ntok)
      case (6_ik:) ; read(buf, *, iostat=ios) y, mo, d, h, mi, se
      case (5_ik)  ; read(buf, *, iostat=ios) y, mo, d, h, mi
      case (4_ik)  ; read(buf, *, iostat=ios) y, mo, d, h
      case (3_ik)  ; read(buf, *, iostat=ios) y, mo, d
      case (2_ik)  ; read(buf, *, iostat=ios) y, mo
      case (1_ik)  ; read(buf, *, iostat=ios) y
      case default ; ios = -1_ik
      end select

      t = ed_time_t(year=y, month=mo, day=d, hour=h, minute=mi, second=se)
      ok = (ios == 0_ik) .and. time_valid(t)
   end subroutine string_to_time

   !----- The base date of a CF time-units string, "<unit> since YYYY-MM-DD[ HH:MM:SS]";   !
   !      the caller checks the unit. ok is false without "since" or with an unreadable    !
   !      date. ---------------------------------------------------------------------------!
   subroutine time_units_base(units, base, ok)
      character(len=*),  intent(in)  :: units
      type(ed_time_t),   intent(out) :: base
      logical,           intent(out) :: ok
      integer :: k
      k = index(units, 'since')
      if (k == 0) then ; ok = .false. ; return ; end if
      call string_to_time(adjustl(units(k + 5:)), base, ok)
   end subroutine time_units_base

   !----- "YYYY-MM-DD HH:MM:SS", for attributes and log lines. ----------------------------!
   pure function time_to_string(t) result(s)
      type(ed_time_t), intent(in) :: t
      character(len=19) :: s
      write(s, '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)')                         &
            t%year, t%month, t%day, t%hour, t%minute, t%second
   end function time_to_string

   !----- "YYYYMMDDHHMMSS", the stamp in a state file's name (<prefix>-S-<stamp>.nc). -----!
   pure function time_to_stamp(t) result(s)
      type(ed_time_t), intent(in) :: t
      character(len=14) :: s
      write(s, '(i4.4,5i2.2)') t%year, t%month, t%day, t%hour, t%minute, t%second
   end function time_to_stamp

   !----- The year plus the elapsed fraction of it: a continuous time axis for plots. -----!
   pure real(wp) function time_to_decimal_year(t) result(yr)
      type(ed_time_t), intent(in) :: t
      real(wp) :: frac
      frac = (real(day_of_year(t) - 1_ik, wp)                                                &
              + real(sec_of_day(t), wp) / day_sec) / real(days_in_year(t%year), wp)
      yr = real(t%year, wp) + frac
   end function time_to_decimal_year

   !----- The number of blank-separated tokens in s. --------------------------------------!
   pure integer(ik) function count_tokens(s) result(nt)
      character(len=*), intent(in) :: s
      integer(ik) :: i
      logical     :: in_tok
      nt = 0_ik ; in_tok = .false.
      do i = 1_ik, int(len_trim(s), ik)
         if (s(i:i) /= ' ') then
            if (.not. in_tok) nt = nt + 1_ik
            in_tok = .true.
         else
            in_tok = .false.
         end if
      end do
   end function count_tokens

end module ed_time
