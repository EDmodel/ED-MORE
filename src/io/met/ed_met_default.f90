! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_met_default -- the ED-MORE forcing file source (format = "default",
! docs/science/forcing.md §1): a multi-grid netCDF file whose records the run needs are
! read into memory at open, so no step reads the file. It opens the file (default_file_open)
! and hands the shared ingest one record's values as stored (default_file_record).
!---~---
module ed_met_default
   use iso_c_binding, only : c_int, c_size_t
   use ed_kinds, only : wp, ik
   use ed_time, only : ed_time_t, time_units_base, seconds_between
   use ed_radiation_lib,    only : SWPART_PASSTHROUGH
   use ed_forcing_params,   only : GRIDMATCH_NEAREST, LW_SYNTHESIZE
   use ed_met_types, only : met_record_t, met_source_t, HUMIDITY_QAIR, HUMIDITY_RHAIR,      &
                                  HUMIDITY_TDEW, N_ED_FIELD, ED_FIELD, FLD_TAIR, FLD_RHAIR,       &
                                  FLD_PSURF, FLD_WIND, FLD_U10, FLD_V10, FLD_RAINF, FLD_LWDOWN,   &
                                  FLD_SWDOWN, FLD_PAR_BEAM, FLD_PAR_DIFFUSE, FLD_NIR_BEAM,        &
                                  FLD_NIR_DIFFUSE
   use ed_met_grid, only : nearest_grid_index
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use ed_netcdf_c,   only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f, nc_get_att_text_f,       &
                             nc_get_vara_double, nc_close, nc_inq_var_chunking, nc_error_text,    &
                             nc_packing_t, nc_packing_of, nc_unpacked, NC_NOERR, NC_NOWRITE,      &
                             NC_CHUNKED
   use ed_met_share,   only : MET_OK, MET_ERR_HUMIDITY, MET_ERR_FILE, validate_recycle_window,    &
                              first_record_after, validate_file_against_config, humidity_field,   &
                              add_line
   implicit none
   private

   public :: default_file_open, default_file_record

   !--- RHair is a fraction; a value above this is a percentage written into a fractional field.
   real(wp), parameter :: RH_FRACTION_MAX = 1.5_wp

contains

   !---~---
   !   OPEN an ED-MORE forcing file ("default"): the grid and time dims, the time axis and its
   ! base time from the `time:units` attribute, the checks of the file against the config and
   ! the declared recycle window, and the records the run needs, read into memory. A
   ! rejection returns its status and the reason; the caller closes the source.
   !---~---
   subroutine default_file_open(src, stat, why, run_start, run_end)
      type(met_source_t),            intent(inout) :: src
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(out)   :: why    ! what stops the run when stat /= MET_OK
      type(ed_time_t), optional,   intent(in)      :: run_start, run_end
      integer(c_int)    :: st, ncid
      integer(c_size_t) :: dlen
      character(len=256):: units
      character(len=320):: buf
      logical           :: ok

      stat = MET_OK
      why = ''
      st = nc_open_f(trim(src%fcfg%path), NC_NOWRITE, ncid)
      call check_read(st, 'open', src, stat, why)
      if (stat /= MET_OK) return
      src%ncid = int(ncid, ik)

      st = nc_inq_dimlen_f(ncid, 'grid', dlen)
      call check_read(st, 'its grid dimension', src, stat, why)
      if (stat /= MET_OK) return
      src%ngrid = int(dlen, ik)
      !--- The grid column: grid_index as given, or the nearest to the [site] location.
      if (src%fcfg%grid_match == GRIDMATCH_NEAREST) then
         call resolve_grid_index(src, stat, why)             ! sets src%grid_index in 1..ngrid
         if (stat /= MET_OK) return
      else                                                   ! GRIDMATCH_EXPLICIT
         if (src%fcfg%grid_index < 1_ik .or. src%fcfg%grid_index > src%ngrid) then
            write(buf,'(a,i0,a,i0)') ' met_open: grid_index ', src%fcfg%grid_index, ' out of range 1..', src%ngrid
            call add_line(why, trim(buf))
            stat = MET_ERR_FILE
            return
         end if
      end if

      st = nc_inq_dimlen_f(ncid, 'time', dlen)
      call check_read(st, 'its time dimension', src, stat, why)
      if (stat /= MET_OK) return
      src%nrec = int(dlen, ik)
      if (src%nrec < 1_ik) then
         call add_line(why, ' met_open: the forcing file has no time records')
         stat = MET_ERR_FILE
         return
      end if

      !--- The base time, from the "seconds since <base>" units attribute.
      block
         integer(c_int) :: vid
         st = nc_inq_varid_f(ncid, 'time', vid)
         call check_read(st, 'its time variable', src, stat, why)
         if (stat /= MET_OK) return
         st = nc_get_att_text_f(ncid, vid, 'units', units)
         call check_read(st, 'time:units', src, stat, why)
         if (stat /= MET_OK) return
         call time_units_base(units, src%base_time, ok)
         if (.not. ok) then
            call add_line(why, ' met_open: could not parse time:units "seconds since <base>"')
            stat = MET_ERR_FILE
            return
         end if
      end block

      !--- The whole time axis (seconds since base_time).
      call read_time_axis(src, stat, why)
      if (stat /= MET_OK) return

      !---~---
      !   The wind vector: a file carrying u10 and v10 supplies it, and the speed comes from
      ! the components; a file with Wind alone gives the speed only.
      !---~---
      block
         integer(c_int) :: vu, vv
         src%has_wind_vector = nc_inq_varid_f(ncid, 'u10', vu) == NC_NOERR .and.                  &
                               nc_inq_varid_f(ncid, 'v10', vv) == NC_NOERR
      end block

      !--- The humidity the file carries: exactly one of RHair, Tdew and Qair.
      call detect_humidity(src, ncid, stat, why)
      if (stat /= MET_OK) return

      !--- The file's record spacing and global attributes against [forcing].
      call validate_file_against_config(src, ncid, stat, why)
      if (stat /= MET_OK) return

      !--- The declared recycle window against the file's record stamps, known only here.
      call validate_recycle_window(src, stat, why)
      if (stat /= MET_OK) return

      !---~---
      !   The records the run can read, into memory: the recycle window, or the run period (the
      ! whole file when neither is known, as in unit tests). Only those, because a file in 1 x 1
      ! chunks makes a whole-file read cost HDF5 metadata for every chunk.
      !---~---
      block
         integer(ik) :: r0, r1
         r0 = 1_ik
         r1 = src%nrec
         if (src%fcfg%recycle .and. src%n_cycle_years >= 1_ik) then
            r0 = max(1_ik, src%irec_cycle_first - 1_ik)
            r1 = min(src%nrec, src%irec_cycle_last + 1_ik)
         else if (present(run_start) .and. present(run_end)) then
            r0 = max(1_ik, first_record_after(src, seconds_between(src%base_time, run_start)) - 2_ik)
            r1 = min(src%nrec, first_record_after(src, seconds_between(src%base_time, run_end)) + 1_ik)
         end if
         if (r1 <= r0) r1 = min(src%nrec, r0 + 1_ik)
         call read_series(src, ncid, r0, r1, stat, why)
         if (stat /= MET_OK) return
         st = nc_close(ncid)
         src%ncid = -1_ik
         src%rec_first = r0                          ! a cursor's first bracket: the first records in range
      end block

      !---~---
      !   RHair is a fraction. Written as a percentage it would clip every record to saturation,
      ! a silent 100 % sky, so a value no fraction can reach stops the run here. A gap is left
      ! to the step that reads it, which names its record.
      !---~---
      if (src%humidity == HUMIDITY_RHAIR) then
         block
            integer(ik) :: j
            real(wp)    :: rh_max
            j = src%series_col(FLD_RHAIR)
            rh_max = maxval(src%series(:, j), mask=.not. ieee_is_nan(src%series(:, j)))
            if (rh_max > RH_FRACTION_MAX) then
               write(buf,'(a,f0.2,a)') ' met_open: RHair reaches ', rh_max,                                &
                                       ', but it is a fraction (units "1"), not a percentage.'
               call add_line(why, trim(buf))
               stat = MET_ERR_HUMIDITY
               return
            end if
         end block
      end if
   end subroutine default_file_open

   !---~---
   !   One record's values at the cell, as the file stores them: the humidity in its own form,
   ! and the shortwave total, or its four streams under sw_partition = "passthrough". The
   ! shared ingest (ed_met's read_record) converts, checks and lapses them.
   !---~---
   subroutine default_file_record(src, irec, rec, humidity_value, sw_total)
      type(met_source_t), intent(in)    :: src
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(inout) :: rec
      real(wp),           intent(out)   :: humidity_value, sw_total
      rec%air_temp   = series_value(src, FLD_TAIR,  irec)
      humidity_value = series_value(src, humidity_field(src%humidity), irec)
      rec%air_pressure = series_value(src, FLD_PSURF, irec)
      if (src%has_wind_vector) then
         rec%wind_u = series_value(src, FLD_U10, irec)
         rec%wind_v = series_value(src, FLD_V10, irec)
         rec%wind   = sqrt(rec%wind_u**2 + rec%wind_v**2)
      else
         rec%wind   = series_value(src, FLD_WIND,  irec)
      end if
      rec%rain    = series_value(src, FLD_RAINF, irec)              ! total rainfall rate [kg/m2/s]
      !---~---
      !   LWdown is optional under lwdown_source = "synthesize", which replaces it
      ! (docs/science/forcing.md §11).
      !---~---
      if (src%fcfg%lwdown_source == LW_SYNTHESIZE) then
         rec%lw_down = series_value(src, FLD_LWDOWN, irec, default=0.0_wp)
      else
         rec%lw_down = series_value(src, FLD_LWDOWN, irec)
      end if
      sw_total = 0.0_wp
      if (src%fcfg%sw_partition == SWPART_PASSTHROUGH) then
         rec%par_beam    = series_value(src, FLD_PAR_BEAM,    irec)
         rec%par_diffuse = series_value(src, FLD_PAR_DIFFUSE, irec)
         rec%nir_beam    = series_value(src, FLD_NIR_BEAM,    irec)
         rec%nir_diffuse = series_value(src, FLD_NIR_DIFFUSE, irec)
      else
         sw_total = series_value(src, FLD_SWDOWN, irec)
      end if
   end subroutine default_file_record

   !---~---
   !   Resolve src%grid_index by nearest great-circle distance from [site] lat/lon to the file's
   ! latitude(grid)/longitude(grid) coordinate vectors.
   !---~---
   subroutine resolve_grid_index(src, stat, why)
      type(met_source_t),            intent(inout) :: src
      integer(ik),                   intent(inout) :: stat
      character(len=:), allocatable, intent(inout) :: why
      integer(c_int)    :: st, vlat, vlon
      integer(c_size_t) :: start1(1), count1(1)
      real(wp), allocatable :: lat_grid(:), lon_grid(:)
      allocate(lat_grid(src%ngrid), lon_grid(src%ngrid))
      start1(1) = 0_c_size_t
      count1(1) = int(src%ngrid, c_size_t)
      st = nc_inq_varid_f(int(src%ncid, c_int), 'latitude', vlat)
      if (st == NC_NOERR) st = nc_get_vara_double(int(src%ncid, c_int), vlat, start1, count1, lat_grid)
      call check_read(st, 'its latitude(grid)', src, stat, why)
      if (stat /= MET_OK) return
      st = nc_inq_varid_f(int(src%ncid, c_int), 'longitude', vlon)
      if (st == NC_NOERR) st = nc_get_vara_double(int(src%ncid, c_int), vlon, start1, count1, lon_grid)
      call check_read(st, 'its longitude(grid)', src, stat, why)
      if (stat /= MET_OK) return
      src%grid_index = nearest_grid_index(src%fcfg%longitude, src%fcfg%latitude, lon_grid, lat_grid)
      write(*,'(a,i0,a,f8.3,a,f8.3,a)') 'met_open: grid_match=nearest resolved grid_index=', src%grid_index, &
            ' (site lon=', src%fcfg%longitude, ' lat=', src%fcfg%latitude, ')'
      deallocate(lat_grid, lon_grid)
   end subroutine resolve_grid_index

   !--- Helpers (NetCDF backend).
   !--- Read the whole time coordinate (seconds since base_time) into src%time_sec.
   subroutine read_time_axis(src, stat, why)
      type(met_source_t),            intent(inout) :: src
      integer(ik),                   intent(inout) :: stat
      character(len=:), allocatable, intent(inout) :: why
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start1(1), count1(1)
      if (allocated(src%time_sec)) deallocate(src%time_sec)
      allocate(src%time_sec(src%nrec))
      start1(1) = 0_c_size_t
      count1(1) = int(src%nrec, c_size_t)
      st = nc_inq_varid_f(int(src%ncid, c_int), 'time', vid)
      if (st == NC_NOERR) st = nc_get_vara_double(int(src%ncid, c_int), vid, start1, count1, src%time_sec)
      call check_read(st, 'its time values', src, stat, why)
   end subroutine read_time_axis

   !---~---
   !   A netCDF call of the open that failed: MET_ERR_FILE, and the reason, naming what was being
   ! read and the library's text.
   !---~---
   subroutine check_read(st, what, src, stat, why)
      integer(c_int),                intent(in)    :: st
      character(len=*),              intent(in)    :: what
      type(met_source_t),            intent(in)    :: src
      integer(ik),                   intent(inout) :: stat
      character(len=:), allocatable, intent(inout) :: why
      if (st == NC_NOERR) return
      call add_line(why, ' met_open: cannot read the forcing file '//trim(src%fcfg%path)//': '//what//        &
                         ' ('//trim(nc_error_text(st))//')')
      stat = MET_ERR_FILE
   end subroutine check_read

   !---~---
   !   The humidity form of an ED-MORE forcing file: exactly one of its three variables
   ! (docs/science/forcing.md §1). A file with none cannot drive the canopy air; a file with
   ! two leaves the reader to pick one, and the two may disagree through the provider's own
   ! saturation curve.
   !---~---
   subroutine detect_humidity(src, ncid, stat, why)
      type(met_source_t),            intent(inout) :: src
      integer(c_int),                intent(in)    :: ncid
      integer(ik),                   intent(out)   :: stat
      character(len=:), allocatable, intent(inout) :: why
      integer(ik), parameter :: FORMS(3) = [HUMIDITY_QAIR, HUMIDITY_RHAIR, HUMIDITY_TDEW]
      character(len=320) :: buf
      integer(c_int) :: vid
      integer(ik)    :: i, n
      stat = MET_OK
      n = 0_ik
      do i = 1_ik, size(FORMS, kind=ik)
         if (nc_inq_varid_f(ncid, trim(ED_FIELD(humidity_field(FORMS(i)))), vid) == NC_NOERR) then
            n = n + 1_ik
            src%humidity = FORMS(i)
         end if
      end do
      if (n /= 1_ik) then
         write(buf,'(a,i0,a)') ' met_open: the forcing file carries ', n,                              &
                               ' of the humidity variables RHair, Tdew and Qair; it must carry exactly one.'
         call add_line(why, trim(buf))
         call add_line(why, '   Store the humidity the source measured (RHair for a tower, Tdew for')
         call add_line(why, '   ERA5-Land); ED-MORE converts it with its own saturation curve.')
         stat = MET_ERR_HUMIDITY
      end if
   end subroutine detect_humidity

   !---~---
   !   An ED-MORE forcing file's fields at the cell, records r0..r1, read at open. Each
   ! field's column in the series is found here, once (series_col, 0 for a field the file
   ! lacks). A value equal to the variable's _FillValue or missing_value is a gap, NaN, which
   ! stops the step that reads it (docs/science/forcing.md §10); a packed variable is
   ! unpacked.
   !---~---
   subroutine read_series(src, ncid, r0, r1, stat, why)
      type(met_source_t),            intent(inout) :: src
      integer(c_int),                intent(in)    :: ncid
      integer(ik),                   intent(in)    :: r0, r1
      integer(ik),                   intent(inout) :: stat
      character(len=:), allocatable, intent(inout) :: why
      integer(c_int)    :: st, vid(N_ED_FIELD)
      integer(c_size_t) :: start2(2), count2(2)
      type(nc_packing_t) :: packing
      logical     :: present_(N_ED_FIELD)
      integer(ik) :: j, n
      do j = 1_ik, N_ED_FIELD
         present_(j) = nc_inq_varid_f(ncid, trim(ED_FIELD(j)), vid(j)) == NC_NOERR
      end do
      if (allocated(src%series)) deallocate(src%series)
      allocate(src%series(r0:r1, count(present_)))                   ! indexed by record number
      start2 = [int(r0 - 1_ik, c_size_t), int(src%grid_index - 1_ik, c_size_t)]    ! [time, grid]
      count2 = [int(r1 - r0 + 1_ik, c_size_t), 1_c_size_t]
      src%series_col = 0_ik
      n = 0_ik
      do j = 1_ik, N_ED_FIELD
         if (.not. present_(j)) cycle
         n = n + 1_ik
         src%series_col(j) = n
         if (n == 1_ik) call note_one_record_chunks(ncid, vid(j), r1 - r0 + 1_ik, src%fcfg%path)
         st = nc_get_vara_double(ncid, vid(j), start2, count2, src%series(:, n))
         call check_read(st, 'its '//trim(ED_FIELD(j)), src, stat, why)
         if (stat /= MET_OK) return
         packing = nc_packing_of(ncid, vid(j))
         src%series(:, n) = nc_unpacked(packing, src%series(:, n))
      end do
   end subroutine read_series

   !---~---
   !   A file stored one time record per chunk is read chunk by chunk at open: about 3 s for
   ! BCI's 90,528 records, most of a short run's start-up. Say so once, with the fix. The
   ! values are the same either way.
   !---~---
   subroutine note_one_record_chunks(ncid, varid, nrec, path)
      integer(c_int),   intent(in) :: ncid, varid
      integer(ik),      intent(in) :: nrec
      character(len=*), intent(in) :: path
      integer(c_int)    :: storage
      integer(c_size_t) :: chunks(8)
      if (nrec < 1000_ik) return
      if (nc_inq_var_chunking(ncid, varid, storage, chunks) /= NC_NOERR) return
      if (storage /= NC_CHUNKED .or. chunks(1) /= 1_c_size_t) return      ! chunks(1) is the time axis
      write(*,'(3a)') ' note: ', trim(path), ' is stored one time record per chunk, so reading it'
      write(*,'(a)')  '       at start-up is slow (about 3 s per 90,000 records). Rewriting it once with'
      write(*,'(a)')  '       "nccopy -c time/8760,grid/1 <file> <new file>" keeps every value and reads fast.'
   end subroutine note_one_record_chunks

   !---~---
   !   Field fld at record irec. A field the file does not carry stops the run, unless the
   ! caller gives a default (LWdown when it is synthesized).
   !---~---
   function series_value(src, fld, irec, default) result(val)
      type(met_source_t), intent(in) :: src
      integer(ik),        intent(in) :: fld, irec
      real(wp), optional, intent(in) :: default
      real(wp)    :: val
      integer(ik) :: j
      j = src%series_col(fld)
      if (j == 0_ik) then
         if (present(default)) then
            val = default
            return
         end if
         write(*,'(2a)') ' met_driver: the forcing file has no variable ', trim(ED_FIELD(fld))
         error stop 'met_driver: a required forcing variable is missing'
      end if
      call check_in_series(src, irec)
      val = src%series(irec, j)
   end function series_value

   !--- A record outside the range read at open was not prefetched: a programming error.
   subroutine check_in_series(src, irec)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: irec
      if (irec < lbound(src%series, 1) .or. irec > ubound(src%series, 1)) then
         write(*,'(a,i0,a,i0,a,i0)') ' met_driver: forcing record ', irec, ' outside the range read, ',  &
               lbound(src%series, 1), '..', ubound(src%series, 1)
         error stop 'met_driver: a forcing record was not read before the step (internal error)'
      end if
   end subroutine check_in_series

end module ed_met_default
