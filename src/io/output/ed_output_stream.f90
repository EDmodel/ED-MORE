! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_stream -- stage [5] of the diagnostic output, the netCDF writer: open a per-tier,
! per-time-chunk file (dims, coordinates, the tier's variables, CF metadata), append a record,
! and roll to a new file at a chunk boundary (write_record). Every file has a polygon axis
! after time; a record is packed over the polygons, with the fill value for a polygon without
! it (one that failed earlier in the month) and past a polygon's live cohort and patch counts.
! The cohort and patch axes grow to the largest live count in the file
! (docs/science/diagnostics.md §1, §3). output_serialize drains the queued records in the
! run's I/O phase, and output_close writes the final partial periods and closes the files. The
! only diagnostic module that calls netCDF.
!---~---
module ed_output_stream
   use iso_c_binding, only : c_int, c_size_t, c_double
   use ed_kinds,      only : wp, ik
   use ed_time,       only : ed_time_t, time_to_decimal_year, time_eq
   use ed_constants,     only : n_soil_layer_max
   use ed_soil_params,       only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use ed_netcdf_c,   only : nc_check, nc_close, nc_create_f, nc_def_dim_f, nc_enddef, nc_put_att_double_f,       &
                             nc_put_att_int_f, nc_put_var1_double, nc_put_vara_double, nc_put_vara_int,           &
                             nc_sync, NC_NOCLOBBER, NC_NOERR, NC_DOUBLE, NC_INT, NC_NETCDF4, NC_UNLIMITED
   use ed_netcdf_write,    only : define_var, put_text_att, put_global_text, define_polygon_axis,               &
                                  write_polygon_axis, define_site_statics, write_site_statics, OUTPUT_DEFLATE
   use ed_output_params,   only : FC_DAY, FC_MONTH, FC_YEAR, SYNC_FLUSH, freq_letter, TIER_FAST, N_FREQ
   use ed_output_integrate, only : close_tier, write_order
   use ed_output_types,    only : stream_file_t, pending_record_t, var_desc_t, slab_col, &
                                  output_files_t, output_buffers_t,                                  &
                                  DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT,            &
                                  DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE, diag_setup_t,               &
                                  XTYPE_DOUBLE, XTYPE_INT, AGG_MEAN, AGG_RATE, AGG_TMEAN,          &
                                  MISSING_VALUE, MISSING_INT
   implicit none
   private

   public :: output_serialize, output_close, output_file_path

   character(len=*), parameter :: TITLE = 'ED-MORE diagnostic aggregation output'

contains

   !---~---
   !   Write record i of tier t -- one closed period -- opening or rolling the tier's file as
   ! needed. `bufs` holds every polygon's buffers. Every polygon that holds record i holds the
   ! same period.
   !---~---
   subroutine write_record(files, bufs, t, i)
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(in)    :: bufs(:)
      integer(ik),            intent(in)    :: t, i
      integer(ik) :: bucket, np, p, p1, fc, n_cohort, n_patch
      np = size(bufs, kind=ik)
      if (np /= files%n_polygon)                                                                    &
         error stop 'write_record: one set of buffers per polygon of the file set'
      p1 = 1_ik
      do while (bufs(p1)%queue(t)%n < i)
         p1 = p1 + 1_ik
         if (p1 > np) error stop 'write_record: no polygon holds the record'
      end do
      associate (r1 => bufs(p1)%queue(t)%rec(i))
         do p = p1 + 1_ik, np
            if (bufs(p)%queue(t)%n < i) cycle
            if (.not. time_eq(bufs(p)%queue(t)%rec(i)%t_open, r1%t_open))                      &
               error stop 'write_record: the polygons'' records are not the same period'
         end do
         !---~---
         !   A file spans what its tier's file_chunk says. Its cohort and patch axes grow with the
         ! stand, which the calendar restructures at a month's turn (docs/science/diagnostics.md §3).
         !---~---
         fc = files%file_chunk(t)
         bucket = bucket_key(r1%t_open, fc)
         if (files%stream(t)%ncid < 0_ik .or. bucket /= files%stream(t)%chunk_bucket) then
            !--- The opening record's cohort and patch counts, the largest among the polygons: the chunk sizes.
            n_cohort = 0_ik
            n_patch = 0_ik
            do p = p1, np
               if (bufs(p)%queue(t)%n < i) cycle
               n_cohort = max(n_cohort, bufs(p)%queue(t)%rec(i)%n_cohort)
               n_patch  = max(n_patch,  bufs(p)%queue(t)%rec(i)%n_patch)
            end do
            call stream_close_file(files%stream(t))
            call open_file(files, t, r1, fc, bucket, n_cohort, n_patch)
         end if
      end associate
      call write_one(files, bufs, t, i, p1)
      !---~---
      !   No per-record nc_sync for the FAST tier: ~n_fast_per_slow records a day would each
      ! fsync. The chunk-boundary nc_close still flushes the file. Coarse tiers honour sync.
      !---~---
      if (files%sync_every == SYNC_FLUSH .and. t /= TIER_FAST)                                            &
         call nc_check(nc_sync(int(files%stream(t)%ncid, c_int)), 'nc_sync')
   end subroutine write_record

   subroutine stream_close_file(stream)
      type(stream_file_t), intent(inout) :: stream
      if (stream%ncid >= 0_ik) call nc_check(nc_close(int(stream%ncid, c_int)), 'nc_close stream')
      stream%ncid = -1_ik
      stream%nrec = 0_ik
      stream%chunk_bucket = -1_ik
      stream%cohort_dim = 0_ik
      stream%patch_dim = 0_ik
   end subroutine stream_close_file

   !--- Integer bucket key identifying the time-chunk a period belongs to.
   pure integer(ik) function bucket_key(t, file_chunk) result(b)
      type(ed_time_t),   intent(in) :: t
      integer(ik),       intent(in) :: file_chunk
      select case (file_chunk)
      case (FC_DAY)
         b = t%year*10000_ik + t%month*100_ik + t%day
      case (FC_MONTH)
         b = t%year*100_ik + t%month
      case (FC_YEAR)
         b = t%year
      case default
         b = 0_ik                    ! FC_ALL: single bucket
      end select
   end function bucket_key

   !--- Filename stamp for the time-chunk (empty for FC_ALL -> stampless -Y.nc).
   pure function chunk_stamp(t, file_chunk) result(s)
      type(ed_time_t),   intent(in) :: t
      integer(ik),       intent(in) :: file_chunk
      character(len=16) :: s
      s = ''
      select case (file_chunk)
      case (FC_DAY)
         write(s,'(i4.4,i2.2,i2.2)') t%year, t%month, t%day
      case (FC_MONTH)
         write(s,'(i4.4,i2.2)')      t%year, t%month
      case (FC_YEAR)
         write(s,'(i4.4)')           t%year
      end select
   end function chunk_stamp

   !---~---
   !   The path of tier `tier`'s file that holds a record stamped `t_open`:
   ! <dir>/<prefix>-<letter>[-<stamp>].nc, the stamp the record's day, month or year as the
   ! tier's file_chunk says, none for one file of all the run's records.
   !---~---
   pure function output_file_path(files, tier, t_open) result(path)
      type(output_files_t), intent(in) :: files
      integer(ik),          intent(in) :: tier
      type(ed_time_t),      intent(in) :: t_open
      character(len=512) :: path
      character(len=16)  :: stamp
      stamp = chunk_stamp(t_open, files%file_chunk(tier))
      if (len_trim(stamp) > 0) then
         path = trim(files%dir)//'/'//trim(files%prefix)//'-'//freq_letter(files%stream(tier)%freq)//'-'//  &
                trim(stamp)//'.nc'
      else
         path = trim(files%dir)//'/'//trim(files%prefix)//'-'//freq_letter(files%stream(tier)%freq)//'.nc'
      end if
   end function output_file_path

   pure function cell_methods_of(agg) result(cm)
      integer(ik), intent(in) :: agg
      character(len=16) :: cm
      select case (agg)
      case (AGG_MEAN, AGG_TMEAN, AGG_RATE)
         cm = 'time: mean'
      case default
         cm = 'time: point'   ! AGG_LAST
      end select
   end function cell_methods_of

   !---~---
   !   Create a file: dims, the polygon and site axes and their coordinates, the tier's
   ! variables and CF metadata. The cohort and patch axes are unlimited, so they grow to the
   ! largest live count of the file's records; the opening record's counts, the largest among
   ! a region's polygons, set their chunk sizes. A file that already exists is never
   ! overwritten: the run checked before it started (ed_region), and creating one stops it.
   !---~---
   subroutine open_file(files, tier, pr, file_chunk, bucket, n_cohort, n_patch)
      type(output_files_t),   intent(inout) :: files
      integer(ik),            intent(in)    :: tier, file_chunk, bucket
      type(pending_record_t), intent(in)    :: pr
      integer(ik),            intent(in)    :: n_cohort, n_patch   ! the opening record's live counts
      character(len=512) :: path
      integer(c_int)     :: ncid, dt, dpo, dsi, dc, dp, ds, dpf, dsz, dnv, vid, dims1(1), dims2(2), dims3(3), st
      integer(c_int)     :: dims4(4)
      integer(c_int)     :: v_row, v_col, v_nsite, v_chave, v_bnds
      integer(c_size_t)  :: st1(1), cn1(1)
      integer(ik)        :: j, k, np, cohort_dim, patch_dim
      logical            :: hasc, hasp, hass, haspf, hassz
      associate (stream => files%stream(tier), reg => files%reg, setup => files%diag)
      np   = files%n_polygon

      !---~---
      !   Which trailing dims does this tier need? A 2-D (soil layer x patch) variable needs
      ! both the patch and soil dims, which is why DIM_SOIL_PATCH sets two flags.
      !---~---
      hasc = .false.
      hasp = .false.
      hass = .false.
      haspf = .false.
      hassz = .false.
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         select case (reg%var(k)%dim)
         case (DIM_COHORT)
            hasc = .true.
         case (DIM_PATCH)
            hasp = .true.
         case (DIM_SOIL)
            hass = .true.
         case (DIM_PFT)
            haspf = .true.
         case (DIM_SIZE)
            hassz = .true.
         case (DIM_SOIL_PATCH)
            hasp = .true.
            hass = .true.
         end select
      end do

      !--- The chunk sizes along the growing axes; max(.,1) for bare ground.
      cohort_dim = max(n_cohort, 1_ik)
      patch_dim  = max(n_patch,  1_ik)

      path = output_file_path(files, tier, pr%t_open)

      st = nc_create_f(trim(path), ior(NC_NETCDF4, NC_NOCLOBBER), ncid)
      if (st /= NC_NOERR) then
         write(*,'(3a)') ' output: cannot create ', trim(path), ' (a file of that name exists, or the folder '// &
                         'cannot be written); an existing output file is never overwritten'
         call nc_check(st, 'stream nc_create')
      end if
      call nc_check(nc_def_dim_f(ncid, 'time', NC_UNLIMITED, dt), 'dim time')
      dpo = -1_c_int
      dc = -1_c_int
      dp = -1_c_int
      ds = -1_c_int
      dpf = -1_c_int
      dsz = -1_c_int
      call define_polygon_axis(ncid, np, dpo)
      call nc_check(nc_def_dim_f(ncid, 'site', int(files%site_max, c_size_t), dsi), 'dim site')
      if (hasc) call nc_check(nc_def_dim_f(ncid, 'cohort', NC_UNLIMITED, dc), 'dim cohort')
      if (hasp) call nc_check(nc_def_dim_f(ncid, 'patch',  NC_UNLIMITED, dp), 'dim patch')
      if (hass) call nc_check(nc_def_dim_f(ncid, 'soil',   int(n_soil_layer_max, c_size_t), ds), 'dim soil')
      if (haspf) call nc_check(nc_def_dim_f(ncid, 'pft',   int(max(setup%n_pft,1_ik), c_size_t), dpf), 'dim pft')
      if (hassz) call nc_check(nc_def_dim_f(ncid, 'dbh_class',                                       &
                               int(max(setup%n_dbh_class,1_ik), c_size_t), dsz), 'dim dbh_class')

      !--- time coordinate + calendar companions (period-start stamp).
      stream%v_time  = int(def_scalar_var(ncid, dt, 'time',  NC_DOUBLE, 'year', 'decimal calendar year (period start)'), ik)
      stream%v_year  = int(def_scalar_var(ncid, dt, 'year',  NC_INT,    '1', 'calendar year (period start)'), ik)
      stream%v_month = int(def_scalar_var(ncid, dt, 'month', NC_INT,    '1', 'calendar month (period start)'), ik)
      stream%v_day   = int(def_scalar_var(ncid, dt, 'day',   NC_INT,    '1', 'calendar day (period start)'), ik)
      !---~---
      !   The FAST tier: clock companions, so a reader can group by hour without decoding the
      ! decimal `time`; the period-start stamp, as the calendar companions above.
      !---~---
      stream%v_hour = -1_ik
      stream%v_minute = -1_ik
      stream%v_second = -1_ik
      if (tier == TIER_FAST) then
         stream%v_hour   = int(def_scalar_var(ncid, dt, 'hour',   NC_INT, '1', 'clock hour (period start)'), ik)
         stream%v_minute = int(def_scalar_var(ncid, dt, 'minute', NC_INT, '1', 'clock minute (period start)'), ik)
         stream%v_second = int(def_scalar_var(ncid, dt, 'second', NC_INT, '1', 'clock second (period start)'), ik)
      end if
      !--- The live counts, per record and polygon.
      if (hasc .or. hasp) then
         dims2 = [dt, dpo]
         if (hasc) then
            call define_var(ncid, 'n_cohort', NC_INT, dims2, vid,                                              &
                            long_name='the polygon''s live cohorts this record (the rest is fill)')
            call nc_check(nc_put_att_int_f(ncid, vid, '_FillValue', NC_INT, int(MISSING_INT, c_int)), 'fill')
            stream%v_ncohort = int(vid, ik)
         end if
         if (hasp) then
            call define_var(ncid, 'n_patch', NC_INT, dims2, vid,                                               &
                            long_name='the polygon''s live patches this record (the rest is fill)')
            call nc_check(nc_put_att_int_f(ncid, vid, '_FillValue', NC_INT, int(MISSING_INT, c_int)), 'fill')
            stream%v_npatch = int(vid, ik)
         end if
      end if

      !---~---
      !   The polygon axis (define_polygon_axis: its id and cell centre) and its grid indices, so a
      ! region maps back onto the lat/lon grid in one step.
      !---~---
      dims1 = [dpo]
      call define_var(ncid, 'row', NC_INT, dims1, v_row, long_name='0-based latitude index on the forcing grid')
      call define_var(ncid, 'col', NC_INT, dims1, v_col, long_name='0-based longitude index on the forcing grid')
      call define_var(ncid, 'n_site', NC_INT, dims1, v_nsite,                                         &
                      long_name='the polygon''s sites (the rest of the site axis is fill)')
      !--- Each polygon's sites, as it runs them: area, terrain and soil column.
      call define_site_statics(ncid, dpo, dsi, [SOIL_BC_DRAIN, SOIL_BC_AQUIFER], 'drain aquifer', by_row=.false.)
      if (allocated(files%polygon_chave_e)) then
         call define_var(ncid, 'chave_e_py', NC_DOUBLE, dims1, v_chave, units='1',                       &
                         long_name='Chave et al. (2014) environmental stress index E')
      end if

      !---~---
      !   Axis coordinates: the pft and dbh_class axes have run-dependent lengths, so each file
      ! carries its own coordinates and stands on its own (docs/science/diagnostics.md §3).
      !---~---
      stream%v_pft = -1_ik
      stream%v_dbh_lower = -1_ik
      stream%v_dbh_upper = -1_ik
      stream%v_soil_z = -1_ik
      if (haspf) then
         dims1 = [dpf]
         call define_var(ncid, 'pft', NC_INT, dims1, vid, long_name='plant functional type index')
         stream%v_pft = int(vid, ik)
      end if
      !---~---
      !   Each site has its own layers: the coordinate is per site, with each layer's top and
      ! bottom as its CF bounds, and the fill value below the site's bottom.
      !---~---
      if (hass) then
         call nc_check(nc_def_dim_f(ncid, 'nv', 2_c_size_t, dnv), 'dim nv')
         dims3 = [dpo, dsi, ds]
         call define_var(ncid, 'soil_z', NC_DOUBLE, dims3, vid, units='m',                             &
                         long_name='the site''s soil layer node depths (negative downward)')
         call put_text_att(ncid, vid, 'bounds', 'soil_z_bnds')
         call nc_check(nc_put_att_double_f(ncid, vid, '_FillValue', NC_DOUBLE, real(MISSING_VALUE, c_double)), 'fill')
         stream%v_soil_z = int(vid, ik)
         dims4 = [dpo, dsi, ds, dnv]
         call define_var(ncid, 'soil_z_bnds', NC_DOUBLE, dims4, v_bnds, units='m',                     &
                         long_name='the site''s soil layer tops and bottoms (negative downward)')
         call nc_check(nc_put_att_double_f(ncid, v_bnds, '_FillValue', NC_DOUBLE, real(MISSING_VALUE, c_double)), &
                       'fill')
      end if
      if (hassz) then
         dims1 = [dsz]
         call define_var(ncid, 'dbh_lower', NC_DOUBLE, dims1, vid, units='cm', long_name='DBH class lower edge (inclusive)')
         stream%v_dbh_lower = int(vid, ik)
         call define_var(ncid, 'dbh_upper', NC_DOUBLE, dims1, vid, units='cm',                         &
                         long_name='DBH class upper edge (exclusive; the last class also holds every larger plant)')
         stream%v_dbh_upper = int(vid, ik)
      end if

      !--- registry-driven variable definitions (dims per DIM_*, chunk+deflate + CF attrs).
      stream%vid = -1_ik
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         vid = def_registry_var(ncid, dt, dpo, dsi, dc, dp, ds, dpf, dsz, reg%var(k), np, files%site_max,  &
                                cohort_dim, patch_dim, setup)
         stream%vid(k) = int(vid, ik)
      end do

      call put_global_text(ncid, 'title', TITLE)
      call put_global_text(ncid, 'Conventions', 'CF-1.10')
      if (len_trim(files%forcing_qair) > 0) call put_global_text(ncid, 'forcing_qair', trim(files%forcing_qair))
      call nc_check(nc_enddef(ncid), 'stream enddef')

      !--- Write the axis coordinates once, right after enddef (they do not vary by record).
      st1 = [0_c_size_t]
      cn1 = [int(np, c_size_t)]
      call write_polygon_axis(ncid, files%polygon_id, files%polygon_lat, files%polygon_lon)
      call nc_check(nc_put_vara_int(ncid, v_row, st1, cn1, int(files%polygon_row, c_int)), 'put row')
      call nc_check(nc_put_vara_int(ncid, v_col, st1, cn1, int(files%polygon_col, c_int)), 'put col')
      call nc_check(nc_put_vara_int(ncid, v_nsite, st1, cn1, int(files%polygon_n_site, c_int)), 'put n_site')
      call write_site_statics(ncid, 0_ik, files%site_area, files%site_elevation, files%site_slope,              &
                              files%site_aspect, files%site_n_soil, files%site_soil_depth, files%site_bottom_bc, &
                              files%site_drain_factor, files%site_deep_temp)
      if (allocated(files%polygon_chave_e))                                                       &
         call nc_check(nc_put_vara_double(ncid, v_chave, st1, cn1, real(files%polygon_chave_e, c_double)), &
                       'put chave_e_py')
      if (hass) call write_site_soil_coord(ncid, int(stream%v_soil_z, c_int), v_bnds, files)
      if (haspf) call write_pft_coord(ncid, stream%v_pft, setup%n_pft)
      if (hassz) call write_size_coord(ncid, stream%v_dbh_lower, stream%v_dbh_upper, setup)

      stream%ncid = int(ncid, ik)
      stream%nrec = 0_ik
      stream%chunk_bucket = bucket
      stream%has_cohort = hasc
      stream%has_patch = hasp
      stream%has_soil = hass
      stream%has_pft = haspf
      stream%has_size = hassz
      stream%cohort_dim = 0_ik
      stream%patch_dim = 0_ik
      stream%d_time = int(dt, ik)
      stream%d_polygon = int(dpo, ik)
      stream%d_cohort = int(dc, ik)
      stream%d_patch = int(dp, ik)
      stream%d_soil = int(ds, ik)
      stream%d_pft = int(dpf, ik)
      stream%d_size = int(dsz, ik)
      stream%d_site = int(dsi, ik)
      write(*,'(2a)') ' output: ', trim(path)
      end associate
   end subroutine open_file

   !---~---
   !   The soil coordinate: each site's own layer node depths and, as their bounds, each
   ! layer's top and bottom [m], negative downward, so its bottom layer ends at its own soil
   ! depth. Below a site's bottom layer, and for a polygon's missing sites, both are the fill
   ! value, as its data are.
   !---~---
   subroutine write_site_soil_coord(ncid, v_z, v_bnds, files)
      integer(c_int),       intent(in) :: ncid, v_z, v_bnds
      type(output_files_t), intent(in) :: files
      real(c_double), allocatable :: z(:,:,:), b(:,:,:,:)
      integer(c_size_t) :: st3(3), cn3(3), st4(4), cn4(4)
      integer(ik)       :: p, is, k, ns
      ns = files%site_max
      allocate(z(n_soil_layer_max, ns, files%n_polygon), b(2, n_soil_layer_max, ns, files%n_polygon))
      z = real(MISSING_VALUE, c_double)
      b = real(MISSING_VALUE, c_double)
      do p = 1_ik, files%n_polygon
         do is = 1_ik, files%polygon_n_site(p)
            do k = 1_ik, files%site_n_soil(is, p)
               z(k, is, p)    = real(files%site_soil_z(k, is, p), c_double)
               b(1, k, is, p) = real(files%site_soil_face(k, is, p), c_double)
               b(2, k, is, p) = real(files%site_soil_face(k + 1_ik, is, p), c_double)
            end do
         end do
      end do
      st3 = 0_c_size_t
      cn3 = [int(files%n_polygon, c_size_t), int(ns, c_size_t), int(n_soil_layer_max, c_size_t)]
      call nc_check(nc_put_vara_double(ncid, v_z, st3, cn3, z), 'put soil_z')
      st4 = 0_c_size_t
      cn4 = [cn3(1), cn3(2), cn3(3), 2_c_size_t]
      call nc_check(nc_put_vara_double(ncid, v_bnds, st4, cn4, b), 'put soil_z_bnds')
   end subroutine write_site_soil_coord

   !--- The pft coordinate: the 1-based PFT indices this run carries.
   subroutine write_pft_coord(ncid, vid_ik, n_pft)
      integer(c_int), intent(in) :: ncid
      integer(ik),    intent(in) :: vid_ik, n_pft
      integer(c_int)    :: iarr(max(n_pft,1_ik))
      integer(c_size_t) :: st(1), cn(1)
      integer(ik)       :: i
      do i = 1_ik, max(n_pft, 1_ik)
         iarr(i) = int(i, c_int)
      end do
      st = [0_c_size_t]
      cn = [int(max(n_pft,1_ik), c_size_t)]
      call nc_check(nc_put_vara_int(ncid, int(vid_ik, c_int), st, cn, iarr), 'put pft coord')
   end subroutine write_pft_coord

   !---~---
   !   The dbh_class coordinates: the lower and upper edge of each class, so a reader can
   ! label the axis without knowing the run's TOML.
   !---~---
   subroutine write_size_coord(ncid, vlo_ik, vhi_ik, setup)
      integer(c_int),      intent(in) :: ncid
      integer(ik),         intent(in) :: vlo_ik, vhi_ik
      type(diag_setup_t),  intent(in) :: setup
      real(c_double)    :: lo(max(setup%n_dbh_class,1_ik)), hi(max(setup%n_dbh_class,1_ik))
      integer(c_size_t) :: st(1), cn(1)
      integer(ik)       :: i, nc
      nc = max(setup%n_dbh_class, 1_ik)
      do i = 1_ik, nc
         lo(i) = real(setup%dbh_edges(i),        c_double)
         hi(i) = real(setup%dbh_edges(i + 1_ik), c_double)
      end do
      st = [0_c_size_t]
      cn = [int(nc, c_size_t)]
      call nc_check(nc_put_vara_double(ncid, int(vlo_ik, c_int), st, cn, lo), 'put dbh_lower')
      call nc_check(nc_put_vara_double(ncid, int(vhi_ik, c_int), st, cn, hi), 'put dbh_upper')
   end subroutine write_size_coord

   !--- Define a 1-D record variable (time) + units/long_name; returns its varid.
   integer(c_int) function def_scalar_var(ncid, dt, name, xtype, units, lname) result(vid)
      integer(c_int),   intent(in) :: ncid, dt, xtype
      character(len=*), intent(in) :: name, units, lname
      call define_var(ncid, name, xtype, [dt], vid, units=units, long_name=lname)
   end function def_scalar_var

   !---~---
   !   Define one registry variable: (time, polygon[, axis...]), chunked one record at a time
   ! and deflated, with CF attributes and _FillValue. A soil-by-patch variable is (time,
   ! polygon, patch, soil) and a site soil column (time, polygon, site, soil): the slab is
   ! patch- or site-major and strided by n_soil_layer_max, netCDF's row-major order.
   !---~---
   integer(c_int) function def_registry_var(ncid, dt, dpo, dsi, dc, dp, ds, dpf, dsz, v, np, ns,    &
                                            cohort_dim, patch_dim, setup) result(vid)
      integer(c_int),      intent(in) :: ncid, dt, dpo, dsi, dc, dp, ds, dpf, dsz
      type(var_desc_t),    intent(in) :: v
      integer(ik),         intent(in) :: np                     ! polygons
      integer(ik),         intent(in) :: ns                     ! the site axis's length
      integer(ik),         intent(in) :: cohort_dim, patch_dim  ! the chunk sizes along the growing axes
      type(diag_setup_t),  intent(in) :: setup
      integer(c_int)    :: xt, dims(4), rank
      integer(c_size_t) :: chunk(4)
      character(len=16) :: cm
      xt = merge(NC_INT, NC_DOUBLE, v%xtype == XTYPE_INT)
      rank = 2_c_int
      dims(1:2)  = [dt, dpo]
      chunk(1:2) = [1_c_size_t, int(np, c_size_t)]
      select case (v%dim)
      case (DIM_COHORT)
         rank = rank + 1_c_int
         dims(rank) = dc
         chunk(rank) = int(cohort_dim, c_size_t)
      case (DIM_PATCH)
         rank = rank + 1_c_int
         dims(rank) = dp
         chunk(rank) = int(patch_dim, c_size_t)
      case (DIM_SOIL)
         rank = rank + 1_c_int
         dims(rank) = dsi
         chunk(rank) = int(ns, c_size_t)
         rank = rank + 1_c_int
         dims(rank) = ds
         chunk(rank) = int(n_soil_layer_max, c_size_t)
      case (DIM_PFT)
         rank = rank + 1_c_int
         dims(rank) = dpf
         chunk(rank) = int(max(setup%n_pft, 1_ik), c_size_t)
      case (DIM_SIZE)
         rank = rank + 1_c_int
         dims(rank) = dsz
         chunk(rank) = int(max(setup%n_dbh_class, 1_ik), c_size_t)
      case (DIM_SOIL_PATCH)
         rank = rank + 1_c_int
         dims(rank) = dp
         chunk(rank) = int(patch_dim, c_size_t)
         rank = rank + 1_c_int
         dims(rank) = ds
         chunk(rank) = int(n_soil_layer_max, c_size_t)
      case (DIM_SITE)
         rank = rank + 1_c_int
         dims(rank) = dsi
         chunk(rank) = int(ns, c_size_t)
      end select
      call define_var(ncid, trim(v%name), xt, dims(1:rank), vid, units=trim(v%units), long_name=trim(v%long_name), &
                      chunk=chunk(1:rank), deflate=OUTPUT_DEFLATE)
      cm = cell_methods_of(v%agg)
      call put_text_att(ncid, vid, 'cell_methods', trim(cm))
      if (v%xtype == XTYPE_INT) then
         call nc_check(nc_put_att_int_f(ncid, vid, '_FillValue', NC_INT, int(MISSING_INT, c_int)), 'fill')
      else
         call nc_check(nc_put_att_double_f(ncid, vid, '_FillValue', NC_DOUBLE, real(MISSING_VALUE, c_double)), 'fill')
      end if
   end function def_registry_var

   !---~---
   !   Append record i: the calendar and the live counts from polygon p1 (the first that holds
   ! it), then each variable, packed over the polygons into one array the length of the file's
   ! axis, the fill value past a polygon's live length and for a polygon without the record.
   !---~---
   subroutine write_one(files, bufs, t, i, p1)
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(in)    :: bufs(:)
      integer(ik),            intent(in)    :: t, i, p1
      integer(c_int)    :: ncid, rank
      integer(c_size_t) :: t0, i1(1), st(4), cn(4)
      integer(ik)       :: j, k, np, p, n, ns, c, m, n_cohort, n_patch
      real(c_double),  allocatable :: x(:,:)
      integer(c_int),  allocatable :: kx(:,:)
      associate (stream => files%stream(t), reg => files%reg, r1 => bufs(p1)%queue(t)%rec(i))
      ncid = int(stream%ncid, c_int)
      t0 = int(stream%nrec, c_size_t)
      i1 = [t0]
      np  = size(bufs, kind=ik)
      !---~---
      !   This record's cohort and patch extent: the largest live count among the polygons. Past
      ! it, and past a polygon's own count, the growing axes read as the fill value.
      !---~---
      n_cohort = 0_ik
      n_patch = 0_ik
      do p = 1_ik, np
         if (bufs(p)%queue(t)%n < i) cycle
         n_cohort = max(n_cohort, bufs(p)%queue(t)%rec(i)%n_cohort)
         n_patch  = max(n_patch,  bufs(p)%queue(t)%rec(i)%n_patch)
      end do
      stream%cohort_dim = max(stream%cohort_dim, n_cohort)
      stream%patch_dim  = max(stream%patch_dim,  n_patch)
      !--- calendar (period start) + live counts.
      call nc_check(nc_put_var1_double(ncid, int(stream%v_time, c_int), i1, time_to_decimal_year(r1%t_open)), 'put time')
      call put_int_rec(ncid, stream%v_year,  t0, int(r1%t_open%year,  c_int))
      call put_int_rec(ncid, stream%v_month, t0, int(r1%t_open%month, c_int))
      call put_int_rec(ncid, stream%v_day,   t0, int(r1%t_open%day,   c_int))
      if (stream%v_hour   >= 0_ik) call put_int_rec(ncid, stream%v_hour,   t0, int(r1%t_open%hour,   c_int))
      if (stream%v_minute >= 0_ik) call put_int_rec(ncid, stream%v_minute, t0, int(r1%t_open%minute, c_int))
      if (stream%v_second >= 0_ik) call put_int_rec(ncid, stream%v_second, t0, int(r1%t_open%second, c_int))
      if (stream%has_cohort) call put_count_rec(stream%v_ncohort, .true.)
      if (stream%has_patch)  call put_count_rec(stream%v_npatch,  .false.)

      do j = 1_ik, reg%nidx(t)
         k = reg%idx_freq(j, t)
         !--- The variable's axis length in this file, and the hyperslab: (time, polygon[, axis]).
         select case (reg%var(k)%dim)
         case (DIM_SCALAR)
            n = 1_ik
         case (DIM_COHORT)
            n = n_cohort
         case (DIM_PATCH)
            n = n_patch
         case (DIM_SOIL)
            n = files%site_max * n_soil_layer_max
         case (DIM_PFT)
            n = max(files%diag%n_pft, 1_ik)
         case (DIM_SIZE)
            n = max(files%diag%n_dbh_class, 1_ik)
         case (DIM_SOIL_PATCH)
            n = n_patch * n_soil_layer_max
         case (DIM_SITE)
            n = files%site_max
         case default
            n = 1_ik
         end select
         if (n <= 0_ik) cycle                       ! no cohort or patch this record: all fill
         rank = 2_c_int
         st(1:2) = [t0, 0_c_size_t]
         cn(1:2) = [1_c_size_t, int(np, c_size_t)]
         if (reg%var(k)%dim == DIM_SOIL_PATCH) then
            rank = rank + 2_c_int
            st(rank-1:rank) = 0_c_size_t
            cn(rank-1:rank) = [int(n_patch, c_size_t), int(n_soil_layer_max, c_size_t)]
         else if (reg%var(k)%dim == DIM_SOIL) then
            rank = rank + 2_c_int
            st(rank-1:rank) = 0_c_size_t
            cn(rank-1:rank) = [int(files%site_max, c_size_t), int(n_soil_layer_max, c_size_t)]
         else if (reg%var(k)%dim /= DIM_SCALAR) then
            rank = rank + 1_c_int
            st(rank) = 0_c_size_t
            cn(rank) = int(n, c_size_t)
         end if
         if (reg%var(k)%xtype == XTYPE_INT) then
            allocate(kx(n, np))
            kx = int(MISSING_INT, c_int)
            do p = 1_ik, np
               if (bufs(p)%queue(t)%n < i) cycle
               associate (r => bufs(p)%queue(t)%rec(i))
                  if (reg%var(k)%dim == DIM_SCALAR) then
                     kx(1, p) = real_to_int(r%sval(k), r%svalid(k))
                  else
                     ns = min(record_length(r, reg%var(k)%dim, k), n)
                     c = slab_col(r, k)
                     do m = 1_ik, ns
                        kx(m, p) = real_to_int(r%slab(m, c), r%slabvalid(m, c))
                     end do
                  end if
               end associate
            end do
            call nc_check(nc_put_vara_int(ncid, int(stream%vid(k), c_int), st(1:rank), cn(1:rank), kx),   &
                          'put '//trim(reg%var(k)%name))
            deallocate(kx)
         else
            allocate(x(n, np))
            x = real(MISSING_VALUE, c_double)
            do p = 1_ik, np
               if (bufs(p)%queue(t)%n < i) cycle
               associate (r => bufs(p)%queue(t)%rec(i))
                  if (reg%var(k)%dim == DIM_SCALAR) then
                     x(1, p) = r%sval(k)
                  else
                     ns = min(record_length(r, reg%var(k)%dim, k), n)
                     c = slab_col(r, k)
                     if (ns > 0_ik) x(1:ns, p) = r%slab(1:ns, c)
                  end if
               end associate
            end do
            call nc_check(nc_put_vara_double(ncid, int(stream%vid(k), c_int), st(1:rank), cn(1:rank), x),  &
                          'put '//trim(reg%var(k)%name))
            deallocate(x)
         end if
      end do
      stream%nrec = stream%nrec + 1_ik
      end associate

   contains
      !--- A record's live counts, one per polygon; the fill value for a polygon without it.
      subroutine put_count_rec(vid_ik, cohorts)
         integer(ik), intent(in) :: vid_ik
         logical,     intent(in) :: cohorts
         integer(c_int)    :: counts(np)
         integer(c_size_t) :: st2(2), cn2(2)
         integer(ik)       :: q
         counts = int(MISSING_INT, c_int)
         do q = 1_ik, np
            if (bufs(q)%queue(t)%n < i) cycle
            if (cohorts) then
               counts(q) = int(bufs(q)%queue(t)%rec(i)%n_cohort, c_int)
            else
               counts(q) = int(bufs(q)%queue(t)%rec(i)%n_patch, c_int)
            end if
         end do
         st2 = [t0, 0_c_size_t]
         cn2 = [1_c_size_t, int(np, c_size_t)]
         call nc_check(nc_put_vara_int(ncid, int(vid_ik, c_int), st2, cn2, counts), 'put live counts')
      end subroutine put_count_rec
   end subroutine write_one

   !---~---
   !   How many slab entries of variable k a record holds: its live length, or for the soil-by-
   ! patch slab (patch-major, strided by n_soil_layer_max) every layer of its live patches.
   !---~---
   pure integer(ik) function record_length(r, dim, k) result(ns)
      type(pending_record_t), intent(in) :: r
      integer(ik),            intent(in) :: dim, k
      if (dim == DIM_SOIL_PATCH) then
         ns = r%n_patch * n_soil_layer_max
      else
         ns = r%nslab(k)
      end if
   end function record_length


   !--- Map a normalized real (with validity) to an int record value (MISSING_INT if invalid).
   pure integer(c_int) function real_to_int(x, valid) result(iv)
      real(wp), intent(in) :: x
      logical,  intent(in) :: valid
      if (valid) then
         iv = int(nint(x, ik), c_int)
      else
         iv = int(MISSING_INT, c_int)
      end if
   end function real_to_int

   subroutine put_int_rec(ncid, vid_ik, t0, val)
      integer(c_int),    intent(in) :: ncid
      integer(ik),       intent(in) :: vid_ik
      integer(c_size_t), intent(in) :: t0
      integer(c_int),    intent(in) :: val
      integer(c_size_t) :: st(1), cn(1)
      integer(c_int)    :: vv(1)
      st = [t0]
      cn = [1_c_size_t]
      vv = [val]
      call nc_check(nc_put_vara_int(ncid, int(vid_ik, c_int), st, cn, vv), 'put int rec')
   end subroutine put_int_rec


   !---~---
   !   Drain every polygon's queued records to the per-tier files, in closing order, then empty
   ! the queues; called from the I/O phase at month boundaries, after a failed step and at the
   ! end of the run, so the queues hold one month.
   ! Record i of tier t is one period across the polygons. A polygon that failed earlier in
   ! the month closed fewer, so the longest queue sets the count and the others get the fill
   ! value: a failure loses no other polygon's month.
   !---~---
   subroutine output_serialize(files, bufs)
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      integer(ik) :: t, i, p, nrec
      if (.not. files%enabled) return
      call write_order(files, bufs)                 ! each cohort record in its month's end order
      do t = 1_ik, N_FREQ
         nrec = 0_ik
         do p = 1_ik, size(bufs, kind=ik)
            nrec = max(nrec, bufs(p)%queue(t)%n)
         end do
         do i = 1_ik, nrec
            call write_record(files, bufs, t, i)
         end do
         do p = 1_ik, size(bufs, kind=ik)
            bufs(p)%queue(t)%n = 0_ik
         end do
      end do
   end subroutine output_serialize

   !---~---
   !   End of run: optionally close every polygon's final partial windows (the FAST tier's too)
   ! and write them, then close the files.
   !---~---
   subroutine output_close(files, bufs, flush_partial)
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      logical, optional,      intent(in)    :: flush_partial
      logical     :: fp
      integer(ik) :: t, p
      if (.not. files%enabled) return
      fp = .true.
      if (present(flush_partial)) fp = flush_partial
      if (fp) then
         do p = 1_ik, size(bufs, kind=ik)
            do t = 1_ik, N_FREQ
               if (bufs(p)%has_data(t)) call close_tier(files, bufs(p), t)
            end do
         end do
         call output_serialize(files, bufs)
      end if
      do t = 1_ik, N_FREQ
         call stream_close_file(files%stream(t))
      end do
   end subroutine output_close

end module ed_output_stream
