! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_init_restart -- the state (restart) file: every polygon's prognostic state at an
! instant, in a netCDF-4 checkpoint a resumed run continues from exactly
! (docs/science/order_of_processes.md §6, The restart file). restart_file_open, _write and
! _close write a file a polygon's row at a time; read_restart_file reads a row, and
! restart_file_slots checks a file before a run starts from it. A file holds every
! variable, and the reader requires each one. Each level's variables are defined, written
! and read in one place, in its state module's section order.
! Diagnostic output is ed_output_*.
!---~---
module ed_init_restart
   use iso_c_binding,       only : c_int, c_size_t
   use ed_kinds,            only : wp, ik
   use ed_biogeochem_params, only : N_SOIL_POOL, SOIL_POOL_NAME
   use ed_constants,        only : n_soil_layer_max, n_snow_layer_max
   use ed_state_defaults,   only : is_unset_positive
   use ed_config,           only : ed_config_t, growth_window_steps
   use ed_time,             only : ed_time_t, time_to_stamp, time_to_decimal_year
   use ed_soil_params,      only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use ed_netcdf_c,         only : nc_check, nc_close, nc_create_f, nc_def_dim_f, nc_enddef, nc_get_att_text_f,   &
                                   nc_get_vara_double, nc_get_vara_int, nc_inq_dimlen, nc_inq_dimlen_f,           &
                                   nc_inq_varid_f, nc_open_f, nc_put_vara_double, nc_put_vara_int, NC_CLOBBER,    &
                                   NC_DOUBLE, NC_FILL_DOUBLE, NC_GLOBAL, NC_INT, NC_NETCDF4, NC_NOERR, NC_NOWRITE
   use ed_netcdf_write,     only : define_var, put_global_text, define_polygon_axis, write_polygon_axis,             &
                                   define_site_statics, write_site_statics, RESTART_DEFLATE
   use ed_polygon_state_vars,           only : polygon_t, alloc_polygon, rebuild_patch_map, rebuild_cohort_map
   use ed_site_state_vars,              only : make_site_surface
   use ed_cohort_state_vars,            only : set_cohort_size_from_carbon
   use ed_slow_update, only : sort_cohorts
   implicit none
   private

   public :: restart_file_t, restart_file_open, restart_file_write, restart_file_close
   public :: read_restart_file, restart_file_slots

   character(len=*), parameter :: STATE_TITLE = 'ED-MORE instantaneous state (restart) file'
   !--- The layout written and read, a global attribute of every file; a change to it counts up.
   integer,          parameter :: N_META_INT = 11, N_META_REAL = 2
   character(len=*), parameter :: SURF_WATER_UNITS  = 'kg/plant'   ! the canopy films' units attribute

   !--- An open state file, its polygon axis ahead of every other.
   type :: restart_file_t
      integer(c_int)     :: ncid   = -1_c_int
      character(len=512) :: path   = ''
   end type restart_file_t

   !--- The dimensions of a file being defined.
   type :: file_dims_t
      integer(c_int) :: ncid = -1_c_int
      integer(c_int) :: polygon, cohort, patch, site, pft, meta_int, meta_real, soil_layer, snow_layer
      integer(c_int) :: growth_window
   end type file_dims_t

   !--- One polygon's row (1-based) of an open file.
   type :: file_row_t
      integer(c_int) :: ncid = -1_c_int
      integer(ik)    :: row  = 1_ik
   end type file_row_t

contains

   !---~---
   !   0. THE FILE: open, write, close and read, the checks before a run starts from one,
   !    and one variable defined, written or read at a polygon's row
   !---~---

   !---~---
   !   Create <dir>/<prefix>-S-YYYYMMDDHHMMSS.nc and define every dimension and variable.
   ! The cohort and patch axes are n_cohort and n_patch long (at least 1), the largest of
   ! the polygons'. A `polygon` axis leads every variable, `polygon_id` names its rows,
   ! and the variables are compressed, since most polygons fill only part of the padded
   ! axes. Why each group persists: docs/science/order_of_processes.md §6, The restart
   ! file.
   !---~---
   subroutine restart_file_open(sf, dir, prefix, now, n_pft, n_window, n_cohort, n_patch, n_site,     &
                              polygon_id, polygon_lat, polygon_lon, restructure_pending, restructure_new_year)
      type(restart_file_t),  intent(out) :: sf
      character(len=*),      intent(in)  :: dir, prefix
      type(ed_time_t),       intent(in)  :: now
      integer(ik),           intent(in)  :: n_pft, n_window, n_cohort, n_patch, n_site
      integer(ik),           intent(in)  :: polygon_id(:)                  ! the polygons' ids ...
      real(wp),              intent(in)  :: polygon_lat(:), polygon_lon(:)  ! ... and locations [deg]
      logical,     optional, intent(in)  :: restructure_pending, restructure_new_year
      type(file_dims_t) :: d
      character(len=5)  :: pending
      integer(c_int)    :: ncid

      !--- The stamp is the simulated date: <prefix>-S-YYYYMMDDHHMMSS.nc.
      sf%path = trim(dir)//'/'//trim(prefix)//'-S-'//time_to_stamp(now)//'.nc'
      call nc_check(nc_create_f(trim(sf%path), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'state nc_create')
      sf%ncid = ncid
      d%ncid = ncid
      call define_polygon_axis(ncid, size(polygon_id, kind=ik), d%polygon)
      call nc_check(nc_def_dim_f(ncid, 'cohort', int(max(n_cohort,1_ik), c_size_t), d%cohort), 'state dim cohort')
      call nc_check(nc_def_dim_f(ncid, 'patch',  int(max(n_patch,1_ik), c_size_t), d%patch),   'state dim patch')
      call nc_check(nc_def_dim_f(ncid, 'site',   int(max(n_site,1_ik), c_size_t),  d%site),    'state dim site')
      call nc_check(nc_def_dim_f(ncid, 'pft',    int(n_pft, c_size_t),             d%pft),     'state dim pft')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_int',  int(N_META_INT, c_size_t), d%meta_int), 'state dim mi')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_real', int(N_META_REAL, c_size_t), d%meta_real), 'state dim mr')
      call nc_check(nc_def_dim_f(ncid, 'soil_layer', int(n_soil_layer_max, c_size_t), d%soil_layer), 'state dim soill')
      call nc_check(nc_def_dim_f(ncid, 'snow_layer', int(n_snow_layer_max, c_size_t), d%snow_layer), 'state dim snowl')
      call nc_check(nc_def_dim_f(ncid, 'growth_window', int(max(n_window,1_ik), c_size_t), d%growth_window),     &
                    'state dim gwin')
      call define_polygon_vars(d)
      call define_cohort_vars(d)
      call define_patch_vars(d)
      call put_global_text(ncid, 'title', STATE_TITLE)
      pending = 'none'
      if (present(restructure_pending)) then
         if (restructure_pending) pending = 'month'
         if (restructure_pending .and. present(restructure_new_year)) then
            if (restructure_new_year) pending = 'year'
         end if
      end if
      call put_global_text(ncid, 'restructure_pending', trim(pending))
      call nc_check(nc_enddef(ncid), 'state enddef')
      call write_polygon_axis(ncid, polygon_id, polygon_lat, polygon_lon)
   end subroutine restart_file_open

   !--- Write one polygon's state into its row (1-based) of an open file.
   subroutine restart_file_write(sf, poly, now, polygon)
      type(restart_file_t), intent(in) :: sf
      type(polygon_t),      intent(in) :: poly
      type(ed_time_t),      intent(in) :: now
      integer(ik),          intent(in) :: polygon
      type(file_row_t) :: r
      r = file_row_t(ncid=sf%ncid, row=polygon)
      call put_polygon_vars(r, poly, now)
      if (poly%cohort%n > 0_ik) call put_cohort_vars(r, poly)
      if (poly%patch%n > 0_ik)  call put_patch_vars(r, poly)
   end subroutine restart_file_write

   subroutine restart_file_close(sf)
      type(restart_file_t), intent(inout) :: sf
      call nc_check(nc_close(sf%ncid), 'state nc_close')
      sf%ncid = -1_c_int
      write(*,'(2a)') ' state : ', trim(sf%path)
   end subroutine restart_file_close

   !---~---
   !   Reconstruct a polygon from row `polygon` (1-based) of a state file written by
   ! restart_file_write. Each site's surface is made from `site_cfg` (else `cfg`) and the
   ! geometry from the carbon pools; the maps and sort order are rebuilt, and the
   ! restructuring owed is the file's restructure_pending
   ! (docs/science/order_of_processes.md §5). found = .false. and `why` say the file
   ! cannot be opened or has a PFT, site or soil-layer count unlike the run's; a file
   ! without one of the variables stops the run, naming it.
   !---~---
   subroutine read_restart_file(poly, cfg, path, restart_time, found, polygon, site_cfg, why)
      type(polygon_t),     intent(out) :: poly
      type(ed_config_t),   intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      type(ed_time_t),     intent(out) :: restart_time
      logical,             intent(out) :: found
      integer(ik),         intent(in)  :: polygon    ! its row (1-based)
      type(ed_config_t), optional, intent(in) :: site_cfg(:)   ! each site's configuration
      character(len=:), allocatable, optional, intent(out) :: why  ! why found = .false.; '' when found
      character(len=:), allocatable :: msg
      type(file_row_t) :: r
      integer(c_int)   :: ncid, st
      logical          :: same_window

      found = .false.
      restart_time = ed_time_t()
      if (present(why)) why = 'cannot be opened'
      if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) return
      r = file_row_t(ncid=ncid, row=polygon)
      call get_polygon_vars(r, poly, cfg, restart_time, msg, site_cfg)
      if (len(msg) > 0) then
         st = nc_close(ncid)
         if (present(why)) why = msg
         return
      end if
      same_window = .true.
      if (poly%cohort%n > 0_ik) call get_cohort_vars(r, poly, cfg, same_window)
      if (poly%patch%n > 0_ik)  call get_patch_vars(r, poly)
      call nc_check(nc_close(ncid), 'state nc_close (read)')
      call derive_after_read(poly, cfg, same_window)
      if (present(why)) why = ''
      found = .true.
   end subroutine read_restart_file

   !---~---
   !   Each polygon's row of the restart file at `path`, by its id. why = '' when the file can
   ! start a run of these polygons, else the reason it cannot: it cannot be opened, or does
   ! not hold one of them.
   !---~---
   subroutine restart_file_slots(path, ids, slot, why)
      character(len=*),              intent(in)  :: path
      integer(ik),                   intent(in)  :: ids(:)
      integer(ik), allocatable,      intent(out) :: slot(:)
      character(len=:), allocatable, intent(out) :: why
      integer(ik), allocatable :: file_ids(:)
      integer(ik)       :: p
      integer(c_int)    :: ncid, st
      character(len=24) :: id_str
      allocate(slot(size(ids)), source=0_ik)
      if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) then
         why = 'cannot be opened'
         return
      end if
      why = ''
      call file_polygon_ids(ncid, file_ids)
      st = nc_close(ncid)
      do p = 1_ik, size(ids, kind=ik)
         slot(p) = int(findloc(file_ids, ids(p), dim=1), ik)
         if (slot(p) == 0_ik) then
            write(id_str, '(i0)') ids(p)
            why = 'does not hold polygon '//trim(id_str)//', and must hold every polygon of the run'
            return
         end if
      end do
   end subroutine restart_file_slots

   !--- An open state file's polygon ids, in its rows' order.
   subroutine file_polygon_ids(ncid, polygon_id)
      integer(c_int),           intent(in)  :: ncid
      integer(ik), allocatable, intent(out) :: polygon_id(:)
      integer(c_size_t) :: n, start1(1), count1(1)
      call nc_check(nc_inq_dimlen_f(ncid, 'polygon', n), 'state inq polygon')
      allocate(polygon_id(n))
      start1 = 0_c_size_t
      count1 = n
      call nc_check(nc_get_vara_int(ncid, var_id(ncid, 'polygon_id'), start1, count1, polygon_id),         &
                    'get polygon_id')
   end subroutine file_polygon_ids

   !--- A per-polygon variable: one value per row, uncompressed, with its long_name.
   subroutine def_polygon_var(d, name, xtype, lname)
      type(file_dims_t), intent(in) :: d
      character(len=*),  intent(in) :: name, lname
      integer(c_int),    intent(in) :: xtype
      integer(c_int) :: v
      call define_var(d%ncid, name, xtype, [d%polygon], v, long_name=lname)
   end subroutine def_polygon_var

   !---~---
   !   A state variable on its own dimensions, behind the polygon axis, with its long_name and,
   ! where the reader needs it, its units; one polygon's row per chunk, since each polygon
   ! writes its own.
   !---~---
   subroutine def_var(d, name, xtype, dimids, lname, units)
      type(file_dims_t), intent(in) :: d
      character(len=*),  intent(in) :: name, lname
      integer(c_int),    intent(in) :: xtype, dimids(:)
      character(len=*),  intent(in), optional :: units
      integer(c_int)    :: v, dims(size(dimids) + 1)
      integer(c_size_t) :: chunk(size(dimids) + 1)
      integer           :: k
      dims(1) = d%polygon
      dims(2:) = dimids
      chunk(1) = 1_c_size_t
      do k = 2, size(dims)
         call nc_check(nc_inq_dimlen(d%ncid, dims(k), chunk(k)), 'state dimlen '//name)
         chunk(k) = max(chunk(k), 1_c_size_t)
      end do
      call define_var(d%ncid, name, xtype, dims, v, units=units, long_name=lname, chunk=chunk,     &
                      deflate=RESTART_DEFLATE)
   end subroutine def_var

   !--- A per-polygon value (n_cohort, n_patch, n_site, n_window), at the polygon's row.
   subroutine put_polygon_int(r, name, n)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: n
      integer(c_size_t) :: st(1), ct(1)
      integer(ik)       :: x(1)
      st = int(r%row - 1_ik, c_size_t)
      ct = 1_c_size_t
      x = n
      call nc_check(nc_put_vara_int(r%ncid, var_id(r%ncid, name), st, ct, x), 'put '//name)
   end subroutine put_polygon_int

   subroutine put_polygon_dbl(r, name, x)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: x
      integer(c_size_t) :: st(1), ct(1)
      real(wp)          :: y(1)
      st = int(r%row - 1_ik, c_size_t)
      ct = 1_c_size_t
      y = x
      call nc_check(nc_put_vara_double(r%ncid, var_id(r%ncid, name), st, ct, y), 'put '//name)
   end subroutine put_polygon_dbl

   real(wp) function get_polygon_dbl(r, name) result(x)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(c_size_t) :: st(1), ct(1)
      real(wp)          :: y(1)
      st = int(r%row - 1_ik, c_size_t)
      ct = 1_c_size_t
      call nc_check(nc_get_vara_double(r%ncid, var_id(r%ncid, name), st, ct, y), 'get '//name)
      x = y(1)
   end function get_polygon_dbl

   integer(ik) function get_polygon_int(r, name) result(n)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(c_size_t) :: st(1), ct(1)
      integer(ik)       :: x(1)
      st = int(r%row - 1_ik, c_size_t)
      ct = 1_c_size_t
      call nc_check(nc_get_vara_int(r%ncid, var_id(r%ncid, name), st, ct, x), 'get '//name)
      n = x(1)
   end function get_polygon_int

   !---~---
   !   One vector into, or out of, the variable `name` at the polygon's row: its first entries
   ! (a per-cohort, per-patch, per-site or meta vector), or with `item`, the layers of that
   ! cohort or patch (1-based) of a variable with a second axis.
   !---~---
   subroutine put_dbl(r, name, x, item)
      type(file_row_t),      intent(in) :: r
      character(len=*),      intent(in) :: name
      real(wp),              intent(in) :: x(:)
      integer(ik), optional, intent(in) :: item
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), item, st, ct, nd)
      call nc_check(nc_put_vara_double(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'put '//name)
   end subroutine put_dbl

   subroutine put_int(r, name, x)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: x(:)
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), st=st, ct=ct, nd=nd)
      call nc_check(nc_put_vara_int(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'put '//name)
   end subroutine put_int

   subroutine get_dbl(r, name, x, item)
      type(file_row_t),      intent(in)  :: r
      character(len=*),      intent(in)  :: name
      real(wp),              intent(out) :: x(:)
      integer(ik), optional, intent(in)  :: item
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), item, st, ct, nd)
      call nc_check(nc_get_vara_double(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'get '//name)
   end subroutine get_dbl

   subroutine get_int(r, name, x)
      type(file_row_t), intent(in)  :: r
      character(len=*), intent(in)  :: name
      integer(ik),      intent(out) :: x(:)
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), st=st, ct=ct, nd=nd)
      call nc_check(nc_get_vara_int(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'get '//name)
   end subroutine get_int

   integer(c_int) function var_id(ncid, name)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name
      call nc_check(nc_inq_varid_f(ncid, name, var_id), 'state inq '//name)
   end function var_id

   !---~---
   !   Where a polygon's n values sit: its row; then entry `item` (1-based) of the cohort or
   ! patch axis, when given; then the first n entries.
   !---~---
   pure subroutine where_in_file(r, n, item, st, ct, nd)
      type(file_row_t),      intent(in)  :: r
      integer(ik),           intent(in)  :: n
      integer(ik), optional, intent(in)  :: item
      integer(c_size_t),     intent(out) :: st(3), ct(3)
      integer,               intent(out) :: nd
      nd = 1
      st(1) = int(r%row - 1_ik, c_size_t)
      ct(1) = 1_c_size_t
      if (present(item)) then
         nd = nd + 1
         st(nd) = int(item - 1_ik, c_size_t)
         ct(nd) = 1_c_size_t
      end if
      nd = nd + 1
      st(nd) = 0_c_size_t
      ct(nd) = int(n, c_size_t)
   end subroutine where_in_file

   !---~---
   !   1. THE POLYGON AND ITS SITES, in ed_polygon_state_vars' and ed_site_state_vars'
   !    section order
   !---~---

   subroutine define_polygon_vars(d)
      type(file_dims_t), intent(in) :: d
      !--- 1. Identity and location: the polygon axis (restart_file_open).
      !--- 2. Dimensions and id counters, with the file's time in meta_int.
      call def_polygon_var(d, 'n_cohort',   NC_INT,    'the polygon''s live cohorts (the rest is fill)')
      call def_polygon_var(d, 'n_patch',    NC_INT,    'the polygon''s live patches (the rest is fill)')
      call def_polygon_var(d, 'n_site',     NC_INT,    'the polygon''s sites (the rest is fill)')
      call def_polygon_var(d, 'n_window',   NC_INT,    'the polygon''s growth-window length')
      call def_var(d, 'meta_int',  NC_INT,    [d%meta_int],                                  &
                   '[n_cohort,n_patch,n_pft,next_cohort_id,next_patch_id,year,month,day,hour,minute,second]')
      call def_var(d, 'meta_real', NC_DOUBLE, [d%meta_real], '[area of the sites (1), decimal_year]')
      !---~---
      !   3. Calendar: the cohorts' ring-buffer slot (the restructuring owed is a global
      ! attribute).
      !---~---
      call def_var(d, 'growth_hist_pos_py', NC_INT, [d%meta_real], 'ring-buffer write slot (element 1)')
      !--- 4. The forcing's memory.
      call def_polygon_var(d, 'kt_last_day_py', NC_DOUBLE, 'the last daytime clearness index [-]')
      !---~---
      !   The sites: 1. each one's area and 3. the surface it was made with, as the diagnostic
      ! files carry them (its layer count is checked on read); 4. thermal acclimation's
      ! running mean (docs/science/leaf_gas_exchange.md §2).
      !---~---
      call define_site_statics(d%ncid, d%polygon, d%site, [SOIL_BC_DRAIN, SOIL_BC_AQUIFER], 'drain aquifer',    &
                               by_row=.true., deflate=RESTART_DEFLATE)
      call def_var(d, 't_growth_avg_si', NC_DOUBLE, [d%site], 'the site''s running-mean growth temperature [K]')
   end subroutine define_polygon_vars

   subroutine put_polygon_vars(r, poly, now)
      type(file_row_t),  intent(in) :: r
      type(polygon_t),   intent(in) :: poly
      type(ed_time_t),   intent(in) :: now
      integer(ik) :: ns, is, nsoil(poly%site%n, 1), bc(poly%site%n, 1)
      real(wp)    :: area(poly%site%n, 1), elev(poly%site%n, 1), slope(poly%site%n, 1), aspect(poly%site%n, 1)
      real(wp)    :: depth(poly%site%n, 1), drain(poly%site%n, 1), tdeep(poly%site%n, 1)
      ns = poly%site%n
      !--- 2. Dimensions and id counters, and the time.
      call put_polygon_int(r, 'n_cohort', poly%cohort%n)
      call put_polygon_int(r, 'n_patch',  poly%patch%n)
      call put_polygon_int(r, 'n_site',   ns)
      call put_polygon_int(r, 'n_window', int(size(poly%cohort%growth_hist, 1), ik))
      call put_int(r, 'meta_int', [poly%cohort%n, poly%patch%n, poly%n_pft, poly%next_cohort_id,              &
                                   poly%next_patch_id, now%year, now%month, now%day, now%hour, now%minute,     &
                                   now%second])
      call put_dbl(r, 'meta_real', [sum(poly%site%area(1:ns)), time_to_decimal_year(now)])
      !--- 3. Calendar.
      call put_int(r, 'growth_hist_pos_py', [poly%growth_hist_pos])
      !--- 4. The forcing's memory.
      call put_polygon_dbl(r, 'kt_last_day_py', poly%kt_last_day)
      !--- The sites, the surface's numbers through contiguous buffers.
      do is = 1_ik, ns
         associate (surf => poly%site%surface(is))
            area(is, 1)   = poly%site%area(is)
            elev(is, 1)   = surf%elevation_m
            slope(is, 1)  = surf%slope_deg
            aspect(is, 1) = surf%aspect_deg
            nsoil(is, 1)  = surf%soil_profile%n_active
            depth(is, 1)  = surf%soil%depth
            bc(is, 1)     = surf%soil%water_bottom_bc
            drain(is, 1)  = merge(real(NC_FILL_DOUBLE, wp), surf%soil%drain_factor, surf%soil%water_bottom_bc == SOIL_BC_AQUIFER)
            tdeep(is, 1)  = surf%soil%deep_temp
         end associate
      end do
      call write_site_statics(r%ncid, r%row - 1_ik, area, elev, slope, aspect, nsoil, depth, bc, drain, tdeep)
      call put_dbl(r, 't_growth_avg_si', poly%site%t_growth_avg(1:ns))
   end subroutine put_polygon_vars

   !---~---
   !   The polygon, allocated to the file's counts, and its sites, each surface made from
   ! site_cfg (else cfg). why = '' unless the file's PFT, site or soil-layer count is unlike
   ! the run's.
   !---~---
   subroutine get_polygon_vars(r, poly, cfg, restart_time, why, site_cfg)
      type(file_row_t),              intent(in)  :: r
      type(polygon_t),               intent(out) :: poly
      type(ed_config_t),           intent(in)    :: cfg
      type(ed_time_t),               intent(out) :: restart_time
      character(len=:), allocatable, intent(out) :: why
      type(ed_config_t), intent(in), optional :: site_cfg(:)
      character(len=16)  :: pending
      character(len=160) :: msg
      integer(ik) :: meta_i(N_META_INT), ns, is

      why = ''
      !--- 2. Dimensions and id counters, and the file's time.
      call get_int(r, 'meta_int',  meta_i)
      if (meta_i(3) /= cfg%pft%n) then
         write(msg, '(a,i0,a,i0)') 'holds ', meta_i(3), ' PFTs and the configuration ', cfg%pft%n
         why = trim(msg)
         return
      end if
      ns = get_polygon_int(r, 'n_site')
      if (present(site_cfg)) then
         if (size(site_cfg, kind=ik) /= ns) then
            write(msg, '(a,i0,a,i0,a)') 'holds ', ns, ' sites and the configuration ', size(site_cfg),        &
                                        ' for the polygon'
            why = trim(msg)
            return
         end if
      end if
      call alloc_polygon(poly, cfg%pft%n, max(meta_i(1), 1_ik), max(meta_i(2), 1_ik), growth_window_steps(cfg), &
                         cfg%pft, n_site=ns)
      poly%cohort%n = meta_i(1)
      poly%kt_last_day = get_polygon_dbl(r, 'kt_last_day_py')
      poly%patch%n = meta_i(2)
      poly%next_cohort_id = meta_i(4)
      poly%next_patch_id = meta_i(5)
      restart_time = ed_time_t(year=meta_i(6), month=meta_i(7), day=meta_i(8),               &
                               hour=meta_i(9), minute=meta_i(10), second=meta_i(11))
      !--- 3. Calendar: the restructuring owed (the ring's slot comes with the rings, section 2).
      call nc_check(nc_get_att_text_f(r%ncid, NC_GLOBAL, 'restructure_pending', pending), 'state restructure_pending')
      poly%restructure_pending  = trim(pending) == 'month' .or. trim(pending) == 'year'
      poly%restructure_new_year = trim(pending) == 'year'
      !---~---
      !   The sites: 3. each surface; 1. each area, and each site's soil state is per layer of
      ! the column it was written on.
      !---~---
      do is = 1_ik, ns
         if (present(site_cfg)) then
            call make_site_surface(poly%site%surface(is), site_cfg(is))
         else
            call make_site_surface(poly%site%surface(is), cfg)
         end if
      end do
      call get_dbl(r, 'area_si', poly%site%area(1:ns))
      block
         integer(ik) :: nsoil(ns)
         call get_int(r, 'n_soil_si', nsoil)
         do is = 1_ik, ns
            if (nsoil(is) /= poly%site%surface(is)%soil_profile%n_active) then
               write(msg, '(a,i0,a,i0,a,i0)') 'has ', nsoil(is), ' soil layers at site ', is,                &
                     ' and the configuration ', poly%site%surface(is)%soil_profile%n_active
               why = trim(msg)
               return
            end if
         end do
      end block
      !--- 4. Growth temperature.
      call get_dbl(r, 't_growth_avg_si', poly%site%t_growth_avg(1:ns))
   end subroutine get_polygon_vars

   !--- 2. THE COHORTS, in ed_cohort_state_vars' section order

   subroutine define_cohort_vars(d)
      type(file_dims_t), intent(in) :: d
      !--- 1. Identity and density.
      call def_var(d, 'pft_co',              NC_INT,    [d%cohort], 'plant functional type index')
      call def_var(d, 'global_cohort_id', NC_INT,    [d%cohort], 'persistent cohort id')
      call def_var(d, 'owner_patch',      NC_INT,    [d%cohort], 'owning patch index (1-based)')
      call def_var(d, 'nplant_co',           NC_DOUBLE, [d%cohort], 'plant number density [plant/m2]')
      !--- 2. The carbon pools, as the run holds them.
      call def_var(d, 'leaf_carbon_co',          NC_DOUBLE, [d%cohort], 'leaf carbon per plant [kgC]')
      call def_var(d, 'fineroot_carbon_co',      NC_DOUBLE, [d%cohort], 'fine-root carbon per plant [kgC]')
      call def_var(d, 'wood_carbon_co',          NC_DOUBLE, [d%cohort], 'wood carbon per plant [kgC]')
      call def_var(d, 'nonstructural_carbon_co', NC_DOUBLE, [d%cohort], 'nonstructural carbon per plant [kgC]')
      !---~---
      !   3. Size: the reader derives it from the pools, so dbh, height, basal_area, agb and
      ! leaf_area are for readers of the file; overtopping_lai is the cohort's light.
      !---~---
      call def_var(d, 'dbh_co',             NC_DOUBLE, [d%cohort], 'diameter at breast height [cm]')
      call def_var(d, 'height_co',          NC_DOUBLE, [d%cohort], 'height [m]')
      call def_var(d, 'basal_area_co',      NC_DOUBLE, [d%cohort], 'basal area per plant [cm2]')
      call def_var(d, 'agb_co',             NC_DOUBLE, [d%cohort], 'aboveground biomass per plant [kgC]')
      call def_var(d, 'leaf_area_co',       NC_DOUBLE, [d%cohort], 'leaf area per plant [m2]')
      call def_var(d, 'overtopping_lai_co', NC_DOUBLE, [d%cohort], 'LAI of the taller cohorts in the patch [m2/m2]')
      !---~---
      !   4. Water: the tissue water; the canopy interception films ([fast].canopy_water_on),
      ! per plant and marked so; the stomatal feedback's predawn potential
      ! (docs/science/leaf_gas_exchange.md §4.2).
      !---~---
      call def_var(d, 'leaf_water_mass_co', NC_DOUBLE, [d%cohort], 'per-cohort internal leaf water [kg/plant]')
      call def_var(d, 'wood_water_mass_co', NC_DOUBLE, [d%cohort], 'per-cohort internal wood water [kg/plant]')
      call def_var(d, 'leaf_surf_water_co', NC_DOUBLE, [d%cohort], 'water on the leaf surfaces [kg/plant]',          &
                   units=SURF_WATER_UNITS)
      call def_var(d, 'wood_surf_water_co', NC_DOUBLE, [d%cohort], 'water on the wood surfaces [kg/plant]',          &
                   units=SURF_WATER_UNITS)
      call def_var(d, 'dmax_psi_leaf_co',   NC_DOUBLE, [d%cohort],                                                 &
                   'yesterday daily-max leaf water potential [MPa] (drives beta_stomata)')
      !--- 5. Energy: the tissue energies, with their temperatures as twins.
      call def_var(d, 'leaf_energy_co', NC_DOUBLE, [d%cohort], 'per-cohort leaf internal energy [J/plant] (prognostic)')
      call def_var(d, 'wood_energy_co', NC_DOUBLE, [d%cohort], 'per-cohort wood internal energy [J/plant] (prognostic)')
      call def_var(d, 'leaf_temp_co',   NC_DOUBLE, [d%cohort], 'per-cohort leaf temperature [K] (diagnosed twin)')
      call def_var(d, 'wood_temp_co',   NC_DOUBLE, [d%cohort], 'per-cohort wood temperature [K] (diagnosed twin)')
      !--- 6. Plastic leaf traits.
      call def_var(d, 'sla_co',                NC_DOUBLE, [d%cohort], 'specific leaf area [m2/kgC] (plastic)')
      call def_var(d, 'vcmax25_co',            NC_DOUBLE, [d%cohort], 'max carboxylation @25C [umol/m2/s] (plastic)')
      call def_var(d, 'rd25_co',               NC_DOUBLE, [d%cohort], 'leaf dark respiration @25C [umol/m2/s] (plastic)')
      call def_var(d, 'leaf_turnover_rate_co', NC_DOUBLE, [d%cohort], 'baseline leaf turnover rate [1/yr] (plastic)')
      !--- 7. Leaf phenology: the season's cue sums and running light mean.
      call def_var(d, 'leaf_flush_tendency_co', NC_DOUBLE, [d%cohort], 'smoothed flush signal [0,1]')
      call def_var(d, 'leaf_shed_tendency_co',  NC_DOUBLE, [d%cohort], 'smoothed senescence signal [0,1]')
      call def_var(d, 'growing_degree_days_co', NC_DOUBLE, [d%cohort], 'warmth sum since midwinter [K day]')
      call def_var(d, 'cold_degree_days_co',    NC_DOUBLE, [d%cohort], 'cold sum since midsummer [K day]')
      call def_var(d, 'dry_psi_sum_co',         NC_DOUBLE, [d%cohort], 'predawn leaf psi below the TLP, summed [MPa day]')
      call def_var(d, 'wet_psi_sum_co',         NC_DOUBLE, [d%cohort], 'predawn leaf psi above the TLP, summed [MPa day]')
      call def_var(d, 'light_hours_avg_co',    NC_DOUBLE, [d%cohort],                                             &
                   'running-mean hours of light at the cohort top [h/day]')
      !--- 8. Growth history: growth_avg (the mortality predictor) and the ring buffer behind it.
      call def_var(d, 'growth_avg_co',   NC_DOUBLE, [d%cohort], 'simple moving-average growth [cm/yr] (mortality predictor)')
      call def_var(d, 'growth_accum_co', NC_DOUBLE, [d%cohort], 'running sum of the growth window [cm/yr]')
      call def_var(d, 'growth_count_co', NC_INT,    [d%cohort], 'samples in the growth window')
      call def_var(d, 'growth_hist_co',  NC_DOUBLE, [d%cohort, d%growth_window], 'growth-rate ring buffer [cm/yr]')
      !--- 9. Fast-to-slow accumulators: the day's running maximum leaf potential.
      call def_var(d, 'dmax_psi_leaf_accum_co', NC_DOUBLE, [d%cohort],                                             &
                   'running daily-max leaf water potential accumulator [MPa]')
   end subroutine define_cohort_vars

   subroutine put_cohort_vars(r, poly)
      type(file_row_t), intent(in) :: r
      type(polygon_t),  intent(in) :: poly
      integer(ik) :: ncoh, nwin, ic
      ncoh = poly%cohort%n
      nwin = int(size(poly%cohort%growth_hist, 1), ik)
      associate (c => poly%cohort)
         !--- 1. Identity and density.
         call put_int(r, 'pft_co',              c%pft(1:ncoh))
         call put_int(r, 'global_cohort_id', c%global_id(1:ncoh))
         call put_int(r, 'owner_patch',      c%owner_patch(1:ncoh))
         call put_dbl(r, 'nplant_co',           c%nplant(1:ncoh))
         !--- 2. Carbon pools.
         call put_dbl(r, 'leaf_carbon_co',          c%leaf_carbon(1:ncoh))
         call put_dbl(r, 'fineroot_carbon_co',      c%fineroot_carbon(1:ncoh))
         call put_dbl(r, 'wood_carbon_co',          c%wood_carbon(1:ncoh))
         call put_dbl(r, 'nonstructural_carbon_co', c%nonstructural_carbon(1:ncoh))
         !--- 3. Size.
         call put_dbl(r, 'dbh_co',             c%dbh(1:ncoh))
         call put_dbl(r, 'height_co',          c%height(1:ncoh))
         call put_dbl(r, 'basal_area_co',      c%basal_area(1:ncoh))
         call put_dbl(r, 'agb_co',             c%agb(1:ncoh))
         call put_dbl(r, 'leaf_area_co',       c%leaf_area(1:ncoh))
         call put_dbl(r, 'overtopping_lai_co', c%overtopping_lai(1:ncoh))
         !--- 4. Water.
         call put_dbl(r, 'leaf_water_mass_co', c%leaf_water_mass(1:ncoh))
         call put_dbl(r, 'wood_water_mass_co', c%wood_water_mass(1:ncoh))
         call put_dbl(r, 'leaf_surf_water_co', c%leaf_surf_water(1:ncoh))
         call put_dbl(r, 'wood_surf_water_co', c%wood_surf_water(1:ncoh))
         call put_dbl(r, 'dmax_psi_leaf_co',   c%dmax_psi_leaf(1:ncoh))
         !--- 5. Energy.
         call put_dbl(r, 'leaf_energy_co', c%leaf_energy(1:ncoh))
         call put_dbl(r, 'wood_energy_co', c%wood_energy(1:ncoh))
         call put_dbl(r, 'leaf_temp_co',   c%leaf_temp(1:ncoh))
         call put_dbl(r, 'wood_temp_co',   c%wood_temp(1:ncoh))
         !--- 6. Plastic leaf traits.
         call put_dbl(r, 'sla_co',                c%sla(1:ncoh))
         call put_dbl(r, 'vcmax25_co',            c%vcmax25(1:ncoh))
         call put_dbl(r, 'rd25_co',               c%rd25(1:ncoh))
         call put_dbl(r, 'leaf_turnover_rate_co', c%leaf_turnover_rate(1:ncoh))
         !--- 7. Leaf phenology.
         call put_dbl(r, 'leaf_flush_tendency_co', c%leaf_flush_tendency(1:ncoh))
         call put_dbl(r, 'leaf_shed_tendency_co',  c%leaf_shed_tendency(1:ncoh))
         call put_dbl(r, 'growing_degree_days_co', c%growing_degree_days(1:ncoh))
         call put_dbl(r, 'cold_degree_days_co',    c%cold_degree_days(1:ncoh))
         call put_dbl(r, 'dry_psi_sum_co',         c%dry_psi_sum(1:ncoh))
         call put_dbl(r, 'wet_psi_sum_co',         c%wet_psi_sum(1:ncoh))
         call put_dbl(r, 'light_hours_avg_co',    c%light_hours_mean(1:ncoh))
         !--- 8. Growth history.
         call put_dbl(r, 'growth_avg_co',   c%growth_avg(1:ncoh))
         call put_dbl(r, 'growth_accum_co', c%growth_accum(1:ncoh))
         call put_int(r, 'growth_count_co', c%growth_count(1:ncoh))
         if (nwin > 0_ik) then
            do ic = 1_ik, ncoh
               call put_dbl(r, 'growth_hist_co', c%growth_hist(1:nwin, ic), item=ic)
            end do
         end if
         !--- 9. Fast-to-slow accumulators.
         call put_dbl(r, 'dmax_psi_leaf_accum_co', c%dmax_psi_leaf_accum(1:ncoh))
      end associate
   end subroutine put_cohort_vars

   !---~---
   !   The cohorts as the file holds them. same_window says the file's growth window is this
   ! run's length, so its rings are read; otherwise section 4 seeds them from growth_avg.
   !---~---
   subroutine get_cohort_vars(r, poly, cfg, same_window)
      type(file_row_t),    intent(in)    :: r
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      logical,             intent(out)   :: same_window
      integer(c_size_t) :: nwin_file
      integer(ik)       :: ncoh, nwin, i, gpos(1)
      ncoh = poly%cohort%n
      associate (c => poly%cohort)
         !--- 1. Identity and density.
         call get_int(r, 'pft_co',              c%pft(1:ncoh))
         call get_int(r, 'global_cohort_id', c%global_id(1:ncoh))
         call get_int(r, 'owner_patch',      c%owner_patch(1:ncoh))
         call get_dbl(r, 'nplant_co',           c%nplant(1:ncoh))
         !--- 2. Carbon pools.
         call get_dbl(r, 'leaf_carbon_co',          c%leaf_carbon(1:ncoh))
         call get_dbl(r, 'fineroot_carbon_co',      c%fineroot_carbon(1:ncoh))
         call get_dbl(r, 'wood_carbon_co',          c%wood_carbon(1:ncoh))
         call get_dbl(r, 'nonstructural_carbon_co', c%nonstructural_carbon(1:ncoh))
         !---~---
         !   3. Size: derived from the pools in section 4; the overtopping LAI is the cohort's
         ! light.
         !---~---
         call get_dbl(r, 'overtopping_lai_co', c%overtopping_lai(1:ncoh))
         !--- 4. Water.
         call get_dbl(r, 'leaf_water_mass_co', c%leaf_water_mass(1:ncoh))
         call get_dbl(r, 'wood_water_mass_co', c%wood_water_mass(1:ncoh))
         call get_dbl(r, 'leaf_surf_water_co', c%leaf_surf_water(1:ncoh))
         call get_dbl(r, 'wood_surf_water_co', c%wood_surf_water(1:ncoh))
         call get_dbl(r, 'dmax_psi_leaf_co',   c%dmax_psi_leaf(1:ncoh))
         !--- 5. Energy.
         call get_dbl(r, 'leaf_energy_co', c%leaf_energy(1:ncoh))
         call get_dbl(r, 'wood_energy_co', c%wood_energy(1:ncoh))
         call get_dbl(r, 'leaf_temp_co',   c%leaf_temp(1:ncoh))
         call get_dbl(r, 'wood_temp_co',   c%wood_temp(1:ncoh))
         !--- 6. Plastic leaf traits.
         call get_dbl(r, 'sla_co',                c%sla(1:ncoh))
         call get_dbl(r, 'vcmax25_co',            c%vcmax25(1:ncoh))
         call get_dbl(r, 'rd25_co',               c%rd25(1:ncoh))
         call get_dbl(r, 'leaf_turnover_rate_co', c%leaf_turnover_rate(1:ncoh))
         !--- 7. Leaf phenology.
         call get_dbl(r, 'leaf_flush_tendency_co', c%leaf_flush_tendency(1:ncoh))
         call get_dbl(r, 'leaf_shed_tendency_co',  c%leaf_shed_tendency(1:ncoh))
         call get_dbl(r, 'growing_degree_days_co', c%growing_degree_days(1:ncoh))
         call get_dbl(r, 'cold_degree_days_co',    c%cold_degree_days(1:ncoh))
         call get_dbl(r, 'dry_psi_sum_co',         c%dry_psi_sum(1:ncoh))
         call get_dbl(r, 'wet_psi_sum_co',         c%wet_psi_sum(1:ncoh))
         call get_dbl(r, 'light_hours_avg_co',    c%light_hours_mean(1:ncoh))
         !---~---
         !   8. Growth history: the ring buffer, with the polygon's write slot, when the file's
         ! window is this run's length.
         !---~---
         call get_dbl(r, 'growth_avg_co', c%growth_avg(1:ncoh))
         nwin = growth_window_steps(cfg)
         call nc_check(nc_inq_dimlen_f(r%ncid, 'growth_window', nwin_file), 'inq growth_window')
         same_window = int(nwin_file, ik) == nwin .and. int(size(c%growth_hist, 1), ik) == nwin
         if (same_window) then
            call get_dbl(r, 'growth_accum_co', c%growth_accum(1:ncoh))
            call get_int(r, 'growth_count_co', c%growth_count(1:ncoh))
            do i = 1_ik, ncoh
               call get_dbl(r, 'growth_hist_co', c%growth_hist(1:nwin, i), item=i)
            end do
            call get_int(r, 'growth_hist_pos_py', gpos)
            poly%growth_hist_pos = gpos(1)
         end if
         !--- 9. Fast-to-slow accumulators.
         call get_dbl(r, 'dmax_psi_leaf_accum_co', c%dmax_psi_leaf_accum(1:ncoh))
      end associate
   end subroutine get_cohort_vars

   !--- 3. THE PATCHES, in ed_patch_state_vars' section order

   subroutine define_patch_vars(d)
      type(file_dims_t), intent(in) :: d
      !--- 1. Identity and area.
      call def_var(d, 'global_patch_id', NC_INT,    [d%patch], 'persistent patch id')
      call def_var(d, 'dist_type_pa',       NC_INT,    [d%patch], 'disturbance type (1=primary,2=treefall)')
      call def_var(d, 'owner_site',      NC_INT,    [d%patch], 'owning site index (1-based)')
      call def_var(d, 'area_pa',      NC_DOUBLE, [d%patch], 'patch area fraction of its site')
      call def_var(d, 'age_pa',       NC_DOUBLE, [d%patch], 'time since last disturbance [yr]')
      !--- 3. Recruit bank.
      call def_var(d, 'recruit_pool_pa', NC_DOUBLE, [d%patch, d%pft], 'carry-forward recruit pool [plant/m2]')
      !--- 4-6. The fast reservoirs (canopy air, snow, soil, pond), with their diagnosed twins.
      call def_var(d, 'cas_can_enthalpy_pa', NC_DOUBLE, [d%patch], 'CAS specific enthalpy [J/kg] (prognostic)')
      call def_var(d, 'cas_can_shv_pa',      NC_DOUBLE, [d%patch], 'CAS specific humidity [kg/kg] (prognostic)')
      call def_var(d, 'cas_can_co2_pa',      NC_DOUBLE, [d%patch], 'CAS CO2 mixing ratio [umol/mol] (prognostic)')
      call def_var(d, 'cas_can_temp_pa',     NC_DOUBLE, [d%patch], 'CAS temperature [K] (diagnosed twin)')
      call def_var(d, 'snow_nlayer_pa', NC_INT,    [d%patch], 'active snow layer count (0 = no snow)')
      call def_var(d, 'snow_swe_layer_pa',    NC_DOUBLE, [d%patch, d%snow_layer], 'snow water-equivalent mass [kg/m2] (prognostic)')
      call def_var(d, 'snow_energy_layer_pa', NC_DOUBLE, [d%patch, d%snow_layer],                                   &
                   'snow extensive internal energy [J/m2] (prognostic)')
      call def_var(d, 'snow_depth_layer_pa',  NC_DOUBLE, [d%patch, d%snow_layer], 'snow geometric depth [m] (prognostic)')
      call def_var(d, 'snow_temp_layer_pa',   NC_DOUBLE, [d%patch, d%snow_layer], 'snow temperature [K] (diagnosed twin)')
      call def_var(d, 'snow_fliq_layer_pa',   NC_DOUBLE, [d%patch, d%snow_layer], 'snow liquid fraction [-] (diagnosed twin)')
      call def_var(d, 'soil_theta_layer_pa',  NC_DOUBLE, [d%patch, d%soil_layer], 'soil volumetric moisture [m3/m3] (prognostic)')
      call def_var(d, 'soil_energy_layer_pa', NC_DOUBLE, [d%patch, d%soil_layer],                                   &
                   'soil volumetric internal energy [J/m3] (prognostic)')
      call def_var(d, 'soil_temp_layer_pa',   NC_DOUBLE, [d%patch, d%soil_layer], 'soil temperature [K] (diagnosed twin)')
      call def_var(d, 'soil_fliq_layer_pa',   NC_DOUBLE, [d%patch, d%soil_layer], 'soil liquid-water fraction [-] (diagnosed twin)')
      call def_var(d, 'soil_w_surface_pa',      NC_DOUBLE, [d%patch], 'ponded surface water [kg/m2]')
      call def_var(d, 'soil_w_surface_enth_pa', NC_DOUBLE, [d%patch], 'ponded-water internal energy [J/m2]')
      !--- 7. Soil carbon (docs/science/soil_carbon.md §2).
      call def_var(d, 'soilc_fast_grnd_pa',   NC_DOUBLE, [d%patch], 'fast/metabolic litter carbon, above-ground [kgC/m2]')
      call def_var(d, 'soilc_fast_soil_pa',   NC_DOUBLE, [d%patch], 'fast/metabolic litter carbon, below-ground [kgC/m2]')
      call def_var(d, 'soilc_struct_grnd_pa', NC_DOUBLE, [d%patch], 'structural litter + CWD carbon, above-ground [kgC/m2]')
      call def_var(d, 'soilc_struct_soil_pa', NC_DOUBLE, [d%patch], 'structural litter + CWD carbon, below-ground [kgC/m2]')
      call def_var(d, 'soilc_microbial_pa',   NC_DOUBLE, [d%patch], 'microbial SOM carbon [kgC/m2]')
      call def_var(d, 'soilc_slow_pa',        NC_DOUBLE, [d%patch], 'slow/humified SOM carbon [kgC/m2]')
      call def_var(d, 'soilc_passive_pa',     NC_DOUBLE, [d%patch], 'passive SOM carbon [kgC/m2]')
      call def_var(d, 'soilc_lignin_grnd_pa', NC_DOUBLE, [d%patch], 'structural litter lignin, above-ground [kgC/m2]')
      call def_var(d, 'soilc_lignin_soil_pa', NC_DOUBLE, [d%patch], 'structural litter lignin, below-ground [kgC/m2]')
      !--- 8. The daily slow-to-fast rates: the turnover shed water and the slow loop's CO2.
      call def_var(d, 'shed_water_rate_pa', NC_DOUBLE, [d%patch], 'turnover shed water to the ground [kg/m2/s]')
      call def_var(d, 'shed_water_enth_pa', NC_DOUBLE, [d%patch], 'the enthalpy that water carries [W/m2]')
      call def_var(d, 'slow_co2_rate_pa',   NC_DOUBLE, [d%patch], 'slow-loop CO2 flux the fast NEE adds [umol/m2/s]')
      !--- 10. Solver memory: the step the fast controller last accepted, its warm start.
      call def_var(d, 'adapt_dt_last_pa', NC_DOUBLE, [d%patch], 'fast-controller warm-start step [s] (0 = cold)')
   end subroutine define_patch_vars

   !---~---
   !   The patches. Each per-patch number of a store goes through a contiguous buffer: a
   ! component section is never an actual argument (docs/building.md, "Compilers").
   !---~---
   subroutine put_patch_vars(r, poly)
      type(file_row_t), intent(in) :: r
      type(polygon_t),  intent(in) :: poly
      real(wp)    :: buf(poly%patch%n)
      integer(ik) :: ibuf(poly%patch%n), npat, npft, ip, kp
      npat = poly%patch%n
      npft = poly%n_pft
      associate (p => poly%patch)
         !--- 1. Identity and area.
         call put_int(r, 'global_patch_id', p%global_id(1:npat))
         call put_int(r, 'dist_type_pa',       p%dist_type(1:npat))
         call put_int(r, 'owner_site',      p%owner_site(1:npat))
         call put_dbl(r, 'area_pa',      p%area(1:npat))
         call put_dbl(r, 'age_pa',       p%age(1:npat))
         !--- 3. Recruit bank.
         do ip = 1_ik, npat
            call put_dbl(r, 'recruit_pool_pa', p%recruit_pool(1:npft, ip), item=ip)
         end do
         !--- 4. Canopy air.
         buf = p%cas(1:npat)%can_enthalpy
         call put_dbl(r, 'cas_can_enthalpy_pa', buf)
         buf = p%cas(1:npat)%can_shv
         call put_dbl(r, 'cas_can_shv_pa', buf)
         buf = p%cas(1:npat)%can_co2
         call put_dbl(r, 'cas_can_co2_pa', buf)
         buf = p%cas(1:npat)%can_temp
         call put_dbl(r, 'cas_can_temp_pa', buf)
         !--- 5. Snow.
         ibuf = p%snow(1:npat)%nlayer
         call put_int(r, 'snow_nlayer_pa', ibuf)
         do ip = 1_ik, npat
            call put_dbl(r, 'snow_swe_layer_pa',    p%snow(ip)%swe(1:n_snow_layer_max),         item=ip)
            call put_dbl(r, 'snow_energy_layer_pa', p%snow(ip)%snow_energy(1:n_snow_layer_max), item=ip)
            call put_dbl(r, 'snow_depth_layer_pa',  p%snow(ip)%snow_depth(1:n_snow_layer_max),  item=ip)
            call put_dbl(r, 'snow_temp_layer_pa',   p%snow(ip)%snow_temp(1:n_snow_layer_max),   item=ip)
            call put_dbl(r, 'snow_fliq_layer_pa',   p%snow(ip)%snow_fliq(1:n_snow_layer_max),   item=ip)
         end do
         !--- 6. Soil water and heat, and the pond.
         do ip = 1_ik, npat
            call put_dbl(r, 'soil_theta_layer_pa',  p%soil_water(ip)%theta(1:n_soil_layer_max),        item=ip)
            call put_dbl(r, 'soil_energy_layer_pa', p%soil_energy(ip)%soil_energy(1:n_soil_layer_max), item=ip)
            call put_dbl(r, 'soil_temp_layer_pa',   p%soil_energy(ip)%soil_temp(1:n_soil_layer_max),   item=ip)
            call put_dbl(r, 'soil_fliq_layer_pa',   p%soil_energy(ip)%soil_fliq(1:n_soil_layer_max),   item=ip)
         end do
         buf = p%soil_water(1:npat)%w_surface
         call put_dbl(r, 'soil_w_surface_pa', buf)
         buf = p%soil_water(1:npat)%w_surface_energy
         call put_dbl(r, 'soil_w_surface_enth_pa', buf)
         !--- 7. Soil carbon: each pool by its name, then the structural pools' lignin.
         do kp = 1_ik, N_SOIL_POOL
            buf = p%soil_carbon(1:npat)%carbon(kp)
            call put_dbl(r, 'soilc_'//trim(SOIL_POOL_NAME(kp))//'_pa', buf)
         end do
         buf = p%soil_carbon(1:npat)%lignin(1)
         call put_dbl(r, 'soilc_lignin_grnd_pa', buf)
         buf = p%soil_carbon(1:npat)%lignin(2)
         call put_dbl(r, 'soilc_lignin_soil_pa', buf)
         !--- 8. Daily slow-to-fast rates.
         call put_dbl(r, 'shed_water_rate_pa', p%shed_water_rate(1:npat))
         call put_dbl(r, 'shed_water_enth_pa', p%shed_water_enth(1:npat))
         call put_dbl(r, 'slow_co2_rate_pa',   p%slow_co2_rate(1:npat))
         !--- 10. Solver memory.
         call put_dbl(r, 'adapt_dt_last_pa', p%adapt_dt_last(1:npat))
      end associate
   end subroutine put_patch_vars

   !---~---
   !   The patches as the file holds them, each per-patch number of a store through a contiguous
   ! buffer.
   !---~---
   subroutine get_patch_vars(r, poly)
      type(file_row_t), intent(in)    :: r
      type(polygon_t),  intent(inout) :: poly
      real(wp)    :: buf(poly%patch%n)
      integer(ik) :: ibuf(poly%patch%n), npat, npft, ip, kp
      npat = poly%patch%n
      npft = poly%n_pft
      associate (p => poly%patch)
         !--- 1. Identity and area.
         call get_int(r, 'global_patch_id', p%global_id(1:npat))
         call get_int(r, 'dist_type_pa',       p%dist_type(1:npat))
         call get_int(r, 'owner_site',      p%owner_site(1:npat))
         call get_dbl(r, 'area_pa',      p%area(1:npat))
         call get_dbl(r, 'age_pa',       p%age(1:npat))
         !--- 3. Recruit bank.
         do ip = 1_ik, npat
            call get_dbl(r, 'recruit_pool_pa', p%recruit_pool(1:npft, ip), item=ip)
         end do
         !--- 4. Canopy air.
         call get_dbl(r, 'cas_can_enthalpy_pa', buf)
         p%cas(1:npat)%can_enthalpy = buf
         call get_dbl(r, 'cas_can_shv_pa', buf)
         p%cas(1:npat)%can_shv = buf
         call get_dbl(r, 'cas_can_co2_pa', buf)
         p%cas(1:npat)%can_co2 = buf
         call get_dbl(r, 'cas_can_temp_pa', buf)
         p%cas(1:npat)%can_temp = buf
         !--- 5. Snow.
         call get_int(r, 'snow_nlayer_pa', ibuf)
         p%snow(1:npat)%nlayer = ibuf
         do ip = 1_ik, npat
            call get_dbl(r, 'snow_swe_layer_pa',    p%snow(ip)%swe(1:n_snow_layer_max),         item=ip)
            call get_dbl(r, 'snow_energy_layer_pa', p%snow(ip)%snow_energy(1:n_snow_layer_max), item=ip)
            call get_dbl(r, 'snow_depth_layer_pa',  p%snow(ip)%snow_depth(1:n_snow_layer_max),  item=ip)
            call get_dbl(r, 'snow_temp_layer_pa',   p%snow(ip)%snow_temp(1:n_snow_layer_max),   item=ip)
            call get_dbl(r, 'snow_fliq_layer_pa',   p%snow(ip)%snow_fliq(1:n_snow_layer_max),   item=ip)
         end do
         !--- 6. Soil water and heat, and the pond.
         do ip = 1_ik, npat
            call get_dbl(r, 'soil_theta_layer_pa',  p%soil_water(ip)%theta(1:n_soil_layer_max),        item=ip)
            call get_dbl(r, 'soil_energy_layer_pa', p%soil_energy(ip)%soil_energy(1:n_soil_layer_max), item=ip)
            call get_dbl(r, 'soil_temp_layer_pa',   p%soil_energy(ip)%soil_temp(1:n_soil_layer_max),   item=ip)
            call get_dbl(r, 'soil_fliq_layer_pa',   p%soil_energy(ip)%soil_fliq(1:n_soil_layer_max),   item=ip)
         end do
         call get_dbl(r, 'soil_w_surface_pa', buf)
         p%soil_water(1:npat)%w_surface = buf
         call get_dbl(r, 'soil_w_surface_enth_pa', buf)
         p%soil_water(1:npat)%w_surface_energy = buf
         !--- 7. Soil carbon: each pool by its name, then the structural pools' lignin.
         do kp = 1_ik, N_SOIL_POOL
            call get_dbl(r, 'soilc_'//trim(SOIL_POOL_NAME(kp))//'_pa', buf)
            p%soil_carbon(1:npat)%carbon(kp) = buf
         end do
         call get_dbl(r, 'soilc_lignin_grnd_pa', buf)
         p%soil_carbon(1:npat)%lignin(1) = buf
         call get_dbl(r, 'soilc_lignin_soil_pa', buf)
         p%soil_carbon(1:npat)%lignin(2) = buf
         !--- 8. Daily slow-to-fast rates.
         call get_dbl(r, 'shed_water_rate_pa', p%shed_water_rate(1:npat))
         call get_dbl(r, 'shed_water_enth_pa', p%shed_water_enth(1:npat))
         call get_dbl(r, 'slow_co2_rate_pa',   p%slow_co2_rate(1:npat))
         !--- 10. Solver memory.
         call get_dbl(r, 'adapt_dt_last_pa', p%adapt_dt_last(1:npat))
      end associate
   end subroutine get_patch_vars

   !---~---
   !   4. AFTER THE READ: the geometry from the pools, the growth rings of another window
   !    length, and the maps and sort order
   !---~---

   subroutine derive_after_read(poly, cfg, same_window)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      logical,             intent(in)    :: same_window
      integer(ik) :: ncoh, nwin, i
      ncoh = poly%cohort%n

      !--- The geometry from the pools, as the run derives it.
      do i = 1_ik, ncoh
         call set_cohort_size_from_carbon(poly%cohort, i)
      end do
      !---~---
      !   The growth rings, when this run's window is another length than the file's: seeded
      ! full of the saved growth_avg.
      !---~---
      if (.not. same_window) then
         nwin = growth_window_steps(cfg)
         do i = 1_ik, ncoh
            if (.not. is_unset_positive(poly%cohort%growth_avg(i))) then
               poly%cohort%growth_count(i)   = nwin
               poly%cohort%growth_accum(i)   = poly%cohort%growth_avg(i) * real(nwin, wp)
               poly%cohort%growth_hist(:, i) = poly%cohort%growth_avg(i)
            end if
         end do
      end if

      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine derive_after_read

end module ed_init_restart
