! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_state_file -- the state (restart) file of ed_init_restart
! (docs/science/order_of_processes.md §6, The restart file):
!   1. a polygon of two sites, three patches and four cohorts, every store stamped with a
!      value of its own, is written (A), read back and written again (B): every variable of
!      B is A's, bit for bit;
!   2. every variable of A differs from what a reader that skipped it would hold (Z, the
!      polygon the reader allocates before it reads a store), so (1) cannot pass on a store
!      left at its default. The few that cannot differ are the polygon's structure and what
!      the configuration gives, listed in MADE_BY_THE_RUN;
!   3. the restructuring owed is recorded: a year's, a month's, none;
!   4. a run of another growth-window length seeds each cohort's rings full of its growth_avg;
!   5. the reader refuses a file it cannot start from, and says why; and a run whose polygons
!      the file does not each hold a row for is refused (restart_file_slots, as region_open
!      matches them): a site run (polygon 0) a region's file, a region a site's file;
!   6. a file from before the ground's surface layers were one record: its snow pack's one layer
!      is the water layer, its pond joins that layer, its organic_* variables are the porous
!      layers, and a density it held is not read; without organic_* variables the porous layers
!      come back empty, and a run that starts from such a file with a moss_depth makes its moss.
! A store a restart drops comes back at its default: stomata reopened mid-drought, phenology
! reborn at flush 1 and shed 0, soil bone dry at 0 K. Section 1 catches it for every variable
! the file holds, and section 2 makes sure a new variable is stamped here before it counts.
!---~---
program test_state_file
   use iso_c_binding,           only : c_int, c_size_t, c_double, c_char, c_null_char
   use iso_fortran_env,         only : int64
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : n_soil_layer_max
   use ed_biogeochem_params,    only : N_SOIL_POOL
   use ed_config,               only : ed_config_t, growth_window_steps, INIT_RESTART
   use ed_therm_lib,            only : temp_to_internal_energy, internal_energy_liquid
   use ed_ground_surface_params, only : N_POROUS, POROUS_MOSS
   use ed_init,                 only : init_polygon
   use ed_polygon_state_vars,   only : polygon_t, alloc_polygon
   use ed_site_state_vars,      only : make_site_surface
   use ed_cohort_state_vars,    only : hold_cohort_tissue_temps, set_cohort_size_from_carbon
   use ed_init_bare_ground,     only : init_bare_sites, add_cohort, finalize_init
   use ed_init_restart,         only : restart_file_t, restart_file_open, restart_file_write, restart_file_close, &
                                       read_restart_file, restart_file_slots
   use ed_time,                 only : ed_time_t
   use ed_netcdf_c,             only : NC_NOERR, NC_NOWRITE, nc_open_f, nc_close, nc_inq_varid_f, nc_inq_dimlen,  &
                                       nc_get_vara_double, NC_DOUBLE, nc_def_dim_f, nc_def_var_f, nc_enddef,         &
                                       nc_put_vara_double, nc_inq_dimid, cstr
   use ed_test_support,         only : build_test_config
   use ed_test_assert,          only : check, check_true, check_close, test_report
   implicit none

   interface
      integer(c_int) function nc_inq_nvars(ncid, nvarsp) bind(c, name='nc_inq_nvars')
         import :: c_int
         integer(c_int), value, intent(in)  :: ncid
         integer(c_int),        intent(out) :: nvarsp
      end function nc_inq_nvars
      integer(c_int) function nc_inq_var(ncid, varid, name, xtypep, ndimsp, dimidsp, nattsp) bind(c, name='nc_inq_var')
         import :: c_int, c_char
         integer(c_int), value,  intent(in)  :: ncid, varid
         character(kind=c_char), intent(out) :: name(*)
         integer(c_int),         intent(out) :: xtypep, ndimsp, dimidsp(*), nattsp
      end function nc_inq_var
      integer(c_int) function nc_inq_attname(ncid, varid, attnum, name) bind(c, name='nc_inq_attname')
         import :: c_int, c_char
         integer(c_int), value,  intent(in)  :: ncid, varid, attnum
         character(kind=c_char), intent(out) :: name(*)
      end function nc_inq_attname
      integer(c_int) function nc_redef(ncid) bind(c, name='nc_redef')
         import :: c_int
         integer(c_int), value, intent(in) :: ncid
      end function nc_redef
      integer(c_int) function nc_rename_var(ncid, varid, name) bind(c, name='nc_rename_var')
         import :: c_int, c_char
         integer(c_int), value,  intent(in) :: ncid, varid
         character(kind=c_char), intent(in) :: name(*)
      end function nc_rename_var
   end interface

   integer(c_int), parameter :: NC_WRITE = 1_c_int   ! netCDF's open mode for writing

   !---~---
   !   The variables a reader that skipped them would still hold: the polygon's structure,
   ! which the reader allocates to, the file's time and the polygon axis, and what each
   ! site's surface takes from the configuration (it is written for readers of the file).
   !---~---
   character(len=*), parameter :: MADE_BY_THE_RUN(*) = [character(len=16) ::                        &
      'polygon_id', 'lat', 'lon', 'n_cohort', 'n_patch', 'n_site', 'n_window',                      &
      'meta_real', 'elevation_si', 'slope_si', 'aspect_si', 'n_soil_si', 'soil_depth_si',             &
      'bottom_bc_si', 'drain_factor_si', 'deep_temp_si']
   character(len=*), parameter :: FILE_A = './test_sf_a-S-20000101000000.nc'
   character(len=*), parameter :: FILE_B = './test_sf_b-S-20000101000000.nc'
   character(len=*), parameter :: FILE_Z = './test_sf_z-S-20000101000000.nc'
   character(len=*), parameter :: FILE_R = './test_sf_r-S-20000101000000.nc'   ! a region's: polygons 35, 18, 53
   character(len=*), parameter :: FILE_O = './test_sf_o-S-20000101000000.nc'   ! an older file, organic_* and all
   character(len=*), parameter :: FILE_N = './test_sf_n-S-20000101000000.nc'   ! an older file with no organic_*

   type(ed_config_t)   :: cfg
   type(polygon_t)     :: poly, back
   type(ed_time_t)     :: now, rt
   logical             :: found

   cfg = build_test_config()
   now = ed_time_t(year=2000, month=1, day=1, hour=0, minute=0, second=0)
   call stamped_polygon(poly)

   call test_write_read_write()
   call test_not_at_defaults()
   call test_restructuring_owed()
   call test_other_window()
   call test_refusals()
   call test_older_file()

   call test_report('test_state_file')

contains

   !--- 1. A, read back and written again, is B, variable for variable.
   subroutine test_write_read_write()
      character(len=64), allocatable :: names(:)
      logical,           allocatable :: same(:)
      integer :: k
      print '(a)', '-- 1. written, read back and written again --'
      call write_state(poly, 'test_sf_a', .true., .true.)
      call read_restart_file(back, cfg, FILE_A, rt, found, polygon=1_ik)
      call check_true('the file is read back', found)
      call check_true('with its cohorts, patches and sites', back%cohort%n == poly%cohort%n .and.           &
                      back%patch%n == poly%patch%n .and. back%site%n == poly%site%n)
      call write_state(back, 'test_sf_b', .true., .true.)
      call compare_files(FILE_A, FILE_B, names, same)
      call check_true('the file holds variables', size(names) > 50)
      do k = 1, size(names)
         call check_true(trim(names(k))//' is read back bit for bit', same(k))
      end do
   end subroutine test_write_read_write

   !--- 2. A reader that skipped a variable would hold Z's: A must differ from it.
   subroutine test_not_at_defaults()
      type(polygon_t) :: z
      character(len=64), allocatable :: names(:)
      logical,           allocatable :: same(:)
      integer(ik) :: is
      integer     :: k
      print '(a)', '-- 2. every variable away from the reader''s default --'
      !--- What get_polygon_vars makes before it reads a store.
      call alloc_polygon(z, cfg%pft%n, poly%cohort%n, poly%patch%n, growth_window_steps(cfg), cfg%pft,      &
                         n_site=poly%site%n)
      z%cohort%n = poly%cohort%n
      z%patch%n  = poly%patch%n
      do is = 1_ik, z%site%n
         call make_site_surface(z%site%surface(is), cfg)
      end do
      call write_state(z, 'test_sf_z', .true., .true.)
      call compare_files(FILE_A, FILE_Z, names, same)
      do k = 1, size(names)
         if (any(names(k) == MADE_BY_THE_RUN)) cycle
         call check_true(trim(names(k))//' is stamped away from its default', .not. same(k))
      end do
   end subroutine test_not_at_defaults

   !--- 3. The restructuring owed: the driver's yearly checkpoint is written with it pending.
   subroutine test_restructuring_owed()
      print '(a)', '-- 3. the restructuring owed --'
      call read_restart_file(back, cfg, FILE_A, rt, found, polygon=1_ik)
      call check_true('a year''s restructuring owed reads so', found .and. back%restructure_pending .and.  &
                      back%restructure_new_year)
      call write_state(poly, 'test_sf_a', .true., .false.)
      call read_restart_file(back, cfg, FILE_A, rt, found, polygon=1_ik)
      call check_true('a month''s reads so', found .and. back%restructure_pending .and.                  &
                      .not. back%restructure_new_year)
      call write_state(poly, 'test_sf_a', .false., .false.)
      call read_restart_file(back, cfg, FILE_A, rt, found, polygon=1_ik)
      call check_true('none owed reads so', found .and. .not. back%restructure_pending .and.             &
                      .not. back%restructure_new_year)
      call check_true('the file''s time is read', rt%year == 2000 .and. rt%month == 1 .and. rt%day == 1)
   end subroutine test_restructuring_owed

   !--- 4. A run of another growth-window length: each cohort's rings full of its average.
   subroutine test_other_window()
      type(ed_config_t) :: cfg_short
      integer(ik) :: nwin, i
      print '(a)', '-- 4. another growth-window length --'
      cfg_short = cfg
      cfg_short%demography%growth_memory_days = cfg%demography%growth_memory_days / 3.0_wp
      nwin = growth_window_steps(cfg_short)
      call check_true('the window is shorter', nwin < growth_window_steps(cfg), real(nwin, wp))
      call read_restart_file(back, cfg_short, FILE_A, rt, found, polygon=1_ik)
      call check_true('the file is read', found)
      call check_true('the rings take the run''s length', size(back%cohort%growth_hist, 1) == nwin)
      do i = 1_ik, back%cohort%n
         call check_true('the ring is full', back%cohort%growth_count(i) == nwin)
         call check_true('of the cohort''s growth_avg', all(back%cohort%growth_hist(1:nwin, i) ==           &
                         back%cohort%growth_avg(i)))
         call check_true('and growth_accum is their sum', abs(back%cohort%growth_accum(i) -                &
                         back%cohort%growth_avg(i) * real(nwin, wp)) <= 1.0e-12_wp * abs(back%cohort%growth_accum(i)))
      end do
   end subroutine test_other_window

   !--- 5. The files a run cannot start from.
   subroutine test_refusals()
      type(ed_config_t) :: cfg_pft
      character(len=:), allocatable :: why
      integer(ik), allocatable :: slot(:)
      print '(a)', '-- 5. the files a run cannot start from --'
      call read_restart_file(back, cfg, './no_such_state.nc', rt, found, polygon=1_ik, why=why)
      call check_true('a missing file: not found, "cannot be opened"', .not. found .and.                  &
                      index(why, 'cannot be opened') > 0)
      cfg_pft = cfg
      cfg_pft%pft%n = cfg%pft%n - 1_ik
      call read_restart_file(back, cfg_pft, FILE_A, rt, found, polygon=1_ik, why=why)
      call check_true('another PFT count: refused, naming both', .not. found .and.                        &
                      index(why, 'holds 3 PFTs and the configuration 2') > 0)
      call read_restart_file(back, cfg, FILE_A, rt, found, polygon=1_ik, site_cfg=[cfg], why=why)
      call check_true('another site count: refused, naming both', .not. found .and.                       &
                      index(why, 'holds 2 sites and the configuration 1') > 0)
      !--- A region's file in a site run, and a site's file (FILE_A, polygon 0) in a region.
      call write_region_state('test_sf_r', [35_ik, 18_ik, 53_ik])
      call restart_file_slots(FILE_R, [0_ik], slot, why)
      call check_true('a site run refuses a region state file', index(why, 'does not hold polygon 0') > 0)
      call restart_file_slots(FILE_A, [35_ik, 18_ik, 53_ik], slot, why)
      call check_true('a region refuses a site state file', index(why, 'does not hold polygon 35') > 0)
   end subroutine test_refusals

   !---~---
   !   6. A's ground-surface variables renamed away and an older file's put in their place: each
   ! patch's snow pack as its one layer, its pond apart, its organic_* variables (or none) and a
   ! density no reader takes. The pack is 10 + ip kg/m2 at 265 K at the first pack's density, the
   ! pond ip / 2 kg/m2 at 274 K, and the organic_* values each patch's own.
   !---~---
   subroutine test_older_file()
      real(wp), parameter :: MOSS_DEPTH = 0.05_wp   ! [m]
      real(c_double)    :: snow(3, 5), pond(3, 2), organic(N_POROUS, 3, 3)
      type(ed_config_t) :: cfg_r
      type(polygon_t)   :: started
      type(ed_time_t)   :: now_r
      logical           :: ok
      integer(ik)       :: ip, k
      character(len=2)  :: tag
      print '(a)', '-- 6. a state file from before the ground''s surface layers were one record --'
      do ip = 1_ik, 3_ik
         snow(ip, 1) = 10.0_wp + real(ip, wp)
         snow(ip, 2) = temp_to_internal_energy(0.0_wp, snow(ip, 1), 265.0_wp, 0.0_wp)
         snow(ip, 3) = snow(ip, 1) / 250.0_wp
         snow(ip, 4) = 265.0_wp
         snow(ip, 5) = 0.0_wp
         pond(ip, 1) = 0.5_wp * real(ip, wp)
         pond(ip, 2) = pond(ip, 1) * internal_energy_liquid(274.0_wp)
         organic(:, 1, ip) = [1.1_wp, 2.0_wp] + 0.1_wp * real(ip, wp)
         organic(:, 2, ip) = [3.0_wp, 15.0_wp] + real(ip, wp)
         organic(:, 3, ip) = [1.0e6_wp, 5.0e6_wp] * real(1_ik + ip, wp)
      end do
      call write_older_state('test_sf_o', FILE_O, snow, pond, organic, .true.)
      call read_restart_file(back, cfg, FILE_O, rt, found, polygon=1_ik)
      call check_true('older file: it is read', found .and. back%patch%n == 3_ik)
      do ip = 1_ik, 3_ik
         write(tag, '(i0)') ip
         associate (gs => back%patch%ground_surface(ip))
            call check_close('older file, patch '//trim(tag)//': the water layer holds the pack and the pond',       &
                             gs%water%water, snow(ip, 1) + pond(ip, 1), 1.0e-12_wp)
            call check_close('older file, patch '//trim(tag)//': and their energy', gs%water%energy,                 &
                             snow(ip, 2) + pond(ip, 2), 1.0e-9_wp * abs(snow(ip, 2)))
            call check_close('older file, patch '//trim(tag)//': as deep as the pack', gs%water%matrix_depth,         &
                             snow(ip, 3), 0.0_wp)
            do k = 1_ik, N_POROUS
               call check_close('older file, patch '//trim(tag)//': a porous layer''s matter, the file''s organic', &
                                gs%porous(k)%matter, organic(k, 1, ip), 0.0_wp)
               call check_close('older file, patch '//trim(tag)//': its water', gs%porous(k)%water,                  &
                                organic(k, 2, ip), 0.0_wp)
               call check_close('older file, patch '//trim(tag)//': its energy', gs%porous(k)%energy,                &
                                organic(k, 3, ip), 0.0_wp)
            end do
         end associate
      end do
      !--- No organic_* variables: the porous layers come back empty.
      call write_older_state('test_sf_n', FILE_N, snow, pond, organic, .false.)
      call read_restart_file(back, cfg, FILE_N, rt, found, polygon=1_ik)
      call check_true('older file without organic_* variables: it is read', found)
      do ip = 1_ik, 3_ik
         do k = 1_ik, N_POROUS
            associate (layer => back%patch%ground_surface(ip)%porous(k))
               call check_true('older file without organic_* variables: the porous layers are empty',              &
                               layer%matter == 0.0_wp .and. layer%water == 0.0_wp .and. layer%energy == 0.0_wp)
            end associate
         end do
      end do
      !--- A run that starts from it with a moss_depth makes its moss on every patch.
      cfg_r = cfg
      cfg_r%init%init_mode    = INIT_RESTART
      cfg_r%init%restart_file = FILE_N
      cfg_r%soil%moss_depth   = MOSS_DEPTH
      now_r = now
      call init_polygon(started, cfg_r, 'an older state file', .false., 1_ik, now_r, ok, 0_ik, 0.0_wp, 0.0_wp,       &
                        0.0_wp)
      call check_true('a run starts from an older file', ok .and. started%patch%n == 3_ik)
      do ip = 1_ik, started%patch%n
         associate (moss => started%patch%ground_surface(ip)%porous(POROUS_MOSS))
            call check_close('a run from an older file makes the moss its moss_depth gives', moss%matter,         &
                             MOSS_DEPTH * cfg%ground_surface%porous(POROUS_MOSS)%bulk_density, 1.0e-12_wp)
            call check_true('and the moss holds water', moss%water > 0.0_wp, moss%water)
         end associate
      end do
   end subroutine test_older_file

   !---~---
   !   The fixture's state as an older file held it, at path, written as prefix: the ground-surface
   ! variables renamed to names no reader takes, the porous layers' to the older file's organic_*
   ! names when with_organic, and the snow pack's layer (snow: water, energy, depth, temperature
   ! and liquid fraction by patch), the pond (pond: water and energy by patch) and the density
   ! the older file held (organic_density_pa) defined and put.
   !---~---
   subroutine write_older_state(prefix, path, snow, pond, organic, with_organic)
      character(len=*), intent(in) :: prefix, path
      real(c_double),   intent(in) :: snow(:,:), pond(:,:), organic(:,:,:)
      logical,          intent(in) :: with_organic
      character(len=*), parameter :: LAYER(5)  = [character(len=12) :: 'water', 'energy', 'matrix_depth', 'temp', &
                                                  'fliq']
      character(len=*), parameter :: SNOW_VAR(5) = [character(len=6) :: 'swe', 'energy', 'depth', 'temp', 'fliq']
      character(len=*), parameter :: POROUS(3)  = [character(len=6) :: 'matter', 'water', 'energy']
      character(len=*), parameter :: ORGANIC_VAR(3) = [character(len=6) :: 'mass', 'water', 'energy']
      real(c_double)    :: density(3), values(N_POROUS * 3)
      integer(c_int)    :: ncid, d_poly, d_patch, d_snow, v_snow(5), v_pond(2), v_dens, varid
      integer(c_size_t) :: n_patch
      integer           :: k
      n_patch = int(size(pond, 1), c_size_t)
      call write_state(poly, prefix, .false., .false.)
      !--- The renames first, on their own: a netCDF-4 file takes them in data mode.
      call nc_ok(nc_open_f(path, NC_WRITE, ncid), 'open '//path)
      call read_every_attribute(ncid)
      do k = 1, size(LAYER)
         call rename_var(ncid, 'ground_water_layer_'//trim(LAYER(k))//'_pa', 'unread_water_layer_'//trim(LAYER(k))//'_pa')
      end do
      do k = 1, size(POROUS)
         if (with_organic) then
            call rename_var(ncid, 'ground_porous_layer_'//trim(POROUS(k))//'_pa', 'organic_'//trim(ORGANIC_VAR(k))//'_pa')
         else
            call rename_var(ncid, 'ground_porous_layer_'//trim(POROUS(k))//'_pa', 'unread_porous_'//trim(POROUS(k))//'_pa')
         end if
      end do
      call nc_ok(nc_close(ncid), 'close after the renames')
      !--- Then the older file's variables, defined and put.
      call nc_ok(nc_open_f(path, NC_WRITE, ncid), 'reopen '//path)
      call nc_ok(nc_redef(ncid), 'redef')
      d_poly  = dim_id(ncid, 'polygon')
      d_patch = dim_id(ncid, 'patch')
      call nc_ok(nc_def_dim_f(ncid, 'snow_layer', 1_c_size_t, d_snow), 'def snow_layer')
      do k = 1, size(SNOW_VAR)
         call nc_ok(nc_def_var_f(ncid, 'snow_'//trim(SNOW_VAR(k))//'_layer_pa', NC_DOUBLE, 3_c_int,                &
                                 [d_poly, d_patch, d_snow], v_snow(k)), 'def snow_'//trim(SNOW_VAR(k)))
      end do
      call nc_ok(nc_def_var_f(ncid, 'soil_w_surface_pa', NC_DOUBLE, 2_c_int, [d_poly, d_patch], v_pond(1)), 'def pond')
      call nc_ok(nc_def_var_f(ncid, 'soil_w_surface_enth_pa', NC_DOUBLE, 2_c_int, [d_poly, d_patch], v_pond(2)),       &
                 'def pond energy')
      call nc_ok(nc_def_var_f(ncid, 'organic_density_pa', NC_DOUBLE, 2_c_int, [d_poly, d_patch], v_dens),             &
                 'def density')
      call nc_ok(nc_enddef(ncid), 'enddef')
      do k = 1, size(SNOW_VAR)
         call nc_ok(nc_put_vara_double(ncid, v_snow(k), [0_c_size_t, 0_c_size_t, 0_c_size_t],                     &
                                       [1_c_size_t, n_patch, 1_c_size_t], snow(:, k)), 'put snow')
      end do
      do k = 1, 2
         call nc_ok(nc_put_vara_double(ncid, v_pond(k), [0_c_size_t, 0_c_size_t], [1_c_size_t, n_patch], pond(:, k)), &
                    'put pond')
      end do
      density = 999.0_wp
      call nc_ok(nc_put_vara_double(ncid, v_dens, [0_c_size_t, 0_c_size_t], [1_c_size_t, n_patch], density),        &
                 'put density')
      if (with_organic) then
         do k = 1, size(ORGANIC_VAR)
            call nc_ok(nc_inq_varid_f(ncid, 'organic_'//trim(ORGANIC_VAR(k))//'_pa', varid), 'inq organic')
            values = reshape(organic(:, k, :), [N_POROUS * 3])
            call nc_ok(nc_put_vara_double(ncid, varid, [0_c_size_t, 0_c_size_t, 0_c_size_t],                      &
                                          [1_c_size_t, n_patch, int(N_POROUS, c_size_t)], values), 'put organic')
         end do
      end if
      call nc_ok(nc_close(ncid), 'close '//path)
   end subroutine write_older_state

   !---~---
   !   Every variable's attributes listed, as a reader that opens a file whole lists them. The
   ! netCDF library reads a netCDF-4 variable's attributes only when asked, and renaming a
   ! variable before they are read makes the file's close fail.
   !---~---
   subroutine read_every_attribute(ncid)
      integer(c_int), intent(in) :: ncid
      character(kind=c_char) :: name(257)
      integer(c_int) :: nvar, varid, xtype, ndims, dimids(8), natts, att
      call nc_ok(nc_inq_nvars(ncid, nvar), 'inq nvars')
      do varid = 0_c_int, nvar - 1_c_int
         call nc_ok(nc_inq_var(ncid, varid, name, xtype, ndims, dimids, natts), 'inq var')
         do att = 0_c_int, natts - 1_c_int
            call nc_ok(nc_inq_attname(ncid, varid, att, name), 'inq attname')
         end do
      end do
   end subroutine read_every_attribute

   !--- Rename variable old_name to new_name in an open file.
   subroutine rename_var(ncid, old_name, new_name)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: old_name, new_name
      character(kind=c_char) :: cname(len_trim(new_name) + 1)
      integer(c_int) :: varid
      call nc_ok(nc_inq_varid_f(ncid, old_name, varid), 'inq '//old_name)
      cname = cstr(new_name)
      call nc_ok(nc_rename_var(ncid, varid, cname), 'rename '//old_name)
   end subroutine rename_var

   !--- A dimension's id in an open file.
   integer(c_int) function dim_id(ncid, name) result(id)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      call nc_ok(nc_inq_dimid(ncid, cname, id), 'inq dim '//name)
   end function dim_id

   !--- A netCDF call that must succeed.
   subroutine nc_ok(status, what)
      integer(c_int),   intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= NC_NOERR) print '(3a,i0)', ' older file: netCDF ', what, ' returned ', status
      call check(status == NC_NOERR, 'older file: netCDF '//what)
   end subroutine nc_ok

   !---~---
   !   The fixture: two sites (0.3 and 0.7 of the polygon) of one and two patches, four
   ! cohorts, and every store the file holds set to a value of its own, away from its
   ! default. The geometry is then derived from the pools, as the reader derives it.
   !---~---
   subroutine stamped_polygon(p)
      type(polygon_t), intent(out) :: p
      integer(ik) :: i, ip, k, nwin
      real(wp)    :: x
      call init_bare_sites(p, cfg, [1_ik, 2_ik], [cfg, cfg], [0.3_wp, 0.7_wp])
      call add_cohort(p, cfg, 1_ik, 1_ik, 0.20_wp, 20.0_wp)
      call add_cohort(p, cfg, 1_ik, 3_ik, 0.50_wp, 4.0_wp)
      call add_cohort(p, cfg, 2_ik, 2_ik, 0.10_wp, 35.0_wp)
      call add_cohort(p, cfg, 3_ik, 3_ik, 0.30_wp, 12.0_wp)
      call finalize_init(p)
      p%kt_last_day = 0.61_wp
      p%growth_hist_pos = 7_ik
      p%site%t_growth_avg(1:2) = [296.5_wp, 297.25_wp]
      nwin = int(size(p%cohort%growth_hist, 1), ik)
      do i = 1_ik, p%cohort%n
         x = 0.01_wp * real(i, wp)
         associate (c => p%cohort)
            c%overtopping_lai(i)     = 0.4_wp + x
            c%leaf_water_mass(i)     = 0.0123_wp + x
            c%wood_water_mass(i)     = 0.0456_wp + x
            c%leaf_surf_water(i)     = 0.020_wp + x
            c%wood_surf_water(i)     = 0.007_wp + x
            !--- Negative and distinct: both sentinels have the other sign.
            c%dmax_psi_leaf(i)       = -0.75_wp - x
            c%dmax_psi_leaf_accum(i) = -0.40_wp - x
            c%leaf_temp(i)           = 291.5_wp + real(i, wp)
            c%wood_temp(i)           = 289.5_wp + real(i, wp)
            c%sla(i)                 = 21.0_wp + x
            c%vcmax25(i)             = 42.0_wp + x
            c%rd25(i)                = 0.65_wp + x
            c%leaf_turnover_rate(i)  = 0.7_wp + x
            c%leaf_flush_tendency(i) = 0.30_wp + x
            c%leaf_shed_tendency(i)  = 0.70_wp - x
            c%growing_degree_days(i) = 123.0_wp + real(i, wp)
            c%cold_degree_days(i)    = 17.0_wp + real(i, wp)
            c%dry_psi_sum(i)         = 5.0_wp + real(i, wp)
            c%wet_psi_sum(i)         = 2.0_wp + real(i, wp)
            c%light_hours_mean(i)    = 10.0_wp + x
            c%growth_hist(1:nwin, i) = [(0.1_wp * real(k, wp) + x, k = 1_ik, nwin)]
            c%growth_count(i)        = nwin - i
            c%growth_accum(i)        = sum(c%growth_hist(1:nwin - i, i))
            c%growth_avg(i)          = c%growth_accum(i) / real(nwin - i, wp)
         end associate
         !--- The energies the stamped temperatures and waters give, and the geometry.
         call hold_cohort_tissue_temps(p%cohort, i)
         call set_cohort_size_from_carbon(p%cohort, i)
      end do
      do ip = 1_ik, p%patch%n
         x = 0.01_wp * real(ip, wp)
         associate (pa => p%patch)
            pa%dist_type(ip) = 2_ik
            pa%age(ip) = 37.0_wp + real(ip, wp)
            pa%recruit_pool(1:p%n_pft, ip) = [(1.0e-3_wp * real(k, wp) + x, k = 1_ik, p%n_pft)]
            pa%cas(ip)%can_enthalpy = 3.1e5_wp + 1.0e3_wp * x
            pa%cas(ip)%can_shv      = 0.011_wp + 1.0e-3_wp * x
            pa%cas(ip)%can_co2      = 415.0_wp + x
            pa%cas(ip)%can_temp     = 291.0_wp + x
            pa%cas(ip)%air_mass_capacity  = 24.0_wp + x
            pa%cas(ip)%air_molar_capacity = 820.0_wp + 10.0_wp * x
            pa%ground_surface(ip)%water%water        = 4.0_wp + x
            pa%ground_surface(ip)%water%energy       = -1.0e5_wp - x
            pa%ground_surface(ip)%water%matrix_depth = 0.02_wp + x
            pa%ground_surface(ip)%water%temp         = 270.0_wp + x
            pa%ground_surface(ip)%water%fliq         = 0.05_wp + x
            pa%ground_surface(ip)%porous%matter      = [1.3_wp + x, 2.1_wp + x]
            pa%ground_surface(ip)%porous%water       = [0.8_wp + x, 9.5_wp + x]
            pa%ground_surface(ip)%porous%energy      = [3.1e6_wp + x, 9.2e6_wp + x]
            pa%soil_water(ip)%theta(1:n_soil_layer_max) = [(0.20_wp + 0.005_wp * real(k, wp) + x, k = 1_ik, n_soil_layer_max)]
            pa%soil_energy(ip)%soil_energy(1:n_soil_layer_max) = [(2.0e8_wp + 1.0e6_wp * real(k, wp) + x,         &
                                                                   k = 1_ik, n_soil_layer_max)]
            pa%soil_energy(ip)%soil_temp(1:n_soil_layer_max) = [(284.0_wp + 0.1_wp * real(k, wp) + x,             &
                                                                 k = 1_ik, n_soil_layer_max)]
            pa%soil_energy(ip)%soil_fliq(1:n_soil_layer_max) = 0.9_wp + x
            !--- The daily slow-to-fast rates.
            pa%shed_water_rate(ip) = 3.0e-8_wp + 1.0e-9_wp * x
            pa%shed_water_enth(ip) = pa%shed_water_rate(ip) * 1.1e5_wp
            pa%slow_co2_rate(ip)   = 0.8_wp + x
            pa%adapt_dt_last(ip)   = 120.0_wp + x
            pa%soil_carbon(ip)%carbon(1:N_SOIL_POOL) = [(0.10_wp + 0.01_wp * real(k, wp) + x, k = 1_ik, N_SOIL_POOL)]
            pa%soil_carbon(ip)%lignin(1:2) = [0.03_wp + x, 0.04_wp + x]
         end associate
      end do
   end subroutine stamped_polygon

   !--- A polygon's state, as the one polygon of a state file (id 0, at 0 N 0 E).
   subroutine write_state(s, prefix, pending, new_year)
      type(polygon_t),  intent(in) :: s
      character(len=*), intent(in) :: prefix
      logical,          intent(in) :: pending, new_year
      type(restart_file_t) :: sf
      call restart_file_open(sf, '.', prefix, now, s%n_pft, int(size(s%cohort%growth_hist, 1), ik),            &
                             s%cohort%n, s%patch%n, s%site%n, polygon_id=[0_ik], polygon_lat=[0.0_wp],        &
                             polygon_lon=[0.0_wp], restructure_pending=pending, restructure_new_year=new_year)
      call restart_file_write(sf, s, now, 1_ik)
      call restart_file_close(sf)
   end subroutine write_state

   !--- A region's state file: the fixture's polygon in each row, the rows the polygons `ids`.
   subroutine write_region_state(prefix, ids)
      character(len=*), intent(in) :: prefix
      integer(ik),      intent(in) :: ids(:)
      type(restart_file_t) :: sf
      integer(ik) :: row
      call restart_file_open(sf, '.', prefix, now, poly%n_pft, int(size(poly%cohort%growth_hist, 1), ik),      &
                             poly%cohort%n, poly%patch%n, poly%site%n, polygon_id=ids,                       &
                             polygon_lat=[(0.0_wp, row = 1_ik, size(ids, kind=ik))],                         &
                             polygon_lon=[(0.0_wp, row = 1_ik, size(ids, kind=ik))])
      do row = 1_ik, size(ids, kind=ik)
         call restart_file_write(sf, poly, now, row)
      end do
      call restart_file_close(sf)
   end subroutine write_region_state

   !--- Every variable of file a, and whether file b holds the same values, bit for bit.
   subroutine compare_files(a, b, names, same)
      character(len=*),               intent(in)  :: a, b
      character(len=64), allocatable, intent(out) :: names(:)
      logical,           allocatable, intent(out) :: same(:)
      real(c_double), allocatable :: xa(:), xb(:)
      integer(c_int) :: nca, ncb, nvar, varid, vb, st
      integer        :: k
      st = nc_open_f(a, NC_NOWRITE, nca)
      st = nc_open_f(b, NC_NOWRITE, ncb)
      st = nc_inq_nvars(nca, nvar)
      allocate(names(nvar), same(nvar))
      do k = 1, nvar
         varid = int(k - 1, c_int)
         call read_var(nca, varid, names(k), xa)
         same(k) = .false.
         if (nc_inq_varid_f(ncb, trim(names(k)), vb) /= NC_NOERR) cycle
         call read_var(ncb, vb, names(k), xb)
         if (size(xa) /= size(xb)) cycle
         same(k) = all(transfer(xa, 0_int64, size(xa)) == transfer(xb, 0_int64, size(xb)))
      end do
      st = nc_close(nca)
      st = nc_close(ncb)
   end subroutine compare_files

   !--- One variable, whole, as doubles (an integer variable converts exactly).
   subroutine read_var(ncid, varid, name, x)
      integer(c_int),              intent(in)    :: ncid, varid
      character(len=64),           intent(inout) :: name
      real(c_double), allocatable, intent(out)   :: x(:)
      character(kind=c_char) :: cname(257)
      integer(c_int)    :: xtype, ndims, dimids(8), natts, st
      integer(c_size_t) :: start(8), count(8)
      integer           :: d, n
      st = nc_inq_var(ncid, varid, cname, xtype, ndims, dimids, natts)
      n = findloc(cname, c_null_char, dim=1) - 1
      name = transfer(cname(1:n), repeat(' ', n))
      do d = 1, ndims
         st = nc_inq_dimlen(ncid, dimids(d), count(d))
      end do
      start = 0_c_size_t
      allocate(x(product(count(1:ndims))))
      st = nc_get_vara_double(ncid, varid, start, count, x)
   end subroutine read_var

end program test_state_file
