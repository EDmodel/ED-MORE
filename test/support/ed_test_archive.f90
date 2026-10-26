! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_test_archive -- a tiny synthetic ED_ERA5land archive for the tests, written                 !
! through the netCDF C API: a 4 x 18 grid (10 deg x 20 deg) spanning the antimeridian and two     !
! 16-cell chunk columns wide, every month of 2021 and January 2022 (so a recycle window may start  !
! after January 1st), two no-data cells, a mostly-water column, and values that are analytic       !
! functions of cell and hour, so every check has a known answer. Shared by test_met_era5land (the  !
! reader) and test_region (a region against site runs). write_surface_maps adds the [maps] maps:  !
! soil depth, texture and chemistry on SoilGrids' six layers, and Chave's E, each a function of the  !
! cell (soil_depth_at, ..., texture_at), empty at the no-data cells.                                  !
!==========================================================================================!
module ed_test_archive
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ed_kinds,             only : wp, sp, ik
   use ed_time,              only : ed_time_t, seconds_between
   use ed_met_era5land,      only : era5land_path, era5land_default_template, era5land_month_hours, &
                                    ERA_NVAR, ERA_VAR_NAME, ERA_VAR_UNITS, ERA_CHUNK
   use ed_met_types,         only : met_grid_t
   use ed_netcdf_c
   use iso_c_binding,        only : c_int, c_size_t, c_double
   use ed_soil_hydr_lib,     only : soil_texture_t
   implicit none
   private

   public :: NLAT, NLON, LAT0, DLAT, LON0, DLON, T0
   public :: is_valid, field, air_temperature, hour_index, write_archive, write_static, write_month
   public :: archive_grid
   public :: write_map
   public :: write_surface_maps, soil_depth_at, bc_class_at, f_leak_at, deep_temp_at, chave_e_at, texture_at
   public :: SG_TOP, SG_BOTTOM

   integer(ik), parameter :: NLAT = 4_ik, NLON = 18_ik
   real(wp),    parameter :: LAT0 = 45.0_wp, DLAT = -10.0_wp, LON0 = -180.0_wp, DLON = 20.0_wp
   real(wp),    parameter :: TWO_PI = 6.283185307179586_wp
   !----- The archive's first stamp; g below counts hours from it (0-based). ---------------------!
   type(ed_time_t), parameter :: T0 = ed_time_t(2021_ik, 1_ik, 1_ik, 1_ik, 0_ik, 0_ik)
   !----- The soil maps' layers: SoilGrids 2.0's depths [m]. ---------------------------------------!
   real(wp), parameter :: SG_TOP(6)    = [0.0_wp, 0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp]
   real(wp), parameter :: SG_BOTTOM(6) = [0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp, 2.00_wp]

contains

   !----- The archive's grid, from one of its static files. -----------------------------------------!
   function archive_grid(static) result(grid)
      character(len=*), intent(in) :: static
      type(met_grid_t) :: grid
      grid = met_grid_t(static_file=static, chunk=ERA_CHUNK)
   end function archive_grid

   !----- No data at (row 0, col 5) and (row 2, col 0); everything else is valid. ----------------!
   pure logical function is_valid(r, c)
      integer(ik), intent(in) :: r, c
      is_valid = .not. ((r == 0_ik .and. c == 5_ik) .or. (r == 2_ik .and. c == 0_ik))
   end function is_valid

   !----- The value of variable v at grid (row r, col c) and hour g, as the archive stores it. ---!
   pure real(wp) function field(v, r, c, g) result(x)
      integer(ik), intent(in) :: v, r, c, g
      real(wp) :: hr, day
      hr = real(g, wp) ; day = real(modulo(g, 24_ik), wp)
      select case (v)
      case (1) ; x = air_temperature(r, c, hr)                                                   ! Tair
      case (2) ; x = air_temperature(r, c, hr) - 4.0_wp - 0.5_wp * sin(TWO_PI * hr / 24.0_wp)    ! Tdew
      case (3) ; x = 90000.0_wp + 100.0_wp * real(r, wp) + 10.0_wp * real(c, wp) + day         ! PSurf
      case (4) ; x = 2.0_wp + 3.0_wp * sin(TWO_PI * hr / 24.0_wp) + 0.1_wp * real(c, wp)       ! u10
      case (5) ; x = -1.0_wp + 2.0_wp * cos(TWO_PI * hr / 17.0_wp) - 0.1_wp * real(r, wp)      ! v10
      case (6) ; x = merge(1.0e-4_wp, 0.0_wp, modulo(g, 7_ik) == 3_ik)                          ! Rainf
      case (7) ; x = max(0.0_wp, 600.0_wp * sin(TWO_PI * (day - 6.0_wp) / 24.0_wp))              ! SWdown
      case default ; x = 300.0_wp + 0.001_wp * hr + real(c, wp)                                  ! LWdown
      end select
      x = real(real(x, sp), wp)                              ! the archive holds float32
   end function field

   pure real(wp) function air_temperature(r, c, hr) result(x)
      integer(ik), intent(in) :: r, c
      real(wp),    intent(in) :: hr
      x = 260.0_wp + 0.002_wp * hr + real(r, wp) + 0.1_wp * real(c, wp)
   end function air_temperature

   !----- Hours from T0 to t. -----------------------------------------------------------------------!
   integer(ik) function hour_index(t) result(g)
      type(ed_time_t), intent(in) :: t
      g = nint(seconds_between(T0, t) / 3600.0_wp, ik)
   end function hour_index

   !=======================================================================================!
   !  The synthetic archive.                                                                     !
   !=======================================================================================!
   subroutine write_archive(dir)
      character(len=*), intent(in) :: dir
      integer(ik) :: m, v
      call execute_command_line('mkdir -p '//dir)
      call write_static(dir//'/ED_ERA5land_static.nc', .false.)
      call write_static(dir//'/static_allvalid.nc', .true.)
      do m = 1_ik, 12_ik
         do v = 1_ik, ERA_NVAR
            call write_month(dir, v, 2021_ik, m)
         end do
      end do
      do v = 1_ik, ERA_NVAR
         call write_month(dir, v, 2022_ik, 1_ik)
      end do
   end subroutine write_archive

   subroutine define_grid(ncid, dlat_id, dlon_id)
      integer(c_int), intent(in)  :: ncid
      integer(c_int), intent(out) :: dlat_id, dlon_id
      integer(c_int) :: st, vla, vlo, dims1(1)
      st = nc_def_dim_f(ncid, 'lat', int(NLAT, c_size_t), dlat_id) ; call nc_check(st, 'def lat dim')
      st = nc_def_dim_f(ncid, 'lon', int(NLON, c_size_t), dlon_id) ; call nc_check(st, 'def lon dim')
      dims1 = dlat_id
      st = nc_def_var_f(ncid, 'lat', NC_DOUBLE, 1, dims1, vla) ; call nc_check(st, 'def lat')
      dims1 = dlon_id
      st = nc_def_var_f(ncid, 'lon', NC_DOUBLE, 1, dims1, vlo) ; call nc_check(st, 'def lon')
   end subroutine define_grid

   subroutine put_grid(ncid)
      integer(c_int), intent(in) :: ncid
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start1(1), count1(1)
      real(c_double)    :: lat(NLAT), lon(NLON)
      integer(ik) :: k
      lat = [(LAT0 + DLAT * real(k, wp), k = 0_ik, NLAT - 1_ik)]
      lon = [(LON0 + DLON * real(k, wp), k = 0_ik, NLON - 1_ik)]
      start1 = 0_c_size_t
      st = nc_inq_varid_f(ncid, 'lat', vid) ; count1 = int(NLAT, c_size_t)
      st = nc_put_vara_double(ncid, vid, start1, count1, lat) ; call nc_check(st, 'put lat')
      st = nc_inq_varid_f(ncid, 'lon', vid) ; count1 = int(NLON, c_size_t)
      st = nc_put_vara_double(ncid, vid, start1, count1, lon) ; call nc_check(st, 'put lon')
   end subroutine put_grid

   subroutine write_static(path, all_valid)
      character(len=*), intent(in) :: path
      logical,          intent(in) :: all_valid
      integer(c_int)    :: st, ncid, dla, dlo, vv, ve, vf, dims2(2)
      integer(c_size_t) :: start2(2), count2(2)
      real(c_double)    :: valid(NLON, NLAT), elev(NLON, NLAT), lfrac(NLON, NLAT)
      integer(ik) :: r, c
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create static')
      call define_grid(ncid, dla, dlo)
      dims2 = [dla, dlo]
      st = nc_def_var_f(ncid, 'valid', NC_BYTE, 2, dims2, vv) ; call nc_check(st, 'def valid')
      st = nc_def_var_f(ncid, 'elevation', NC_FLOAT, 2, dims2, ve) ; call nc_check(st, 'def elevation')
      st = nc_def_var_f(ncid, 'land_fraction', NC_FLOAT, 2, dims2, vf) ; call nc_check(st, 'def land_fraction')
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef static')
      call put_grid(ncid)
      do r = 0_ik, NLAT - 1_ik
         do c = 0_ik, NLON - 1_ik
            valid(c + 1_ik, r + 1_ik) = merge(1.0_c_double, 0.0_c_double, all_valid .or. is_valid(r, c))
            elev(c + 1_ik, r + 1_ik)  = 100.0_c_double * real(r, c_double) + real(c, c_double)
            lfrac(c + 1_ik, r + 1_ik) = merge(0.3_c_double, 1.0_c_double, c == 1_ik)   ! column 1: mostly water
         end do
      end do
      start2 = 0_c_size_t ; count2 = [int(NLAT, c_size_t), int(NLON, c_size_t)]
      st = nc_put_vara_double(ncid, vv, start2, count2, valid) ; call nc_check(st, 'put valid')
      st = nc_put_vara_double(ncid, ve, start2, count2, elev) ; call nc_check(st, 'put elevation')
      st = nc_put_vara_double(ncid, vf, start2, count2, lfrac) ; call nc_check(st, 'put land_fraction')
      st = nc_close(ncid) ; call nc_check(st, 'close static')
   end subroutine write_static

   !----- A static map on the archive's grid (ed_maps): values(lon, lat, layer) as a 2-D      !
   !      (lat, lon) variable, or, with layer_dim, a layered (layer, lat, lon) one, with a _FillValue.  !
   !      lat_shift moves the file's lat axis, to make a map that is off the grid. -----------------------!
   subroutine write_map(path, name, values, fill, layer_dim, lat_shift)
      character(len=*),           intent(in) :: path, name
      real(wp),                   intent(in) :: values(:,:,:)
      real(wp),                   intent(in) :: fill
      character(len=*), optional, intent(in) :: layer_dim
      real(wp),         optional, intent(in) :: lat_shift
      integer(c_int)    :: st, ncid, dla, dlo, dly, vid, vla, dims2(2), dims3(3)
      integer(c_size_t) :: start3(3), count3(3), start2(2), count2(2), start1(1), count1(1)
      real(c_double)    :: buf(NLON, NLAT, size(values, 3)), lat(NLAT)
      integer(ik)       :: k
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create map')
      call define_grid(ncid, dla, dlo)
      if (present(layer_dim)) then
         st = nc_def_dim_f(ncid, layer_dim, int(size(values, 3), c_size_t), dly) ; call nc_check(st, 'def layer')
         dims3 = [dly, dla, dlo]
         st = nc_def_var_f(ncid, name, NC_DOUBLE, 3, dims3, vid) ; call nc_check(st, 'def map')
      else
         dims2 = [dla, dlo]
         st = nc_def_var_f(ncid, name, NC_DOUBLE, 2, dims2, vid) ; call nc_check(st, 'def map')
      end if
      st = nc_put_att_double_f(ncid, vid, '_FillValue', NC_DOUBLE, real(fill, c_double))
      call nc_check(st, 'put _FillValue')
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef map')
      call put_grid(ncid)
      if (present(lat_shift)) then
         lat = [(LAT0 + DLAT * real(k, wp) + lat_shift, k = 0_ik, NLAT - 1_ik)]
         st = nc_inq_varid_f(ncid, 'lat', vla)
         start1 = 0_c_size_t ; count1 = int(NLAT, c_size_t)
         st = nc_put_vara_double(ncid, vla, start1, count1, lat) ; call nc_check(st, 'shift lat')
      end if
      buf = real(values, c_double)
      if (present(layer_dim)) then
         start3 = 0_c_size_t
         count3 = [int(size(values, 3), c_size_t), int(NLAT, c_size_t), int(NLON, c_size_t)]
         st = nc_put_vara_double(ncid, vid, start3, count3, buf)
      else
         start2 = 0_c_size_t ; count2 = [int(NLAT, c_size_t), int(NLON, c_size_t)]
         st = nc_put_vara_double(ncid, vid, start2, count2, buf)
      end if
      call nc_check(st, 'put map')
      st = nc_close(ncid) ; call nc_check(st, 'close map')
   end subroutine write_map

   !=======================================================================================!
   !  The [maps] files, each value a function of its cell (and layer):                         !
   !  dir/soil_depth.nc (soil_depth, bc_class, f_leak, deep_temp), dir/soil_texture.nc and       !
   !  dir/soil_chemistry.nc (on SoilGrids' layers, with depth_bnds), dir/chave_e.nc. The no-data  !
   !  cells are empty; gap = [row, col] also leaves sand empty at that cell, and f_leak sets     !
   !  every cell's leaky drainage factor to one value.                                           !
   !=======================================================================================!
   pure real(wp) function soil_depth_at(r, c) result(x)
      integer(ik), intent(in) :: r, c
      x = 0.8_wp + 0.3_wp * real(r, wp) + 0.05_wp * real(c, wp)
   end function soil_depth_at

   pure integer(ik) function bc_class_at(r, c) result(k)
      integer(ik), intent(in) :: r, c
      k = 1_ik + modulo(r + c, 4_ik)
   end function bc_class_at

   !----- A leaky base's drainage factor, in (0, 1] and different between cells (the ED soil maps   !
   !      hold 0.05 everywhere). --------------------------------------------------------------------!
   pure real(wp) function f_leak_at(r, c) result(x)
      integer(ik), intent(in) :: r, c
      x = 0.02_wp + 0.01_wp * real(modulo(r + 2_ik * c, 5_ik), wp)
   end function f_leak_at

   pure real(wp) function deep_temp_at(r, c) result(x)
      integer(ik), intent(in) :: r, c
      x = 280.0_wp + real(r, wp) + 0.1_wp * real(c, wp)
   end function deep_temp_at

   pure real(wp) function chave_e_at(r, c) result(x)
      integer(ik), intent(in) :: r, c
      x = 0.05_wp + 0.01_wp * real(c, wp) - 0.02_wp * real(r, wp)
   end function chave_e_at

   pure function texture_at(l, r, c) result(t)
      integer(ik), intent(in) :: l, r, c
      type(soil_texture_t) :: t
      t%sand = 25.0_wp + 2.0_wp * real(l, wp) + 0.5_wp * real(c, wp)
      t%clay = 45.0_wp - 3.0_wp * real(l, wp) + real(r, wp)
      t%silt = 100.0_wp - t%sand - t%clay
      t%bulk_density = 1.15_wp + 0.03_wp * real(l, wp)
      t%coarse_frag  = 2.0_wp + real(l, wp)
      t%org_carbon   = 30.0_wp / real(l, wp)
      t%cec          = 15.0_wp + real(r, wp)
      t%ph           = 5.5_wp + 0.1_wp * real(l, wp)
   end function texture_at

   subroutine write_surface_maps(dir, gap, f_leak)
      character(len=*),      intent(in) :: dir
      integer(ik), optional, intent(in) :: gap(2)
      real(wp),    optional, intent(in) :: f_leak
      real(wp), parameter :: FILL = -9999.0_wp
      real(wp) :: v2(NLON, NLAT, 1, 4), v3(NLON, NLAT, 6, 5), vc(NLON, NLAT, 6, 3), ve(NLON, NLAT, 1, 1)
      type(soil_texture_t) :: t
      integer(ik) :: r, c, l
      v2 = FILL ; v3 = FILL ; vc = FILL ; ve = FILL
      do r = 0_ik, NLAT - 1_ik
         do c = 0_ik, NLON - 1_ik
            if (.not. is_valid(r, c)) cycle
            v2(c + 1_ik, r + 1_ik, 1, :) = [soil_depth_at(r, c), real(bc_class_at(r, c), wp), f_leak_at(r, c), &
                                            deep_temp_at(r, c)]
            ve(c + 1_ik, r + 1_ik, 1, 1) = chave_e_at(r, c)
            do l = 1_ik, 6_ik
               t = texture_at(l, r, c)
               v3(c + 1_ik, r + 1_ik, l, :) = [t%sand, t%silt, t%clay, t%bulk_density, t%coarse_frag]
               vc(c + 1_ik, r + 1_ik, l, :) = [t%org_carbon, t%cec, t%ph]
            end do
         end do
      end do
      if (present(gap)) v3(gap(2) + 1_ik, gap(1) + 1_ik, :, 1) = FILL
      if (present(f_leak)) where (v2(:, :, 1, 3) /= FILL) v2(:, :, 1, 3) = f_leak
      call execute_command_line('mkdir -p '//dir)
      call write_fields(dir//'/soil_depth.nc', [character(len=16) :: 'soil_depth', 'bc_class', 'f_leak',     &
                        'deep_temp'], v2, FILL, .false.)
      call write_fields(dir//'/soil_texture.nc', [character(len=16) :: 'sand', 'silt', 'clay', 'bulk_density', &
                        'coarse_fragments'], v3, FILL, .true.)
      call write_fields(dir//'/soil_chemistry.nc', [character(len=16) :: 'organic_carbon', 'cec', 'ph_h2o'], &
                        vc, FILL, .true.)
      call write_fields(dir//'/chave_e.nc', [character(len=16) :: 'chave_e'], ve, FILL, .false.)
   end subroutine write_surface_maps

   !----- Several map variables in one file: values(lon, lat, layer, variable), 2-D (lat, lon) or, !
   !      layered, (depth, lat, lon) with SoilGrids' depth_bnds. ------------------------------------!
   subroutine write_fields(path, names, values, fill, layered)
      character(len=*), intent(in) :: path, names(:)
      real(wp),         intent(in) :: values(:,:,:,:)
      real(wp),         intent(in) :: fill
      logical,          intent(in) :: layered
      integer(c_int)    :: st, ncid, dla, dlo, dly, dnv, vb, dims2(2), dims3(3)
      integer(c_int)    :: vid(size(names))
      integer(c_size_t) :: start2(2), count2(2), start3(3), count3(3)
      real(c_double)    :: bnds(2, 6), buf(NLON, NLAT, size(values, 3))
      integer(ik)       :: f
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create '//path)
      call define_grid(ncid, dla, dlo)
      if (layered) then
         st = nc_def_dim_f(ncid, 'depth', int(size(values, 3), c_size_t), dly) ; call nc_check(st, 'def depth')
         st = nc_def_dim_f(ncid, 'nv', 2_c_size_t, dnv) ; call nc_check(st, 'def nv')
         dims2 = [dly, dnv]
         st = nc_def_var_f(ncid, 'depth_bnds', NC_DOUBLE, 2, dims2, vb) ; call nc_check(st, 'def depth_bnds')
      end if
      do f = 1_ik, size(names, kind=ik)
         if (layered) then
            dims3 = [dly, dla, dlo]
            st = nc_def_var_f(ncid, trim(names(f)), NC_DOUBLE, 3, dims3, vid(f))
         else
            dims2 = [dla, dlo]
            st = nc_def_var_f(ncid, trim(names(f)), NC_DOUBLE, 2, dims2, vid(f))
         end if
         call nc_check(st, 'def '//trim(names(f)))
         st = nc_put_att_double_f(ncid, vid(f), '_FillValue', NC_DOUBLE, real(fill, c_double))
         call nc_check(st, 'put _FillValue')
      end do
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef '//path)
      call put_grid(ncid)
      if (layered) then
         bnds(1, :) = SG_TOP ; bnds(2, :) = SG_BOTTOM
         start2 = 0_c_size_t ; count2 = [6_c_size_t, 2_c_size_t]
         st = nc_put_vara_double(ncid, vb, start2, count2, bnds) ; call nc_check(st, 'put depth_bnds')
      end if
      do f = 1_ik, size(names, kind=ik)
         buf = real(values(:, :, :, f), c_double)
         if (layered) then
            start3 = 0_c_size_t
            count3 = [int(size(values, 3), c_size_t), int(NLAT, c_size_t), int(NLON, c_size_t)]
            st = nc_put_vara_double(ncid, vid(f), start3, count3, buf)
         else
            start2 = 0_c_size_t ; count2 = [int(NLAT, c_size_t), int(NLON, c_size_t)]
            st = nc_put_vara_double(ncid, vid(f), start2, count2, buf)
         end if
         call nc_check(st, 'put '//trim(names(f)))
      end do
      st = nc_close(ncid) ; call nc_check(st, 'close '//path)
   end subroutine write_fields

   subroutine write_month(dir, v, y, m, fill)
      character(len=*), intent(in) :: dir
      integer(ik),      intent(in) :: v, y, m
      real(c_double), optional, intent(in) :: fill   !< a numeric _FillValue to declare (a bad file)
      character(len=:), allocatable :: path
      integer(c_int)    :: st, ncid, dt, dla, dlo, vt, vx, dims1(1), dims3(3)
      integer(c_size_t) :: start1(1), count1(1), start3(3), count3(3)
      real(c_double), allocatable :: tsec(:), dat(:,:,:)
      real(c_double) :: nan
      integer(ik) :: nt, g0, h, r, c
      nan = ieee_value(0.0_c_double, ieee_quiet_nan)
      nt = era5land_month_hours(y, m)
      g0 = hour_index(ed_time_t(y, m, 1_ik, 1_ik))
      path = trim(era5land_path(era5land_default_template(), dir, ERA_VAR_NAME(v), y, m))
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create '//path)
      st = nc_def_dim_f(ncid, 'time', int(nt, c_size_t), dt) ; call nc_check(st, 'def time dim')
      call define_grid(ncid, dla, dlo)
      dims1 = dt
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt) ; call nc_check(st, 'def time')
      st = nc_put_att_text_f(ncid, vt, 'units', 33_c_size_t, 'seconds since 1970-01-01 00:00:00')
      call nc_check(st, 'time units')
      dims3 = [dt, dla, dlo]
      st = nc_def_var_f(ncid, trim(ERA_VAR_NAME(v)), NC_FLOAT, 3, dims3, vx) ; call nc_check(st, 'def var')
      st = nc_put_att_text_f(ncid, vx, 'units', int(len_trim(ERA_VAR_UNITS(v)), c_size_t), trim(ERA_VAR_UNITS(v)))
      call nc_check(st, 'var units')
      if (present(fill)) then
         st = nc_put_att_double_f(ncid, vx, '_FillValue', NC_FLOAT, fill) ; call nc_check(st, 'var fill')
      end if
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'avg_convention', 3_c_size_t, 'end') ; call nc_check(st, 'avg_convention')
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef month')
      call put_grid(ncid)
      allocate(tsec(nt), dat(NLON, NLAT, nt))
      do h = 1_ik, nt
         tsec(h) = seconds_between(ed_time_t(1970_ik, 1_ik, 1_ik), T0) + 3600.0_wp * real(g0 + h - 1_ik, wp)
         do r = 0_ik, NLAT - 1_ik
            do c = 0_ik, NLON - 1_ik
               dat(c + 1_ik, r + 1_ik, h) = merge(field(v, r, c, g0 + h - 1_ik), nan, is_valid(r, c))
            end do
         end do
      end do
      start1 = 0_c_size_t ; count1 = int(nt, c_size_t)
      st = nc_put_vara_double(ncid, vt, start1, count1, tsec) ; call nc_check(st, 'put time')
      start3 = 0_c_size_t ; count3 = [int(nt, c_size_t), int(NLAT, c_size_t), int(NLON, c_size_t)]
      st = nc_put_vara_double(ncid, vx, start3, count3, dat)
      call nc_check(st, 'put '//path)
      st = nc_close(ncid) ; call nc_check(st, 'close '//path)
   end subroutine write_month

end module ed_test_archive
