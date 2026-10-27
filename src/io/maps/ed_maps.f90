! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_maps -- every static per-cell map a run reads; no other module opens a map
! (docs/configuration.md, Input maps). A map is on the forcing grid, so a cell's value
! sits at its own (row, col) and ED-MORE never regrids; a value the map does not hold
! reads as NaN. map_values_at_cells and map_layers_at_cells read one variable at a list of
! cells; apply_cell_mask keeps a region's cells; read_map_inputs reads every [maps] map at
! a run's cells, which ed_config's apply_map_inputs turns into each polygon's soil column,
! bottom boundary and Chave E (docs/science/soil_biophysics.md §The column from the soil
! maps).
!---~---
module ed_maps
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use iso_c_binding,         only : c_int, c_size_t, c_double
   use ed_kinds,              only : wp, ik
   use ed_netcdf_c,           only : nc_close, nc_get_vara_double, nc_inq_dimlen_f, nc_inq_varid_f, nc_open_f,   &
                                     nc_packing_t, nc_packing_of, nc_unpacked, NC_NOERR, NC_NOWRITE
   use ed_met_types,          only : met_grid_t, met_cells_t
   use ed_met_grid,           only : same_grid, keep_cells, GRID_OK, GRID_ERR_OPEN
   use ed_maps_params,        only : maps_params_t, soil_maps_on, map_inputs_t, MAP_BC_FREE_DRAIN, MAP_BC_LEAKY, &
                                     MAP_BC_SEALED, MAP_BC_AQUIFER
   implicit none
   private

   public :: map_values_at_cells, map_layers_at_cells, apply_cell_mask, map_status_message
   public :: MAP_OK, MAP_ERR_OPEN, MAP_ERR_GRID, MAP_ERR_VAR
   public :: read_map_inputs

   integer(ik), parameter :: MAP_OK       = 0_ik
   integer(ik), parameter :: MAP_ERR_OPEN = 1_ik   !--- the file (or the grid's static file) cannot be opened
   integer(ik), parameter :: MAP_ERR_GRID = 2_ik   !--- the file is not on the forcing grid
   integer(ik), parameter :: MAP_ERR_VAR  = 3_ik   !--- the variable (or its layer dimension) is not there

   !--- The soil maps' layers (ED-SoilMaps v2: SoilGrids 2.0's six depths, to 2 m).
   character(len=*), parameter :: SOIL_LAYER_DIM = 'depth', SOIL_LAYER_BOUNDS = 'depth_bnds'

contains

   !--- A 2-D (lat, lon) map variable at every cell: values(cell), NaN where the map is empty.
   subroutine map_values_at_cells(path, var, grid, cells, values, stat)
      character(len=*),      intent(in)  :: path, var
      type(met_grid_t),      intent(in)  :: grid
      type(met_cells_t),     intent(in)  :: cells
      real(wp), allocatable, intent(out) :: values(:)
      integer(ik),           intent(out) :: stat
      real(wp), allocatable :: v(:,:)
      call read_cells(path, var, '', grid, cells, v, stat)
      if (stat == MAP_OK) values = v(1, :)
   end subroutine map_values_at_cells

   !---~---
   !   A layered (layer, lat, lon) map variable at every cell: values(layer, cell), NaN where the
   ! map is empty. `layer_dim` names the variable's first dimension.
   !---~---
   subroutine map_layers_at_cells(path, var, layer_dim, grid, cells, values, stat)
      character(len=*),      intent(in)  :: path, var, layer_dim
      type(met_grid_t),      intent(in)  :: grid
      type(met_cells_t),     intent(in)  :: cells
      real(wp), allocatable, intent(out) :: values(:,:)
      integer(ik),           intent(out) :: stat
      call read_cells(path, var, layer_dim, grid, cells, values, stat)
   end subroutine map_layers_at_cells

   !---~---
   !   The cell mask of a region (docs/configuration.md, Regional runs): keep the cells whose
   ! mask value is at least mask_min, not those it leaves empty; n_dropped counts the rest.
   !---~---
   subroutine apply_cell_mask(cells, mask_path, mask_var, mask_min, grid, n_dropped, stat)
      type(met_cells_t), intent(inout) :: cells
      character(len=*),  intent(in)    :: mask_path, mask_var
      type(met_grid_t),  intent(in)    :: grid
      real(wp),          intent(in)    :: mask_min
      integer(ik),       intent(out)   :: n_dropped, stat
      real(wp), allocatable :: v(:)
      logical,  allocatable :: keep(:)
      n_dropped = 0_ik
      call map_values_at_cells(mask_path, mask_var, grid, cells, v, stat)
      if (stat /= MAP_OK) return
      keep = .not. ieee_is_nan(v)
      where (keep) keep = v >= mask_min
      n_dropped = cells%ncell - int(count(keep), ik)
      call keep_cells(cells, keep)
   end subroutine apply_cell_mask

   !--- What a status code means, for the caller's message.
   function map_status_message(stat) result(msg)
      integer(ik), intent(in)       :: stat
      character(len=:), allocatable :: msg
      select case (stat)
      case (MAP_OK)
         msg = 'read'
      case (MAP_ERR_OPEN)
         msg = 'cannot be opened'
      case (MAP_ERR_GRID)
         msg = 'is not on the forcing grid (its lat and lon must be the static file''s)'
      case (MAP_ERR_VAR)
         msg = 'does not hold the variable (or its layer dimension)'
      case default
         msg = 'cannot be read'
      end select
   end function map_status_message

   !---~---
   !   read_map_inputs -- every [maps] map at every cell of a run: the soil depth map's four
   ! variables, the texture and chemistry maps' layers, and the Chave map's E. msg says
   ! what fails, naming the map, the variable and any cell: a failed read, a gap, a
   ! bc_class not 1-4, or a leaky base's f_leak outside (0, 1]; it is '' when every value
   ! is there.
   !---~---
   subroutine read_map_inputs(opts, grid, cells, inp, msg)
      type(maps_params_t),           intent(in)  :: opts
      type(met_grid_t),              intent(in)  :: grid
      type(met_cells_t),             intent(in)  :: cells
      type(map_inputs_t),            intent(out) :: inp
      character(len=:), allocatable, intent(out) :: msg
      real(wp), allocatable :: v(:,:), top_chem(:), bottom_chem(:)
      integer(ik) :: k
      msg = ''
      if (soil_maps_on(opts)) then
         inp%soil = .true.
         call need(opts%soil_depth_map, 'soil_depth', '', v)
         if (len(msg) > 0) return
         inp%soil_depth = v(1, :)
         call need(opts%soil_depth_map, 'bc_class', '', v)
         if (len(msg) > 0) return
         inp%bc_class = v(1, :)
         do k = 1_ik, cells%ncell
            if (all(nint(inp%bc_class(k), ik) /= [MAP_BC_FREE_DRAIN, MAP_BC_LEAKY, MAP_BC_SEALED, MAP_BC_AQUIFER])) then
               msg = 'soil map '//trim(opts%soil_depth_map)//': bc_class is not 1-4 at '//cell_name(k)
               return
            end if
         end do
         call need(opts%soil_depth_map, 'f_leak', '', v)
         if (len(msg) > 0) return
         inp%f_leak = v(1, :)
         do k = 1_ik, cells%ncell
            if (nint(inp%bc_class(k), ik) == MAP_BC_LEAKY .and. .not. (inp%f_leak(k) > 0.0_wp .and.           &
                                                                  inp%f_leak(k) <= 1.0_wp)) then
               msg = 'soil map '//trim(opts%soil_depth_map)//': f_leak is not in (0, 1] at '//cell_name(k)// &
                     ', whose base is leaky'
               return
            end if
         end do
         call need(opts%soil_depth_map, 'deep_temp', '', v)
         if (len(msg) > 0) return
         inp%deep_temp = v(1, :)

         !--- The layers: the texture and the chemistry maps must share them.
         call layer_bounds(opts%soil_texture_map, inp%layer_top, inp%layer_bottom)
         if (len(msg) > 0) return
         call layer_bounds(opts%soil_chemistry_map, top_chem, bottom_chem)
         if (len(msg) > 0) return
         if (size(top_chem) /= size(inp%layer_top)) then
            msg = 'soil maps '//trim(opts%soil_texture_map)//' and '//trim(opts%soil_chemistry_map)// &
                  ' do not have the same layers'
            return
         end if
         if (any(top_chem /= inp%layer_top) .or. any(bottom_chem /= inp%layer_bottom)) then
            msg = 'soil maps '//trim(opts%soil_texture_map)//' and '//trim(opts%soil_chemistry_map)// &
                  ' do not have the same layers'
            return
         end if
         allocate(inp%texture(size(inp%layer_top), cells%ncell))
         call need(opts%soil_texture_map, 'sand', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%sand = v
         call need(opts%soil_texture_map, 'silt', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%silt = v
         call need(opts%soil_texture_map, 'clay', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%clay = v
         call need(opts%soil_texture_map, 'bulk_density', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%bulk_density = v
         call need(opts%soil_texture_map, 'coarse_fragments', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%coarse_frag = v
         call need(opts%soil_chemistry_map, 'organic_carbon', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%org_carbon = v
         call need(opts%soil_chemistry_map, 'cec', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%cec = v
         call need(opts%soil_chemistry_map, 'ph_h2o', SOIL_LAYER_DIM, v)
         if (len(msg) > 0) return
         inp%texture%ph = v
      end if
      if (len_trim(opts%chave_e_map) > 0) then
         inp%chave = .true.
         call need(opts%chave_e_map, 'chave_e', '', v)
         if (len(msg) > 0) return
         inp%chave_e = v(1, :)
      end if

   contains

      !--- One variable at every cell, and every layer of it; msg on a failed read or a gap.
      subroutine need(path, var, layer_dim, values)
         character(len=*),      intent(in)  :: path, var, layer_dim
         real(wp), allocatable, intent(out) :: values(:,:)
         integer(ik) :: st, c
         call read_cells(trim(path), var, layer_dim, grid, cells, values, st)
         if (st /= MAP_OK) then
            msg = 'map '//trim(path)//' '//map_status_message(st)//': '//var
            return
         end if
         do c = 1_ik, cells%ncell
            if (any(ieee_is_nan(values(:, c)))) then
               msg = 'map '//trim(path)//': '//var//' has no value at '//cell_name(c)
               return
            end if
         end do
      end subroutine need

      !--- A soil map's layers: the tops and bottoms [m] of its depth_bnds(depth, nv).
      subroutine layer_bounds(path, top, bottom)
         character(len=*),      intent(in)  :: path
         real(wp), allocatable, intent(out) :: top(:), bottom(:)
         integer(c_int)    :: ncid, vid, st, st_close
         integer(c_size_t) :: n, start2(2), count2(2)
         real(c_double), allocatable :: b(:,:)
         msg = 'map '//trim(path)//' '//map_status_message(MAP_ERR_OPEN)
         if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) return
         msg = 'map '//trim(path)//' '//map_status_message(MAP_ERR_VAR)//': '//SOIL_LAYER_BOUNDS
         st = nc_inq_dimlen_f(ncid, SOIL_LAYER_DIM, n)
         if (st == NC_NOERR) st = nc_inq_varid_f(ncid, SOIL_LAYER_BOUNDS, vid)
         if (st == NC_NOERR) then
            allocate(b(2, n))
            start2 = 0_c_size_t
            count2 = [n, 2_c_size_t]
            st = nc_get_vara_double(ncid, vid, start2, count2, b)
         end if
         st_close = nc_close(ncid)
         if (st /= NC_NOERR .or. st_close /= NC_NOERR) return
         top = real(b(1, :), wp)
         bottom = real(b(2, :), wp)
         msg = ''
      end subroutine layer_bounds

      !--- A cell, as a message names it.
      function cell_name(c) result(name)
         integer(ik), intent(in) :: c
         character(len=:), allocatable :: name
         character(len=96) :: buf
         write(buf, '(a,f0.2,a,f0.2,a,i0,a,i0,a)') 'the cell (', cells%lat(c), ', ', cells%lon(c),      &
                                                  ') (row ', cells%row(c), ', col ', cells%col(c), ')'
         name = trim(buf)
      end function cell_name

   end subroutine read_map_inputs

   !---~---
   !   The read: check the grid, read the band of rows the cells span (every column and
   ! layer, in one read per variable, which also serves a box across 180 degrees), then
   ! pick each cell's values. _FillValue, missing_value and NaN read as NaN, and packed
   ! values are unpacked. layer_dim = '' reads a 2-D variable as one layer.
   !---~---
   subroutine read_cells(path, var, layer_dim, grid, cells, values, stat)
      character(len=*),      intent(in)  :: path, var, layer_dim
      type(met_grid_t),      intent(in)  :: grid
      type(met_cells_t),     intent(in)  :: cells
      real(wp), allocatable, intent(out) :: values(:,:)
      integer(ik),           intent(out) :: stat
      integer(c_int)    :: ncid, vid, st
      integer(c_size_t) :: nlayer_c, start3(3), count3(3), start2(2), count2(2)
      integer(ik)       :: nlayer, r0, nr, k, l
      real(c_double), allocatable :: band(:,:,:)
      type(nc_packing_t) :: packing

      allocate(values(1, cells%ncell))
      values = ieee_value(1.0_wp, ieee_quiet_nan)
      call same_grid(path, grid, stat)
      if (stat /= GRID_OK) then
         stat = merge(MAP_ERR_OPEN, MAP_ERR_GRID, stat == GRID_ERR_OPEN)
         return
      end if
      stat = MAP_ERR_OPEN
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) return
      stat = MAP_ERR_VAR
      st = nc_inq_varid_f(ncid, trim(var), vid)
      if (st /= NC_NOERR) then
         st = nc_close(ncid)
         return
      end if
      nlayer = 1_ik
      if (len_trim(layer_dim) > 0) then
         st = nc_inq_dimlen_f(ncid, trim(layer_dim), nlayer_c)
         if (st /= NC_NOERR) then
            st = nc_close(ncid)
            return
         end if
         nlayer = int(nlayer_c, ik)
      end if
      if (cells%ncell == 0_ik) then
         st = nc_close(ncid)
         deallocate(values)
         allocate(values(nlayer, 0))
         stat = MAP_OK
         return
      end if

      !--- The band of rows the cells span, every column: one hyperslab, x(lon, lat, layer).
      r0 = minval(cells%row)
      nr = maxval(cells%row) - r0 + 1_ik
      allocate(band(cells%nlon, nr, nlayer))
      if (len_trim(layer_dim) > 0) then
         start3 = [0_c_size_t, int(r0, c_size_t), 0_c_size_t]
         count3 = [int(nlayer, c_size_t), int(nr, c_size_t), int(cells%nlon, c_size_t)]
         st = nc_get_vara_double(ncid, vid, start3, count3, band)
      else
         start2 = [int(r0, c_size_t), 0_c_size_t]
         count2 = [int(nr, c_size_t), int(cells%nlon, c_size_t)]
         st = nc_get_vara_double(ncid, vid, start2, count2, band)
      end if
      if (st /= NC_NOERR) then
         st = nc_close(ncid)
         return
      end if
      packing = nc_packing_of(ncid, vid)
      st = nc_close(ncid)

      deallocate(values)
      allocate(values(nlayer, cells%ncell))
      !--- Empty where the map says so; what is there, unpacked.
      do k = 1_ik, cells%ncell
         do l = 1_ik, nlayer
            values(l, k) = nc_unpacked(packing, band(cells%col(k) + 1_ik, cells%row(k) - r0 + 1_ik, l))
         end do
      end do
      stat = MAP_OK
   end subroutine read_cells

end module ed_maps
