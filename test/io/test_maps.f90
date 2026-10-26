! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_maps -- meds_maps, the one module that reads static per-cell maps (docs/configuration.md §Input !
! maps), on the synthetic archive's 4 x 18 grid.                                                      !
!                                                                                          !
! A 2-D map is read at a box's cells, each at its own (row, col); a cell the map leaves empty reads  !
! as NaN. A layered map gives every layer at every cell. A map off the forcing grid, a missing       !
! variable and a missing file are refused, each with its own status. The cell mask keeps the cells   !
! at or above its threshold, drops the empty one, and leaves the kept cells grouped by chunk the way !
! the forcing reader needs them.                                                                    !
!                                                                                          !
! The [maps] inputs: a run cell the soil maps leave empty is named, in a 2-D map and in a layered    !
! one; every value is read at its own cell; a column layer takes the map layers it overlaps,         !
! weighted by the overlap, and below the deepest map layer that layer's soil; the bottom boundary    !
! and the thermal anchor follow the map (the anchor's distance is docs/science/soil_biophysics.md    !
! eq. (13)): a leaky base drains at its cell's f_leak, a free one at 1 and a sealed one at 0, and a  !
! leaky f_leak outside (0, 1] is refused.                                                            !
! The work directory is the first argument.                                                          !
!==========================================================================================!
program test_maps
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use meds_kinds,                 only : wp, ik
   use meds_test_assert,           only : check_true, test_report
   use meds_test_archive,          only : NLAT, NLON, write_static, write_map, write_surface_maps, archive_grid, &
                                          soil_depth_at, bc_class_at, f_leak_at, deep_temp_at, chave_e_at, &
                                          texture_at, SG_TOP, SG_BOTTOM
   use meds_met_types,             only : met_grid_t, met_cells_t
   use meds_met_grid,              only : select_box_cells, keep_cells, GRID_OK
   use meds_maps,                  only : map_values_at_cells, map_layers_at_cells, apply_cell_mask,  &
                                          MAP_OK, MAP_ERR_OPEN, MAP_ERR_GRID, MAP_ERR_VAR,              &
                                          read_map_inputs
   use meds_config,                only : meds_config_t, alloc_pft_table, apply_map_inputs
   use meds_maps_params,           only : maps_params_t, map_inputs_t
   use meds_soil_hydr_lib,         only : soil_texture_t, soil_texture_to_hydraulics
   use meds_therm_lib,             only : soil_texture_to_thermal
   use meds_soil_params,           only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER, ANNUAL_DAMPING_DEPTH
   implicit none

   real(wp), parameter :: FILL = -9999.0_wp
   integer(ik), parameter :: NL = 3_ik
   character(len=512) :: work
   character(len=:), allocatable :: static
   type(met_grid_t) :: grid
   type(met_cells_t)  :: cells, masked
   real(wp) :: grid2(NLON, NLAT, 1), grid3(NLON, NLAT, NL)
   real(wp), allocatable :: v(:), v3(:,:)
   integer(ik) :: st, r, c, l, k, n_dropped
   logical :: ok
   type(maps_params_t)    :: opts
   type(map_inputs_t)     :: inp
   type(meds_config_t)    :: cfg, cfg2
   type(soil_texture_t)   :: t5, t6, tex
   character(len=:), allocatable :: msg
   real(wp) :: ts, tr, al, vn, ks, kso, kd, cv

   call get_command_argument(1, work)
   if (len_trim(work) == 0) error stop 'test_maps: pass the work directory'
   call execute_command_line('rm -rf "'//trim(work)//'" && mkdir -p "'//trim(work)//'"')
   static = trim(work)//'/static.nc'
   call write_static(static, .true.)
   grid = archive_grid(static)

   !----- The box: rows 0-3 (45 .. 15 N) and columns 0-3 (180 .. 120 W), every cell valid. ---------!
   call select_box_cells(grid, [45.0_wp, -180.0_wp, 15.0_wp, -120.0_wp], cells, st)
   call check_true('the box selects its 16 cells', st == GRID_OK .and. cells%ncell == 16_ik)

   !----- A 2-D map: 100 row + col, with (row 1, col 2) left empty. ----------------------------------!
   do r = 0_ik, NLAT - 1_ik
      do c = 0_ik, NLON - 1_ik
         grid2(c + 1_ik, r + 1_ik, 1) = 100.0_wp * real(r, wp) + real(c, wp)
      end do
   end do
   grid2(3, 2, 1) = FILL
   call write_map(trim(work)//'/map2.nc', 'value', grid2, FILL)
   call map_values_at_cells(trim(work)//'/map2.nc', 'value', grid, cells, v, st)
   call check_true('a 2-D map is read', st == MAP_OK .and. size(v) == cells%ncell)
   ok = .true.
   do k = 1_ik, cells%ncell
      if (cells%row(k) == 1_ik .and. cells%col(k) == 2_ik) then
         ok = ok .and. ieee_is_nan(v(k))
      else
         ok = ok .and. v(k) == 100.0_wp * real(cells%row(k), wp) + real(cells%col(k), wp)
      end if
   end do
   call check_true('each cell reads its own (row, col); the empty cell reads NaN', ok)

   !----- A layered map: 1000 layer + 100 row + col. -----------------------------------------------!
   do l = 1_ik, NL
      grid3(:, :, l) = grid2(:, :, 1) + 1000.0_wp * real(l, wp)
   end do
   grid3(3, 2, :) = FILL
   call write_map(trim(work)//'/map3.nc', 'value', grid3, FILL, layer_dim='depth')
   call map_layers_at_cells(trim(work)//'/map3.nc', 'value', 'depth', grid, cells, v3, st)
   call check_true('a layered map is read', st == MAP_OK .and. size(v3, 1) == NL .and. size(v3, 2) == cells%ncell)
   ok = .true.
   do k = 1_ik, cells%ncell
      do l = 1_ik, NL
         if (cells%row(k) == 1_ik .and. cells%col(k) == 2_ik) then
            ok = ok .and. ieee_is_nan(v3(l, k))
         else
            ok = ok .and. v3(l, k) == 1000.0_wp * real(l, wp) + 100.0_wp * real(cells%row(k), wp)       &
                                    + real(cells%col(k), wp)
         end if
      end do
   end do
   call check_true('every layer of every cell, NaN where empty', ok)

   !----- Refusals. ---------------------------------------------------------------------------------!
   call write_map(trim(work)//'/offgrid.nc', 'value', grid2, FILL, lat_shift=0.05_wp)
   call map_values_at_cells(trim(work)//'/offgrid.nc', 'value', grid, cells, v, st)
   call check_true('a map off the forcing grid is refused', st == MAP_ERR_GRID)
   call map_values_at_cells(trim(work)//'/map2.nc', 'other', grid, cells, v, st)
   call check_true('a missing variable is refused', st == MAP_ERR_VAR)
   call map_layers_at_cells(trim(work)//'/map3.nc', 'value', 'layer', grid, cells, v3, st)
   call check_true('a missing layer dimension is refused', st == MAP_ERR_VAR)
   call map_values_at_cells(trim(work)//'/none.nc', 'value', grid, cells, v, st)
   call check_true('a missing file is refused', st == MAP_ERR_OPEN)

   !----- The cell mask: keep the even columns at a 0.5 threshold; the empty cell is dropped. --------!
   do c = 0_ik, NLON - 1_ik
      grid2(c + 1_ik, :, 1) = merge(1.0_wp, 0.0_wp, modulo(c, 2_ik) == 0_ik)
   end do
   grid2(3, 2, 1) = FILL                                    ! (row 1, col 2): an even column, but empty
   call write_map(trim(work)//'/mask.nc', 'keep', grid2, FILL)
   masked = cells
   call apply_cell_mask(masked, trim(work)//'/mask.nc', 'keep', 0.5_wp, grid, n_dropped, st)
   call check_true('the mask is read', st == MAP_OK)
   call check_true('it keeps the 7 even-column cells that are not empty', masked%ncell == 7_ik            &
                   .and. n_dropped == 9_ik)
   call check_true('every kept cell is in an even column, none is the empty one',                         &
                   all(modulo(masked%col, 2_ik) == 0_ik) .and.                                           &
                   .not. any(masked%row == 1_ik .and. masked%col == 2_ik))
   call check_true('the kept cells keep their order', all(masked%row(2:) >= masked%row(:masked%ncell - 1_ik)))
   call check_true('the kept cells are regrouped by chunk', masked%chunk_first(masked%nchunk + 1_ik) - 1_ik  &
                   == masked%ncell .and. size(masked%by_chunk) == masked%ncell)
   !----- A fraction and a threshold: 0.1 col, kept at >= 0.25, so columns 3 and up. -----------------!
   do c = 0_ik, NLON - 1_ik
      grid2(c + 1_ik, :, 1) = 0.1_wp * real(c, wp)
   end do
   call write_map(trim(work)//'/frac.nc', 'natural', grid2, FILL)
   masked = cells
   call apply_cell_mask(masked, trim(work)//'/frac.nc', 'natural', 0.25_wp, grid, n_dropped, st)
   call check_true('a fraction map with a threshold keeps column 3 only', st == MAP_OK .and.               &
                   masked%ncell == 4_ik .and. all(masked%col == 3_ik))

   !----- The [maps] inputs. The box holds (row 2, col 0), which the maps leave empty: named. ------!
   call write_surface_maps(trim(work)//'/surface')
   opts%soil_depth_map     = trim(work)//'/surface/soil_depth.nc'
   opts%soil_texture_map   = trim(work)//'/surface/soil_texture.nc'
   opts%soil_chemistry_map = trim(work)//'/surface/soil_chemistry.nc'
   opts%chave_e_map        = trim(work)//'/surface/chave_e.nc'
   call read_map_inputs(opts, grid, cells, inp, msg)
   call check_true('a run cell the maps leave empty is named, with the variable',                        &
                   index(msg, 'soil_depth has no value at the cell (25.00, -180.00) (row 2, col 0)') > 0)
   masked = cells
   call keep_cells(masked, .not. (cells%row == 2_ik .and. cells%col == 0_ik))
   call read_map_inputs(opts, grid, masked, inp, msg)
   call check_true('the maps are read at the other 15 cells', len(msg) == 0 .and. inp%soil .and. inp%chave)
   call check_true('the soil maps'' layers are SoilGrids''', all(inp%layer_top == SG_TOP) .and.           &
                   all(inp%layer_bottom == SG_BOTTOM))
   ok = .true.
   do k = 1_ik, masked%ncell
      r = masked%row(k) ; c = masked%col(k)
      t5 = texture_at(3_ik, r, c) ; t6 = texture_at(6_ik, r, c)
      ok = ok .and. close(inp%soil_depth(k), soil_depth_at(r, c)) .and. nint(inp%bc_class(k)) == bc_class_at(r, c) &
              .and. close(inp%f_leak(k), f_leak_at(r, c))                                                         &
              .and. close(inp%deep_temp(k), deep_temp_at(r, c)) .and. close(inp%chave_e(k), chave_e_at(r, c))       &
              .and. close(inp%texture(3, k)%clay, t5%clay) .and. close(inp%texture(6, k)%ph, t6%ph)
   end do
   call check_true('every value is its own cell''s', ok)

   !----- Cell (1, 0): 1.1 m deep, so 5 layers, the last 0.6-1.1 m: 0.4 m of map layer 5 and 0.1 m of 6. !
   do k = 1_ik, masked%ncell
      if (masked%row(k) == 1_ik .and. masked%col(k) == 0_ik) exit
   end do
   call alloc_pft_table(cfg%pft, 1_ik)      ! the cell's E goes into every PFT's allometry
   call apply_map_inputs(inp, k, cfg)
   t5 = texture_at(5_ik, 1_ik, 0_ik) ; t6 = texture_at(6_ik, 1_ik, 0_ik)
   tex%sand = (0.4_wp*t5%sand + 0.1_wp*t6%sand) / 0.5_wp ; tex%silt = (0.4_wp*t5%silt + 0.1_wp*t6%silt) / 0.5_wp
   tex%clay = (0.4_wp*t5%clay + 0.1_wp*t6%clay) / 0.5_wp
   tex%org_carbon = (0.4_wp*t5%org_carbon + 0.1_wp*t6%org_carbon) / 0.5_wp
   tex%cec = (0.4_wp*t5%cec + 0.1_wp*t6%cec) / 0.5_wp ; tex%ph = (0.4_wp*t5%ph + 0.1_wp*t6%ph) / 0.5_wp
   tex%bulk_density = (0.4_wp*t5%bulk_density + 0.1_wp*t6%bulk_density) / 0.5_wp
   tex%coarse_frag  = (0.4_wp*t5%coarse_frag + 0.1_wp*t6%coarse_frag) / 0.5_wp
   call soil_texture_to_hydraulics(tex, ts, tr, al, vn, ks)
   call soil_texture_to_thermal(tex%sand, tex%clay, ts, kso, kd, cv)
   associate (sc => cfg%soil)
      call check_true('the column ends at the map''s depth, on the default layers',                      &
                      close(sc%depth, soil_depth_at(1_ik, 0_ik)) .and. sc%n_layer_values == 5_ik)
      call check_true('a layer across two map layers takes their overlap-weighted soil',                 &
                      close(sc%layer_theta_sat(5), ts) .and. close(sc%layer_theta_res(5), tr) .and.       &
                      close(sc%layer_curve_par_a(5), al) .and. close(sc%layer_curve_par_n(5), vn) .and.   &
                      close(sc%layer_ksat(5), ks))
      call check_true('and its thermal properties from the same soil', close(sc%layer_solid_conductivity(5), kso) &
                      .and. close(sc%layer_dry_conductivity(5), kd) .and. close(sc%layer_dry_heat_capacity(5), cv))
      call soil_texture_to_hydraulics(texture_at(2_ik, 1_ik, 0_ik), ts, tr, al, vn, ks)
      call check_true('a layer on one map layer takes that layer''s soil', close(sc%layer_theta_sat(2), ts) &
                      .and. close(sc%layer_ksat(2), ks))
   end associate
   call check_true('a leaky bottom drains at its cell''s f_leak', bc_class_at(1_ik, 0_ik) == 2_ik .and. &
                   cfg%soil%water_bottom_bc == SOIL_BC_DRAIN .and. cfg%soil%drain_factor == f_leak_at(1_ik, 0_ik) &
                   .and. f_leak_at(1_ik, 0_ik) < 1.0_wp)
   call check_true('the thermal anchor: the map''s temperature, d/sqrt(2) below the bottom node',          &
                   close(cfg%soil%deep_temp, deep_temp_at(1_ik, 0_ik)) .and.                              &
                   close(cfg%soil%deep_depth, 0.85_wp + ANNUAL_DAMPING_DEPTH / sqrt(2.0_wp)))
   call check_true('Chave''s E is the cell''s', close(cfg%pft%allom(1)%chave_e, chave_e_at(1_ik, 0_ik)))

   !----- Rows 2-3, columns 14-17: cell (3, 17) is 2.55 m deep, 8 layers, the last 2.0-2.55 m, below    !
   !      the maps' 2 m: it takes the deepest map layer's soil. Its bottom is (3 + 17) mod 4 + 1 = 1. ---!
   call select_box_cells(grid, [26.0_wp, 99.0_wp, 14.0_wp, 161.0_wp], cells, st)
   call read_map_inputs(opts, grid, cells, inp, msg)
   do k = 1_ik, cells%ncell
      if (cells%row(k) == 3_ik .and. cells%col(k) == 17_ik) exit
   end do
   call check_true('the second box is read', len(msg) == 0 .and. k <= cells%ncell)
   call alloc_pft_table(cfg2%pft, 1_ik)
   call apply_map_inputs(inp, k, cfg2)
   call soil_texture_to_hydraulics(texture_at(6_ik, 3_ik, 17_ik), ts, tr, al, vn, ks)
   call check_true('a layer below the maps'' 2 m takes the deepest map layer''s soil',                     &
                   cfg2%soil%n_layer_values == 8_ik .and. close(cfg2%soil%layer_theta_sat(8), ts)               &
                   .and. close(cfg2%soil%layer_curve_par_n(8), vn))
   call check_true('a free-draining bottom drains at f = 1', bc_class_at(3_ik, 17_ik) == 1_ik .and.     &
                   cfg2%soil%water_bottom_bc == SOIL_BC_DRAIN .and. cfg2%soil%drain_factor == 1.0_wp)
   do k = 1_ik, cells%ncell
      if (bc_class_at(cells%row(k), cells%col(k)) == 3_ik) exit
   end do
   call apply_map_inputs(inp, k, cfg2)
   call check_true('a sealed bottom drains at f = 0, as bedrock',                                       &
                   cfg2%soil%water_bottom_bc == SOIL_BC_DRAIN .and. cfg2%soil%drain_factor == 0.0_wp)
   do k = 1_ik, cells%ncell
      if (bc_class_at(cells%row(k), cells%col(k)) == 4_ik) exit
   end do
   call apply_map_inputs(inp, k, cfg2)
   call check_true('an aquifer bottom is an aquifer', cfg2%soil%water_bottom_bc == SOIL_BC_AQUIFER)

   !----- A leaky base's factor must be a fraction of free drainage: f_leak = 0 is refused, naming   !
   !      the cell. ------------------------------------------------------------------------------------!
   call write_surface_maps(trim(work)//'/surface_f0', f_leak=0.0_wp)
   opts%soil_depth_map     = trim(work)//'/surface_f0/soil_depth.nc'
   call read_map_inputs(opts, grid, cells, inp, msg)
   call check_true('a leaky base with f_leak = 0 is refused, naming the cell',                          &
                   index(msg, 'f_leak is not in (0, 1] at the cell (25.00, 120.00) (row 2, col 15), '//   &
                         'whose base is leaky') > 0)

   !----- A run cell a layered map leaves empty is named too: a texture map without sand at       !
   !      (row 1, col 0), read at the first box's 15 cells, whose depths the depth map gives. -----!
   call write_surface_maps(trim(work)//'/surface_gap', gap=[1_ik, 0_ik])
   opts%soil_depth_map   = trim(work)//'/surface/soil_depth.nc'
   opts%soil_texture_map = trim(work)//'/surface_gap/soil_texture.nc'
   call read_map_inputs(opts, grid, masked, inp, msg)
   call check_true('a run cell a layered map leaves empty is named, with the variable',                  &
                   index(msg, 'sand has no value at the cell (35.00, -180.00) (row 1, col 0)') > 0)

   call test_report('test_maps')

contains

   pure logical function close(a, b)
      real(wp), intent(in) :: a, b
      close = abs(a - b) <= 1.0e-12_wp * max(abs(a), abs(b), 1.0_wp)
   end function close

end program test_maps
