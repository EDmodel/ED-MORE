! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_met_grid -- where on the Earth the met data are: the cells of a gridded product a
! run reads, from the product's static file (met_grid_t), and whether another file lies on
! that grid (docs/science/forcing.md §1). select_site_cell takes a site's cell,
! select_box_cells a region's, keep_cells drops the cells a mask leaves out, and same_grid
! checks a map's axes, since ED-MORE never regrids. It reads the static file, never
! weather.
!---~---
module ed_met_grid
   use iso_c_binding,  only : c_int, c_size_t
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : deg_to_rad
   use ed_met_types,   only : met_grid_t, met_cells_t
   use ed_netcdf_c,    only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f, nc_get_vara_double, nc_close,   &
                              NC_NOERR, NC_NOWRITE
   implicit none
   private

   public :: select_site_cell, select_box_cells, keep_cells, same_grid
   public :: great_circle_distance, nearest_grid_index
   public :: GRID_OK, GRID_ERR_OPEN, GRID_ERR_LAYOUT, GRID_ERR_NO_CELL

   real(wp), parameter :: EARTH_RADIUS_M = 6.371e6_wp  ! [m] mean Earth radius
   real(wp), parameter :: GRID_TOL = 1.0e-4_wp         ! [deg] slack for "these coordinates are equal"

   !--- Status codes.
   integer(ik), parameter :: GRID_OK          = 0_ik
   integer(ik), parameter :: GRID_ERR_OPEN    = 1_ik   ! a file is missing or unreadable
   integer(ik), parameter :: GRID_ERR_LAYOUT  = 2_ik   ! no regular lat/lon axes, another grid's, or a field missing
   integer(ik), parameter :: GRID_ERR_NO_CELL = 3_ik   ! no valid cell within the distance limit / in the box

contains

   !---~---
   !   A SITE'S CELL (docs/science/forcing.md §1): the cell containing the site by
   ! regular-grid arithmetic; if it has no data, the nearest valid cell by great-circle
   ! distance, lowest (row-major) index on a tie, up to max_km. Rows farther than max_km
   ! in latitude alone cannot qualify, so only the rows within that band are read, across
   ! all longitudes (which handles the antimeridian for free).
   !---~---
   subroutine select_site_cell(grid, site_lat, site_lon, max_km, cells, distance_km, stat)
      type(met_grid_t),  intent(in)  :: grid
      real(wp),          intent(in)  :: site_lat, site_lon, max_km
      type(met_cells_t), intent(out) :: cells
      real(wp),          intent(out) :: distance_km
      integer(ik),       intent(out) :: stat
      integer(c_int) :: ncid, st
      real(wp), allocatable :: lat(:), lon(:), valid(:,:), elev(:,:)
      real(wp)    :: row_km, d, dlat
      integer(ik) :: j0, rj, r0, r1, j, i, best_j, best_i

      distance_km = huge(1.0_wp)
      call open_grid_file(grid%static_file, grid, ncid, lat, lon, stat)
      if (stat /= GRID_OK) return
      dlat = lat(2) - lat(1)
      j0 = min(size(lat, kind=ik) - 1_ik, max(0_ik, nint((site_lat - lat(1)) / dlat, ik)))
      row_km = great_circle_distance(0.0_wp, lat(1), 0.0_wp, lat(2)) / 1000.0_wp
      rj = ceiling(max_km / row_km, ik) + 1_ik
      r0 = max(0_ik, j0 - rj)
      r1 = min(size(lat, kind=ik) - 1_ik, j0 + rj)
      call read_2d(ncid, grid%valid_var, r0, r1 - r0 + 1_ik, 0_ik, size(lon, kind=ik), valid, stat)
      if (stat /= GRID_OK) then
         st = nc_close(ncid)
         return
      end if

      best_j = -1_ik
      best_i = -1_ik
      do j = r0, r1
         do i = 0_ik, size(lon, kind=ik) - 1_ik
            if (valid(i + 1_ik, j - r0 + 1_ik) < 0.5_wp) cycle
            d = great_circle_distance(site_lon, site_lat, lon(i + 1_ik), lat(j + 1_ik)) / 1000.0_wp
            if (d < distance_km) then
               distance_km = d
               best_j = j
               best_i = i
            end if
         end do
      end do
      if (best_j < 0_ik .or. distance_km > max_km) then
         st = nc_close(ncid)
         stat = GRID_ERR_NO_CELL
         return
      end if
      call read_2d(ncid, grid%elevation_var, best_j, 1_ik, best_i, 1_ik, elev, stat)
      st = nc_close(ncid)
      if (stat /= GRID_OK) return

      cells%nlat = size(lat, kind=ik)
      cells%nlon = size(lon, kind=ik)
      cells%ncell = 1_ik
      cells%row = [best_j]
      cells%col = [best_i]
      cells%lat = [lat(best_j + 1_ik)]
      cells%lon = [lon(best_i + 1_ik)]
      cells%elevation = [elev(1, 1)]
      cells%chunk = grid%chunk
      call group_by_chunk(cells)
   end subroutine select_site_cell

   !---~---
   !   A REGION'S CELLS (docs/configuration.md, Regional runs): the valid cells whose
   ! centres lie in [south, north] x [west, east], in row-major order north to south and
   ! west to east. A box with west > east crosses the antimeridian: its columns run from
   ! west to the grid's east edge, then on from the west edge to east. nwse = [north,
   ! west, south, east] in degrees.
   !---~---
   subroutine select_box_cells(grid, nwse, cells, stat, land_fraction_min)
      type(met_grid_t),  intent(in)  :: grid
      real(wp),          intent(in)  :: nwse(4)
      type(met_cells_t), intent(out) :: cells
      integer(ik),       intent(out) :: stat
      !---~---
      !   Keep only cells whose static land fraction is at least this (a region's selection
      ! rule, docs/configuration.md, Regional runs); absent, every valid cell is kept.
      !---~---
      real(wp), optional, intent(in)  :: land_fraction_min
      integer(c_int) :: ncid, st
      real(wp), allocatable :: lat(:), lon(:), valid(:,:), elev(:,:), lfrac(:,:)
      real(wp)    :: lf_min
      integer(ik), allocatable :: cols(:), col_index(:)
      real(wp)    :: west, east
      integer(ik) :: nlat, nlon, r0, r1, j, i, k, n

      call open_grid_file(grid%static_file, grid, ncid, lat, lon, stat)
      if (stat /= GRID_OK) return
      nlat = size(lat, kind=ik)
      nlon = size(lon, kind=ik)
      r0 = -1_ik
      r1 = -1_ik
      do j = 0_ik, nlat - 1_ik
         if (lat(j + 1_ik) <= nwse(1) + GRID_TOL .and. lat(j + 1_ik) >= nwse(3) - GRID_TOL) then
            if (r0 < 0_ik) r0 = j
            r1 = j
         end if
      end do
      west = wrap_longitude(nwse(2), lon(1))
      east = wrap_longitude(nwse(4), lon(1))
      !--- A named index array, not an implied-do inside pack (compiler trap, docs/building.md).
      col_index = [(i, i = 0_ik, nlon - 1_ik)]
      if (west <= east) then
         cols = pack(col_index, lon >= west - GRID_TOL .and. lon <= east + GRID_TOL)
      else
         cols = [pack(col_index, lon >= west - GRID_TOL), pack(col_index, lon <= east + GRID_TOL)]
      end if
      if (r0 < 0_ik .or. size(cols) == 0) then
         st = nc_close(ncid)
         stat = GRID_ERR_NO_CELL
         return
      end if
      call read_2d(ncid, grid%valid_var, r0, r1 - r0 + 1_ik, 0_ik, nlon, valid, stat)
      if (stat == GRID_OK) call read_2d(ncid, grid%elevation_var, r0, r1 - r0 + 1_ik, 0_ik, nlon, elev, stat)
      if (stat == GRID_OK) call read_2d(ncid, grid%land_fraction_var, r0, r1 - r0 + 1_ik, 0_ik, nlon, lfrac, stat)
      st = nc_close(ncid)
      if (stat /= GRID_OK) return

      lf_min = -1.0_wp
      if (present(land_fraction_min)) lf_min = land_fraction_min
      where (lfrac < lf_min) valid = 0.0_wp
      n = int(count(valid(cols + 1_ik, :) > 0.5_wp), ik)
      if (n == 0_ik) then
         stat = GRID_ERR_NO_CELL
         return
      end if
      cells%nlat = nlat
      cells%nlon = nlon
      cells%ncell = n
      allocate(cells%row(n), cells%col(n), cells%lat(n), cells%lon(n), cells%elevation(n), cells%land_fraction(n))
      k = 0_ik
      do j = r0, r1
         do i = 1_ik, size(cols, kind=ik)
            if (valid(cols(i) + 1_ik, j - r0 + 1_ik) < 0.5_wp) cycle
            k = k + 1_ik
            cells%row(k) = j
            cells%col(k) = cols(i)
            cells%lat(k) = lat(j + 1_ik)
            cells%lon(k) = lon(cols(i) + 1_ik)
            cells%elevation(k) = elev(cols(i) + 1_ik, j - r0 + 1_ik)
            cells%land_fraction(k) = lfrac(cols(i) + 1_ik, j - r0 + 1_ik)
         end do
      end do
      cells%chunk = grid%chunk
      call group_by_chunk(cells)
   end subroutine select_box_cells

   !---~---
   !   Keep the cells where `keep` is true, in their order, and regroup them by chunk (a
   ! region's cell mask, ed_maps).
   !---~---
   subroutine keep_cells(cells, keep)
      type(met_cells_t), intent(inout) :: cells
      logical,           intent(in)    :: keep(:)
      cells%row = pack(cells%row, keep)
      cells%col = pack(cells%col, keep)
      cells%lat = pack(cells%lat, keep)
      cells%lon = pack(cells%lon, keep)
      cells%elevation = pack(cells%elevation, keep)
      if (allocated(cells%land_fraction)) cells%land_fraction = pack(cells%land_fraction, keep)
      cells%ncell = int(count(keep), ik)
      deallocate(cells%chunk_row, cells%chunk_col, cells%chunk_first, cells%chunk_nrow, cells%chunk_ncol,     &
                 cells%by_chunk)
      call group_by_chunk(cells)
   end subroutine keep_cells

   !---~---
   !   ON THE GRID: GRID_OK when the file at `path` has the static file's lat and lon axes
   ! (within GRID_TOL), so a cell's (row, col) indexes both. The static maps a run reads
   ! (ed_maps) must be on this grid: ED-MORE never regrids.
   !---~---
   subroutine same_grid(path, grid, stat)
      character(len=*), intent(in)  :: path
      type(met_grid_t), intent(in)  :: grid
      integer(ik),      intent(out) :: stat
      integer(c_int) :: ncid, st
      real(wp), allocatable :: lat_s(:), lon_s(:), lat_m(:), lon_m(:)
      call open_grid_file(grid%static_file, grid, ncid, lat_s, lon_s, stat)
      if (stat /= GRID_OK) return
      st = nc_close(ncid)
      call open_grid_file(path, grid, ncid, lat_m, lon_m, stat)
      if (stat /= GRID_OK) return
      st = nc_close(ncid)
      stat = GRID_ERR_LAYOUT
      if (size(lat_m) /= size(lat_s) .or. size(lon_m) /= size(lon_s)) return
      if (any(abs(lat_m - lat_s) > GRID_TOL) .or. any(abs(lon_m - lon_s) > GRID_TOL)) return
      stat = GRID_OK
   end subroutine same_grid

   !---~---
   !   DISTANCE (docs/science/forcing.md §1): the great-circle distance, and the nearest of
   ! a list of points with a strict '<', so a tie keeps the lowest index. A site's cell
   ! is held to [forcing].max_distance_km; the mean Earth radius keeps it within 0.5 %.
   !---~---
   pure function great_circle_distance(lon1_deg, lat1_deg, lon2_deg, lat2_deg) result(dist_m)
      real(wp), intent(in) :: lon1_deg, lat1_deg, lon2_deg, lat2_deg
      real(wp) :: dist_m, lat1, lat2, dlon, x, y
      lat1 = lat1_deg * deg_to_rad
      lat2 = lat2_deg * deg_to_rad
      dlon = (lon2_deg - lon1_deg) * deg_to_rad
      x = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon)
      y = sqrt( (cos(lat2) * sin(dlon)) ** 2                                                    &
              + (cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)) ** 2 )
      dist_m = EARTH_RADIUS_M * atan2(y, x)
   end function great_circle_distance

   pure function nearest_grid_index(site_lon, site_lat, lon, lat) result(idx)
      real(wp), intent(in) :: site_lon, site_lat, lon(:), lat(:)
      integer(ik) :: idx, i
      real(wp) :: dmin, d
      idx = 1_ik
      dmin = huge(1.0_wp)
      do i = 1_ik, int(size(lon), ik)
         d = great_circle_distance(site_lon, site_lat, lon(i), lat(i))
         if (d < dmin) then
            dmin = d
            idx = i
         end if
      end do
   end function nearest_grid_index

   !--- A longitude moved into [lon0, lon0 + 360), the grid's own range.
   pure real(wp) function wrap_longitude(x, lon0) result(y)
      real(wp), intent(in) :: x, lon0
      y = lon0 + modulo(x - lon0, 360.0_wp)
   end function wrap_longitude

   !---~---
   !   Group the cells by the grid's square chunk with a stable counting sort, so each touched
   ! chunk column is read once per variable-month and the cells keep their order. Each chunk's
   ! read box is the rows and columns its cells occupy: HDF5 decompresses the whole chunk
   ! either way, but a site then copies one cell instead of the whole chunk.
   !---~---
   subroutine group_by_chunk(cells)
      type(met_cells_t), intent(inout) :: cells
      integer(ik), allocatable :: id(:), cnt(:), next(:)
      integer(ik) :: ncx, ncy, c, k, m
      ncx = (cells%nlon + cells%chunk - 1_ik) / cells%chunk
      ncy = (cells%nlat + cells%chunk - 1_ik) / cells%chunk
      allocate(id(cells%ncell), cnt(0:ncx*ncy))
      id = (cells%row / cells%chunk) * ncx + cells%col / cells%chunk
      cnt = 0_ik
      do c = 1_ik, cells%ncell
         cnt(id(c)) = cnt(id(c)) + 1_ik
      end do
      cells%nchunk = int(count(cnt(0:ncx*ncy - 1_ik) > 0_ik), ik)
      allocate(cells%chunk_row(cells%nchunk), cells%chunk_col(cells%nchunk), cells%chunk_first(cells%nchunk + 1_ik))
      allocate(cells%chunk_nrow(cells%nchunk), cells%chunk_ncol(cells%nchunk))
      allocate(cells%by_chunk(cells%ncell), next(0:ncx*ncy - 1_ik))
      k = 0_ik
      m = 1_ik
      do c = 0_ik, ncx*ncy - 1_ik
         if (cnt(c) == 0_ik) cycle
         k = k + 1_ik
         cells%chunk_first(k) = m
         next(c) = m
         m = m + cnt(c)
      end do
      cells%chunk_first(cells%nchunk + 1_ik) = m
      do c = 1_ik, cells%ncell
         cells%by_chunk(next(id(c))) = c
         next(id(c)) = next(id(c)) + 1_ik
      end do
      do k = 1_ik, cells%nchunk
         associate (members => cells%by_chunk(cells%chunk_first(k):cells%chunk_first(k + 1_ik) - 1_ik))
            cells%chunk_row(k)  = minval(cells%row(members))
            cells%chunk_col(k)  = minval(cells%col(members))
            cells%chunk_nrow(k) = maxval(cells%row(members)) - cells%chunk_row(k) + 1_ik
            cells%chunk_ncol(k) = maxval(cells%col(members)) - cells%chunk_col(k) + 1_ik
         end associate
      end do
   end subroutine group_by_chunk

   !---~---
   !   Open a file on the grid and read its (regular, 1-D) lat and lon axes; closed again unless
   ! stat = GRID_OK.
   !---~---
   subroutine open_grid_file(path, grid, ncid, lat, lon, stat)
      character(len=*),      intent(in)  :: path
      type(met_grid_t),      intent(in)  :: grid
      integer(c_int),        intent(out) :: ncid
      real(wp), allocatable, intent(out) :: lat(:), lon(:)
      integer(ik),           intent(out) :: stat
      integer(c_int) :: st
      stat = GRID_ERR_OPEN
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) return
      stat = GRID_ERR_LAYOUT
      call read_axis(ncid, trim(grid%lat_var), lat, st)
      if (st == NC_NOERR) call read_axis(ncid, trim(grid%lon_var), lon, st)
      if (st == NC_NOERR) then
         if (size(lat) >= 2 .and. size(lon) >= 2) then            ! a regular grid
            if (all(abs((lat(2:) - lat(:size(lat) - 1)) - (lat(2) - lat(1))) <= GRID_TOL) .and.        &
                all(abs((lon(2:) - lon(:size(lon) - 1)) - (lon(2) - lon(1))) <= GRID_TOL)) stat = GRID_OK
         end if
      end if
      if (stat /= GRID_OK) st = nc_close(ncid)
   end subroutine open_grid_file

   subroutine read_axis(ncid, name, x, st)
      integer(c_int),        intent(in)  :: ncid
      character(len=*),      intent(in)  :: name
      real(wp), allocatable, intent(out) :: x(:)
      integer(c_int),        intent(out) :: st
      integer(c_int)    :: vid
      integer(c_size_t) :: n, start1(1), count1(1)
      st = nc_inq_dimlen_f(ncid, name, n)
      if (st /= NC_NOERR) return
      st = nc_inq_varid_f(ncid, name, vid)
      if (st /= NC_NOERR) return
      allocate(x(n))
      start1 = 0_c_size_t
      count1 = n
      st = nc_get_vara_double(ncid, vid, start1, count1, x)
   end subroutine read_axis

   !--- rows r0 .. r0+nr-1 and columns c0 .. c0+nc-1 (0-based) of a (lat, lon) variable, as x(lon, lat).
   subroutine read_2d(ncid, name, r0, nr, c0, nc, x, stat)
      integer(c_int),        intent(in)  :: ncid
      character(len=*),      intent(in)  :: name
      integer(ik),           intent(in)  :: r0, nr, c0, nc
      real(wp), allocatable, intent(out) :: x(:,:)
      integer(ik),           intent(out) :: stat
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start2(2), count2(2)
      stat = GRID_ERR_LAYOUT
      st = nc_inq_varid_f(ncid, trim(name), vid)
      if (st /= NC_NOERR) return
      allocate(x(nc, nr))
      start2 = [int(r0, c_size_t), int(c0, c_size_t)]
      count2 = [int(nr, c_size_t), int(nc, c_size_t)]
      st = nc_get_vara_double(ncid, vid, start2, count2, x)
      if (st /= NC_NOERR) return
      stat = GRID_OK
   end subroutine read_2d

end module ed_met_grid
