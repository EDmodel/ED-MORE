! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_netcdf_write -- what the two netCDF writers share, the restart file                   !
! (ed_init_restart) and the diagnostic files (ed_output_stream): a variable defined        !
! with its attributes, chunks and compression in one call, text and integer attributes,    !
! the zlib levels, the polygon axis both files lead with, its id and cell centre, and the  !
! per-site statics both carry, under one set of names.                                     !
!==========================================================================================!
module ed_netcdf_write
   use iso_c_binding, only : c_int, c_size_t, c_double, c_char
   use ed_kinds,      only : wp, ik
   use ed_netcdf_c,   only : cstr, nc_check, nc_def_dim_f, nc_def_var_f, nc_def_var_chunking, nc_def_var_deflate, &
                             nc_inq_dimlen, nc_inq_varid_f, nc_put_att_double_f, nc_put_att_int, nc_put_att_int_f, &
                             nc_put_att_text_f, nc_put_vara_double, nc_put_vara_int, NC_CHUNKED, NC_DOUBLE,      &
                             NC_FILL_DOUBLE, NC_FILL_INT, NC_GLOBAL, NC_INT
   implicit none
   private

   public :: define_var, put_text_att, put_int_atts, put_global_text
   public :: define_polygon_axis, write_polygon_axis, define_site_statics, write_site_statics
   public :: RESTART_DEFLATE, OUTPUT_DEFLATE

   !----- The zlib levels: the restart file's variables, and the diagnostic files'. ------------!
   integer(c_int), parameter :: RESTART_DEFLATE = 1_c_int
   integer(c_int), parameter :: OUTPUT_DEFLATE  = 4_c_int

contains

   !----- Define a variable on `dims` (slowest first), with its units and long_name when given,  !
   !      its chunks, and its zlib level. -----------------------------------------------------!
   subroutine define_var(ncid, name, xtype, dims, vid, units, long_name, chunk, deflate)
      integer(c_int),    intent(in)  :: ncid, xtype, dims(:)
      character(len=*),  intent(in)  :: name
      integer(c_int),    intent(out) :: vid
      character(len=*),  intent(in), optional :: units, long_name
      integer(c_size_t), intent(in), optional :: chunk(:)
      integer(c_int),    intent(in), optional :: deflate
      call nc_check(nc_def_var_f(ncid, name, xtype, size(dims, kind=c_int), dims, vid), 'def '//name)
      if (present(chunk))   call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED, chunk), 'chunk '//name)
      if (present(deflate)) call nc_check(nc_def_var_deflate(ncid, vid, 1_c_int, 1_c_int, deflate), 'deflate '//name)
      if (present(units))     call put_text_att(ncid, vid, 'units', units)
      if (present(long_name)) call put_text_att(ncid, vid, 'long_name', long_name)
   end subroutine define_var

   !----- One text attribute on a defined variable. -------------------------------------------!
   subroutine put_text_att(ncid, vid, name, text)
      integer(c_int),   intent(in) :: ncid, vid
      character(len=*), intent(in) :: name, text
      call nc_check(nc_put_att_text_f(ncid, vid, name, int(len_trim(text), c_size_t), text), 'att '//name)
   end subroutine put_text_att

   !----- One integer-array attribute on a defined variable. ----------------------------------!
   subroutine put_int_atts(ncid, vid, name, vals)
      integer(c_int),   intent(in) :: ncid, vid
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: vals(:)
      character(kind=c_char) :: cname(len_trim(name) + 1)
      integer(c_int)         :: iv(size(vals))
      cname = cstr(name)
      iv = int(vals, c_int)
      call nc_check(nc_put_att_int(ncid, vid, cname, NC_INT, int(size(iv), c_size_t), iv), 'att '//name)
   end subroutine put_int_atts

   !----- One global text attribute. ------------------------------------------------------------!
   subroutine put_global_text(ncid, name, text)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name, text
      call nc_check(nc_put_att_text_f(ncid, NC_GLOBAL, name, int(len_trim(text), c_size_t), text), 'global '//name)
   end subroutine put_global_text

   !----- The polygon axis: the `polygon` dimension, n_polygon long, and its id and cell centre, one !
   !      value per polygon; d_polygon is the dimension. ---------------------------------------!
   subroutine define_polygon_axis(ncid, n_polygon, d_polygon)
      integer(c_int), intent(in)  :: ncid
      integer(ik),    intent(in)  :: n_polygon
      integer(c_int), intent(out) :: d_polygon
      integer(c_int) :: v
      call nc_check(nc_def_dim_f(ncid, 'polygon', int(n_polygon, c_size_t), d_polygon), 'dim polygon')
      call define_var(ncid, 'polygon_id', NC_INT, [d_polygon], v,                                  &
                      long_name='forcing cell row-major index: row * n_lon + col')
      call define_var(ncid, 'lat', NC_DOUBLE, [d_polygon], v, units='degrees_north')
      call put_text_att(ncid, v, 'standard_name', 'latitude')
      call define_var(ncid, 'lon', NC_DOUBLE, [d_polygon], v, units='degrees_east')
      call put_text_att(ncid, v, 'standard_name', 'longitude')
   end subroutine define_polygon_axis

   !----- The polygon axis's values, every polygon at once, after the definitions end. ----------!
   subroutine write_polygon_axis(ncid, polygon_id, lat, lon)
      integer(c_int), intent(in) :: ncid
      integer(ik),    intent(in) :: polygon_id(:)
      real(wp),       intent(in) :: lat(:), lon(:)   !< [deg]
      integer(c_size_t) :: start1(1), count1(1)
      integer(c_int)    :: v, ids(size(polygon_id))
      real(c_double)    :: x(size(lat))
      start1 = 0_c_size_t
      count1 = int(size(polygon_id), c_size_t)
      !----- Each through a named array: no function result goes straight into a call (CLAUDE.md). -!
      ids = int(polygon_id, c_int)
      call nc_check(nc_inq_varid_f(ncid, 'polygon_id', v), 'inq polygon_id')
      call nc_check(nc_put_vara_int(ncid, v, start1, count1, ids), 'put polygon_id')
      x = real(lat, c_double)
      call nc_check(nc_inq_varid_f(ncid, 'lat', v), 'inq lat')
      call nc_check(nc_put_vara_double(ncid, v, start1, count1, x), 'put lat')
      x = real(lon, c_double)
      call nc_check(nc_inq_varid_f(ncid, 'lon', v), 'inq lon')
      call nc_check(nc_put_vara_double(ncid, v, start1, count1, x), 'put lon')
   end subroutine write_polygon_axis

   !---------------------------------------------------------------------------------------!
   ! The per-site statics, on (polygon, site): each site's area, terrain, soil column and   !
   ! its bottom, as the run runs them, with netCDF's default fill where a polygon has no   !
   ! such site and for the drainage factor of an aquifer bottom. bc_codes and bc_names are  !
   ! site_bottom_bc's flag values and meanings. by_row chunks one polygon's row, for a      !
   ! writer that writes a row at a time.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine define_site_statics(ncid, d_polygon, d_site, bc_codes, bc_names, by_row, deflate)
      integer(c_int),           intent(in) :: ncid, d_polygon, d_site
      integer(ik),              intent(in) :: bc_codes(:)
      character(len=*),         intent(in) :: bc_names
      logical,                  intent(in) :: by_row
      integer(c_int), optional, intent(in) :: deflate
      integer(c_size_t) :: chunk(2)
      integer(c_int)    :: v
      chunk(1) = 1_c_size_t
      call nc_check(nc_inq_dimlen(ncid, d_site, chunk(2)), 'dimlen site')
      call def_static('area_si',       NC_DOUBLE, '1', 'site area fraction of the polygon')
      call def_static('elevation_si',  NC_DOUBLE, 'm', 'site elevation')
      call def_static('slope_si',      NC_DOUBLE, 'degree', 'site slope')
      call def_static('aspect_si',     NC_DOUBLE, 'degree', 'site aspect, clockwise from north')
      call def_static('n_soil_si',     NC_INT,    '1', 'the site''s active soil layers')
      call def_static('soil_depth_si', NC_DOUBLE, 'm', 'soil column depth, the bottom of the deepest layer')
      call def_static('bottom_bc_si',  NC_INT,    '1',                                                    &
                      'bottom water boundary as run: drainage at site_drain_factor, or aquifer')
      call put_int_atts(ncid, v, 'flag_values', bc_codes)
      call put_text_att(ncid, v, 'flag_meanings', bc_names)
      call def_static('drain_factor_si', NC_DOUBLE, '1', 'drainage factor f of the bottom flux q = f K: 1 free '// &
                      'drainage, the map''s f_leak on a leaky base, 0 sealed (fill on an aquifer)')
      call def_static('deep_temp_si',  NC_DOUBLE, 'K', 'deep soil temperature, the thermal anchor below the column')
   contains
      subroutine def_static(name, xtype, units, lname)
         character(len=*), intent(in) :: name, units, lname
         integer(c_int),   intent(in) :: xtype
         if (by_row) then
            call define_var(ncid, name, xtype, [d_polygon, d_site], v, units=units, long_name=lname,         &
                            chunk=chunk, deflate=deflate)
         else
            call define_var(ncid, name, xtype, [d_polygon, d_site], v, units=units, long_name=lname,         &
                            deflate=deflate)
         end if
         if (xtype == NC_INT) then
            call nc_check(nc_put_att_int_f(ncid, v, '_FillValue', NC_INT, NC_FILL_INT), 'fill '//name)
         else
            call nc_check(nc_put_att_double_f(ncid, v, '_FillValue', NC_DOUBLE, NC_FILL_DOUBLE), 'fill '//name)
         end if
      end subroutine def_static
   end subroutine define_site_statics

   !----- The per-site statics of the polygons from row row0 (0-based) on, each array (site,     !
   !      polygon); after the definitions end. -------------------------------------------------!
   subroutine write_site_statics(ncid, row0, area, elevation, slope, aspect, n_soil, soil_depth, bottom_bc, &
                                 drain_factor, deep_temp)
      integer(c_int), intent(in) :: ncid
      integer(ik),    intent(in) :: row0
      real(wp),       intent(in) :: area(:,:), elevation(:,:), slope(:,:), aspect(:,:), soil_depth(:,:)
      real(wp),       intent(in) :: drain_factor(:,:), deep_temp(:,:)
      integer(ik),    intent(in) :: n_soil(:,:), bottom_bc(:,:)
      integer(c_size_t) :: st(2), ct(2)
      st = [int(row0, c_size_t), 0_c_size_t]
      ct = [int(size(area, 2), c_size_t), int(size(area, 1), c_size_t)]
      call put_dbl('area_si', area)
      call put_dbl('elevation_si', elevation)
      call put_dbl('slope_si', slope)
      call put_dbl('aspect_si', aspect)
      call put_int('n_soil_si', n_soil)
      call put_dbl('soil_depth_si', soil_depth)
      call put_int('bottom_bc_si', bottom_bc)
      call put_dbl('drain_factor_si', drain_factor)
      call put_dbl('deep_temp_si', deep_temp)
   contains
      !----- Each through a named array: no function result goes straight into a call (CLAUDE.md). -!
      subroutine put_dbl(name, x)
         character(len=*), intent(in) :: name
         real(wp),         intent(in) :: x(:,:)
         real(c_double) :: buf(size(x, 1), size(x, 2))
         integer(c_int) :: v
         buf = real(x, c_double)
         call nc_check(nc_inq_varid_f(ncid, name, v), 'inq '//name)
         call nc_check(nc_put_vara_double(ncid, v, st, ct, buf), 'put '//name)
      end subroutine put_dbl
      subroutine put_int(name, x)
         character(len=*), intent(in) :: name
         integer(ik),      intent(in) :: x(:,:)
         integer(c_int) :: buf(size(x, 1), size(x, 2))
         integer(c_int) :: v
         buf = int(x, c_int)
         call nc_check(nc_inq_varid_f(ncid, name, v), 'inq '//name)
         call nc_check(nc_put_vara_int(ncid, v, st, ct, buf), 'put '//name)
      end subroutine put_int
   end subroutine write_site_statics

end module ed_netcdf_write
