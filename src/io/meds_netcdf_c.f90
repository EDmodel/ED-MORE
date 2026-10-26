! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_netcdf_c -- a thin ISO_C_BINDING layer over the netCDF C library, for the state     !
! file, the forcing and map readers and the output streams. netCDF-Fortran's .mod files    !
! are compiler-specific, while the C library links from any compiler, so MEDS binds the C  !
! entry points it uses directly. C counts from 0 (callers pass 0-based starts), lists      !
! dimensions slowest-varying first ([time, cohort]), and takes sizes as size_t and         !
! NUL-terminated strings (cstr); the `_f` wrappers take Fortran strings. The readers share !
! one rule for a variable's gaps and packing (nc_packing_of, nc_unpacked).                 !
!==========================================================================================!
module meds_netcdf_c
   use iso_c_binding, only : c_int, c_size_t, c_double, c_float, c_char, c_ptr, c_null_char,           &
                             c_f_pointer, c_associated
   use meds_kinds,    only : ik
   implicit none
   public

   !----- netCDF constants (values from netcdf.h, v4.9.2). --------------------------------!
   integer(c_int), parameter :: NC_NOERR     = 0_c_int
   integer(c_int), parameter :: NC_NETCDF4   = 4096_c_int   ! 0x1000: HDF5-backed format
   integer(c_int), parameter :: NC_CLOBBER   = 0_c_int
   integer(c_int), parameter :: NC_NOWRITE   = 0_c_int      ! read-only open mode
   integer(c_int), parameter :: NC_GLOBAL    = -1_c_int
   integer(c_int), parameter :: NC_BYTE      = 1_c_int
   integer(c_int), parameter :: NC_INT       = 4_c_int
   integer(c_int), parameter :: NC_FLOAT     = 5_c_int
   integer(c_int), parameter :: NC_DOUBLE    = 6_c_int
   integer(c_int), parameter :: NC_CHUNKED   = 0_c_int      ! storage mode, nc_def/inq_var_chunking
   integer(c_size_t), parameter :: NC_UNLIMITED = 0_c_size_t
   !----- netCDF's default fill values: what a variable holds where nothing was written. --------!
   real(c_double),    parameter :: NC_FILL_DOUBLE = 9.9692099683868690e+36_c_double
   integer(c_int),    parameter :: NC_FILL_INT    = -2147483647_c_int

   !----- A variable's gap and packing attributes, NaN where it has none (nc_packing_of). ----------!
   type :: nc_packing_t
      real(c_double) :: fill     !< _FillValue
      real(c_double) :: missing  !< missing_value
      real(c_double) :: scale    !< scale_factor
      real(c_double) :: offset   !< add_offset
   end type nc_packing_t

   interface
      integer(c_int) function nc_create(path, cmode, ncidp) bind(c, name="nc_create")
         import :: c_char, c_int
         character(kind=c_char), intent(in)  :: path(*)
         integer(c_int), value, intent(in)   :: cmode
         integer(c_int),        intent(out)  :: ncidp
      end function nc_create

      integer(c_int) function nc_def_dim(ncid, name, length, idp) bind(c, name="nc_def_dim")
         import :: c_char, c_int, c_size_t
         integer(c_int), value,    intent(in)  :: ncid
         character(kind=c_char),   intent(in)  :: name(*)
         integer(c_size_t), value, intent(in)  :: length
         integer(c_int),           intent(out) :: idp
      end function nc_def_dim

      integer(c_int) function nc_def_var(ncid, name, xtype, ndims, dimids, varidp)           &
                              bind(c, name="nc_def_var")
         import :: c_char, c_int
         integer(c_int), value,  intent(in)  :: ncid
         character(kind=c_char), intent(in)  :: name(*)
         integer(c_int), value,  intent(in)  :: xtype
         integer(c_int), value,  intent(in)  :: ndims
         integer(c_int),         intent(in)  :: dimids(*)
         integer(c_int),         intent(out) :: varidp
      end function nc_def_var

      integer(c_int) function nc_put_att_text(ncid, varid, name, length, tp)                 &
                              bind(c, name="nc_put_att_text")
         import :: c_char, c_int, c_size_t
         integer(c_int), value,    intent(in) :: ncid, varid
         character(kind=c_char),   intent(in) :: name(*)
         integer(c_size_t), value, intent(in) :: length
         character(kind=c_char),   intent(in) :: tp(*)
      end function nc_put_att_text

      integer(c_int) function nc_put_att_double(ncid, varid, name, xtype, length, op)        &
                              bind(c, name="nc_put_att_double")
         import :: c_char, c_int, c_size_t, c_double
         integer(c_int), value,    intent(in) :: ncid, varid, xtype
         character(kind=c_char),   intent(in) :: name(*)
         integer(c_size_t), value, intent(in) :: length
         real(c_double),           intent(in) :: op(*)
      end function nc_put_att_double

      integer(c_int) function nc_put_att_int(ncid, varid, name, xtype, length, op)           &
                              bind(c, name="nc_put_att_int")
         import :: c_char, c_int, c_size_t
         integer(c_int), value,    intent(in) :: ncid, varid, xtype
         character(kind=c_char),   intent(in) :: name(*)
         integer(c_size_t), value, intent(in) :: length
         integer(c_int),           intent(in) :: op(*)
      end function nc_put_att_int

      integer(c_int) function nc_sync(ncid) bind(c, name="nc_sync")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid
      end function nc_sync

      integer(c_int) function nc_def_var_deflate(ncid, varid, shuffle, deflate, level)       &
                              bind(c, name="nc_def_var_deflate")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid, varid, shuffle, deflate, level
      end function nc_def_var_deflate

      integer(c_int) function nc_def_var_chunking(ncid, varid, storage, chunks)              &
                              bind(c, name="nc_def_var_chunking")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in) :: ncid, varid, storage
         integer(c_size_t),     intent(in) :: chunks(*)
      end function nc_def_var_chunking

      integer(c_int) function nc_inq_var_chunking(ncid, varid, storagep, chunksizesp)        &
                              bind(c, name="nc_inq_var_chunking")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in)  :: ncid, varid
         integer(c_int),        intent(out) :: storagep
         integer(c_size_t),     intent(out) :: chunksizesp(*)
      end function nc_inq_var_chunking

      integer(c_int) function nc_enddef(ncid) bind(c, name="nc_enddef")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid
      end function nc_enddef

      integer(c_int) function nc_put_vara_double(ncid, varid, startp, countp, op)            &
                              bind(c, name="nc_put_vara_double")
         import :: c_int, c_size_t, c_double
         integer(c_int), value, intent(in) :: ncid, varid
         integer(c_size_t),     intent(in) :: startp(*), countp(*)
         real(c_double),        intent(in) :: op(*)
      end function nc_put_vara_double

      integer(c_int) function nc_put_vara_int(ncid, varid, startp, countp, op)               &
                              bind(c, name="nc_put_vara_int")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in) :: ncid, varid
         integer(c_size_t),     intent(in) :: startp(*), countp(*)
         integer(c_int),        intent(in) :: op(*)
      end function nc_put_vara_int

      integer(c_int) function nc_put_var1_double(ncid, varid, indexp, op)                     &
                              bind(c, name="nc_put_var1_double")
         import :: c_int, c_size_t, c_double
         integer(c_int), value, intent(in) :: ncid, varid
         integer(c_size_t),     intent(in) :: indexp(*)
         real(c_double),        intent(in) :: op
      end function nc_put_var1_double

      integer(c_int) function nc_close(ncid) bind(c, name="nc_close")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid
      end function nc_close

      !----- Read side (the state reader, the forcing and map readers). --------------------!
      integer(c_int) function nc_open(path, omode, ncidp) bind(c, name="nc_open")
         import :: c_char, c_int
         character(kind=c_char), intent(in)  :: path(*)
         integer(c_int), value,  intent(in)  :: omode
         integer(c_int),         intent(out) :: ncidp
      end function nc_open

      integer(c_int) function nc_inq_varid(ncid, name, varidp) bind(c, name="nc_inq_varid")
         import :: c_char, c_int
         integer(c_int), value,  intent(in)  :: ncid
         character(kind=c_char), intent(in)  :: name(*)
         integer(c_int),         intent(out) :: varidp
      end function nc_inq_varid

      integer(c_int) function nc_get_vara_int(ncid, varid, startp, countp, vals)             &
                              bind(c, name="nc_get_vara_int")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in)  :: ncid, varid
         integer(c_size_t),     intent(in)  :: startp(*), countp(*)
         integer(c_int),        intent(out) :: vals(*)
      end function nc_get_vara_int

      integer(c_int) function nc_get_vara_double(ncid, varid, startp, countp, vals)          &
                              bind(c, name="nc_get_vara_double")
         import :: c_int, c_size_t, c_double
         integer(c_int), value, intent(in)  :: ncid, varid
         integer(c_size_t),     intent(in)  :: startp(*), countp(*)
         real(c_double),        intent(out) :: vals(*)
      end function nc_get_vara_double

      !----- Float reads, for files that store float32 (the ED_ERA5land archive): no conversion. !
      integer(c_int) function nc_get_vara_float(ncid, varid, startp, countp, vals)           &
                              bind(c, name="nc_get_vara_float")
         import :: c_int, c_size_t, c_float
         integer(c_int), value, intent(in)  :: ncid, varid
         integer(c_size_t),     intent(in)  :: startp(*), countp(*)
         real(c_float),         intent(out) :: vals(*)
      end function nc_get_vara_float

      !----- Dimension + attribute inquiry (the forcing, map and state readers). -----------!
      integer(c_int) function nc_inq_dimid(ncid, name, idp) bind(c, name="nc_inq_dimid")
         import :: c_char, c_int
         integer(c_int), value,  intent(in)  :: ncid
         character(kind=c_char), intent(in)  :: name(*)
         integer(c_int),         intent(out) :: idp
      end function nc_inq_dimid

      integer(c_int) function nc_inq_dimlen(ncid, dimid, lenp) bind(c, name="nc_inq_dimlen")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in)  :: ncid, dimid
         integer(c_size_t),     intent(out) :: lenp
      end function nc_inq_dimlen

      integer(c_int) function nc_inq_attlen(ncid, varid, name, lenp) bind(c, name="nc_inq_attlen")
         import :: c_char, c_int, c_size_t
         integer(c_int), value,  intent(in)  :: ncid, varid
         character(kind=c_char), intent(in)  :: name(*)
         integer(c_size_t),      intent(out) :: lenp
      end function nc_inq_attlen

      integer(c_int) function nc_get_att_text(ncid, varid, name, tp) bind(c, name="nc_get_att_text")
         import :: c_char, c_int
         integer(c_int), value,  intent(in)  :: ncid, varid
         character(kind=c_char), intent(in)  :: name(*)
         character(kind=c_char), intent(out) :: tp(*)
      end function nc_get_att_text

      integer(c_int) function nc_get_att_double(ncid, varid, name, ip) bind(c, name="nc_get_att_double")
         import :: c_char, c_int, c_double
         integer(c_int), value,  intent(in)  :: ncid, varid
         character(kind=c_char), intent(in)  :: name(*)
         real(c_double),         intent(out) :: ip(*)
      end function nc_get_att_double

      type(c_ptr) function nc_strerror(ncerr) bind(c, name="nc_strerror")
         import :: c_int, c_ptr
         integer(c_int), value, intent(in) :: ncerr
      end function nc_strerror
   end interface

contains

   !----- Build a NUL-terminated C string from a Fortran string. --------------------------!
   pure function cstr(f) result(c)
      character(len=*), intent(in)  :: f
      character(kind=c_char) :: c(len_trim(f) + 1)
      integer(ik) :: i, n
      n = len_trim(f)
      do i = 1_ik, n
         c(i) = f(i:i)
      end do
      c(n + 1) = c_null_char
   end function cstr

   !---------------------------------------------------------------------------------------!
   ! String-argument wrappers, for plain character(len=*) arguments; the raw bind(c) forms !
   ! stay for the numeric calls. Each binds cstr() to a named local before the C call:     !
   ! no function result goes straight into a call (CLAUDE.md, compiler traps).             !
   !---------------------------------------------------------------------------------------!
   integer(c_int) function nc_create_f(path, cmode, ncidp) result(st)
      character(len=*), intent(in)  :: path
      integer(c_int),   intent(in)  :: cmode
      integer(c_int),   intent(out) :: ncidp
      character(kind=c_char) :: cpath(len_trim(path) + 1)
      cpath = cstr(path)
      st = nc_create(cpath, cmode, ncidp)
   end function nc_create_f

   integer(c_int) function nc_def_dim_f(ncid, name, length, idp) result(st)
      integer(c_int),    intent(in)  :: ncid
      character(len=*),  intent(in)  :: name
      integer(c_size_t), intent(in)  :: length
      integer(c_int),    intent(out) :: idp
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      st = nc_def_dim(ncid, cname, length, idp)
   end function nc_def_dim_f

   integer(c_int) function nc_def_var_f(ncid, name, xtype, ndims, dimids, varidp) result(st)
      integer(c_int),   intent(in)  :: ncid, xtype, ndims
      character(len=*), intent(in)  :: name
      integer(c_int),   intent(in)  :: dimids(*)
      integer(c_int),   intent(out) :: varidp
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      st = nc_def_var(ncid, cname, xtype, ndims, dimids, varidp)
   end function nc_def_var_f

   integer(c_int) function nc_put_att_text_f(ncid, varid, name, length, tp) result(st)
      integer(c_int),    intent(in) :: ncid, varid
      character(len=*),  intent(in) :: name, tp
      integer(c_size_t), intent(in) :: length
      character(kind=c_char) :: cname(len_trim(name) + 1)
      character(kind=c_char) :: ctp(len_trim(tp) + 1)
      cname = cstr(name)
      ctp = cstr(tp)
      st = nc_put_att_text(ncid, varid, cname, length, ctp)
   end function nc_put_att_text_f

   !----- A one-value numeric attribute (e.g. _FillValue), through named locals as above. --!
   integer(c_int) function nc_put_att_double_f(ncid, varid, name, xtype, val) result(st)
      integer(c_int),   intent(in) :: ncid, varid, xtype
      character(len=*), intent(in) :: name
      real(c_double),   intent(in) :: val
      character(kind=c_char) :: cname(len_trim(name) + 1)
      real(c_double)         :: opv(1)
      cname = cstr(name)
      opv(1) = val
      st = nc_put_att_double(ncid, varid, cname, xtype, 1_c_size_t, opv)
   end function nc_put_att_double_f

   integer(c_int) function nc_put_att_int_f(ncid, varid, name, xtype, val) result(st)
      integer(c_int),   intent(in) :: ncid, varid, xtype
      character(len=*), intent(in) :: name
      integer(c_int),   intent(in) :: val
      character(kind=c_char) :: cname(len_trim(name) + 1)
      integer(c_int)         :: opv(1)
      cname = cstr(name)
      opv(1) = val
      st = nc_put_att_int(ncid, varid, cname, xtype, 1_c_size_t, opv)
   end function nc_put_att_int_f

   integer(c_int) function nc_open_f(path, omode, ncidp) result(st)
      character(len=*), intent(in)  :: path
      integer(c_int),   intent(in)  :: omode
      integer(c_int),   intent(out) :: ncidp
      character(kind=c_char) :: cpath(len_trim(path) + 1)
      cpath = cstr(path)
      st = nc_open(cpath, omode, ncidp)
   end function nc_open_f

   integer(c_int) function nc_inq_varid_f(ncid, name, varidp) result(st)
      integer(c_int),   intent(in)  :: ncid
      character(len=*), intent(in)  :: name
      integer(c_int),   intent(out) :: varidp
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      st = nc_inq_varid(ncid, cname, varidp)
   end function nc_inq_varid_f

   integer(c_int) function nc_inq_dimlen_f(ncid, name, length) result(st)
      integer(c_int),    intent(in)  :: ncid
      character(len=*),  intent(in)  :: name
      integer(c_size_t), intent(out) :: length
      integer(c_int) :: dimid
      character(kind=c_char) :: cname(len_trim(name) + 1)
      length = 0_c_size_t
      cname  = cstr(name)
      st = nc_inq_dimid(ncid, cname, dimid)
      if (st /= NC_NOERR) return
      st = nc_inq_dimlen(ncid, dimid, length)
   end function nc_inq_dimlen_f

   !----- Read a text attribute into a Fortran string (netCDF text attrs are not NUL-       !
   !      terminated, so inquire the length first, then copy exactly that many chars).      !
   integer(c_int) function nc_get_att_text_f(ncid, varid, name, text) result(st)
      integer(c_int),   intent(in)  :: ncid, varid
      character(len=*), intent(in)  :: name
      character(len=*), intent(out) :: text
      character(kind=c_char) :: cname(len_trim(name) + 1)
      character(kind=c_char), allocatable :: cbuf(:)
      integer(c_size_t) :: alen
      integer(ik) :: i, n
      text  = ''
      cname = cstr(name)
      st = nc_inq_attlen(ncid, varid, cname, alen)
      if (st /= NC_NOERR) return
      allocate(cbuf(max(int(alen, ik), 1_ik)))
      st = nc_get_att_text(ncid, varid, cname, cbuf)
      if (st /= NC_NOERR) return
      n = min(int(alen, ik), int(len(text), ik))
      do i = 1_ik, n
         text(i:i) = cbuf(i)
      end do
   end function nc_get_att_text_f

   !----- A one-value numeric attribute, converted to double by the library whatever its stored  !
   !      type. A missing or longer attribute returns the error status and leaves `value` at 0.    !
   integer(c_int) function nc_get_att_double_f(ncid, varid, name, value) result(st)
      integer(c_int),   intent(in)  :: ncid, varid
      character(len=*), intent(in)  :: name
      real(c_double),   intent(out) :: value
      character(kind=c_char) :: cname(len_trim(name) + 1)
      integer(c_size_t) :: alen
      real(c_double)    :: buf(1)
      value = 0.0_c_double
      cname = cstr(name)
      st = nc_inq_attlen(ncid, varid, cname, alen)
      if (st /= NC_NOERR) return
      if (alen /= 1_c_size_t) then
         st = -1_c_int
         return
      end if
      st = nc_get_att_double(ncid, varid, cname, buf)
      if (st == NC_NOERR) value = buf(1)
   end function nc_get_att_double_f

   !----- A variable's gaps and packing, as the netCDF attribute conventions declare them: a value !
   !      equal to _FillValue or missing_value is a gap, and scale_factor and add_offset unpack the !
   !      rest. Read once per variable here; nc_unpacked applies it value by value. -------------!
   function nc_packing_of(ncid, varid) result(packing)
      integer(c_int), intent(in) :: ncid, varid
      type(nc_packing_t)         :: packing
      packing%fill    = attribute_or_nan('_FillValue')
      packing%missing = attribute_or_nan('missing_value')
      packing%scale   = attribute_or_nan('scale_factor')
      packing%offset  = attribute_or_nan('add_offset')
   contains
      real(c_double) function attribute_or_nan(name) result(x)
         use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
         character(len=*), intent(in) :: name
         real(c_double) :: a
         x = ieee_value(1.0_c_double, ieee_quiet_nan)
         if (nc_get_att_double_f(ncid, varid, name, a) == NC_NOERR) x = a
      end function attribute_or_nan
   end function nc_packing_of

   !----- One stored value as data: NaN for a gap, otherwise unpacked. An attribute the variable  !
   !      lacks (NaN in the record) does nothing, so a plain variable reads as stored. ------------!
   elemental real(c_double) function nc_unpacked(packing, raw) result(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
      type(nc_packing_t), intent(in) :: packing
      real(c_double),     intent(in) :: raw
      x = raw
      if (.not. ieee_is_nan(packing%fill)) then
         if (raw == packing%fill) x = ieee_value(x, ieee_quiet_nan)
      end if
      if (.not. ieee_is_nan(packing%missing)) then
         if (raw == packing%missing) x = ieee_value(x, ieee_quiet_nan)
      end if
      if (.not. ieee_is_nan(packing%scale))  x = x * packing%scale
      if (.not. ieee_is_nan(packing%offset)) x = x + packing%offset
   end function nc_unpacked

   !----- The netCDF library's text for a status ('' when it has none). ---------------------!
   function nc_error_text(status) result(text)
      integer(c_int), intent(in) :: status
      character(len=256)         :: text
      character(kind=c_char), pointer :: cbuf(:)
      type(c_ptr) :: cp
      integer     :: i
      text = ''
      cp = nc_strerror(status)
      if (.not. c_associated(cp)) return
      call c_f_pointer(cp, cbuf, [256])
      do i = 1, 256
         if (cbuf(i) == c_null_char) exit
         text(i:i) = cbuf(i)
      end do
   end function nc_error_text

   !----- Abort with the netCDF error message if `status` is not NC_NOERR. ----------------!
   subroutine nc_check(status, context)
      integer(c_int),   intent(in) :: status
      character(len=*), intent(in) :: context
      character(len=256) :: text
      if (status == NC_NOERR) return
      write(*,'(3a,i0)') 'meds_netcdf_c: netCDF error in ', trim(context), ', status=', status
      text = nc_error_text(status)
      if (len_trim(text) > 0) write(*,'(2a)') '   ', trim(text)
      error stop 'meds_netcdf_c: netCDF call failed'
   end subroutine nc_check

end module meds_netcdf_c
