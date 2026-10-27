! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_setup -- a run's output file set, laid out once: open_output_files builds
! the registry from the config, applies [output].io_config's per-variable overrides, lays
! out the polygon axis and sizes the records; attach_output gives each polygon its
! buffers, switches on the diagnostic blocks the set's variables read
! (activate_polygon_diag) and fills its row of the per-polygon variables. netCDF-free
! (docs/science/diagnostics.md §5).
!---~---
module ed_output_setup
   use ed_kinds,            only : wp, ik
   use ed_config,           only : ed_config_t
   use ed_forcing_params,   only : MET_BACKEND_ERA5LAND
   use ed_polygon_state_vars,   only : polygon_t
   use ed_diag_vars,            only : N_CO_FAST_ROWS, N_CO_ROWS, N_PA_ROWS, N_PY_ROWS, N_SI_ROWS,               &
                                       alloc_cohort_diag, alloc_patch_diag, alloc_site_diag
   use ed_output_params,    only : N_FREQ, TIER_FAST, N_DBH_CLASS_DEFAULT, DBH_EDGES_DEFAULT, FREQ_NONE
   use ed_output_types,     only : output_files_t, output_buffers_t, MAX_DBH_CLASS, DIM_COHORT, DIM_PATCH,        &
                                   DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE, MISSING_VALUE,          &
                                   MISSING_INT, DIM_SCALAR
   use ed_constants,         only : n_soil_layer_max
   use ed_output_registry,    only : build_output_registry, freq_bit, parse_stream_mask, selection_t,           &
                                     apply_selections, resolve_streams, OVR_TRUE, OVR_FALSE, OVR_MASK
   use ed_output_integrate,   only : alloc_integ_buffer
   use ed_output_extract,     only : FLD_C_DIAG0, FLD_P_DIAG0, FLD_PY_DIAG0, FLD_SI_DIAG0, SRC_F_CO0, SRC_NONE
   use ed_allometry_lib,      only : HEIGHT_CHAVE
   use ed_soil_params,        only : SOIL_BC_AQUIFER
   use ed_toml,               only : toml_table_t, toml_parse_file
   implicit none
   private

   public :: open_output_files, attach_output
   public :: setup_output_files, set_output_polygons, finalize_output_files, alloc_output_buffers
   public :: alloc_output_set, activate_polygon_diag

contains

   !---~---
   !   open_output_files -- lay out one output file set: the registry from the config, the
   ! per-variable overrides of [output].io_config, and the polygon axis with each
   ! polygon's cell and site count. Each polygon that writes into the set then calls
   ! attach_output, which gives the files that polygon's own sites.
   !---~---
   subroutine open_output_files(cfg, files, verbose, polygon_id, n_site, lat, lon, row, col)
      type(ed_config_t),  intent(in)      :: cfg
      type(output_files_t), intent(inout) :: files
      logical,              intent(in)    :: verbose
      integer(ik),          intent(in)    :: polygon_id(:)   !--- the polygons' ids ...
      integer(ik),          intent(in)    :: n_site(:)       !--- ... their site counts ...
      real(wp),             intent(in)    :: lat(:), lon(:)  !--- ... their cells' centres [deg] ...
      integer(ik),          intent(in)    :: row(:), col(:)  !--- ... and grid indices (0-based)
      call setup_output_files(files, cfg)
      if (len_trim(cfg%output%io_config) > 0)                                                     &
         call apply_io_overrides(files, trim(cfg%output%io_config), verbose)
      call set_output_polygons(files, polygon_id, n_site, lat, lon, row, col,                      &
                               any(cfg%pft%allom(1:cfg%pft%n)%height_allometry == HEIGHT_CHAVE))
      call finalize_output_files(files)
   end subroutine open_output_files

   !---~---
   !   A polygon's share of an output file set: its buffers, which read the soil diagnostics
   ! when its fast loop runs, the diagnostic blocks it must accumulate for the set's
   ! variables, and its row of the per-polygon variables: each site's area, terrain, layers
   ! and soil column, and the height allometry's site index its cell's map gave it.
   !---~---
   subroutine attach_output(cfg, files, bufs, poly, p)
      type(ed_config_t),    intent(in)      :: cfg
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(inout) :: poly
      integer(ik),            intent(in)    :: p      !--- the polygon's row of the file set
      integer(ik) :: is
      call alloc_output_buffers(files, bufs)
      bufs%diag%soil_ready = cfg%fast_biophysics_on
      call activate_polygon_diag(files, poly)
      do is = 1_ik, poly%site%n
         associate (surf => poly%site%surface(is))
            files%site_area(is, p)         = poly%site%area(is)
            files%site_elevation(is, p)    = surf%elevation_m
            files%site_slope(is, p)        = surf%slope_deg
            files%site_aspect(is, p)       = surf%aspect_deg
            files%site_n_soil(is, p)       = surf%soil_profile%n_active
            files%site_soil_z(:, is, p)    = surf%soil_profile%z_node
            files%site_soil_face(:, is, p) = surf%soil_profile%soil_layer_z
            files%site_soil_depth(is, p)   = surf%soil%depth
            files%site_bottom_bc(is, p)    = surf%soil%water_bottom_bc
            files%site_drain_factor(is, p) = merge(MISSING_VALUE, surf%soil%drain_factor,           &
                                                   surf%soil%water_bottom_bc == SOIL_BC_AQUIFER)
            files%site_deep_temp(is, p)    = surf%soil%deep_temp
         end associate
      end do
      !--- E is the polygon's, the same in every PFT's record.
      if (allocated(files%polygon_chave_e)) files%polygon_chave_e(p) = poly%cohort%pft_table%allom(1)%chave_e
   end subroutine attach_output

   !---~---
   !   Apply [output].io_config's [variables] lines to the file set's registry (the grammar and
   ! the rules: docs/science/diagnostics.md §5). A key is a variable's name or a pattern with
   ! `*`, quoted or not; a value true, false or a quoted "F D M Y" string.
   !---~---
   subroutine apply_io_overrides(files, tomlpath, verbose)
      type(output_files_t),   intent(inout) :: files
      character(len=*),       intent(in)    :: tomlpath
      logical,                intent(in)    :: verbose
      type(toml_table_t) :: tt
      type(selection_t), allocatable :: sel(:)
      logical            :: ok
      integer(ik)        :: i, n, mask, status
      character(len=256) :: raw, sval
      character(len=64)  :: key
      character(len=8)   :: bad
      call toml_parse_file(tomlpath, tt, ok)
      if (.not. ok) error stop 'ed_output_setup: cannot read [output].io_config = '//trim(tomlpath)
      allocate(sel(max(tt%n, 1_ik)))
      n = 0_ik
      do i = 1_ik, tt%n
         if (len_trim(tt%key(i)) <= 10) cycle
         if (tt%key(i)(1:10) /= 'variables.') cycle
         key = adjustl(tt%key(i)(11:))
         if (key(1:1) == '"') key = key(2:index(key(2:), '"'))     ! a quoted key: "mean_*_py"
         raw = adjustl(tt%val(i))
         n = n + 1_ik
         sel(n)%key = key
         select case (trim(raw))
         case ('true')
            sel(n)%kind = OVR_TRUE
         case ('false')
            sel(n)%kind = OVR_FALSE
         case default
            if (raw(1:1) /= '"')                                                                  &
               error stop 'ed_output_setup: io_config bad value for '//trim(key)//                     &
                          ' (expected true|false or a quoted "F D M Y" string)'
            sval = raw(2:index(raw(2:), '"'))
            call parse_stream_mask(trim(sval), mask, status, bad)
            if (status /= 0_ik)                                                                   &
               error stop 'ed_output_setup: io_config unknown stream token "'//trim(bad)//              &
                          '" for '//trim(key)
            sel(n)%kind = OVR_MASK
            sel(n)%mask = mask
         end select
      end do
      call apply_selections(files%reg, sel(1:n))
      call resolve_streams(files%reg)
      if (verbose) write(*,'(3a)') ' output: applied the selection in ', trim(tomlpath), ''
   end subroutine apply_io_overrides

   !---~---
   !   Build a file set's registry and settings, with no buffers yet, so a caller can apply the
   ! ed_io_config.toml overrides (which change the buffers needed) before
   ! finalize_output_files and alloc_output_buffers.
   !---~---
   subroutine setup_output_files(files, cfg)
      type(output_files_t),   intent(out) :: files
      type(ed_config_t),    intent(in)    :: cfg
      integer(ik) :: ne
      files%enabled = cfg%output%enabled
      files%fast_loop_on = cfg%fast_biophysics_on
      call build_output_registry(files%reg, cfg)
      files%cohort_max = cfg%output%cohort_max
      files%patch_max  = cfg%output%patch_max
      files%dir        = cfg%output%dir
      files%prefix     = cfg%output%prefix
      files%file_chunk = cfg%output%file_chunk
      files%sync_every = cfg%output%sync_every
      files%fast_interval_steps = cfg%output%fast_interval_steps
      if (cfg%forcing%forcing_on .and. cfg%forcing%format == MET_BACKEND_ERA5LAND)                   &
         files%forcing_qair = 'computed by ED-MORE from the ED_ERA5land 2 m dewpoint Td and surface '//         &
                            'pressure P: q = 0.62196 e / (P - 0.37804 e), e = e_sat(Td), the Bolton (1980) '// &
                            'liquid-water saturation vapour pressure of ed_therm_lib'

      !---~---
      !   The derived diagnostics' parameters (diag_setup_t). n_pft is the run-time PFT count,
      ! so the files carry a pft coordinate (docs/science/diagnostics.md §3).
      !---~---
      files%diag%n_pft  = cfg%pft%n
      !--- DBH size classes: TOML edges if given, else the built-in inventory set.
      if (cfg%output%n_dbh_class > 0_ik) then
         files%diag%n_dbh_class = min(cfg%output%n_dbh_class, MAX_DBH_CLASS)
         ne = files%diag%n_dbh_class + 1_ik
         files%diag%dbh_edges(1:ne) = cfg%output%dbh_class_edges(1:ne)
      else
         files%diag%n_dbh_class = N_DBH_CLASS_DEFAULT
         files%diag%dbh_edges(1:N_DBH_CLASS_DEFAULT + 1_ik) = DBH_EDGES_DEFAULT
      end if
      call check_dbh_edges(files%diag%dbh_edges, files%diag%n_dbh_class)
      !---~---
      !   attach_output sets soil_ready on a polygon's buffers when its fast loop runs; until
      ! then the psi/wetness diagnostics emit _FillValue.
      !---~---
      files%diag%soil_ready = .false.
   end subroutine setup_output_files

   !--- Largest slab length any live variable can produce (see finalize_output_files).
   pure integer(ik) function live_max_slab(files) result(cap)
      type(output_files_t),   intent(in) :: files
      integer(ik) :: k
      cap = 1_ik
      do k = 1_ik, files%reg%nvar
         if (.not. files%reg%var(k)%enabled)      cycle
         if (files%reg%var(k)%streams == FREQ_NONE) cycle
         cap = max(cap, dim_capacity(files, files%reg%var(k)%dim))
      end do
   end function live_max_slab

   !--- Slab capacity of one axis.
   pure integer(ik) function dim_capacity(files, dm) result(cap)
      type(output_files_t),   intent(in) :: files
      integer(ik),            intent(in) :: dm
      select case (dm)
      case (DIM_COHORT)
         cap = files%cohort_max
      case (DIM_PATCH)
         cap = files%patch_max
      case (DIM_SOIL)
         cap = files%site_max * n_soil_layer_max
      case (DIM_PFT)
         cap = max(files%diag%n_pft, 1_ik)
      case (DIM_SIZE)
         cap = max(files%diag%n_dbh_class, 1_ik)
      case (DIM_SOIL_PATCH)
         cap = files%patch_max * n_soil_layer_max
      case (DIM_SITE)
         cap = files%site_max
      case default
         cap = 0_ik
      end select
   end function dim_capacity

   !---~---
   !   The size-class edges must be strictly ascending, or dbh_class_index mis-bins without a
   ! word and the closure of docs/science/diagnostics.md eq. (3) fails; stop at start-up.
   !---~---
   subroutine check_dbh_edges(edges, n_class)
      real(wp),    intent(in) :: edges(:)
      integer(ik), intent(in) :: n_class
      integer(ik) :: j
      do j = 1_ik, n_class
         if (.not. (edges(j+1_ik) > edges(j)))                                                   &
            error stop 'ed_output_setup: [output].dbh_class_edges must be strictly ascending'
      end do
   end subroutine check_dbh_edges

   !---~---
   !   Finish the file set once the registry is final (after setup_output_files and any
   ! per-variable overrides): the per-tier stream handles and the slab size.
   !---~---
   subroutine finalize_output_files(files)
      type(output_files_t),  intent(inout) :: files
      integer(ik) :: t
      do t = 1_ik, N_FREQ
         files%stream(t)%freq  = freq_bit(t)
         allocate(files%stream(t)%vid(files%reg%nvar))
         files%stream(t)%vid = -1_ik
      end do
      !---~---
      !   Size the records' slab to the largest live axis, not the largest that exists, so a
      ! run pays no memory for an axis it switched off. It is sized here, after every
      ! override, since an override can switch on a longer slab variable.
      !---~---
      files%max_slab = live_max_slab(files)
      !---~---
      !   Whether the cohort and patch caps bind: only a live variable on that axis fills a slab of
      ! that size.
      !---~---
      files%cohort_axis = live_axis(files, [DIM_COHORT])
      files%patch_axis  = live_axis(files, [DIM_PATCH, DIM_SOIL_PATCH])
   end subroutine finalize_output_files

   !--- Does any live variable have one of these axes?
   pure logical function live_axis(files, dims) result(yes)
      type(output_files_t), intent(in) :: files
      integer(ik),          intent(in) :: dims(:)
      integer(ik) :: k
      yes = .false.
      do k = 1_ik, files%reg%nvar
         if (.not. files%reg%var(k)%enabled)        cycle
         if (files%reg%var(k)%streams == FREQ_NONE) cycle
         if (any(dims == files%reg%var(k)%dim)) then
            yes = .true.
            return
         end if
      end do
   end function live_axis

   !---~---
   !   Allocate one polygon's buffers for a finalized file set: the integrator buffers of
   ! every live (variable, tier) pair and the per-tier scratch records.
   !---~---
   subroutine alloc_output_buffers(files, bufs)
      type(output_files_t),  intent(in)  :: files
      type(output_buffers_t), intent(out) :: bufs
      integer(ik) :: t, k, nv, cap, bit
      nv = files%reg%nvar
      allocate(bufs%buf(nv, N_FREQ))
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         do k = 1_ik, nv
            if (files%reg%var(k)%enabled .and. iand(files%reg%var(k)%streams, bit) /= 0_ik) then
               cap = dim_capacity(files, files%reg%var(k)%dim)
               call alloc_integ_buffer(bufs%buf(k,t), files%reg%var(k), bit, cap)
               bufs%buf(k,t)%var_id = k
            end if
         end do
      end do
      !---~---
      !   Each tier's scratch record: every variable's scalar, and a slab column for each of the
      ! tier's live slab variables only, so a variable nobody writes costs nothing here.
      !---~---
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         allocate(bufs%pending(t)%sval(nv), bufs%pending(t)%svalid(nv), bufs%pending(t)%nslab(nv))
         allocate(bufs%pending(t)%col(nv), source=0_ik)
         cap = 0_ik
         do k = 1_ik, nv
            if (files%reg%var(k)%enabled .and. iand(files%reg%var(k)%streams, bit) /= 0_ik .and.     &
                files%reg%var(k)%dim /= DIM_SCALAR) then
               cap = cap + 1_ik
               bufs%pending(t)%col(k) = cap
            end if
         end do
         allocate(bufs%pending(t)%slab(files%max_slab, max(cap, 1_ik)),                              &
                  bufs%pending(t)%slabvalid(files%max_slab, max(cap, 1_ik)))
         bufs%pending(t)%used = .false.
      end do
      !---~---
      !   The fast loop stages its sub-steps for the F tier and for the snapshots the coarse
      ! tiers take from the last one.
      !---~---
      bufs%fast_on = files%reg%nidx(TIER_FAST) > 0_ik
      do t = TIER_FAST + 1_ik, N_FREQ
         do k = 1_ik, files%reg%nidx(t)
            if (files%reg%var(files%reg%idx_freq(k, t))%snap_staged) bufs%fast_on = .true.
         end do
      end do
      bufs%fast_on = bufs%fast_on .and. files%enabled
      bufs%diag = files%diag
   end subroutine alloc_output_buffers

   !---~---
   !   Convenience: registry, stream handles and one polygon's buffers in one call, with an
   ! optional selection in place of an io_config file.
   !---~---
   subroutine alloc_output_set(files, bufs, cfg, sel)
      type(output_files_t),  intent(out) :: files
      type(output_buffers_t), intent(out) :: bufs
      type(ed_config_t),   intent(in)  :: cfg
      type(selection_t), optional, intent(in) :: sel(:)   !--- io_config lines, as apply_io_overrides reads them
      call setup_output_files(files, cfg)
      if (present(sel)) then
         call apply_selections(files%reg, sel)
         call resolve_streams(files%reg)
      end if
      call set_output_polygons(files, [0_ik], [1_ik], [0.0_wp], [0.0_wp], [0_ik], [0_ik], .false.)
      call finalize_output_files(files)
      call alloc_output_buffers(files, bufs)
   end subroutine alloc_output_set

   !---~---
   !   The file set's polygon and site axes: one row per polygon, its id, location, grid indices
   ! and site count; the site axis is the largest polygon's. The rest of a row -- its sites'
   ! areas, terrain and soil columns -- is filled when the polygon attaches to the set
   ! (attach_output); a polygon's missing sites keep the fill value.
   !---~---
   subroutine set_output_polygons(files, polygon_id, n_site, lat, lon, row, col, with_chave_e)
      type(output_files_t), intent(inout) :: files
      integer(ik),          intent(in)    :: polygon_id(:), n_site(:), row(:), col(:)
      real(wp),             intent(in)    :: lat(:), lon(:)
      logical,              intent(in)    :: with_chave_e   !--- the Chave height allometry's index E, per polygon
      integer(ik) :: n, ns
      n  = size(polygon_id, kind=ik)
      ns = max(maxval(n_site), 1_ik)
      files%n_polygon      = n
      files%site_max       = ns
      files%polygon_id     = polygon_id
      files%polygon_n_site = n_site
      files%polygon_row    = row
      files%polygon_col    = col
      files%polygon_lat    = lat
      files%polygon_lon    = lon
      allocate(files%site_area(ns, n), files%site_elevation(ns, n), files%site_slope(ns, n),        &
               files%site_aspect(ns, n), files%site_soil_depth(ns, n), files%site_drain_factor(ns, n), &
               files%site_deep_temp(ns, n), source=MISSING_VALUE)
      allocate(files%site_n_soil(ns, n), source=MISSING_INT)
      allocate(files%site_bottom_bc(ns, n), source=MISSING_INT)
      allocate(files%site_soil_z(n_soil_layer_max, ns, n), source=0.0_wp)
      allocate(files%site_soil_face(n_soil_layer_max + 1_ik, ns, n), source=0.0_wp)
      if (with_chave_e) allocate(files%polygon_chave_e(n), source=0.0_wp)
   end subroutine set_output_polygons

   !---~---
   !   Switch each diagnostic block (cohort, slow cohort, patch, polygon, site) on iff a live
   ! variable of the final registry reads it. A block that is off stays unallocated and its
   ! ed_diag_vars entry points do nothing, so a run pays nothing for what it does not report.
   !---~---
   subroutine activate_polygon_diag(files, poly)
      type(output_files_t),   intent(in)    :: files
      type(polygon_t),        intent(inout) :: poly
      logical     :: need_c, need_p, need_s, need_y, need_si
      integer(ik) :: t, j, k, src
      need_c = .false.
      need_p = .false.
      need_s = .false.
      need_y = .false.
      need_si = .false.
      if (files%enabled) then
         do t = 1_ik, N_FREQ
            do j = 1_ik, files%reg%nidx(t)
               k   = files%reg%idx_freq(j, t)
               !---~---
               !   The F tier, and a snapshot that reads the staging, need the block that stages
               ! their source, not the one the coarse means accumulate in.
               !---~---
               src = files%reg%var(k)%source_id
               if (t == TIER_FAST .or. (files%reg%var(k)%snap_staged .and.                           &
                                        files%reg%var(k)%source_fast /= SRC_NONE))                   &
                  src = files%reg%var(k)%source_fast
               if (src > FLD_C_DIAG0 .and. src <= FLD_C_DIAG0 + N_CO_FAST_ROWS) need_c = .true.
               if (src > SRC_F_CO0 .and. src <= SRC_F_CO0 + N_CO_FAST_ROWS) need_c = .true.
               if (src > FLD_P_DIAG0 .and. src <= FLD_P_DIAG0 + N_PA_ROWS) need_p = .true.
               if (src > FLD_C_DIAG0 + N_CO_FAST_ROWS .and. src <= FLD_C_DIAG0 + N_CO_ROWS) need_s = .true.
               if (src > FLD_PY_DIAG0 .and. src <= FLD_PY_DIAG0 + N_PY_ROWS) need_y = .true.
               if (src > FLD_SI_DIAG0 .and. src <= FLD_SI_DIAG0 + N_SI_ROWS) need_si = .true.
            end do
         end do
      end if
      !---~---
      !   Only the fast loop fills the per-cohort and polygon blocks, and the fast rows of the
      ! patch block: in a slow-only run they stay off, and their variables read as missing.
      !---~---
      need_c = need_c .and. files%fast_loop_on
      need_y = need_y .and. files%fast_loop_on
      need_si = (need_si .and. files%fast_loop_on) .or. poly%site%diag%active
      poly%patch%diag%fast_rows = files%fast_loop_on
      !--- Additive: a block already on stays on.
      need_c = need_c .or. poly%cohort%diag%active
      need_s = need_s .or. poly%cohort%sdiag%active
      need_p = need_p .or. poly%patch%diag%active
      need_y = need_y .or. poly%diag%active
      call alloc_cohort_diag(poly%cohort%diag,  max(poly%cohort%cap, 1_ik), need_c)
      call alloc_cohort_diag(poly%cohort%sdiag, max(poly%cohort%cap, 1_ik), need_s, slow_rows=.true.)
      call alloc_patch_diag (poly%patch%diag,   max(poly%patch%cap,  1_ik), need_p)
      poly%diag%active = need_y                    ! the polygon block: fixed size, nothing to allocate
      call alloc_site_diag(poly%site%diag, poly%site%n, need_si)
   end subroutine activate_polygon_diag

end module ed_output_setup
