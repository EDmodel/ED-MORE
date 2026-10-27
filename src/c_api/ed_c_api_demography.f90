! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_c_api_demography -- the ISO_C_BINDING shim that exposes the DEMOGRAPHIC model to C /
! Python (ctypes). It goes into the optional shared library libedmore.so (ED_BUILD_PYLIB) and
! into a mandatory ctest target, so a signature change is a build failure in a default build.
!
!   `polygon_t` (allocatable components, so not f2py-able) crosses as an OPAQUE HANDLE: a
! module-`save` registry of live sites indexed by a small integer the caller holds. Config is
! likewise an opaque handle loaded from the same TOML `edmore` uses. Only flat arrays
! (copied into caller buffers) and scalars cross the boundary; getters always COPY.
!
!   WHY THIS IS ONE FILE, not a `site` / `demography` pair: that registry (g_poly / g_cfg /
! site_used / cfg_used / g_generation) is module state
! EVERY entry point indexes. Splitting lifecycle from verbs would have to hoist it into a third
! module both `use`, which buys nothing at this size and gives the registry two places to be
! got wrong.
!
!   It drives the slow model from Python: load config -> create site -> init (bare ground or
! census) -> advance_slow (the driver's carbon orchestration) or apply_rates (caller-supplied
! demographic rates) -> read state -> free. A per-handle GENERATION counter (bumped by every
! init and advance; fuse/fiss reorder the SoA) lets the caller detect a stale positional
! snapshot -- global_id is the only stable key.
!---~---
module ed_c_api_demography
   use iso_c_binding
   use ed_kinds,                    only : wp, ik
   use ed_constants,                only : tiny_num
   use ed_config,                   only : ed_config_t, growth_window_steps
   use ed_config_io,                only : load_ed_config
   use ed_cohort_state_vars,         only : cohort_deriv_alloc
   use ed_polygon_state_vars,        only : polygon_t, free_polygon
   use ed_init_bare_ground, only : init_bare_ground
   use ed_init_census, only : init_from_census, restructure_census_stand
   use ed_slow_dynamics,   only : vegetation_dynamics, restructure_stand
   use ed_slow_ledger,     only : slow_ledger_t
   use ed_slow_state,      only : clear_litter
   use ed_slow_update,   only : accumulate_recruit_pool, update_cohort_states, fill_cohort_deriv,        &
                                age_patches, update_overtopping_lai, sort_cohorts
   use ed_output_reduce,   only : total_agb, total_lai, total_nplant, total_basal_area
   use ed_allometry_lib,            only : dbh_to_height, size_to_leaf_carbon, size_to_wood_carbon,  &
                                           storage_carbon_target
   implicit none
   private

   !---~---
   !   The shim's entry points. PUBLIC so a ctest target can `use` this module and call them
   ! as Fortran; the bind(c) names are exported at link level either way, but a private
   ! procedure is unreachable from a test, and a shim no test compiles can stop compiling
   ! without anyone noticing (see the header).
   !---~---
   public :: ed_config_load, ed_config_n_pft, ed_config_dt_years
   public :: ed_site_create, ed_site_init_bare, ed_site_init_census, ed_site_free
   public :: ed_advance_slow, ed_apply_rates
   public :: ed_site_generation, ed_site_n_patch, ed_site_n_cohort
   public :: ed_site_total_agb, ed_site_total_lai, ed_site_total_nplant
   public :: ed_site_total_basal_area, ed_site_get_real, ed_site_get_int
   public :: ed_site_get_patch_real, ed_site_get_patch_int

   integer, parameter :: MAXH = 8
   type(polygon_t),     target, save :: g_poly(MAXH)
   type(ed_config_t),   target, save :: g_cfg(MAXH)
   logical,             save :: site_used(MAXH) = .false.
   logical,             save :: cfg_used(MAXH)  = .false.
   integer(c_long),     save :: g_generation(MAXH) = 0_c_long

contains

   !--- Build a Fortran string from an explicit-length C char buffer (no strlen needed).
   pure function f_string(cbuf, n) result(s)
      character(kind=c_char), intent(in) :: cbuf(*)
      integer(c_int),         intent(in) :: n
      character(len=n) :: s
      integer :: i
      do i = 1, n
         s(i:i) = cbuf(i)
      end do
   end function f_string

   !--- Load a config from a TOML path -> opaque cfg handle (>=1), or -1 if the registry is full.
   function ed_config_load(path, path_len) result(h) bind(c, name="ed_config_load")
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value, intent(in)  :: path_len
      integer(c_int)                     :: h
      integer :: i
      h = -1_c_int
      do i = 1, MAXH
         if (.not. cfg_used(i)) then
            h = int(i, c_int)
            exit
         end if
      end do
      if (h < 0_c_int) return
      call load_ed_config(f_string(path, path_len), g_cfg(h))
      cfg_used(h) = .true.
   end function ed_config_load

   !--- Create an empty site slot -> opaque site handle (>=1), or -1 if full.
   function ed_site_create() result(h) bind(c, name="ed_site_create")
      integer(c_int) :: h
      integer :: i
      h = -1_c_int
      do i = 1, MAXH
         if (.not. site_used(i)) then
            h = int(i, c_int)
            site_used(i) = .true.
            g_generation(i) = 0_c_long
            exit
         end if
      end do
   end function ed_site_create

   subroutine ed_site_init_bare(sh, ch, n_patch) bind(c, name="ed_site_init_bare")
      integer(c_int), value, intent(in) :: sh, ch, n_patch
      call init_bare_ground(g_poly(sh), g_cfg(ch), int(n_patch, ik))
      g_generation(sh) = g_generation(sh) + 1_c_long
   end subroutine ed_site_init_bare

   !---~---
   !   Start the site from the census file named by the config's [init].census_file, and let
   ! the slow step's own operators restructure it, as edmore does: a census arrives with a
   ! cohort per measured size and a patch per plot cell. Returns 1 when the file was read, 0
   ! when it could not be (the site is then left unbuilt). The overtopping LAI is computed
   ! for the restructured stand, so rates read on the first step see its canopy.
   !---~---
   function ed_site_init_census(sh, ch) result(ok) bind(c, name="ed_site_init_census")
      integer(c_int), value, intent(in) :: sh, ch
      integer(c_int) :: ok
      logical :: found
      call init_from_census(g_poly(sh), g_cfg(ch), trim(g_cfg(ch)%init%census_file), found)
      if (found) then
         call restructure_census_stand(g_poly(sh), g_cfg(ch))
         call update_overtopping_lai(g_poly(sh))
      end if
      ok = merge(1_c_int, 0_c_int, found)
      g_generation(sh) = g_generation(sh) + 1_c_long
   end function ed_site_init_census

   !---~---
   !   Advance one slow step (the driver's CARBON orchestration + cadence). Bumps the
   ! generation (the SoA is reordered by fuse/fission). is_new_month/is_new_year are 0/1:
   ! the step's END crosses that boundary, and the stand is restructured (restructure_stand)
   ! right after the step. The driver restructures at the start of the NEXT step instead,
   ! after its output has read this one, and then resizes the canopy air; with no fast loop
   ! and no output on this path, both orders leave the same stand.
   !---~---
   subroutine ed_advance_slow(sh, ch, is_new_month, is_new_year) bind(c, name="ed_advance_slow")
      integer(c_int), value, intent(in) :: sh, ch, is_new_month, is_new_year
      type(slow_ledger_t) :: no_ledger   ! inactive: this path books no ledger
      !---~---
      !   The litter this step accumulates lands in poly%patch%litter_in and is left there: this
      ! path steps the stand's demography alone and runs no soil step. The accumulator starts
      ! each step empty.
      !---~---
      call clear_litter(g_poly(sh))
      call vegetation_dynamics(g_poly(sh), g_cfg(ch), no_ledger)
      call age_patches(g_poly(sh)%patch, g_cfg(ch)%dt_years)
      call restructure_stand(g_poly(sh), g_cfg(ch), is_new_month /= 0_c_int, is_new_year /= 0_c_int, no_ledger)
      g_generation(sh) = g_generation(sh) + 1_c_long
   end subroutine ed_advance_slow

   !---~---
   !   Apply CALLER-SUPPLIED demographic rates (the empirical / Python path): grow by the
   ! per-cohort growth [cm/yr], die by the per-cohort mortality [1/yr], recruit from the
   ! per-(PFT,patch) recruitment [plant/m2/yr], then sort, age the patches and restructure the
   ! stand with the carbon path's own restructure_stand: monthly recruitment, cohort fusion and
   ! fission; yearly treefall disturbance and patch fusion. `recr` is a flat column-major
   ! (PFT-fastest) npft*npatch array. is_new_month / is_new_year are 0/1.
   !---~---
   subroutine ed_apply_rates(sh, ch, growth, mortality, recr, is_new_month, is_new_year) &
                             bind(c, name="ed_apply_rates")
      integer(c_int), value, intent(in) :: sh, ch, is_new_month, is_new_year
      real(c_double),        intent(in) :: growth(*), mortality(*), recr(*)
      integer(ik) :: n, np, npft, ip, pf, i, n_window
      real(wp), allocatable :: g(:), m(:), rec(:,:)
      real(wp)              :: seed_c
      real(wp) :: dbh_new, height_new
      real(wp) :: lc_new, fc_new, wc_new, nc_new
      type(slow_ledger_t) :: no_ledger   ! inactive: this path books no ledger

      associate (poly => g_poly(sh), cfg => g_cfg(ch))
         n = poly%cohort%n
         np = poly%patch%n
         npft = cfg%pft%n
         allocate(g(max(n,1)), m(max(n,1)), rec(npft, max(np,1)))
         if (n > 0_ik) then
            g(1:n) = real(growth(1:n),    wp)
            m(1:n) = real(mortality(1:n), wp)
         end if
         do ip = 1_ik, np
            do pf = 1_ik, npft
               rec(pf, ip) = real(recr((ip - 1_ik) * npft + pf), wp)
            end do
         end do

         call clear_litter(poly)   ! the accumulator starts each step empty, as on the carbon path
         n_window = growth_window_steps(cfg)
         poly%growth_hist_pos = mod(poly%growth_hist_pos, n_window) + 1_ik
         !---~---
         !   Empirical growth + mortality via the TENDENCY seam: per cohort, forward-derive the
         ! on-allometry pools of the new size from the SUPPLIED growth rate [cm/yr] (mirrors
         ! set_cohort_size), back out each pool's time-
         ! derivative, and advance the moving-average ring buffer with the supplied rate (the
         ! empirical convention -- records the request even when dbh is capped); then let the
         ! core engine's pure applier advance the state. There is deliberately NO empirical
         ! twin of the driver's carbon update_cohort_derivatives -- this is thin C-API glue, not
         ! a named peer computer.
         !---~---
         call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
         associate (cohort => poly%cohort, pft => cfg%pft)
            do i = 1_ik, cohort%n
               pf = cohort%pft(i)
               !---~---
               !   Forward allometry from the supplied dbh growth. It CALLS
               ! ed_allometry_lib rather than inlining it, so a coefficient change
               ! reaches the model and this shim alike.
               !---~---
               dbh_new    = cohort%dbh(i) + g(i) * cfg%dt_years
               height_new = dbh_to_height(cohort%pft_table%allom(pf), dbh_new)
               lc_new     = size_to_leaf_carbon(cohort%pft_table%allom(pf), dbh_new, height_new, cohort%sla(i))
               wc_new     = size_to_wood_carbon(cohort%pft_table%allom(pf), dbh_new, height_new, pft%wood_density(pf))
               fc_new     = pft%root_to_leaf_ratio(pf) * lc_new
               nc_new     = storage_carbon_target(pft%storage_cushion(pf), lc_new, fc_new)
               !---~---
               !   Back out the tendencies + advance the ring buffer with the SUPPLIED growth
               ! rate (empirical behaviour -- records the request even when dbh is capped).
               !---~---
               call fill_cohort_deriv(cohort, i, poly%deriv, cfg%dt_years, m(i), g(i),                &
                                      lc_new, fc_new, wc_new, nc_new, n_window, poly%growth_hist_pos)
            end do
         end associate
         !---~---
         !   Credit the recruit pool with THIS step's supplied rate, then let the monthly
         ! apply_recruitment spawn from it. The carbon driver credits the pool the same way,
         ! every step, which keeps it an exact carbon quantity; this path drives the same
         ! operators, so it does the same. `seed_c` is the ledger's external-import report and
         ! has no consumer here.
         !---~---
         call accumulate_recruit_pool(poly, cfg, rec, cfg%dt_years, seed_c)
         call update_cohort_states(poly%cohort, poly%deriv, cfg%dt_years, cfg%demography%negligible_nplant)
         !---~---
         !   From here on the step is the carbon path's: sort every step (as vegetation_dynamics
         ! does), age the patches, and restructure on the calendar's boundaries with the same
         ! restructure_stand that ed_advance_slow calls, which ends by refreshing the
         ! overtopping LAI the caller's next rates read.
         !---~---
         call sort_cohorts(poly)
         call age_patches(poly%patch, cfg%dt_years)
         call restructure_stand(poly, cfg, is_new_month /= 0_c_int, is_new_year /= 0_c_int, no_ledger)
      end associate
      g_generation(sh) = g_generation(sh) + 1_c_long
   end subroutine ed_apply_rates

   function ed_site_n_patch(sh) result(np) bind(c, name="ed_site_n_patch")
      integer(c_int), value, intent(in) :: sh
      integer(c_int)                    :: np
      np = int(g_poly(sh)%patch%n, c_int)
   end function ed_site_n_patch

   function ed_config_n_pft(ch) result(npft) bind(c, name="ed_config_n_pft")
      integer(c_int), value, intent(in) :: ch
      integer(c_int)                    :: npft
      npft = int(g_cfg(ch)%pft%n, c_int)
   end function ed_config_n_pft

   !--- Read accessors (scalars). All COPY out; none alias the SoA.
   function ed_config_dt_years(ch) result(v) bind(c, name="ed_config_dt_years")
      integer(c_int), value, intent(in) :: ch
      real(c_double)                    :: v
      v = real(g_cfg(ch)%dt_years, c_double)
   end function ed_config_dt_years

   function ed_site_generation(sh) result(g) bind(c, name="ed_site_generation")
      integer(c_int), value, intent(in) :: sh
      integer(c_long)                   :: g
      g = g_generation(sh)
   end function ed_site_generation

   function ed_site_n_cohort(sh) result(n) bind(c, name="ed_site_n_cohort")
      integer(c_int), value, intent(in) :: sh
      integer(c_int)                    :: n
      n = int(g_poly(sh)%cohort%n, c_int)
   end function ed_site_n_cohort

   function ed_site_total_agb(sh) result(v) bind(c, name="ed_site_total_agb")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_agb(g_poly(sh)), c_double)
   end function ed_site_total_agb

   function ed_site_total_lai(sh) result(v) bind(c, name="ed_site_total_lai")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_lai(g_poly(sh)), c_double)
   end function ed_site_total_lai

   function ed_site_total_nplant(sh) result(v) bind(c, name="ed_site_total_nplant")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_nplant(g_poly(sh)), c_double)
   end function ed_site_total_nplant

   function ed_site_total_basal_area(sh) result(v) bind(c, name="ed_site_total_basal_area")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_basal_area(g_poly(sh)), c_double)
   end function ed_site_total_basal_area

   !--- Copy-out per-cohort field getters (caller allocates buf of length n_cohort).
   subroutine ed_site_get_real(sh, field_id, buf) bind(c, name="ed_site_get_real")
      integer(c_int), value, intent(in)  :: sh, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n
      associate (c => g_poly(sh)%cohort)
         n = c%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int)
            buf(1:n) = real(c%dbh(1:n),             c_double)
         case (1_c_int)
            buf(1:n) = real(c%height(1:n),          c_double)
         case (2_c_int)
            buf(1:n) = real(c%nplant(1:n),          c_double)
         case (3_c_int)
            buf(1:n) = real(c%agb(1:n),             c_double)
         case (4_c_int)
            buf(1:n) = real(c%leaf_area(1:n),       c_double)
         case (5_c_int)
            buf(1:n) = real(c%overtopping_lai(1:n), c_double)
         case (6_c_int)
            buf(1:n) = real(c%growth_avg(1:n),      c_double)
         case (7_c_int)
            buf(1:n) = real(c%wood_carbon(1:n),     c_double)
         end select
      end associate
   end subroutine ed_site_get_real

   subroutine ed_site_get_int(sh, field_id, buf) bind(c, name="ed_site_get_int")
      integer(c_int), value, intent(in)  :: sh, field_id
      integer(c_int),        intent(out) :: buf(*)
      integer(ik) :: n
      associate (c => g_poly(sh)%cohort)
         n = c%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int)
            buf(1:n) = int(c%pft(1:n),         c_int)
         case (1_c_int)
            buf(1:n) = int(c%owner_patch(1:n), c_int)
         case (2_c_int)
            buf(1:n) = int(c%global_id(1:n),   c_int)
         end select
      end associate
   end subroutine ed_site_get_int

   !---~---
   !   PATCH-level copy-out, the exact analogue of the two above. The demography
   ! example drives this engine from Python with its own vital-rate laws and writes its own
   ! output; the cohort getters alone cannot describe a stand, because every extensive
   ! quantity is per m2 of ITS OWN patch and needs the patch areas to reach the site, and
   ! the cohort->patch CSR map is what makes the file readable at all.
   !---~---
   subroutine ed_site_get_patch_real(sh, field_id, buf) bind(c, name="ed_site_get_patch_real")
      integer(c_int), value, intent(in)  :: sh, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n
      associate (p => g_poly(sh)%patch)
         n = p%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int)
            buf(1:n) = real(p%area(1:n), c_double)
         case (1_c_int)
            buf(1:n) = real(p%age(1:n),  c_double)
         end select
      end associate
   end subroutine ed_site_get_patch_real

   subroutine ed_site_get_patch_int(sh, field_id, buf) bind(c, name="ed_site_get_patch_int")
      integer(c_int), value, intent(in)  :: sh, field_id
      integer(c_int),        intent(out) :: buf(*)
      integer(ik) :: n
      associate (p => g_poly(sh)%patch)
         n = p%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int)
            buf(1:n) = int(p%dist_type(1:n),     c_int)
         case (1_c_int)
            buf(1:n) = int(p%global_id(1:n),     c_int)
         case (2_c_int)
            buf(1:n) = int(p%cohort_offset(1:n), c_int)
         case (3_c_int)
            buf(1:n) = int(p%cohort_count(1:n),  c_int)
         end select
      end associate
   end subroutine ed_site_get_patch_int

   subroutine ed_site_free(sh) bind(c, name="ed_site_free")
      integer(c_int), value, intent(in) :: sh
      call free_polygon(g_poly(sh))
      site_used(sh) = .false.
   end subroutine ed_site_free

end module ed_c_api_demography
