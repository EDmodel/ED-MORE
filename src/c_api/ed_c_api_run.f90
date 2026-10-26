! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_c_api_run -- the ISO_C_BINDING shim that exposes the FULL COUPLED MODEL (`ed.model`).     !
!                                                                                          !
! The other shims expose pieces: `ed_c_api_leaf`, `ed_c_api_phenology` and                        !
! `ed_c_api_canopy` the kernels, and `ed_c_api_demography` the slow loop alone (no                !
! biophysics, no forcing, no output streams). This one drives what the `ed_main` EXECUTABLE       !
! drives -- coupled fast + slow, live met forcing, netCDF output, both conservation ledgers -- by !
! holding a run (`ed_region_t`) and handing the caller the OPEN / STEP / FINALIZE seam. The       !
! getters read the run's first polygon: the site of a one-polygon run.                            !
!                                                                                          !
! The caller owns the time loop: a Python process opens a run, steps it, reads the site state as   !
! it goes, and finalizes it. `python -m ed.model CONFIG` is that loop, and                         !
! scripts/calibrate_fast runs its trials through it. Nothing here re-implements the model --       !
! `region_step` takes one slow step of what the executable steps a month at a time, so a          !
! Python-driven run and a `ed_main` run of the same config produce byte-identical output.          !
!                                                                                          !
! Conventions match the other shims: `polygon_t` and friends never cross the boundary. A run is an  !
! opaque integer handle into a module-`save` registry; scalars return by value; arrays COPY into a  !
! caller-allocated buffer. The bind(c) entry points are also PUBLIC Fortran so `test_c_api_run`     !
! can call them -- a shim no test calls can stop compiling without anyone noticing.                 !
!                                                                                          !
! `region_step` RETURNS a NaN status rather than `error stop`ping, which is what makes this safe to  !
! call from an interpreter: a bad state surfaces as a Python exception, not a dead process.           !
!==========================================================================================!
module ed_c_api_run
   use iso_c_binding
   use ed_kinds,             only :   wp, ik
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ed_diag_vars,         only :   CO_NPP_LEAF, CO_NPP_FINEROOT, CO_NPP_WOOD, CO_NPP_STORAGE, CO_NPP_REPRO, CO_GROWTH_RESP, &
                                      CO_EXUDATE
   use ed_region,            only :   ed_region_t, region_open, region_step, region_finalize,  &
                                      region_free, region_done
   use ed_polygon_state_vars,   only : patch_polygon_area
   use ed_output_reduce, only   : total_agb, total_lai, total_nplant, total_basal_area
   implicit none
   private

   public :: ed_run_open,   ed_run_step, ed_run_is_done, ed_run_finalize, ed_run_free
   public :: ed_run_year,   ed_run_month, ed_run_day, ed_run_istep, ed_run_iyear
   public :: ed_run_n_patch, ed_run_n_cohort, ed_run_n_site
   public :: ed_run_total_agb, ed_run_total_lai, ed_run_total_nplant
   public :: ed_run_total_basal_area, ed_run_soil_carbon,     ed_run_get_real, ed_run_get_int
   public :: ed_run_get_patch_real, ed_run_get_patch_int,     ed_run_get_site_real

   !----- Small registry: a handful of concurrent runs is plenty (the example opens two, one per   !
   !      stage, and could hold both at once). Fixed-size and module-`save` for the same reason    !
   !      ed_c_api_demography's is: the handle has to survive across ctypes calls.                 !
   integer, parameter :: MAXR = 4
   type(ed_region_t), target,   save :: g_run(MAXR)
   logical,                  save :: run_used(MAXR) = .false.

contains

   !----- Build a Fortran string from an explicit-length C char buffer (no strlen needed). -------!
   pure function f_string(cbuf, n) result(s)
      character(kind=c_char), intent(in) :: cbuf(*)
      integer(c_int),         intent(in) :: n
      character(len=n) :: s
      integer :: i
      do i = 1, n
         s(i:i) = cbuf(i)
      end do
   end function f_string

   !----- Guard every accessor: a handle the caller never opened (or already freed) must return a  !
   !      sentinel, not index a stale registry slot.  --------------------------------------------!
   pure logical function live(h)
      integer(c_int), intent(in) :: h
      live = (h >= 1_c_int .and. h <= int(MAXR, c_int))
      if (live) live = run_used(h)
   end function live

   !=======================================================================================!
   !  Lifecycle                                                                             !
   !=======================================================================================!

   !----- Open a run from a TOML config -> opaque handle (>=1); -1 registry full, -2 open failed. -!
   !      `verbose` /= 0 keeps the progress lines ed_main prints on stdout.                        !
   function ed_run_open(path, path_len, verbose) result(h) bind(c, name="ed_run_open")
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value,  intent(in) :: path_len, verbose
      integer(c_int)                     :: h
      integer :: i
      logical :: ok
      h = -1_c_int
      do i = 1, MAXR
         if (.not. run_used(i)) then
            h = int(i, c_int)
            exit
         end if
      end do
      if (h < 0_c_int) return
      call region_open(f_string(path, path_len), g_run(h), ok, verbose=(verbose /= 0_c_int))
      if (.not. ok) then
         h = -2_c_int
         return
      end if
      run_used(h) = .true.
   end function ed_run_open

   !----- Advance ONE slow step. Returns the run status (RUN_* in ed_polygon): 0 stepped,       !
   !      1 already finished, 2 NaN in the state, 4 a cohort or soil-carbon state out of bounds. -!
   function ed_run_step(h)   result(status) bind(c, name="ed_run_step")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: status
      integer(ik) :: st
      if (.not. live(h)) then
         status = -1_c_int
         return
      end if
      call region_step(g_run(h), st)
      status = int(st, c_int)
   end function ed_run_step

   !----- 1 when the calendar has reached end_time. The loop predicate, so the CALLER writes the   !
   !      `while` -- Python owning the time loop is the point of this shim.  ---------------------!
   function ed_run_is_done(h) result(d) bind(c, name="ed_run_is_done")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: d
      if (.not. live(h)) then
         d = 1_c_int
         return
      end if
      d = merge(1_c_int, 0_c_int, region_done(g_run(h)))
   end function ed_run_is_done

   !----- Terminal checkpoint + conservation reports + close the streams. The state stays readable !
   !      afterwards (that is why free is separate).  --------------------------------------------!
   function ed_run_finalize(h) result(status) bind(c, name="ed_run_finalize")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: status
      integer(ik) :: st
      if (.not. live(h)) then
         status = -1_c_int
         return
      end if
      call region_finalize(g_run(h), st)
      status = int(st, c_int)
   end function ed_run_finalize

   subroutine ed_run_free(h) bind(c, name="ed_run_free")
      integer(c_int), value, intent(in) :: h
      if (.not. live(h)) return
      call region_free(g_run(h))
      run_used(h) = .false.
   end subroutine ed_run_free

   !=======================================================================================!
   !  Where the calendar is                                                                 !
   !=======================================================================================!
   function ed_run_year(h)   result(v) bind(c, name="ed_run_year")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%now%year, c_int)
   end function ed_run_year

   function ed_run_month(h) result(v)   bind(c, name="ed_run_month")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%now%month, c_int)
   end function ed_run_month

   function ed_run_day(h) result(v) bind(c, name="ed_run_day")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%now%day, c_int)
   end function ed_run_day

   function ed_run_istep(h) result(v)   bind(c, name="ed_run_istep")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%istep, c_int)
   end function ed_run_istep

   function ed_run_iyear(h) result(v)   bind(c, name="ed_run_iyear")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%iyear, c_int)
   end function ed_run_iyear

   !=======================================================================================!
   !  Site aggregates -- the same reducers print_summary uses, so a Python progress line and !
   !  the executable's own summary line report the same numbers.                             !
   !=======================================================================================!
   function ed_run_n_patch(h) result(v) bind(c, name="ed_run_n_patch")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%poly(1)%patch%n, c_int)
   end function ed_run_n_patch

   function ed_run_n_site(h) result(v) bind(c, name="ed_run_n_site")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%poly(1)%site%n, c_int)
   end function ed_run_n_site

   function ed_run_n_cohort(h) result(v)   bind(c, name="ed_run_n_cohort")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int
      if (live(h)) v = int(g_run(h)%poly(1)%cohort%n, c_int)
   end function ed_run_n_cohort

   function ed_run_total_agb(h) result(v) bind(c, name="ed_run_total_agb")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double
      if (live(h)) v = real(total_agb(g_run(h)%poly(1)), c_double)
   end function ed_run_total_agb

   function ed_run_total_lai(h) result(v) bind(c, name="ed_run_total_lai")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double
      if (live(h)) v = real(total_lai(g_run(h)%poly(1)), c_double)
   end function ed_run_total_lai

   function ed_run_total_nplant(h) result(v) bind(c, name="ed_run_total_nplant")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double
      if (live(h)) v = real(total_nplant(g_run(h)%poly(1)), c_double)
   end function ed_run_total_nplant

   function ed_run_total_basal_area(h) result(v) bind(c, name="ed_run_total_basal_area")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double
      if (live(h)) v = real(total_basal_area(g_run(h)%poly(1)), c_double)
   end function ed_run_total_basal_area

   !----- Site total soil carbon [kgC/m2], area-weighted over patches: all seven CENTURY pools.    !
   !      Identically zero unless [soil_carbon].soil_carbon_on -- which is exactly what makes it    !
   !      worth watching from the driver loop, since with the feature off litter is DISCARDED       !
   !      rather than stored, and this number stays flat at zero while the stand grows.             !
   function ed_run_soil_carbon(h) result(v)   bind(c, name="ed_run_soil_carbon")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      real(wp)    :: tot
      integer(ik) :: ip
      v = 0.0_c_double
      if (.not. live(h)) return
      tot = 0.0_wp
      associate (p => g_run(h)%poly(1)%patch)
         do ip = 1_ik, p%n
            tot = tot + patch_polygon_area(g_run(h)%poly(1), ip) * sum(p%soil_carbon(ip)%carbon)
         end do
      end associate
      v = real(tot, c_double)
   end function ed_run_soil_carbon

   !=======================================================================================!
   !  Copy-out per-cohort getters (caller allocates buf of length n_cohort). Field ids match !
   !  ed_c_api_demography's, so the two shims read the same way from Python.                 !
   !=======================================================================================!
   subroutine ed_run_get_real(h, field_id, buf) bind(c, name="ed_run_get_real")
      integer(c_int), value, intent(in)  :: h, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n
      if (.not. live(h)) return
      associate (c => g_run(h)%poly(1)%cohort)
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
         !----- Run only (a demography site has no fast loop): gross GPP over the last slow step. --!
         case (8_c_int)
            buf(1:n) = real(c%gpp_accum(1:n),       c_double)
         !----- Run only: the maintenance respiration over the last slow step [kgC/plant] the        !
         !      allocator charged against that GPP -- leaves, stems (coarse roots included), fine     !
         !      roots. ----------------------------------------------------------------------------!
         case (11_c_int)
            buf(1:n) = real(c%leaf_resp_accum(1:n), c_double)
         case (12_c_int)
            buf(1:n) = real(c%stem_resp_accum(1:n), c_double)
         case (13_c_int)
            buf(1:n) = real(c%root_resp_accum(1:n), c_double)
         !----- The live pools [kgC/plant] and the baseline leaf turnover rate [1/yr]. ----------!
         case (14_c_int)
            buf(1:n) = real(c%leaf_carbon(1:n),          c_double)
         case (15_c_int)
            buf(1:n) = real(c%fineroot_carbon(1:n),      c_double)
         case (16_c_int)
            buf(1:n) = real(c%nonstructural_carbon(1:n), c_double)
         case (17_c_int)
            buf(1:n) = real(c%leaf_turnover_rate(1:n),   c_double)
         !----- The last slow step's allocation [kgC/plant/yr]: NPP (leaf, fine-root, wood,        !
         !      reproduction and storage growth) and its parts, the growth respiration and the root  !
         !      exudate; NaN unless the output needs the slow cohort diagnostics that hold them (a   !
         !      site NPP total in the output list switches them on). --------------------------------!
         case (9_c_int, 10_c_int, 18_c_int:23_c_int)
            if (.not. c%sdiag%active) then
               buf(1:n) = ieee_value(1.0_c_double, ieee_quiet_nan)
               return
            end if
            select case (field_id)
            case (9_c_int)
               buf(1:n) = real(c%sdiag%v(CO_NPP_LEAF, 1:n) + c%sdiag%v(CO_NPP_FINEROOT, 1:n)          &
                             + c%sdiag%v(CO_NPP_WOOD, 1:n) + c%sdiag%v(CO_NPP_STORAGE, 1:n)          &
                             + c%sdiag%v(CO_NPP_REPRO, 1:n), c_double)
            case (10_c_int)
               buf(1:n) = real(c%sdiag%v(CO_NPP_WOOD,     1:n), c_double)
            case (18_c_int)
               buf(1:n) = real(c%sdiag%v(CO_NPP_LEAF,     1:n), c_double)
            case (19_c_int)
               buf(1:n) = real(c%sdiag%v(CO_NPP_FINEROOT, 1:n), c_double)
            case (20_c_int)
               buf(1:n) = real(c%sdiag%v(CO_NPP_STORAGE,  1:n), c_double)
            case (21_c_int)
               buf(1:n) = real(c%sdiag%v(CO_NPP_REPRO,    1:n), c_double)
            case (22_c_int)
               buf(1:n) = real(c%sdiag%v(CO_GROWTH_RESP,  1:n), c_double)
            case (23_c_int)
               buf(1:n) = real(c%sdiag%v(CO_EXUDATE,      1:n), c_double)
            end select
         end select
      end associate
   end subroutine ed_run_get_real

   subroutine ed_run_get_int(h, field_id,   buf) bind(c, name="ed_run_get_int")
      integer(c_int), value, intent(in)  :: h, field_id
      integer(c_int),        intent(out) :: buf(*)
      integer(ik) :: n
      if (.not. live(h)) return
      associate (c => g_run(h)%poly(1)%cohort)
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
   end subroutine ed_run_get_int

   !----- Copy-out per-patch getters (caller allocates buf of length n_patch). A cohort's     !
   !      owner_patch indexes them, so a polygon total of a cohort quantity is                 !
   !      sum(share(owner_patch) * nplant * x), the share being the patch's fraction of the     !
   !      polygon: its area in its site times the site's area. ----------------------------------!
   subroutine ed_run_get_patch_real(h,   field_id, buf) bind(c, name="ed_run_get_patch_real")
      integer(c_int), value, intent(in)  :: h, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n, ip
      if (.not. live(h)) return
      associate (p => g_run(h)%poly(1)%patch)
         n = p%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int)
            buf(1:n) = real(p%area(1:n), c_double)                 ! fraction of its site
         case (1_c_int)
            buf(1:n) = real(p%age(1:n),  c_double)                 ! [yr] since disturbance
         case (2_c_int)                                            ! fraction of the polygon
            do ip = 1_ik, n
               buf(ip) = real(patch_polygon_area(g_run(h)%poly(1), ip), c_double)
            end do
         end select
      end associate
   end subroutine ed_run_get_patch_real

   subroutine ed_run_get_patch_int(h, field_id, buf) bind(c, name="ed_run_get_patch_int")
      integer(c_int), value, intent(in)  :: h, field_id
      integer(c_int),        intent(out) :: buf(*)
      integer(ik) :: n
      if (.not. live(h)) return
      associate (p => g_run(h)%poly(1)%patch)
         n = p%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int)                                            ! its site (1-based)
            buf(1:n) = int(p%owner_site(1:n), c_int)
         end select
      end associate
   end subroutine ed_run_get_patch_int

   !----- Copy-out per-site getter (caller allocates buf of length n_site). ------------------!
   subroutine ed_run_get_site_real(h, field_id, buf) bind(c, name="ed_run_get_site_real")
      integer(c_int), value, intent(in)  :: h, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n
      if (.not. live(h)) return
      associate (s => g_run(h)%poly(1)%site)
         n = s%n
         select case (field_id)
         case (0_c_int)                                            ! fraction of the polygon
            buf(1:n) = real(s%area(1:n), c_double)
         end select
      end associate
   end subroutine ed_run_get_site_real

end module ed_c_api_run
