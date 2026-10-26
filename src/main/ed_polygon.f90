! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_polygon -- one polygon of a run and the step that advances it                         !
! (docs/science/order_of_processes.md §2).                                                 !
!                                                                                          !
! A polygon is two records side by side in the run: its state (polygon_t) and what the run !
! keeps to step and check it (polygon_run_t: its forcing cursor, its conservation ledgers, !
! its status). The second is not part of the state because its types come from the fast    !
! and slow libraries, which use the state library. A run is one polygon or many, each      !
! stepped by polygon_step, so a polygon of a region computes what a one-polygon run at its !
! cell computes.                                                                           !
!                                                                                          !
! polygon_step does no file work: the caller loads the step's forcing before it            !
! (met_prefetch) and writes the queued output records after a month closes.                !
!==========================================================================================!
module ed_polygon
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ed_kinds,                      only : wp, ik
   use ed_config,                     only : ed_config_t
   use ed_time,                       only : ed_time_t, calendar_step_t, time_to_string
   use ed_budget_check,               only : budget_report
   use ed_cohort_state_vars,          only : cohort_bad_field, cohort_field_finite, cohort_field_name
   use ed_patch_state_vars,           only : canopy_air_bad_field, canopy_air_field_name
   use ed_polygon_state_vars,         only : polygon_t, reset_step_diagnostics
   use ed_stepper,                    only : run_ledgers_t, advance_one_step
   use ed_slow_update,                only : max_cohort_count
   use ed_slow_ledger,                only : slow_ledger_report
   use ed_soil_carbon,                only : soil_carbon_bad_pool
   use ed_biogeochem_params,          only : SOIL_POOL_NAME
   use ed_met_types,                  only : met_source_t, met_cursor_t
   use ed_met,                        only : met_cursor_init
   use ed_output_types,               only : output_files_t, output_buffers_t
   use ed_output_integrate,           only : output_tick
   use ed_output_reduce,              only : total_area
   implicit none
   private

   public :: polygon_run_t, polygon_prepare, polygon_step, polygon_report
   public :: RUN_OK, RUN_FINISHED, RUN_ERR_NAN, RUN_ERR_AREA, RUN_ERR_BOUNDS

   !----- A step's or a run's status. A failure is a status, not an `error stop`, so a library     !
   !      caller (Python) survives it. ------------------------------------------------------------!
   integer(ik), parameter :: RUN_OK         = 0_ik   !< a slow step was taken
   integer(ik), parameter :: RUN_FINISHED   = 1_ik   !< the calendar had reached end_time; nothing done
   integer(ik), parameter :: RUN_ERR_NAN    = 2_ik   !< NaN in the state or the energy/water ledger
   integer(ik), parameter :: RUN_ERR_AREA   = 3_ik   !< patch areas do not sum to 1 (finalize only)
   integer(ik), parameter :: RUN_ERR_BOUNDS = 4_ik   !< a cohort or soil-carbon state outside its bounds

   !----- What the run keeps to step and check one polygon; its state is the polygon_t beside      !
   !      it. -------------------------------------------------------------------------------------!
   type :: polygon_run_t
      character(len=40)   :: label = ''                   !< how messages name it ('' for a site run)
      integer(ik)         :: patch_threads = 1_ik         !< the patch loop's threads ([run].n_threads, or 1)
      type(met_cursor_t)  :: met_cur                      !< set only if forcing_on
      type(run_ledgers_t) :: ledgers                      !< its conservation ledgers over the run
      !----- The most patches, and the most cohorts in one patch, a restructuring left, and when:  !
      !      above max_patch or max_cohort, the fusion tolerance ceilings kept dissimilar ones     !
      !      apart. -------------------------------------------------------------------------------!
      integer(ik)         :: most_patches = 0_ik, most_cohorts = 0_ik
      character(len=19)   :: most_patches_when = '', most_cohorts_when = ''
      real(wp)            :: area_start = 0.0_wp          !< the patch areas' sum at the start
      integer(ik)         :: status = RUN_OK           !< the last step's status
   end type polygon_run_t

contains

   !---------------------------------------------------------------------------------------!
   ! polygon_prepare -- the run record of a polygon whose start is complete                !
   ! (init_polygon), on a freshly allocated record: its starting area, its patch threads,  !
   ! its forcing cursor at the polygon's cell and location, and whether its slow ledger    !
   ! books.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_prepare(cfg, met_src, cell, poly, run, patch_threads, verbose)
      type(ed_config_t),   intent(in)    :: cfg
      type(met_source_t),  intent(in)    :: met_src
      integer(ik),         intent(in)    :: cell            !< the polygon's cell in the forcing source's list
      type(polygon_t),     intent(in)    :: poly
      type(polygon_run_t), intent(inout) :: run
      integer(ik),         intent(in)    :: patch_threads   !< [run].n_threads for one polygon, 1 for a region's
      logical,             intent(in)    :: verbose
      run%area_start    = total_area(poly)
      run%patch_threads = patch_threads
      if (cfg%fast_biophysics_on) then
         if (cfg%forcing%forcing_on) then
            call met_cursor_init(met_src, run%met_cur, cell, poly%latitude_deg, poly%longitude_deg,      &
                                 poly%elevation_m)
         end if
         if (verbose) write(*,'(a)') ' fast  : sub-daily biophysics ON'
      end if
      run%ledgers%slow%active = cfg%slow_on .and. cfg%slow_ledger_on
   end subroutine polygon_prepare

   !---------------------------------------------------------------------------------------!
   ! polygon_step -- one slow step of one polygon, `step`: the coupled step (ed_stepper,   !
   ! in the order of processes), the output tick, the step's diagnostics zeroed and the    !
   ! restructuring a month boundary owes, and the state guard. The step's forcing is       !
   ! already loaded into `met_src`; closed output records are queued in the polygon's      !
   ! buffers, `out_bufs`.                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_step(cfg, met_src, out_files, out_bufs, poly, run, step, status)
      type(ed_config_t),      intent(in)    :: cfg
      type(met_source_t),     intent(in)    :: met_src
      type(output_files_t),   intent(in)    :: out_files
      type(output_buffers_t), intent(inout) :: out_bufs   !< this polygon's share of out_files
      type(polygon_t),        intent(inout) :: poly
      type(polygon_run_t),    intent(inout) :: run
      type(calendar_step_t),  intent(in)    :: step
      integer(ik),            intent(out)   :: status
      integer(ik) :: n_cohort_max
      logical     :: restructured

      !----- The coupled step (ed_stepper): the restructuring the last step's boundary owes, the   !
      !      thermal acclimation, the fast loop, which stages its sub-steps in `out_bufs` when the !
      !      output reads them, and the slow dynamics. The daily step neither adds nor removes     !
      !      patches or cohorts, so the most a restructuring left are counted after it. -----------!
      restructured = poly%restructure_pending
      call advance_one_step(poly, cfg, run%ledgers, met_src=met_src, met_cur=run%met_cur,               &
                            step_start=step%start, out_bufs=out_bufs, patch_threads=run%patch_threads)
      if (restructured) then
         if (poly%patch%n > run%most_patches) then
            run%most_patches      = poly%patch%n
            run%most_patches_when = time_to_string(step%start)
         end if
         n_cohort_max = max_cohort_count(poly)
         if (n_cohort_max > run%most_cohorts) then
            run%most_cohorts      = n_cohort_max
            run%most_cohorts_when = time_to_string(step%start)
         end if
      end if

      !----- Output: replay the staged FAST samples, then fold this step -- its fluxes over the    !
      !      step and its state at its end -- into the slower tiers' windows that hold its start,  !
      !      closing each period its end has left and queueing it for the I/O phase. --------------!
      if (out_files%enabled) call output_tick(out_files, out_bufs, poly, step, cfg%dt_slow, cfg%dt_fast)

      !----- The output has read the step: zero its diagnostics. A step that ends on a month       !
      !      boundary leaves the restructuring owed to the next step, so the ending period's       !
      !      records describe one stand and what the restructuring does belongs to the period that !
      !      begins (docs/science/diagnostics.md §4). ---------------------------------------------!
      call reset_step_diagnostics(poly)
      poly%restructure_pending  = step%new_month
      poly%restructure_new_year = step%new_year

      status = check_polygon_state(cfg, poly, run, step%end)
      run%status = status

   end subroutine polygon_step

   !---------------------------------------------------------------------------------------!
   ! check_polygon_state -- the guard every step ends with: the run's energy and water     !
   ! ledgers, each cohort's state, each patch's canopy air and, with soil carbon on, each  !
   ! soil-carbon pool (docs/science/order_of_processes.md §2). The first value found that  !
   ! is not finite is RUN_ERR_NAN, the first outside its bounds RUN_ERR_BOUNDS, reported   !
   ! with what and where. An operator keeps these in bounds by construction, so either is  !
   ! a bug, and stopping on the step it happens is what makes it findable.                 !
   !---------------------------------------------------------------------------------------!
   integer(ik) function check_polygon_state(cfg, poly, run, now) result(status)
      type(ed_config_t),   intent(in) :: cfg
      type(polygon_t),     intent(in) :: poly
      type(polygon_run_t), intent(in) :: run
      type(ed_time_t),     intent(in) :: now
      character(len=160) :: where
      integer(ik)        :: i, k
      status = RUN_OK
      !----- One NaN flux or store turns a ledger's running residual NaN for good. ----------------!
      if (cfg%fast_biophysics_on) then
         if (.not. (ieee_is_finite(run%ledgers%energy%resid_sum) .and.                             &
                    ieee_is_finite(run%ledgers%water%resid_sum))) then
            call report_failure(run, 'a NaN reached the energy/water ledger.', '',                   &
                                'in the step ending '//time_to_string(now))
            status = RUN_ERR_NAN
            return
         end if
      end if
      do i = 1_ik, poly%cohort%n
         k = cohort_bad_field(poly%cohort, i)
         if (k == 0_ik) cycle
         write(where,'(3a,i0,a,i0,a,i0)') 'field ', trim(cohort_field_name(k)), ', cohort id ',     &
               poly%cohort%global_id(i), ' (pft ', poly%cohort%pft(i), ') in patch ', poly%cohort%owner_patch(i)
         if (cohort_field_finite(poly%cohort, i, k)) then
            call report_failure(run, 'a cohort''s state is outside its physical bounds.', where,       &
                                'date '//time_to_string(now),                                          &
                                [character(len=64) :: 'Density, carbon and water are amounts: never negative, never', &
                                                      'infinite. Every operator keeps them so by construction.'])
            status = RUN_ERR_BOUNDS
         else
            call report_failure(run, 'a cohort''s state is not finite.', where, 'date '//time_to_string(now))
            status = RUN_ERR_NAN
         end if
         return
      end do
      if (allocated(poly%patch%cas)) then
         do i = 1_ik, poly%patch%n
            k = canopy_air_bad_field(poly%patch%cas(i))
            if (k == 0_ik) cycle
            write(where,'(3a,i0,a,i0)') 'field ', trim(canopy_air_field_name(k)), ', patch ', i, ' of ',  &
                  poly%patch%n
            call report_failure(run, 'a patch''s canopy air is not finite.', where, 'date '//time_to_string(now))
            status = RUN_ERR_NAN
            return
         end do
      end if
      if (cfg%soil_carbon%soil_carbon_on) then
         do i = 1_ik, poly%patch%n
            k = soil_carbon_bad_pool(poly%patch%soil_carbon(i))
            if (k == 0_ik) cycle
            write(where,'(3a,i0,a,i0)') 'pool ', trim(SOIL_POOL_NAME(k)), ', patch ', i, ' of ',   &
                  poly%patch%n
            call report_failure(run, 'a CENTURY soil-carbon pool is physically impossible.', where,       &
                                'date '//time_to_string(now),                                          &
                                [character(len=64) :: 'A carbon pool is a mass: it cannot be negative, and the',      &
                                                      'ceiling is ~500x the richest real soil. Conservation can',     &
                                                      'hold perfectly while this is true -- see the slow ledger.'])
            status = RUN_ERR_BOUNDS
            return
         end do
      end if
   end function check_polygon_state

   !----- The guard's report: what failed, which polygon (a region's), where in it, when, and a    !
   !      note. -----------------------------------------------------------------------------------!
   subroutine report_failure(run, what, where, when, note)
      type(polygon_run_t), intent(in) :: run
      character(len=*),    intent(in) :: what, where, when
      character(len=*),    intent(in), optional :: note(:)
      integer :: i
      write(*,'(2a)') ' ERROR: ', what
      if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
      if (len_trim(where) > 0) write(*,'(2a)') '        ', trim(where)
      write(*,'(2a)') '        ', when
      if (present(note)) then
         do i = 1, size(note)
            write(*,'(2a)') '        ', trim(note(i))
         end do
      end if
   end subroutine report_failure

   !---------------------------------------------------------------------------------------!
   ! polygon_report -- one polygon's end-of-run conservation reports: the two whole-column !
   ! budgets and the face closure, the fusion ceilings, and the slow ledger with the       !
   ! soil-carbon seam. Each is printed whether or not it ever breached: the run's signed   !
   ! cumulative residual is the number a per-step tolerance cannot see.                    !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_report(cfg, run)
      type(ed_config_t),   intent(in) :: cfg
      type(polygon_run_t), intent(in) :: run
      if (cfg%fast_biophysics_on) then
         call budget_report(run%ledgers%energy, 'whole_energy', 'J/m2',  'W/m2')
         call budget_report(run%ledgers%water,  'whole_water',  'kg/m2', 'kg/m2/s')
         !----- The face closure should be machine zero: a nonzero value means soil enthalpy moved !
         !      on a mass flux the committed water never moved, a misplacement between layers that !
         !      both whole-column budgets sum away. -----------------------------------------------!
         write(*,'(a,es12.3,a,es12.3,a,i0)') ' faces[soil_layer_mass]  worst = ',                  &
               run%ledgers%face%worst, ' kg/m2   mean |resid| = ',                                   &
               run%ledgers%face%abs_sum / max(real(run%ledgers%face%n_check, wp), 1.0_wp),            &
               ' kg/m2   checks = ', run%ledgers%face%n_check
         if (run%ledgers%face%n_fail > 0_ik)                                                        &
            write(*,'(a,i0,a)') ' WARNING: ', run%ledgers%face%n_fail,                                &
               ' soil-water solves did not converge (see [fast].debug_error to make this fatal)'
         if (run%ledgers%energy%n_fail + run%ledgers%water%n_fail > 0_ik)                          &
            write(*,'(a,i0,a)') ' WARNING: ', run%ledgers%energy%n_fail + run%ledgers%water%n_fail,   &
               ' whole-column budget checks breached tolerance (see [fast].debug_error to make this fatal)'
      end if
      !----- Fusion targets the tolerance ceilings overrode. Silent when every restructuring      !
      !      reached them, which is the usual case.  -------------------------------------------------!
      if (cfg%demography%patch_fusion_on .and. run%most_patches > cfg%demography%max_patch)               &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: patch fusion left ', run%most_patches, ' patches on ',         &
               trim(run%most_patches_when), ' (max_patch = ', cfg%demography%max_patch,                            &
               '): patch_light_tol_max kept dissimilar patches apart'
      if (cfg%demography%cohort_fusion_on .and. run%most_cohorts > cfg%demography%max_cohort)             &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: cohort fusion left ', run%most_cohorts, ' cohorts in one patch on ',  &
               trim(run%most_cohorts_when), ' (max_cohort = ', cfg%demography%max_cohort,                          &
               '): cohort_size_tol_max kept dissimilar cohorts apart'
      call slow_ledger_report(run%ledgers%slow)
   end subroutine polygon_report

end module ed_polygon
