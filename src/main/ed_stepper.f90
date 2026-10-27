! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_stepper -- one coupled slow step of one polygon, in the order of processes
! (docs/science/order_of_processes.md §2): the restructuring the last step's boundary
! owes, the thermal acclimation, the fast loop, the slow dynamics. A library of its own,
! so the tests drive the same step the run does.
!---~---
module ed_stepper
   use ed_kinds,                only : ik
   use ed_constants,            only : day_sec
   use ed_config,               only : ed_config_t
   use ed_time,                 only : ed_time_t
   use ed_budget_check,         only : budget_t
   use ed_polygon_state_vars,   only : polygon_t
   use ed_slow_dynamics,        only : advance_slow_dynamics, advance_boundary_dynamics, settle_tissue_water
   use ed_slow_state,           only : acclimate_site_leaf_tables
   use ed_slow_ledger,          only : slow_ledger_t
   use ed_fast_dynamics,        only : fast_dynamics
   use ed_met_types,            only : met_source_t, met_cursor_t
   use ed_output_types,         only : output_buffers_t
   implicit none
   private

   public :: run_ledgers_t, advance_one_step

   !---~---
   !   A polygon's conservation ledgers over the run: the fast loop's whole-column energy and
   ! water budgets and each soil layer face's water (did the heat move with it), and the slow
   ! step's ledger of the site's stores, which also keeps the soil-carbon seam's worst
   ! values. One never activated books nothing, so a test can pass a fresh one.
   !---~---
   type :: run_ledgers_t
      type(budget_t)      :: energy, water, face
      type(slow_ledger_t) :: slow
   end type run_ledgers_t

contains

   !---~---
   !   Advance one slow step of one polygon. The fast loop takes the forcing when forcing_on
   ! and the source, cursor and step_start are all given, and otherwise the reference
   ! climate; leaf phenology advances when step_start is given. A test passes what it
   ! needs; the run passes all.
   !---~---
   subroutine advance_one_step(poly, cfg, ledgers, met_src, met_cur, step_start, out_bufs, patch_threads)
      type(polygon_t),        intent(inout) :: poly
      type(ed_config_t),    intent(in)      :: cfg
      type(run_ledgers_t),    intent(inout) :: ledgers
      type(met_source_t),     intent(in),    optional :: met_src        !--- the run's forcing source
      type(met_cursor_t),     intent(inout), optional :: met_cur        !--- this polygon's forcing cursor
      type(ed_time_t),      intent(in),      optional :: step_start     !--- the step's start date
      type(output_buffers_t), intent(inout), optional :: out_bufs       !--- this polygon's FAST-tier staging
      integer(ik),            intent(in),    optional :: patch_threads  !--- the patch loop's threads (default 1)

      !---~---
      !   The restructuring the last step's boundary owes (docs/science/order_of_processes.md
      ! §5), before this step's fast loop sees the stand, unless the master slow_on freeze
      ! holds the stand.
      !---~---
      if (poly%restructure_pending) then
         if (cfg%slow_on) call advance_boundary_dynamics(poly, cfg, .true., poly%restructure_new_year,   &
                                                         ledgers%slow)
         poly%restructure_pending  = .false.
         poly%restructure_new_year = .false.
      end if

      !---~---
      !   Thermal acclimation (docs/science/leaf_gas_exchange.md §2): each site's
      ! growth-temperature mean takes the daily mean the last step's fast loop accumulated,
      ! and its leaf table follows. The mean tracks the preceding weeks, so one step's delay
      ! is a thirtieth of its window.
      !---~---
      if (cfg%leaf_physiology%thermal_acclimation)                                                  &
         call acclimate_site_leaf_tables(poly, cfg%leaf_physiology, max(1_ik, nint(cfg%dt_slow / day_sec, ik)))

      !---~---
      !   The fast loop. The tissue water above what each cohort's biomass now holds is shed
      ! first, and the slow ledger books it (docs/science/order_of_processes.md §2).
      !---~---
      if (cfg%fast_biophysics_on) then
         call settle_tissue_water(poly, cfg, ledgers%slow)
         call fast_dynamics(poly, cfg, met_src=met_src, met_cur=met_cur, step_start=step_start,          &
                            out_bufs=out_bufs, run_energy_budget=ledgers%energy,                       &
                            run_water_budget=ledgers%water, run_face_budget=ledgers%face,               &
                            patch_threads=patch_threads)
      end if

      !---~---
      !   The slow dynamics, unless slow_on freezes the stand while the fast loop runs: leaf
      ! phenology, the vegetation dynamics and the soil carbon.
      !---~---
      if (cfg%slow_on) call advance_slow_dynamics(poly, cfg, ledgers%slow, step_start)
   end subroutine advance_one_step

end module ed_stepper
