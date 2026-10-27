! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_slow_dynamics -- the order of the slow loop (docs/science/order_of_processes.md
! §4-§5), and where the slow ledger marks its phases (conservation.md §3). It loops over
! nothing: ed_slow_state reads the stand and calls the processes, ed_slow_update
! writes the stand, and ed_slow_diag writes the slow rows.
!   * advance_slow_dynamics     -- one slow step: vegetation_dynamics, patch ageing, the
!                                  canopy-air depth, then the soil-carbon step.
!   * vegetation_dynamics       -- phenology, plasticity, allocation, the vital rates, and
!                                  the commit of growth and mortality.
!   * advance_boundary_dynamics -- a calendar boundary: restructure_stand (the monthly
!                                  cohort operators, the yearly disturbance and patch
!                                  operators), then the canopy-air depth.
!   * settle_tissue_water       -- before the fast loop, the tissue water above capacity.
!---~---
module ed_slow_dynamics
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_time,                 only : ed_time_t, day_of_year
   use ed_polygon_state_vars,   only : polygon_t
   use ed_patch_state_vars,     only : soil_carbon_state_t
   use ed_slow_ledger,          only : slow_ledger_t, slow_ledger_open, slow_ledger_heat_open,           &
                                       slow_ledger_heat_close, slow_ledger_mark, slow_ledger_rebase,     &
                                       slow_ledger_declare, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC,         &
                                       SLOW_PHASE_TISSUE_WATER, slow_fast_carbon_handover,               &
                                       SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW, SLOW_PHASE_RECRUIT,         &
                                       SLOW_PHASE_COHORT, SLOW_PHASE_DISTURB, SLOW_PHASE_PATCH,          &
                                       slow_co2_handoff, litter_carbon_total
   use ed_slow_state,           only : slow_step_data_t, decide_recruitment, decide_disturbance,         &
                                       kill_report_t, kill_litter, clear_litter,                         &
                                       advance_leaf_phenology, advance_plant_traits,                     &
                                       compute_carbon_allocation, shed_turnover_water,                   &
                                       compute_vital_rates, fill_new_leaf_water,                         &
                                       accumulate_mortality_litter, shed_mortality_water,                &
                                       shed_excess_tissue_water, advance_soil_carbon
   use ed_slow_update,          only : accumulate_recruit_pool, commit_vegetation_step, sort_cohorts,    &
                                       update_overtopping_lai, age_patches, refresh_canopy_depth,        &
                                       apply_recruitment, fuse_cohorts, terminate_cohorts,               &
                                       split_cohorts, apply_patch_disturbance, sort_patches,             &
                                       fuse_patches, terminate_patches, commit_trait_carbon,             &
                                       commit_soil_carbon
   use ed_slow_diag,            only : open_slow_cohort_rows, write_slow_patch_rows
   implicit none
   private

   public :: advance_slow_dynamics, advance_boundary_dynamics, settle_tissue_water
   public :: vegetation_dynamics, restructure_stand

   !--- Patch disturbance runs once a year, so it integrates its yearly rate over this interval.
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   !--- [yr]

contains

   !---~---
   !   Advance one slow step (order_of_processes.md §4): vegetation dynamics, patch ageing,
   ! the canopy-air depth, then the daily soil-carbon step. Without `step_start` (a bare
   ! test call) phenology is not advanced.
   !---~---
   subroutine advance_slow_dynamics(poly, cfg, ledger, step_start)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_ledger_t), intent(inout) :: ledger   !--- the slow ledger; an inactive one books nothing
      type(ed_time_t),   intent(in), optional :: step_start   !--- the slow step's start (its day of year)
      type(soil_carbon_state_t), allocatable :: soil_pools(:)

      call slow_ledger_open(ledger, poly, cfg, step_start)

      if (present(step_start)) then
         call vegetation_dynamics(poly, cfg, ledger, day_of_year(step_start))
      else
         call vegetation_dynamics(poly, cfg, ledger)
      end if

      call age_patches(poly%patch, cfg%dt_years)

      !---~---
      !   Resize each canopy-air volume to the stand growth and mortality left, so the
      ! fast loop sees one depth all day.
      !---~---
      call refresh_canopy_depth(poly, cfg, ledger)
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)

      if (cfg%soil_carbon%soil_carbon_on) then
         call advance_soil_carbon(poly, cfg, ledger, soil_pools)
         call commit_soil_carbon(poly%patch, soil_pools)
      end if
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_SOILC)
      call clear_litter(poly)
   end subroutine advance_slow_dynamics

   !---~---
   !   Advance the vegetation one slow step: phenology, trait plasticity, carbon allocation, the
   ! vital rates, and the growth and mortality commit. The calendar's restructuring runs at the
   ! boundary, in restructure_stand (docs/science/order_of_processes.md §4).
   !---~---
   subroutine vegetation_dynamics(poly, cfg, ledger, doy)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_ledger_t), intent(inout) :: ledger   !--- the slow ledger; an inactive one books nothing
      integer(ik),         intent(in), optional :: doy   !--- day-of-year at the step start (drives phenology)
      type(slow_step_data_t) :: step
      real(wp)               :: mort_water, seed_rain_carbon, handover, litter_open
      real(wp), allocatable  :: leaf_carbon_new(:), storage_new(:)

      !---~---
      !   0. Leaf phenology: advance the leaf tendencies one daily step, before the allocation
      ! reads them (docs/science/plant_phenology.md). A caller with no calendar (the C API's
      ! slow-only path) omits doy, and the tendencies are not advanced.
      !---~---
      if (present(doy)) call advance_leaf_phenology(poly, cfg, doy)

      !--- The litter the step opens with, the last boundary's kills, was declared there.
      litter_open = 0.0_wp
      if (ledger%active) litter_open = litter_carbon_total(poly, poly%patch%litter_in)

      !---~---
      !   0b. Light trait plasticity (opt-in): acclimate the leaf traits toward their shaded
      ! targets from last step's overtopping LAI, before the allocation reads them. The
      ! leaf carbon an SLA rise resorbs to storage (no heat capacity) changes the tissue
      ! energy at unchanged temperature: declared (docs/science/conservation.md §4).
      !---~---
      if (cfg%trait_dynamics%trait_plasticity_on) then
         call advance_plant_traits(poly, cfg%dt_years, leaf_carbon_new, storage_new)
         call slow_ledger_heat_open(ledger, poly)
         call commit_trait_carbon(poly, leaf_carbon_new, storage_new)
         call slow_ledger_heat_close(ledger, poly)
      end if

      !--- 1. Carbon allocation, which adds the turnover litter to the per-patch accumulator.
      call open_slow_cohort_rows(poly)
      call compute_carbon_allocation(poly, cfg, step)

      !---~---
      !   1b. Shed the turnover's tissue water before the commit, which changes the pools the
      ! loss is measured against (docs/science/soil_plant_water.md §2).
      !---~---
      call shed_turnover_water(poly, cfg, step, ledger)
      !---~---
      !   Every movement of the allocation phase is declared, so a residual here is a movement
      ! nobody declared (docs/science/conservation.md §3).
      !---~---
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_ALLOCATE)

      !---~---
      !   2. The vital rates, before the commit: mortality reads the growth average of the steps
      ! before this one (docs/science/demography.md §1).
      !---~---
      call compute_vital_rates(poly, cfg, step)

      !---~---
      !   2b. Credit the recruit pool with this step's recruitment, and declare the seed rain
      ! as carbon in: no parent at the site paid for it (docs/science/demography.md §2).
      !---~---
      call accumulate_recruit_pool(poly, cfg, step%recruitment, cfg%dt_years, seed_rain_carbon)
      call slow_ledger_declare(ledger, carbon_in = seed_rain_carbon)

      !---~---
      !   3. Growth and mortality. The commit and the new leaves' water leave the temperatures
      ! alone, so the energy change is the tissues' thermal mass: one exchange, declared,
      ! not a leak (docs/science/conservation.md §4). The fast loop's net carbon is valued
      ! at the pre-commit density it accumulated on; after the commit it would miss the
      ! dying plants' share (conservation.md §3).
      !---~---
      handover = 0.0_wp
      if (ledger%active) handover = slow_fast_carbon_handover(poly, cfg)
      call slow_ledger_heat_open(ledger, poly)
      call commit_vegetation_step(poly, cfg, step)
      call fill_new_leaf_water(poly, cfg, step)
      call slow_ledger_heat_close(ledger, poly)

      !---~---
      !   Mortality's litter and water, on the density the commit removed and the pools it left
      ! (docs/science/demography.md eq. (9)); recorded whether or not soil carbon is on.
      !---~---
      call accumulate_mortality_litter(poly, cfg, step)
      call shed_mortality_water(poly, cfg, step, mort_water)

      !--- Growth changed heights: re-sort tallest first (docs/science/demography.md §3).
      call sort_cohorts(poly)

      !---~---
      !   The growth phase declares the fast loop's net carbon (in), the litter and the CO2 owed
      ! to the canopy air (out) and the dying plants' water, so it closes to round-off
      ! (docs/science/conservation.md §3).
      !---~---
      if (ledger%active) then
         call slow_ledger_declare(ledger, carbon_in  = handover,                                    &
                                          carbon_out = litter_carbon_total(poly, poly%patch%litter_in)  &
                                                     - litter_open                                    &
                                                     + slow_co2_handoff(poly, cfg),                &
                                          water_out  = mort_water)
      end if
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_GROW)

      !--- The slow patch rows, before the restructuring permutes the patch axis.
      call write_slow_patch_rows(poly, cfg, step%recruitment)

      !--- 4. Refresh the overtopping LAI on the sorted stand (docs/science/demography.md §6).
      call update_overtopping_lai(poly)
   end subroutine vegetation_dynamics

   !---~---
   !   advance_boundary_dynamics -- the calendar boundary (order_of_processes.md §5):
   ! restructure the stand, then resize the canopy-air volumes; a no-op off a month
   ! boundary. Its marks continue the last step's ledger window.
   !---~---
   subroutine advance_boundary_dynamics(poly, cfg, is_new_month, is_new_year, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      type(slow_ledger_t), intent(inout) :: ledger   !--- the slow ledger; an inactive one books nothing
      if (.not. is_new_month) return
      if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
      call restructure_stand(poly, cfg, is_new_month, is_new_year, ledger)
      call refresh_canopy_depth(poly, cfg, ledger)
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)
   end subroutine advance_boundary_dynamics

   !---~---
   !   restructure_stand -- the calendar's restructuring: monthly recruitment, cohort
   ! fusion, cull, fission and sort; yearly patch disturbance, then patch fusion and
   ! termination. It runs at the boundary, after the output has read the period's last
   ! step, so each record describes one stand and the boundary's events belong to the
   ! period it opens (docs/science/order_of_processes.md §5). The operators change the
   ! tissues' thermal mass (fusion and fission re-derive the sapwood ring, nonlinear in
   ! dbh; the cull and the treefall remove cohorts): one exchange per phase, declared
   ! (conservation.md §4).
   !---~---
   subroutine restructure_stand(poly, cfg, is_new_month, is_new_year, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year   !--- which boundary this is
      type(slow_ledger_t), intent(inout) :: ledger   !--- the slow ledger; an inactive one books nothing
      real(wp) :: cull_water, rec_carbon, rec_heat, rec_water, dist_water, frac, kill_carbon
      logical  :: monthly, yearly
      logical, allocatable :: found(:,:), survives(:)
      type(kill_report_t)  :: kills

      !--- The calendar cadence; each operator also has its own switch.
      monthly = is_new_month .and. cfg%demography%demography_on
      yearly  = is_new_year  .and. cfg%demography%demography_on

      !---~---
      !   Cohort restructuring (monthly): recruitment, fusion, cull, fission and sort.
      ! Recruitment is marked on its own: it alone draws carbon, heat and water from
      ! outside the model, checkable only if no other operator shares its phase
      ! (docs/science/conservation.md §3).
      !---~---
      if (monthly) then
         if (cfg%demography%recruitment_on) then
            call decide_recruitment(poly, cfg, found)
            call apply_recruitment(poly, cfg, found, rec_carbon, rec_heat, rec_water)
            call slow_ledger_declare(ledger, carbon_in = rec_carbon, energy_in = rec_heat, water_in = rec_water)
         end if
         call slow_ledger_heat_open(ledger, poly)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_RECRUIT)
         if (cfg%demography%cohort_fusion_on) call fuse_cohorts(poly, cfg)
         call terminate_cohorts(poly, cfg, cull_water, kills)
         call kill_litter(poly, kills, kill_carbon)
         if (cfg%demography%cohort_fission_on) call split_cohorts(poly, cfg)
         call sort_cohorts(poly)
         call slow_ledger_declare(ledger, carbon_out = kill_carbon, water_out = cull_water)
         call slow_ledger_heat_close(ledger, poly)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_COHORT)
      end if

      !--- Treefall gaps, then patch fusion and termination (yearly).
      if (yearly .and. cfg%demography%patch_disturbance_on) then
         call decide_disturbance(poly, cfg, PATCH_DYNAMICS_INTERVAL, frac, survives)
         call slow_ledger_heat_open(ledger, poly)
         call apply_patch_disturbance(poly, cfg, frac, survives, dist_water, kills)
         call kill_litter(poly, kills, kill_carbon)
         call slow_ledger_declare(ledger, carbon_out = kill_carbon, water_out = dist_water)
         call slow_ledger_heat_close(ledger, poly)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_DISTURB)
      end if
      if (yearly) then
         call slow_ledger_heat_open(ledger, poly)
         call sort_patches(poly)
         if (cfg%demography%patch_fusion_on) call fuse_patches(poly, cfg)
         call terminate_patches(poly, cfg)
         if (cfg%demography%cohort_fusion_on) call fuse_cohorts(poly, cfg)
         call terminate_cohorts(poly, cfg, cull_water, kills)
         call kill_litter(poly, kills, kill_carbon)
         call sort_cohorts(poly)
         call slow_ledger_declare(ledger, carbon_out = kill_carbon, water_out = cull_water)
         call slow_ledger_heat_close(ledger, poly)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_PATCH)
      end if

      !--- Refresh the overtopping LAI for the restructured stand.
      call update_overtopping_lai(poly)
   end subroutine restructure_stand

   !---~---
   !   settle_tissue_water -- before the fast loop, shed the tissue water above what the
   ! biomass holds and declare it: the water as a handoff to the ground, the heat as a
   ! thermal-mass change (conservation.md §3-§4), in a phase of its own.
   !---~---
   subroutine settle_tissue_water(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_ledger_t), intent(inout) :: ledger   !--- the slow ledger; an inactive one books nothing
      real(wp) :: water_shed
      if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
      call slow_ledger_heat_open(ledger, poly)
      call shed_excess_tissue_water(poly, cfg, water_shed)
      call slow_ledger_declare(ledger, water_out = water_shed)
      call slow_ledger_heat_close(ledger, poly)
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_TISSUE_WATER)
   end subroutine settle_tissue_water

end module ed_slow_dynamics
