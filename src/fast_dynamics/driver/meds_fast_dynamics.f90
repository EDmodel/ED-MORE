! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_dynamics -- the polygon's fast loop over one slow step, the window             !
! (docs/science/order_of_processes.md §3). fast_dynamics samples the forcing, then for     !
! each patch gathers its state (meds_fast_state), runs n_fast_per_slow sub-steps of        !
! fast_step on its forcing (meds_fast_forcing), records each sub-step for the output       !
! (meds_fast_diag) and for the slow step, and writes the state back. The polygon's sums    !
! are folded in patch order after the loop.                                                !
!==========================================================================================!
module meds_fast_dynamics
   use meds_kinds,              only : wp, ik, lk
   use meds_numerics,           only : ascending_order
   use meds_constants,          only : latent_heat_vap
   use meds_config,             only : meds_config_t
   use meds_budget_check,       only : budget_t, budget_merge
   use meds_time,               only : meds_time_t
   use meds_output_types,       only : output_buffers_t
   use meds_met_types,          only : met_source_t, met_cursor_t, met_forcing_t
   use meds_met,                only : met_co2
   use meds_forcing_params,     only : CO2_SOURCE_FILE
   use meds_polygon_state_vars, only : polygon_t, patch_polygon_area
   use meds_fast_integ_types,         only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,   &
                                       fast_workspace_t, step_report_t, ensure_static_capacity,                  &
                                       ensure_dynamic_capacity, ensure_forcing_capacity, ensure_output_capacity, &
                                       N_WORK, WORK_INTEG_STEPS, WORK_INTEG_REJECT, WORK_WATER_NSUB,             &
                                       WORK_WATER_NEWTON, WORK_WATER_REJECT, WORK_WATER_NFAIL, WORK_RK45_RESCUE, &
                                       WORK_CLAMP_STAGE
   use meds_fast_integ_stepper,          only : fast_step
   use meds_fast_state,         only : reset_fast_to_slow, gather_patch, scatter_patch, accumulate_fast_to_slow,  &
                                       seed_start_canopy_air
   use meds_fast_forcing,       only : sample_window, site_forcing, patch_forcing
   use meds_fast_diag,          only : size_staging, record_substep, fold_window
   !$ use omp_lib,               only : omp_get_thread_num
   implicit none
   private

   public :: fast_dynamics

   !----- Slots of the reduction staging: each (patch, sub-step) writes its own cells, folded in patch !
   !      order after the loop, so the sums do not move with the thread count                       !
   !      (docs/science/numerical_scheme.md §6a). ---------------------------------------------------!
   integer(ik), parameter :: RED_ET          = 1_ik   !< polygon ET [kg/m2], area-weighted
   integer(ik), parameter :: RED_PHENO_TAIR  = 2_ik   !< daily-mean air-temperature accumulator [K]
   !----- The step report's work counters, area-weighted: slot RED_WORK + WORK_*. ------------------!
   integer(ik), parameter :: RED_WORK        = 2_ik
   !----- Rain and snow as split at each patch's canopy-air top, area-weighted (docs/science/forcing.md §7), !
   !      and the water leaving the soil column. ----------------------------------------------------------!
   integer(ik), parameter :: RED_RAINF       = RED_WORK + N_WORK + 1_ik   !< rain [kg/m2/s]
   integer(ik), parameter :: RED_SNOWFALL    = RED_WORK + N_WORK + 2_ik   !< snow [kg/m2/s]
   integer(ik), parameter :: RED_RUNOFF_SURF = RED_WORK + N_WORK + 3_ik   !< surface runoff [kg/m2]
   integer(ik), parameter :: RED_DRAINAGE    = RED_WORK + N_WORK + 4_ik   !< drainage at the column's bottom [kg/m2]
   integer(ik), parameter :: N_RED           = RED_DRAINAGE

   !----- Everything one fast step reads and writes besides cfg and the site, one per thread: the  !
   !      patch's static and dynamic state, its forcing, the step's output and its workspace.     !
   !      Each is sized once per window, to the largest patch, and every reader loops to the       !
   !      active count, never to size(). -------------------------------------------------------!
   type :: fast_step_data_t
      type(fast_static_state_t)  :: static
      type(fast_dynamic_state_t) :: dynamic
      type(fast_forcing_t)       :: forcing
      type(fast_output_t)        :: out
      !----- Allocated by ensure_step_data: as a plain component, the record's default             !
      !      initialization grows large enough that gfortran 11 warns, falsely, that the array's   !
      !      allocation reads uninitialized descriptors. -------------------------------------------!
      type(fast_workspace_t), allocatable :: work
   end type fast_step_data_t

contains

   !=======================================================================================!
   !  Advance every patch of the polygon by one slow step of fast biophysics:              !
   !  n_fast_per_slow sub-steps on the forcing source's samples, or on one fixed record    !
   !  without one. The optional outputs report the worst whole-column budget residuals and !
   !  the failure count, for a caller that asserts conservation. A start's first window     !
   !  (poly%seed_canopy_air) puts the canopy air at the window's first sample and clears    !
   !  the flag.                                                                            !
   !=======================================================================================!
   subroutine fast_dynamics(poly, cfg, met_src, met_cur, step_start, worst_energy, worst_water,  &
                            n_budget_fail, out_bufs, run_energy_budget, run_water_budget,      &
                            run_face_budget, patch_threads)
      type(polygon_t),      intent(inout) :: poly
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t), optional, intent(in)    :: met_src       !< the run's forcing source (read-only here)
      type(met_cursor_t), optional, intent(inout) :: met_cur       !< this polygon's cursor (per-sub-step met)
      type(meds_time_t),  optional, intent(in)    :: step_start    !< calendar time at the start of this slow step
      real(wp),    optional, intent(out)  :: worst_energy, worst_water
      integer(ik), optional, intent(out)  :: n_budget_fail
      type(output_buffers_t), optional, intent(inout) :: out_bufs  !< this polygon's FAST-tier staging (when on)
      !----- Run-level whole-column ledgers, which take the step's area-weighted ledgers: their signed  !
      !      cumulative residual shows a one-signed bias below the per-step tolerance. ----------------!
      type(budget_t), optional, intent(inout) :: run_energy_budget, run_water_budget
      !----- The per-layer face-closure ledger, folded the same way: it sees a vertical misplacement. -!
      type(budget_t), optional, intent(inout) :: run_face_budget
      !----- Threads for the patch loop: [run].n_threads in a one-polygon run, else 1 (the default). --!
      integer(ik),         optional, intent(in) :: patch_threads

      !----- The window's samples, each site's after its terrain, and their times: the same for every  !
      !      patch, made before the patch loop. Without forcing, every one is the fixed record. -------!
      type(met_forcing_t), allocatable :: met_sample(:)     !< (sub-step)
      type(met_forcing_t), allocatable :: met_site(:,:)     !< (sub-step, site)
      type(meds_time_t),   allocatable :: t_sample(:)       !< (sub-step)
      type(met_forcing_t) :: met_ref                        !< the air every sub-step reads without forcing
      !----- Reduction staging: one cell per (slot, sub-step, patch), assigned once in the patch loop  !
      !      (so never zeroed) and folded into the polygon in patch order after it. ----------------!
      real(wp),       allocatable :: red_site(:,:,:)                           !< (N_RED, sub-step, patch)
      real(wp),       allocatable :: red_worst_energy(:), red_worst_water(:)   !< (patch); max-folded
      type(budget_t), allocatable :: red_budget_energy(:), red_budget_water(:) !< (patch); area-merged
      type(budget_t), allocatable :: red_budget_face(:)                        !< (patch); area-merged
      integer(ik),    allocatable :: red_nfail(:)                              !< (patch) budget-failure counts
      real(wp),       allocatable :: rain_sub(:), snow_sub(:)                  !< (sub-step) the patches' split
      type(budget_t) :: site_energy_budget, site_water_budget, site_face_budget
      !=========================================================================================!
      !  The step data: one per thread, indexed by thread number and aliased with `associate`  !
      !  in the loop. An array, not OpenMP private copies or BLOCK locals, which ifx builds    !
      !  through shared compiler-made storage (CLAUDE.md, compiler traps; docs/building.md).   !
      !=========================================================================================!
      type(fast_step_data_t), allocatable :: step_data(:)
      real(wp)    :: area
      integer(ik) :: ip, isub, npatch, ncoh_max, nsub, n_thread, ith, is_p, k, iv
      integer(ik), allocatable :: visit(:)          !< the patch loop's order, costliest first
      real(wp),    allocatable :: neg_cost(:)
      integer(lk) :: tick0, tick1, tick_rate
      logical     :: do_forcing, do_fast

      !----- Live forcing drives the sub-steps when it is on and a source, cursor and step time are   !
      !      given; otherwise every sub-step reads met_forcing_t's defaults, the reference climate,   !
      !      at the run's CO2 (docs/science/forcing.md §2, §12). -------------------------------------!
      do_forcing = cfg%forcing%forcing_on .and. present(met_src) .and. present(met_cur)          &
                   .and. present(step_start)
      met_ref = met_forcing_t()
      met_ref%co2 = cfg%forcing%co2_const
      if (cfg%forcing%co2_source == CO2_SOURCE_FILE .and. present(met_src) .and. present(step_start))    &
         met_ref%co2 = met_co2(met_src, step_start)
      npatch   = poly%patch%n
      nsub     = cfg%n_fast_per_slow
      n_thread = 1_ik
      if (present(patch_threads)) n_thread = max(1_ik, patch_threads)

      call reset_fast_to_slow(poly, cfg)
      !----- FAST-tier staging, filled when the F tier or a snapshot of a fast-loop value reads it. --!
      do_fast = present(out_bufs)
      if (do_fast) do_fast = out_bufs%fast_on
      call size_staging(poly, nsub, cfg%dt_fast, do_fast, step_start, out_bufs)

      !----- The window's forcing. ----------------------------------------------------------------!
      allocate(met_sample(nsub), t_sample(nsub), met_site(nsub, poly%site%n))
      if (do_forcing) then
         call sample_window(met_src, met_cur, poly, cfg, step_start, met_sample, t_sample)
         call site_forcing(poly, cfg, met_cur, met_sample, t_sample, met_site)
      else
         met_sample = met_ref
         met_site   = met_ref
      end if
      !----- A start's canopy air is the air of the window's first sample (docs/science/            !
      !      canopy_air_space_biophysics.md, the start). -------------------------------------------!
      if (poly%seed_canopy_air) then
         call seed_start_canopy_air(poly, met_site(1, :))
         poly%seed_canopy_air = .false.
      end if

      !----- The per-(sub-step, patch) reduction staging; every cell is assigned, so none is zeroed. -!
      allocate(red_site(N_RED, nsub, max(npatch,1_ik)))
      allocate(red_budget_energy(max(npatch,1_ik)), red_budget_water(max(npatch,1_ik)))
      allocate(red_budget_face(max(npatch,1_ik)))
      allocate(red_worst_energy(max(npatch,1_ik)), red_worst_water(max(npatch,1_ik)),               &
               red_nfail(max(npatch,1_ik)))

      !----- The step data, sized to the largest patch. -------------------------------------------!
      ncoh_max = 0_ik
      do ip = 1_ik, npatch
         ncoh_max = max(ncoh_max, poly%patch%cohort_count(ip))
      end do
      allocate(step_data(n_thread))
      do ith = 1_ik, n_thread
         call ensure_step_data(step_data(ith), ncoh_max, poly%cohort%diag%active)
      end do

      !----- The costliest patches first, by each one's wall time in the last window, so the last     !
      !      patches the threads pick up are short ones. The order decides only which thread runs a  !
      !      patch and when, so every result is the same. Before any measurement (a start or a       !
      !      restart) the order is the patches' own. ----------------------------------------------!
      allocate(visit(npatch), neg_cost(npatch))
      visit = [(iv, iv = 1_ik, npatch)]
      if (npatch > 1_ik) then
         if (maxval(poly%patch%fast_cost(1:npatch)) > 0.0_wp) then
            neg_cost = -poly%patch%fast_cost(1:npatch)
            call ascending_order(neg_cost, npatch, visit)
         end if
      end if

      !=========================================================================================!
      !  The parallel patch loop (docs/science/numerical_scheme.md §6a): patches are            !
      !  independent within a dt_fast. schedule(dynamic, 1), because the adaptive march makes   !
      !  their costs very unequal. Without OpenMP the directives are comments and there is one  !
      !  step data.                                                                             !
      !=========================================================================================!
      !$omp parallel do default(shared) schedule(dynamic, 1) num_threads(n_thread)                  &
      !$omp    private(iv, ip, ith, isub, is_p, area, k, tick0, tick1, tick_rate)
      do iv = 1_ik, npatch
         ip = visit(iv)
         call system_clock(tick0, tick_rate)
         !----- This thread's step data; without OpenMP the `!$` line is a comment: slot 1. --------!
         ith = 1_ik
         !$ ith = int(omp_get_thread_num(), ik) + 1_ik
         associate (sd => step_data(ith))
         is_p = poly%patch%owner_site(ip)
         area = patch_polygon_area(poly, ip)
         call gather_patch(poly, cfg, ip, sd%static, sd%dynamic)
         sd%out%report = step_report_t()
         do isub = 1_ik, nsub
            !----- The patch's forcing, the phenology's air temperature and the precipitation split. -!
            call patch_forcing(cfg, sd%static, sd%dynamic, met_site(isub, is_p), do_forcing, sd%forcing)
            red_site(RED_PHENO_TAIR, isub, ip) = met_site(isub, is_p)%air_temp
            red_site(RED_RAINF,      isub, ip) = area * sd%forcing%rainfall
            red_site(RED_SNOWFALL,   isub, ip) = area * sd%forcing%snowfall

            call fast_step(cfg%dt_fast, cfg, poly%site%surface(is_p), sd%static, sd%forcing, sd%dynamic,  &
                           sd%out, sd%work)

            !----- The CAS -> atmosphere latent flux integrated into ET [kg/m2 = mm], the water leaving !
            !      the soil column and the step's work, area-weighted. -----------------------------!
            red_site(RED_ET,          isub, ip) = area * (sd%out%le_flux / latent_heat_vap) * cfg%dt_fast
            red_site(RED_RUNOFF_SURF, isub, ip) = area * sd%out%report%exports%runoff_surf_export
            red_site(RED_DRAINAGE,    isub, ip) = area * sd%out%report%exports%drainage_export
            do k = 1_ik, N_WORK
               red_site(RED_WORK + k, isub, ip) = area * real(sd%out%report%work(k), wp)
            end do
            call record_substep(poly, ip, isub, cfg%dt_fast, sd%static, sd%dynamic, sd%forcing, sd%out,  &
                                do_fast, out_bufs)
            call accumulate_fast_to_slow(poly, cfg, ip, sd%forcing, sd%out)
         end do
         call scatter_patch(sd%dynamic, poly, ip)

         associate (ledgers => sd%out%report%ledgers)
            red_worst_energy(ip)  = ledgers%whole_energy%worst
            red_worst_water(ip)   = ledgers%whole_water%worst
            red_budget_energy(ip) = ledgers%whole_energy
            red_budget_water(ip)  = ledgers%whole_water
            red_budget_face(ip)   = ledgers%soil_face_mass
            red_nfail(ip)         = ledgers%whole_energy%n_fail + ledgers%whole_water%n_fail
         end associate
         end associate
         call system_clock(tick1)
         poly%patch%fast_cost(ip) = real(tick1 - tick0, wp) / real(max(tick_rate, 1_lk), wp)
      end do
      !$omp end parallel do

      !=========================================================================================!
      !  Fold the staged cells into the polygon in patch order, patch outer and sub-step inner, !
      !  as a serial loop sums them, so the totals do not move with the thread count            !
      !  (docs/science/numerical_scheme.md §6a). A site's patches share its air temperature,   !
      !  so the mean over (patch, sub-step) samples is the daily mean.                          !
      !=========================================================================================!
      poly%site%pheno_tair_sum = 0.0_wp
      poly%site%pheno_tair_n   = 0_ik
      poly%et_accum            = 0.0_wp
      poly%runoff_surf_accum   = 0.0_wp
      poly%drainage_accum      = 0.0_wp
      poly%work_integ_steps    = 0.0_wp
      poly%work_integ_rej      = 0.0_wp
      poly%work_water_nsub     = 0.0_wp
      poly%work_water_newton   = 0.0_wp
      poly%work_water_nfail    = 0.0_wp
      poly%work_water_reject   = 0.0_wp
      poly%work_rk45_rescue    = 0.0_wp
      poly%work_clamp_stage    = 0.0_wp
      do ip = 1_ik, npatch
         do isub = 1_ik, nsub
            poly%et_accum          = poly%et_accum          + red_site(RED_ET,          isub, ip)
            poly%runoff_surf_accum = poly%runoff_surf_accum + red_site(RED_RUNOFF_SURF, isub, ip)
            poly%drainage_accum    = poly%drainage_accum    + red_site(RED_DRAINAGE,    isub, ip)
            associate (is_s => poly%patch%owner_site(ip))
               poly%site%pheno_tair_sum(is_s) = poly%site%pheno_tair_sum(is_s) + red_site(RED_PHENO_TAIR, isub, ip)
               poly%site%pheno_tair_n(is_s)   = poly%site%pheno_tair_n(is_s) + 1_ik
            end associate
            associate (w => red_site(RED_WORK + 1_ik:RED_WORK + N_WORK, isub, ip))
               poly%work_integ_steps  = poly%work_integ_steps  + w(WORK_INTEG_STEPS)
               poly%work_integ_rej    = poly%work_integ_rej    + w(WORK_INTEG_REJECT)
               poly%work_water_nsub   = poly%work_water_nsub   + w(WORK_WATER_NSUB)
               poly%work_water_newton = poly%work_water_newton + w(WORK_WATER_NEWTON)
               poly%work_water_reject = poly%work_water_reject + w(WORK_WATER_REJECT)
               poly%work_water_nfail  = poly%work_water_nfail  + w(WORK_WATER_NFAIL)
               poly%work_rk45_rescue  = poly%work_rk45_rescue  + w(WORK_RK45_RESCUE)
               poly%work_clamp_stage  = poly%work_clamp_stage  + w(WORK_CLAMP_STAGE)
            end associate
         end do
      end do
      !----- The precipitation split each sub-step's patches received, for the forcing echo. -------!
      allocate(rain_sub(nsub), snow_sub(nsub))
      do isub = 1_ik, nsub
         rain_sub(isub) = 0.0_wp
         snow_sub(isub) = 0.0_wp
         do ip = 1_ik, npatch
            rain_sub(isub) = rain_sub(isub) + red_site(RED_RAINF,    isub, ip)
            snow_sub(isub) = snow_sub(isub) + red_site(RED_SNOWFALL, isub, ip)
         end do
      end do
      call fold_window(poly, cfg%dt_fast, met_sample, met_site, rain_sub, snow_sub, do_fast, out_bufs)

      if (present(worst_energy)) then
         worst_energy = 0.0_wp
         do ip = 1_ik, npatch
            worst_energy = max(worst_energy, red_worst_energy(ip))
         end do
      end if
      if (present(worst_water)) then
         worst_water = 0.0_wp
         do ip = 1_ik, npatch
            worst_water = max(worst_water, red_worst_water(ip))
         end do
      end if
      if (present(n_budget_fail)) n_budget_fail = sum(red_nfail(1:npatch))
      if (present(run_energy_budget) .or. present(run_water_budget) .or. present(run_face_budget)) then
         site_energy_budget = budget_t()
         site_water_budget = budget_t()
         site_face_budget = budget_t()
         do ip = 1_ik, npatch
            call budget_merge(site_energy_budget, red_budget_energy(ip), patch_polygon_area(poly, ip))
            call budget_merge(site_water_budget,  red_budget_water(ip),  patch_polygon_area(poly, ip))
            call budget_merge(site_face_budget,   red_budget_face(ip),   patch_polygon_area(poly, ip))
         end do
         if (present(run_energy_budget)) call budget_merge(run_energy_budget, site_energy_budget, 1.0_wp)
         if (present(run_water_budget))  call budget_merge(run_water_budget,  site_water_budget,  1.0_wp)
         if (present(run_face_budget))   call budget_merge(run_face_budget,   site_face_budget,   1.0_wp)
      end if
   end subroutine fast_dynamics

   !----- Size one thread's step data for patches of up to n cohorts; capture turns the per-cohort !
   !      diagnostic capture on. --------------------------------------------------------------------!
   subroutine ensure_step_data(sd, n, capture)
      type(fast_step_data_t), intent(inout) :: sd
      integer(ik),            intent(in)    :: n
      logical,                intent(in)    :: capture
      call ensure_static_capacity(sd%static, n)
      call ensure_dynamic_capacity(sd%dynamic, n)
      call ensure_forcing_capacity(sd%forcing, n)
      call ensure_output_capacity(sd%out, n, capture)
      if (.not. allocated(sd%work)) allocate(sd%work)
   end subroutine ensure_step_data

end module meds_fast_dynamics
