! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_polygon_state_vars -- the polygon's ecological state: its site, patch and cohort      !
! blocks, the polygon-level numbers, and the routines that need more than one block (the   !
! patch and cohort maps, the persistent ids, the step diagnostic reset).                   !
!                                                                                          !
! The run keeps beside each polygon_t a polygon_run_t (src/main) with its forcing, budgets !
! and ledgers.                                                                             !
!==========================================================================================!
module ed_polygon_state_vars
   use ed_diag_vars,          only :   polygon_diag_block, reset_cohort_diag, reset_patch_diag,          &
                                       reset_polygon_diag, reset_site_diag
   use ed_kinds,                only : wp, ik
   use ed_pft_params,           only : pft_table_t
   use ed_pft_derive,           only : derive_recruit_carbon
   use ed_cohort_state_vars,    only : cohort_block, cohort_deriv_block, alloc_cohort_block, free_cohort_block, reorder_cohort_block
   use ed_patch_state_vars,     only : patch_block, alloc_patch_block, free_patch_block
   use ed_site_state_vars,      only : site_block, alloc_site_block, free_site_block
   implicit none
   private

   public :: polygon_t, alloc_polygon, free_polygon, rebuild_patch_map, rebuild_cohort_map
   public :: patch_polygon_area, patch_polygon_areas
   public :: assign_cohort_id, assign_patch_id, reset_step_diagnostics

   type :: polygon_t
      !==========================================================================================!
      ! 0. THE BLOCKS IT HOLDS                                                                   !
      !==========================================================================================!
      type(site_block)   :: site
      type(patch_block)  :: patch
      type(cohort_block) :: cohort
      !==========================================================================================!
      ! 1. IDENTITY AND LOCATION                                                                 !
      !==========================================================================================!
      integer(ik)        :: id = 0_ik       !< global id: the cell's row-major index (0 for a site run)
      !----- Where it is: the cell centre and orography of a region's polygon, the [site] of a site   !
      !      run. The forcing cursor, the leaf phenology's day length and the output all read it. ---!
      real(wp)           :: latitude_deg = 0.0_wp, longitude_deg = 0.0_wp, elevation_m = 0.0_wp
      !==========================================================================================!
      ! 2. DIMENSIONS AND ID COUNTERS                                                            !
      !==========================================================================================!
      integer(ik)        :: n_pft     = 0_ik
      !----- Monotonic counters for the persistent global ids (never reused within a run). ---!
      integer(ik)        :: next_cohort_id = 1_ik
      integer(ik)        :: next_patch_id  = 1_ik
      !==========================================================================================!
      ! 3. CALENDAR: the cohorts' shared sample clock and the restructuring owed                 !
      !==========================================================================================!
      !----- Global ring-buffer write slot for the moving-average growth (same for all cohorts !
      !      since they all take one sample per step); advanced by commit_vegetation_step.      !
      integer(ik)        :: growth_hist_pos = 0_ik
      !----- The calendar restructuring the last step's end owes: set when a step ends on a month    !
      !      (and, with restructure_new_year, a year) boundary, run at the start of the next step,   !
      !      before its fast loop. Between the two sit the output tick, the I/O phase and any        !
      !      checkpoint, so a checkpoint holds the stand before the restructuring together with      !
      !      this flag, and a run resumed from it restructures first, exactly as the continuous run  !
      !      does. ----------------------------------------------------------------------------------!
      logical            :: restructure_pending  = .false.
      logical            :: restructure_new_year = .false.
      !----- A bare-ground or census start owes its canopy air the first forcing sample: the fast    !
      !      loop's first window seeds it and clears this flag. A restart's state file holds the     !
      !      canopy air, so a restart owes nothing, and no state file is written before that window. !
      logical            :: seed_canopy_air = .false.
      !==========================================================================================!
      ! 4. THE FORCING'S MEMORY                                                                  !
      !==========================================================================================!
      !----- The last daytime clearness index, which the longwave synthesis holds through the night  !
      !      (docs/science/forcing.md §11). The forcing cursor works on a copy during each fast loop  !
      !      (fast_dynamics), so a state file carries it and a restart resumes the same night. ---------!
      real(wp)           :: kt_last_day = 1.0_wp
      !==========================================================================================!
      ! 5. TENDENCIES (scratch)                                                                     !
      !==========================================================================================!
      !----- TRANSIENT slow-loop scratch: the per-cohort tendency bundle the driver fills each     !
      !      step and update_cohort_states applies. Held here so it allocates once (resize-on-      !
      !      grow) instead of per step; refilled every step, so it is NOT lockstep-reordered and    !
      !      NOT serialized (a scratch buffer, not prognostic state).                               !
      type(cohort_deriv_block) :: deriv
      !==========================================================================================!
      ! 6. OUTPUT ACCUMULATORS (only output reads them; never restarted)                         !
      !==========================================================================================!
      !----- TRANSIENT fast-loop diagnostics that are the same across the polygon -- its forcing   !
      !      (PY_* in ed_diag_vars): one slot, the cohort and patch blocks' polygon twin. ---        !
      type(polygon_diag_block) :: diag
      !----- The polygon's evapotranspiration accumulator [kg/m2 = mm]: the fast loop sums the       !
      !      area-weighted canopy-air -> atmosphere water-vapour flux * dt_fast over the              !
      !      slow step (reset each step, mirrors the gpp_accum/pheno lifecycle); read as a diagnostic   !
      !      (SRC_PY_ET) at the output tick. Transient (never restarted).                                !
      real(wp)           :: et_accum       = 0.0_wp
      !----- And the water the soil column lost the same way [kg/m2 = mm]: the pond's spill (surface   !
      !      runoff) and the bottom flux (drainage; negative when an aquifer feeds the column). --------!
      real(wp)           :: runoff_surf_accum = 0.0_wp
      real(wp)           :: drainage_accum    = 0.0_wp
      !----- Integrator WORK accumulators: summed over patches and fast sub-steps within a slow      !
      !      step, reset alongside et_accum. These are the COST axis of a scheme comparison --       !
      !      wall-clock is too coarse and too machine-dependent to rank integrators, whereas         !
      !      step/sub-solver counts are exact and reproducible. Transient (diagnostic only,          !
      !      never restarted). Real-valued so the output integrators can average them over a         !
      !      period like any other diagnostic. ------------------------------------------------------!
      real(wp)           :: work_integ_steps = 0.0_wp   !< accepted integrator sub-steps
      real(wp)           :: work_integ_rej   = 0.0_wp   !< rejected integrator steps
      !----- The soil-plant water network's work. A nonzero work_water_nfail is also a correctness     !
      !      flag: a sub-step taken at the shortest length although Newton had not converged keeps      !
      !      water exactly but leaves the stores off their curves. -------------------------------------!
      real(wp)           :: work_water_nsub   = 0.0_wp   !< accepted sub-steps
      real(wp)           :: work_water_newton = 0.0_wp   !< Newton iterations
      real(wp)           :: work_water_reject = 0.0_wp   !< rejected sub-steps
      real(wp)           :: work_water_nfail  = 0.0_wp   !< sub-steps taken although unconverged
      !----- INTEGRATOR-HEALTH accumulators. Same cadence and area weighting as the work             !
      !      counters above, but these measure CORRECTNESS rather than cost: how often a scheme      !
      !      had to leave its own model to stay standing, and how much mass/energy that cost         !
      !      outside the conservation ledger. Without them a run cannot report whether it            !
      !      actually ran the integrator it was configured for (docs/science/numerical_scheme.md     !
      !      §3, `rk45`). ---------------------------------------------------------------------------!
      real(wp)           :: work_rk45_rescue = 0.0_wp   !< dt_fast steps RK45 abandoned and redid with the ESDIRK2
      real(wp)           :: work_clamp_stage = 0.0_wp   !< stage-input clamp activations (throwaway states)
   end type polygon_t

contains

   !=======================================================================================!
   !  LIFECYCLE                                                                            !
   !=======================================================================================!

   !=======================================================================================!
   !  Allocation                                                                            !
   !=======================================================================================!
   subroutine alloc_polygon(poly, n_pft, coh_cap, pat_cap, n_growth_window, pft, n_site)
      type(polygon_t),   intent(out) :: poly
      integer(ik),       intent(in)  :: n_pft, coh_cap, pat_cap, n_growth_window
      type(pft_table_t), intent(in)  :: pft     !< the polygon's PFT traits and allometry (cfg%pft of its config)
      integer(ik), optional, intent(in) :: n_site   !< its sites (1 when absent)
      integer(ik) :: ns
      poly%n_pft = n_pft
      poly%cohort%pft_table = pft
      call derive_recruit_carbon(poly%cohort%pft_table)
      !----- Its sites, sharing the polygon equally until the caller gives their areas; every   !
      !      patch starts in the first. --------------------------------------------------------!
      ns = 1_ik
      if (present(n_site)) ns = max(n_site, 1_ik)
      call alloc_site_block(poly%site, ns, pft%leaf_photo)
      poly%site%n          = ns
      poly%site%area(1:ns) = 1.0_wp / real(ns, wp)
      call alloc_patch_block(poly%patch, max(pat_cap, 1_ik), n_pft)
      call alloc_cohort_block(poly%cohort, max(coh_cap, 1_ik), max(n_growth_window, 1_ik))
   end subroutine alloc_polygon

   !----- Each level frees its own fields: the site, patch and cohort files list them. -------!
   subroutine free_polygon(poly)
      type(polygon_t), intent(inout) :: poly
      call free_site_block(poly%site)
      call free_patch_block(poly%patch)
      call free_cohort_block(poly%cohort)
   end subroutine free_polygon

   !=======================================================================================!
   !  STRUCTURE                                                                            !
   !=======================================================================================!

   !=======================================================================================!
   !  Rebuild the site -> patch map from owner_site. The patch block is sorted by site      !
   !  (sort_patches), so each site's patches are one slice; a block out of that order stops. !
   !=======================================================================================!
   subroutine rebuild_patch_map(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik) :: ip, is
      associate (site => poly%site, patch => poly%patch)
         site%patch_count(1:site%n) = 0_ik
         do ip = 1_ik, patch%n
            is = patch%owner_site(ip)
            if (is < 1_ik .or. is > site%n) error stop 'rebuild_patch_map: owner_site out of range'
            if (ip > 1_ik) then
               if (is < patch%owner_site(ip - 1_ik)) error stop 'rebuild_patch_map: patches not sorted by site'
            end if
            site%patch_count(is) = site%patch_count(is) + 1_ik
         end do
         site%patch_offset(1) = 1_ik
         do is = 2_ik, site%n
            site%patch_offset(is) = site%patch_offset(is - 1_ik) + site%patch_count(is - 1_ik)
         end do
      end associate
   end subroutine rebuild_patch_map

   !=======================================================================================!
   !  Rebuild the CSR patch map by stably regrouping cohorts by owner_patch (counting sort).!
   !=======================================================================================!
   subroutine rebuild_cohort_map(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik), allocatable :: perm(:), pos(:), slot(:)
      integer(ik)              :: i, ip, np, nc
      np = poly%patch%n
      nc = poly%cohort%n
      if (np < 1_ik) then                                  ! no patches: nothing to group
         if (nc > 0_ik) error stop 'rebuild_cohort_map: cohorts present with no patches'
         return
      end if
      poly%patch%cohort_count(1:np) = 0_ik
      do i = 1_ik, nc
         ip = poly%cohort%owner_patch(i)
         if (ip < 1_ik .or. ip > np) error stop 'rebuild_cohort_map: owner_patch out of range'
         poly%patch%cohort_count(ip) = poly%patch%cohort_count(ip) + 1_ik
      end do
      poly%patch%cohort_offset(1) = 1_ik
      do ip = 2_ik, np
         poly%patch%cohort_offset(ip) = poly%patch%cohort_offset(ip-1) + poly%patch%cohort_count(ip-1)
      end do
      !----- Stable placement: preserve within-patch order. -------------------------------!
      allocate(pos(np), slot(nc), perm(nc))
      pos(1:np) = poly%patch%cohort_offset(1:np)
      do i = 1_ik, nc
         ip = poly%cohort%owner_patch(i)
         slot(i) = pos(ip)
         pos(ip) = pos(ip) + 1_ik
      end do
      do i = 1_ik, nc
         perm(slot(i)) = i
      end do
      call reorder_cohort_block(poly%cohort, perm, nc)
   end subroutine rebuild_cohort_map

   !=======================================================================================!
   !  AREA                                                                                 !
   !=======================================================================================!

   !----- A patch's share of the polygon: its area is a fraction of its site, and its site's   !
   !      a fraction of the polygon. Every polygon total weights its patches by this; anything !
   !      within one site (fusion, a disturbance gap, the patch thresholds) by the patch area. !
   pure real(wp) function patch_polygon_area(poly, ip) result(a)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: ip
      a = poly%site%area(poly%patch%owner_site(ip)) * poly%patch%area(ip)
   end function patch_polygon_area

   !----- Every patch's share, in patch order: a contiguous weight for a loop over the patches. !
   pure subroutine patch_polygon_areas(poly, share)
      type(polygon_t), intent(in)  :: poly
      real(wp),        intent(out) :: share(:)
      integer(ik) :: ip
      do ip = 1_ik, poly%patch%n
         share(ip) = patch_polygon_area(poly, ip)
      end do
   end subroutine patch_polygon_areas

   !=======================================================================================!
   !  IDENTITY                                                                             !
   !=======================================================================================!

   !=======================================================================================!
   !  Persistent identity: stamp a slot with the next unused global id and bump the         !
   !  monotonic counter. Ids are never reused within a run, so a terminated/fused-away id    !
   !  simply never reappears -- exactly what an external tracker needs.                      !
   !=======================================================================================!
   subroutine assign_cohort_id(poly, slot)
      type(polygon_t), intent(inout) :: poly
      integer(ik),  intent(in)    :: slot
      poly%cohort%global_id(slot) = poly%next_cohort_id
      poly%next_cohort_id = poly%next_cohort_id + 1_ik
   end subroutine assign_cohort_id

   subroutine assign_patch_id(poly, slot)
      type(polygon_t), intent(inout) :: poly
      integer(ik),  intent(in)    :: slot
      poly%patch%global_id(slot) = poly%next_patch_id
      poly%next_patch_id = poly%next_patch_id + 1_ik
   end subroutine assign_patch_id

   !=======================================================================================!
   !  PER STEP                                                                             !
   !=======================================================================================!

   !----- Zero the step diagnostic blocks -- per cohort, per patch, per polygon -- after the       !
   !      output has read them. They accumulate from one output tick to the next: the fast loop    !
   !      adds each sub-step, the slow step its rows, and a calendar boundary between two steps    !
   !      its events (disturbance, culls), which therefore reach the record of the step that       !
   !      follows. The slow cohort block is SET each step instead and needs no reset here. --------!
   subroutine reset_step_diagnostics(poly)
      type(polygon_t), intent(inout) :: poly
      call reset_cohort_diag(poly%cohort%diag)
      call reset_patch_diag(poly%patch%diag)
      call reset_polygon_diag(poly%diag)
      call reset_site_diag(poly%site%diag)
   end subroutine reset_step_diagnostics

end module ed_polygon_state_vars
