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
   use ed_constants,            only : c2b, rho_h2o, cp_ice
   use ed_therm_lib,            only : internal_energy_liquid
   use ed_biogeochem_params,    only : IP_FAST_GRND
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area
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
                                       shed_excess_tissue_water, advance_soil_carbon, hold_frozen_growth
   use ed_slow_update,          only : accumulate_recruit_pool, commit_vegetation_step, sort_cohorts,    &
                                       update_overtopping_lai, age_patches, refresh_canopy_depth,        &
                                       apply_recruitment, fuse_cohorts, terminate_cohorts,               &
                                       split_cohorts, apply_patch_disturbance, sort_patches,             &
                                       fuse_patches, terminate_patches, commit_trait_carbon,             &
                                       commit_soil_carbon, follow_porous_layers
   use ed_slow_diag,            only : open_slow_cohort_rows, write_slow_patch_rows
   use ed_soil_hydr_lib,          only : soil_theta_to_head, soil_head_to_theta, soil_layer_curve_t, PSI_AIR_DRY
   use ed_ground_surface_lib,     only : porous_layer_matter, porous_layer_temp
   use ed_ground_surface_params,  only : N_POROUS, MASS_FROM_METABOLIC_LITTER, MASS_FROM_SITE_DEPTH
   implicit none
   private

   public :: advance_slow_dynamics, advance_boundary_dynamics, settle_tissue_water
   public :: vegetation_dynamics, restructure_stand, follow_porous_sources

   !--- Patch disturbance runs once a year, so it integrates its yearly rate over this interval.
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   ! [yr]
   !---~---
   !   The search for the head the top soil layer and new porous matter share (common_head): its
   ! first step down [m], its most tenfold steps (past the air-dry head, where it stops), and its
   ! most halvings. Numerical choices: the root lies far inside them.
   !---~---
   real(wp),    parameter :: HEAD_STEP_MIN    = 1.0e-3_wp
   integer(ik), parameter :: HEAD_BRACKET_MAX = 12_ik
   integer(ik), parameter :: HEAD_BISECT_MAX  = 60_ik

contains

   !---~---
   !   Advance one slow step (order_of_processes.md §4): vegetation dynamics, patch ageing,
   ! the canopy-air depth, then the daily soil-carbon step. Without `step_start` (a bare
   ! test call) phenology is not advanced.
   !---~---
   subroutine advance_slow_dynamics(poly, cfg, ledger, step_start)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_ledger_t), intent(inout) :: ledger   ! the slow ledger; an inactive one books nothing
      type(ed_time_t),   intent(in), optional :: step_start   ! the slow step's start (its day of year)
      type(soil_carbon_state_t), allocatable :: soil_pools(:)
      real(wp),                  allocatable :: heat(:), water(:)
      integer(ik) :: ip

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
      !--- The porous layers follow their matter, the water it brings drawn from the top soil layer; the dry heat declared.
      allocate(heat(poly%patch%n), water(poly%patch%n))
      call follow_porous_sources(poly, cfg, .true., heat, water)
      call slow_ledger_declare(ledger, energy_in = sum([(patch_polygon_area(poly, ip) * heat(ip),                 &
                                                         ip = 1_ik, poly%patch%n)]))
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_SOILC)
      call clear_litter(poly)
   end subroutine advance_slow_dynamics

   !---~---
   !   follow_porous_sources -- each patch's porous layers follow the matter their mass_source gives
   ! (docs/science/surface_organic_layer.md): the above-ground metabolic litter pool's dry matter
   ! (the structural pool, which holds the coarse woody debris, waits for a fine-litter split), a
   ! site's [soil].moss_depth at the material's bulk density, or none, each a layer only from
   ! min_depth (porous_layer_matter), at a run's start, at a restart and at every slow step
   ! (follow_porous_layers). The matter that arrives brings the water its own volume holds at a
   ! matric head, with the energy that water has in the top soil layer, at that layer's
   ! temperature and liquid fraction. During a run it is drawn from the top soil layer (draw), at
   ! the head the two share once it has left the layer (common_head), so a large fall of litter
   ! on dry soil takes what the soil can give, never drawing it below its air-dry head; the
   ! layer's temperature and liquid fraction are then read again from its energy. At a start it arrives
   ! from outside, at the head the top soil layer's liquid gives (on its freezing curve where it
   ! is frozen). heat and water return what each patch took from outside, which the caller
   ! declares.
   !---~---
   subroutine follow_porous_sources(poly, cfg, draw, heat, water)
      type(polygon_t),   intent(inout) :: poly
      type(ed_config_t), intent(in)    :: cfg
      logical,           intent(in)    :: draw
      real(wp),          intent(out)   :: heat(:), water(:)   ! [J/m2], [kg/m2] (patch)
      real(wp)    :: mass(N_POROUS, poly%patch%n), theta_in(N_POROUS, poly%patch%n)
      real(wp)    :: t_new(poly%patch%n), u_water(poly%patch%n), spare(poly%patch%n), dz(poly%patch%n)
      real(wp)    :: head, t1, fl1
      real(wp)    :: th_liq_top, dz_top, arriving(N_POROUS)   ! common_head's: the top soil layer, the matter arriving
      type(soil_layer_curve_t) :: c_top
      integer(ik) :: ip, k
      do ip = 1_ik, poly%patch%n
         associate (surf => poly%site%surface(poly%patch%owner_site(ip)), sc => poly%patch%soil_carbon(ip),     &
                    c => poly%site%surface(poly%patch%owner_site(ip))%soil_profile%layer_curve(1),           &
                    theta1 => poly%patch%soil_water(ip)%theta(1))
            t1  = poly%patch%soil_energy(ip)%soil_temp(1)
            fl1 = poly%patch%soil_energy(ip)%soil_fliq(1)
            head = soil_theta_to_head(c%retention, fl1 * theta1, c%theta_sat, c%theta_res, c%par_a, c%par_n)
            do k = 1_ik, N_POROUS
               associate (m => cfg%ground_surface%porous(k))
                  mass(k, ip) = 0.0_wp
                  if (m%mass_source == MASS_FROM_METABOLIC_LITTER)                                          &
                     mass(k, ip) = porous_layer_matter(c2b * sc%carbon(IP_FAST_GRND), m%bulk_density,           &
                                                       cfg%ground_surface%min_depth)
                  if (m%mass_source == MASS_FROM_SITE_DEPTH)                                                &
                     mass(k, ip) = porous_layer_matter(surf%soil%moss_depth * m%bulk_density, m%bulk_density,   &
                                                       cfg%ground_surface%min_depth)
               end associate
            end do
            if (draw) then
               c_top      = c
               th_liq_top = fl1 * theta1
               dz_top     = surf%soil_profile%dz(1)
               arriving   = max(0.0_wp, mass(:, ip) - poly%patch%ground_surface(ip)%porous%matter)
               head       = common_head(head)
            end if
            do k = 1_ik, N_POROUS
               associate (m => cfg%ground_surface%porous(k))
                  theta_in(k, ip) = soil_head_to_theta(m%retention, head, m%theta_sat, m%theta_res, m%curve_par_a, &
                                                       m%curve_par_n)
               end associate
            end do
            t_new(ip)   = t1
            u_water(ip) = fl1 * internal_energy_liquid(t1) + (1.0_wp - fl1) * cp_ice * t1
            dz(ip)      = surf%soil_profile%dz(1)
            spare(ip)   = max(0.0_wp, fl1 * theta1 - soil_head_to_theta(c%retention, PSI_AIR_DRY, c%theta_sat,       &
                                                                    c%theta_res, c%par_a, c%par_n)) * rho_h2o * dz(ip)
         end associate
      end do
      call follow_porous_layers(poly%patch, cfg%ground_surface, mass, t_new, theta_in, u_water, spare, dz, draw,   &
                                heat, water)
      if (.not. draw) return
      do ip = 1_ik, poly%patch%n
         associate (prof => poly%site%surface(poly%patch%owner_site(ip))%soil_profile,                           &
                    se => poly%patch%soil_energy(ip))
            call porous_layer_temp(se%soil_energy(1), poly%patch%soil_water(ip)%theta(1), prof%dry_heat_capacity(1), &
                                   prof%layer_curve(1), se%soil_temp(1), se%soil_fliq(1))
         end associate
      end do
   contains
      !---~---
      !   [m] The matric head the top soil layer (c_top, its liquid th_liq_top over its depth dz_top)
      ! and the matter arriving (arriving [kg/m2], per layer) share once the water the matter takes
      ! has left the layer, at head h1 before: the root of what the matter holds at a head less what
      ! the layer gives down to it (excess), which falls as the head falls, so the layer is drawn
      ! down only to where the two agree, and never below its air-dry head (PSI_AIR_DRY), where the
      ! search stops and the layer gives what it holds above that (follow_porous_layers' spare).
      ! Bracketed by tenfold steps down from h1 and halved geometrically.
      !---~---
      real(wp) function common_head(h1) result(h)
         real(wp), intent(in) :: h1
         real(wp)    :: lo, hi, mid
         integer(ik) :: it
         h = h1
         if (.not. any(arriving > 0.0_wp) .or. excess(h1) <= 0.0_wp .or. h1 <= PSI_AIR_DRY) return
         hi = h1
         lo = min(h1, -HEAD_STEP_MIN)
         do it = 1_ik, HEAD_BRACKET_MAX
            if (excess(lo) <= 0.0_wp .or. lo <= PSI_AIR_DRY) exit
            hi = lo
            lo = max(10.0_wp * lo, PSI_AIR_DRY)
         end do
         if (excess(lo) > 0.0_wp) then
            h = lo
            return
         end if
         do it = 1_ik, HEAD_BISECT_MAX
            mid = -sqrt(lo * min(hi, -HEAD_STEP_MIN))
            if (excess(mid) > 0.0_wp) then
               hi = mid
            else
               lo = mid
            end if
         end do
         h = lo
      end function common_head
      !--- [kg/m2] At head hh, what the arriving matter holds less what the top soil layer gives down to hh.
      real(wp) function excess(hh)
         real(wp), intent(in) :: hh
         integer(ik) :: kk
         excess = -rho_h2o * dz_top * (th_liq_top - soil_head_to_theta(c_top%retention, hh, c_top%theta_sat,        &
                                                                       c_top%theta_res, c_top%par_a, c_top%par_n))
         do kk = 1_ik, N_POROUS
            associate (m => cfg%ground_surface%porous(kk))
               excess = excess + soil_head_to_theta(m%retention, hh, m%theta_sat, m%theta_res, m%curve_par_a,          &
                                                    m%curve_par_n) * rho_h2o * arriving(kk) / m%bulk_density
            end associate
         end do
      end function excess
   end subroutine follow_porous_sources

   !---~---
   !   Advance the vegetation one slow step: phenology, trait plasticity, carbon allocation, the
   ! vital rates, and the growth and mortality commit. The calendar's restructuring runs at the
   ! boundary, in restructure_stand (docs/science/order_of_processes.md §4).
   !---~---
   subroutine vegetation_dynamics(poly, cfg, ledger, doy)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_ledger_t), intent(inout) :: ledger   ! the slow ledger; an inactive one books nothing
      integer(ik),         intent(in), optional :: doy   ! day-of-year at the step start (drives phenology)
      type(slow_step_data_t) :: step
      real(wp)               :: mort_water, seed_rain_carbon, handover, litter_open, frozen_carbon
      real(wp), allocatable  :: leaf_carbon_new(:), storage_new(:)

      !--- The step's slow cohort rows, opened before the phenology writes its cues to them.
      call open_slow_cohort_rows(poly)

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
      ! as carbon in: no parent at the site paid for it (docs/science/demography.md §2). With
      ! the stand frozen ([demography].demography_on = false) nothing recruits, grows or dies:
      ! the wood growth and the recruits' share of the reproduction are held out of the stand
      ! as frozen growth carbon, declared out below (demography.md §7).
      !---~---
      frozen_carbon = 0.0_wp
      if (cfg%demography%demography_on) then
         call accumulate_recruit_pool(poly, cfg, step%recruitment, cfg%dt_years, seed_rain_carbon)
         call slow_ledger_declare(ledger, carbon_in = seed_rain_carbon)
      else
         call hold_frozen_growth(poly, cfg, step, frozen_carbon)
      end if

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
      !   The growth phase declares the fast loop's net carbon (in), the litter, the CO2 owed
      ! to the canopy air and a frozen stand's unbuilt growth (out) and the dying plants' water,
      ! so it closes to round-off (docs/science/conservation.md §3).
      !---~---
      if (ledger%active) then
         call slow_ledger_declare(ledger, carbon_in  = handover,                                    &
                                          carbon_out = litter_carbon_total(poly, poly%patch%litter_in)  &
                                                     - litter_open                                    &
                                                     + slow_co2_handoff(poly, cfg)                 &
                                                     + frozen_carbon,                              &
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
      type(slow_ledger_t), intent(inout) :: ledger   ! the slow ledger; an inactive one books nothing
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
   ! diameter; the cull and the treefall remove cohorts): one exchange per phase, declared
   ! (conservation.md §4).
   !---~---
   subroutine restructure_stand(poly, cfg, is_new_month, is_new_year, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year   ! which boundary this is
      type(slow_ledger_t), intent(inout) :: ledger   ! the slow ledger; an inactive one books nothing
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
      type(slow_ledger_t), intent(inout) :: ledger   ! the slow ledger; an inactive one books nothing
      real(wp) :: water_shed
      if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
      call slow_ledger_heat_open(ledger, poly)
      call shed_excess_tissue_water(poly, cfg, water_shed)
      call slow_ledger_declare(ledger, water_out = water_shed)
      call slow_ledger_heat_close(ledger, poly)
      call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_TISSUE_WATER)
   end subroutine settle_tissue_water

end module ed_slow_dynamics
