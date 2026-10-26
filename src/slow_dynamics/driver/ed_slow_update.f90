! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_slow_update -- the one place a slow step writes the stand's carbon and structure      !
! (docs/science/demography.md §2-§6). The laws of ed_demography decide; this module        !
! applies what it is given and never uses them (src/README.md rule 4). In four sections:   !
!   1. The step's commit: each cohort's tendencies and growth history, its pools (eq. (8)) !
!      and density, the recruit pool's credit, patch ageing, the overtopping LAI           !
!      (eq. (11)) and the canopy-air depth.                                                !
!   2. Recruitment: apply_recruitment founds a cohort from each full recruit pool (§3).    !
!   3. Disturbance: apply_patch_disturbance carves each site's treefall gap (§4).          !
!   4. Resolution: sort, fuse, split and cull the cohorts (§3); sort, fuse and terminate   !
!      the patches on their light profiles (§4, eq. (10)).                                 !
! What each operator conserves and hands to the slow ledger is demography.md §5.           !
!==========================================================================================!
module ed_slow_update
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : tiny_num, lnexp_min, lnexp_max, almost_one, yr_sec, rho_air_ref, umol_2_kgc
   use ed_numerics,             only : equal_height_layers, conserved, pair_mean
   use ed_cohort_state_vars,    only : cohort_block, cohort_deriv_block, set_cohort_size_from_carbon,    &
                                       hold_cohort_tissue_temps, reorder_cohort_block,                   &
                                       compact_cohort_block, expand_cohort_block, copy_cohort,           &
                                       init_cohort, cohort_tissue_heat_capacity, cohort_film_water,      &
                                       cohort_water_enthalpy, cohort_tissue_water,                       &
                                       fuse_cohort_fast_state, set_cohort_tissue_temps,                  &
                                       cohort_deriv_alloc
   use ed_patch_state_vars,     only : patch_block, soil_carbon_state_t, blend_patch_stores,             &
                                       add_shed_water, patch_physical_state, DIST_TREEFALL,              &
                                       expand_patch_block, reorder_patch_block, cas_set_depth
   use ed_polygon_state_vars,   only : polygon_t, rebuild_cohort_map, assign_cohort_id,                  &
                                       patch_polygon_area, rebuild_patch_map, assign_patch_id
   use ed_allometry_lib,        only : height_to_dbh, carbon_to_structure
   use ed_config,               only : ed_config_t, growth_window_steps
   use ed_diag_vars,            only : fuse_cohort_diag, PA_MORT_C_CULL, blend_patch_diag,               &
                                       inherit_patch_diag, PA_DISTURB_AREA, PA_MORT_C_DISTURB,           &
                                       CO_DDBH_DT, CO_DAGB_DT, CO_MORT_RATE
   use ed_slow_ledger,          only : slow_ledger_t, slow_ledger_declare
   use ed_slow_state,           only : slow_step_data_t, kill_report_t, open_kill_report, report_kill
   implicit none
   private

   public :: commit_vegetation_step, update_cohort_derivatives, fill_cohort_deriv, update_cohort_states
   public :: advance_pool, accumulate_recruit_pool, age_patches, update_overtopping_lai
   public :: commit_trait_carbon, commit_soil_carbon
   public :: refresh_canopy_depth, apply_recruitment, apply_patch_disturbance
   public :: sort_cohorts, max_cohort_count, fuse_cohorts, fuse_2_cohorts, split_cohorts, terminate_cohorts
   public :: sort_patches, patch_light_profile, fuse_patches, fuse_2_patches, terminate_patches

   !----- How far past its bound a quantity may land from round-off alone, in epsilons of the magnitudes !
   !      it touches (demography.md §2). An emptied pool, or a density decline aimed at the floor, gets  !
   !      there through four roundings of half an epsilon each, so 4 leaves a factor of two to spare. ---!
   real(wp), parameter :: ROUNDOFF_EPSILONS = 4.0_wp

contains

   !=======================================================================================!
   ! 1. The step's commit: tendencies, pools and density, recruit pool, ageing, canopy     !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Commit one step's growth and mortality (docs/science/demography.md §2): take the growth !
   ! history's slot for this step's sample, build each cohort's tendencies, keep the density !
   ! and leaf carbon the commit starts from in the step record, then advance the pools and   !
   ! the density.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine commit_vegetation_step(poly, cfg, step)
      type(polygon_t),        intent(inout) :: poly
      type(ed_config_t),    intent(in)      :: cfg
      type(slow_step_data_t), intent(inout) :: step   !< reads the pools' changes and mortality; fills the *_before
      integer(ik) :: n_window

      n_window = growth_window_steps(cfg)
      poly%growth_hist_pos = mod(poly%growth_hist_pos, n_window) + 1_ik
      call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
      call update_cohort_derivatives(poly, step, cfg%dt_years, n_window, poly%growth_hist_pos)
      step%nplant_before      = poly%cohort%nplant(1:poly%cohort%n)
      step%leaf_carbon_before = poly%cohort%leaf_carbon(1:poly%cohort%n)
      call update_cohort_states(poly%cohort, poly%deriv, cfg%dt_years, cfg%demography%negligible_nplant)
   end subroutine commit_vegetation_step

   !---------------------------------------------------------------------------------------!
   ! Each cohort's tendencies for the applier, which commits them: the four pools' rates of     !
   ! change, the log-space mortality rate, and the realized diameter rate, entered in the     !
   ! growth history that mortality averages. Nothing is committed here; the diameter comes from  !
   ! the allometry's inverse at the tentative pools (docs/science/demography.md §1-§2).          !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_derivatives(poly, step, dt_yr, n_window, hist_pos)
      type(polygon_t),          intent(inout) :: poly       !< fills poly%deriv and the growth history
      type(slow_step_data_t),   intent(in)    :: step       !< reads the pools' changes and mortality
      real(wp),                 intent(in)    :: dt_yr
      integer(ik),              intent(in)    :: n_window, hist_pos
      integer(ik) :: i
      real(wp)    :: lc_new, fc_new, wc_new, nc_new
      real(wp)    :: dbh_new, height_new, ba_new, agb_new, la_new, dbh_rate

      associate (cohort => poly%cohort, mortality => step%mortality)
         do i = 1_ik, cohort%n
            !----- Tentative pools by the applier's own rule (advance_pool): an overdraw stays       !
            !      negative, for the bounds guard to name (docs/science/demography.md eq. (8)). -----!
            lc_new = advance_pool(cohort%leaf_carbon(i),          step%d_leaf(i))
            fc_new = advance_pool(cohort%fineroot_carbon(i),      step%d_fineroot(i))
            wc_new = advance_pool(cohort%wood_carbon(i),          step%d_wood(i))
            nc_new = advance_pool(cohort%nonstructural_carbon(i), step%d_nonstructural(i))
            !----- The tentative geometry, for the reported growth rates; update_cohort_states     !
            !      derives the same geometry from the committed pools. -----------------------------!
            call carbon_to_structure(cohort%pft_table%allom(cohort%pft(i)), wc_new, lc_new,          &
                                     cohort%pft_table%wood_density(cohort%pft(i)), cohort%sla(i),   &
                                     dbh_new, height_new, ba_new, agb_new, la_new)
            !----- The pools' tendencies, and the realized diameter rate into the growth history. -!
            dbh_rate = (dbh_new - cohort%dbh(i)) / dt_yr
            call fill_cohort_deriv(cohort, i, poly%deriv, dt_yr, mortality(i), dbh_rate,            &
                                   lc_new, fc_new, wc_new, nc_new, n_window, hist_pos)
            !----- Copy the reported rates into the slow cohort block, which follows the cohort     !
            !      reordering; poly%deriv does not (docs/science/diagnostics.md §8, "The one trap"). !
            if (cohort%sdiag%active) then
               cohort%sdiag%v(CO_DDBH_DT,   i) = dbh_rate
               cohort%sdiag%v(CO_DAGB_DT,   i) = (agb_new - cohort%agb(i)) / dt_yr
               cohort%sdiag%v(CO_MORT_RATE, i) = mortality(i)
               cohort%sdiag%w(i)               = 1.0_wp
            end if
         end do
      end associate
   end subroutine update_cohort_derivatives

   !---------------------------------------------------------------------------------------!
   ! Cohort i's tendency from the new carbon pools the caller computed (carbon allocation  !
   ! or empirical growth law): each pool's rate is (new - old)/dt, the density's is        !
   ! -mortality in log space. It also adds ring_sample to the growth history (eq. (2)).    !
   !---------------------------------------------------------------------------------------!
   subroutine fill_cohort_deriv(cohort, i, deriv, dt_yr, mortality, ring_sample,                   &
                                leaf_carbon_new, fineroot_carbon_new, wood_carbon_new,              &
                                nonstructural_carbon_new, n_window, hist_pos)
      type(cohort_block),       intent(inout) :: cohort
      integer(ik),              intent(in)    :: i
      type(cohort_deriv_block), intent(inout) :: deriv
      real(wp),                 intent(in)    :: dt_yr, mortality, ring_sample
      real(wp),                 intent(in)    :: leaf_carbon_new, fineroot_carbon_new, wood_carbon_new
      real(wp),                 intent(in)    :: nonstructural_carbon_new
      integer(ik),              intent(in)    :: n_window, hist_pos

      deriv%d_leaf_carbon_dt(i)          = (leaf_carbon_new          - cohort%leaf_carbon(i))          / dt_yr
      deriv%d_fineroot_carbon_dt(i)      = (fineroot_carbon_new      - cohort%fineroot_carbon(i))      / dt_yr
      deriv%d_wood_carbon_dt(i)          = (wood_carbon_new          - cohort%wood_carbon(i))          / dt_yr
      deriv%d_nonstructural_carbon_dt(i) = (nonstructural_carbon_new - cohort%nonstructural_carbon(i)) / dt_yr
      deriv%dln_nplant_dt(i)             = -mortality
      !----- Write the sample over the oldest and sum the window from the buffer itself, not as !
      !      a running total (demography.md §1); the slots not yet written are 0, cleared by    !
      !      init_cohort. ----------------------------------------------------------------------!
      if (cohort%growth_count(i) < n_window) cohort%growth_count(i) = cohort%growth_count(i) + 1_ik
      cohort%growth_hist(hist_pos, i) = ring_sample
      cohort%growth_accum(i) = sum(cohort%growth_hist(1:n_window, i))
      cohort%growth_avg(i) = cohort%growth_accum(i) / real(cohort%growth_count(i), wp)
   end subroutine fill_cohort_deriv

   !---------------------------------------------------------------------------------------!
   ! Advance every cohort's pools (advance_pool, eq. (8)) and density (advance_density, in !
   ! log space, stopped at the floor `negligible`) by the tendencies over dt_yr, then      !
   ! derive the geometry from the committed carbon, as fusion and fission derive it. The   !
   ! tissue made or lost came or went at the tissue's temperature, which its energy is     !
   ! held to (demography.md §2; conservation.md §4).                                       !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_states(cohort, deriv, dt_yr, negligible)
      type(cohort_block),       intent(inout) :: cohort
      type(cohort_deriv_block), intent(in)    :: deriv
      real(wp),                 intent(in)    :: dt_yr, negligible
      integer(ik) :: i
      do i = 1_ik, cohort%n
         cohort%leaf_carbon(i)          = advance_pool(cohort%leaf_carbon(i),          deriv%d_leaf_carbon_dt(i) * dt_yr)
         cohort%fineroot_carbon(i)      = advance_pool(cohort%fineroot_carbon(i),      deriv%d_fineroot_carbon_dt(i) * dt_yr)
         cohort%wood_carbon(i)          = advance_pool(cohort%wood_carbon(i),          deriv%d_wood_carbon_dt(i) * dt_yr)
         cohort%nonstructural_carbon(i) = advance_pool(cohort%nonstructural_carbon(i),                    &
                                                       deriv%d_nonstructural_carbon_dt(i) * dt_yr)
         cohort%nplant(i)               = advance_density(cohort%nplant(i), deriv%dln_nplant_dt(i) * dt_yr, negligible)
         call set_cohort_size_from_carbon(cohort, i)
         call hold_cohort_tissue_temps(cohort, i)
      end do
   end subroutine update_cohort_states

   !---------------------------------------------------------------------------------------!
   ! The density after the log-space change dln, stopped at the floor (demography.md §2).  !
   ! A decline that reaches the floor, to within the round-off of comparing the two in log !
   ! space, ends exactly on it, so the cull finds it there; a cohort already at or below   !
   ! the floor is not taken lower. The `1 +` counts the quotient's own rounding.           !
   !---------------------------------------------------------------------------------------!
   elemental pure real(wp) function advance_density(nplant, dln, floor) result(nplant_new)
      real(wp), intent(in) :: nplant   !< [plant/m2] the density
      real(wp), intent(in) :: dln      !< [-] its change in log space over the step, -mortality*dt
      real(wp), intent(in) :: floor    !< [plant/m2] the density floor, negligible_nplant
      real(wp) :: floor_dln
      floor_dln = min(0.0_wp, log(max(floor, tiny_num) / max(nplant, tiny_num)))
      if (floor_dln < 0.0_wp .and. dln - floor_dln <= ROUNDOFF_EPSILONS * epsilon(dln)              &
                                                     * (1.0_wp + abs(floor_dln) + abs(dln))) then
         nplant_new = floor
      else
         nplant_new = nplant * exp(min(max(max(dln, floor_dln), lnexp_min), lnexp_max))
      end if
   end function advance_density

   !---------------------------------------------------------------------------------------!
   ! A carbon pool advanced by its change over the step, eq. (8). A result below zero by   !
   ! no more than ROUNDOFF_EPSILONS of |pool| + |change| is round-off and set to 0. An     !
   ! overdraw further below is left as it is, so the bounds check (cohort_bad_field) stops !
   ! the run and names the pool (demography.md §2).                                        !
   !---------------------------------------------------------------------------------------!
   elemental pure real(wp) function advance_pool(pool, change) result(pool_new)
      real(wp), intent(in) :: pool, change    !< [kgC/plant] the pool, and its change over the step
      pool_new = pool + change
      if (pool_new < 0.0_wp .and.                                                                  &
          -pool_new <= ROUNDOFF_EPSILONS * epsilon(pool_new) * (abs(pool) + abs(change))) pool_new = 0.0_wp
   end function advance_pool

   !---------------------------------------------------------------------------------------!
   ! Credit one step's recruitment into each (PFT, patch) recruit pool, and report the seed     !
   ! rain, the part no parent at the site paid for. The parents' share is debited and credited  !
   ! in the same phase, so it is internal. include_pft gates the seed rain as it gates the rate  !
   ! (docs/science/demography.md §2; conservation.md §3).                                        !
   !---------------------------------------------------------------------------------------!
   subroutine accumulate_recruit_pool(poly, cfg, recruitment, dt_yr, seed_rain_carbon)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      real(wp),            intent(in)    :: dt_yr
      real(wp),            intent(out)   :: seed_rain_carbon  !< [kgC/m2 site] the external part
      integer(ik) :: ip, pf

      seed_rain_carbon = 0.0_wp
      associate (patch => poly%patch, pft => poly%cohort%pft_table)
         do ip = 1_ik, patch%n
            do pf = 1_ik, poly%n_pft
               patch%recruit_pool(pf, ip) = patch%recruit_pool(pf, ip) + recruitment(pf, ip) * dt_yr
               if (pft%include_pft(pf) /= 1_ik) cycle
               seed_rain_carbon = seed_rain_carbon + patch_polygon_area(poly, ip)                  &
                                * pft%seed_rain_recruits(pf) * dt_yr * pft%recruit_carbon(pf)
            end do
         end do
      end associate
   end subroutine accumulate_recruit_pool

   !---------------------------------------------------------------------------------------!
   ! Age every patch by dt_yr.                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine age_patches(patch, dt_yr)
      type(patch_block), intent(inout) :: patch
      real(wp),          intent(in)    :: dt_yr
      integer(ik) :: ip
      do ip = 1_ik, patch%n
         patch%age(ip) = patch%age(ip) + dt_yr
      end do
   end subroutine age_patches

   !---------------------------------------------------------------------------------------!
   ! The soil-carbon step's new pools and lignin, per patch, from ed_slow_state's          !
   ! advance_soil_carbon (docs/science/soil_carbon.md §6).                                 !
   !---------------------------------------------------------------------------------------!
   subroutine commit_soil_carbon(patch, pools_new)
      type(patch_block),         intent(inout) :: patch
      type(soil_carbon_state_t), intent(in)    :: pools_new(:)   !< (patch)
      patch%soil_carbon(1:patch%n) = pools_new(1:patch%n)
   end subroutine commit_soil_carbon

   !---------------------------------------------------------------------------------------!
   ! The leaf and storage carbon a trait change leaves, at constant total carbon: the      !
   ! light plasticity's resorption and a restart's re-acclimation, from ed_slow_state      !
   ! (docs/science/plant_traits.md §4-§5). Leaf area follows the new SLA; the tissues      !
   ! keep their temperatures.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine commit_trait_carbon(poly, leaf_carbon_new, storage_new)
      type(polygon_t), intent(inout) :: poly
      real(wp),        intent(in)    :: leaf_carbon_new(:)   !< [kgC/plant] (cohort)
      real(wp),        intent(in)    :: storage_new(:)       !< [kgC/plant] (cohort)
      integer(ik) :: j
      do j = 1_ik, poly%cohort%n
         poly%cohort%leaf_carbon(j)          = leaf_carbon_new(j)
         poly%cohort%nonstructural_carbon(j) = storage_new(j)
         call set_cohort_size_from_carbon(poly%cohort, j)   ! leaf_area = sla*leaf_carbon
         call hold_cohort_tissue_temps(poly%cohort, j)      ! the moved leaf carbon left at the leaf's temperature
      end do
   end subroutine commit_trait_carbon

   !---------------------------------------------------------------------------------------!
   ! Each cohort's overtopping LAI, eq. (11) of demography.md §6: the LAI of the taller    !
   ! cohorts in its patch, in one top-down pass. Equal-height cohorts share a layer and do !
   ! not shade each other, as in the radiation (equal_height_layers).                      !
   !---------------------------------------------------------------------------------------!
   subroutine update_overtopping_lai(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik) :: ip, i0, i1, i, j, k, il, n_layer
      real(wp)    :: cum, over_lai, layer_lai

      associate (cohort => poly%cohort, patch => poly%patch)
         do ip = 1_ik, patch%n
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            block
               integer(ik) :: layer(i0:i1)
               call equal_height_layers(cohort%height(i0:i1), i1 - i0 + 1_ik, layer, n_layer)
               cum = 0.0_wp                                ! overtopping LAI accumulator
               i   = i0
               do il = 1_ik, n_layer
                  k = i
                  do while (k < i1)
                     if (layer(k + 1_ik) /= il) exit
                     k = k + 1_ik
                  end do
                  over_lai  = cum
                  layer_lai = 0.0_wp
                  do j = i, k
                     cohort%overtopping_lai(j) = over_lai
                     layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
                  end do
                  cum = cum + layer_lai                     ! whole layer shades those below
                  i   = k + 1_ik
               end do
            end block
         end do
      end associate
   end subroutine update_overtopping_lai

   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_depth -- each patch's canopy-air depth, the tallest cohort plus the    !
   ! freeboard, floored. The only writer of can_depth, which sets the canopy air's         !
   ! capacity and the height the forcing is moved to (canopy_air_space_biophysics.md;      !
   ! forcing.md Step 2).                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_depth(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_ledger_t), intent(inout) :: ledger   !< an inactive one books nothing
      real(wp)    :: h_top, depth_new, de, dw, dc
      integer(ik) :: ip, i0, i1, i
      logical     :: book
      if (.not. allocated(poly%patch%cas)) return
      book = ledger%active   ! without an active ledger the resize is a plain geometry update
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         h_top = 0.0_wp
         do i = i0, i1
            h_top = max(h_top, poly%cohort%height(i))
         end do
         depth_new = max(cfg%aero%min_canopy_depth, h_top + cfg%aero%canopy_freeboard)
         if (book) then
            !----- Air entrained or detrained is declared (conservation.md §3). -----------!
            call cas_set_depth(poly%patch%cas(ip), depth_new, rho_air=rho_air_ref,                 &
                               de_open=de, dw_open=dw, dc_open=dc)
            call slow_ledger_declare(ledger,                                                       &
                     energy_in=patch_polygon_area(poly, ip) * de, water_in=patch_polygon_area(poly, ip) * dw, &
                     carbon_in=patch_polygon_area(poly, ip) * dc * umol_2_kgc)
         else
            call cas_set_depth(poly%patch%cas(ip), depth_new)
         end if
      end do
   end subroutine refresh_canopy_depth

   !=======================================================================================!
   ! 2. Recruitment: found a cohort from each full recruit pool                            !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Recruitment: a recruit pool the founding law releases founds one cohort at the        !
   ! smallest tracked height and is reset; a smaller pool carries over. The driver credits !
   ! the pool daily (accumulate_recruit_pool) and this consumes it monthly (demography.md  !
   ! §2-§3). It changes the cohort count, so it lives with the cohort operators.           !
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(poly, cfg, found, carbon_drawn, heat_drawn, water_drawn)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      logical,             intent(in)    :: found(:,:)   !< (pft, patch) the pools that found a cohort
      !----- What the recruits drew from outside the model, for the caller to declare to the   !
      !      slow ledger: their carbon beyond the carbon_min per plant the pool paid for, and  !
      !      their tissues' heat and water. Recruitment stands in for the sizes below          !
      !      min_cohort_height, which the model does not represent, so the draw comes from the !
      !      free atmosphere, not the canopy air (demography.md §3; conservation.md §3). ------!
      real(wp),            intent(out)   :: carbon_drawn  !< [kgC/m2 site]
      real(wp),            intent(out)   :: heat_drawn    !< [J/m2 site]
      real(wp),            intent(out)   :: water_drawn   !< [kg/m2 site]
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh, endow, t_birth, psi_birth, c_tot, e_tot, w_tot

      c_tot = 0.0_wp
      e_tot = 0.0_wp
      w_tot = 0.0_wp
      carbon_drawn = 0.0_wp
      heat_drawn   = 0.0_wp
      water_drawn  = 0.0_wp
      np = poly%patch%n
      if (np < 1_ik) return

      !----- Count the pools that found a cohort. -------------------------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, poly%n_pft
            if (found(pf, ip)) nspawn = nspawn + 1_ik
         end do
      end do
      !----- No pool reached it: nothing is born, so nothing is drawn (the seed rain is declared !
      !      daily, as it arrives, by accumulate_recruit_pool). ---------------------------------!
      if (nspawn == 0_ik) return

      call expand_cohort_block(poly%cohort, poly%cohort%n + nspawn)
      n_before = poly%cohort%n
      m = poly%cohort%n
      associate (cohort => poly%cohort, patch => poly%patch, pft => poly%cohort%pft_table)
         do ip = 1_ik, np
            do pf = 1_ik, poly%n_pft
               if (.not. found(pf, ip)) cycle
               m = m + 1_ik
               !----- Every PFT recruits at the smallest tracked height, on its own curve. -------!
               recruit_dbh = height_to_dbh(cohort%pft_table%allom(pf), pft%min_cohort_height)
               !----- Born in its patch's state: at the canopy-air temperature, holding the water     !
               !      its tissues hold at the top soil layer's potential (patch_physical_state). -------!
               call patch_physical_state(patch%cas(ip), patch%soil_water(ip),                          &
                                         poly%site%surface(patch%owner_site(ip)), t_birth, psi_birth)
               call init_cohort(cohort, m, pft, pf, ip, patch%recruit_pool(pf, ip), recruit_dbh,   &
                                t_birth, psi_birth)
               !----- The draw: everything the new cohort holds, less the reproduction carbon the !
               !      pool paid for at carbon_min per plant. -------------------------------------!
               endow = cohort%leaf_carbon(m) + cohort%fineroot_carbon(m)                           &
                     + cohort%wood_carbon(m) + cohort%nonstructural_carbon(m)
               c_tot = c_tot + patch_polygon_area(poly, ip) * patch%recruit_pool(pf, ip)            &
                             * (endow - pft%recruit_carbon(pf))
               e_tot = e_tot + patch_polygon_area(poly, ip) * cohort%nplant(m)                       &
                             * (cohort%leaf_energy(m) + cohort%wood_energy(m))
               w_tot = w_tot + patch_polygon_area(poly, ip) * cohort_tissue_water(cohort, m)
               patch%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         cohort%n = m
      end associate

      carbon_drawn = c_tot
      heat_drawn   = e_tot
      water_drawn  = w_tot

      !----- Each new cohort takes a new id. ----------------------------------------------!
      do k = n_before + 1_ik, poly%cohort%n
         call assign_cohort_id(poly, k)
      end do

      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine apply_recruitment

   !=======================================================================================!
   ! 3. Disturbance: each site's treefall gap                                              !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Treefall gaps (demography.md §4): the share `frac` of every patch's area the law      !
   ! decided becomes one age-0 gap per site. The cohorts it says survive move into it, the !
   ! rest die in it, and each donor keeps (1-frac) of its area. Site                       !
   ! area is conserved; plant number is not, since the killed canopy is the disturbance.   !
   !---------------------------------------------------------------------------------------!
   subroutine apply_patch_disturbance(poly, cfg, frac, survives, water_shed, kills)
      type(polygon_t),       intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      real(wp),            intent(in)    :: frac          !< [-] the disturbed share (decide_disturbance)
      logical,             intent(in)    :: survives(:)   !< (cohort) survives into the gap
      !----- What the killed canopy handed on, for the caller to declare to the slow ledger. ---!
      real(wp),            intent(out)   :: water_shed  !< [kg/m2 polygon] tissue + film water -> the gaps' ground
      type(kill_report_t), intent(out)   :: kills       !< the killed canopy's carbon (kill_litter)
      integer(ik) :: np0, newp, d, i, i0, i1, m, m0, nsurv, is, p0, p1
      real(wp)    :: new_area, atot
      real(wp), allocatable :: donor_area(:)
      real(wp)    :: lost_density
      real(wp)    :: kill_frac, w_kill, w_tot

      w_tot = 0.0_wp
      water_shed = 0.0_wp
      call open_kill_report(kills, max(poly%cohort%n, 1_ik))

      np0 = poly%patch%n
      if (np0 < 1_ik .or. frac <= tiny_num) return

      !----- Disturbed area fraction this step, and the aggregate new-patch area. ----------!
      new_area = frac * sum(poly%patch%area(1:np0))
      if (new_area <= tiny_num) return

      !----- The disturbed area, written on the donors before the gaps are added; each gap inherits !
      !      its donors' rows (inherit_patch_diag, below). The rate frac/dt_yr [1/yr] times the     !
      !      weight dt_yr x yr_sec [s] is frac x yr_sec (diagnostics.md §4, "The slow rows"). ------!
      if (poly%patch%diag%active) then
         do i = 1_ik, np0
            poly%patch%diag%v(PA_DISTURB_AREA, i) = poly%patch%diag%v(PA_DISTURB_AREA, i)          &
                                                  + frac * yr_sec
         end do
      end if

      !----- Count the survivors. The other branch writes the disturbance's mortality carbon on the    !
      !      donors, per m2 of their original ground: the shrunk donors and the gap that inherits      !
      !      their rows weigh it by that area, so the site total is exact (dt_yr cancels as above). ---!
      nsurv = 0_ik
      do i = 1_ik, poly%cohort%n
         if (survives(i)) then
            nsurv = nsurv + 1_ik
         else if (poly%patch%diag%active) then
            d = poly%cohort%owner_patch(i)
            poly%patch%diag%v(PA_MORT_C_DISTURB, d) = poly%patch%diag%v(PA_MORT_C_DISTURB, d)      &
               + poly%cohort%nplant(i) * frac * yr_sec                                             &
                 * (poly%cohort%leaf_carbon(i)     + poly%cohort%fineroot_carbon(i)                &
                  + poly%cohort%wood_carbon(i)     + poly%cohort%nonstructural_carbon(i))
         end if
      end do

      call expand_patch_block(poly%patch, np0 + poly%site%n, poly%n_pft)
      call expand_cohort_block(poly%cohort, poly%cohort%n + nsurv)
      m0 = poly%cohort%n                         ! cohorts before the survivors are moved in
      allocate(donor_area(np0))

      !----- One gap per site, cut from that site's patches: their areas are fractions of it. The  !
      !      gaps go at the end of the block, and sort_patches below puts each in its site. -------!
      newp = np0
      do is = 1_ik, poly%site%n
         p0 = poly%site%patch_offset(is)
         p1 = p0 + poly%site%patch_count(is) - 1_ik
         if (p1 < p0) cycle
         new_area = frac * sum(poly%patch%area(p0:p1))
         if (new_area <= tiny_num) cycle
         newp = newp + 1_ik
         call assign_patch_id(poly, newp)   ! the gap's id, which the kill report names it by

         associate (cohort => poly%cohort, patch => poly%patch)
            !----- A new age-0 treefall gap. ----------------------------------------------!
            patch%area(newp)           = new_area
            patch%age(newp)            = 0.0_wp
            patch%dist_type(newp)      = DIST_TREEFALL
            patch%owner_site(newp)     = is
            !----- A fresh gap has no cohorts of its own yet, so none of its own rates (demography.md §4). ---!
            patch%shed_water_rate(newp) = 0.0_wp
            patch%shed_water_enth(newp) = 0.0_wp
            patch%slow_co2_rate(newp)   = 0.0_wp   ! fresh gap: no cohort of its own respired today
            patch%adapt_dt_last(newp)   = 0.0_wp   ! a fresh gap cold-starts the controller
            patch%fast_cost(newp)       = 0.0_wp   ! and is not yet measured
            patch%n = newp

            !----- The gap inherits the area-weighted mean of its donors' stores, the recruit pool among !
            !      them: a gap does not sterilise the ground it opens (demography.md §4). ---------------!
            atot = sum(patch%area(p0:p1))                       ! donors still hold pre-disturbance area here
            if (atot > tiny_num) then
               call blend_patch_stores(patch, newp, patch%area(p0) / atot, p0, 0.0_wp, p0)
               do d = p0 + 1_ik, p1
                  call blend_patch_stores(patch, newp, 1.0_wp, newp, patch%area(d) / atot, d)
               end do
            end if

            !----- Survivors move in at the density the disturbed share of their patch gives over the gap's !
            !      area; the killed canopy's carbon, at that density, is reported for the gap's litter     !
            !      (soil_carbon.md §5). --------------------------------------------------------------------!
            m = cohort%n
            do d = p0, p1
               i0 = patch%cohort_offset(d)
               i1 = i0 + patch%cohort_count(d) - 1_ik
               do i = i0, i1
                  if (.not. survives(i)) then                                     ! canopy dies in gap
                     !----- The killed canopy's water goes to the gap's ground down the channel turnover !
                     !      shedding uses (conservation.md §3); its tissue heat is a thermal-mass change !
                     !      the caller declares (conservation.md §4). kill_frac is the share of the      !
                     !      cohort's area density the gap takes, the factor the survivor branch uses. ---!
                     kill_frac = frac * patch%area(d)
                     w_kill = kill_frac * (cohort_tissue_water(cohort, i) + cohort_film_water(cohort, i))
                     call add_shed_water(patch, newp, w_kill / new_area,                               &
                                         kill_frac * cohort%nplant(i) * cohort_water_enthalpy(cohort, i)   &
                                         / new_area, cfg%dt_slow)
                     w_tot = w_tot + poly%site%area(is) * w_kill   ! w_kill is per m2 of the site
                     lost_density = cohort%nplant(i) * (frac * patch%area(d) / new_area)
                     call report_kill(kills, patch%global_id(newp), cohort%pft(i),                    &
                              lost_density * cohort%leaf_carbon(i), lost_density * cohort%fineroot_carbon(i), &
                              lost_density * cohort%wood_carbon(i),                                    &
                              lost_density * cohort%nonstructural_carbon(i))
                     cycle
                  end if
                  m = m + 1_ik
                  call copy_cohort(cohort, m, i)
                  cohort%nplant(m)      = cohort%nplant(i) * (frac * patch%area(d) / new_area)
                  cohort%owner_patch(m) = newp
               end do
            end do
            cohort%n = m

            !----- Donors keep the undisturbed remainder (all their cohorts intact). ----------!
            do d = p0, p1
               patch%area(d) = (1.0_wp - frac) * patch%area(d)
            end do
         end associate

         !----- The gap takes its donors' diagnostic rows blended on area, so the stand's         !
         !      Sum(area x value) is unchanged by the carving (diagnostics.md §4). Every donor     !
         !      lost the same share, so their current areas blend as their original ones. --------!
         donor_area = 0.0_wp
         donor_area(p0:p1) = poly%patch%area(p0:p1)
         call inherit_patch_diag(poly%patch%diag, newp, donor_area)
      end do

      water_shed = w_tot

      !----- The survivors copied into the gaps take new ids; the cohorts they came from keep theirs. ---!
      do i = m0 + 1_ik, poly%cohort%n
         call assign_cohort_id(poly, i)
      end do

      call sort_patches(poly)                    ! each gap into its site's slice; both maps rebuilt
      call sort_cohorts(poly)
   end subroutine apply_patch_disturbance

   !=======================================================================================!
   ! 4. Resolution: sort, fuse, split and cull cohorts; sort, fuse and terminate patches   !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Sort each patch's cohorts tallest first: height descending, diameter descending on    !
   ! ties. The overtopping LAI and the fusion scans rely on the order (demography.md §3).  !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik), allocatable :: perm(:)
      integer(ik) :: ip, i0, i1, i, j, key
      integer(ik) :: n

      n = poly%cohort%n
      if (n < 1_ik) return
      allocate(perm(n))
      do i = 1_ik, n
         perm(i) = i
      end do
      associate (cohort => poly%cohort, patch => poly%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            !----- Insertion sort of the index window perm(i0:i1), descending. ------------!
            do i = i0 + 1_ik, i1
               key = perm(i)
               j   = i - 1_ik
               do while (j >= i0)
                  if (higher(cohort%height(perm(j)), cohort%dbh(perm(j)), cohort%height(key), cohort%dbh(key))) exit
                  perm(j + 1_ik) = perm(j)
                  j = j - 1_ik
               end do
               perm(j + 1_ik) = key
            end do
         end do
      end associate
      call reorder_cohort_block(poly%cohort, perm, n)
   end subroutine sort_cohorts

   !----- True if (h1,d1) ranks at least as high as (h2,d2) under the sort order. ----------!
   pure logical function higher(h1, d1, h2, d2)
      real(wp), intent(in) :: h1, d1, h2, d2
      if (h1 > h2) then
         higher = .true.
      else if (h1 < h2) then
         higher = .false.
      else
         higher = (d1 >= d2)
      end if
   end function higher

   integer(ik) function max_cohort_count(poly)
      type(polygon_t), intent(in) :: poly
      if (poly%patch%n < 1_ik) then
         max_cohort_count = 0_ik
      else
         max_cohort_count = maxval(poly%patch%cohort_count(1:poly%patch%n))
      end if
   end function max_cohort_count

   !---------------------------------------------------------------------------------------!
   ! Cohort fusion in passes, the tolerance growing geometrically from                     !
   ! cohort_size_tol_min, until no patch holds more than max_cohort cohorts              !
   ! (demography.md §3).                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_cohorts(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: ifus

      tol = cfg%demography%cohort_size_tol_min
      do ifus = 1_ik, cfg%demography%n_cohort_fusion_iter
         call fuse_pass(poly, cfg, tol)
         if (max_cohort_count(poly) <= cfg%demography%max_cohort) exit
         tol = tol * cfg%demography%cohort_size_tol_mult
      end do
   end subroutine fuse_cohorts

   !----- One fusion sweep over all patches at a fixed tolerance. --------------------------!
   subroutine fuse_pass(poly, cfg, tol)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical, allocatable :: alive(:)
      integer(ik) :: ip, i0, i1, recc, donc
      real(wp)    :: diff_height, lai_comb, h_top
      logical     :: did_fuse

      did_fuse = .false.
      allocate(alive(poly%cohort%n))
      alive = .true.
      associate (cohort => poly%cohort, patch => poly%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            if (i1 < i0) cycle
            !----- The tolerance is a share of the patch's top height (its first cohort), so it !
            !      scales with the stand: a cohort that pulls ahead in a young stand keeps its  !
            !      place in the light, and a dense stand can form the size hierarchy it thins   !
            !      by (demography.md §3). ------------------------------------------------------!
            h_top = cohort%height(i0)
            do recc = i0, i1 - 1_ik
               if (.not. alive(recc)) cycle
               do donc = recc + 1_ik, i1
                  if (.not. alive(donc)) cycle
                  if (cohort%pft(donc) /= cohort%pft(recc)) cycle
                  lai_comb = cohort%nplant(recc) * cohort%leaf_area(recc) + cohort%nplant(donc) * cohort%leaf_area(donc)
                  if (lai_comb >= cfg%demography%cohort_lai_cap) cycle
                  diff_height = abs(cohort%height(donc) - cohort%height(recc))
                  if (diff_height < h_top * tol) then
                     call fuse_2_cohorts(poly, recc, donc, cfg%demography%conservation_tol)
                     alive(donc) = .false.
                     did_fuse    = .true.
                  end if
               end do
            end do
         end do
      end associate
      if (did_fuse) then
         call compact_cohort_block(poly%cohort, alive)
         call rebuild_cohort_map(poly)
         call sort_cohorts(poly)
      end if
   end subroutine fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge cohort donc into recc. Plant number and the four carbon pools are conserved     !
   ! (each pool a density-weighted mean per plant), and with them biomass and leaf area;   !
   ! the geometry is re-derived from the merged carbon, never averaged (demography.md §3). !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(poly, recc, donc, conservation_tol)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: conservation_tol
      real(wp) :: nr, nd, ntot, agb_tot, agb_new, wr, wd, wtot
      real(wp) :: lc_tot, fc_tot, wc_tot, nc_tot                !< pre-merge pool totals [kgC/m2]
      real(wp) :: cr, cd                                        !< the two cohorts' leaf carbon [kgC/m2]
      real(wp) :: hlr, hld, hwr, hwd                            !< their leaf and wood heat capacities [J/m2/K]

      associate (cohort => poly%cohort)
         nr      = cohort%nplant(recc)
         nd      = cohort%nplant(donc)
         ntot    = nr + nd
         agb_tot = nr * cohort%agb(recc) + nd * cohort%agb(donc)     ! [kgC/m2] conserved
         !----- The fast-loop state fuses by the policy fuse_cohort_fast_state declares per field   !
         !      (each per-plant amount, the tissue energies among them, density-weighted). The leaf !
         !      traits are set before the geometry is re-derived, since sla maps leaf carbon to     !
         !      leaf area; each blends on the quantity it is per (plant_traits.md §5). -------------!
         wr = nr * cohort%leaf_area(recc)
         wd = nd * cohort%leaf_area(donc)
         wtot = wr + wd
         cr = nr * cohort%leaf_carbon(recc)
         cd = nd * cohort%leaf_carbon(donc)
         !----- The tissue heat capacities before fuse_cohort_fast_state blends the water into   !
         !      them: the diagnostic temperatures fuse on the weights the prognostic ones do. ---!
         call cohort_tissue_heat_capacity(cohort, recc, hlr, hwr)
         call cohort_tissue_heat_capacity(cohort, donc, hld, hwd)
         call fuse_cohort_fast_state(cohort, recc, donc, nr, nd)
         cohort%sla(recc)     = pair_mean(cohort%sla(recc),     cr, cohort%sla(donc),     cd)
         cohort%vcmax25(recc) = pair_mean(cohort%vcmax25(recc), wr, cohort%vcmax25(donc), wd)
         cohort%rd25(recc)    = pair_mean(cohort%rd25(recc),    wr, cohort%rd25(donc),    wd)
         cohort%leaf_turnover_rate(recc) = pair_mean(cohort%leaf_turnover_rate(recc), cr,                &
                                                     cohort%leaf_turnover_rate(donc), cd)
         !----- The diagnostic rows fuse on leaf area, density or heat capacity, as ed_diag_vars   !
         !      declares per field (CO_FUSION_KIND). ----------------------------------------------!
         call fuse_cohort_diag(cohort%diag,  recc, donc, wr, wd, nr, nd, hlr, hld, hwr, hwd)
         call fuse_cohort_diag(cohort%sdiag, recc, donc, wr, wd, nr, nd, hlr, hld, hwr, hwd)
         !----- The survivor keeps its own growth history and predawn water status: the two !
         !      cohorts already match in height and LAI, so these differ little. -----------!
         cohort%nplant(recc) = ntot
         !----- The four carbon pools, density-weighted, then the geometry from them;       !
         !      set_cohort_size_from_carbon reads the pools and leaves them as they are. ---!
         lc_tot = nr * cohort%leaf_carbon(recc)          + nd * cohort%leaf_carbon(donc)
         fc_tot = nr * cohort%fineroot_carbon(recc)      + nd * cohort%fineroot_carbon(donc)
         wc_tot = nr * cohort%wood_carbon(recc)          + nd * cohort%wood_carbon(donc)
         nc_tot = nr * cohort%nonstructural_carbon(recc) + nd * cohort%nonstructural_carbon(donc)
         cohort%leaf_carbon(recc)          = lc_tot / ntot
         cohort%fineroot_carbon(recc)      = fc_tot / ntot
         cohort%wood_carbon(recc)          = wc_tot / ntot
         cohort%nonstructural_carbon(recc) = nc_tot / ntot
         call set_cohort_size_from_carbon(poly%cohort, recc)   ! dbh=wood_to_dbh(wood_carbon); pools kept
         !----- The merged tissues' temperatures, from the energies fused above and the merged masses. !
         call set_cohort_tissue_temps(poly%cohort, recc)
         !----- Each pool's conservation, asserted; the pools are non-negative, so each total is !
         !      also its magnitude. -------------------------------------------------------------!
         if (.not. (conserved(lc_tot, ntot*cohort%leaf_carbon(recc),          lc_tot, conservation_tol) .and. &
                    conserved(fc_tot, ntot*cohort%fineroot_carbon(recc),      fc_tot, conservation_tol) .and. &
                    conserved(wc_tot, ntot*cohort%wood_carbon(recc),          wc_tot, conservation_tol) .and. &
                    conserved(nc_tot, ntot*cohort%nonstructural_carbon(recc), nc_tot, conservation_tol)))    &
            error stop 'fuse_2_cohorts: carbon-pool conservation violated'
         !----- Biomass (a share of the wood carbon) and leaf area (sla x leaf carbon) follow from !
         !      the pools and the leaf-carbon-weighted sla, conserved above; asserted too. --------!
         agb_new = cohort%nplant(recc) * cohort%agb(recc)
         if (.not. conserved(agb_tot, agb_new, agb_tot, conservation_tol)) then
            write(*,'(a,2es24.16)') ' fuse_2_cohorts: AGB before, after ', agb_tot, agb_new
            error stop 'fuse_2_cohorts: AGB conservation violated'
         end if
         if (.not. conserved(wtot, ntot*cohort%leaf_area(recc), wtot, conservation_tol)) then
            write(*,'(a,2es24.16)') ' fuse_2_cohorts: leaf area before, after ', wtot,          &
                                    ntot*cohort%leaf_area(recc)
            error stop 'fuse_2_cohorts: leaf-area conservation violated'
         end if
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose LAI exceeds cohort_lai_cap into two of half the density,     !
   ! repeating until none does (demography.md §3).                                         !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: eps, agb_before, agb_after, wc0

      eps = cfg%demography%split_eps

      do iter = 1_ik, cfg%demography%n_cohort_fusion_iter
         n0 = poly%cohort%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (poly%cohort%nplant(i) * poly%cohort%leaf_area(i) > cfg%demography%cohort_lai_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         agb_before = sum(poly%cohort%nplant(1:n0) * poly%cohort%agb(1:n0))
         call expand_cohort_block(poly%cohort, n0 + nsplit)
         m = n0
         associate (cohort => poly%cohort)
            do i = 1_ik, n0
               if (cohort%nplant(i) * cohort%leaf_area(i) <= cfg%demography%cohort_lai_cap) cycle
               cohort%nplant(i) = 0.5_wp * cohort%nplant(i)
               m = m + 1_ik
               !----- The wood carbon, not the diameter, moves by +/-eps and the other pools stay !
               !      equal, so every pool and the biomass are conserved exactly, even where the  !
               !      height is capped (demography.md §3). ---------------------------------------!
               wc0 = cohort%wood_carbon(i)
               call copy_cohort(cohort, m, i)           ! copies halved nplant + params + pools
               cohort%wood_carbon(i) = wc0 * (1.0_wp + eps)
               cohort%wood_carbon(m) = wc0 * (1.0_wp - eps)
               call set_cohort_size_from_carbon(cohort, i)
               call set_cohort_size_from_carbon(cohort, m)
               !----- Both daughters keep their parent's temperatures: the wood carbon each gained  !
               !      or lost is held at them, and the two together hold the parent's energy. ----!
               call hold_cohort_tissue_temps(cohort, i)
               call hold_cohort_tissue_temps(cohort, m)
            end do
            cohort%n = m
            agb_after  = sum(cohort%nplant(1:m) * cohort%agb(1:m))
         end associate
         !----- The '-eps' daughters (slots n0+1..m) are new cohorts and take new ids; the !
         !      '+eps' half keeps slot i and its parent's id. -----------------------------!
         do i = n0 + 1_ik, m
            call assign_cohort_id(poly, i)
         end do
         if (.not. conserved(agb_before, agb_after, agb_before, cfg%demography%conservation_tol))             &
            error stop 'split_cohorts: AGB conservation violated'
         call rebuild_cohort_map(poly)
         call sort_cohorts(poly)
      end do
   end subroutine split_cohorts

   !---------------------------------------------------------------------------------------!
   ! The cull: remove the cohorts whose n x AGB is below min_cohort_agb, and those at or   !
   ! below the density floor; mortality stops a cohort exactly on the floor, so the test   !
   ! includes equality. A culled cohort's water goes to the ground; its carbon is reported !
   ! in `kills`, for its patch's litter (demography.md §3; soil_carbon.md §5).             !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(poly, cfg, water_shed, kills)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t), intent(in)    :: cfg
      !----- What the cull handed on, for the caller to declare to the slow ledger            !
      !      (conservation.md §3). ------------------------------------------------------------!
      real(wp),           intent(out) :: water_shed  !< [kg/m2 site] tissue + film water -> the ground
      type(kill_report_t), intent(out) :: kills      !< the culled cohorts' carbon (kill_litter)
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n, ip
      real(wp)             :: w_tot, w_cohort

      w_tot = 0.0_wp
      water_shed = 0.0_wp
      n = poly%cohort%n
      call open_kill_report(kills, max(n, 1_ik))
      if (n < 1_ik) return
      allocate(keep(n))
      associate (cohort => poly%cohort, pft => poly%cohort%pft_table, patch => poly%patch)
         do i = 1_ik, n
            keep(i) = (cohort%nplant(i) * cohort%agb(i) >= cfg%demography%min_cohort_agb) .and.               &
                      (cohort%nplant(i) > cfg%demography%negligible_nplant)
            if (keep(i)) cycle
            ip = cohort%owner_patch(i)
            !----- The culled cohort's water goes to the ground down the channel turnover        !
            !      shedding uses (conservation.md §3). The heat of its tissues is a thermal-mass !
            !      change the caller declares with the phase's (conservation.md §4). ------------!
            w_cohort = cohort_tissue_water(cohort, i) + cohort_film_water(cohort, i)
            call add_shed_water(patch, ip, w_cohort, cohort%nplant(i) * cohort_water_enthalpy(cohort, i), &
                                cfg%dt_slow)
            w_tot = w_tot + patch_polygon_area(poly, ip) * w_cohort
            !----- The cull's mortality carbon, a rate [kgC/m2/yr] times the step's seconds, on    !
            !      the cohort's patch, which the cohort restructuring leaves in place              !
            !      (diagnostics.md §4, "The slow rows"). ------------------------------------------!
            if (patch%diag%active)                                                                 &
               patch%diag%v(PA_MORT_C_CULL, ip) = patch%diag%v(PA_MORT_C_CULL, ip)                 &
                  + (cohort%nplant(i) * (cohort%leaf_carbon(i) + cohort%fineroot_carbon(i)         &
                                       + cohort%wood_carbon(i) + cohort%nonstructural_carbon(i))   &
                     / cfg%dt_years) * cfg%dt_slow
            call report_kill(kills, patch%global_id(ip), cohort%pft(i),                             &
                     cohort%nplant(i) * cohort%leaf_carbon(i), cohort%nplant(i) * cohort%fineroot_carbon(i), &
                     cohort%nplant(i) * cohort%wood_carbon(i),                                     &
                     cohort%nplant(i) * cohort%nonstructural_carbon(i))
         end do
      end associate
      water_shed = w_tot
      if (all(keep)) return
      call compact_cohort_block(poly%cohort, keep)
      call rebuild_cohort_map(poly)
   end subroutine terminate_cohorts

   !---------------------------------------------------------------------------------------!
   ! Sort by site, then oldest first; remap owner_patch and rebuild both maps.             !
   !---------------------------------------------------------------------------------------!
   subroutine sort_patches(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik), allocatable :: pperm(:), inv(:)
      integer(ik) :: np, i, j, key, k
      np = poly%patch%n
      if (np < 2_ik) return
      allocate(pperm(np), inv(np))
      do i = 1_ik, np
         pperm(i) = i
      end do
      associate (patch => poly%patch, cohort => poly%cohort)
         !----- Insertion sort patch index by site, then by age descending. ----------------!
         do i = 2_ik, np
            key = pperm(i)
            j   = i - 1_ik
            do while (j >= 1_ik)
               if (patch%owner_site(pperm(j)) < patch%owner_site(key)) exit
               if (patch%owner_site(pperm(j)) == patch%owner_site(key) .and.                        &
                   patch%age(pperm(j)) >= patch%age(key)) exit
               pperm(j + 1_ik) = pperm(j)
               j = j - 1_ik
            end do
            pperm(j + 1_ik) = key
         end do
         !----- Apply the permutation to every patch field, the diagnostics included. ---------!
         call reorder_patch_block(patch, pperm, np)
         !----- Remap owner_patch: old index -> new position. -----------------------------!
         do k = 1_ik, np
            inv(pperm(k)) = k
         end do
         do i = 1_ik, cohort%n
            cohort%owner_patch(i) = inv(cohort%owner_patch(i))
         end do
      end associate
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine sort_patches

   !---------------------------------------------------------------------------------------!
   ! A patch's light profile, eq. (10) of demography.md §4: the fraction of full sunlight  !
   ! left at each height layer, exp(-k x the LAI at and above it), with k =                !
   ! light_ext_coeff_fusefiss. Layer 1 is the ground, layer n_height_layers the top.       !
   !---------------------------------------------------------------------------------------!
   subroutine patch_light_profile(poly, cfg, ip, light)
      type(polygon_t),       intent(in)  :: poly
      type(ed_config_t),   intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip
      real(wp),            intent(out) :: light(:)         ! size cfg%demography%n_height_layers
      call light_profile_of(poly, cfg, ip, poly%patch%cohort_offset(ip),                                   &
                            poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik, light)
   end subroutine patch_light_profile

   !----- The light profile of patch ip's cohorts among positions i0..i1 of the cohort block: its !
   !      slice, or in a fusion pass the slice widened to the patches already fused into it. -----!
   subroutine light_profile_of(poly, cfg, ip, i0, i1, light)
      type(polygon_t),     intent(in)  :: poly
      type(ed_config_t),   intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip, i0, i1
      real(wp),            intent(out) :: light(:)         ! size cfg%demography%n_height_layers
      real(wp)    :: cum_lai(cfg%demography%n_height_layers)
      integer(ik) :: i, ihgt, nl
      nl = cfg%demography%n_height_layers
      cum_lai = 0.0_wp
      do i = i0, i1
         if (poly%cohort%owner_patch(i) /= ip) cycle
         ihgt = min(nl, max(1_ik, count(cfg%demography%height_edges < poly%cohort%height(i)) + 1_ik))
         cum_lai(ihgt) = cum_lai(ihgt) + poly%cohort%nplant(i) * poly%cohort%leaf_area(i)
      end do
      !----- Integrate top-down: cum_lai(ihgt) = LAI in this layer and every layer above. --!
      do ihgt = nl - 1_ik, 1_ik, -1_ik
         cum_lai(ihgt) = cum_lai(ihgt) + cum_lai(ihgt + 1_ik)
      end do
      do ihgt = 1_ik, nl
         light(ihgt) = exp(-cfg%demography%light_ext_coeff_fusefiss * cum_lai(ihgt))
      end do
   end subroutine light_profile_of

   !---------------------------------------------------------------------------------------!
   ! Patch fusion (demography.md §4) in passes, the tolerance growing from                 !
   ! patch_light_tol to patch_light_tol_max, until no site holds more than max_patch       !
   ! patches. Patches less alike than the last tolerance stay apart even then. Only        !
   ! patches of one site fuse.                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_patches(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: it

      tol = cfg%demography%patch_light_tol
      do it = 1_ik, cfg%demography%n_patch_fusion_iter
         call patch_fuse_pass(poly, cfg, tol)
         if (maxval(poly%site%patch_count(1:poly%site%n)) <= cfg%demography%max_patch) exit
         tol = tol * cfg%demography%patch_light_tol_mult
      end do
   end subroutine fuse_patches

   !----- One fusion pass at a fixed tolerance: a receptor only absorbs patches after it, and an !
   !      absorbed patch absorbs nothing, so the cohort block can stay as it is until the pass   !
   !      ends. patch_compact then regroups it once. --------------------------------------------!
   subroutine patch_fuse_pass(poly, cfg, tol)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,     allocatable :: alive(:)
      real(wp),    allocatable :: prof(:,:), lr(:), ld(:)
      integer(ik), allocatable :: last(:), ncoh(:)   !< each patch's last cohort position and its count
      real(wp)    :: davg, dmax, dif
      integer(ik) :: np, recp, donp, ihgt, npop, nl
      logical     :: similar, both_empty

      np = poly%patch%n
      if (np < 2_ik) return
      nl = cfg%demography%n_height_layers
      allocate(alive(np), prof(nl, np), lr(nl), ld(nl), last(np), ncoh(np))
      alive = .true.
      do recp = 1_ik, np
         ncoh(recp) = poly%patch%cohort_count(recp)
         last(recp) = poly%patch%cohort_offset(recp) + ncoh(recp) - 1_ik
         call light_profile_of(poly, cfg, recp, poly%patch%cohort_offset(recp), last(recp), prof(:, recp))
      end do

      do recp = 1_ik, np - 1_ik
         if (.not. alive(recp)) cycle
         do donp = recp + 1_ik, np
            if (.not. alive(donp)) cycle
            if (poly%patch%owner_site(donp) /= poly%patch%owner_site(recp)) cycle
            if (poly%patch%dist_type(donp) /= poly%patch%dist_type(recp)) cycle
            both_empty = (ncoh(recp) == 0_ik .and. ncoh(donp) == 0_ik)
            if (both_empty) then
               similar = .true.
            else
               lr = prof(:, recp)
               ld = prof(:, donp)
               davg = 0.0_wp
               dmax = 0.0_wp
               npop = 0_ik
               do ihgt = 1_ik, nl
                  !----- Skip layers above both canopies (both at full light). -----------!
                  if (lr(ihgt) >= almost_one .and. ld(ihgt) >= almost_one) cycle
                  dif  = abs(lr(ihgt) - ld(ihgt))
                  davg = davg + dif
                  dmax = max(dmax, dif)
                  npop = npop + 1_ik
               end do
               if (npop > 0_ik) davg = davg / real(npop, wp)
               similar = (davg <= tol) .and. (dmax <= tol * cfg%demography%patch_light_maxdev_factor)
            end if
            if (similar) then
               call fuse_patch_into(poly, recp, donp, poly%patch%cohort_offset(recp), last(recp))
               last(recp) = max(last(recp), last(donp))
               ncoh(recp) = ncoh(recp) + ncoh(donp)
               call light_profile_of(poly, cfg, recp, poly%patch%cohort_offset(recp), last(recp), prof(:, recp))
               alive(donp) = .false.
            end if
         end do
      end do

      call patch_compact(poly, alive)
      call sort_cohorts(poly)
   end subroutine patch_fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge patch donp into recp: the stores blend on area and the densities scale by each  !
   ! patch's share of the fused area, which conserves plant number (demography.md §4).     !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(poly, recp, donp)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recp, donp
      call fuse_patch_into(poly, recp, donp, poly%patch%cohort_offset(recp),                               &
                           poly%patch%cohort_offset(recp) + poly%patch%cohort_count(recp) - 1_ik)
   end subroutine fuse_2_patches

   !----- fuse_2_patches, with the receptor's cohorts those it owns among positions r0..r1 (its slice, or a  !
   !      fusion pass's widened one); the donor has absorbed nothing, so its own slice holds its cohorts. ---!
   subroutine fuse_patch_into(poly, recp, donp, r0, r1)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recp, donp, r0, r1
      real(wp)    :: ar, ad, anew, rawgt, dawgt
      integer(ik) :: i, i0, i1

      associate (patch => poly%patch, cohort => poly%cohort)
         ar = patch%area(recp)
         ad = patch%area(donp)
         anew = ar + ad
         if (anew <= tiny_num) then                     ! degenerate: still hand donor cohorts to recp
            i0 = patch%cohort_offset(donp)
            i1 = i0 + patch%cohort_count(donp) - 1_ik
            do i = i0, i1
               cohort%owner_patch(i) = recp
            end do
            return
         end if
         rawgt = ar / anew
         dawgt = ad / anew
         !----- The fields in the patch file's section order; the area, the weight, is set last. ---!
         !----- 1. The age, area-weighted. ---------------------------------------------------------!
         patch%age(recp)            = rawgt * patch%age(recp)            + dawgt * patch%age(donp)
         !----- 3-7, 9. The stores, area-weighted, which conserves each; the temperatures and liquid !
         !      fractions are re-diagnosed in the fast loop. -------------------------------------!
         call blend_patch_stores(patch, recp, rawgt, recp, dawgt, donp)
         !----- 8. The shed water, its enthalpy and the CO2 owed are per-area rates: area-weighted. ---!
         patch%shed_water_rate(recp) = rawgt*patch%shed_water_rate(recp) + dawgt*patch%shed_water_rate(donp)
         patch%shed_water_enth(recp) = rawgt*patch%shed_water_enth(recp) + dawgt*patch%shed_water_enth(donp)
         patch%slow_co2_rate(recp)   = rawgt*patch%slow_co2_rate(recp)   + dawgt*patch%slow_co2_rate(donp)
         !----- 10. adapt_dt_last only seeds the step controller; area-weighted to be deterministic. ---!
         patch%adapt_dt_last(recp)   = rawgt*patch%adapt_dt_last(recp)   + dawgt*patch%adapt_dt_last(donp)
         !----- fast_cost orders the fast loop only: the fused patch carries both patches' work. --------!
         patch%fast_cost(recp)       = patch%fast_cost(recp) + patch%fast_cost(donp)
         !----- 11. The diagnostic rows, each per m2 of its own patch's ground: area-weighted. ---!
         call blend_patch_diag(patch%diag, recp, donp, ar, ad)
         !----- Scale the receptor's densities by its share of the fused area (per-plant amounts stay). ---!
         do i = r0, r1
            if (cohort%owner_patch(i) /= recp) cycle
            cohort%nplant(i) = cohort%nplant(i) * rawgt
         end do
         !----- Scale the donor's densities likewise and hand its cohorts to the receptor. ---!
         i0 = patch%cohort_offset(donp)
         i1 = i0 + patch%cohort_count(donp) - 1_ik
         do i = i0, i1
            cohort%nplant(i)      = cohort%nplant(i) * dawgt
            cohort%owner_patch(i) = recp
         end do
         patch%area(recp) = anew
         patch%area(donp) = 0.0_wp
      end associate
   end subroutine fuse_patch_into

   !---------------------------------------------------------------------------------------!
   ! Termination: merge each patch below min_patch_area into its site's largest surviving  !
   ! patch, then renormalize each site's areas to sum to 1 (demography.md §4).             !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(ed_config_t), intent(in)    :: cfg
      logical, allocatable :: pkeep(:)
      integer(ik)          :: np, ip, imax, is, p0, p1, r0, r1
      real(wp)             :: s

      np = poly%patch%n
      if (np < 1_ik) return
      allocate(pkeep(np))
      pkeep = poly%patch%area(1:np) >= cfg%demography%min_patch_area
      !----- Every site keeps at least its largest patch. -------------------------------------!
      do is = 1_ik, poly%site%n
         p0 = poly%site%patch_offset(is)
         p1 = p0 + poly%site%patch_count(is) - 1_ik
         if (p1 < p0) cycle
         if (.not. any(pkeep(p0:p1))) then
            imax = p0 - 1_ik + maxloc(poly%patch%area(p0:p1), dim=1)
            pkeep(imax) = .true.
         end if
      end do

      if (.not. all(pkeep)) then
         !----- Merged, not deleted: dropping a sliver and renormalizing would redistribute every conserved !
         !      quantity across the site, while merging conserves exactly, so the renormalization below     !
         !      corrects round-off only (demography.md §4). The survivor is in the sliver's own site, and a  !
         !      sliver absorbs nothing, so each merge reads the survivor's cohorts over the site's cohort    !
         !      range and patch_compact regroups the block once, as patch_fuse_pass does. -----------------!
         do is = 1_ik, poly%site%n
            p0 = poly%site%patch_offset(is)
            p1 = p0 + poly%site%patch_count(is) - 1_ik
            if (p1 < p0) cycle
            if (all(pkeep(p0:p1))) cycle
            imax = p0 - 1_ik + maxloc(poly%patch%area(p0:p1), dim=1)
            if (.not. pkeep(imax)) imax = p0 - 1_ik + maxloc(poly%patch%area(p0:p1), dim=1, mask=pkeep(p0:p1))
            r0 = poly%patch%cohort_offset(p0)
            r1 = poly%patch%cohort_offset(p1) + poly%patch%cohort_count(p1) - 1_ik
            do ip = p0, p1
               if (pkeep(ip) .or. ip == imax) cycle
               call fuse_patch_into(poly, imax, ip, r0, r1)
            end do
         end do
         call patch_compact(poly, pkeep)
         call sort_cohorts(poly)   ! tallest first again, which the cohort fusion after it scans by
      end if

      !----- Renormalize each site's areas to sum to 1; after the merges, round-off only. ---!
      do is = 1_ik, poly%site%n
         p0 = poly%site%patch_offset(is)
         p1 = p0 + poly%site%patch_count(is) - 1_ik
         if (p1 < p0) cycle
         s = sum(poly%patch%area(p0:p1))
         if (s > tiny_num) then
            do ip = p0, p1
               poly%patch%area(ip) = poly%patch%area(ip) / s
            end do
         end if
      end do
   end subroutine terminate_patches

   !---------------------------------------------------------------------------------------!
   ! Drop the patches pkeep leaves out, with their cohorts; remap owner_patch to the new   !
   ! patch indices and rebuild both maps.                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine patch_compact(poly, pkeep)
      type(polygon_t), intent(inout) :: poly
      logical,         intent(in)    :: pkeep(:)
      integer(ik), allocatable :: newidx(:), keep_idx(:)
      logical,     allocatable :: ckeep(:)
      integer(ik) :: np, ip, k, i, nc

      np = poly%patch%n
      allocate(newidx(np), keep_idx(max(np, 1_ik)))
      k = 0_ik
      do ip = 1_ik, np
         if (pkeep(ip)) then
            k = k + 1_ik
            newidx(ip) = k
            keep_idx(k) = ip
         else
            newidx(ip) = 0_ik
         end if
      end do

      !----- Cohorts: keep only those whose patch survives; remap owner. ------------------!
      nc = poly%cohort%n
      allocate(ckeep(nc))
      do i = 1_ik, nc
         ckeep(i) = pkeep(poly%cohort%owner_patch(i))
         if (ckeep(i)) poly%cohort%owner_patch(i) = newidx(poly%cohort%owner_patch(i))
      end do
      call compact_cohort_block(poly%cohort, ckeep)

      !----- Patches: keep the survivors, in order; their indices are reorder_patch_block's permutation. ---!
      call reorder_patch_block(poly%patch, keep_idx, k)
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine patch_compact

end module ed_slow_update
