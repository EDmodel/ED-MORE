! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_slow_state -- what the slow processes read from the stand and write back
! (docs/science/order_of_processes.md §4). Each routine loops over the cohorts or patches,
! calls a law of slow_dynamics/processes and writes that process's own fields; the step's
! tendencies and decisions go into slow_step_data_t, which ed_slow_update commits.
!   * leaf: advance_leaf_phenology, advance_plant_traits, reacclimate_plant_traits,
!     acclimate_site_leaf_tables;
!   * carbon: compute_carbon_allocation, compute_vital_rates, accumulate_mortality_litter;
!   * tissue water: shed_turnover_water, fill_new_leaf_water, shed_mortality_water,
!     shed_excess_tissue_water;
!   * soil carbon: advance_soil_carbon.
!---~---
module ed_slow_state
   use ed_kinds,                     only : wp, ik
   use ed_constants,                 only : day_sec, tiny_num, grav_head, kgc_2_umol
   use ed_config,                    only : ed_config_t, growth_window_steps
   use ed_plant_hydr_lib,            only : water_curve_t, tissue_psi_to_water, tissue_water_to_psi,     &
                                            tissue_capped_water
   use ed_therm_lib,                 only : internal_energy_liquid
   use ed_allometry_lib,             only : size_to_leaf_carbon, storage_carbon_target
   use ed_cohort_state_vars,         only : cohort_carbon_supply, cohort_wood_store_carbon,              &
                                            hold_cohort_tissue_temps, cohort_water_enthalpy
   use ed_polygon_state_vars,        only : polygon_t, patch_polygon_area
   use ed_slow_ledger,               only : slow_ledger_t, slow_ledger_declare, slow_ledger_soil_check
   use ed_plant_carbon_allocation,   only : plant_carbon_allocation, cohort_carbon_demand,               &
                                            growth_window_so_far, growth_sink_limitation
   use ed_diag_vars,                 only : CO_NPP_LEAF, CO_NPP_FINEROOT, CO_NPP_WOOD, CO_NPP_STORAGE,   &
                                            CO_NPP_REPRO, CO_GROWTH_RESP, CO_EXUDATE,   &
                                            PA_MORT_C_BACKGROUND
   use ed_leaf_params,               only : leaf_physiology_params_t
   use ed_patch_state_vars,          only : litter_input_t, add_shed_water, soil_carbon_state_t
   use ed_state_defaults,            only : is_unset_positive
   use ed_biogeochem_params,         only : N_SOIL_POOL
   use ed_plant_phenology,           only : pheno_env_t, pheno_state_t, pheno_out_t, phenology_kernel
   use ed_plant_traits,              only : light_plastic_traits, update_plastic_trait,                  &
                                            acclimate_leaf_photo_table
   use ed_soil_carbon,               only : add_necromass, build_litter_input, soil_carbon_step,         &
                                            soil_carbon_out_t
   use ed_demography,                only : founds_cohort, disturbed_fraction, survives_treefall,        &
                                            recruit_efficiency, npp_to_growth, camac_mortality,          &
                                            decline_to_floor_mortality, npp_to_recruitment
   implicit none
   private

   public :: slow_step_data_t
   public :: advance_leaf_phenology, advance_plant_traits, reacclimate_plant_traits, acclimate_site_leaf_tables
   public :: compute_carbon_allocation, compute_vital_rates, accumulate_mortality_litter
   public :: shed_turnover_water, fill_new_leaf_water, shed_mortality_water, shed_excess_tissue_water
   public :: advance_soil_carbon, decide_recruitment, decide_disturbance
   public :: kill_report_t, open_kill_report, report_kill, kill_litter, clear_litter

   !---~---
   !   One slow step's working record: what the processes decided, per cohort unless noted,
   ! for ed_slow_update to commit (docs/science/order_of_processes.md §4).
   !---~---
   type :: slow_step_data_t
      !---~---
      !   The net change of the four carbon pools [kgC/plant]: growth less shed for leaf
      ! and fine root, growth for wood, allocation plus resorption
      ! for storage. Not the reported NPP (the gross CO_NPP_* rows): a canopy that sheds
      ! more than it grows has a negative leaf change, the sign shed_turnover_water reads.
      !---~---
      real(wp), allocatable :: d_leaf(:), d_fineroot(:), d_wood(:), d_nonstructural(:)
      real(wp), allocatable :: npp_repro(:)          !--- [kgC/plant] reproduction growth
      logical,  allocatable :: exhausted(:)          !--- no leaves and no storage (plant_carbon_allocation.md §2)
      real(wp), allocatable :: mortality(:)          !--- [1/yr]
      real(wp), allocatable :: recruitment(:,:)      !--- [plant/m2/yr] (pft, patch)
      !---~---
      !   The density and leaf carbon before the commit, which the mortality litter and
      ! the new leaves' water are measured against.
      !---~---
      real(wp), allocatable :: nplant_before(:)      !--- [plant/m2]
      real(wp), allocatable :: leaf_carbon_before(:) !--- [kgC/plant]
      !---~---
      !   The share of the reproduction carbon that becomes recruits, under the patch's
      ! shade (demography.md eq. (5)); the recruits and the seed lost to litter both read
      ! it, so the two add up to the reproduction carbon.
      !---~---
      real(wp), allocatable :: recruit_efficiency(:) !--- [-]
   end type slow_step_data_t

   !---~---
   !   The plants a stand operator killed (the cull, a treefall gap), one row per cohort:
   ! the patch its litter goes to, by global id, which a re-sort keeps; its PFT; and
   ! its carbon by tissue [kgC/m2 of that patch]. kill_litter makes it litter.
   !---~---
   type :: kill_report_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: patch_id(:), pft(:)
      real(wp),    allocatable :: leaf(:), fineroot(:), wood(:), storage(:)
   end type kill_report_t

contains

   !---~---
   !   Advance every cohort's leaf phenology one slow step (docs/science/plant_phenology.md): build
   ! the daily drivers from the fast loop's accumulators (each site's air temperature, each
   ! cohort's PAR hours) and the cohort's predawn leaf psi, run the kernel, and write back the
   ! tendencies and cue memory. No carbon moves here. A step with no fast sub-steps is skipped,
   ! so the memory never advances on a 0/0 mean.
   !---~---
   subroutine advance_leaf_phenology(poly, cfg, doy)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      integer(ik),         intent(in)    :: doy
      type(pheno_env_t)    :: env
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: i, pf, is
      real(wp)    :: dt_days, temp_day(poly%site%n)
      logical     :: north

      if (all(poly%site%pheno_tair_n(1:poly%site%n) < 1_ik)) return   ! no fast sub-steps -> no drivers
      dt_days  = cfg%dt_slow / day_sec
      do is = 1_ik, poly%site%n
         temp_day(is) = site_air_temp_day(poly, is)
      end do
      north = poly%latitude_deg >= 0.0_wp   ! the season windows' hemisphere

      do i = 1_ik, poly%cohort%n
         pf = poly%cohort%pft(i)

         !---~---
         !   The daily drivers. dmax_psi_leaf is the last completed predawn value; a cohort born
         ! today reads 0, well watered, for one day, as in the leaf kernel.
         !---~---
         env%temp_day         = temp_day(poly%patch%owner_site(poly%cohort%owner_patch(i)))
         !--- Hours the PAR at the cohort's top exceeded par_min over the step, per day.
         env%par_hours        = poly%cohort%light_hours_accum(i) / dt_days
         env%predawn_leaf_psi = min(0.0_wp, poly%cohort%dmax_psi_leaf(i))
         env%doy              = doy
         env%hemis_north      = north

         !--- Pack the cohort's tendencies and cue memory, advance, unpack.
         state%leaf_flush_tendency = poly%cohort%leaf_flush_tendency(i)
         state%leaf_shed_tendency  = poly%cohort%leaf_shed_tendency(i)
         state%growing_degree_days = poly%cohort%growing_degree_days(i)
         state%cold_degree_days    = poly%cohort%cold_degree_days(i)
         state%wet_psi_sum         = poly%cohort%wet_psi_sum(i)
         state%dry_psi_sum         = poly%cohort%dry_psi_sum(i)
         state%light_hours_mean    = poly%cohort%light_hours_mean(i)
         call phenology_kernel(env, poly%cohort%pft_table%pheno(pf), dt_days, state, out)
         poly%cohort%leaf_flush_tendency(i) = state%leaf_flush_tendency
         poly%cohort%leaf_shed_tendency(i)  = state%leaf_shed_tendency
         poly%cohort%growing_degree_days(i) = state%growing_degree_days
         poly%cohort%cold_degree_days(i)    = state%cold_degree_days
         poly%cohort%wet_psi_sum(i)         = state%wet_psi_sum
         poly%cohort%dry_psi_sum(i)         = state%dry_psi_sum
         poly%cohort%light_hours_mean(i)    = state%light_hours_mean
      end do
   end subroutine advance_leaf_phenology

   !---~---
   !   Light trait plasticity over one slow step: each cohort's sla, vcmax25, rd25 and leaf
   ! turnover rate relax toward their targets at the LAI above, as fast as leaves are
   ! replaced. Leaf carbon above the new leaf target is resorbed to storage:
   ! leaf_carbon_new and storage_new are the pools ed_slow_update's commit_trait_carbon
   ! writes. Called only with trait_plasticity_on (docs/science/plant_traits.md eqs.
   ! (1)-(2), §4).
   !---~---
   subroutine advance_plant_traits(poly, dt_yr, leaf_carbon_new, storage_new, instantaneous)
      type(polygon_t),       intent(inout) :: poly
      real(wp),              intent(in)    :: dt_yr
      real(wp), allocatable, intent(out)   :: leaf_carbon_new(:)   !--- [kgC/plant] (cohort)
      real(wp), allocatable, intent(out)   :: storage_new(:)       !--- [kgC/plant] (cohort)
      logical, optional,     intent(in)    :: instantaneous   !--- census restart: jump straight to the target
      integer(ik) :: j, pf
      real(wp)    :: sla_t, vcmax_t, rd_t, turnover_t, turnover_now, leaf_target, excess
      logical     :: instant
      instant = .false.
      if (present(instantaneous)) instant = instantaneous
      allocate(leaf_carbon_new(poly%cohort%n), storage_new(poly%cohort%n))
      associate (cohort => poly%cohort, pft => poly%cohort%pft_table)
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            !--- Light-acclimated targets from this cohort's cumulative LAI above.
            call light_plastic_traits(cohort%overtopping_lai(j),                                 &
                     pft%sla(pf), pft%vcmax25(pf), pft%rd25(pf), pft%leaf_turnover_rate_toc(pf),  &
                     pft%kplastic_sla(pf), pft%kplastic_vcmax(pf), pft%kplastic_rd(pf),             &
                     pft%kplastic_leaf_turnover(pf), sla_t, vcmax_t, rd_t, turnover_t)
            !--- Each trait toward its target, or straight to it on a census start.
            turnover_now = cohort%leaf_turnover_rate(j)
            cohort%sla(j)     = update_plastic_trait(cohort%sla(j),     sla_t,   turnover_now, dt_yr, instant)
            cohort%vcmax25(j) = update_plastic_trait(cohort%vcmax25(j), vcmax_t, turnover_now, dt_yr, instant)
            cohort%rd25(j)    = update_plastic_trait(cohort%rd25(j),    rd_t,    turnover_now, dt_yr, instant)
            cohort%leaf_turnover_rate(j) = update_plastic_trait(cohort%leaf_turnover_rate(j), turnover_t,      &
                                                                turnover_now, dt_yr, instant)
            !---~---
            !   A higher SLA lowers the leaf target: the excess leaf carbon goes to storage, so
            ! leaf area stays at allometry; a lower SLA leaves a deficit growth fills.
            !---~---
            leaf_target = size_to_leaf_carbon(cohort%pft_table%allom(cohort%pft(j)), cohort%dbh(j), cohort%height(j), &
                                           cohort%sla(j))
            leaf_carbon_new(j) = cohort%leaf_carbon(j)
            storage_new(j)     = cohort%nonstructural_carbon(j)
            if (cohort%leaf_carbon(j) > leaf_target) then
               excess = cohort%leaf_carbon(j) - leaf_target
               leaf_carbon_new(j) = leaf_target
               storage_new(j)     = cohort%nonstructural_carbon(j) + excess
            end if
         end do
      end associate
   end subroutine advance_plant_traits

   !---~---
   !   reacclimate_plant_traits -- a restart's plastic leaf traits from this run's PFT file
   ! ([init].reacclimate_traits): the census start's instant acclimation with plasticity on,
   ! the top-of-canopy values without. Leaf area is kept: leaf carbon scales by the SLA's
   ! change, storage taking up the difference, and a lower SLA takes at most the storage there
   ! is, so no carbon is made. leaf_carbon_new and storage_new are the pools
   ! ed_slow_update's commit_trait_carbon writes (docs/configuration.md, "How a run
   ! starts"; plant_traits.md §5).
   !---~---
   subroutine reacclimate_plant_traits(poly, cfg, leaf_carbon_new, storage_new)
      type(polygon_t),       intent(inout) :: poly
      type(ed_config_t),   intent(in)      :: cfg
      real(wp), allocatable, intent(out)   :: leaf_carbon_new(:)   !--- [kgC/plant] (cohort)
      real(wp), allocatable, intent(out)   :: storage_new(:)       !--- [kgC/plant] (cohort)
      integer(ik) :: j, pf
      real(wp)    :: sla_t, vcmax_t, rd_t, turnover_t, sla_old
      allocate(leaf_carbon_new(poly%cohort%n), storage_new(poly%cohort%n))
      associate (cohort => poly%cohort, pft => poly%cohort%pft_table)
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            if (cfg%trait_dynamics%trait_plasticity_on) then
               call light_plastic_traits(cohort%overtopping_lai(j),                              &
                        pft%sla(pf), pft%vcmax25(pf), pft%rd25(pf), pft%leaf_turnover_rate_toc(pf), &
                        pft%kplastic_sla(pf), pft%kplastic_vcmax(pf), pft%kplastic_rd(pf),          &
                        pft%kplastic_leaf_turnover(pf), sla_t, vcmax_t, rd_t, turnover_t)
            else
               sla_t = pft%sla(pf)
               vcmax_t = pft%vcmax25(pf)
               rd_t = pft%rd25(pf)
               turnover_t = pft%leaf_turnover_rate_toc(pf)
            end if
            sla_old           = cohort%sla(j)
            cohort%sla(j)     = sla_t
            cohort%vcmax25(j) = vcmax_t
            cohort%rd25(j)    = rd_t
            cohort%leaf_turnover_rate(j) = turnover_t
            leaf_carbon_new(j) = cohort%leaf_carbon(j)
            storage_new(j)     = cohort%nonstructural_carbon(j)
            if (cohort%sla(j) == sla_old) cycle
            leaf_carbon_new(j) = cohort%leaf_carbon(j) * sla_old / max(cohort%sla(j), tiny_num)
            leaf_carbon_new(j) = min(leaf_carbon_new(j), cohort%leaf_carbon(j) + cohort%nonstructural_carbon(j))
            storage_new(j)     = cohort%nonstructural_carbon(j) + cohort%leaf_carbon(j) - leaf_carbon_new(j)
         end do
      end associate
   end subroutine reacclimate_plant_traits

   !---~---
   !   acclimate_site_leaf_tables -- advance each site's growth-temperature running mean
   ! (weight dt/window, on the last fast loop's daily mean) and acclimate its leaf table
   ! (leaf_gas_exchange.md §2 "Thermal acclimation", eq. (4)). An unset mean takes the
   ! first day's value, so the optimum starts from the site's own climate.
   !---~---
   subroutine acclimate_site_leaf_tables(poly, lp, step_days)
      type(polygon_t),                intent(inout) :: poly
      type(leaf_physiology_params_t), intent(in)    :: lp          !--- [leaf_physiology]
      integer(ik),                    intent(in)    :: step_days
      real(wp)    :: t_day, w
      integer(ik) :: is
      do is = 1_ik, poly%site%n
         if (poly%site%pheno_tair_n(is) < 1_ik) cycle      ! no fast sub-steps ran: nothing to average
         t_day = site_air_temp_day(poly, is)
         associate (t_growth => poly%site%t_growth_avg(is))
            if (is_unset_positive(t_growth)) then
               t_growth = t_day
            else
               w = real(step_days, wp) / max(lp%acclim_window_days, real(step_days, wp))
               t_growth = t_growth + w * (t_day - t_growth)
            end if
            call acclimate_leaf_photo_table(t_growth, lp%acclim_ds_vcmax_a, lp%acclim_ds_vcmax_b,            &
                                            lp%acclim_ds_jmax_a, lp%acclim_ds_jmax_b, lp%acclim_jv_a,        &
                                            lp%acclim_jv_b, poly%site%leaf_photo(is))
         end associate
      end do
   end subroutine acclimate_site_leaf_tables

   !---~---
   !   A site's daily-mean air temperature [K] over the fast loop's sub-steps: its own forcing,
   ! after its terrain. Phenology and thermal acclimation both read it.
   !---~---
   pure real(wp) function site_air_temp_day(poly, is) result(t_day)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      t_day = poly%site%pheno_tair_sum(is) / real(max(poly%site%pheno_tair_n(is), 1_ik), wp)
   end function site_air_temp_day

   !---~---
   !   Each cohort's carbon for the step: gather the day's GPP and maintenance respiration, apply
   ! turnover first (this step's litter) so the growth demands are formed against the post-shed
   ! pools, allocate with plant_carbon_allocation, and form the net per-pool change. Turnover
   ! first lets a canopy that is not senescing replace its loss within the step and hold its
   ! target. The litter and the CO2 owed go to the patch; the pools' changes, the reproduction
   ! carbon and the exhausted flags to the step record
   ! (docs/science/plant_carbon_allocation.md §1-§2, §4-§5).
   !---~---
   subroutine compute_carbon_allocation(poly, cfg, step)
      type(polygon_t),         intent(inout) :: poly   !--- no pool moves here; inout for the litter, rows and patch rates
      type(ed_config_t),     intent(in)      :: cfg
      type(slow_step_data_t),  intent(inout) :: step
      integer(ik) :: n, j, pf, ip
      real(wp)    :: dt_yr
      real(wp)    :: leaf_target, gross_gpp, resp_maint, dt_day, r2l, leaf_turn, root_turn
      real(wp)    :: leaf_demand, fineroot_demand, store_demand, repro_frac
      real(wp)    :: leaf_shed_c, fineroot_shed_c
      real(wp)    :: g_leaf, g_fineroot, g_wood, npp_store, g_repro, growth_resp, deficit
      real(wp)    :: leaf_retained_c, fineroot_retained_c
      real(wp)    :: seed_lost, exudate, window_yr, grown_in_window
      real(wp), allocatable :: co2_owed(:)
      integer(ik) :: n_window, next_pos

      n = poly%cohort%n
      dt_yr  = cfg%dt_years
      dt_day = cfg%dt_slow / day_sec            ! days in the slow step (phenology rates are [1/day])
      allocate(step%d_leaf(n), step%d_fineroot(n), step%d_wood(n), step%d_nonstructural(n),          &
               step%npp_repro(n), step%exhausted(n))
      allocate(co2_owed(poly%patch%n))
      co2_owed = 0.0_wp
      step%recruit_efficiency = cohort_recruit_efficiency(poly)
      !---~---
      !   The growth-history slot this step's sample will fill (commit_vegetation_step advances
      ! growth_hist_pos the same way): the sink limit's window leaves out its old sample.
      !---~---
      n_window = growth_window_steps(cfg)
      next_pos = mod(poly%growth_hist_pos, n_window) + 1_ik
      associate (cohort => poly%cohort, pft => poly%cohort%pft_table, lit => poly%patch%litter_in)
         do j = 1_ik, n
            pf  = cohort%pft(j)
            r2l = pft%root_to_leaf_ratio(pf)
            !--- The leaf target at the cohort's own SLA conserves leaf area (docs/science/plant_traits.md §4).
            leaf_target  = size_to_leaf_carbon(cohort%pft_table%allom(pf), cohort%dbh(j), cohort%height(j), &
                                            cohort%sla(j))
            leaf_turn    = cohort%leaf_turnover_rate(j)   ! baseline leaf turnover [1/yr]
            store_demand = max(0.0_wp, storage_carbon_target(pft%storage_cushion(pf), leaf_target,      &
                                                             r2l * leaf_target) - cohort%nonstructural_carbon(j))
            !--- Reproduction starts at min_repro_height (docs/science/plant_carbon_allocation.md §2).
            repro_frac   = merge(pft%reproduction_investment_fraction(pf), 0.0_wp,                  &
                                 cohort%height(j) >= pft%min_repro_height(pf))
            !---~---
            !   Fine-root turnover, slowed with the leaves' when fineroot_turnover_plastic
            ! (docs/science/plant_traits.md §4); then the turnover and the growth demands.
            !---~---
            root_turn = pft%fineroot_turnover_rate(pf)
            if (cfg%trait_dynamics%fineroot_turnover_plastic)                                                     &
               root_turn = root_turn * cohort%leaf_turnover_rate(j) / pft%leaf_turnover_rate_toc(pf)
            call cohort_carbon_supply(cohort, j, cfg%fast_biophysics_on, cfg%gpp_ref, dt_yr,       &
                                      gross_gpp, resp_maint)
            call cohort_carbon_demand(cohort%leaf_flush_tendency(j), cohort%leaf_shed_tendency(j), &
                     pft%pheno(pf)%flush_rate_max, pft%pheno(pf)%shed_rate_max, leaf_turn,         &
                     root_turn,                                                                    &
                     pft%pheno(pf)%min_leaf_cover, pft%pheno(pf)%bare_leaf_cover,                  &
                     cohort%leaf_carbon(j), cohort%fineroot_carbon(j), leaf_target, dt_day, r2l,   &
                     leaf_shed_c, fineroot_shed_c, leaf_demand,                                    &
                     fineroot_demand, pft%retained_carbon_fraction(pf), leaf_retained_c,           &
                     fineroot_retained_c)
            call growth_window_so_far(cohort%growth_count(j), cohort%growth_accum(j),              &
                                      cohort%growth_hist(next_pos, j), n_window, dt_yr,             &
                                      window_yr, grown_in_window)
            !---~---
            !   Allocate the day's carbon: growth respiration on realized growth, wood capped by
            ! the sink limit (docs/science/plant_carbon_allocation.md §2-§3).
            !---~---
            call plant_carbon_allocation(gpp=gross_gpp, resp_maint=resp_maint,                     &
                     growth_resp_frac=pft%growth_resp_factor(pf), storage=cohort%nonstructural_carbon(j),                    &
                     leaf_demand=leaf_demand, fineroot_demand=fineroot_demand,                     &
                     storage_demand=store_demand, repro_frac=repro_frac,                           &
                     growth_leaf=g_leaf, growth_fineroot=g_fineroot, growth_wood=g_wood,          &
                     npp_store=npp_store, growth_repro=g_repro, growth_resp=growth_resp,          &
                     deficit=deficit,                                                              &
                     wood_growth_max=growth_sink_limitation(                                       &
                        allom=cohort%pft_table%allom(pf),                                            &
                        wood_carbon=cohort%wood_carbon(j), dbh=cohort%dbh(j),                       &
                        wood_density=pft%wood_density(pf),                                          &
                        max_relative_growth_rate=pft%max_relative_growth_rate(pf),                  &
                        max_absolute_growth_rate=pft%max_absolute_growth_rate(pf),                  &
                        max_absolute_growth_exponent=pft%max_absolute_growth_exponent(pf),          &
                        dt_yr=dt_yr, window_yr=window_yr, grown_in_window=grown_in_window),         &
                     exudate=exudate)
            !--- The net per-pool change: growth less the shed for leaf and fine root.
            step%d_leaf(j)          = g_leaf     - leaf_shed_c
            step%d_fineroot(j)      = g_fineroot - fineroot_shed_c
            step%d_wood(j)          = g_wood
            !---~---
            !   Storage: the allocation, plus the resorbed carbon, which d_leaf and d_fineroot
            ! remove in full (docs/science/plant_carbon_allocation.md eq. (2)).
            !---~---
            step%d_nonstructural(j) = npp_store + leaf_retained_c + fineroot_retained_c
            step%npp_repro(j)       = g_repro
            !---~---
            !   No leaves and no storage: the maintenance left unpaid once storage is spent exceeds
            ! all the leaf carbon, so no flush can be paid for
            ! (docs/science/plant_carbon_allocation.md §2, "A cohort with no leaves and no storage").
            !---~---
            step%exhausted(j)       = deficit > cohort%leaf_carbon(j)
            !---~---
            !   Slow cohort rows: the allocation split, growth respiration and exudate, as rates
            ! [kgC/plant/yr], so the carbon budget closes from the output.
            !---~---
            if (cohort%sdiag%active) then
               cohort%sdiag%v(CO_NPP_LEAF,     j) = g_leaf     / dt_yr
               cohort%sdiag%v(CO_NPP_FINEROOT, j) = g_fineroot / dt_yr
               cohort%sdiag%v(CO_NPP_WOOD,     j) = g_wood     / dt_yr
               cohort%sdiag%v(CO_NPP_STORAGE,  j) = npp_store  / dt_yr
               cohort%sdiag%v(CO_NPP_REPRO,    j) = g_repro    / dt_yr
               cohort%sdiag%v(CO_GROWTH_RESP,  j) = growth_resp / dt_yr
               cohort%sdiag%v(CO_EXUDATE,      j) = exudate / dt_yr
               cohort%sdiag%w(j)                  = 1.0_wp
            end if
            !---~---
            !   Growth respiration, less the starvation deficit the fast loop over-reported, is
            ! owed to the canopy air (docs/science/plant_respiration.md §6).
            !---~---
            ip = cohort%owner_patch(j)
            co2_owed(ip) = co2_owed(ip) + cohort%nplant(j) * (growth_resp - deficit)

            !---~---
            !   The seed that does not establish is litter that falls with the canopy
            ! (docs/science/demography.md §1, "Recruitment").
            !---~---
            seed_lost = g_repro * cohort%nplant(j)                                              &
                      * (1.0_wp - min(max(step%recruit_efficiency(j), 0.0_wp), 1.0_wp))
            !---~---
            !   The unresorbed share of the turnover and the seed lost go to the patch's litter
            ! (docs/science/plant_carbon_allocation.md eq. (2)).
            !---~---
            call add_necromass(lit(ip), (leaf_shed_c - leaf_retained_c) * cohort%nplant(j),       &
                     (fineroot_shed_c - fineroot_retained_c) * cohort%nplant(j), 0.0_wp, 0.0_wp,   &
                     seed_lost, pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                       &
                     cohort%pft_table%allom(pf)%aboveground_frac, pft%struct_lignin_frac(pf))
            !--- Root exudate, the carbon the sink limit left, is below-ground labile litter.
            lit(ip)%labile_soil = lit(ip)%labile_soil + exudate * cohort%nplant(j)
            lit(ip)%exudate_c   = lit(ip)%exudate_c   + exudate * cohort%nplant(j)
         end do
      end associate

      !--- The CO2 owed as a molar rate [umol/m2/s] the next fast window emits.
      do ip = 1_ik, poly%patch%n
         poly%patch%slow_co2_rate(ip) = co2_owed(ip) / cfg%dt_slow * kgc_2_umol
      end do
   end subroutine compute_carbon_allocation

   !---~---
   !   Each cohort's recruit efficiency this step, by the law of ed_demography on its
   ! patch's LAI (docs/science/demography.md eq. (5)).
   !---~---
   function cohort_recruit_efficiency(poly) result(efficiency)
      type(polygon_t), intent(in) :: poly
      real(wp), allocatable :: efficiency(:)
      real(wp), allocatable :: patch_lai(:)
      integer(ik) :: j, pf, ip
      allocate(efficiency(poly%cohort%n), patch_lai(max(poly%patch%n, 1_ik)))
      patch_lai = 0.0_wp
      associate (cohort => poly%cohort, pft => poly%cohort%pft_table)
         do j = 1_ik, cohort%n
            ip = cohort%owner_patch(j)
            patch_lai(ip) = patch_lai(ip) + cohort%nplant(j) * cohort%leaf_area(j)
         end do
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            efficiency(j) = recruit_efficiency(pft%recruit_carbon_efficiency(pf), pft%recruit_shade_decay(pf), &
                                               patch_lai(cohort%owner_patch(j)))
         end do
      end associate
   end function cohort_recruit_efficiency

   !---~---
   !   The vital rates from the step's carbon: each cohort's mortality [1/yr] from its wood growth,
   ! and each PFT's recruitment per patch [plant/m2/yr] from the reproduction carbon plus the
   ! seed rain, by the laws of ed_demography (docs/science/demography.md §1), from the wood
   ! and reproduction carbon compute_carbon_allocation put in the step record. The sweep over the
   ! stand stays here, so the laws keep no state.
   !---~---
   subroutine compute_vital_rates(poly, cfg, step)
      type(polygon_t),        intent(in)    :: poly
      type(ed_config_t),      intent(in)    :: cfg
      type(slow_step_data_t), intent(inout) :: step   !--- reads d_wood, npp_repro, exhausted, recruit_efficiency
      integer(ik) :: n, np, npft, j, pf, ip
      real(wp)    :: dbh_rate, dt_yr

      n     = poly%cohort%n
      np    = poly%patch%n
      npft  = poly%cohort%pft_table%n
      dt_yr = cfg%dt_years
      if (allocated(step%mortality))   deallocate(step%mortality)
      if (allocated(step%recruitment)) deallocate(step%recruitment)
      allocate(step%mortality(n), step%recruitment(npft, max(np, 1_ik)))
      step%recruitment = 0.0_wp

      associate (cohort => poly%cohort, pft => poly%cohort%pft_table, mortality => step%mortality,            &
                 recruitment => step%recruitment, npp_wood => step%d_wood, npp_repro => step%npp_repro, &
                 exhausted => step%exhausted)
         !--- The baseline seed rain, from outside the site.
         do ip = 1_ik, np
            do pf = 1_ik, npft
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%seed_rain_recruits(pf)
            end do
         end do
         do j = 1_ik, n
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !--- Prospective diameter rate and the Camac hazard (docs/science/demography.md eqs. (1), (3)).
            dbh_rate = npp_to_growth(cohort%pft_table%allom(pf), cohort%wood_carbon(j), npp_wood(j),  &
                                     cohort%dbh(j), dt_yr, pft%wood_density(pf))
            mortality(j) = camac_mortality(cohort%growth_avg(j), dbh_rate,                          &
                                           cohort%growth_count(j) > 0_ik,                           &
                                           pft%mort_gamma(pf), pft%mort_alpha(pf), pft%mort_beta(pf))
            !---~---
            !   A cohort with no leaves and no storage cannot flush again: its density goes to the
            ! floor this step (docs/science/demography.md eq. (4); plant_carbon_allocation.md §2).
            !---~---
            if (exhausted(j)) mortality(j) = max(mortality(j),                                      &
                    decline_to_floor_mortality(cohort%nplant(j), cfg%demography%negligible_nplant, dt_yr))
            !---~---
            !   Reproduction carbon to recruits, a yearly rate the recruit pool takes daily; a
            ! disabled PFT recruits nothing (docs/science/demography.md eq. (6)).
            !---~---
            if (pft%include_pft(pf) == 1_ik)                                                       &
            recruitment(pf, ip) = recruitment(pf, ip)                                              &
                 + npp_to_recruitment(cohort%nplant(j), npp_repro(j), dt_yr,                       &
                                      step%recruit_efficiency(j), pft%recruit_carbon(pf))
         end do
      end associate
   end subroutine compute_vital_rates

   !---~---
   !   Mortality litter: every pool of the plants that died this step, (n0 - n1).p1, on the pools
   ! and the density drop of the commit it follows, so litter is valued on the new pools and a
   ! cohort the floor stopped gets none (docs/science/demography.md eq. (9)). It goes to the same
   ! per-patch accumulator as the turnover litter.
   !---~---
   subroutine accumulate_mortality_litter(poly, cfg, step)
      type(polygon_t),        intent(inout) :: poly   !--- the litter accumulator and the mortality-carbon row
      type(ed_config_t),      intent(in)    :: cfg
      type(slow_step_data_t), intent(in)    :: step   !--- reads nplant_before
      integer(ik) :: j, pf, ip
      real(wp)    :: died_nplant

      associate (cohort => poly%cohort, pft => poly%cohort%pft_table, lit => poly%patch%litter_in,            &
                 nplant_before => step%nplant_before)
         do j = 1_ik, cohort%n
            died_nplant = nplant_before(j) - cohort%nplant(j)
            if (died_nplant <= 0.0_wp) cycle
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !---~---
            !   The background mortality's carbon, from the same numbers as the litter below, as
            ! a rate [kgC/m2/yr] times the step's seconds (docs/science/diagnostics.md §4).
            !---~---
            if (poly%patch%diag%active)                                                            &
               poly%patch%diag%v(PA_MORT_C_BACKGROUND, ip) =                                       &
                  poly%patch%diag%v(PA_MORT_C_BACKGROUND, ip)                                      &
                  + (died_nplant * (cohort%leaf_carbon(j) + cohort%fineroot_carbon(j)              &
                                  + cohort%wood_carbon(j) + cohort%nonstructural_carbon(j))        &
                     / cfg%dt_years) * cfg%dt_slow
            call add_necromass(lit(ip), died_nplant * cohort%leaf_carbon(j),                      &
                     died_nplant * cohort%fineroot_carbon(j), died_nplant * cohort%wood_carbon(j),  &
                     died_nplant * cohort%nonstructural_carbon(j), 0.0_wp,                         &
                     pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                                 &
                     cohort%pft_table%allom(pf)%aboveground_frac, pft%struct_lignin_frac(pf))
         end do
      end associate
   end subroutine accumulate_mortality_litter

   !---~---
   !   The turnover's tissue water. A net leaf or fine-root loss sheds water in proportion to the
   ! carbon lost, so the tissue left keeps its water per unit carbon; fine roots shed their
   ! share of the wood store. A net gain moves no water: the store keeps its mass. It runs
   ! before the commit, against the pre-step pools. The water and its enthalpy go to the ground
   ! as patch rates (shed_water_rate, shed_water_enth) the next fast window adds
   ! (docs/science/soil_plant_water.md §2; conservation.md §3).
   !---~---
   subroutine shed_turnover_water(poly, cfg, step, ledger)
      type(polygon_t),        intent(inout) :: poly
      type(ed_config_t),      intent(in)    :: cfg
      type(slow_step_data_t), intent(in)    :: step   !--- reads d_leaf, d_fineroot
      type(slow_ledger_t),    intent(inout) :: ledger
      real(wp), allocatable :: patch_total(:), patch_enth(:)
      real(wp) :: shed_frac, shed_amount
      integer(ik) :: j, ip

      allocate(patch_total(poly%patch%n), patch_enth(poly%patch%n))
      patch_total = 0.0_wp
      patch_enth  = 0.0_wp

      associate (cohort => poly%cohort)
         do j = 1_ik, cohort%n
            ip = cohort%owner_patch(j)
            if (step%d_leaf(j) < 0.0_wp) then
               shed_frac   = min(1.0_wp, -step%d_leaf(j) / max(cohort%leaf_carbon(j), tiny_num))
               shed_amount = cohort%leaf_water_mass(j) * shed_frac
               cohort%leaf_water_mass(j) = cohort%leaf_water_mass(j) - shed_amount
               patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
               patch_enth(ip)  = patch_enth(ip)  + cohort%nplant(j) * shed_amount                        &
                                                 * internal_energy_liquid(cohort%leaf_temp(j))
            end if
            if (step%d_fineroot(j) < 0.0_wp) then
               shed_frac   = min(1.0_wp, -step%d_fineroot(j) / max(cohort%fineroot_carbon(j), tiny_num))
               !--- The fine roots' own share of the store, in proportion to their carbon.
               shed_amount = max(cohort%wood_water_mass(j), 0.0_wp) * cohort%fineroot_carbon(j)   &
                           / max(cohort_wood_store_carbon(cohort%sapwood_carbon(j), cohort%fineroot_carbon(j)), &
                                 tiny_num) * shed_frac
               cohort%wood_water_mass(j) = cohort%wood_water_mass(j) - shed_amount
               patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
               patch_enth(ip)  = patch_enth(ip)  + cohort%nplant(j) * shed_amount                        &
                                                 * internal_energy_liquid(cohort%wood_temp(j))
            end if
            !--- The water leaves at its tissue's temperature: the temperature holds, the energy falls.
            call hold_cohort_tissue_temps(cohort, j)
         end do
      end associate

      do ip = 1_ik, poly%patch%n
         poly%patch%shed_water_rate(ip) = patch_total(ip) / cfg%dt_slow
         poly%patch%shed_water_enth(ip) = patch_enth(ip)  / cfg%dt_slow
      end do

      !--- The shed water and its enthalpy: a declared handoff (docs/science/conservation.md §3).
      if (ledger%active) then
         do ip = 1_ik, poly%patch%n
            call slow_ledger_declare(ledger, water_out  = patch_polygon_area(poly, ip) * patch_total(ip),     &
                                             energy_out = patch_polygon_area(poly, ip) * patch_enth(ip))
         end do
      end if
   end subroutine shed_turnover_water

   !---~---
   !   New leaves take their water from the wood: the leaf tissue the commit added is filled,
   ! with at most the water the wood holds, to the wood's potential less the crown's gravity
   ! head. The plant's water is unchanged, and a canopy re-flushing from bare does not start
   ! at the dry end of its pressure-volume curve. The wood ends a little drier than the
   ! leaves at rest; the next fast step evens them out (docs/science/soil_plant_water.md §2).
   !---~---
   subroutine fill_new_leaf_water(poly, cfg, step)
      type(polygon_t),        intent(inout) :: poly
      type(ed_config_t),      intent(in)    :: cfg
      type(slow_step_data_t), intent(in)    :: step   !--- reads leaf_carbon_before
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp)    :: grown, head, psi_wood, moved
      integer(ik) :: i

      associate (cohort => poly%cohort)
         do i = 1_ik, cohort%n
            grown = cohort%leaf_carbon(i) - step%leaf_carbon_before(i)
            if (grown <= 0.0_wp) cycle
            leaf_curve = cohort%pft_table%hydro(cohort%pft(i))%leaf_curve
            wood_curve = cohort%pft_table%hydro(cohort%pft(i))%wood_curve
            head     = grav_head * cohort%height(i)          ! the crown above the wood datum [MPa]
            psi_wood = tissue_water_to_psi(cohort%wood_water_mass(i), wood_curve,                           &
                          cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
            moved    = min(tissue_psi_to_water(psi_wood - head, leaf_curve, grown), cohort%wood_water_mass(i))
            cohort%leaf_water_mass(i) = cohort%leaf_water_mass(i) + moved
            cohort%wood_water_mass(i) = cohort%wood_water_mass(i) - moved
            call hold_cohort_tissue_temps(cohort, i)
         end do
      end associate
   end subroutine fill_new_leaf_water

   !---~---
   !   The water of the plants that die this step: their tissue water and canopy films go to the
   ! ground down the turnover shed's channel. Runs after the commit, with the dying number
   ! n0 - n1 taken as accumulate_mortality_litter takes it (docs/science/demography.md §2).
   !---~---
   subroutine shed_mortality_water(poly, cfg, step, water_shed)
      type(polygon_t),        intent(inout) :: poly
      type(ed_config_t),      intent(in)    :: cfg
      type(slow_step_data_t), intent(in)    :: step         !--- reads nplant_before
      real(wp),               intent(out)   :: water_shed   !--- [kg/m2 site]
      integer(ik) :: j, ip
      real(wp)    :: died_nplant, w_j

      water_shed = 0.0_wp
      associate (cohort => poly%cohort, patch => poly%patch)
         do j = 1_ik, cohort%n
            died_nplant = step%nplant_before(j) - cohort%nplant(j)
            if (died_nplant <= 0.0_wp) cycle
            ip  = cohort%owner_patch(j)
            !---~---
            !   The density drop times the per-plant water (cohort_tissue_water would carry the
            ! post-commit density and scale the loss by n1/n0).
            !---~---
            w_j = died_nplant * (max(cohort%leaf_water_mass(j), 0.0_wp)                            &
                               + max(cohort%wood_water_mass(j), 0.0_wp)                            &
                               + cohort%leaf_surf_water(j) + cohort%wood_surf_water(j))
            call add_shed_water(patch, ip, w_j, died_nplant * cohort_water_enthalpy(cohort, j), cfg%dt_slow)
            water_shed = water_shed + patch_polygon_area(poly, ip) * w_j
         end do
      end associate
   end subroutine shed_mortality_water

   !---~---
   !   Tissue water above the capacity a cohort's biomass now allows goes to the ground, down the
   ! turnover shed's channel (patch%shed_water_rate). An operator that shrinks a store without
   ! touching its water (an SLA rise, a fusion or split) can leave it above saturation; the
   ! caller declares the water shed (docs/science/soil_plant_water.md §2).
   !---~---
   subroutine shed_excess_tissue_water(poly, cfg, water_shed)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      real(wp),            intent(out)   :: water_shed     !--- [kg/m2 site] the water shed
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp)    :: leaf_kept, wood_kept, excess, excess_enth
      integer(ik) :: i, ip

      water_shed = 0.0_wp
      associate (cohort => poly%cohort, patch => poly%patch)
         do i = 1_ik, cohort%n
            leaf_curve = cohort%pft_table%hydro(cohort%pft(i))%leaf_curve
            wood_curve = cohort%pft_table%hydro(cohort%pft(i))%wood_curve
            leaf_kept = tissue_capped_water(cohort%leaf_water_mass(i), leaf_curve, cohort%leaf_carbon(i))
            wood_kept = tissue_capped_water(cohort%wood_water_mass(i), wood_curve,                &
                           cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
            excess = (cohort%leaf_water_mass(i) - leaf_kept) + (cohort%wood_water_mass(i) - wood_kept)
            if (excess <= 0.0_wp) cycle
            excess_enth = (cohort%leaf_water_mass(i) - leaf_kept) * internal_energy_liquid(cohort%leaf_temp(i)) &
                        + (cohort%wood_water_mass(i) - wood_kept) * internal_energy_liquid(cohort%wood_temp(i))
            cohort%leaf_water_mass(i) = leaf_kept
            cohort%wood_water_mass(i) = wood_kept
            call hold_cohort_tissue_temps(cohort, i)          ! the excess leaves at the tissue's temperature
            ip = cohort%owner_patch(i)
            call add_shed_water(patch, ip, cohort%nplant(i) * excess, cohort%nplant(i) * excess_enth,  &
                                cfg%dt_slow)
            water_shed = water_shed + patch_polygon_area(poly, ip) * cohort%nplant(i) * excess
         end do
      end associate
   end subroutine shed_excess_tissue_water

   !---~---
   !   Step each patch's soil-carbon pools one day on the step's litter,
   ! poly%patch%litter_in (turnover, mortality, and the kills of the boundary before it),
   ! into pools_new, which ed_slow_update's commit_soil_carbon
   ! writes. The ledger keeps the run's worst Rh gap, lignin residual and lambda
   ! (soil_carbon.md §6-§7, §9). Without the fast loop there is no soil climate to
   ! integrate: the pools decay at the fixed scalar spinup_xi, as GPP runs at gpp_ref, and
   ! there is no seam (soil_carbon.md §7).
   !---~---
   subroutine advance_soil_carbon(poly, cfg, ledger, pools_new)
      type(polygon_t),       intent(in)    :: poly
      type(ed_config_t),     intent(in)    :: cfg
      !---~---
      !   The slow ledger: the litter declared leaving the live pools is declared
      ! arriving here, and the day's Rh is declared leaving (conservation.md §3).
      !---~---
      type(slow_ledger_t),   intent(inout) :: ledger
      type(soil_carbon_state_t), allocatable, intent(out) :: pools_new(:)   !--- (patch) the pools after the step
      real(wp)             :: u(N_SOIL_POOL), lignin_in(2), rh_today
      real(wp)             :: xi_int(N_SOIL_POOL)   !--- [day] the scalar integrated over the step
      type(soil_carbon_out_t) :: out
      integer(ik)          :: ip

      allocate(pools_new(poly%patch%n))
      if (ledger%active) then
         do ip = 1_ik, poly%patch%n
            associate (lit => poly%patch%litter_in(ip))
               call slow_ledger_declare(ledger, carbon_in = patch_polygon_area(poly, ip)           &
                        * (lit%labile_grnd + lit%labile_soil + lit%struct_grnd + lit%struct_soil))
            end associate
         end do
      end if
      do ip = 1_ik, poly%patch%n
         call build_litter_input(poly%patch%litter_in(ip), u, lignin_in)
         associate (xa => poly%patch%decomp_xi(ip))
            if (cfg%fast_biophysics_on) then
               xi_int = xa%xi_days
            else
               xi_int = cfg%soil_carbon%spinup_xi * cfg%dt_slow / day_sec
            end if
            call soil_carbon_step(poly%patch%soil_carbon(ip), u, lignin_in, xi_int, cfg%soil_carbon,     &
                                  rh_today, out)
            pools_new(ip) = out%pools
            !---~---
            !   The assertion guard (soil_carbon.md §7): rh_today, the matrix's own
            ! debit, should equal the fast loop's accumulated Rh.
            !---~---
            out%rh_seam_gap = 0.0_wp
            if (cfg%fast_biophysics_on) out%rh_seam_gap = rh_today - xa%rh_fast
            call slow_ledger_soil_check(ledger, out, poly%patch%n)
            !---~---
            !   Declare the matrix's own Rh debit (conservation.md §3): what is left in
            ! the residual is the fast/slow seam gap.
            !---~---
            call slow_ledger_declare(ledger,                                  &
                                          carbon_out = patch_polygon_area(poly, ip) * rh_today)
         end associate
      end do
   end subroutine advance_soil_carbon

   !---~---
   !   The month's founding decision for each (PFT, patch) recruit pool, by the law of
   ! ed_demography (docs/science/demography.md §3).
   !---~---
   subroutine decide_recruitment(poly, cfg, found)
      type(polygon_t),      intent(in)  :: poly
      type(ed_config_t),    intent(in)  :: cfg
      logical, allocatable, intent(out) :: found(:,:)   !--- (pft, patch)
      integer(ik) :: np
      np = poly%patch%n
      allocate(found(poly%n_pft, max(np, 1_ik)))
      found = .false.
      found(:, 1:np) = founds_cohort(poly%patch%recruit_pool(1:poly%n_pft, 1:np), cfg%demography%min_recruit_size)
   end subroutine decide_recruitment

   !---~---
   !   The year's treefall decision, by the laws of ed_demography: the share of every
   ! patch's area the gaps take, and which cohorts survive into them
   ! (docs/science/demography.md §4, eq. (7)).
   !---~---
   subroutine decide_disturbance(poly, cfg, dt_yr, frac, survives)
      type(polygon_t),      intent(in)  :: poly
      type(ed_config_t),    intent(in)  :: cfg
      real(wp),             intent(in)  :: dt_yr      !--- [yr] the interval the yearly rate covers
      real(wp),             intent(out) :: frac       !--- [-] the disturbed share
      logical, allocatable, intent(out) :: survives(:)   !--- (cohort)
      integer(ik) :: n
      frac = 0.0_wp
      if (cfg%demography%patch_disturbance_rate > 0.0_wp) frac = disturbed_fraction(cfg%demography%patch_disturbance_rate, dt_yr)
      n = poly%cohort%n
      allocate(survives(max(n, 1_ik)))
      survives = .true.
      survives(1:n) = survives_treefall(poly%cohort%height(1:n), cfg%demography%disturbance_survive_height)
   end subroutine decide_disturbance

   !---~---
   !   An empty kill report with room for `capacity` rows: an operator kills each cohort
   ! at most once.
   !---~---
   subroutine open_kill_report(kills, capacity)
      type(kill_report_t), intent(out) :: kills
      integer(ik),         intent(in)  :: capacity
      allocate(kills%patch_id(capacity), kills%pft(capacity), kills%leaf(capacity),             &
               kills%fineroot(capacity), kills%wood(capacity), kills%storage(capacity))
      kills%n = 0_ik
   end subroutine open_kill_report

   !--- One killed cohort's row of the report [kgC/m2 of the patch `patch_id`].
   subroutine report_kill(kills, patch_id, pft, leaf, fineroot, wood, storage)
      type(kill_report_t), intent(inout) :: kills
      integer(ik),         intent(in)    :: patch_id, pft
      real(wp),            intent(in)    :: leaf, fineroot, wood, storage
      kills%n = kills%n + 1_ik
      kills%patch_id(kills%n) = patch_id
      kills%pft(kills%n)      = pft
      kills%leaf(kills%n)     = leaf
      kills%fineroot(kills%n) = fineroot
      kills%wood(kills%n)     = wood
      kills%storage(kills%n)  = storage
   end subroutine report_kill

   !---~---
   !   The killed plants' carbon into their patches' litter, split as background mortality's
   ! is (soil_carbon.md §5); the next soil step takes it up. `carbon` [kgC/m2 polygon] is
   ! what left the plants, for the ledger to declare (conservation.md §3).
   !---~---
   subroutine kill_litter(poly, kills, carbon)
      type(polygon_t),     intent(inout) :: poly
      type(kill_report_t), intent(in)    :: kills
      real(wp),            intent(out)   :: carbon
      integer(ik) :: k, ip, pf
      carbon = 0.0_wp
      associate (pft => poly%cohort%pft_table)
         do k = 1_ik, kills%n
            ip = findloc(poly%patch%global_id(1:poly%patch%n), kills%patch_id(k), dim=1)
            if (ip < 1_ik) error stop 'kill_litter: the patch a kill was reported on is gone'
            pf = kills%pft(k)
            call add_necromass(poly%patch%litter_in(ip), kills%leaf(k), kills%fineroot(k),           &
                     kills%wood(k), kills%storage(k), 0.0_wp, pft%f_labile_leaf(pf),                 &
                     pft%f_labile_stem(pf), pft%allom(pf)%aboveground_frac, pft%struct_lignin_frac(pf))
            carbon = carbon + patch_polygon_area(poly, ip)                                          &
                              * (kills%leaf(k) + kills%fineroot(k) + kills%wood(k) + kills%storage(k))
         end do
      end associate
   end subroutine kill_litter

   !---~---
   !   The step's litter has gone to the soil pools or, with soil carbon off, left the
   ! model: the accumulator starts the next step empty (soil_carbon.md §5).
   !---~---
   subroutine clear_litter(poly)
      type(polygon_t), intent(inout) :: poly
      poly%patch%litter_in(1:poly%patch%n) = litter_input_t()
   end subroutine clear_litter

end module ed_slow_state
