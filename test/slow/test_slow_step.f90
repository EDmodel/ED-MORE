! SPDX-License-Identifier: Apache-2.0
!---~---
!   Leaf phenology's wiring (ed_slow_state's advance_leaf_phenology): a two-cohort site
! through one synthetic Ithaca year; a step with no air temperature is skipped; each cue
! driver against a hand-computed value; a drought-deciduous and a light-exchanging PFT end
! to end.
!---~---
module test_slow_step_phenology
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_pft_derive,           only : pft_leaf_psi_tlp
   use ed_polygon_state_vars,   only : polygon_t
   use ed_init_bare_ground,     only : init_bare_ground, add_cohort
   use ed_slow_state,           only : advance_leaf_phenology
   use ed_constants,            only : day_sec
   use ed_time,                 only : daylength
   use ed_pft_params,           only : CUE_TEMP, CUE_NONE, CUE_WATER, CUE_LIGHT
   use ed_test_support,         only : build_test_config, check_close, check_int, check_true,            &
                                       test_report
   implicit none
   private

   public :: run_phenology_wiring

   real(wp), parameter :: TWOPI = 6.283185307179586_wp
   type(ed_config_t)   :: cfg
   type(polygon_t)     :: poly
   integer(ik) :: doy
   real(wp)    :: fl_decid_200, sh_decid_200, fl_decid_340, sh_decid_340
   real(wp)    :: fl_none_200, sh_none_200, fl_none_340, sh_none_340, gdd_summer

contains

   subroutine run_phenology_wiring()
      !---~---
      !   Config: PFT 1 temperate deciduous (warmth x light, cold x short days), the rest
      ! with no cues. This test calls advance_leaf_phenology directly.
      !---~---
      cfg = build_test_config()
      cfg%forcing%latitude               = 42.44_wp           ! Ithaca NY (northern hemisphere)
      cfg%pft%pheno(1)%flush_cue_mask        = CUE_TEMP + CUE_LIGHT
      cfg%pft%pheno(1)%shed_cue_mask         = CUE_TEMP + CUE_LIGHT
      cfg%pft%pheno(1)%flush_degree_days     = 92.0_wp
      cfg%pft%pheno(1)%shed_base_temp        = 290.37_wp
      cfg%pft%pheno(1)%shed_degree_days      = 48.0_wp
      cfg%pft%pheno(1)%flush_light_hours     = 10.35_wp
      cfg%pft%pheno(1)%shed_light_hours      = 9.83_wp

      !--- A site with two cohorts: cohort 1 = PFT 1 (deciduous), cohort 2 = PFT 2 (no cues).
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 10.0_wp)
      call check_int('two cohorts created', int(poly%cohort%n, ik), 2_ik)
      call check_close('cohort 1 born flushing',      poly%cohort%leaf_flush_tendency(1), 1.0_wp, 1.0e-9_wp)
      call check_close('cohort 1 born not senescing', poly%cohort%leaf_shed_tendency(1),  0.0_wp, 1.0e-9_wp)

      !--- Drive one synthetic year (dt_slow = 1 day).
      do doy = 1_ik, 365_ik
         call set_daily_drivers(poly, daily_tair(doy), daylength(cfg%forcing%latitude, doy))
         call advance_leaf_phenology(poly, cfg, doy)
         if (doy == 200_ik) then
            fl_decid_200 = poly%cohort%leaf_flush_tendency(1)
            sh_decid_200 = poly%cohort%leaf_shed_tendency(1)
            fl_none_200  = poly%cohort%leaf_flush_tendency(2)
            sh_none_200  = poly%cohort%leaf_shed_tendency(2)
            gdd_summer   = poly%cohort%growing_degree_days(1)
         end if
         if (doy == 340_ik) then
            fl_decid_340 = poly%cohort%leaf_flush_tendency(1)
            sh_decid_340 = poly%cohort%leaf_shed_tendency(1)
            fl_none_340  = poly%cohort%leaf_flush_tendency(2)
            sh_none_340  = poly%cohort%leaf_shed_tendency(2)
         end if
      end do


      !--- 3. A no-temperature step is skipped (tendencies + memory unchanged).
      block
         real(wp) :: fl_before, sh_before, gdd_before
         fl_before  = poly%cohort%leaf_flush_tendency(1)
         sh_before  = poly%cohort%leaf_shed_tendency(1)
         gdd_before = poly%cohort%growing_degree_days(1)
         call set_daily_drivers(poly, 0.0_wp, 0.0_wp)
         poly%site%pheno_tair_n(1) = 0_ik                            ! no fast sub-steps ran
         call advance_leaf_phenology(poly, cfg, 1_ik)
         call check_true('no-temperature step leaves the flush tendency unchanged', &
                         abs(poly%cohort%leaf_flush_tendency(1) - fl_before) < tiny(1.0_wp))
         call check_true('no-temperature step leaves the shed tendency unchanged', &
                         abs(poly%cohort%leaf_shed_tendency(1) - sh_before) < tiny(1.0_wp))
         call check_true('no-temperature step leaves the warmth sum unchanged', &
                         abs(poly%cohort%growing_degree_days(1) - gdd_before) < tiny(1.0_wp))
      end block

      !--- 4. THE CUE DRIVERS, each against a hand-computed value.

      !---~---
      !   4a. Temperature: 10 days at 288.15 K add 10 x (288.15 - 278.15) to the warmth sum;
      ! before midsummer the cold sum stays 0, after it 10 days at 280.15 K add
      ! 10 x (290.37 - 280.15).
      !---~---
      block
         call reset_pheno_memory(poly)
         do doy = 100_ik, 109_ik
            call set_daily_drivers(poly, 288.15_wp, 12.0_wp)
            call advance_leaf_phenology(poly, cfg, doy)
         end do
         call check_close('warmth sum after 10 days 10 K above base', poly%cohort%growing_degree_days(1), &
                          100.0_wp, 1.0e-9_wp)
         call check_close('no cold sum before midsummer', poly%cohort%cold_degree_days(1), 0.0_wp, 0.0_wp)
         call reset_pheno_memory(poly)
         do doy = 250_ik, 259_ik
            call set_daily_drivers(poly, 280.15_wp, 12.0_wp)
            call advance_leaf_phenology(poly, cfg, doy)
         end do
         call check_close('cold sum after 10 days below the shed base', poly%cohort%cold_degree_days(1), &
                          10.0_wp * (290.37_wp - 280.15_wp), 1.0e-9_wp)
      end block

      !---~---
      !   4b. Water: the dry sum adds (psi_tlp - psi_pd) per day, against the PFT's DERIVED
      ! turgor-loss point; one wet day adds wet credit without wiping it.
      !---~---
      block
         real(wp) :: tlp, dry7
         poly%cohort%pft_table%pheno(1)%flush_cue_mask = CUE_NONE
         poly%cohort%pft_table%pheno(1)%shed_cue_mask  = CUE_WATER
         tlp = pft_leaf_psi_tlp(cfg%pft, 1_ik)
         call reset_pheno_memory(poly)
         poly%cohort%dmax_psi_leaf(1) = tlp - 1.0_wp           ! 1 MPa below the TLP
         do doy = 1_ik, 7_ik
            call set_daily_drivers(poly, 290.0_wp, 12.0_wp)
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         dry7 = poly%cohort%dry_psi_sum(1)
         call check_close('dry sum after 7 days 1 MPa below the TLP', dry7, 7.0_wp, 1.0e-12_wp)
         poly%cohort%dmax_psi_leaf(1) = tlp + 0.5_wp           ! one wet day
         call set_daily_drivers(poly, 290.0_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, 200_ik)
         call check_close('one wet day adds wet credit', poly%cohort%wet_psi_sum(1), 0.5_wp, 1.0e-12_wp)
         call check_close('one wet day does not wipe the dry sum', poly%cohort%dry_psi_sum(1), dry7, 0.0_wp)
      end block

      !---~---
      !   4c. Light: x += w*(h - x) with w = dt/window. From x = 0 with a constant input h
      ! and w = 0.1, after n days x = h*(1 - 0.9^n) EXACTLY. Each cohort reads its OWN
      ! hours of light (light_hours_accum, counted at its top): a shaded cohort's mean
      ! stays lower. A cohort with no light memory starts from its first day's hours.
      !---~---
      block
         real(wp) :: expect
         poly%cohort%pft_table%pheno(1:2)%shed_cue_mask = CUE_LIGHT
         poly%cohort%pft_table%pheno(1:2)%light_window  = 10.0_wp
         call reset_pheno_memory(poly)
         call set_daily_drivers(poly, 290.0_wp, 9.0_wp)
         call advance_leaf_phenology(poly, cfg, 200_ik)
         call check_close('a cohort with no light memory starts from its first day', &
                          poly%cohort%light_hours_mean(1), 9.0_wp, 0.0_wp)
         poly%cohort%light_hours_mean(1:2) = 0.0_wp
         do doy = 1_ik, 5_ik
            call set_daily_drivers(poly, 290.0_wp, 12.0_wp)
            poly%cohort%light_hours_accum(2) = 3.0_wp * cfg%dt_slow / day_sec    ! cohort 2 in shade
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         expect = 12.0_wp * (1.0_wp - 0.9_wp**5)
         call check_close('light running mean after 5 days of 12 h',                                &
                          poly%cohort%light_hours_mean(1), expect, 1.0e-9_wp)
         call check_close('a shaded cohort reads its own hours of light (3 h)',                     &
                          poly%cohort%light_hours_mean(2), expect / 4.0_wp, 1.0e-9_wp)
         poly%cohort%pft_table%pheno(2)%shed_cue_mask = CUE_NONE
      end block

      !--- 5. ACCEPTANCE: the drought-deciduous and light-exchanging habits, end to end.

      !---~---
      !   Drought deciduous (flush and shed both on water): full when watered, senescing under
      ! sustained drought, flushing again on rewet.
      !---~---
      block
         real(wp) :: shed_wet, shed_dry, flush_dry, flush_rewet, tlp
         integer(ik) :: d
         poly%cohort%pft_table%pheno(1)%flush_cue_mask = CUE_WATER
         poly%cohort%pft_table%pheno(1)%shed_cue_mask  = CUE_WATER
         tlp = pft_leaf_psi_tlp(cfg%pft, 1_ik)
         call reset_pheno_memory(poly)
         poly%cohort%dmax_psi_leaf(1) = 0.5_wp * tlp               ! well watered
         do d = 1_ik, 60_ik
            call set_daily_drivers(poly, 295.0_wp, 12.0_wp)
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         shed_wet = poly%cohort%leaf_shed_tendency(1)
         poly%cohort%dmax_psi_leaf(1) = tlp - 2.0_wp               ! sustained drought, past the TLP
         do d = 1_ik, 30_ik
            call set_daily_drivers(poly, 295.0_wp, 12.0_wp)
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         shed_dry  = poly%cohort%leaf_shed_tendency(1)
         flush_dry = poly%cohort%leaf_flush_tendency(1)
         poly%cohort%dmax_psi_leaf(1) = 0.5_wp * tlp               ! rewet
         do d = 1_ik, 30_ik
            call set_daily_drivers(poly, 295.0_wp, 12.0_wp)
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         flush_rewet = poly%cohort%leaf_flush_tendency(1)
         call check_true('drought-deciduous: no senescence when watered',  shed_wet < 0.05_wp,   shed_wet)
         call check_true('drought-deciduous: senesces under drought',      shed_dry > 0.8_wp,    shed_dry)
         call check_true('drought-deciduous: flush suppressed in drought', flush_dry < 0.2_wp,   flush_dry)
         call check_true('drought-deciduous: flushes again on rewet',      flush_rewet > 0.8_wp, flush_rewet)
      end block

      !--- Light-driven leaf exchange: flushing throughout, senescence rising with light.
      block
         real(wp) :: shed_dim, shed_bright, flush_bright
         integer(ik) :: d
         poly%cohort%pft_table%pheno(1)%flush_cue_mask       = CUE_NONE
         poly%cohort%pft_table%pheno(1)%shed_cue_mask        = CUE_LIGHT
         poly%cohort%pft_table%pheno(1)%shed_light_hours     = 6.0_wp
         poly%cohort%pft_table%pheno(1)%shed_light_sharpness = 2.0_wp         ! > 0: many bright hours trigger
         call reset_pheno_memory(poly)
         do d = 1_ik, 60_ik
            call set_daily_drivers(poly, 295.0_wp, 2.0_wp)       ! few bright hours
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         shed_dim = poly%cohort%leaf_shed_tendency(1)
         call reset_pheno_memory(poly)
         do d = 1_ik, 60_ik
            call set_daily_drivers(poly, 295.0_wp, 10.0_wp)      ! many bright hours
            call advance_leaf_phenology(poly, cfg, 200_ik)
         end do
         shed_bright  = poly%cohort%leaf_shed_tendency(1)
         flush_bright = poly%cohort%leaf_flush_tendency(1)
         call check_true('leaf exchange: little senescence under dim light', shed_dim < 0.05_wp,     shed_dim)
         call check_true('leaf exchange: senescence rises with light',       shed_bright > 0.8_wp,   shed_bright)
         call check_true('leaf exchange: still flushing while exchanging',   flush_bright > 0.99_wp, flush_bright)
      end block
   end subroutine run_phenology_wiring

   !---~---
   !   Fill what the fast loop fills: the site air-temperature sum and, per cohort, the time
   ! [h] its top's PAR exceeded par_min over the step -- here `hours` a day at every cohort,
   ! so the driver's per-day value recovers it exactly.
   !---~---
   subroutine set_daily_drivers(poly, tair, hours)
      type(polygon_t), intent(inout) :: poly
      real(wp),     intent(in)    :: tair, hours
      poly%site%pheno_tair_sum(1) = tair
      poly%site%pheno_tair_n(1) = 1_ik
      poly%cohort%light_hours_accum(1:poly%cohort%n) = hours * cfg%dt_slow / day_sec
   end subroutine set_daily_drivers
   !--- Clear every phenology memory so each cue block starts from a known state.
   subroutine reset_pheno_memory(poly)
      type(polygon_t), intent(inout) :: poly
      associate (n => poly%cohort%n)
         poly%cohort%leaf_flush_tendency(1:n) = 1.0_wp
         poly%cohort%leaf_shed_tendency(1:n)  = 0.0_wp
         poly%cohort%growing_degree_days(1:n) = 0.0_wp
         poly%cohort%cold_degree_days(1:n)    = 0.0_wp
         poly%cohort%dry_psi_sum(1:n)         = 0.0_wp
         poly%cohort%wet_psi_sum(1:n)         = 0.0_wp
         poly%cohort%light_hours_mean(1:n)    = -1.0_wp      ! no light memory
      end associate
   end subroutine reset_pheno_memory
   !--- Ithaca-like daily-mean air temperature [K]: ~270 in winter, ~297 in summer.
   pure real(wp) function daily_tair(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 14.0_wp * sin(TWOPI * (real(doy, wp) - 110.0_wp) / 365.0_wp)
   end function daily_tair

end module test_slow_step_phenology

!---~---
!   Thermal acclimation's wiring (ed_slow_state's acclimate_site_leaf_tables, and its
! switch in advance_one_step; docs/science/leaf_gas_exchange.md §2): each site's growth
! temperature starts from its own first day and follows its own days as a running mean,
! and each site's leaf table is the PFT table acclimated at its own mean; a site whose
! fast loop did not run keeps both; with acclimation off every table stays the PFT's own.
!---~---
module test_slow_step_acclimation
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : t_kelvin
   use ed_config,               only : ed_config_t
   use ed_pft_params,           only : leaf_photo_table_t
   use ed_polygon_state_vars,   only : polygon_t
   use ed_state_defaults,       only : is_unset_positive
   use ed_init_bare_ground,     only : init_bare_sites, finalize_init
   use ed_plant_traits,         only : acclimate_leaf_photo_table
   use ed_slow_state,           only : acclimate_site_leaf_tables
   use ed_stepper,              only : advance_one_step, run_ledgers_t
   use ed_test_support,         only : banner, build_test_config, check_close, check_true
   implicit none
   private

   public :: run_thermal_acclimation

   real(wp), parameter :: T_COOL = t_kelvin + 10.0_wp   ! [K] the cool site's daily mean
   real(wp), parameter :: T_WARM = t_kelvin + 30.0_wp   ! [K] the warm site's
   real(wp), parameter :: WINDOW = 30.0_wp              ! [day] the running mean's window

contains

   subroutine run_thermal_acclimation()
      type(ed_config_t)        :: cfg
      type(polygon_t)          :: poly
      type(run_ledgers_t)      :: ledgers   ! never activated: the steps book nothing
      type(leaf_photo_table_t) :: warm_table
      real(wp)    :: t_new, gap0
      integer(ik) :: iday, is
      call banner('thermal acclimation: each site''s growth temperature and leaf table')
      cfg = build_test_config()
      cfg%leaf_physiology%thermal_acclimation = .true.
      cfg%leaf_physiology%acclim_window_days  = WINDOW

      !---~---
      !   1. A first day, 10 C at one site and 30 C at the other: each unset mean takes its own
      ! site's day, and each site's table is the PFT table acclimated there, so the warm site's
      ! entropy term is the lower.
      !---~---
      call two_sites(cfg, poly)
      call set_day(poly, [T_COOL, T_WARM])
      call acclimate_site_leaf_tables(poly, cfg%leaf_physiology, 1_ik)
      call check_close('an unset mean takes the cool site''s first day', poly%site%t_growth_avg(1), T_COOL, 0.0_wp)
      call check_close('an unset mean takes the warm site''s first day', poly%site%t_growth_avg(2), T_WARM, 0.0_wp)
      call check_tables_follow(poly, cfg, 'first day')
      call check_true('the warm site''s Vcmax entropy term is the lower',                                    &
                      poly%site%leaf_photo(2)%pft(1)%ds_vcmax < poly%site%leaf_photo(1)%pft(1)%ds_vcmax,     &
                      poly%site%leaf_photo(1)%pft(1)%ds_vcmax - poly%site%leaf_photo(2)%pft(1)%ds_vcmax)

      !---~---
      !   2. Thirty days at 20 C at the cool site, on a 30-day window: each day closes 1/30 of the
      ! gap to the day's temperature, so (29/30)**30 of the first gap is left. The warm site's
      ! fast loop runs no sub-step: its mean and its table stay.
      !---~---
      t_new = t_kelvin + 20.0_wp
      gap0  = poly%site%t_growth_avg(1) - t_new
      warm_table = poly%site%leaf_photo(2)
      do iday = 1_ik, 30_ik
         call set_day(poly, [t_new, T_COOL])
         poly%site%pheno_tair_n(2) = 0_ik
         call acclimate_site_leaf_tables(poly, cfg%leaf_physiology, 1_ik)
      end do
      call check_close('a month later (29/30)**30 of the gap is left', poly%site%t_growth_avg(1) - t_new,       &
                       gap0 * (1.0_wp - 1.0_wp / WINDOW)**30, 1.0e-9_wp)
      call check_tables_follow(poly, cfg, 'a month later')
      call check_close('a site with no sub-step keeps its mean', poly%site%t_growth_avg(2), T_WARM, 0.0_wp)
      call check_true('and its table', same_acclimated(poly%site%leaf_photo(2), warm_table))

      !--- 3. A step longer than the window takes the day's mean outright.
      call set_day(poly, [T_WARM, T_WARM])
      call acclimate_site_leaf_tables(poly, cfg%leaf_physiology, 2_ik * nint(WINDOW, ik))
      call check_close('a step longer than the window takes the day''s mean', poly%site%t_growth_avg(1), T_WARM,  &
                       1.0e-12_wp)

      !---~---
      !   4. Through advance_one_step, after a fast loop that saw 10 C and 30 C: with acclimation
      ! off every site's mean stays unset and its table the PFT's own; on, the same step seeds
      ! the means and acclimates the tables.
      !---~---
      cfg%leaf_physiology%thermal_acclimation = .false.
      call two_sites(cfg, poly)
      call set_day(poly, [T_COOL, T_WARM])
      call advance_one_step(poly, cfg, ledgers)
      do is = 1_ik, poly%site%n
         call check_true('acclimation off: the mean stays unset', is_unset_positive(poly%site%t_growth_avg(is)),  &
                         poly%site%t_growth_avg(is))
         call check_true('acclimation off: the table is the PFT''s own',                                          &
                         same_acclimated(poly%site%leaf_photo(is), cfg%pft%leaf_photo))
      end do
      cfg%leaf_physiology%thermal_acclimation = .true.
      call two_sites(cfg, poly)
      call set_day(poly, [T_COOL, T_WARM])
      call advance_one_step(poly, cfg, ledgers)
      call check_close('acclimation on: the step seeds the cool site''s mean', poly%site%t_growth_avg(1), T_COOL,  &
                       0.0_wp)
      call check_tables_follow(poly, cfg, 'through the step')
   end subroutine run_thermal_acclimation

   !--- A polygon of two bare sites, each half of it, on the run's config.
   subroutine two_sites(cfg, poly)
      type(ed_config_t), intent(in)    :: cfg
      type(polygon_t),     intent(out) :: poly
      call init_bare_sites(poly, cfg, [1_ik, 1_ik], [cfg, cfg], [0.5_wp, 0.5_wp])
      call finalize_init(poly)
   end subroutine two_sites

   !--- Each site's daily-mean air temperature [K], as the last fast loop's one sub-step.
   subroutine set_day(poly, t_air)
      type(polygon_t), intent(inout) :: poly
      real(wp),        intent(in)    :: t_air(:)
      poly%site%pheno_tair_sum(1:2) = t_air
      poly%site%pheno_tair_n(1:2)   = 1_ik
   end subroutine set_day

   !--- Each site's table is the PFT table acclimated at that site's running mean.
   subroutine check_tables_follow(poly, cfg, tag)
      type(polygon_t),     intent(in) :: poly
      type(ed_config_t),   intent(in) :: cfg
      character(len=*),    intent(in) :: tag
      type(leaf_photo_table_t) :: expected
      integer(ik) :: is
      do is = 1_ik, poly%site%n
         expected = cfg%pft%leaf_photo
         associate (lp => cfg%leaf_physiology)
            call acclimate_leaf_photo_table(poly%site%t_growth_avg(is), lp%acclim_ds_vcmax_a, lp%acclim_ds_vcmax_b, &
                                            lp%acclim_ds_jmax_a, lp%acclim_ds_jmax_b, lp%acclim_jv_a,               &
                                            lp%acclim_jv_b, expected)
         end associate
         call check_true(tag//': the site''s table is acclimated at its own mean',                               &
                         same_acclimated(poly%site%leaf_photo(is), expected), poly%site%t_growth_avg(is))
      end do
   end subroutine check_tables_follow

   !--- Whether two tables hold the same entropy terms and Jmax:Vcmax ratios, every PFT.
   pure logical function same_acclimated(a, b)
      type(leaf_photo_table_t), intent(in) :: a, b
      same_acclimated = all(a%pft%ds_vcmax == b%pft%ds_vcmax) .and. all(a%pft%ds_jmax == b%pft%ds_jmax)       &
                        .and. all(a%jmax_vcmax_ratio == b%jmax_vcmax_ratio)
   end function same_acclimated

end module test_slow_step_acclimation
!---~---
!   The carbon step: a step without the fast loop grows wood and diameter, the driver's
! guards see the fast loop's states, the cached wood geometry follows the growth, a shaded
! cohort acclimates and resorbs, its fine roots turn over with its leaves when
! fineroot_turnover_plastic says so, the turnover sheds its water, a re-flushing cohort
! fills its new leaves on its own PFT's water curves, and a restart's re-acclimation takes
! only the storage there is.
!---~---
module test_slow_step_carbon
   use ed_kinds,                    only : wp, ik
   use ed_constants,                only : grav_head
   use ed_therm_lib,                only : internal_energy_liquid
   use ed_plant_hydr_lib,           only : water_curve_t, tissue_water_to_psi, tissue_psi_to_water
   use ed_allometry_lib, only : dbh_to_wai, sapwood_fraction, carbon_to_structure
   use ed_config,                 only : ed_config_t, growth_window_steps
   use ed_cohort_state_vars,   only : cohort_deriv_alloc, cohort_bad_field, cohort_field_finite, cohort_field_name
   use ed_patch_state_vars,    only : canopy_air_bad_field
   use ed_polygon_state_vars, only : polygon_t
   use ed_slow_update, only : update_cohort_states, terminate_cohorts, update_cohort_derivatives, commit_trait_carbon
   use ed_plant_carbon_allocation, only : plant_carbon_allocation
   use ed_slow_state,   only : kill_report_t, advance_plant_traits, shed_turnover_water, reacclimate_plant_traits,      &
                               slow_step_data_t, compute_carbon_allocation, clear_litter, fill_new_leaf_water
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_stepper,                only : advance_one_step, run_ledgers_t
   use ed_test_support, only : banner, build_test_config, check, check_close
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use ed_slow_ledger,  only : slow_ledger_t
   implicit none
   private

   public :: run_carbon_step
   real(wp), allocatable :: leaf_new(:), storage_new(:)   ! the traits' leaf and storage carbon
   type(slow_ledger_t) :: no_ledger    ! inactive: a test step books nothing
   type(run_ledgers_t) :: ledgers      ! never activated: the coupled steps book nothing

   type(ed_config_t)   :: cfg
   type(polygon_t)     :: poly
   real(wp)    :: wood0, dbh0
   integer(ik) :: istep

contains

   subroutine run_carbon_step()
      call banner('the carbon step: growth, its geometry, plasticity and the turnover shed')
      cfg = build_test_config()

      !--- 2. A step without the fast loop grows wood_carbon -> dbh; leaf_area stays leaf_carbon*sla.
      cfg%gpp_ref = 0.5_wp
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 3_ik, 0.1_wp, 20.0_wp)       ! climax cohort
      call finalize_init(poly)
      wood0 = poly%cohort%wood_carbon(1)
      dbh0 = poly%cohort%dbh(1)
      do istep = 1_ik, 30_ik
         call advance_one_step(poly, cfg, ledgers)                  ! growth + mortality only
      end do
      call check(poly%cohort%wood_carbon(1) > wood0, 'carbon growth did not add wood_carbon')
      call check(poly%cohort%dbh(1) > dbh0,          'dbh did not increase with wood_carbon (the flip)')
      call check_close(poly%cohort%leaf_area(1), poly%cohort%leaf_carbon(1) * cfg%pft%sla(3_ik),   &
                       1.0e-9_wp, 'leaf_area /= leaf_carbon*sla after a slow-only step')
      call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'a slow-only step left the cohort out of bounds')

      !---~---
      !   The driver's checks see the fast loop's states too: a non-finite canopy air or tissue
      ! energy is what the integrators commit when even their smallest step fails.
      !---~---
      block
         real(wp)    :: keep
         integer(ik) :: k
         if (allocated(poly%patch%cas)) then
            keep = poly%patch%cas(1)%can_shv
            poly%patch%cas(1)%can_shv = ieee_value(keep, ieee_quiet_nan)
            call check(canopy_air_bad_field(poly%patch%cas(1)) == 2_ik, 'a NaN canopy-air humidity is caught')
            poly%patch%cas(1)%can_shv = keep
            call check(canopy_air_bad_field(poly%patch%cas(1)) == 0_ik, 'the restored canopy air passes')
         end if
         keep = poly%cohort%leaf_energy(1)
         poly%cohort%leaf_energy(1) = ieee_value(keep, ieee_positive_inf)
         k = cohort_bad_field(poly%cohort, 1_ik)
         call check(trim(cohort_field_name(k)) == 'leaf_energy' .and. .not. cohort_field_finite(poly%cohort, 1_ik, k), &
                    'an infinite leaf energy is caught, as not finite')
         poly%cohort%leaf_energy(1) = keep
         call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'the restored cohort passes')
      end block

      !---~---
      !   2b. The cached WOOD geometry must describe TODAY's tree, not the one 30 steps ago.
      ! wood_area / sapwood_carbon / sapwood_area are derived once per size change rather
      ! than recomputed every dt_fast, so the fast loop reads them in place. That is only
      ! safe while every path that moves dbh, basal_area or wood_carbon re-derives them,
      ! growth included. Assert it directly: a stale cache is invisible in a conservation
      ! ledger and shows up only as a slowly wrong tree.
      !---~---
      block
         real(wp) :: wa_fresh, fs_fresh
         associate (pft => poly%cohort%pft_table, pf => poly%cohort%pft(1))
            wa_fresh = dbh_to_wai(pft%allom(pf), poly%cohort%dbh(1), 1.0_wp)
            fs_fresh = sapwood_fraction(pft%allom(pf), poly%cohort%dbh(1))
         end associate
         call check(wa_fresh > 0.0_wp, 'wood-geometry fixture is degenerate (zero wood area)')
         call check_close(poly%cohort%wood_area(1), wa_fresh, 1.0e-13_wp,                           &
                          'cached wood_area is stale after growth')
         call check_close(poly%cohort%sapwood_carbon(1), fs_fresh * poly%cohort%wood_carbon(1),     &
                          1.0e-13_wp, 'cached sapwood_carbon is stale after growth')
         call check_close(poly%cohort%sapwood_area(1), fs_fresh * poly%cohort%basal_area(1),        &
                          1.0e-13_wp, 'cached sapwood_area is stale after growth')
      end block

      !--- 3. Trait plasticity: a SHADED cohort acclimates -- SLA up, Vcmax down, turnover down.
      block
         real(wp) :: sla0, vc0, ll0, ns0
         cfg%trait_dynamics%trait_plasticity_on = .true.
         poly%cohort%overtopping_lai(1) = 4.0_wp                  ! deep shade
         sla0 = poly%cohort%sla(1)
         vc0 = poly%cohort%vcmax25(1)
         ll0 = poly%cohort%leaf_turnover_rate(1)
         ns0  = poly%cohort%nonstructural_carbon(1)
         do istep = 1_ik, 200_ik                                  ! many steps so the slow relaxation shows
            call advance_plant_traits(poly, 1.0_wp, leaf_new, storage_new)
            call commit_trait_carbon(poly, leaf_new, storage_new)
         end do
         call check(poly%cohort%sla(1)     > sla0, 'shaded cohort: SLA acclimates upward')
         call check(poly%cohort%vcmax25(1) < vc0,  'shaded cohort: Vcmax acclimates downward')
         !---~---
         !   PFT 3 (climax) turns its leaves over at 0.33 /yr, below the ~0.39 /yr ED2 crossover, so
         ! kplastic_leaf_turnover is capped at 0 (docs/science/plant_traits.md eq. (4)) => a
         ! slow-turnover PFT keeps ONE rate through the canopy (never faster in shade).
         !---~---
         call check_close(poly%cohort%leaf_turnover_rate(1), ll0, 1.0e-12_wp,                       &
                          'slow-turnover PFT: leaf turnover rate is canopy-invariant (capped)')
         !--- SLA rose => the leaf pool overshot allometry; the excess was resorbed to storage.
         call check(poly%cohort%nonstructural_carbon(1) > ns0, 'SLA overshoot resorbed to storage')
         call check_close(poly%cohort%leaf_area(1), poly%cohort%leaf_carbon(1) * poly%cohort%sla(1),   &
                          1.0e-9_wp, 'leaf_area stays leaf_carbon*sla after resorption')
      end block

      !---~---
      !   3b. Plastic fine-root turnover (compute_carbon_allocation; docs/science/plant_traits.md
      ! §4). With fineroot_turnover_plastic the fine roots turn over at the PFT's rate times the
      ! cohort's leaf turnover over its top-of-canopy value; without it at the PFT's own rate,
      ! whatever the leaves do. Read from the fine-root litter one carbon step sheds, for a cohort
      ! at its top-of-canopy leaf turnover and at half of it.
      !---~---
      block
         real(wp) :: gamma_toc, root_off_half, root_off_toc, root_on_half, root_on_toc
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
         call finalize_init(poly)
         gamma_toc = poly%cohort%pft_table%leaf_turnover_rate_toc(poly%cohort%pft(1))
         cfg%trait_dynamics%trait_plasticity_on = .true.

         cfg%trait_dynamics%fineroot_turnover_plastic = .false.
         poly%cohort%leaf_turnover_rate(1) = gamma_toc
         root_off_toc  = fineroot_litter(poly, cfg)
         poly%cohort%leaf_turnover_rate(1) = 0.5_wp * gamma_toc
         root_off_half = fineroot_litter(poly, cfg)

         cfg%trait_dynamics%fineroot_turnover_plastic = .true.
         poly%cohort%leaf_turnover_rate(1) = gamma_toc
         root_on_toc  = fineroot_litter(poly, cfg)
         poly%cohort%leaf_turnover_rate(1) = 0.5_wp * gamma_toc
         root_on_half = fineroot_litter(poly, cfg)
         cfg%trait_dynamics%fineroot_turnover_plastic = .false.

         call check(root_off_toc > 0.0_wp, 'the fixture''s fine roots shed no litter')
         call check_close(root_off_half, root_off_toc, 1.0e-14_wp,                                        &
                          'not plastic: the fine roots turn over at the PFT''s rate, whatever the leaves do')
         call check_close(root_on_toc, root_off_toc, 1.0e-14_wp,                                          &
                          'plastic, at the top-of-canopy leaf turnover: the PFT''s own fine-root rate')
         call check_close(root_on_half, 0.5_wp * root_off_toc, 1.0e-12_wp,                                &
                          'plastic, at half the top-of-canopy leaf turnover: half the fine-root litter')
      end block

      !---~---
      !   4. shed_turnover_water: a NET leaf/root LOSS sheds water in the SAME proportion as the
      ! carbon loss (remaining tissue's rwc unchanged), credited to the patch's shed_water_rate;
      ! a NET GAIN touches nothing (the tissue keeps its water mass; its capacity grows).
      !---~---
      block
         type(slow_step_data_t) :: step
         real(wp) :: leaf_water0, wood_water0, leaf_carbon0, fineroot_carbon0, nplant0
         real(wp) :: expect_leaf_shed, expect_wood_shed

         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
         call finalize_init(poly)
         poly%cohort%leaf_water_mass(1) = 4.0_wp
         poly%cohort%wood_water_mass(1) = 9.0_wp
         leaf_carbon0     = poly%cohort%leaf_carbon(1)
         fineroot_carbon0 = poly%cohort%fineroot_carbon(1)
         nplant0          = poly%cohort%nplant(1)
         leaf_water0      = poly%cohort%leaf_water_mass(1)
         wood_water0      = poly%cohort%wood_water_mass(1)

         allocate(step%d_leaf(1), step%d_fineroot(1), step%d_wood(1), step%d_nonstructural(1))
         step%d_wood = 0.0_wp
         step%d_nonstructural = 0.0_wp

         !--- (a) Net LOSS on both leaf and fine root: shed proportional to the carbon fraction lost.
         step%d_leaf(1)     = -0.25_wp * leaf_carbon0        ! sheds 25% of the leaf pool this step
         step%d_fineroot(1) = -0.10_wp * fineroot_carbon0    ! sheds 10% of the fine-root pool
         call shed_turnover_water(poly, cfg, step, no_ledger)
         expect_leaf_shed = leaf_water0 * 0.25_wp
         !---~---
         !   The fine roots hold their carbon's share of the store they share with the sapwood, and
         ! shed a tenth of that: the store's water per unit of its carbon is unchanged once the
         ! fine roots have lost their tenth.
         !---~---
         expect_wood_shed = wood_water0 * 0.10_wp * fineroot_carbon0                                    &
                          / (poly%cohort%sapwood_carbon(1) + fineroot_carbon0)
         call check_close(poly%cohort%leaf_water_mass(1), leaf_water0 - expect_leaf_shed, 1.0e-9_wp, &
                          'shed_turnover_water: leaf water sheds proportional to the carbon loss')
         call check_close(poly%cohort%wood_water_mass(1), wood_water0 - expect_wood_shed, 1.0e-9_wp, &
                          'shed_turnover_water: the fine roots shed their own share of the wood store')
         call check_close(poly%cohort%wood_water_mass(1)                                               &
                          / (poly%cohort%sapwood_carbon(1) + 0.9_wp * fineroot_carbon0),               &
                          wood_water0 / (poly%cohort%sapwood_carbon(1) + fineroot_carbon0), 1.0e-12_wp,   &
                          'shed_turnover_water: the wood store keeps its water per unit carbon')
         call check_close(poly%patch%shed_water_rate(1),                                              &
                          nplant0 * (expect_leaf_shed + expect_wood_shed) / cfg%dt_slow, 1.0e-9_wp,    &
                          'shed_turnover_water: patch rate = nplant-weighted shed amount / dt_slow')
         call check_close(poly%patch%shed_water_enth(1),                                              &
                          nplant0 * (expect_leaf_shed * internal_energy_liquid(poly%cohort%leaf_temp(1))  &
                                   + expect_wood_shed * internal_energy_liquid(poly%cohort%wood_temp(1))) &
                          / cfg%dt_slow, 1.0e-9_wp,                                                   &
                          'shed_turnover_water: the water leaves with its enthalpy at its tissue''s temperature')

         !--- (b) Net GAIN: touches nothing (mass conserved; the larger capacity reads a lower rwc).
         poly%cohort%leaf_water_mass(1) = leaf_water0    ! reset
         poly%cohort%wood_water_mass(1) = wood_water0
         step%d_leaf(1)     = 0.01_wp * leaf_carbon0
         step%d_fineroot(1) = 0.01_wp * fineroot_carbon0
         call shed_turnover_water(poly, cfg, step, no_ledger)
         call check_close(poly%cohort%leaf_water_mass(1), leaf_water0, 1.0e-12_wp,                    &
                          'shed_turnover_water: net GAIN leaves leaf water mass untouched')
         call check_close(poly%cohort%wood_water_mass(1), wood_water0, 1.0e-12_wp,                    &
                          'shed_turnover_water: net GAIN leaves wood water mass untouched')
         call check_close(poly%patch%shed_water_rate(1), 0.0_wp, 1.0e-12_wp,                           &
                          'shed_turnover_water: net GAIN sets patch rate back to 0 (SET, not accumulated)')

         !--- (c) Full snap-to-bare (net loss == the entire pool): sheds to EXACTLY 0.
         poly%cohort%leaf_water_mass(1) = leaf_water0
         step%d_leaf(1)     = -leaf_carbon0
         step%d_fineroot(1) = 0.0_wp
         call shed_turnover_water(poly, cfg, step, no_ledger)
         call check_close(poly%cohort%leaf_water_mass(1), 0.0_wp, 1.0e-9_wp,                          &
                          'shed_turnover_water: full leaf-carbon loss sheds ALL leaf water')
      end block

      !---~---
      !   4b. fill_new_leaf_water reads the cohort's OWN water curves. A bare PFT 2 cohort re-flushes
      ! from a stem a tenth drier than at birth: its new leaves take the water that puts them at
      ! the wood's potential less the crown's gravity head, the wood's potential and the leaves'
      ! water both read on PFT 2's curves (docs/science/soil_plant_water.md §2). The fixture's PFTs
      ! differ in every curve parameter, so PFT 1's curves give the leaves another amount.
      !---~---
      block
         integer(ik), parameter :: PFT = 2_ik
         type(slow_step_data_t) :: step
         type(water_curve_t)    :: leaf_own, wood_own, leaf_pft1, wood_pft1
         real(wp) :: wood_water0, wood_store, head, psi_wood0, psi_leaf, leaf_water, leaf_water_pft1
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, PFT, 0.3_wp, 16.0_wp)
         call finalize_init(poly)
         leaf_own  = cfg%pft%hydro(PFT)%leaf_curve
         wood_own  = cfg%pft%hydro(PFT)%wood_curve
         leaf_pft1 = cfg%pft%hydro(1_ik)%leaf_curve
         wood_pft1 = cfg%pft%hydro(1_ik)%wood_curve
         !--- Bare before the commit, so all its leaf carbon is new and holds no water yet.
         allocate(step%leaf_carbon_before(1))
         step%leaf_carbon_before(1)     = 0.0_wp
         poly%cohort%leaf_water_mass(1) = 0.0_wp
         poly%cohort%wood_water_mass(1) = 0.9_wp * poly%cohort%wood_water_mass(1)
         wood_water0 = poly%cohort%wood_water_mass(1)
         wood_store  = poly%cohort%sapwood_carbon(1) + poly%cohort%fineroot_carbon(1)
         head        = grav_head * poly%cohort%height(1)
         call fill_new_leaf_water(poly, cfg, step)
         leaf_water = poly%cohort%leaf_water_mass(1)
         call check(leaf_water > 0.0_wp, 'fill_new_leaf_water: the PFT 2 cohort''s new leaves take water')
         psi_wood0 = tissue_water_to_psi(wood_water0, wood_own, wood_store)
         psi_leaf  = tissue_water_to_psi(leaf_water, leaf_own, poly%cohort%leaf_carbon(1))
         call check(abs(psi_leaf - (psi_wood0 - head)) <= 1.0e-9_wp,                                         &
                    'fill_new_leaf_water: PFT 2''s new leaves sit at its wood''s potential less the head, on its curves')
         !--- What PFT 1's curves would put in the same leaves: the check above tells the two apart.
         leaf_water_pft1 = tissue_psi_to_water(tissue_water_to_psi(wood_water0, wood_pft1, wood_store) - head,    &
                                               leaf_pft1, poly%cohort%leaf_carbon(1))
         call check(abs(leaf_water_pft1 - leaf_water) > 0.01_wp * leaf_water,                                 &
                    'fill_new_leaf_water: PFT 1''s curves would give the leaves another amount')
      end block

      !---~---
      !   Re-acclimation to a lower SLA asks storage for more leaf carbon: it takes only what
      ! storage holds, so the cohort's carbon is unchanged and its leaf area falls short.
      !---~---
      block
         real(wp) :: carbon0, leaf0, store0
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
         call finalize_init(poly)
         poly%cohort%sla(1) = 2.0_wp * cfg%pft%sla(2_ik)       ! the file's SLA is half the cohort's
         leaf0  = poly%cohort%leaf_carbon(1)
         store0 = 0.25_wp * leaf0                              ! storage covers a quarter of the doubling
         poly%cohort%nonstructural_carbon(1) = store0
         carbon0 = leaf0 + store0
         call reacclimate_plant_traits(poly, cfg, leaf_new, storage_new)
         call commit_trait_carbon(poly, leaf_new, storage_new)
         call check(poly%cohort%sla(1) == cfg%pft%sla(2_ik), 'reacclimation takes the file''s SLA')
         call check(poly%cohort%nonstructural_carbon(1) == 0.0_wp, 'reacclimation empties storage, no further')
         call check_close(poly%cohort%leaf_carbon(1) + poly%cohort%nonstructural_carbon(1), carbon0,      &
                          1.0e-12_wp, 'reacclimation conserves the cohort''s leaf and storage carbon')
         call check(poly%cohort%leaf_carbon(1) < 2.0_wp * leaf0, 'the leaves take only what storage has')
      end block
   end subroutine run_carbon_step

   !---~---
   !   The fine-root carbon [kgC/m2] one compute_carbon_allocation sheds to patch 1's litter,
   ! from empty litter. The step moves no pool, so repeated calls see the same cohort.
   !---~---
   real(wp) function fineroot_litter(poly, cfg) result(shed)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(slow_step_data_t) :: step
      call clear_litter(poly)
      call compute_carbon_allocation(poly, cfg, step)
      shed = poly%patch%litter_in(1)%fineroot_c
   end function fineroot_litter

end module test_slow_step_carbon

!---~---
!   The slow step's paths the ledger cannot see: what the allocator routes out of the slow
! tier, mortality's water, the litter in lockstep with its patch and with soil carbon off,
! and reproduction reaching the recruit pool.
!---~---
module test_slow_step_ledger_paths
   use ed_kinds,                only : wp, ik
   use ed_biogeochem_params,    only : IP_SLOW
   use ed_config,               only : ed_config_t
   use ed_soil_derive,          only : build_profile_hydraulics, soil_layers, build_soil_profile,           &
                                       soil_start_theta, top_soil_psi
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area
   use ed_patch_state_vars, only : cas_set_depth, patch_physical_state, litter_input_t
   use ed_soil_params,        only : soil_profile_t
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_cohort_state_vars,  only : init_cohort, set_cohort_size_from_carbon, hold_cohort_tissue_temps
   use ed_slow_update, only : terminate_cohorts, apply_patch_disturbance, fuse_2_patches, refresh_canopy_depth
   use ed_slow_dynamics, only : vegetation_dynamics, restructure_stand, settle_tissue_water
   use ed_slow_state, only : kill_report_t, kill_litter, shed_turnover_water, slow_step_data_t
   use ed_slow_ledger,          only : slow_store_t, slow_ledger_t, slow_site_store,               &
                                       slow_ledger_open, slow_ledger_declare, slow_ledger_mark,    &
                                       slow_fast_carbon_handover, slow_tissue_energy, SLOW_PHASE_GROW, &
                                       SLOW_PHASE_ALLOCATE, SLOW_PHASE_RECRUIT, SLOW_PHASE_TISSUE_WATER, &
                                       litter_carbon_total
   use ed_plant_hydr_lib,       only : water_curve_t, tissue_water_to_psi
   use ed_constants,            only : grav_head, n_soil_layer_max, umol_2_kgc
   use ed_test_support, only : banner, build_test_config, check, check_close
   use ed_constants,          only : rho_air_ref
   implicit none
   private

   public :: run_ledger_paths
   type(slow_ledger_t) :: no_ledger    ! inactive: a test step books nothing
   real(wp)    :: soil_bottoms(n_soil_layer_max)   ! [m] the column's layer bottoms, as the ledger builds them
   integer(ik) :: n_soil_layer

   real(wp), parameter  :: RHO = rho_air_ref    ! the density the canopy-air store is valued at
   type(ed_config_t)    :: cfg
   type(polygon_t)      :: poly
   type(soil_profile_t) :: soil

contains

   subroutine run_ledger_paths()
      call banner('the slow step: routing, mortality water, litter, recruit pool')

      cfg = build_test_config()
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.3_wp, 18.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.8_wp,  6.0_wp)
      call finalize_init(poly)
      call check(poly%cohort%n == 2_ik .and. poly%patch%n == 1_ik, 'fixture: 2 cohorts in 1 patch')

      associate (sc => cfg%soil)
         call soil_layers(sc, soil_bottoms, n_soil_layer)
         call build_profile_hydraulics(n_soil_layer, sc%retention, soil_bottoms, sc%theta_sat, sc%theta_res, &
                                     sc%ksat, sc%curve_par_a, sc%curve_par_n, sc%psi_fc, soil)
      end associate

      !---~---
      !   Give every reservoir a non-zero value, so a perturbation below is a change to a store
      ! that already exists rather than the creation of one.
      !---~---
      poly%patch%soil_water(1)%theta(1:soil%n_active)      = 0.30_wp
      poly%patch%soil_energy(1)%soil_energy(1:soil%n_active) = 5.0e7_wp
      poly%patch%soil_water(1)%w_surface      = 2.0_wp
      poly%patch%soil_water(1)%w_surface_energy = 8.0e5_wp
      poly%patch%snow(1)%swe(1)           = 5.0_wp
      poly%patch%snow(1)%snow_energy(1)   = -1.0e6_wp
      poly%patch%cas(1)%can_shv      = 0.008_wp
      poly%patch%cas(1)%can_enthalpy = 2.9e5_wp
      poly%patch%cas(1)%can_co2      = 400.0_wp
      poly%patch%cas(1)%can_depth    = 20.0_wp
      poly%patch%soil_carbon(1)%carbon(IP_SLOW) = 3.0_wp
      poly%patch%recruit_pool(:,1)          = 0.0_wp
      poly%cohort%leaf_water_mass(1:2) = 0.02_wp
      poly%cohort%wood_water_mass(1:2) = 0.05_wp
      poly%cohort%leaf_surf_water(1:2) = 0.01_wp
      poly%cohort%wood_surf_water(1:2) = 0.03_wp
      call hold_cohort_tissue_temps(poly%cohort, 1_ik)   ! the tissues' energy with that internal and surface water
      call hold_cohort_tissue_temps(poly%cohort, 2_ik)

      call check_routing()
      call check_mortality_water()
      call check_litter_patch_lockstep()
      call check_litter_without_soil_carbon()
      call check_recruit_pool()
   end subroutine run_ledger_paths

   subroutine check_routing()
      real(wp) :: co2_on, co2_off, lit_lossy, lit_perfect, dummy
      real(wp) :: lab0, lab1, wood0, wood1, n0, dd, dty
      !---~---
      !   Each variant runs on a FRESH stand. vegetation_dynamics COMMITS growth, so calling it
      ! twice on one site compares two different forests and the second comparison passes for
      ! the wrong reason (a mutation that should fail it does not).
      !---~---
      call run_variant(0.3_wp, 1.0_wp, co2_on,  dummy)
      call run_variant(0.0_wp, 1.0_wp, co2_off, dummy)
      call check(co2_on > 0.0_wp, 'routing: growth respiration reaches patch%slow_co2_rate')
      call check(abs(co2_off) < abs(co2_on),                                                       &
                 'routing: zero construction cost => nothing owed (the channel is not a constant)')

      call run_variant(0.3_wp, 0.0_wp, dummy, lit_lossy)     ! every seed dies -> all necromass
      call run_variant(0.3_wp, 1.0_wp, dummy, lit_perfect)   ! every seed establishes -> none
      call check(lit_lossy > lit_perfect,                                                          &
                 'routing: the unestablished seed fraction becomes litter, not nothing')

      !---~---
      !   STARVATION rides the same channel with the OPPOSITE sign. The fast loop has already
      ! exhaled the full maintenance respiration; a stand that could not fund it means the
      ! loop OVER-reported, so the correction owed to the canopy air is NEGATIVE. A sign
      ! error here would quietly turn a starving forest into a CO2 source.
      !---~---
      call run_variant(0.3_wp, 1.0_wp, co2_on, dummy, starve=.true.)
      call check(co2_on < 0.0_wp, 'routing: a starving stand OWES the canopy air a negative flux')

      !---~---
      !   SINK LIMIT: the wood a binding limit leaves unbuilt reaches the soil as exudate, one for
      ! one, and since exudate costs (1+g) like wood the CO2 owed does not change. Mortality is
      ! off in both runs: it rises as growth falls, so the limited stand would also lose more
      ! trees to litter.
      !---~---
      call run_variant(0.3_wp, 1.0_wp, co2_on, dummy, lab_soil=lab0, wood=wood0, nplant0=n0,       &
                       immortal=.true.)
      call run_variant(0.3_wp, 1.0_wp, co2_off, dummy, max_rgr=1.0e-3_wp, lab_soil=lab1, wood=wood1, &
                       immortal=.true.)
      call check(wood1 < wood0, 'sink limit: a binding limit builds less wood')
      call check_close(lab1 - lab0, (wood0 - wood1) * n0, 1.0e-9_wp,                               &
                       'sink limit: the unbuilt wood arrives as below-ground labile litter')
      call check_close(co2_off, co2_on, 1.0e-9_wp,                                                 &
                       'sink limit: the exudate pays the growth respiration the wood would have')
      !--- An ABSOLUTE limit that binds grows the stem by exactly the limit: 0.2 cm/yr for a step.
      call run_variant(0.3_wp, 1.0_wp, dummy, dummy, max_agr=0.2_wp, ddbh=dd, dt_yr=dty, immortal=.true.)
      call check_close(dd, 0.2_wp * dty, 1.0e-9_wp, 'sink limit: a binding absolute limit grows dbh by rate x dt')
   end subroutine check_routing
   !--- One slow step on a FRESH stand, returning the two routed quantities (CO2 rate, litter).
   subroutine run_variant(g_resp, repro_eff, co2_rate, litter_total, starve, max_rgr,              &
                          lab_soil, wood, nplant0, immortal, max_agr, ddbh, dt_yr)
      real(wp), intent(in)  :: g_resp, repro_eff
      real(wp), intent(out) :: co2_rate, litter_total
      logical, optional, intent(in) :: starve
      real(wp), optional, intent(in)  :: max_rgr    ! pft.max_relative_growth_rate (default 0, no limit)
      real(wp), optional, intent(out) :: lab_soil   ! patch 1's below-ground labile litter
      real(wp), optional, intent(out) :: wood       ! the cohort's wood carbon after the step
      real(wp), optional, intent(out) :: nplant0    ! the cohort's density before the step
      logical,  optional, intent(in)  :: immortal   ! .true. => no background or low-growth mortality
      real(wp), optional, intent(in)  :: max_agr    ! pft.max_absolute_growth_rate (default 0, no limit)
      real(wp), optional, intent(out) :: ddbh       ! the cohort's dbh growth over the step
      real(wp), optional, intent(out) :: dt_yr      ! the step's length [yr]
      real(wp) :: dbh_start
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography%demography_on      = .false.        ! isolate allocation from the structural operators
      c%soil_carbon%soil_carbon_on     = .true.
      c%pft%growth_resp_factor(:)      = g_resp
      c%pft%recruit_carbon_efficiency(:) = repro_eff
      !---~---
      !   A SLOW leaf turnover, so turnover does not eat the whole supply before reproduction
      ! is reached. The allocator's priority order is leaf/root growth, then storage, then
      ! reproduction, then wood -- with the default ~1 /yr turnover this cohort's daily leaf
      ! replacement is 14x the photosynthate below and `avail` reaches zero before any seed
      ! is made, so the seed-loss assertion would pass while asserting nothing.
      !---~---
      c%pft%leaf_turnover_rate_toc(:) = 0.01_wp
      if (present(max_rgr)) c%pft%max_relative_growth_rate(:) = max_rgr
      if (present(max_agr)) c%pft%max_absolute_growth_rate(:) = max_agr
      if (present(dt_yr))   dt_yr = c%dt_years
      if (present(immortal)) then
         if (immortal) then
            c%pft%mort_gamma(:) = 0.0_wp
            c%pft%mort_alpha(:) = 0.0_wp
         end if
      end if
      call init_bare_ground(st, c, 1_ik)
      !---~---
      !   40 cm => ~25 m, clear of the 13.9 m reproduction height: a cohort below it allocates
      ! NOTHING to reproduction, and the seed-loss assertion would then be vacuously true.
      !---~---
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp      ! generous: a clear surplus to allocate
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      if (present(starve)) then
         if (starve) then
            !--- No photosynthate, a real maintenance bill, and no reserves to pay it with.
            st%cohort%gpp_accum(1:st%cohort%n)            = 0.0_wp
            st%cohort%leaf_resp_accum(1:st%cohort%n)      = 1.0_wp
            st%cohort%nonstructural_carbon(1:st%cohort%n) = 0.0_wp
         end if
      end if
      if (present(nplant0)) nplant0 = st%cohort%nplant(1)
      dbh_start = st%cohort%dbh(1)
      call vegetation_dynamics(st, c, no_ledger)
      if (present(ddbh)) ddbh = st%cohort%dbh(1) - dbh_start
      if (present(lab_soil)) lab_soil = st%patch%litter_in(1)%labile_soil
      if (present(wood))     wood     = st%cohort%wood_carbon(1)
      co2_rate     = st%patch%slow_co2_rate(1)
      !---~---
      !   The litter accumulator is patch state, so it is read off the site rather than
      ! returned: st%patch%litter_in(1), in lockstep with patch 1.
      !---~---
      litter_total = st%patch%litter_in(1)%labile_grnd + st%patch%litter_in(1)%labile_soil          &
                   + st%patch%litter_in(1)%struct_grnd + st%patch%litter_in(1)%struct_soil
   end subroutine run_variant
   subroutine check_mortality_water()
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      real(wp) :: tis0, tis1, shed0, shed1
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography%demography_on      = .false.        ! isolate the commit from the structural operators
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp     ! grow, so p1 /= p0 and the basis matters
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%leaf_water_mass(1:st%cohort%n) = 0.7_wp
      st%cohort%wood_water_mass(1:st%cohort%n) = 2.3_wp
      tis0  = tissue_water_total(st)
      shed0 = shed_total(st, c)
      call vegetation_dynamics(st, c, no_ledger)
      tis1  = tissue_water_total(st)
      shed1 = shed_total(st, c)
      call check(tis1 < tis0, 'mortality water: the tissue store actually falls')
      call check_close(shed1 - shed0, tis0 - tis1, 1.0e-12_wp,                                     &
                       'every kg that leaves tissue must arrive in the patch shed channel')
   end subroutine check_mortality_water
   pure function tissue_water_total(st) result(w)
      type(polygon_t), intent(in) :: st
      real(wp)    :: w
      integer(ik) :: ip, i, i0, i1
      w = 0.0_wp
      do ip = 1_ik, st%patch%n
         i0 = st%patch%cohort_offset(ip)
         i1 = i0 + st%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            w = w + st%patch%area(ip) * st%cohort%nplant(i)                                        &
                  * (st%cohort%leaf_water_mass(i) + st%cohort%wood_water_mass(i))
         end do
      end do
   end function tissue_water_total
   pure function shed_total(st, c) result(w)
      type(polygon_t),     intent(in) :: st
      type(ed_config_t),   intent(in) :: c
      real(wp)    :: w
      integer(ik) :: ip
      w = 0.0_wp
      do ip = 1_ik, st%patch%n
         w = w + st%patch%area(ip) * st%patch%shed_water_rate(ip) * c%dt_slow
      end do
   end function shed_total
   !---~---
   !   THE LITTER ACCUMULATOR MUST RIDE THE PATCH ARRAY.
   !
   !   restructure_stand (apply_patch_disturbance) CREATES patches after vegetation_dynamics has
   ! filled the step's litter, and the consumer (advance_soil_carbon) loops to the NEW
   ! patch count. An accumulator sized to the patch count when it was filled is read past its
   ! end, feeding undefined memory into the CENTURY litter input (1.8e9 kgC/m2 on a gap patch
   ! born that step, over a 50-year spin-up).
   !
   !   The slow ledger cannot catch that: it DECLARES the same out-of-bounds read as a boundary
   ! term and so balances perfectly against the garbage. Only `-check bounds` sees it.
   !
   !   So the assertion is the INVARIANT -- the accumulator is patch state, at least as long as
   ! the patch array after the operators that change the patch count have run -- plus the
   ! quieter half: a FUSION must carry the litter across area-weighted rather than dropping or
   ! duplicating it (a local array misattributes one patch's litter to another).
   !---~---
   subroutine check_litter_patch_lockstep()
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      real(wp)    :: tot0, tot1
      integer(ik) :: np0

      c = build_test_config()
      c%fast_biophysics_on   = .true.
      c%demography%demography_on        = .true.        ! the structural operators are the point here
      c%demography%patch_disturbance_on = .true.
      c%soil_carbon%soil_carbon_on       = .true.        ! the only consumer of the accumulator
      call init_bare_ground(st, c, 2_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call add_cohort(st, c, 2_ik, 1_ik, 0.3_wp, 35.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      np0 = st%patch%n

      call vegetation_dynamics(st, c, no_ledger)
      call restructure_stand(st, c, .true., .true., no_ledger)   ! monthly + annual: every operator fires

      !---~---
      !   THE invariant. `>=` not `==`: the patch block is capacity-allocated, so the array is
      ! allowed to be longer than the live patch count -- it must never be SHORTER.
      !---~---
      call check(size(st%patch%litter_in) >= st%patch%n,                                           &
                 'the litter accumulator is at least as long as the patch array')
      call check(allocated(st%patch%litter_in), 'the litter accumulator is patch state, not a local')

      !---~---
      !   A fusion must MOVE the litter, not lose it. Fusing two patches area-weighted leaves
      ! the site-wide litter flux unchanged, which is the same rule blend_soil_carbon obeys
      ! -- and has to, since this litter is on its way into those pools.
      !---~---
      if (st%patch%n >= 2_ik) then
         tot0 = litter_site_total(st)
         call fuse_2_patches(st, 1_ik, 2_ik)
         tot1 = litter_site_total(st)
         call check_close(tot1, tot0, 1.0e-12_wp * max(abs(tot0), 1.0_wp),                         &
                          'patch fusion conserves the site-wide litter flux')
      end if
   end subroutine check_litter_patch_lockstep
   !--- Area-weighted site total of today's litter [kgC/m2/day].
   pure function litter_site_total(st) result(t)
      type(polygon_t), intent(in) :: st
      real(wp)    :: t
      integer(ik) :: ip
      t = 0.0_wp
      do ip = 1_ik, st%patch%n
         t = t + st%patch%area(ip) * (st%patch%litter_in(ip)%labile_grnd                           &
                                    + st%patch%litter_in(ip)%labile_soil                           &
                                    + st%patch%litter_in(ip)%struct_grnd                           &
                                    + st%patch%litter_in(ip)%struct_soil)
      end do
   end function litter_site_total
   !---~---
   !   With the soil-carbon pools off, the litter is still recorded and still declared leaving at
   ! the growth phase: the same litter, and the same growth-phase residual, as with them on.
   ! Undeclared, the phase would read it as a carbon leak.
   !---~---
   subroutine check_litter_without_soil_carbon()
      real(wp) :: lit_on, lit_off, res_on, res_off
      call grow_one_step(.true.,  lit_on,  res_on)
      call grow_one_step(.false., lit_off, res_off)
      call check(lit_off > 0.0_wp, 'soil carbon off: the litter is still recorded')
      call check_close(lit_off, lit_on, 1.0e-12_wp * lit_on, 'soil carbon off: the same litter as with it on')
      call check_close(res_off, res_on, 1.0e-12_wp * max(lit_on, 1.0e-12_wp),                        &
                       'soil carbon off: the growth phase closes as it does with it on')
   end subroutine check_litter_without_soil_carbon
   subroutine grow_one_step(soil_carbon_on, litter, grow_residual)
      logical,  intent(in)  :: soil_carbon_on
      real(wp), intent(out) :: litter, grow_residual
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      type(slow_ledger_t) :: lg
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography%demography_on      = .false.        ! one step of growth, turnover and mortality, no restructuring
      c%soil_carbon%soil_carbon_on     = soil_carbon_on
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      lg%active = .true.
      call slow_ledger_open(lg, st, c)
      call vegetation_dynamics(st, c, lg)
      litter        = litter_site_total(st)
      grow_residual = lg%resid_sum(SLOW_PHASE_GROW)%carbon
   end subroutine grow_one_step
   !---~---
   !   The recruit pool is credited EVERY step, and what it holds is the reproduction carbon
   ! the parents were debited for. Two assertions, one per half of the contract: the
   ! CADENCE (credited every step, not sampled on month boundaries from one day's rate
   ! scaled up to stand for the month) and the CARBON LINK (the credit must follow
   ! recruit_carbon_efficiency, which is what makes pool x carbon_min the debited carbon).
   ! A third: under recruit_shade_decay the credit falls as exp(-k * LAI) of the patch.
   !---~---
   subroutine check_recruit_pool()
      real(wp) :: p1, p2, p_lo, p_hi, p_one, p_shade, lai
      call pool_after_steps(3_ik, 1.0e-3_wp, p1)
      call pool_after_steps(6_ik, 1.0e-3_wp, p2)
      call check(p1 > 0.0_wp, 'recruit pool: credited on an ordinary step, not only at month end')
      !---~---
      !   Six steps credit about twice what three do. NOT exactly twice: the stand grows
      ! between steps, so later steps produce a little more reproduction carbon. The point
      ! is that every step contributes -- a monthly sampling would leave p1 and p2 both
      ! identically zero, since neither run crosses a month boundary.
      !---~---
      call check(abs(p2 - 2.0_wp * p1) <= 0.01_wp * p2,                                            &
                 'recruit pool: six steps credit ~twice what three do (it integrates, not samples)')
      call pool_after_steps(3_ik, 1.0e-3_wp, p_lo)
      call pool_after_steps(3_ik, 2.0e-3_wp, p_hi)
      call check_close(p_hi, 2.0_wp * p_lo, 1.0e-9_wp,                                             &
                       'recruit pool: the credit follows recruit_carbon_efficiency')
      !--- One step, so both runs see the same stand and its LAI exactly.
      call pool_after_steps(1_ik, 1.0e-3_wp, p_one, lai=lai)
      call pool_after_steps(1_ik, 1.0e-3_wp, p_shade, shade_decay=0.5_wp)
      call check(lai > 0.1_wp, 'recruit pool: the shade test stand has leaves')
      call check_close(p_shade, p_one * exp(-0.5_wp * lai), 1.0e-9_wp,                             &
                       'recruit pool: the credit falls as exp(-recruit_shade_decay * patch LAI)')
   end subroutine check_recruit_pool
   !---~---
   !   Run `nstep` ordinary (non-month-boundary) slow steps on a fresh stand and return the
   ! accumulated recruit pool. Seed rain is off, so what accumulates is reproduction alone.
   !---~---
   subroutine pool_after_steps(nstep, repro_eff, pool, shade_decay, lai)
      integer(ik), intent(in)  :: nstep
      real(wp),    intent(in)  :: repro_eff
      real(wp),    intent(out) :: pool
      real(wp), optional, intent(in)  :: shade_decay   ! pft.recruit_shade_decay (default 0)
      real(wp), optional, intent(out) :: lai           ! the stand's LAI before the first step
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      integer(ik)         :: k
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography%demography_on      = .false.
      c%pft%recruit_carbon_efficiency(:) = repro_eff
      if (present(shade_decay)) c%pft%recruit_shade_decay(:) = shade_decay
      c%pft%seed_rain_recruits(:)      = 0.0_wp
      c%pft%leaf_turnover_rate_toc(:)  = 0.01_wp    ! [1/yr] keep turnover from eating the supply
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)   ! above the reproduction height
      call finalize_init(st)
      if (present(lai)) lai = sum(st%cohort%nplant(1:st%cohort%n) * st%cohort%leaf_area(1:st%cohort%n))
      do k = 1_ik, nstep
         st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
         st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
         st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
         st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
         call vegetation_dynamics(st, c, no_ledger)
      end do
      pool = st%patch%recruit_pool(1, 1)
   end subroutine pool_after_steps

end module test_slow_step_ledger_paths

!---~---
!   Mortality's carbon by pathway: the cull, the treefall, background mortality and the
! death of a cohort with no leaves and no storage, each on its own row.
!---~---
module test_slow_step_mortality
   use ed_kinds,                 only : wp, ik
   use ed_constants,             only : yr_sec
   use ed_config,                only : ed_config_t
   use ed_cohort_state_vars,     only : set_cohort_size_from_carbon
   use ed_polygon_state_vars,    only : polygon_t
   use ed_diag_vars,             only : alloc_patch_diag, patch_diag_value, PA_MORT_C_BACKGROUND, PA_MORT_C_CULL, &
                                        PA_MORT_C_DISTURB
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_slow_update, only : terminate_cohorts, apply_patch_disturbance
   use ed_slow_dynamics, only : advance_slow_dynamics, advance_boundary_dynamics
   use ed_test_support,        only : banner, build_test_config, check, check_close, test_report
   use ed_slow_ledger,  only : slow_ledger_t
   use ed_slow_state,  only : kill_report_t, decide_disturbance
   implicit none
   private

   public :: run_mortality_pathways
   type(kill_report_t)  :: kills   ! what the cull and the treefall killed
   real(wp)             :: dist_frac          ! the treefall decision's share
   logical, allocatable :: dist_survives(:)   ! and its survivors
   type(slow_ledger_t) :: no_ledger    ! inactive: a test step books nothing
   real(wp)            :: shed_water   ! an operator's shed water, unread here

contains

   subroutine run_mortality_pathways()
      call banner('mortality carbon by pathway')
      call check_cull()
      call check_disturbance()
      call check_background()
      call check_exhausted()
   end subroutine run_mortality_pathways

   !---~---
   !   Arm the per-patch diagnostic block the way a real run does, minus the registry: active,
   ! sized, `n` set to the patch count, and carrying the weight ONE SLOW STEP of fast
   ! sub-steps accumulates. The weight is what patch_diag_value divides out, so a fixture
   ! that leaves it at zero reads every slot as 0 and passes nothing.
   !---~---
   subroutine arm_patch_diag(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      call alloc_patch_diag(poly%patch%diag, poly%patch%n + 2_ik, .true.)
      poly%patch%diag%n = poly%patch%n
      poly%patch%diag%w(1:poly%patch%diag%cap) = cfg%dt_slow
   end subroutine arm_patch_diag
   !---~---
   !   Total carbon one individual of cohort `i` carries. A whole-individual death takes every
   ! pool, which is what separates it from a turnover shed.
   !---~---
   pure real(wp) function plant_carbon(poly, i) result(c)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in) :: i
      c = poly%cohort%leaf_carbon(i) + poly%cohort%fineroot_carbon(i)                              &
        + poly%cohort%wood_carbon(i) + poly%cohort%nonstructural_carbon(i)
   end function plant_carbon
   !---~---
   !   Site total of a patch diagnostic: Sum(area * value), the aggregation the output layer
   ! performs, times the elapsed years to turn the declared [kgC/m2/yr] rate back into the
   ! amount that flowed.
   !---~---
   real(wp) function site_amount(poly, cfg, slot) result(a)
      type(polygon_t),     intent(in) :: poly
      type(ed_config_t),   intent(in) :: cfg
      integer(ik),         intent(in) :: slot
      real(wp)    :: x(16)
      integer(ik) :: n, p
      call patch_diag_value(poly%patch%diag, slot, x, n)
      a = 0.0_wp
      do p = 1_ik, n
         a = a + poly%patch%area(p) * x(p)
      end do
      a = a * (cfg%dt_slow / yr_sec)
   end function site_amount
   !--- CULL: a cohort below the tracked size floor.
   subroutine check_cull()
      type(ed_config_t)   :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: expect

      cfg = build_test_config()
      !---~---
      !   soil_carbon_on OFF deliberately. Both this pathway's litter and the disturbance
      ! pathway's sit behind that switch, and the diagnostic must not: how much biomass died
      ! is demography. With the default (.true.) this part would pass even with the write on
      ! the wrong side of the guard, which is the one placement error it exists to catch.
      !---~---
      cfg%soil_carbon%soil_carbon_on    = .false.
      !---~---
      !   Raise the density floor so a NORMAL cohort is culled: the alternative is a cohort so
      ! small its carbon is at the rounding floor, which would make the assertion vacuous.
      !---~---
      cfg%demography%negligible_nplant = 0.30_wp

      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.20_wp, 20.0_wp)   ! below the floor -> culled
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp, 20.0_wp)   ! above it        -> kept
      call finalize_init(poly)
      call arm_patch_diag(poly, cfg)

      expect = 0.20_wp * plant_carbon(poly, 1_ik)
      call check(expect > 0.0_wp, 'cull fixture: the doomed cohort carries carbon')

      call terminate_cohorts(poly, cfg, shed_water, kills)
      call check(poly%cohort%n == 1_ik, 'cull fixture: exactly one cohort was culled')
      call check_close(site_amount(poly, cfg, PA_MORT_C_CULL), expect, 1.0e-12_wp,                 &
                       'PA_MORT_C_CULL == the culled cohort''s every pool x its density')
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), 0.0_wp, 1.0e-12_wp,           &
                       'a cull is not also counted as background mortality')
      call check_close(site_amount(poly, cfg, PA_MORT_C_DISTURB), 0.0_wp, 1.0e-12_wp,              &
                       'a cull is not also counted as disturbance mortality')
   end subroutine check_cull
   !--- DISTURBANCE: the canopy killed when a gap opens.
   subroutine check_disturbance()
      type(ed_config_t)   :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: frac, expect

      cfg = build_test_config()
      cfg%soil_carbon%soil_carbon_on = .false.          ! the split is demography -- see check_cull
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall  -> dies in the gap
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short -> survives into it
      call finalize_init(poly)
      call check(poly%cohort%height(1) >= cfg%demography%disturbance_survive_height,                          &
                 'disturbance fixture: cohort 1 is canopy')
      call check(poly%cohort%height(2) <  cfg%demography%disturbance_survive_height,                          &
                 'disturbance fixture: cohort 2 is understorey')
      call arm_patch_diag(poly, cfg)

      !---~---
      !   The oracle, on the state BEFORE the call, since the call removes the canopy cohort.
      ! One patch of unit area, so the site amount is just frac x its density x its pools.
      !---~---
      frac   = 1.0_wp - exp(-cfg%demography%patch_disturbance_rate * 1.0_wp)
      expect = poly%patch%area(1) * frac * poly%cohort%nplant(1) * plant_carbon(poly, 1_ik)
      call check(expect > 0.0_wp, 'disturbance fixture: the canopy carries carbon')

      call decide_disturbance(poly, cfg, 1.0_wp, dist_frac, dist_survives)
      call apply_patch_disturbance(poly, cfg, dist_frac, dist_survives, shed_water, kills)

      !---~---
      !   Asserted on the SITE AGGREGATE, not on the donor slot, because that is where the
      ! donor's area shrink lands: the carbon died on the frac of the donor that became the
      ! gap, and only Sum(area * value) after the shrink can show whether it was valued on
      ! the right ground.
      !---~---
      call check_close(site_amount(poly, cfg, PA_MORT_C_DISTURB), expect, 1.0e-12_wp,              &
                       'PA_MORT_C_DISTURB site total == the canopy carbon the gap consumed')
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), 0.0_wp, 1.0e-12_wp,           &
                       'disturbance is not also counted as background mortality')
      call check_close(site_amount(poly, cfg, PA_MORT_C_CULL), 0.0_wp, 1.0e-12_wp,                 &
                       'disturbance is not also counted as a cull')
   end subroutine check_disturbance
   !--- BACKGROUND: the continuous hazard, through the real slow driver.
   subroutine check_background()
      type(ed_config_t)   :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: n0, expect

      cfg = build_test_config()
      cfg%fast_biophysics_on = .false.       ! the hazard is a state question, not a GPP one
      cfg%soil_carbon%soil_carbon_on     = .false.       ! and the split is emitted whatever the soil pools do

      !---~---
      !   ONE cohort, so the per-step sort_cohorts cannot permute the axis between the density
      ! captured before the step and the pools read after it.
      !---~---
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp, 20.0_wp)
      call finalize_init(poly)
      call arm_patch_diag(poly, cfg)
      !---~---
      !   With the fast loop off, the slow step weights the block itself, so the fixture's
      ! stand-in for the fast sub-steps' weight would count the step twice.
      !---~---
      poly%patch%diag%w = 0.0_wp

      n0 = poly%cohort%nplant(1)
      call advance_slow_dynamics(poly, cfg, no_ledger)                      ! no month/year boundary
      call check(poly%cohort%n == 1_ik, 'background fixture: the cohort survived the step')
      call check(n0 - poly%cohort%nplant(1) > 0.0_wp, 'the step actually killed some density')

      !---~---
      !   Valued on the density the applier REMOVED and the pools it LEFT: growth happens at the
      ! old density and mortality at the new pools, which is the decomposition
      ! accumulate_mortality_litter's header sets out. Reading the pools after the step is
      ! therefore not a convenience -- it is the definition.
      !---~---
      expect = (n0 - poly%cohort%nplant(1)) * plant_carbon(poly, 1_ik)
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), expect, 1.0e-12_wp,           &
                       'PA_MORT_C_BACKGROUND == the removed density x the pools it left')
      call check_close(site_amount(poly, cfg, PA_MORT_C_CULL), 0.0_wp, 1.0e-12_wp,                 &
                       'background mortality is not also counted as a cull')
      call check_close(site_amount(poly, cfg, PA_MORT_C_DISTURB), 0.0_wp, 1.0e-12_wp,              &
                       'background mortality is not also counted as disturbance')
   end subroutine check_background
   !---~---
   !   A cohort with no leaves and no storage dies; one living on its storage does not.
   ! Both are leafless and dormant (no flush, so nothing grows) with a day's stem respiration
   ! and no photosynthesis, as in a dry season; only their storage differs. The first cannot pay
   ! its upkeep, so its density goes to the floor on its first step, through the background
   ! pathway, and the next month boundary culls it. The second pays its upkeep from storage
   ! through a 150-day dry season and loses density only at the growth-driven hazard.
   !---~---
   subroutine check_exhausted()
      type(ed_config_t)   :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: n0, expect
      integer(ik)         :: day

      cfg = build_test_config()
      cfg%fast_biophysics_on = .true.        ! the allocator reads the accumulators set below
      cfg%soil_carbon%soil_carbon_on     = .false.

      !--- No storage: dead on the first step, at the floor, its carbon down the background path.
      call leafless_stand(poly, cfg, 0.0_wp)
      call arm_patch_diag(poly, cfg)        ! fast loop on: the fixture's weight stands in for it
      n0     = poly%cohort%nplant(1)
      expect = (n0 - cfg%demography%negligible_nplant) * plant_carbon(poly, 1_ik)
      call advance_slow_dynamics(poly, cfg, no_ledger)
      call check(poly%cohort%nplant(1) == cfg%demography%negligible_nplant,                                   &
                 'no leaves, no storage: the density reaches the floor on the first step')
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), expect, 1.0e-9_wp,            &
                       'no leaves, no storage: its carbon dies down the background pathway')
      do day = 2_ik, 30_ik
         call advance_slow_dynamics(poly, cfg, no_ledger)
      end do
      call advance_boundary_dynamics(poly, cfg, .true., .false., no_ledger)
      call check(poly%cohort%n == 0_ik, 'no leaves, no storage: culled at the month boundary')

      !--- Storage for the season: the cohort lives on it, leafless, for 150 days.
      call leafless_stand(poly, cfg, 0.5_wp)
      n0 = poly%cohort%nplant(1)
      do day = 1_ik, 150_ik
         call advance_slow_dynamics(poly, cfg, no_ledger)
         if (mod(day, 30_ik) == 0_ik) call advance_boundary_dynamics(poly, cfg, .true., .false., no_ledger)
      end do
      call check(poly%cohort%n == 1_ik, 'leafless with storage: survives a 150-day dry season')
      call check(poly%cohort%nonstructural_carbon(1) > 0.0_wp .and. poly%cohort%leaf_carbon(1) == 0.0_wp, &
                 'leafless with storage: it paid its upkeep from storage and stayed leafless')
      call check(poly%cohort%nplant(1) > 0.9_wp * n0,                                              &
                 'leafless with storage: only the growth-driven hazard thinned it')
   end subroutine check_exhausted
   !---~---
   !   One leafless, dormant cohort holding `storage` [kgC/plant], with a day's stem respiration
   ! and no photosynthesis on every step (there is no fast loop to overwrite the accumulators).
   !---~---
   subroutine leafless_stand(poly, cfg, storage)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      real(wp),            intent(in)    :: storage
      real(wp), parameter :: STEM_RESP_DAY = 1.0e-3_wp   ! [kgC/plant] one day's stem respiration
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.05_wp, 20.0_wp)
      call finalize_init(poly)
      poly%cohort%leaf_carbon(1)          = 0.0_wp
      poly%cohort%nonstructural_carbon(1) = storage
      call set_cohort_size_from_carbon(poly%cohort, 1_ik)
      poly%cohort%leaf_flush_tendency(1) = 0.0_wp      ! dormant: no flush, no turnover
      poly%cohort%gpp_accum(1)           = 0.0_wp
      poly%cohort%leaf_resp_accum(1)     = 0.0_wp
      poly%cohort%stem_resp_accum(1)     = STEM_RESP_DAY
      poly%cohort%root_resp_accum(1)     = 0.0_wp
   end subroutine leafless_stand

end module test_slow_step_mortality

!---~---
!   The soil-carbon seam through the real driver: the daily step's Rh is the fast loop's,
! the pool moves and stays non-negative, and soil carbon off leaves it untouched; without
! the fast loop the pools decay at spinup_xi and there is no seam.
!---~---
module test_slow_step_seam
   use ed_kinds,                 only : wp, ik
   use ed_biogeochem_params,     only : IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, &
                                        IP_SLOW, N_SOIL_POOL
   use ed_respiration,           only : decomposition_transfer_matrix
   use ed_config,                only : ed_config_t
   use ed_soil_derive,           only : build_profile_hydraulics, build_profile_thermal
   use ed_soil_params,           only : soil_profile_t
   use ed_polygon_state_vars,    only : polygon_t
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_test_soil_grid, only : geometric_layer_bottoms
   use ed_constants,     only : n_soil_layer_max, day_sec
   use ed_soil_hydr_lib,   only : SOIL_RETENTION_VG
   use ed_fast_dynamics,         only : fast_dynamics
   use ed_init,                  only : init_fast_state
   use ed_slow_dynamics, only : advance_slow_dynamics
   use ed_patch_state_vars, only : litter_input_t
   use ed_slow_ledger, only : slow_ledger_t
   use ed_test_support, only : banner, build_test_config, check, check_close, check_true, test_report
   implicit none
   private

   public :: run_soil_seam
   type(slow_ledger_t) :: no_ledger    ! inactive: a test step books nothing
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks are computed on

   integer(ik)            :: n_budget_fail
   integer(ik), parameter :: NSL = 10_ik
   type(ed_config_t)    :: cfg
   type(polygon_t)      :: poly
   type(soil_profile_t) :: test_profile   ! the site's soil column (use_test_column)
   real(wp) :: total0, total1
   type(slow_ledger_t) :: seam   ! an inactive ledger, which keeps the run's seam worsts
   real(wp) :: we, ww

contains

   subroutine run_soil_seam()
      call banner('slow soil-carbon biogeochemistry (double-count gate)')

      cfg = build_test_config()
      cfg%fast_biophysics_on = .true.
      cfg%dt_fast            = 150.0_wp
      cfg%n_fast_per_slow    = 48_ik   ! 48 x 150 s = a 2 h fast window
      cfg%soil_carbon%soil_carbon_on     = .true.

      call geometric_layer_bottoms(2.0_wp, NSL, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(NSL, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,          &
                             2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, test_profile)
      call build_profile_thermal(NSL, 3.0_wp, 0.15_wp, 2.0e6_wp, test_profile)
      cfg%init%soil_temp = 295.0_wp

      call init_bare_ground(poly, cfg, 1_ik)

      call use_test_column(poly)
      call add_cohort(poly, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)     ! one climax cohort (turnover -> litter)
      call finalize_init(poly)
      call init_fast_state(poly, cfg)

      !--- Seed a nonzero starting pool so decomposition (and hence Rh) is actually active.
      poly%patch%soil_carbon(1)%carbon(IP_FAST_GRND)   = 2.0_wp
      poly%patch%soil_carbon(1)%carbon(IP_FAST_SOIL)   = 1.0_wp
      poly%patch%soil_carbon(1)%carbon(IP_STRUCT_GRND) = 20.0_wp
      poly%patch%soil_carbon(1)%carbon(IP_STRUCT_SOIL) = 10.0_wp
      poly%patch%soil_carbon(1)%carbon(IP_SLOW)        = 500.0_wp

      total0 = sum(poly%patch%soil_carbon(1)%carbon)

      !---~---
      !   1. Run one real slow step: fast loop (accumulates decomp_xi against the FROZEN
      ! pool seeded above) THEN advance_slow_dynamics (vegetation turnover -> litter, then the
      ! daily soil_carbon_step consuming that litter + xi_int).
      !---~---
      call fast_dynamics(poly, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
      call check(n_budget_fail == 0_ik, 'fast loop whole-column budgets closed (n_fail == 0)')

      block
         real(wp) :: rh_fast_accum_ref
         rh_fast_accum_ref = poly%patch%decomp_xi(1)%rh_fast
         call check(rh_fast_accum_ref > 0.0_wp, 'fast loop accumulated positive matrix Rh over the day')

         call advance_slow_dynamics(poly, cfg, seam)

         call check_true('rh_seam_gap ~ 0 (the double-count gate)',                                   &
                         seam%worst_rh_gap < 1.0e-9_wp * max(rh_fast_accum_ref, 1.0_wp), seam%worst_rh_gap)
         !---~---
         !   The other two seam numbers the same call accumulates. Lambda is the freeze number:
         ! the fraction of a pool one slow step withdraws. It must be nonzero (the step DID
         ! decompose) and far below 1 (or freezing the pool across the step is unsound).
         !---~---
         call check_true('lignin passive-tracer residual ~ 0', seam%worst_lignin < 1.0e-12_wp,        &
                         seam%worst_lignin)
         call check_true('lambda is nonzero -- the step actually decomposed', seam%worst_lambda > 0.0_wp, &
                         seam%worst_lambda)
         call check_true('lambda << 1 -- the frozen-pool approximation is sound at a daily step',      &
                         seam%worst_lambda < 0.1_wp, seam%worst_lambda)
      end block

      total1 = sum(poly%patch%soil_carbon(1)%carbon)

      !---~---
      !   Direction is not asserted (this cohort's turnover litter can outweigh or be outweighed by
      ! Rh depending on the seeded pool -- test_soil_carbon's INDEPENDENT mass check (Rh
      ! recomputed from er*xi_int*K*X0, not from the closure definition) already guarantees
      ! dC_pool == litter_in - rh_today to machine precision); this just checks the driver seam
      ! actually moved the pool and stayed physical.
      !---~---
      call check(abs(total1 - total0) > 1.0e-6_wp, 'pool total changed (the daily step actually ran)')
      call check(total1 > 0.0_wp, 'pool total stayed physical (non-negative)')
      write(*,'(a,es12.5,a,es12.5,a)') '   (pool total: ', total0, ' -> ', total1, ' kgC/m2)'
      write(*,'(a,es12.5,a,es12.5,a)') '   (rh_seam_gap: ', seam%worst_rh_gap,                       &
                                       '   lambda: ', seam%worst_lambda, ')'

      !--- 2. OFF-PATH: soil_carbon_on = .false. must leave the pool exactly untouched.
      block
         type(polygon_t) :: poly2
         real(wp)     :: totalA, totalB
         cfg%soil_carbon%soil_carbon_on = .false.
         call init_bare_ground(poly2, cfg, 1_ik)
         call use_test_column(poly2)
         call add_cohort(poly2, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)
         call finalize_init(poly2)
         call init_fast_state(poly2, cfg)
         poly2%patch%soil_carbon(1)%carbon(IP_SLOW) = 500.0_wp
         totalA = poly2%patch%soil_carbon(1)%carbon(IP_SLOW)
         call fast_dynamics(poly2, cfg)
         call advance_slow_dynamics(poly2, cfg, no_ledger)
         totalB = sum(poly2%patch%soil_carbon(1)%carbon)
         call check_close(totalB, totalA, 1.0e-12_wp, 'soil_carbon_on=.false.: pool exactly untouched')
      end block

      !---~---
      !   3. WITHOUT THE FAST LOOP the pools decay at spinup_xi over the step and there is no seam
      ! (soil_carbon.md §7). Lambda, the step's largest xi_int * k, is then spinup_xi times the
      ! step's days times the largest decay rate; with xi_int left at zero it was 0.
      !---~---
      block
         type(polygon_t)     :: poly3
         type(slow_ledger_t) :: seam3
         real(wp) :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
         real(wp), parameter :: XI = 0.4_wp
         cfg%soil_carbon%soil_carbon_on = .true.
         cfg%fast_biophysics_on         = .false.
         cfg%soil_carbon%spinup_xi      = XI
         call init_bare_ground(poly3, cfg, 1_ik)
         call use_test_column(poly3)
         call add_cohort(poly3, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)
         call finalize_init(poly3)
         call init_fast_state(poly3, cfg)
         poly3%patch%soil_carbon(1)%carbon(IP_FAST_GRND) = 2.0_wp
         poly3%patch%soil_carbon(1)%carbon(IP_SLOW)      = 500.0_wp
         call decomposition_transfer_matrix(poly3%patch%soil_carbon(1), cfg%soil_carbon, a_mat, k_diag, er)
         call advance_slow_dynamics(poly3, cfg, seam3)
         call check_close(seam3%worst_lambda, XI * cfg%dt_slow / day_sec * maxval(k_diag), 1.0e-12_wp,  &
                          'no fast loop: the pools decay at spinup_xi over the step')
         call check_true('no fast loop: no seam gap', seam3%worst_rh_gap == 0.0_wp, seam3%worst_rh_gap)
      end block
   end subroutine run_soil_seam

   !---~---
   !   The test's hand-built soil column is the site's, which the fast loop runs on, starting
   ! at 0.30 in every layer.
   !---~---
   subroutine use_test_column(p)
      type(polygon_t), intent(inout) :: p
      p%site%surface(1)%soil_profile = test_profile
      p%site%surface(1)%theta_start  = 0.30_wp
   end subroutine use_test_column

end module test_slow_step_seam

!---~---
!   The slow rows' units: each reported rate times the elapsed time is the amount the step
! moved, with and without a fast loop.
!---~---
module test_slow_step_diag
   use ed_kinds,                 only : wp, ik
   use ed_constants,             only : yr_sec, n_soil_layer_max
   use ed_config,                only : ed_config_t
   use ed_soil_derive,           only : build_profile_hydraulics, build_profile_thermal
   use ed_soil_params,           only : soil_profile_t
   use ed_polygon_state_vars,    only : polygon_t, reset_step_diagnostics
   use ed_diag_vars,             only : alloc_patch_diag, patch_diag_value, patch_diag_slow_row, PA_LITTER_LEAF, &
                                        PA_LITTER_FINEROOT, PA_LITTER_STRUCT, PA_RECRUIT_NPLANT, PA_DISTURB_AREA, PA_LE, PA_GPP
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_test_soil_grid, only : geometric_layer_bottoms
   use ed_soil_hydr_lib,         only : SOIL_RETENTION_VG
   use ed_fast_dynamics,         only : fast_dynamics
   use ed_init,                  only : init_fast_state
   use ed_slow_dynamics, only : advance_slow_dynamics, advance_boundary_dynamics, vegetation_dynamics
   use ed_test_support,          only : banner, build_test_config, check, check_close, check_true, &
                                        test_report
   use ed_slow_ledger,  only : slow_ledger_t
   implicit none
   private

   public :: run_slow_rows
   type(slow_ledger_t) :: no_ledger    ! inactive: a test step books nothing
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks are computed on

   integer(ik), parameter :: NSL = 10_ik
   type(ed_config_t)    :: cfg
   type(polygon_t)      :: poly
   type(soil_profile_t) :: test_profile   ! the site's soil column (use_test_column)
   real(wp)             :: x(8), pool0, pool1, elapsed_yr, frac_expected
   integer(ik)          :: np

contains

   subroutine run_slow_rows()
      call banner('slow per-patch diagnostic units')

      cfg = build_test_config()
      cfg%fast_biophysics_on = .true.
      !---~---
      !   ONE full slow step of fast sub-steps: 96 x 900 s = 86400 s = cfg%dt_slow. This is the
      ! property check 1 asserts; it is set here rather than inherited so the mirror is explicit.
      !---~---
      cfg%dt_fast            = 900.0_wp
      cfg%n_fast_per_slow    = 96_ik
      cfg%soil_carbon%soil_carbon_on     = .true.

      call geometric_layer_bottoms(2.0_wp, NSL, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(NSL, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,        &
                             2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, test_profile)
      call build_profile_thermal(NSL, 3.0_wp, 0.15_wp, 2.0e6_wp, test_profile)
      cfg%init%soil_temp = 295.0_wp

      call init_bare_ground(poly, cfg, 1_ik)

      call use_test_column(poly)
      call add_cohort(poly, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)     ! one climax cohort (turnover -> litter)
      call finalize_init(poly)
      call init_fast_state(poly, cfg)

      !---~---
      !   Turn the patch diagnostic block on by hand. In production ed_output_registry does this
      ! from the enabled-variable set; a test has no registry, and the block is a no-op at every
      ! entry point until it is active.
      !---~---
      call alloc_patch_diag(poly%patch%diag, max(poly%patch%cap, 1_ik), .true.)

      !--- 1. THE FIXTURE MIRRORS THE DRIVER: one fast window weighs exactly one slow step.
      call fast_dynamics(poly, cfg)
      call check_close(poly%patch%diag%w(1), cfg%dt_slow, 1.0e-9_wp,                                  &
                       'fast window accumulates exactly dt_slow of weight (the fixture mirrors the driver)')
      elapsed_yr = poly%patch%diag%w(1) / yr_sec

      !---~---
      !   2. LITTER: reported rate x elapsed == the amount biogeochemistry consumed.
      ! The cohorts' day writes the rows and fills the accumulator; the soil step would take
      ! the litter up and clear it, so the check reads it before.
      !---~---
      pool0 = sum(poly%patch%recruit_pool(:, 1))
      call vegetation_dynamics(poly, cfg, no_ledger)
      pool1 = sum(poly%patch%recruit_pool(:, 1))

      call check(poly%patch%litter_in(1)%leaf_c > 0.0_wp, 'the step actually shed leaf litter')

      call patch_diag_value(poly%patch%diag, PA_LITTER_LEAF, x, np)
      call check(np >= 1_ik, 'the patch diagnostic block reports at least one slot')
      call check_close(x(1) * elapsed_yr, poly%patch%litter_in(1)%leaf_c, 1.0e-10_wp,                 &
                       'PA_LITTER_LEAF rate x elapsed == the leaf litter amount')
      write(*,'(a,es12.5,a)') '   (leaf litterfall: ', x(1), ' kgC/m2/yr)'

      call patch_diag_value(poly%patch%diag, PA_LITTER_FINEROOT, x, np)
      call check_close(x(1) * elapsed_yr, poly%patch%litter_in(1)%fineroot_c, 1.0e-10_wp,             &
                       'PA_LITTER_FINEROOT rate x elapsed == the fine-root litter amount')

      call patch_diag_value(poly%patch%diag, PA_LITTER_STRUCT, x, np)
      call check_close(x(1) * elapsed_yr,                                                             &
                       poly%patch%litter_in(1)%struct_grnd + poly%patch%litter_in(1)%struct_soil,     &
                       1.0e-10_wp, 'PA_LITTER_STRUCT rate x elapsed == the structural litter amount')

      !---~---
      !   3. RECRUITMENT: the same identity against the pool the SAME step credits. The pool is
      ! consumed only on a month boundary, and this step is not one, so its change IS the amount.
      !---~---
      call patch_diag_value(poly%patch%diag, PA_RECRUIT_NPLANT, x, np)
      call check(pool1 - pool0 > 0.0_wp, 'the step actually credited the recruit pool')
      call check_close(x(1) * elapsed_yr, pool1 - pool0, 1.0e-10_wp,                                  &
                       'PA_RECRUIT_NPLANT rate x elapsed == the recruit-pool credit')

      !---~---
      !   4. DISTURBANCE: fires on a year boundary, in a different file with its own expression, so
      ! it gets its own check. The amount is the area fraction its hazard removes over the one-
      ! year PATCH_DYNAMICS_INTERVAL it is called with. As in the driver, the boundary comes after
      ! the output has read the step and zeroed the block, so its event lands in the NEXT window,
      ! which the next fast loop then weighs.
      !---~---
      block
         call reset_step_diagnostics(poly)                          ! the output has read step 1
         call advance_boundary_dynamics(poly, cfg, .true., .true., no_ledger)  ! new month AND new year
         call fast_dynamics(poly, cfg)                              ! the next step's window
         elapsed_yr = poly%patch%diag%w(1) / yr_sec
         call check_close(poly%patch%diag%w(1), cfg%dt_slow, 1.0e-9_wp,                               &
                          'the window after the boundary weighs exactly one slow step')
         frac_expected = 1.0_wp - exp(-cfg%demography%patch_disturbance_rate * 1.0_wp)   ! PATCH_DYNAMICS_INTERVAL
         call check(frac_expected > 0.0_wp, 'the test config actually disturbs')
         call patch_diag_value(poly%patch%diag, PA_DISTURB_AREA, x, np)
         call check_close(x(1) * elapsed_yr, frac_expected, 1.0e-9_wp,                                &
                          'PA_DISTURB_AREA rate x elapsed == the area fraction disturbed')
         write(*,'(a,es12.5,a,es12.5,a)') '   (disturbance: reported ', x(1),                         &
                                          ' 1/yr  vs configured hazard ', cfg%demography%patch_disturbance_rate, ' 1/yr)'
      end block

      !---~---
      !   5. A SLOW-ONLY RUN. No fast loop weighs the block or sets its patch count, so the
      ! slow step does both: the slow rows read their rates exactly as with a fast loop, and the
      ! fast rows are not slow rows (the output layer reports those as missing).
      !---~---
      block
         type(ed_config_t) :: cfg_slow
         cfg_slow = cfg
         cfg_slow%fast_biophysics_on = .false.
         call reset_step_diagnostics(poly)
         poly%patch%diag%n = 0_ik                                   ! as if no fast loop had ever run
         pool0 = sum(poly%patch%recruit_pool(:, 1))
         call advance_slow_dynamics(poly, cfg_slow, no_ledger)
         pool1 = sum(poly%patch%recruit_pool(:, 1))
         call check_close(poly%patch%diag%w(1), cfg%dt_slow, 1.0e-9_wp,                               &
                          'slow-only: the slow step weighs the block by dt_slow')
         call check(poly%patch%diag%n == poly%patch%n, 'slow-only: the slow step sets the patch count')
         elapsed_yr = poly%patch%diag%w(1) / yr_sec
         call patch_diag_value(poly%patch%diag, PA_RECRUIT_NPLANT, x, np)
         call check(np >= 1_ik .and. pool1 - pool0 > 0.0_wp, 'slow-only: the recruit row has a slot')
         call check_close(x(1) * elapsed_yr, pool1 - pool0, 1.0e-10_wp,                               &
                          'slow-only: PA_RECRUIT_NPLANT rate x elapsed == the recruit-pool credit')
         call check_true('slow-only: litter, recruitment and disturbance are slow rows',              &
                         patch_diag_slow_row(PA_LITTER_LEAF) .and. patch_diag_slow_row(PA_RECRUIT_NPLANT) &
                         .and. patch_diag_slow_row(PA_DISTURB_AREA))
         call check_true('slow-only: latent heat and GPP are fast rows',                              &
                         .not. (patch_diag_slow_row(PA_LE) .or. patch_diag_slow_row(PA_GPP)))
      end block
   end subroutine run_slow_rows

   !---~---
   !   The test's hand-built soil column is the site's, which the fast loop runs on, starting
   ! at 0.30 in every layer.
   !---~---
   subroutine use_test_column(p)
      type(polygon_t), intent(inout) :: p
      p%site%surface(1)%soil_profile = test_profile
      p%site%surface(1)%theta_start  = 0.30_wp
   end subroutine use_test_column

end module test_slow_step_diag

!---~---
!   test_slow_step -- the slow step through ed_slow_state and ed_slow_dynamics
! (docs/science/order_of_processes.md §4): leaf phenology's wiring, thermal acclimation's
! wiring, the carbon step, the paths the ledger cannot see, mortality by pathway, the
! soil-carbon seam and the slow rows' units. One module per section.
!---~---
program test_slow_step
   use test_slow_step_phenology, only : run_phenology_wiring
   use test_slow_step_acclimation, only : run_thermal_acclimation
   use test_slow_step_carbon, only : run_carbon_step
   use test_slow_step_ledger_paths, only : run_ledger_paths
   use test_slow_step_mortality, only : run_mortality_pathways
   use test_slow_step_seam, only : run_soil_seam
   use test_slow_step_diag, only : run_slow_rows
   use ed_test_support, only : test_report
   implicit none

   call run_phenology_wiring()
   call run_thermal_acclimation()
   call run_carbon_step()
   call run_ledger_paths()
   call run_mortality_pathways()
   call run_soil_seam()
   call run_slow_rows()

   call test_report('test_slow_step')
end program test_slow_step
