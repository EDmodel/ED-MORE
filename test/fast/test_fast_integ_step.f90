! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_fast_integ_step -- one fast step of a column, end to end (ed_fast_integ_stepper:    !
! fast_step; docs/science/numerical_scheme.md §2-§5), on the shared fast fixture: the      !
! frozen record, the water network's step, the march under ESDIRK2 or RK45, the update and !
! the seven ledgers. A scenario that both schemes take runs on both.                       !
!   1. THE COLUMN (static_state_fixture, build_frozen): the gathered view brings the       !
!      sapwood area in m2, the PFT's aboveground fraction and its leaf and wood            !
!      emissivities to the step.                                                           !
!   2. A DIURNAL DAY with a morning rain (ESDIRK2): every ledger closes; the canopy air,   !
!      the leaf, the soil and the leaf potential follow the sun; midday is unstable and    !
!      photosynthesis draws the canopy CO2 down; the soil's swing damps with depth.        !
!   3. THE REFERENCE LEVEL: a dark column at the air's temperature reports almost no       !
!      sensible heat, on both schemes.                                                     !
!   4. ESDIRK2 AND RK45 AGREE on a forced day at the shipped 900 s: both close every       !
!      ledger, and the day's sensible and latent heat, net radiation, GPP, NEE, ground     !
!      heat and transpiration agree within what their error control allows.                 !
!   5. THE LEDGERS ON EVERY PATH, on both schemes: a dry day, steady rain, shed water,     !
!      canopy films; a sealed column saturating and ponding; a free-draining one          !
!      saturating under a storm (ESDIRK2); an aquifer wetting a dry column from below.     !
!   6. SNOW: a melting pack, shed water under a pack and snowfall each close the column's  !
!      water and energy on both schemes; snowfall arrives in the column; snowfall too light !
!      to start a pack reaches the ground as ice.                                           !
!   7. PLANT WATER AND HEAT: the leaf potential converges in dt_fast; the wood's heat      !
!      store lags the canopy air on both schemes; refilling wood takes no heat from water  !
!      at its own temperature; the wood's heat capacity follows its water; the canopy      !
!      films stay out of the tissues' heat capacity for the step; a new cohort's predawn   !
!      potential follows a Campbell soil's curve.                                          !
!   8. RK45 AT ITS EDGES: a drydown to theta_res, a just-recruited cohort, a dense canopy   !
!      in a cold snap, and the dry soil's uptake seam.                                      !
!   9. THE REPORTED NET RADIATION (fast_step's out%rnet) counts the wood: with the leaves  !
!      and the ground dark, shortwave on the wood alone raises it by what the wood         !
!      absorbs, less only the warmer wood's extra emission, on both schemes.               !
! The building blocks are test_fast_integ_kernels'; the schemes' orders and the adaptive  !
! march are test_fast_integ_schemes'.                                                       !
!==========================================================================================!
program test_fast_integ_step
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : latent_heat_fusion, rho_h2o, grav, pio4, cp_liq, grav_head, mmh2o,          &
                                       tissue_c_wood_dry, tissue_c_leaf
   use ed_fast_params,          only : INTEG_ESDIRK2, INTEG_RK45, RK45_RTOL
   use ed_time,                 only : solar_cosz
   use ed_therm_lib,            only : temp_to_internal_energy, internal_energy_liquid, sat_specific_humidity,     &
                                       tissue_temp_to_energy, tissue_dry_heat_capacity
   use ed_soil_hydr_lib,        only : SOIL_RETENTION_CAMPBELL, SOIL_SPECIFIC_STORAGE, soil_theta_to_head
   use ed_soil_params,          only : soil_params_t, SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use ed_site_state_vars,      only : make_site_surface
   use ed_diag_vars,            only : CO_PSI_PREDAWN, CO_TRANSP
   use ed_plant_hydr_lib,       only : tissue_water_to_psi
   use ed_fast_state,           only : derive_cohort_terms
   use ed_fast_integ_types,     only : frozen_t, march_state_t, ensure_output_capacity, WORK_RK45_RESCUE,          &
                                       WORK_INTEG_STEPS
   use ed_test_fast_fixture,    only : fast_column_t, make_fast_column, gather_column, reset_column, step_column,  &
                                       set_diurnal_forcing, set_noon_forcing, set_coldsnap_forcing, set_air,        &
                                       set_ground_tiles, frozen_record, all_ledgers_closed, NSL, DT_FAST,           &
                                       STEPS_PER_HOUR, STEPS_PER_DAY, T0, THETA0, LAT, SIM_DATE, SHV_AIR
   use ed_test_assert,          only : check, check_true, test_report
   implicit none

   integer(ik), parameter :: SCHEMES(2) = [INTEG_ESDIRK2, INTEG_RK45]

   !----- One day of the column (march_day): what drives it beyond the fixture's diurnal forcing, !
   !      and what it records. ----------------------------------------------------------------!
   type :: day_t
      real(wp)    :: theta_seed   = THETA0   !< [m3/m3] the soil's starting water
      integer(ik) :: hours        = 24_ik    !< the day's length
      real(wp)    :: rain_pulse   = 0.0_wp   !< [kg/m2/s] rain over 03:00-07:00
      real(wp)    :: rain_all_day = 0.0_wp   !< [kg/m2/s] rain all day
      real(wp)    :: snowfall     = 0.0_wp   !< [kg/m2/s] snowfall all day, under sub-freezing air
      logical     :: cold_air     = .false.  !< the snowfall's sub-freezing air, without the snowfall
      real(wp)    :: snow_seed    = 0.0_wp   !< [kg/m2] a pack at the start
      real(wp)    :: shed         = 0.0_wp   !< [kg/m2/s] the day's turnover shed water
      real(wp)    :: shv_air      = SHV_AIR  !< [kg/kg] the air's humidity
      logical     :: sealed       = .false.  !< a bedrock bottom: no drainage
      logical     :: aquifer      = .false.  !< an aquifer bottom
      logical     :: dark_step    = .false.  !< one dark step with the air at the column's temperature
      !----- Recorded. ----------------------------------------------------------------------!
      real(wp) :: ct_noon = 0.0_wp, ct_night = 0.0_wp, co2_noon = 0.0_wp, co2_night = 0.0_wp
      real(wp) :: tleaf_noon = 0.0_wp, tleaf_night = 0.0_wp, ustar_noon = 0.0_wp, dtheta_noon = 0.0_wp
      real(wp) :: gpp_noon = 0.0_wp, nee_noon = 0.0_wp, psileaf_noon = 0.0_wp, psileaf_night = 0.0_wp
      real(wp) :: ss_min = 1.0e9_wp, ss_max = -1.0e9_wp, sd_min = 1.0e9_wp, sd_max = -1.0e9_wp
      real(wp) :: th_min = 1.0e9_wp, th_max = -1.0e9_wp, theta_min = 1.0e9_wp
      real(wp) :: film_peak = 0.0_wp, pond_peak = 0.0_wp, theta_peak = 0.0_wp, ood_peak = 0.0_wp
      real(wp) :: wood_lag = 0.0_wp, worst_super = -1.0_wp
      real(wp) :: h_first = 0.0_wp, gah_first = 0.0_wp
      logical  :: physical = .true., snow_physical = .true.
   end type day_t

   type(fast_column_t) :: col

   call make_fast_column(col)

   call test_the_column()
   call test_diurnal_day()
   call test_reference_level()
   call test_schemes_agree()
   call test_dry_day()
   call test_wet_day()
   call test_shed_water()
   call test_canopy_films()
   call test_sealed_saturated()
   call test_free_drain_saturated()
   call test_aquifer()
   call test_snow_melt()
   call test_shed_under_pack()
   call test_snowfall_conserved()
   call test_light_snowfall()
   call test_psi_dt_convergence()
   call test_prognostic_wood()
   call test_refilling_wood()
   call test_wood_capacity()
   call test_film_outside_capacity()
   call test_campbell_seed()
   call test_rk45_drydown()
   call test_rk45_tiny_cohort()
   call test_rk45_dense_cold()
   call test_rk45_dry_seam()
   call test_rnet_counts_wood()

   call test_report('test_fast_integ_step')

contains

   !=======================================================================================!
   ! Helpers.                                                                              !
   !=======================================================================================!

   !----- Run the column on a scheme. --------------------------------------------------------------!
   subroutine use_scheme(scheme)
      integer(ik), intent(in) :: scheme
      col%cfg%time_integrator = scheme
   end subroutine use_scheme

   character(len=7) function scheme_name(scheme) result(nm)
      integer(ik), intent(in) :: scheme
      nm = merge('ESDIRK2', 'RK45   ', scheme == INTEG_ESDIRK2)
   end function scheme_name

   !----- One day of the column from a fresh state: the fixture's 21 June forcing, step i at          !
   !      (i - 1/2)*DT_FAST after midnight, with the day's rain, snow, cold air, shed water and bottom  !
   !      boundary, recording the noon (13:30) and night (00:30) states and the day's extremes. --------!
   subroutine march_day(d)
      type(day_t), intent(inout) :: d
      real(wp)    :: t_sec, cosz, qsat, psi, t_air
      integer(ik) :: istep
      call reset_column(col, d%theta_seed)
      if (d%sealed) then
         col%site%soil%water_bottom_bc = SOIL_BC_DRAIN
         col%site%soil%drain_factor = 0.0_wp
      end if
      if (d%aquifer) col%site%soil%water_bottom_bc = SOIL_BC_AQUIFER
      col%static%shed_water_rate = d%shed
      col%static%shed_water_enth = d%shed * internal_energy_liquid(T0)   ! shed by tissues at T0
      !----- A pack, seeded as the model seeds its first one (250 kg/m3), at 270 K. ----------------!
      if (d%snow_seed > 0.0_wp) then
         col%dynamic%snow%swe(1)         = d%snow_seed
         col%dynamic%snow%snow_energy(1) = temp_to_internal_energy(0.0_wp, d%snow_seed, 270.0_wp, 0.0_wp)
         col%dynamic%snow%snow_depth(1)  = d%snow_seed / 250.0_wp
         col%dynamic%snow%nlayer         = 1_ik
      end if
      do istep = 1_ik, d%hours * STEPS_PER_HOUR
         t_sec = (real(istep, wp) - 0.5_wp) * DT_FAST
         cosz  = solar_cosz(SIM_DATE, t_sec, LAT)
         call set_diurnal_forcing(col, t_sec)
         t_air = col%forcing%air_temp
         call set_air(col%forcing, t_air, d%shv_air)
         if (istep >= 3_ik * STEPS_PER_HOUR .and. istep <= 7_ik * STEPS_PER_HOUR) col%forcing%rainfall = d%rain_pulse
         if (d%rain_all_day > 0.0_wp) col%forcing%rainfall = d%rain_all_day
         col%forcing%snowfall = d%snowfall
         if (d%snowfall > 0.0_wp .or. d%cold_air) call set_air(col%forcing, 268.0_wp + 3.0_wp * (cosz - 0.3_wp), d%shv_air)
         if (d%dark_step) then
            call set_air(col%forcing, T0, d%shv_air)
            col%forcing%abs_sw = 0.0_wp
            col%forcing%abs_par = 0.0_wp
            col%forcing%abs_sw_ground = 0.0_wp
            col%forcing%rainfall = 0.0_wp
         end if
         call set_ground_tiles(col)
         call step_column(col, DT_FAST)
         if (d%dark_step) then
            d%h_first   = col%out%h_flux
            d%gah_first = col%forcing%rho_air * col%out%aero%ustar * col%out%aero%heat_transfer_factor
            exit
         end if
         associate (dy => col%dynamic)
            d%ss_min = min(d%ss_min, dy%soil_energy%soil_temp(1))
            d%ss_max = max(d%ss_max, dy%soil_energy%soil_temp(1))
            d%sd_min = min(d%sd_min, dy%soil_energy%soil_temp(NSL))
            d%sd_max = max(d%sd_max, dy%soil_energy%soil_temp(NSL))
            d%th_min = min(d%th_min, dy%soil_water%theta(1))
            d%th_max = max(d%th_max, dy%soil_water%theta(1))
            d%theta_min  = min(d%theta_min, minval(dy%soil_water%theta(1:NSL)))
            d%theta_peak = max(d%theta_peak, maxval(dy%soil_water%theta(1:NSL)))
            d%film_peak  = max(d%film_peak, dy%leaf_surf_water(1) + dy%wood_surf_water(1))
            d%pond_peak  = max(d%pond_peak, dy%soil_water%w_surface)
            d%ood_peak   = max(d%ood_peak, col%out%report%theta_ood_max)
            d%wood_lag   = max(d%wood_lag, abs(dy%wood_temp(1) - dy%cas%can_temp))
            qsat = sat_specific_humidity(dy%cas%can_temp, col%forcing%press)
            d%worst_super = max(d%worst_super, dy%cas%can_shv - qsat)
            psi = tissue_water_to_psi(dy%leaf_water_mass(1), col%cfg%pft%hydraulics%pft(1)%leaf_curve,             &
                                      col%static%bleaf(1))
            d%physical = d%physical .and. dy%cas%can_temp > 270.0_wp .and. dy%cas%can_temp < 325.0_wp              &
                         .and. dy%cas%can_shv > 0.0_wp .and. dy%leaf_temp(1) > 250.0_wp .and. dy%leaf_temp(1) < 340.0_wp &
                         .and. all(dy%soil_water%theta(1:NSL) > 0.0_wp) .and. all(dy%soil_water%theta(1:NSL) < 0.6_wp) &
                         .and. all(dy%soil_energy%soil_temp(1:NSL) > 260.0_wp)                                   &
                         .and. all(dy%soil_energy%soil_temp(1:NSL) < 340.0_wp)                                   &
                         .and. dy%leaf_water_mass(1) > 0.0_wp .and. dy%wood_water_mass(1) > 0.0_wp               &
                         .and. psi < 0.5_wp .and. psi > -12.0_wp
            if (dy%snow%nlayer >= 1_ik) then
               d%snow_physical = d%snow_physical .and. dy%snow%swe(1) >= 0.0_wp .and. dy%snow%snow_temp(1) > 200.0_wp &
                                 .and. dy%snow%snow_temp(1) < 290.0_wp                                               &
                                 .and. all(dy%soil_energy%soil_temp(1:NSL) > 230.0_wp)                              &
                                 .and. all(dy%soil_energy%soil_temp(1:NSL) < 340.0_wp)
            end if
            if (istep == 13_ik * STEPS_PER_HOUR + STEPS_PER_HOUR / 2_ik) then          ! 13:30
               d%ct_noon     = dy%cas%can_temp
               d%tleaf_noon  = dy%leaf_temp(1)
               d%co2_noon    = dy%cas%can_co2
               d%ustar_noon  = col%out%aero%ustar
               d%dtheta_noon = dy%cas%can_temp - col%forcing%air_temp
               d%gpp_noon    = col%out%report%exports%gpp_last
               d%nee_noon    = col%out%report%exports%nee_last
               d%psileaf_noon = psi
            end if
            if (istep == STEPS_PER_HOUR / 2_ik) then                                     ! 00:30
               d%ct_night    = dy%cas%can_temp
               d%tleaf_night = dy%leaf_temp(1)
               d%co2_night   = dy%cas%can_co2
               d%psileaf_night = psi
            end if
         end associate
      end do
   end subroutine march_day

   !----- The column's water at the day's end: the soil layers, the pond and the pack [kg/m2]. -------!
   real(wp) function column_water() result(w)
      w = col%dynamic%soil_water%w_surface + col%dynamic%snow%swe(1)                                            &
          + sum(col%dynamic%soil_water%theta(1:NSL) * col%site%soil_profile%dz(1:NSL)) * rho_h2o
   end function column_water

   !----- The soil's and the pond's energy [J/m2]. -------------------------------------------------!
   real(wp) function ground_energy() result(e)
      e = sum(col%dynamic%soil_energy%soil_energy(1:NSL) * col%site%soil_profile%dz(1:NSL))                    &
          + col%dynamic%soil_water%w_surface_energy
   end function ground_energy

   !----- The heat a tissue gains per kelvin at T0 [J/m2/K], with dry heat capacity dry_hcap and     !
   !      holding wmass of water: the slope of its energy in temperature (tissue_temp_to_energy). ---!
   real(wp) function heat_per_kelvin(dry_hcap, wmass) result(c)
      real(wp), intent(in) :: dry_hcap   !< [J/m2/K]
      real(wp), intent(in) :: wmass      !< [kg/m2]
      c = tissue_temp_to_energy(dry_hcap, wmass, T0 + 1.0_wp) - tissue_temp_to_energy(dry_hcap, wmass, T0)
   end function heat_per_kelvin

   !=======================================================================================!
   ! 1. THE COLUMN: what the gathered view brings to the step.                               !
   !=======================================================================================!
   subroutine test_the_column()
      type(frozen_t)      :: fz
      type(march_state_t) :: ys
      real(wp) :: eps_leaf0(col%cfg%pft%n), eps_wood0(col%cfg%pft%n)
      print '(a)', 'test_the_column:'
      !----- The hydraulics take the sapwood area in m2; the cohort keeps it in cm2. Passed through     !
      !      unconverted it would be 10^4 times the stem's own cross-section. ------------------------!
      call check_true('the sapwood area reaches the hydraulics in m2, within the stem',                          &
                      col%static%sap_area(1) > 0.0_wp .and. col%static%sap_area(1) <= pio4 * 0.2_wp ** 2,         &
                      col%static%sap_area(1))
      !----- The PFT's aboveground fraction (stem respiration) reaches the view: PFT table ->          !
      !      init_cohort -> the cohort block -> the gather. -----------------------------------------!
      call check('the PFT''s aboveground fraction reaches the view', col%static%aboveground_frac(1),             &
                 col%cfg%pft%allom(1)%aboveground_frac, 1.0e-14_wp)
      !----- The energy balance's emission slope 4 eps sigma T^3 takes each cohort's PFT leaf and wood   !
      !      emissivity, the ones the radiation solver emits with. Distinct values catch a swap. -------!
      eps_leaf0 = col%cfg%pft%leaf_emissivity
      eps_wood0 = col%cfg%pft%wood_emissivity
      col%cfg%pft%leaf_emissivity = 0.93_wp
      col%cfg%pft%wood_emissivity = 0.87_wp
      call gather_column(col)
      call reset_column(col)
      call set_noon_forcing(col)
      call frozen_record(col, fz, ys)
      call check('build_frozen: the energy balance takes the PFT''s leaf emissivity', fz%tissue%leaf_emiss(1),   &
                 0.93_wp, 0.0_wp)
      call check('build_frozen: the energy balance takes the PFT''s wood emissivity', fz%tissue%wood_emiss(1),   &
                 0.87_wp, 0.0_wp)
      col%cfg%pft%leaf_emissivity = eps_leaf0
      col%cfg%pft%wood_emissivity = eps_wood0
      call gather_column(col)
   end subroutine test_the_column

   !=======================================================================================!
   ! 2. A DIURNAL DAY with a morning rain pulse (03:00-07:00), on ESDIRK2.                   !
   !=======================================================================================!
   subroutine test_diurnal_day()
      type(day_t) :: d
      print '(a)', 'test_diurnal_day:'
      call use_scheme(INTEG_ESDIRK2)
      d = day_t(rain_pulse = 1.5e-4_wp)
      call march_day(d)
      call check_true('DAY: all seven ledgers close on every step', all_ledgers_closed(col%out%report))
      call check_true('DAY: canopy air warmer near noon than at night', d%ct_noon > d%ct_night, d%ct_noon - d%ct_night)
      !----- The turbulence regime itself, which every invariant above passes in any regime: with the   !
      !      air at its type default (298.15 K) instead of the forcing's 282-294 K, the surface layer    !
      !      is stable all day. Correctly driven, midday is unstable (this fails by -0.56 K when the     !
      !      frozen build stops handing the aerodynamics the forcing's air temperature) and u* sits     !
      !      off its floor, which bounds the decoupled regime a fixed reference gives. ----------------!
      call check_true('DAY: midday is unstable (canopy air warmer than the air)', d%dtheta_noon > 0.0_wp,          &
                      d%dtheta_noon)
      call check_true('DAY: midday u* is off the surface layer''s floor', d%ustar_noon > 1.5_wp * col%cfg%aero%ustmin, &
                      d%ustar_noon)
      call check_true('DAY: the leaf warms under the sun', d%tleaf_noon > d%tleaf_night, d%tleaf_noon - d%tleaf_night)
      call check_true('DAY: the soil''s diurnal swing damps with depth', (d%ss_max - d%ss_min) > (d%sd_max - d%sd_min), &
                      (d%ss_max - d%ss_min) - (d%sd_max - d%sd_min))
      call check_true('DAY: the soil moisture responds to the rain', d%th_max - d%th_min > 1.0e-4_wp, d%th_max - d%th_min)
      call check_true('DAY: the soil moisture stays physical', d%th_min > 0.05_wp .and. d%th_max < 0.43_wp, d%th_max)
      call check_true('DAY: photosynthesis by day (GPP > 1)', d%gpp_noon > 1.0_wp, d%gpp_noon)
      call check_true('DAY: net uptake by day (NEE < 0)', d%nee_noon < 0.0_wp, d%nee_noon)
      call check_true('DAY: canopy CO2 lower at midday than at night', d%co2_noon < d%co2_night, d%co2_night - d%co2_noon)
      call check_true('DAY: the leaf water is under tension by day', d%psileaf_noon < 0.0_wp, d%psileaf_noon)
      call check_true('DAY: the leaf is under more tension at midday than at night', d%psileaf_noon < d%psileaf_night, &
                      d%psileaf_noon - d%psileaf_night)
      !----- With the interior faces advecting the water's enthalpy, the soil stays bounded at both     !
      !      ends; with the faces lumped this day reaches a 361 K soil surface while the whole-column   !
      !      ledger still closes, the error being purely vertical (docs/science/numerical_scheme.md §5b, !
      !      "Every soil-water flow comes from one solve"). -------------------------------------------!
      call check_true('DAY: the soil surface stays bounded', d%ss_min > 270.0_wp .and. d%ss_max < 330.0_wp, d%ss_max)
      call check_true('DAY: the deep soil stays bounded', d%sd_min > 270.0_wp .and. d%sd_max < 330.0_wp, d%sd_max)
   end subroutine test_diurnal_day

   !=======================================================================================!
   ! 3. THE REFERENCE LEVEL. The surface layer compares the canopy air with the air at the   !
   ! canopy-air top (zref), both referenced to that level, and the reported sensible heat is !
   ! the flux the whole-column ledger exchanges. A dark column at the air's own temperature  !
   ! has nothing driving a sensible flux, so its first step reports |H| far below            !
   ! g_atm_heat*g*zref; referencing only the air's potential temperature to the ground would !
   ! report about -g_atm_heat*g*zref (docs/science/canopy_aerodynamics.md §2).               !
   !=======================================================================================!
   subroutine test_reference_level()
      type(day_t) :: d
      integer(ik) :: is
      print '(a)', 'test_reference_level:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         d = day_t(dark_step = .true.)
         call march_day(d)
         print '(3a,f8.3,a,f8.3,a)', '   ', scheme_name(SCHEMES(is)), ' H = ', d%h_first, ' W/m2; g_atm_heat g zref = ', &
               d%gah_first * grav * col%static%z_top, ' W/m2'
         call check_true(trim(scheme_name(SCHEMES(is)))//': a dark isothermal column reports |H| < 5% of '//          &
                         'g_atm_heat*g*zref', abs(d%h_first) < 0.05_wp * d%gah_first * grav * col%static%z_top, d%h_first)
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_reference_level

   !=======================================================================================!
   ! 4. ESDIRK2 AND RK45 AGREE on a forced day: the fixture's 21 June with a morning rain on  !
   ! the canopy films, at the shipped dt_fast of 900 s, each scheme at RK45's shipped         !
   ! tolerance. Both close every ledger, and the day's totals agree.                          !
   !                                                                                          !
   ! THE TOLERANCE. Both march the same frozen record, water network and update; they differ !
   ! only in how the canopy air, the soil heat, the tissues and the films are integrated      !
   ! within a step. Each holds every sub-step's local error in the canopy air and the soil    !
   ! heat below one tolerance unit, atol + rtol|y| (eq. (3)): at rtol = 1e-4 that is          !
   ! 50 + 1e-4 * 2.9e5 ~ 80 J/kg of canopy-air enthalpy, 0.08 K. The canopy air forgets an    !
   ! error within its relaxation time (minutes), so the two trajectories stay within about    !
   ! one unit of each other instead of adding up over the day's 96 steps. One unit moves the  !
   ! sensible heat by g_atm cp * 0.08 K, ~7 W/m2 at this column's night-time coupling of      !
   ! ~90 W/m2/K (more by day), 5-6 % of the day's mean |H| (~120 W/m2), and the other fluxes  !
   ! in proportion. So each day's total may differ by up to 5 % of the day's integrated      !
   ! |flux|, the bound asserted. The schemes in fact agree far closer (a few 1e-4), while a   !
   ! wiring error in either (a dropped flux term, a stage fed the wrong input) moves a total  !
   ! by tens of percent. ESDIRK2's films are advanced by a map and RK45's in its tableau, so  !
   ! the latent heat and the transpiration are where the two differ most.                     !
   !=======================================================================================!
   subroutine test_schemes_agree()
      real(wp),    parameter :: H = 900.0_wp, AGREE = 0.05_wp
      integer(ik), parameter :: NSTEP = nint(86400.0_wp / H, ik), NQ = 7_ik
      character(len=*), parameter :: QNAME(NQ) = [character(len=24) :: 'sensible heat', 'latent heat',             &
                                                   'net radiation', 'GPP', 'NEE', 'ground heat', 'transpiration']
      real(wp)    :: total(NQ, 2), gross(NQ, 2), amount(NQ), e_ground, t_sec, rtol0
      integer(ik) :: is, istep, q, n
      print '(a)', 'test_schemes_agree:'
      n = col%static%n
      rtol0 = col%cfg%rtol
      col%cfg%canopy_water_on = .true.
      col%cfg%rtol = RK45_RTOL
      call ensure_output_capacity(col%out, n, .true.)
      total = 0.0_wp
      gross = 0.0_wp
      do is = 1_ik, 2_ik
         call use_scheme(SCHEMES(is))
         call reset_column(col)
         do istep = 1_ik, NSTEP
            t_sec = (real(istep, wp) - 0.5_wp) * H
            call set_diurnal_forcing(col, t_sec)
            if (t_sec >= 3.0_wp * 3600.0_wp .and. t_sec <= 7.0_wp * 3600.0_wp) col%forcing%rainfall = 1.5e-4_wp
            e_ground = ground_energy()
            call step_column(col, H)
            !----- The step's amounts: sensible and latent heat and net radiation [J/m2], GPP and NEE     !
            !      [umol/m2], the soil and pond's heat gain [J/m2], transpiration [kg/m2]. -----------!
            amount(1) = col%out%h_flux * H
            amount(2) = col%out%le_flux * H
            amount(3) = col%out%rnet * H
            amount(4) = col%out%report%exports%gpp_last * H
            amount(5) = col%out%report%exports%nee_last * H
            amount(6) = ground_energy() - e_ground
            amount(7) = sum(col%out%cdiag(CO_TRANSP, 1:n) * col%static%nplant(1:n) * col%static%leaf_area(1:n))   &
                        * mmh2o * H
            total(:, is) = total(:, is) + amount
            gross(:, is) = gross(:, is) + abs(amount)
         end do
         call check_true(trim(scheme_name(SCHEMES(is)))//' at 900 s: all seven ledgers close on every step',          &
                         all_ledgers_closed(col%out%report))
      end do
      do q = 1_ik, NQ
         print '(a,a24,a,es12.4,a,es12.4,a,f8.4)', '   ', QNAME(q), ' ESDIRK2 ', total(q, 1), ' RK45 ', total(q, 2),   &
               ' diff / |flux| ', abs(total(q, 1) - total(q, 2)) / max(gross(q, 1), gross(q, 2), tiny(1.0_wp))
         call check_true('ESDIRK2 ~ RK45 over a forced day: '//trim(QNAME(q))//' within 5% of the day''s |flux|',     &
                         abs(total(q, 1) - total(q, 2)) <= AGREE * max(gross(q, 1), gross(q, 2)),                  &
                         (total(q, 1) - total(q, 2)) / max(gross(q, 1), gross(q, 2), tiny(1.0_wp)))
      end do
      call ensure_output_capacity(col%out, n, .false.)
      col%cfg%canopy_water_on = .false.
      col%cfg%rtol = rtol0
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_schemes_agree

   !=======================================================================================!
   ! 5. THE LEDGERS ON EVERY PATH.                                                            !
   !=======================================================================================!

   !----- A dry day over moist free-draining soil, at RK45's shipped tolerance, on ESDIRK2 and     !
   !      RK45: all seven ledgers close on every dispatched step, the canopy air's energy and the   !
   !      whole column's to round-off, and the state stays physical and the canopy air below       !
   !      saturation all day. -------------------------------------------------------------------!
   subroutine test_dry_day()
      type(day_t) :: d
      character(len=16) :: tag
      real(wp)    :: rtol0
      integer(ik) :: ir
      print '(a)', 'test_dry_day:'
      rtol0 = col%cfg%rtol
      col%cfg%rtol = RK45_RTOL
      do ir = 1_ik, 2_ik
         select case (ir)
         case (1_ik)
            call use_scheme(INTEG_ESDIRK2)
            tag = 'ESDIRK2'
         case default
            call use_scheme(INTEG_RK45)
            tag = 'RK45'
         end select
         d = day_t()
         call march_day(d)
         call check_true(trim(tag)//' dry day: all seven ledgers close', all_ledgers_closed(col%out%report))
         call check_true(trim(tag)//' dry day: the ledgers ran on every dispatched dt_fast',                       &
                         col%out%report%ledgers%cas_energy%n_check == STEPS_PER_DAY,                                &
                         real(col%out%report%ledgers%cas_energy%n_check, wp))
         call check_true(trim(tag)//' dry day: canopy-air energy closes to round-off',                            &
                         col%out%report%ledgers%cas_energy%worst < 1.0e-3_wp, col%out%report%ledgers%cas_energy%worst)
         call check_true(trim(tag)//' dry day: whole-column energy closes within 1 J',                            &
                         col%out%report%ledgers%whole_energy%worst < 1.0_wp, col%out%report%ledgers%whole_energy%worst)
         call check_true(trim(tag)//' dry day: the state stays physical (canopy air, leaf, soil, plant water)',    &
                         d%physical)
         call check_true(trim(tag)//' dry day: the canopy air stays below saturation', d%worst_super <= 1.0e-4_wp,  &
                         d%worst_super)
         print '(3a,i0,a,es10.3,a,es10.3)', '   ', trim(tag), ' last step: sub-steps = ',                            &
               col%out%report%work(WORK_INTEG_STEPS), '; worst whole-column energy ',                                &
               col%out%report%ledgers%whole_energy%worst, ' J/m2, water ', col%out%report%ledgers%whole_water%worst
      end do
      col%cfg%rtol = rtol0
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_dry_day

   !----- Steady rain all day (8e-5 kg/m2/s, ~0.29 mm/h) over free-draining soil, on both schemes:    !
   !      the boundary water-enthalpy advection under a real rainfall rate. Free drainage empties a   !
   !      wet column whatever falls on it, so the rain wets the soil when the column ends wetter than !
   !      the same day without it. ---------------------------------------------------------------!
   subroutine test_wet_day()
      type(day_t) :: d
      real(wp)    :: theta_dry
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_wet_day:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t()
         call march_day(d)
         theta_dry = sum(col%dynamic%soil_water%theta(1:NSL))
         d = day_t(rain_all_day = 8.0e-5_wp)
         call march_day(d)
         call check_true(trim(nm)//' wet day: the ledgers ran on every step (no guard stopped it)',                &
                         col%out%report%ledgers%whole_water%n_check == STEPS_PER_DAY,                               &
                         real(col%out%report%ledgers%whole_water%n_check, wp))
         call check_true(trim(nm)//' wet day: the rain wetted the soil (wetter than the dry day)',                 &
                         sum(col%dynamic%soil_water%theta(1:NSL)) > theta_dry,                                     &
                         sum(col%dynamic%soil_water%theta(1:NSL)) - theta_dry)
         call check_true(trim(nm)//' wet day: all seven ledgers close', all_ledgers_closed(col%out%report))
         call check_true(trim(nm)//' wet day: whole-column energy closes within 1 J',                             &
                         col%out%report%ledgers%whole_energy%worst < 1.0_wp, col%out%report%ledgers%whole_energy%worst)
         call check_true(trim(nm)//' wet day: soil energy closes to round-off',                                   &
                         col%out%report%ledgers%soil_energy%worst < 1.0e-3_wp, col%out%report%ledgers%soil_energy%worst)
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_wet_day

   !----- The day's leaf and root turnover shed water (8e-5 kg/m2/s), a patch input distinct from     !
   !      the rain, which stays 0: it alone wets the soil on both schemes, and both whole-column       !
   !      ledgers close with no wiring of its own, its enthalpy riding the infiltration's. ---------!
   subroutine test_shed_water()
      type(day_t) :: d
      real(wp)    :: theta_dry
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_shed_water:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t()
         call march_day(d)
         theta_dry = sum(col%dynamic%soil_water%theta(1:NSL))
         d = day_t(shed = 8.0e-5_wp)
         call march_day(d)
         call check_true(trim(nm)//' shed water: it alone wetted the soil (wetter than without it)',               &
                         sum(col%dynamic%soil_water%theta(1:NSL)) > theta_dry,                                     &
                         sum(col%dynamic%soil_water%theta(1:NSL)) - theta_dry)
         call check_true(trim(nm)//' shed water: whole-column water closes',                                      &
                         col%out%report%ledgers%whole_water%n_fail == 0_ik, real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' shed water: whole-column energy closes',                                     &
                         col%out%report%ledgers%whole_energy%n_fail == 0_ik,                                       &
                         real(col%out%report%ledgers%whole_energy%n_fail, wp))
         call check_true(trim(nm)//' shed water: whole-column energy closes within 1 J',                          &
                         col%out%report%ledgers%whole_energy%worst < 1.0_wp, col%out%report%ledgers%whole_energy%worst)
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_shed_water

   !----- The canopy films ([fast].canopy_water_on) under the morning rain, on both schemes: the      !
   !      canopy intercepts some of it, and the whole column's water and energy close, each film     !
   !      being part of its tissue's energy at the tissue's temperature, its capture brought in at   !
   !      the rain's liquid enthalpy and its evaporation paying the vapour's. The wiring, from the    !
   !      interception through the films to the canopy air and the ledgers; the films' own laws are  !
   !      test_canopy_water's and test_plant_energy's. ------------------------------------------!
   subroutine test_canopy_films()
      type(day_t) :: d
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_canopy_films:'
      col%cfg%canopy_water_on = .true.
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t(rain_pulse = 1.5e-4_wp)
         call march_day(d)
         print '(3a,es10.3,a)', '   ', trim(nm), ' peak film ', d%film_peak, ' kg/m2'
         call check_true(trim(nm)//' films: whole-column water closes', col%out%report%ledgers%whole_water%n_fail == 0_ik, &
                         real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' films: the morning rain was intercepted', d%film_peak > 0.0_wp, d%film_peak)
         call check_true(trim(nm)//' films: whole-column energy closes to round-off',                             &
                         col%out%report%ledgers%whole_energy%worst < 1.0e-3_wp, col%out%report%ledgers%whole_energy%worst)
      end do
      col%cfg%canopy_water_on = .false.
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_canopy_films

   !----- A sealed (bedrock) column seeded just below theta_sat under ~29 mm/h of rain, on both      !
   !      schemes: it saturates (the water network raises the layer's potential and drives the       !
   !      excess to the pond and the next layer) and overflows its pond (surface runoff), paths no    !
   !      other scenario takes. Both move mass and enthalpy as a pair: the pond holds its own        !
   !      enthalpy and runoff leaves at its temperature (docs/science/numerical_scheme.md §5b, "The  !
   !      pond holds heat"). -----------------------------------------------------------------------!
   subroutine test_sealed_saturated()
      type(day_t) :: d
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_sealed_saturated:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t(theta_seed = 0.428_wp, rain_all_day = 8.0e-3_wp, sealed = .true.)
         call march_day(d)
         call check_true(trim(nm)//' sealed: the column reached theta_sat (the saturation path is live)',          &
                         d%theta_peak >= 0.43_wp - 1.0e-12_wp, d%theta_peak)
         call check_true(trim(nm)//' sealed: the pond overflowed (the runoff path is live)',                       &
                         d%pond_peak >= col%site%soil%w_pond_max - 1.0e-9_wp, d%pond_peak)
         !----- Energy as a bound: nothing clamps the committed state, but whether a stage clamp bites  !
         !      depends on the trajectory and the compiler; a clamp that rebuilt the committed soil      !
         !      energy with no ledger term would show here (1.5e5 J/m2 on nvfortran). ------------------!
         call check_true(trim(nm)//' sealed: whole-column energy closes (no unbooked clamp on the commit)',        &
                         col%out%report%ledgers%whole_energy%worst < 1.0e-3_wp, col%out%report%ledgers%whole_energy%worst)
         !----- Water to machine precision: drainage, ponding and runoff come from the water network solve  !
         !      whose theta both schemes commit, so every flow and the state it moves agree. ----------!
         call check_true(trim(nm)//' sealed: whole-column water closes', col%out%report%ledgers%whole_water%n_fail == 0_ik, &
                         real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' sealed: whole-column water closes to machine precision',                     &
                         col%out%report%ledgers%whole_water%worst < 1.0e-6_wp, col%out%report%ledgers%whole_water%worst)
         !----- The vertical check the whole-column ledgers cannot make: each layer's water change against !
         !      the faces the soil energy is charged for, the network's own per-layer closure at the       !
         !      saturated corner where its pond and runoff paths are live. A misplacement between layers  !
         !      still sums correctly against the boundary, even at a 345 K soil surface. ---------------!
         call check_true(trim(nm)//' sealed: per-layer face closure holds',                                       &
                         col%out%report%ledgers%soil_face_mass%worst < 1.0e-8_wp,                                    &
                         col%out%report%ledgers%soil_face_mass%worst)
         !----- The soil surface sees heat advected on faces that are not the scheme's own, which no ledger !
         !      sees: one such error moves it by ~50 K; with the faces right it sits near 295 K. A 29 mm/h   !
         !      rain cannot warm a soil surface past ~300 K, so 310 K has diagnostic power. --------------!
         call check_true(trim(nm)//' sealed: the soil surface stays physical (the faces on the committed theta)',    &
                         d%ss_min > 250.0_wp .and. d%ss_max < 310.0_wp, d%ss_max)
         !----- A saturated layer holds more only through its specific storage, at a positive head, at most !
         !      the hydrostatic one in a sealed column: the column's depth plus the pond (docs/science/      !
         !      numerical_scheme.md §5b, "Saturation is the water network's"). ---------------------------!
         call check_true(trim(nm)//' sealed: theta above theta_sat only through specific storage (hydrostatic head)', &
                         d%theta_peak <= 0.43_wp + SOIL_SPECIFIC_STORAGE                                             &
                                       * (-col%site%soil_profile%soil_layer_z(NSL + 1_ik)                            &
                                          + col%site%soil%w_pond_max / rho_h2o), d%theta_peak)
         !----- How far an RK45 stage's starting theta lay outside [theta_res, theta_sat] (ESDIRK2 has no    !
         !      explicit stage): under the water network only above saturation, by the specific storage,  !
         !      so the bound is far looser than what is reached (below 1e-4). --------------------------!
         if (SCHEMES(is) == INTEG_RK45) then
            print '(a,es10.3)', '   RK45 sealed: worst stage excursion of theta ', d%ood_peak
            call check_true('RK45 sealed: the constitutive-domain excursion stays small', d%ood_peak < 2.0e-2_wp,       &
                            d%ood_peak)
         end if
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_sealed_saturated

   !----- A free-draining column seeded just below theta_sat under ~29 mm/h of rain, far above its   !
   !      drainage (ksat 2.89e-6 m/s, ~2.9e-3 kg/m2/s), saturates and overflows its pond too (on      !
   !      ESDIRK2, the water network's route being the same on both). ----------------------------!
   subroutine test_free_drain_saturated()
      type(day_t) :: d
      print '(a)', 'test_free_drain_saturated:'
      call use_scheme(INTEG_ESDIRK2)
      d = day_t(theta_seed = 0.428_wp, rain_all_day = 8.0e-3_wp)
      call march_day(d)
      call check_true('free drain: the column reached theta_sat', d%theta_peak >= 0.43_wp - 1.0e-12_wp, d%theta_peak)
      call check_true('free drain: the pond overflowed', d%pond_peak >= col%site%soil%w_pond_max - 1.0e-9_wp, d%pond_peak)
      call check_true('free drain: whole-column energy closes through saturation and runoff',                      &
                      col%out%report%ledgers%whole_energy%n_fail == 0_ik, col%out%report%ledgers%whole_energy%worst)
      call check_true('free drain: whole-column water closes', col%out%report%ledgers%whole_water%n_fail == 0_ik,     &
                      col%out%report%ledgers%whole_water%worst)
      call check_true('free drain: the soil surface stays physical', d%ss_min > 250.0_wp .and. d%ss_max < 340.0_wp,  &
                      d%ss_max)
   end subroutine test_free_drain_saturated

   !----- The aquifer bottom is a head-driven, two-way boundary with no state behind it, and both     !
   !      schemes commit the water network's theta, so each takes it like any other face: a column    !
   !      started dry (0.15) wets from below over 12 h with both whole-column ledgers closed. -------!
   subroutine test_aquifer()
      type(day_t) :: d
      integer(ik) :: is, nb
      character(len=7) :: nm
      print '(a)', 'test_aquifer:'
      nb = col%site%soil_profile%n_active
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t(theta_seed = 0.15_wp, aquifer = .true., hours = 12_ik)
         call march_day(d)
         call check_true(trim(nm)//' aquifer: whole-column water closes', col%out%report%ledgers%whole_water%n_fail == 0_ik, &
                         real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' aquifer: whole-column energy closes',                                         &
                         col%out%report%ledgers%whole_energy%n_fail == 0_ik,                                        &
                         real(col%out%report%ledgers%whole_energy%n_fail, wp))
         call check_true(trim(nm)//' aquifer: the dry column wets from below', col%dynamic%soil_water%theta(nb) > 0.15_wp, &
                         col%dynamic%soil_water%theta(nb) - 0.15_wp)
         print '(3a,i0,a,i0)', '   ', trim(nm), ' aquifer, last step: sub-steps = ', col%out%report%work(WORK_INTEG_STEPS), &
               ', rescues = ', col%out%report%work(WORK_RK45_RESCUE)
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_aquifer

   !=======================================================================================!
   ! 6. SNOW. The pack's coupling, integrated: the stage both schemes share                  !
   ! (prepare_snow_stage: accumulate, then drain the meltwater), the surface balance against !
   ! the live canopy air in the march, the snowfac-blended ground, sublimation into the       !
   ! canopy air, the commit, and the pack's terms in the whole-column ledgers. The snow       !
   ! kernels alone are test_surface's.                                                        !
   !=======================================================================================!

   !----- A 60 kg/m2 pack under 2e-5 kg/m2/s of snowfall in sub-freezing air, no rain: it melts over   !
   !      the day (75 W/m2 of ground shortwave and no longwave loss beat the snowfall) without        !
   !      vanishing, so the snowfac-blended surface is live all day and melt routes a paired mass and  !
   !      enthalpy into the soil top. A scheme that imports the snow kernels but never calls them runs !
   !      snow-free and passes its ledgers trivially; RK45's melt total agrees with ESDIRK2's, within  !
   !      10 % for the schemes driving the surface along different canopy-air paths, while a stage    !
   !      fed the wrong inputs is off by order 100 %. -----------------------------------------------!
   subroutine test_snow_melt()
      type(day_t) :: d
      real(wp)    :: melt_esdirk2
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_snow_melt:'
      melt_esdirk2 = 0.0_wp
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t(snow_seed = 60.0_wp, snowfall = 2.0e-5_wp)
         call march_day(d)
         print '(3a,f7.3,a,es10.3)', '   ', trim(nm), ' swe = ', col%dynamic%snow%swe(1),                           &
               ' kg/m2; worst water residual ', col%out%report%ledgers%whole_water%worst
         call check_true(trim(nm)//' snow: the pack and the soil stay physical', d%snow_physical,                   &
                         col%dynamic%snow%snow_temp(1))
         call check_true(trim(nm)//' snow: melt and sublimation drained the pack without exhausting it',            &
                         col%dynamic%snow%swe(1) < 60.0_wp .and. col%dynamic%snow%swe(1) > 0.0_wp, col%dynamic%snow%swe(1))
         call check_true(trim(nm)//' snow: whole-column water closes with a pack',                                 &
                         col%out%report%ledgers%whole_water%n_fail == 0_ik, real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' snow: whole-column energy closes to round-off with a melting pack',           &
                         col%out%report%ledgers%whole_energy%worst < 1.0e-3_wp, col%out%report%ledgers%whole_energy%worst)
         if (SCHEMES(is) == INTEG_ESDIRK2) then
            melt_esdirk2 = 60.0_wp - col%dynamic%snow%swe(1)
         else
            call check_true('RK45 snow: the melt total agrees with ESDIRK2''s (one shared stage, the same inputs)',  &
                            abs((60.0_wp - col%dynamic%snow%swe(1)) - melt_esdirk2) < 0.10_wp * melt_esdirk2,       &
                            (60.0_wp - col%dynamic%snow%swe(1)) - melt_esdirk2)
         end if
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_snow_melt

   !----- The day's turnover shed water reaches the pond with the meltwater, at the melt temperature,  !
   !      whenever a pack exists, and both whole-column ledgers book that inflow's enthalpy under the  !
   !      pack too. Left unbooked it is a constant residual of shed*dt*u_liq(t_melt) every winter step  !
   !      (~1.8 J/m2 at 150 s in the Ithaca example, +6.75 MJ/m2 over a 50-year spin-up), small enough  !
   !      to pass the tolerance, so closure is asked to round-off, on both schemes. -----------------!
   subroutine test_shed_under_pack()
      type(day_t) :: d
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_shed_under_pack:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t(snow_seed = 60.0_wp, shed = 2.0e-8_wp)
         call march_day(d)
         print '(3a,es10.3,a,f7.3,a)', '   ', trim(nm), ' worst energy residual ',                                   &
               col%out%report%ledgers%whole_energy%worst, ' J/m2; swe ', col%dynamic%snow%swe(1), ' kg/m2'
         call check_true(trim(nm)//' shed under a pack: whole-column water closes',                                &
                         col%out%report%ledgers%whole_water%n_fail == 0_ik, real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' shed under a pack: whole-column energy closes to round-off',                  &
                         col%out%report%ledgers%whole_energy%worst < 1.0e-3_wp, col%out%report%ledgers%whole_energy%worst)
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_shed_under_pack

   !----- Snowfall is a boundary water input like rain, and every kilogram of it must end in the      !
   !      column: the pack, the soil or the pond. A frozen build that routes only the rainfall loses   !
   !      snowfall*dt every sub-freezing step while the whole-column ledger counts it as an input, and  !
   !      per step that sits below the ledger's tolerance even as it adds up to a large seasonal error. !
   !      So the assertion is a boundary-input identity, not a ledger residual: the same cold day with  !
   !      a seeded pack, snowfall on and off, so melt, sublimation and drainage cancel in the           !
   !      difference, and the column must gain what fell. A scheme that discards it gains ~0; 10 %      !
   !      allows for the second-order response of melt and drainage to the deeper pack. RK45's gain    !
   !      agrees with ESDIRK2's. -------------------------------------------------------------------!
   subroutine test_snowfall_conserved()
      real(wp), parameter :: SNOWFALL = 2.0e-5_wp
      type(day_t) :: d
      real(wp)    :: cw_on, cw_off, gain, gain_esdirk2, total
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_snowfall_conserved:'
      total = SNOWFALL * real(STEPS_PER_DAY, wp) * DT_FAST
      gain_esdirk2 = 0.0_wp
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t(snow_seed = 60.0_wp, snowfall = SNOWFALL)
         call march_day(d)
         cw_on = column_water()
         call check_true(trim(nm)//' snowfall: whole-column water closes',                                         &
                         col%out%report%ledgers%whole_water%n_fail == 0_ik, real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' snowfall: whole-column energy closes',                                        &
                         col%out%report%ledgers%whole_energy%n_fail == 0_ik,                                        &
                         real(col%out%report%ledgers%whole_energy%n_fail, wp))
         d = day_t(snow_seed = 60.0_wp, cold_air = .true.)
         call march_day(d)
         cw_off = column_water()
         gain = cw_on - cw_off
         print '(3a,f8.4,a,f8.4,a)', '   ', trim(nm), ' column water gain ', gain, ' kg/m2; snowfall ', total, ' kg/m2'
         call check_true(trim(nm)//' snowfall: what fell arrived in the column', abs(gain - total) < 0.10_wp * total, gain)
         if (SCHEMES(is) == INTEG_ESDIRK2) then
            gain_esdirk2 = gain
         else
            call check_true('RK45 snowfall: the gain agrees with ESDIRK2''s', abs(gain - gain_esdirk2) < 0.10_wp * total, &
                            gain - gain_esdirk2)
         end if
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_snowfall_conserved

   !----- Snowfall too light to start a pack (snowfall*dt below min_new_snow_mass) reaches the ground   !
   !      with the rain, as ice at the air's temperature. Valued as liquid at the canopy air's, the   !
   !      column would get the snow's fusion enthalpy from nowhere: consistent with its ledger (the    !
   !      boundary term uses the same number), and wrong. The same cold day twice, the mass as light   !
   !      snow and as rain: the soil and pond end the snow day poorer by about L_f per kilogram, the    !
   !      ground paying for the melt; liquid valuation gives about 0. On this column the ratio is     !
   !      ~0.5: the colder, wetter surface loses less to the air and recovers part of the melt's      !
   !      cost. ---------------------------------------------------------------------------------!
   subroutine test_light_snowfall()
      real(wp), parameter :: SNOWFALL = 5.0e-6_wp   ! 7.5e-4 kg/m2 per 150 s step, below min_new_snow_mass (1e-3)
      type(day_t) :: d
      real(wp)    :: e_snow, e_rain, fusion
      print '(a)', 'test_light_snowfall:'
      call use_scheme(INTEG_ESDIRK2)
      d = day_t(snowfall = SNOWFALL)
      call march_day(d)
      e_snow = ground_energy()
      call check_true('light snowfall: it never started a pack', col%dynamic%snow%nlayer == 0_ik,                    &
                      real(col%dynamic%snow%nlayer, wp))
      call check_true('light snowfall: whole-column energy closes to round-off',                                   &
                      col%out%report%ledgers%whole_energy%worst < 1.0e-3_wp, col%out%report%ledgers%whole_energy%worst)
      d = day_t(cold_air = .true., rain_all_day = SNOWFALL)
      call march_day(d)
      e_rain = ground_energy()
      fusion = latent_heat_fusion * SNOWFALL * real(STEPS_PER_DAY, wp) * DT_FAST
      print '(a,es10.3,a,es10.3,a)', '   soil and pond energy, rain day - snow day = ', e_rain - e_snow,              &
            ' J/m2; L_f * snow = ', fusion, ' J/m2'
      call check_true('light snowfall: the snow day leaves the ground poorer than the rain day by ~L_f per kg',     &
                      e_rain - e_snow > 0.3_wp * fusion .and. e_rain - e_snow < 1.5_wp * fusion,                    &
                      (e_rain - e_snow) / fusion)
   end subroutine test_light_snowfall

   !=======================================================================================!
   ! 7. PLANT WATER AND HEAT.                                                                !
   !=======================================================================================!

   !----- The leaf potential converges in dt_fast, and the plant's total water stays dt-insensitive.  !
   !      The leaf is the column's smallest store, so any mismatch between the transpiration charged   !
   !      to it and the sapflow credited to it lands there and grows with dt: an explicit plant-water   !
   !      step charging the end-of-step transpiration against the start-of-step sapflow takes the      !
   !      midday leaf from -0.33 MPa at 5 s to -1.79 MPa at 900 s while every other state converges.   !
   !      The water network solves the leaf implicitly with the soil, and the stages hold              !
   !      transpiration at the network's, so the leaf is charged and credited the same step            !
   !      (docs/science/numerical_scheme.md §5a). No ledger sees the defect: a plant store absorbing a  !
   !      flux mismatch is a closed budget. The same 3 h from solar noon at 75 s and 900 s, each step   !
   !      forced at its own mid-time so both see the same forcing: the defect parts them by ~1 MPa,    !
   !      0.05 MPa is far looser than the ~0.005 achieved, and the window is in daylight, where the     !
   !      transpiration that drives the defect is. --------------------------------------------------!
   subroutine test_psi_dt_convergence()
      real(wp) :: psi_short, psi_long, w_short, w_long
      print '(a)', 'test_psi_dt_convergence:'
      call use_scheme(INTEG_ESDIRK2)
      call run_window(75.0_wp,  psi_short, w_short)
      call run_window(900.0_wp, psi_long,  w_long)
      call check_true('psi: the leaf potential converges in dt_fast (75 s vs 900 s over 3 h)',                      &
                      abs(psi_long - psi_short) < 0.05_wp, abs(psi_long - psi_short))
      call check_true('psi: the plant''s total water is dt-insensitive (set by uptake and transpiration)',          &
                      abs(w_long - w_short) / max(w_short, tiny(1.0_wp)) < 0.05_wp,                                 &
                      abs(w_long - w_short) / max(w_short, tiny(1.0_wp)))
   end subroutine test_psi_dt_convergence

   subroutine run_window(dt, psi_out, w_out)
      real(wp), intent(in)  :: dt
      real(wp), intent(out) :: psi_out, w_out
      real(wp), parameter   :: T_START = 43200.0_wp, T_WINDOW = 10800.0_wp   ! solar noon, 3 h
      real(wp) :: t
      call reset_column(col)
      t = 0.0_wp
      do while (t < T_WINDOW - 1.0e-6_wp)
         call set_diurnal_forcing(col, T_START + t + 0.5_wp * dt)
         call step_column(col, dt)
         t = t + dt
      end do
      psi_out = tissue_water_to_psi(col%dynamic%leaf_water_mass(1), col%cfg%pft%hydraulics%pft(1)%leaf_curve,      &
                                    col%static%bleaf(1))
      w_out = sum(col%dynamic%leaf_water_mass(1:col%static%n) + col%dynamic%wood_water_mass(1:col%static%n))
   end subroutine run_window

   !----- The wood's heat store: tissue_energy_balance relaxes it exactly (docs/science/              !
   !      vegetation_energy_dynamics.md eqs. (3)-(4)) on both schemes; the whole-column energy closes   !
   !      with its store change and its net radiation; the wood lags the canopy air through the day,    !
   !      which is the feature; and the network's per-layer face closure holds on the faces both        !
   !      schemes advect soil heat on. -------------------------------------------------------------!
   subroutine test_prognostic_wood()
      type(day_t) :: d
      integer(ik) :: is
      character(len=7) :: nm
      print '(a)', 'test_prognostic_wood:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         nm = scheme_name(SCHEMES(is))
         d = day_t()
         call march_day(d)
         call check_true(trim(nm)//' wood store: whole-column energy closes',                                       &
                         col%out%report%ledgers%whole_energy%n_fail == 0_ik,                                        &
                         real(col%out%report%ledgers%whole_energy%n_fail, wp))
         call check_true(trim(nm)//' wood store: whole-column water closes',                                        &
                         col%out%report%ledgers%whole_water%n_fail == 0_ik, real(col%out%report%ledgers%whole_water%n_fail, wp))
         call check_true(trim(nm)//' wood store: per-layer face closure holds',                                     &
                         col%out%report%ledgers%soil_face_mass%worst < 1.0e-8_wp,                                    &
                         col%out%report%ledgers%soil_face_mass%worst)
         call check_true(trim(nm)//' wood store: the wood lags the canopy air', d%wood_lag > 1.0e-3_wp, d%wood_lag)
         call check_true(trim(nm)//' wood store: the wood temperature stays physical',                              &
                         col%dynamic%wood_temp(1) > 200.0_wp .and. col%dynamic%wood_temp(1) < 350.0_wp,              &
                         col%dynamic%wood_temp(1))
         print '(3a,i0,a,i0)', '   ', trim(nm), ' wood store, last step: sub-steps = ',                                &
               col%out%report%work(WORK_INTEG_STEPS), ', rescues = ', col%out%report%work(WORK_RK45_RESCUE)
      end do
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_prognostic_wood

   !----- Wood refilling from the soil does not heat itself. Dried wood over moist soil at night, with !
   !      the soil, the wood, the leaf and the air at one temperature, refills, and the water it takes !
   !      arrives at its own temperature, so the wood's heat input is zero. Counting the water's whole  !
   !      enthalpy (~4.4e5 J/kg from the liquid datum) as heat put the wood 13 K above the canopy air at !
   !      a storm front (BCI, 2017-04-17). The water stores take that enthalpy instead, and what the    !
   !      soil gives up is what the tissues and their water gain (build_frozen, record_water_flows). ---!
   subroutine test_refilling_wood()
      type(frozen_t)      :: fz
      type(march_state_t) :: ys
      real(wp) :: refill, gained, given
      integer(ik) :: n
      print '(a)', 'test_refilling_wood:'
      n = col%static%n
      call reset_column(col)
      call set_noon_forcing(col)
      col%forcing%abs_sw = 0.0_wp
      col%forcing%abs_par = 0.0_wp
      col%forcing%abs_sw_ground = 0.0_wp
      call set_air(col%forcing, T0, SHV_AIR)
      call set_ground_tiles(col)
      col%dynamic%wood_water_mass(1:n) = 0.8_wp * col%dynamic%wood_water_mass(1:n)
      call frozen_record(col, fz, ys)
      refill = sum((fz%plant%uptake_frozen(1:n) - fz%plant%sapflow_frozen(1:n)) * col%static%nplant(1:n))
      call check_true('refill: dried wood over moist soil refills', refill > 0.0_wp, refill)
      call check_true('refill: it brings the wood no heat (soil, wood and leaf at one temperature)',                &
                      maxval(abs(fz%tissue%q_wood_net(1:n))) < 1.0e-9_wp * refill * internal_energy_liquid(T0),      &
                      maxval(abs(fz%tissue%q_wood_net(1:n))))
      gained = sum(fz%tissue%q_wood_net(1:n)) + sum(fz%tissue%qwflux_wl(1:n))                                      &
               + sum(fz%tissue%leaf_water_enth(1:n)) + sum(fz%tissue%wood_water_enth(1:n))
      given  = sum(fz%roots%qloss_frozen(1:n))
      call check('refill: what the soil gives up is what the tissues and their water gain', gained, given,          &
                 1.0e-12_wp * abs(given))
   end subroutine test_refilling_wood

   !----- The wood's heat capacity is its dry wood at the dry specific heat plus the water network's   !
   !      wood store, its fine roots' share included (docs/science/vegetation_energy_dynamics.md eq.   !
   !      (8)), so a store drained by a fifth holds a fifth less of that water's heat (build_frozen). -!
   subroutine test_wood_capacity()
      type(frozen_t)      :: fz
      type(march_state_t) :: ys
      real(wp) :: w_store, cap_dry, cap_full, cap_drained
      print '(a)', 'test_wood_capacity:'
      call reset_column(col)
      call set_noon_forcing(col)
      w_store = col%dynamic%wood_water_mass(1)
      cap_dry = col%static%bwood(1) * col%static%nplant(1) * 2.0_wp * tissue_c_wood_dry
      call frozen_record(col, fz, ys)
      cap_full = fz%tissue%wood_hcap_per_dt(1) * DT_FAST
      call check('wood capacity: dry wood plus the wood store''s water', cap_full,                                 &
                 cap_dry + col%static%nplant(1) * w_store * cp_liq, 1.0e-12_wp * cap_full)
      call check_true('wood capacity: the fine roots hold part of the store', col%static%broot(1) > 0.0_wp,          &
                      col%static%broot(1))
      col%dynamic%wood_water_mass(1) = 0.8_wp * w_store
      call frozen_record(col, fz, ys)
      cap_drained = fz%tissue%wood_hcap_per_dt(1) * DT_FAST
      call check('wood capacity: a store drained by a fifth holds a fifth less water heat', cap_full - cap_drained,  &
                 0.2_wp * col%static%nplant(1) * w_store * cp_liq, 1.0e-9_wp * cap_full)
   end subroutine test_wood_capacity

   !----- The canopy films stay out of the tissues' heat capacity for the step (build_frozen; docs/     !
   !      science/vegetation_energy_dynamics.md, "The canopy's interception film"). A film is part of   !
   !      its tissue's energy, but within a step the tissue's balance runs on the tissue's own heat     !
   !      capacity, its dry matter and its internal water, and the commit reads the temperature back    !
   !      with the end-of-step film. Counted in the step's capacity, a film far heavier than its tissue !
   !      that evaporates within the step would leave its heat to the tissue alone: at Palo Verde a     !
   !      nearly leafless cohort's leaf read back at 107 K. The column frozen dry and frozen with its   !
   !      leaf and wood films half full gives the same two capacities, each the heat its tissue's dry   !
   !      matter and internal water gain per kelvin; each film holds enough heat that counting it would !
   !      show. ----------------------------------------------------------------------------------------!
   subroutine test_film_outside_capacity()
      real(wp), parameter :: FILL = 0.5_wp   ! the films' share of their capacity
      type(frozen_t)      :: fz
      type(march_state_t) :: ys
      real(wp) :: leaf_cap_no_film, wood_cap_no_film, leaf_film, wood_film, leaf_own, wood_own
      print '(a)', 'test_film_outside_capacity:'
      col%cfg%canopy_water_on = .true.
      call reset_column(col)
      call set_noon_forcing(col)
      call frozen_record(col, fz, ys)
      leaf_cap_no_film = fz%tissue%leaf_hcap_per_dt(1)
      wood_cap_no_film = fz%tissue%wood_hcap_per_dt(1)
      leaf_film = FILL * col%cfg%pft%leaf_surf_water_max(col%static%pft(1)) * col%static%lai(1)
      wood_film = FILL * col%cfg%pft%wood_surf_water_max(col%static%pft(1)) * col%static%wai(1)
      call reset_column(col)
      col%dynamic%leaf_surf_water(1) = leaf_film
      col%dynamic%wood_surf_water(1) = wood_film
      call frozen_record(col, fz, ys)
      call check_true('film: the wet column''s frozen build sees its films', fz%film%f_wet_c(1) > 0.0_wp,          &
                      fz%film%f_wet_c(1))
      call check('film: the leaves'' capacity for the step leaves their film out', fz%tissue%leaf_hcap_per_dt(1),    &
                 leaf_cap_no_film, 1.0e-12_wp * leaf_cap_no_film)
      call check('film: the wood''s capacity for the step leaves its film out', fz%tissue%wood_hcap_per_dt(1),       &
                 wood_cap_no_film, 1.0e-12_wp * wood_cap_no_film)
      !----- Each is its tissue's own, the dry matter and the internal water, as the column made them. -!
      leaf_own = heat_per_kelvin(tissue_dry_heat_capacity(col%static%nplant(1), col%static%bleaf(1), tissue_c_leaf),  &
                                 col%dynamic%leaf_water_mass(1) * col%static%nplant(1))
      wood_own = heat_per_kelvin(tissue_dry_heat_capacity(col%static%nplant(1), col%static%bwood(1),                 &
                                                          tissue_c_wood_dry),                                       &
                                 col%dynamic%wood_water_mass(1) * col%static%nplant(1))
      call check('film: the leaves'' capacity is their dry matter''s and internal water''s',                       &
                 fz%tissue%leaf_hcap_per_dt(1) * DT_FAST, leaf_own, 1.0e-12_wp * leaf_own)
      call check('film: the wood''s capacity is its dry matter''s and internal water''s',                          &
                 fz%tissue%wood_hcap_per_dt(1) * DT_FAST, wood_own, 1.0e-12_wp * wood_own)
      !----- Counted, either film would move its capacity far beyond the tolerances above. ------------!
      call check_true('film: the leaf film holds heat a capacity counting it would show',                         &
                      heat_per_kelvin(0.0_wp, leaf_film) > 1.0e-6_wp * leaf_own,                                    &
                      heat_per_kelvin(0.0_wp, leaf_film) / leaf_own)
      call check_true('film: the wood film holds heat a capacity counting it would show',                         &
                      heat_per_kelvin(0.0_wp, wood_film) > 1.0e-6_wp * wood_own,                                    &
                      heat_per_kelvin(0.0_wp, wood_film) / wood_own)
      col%cfg%canopy_water_on = .false.
   end subroutine test_film_outside_capacity

   !----- A new cohort's predawn potential is seeded from the top soil layer's, which the leaf's        !
   !      stomatal limb reads until the cohort has a day of its own. On a Campbell soil the seed comes  !
   !      from the Campbell curve: the column's one parameter pair is (psi_sat, b) there, and read as   !
   !      van Genuchten's it gives a wrong potential, or none. The soil is made through the site's      !
   !      surface, as a run makes it. -----------------------------------------------------------!
   subroutine test_campbell_seed()
      real(wp), parameter :: THETA_SAT = 0.43_wp, THETA_RES = 0.0_wp, PSI_SAT = -0.26_wp, B = 5.65_wp
      type(soil_params_t) :: soil0
      real(wp) :: seed, predawn
      print '(a)', 'test_campbell_seed:'
      soil0 = col%cfg%soil
      col%cfg%soil%retention   = SOIL_RETENTION_CAMPBELL
      col%cfg%soil%theta_sat   = THETA_SAT
      col%cfg%soil%theta_res   = THETA_RES
      col%cfg%soil%curve_par_a = PSI_SAT
      col%cfg%soil%curve_par_n = B
      call make_site_surface(col%site, col%cfg)
      call use_scheme(INTEG_RK45)
      call reset_column(col)
      call set_noon_forcing(col)
      call ensure_output_capacity(col%out, col%static%n, .true.)   ! the per-cohort capture, for the seed
      call step_column(col, DT_FAST)
      predawn = col%out%cdiag(CO_PSI_PREDAWN, 1)
      call ensure_output_capacity(col%out, col%static%n, .false.)
      col%cfg%soil = soil0
      call make_site_surface(col%site, col%cfg)
      call use_scheme(INTEG_ESDIRK2)
      seed = grav_head * soil_theta_to_head(SOIL_RETENTION_CAMPBELL, THETA0, THETA_SAT, THETA_RES, PSI_SAT, B)
      call check_true('Campbell seed: the cohort has no predawn potential of its own yet',                          &
                      col%static%dmax_psi_leaf(1) > 0.0_wp, col%static%dmax_psi_leaf(1))
      call check_true('Campbell seed: the top layer is unsaturated (a seed below zero)', seed < 0.0_wp, seed)
      call check_true('Campbell seed: the predawn seed is the Campbell curve''s top-layer potential',               &
                      abs(predawn - seed) <= 1.0e-12_wp * abs(seed), predawn)
   end subroutine test_campbell_seed

   !=======================================================================================!
   ! 8. RK45 AT ITS EDGES.                                                                   !
   !=======================================================================================!

   !----- A drydown toward theta_res, the domain's other edge: two days from just above it (0.085      !
   !      against 0.078) under air at half the humidity, so ground evaporation and transpiration both  !
   !      pull hard. The water network solves each dry layer in effective saturation, which cannot go  !
   !      below zero, and both sinks fall off before the floor (the root sink with the soil-root       !
   !      conductance, evaporation with the dry surface layer), so theta never commits below          !
   !      theta_res and nothing corrects it after the step; both whole-column ledgers close. ---------!
   subroutine test_rk45_drydown()
      type(day_t) :: d
      real(wp)    :: res_min
      print '(a)', 'test_rk45_drydown:'
      call use_scheme(INTEG_RK45)
      d = day_t(theta_seed = 0.085_wp, hours = 48_ik, shv_air = 0.5_wp * SHV_AIR)
      call march_day(d)
      res_min = minval(col%site%soil_profile%theta_res(1:NSL))
      print '(a,es10.3,a,es10.3)', '   theta_min = ', d%theta_min, ', theta_res = ', res_min
      call check_true('RK45 drydown: theta never commits below theta_res', d%theta_min >= res_min - 1.0e-12_wp,      &
                      d%theta_min)
      call check_true('RK45 drydown: whole-column water closes', col%out%report%ledgers%whole_water%worst < 1.0e-6_wp, &
                      col%out%report%ledgers%whole_water%worst)
      call check_true('RK45 drydown: whole-column energy closes', col%out%report%ledgers%whole_energy%worst < 1.0e-3_wp, &
                      col%out%report%ledgers%whole_energy%worst)
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_rk45_drydown

   !----- A just-recruited cohort (leaf and wood area, density and carbon at the magnitude of a 30-year !
   !      Ithaca run's first recruitment) over 00:00-12:00, guarding two failures a mature canopy      !
   !      never meets: the tissue temperature overflowing through a near-zero balance denominator,     !
   !      which cannot happen because the forcing, the coupling and the heat capacity all scale with    !
   !      the tissue's area (docs/science/vegetation_energy_dynamics.md eqs. (1), (8)); and the canopy  !
   !      air or the soil leaving its domain within one explicit stage over a near-bare patch, which    !
   !      the stage clamp and the frozen ground evaporation keep evaluable. The fixture hands every     !
   !      cohort a full canopy's load; a 0.0027-LAI seedling intercepts ~0.1 % of it, so the radiation  !
   !      is scaled by Beer's law, as the two-stream would. -----------------------------------------!
   subroutine test_rk45_tiny_cohort()
      real(wp)    :: rad_frac
      integer(ik) :: istep
      logical     :: physical
      print '(a)', 'test_rk45_tiny_cohort:'
      associate (s => col%static)
         s%lai(1) = 0.00265_wp
         s%wai(1) = 0.000529_wp
         s%leaf_area(1) = 0.00265_wp
         s%bleaf(1) = 0.0203507_wp
         s%bsap(1)  = 0.00212943_wp
         s%broot(1) = 0.0203507_wp
         s%sap_area(1) = 1.0e-4_wp
         s%nplant(1) = 0.01_wp
         call derive_cohort_terms(s, col%cfg%pft%hydraulics, col%site%soil_profile)
      end associate
      call use_scheme(INTEG_RK45)
      call reset_column(col)
      physical = .true.
      rad_frac = 1.0_wp - exp(-0.5_wp * col%static%lai(1))
      do istep = 1_ik, 12_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(col, (real(istep, wp) - 0.5_wp) * DT_FAST)
         col%forcing%abs_sw  = rad_frac * col%forcing%abs_sw
         col%forcing%abs_par = col%forcing%abs_sw
         col%forcing%abs_lw  = rad_frac * col%forcing%abs_lw
         call step_column(col, DT_FAST)
         associate (dy => col%dynamic)
            physical = physical .and. dy%leaf_temp(1) > 250.0_wp .and. dy%leaf_temp(1) < 340.0_wp                   &
                       .and. dy%wood_temp(1) > 250.0_wp .and. dy%wood_temp(1) < 340.0_wp                            &
                       .and. dy%cas%can_temp > 250.0_wp .and. dy%cas%can_temp < 340.0_wp                            &
                       .and. all(dy%soil_water%theta(1:NSL) >= 0.0_wp) .and. all(dy%soil_water%theta(1:NSL) <= 0.6_wp) &
                       .and. all(dy%soil_energy%soil_temp(1:NSL) > 250.0_wp)                                       &
                       .and. all(dy%soil_energy%soil_temp(1:NSL) < 340.0_wp)                                       &
                       .and. col%out%gpp_coh(1) == col%out%gpp_coh(1)
         end associate
      end do
      call check_true('RK45 tiny cohort: stays physical and bounded (00:00-12:00)', physical, col%dynamic%leaf_temp(1))
      call check_true('RK45 tiny cohort: its GPP is a number', col%out%gpp_coh(1) == col%out%gpp_coh(1), col%out%gpp_coh(1))
      call gather_column(col)
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_rk45_tiny_cohort

   !----- A very dense canopy (LAI 8) in a hard cold snap at twice dt_fast: the operating point that  !
   !      rails the explicit RK45 surface at scale. Whether the rescue onto ESDIRK2 fires here or the   !
   !      stage clamps alone hold it, the committed state stays physical (canopy air, leaf and soil in   !
   !      250-340 K, GPP a number). A 30-year Ithaca run collapses at LAI ~1 in its winter of 2050      !
   !      without the rescue, a point no short march reproduces, so the rescue count is printed, not    !
   !      asserted. ---------------------------------------------------------------------------!
   subroutine test_rk45_dense_cold()
      integer(ik) :: istep, rescues
      logical     :: physical
      print '(a)', 'test_rk45_dense_cold:'
      col%static%lai(1) = 8.0_wp
      col%static%wai(1) = 1.5_wp
      col%static%leaf_area(1) = 26.0_wp
      call use_scheme(INTEG_RK45)
      call reset_column(col)
      physical = .true.
      rescues = 0_ik
      do istep = 1_ik, STEPS_PER_DAY
         call set_coldsnap_forcing(col, (real(istep, wp) - 0.5_wp) * DT_FAST)
         call step_column(col, 2.0_wp * DT_FAST)
         rescues = rescues + col%out%report%work(WORK_RK45_RESCUE)
         associate (dy => col%dynamic)
            physical = physical .and. dy%cas%can_temp > 250.0_wp .and. dy%cas%can_temp < 340.0_wp                   &
                       .and. dy%leaf_temp(1) > 250.0_wp .and. dy%leaf_temp(1) < 340.0_wp                            &
                       .and. all(dy%soil_energy%soil_temp(1:NSL) > 250.0_wp)                                       &
                       .and. all(dy%soil_energy%soil_temp(1:NSL) < 340.0_wp)                                       &
                       .and. col%out%gpp_coh(1) == col%out%gpp_coh(1)
         end associate
      end do
      print '(a,i0,a)', '   dense cold snap: ', rescues, ' RK45 -> ESDIRK2 rescues over the run'
      call check_true('RK45 dense cold snap: the canopy stays physical (surface stability)', physical,               &
                      col%dynamic%cas%can_temp)
      call gather_column(col)
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_rk45_dense_cold

   !----- On moderately dry soil (theta 0.14, psi ~ -10 m) the wood-soil interface cancels to machine  !
   !      precision: the soil loses exactly the uptake the wood gains, asked at the ledger's round-off  !
   !      rather than at its relative tolerance. The drydown sits near theta_res and the other runs in  !
   !      wet soil; this one sits between, on the RK45 path itself (no rescue). ---------------------!
   subroutine test_rk45_dry_seam()
      type(day_t) :: d
      real(wp)    :: psi_top
      print '(a)', 'test_rk45_dry_seam:'
      call use_scheme(INTEG_RK45)
      d = day_t(theta_seed = 0.14_wp, hours = 12_ik)
      call march_day(d)
      associate (sp => col%site%soil_profile)
         psi_top = soil_theta_to_head(sp%retention, col%dynamic%soil_water%theta(1), sp%theta_sat(1), sp%theta_res(1), &
                                      sp%curve_par_a(1), sp%curve_par_n(1))
      end associate
      call check_true('RK45 dry seam: the soil really is moderately dry (-150 m < psi < -3 m)',                      &
                      psi_top < -3.0_wp .and. psi_top > -150.0_wp, psi_top)
      call check_true('RK45 dry seam: the RK45 path itself ran (no rescue)',                                         &
                      col%out%report%work(WORK_RK45_RESCUE) == 0_ik, real(col%out%report%work(WORK_RK45_RESCUE), wp))
      call check_true('RK45 dry seam: whole-column water closes to round-off (soil debit = wood credit)',            &
                      col%out%report%ledgers%whole_water%worst < 1.0e-8_wp, col%out%report%ledgers%whole_water%worst)
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_rk45_dry_seam

   !=======================================================================================!
   ! 9. THE REPORTED NET RADIATION. The step reports the net radiation the ground, the leaves  !
   ! and the wood absorbed, as their balances booked it (out%rnet, the output's net radiation; !
   ! docs/science/vegetation_energy_dynamics.md, "The emission base"). One step from the fresh !
   ! column under the noon air, its leaves and ground given no radiation, twice: the wood      !
   ! absorbing nothing, then 40 W/m2 of shortwave. The second reports 40 W/m2 more, less only  !
   ! what the wood, warmed over the step, emits beyond the first (well under 2 %). Leaving the  !
   ! wood out of the sum would report no change at all.                                       !
   !=======================================================================================!
   subroutine test_rnet_counts_wood()
      real(wp), parameter :: SW_WOOD = 40.0_wp, AGREE = 0.02_wp   ! [W/m2], [-]
      real(wp)    :: rnet(2)
      integer(ik) :: is, k
      print '(a)', 'test_rnet_counts_wood:'
      do is = 1_ik, size(SCHEMES, kind=ik)
         call use_scheme(SCHEMES(is))
         do k = 1_ik, 2_ik
            call reset_column(col)
            call set_noon_forcing(col)
            col%forcing%abs_sw        = 0.0_wp
            col%forcing%abs_par       = 0.0_wp
            col%forcing%abs_lw        = 0.0_wp
            col%forcing%abs_sw_ground = 0.0_wp
            col%forcing%abs_lw_ground = 0.0_wp
            col%forcing%abs_sw_wood   = merge(SW_WOOD, 0.0_wp, k == 2_ik)
            col%forcing%abs_lw_wood   = 0.0_wp
            call set_ground_tiles(col)
            call step_column(col, DT_FAST)
            rnet(k) = col%out%rnet
         end do
         print '(3a,2f10.4,a)', '   ', scheme_name(SCHEMES(is)), ': net radiation, the wood dark and lit ', rnet, ' W/m2'
         call check_true(trim(scheme_name(SCHEMES(is)))//': the reported net radiation counts what the wood absorbs', &
                         rnet(2) - rnet(1) <= SW_WOOD * (1.0_wp + 1.0e-12_wp)                                       &
                         .and. rnet(2) - rnet(1) >= (1.0_wp - AGREE) * SW_WOOD, rnet(2) - rnet(1))
      end do
      col%forcing%abs_sw_wood = 0.0_wp
      call use_scheme(INTEG_ESDIRK2)
   end subroutine test_rnet_counts_wood

end program test_fast_integ_step
