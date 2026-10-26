! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_test_fast_fixture -- the column the fast integrator's tests run on, assembled the  !
! way a run assembles one, and the helpers they share.                                     !
!   * make_fast_column: the in-code config (build_test_config) with the test soil, a 2 m   !
!     loam in 10 geometric layers; the site's surface made from it (make_site_surface);    !
!     one 20 cm tree at 224 stems/ha born and gathered by the production path              !
!     (static_state_fixture, LAI ~3; gather_column); the step's output and forcing sized.  !
!   * reset_column: a fresh dynamic state (alloc_dynamic_state) at 288 K, the tissues      !
!     holding water at a mild tension and their energy at that temperature, the soil at a  !
!     given water content; no shed water; the site's soil and the step report as made.     !
!   * The forcing: set_diurnal_forcing (21 June at 40 N, the air following the sun),       !
!     set_noon_forcing (a constant noon), set_coldsnap_forcing (a -20 C day under a weak    !
!     sun); set_air and set_ground_tiles keep the air and the ground's radiation           !
!     consistent after a test changes them.                                                !
!   * step_column: one fast_step of the column.                                            !
!   * frozen_record: the frozen record a step marches against, in the stepper's order: the !
!     frozen build, the water network's step at the run's tolerances, and what it moved.   !
!   * all_ledgers_closed: the seven conservation ledgers closed on every step.             !
!   * kernel_record: the one hand-set record, for the kernel and scheme tests. Its soil is !
!     the fixture's, but its coefficients are set by hand: zero tissue heat capacity, the   !
!     canopy-air conductances held, no water flowing, so a closed form or an observed order !
!     can be read from the right-hand side.                                                !
!==========================================================================================!
module meds_test_fast_fixture
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : rho_h2o, cp_air, stefan_boltzmann, tissue_c_leaf, tissue_c_wood_dry
   use meds_config,            only : meds_config_t
   use meds_time,              only : meds_time_t, solar_cosz
   use meds_therm_lib,         only : cas_temp_to_enthalpy, temp_to_internal_energy, tissue_temp_to_energy,       &
                                      tissue_dry_heat_capacity
   use meds_soil_params,       only : soil_params_t
   use meds_soil_hydr_lib,     only : SOIL_RETENTION_VG
   use meds_site_state_vars,   only : site_surface_t, make_site_surface
   use meds_plant_hydr_lib,    only : tissue_psi_to_water
   use meds_pft_params,        only : hydraulics_table_t
   use meds_fast_state,        only : static_state_fixture
   use meds_fast_integ_types,  only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,   &
                                      fast_workspace_t, step_report_t, march_state_t, frozen_t, error_control_t,  &
                                      alloc_dynamic_state, ensure_forcing_capacity, ensure_output_capacity
   use meds_fast_integ_control, only : build_error_control
   use meds_fast_integ_frozen, only : build_frozen, record_water_flows
   use meds_fast_integ_backward_euler, only : water_step
   use meds_fast_integ_stepper, only : fast_step
   use meds_test_soil_grid,    only : geometric_layer_bottoms
   use meds_test_support,      only : build_test_config
   implicit none
   private

   public :: fast_column_t
   public :: make_fast_column, gather_column, reset_column, step_column
   public :: set_diurnal_forcing, set_noon_forcing, set_coldsnap_forcing, set_air, set_ground_tiles
   public :: frozen_record, all_ledgers_closed, kernel_record
   public :: N_COH, NSL, DT_FAST, STEPS_PER_HOUR, STEPS_PER_DAY, T0, THETA0, LAT, SIM_DATE, SHV_AIR

   !----- The column: one cohort over a 10-layer soil, stepped at 150 s, the step its checks were  !
   !      computed at (the shipped dt_fast is 900 s; a check that needs it passes it itself).      !
   !      Durations and times of day are written in hours and turned into steps here, so a change  !
   !      of DT_FAST keeps what each test covers: step i of a day runs at (i - 1/2)*DT_FAST after  !
   !      midnight, and the sun follows the clock. ------------------------------------------------!
   integer(ik), parameter :: N_COH = 1_ik
   integer(ik), parameter :: NSL   = 10_ik
   real(wp),    parameter :: DT_FAST = 150.0_wp           !< [s]
   integer(ik), parameter :: STEPS_PER_HOUR = nint(3600.0_wp / DT_FAST, ik)
   integer(ik), parameter :: STEPS_PER_DAY  = 24_ik * STEPS_PER_HOUR
   real(wp),    parameter :: T0     = 288.0_wp            !< [K] the column's starting temperature
   real(wp),    parameter :: THETA0 = 0.30_wp             !< [m3/m3] the soil's starting water content
   real(wp),    parameter :: LAT    = 40.0_wp             !< [deg N]
   real(wp),    parameter :: SHV_AIR = 0.008_wp           !< [kg/kg] the air's humidity
   type(meds_time_t), parameter :: SIM_DATE = meds_time_t(2001_ik, 6_ik, 21_ik)
   real(wp),    parameter :: PSI_START = -0.1_wp          !< [MPa] the tissues' starting potential, a mild tension

   !----- Everything one fast step takes, as a run holds it for a patch. --------------------------!
   type :: fast_column_t
      type(meds_config_t)        :: cfg
      type(site_surface_t)       :: site      !< the column's soil and its boundaries
      type(fast_static_state_t)  :: static    !< the cohort, the canopy geometry and its top
      type(fast_dynamic_state_t) :: dynamic
      type(fast_forcing_t)       :: forcing
      type(fast_output_t)        :: out
      type(fast_workspace_t)     :: ws        !< the fast step's storage, reused step to step
   end type fast_column_t

contains

   !----- The column, assembled through production: the config and its [soil], the site's surface, !
   !      the gathered cohort, and the step's output and forcing sized for it. ----------------------!
   subroutine make_fast_column(col)
      type(fast_column_t), intent(out) :: col
      col%cfg = build_test_config()
      col%cfg%dt_fast         = DT_FAST          ! the column's step, not the fixture's shipped 900 s
      col%cfg%n_fast_per_slow = STEPS_PER_DAY
      !----- The test soil: the loam these checks were computed on, in geometric layers that thicken !
      !      with depth (meds_test_soil_grid), given as an explicit list of bottoms. ----------------!
      col%cfg%soil = soil_params_t()
      col%cfg%soil%depth       = 2.0_wp
      col%cfg%soil%retention   = SOIL_RETENTION_VG
      col%cfg%soil%theta_sat   = 0.43_wp
      col%cfg%soil%theta_res   = 0.078_wp
      col%cfg%soil%ksat        = 2.89e-6_wp
      col%cfg%soil%curve_par_a = 3.6_wp
      col%cfg%soil%curve_par_n = 1.56_wp
      col%cfg%soil%psi_fc      = -3.37_wp
      col%cfg%soil%solid_conductivity = 3.0_wp
      col%cfg%soil%dry_conductivity   = 0.15_wp
      col%cfg%soil%dry_heat_capacity  = 2.0e6_wp
      col%cfg%soil%n_layer_bottoms = NSL
      call geometric_layer_bottoms(col%cfg%soil%depth, NSL, 3.0_wp, col%cfg%soil%layer_bottoms)
      call make_site_surface(col%site, col%cfg)
      call gather_column(col)
      call ensure_output_capacity(col%out, N_COH, .false.)
      call ensure_forcing_capacity(col%forcing, N_COH)
      col%forcing%abs_sw_wood = 0.0_wp
      col%forcing%abs_lw_wood = 0.0_wp
      call reset_column(col)
   end subroutine make_fast_column

   !----- The cohort, born through the canonical path (init_cohort + set_cohort_size) and gathered by !
   !      the production gather, so its height, leaf area and the wood carbon its heat capacity is    !
   !      sized from agree; then the canopy geometry and the forcing's reference height. A test that  !
   !      changes the PFT table or the cohort calls it again to gather the column afresh. ----------!
   subroutine gather_column(col)
      type(fast_column_t), intent(inout) :: col
      call static_state_fixture(col%static, col%cfg%pft, [1_ik], [20.0_wp], [0.0224_wp], col%site%soil_profile)
      col%static%geom%veg_height   = 18.0_wp
      col%static%geom%opencan_frac = 0.0_wp
      col%static%geom%snowfac      = 0.0_wp
      col%static%z_top = 30.0_wp
   end subroutine gather_column

   !----- A fresh column state: the canopy air and tissues at T0, the tissue water at PSI_START     !
   !      (alloc_dynamic_state leaves a scratch zero, which the gather overwrites in a run), each    !
   !      tissue's energy at T0 with that water, the soil at theta_seed and T0, no shed water, the   !
   !      site's soil as made, and the step report cleared. ----------------------------------------!
   subroutine reset_column(col, theta_seed)
      type(fast_column_t), intent(inout) :: col
      real(wp), optional,  intent(in)    :: theta_seed   !< [m3/m3] default THETA0
      real(wp)    :: theta
      integer(ik) :: i, k
      theta = THETA0
      if (present(theta_seed)) theta = theta_seed
      call alloc_dynamic_state(col%dynamic, col%static%n, T0, SHV_AIR, 400.0_wp, T0)
      associate (s => col%static, d => col%dynamic, hyd => col%cfg%pft%hydraulics)
         do i = 1_ik, s%n
            d%leaf_water_mass(i) = tissue_psi_to_water(PSI_START, hyd%pft(s%pft(i))%leaf_curve, s%bleaf(i))
            d%wood_water_mass(i) = tissue_psi_to_water(PSI_START, hyd%pft(s%pft(i))%wood_curve, s%bsap(i) + s%broot(i))
            d%leaf_energy(i) = tissue_temp_to_energy(tissue_dry_heat_capacity(s%nplant(i), s%bleaf(i), tissue_c_leaf), &
                                                     d%leaf_water_mass(i) * s%nplant(i), T0)
            d%wood_energy(i) = tissue_temp_to_energy(tissue_dry_heat_capacity(s%nplant(i), s%bwood(i),               &
                                                                              tissue_c_wood_dry),                   &
                                                     d%wood_water_mass(i) * s%nplant(i), T0)
         end do
         d%soil_water%theta(1:NSL) = theta
         do k = 1_ik, NSL
            d%soil_energy%soil_energy(k) = temp_to_internal_energy(col%site%soil_profile%dry_heat_capacity(k),      &
                                                                   theta * rho_h2o, T0, 1.0_wp)
            d%soil_energy%soil_temp(k) = T0
         end do
      end associate
      col%static%shed_water_rate = 0.0_wp
      col%static%shed_water_enth = 0.0_wp
      col%site%soil = col%cfg%soil
      col%out%report = step_report_t()
   end subroutine reset_column

   !----- One fast step of the column. -----------------------------------------------------------!
   subroutine step_column(col, dt)
      type(fast_column_t), intent(inout) :: col
      real(wp),            intent(in)    :: dt   !< [s]
      call fast_step(dt, col%cfg, col%site, col%static, col%forcing, col%dynamic, col%out, col%ws)
   end subroutine step_column

   !----- 21 June at 40 N, t_sec seconds after midnight: 500 W/m2 of shortwave on the leaves and 75 !
   !      on the ground at the zenith, scaled by cos(zenith); the air 288 K +- 6 K with the sun, at  !
   !      SHV_AIR and 400 ppm; no longwave, rain or snow. -------------------------------------------!
   subroutine set_diurnal_forcing(col, t_sec)
      type(fast_column_t), intent(inout) :: col
      real(wp),            intent(in)    :: t_sec   !< [s] after midnight
      real(wp) :: cosz
      cosz = solar_cosz(SIM_DATE, t_sec, LAT)
      call set_light(col%forcing, 500.0_wp * cosz, 75.0_wp * cosz)
      call set_air(col%forcing, 288.0_wp + 6.0_wp * (cosz - 0.3_wp), SHV_AIR)
      call set_ground_tiles(col)
   end subroutine set_diurnal_forcing

   !----- A constant noon: 450 W/m2 on the leaves, 70 on the ground, the air at 295 K. --------------!
   subroutine set_noon_forcing(col)
      type(fast_column_t), intent(inout) :: col
      call set_light(col%forcing, 450.0_wp, 70.0_wp)
      call set_air(col%forcing, 295.0_wp, SHV_AIR)
      call set_ground_tiles(col)
   end subroutine set_noon_forcing

   !----- A hard cold snap over the warm column reset_column seeds: very cold, very dry air (~248-256 K, !
   !      1 g/kg) under a weak winter sun, a large gradient through the canopy. ---------------------!
   subroutine set_coldsnap_forcing(col, t_sec)
      type(fast_column_t), intent(inout) :: col
      real(wp),            intent(in)    :: t_sec   !< [s] after midnight
      real(wp) :: cosz
      cosz = solar_cosz(SIM_DATE, t_sec, LAT)
      call set_light(col%forcing, 120.0_wp * cosz, 20.0_wp * cosz)
      call set_air(col%forcing, 252.0_wp + 4.0_wp * (cosz - 0.3_wp), 0.001_wp)
      call set_ground_tiles(col)
   end subroutine set_coldsnap_forcing

   !----- The radiation on the leaves (shortwave; PAR is taken equal to it) and on the ground, with no !
   !      longwave and no precipitation. ------------------------------------------------------------!
   subroutine set_light(forcing, sw_leaf, sw_ground)
      type(fast_forcing_t), intent(inout) :: forcing
      real(wp),             intent(in)    :: sw_leaf, sw_ground   !< [W/m2]
      forcing%abs_sw        = sw_leaf
      forcing%abs_par       = sw_leaf
      forcing%abs_lw        = 0.0_wp
      forcing%abs_sw_ground = sw_ground
      forcing%abs_lw_ground = 0.0_wp
      forcing%rainfall      = 0.0_wp
      forcing%snowfall      = 0.0_wp
   end subroutine set_light

   !----- The air at the canopy-air top: its temperature, which the surface layer measures the canopy  !
   !      against and which values frozen precipitation, its humidity, its enthalpy and 400 ppm. ------!
   subroutine set_air(forcing, t_air, shv)
      type(fast_forcing_t), intent(inout) :: forcing
      real(wp),             intent(in)    :: t_air   !< [K]
      real(wp),             intent(in)    :: shv     !< [kg/kg]
      forcing%air_temp     = t_air
      forcing%shv_atm      = shv
      forcing%enthalpy_atm = cas_temp_to_enthalpy(t_air, shv)
      forcing%co2_atm      = 400.0_wp
   end subroutine set_air

   !----- Each ground tile's radiation from the ground's: the bare soil nets what the ground absorbs, !
   !      and a pack absorbs its shortwave under a sky that emits what the pack does, so the pack has  !
   !      no longwave loss at the step's start. --------------------------------------------------!
   subroutine set_ground_tiles(col)
      type(fast_column_t), intent(inout) :: col
      col%forcing%rnet_soil      = col%forcing%abs_sw_ground + col%forcing%abs_lw_ground
      col%forcing%sw_snow        = col%forcing%abs_sw_ground
      col%forcing%lw_down_ground = stefan_boltzmann * col%dynamic%snow%snow_temp(1) ** 4
   end subroutine set_ground_tiles

   !----- The frozen record a step of dt (default DT_FAST) marches against and the march's start, in  !
   !      fast_step's order: the frozen build, the water network's step at the run's tolerances, and  !
   !      what it moved recorded into the record. ---------------------------------------------------!
   subroutine frozen_record(col, fz, ys, dt)
      type(fast_column_t), intent(inout) :: col
      type(frozen_t),      intent(out)   :: fz
      type(march_state_t), intent(out)   :: ys
      real(wp), optional,  intent(in)    :: dt   !< [s]
      type(error_control_t) :: ec
      real(wp) :: h
      h = DT_FAST
      if (present(dt)) h = dt
      ec = build_error_control(col%cfg)
      call build_frozen(h, col%cfg, col%site, col%static, col%forcing, col%dynamic, col%out, col%ws%water, fz, ys)
      call water_step(h, col%static, col%site%soil_profile, col%site%soil, col%cfg%pft%hydraulics, fz%water_drive,  &
                      ec%tols, col%ws%water, col%out%report%work)
      call record_water_flows(h, col%cfg, col%site, col%static, col%dynamic, col%ws%water, col%out, fz)
   end subroutine frozen_record

   !----- Did every step close all seven ledgers: the canopy air's energy, water and CO2, the soil's  !
   !      energy and water, and the whole column's energy and water? Those that did not are named. --!
   logical function all_ledgers_closed(report) result(ok)
      type(step_report_t), intent(in) :: report
      associate (l => report%ledgers)
         ok = l%cas_energy%n_fail == 0_ik .and. l%cas_water%n_fail == 0_ik .and. l%cas_co2%n_fail == 0_ik          &
              .and. l%soil_energy%n_fail == 0_ik .and. l%soil_water%n_fail == 0_ik                                &
              .and. l%whole_energy%n_fail == 0_ik .and. l%whole_water%n_fail == 0_ik
         if (l%cas_energy%n_fail > 0_ik)   print '(a,i0)', '     cas_energy fails: ', l%cas_energy%n_fail
         if (l%cas_water%n_fail > 0_ik)    print '(a,i0)', '     cas_water fails: ', l%cas_water%n_fail
         if (l%cas_co2%n_fail > 0_ik)      print '(a,i0)', '     cas_co2 fails: ', l%cas_co2%n_fail
         if (l%soil_energy%n_fail > 0_ik)  print '(a,i0)', '     soil_energy fails: ', l%soil_energy%n_fail
         if (l%soil_water%n_fail > 0_ik)   print '(a,i0)', '     soil_water fails: ', l%soil_water%n_fail
         if (l%whole_energy%n_fail > 0_ik) print '(a,i0)', '     whole_energy fails: ', l%whole_energy%n_fail
         if (l%whole_water%n_fail > 0_ik)  print '(a,i0)', '     whole_water fails: ', l%whole_water%n_fail
      end associate
   end function all_ledgers_closed

   !----- The kernel and scheme tests' record: n cohorts over the site's soil, a representative       !
   !      daytime state, every coefficient set by hand. The tissues hold no heat (zero capacity) and   !
   !      take their longwave about 298.5 K; the canopy-air box has a 20 m depth and its atmospheric   !
   !      conductances held (mo_live off), so the CO2 twin is a decoupled affine equation; the leaves  !
   !      transpire at the frozen sapflow, held as the water network holds it; no water flows through  !
   !      the soil (theta_rate 0), because these marches run for hours on one record, which a run     !
   !      never does. The tissue water starts at -1 MPa in the leaves and -0.5 MPa in the wood, on     !
   !      the PFT's curves. ------------------------------------------------------------------------!
   subroutine kernel_record(site, hyd, n, y, frozen)
      type(site_surface_t),     intent(in)  :: site
      type(hydraulics_table_t), intent(in)  :: hyd
      integer(ik),              intent(in)  :: n
      type(march_state_t),      intent(out) :: y
      type(frozen_t),           intent(out) :: frozen
      real(wp), parameter :: BOX_DEPTH = 20.0_wp      !< [m] the canopy air's depth
      real(wp), parameter :: G_EXCHANGE = 0.3_wp * 0.02_wp   !< [m/s] u* times the heat transfer factor
      real(wp) :: rho, molar_density, root(site%soil_profile%n_active)
      integer(ik) :: i, k, nsl_s
      nsl_s = site%soil_profile%n_active
      rho = 1.2_wp
      molar_density = rho * (1.0_wp - 0.012_wp) / 0.0289655_wp   ! [mol dry air/m3]
      !----- The plants: 0.3 per m2, each taking up and passing 1e-4 kg/s, roots e^(-2d) by layer. ------!
      allocate(frozen%plant%nplant(n), frozen%plant%sapflow_frozen(n), frozen%plant%uptake_frozen(n))
      frozen%plant%nplant = 0.3_wp
      frozen%plant%sapflow_frozen = 1.0e-4_wp
      frozen%plant%uptake_frozen  = 1.0e-4_wp
      do k = 1_ik, nsl_s
         root(k) = exp(2.0_wp * site%soil_profile%soil_layer_z(k)) - exp(2.0_wp * site%soil_profile%soil_layer_z(k+1))
      end do
      allocate(frozen%roots%root_share(nsl_s), frozen%roots%qloss_frozen(n), frozen%roots%root_heat(nsl_s))
      frozen%roots%root_share = root / sum(root)
      frozen%roots%uptake = frozen%plant%uptake_frozen(1) * sum(frozen%plant%nplant)
      frozen%roots%qloss_frozen = 0.0_wp
      frozen%roots%root_heat = 0.0_wp
      frozen%hydrology%q_top  = 1.0e-6_wp
      frozen%hydrology%t_infil = 296.0_wp
      frozen%hydrology%t_bot   = 292.0_wp
      !----- The tissues: leaf area 2.0, 1.5, ... from the top, 250, 200, ... W/m2 of shortwave, -30 of  !
      !      longwave; a series transpiration conductance of 4 mm/s per unit leaf area; no wood area. ---!
      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n),          &
               frozen%tissue%abs_lw(n), frozen%tissue%lai(n), frozen%tissue%h_coeff_w(n),                       &
               frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n),                &
               frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n), frozen%tissue%t_leaf0(n),  &
               frozen%tissue%t_wood0(n), frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n),               &
               frozen%tissue%leaf_emiss(n), frozen%tissue%wood_emiss(n), frozen%tissue%transp_net(n))
      do i = 1_ik, n
         frozen%tissue%lai(i)    = 2.0_wp - 0.5_wp * real(i - 1_ik, wp)
         frozen%tissue%abs_sw(i) = 250.0_wp - 50.0_wp * real(i - 1_ik, wp)
      end do
      frozen%tissue%abs_lw = -30.0_wp
      frozen%tissue%h_coeff_leaf  = 2.0_wp * frozen%tissue%lai * 0.03_wp * rho * cp_air
      frozen%tissue%g_transp_leaf = 0.004_wp * frozen%tissue%lai
      frozen%tissue%h_coeff_w   = 0.0_wp
      frozen%tissue%abs_sw_wood = 0.0_wp
      frozen%tissue%abs_lw_wood = 0.0_wp
      frozen%tissue%wai         = 0.0_wp
      frozen%tissue%leaf_hcap_per_dt = 0.0_wp
      frozen%tissue%wood_hcap_per_dt = 0.0_wp
      frozen%tissue%t_leaf0 = 298.5_wp
      frozen%tissue%t_wood0 = 298.5_wp
      frozen%tissue%qwflux_wl  = 0.0_wp
      frozen%tissue%q_wood_net = 0.0_wp
      frozen%tissue%leaf_emiss = 0.95_wp
      frozen%tissue%wood_emiss = 0.95_wp
      frozen%tissue%transp_net = frozen%plant%sapflow_frozen * frozen%plant%nplant
      frozen%tissue%transp_held = .true.
      !----- No films. -----------------------------------------------------------------------------!
      allocate(frozen%film%f_wet_c(n), frozen%film%g_film_leaf(n), frozen%film%g_film_w(n),                     &
               frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n), frozen%film%transfer_lw(n))
      frozen%film%f_wet_c = 0.0_wp
      frozen%film%g_film_leaf = 0.0_wp
      frozen%film%g_film_w = 0.0_wp
      frozen%film%intercept_leaf = 0.0_wp
      frozen%film%intercept_wood = 0.0_wp
      frozen%film%transfer_lw = 0.0_wp
      !----- The canopy-air box, venting to 300 K, 11 g/kg and 400 ppm air; a net CO2 uptake of 5. ------!
      frozen%cas%rho   = rho
      frozen%cas%press = 101325.0_wp
      frozen%cas%nee_biotic = -5.0_wp
      frozen%cas%box%air_mass_capacity  = rho * BOX_DEPTH
      frozen%cas%box%air_molar_capacity = molar_density * BOX_DEPTH
      frozen%cas%box%atm_conductance_enthalpy = rho * G_EXCHANGE
      frozen%cas%box%atm_conductance_vapor    = rho * G_EXCHANGE
      frozen%cas%box%atm_conductance_co2      = molar_density * G_EXCHANGE
      frozen%cas%box%atm_enthalpy          = cas_temp_to_enthalpy(300.0_wp, 0.011_wp)
      frozen%cas%box%atm_specific_humidity = 0.011_wp
      frozen%cas%box%atm_co2               = 400.0_wp
      frozen%ground%rnet_soil = 50.0_wp
      frozen%ground%ggnet     = 0.02_wp
      frozen%ground%soil_evap = 2.0e-5_wp
      !----- The state: canopy air at 297 K, 12 g/kg and 410 ppm; the soil at 0.30 and 296 K cooling by  !
      !      0.4 K a layer. ----------------------------------------------------------------------------!
      y%cas_enthalpy = cas_temp_to_enthalpy(297.0_wp, 0.012_wp)
      y%cas_shv      = 0.012_wp
      y%cas_co2      = 410.0_wp
      do k = 1_ik, nsl_s
         y%theta(k) = 0.30_wp
         y%soil_energy(k) = temp_to_internal_energy(site%soil_profile%dry_heat_capacity(k), 0.30_wp * rho_h2o,       &
                                                    296.0_wp - 0.4_wp * real(k - 1_ik, wp), 1.0_wp)
      end do
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n), y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water = 0.0_wp
      y%wood_surf_water = 0.0_wp
      do i = 1_ik, n
         y%leaf_water_mass(i) = tissue_psi_to_water(-1.0_wp, hyd%pft(1)%leaf_curve, 0.5_wp)   ! 0.5 kgC of leaf
         y%wood_water_mass(i) = tissue_psi_to_water(-0.5_wp, hyd%pft(1)%wood_curve, 7.0_wp)   ! 5 sapwood, 2 root
      end do
   end subroutine kernel_record

end module meds_test_fast_fixture
