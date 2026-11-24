! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_surface -- the ground's surface and its water layer (ed_surface;
! docs/science/snow_biophysics.md). The layer steps as the integrator steps it: its stage
! (water_layer_stage: what falls joins it and its free water leaves it), the ground surface's
! exchange at its state (ground_surface_exchange), then it takes the step's amounts
! (water_layer_commit); its energy ledger is what the surroundings give up by the ground
! surface's balance.
!   1. COVER (water_layer_cover_of): 0 without snow or for a film, monotone in [0, 1], and
!      Niu & Yang (2007)'s fraction for 1 cm of snow.
!   2. ACCUMULATION: a cold snap's pack grows below t_3ple, all ice; its mass and energy
!      close against what fell, what it exchanged and what drained.
!   3. MELT: the pack pins at t_3ple, its liquid rises and drains as free water, and it
!      melts out; both ledgers close through the plateau.
!   4. ISOTHERMAL: a pack at the temperature of its surroundings, under canopy air
!      saturated over ice, stays as it is with no heat into the soil; under canopy air
!      saturated over liquid it gains mass.
!   5. RAIN ON SNOW (water_layer_stage): warm rain on a cold pack adds its mass and exactly
!      its enthalpy, and no depth.
!   6. STIFF COUPLING: a thin pack under an open patch's bare-ground conductance relaxes
!      toward its balance without overshooting it, over a step and a day of steps; both
!      ledgers close.
!   7. TILE RADIATION (ground_tile_radiation, snow_blend_ground_optics): on half-covered
!      ground the snow and the soil each absorb at their own optics, and weighted by the
!      cover they sum to what the blended ground absorbs.
!   8. BARE GROUND (ground_surface_fluxes): no sensible flux without a temperature
!      difference and no latent flux without evaporation; warm ground heats the canopy
!      air, in proportion to the difference and at the canopy air's specific heat; the
!      ground surface's exchange takes that specific heat from its surroundings.
!   9. THE STAGE (water_layer_stage): light snowfall on frozen ground gathers as a film at the
!      ground's temperature until it is a layer; a dusting on warm ground melts, the ground
!      paying, and drains; a pond with a round-off share of ice is free water and covers nothing;
!      a cold pack freezes a pond within its pores; a warm pond is all free water.
!  10. THE POROUS LAYERS (their helpers in ed_patch_state_vars, and surface_node_balance): a layer
!      holds its theta_sat over its depth, and passes what it holds beyond to the water layer at its
!      temperature; each its own temperature, from its own energy; heat
!      crosses a layer's dried top at its dry conductivity, in series with its wet part. Each
!      layer's curve is its material's at its depth (the moss's at -10 kPa as Zhang et al. 2010
!      measured, the litter's near what drained litter holds), whether frozen or empty. Its
!      balance closes on its energy, a thin layer of little heat capacity relaxes over an hour
!      without passing its surroundings, a dry layer does not evaporate, and a wet one gives no
!      more than it holds; held at the network's rate it gives that rate, evaporation or dew, and
!      takes its root; the soil takes what it conducts down, less the soil's own evaporation's
!      latent heat, and with a pack on it the ground's energy closes on what the pack and the
!      layer keep.
!  11. INSULATION: a deep pack colder than the soil, in balance with the air, passes the soil
!      the heat its whole depth and the top soil node conduct in series, since its one
!      temperature is its surface's; twice as deep, it passes less.
!  12. THE STACK (surface_stack_balance): a thin melting pack over a thin litter over warm soil,
!      solved together: what the nodes take is what they absorb less what they give the air and
!      the soil; Newton's root is the Gauss-Seidel sweeps' root; a single node's is its own
!      closed form's.
!  13. THE FREEZING ONSET: a saturated 3.4 mm litter just above the triple point, evaporating at
!      the network's rate, between a thin wet pack over a ninth of the ground and wet moss, under
!      canopy air swept across the temperature at which the litter's own ends at the triple point:
!      the ground's latent and sensible heat change continuously, with no jump of the litter's
!      vapour times the latent heat of fusion (3.2 W/m2). The canopy air's balance over such a
!      jump has no root (a tundra's late-summer snowfall on wet litter).
! The energy datum the pack is held in is test_therm_lib's.
!---~---
program test_surface
   use ed_test_assert,        only : check, check_true, test_report
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : t_3ple, stefan_boltzmann, rho_ice, cp_organic, cp_ice, latent_heat_fusion
   use ed_therm_lib,          only : temp_to_internal_energy, internal_energy_to_temp, internal_energy_ice,      &
                                     internal_energy_liquid, sat_specific_humidity, cas_heat_capacity,     &
                                     enthalpy_vapor, sat_specific_humidity_temp_deriv
   use ed_patch_state_vars,   only : ground_store_t, water_layer_state_t, porous_layer_state_t, ground_store_present, &
                                     porous_layer_curve, porous_layers_diagnose, porous_layers_spill,             &
                                     porous_layer_energy_at
   use ed_ground_surface_lib, only : water_layer_bulk_density, porous_layer_capacity, porous_layer_resistance,  &
                                     porous_layer_temp
   use ed_ground_surface_params, only : ground_surface_layer_params_t, water_material_t, POROUS_LITTER, POROUS_MOSS, &
                                        N_POROUS
   use ed_soil_hydr_lib,      only : soil_layer_curve_t, soil_head_to_theta
   use ed_surface,            only : surface_env_t, porous_node_env_t, water_stage_out_t, surface_node_t,      &
                                     node_flux_t, surface_node_balance, surface_stack_balance, water_layer_stage,  &
                                     water_layer_commit,                                                        &
                                     water_layer_cover_of, water_layer_ice, ground_surface_exchange, water_layer_exposed, &
                                     store_exchange_t, STORE_WATER, N_SURFACE_STORE,                            &
                                     snow_blend_ground_optics, ground_tile_radiation, ground_surface_fluxes,    &
                                     conductance_in_series
   implicit none

   real(wp), parameter :: MIN_DEPTH = 1.0e-4_wp   ! [m] the shipped [ground_surface_layer].min_depth

   call test_snow_cover()
   call test_accumulate()
   call test_melt_event()
   call test_isothermal()
   call test_rain_on_snow()
   call test_stiff_coupling()
   call test_tile_radiation()
   call test_ground_balance()
   call test_water_layer_stage()
   call test_porous_layers()
   call test_porous_balance()
   call test_insulation()
   call test_stack()
   call test_stack_water()
   call test_freezing_onset()

   call test_report('test_surface')

contains

   !--- Cold-night forcing (accumulation): no SW, net LW loss, cold dry CAS, warmer soil below.
   subroutine cold_env(env)
      type(surface_env_t), intent(out) :: env
      env%sw_snow     = 0.0_wp
      env%lw_down     = stefan_boltzmann * 265.0_wp ** 4 - 40.0_wp
      env%can_temp    = 265.0_wp
      env%can_shv     = 1.0e-3_wp
      env%ggnet       = 0.01_wp
      env%rho_air     = 1.3_wp
      env%press       = 101325.0_wp
      env%t_soil_top  = 273.0_wp
      env%r_soil_top  = 0.05_wp / 1.5_wp
   end subroutine cold_env

   !--- A pack of swe [kg/m2], all ice at temp, at the default snow density.
   function cold_pack(swe, temp, p) result(snow)
      real(wp),            intent(in) :: swe, temp   ! [kg/m2], [K]
      type(water_material_t), intent(in) :: p
      type(water_layer_state_t) :: snow
      snow = water_layer_state_t()
      snow%water = swe
      snow%energy = temp_to_internal_energy(0.0_wp, swe, temp, 0.0_wp)
      snow%matrix_depth  = swe / p%fresh_density
   end function cold_pack

   !---~---
   !   The pack's stage with nothing falling (water_layer_stage at the default min_depth, on soil
   ! at t_soil): its free water out, counted into the mass and energy the pack gave up, and the
   ! heat a film took from the soil into what it gained.
   !---~---
   subroutine drain(snow, p, t_soil, mass_out, e_out, e_in)
      type(water_layer_state_t), intent(inout) :: snow
      type(water_material_t),    intent(in)    :: p
      real(wp),                  intent(in)    :: t_soil
      real(wp),                  intent(inout) :: mass_out, e_out, e_in
      type(water_stage_out_t) :: st
      call water_layer_stage(0.0_wp, 0.0_wp, 273.0_wp, 0.0_wp, 0.0_wp, 900.0_wp, p, MIN_DEPTH, t_soil, 1.0_wp,       &
                             snow, st)
      mass_out = mass_out + st%free%water
      e_out    = e_out    + st%free%energy
      e_in     = e_in     + st%film_heat
   end subroutine drain

   !---~---
   !   One step of the pack as the integrator takes it: the ground surface's exchange at the
   ! pack's state and the given cover (ground_surface_exchange, the call the derivatives make),
   ! then the pack takes the step's amounts (water_layer_commit, the update's). Returns the
   ! energy the surroundings gave the pack [J/m2], by the ground surface's own balance -- its
   ! net radiation less the sensible and latent heat to the canopy air and the heat into the
   ! soil -- the vapour the pack lost [kg/m2], and the heat flux into the soil top [W/m2].
   !---~---
   subroutine snow_step(snow, env, cover, dt, energy_gain, vapour_loss, g_top)
      type(water_layer_state_t),  intent(inout) :: snow          ! the pack
      type(surface_env_t), intent(in)    :: env           ! its surroundings; the pack and cover are set here
      real(wp),            intent(in)    :: cover, dt     ! [-], [s]
      real(wp),            intent(out)   :: energy_gain, vapour_loss, g_top   ! [J/m2], [kg/m2], [W/m2]
      type(surface_env_t) :: e
      real(wp) :: h_ground, le_ground, ground_rad, snow_energy_in, snow_vapour_out
      type(store_exchange_t) :: stores
      integer(ik) :: it_stack, fb_stack
      e = env
      call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)
      e%pack    = snow
      e%snowfac = cover
      e%dt      = dt
      call ground_surface_exchange(e, h_ground, le_ground, g_top, ground_rad, stores, it_stack, fb_stack)
      snow_energy_in  = stores%energy_in(STORE_WATER)
      snow_vapour_out = stores%vapour_out(STORE_WATER)
      call water_layer_commit(snow, dt * snow_energy_in, dt * snow_vapour_out, ground_store_t(), 0.0_wp,           &
                              ground_store_t(), water_layer_ice(snow), e%params)
      energy_gain = dt * (ground_rad - h_ground - le_ground - g_top)
      vapour_loss = dt * snow_vapour_out
   end subroutine snow_step

   subroutine test_snow_cover()
      real(wp), parameter :: GROUND_ROUGH = 0.01_wp   ! [m] the [aerodynamics] soil_rough default
      type(water_material_t) :: p
      real(wp) :: f0, f1, f2
      print '(a)', 'test_snow_cover:'
      f0 = water_layer_cover_of(cold_pack(0.0_wp, 265.0_wp, p), .true., p, GROUND_ROUGH)
      f1 = water_layer_cover_of(cold_pack(1.0_wp, 265.0_wp, p), .true., p, GROUND_ROUGH)
      f2 = water_layer_cover_of(cold_pack(50.0_wp, 265.0_wp, p), .true., p, GROUND_ROUGH)
      call check('snowfac(0)=0', f0, 0.0_wp, 1.0e-12_wp)
      call check_true('snowfac monotone 0<f1<f2<=1', f1 > 0.0_wp .and. f1 < f2 .and. f2 <= 1.0_wp, f2)
      call check_true('deep pack fully covered (~1)', f2 > 0.99_wp, f2)
      call check('a film covers nothing', water_layer_cover_of(cold_pack(1.0_wp, 265.0_wp, p), .false., p,          &
                 GROUND_ROUGH), 0.0_wp, 0.0_wp)
      !--- 1 cm of snow on ground of 1 cm roughness covers about a sixth of it (Niu & Yang 2007).
      f1 = water_layer_cover_of(cold_pack(2.5_wp, 265.0_wp, p), .true., p, GROUND_ROUGH)
      call check('1 cm of snow covers tanh(0.16) of the ground', f1, tanh(0.16_wp), 1.0e-12_wp)
   end subroutine test_snow_cover

   !--- Build a cold pack by accumulation; assert growth, sub-freezing, and closed ledgers.
   subroutine test_accumulate()
      type(water_layer_state_t)  :: snow
      type(surface_env_t) :: env
      type(water_material_t) :: p
      type(water_stage_out_t) :: st
      real(wp) :: dt, snowfall, mass_in, mass_out, e_in, e_out, gain, vapour, g_top
      integer(ik) :: k
      print '(a)', 'test_accumulate:'
      call cold_env(env)
      snow = water_layer_state_t()
      dt = 900.0_wp
      snowfall = 1.0e-3_wp                                   ! [kg/m2/s] ~ 3.6 mm SWE/hr
      mass_in  = 0.0_wp
      mass_out = 0.0_wp
      e_in     = 0.0_wp
      e_out    = 0.0_wp
      do k = 1_ik, 40_ik
         call water_layer_stage(0.0_wp, snowfall, 264.0_wp, 0.0_wp, 0.0_wp, dt, p, MIN_DEPTH, env%t_soil_top,       &
                                1.0_wp, snow, st)
         mass_in  = mass_in + snowfall * dt
         e_in     = e_in + st%inflow_enth + st%film_heat
         mass_out = mass_out + st%free%water
         e_out    = e_out    + st%free%energy
         call snow_step(snow, env, 1.0_wp, dt, gain, vapour, g_top)
         mass_out = mass_out + vapour
         e_in     = e_in     + gain
      end do
      call check_true('pack grew (swe > 30 kg/m2)', snow%water > 30.0_wp, snow%water)
      call check_true('snow sub-freezing (temp < t_3ple)', snow%temp < t_3ple, snow%temp)
      call check_true('fully frozen (fliq = 0)', snow%fliq < 1.0e-9_wp, snow%fliq)
      call check('MASS closure: swe = in - out', snow%water, mass_in - mass_out, 1.0e-9_wp)
      call check('ENERGY closure: energy = in - out', snow%energy, e_in - e_out, 1.0e-4_wp)
   end subroutine test_accumulate

   !--- Warm a cold pack: temp pins at t_3ple, fliq rises, melt drains, swe -> melt-out.
   subroutine test_melt_event()
      type(water_layer_state_t)  :: snow
      type(surface_env_t) :: env
      type(water_material_t) :: p
      real(wp) :: dt, swe0, mass_out, e_in, e_out, plateau_temp, tot_melt, gain, vapour, g_top, w_before
      integer(ik) :: k
      logical :: saw_plateau
      print '(a)', 'test_melt_event:'
      call cold_env(env)
      snow = cold_pack(20.0_wp, 270.0_wp, p)                 ! a 20 kg/m2 cold pack, all ice at 270 K
      swe0 = snow%water
      !--- Strong warm radiative and canopy-air forcing to drive melt.
      env%sw_snow    = 400.0_wp
      env%lw_down    = stefan_boltzmann * 270.0_wp ** 4 + 20.0_wp
      env%can_temp   = 280.0_wp
      env%can_shv    = 5.0e-3_wp
      env%t_soil_top = 275.0_wp
      dt = 900.0_wp
      mass_out = 0.0_wp
      e_out    = 0.0_wp
      tot_melt = 0.0_wp
      e_in = snow%energy                             ! initial store folded into the ledger baseline
      saw_plateau  = .false.
      plateau_temp = 0.0_wp
      do k = 1_ik, 120_ik
         call snow_step(snow, env, 1.0_wp, dt, gain, vapour, g_top)
         mass_out = mass_out + vapour
         e_in     = e_in     + gain
         if (ground_store_present(snow) .and. snow%fliq > 0.05_wp .and. snow%fliq < 0.95_wp) then
            saw_plateau  = .true.
            plateau_temp = snow%temp
         end if
         w_before = mass_out
         call drain(snow, p, env%t_soil_top, mass_out, e_out, e_in)
         tot_melt = tot_melt + (mass_out - w_before)
         if (.not. ground_store_present(snow)) exit
      end do
      call check_true('melted out (the store is empty)', .not. ground_store_present(snow),                     &
                      merge(1.0_wp, 0.0_wp, ground_store_present(snow)))
      call check_true('saw the melt plateau (0<fliq<1)', saw_plateau, merge(1.0_wp, 0.0_wp, saw_plateau))
      call check('plateau temperature = t_3ple', plateau_temp, t_3ple, 1.0e-6_wp)
      call check_true('meltwater drained to soil (> half of swe0)', tot_melt > 0.5_wp * swe0, tot_melt)
      call check('MASS closure: 0 = swe0 - out', 0.0_wp, swe0 - mass_out, 1.0e-7_wp)
      call check('ENERGY closure: 0 = in - out', 0.0_wp, e_in - e_out, 1.0e-2_wp)
   end subroutine test_melt_event

   !---~---
   !   A cold pack at the temperature of its surroundings: the sky emits what it does, the soil
   ! is at its temperature, and the canopy air is saturated over ICE (fliq = 0), as the pack's
   ! surface saturates (docs/science/snow_biophysics.md eq. (1)); so nothing flows and it
   ! stays as it is. Saturated over LIQUID instead, the canopy air holds more vapour than the
   ! frozen surface, 3-8 % more at 268 K, and the pack gains mass by deposition.
   !---~---
   subroutine test_isothermal()
      type(water_layer_state_t)  :: snow, probe
      type(surface_env_t) :: env
      type(water_material_t) :: p
      real(wp) :: swe0, tsnow0, gain, vapour, g_top, m_out, e_out, e_in
      integer(ik) :: k
      print '(a)', 'test_isothermal:'
      snow = cold_pack(15.0_wp, 268.0_wp, p)
      env%sw_snow     = 0.0_wp
      env%lw_down     = stefan_boltzmann * 268.0_wp ** 4      ! the sky emits what the pack does
      env%can_temp    = 268.0_wp
      env%ggnet       = 0.01_wp
      env%rho_air     = 1.3_wp
      env%press       = 101325.0_wp
      env%t_soil_top  = 268.0_wp
      env%r_soil_top  = 0.05_wp / 1.5_wp
      env%can_shv     = sat_specific_humidity(268.0_wp, env%press, 0.0_wp)
      swe0   = snow%water
      tsnow0 = 268.0_wp
      m_out  = 0.0_wp
      e_out  = 0.0_wp
      e_in   = 0.0_wp
      do k = 1_ik, 50_ik
         call snow_step(snow, env, 1.0_wp, 900.0_wp, gain, vapour, g_top)
         call drain(snow, p, env%t_soil_top, m_out, e_out, e_in)
      end do
      call check_true('no heat into the soil (no soil gradient)', abs(g_top) < 1.0e-6_wp, g_top)
      call check_true('swe steady (drift < 0.5 kg/m2)', abs(snow%water - swe0) < 0.5_wp, snow%water - swe0)
      call check_true('snow_temp steady (within 1 K)', abs(snow%temp - tsnow0) < 1.0_wp, snow%temp)
      call check_true('stayed frozen (no spurious melt)', snow%fliq < 1.0e-6_wp, snow%fliq)
      !--- The same frozen pack under canopy air saturated over liquid: deposition.
      probe = cold_pack(15.0_wp, 268.0_wp, p)
      env%can_shv = sat_specific_humidity(268.0_wp, env%press)
      call snow_step(probe, env, 1.0_wp, 900.0_wp, gain, vapour, g_top)
      call check_true('frozen pack under a liquid-saturated CAS gains mass (deposition)', vapour < 0.0_wp, vapour)
   end subroutine test_isothermal

   !---~---
   !   A deep pack at 255 K on soil at 270 K, the sky and the canopy air in balance with it, so
   ! only the soil exchanges with it. Over one short step its temperature hardly moves, and the
   ! soil loses (T_pack - T_soil) / (depth / k_snow + dz / k_soil): the whole pack and the top
   ! soil node in series. A pack twice as deep passes the soil less.
   !---~---
   subroutine test_insulation()
      type(water_layer_state_t)  :: snow
      type(surface_env_t) :: env
      type(water_material_t) :: p
      real(wp) :: gain, vapour, g_top, g_deep, g_expect
      real(wp), parameter :: T_PACK = 255.0_wp, T_SOIL = 270.0_wp, SWE = 100.0_wp   ! [K], [K], [kg/m2]
      print '(a)', 'test_insulation:'
      env%sw_snow     = 0.0_wp
      env%lw_down     = stefan_boltzmann * T_PACK ** 4
      env%can_temp    = T_PACK
      env%ggnet       = 0.01_wp
      env%rho_air     = 1.3_wp
      env%press       = 101325.0_wp
      env%can_shv     = sat_specific_humidity(T_PACK, env%press, 0.0_wp)
      env%t_soil_top  = T_SOIL
      env%r_soil_top  = 0.05_wp / 1.5_wp
      snow = cold_pack(SWE, T_PACK, p)
      g_expect = (T_PACK - T_SOIL) / (snow%matrix_depth / p%conductivity_at_fresh + env%r_soil_top)
      call snow_step(snow, env, 1.0_wp, 60.0_wp, gain, vapour, g_top)
      call check('the soil loses what the whole pack and its top node conduct', g_top, g_expect,           &
                 1.0e-3_wp * abs(g_expect))
      snow = cold_pack(2.0_wp * SWE, T_PACK, p)
      call snow_step(snow, env, 1.0_wp, 60.0_wp, gain, vapour, g_deep)
      call check_true('a pack twice as deep passes the soil less', abs(g_deep) < abs(g_top), g_deep)
   end subroutine test_insulation

   !---~---
   !   Rain on a sub-freezing pack: it adds its mass and exactly its enthalpy, the pack and the
   ! free water it gives together, and no depth: what freezes fills the pores.
   !---~---
   subroutine test_rain_on_snow()
      type(water_layer_state_t)  :: snow
      type(water_material_t) :: p
      type(water_stage_out_t) :: st
      real(wp) :: e0, swe0, d0
      print '(a)', 'test_rain_on_snow:'
      snow = cold_pack(10.0_wp, 265.0_wp, p)
      e0   = snow%energy
      swe0 = snow%water
      d0   = snow%matrix_depth
      call water_layer_stage(2.0e-3_wp, 0.0_wp, 278.0_wp, 0.0_wp, 0.0_wp, 900.0_wp, p, MIN_DEPTH, 268.0_wp, 0.0_wp, &
                             snow, st)                                       ! warm rain (2 mm/900s)
      call check('rain adds mass', snow%water + st%free%water, swe0 + 2.0e-3_wp * 900.0_wp, 1.0e-9_wp)
      call check_true('rain adds energy (warms/refreezes pack)', snow%energy > e0, snow%energy - e0)
      call check('rain enthalpy exact', snow%energy + st%free%energy - e0,                                       &
                 2.0e-3_wp * 900.0_wp * internal_energy_liquid(278.0_wp), 1.0e-6_wp)
      call check('rain adds no depth', snow%matrix_depth, d0, 0.0_wp)
      call check('and the pack is denser by it', water_layer_bulk_density(snow%water, snow%matrix_depth),        &
                 snow%water / d0, 1.0e-12_wp)
   end subroutine test_rain_on_snow

   !---~---
   !   A thin pack (0.3 kg/m2, a fifth of the ground) under the bare-ground conductance of an
   ! open patch: an explicit flux at t_3ple, or the nonlinear sublimation at a linearized
   ! temperature, would drive it hundreds of kelvin past its balance.
   !---~---
   subroutine test_stiff_coupling()
      type(water_layer_state_t)  :: snow
      type(surface_env_t) :: env
      type(water_material_t) :: p
      real(wp), parameter :: SWE = 0.3_wp, DT = 900.0_wp, AF = 0.2_wp
      real(wp) :: e0, m_in, m_out, e_in, e_out, t_min, gain, vapour, g_top
      integer(ik) :: k
      print '(a)', 'test_stiff_coupling:'
      env%ggnet       = 0.065_wp
      env%rho_air     = 1.19_wp
      env%press       = 101325.0_wp
      env%t_soil_top  = t_3ple
      env%r_soil_top  = 0.025_wp / 1.5_wp

      !--- (a) A cold pack (243 K) under warm, moist canopy air: it warms toward the air, not past it.
      snow = cold_pack(SWE, 243.0_wp, p)
      env%sw_snow  = 0.0_wp
      env%lw_down  = stefan_boltzmann * 243.0_wp ** 4 + 20.0_wp
      env%can_temp = 268.8_wp
      env%can_shv  = 1.79e-3_wp
      e0 = snow%energy
      call snow_step(snow, env, AF, DT, gain, vapour, g_top)
      call check_true('a cold pack under warm air and on warm soil warms, and stays frozen',             &
                      snow%temp > 243.0_wp .and. snow%temp < t_3ple .and.                &
                      snow%fliq == 0.0_wp, snow%temp)
      call check('its energy change is what its surroundings gave up', snow%energy - e0, gain,       &
                 1.0e-9_wp * abs(e0))

      !---~---
      !   (b) A melting pack (10 % liquid) under cold, dry canopy air: it freezes through the plateau and
      ! cools toward the air; over a day it stays near the air, and both ledgers close.
      !---~---
      snow = water_layer_state_t()
      snow%water = SWE
      snow%energy = temp_to_internal_energy(0.0_wp, SWE, t_3ple, 0.1_wp)
      snow%matrix_depth  = SWE / p%fresh_density
      env%lw_down  = stefan_boltzmann * t_3ple ** 4 - 60.0_wp
      env%can_temp = 255.0_wp
      env%can_shv  = 3.0e-4_wp
      m_in  = SWE
      m_out = 0.0_wp
      e_in  = snow%energy
      e_out = 0.0_wp
      t_min = t_3ple
      do k = 1_ik, 96_ik
         call snow_step(snow, env, AF, DT, gain, vapour, g_top)
         m_out = m_out + vapour
         e_in  = e_in  + gain
         if (k == 1_ik) call check_true('one step freezes the liquid without passing the air temperature', &
                                        snow%fliq == 0.0_wp .and. snow%temp > env%can_temp - 5.0_wp, &
                                        snow%temp)
         call drain(snow, p, env%t_soil_top, m_out, e_out, e_in)
         if (.not. ground_store_present(snow)) exit
         t_min = min(t_min, snow%temp)
      end do
      call check_true('over a day the pack never falls more than 5 K below the air', t_min > env%can_temp - 5.0_wp, &
                      t_min)
      call check('MASS closure', snow%water, m_in - m_out, 1.0e-12_wp)
      call check('ENERGY closure', snow%energy, e_in - e_out, 1.0e-6_wp * abs(e0))
   end subroutine test_stiff_coupling

   !--- 7. Each tile's radiation at its own optics; the cover-weighted sum is the blended ground's.
   subroutine test_tile_radiation()
      type(water_layer_state_t)  :: snow
      type(water_material_t) :: p
      real(wp), parameter :: DN_VIS = 150.0_wp, DN_NIR = 150.0_wp, DN_LW = 280.0_wp   ! [W/m2] reaching the ground
      real(wp), parameter :: SOIL_VIS = 0.10_wp, SOIL_NIR = 0.20_wp, SOIL_EMISS = 0.95_wp, T_SOIL = 275.0_wp
      real(wp) :: fc, a_vis, a_nir, emiss, temp, rnet_soil, sw_snow, sw_ground, rnet_snow, blended
      print '(a)', 'test_tile_radiation:'
      snow = cold_pack(8.0_wp, 263.0_wp, p)                  ! 3.2 cm deep: cover ~0.47
      call internal_energy_to_temp(snow%energy, snow%water, 0.0_wp, snow%temp, snow%fliq)
      fc = water_layer_cover_of(snow, .true., p, 0.01_wp)
      call check_true('a patchy pack', fc > 0.05_wp .and. fc < 0.95_wp, fc)
      call ground_tile_radiation(DN_VIS, DN_NIR, DN_LW, SOIL_VIS, SOIL_NIR, SOIL_EMISS, T_SOIL, snow, p,          &
                                 rnet_soil, sw_snow, sw_ground)
      rnet_snow = sw_snow + p%emissivity * (DN_LW - stefan_boltzmann * snow%temp ** 4)
      a_vis = SOIL_VIS
      a_nir = SOIL_NIR
      emiss = SOIL_EMISS
      temp  = T_SOIL
      call snow_blend_ground_optics(snow, fc, p, a_vis, a_nir, emiss, temp)
      blended = DN_VIS * (1.0_wp - a_vis) + DN_NIR * (1.0_wp - a_nir) + emiss * DN_LW                        &
                - emiss * stefan_boltzmann * temp ** 4
      call check('the tiles sum to the blended ground', fc * rnet_snow + (1.0_wp - fc) * rnet_soil, blended,   &
                 1.0e-9_wp)
   end subroutine test_tile_radiation

   !---~---
   !   8. The bare ground's skin fluxes to the canopy air, at 15 g/kg; and the exchange's bare
   ! ground, with no pack, under canopy air whose specific heat is set.
   !---~---
   subroutine test_ground_balance()
      real(wp), parameter :: GGNET = 0.02_wp, RHO = 1.2_wp, T_GROUND = 300.0_wp, Q_CAS = 0.015_wp
      type(surface_env_t) :: env
      real(wp) :: cp_cas, h_bare, le_soil, h_2k, h_ground, le_ground, g_top, ground_rad
      type(store_exchange_t) :: stores
      integer(ik) :: it_stack, fb_stack
      print '(a)', 'test_ground_balance:'
      cp_cas = cas_heat_capacity(Q_CAS)
      call ground_surface_fluxes(T_GROUND, T_GROUND, GGNET, RHO, cp_cas, 0.0_wp, h_bare, le_soil)
      call check_true('H_bare = 0 (no gradient)', abs(h_bare) < 1.0e-9_wp, h_bare)
      call check_true('LE_soil = 0 (no evap)',    abs(le_soil) < 1.0e-30_wp, le_soil)
      call ground_surface_fluxes(T_GROUND, T_GROUND - 1.0_wp, GGNET, RHO, cp_cas, 0.0_wp, h_bare, le_soil)
      call ground_surface_fluxes(T_GROUND, T_GROUND - 2.0_wp, GGNET, RHO, cp_cas, 0.0_wp, h_2k, le_soil)
      call check_true('warm ground over cooler canopy air heats it', h_bare > 0.0_wp, h_bare)
      call check('twice the difference, twice the sensible flux', h_2k, 2.0_wp * h_bare, 1.0e-12_wp * h_2k)
      call check('the bare ground''s sensible flux over g rho dT is the canopy air''s specific heat',          &
                 h_bare / (GGNET * RHO), cp_cas, 1.0e-12_wp * cp_cas)
      env%can_temp   = T_GROUND - 1.0_wp
      env%can_shv    = Q_CAS
      env%ggnet      = GGNET
      env%rho_air    = RHO
      env%cp_cas     = cp_cas
      env%t_soil_top = T_GROUND
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores, it_stack, fb_stack)
      call check('the ground surface''s exchange takes the canopy air''s specific heat', h_ground, h_bare,        &
                 1.0e-12_wp * h_bare)
   end subroutine test_ground_balance

   !--- 9. The water layer's stage, before a step.
   subroutine test_water_layer_stage()
      real(wp), parameter :: DT = 150.0_wp, LIGHT = 5.0e-6_wp   ! [s], [kg/m2/s] 7.5e-4 kg/m2 a step
      type(water_material_t)    :: p
      type(water_layer_state_t) :: layer
      type(water_stage_out_t)   :: st
      type(ground_surface_layer_params_t) :: o
      type(porous_layer_state_t) :: thin
      real(wp) :: fallen, heat, e_in
      integer(ik) :: k
      logical  :: film_all_along
      print '(a)', 'test_water_layer_stage:'

      !--- (a) Light snowfall on frozen ground: a film at the ground's temperature until it is a layer.
      layer = water_layer_state_t()
      fallen = 0.0_wp
      heat   = 0.0_wp
      e_in   = 0.0_wp
      film_all_along = .true.
      do k = 1_ik, 60_ik
         call water_layer_stage(0.0_wp, LIGHT, 268.0_wp, 0.0_wp, 0.0_wp, DT, p, MIN_DEPTH, 265.0_wp, 0.0_wp, layer, st)
         fallen = fallen + LIGHT * DT
         heat   = heat + st%film_heat
         e_in   = e_in + st%inflow_enth
         if (k < 30_ik) film_all_along = film_all_along .and. .not. st%is_layer .and. abs(layer%temp - 265.0_wp) < 1.0e-9_wp
         call check_true('light snow on frozen ground: nothing drains', st%free%water == 0.0_wp, st%free%water)
      end do
      call check_true('light snow on frozen ground: a film at the ground''s temperature', film_all_along)
      call check_true('light snow on frozen ground: then a layer', st%is_layer, layer%matrix_depth)
      call check('light snow on frozen ground: it keeps all that fell', layer%water, fallen, 1.0e-12_wp)
      call check('light snow on frozen ground: its energy is what fell and what the ground gave', layer%energy,     &
                 e_in + heat, 1.0e-9_wp * abs(e_in))

      !--- (b) A dusting on warm ground: the film melts, the ground paying for it, and drains.
      layer = water_layer_state_t()
      call water_layer_stage(0.0_wp, 0.02_wp / DT, 270.0_wp, 0.0_wp, 0.0_wp, DT, p, MIN_DEPTH, 280.0_wp, 1.0_wp,      &
                             layer, st)
      call check_true('a dusting on warm ground is a film', .not. st%is_layer)
      call check('it drains whole', st%free%water, 0.02_wp, 1.0e-15_wp)
      call check_true('the ground pays for its melt', st%film_heat > 0.02_wp * latent_heat_fusion, st%film_heat)
      call check('and it leaves at the ground''s temperature', st%free%energy, 0.02_wp * internal_energy_liquid(280.0_wp), &
                 1.0e-9_wp)
      call check_true('nothing is left', .not. ground_store_present(layer))

      !---~---
      !   (c) A 3 mm pond holding a round-off share of ice: its water free, the pond a layer that covers
      ! what a 3 mm layer of its density covers (water_layer_exposed: the free water with the rest).
      !---~---
      layer = water_layer_state_t()
      layer%water  = 3.0_wp
      layer%energy = 3.0_wp * internal_energy_liquid(t_3ple) - 1.0e-6_wp
      call water_layer_stage(0.0_wp, 0.0_wp, 280.0_wp, 0.0_wp, 0.0_wp, DT, p, MIN_DEPTH, 275.0_wp, 1.0_wp, layer, st)
      call check_true('a 3 mm pond with a trace of ice is a layer', st%is_layer)
      call check_true('its water is free, but for the trace and what it holds', st%free%water > 3.0_wp - 1.0e-9_wp, &
                      st%free%water)
      call check('and it covers what 3 mm of water does', water_layer_cover_of(water_layer_exposed(layer, st%free),  &
                 st%is_layer, p, 0.01_wp), tanh(3.0e-3_wp / (p%cover_scale_factor * 0.01_wp * 10.0_wp)), 1.0e-6_wp)

      !--- (d) A cold deep pack and a pond at the triple point, one store: the pack's cold freezes it in its pores.
      layer = cold_pack(40.0_wp, 255.0_wp, p)
      layer%water  = layer%water + 2.0_wp
      layer%energy = layer%energy + 2.0_wp * internal_energy_liquid(t_3ple)
      call water_layer_stage(0.0_wp, 0.0_wp, 255.0_wp, 0.0_wp, 0.0_wp, DT, p, MIN_DEPTH, 268.0_wp, 0.0_wp, layer, st)
      call check('the pack freezes the pond', st%free%water, 0.0_wp, 0.0_wp)
      call check_true('all ice', layer%fliq == 0.0_wp .and. layer%temp < t_3ple, layer%temp)
      call check('within its pores: no deeper', layer%matrix_depth, 40.0_wp / p%fresh_density, 1.0e-15_wp)

      !--- (e) A warm pond: all free water, the layer empty.
      layer = water_layer_state_t()
      layer%water  = 2.0_wp
      layer%energy = 2.0_wp * internal_energy_liquid(285.0_wp)
      call water_layer_stage(0.0_wp, 0.0_wp, 285.0_wp, 0.0_wp, 0.0_wp, DT, p, MIN_DEPTH, 285.0_wp, 1.0_wp, layer, st)
      call check_true('a warm pond is all free water', st%free%water == 2.0_wp .and.                              &
                      st%free%energy == 2.0_wp * internal_energy_liquid(285.0_wp) .and.                           &
                      .not. ground_store_present(layer))

      !--- (f) A dusting on a thin litter (40 J/K dry): the two at one temperature, the litter too small to melt it.
      layer = water_layer_state_t()
      thin%matter = 40.0_wp / cp_organic
      thin%water  = 0.005_wp
      call water_layer_stage(0.0_wp, 0.02_wp / DT, 270.0_wp, 0.0_wp, 0.0_wp, DT, p, MIN_DEPTH, 275.0_wp, 1.0_wp,      &
                             layer, st, below=ground_store_t(water=0.005_wp,                                          &
                                                             energy=temp_to_internal_energy(40.0_wp, 0.005_wp,        &
                                                                                            275.0_wp, 1.0_wp)),       &
                             below_dry_cap=40.0_wp,                                                                   &
                             below_depth=thin%matter / o%porous(POROUS_LITTER)%bulk_density,                          &
                             below_curve=porous_layer_curve(o%porous(POROUS_LITTER), thin))
      call check_true('a dusting on thin litter: a film near the triple point with it, the litter''s water on its curve', &
                      .not. st%is_layer .and. abs(layer%temp - t_3ple) < 0.5_wp, layer%temp)
      call check_true('the litter gave it no more than it held above its all-ice state at the triple point',     &
                      st%film_heat > 0.0_wp .and. st%film_heat < temp_to_internal_energy(40.0_wp, 0.005_wp, 275.0_wp,  &
                      1.0_wp) - temp_to_internal_energy(40.0_wp, 0.005_wp, t_3ple, 0.0_wp), st%film_heat)
   end subroutine test_water_layer_stage

   !--- 10a. The porous layers: 1.35 kg of litter (3 cm) over 2.15 kg of moss (5 cm), dry, each at its temperature.
   function two_layers(o, t_litter, t_moss) result(layer)
      type(ground_surface_layer_params_t), intent(in) :: o
      real(wp),                            intent(in) :: t_litter, t_moss   ! [K]
      type(porous_layer_state_t) :: layer(N_POROUS)
      layer = porous_layer_state_t()
      layer%matter = [0.03_wp * o%porous(POROUS_LITTER)%bulk_density, 0.05_wp * o%porous(POROUS_MOSS)%bulk_density]
      call set_layer_temps(layer, o, [t_litter, t_moss])
   end function two_layers

   !--- Each layer's energy at a temperature on its freezing curve, and its temperature read back.
   subroutine set_layer_temps(layer, o, temps)
      type(porous_layer_state_t),          intent(inout) :: layer(N_POROUS)
      type(ground_surface_layer_params_t), intent(in)    :: o
      real(wp),                            intent(in)    :: temps(N_POROUS)
      integer(ik) :: k
      do k = 1_ik, N_POROUS
         layer(k)%energy = porous_layer_energy_at(layer(k), o%porous(k), temps(k))
      end do
      call porous_layers_diagnose(layer, o)
   end subroutine set_layer_temps

   subroutine test_porous_layers()
      real(wp), parameter :: HEAD_10KPA = -1.0e4_wp / (1000.0_wp * 9.80665_wp)   ! [m]
      type(ground_surface_layer_params_t)   :: o
      type(porous_layer_state_t)    :: layer(N_POROUS)
      type(soil_layer_curve_t) :: c
      type(ground_store_t)     :: pond
      real(wp) :: cap_litter, cap_moss, fill, e0, e_moss
      print '(a)', 'test_porous_layers:'
      layer = two_layers(o, 280.0_wp, 283.0_wp)
      cap_litter = o%porous(POROUS_LITTER)%theta_sat * 1000.0_wp * 0.03_wp            ! 29.4 kg/m2
      cap_moss   = o%porous(POROUS_MOSS)%theta_sat   * 1000.0_wp * 0.05_wp            ! 42.5 kg/m2
      call check('a layer holds its theta_sat over its depth',                                                   &
                 porous_layer_capacity(layer(POROUS_MOSS)%matter, o%porous(POROUS_MOSS)%bulk_density,           &
                                       o%porous(POROUS_MOSS)%theta_sat), cap_moss,                              &
                 1.0e-12_wp)
      call check('each layer its own temperature: the litter''s', layer(POROUS_LITTER)%temp, 280.0_wp, 1.0e-9_wp)
      call check('each layer its own temperature: the moss''s', layer(POROUS_MOSS)%temp, 283.0_wp, 1.0e-9_wp)
      layer%water  = [cap_litter, 2.0_wp]
      call set_layer_temps(layer, o, [270.0_wp, 283.0_wp])
      call check('a frozen litter: its temperature from its energy on its curve', layer(POROUS_LITTER)%temp,      &
                 270.0_wp, 1.0e-7_wp)
      call check_true('a frozen litter over a thawed moss, the litter''s water partly liquid on its curve',         &
                      layer(POROUS_LITTER)%fliq > 0.0_wp .and. layer(POROUS_LITTER)%fliq < 1.0_wp                   &
                      .and. layer(POROUS_MOSS)%fliq == 1.0_wp, layer(POROUS_LITTER)%fliq)
      call check('a saturated layer conducts at its wet conductivity',                                           &
                 porous_layer_resistance(0.03_wp, 1.0_wp, o%porous(POROUS_LITTER)%k_dry, o%porous(POROUS_LITTER)%k_wet), &
                 0.03_wp / o%porous(POROUS_LITTER)%k_wet, 1.0e-12_wp)
      fill = 2.0_wp / cap_moss
      call check('a layer''s dried top at its dry conductivity, the wet rest at its wet one',                    &
                 porous_layer_resistance(0.05_wp, fill, o%porous(POROUS_MOSS)%k_dry, o%porous(POROUS_MOSS)%k_wet),  &
                 0.05_wp * ((1.0_wp - fill) / o%porous(POROUS_MOSS)%k_dry + fill / o%porous(POROUS_MOSS)%k_wet),   &
                 1.0e-12_wp)
      !--- Each layer's curve: its material's at its depth, the residual no more than the liquid it holds.
      c = porous_layer_curve(o%porous(POROUS_MOSS), layer(POROUS_MOSS))
      call check('the moss''s curve: its depth', c%dz, 0.05_wp, 1.0e-12_wp)
      call check('the moss at -10 kPa holds 0.29 (Zhang et al. 2010)',                                             &
                 soil_head_to_theta(c%retention, HEAD_10KPA, c%theta_sat, c%theta_res, c%par_a, c%par_n), 0.294_wp,  &
                 1.0e-3_wp)
      call check('a moss above its residual keeps it', c%theta_res, o%porous(POROUS_MOSS)%theta_res, 0.0_wp)
      layer(POROUS_LITTER)%water = 0.0_wp
      c = porous_layer_curve(o%porous(POROUS_LITTER), layer(POROUS_LITTER))
      call check('an empty litter keeps its material''s residual', c%theta_res, o%porous(POROUS_LITTER)%theta_res, 0.0_wp)
      layer(POROUS_LITTER)%water = cap_litter
      c = porous_layer_curve(o%porous(POROUS_LITTER), layer(POROUS_LITTER))
      call check('the litter at -10 kPa holds 0.14, near drained litter''s 0.12 (Putuhena & Cordery 1996)',        &
                 soil_head_to_theta(c%retention, HEAD_10KPA, c%theta_sat, c%theta_res, c%par_a, c%par_n), 0.1415_wp, &
                 1.0e-3_wp)
      layer%fliq = 0.0_wp
      c = porous_layer_curve(o%porous(POROUS_MOSS), layer(POROUS_MOSS))
      call check('a frozen layer keeps its material''s curve', c%theta_res, o%porous(POROUS_MOSS)%theta_res, 0.0_wp)
      !--- A frozen litter 2 kg beyond its saturation: the excess passes to the water layer, the litter keeps its temperature.
      layer = two_layers(o, 265.0_wp, 283.0_wp)
      layer%water  = [cap_litter + 2.0_wp, 2.0_wp]
      call set_layer_temps(layer, o, [265.0_wp, 283.0_wp])
      e0 = sum(layer%energy)
      e_moss = layer(POROUS_MOSS)%energy
      pond = ground_store_t()
      call porous_layers_spill(layer, o, pond)
      call check('beyond saturation: the excess passes to the water layer', pond%water, 2.0_wp, 1.0e-12_wp)
      call check('beyond saturation: the litter keeps its saturation', layer(POROUS_LITTER)%water, cap_litter,     &
                 1.0e-12_wp)
      call check('beyond saturation: the litter keeps its temperature', layer(POROUS_LITTER)%temp, 265.0_wp,      &
                 1.0e-7_wp)
      call check('beyond saturation: the energy goes with the water', sum(layer%energy) + pond%energy, e0,        &
                 1.0e-9_wp * abs(e0))
      call check('beyond saturation: a layer within its pores keeps all it holds', layer(POROUS_MOSS)%energy,     &
                 e_moss, 0.0_wp)
   end subroutine test_porous_layers

   !--- A porous layer as the stack reads it: water [kg/m2] at temp, its matter's heat, its resistance.
   function porous_env(o, k, matter, water, temp) result(pe)
      type(ground_surface_layer_params_t), intent(in) :: o
      integer(ik),                         intent(in) :: k
      real(wp),                            intent(in) :: matter, water, temp   ! [kg/m2], [kg/m2], [K]
      type(porous_node_env_t) :: pe
      type(porous_layer_state_t) :: layer
      real(wp) :: depth
      layer%matter = matter
      layer%water  = water
      pe%dry_cap = cp_organic * matter
      pe%water   = water
      pe%energy  = porous_layer_energy_at(layer, o%porous(k), temp)
      depth = matter / o%porous(k)%bulk_density
      pe%depth = depth
      pe%curve = porous_layer_curve(o%porous(k), layer)
      call porous_layer_temp(pe%energy / depth, water / (1000.0_wp * depth), pe%dry_cap / depth, pe%curve,         &
                             pe%temp, pe%fliq)
      pe%resistance = porous_layer_resistance(depth, min(1.0_wp, water / porous_layer_capacity(matter,               &
                                              o%porous(k)%bulk_density, o%porous(k)%theta_sat)),                   &
                                              o%porous(k)%k_dry, o%porous(k)%k_wet)
   end function porous_env

   !---~---
   !   10b. One node's balance (surface_node_balance), and the ground surface with the porous
   ! layers: a node held at the network's evaporation or dew gives that rate per unit ground area,
   ! even where the pack covers it all; its energy closes and the exchange it takes gives back the
   ! temperature it solved for; a node whose vapour is the canopy air's evaporates no more than it
   ! holds; and a pack on the ice line is the closed form of its own balance.
   !---~---
   subroutine test_porous_balance()
      real(wp), parameter :: DT = 3600.0_wp
      type(ground_surface_layer_params_t) :: o
      type(surface_env_t)    :: env
      type(surface_node_t)   :: node
      type(node_flux_t)      :: ofx
      type(water_material_t)    :: p
      real(wp) :: h_ground, le_ground, g_top, ground_rad, s_in, o_in
      type(store_exchange_t) :: stores
      integer(ik) :: it_stack, fb_stack
      print '(a)', 'test_porous_balance:'
      call cold_env(env)
      env%dt           = DT
      env%can_temp     = 270.0_wp
      env%t_soil_top   = 285.0_wp
      env%ggnet        = 0.05_wp
      env%rnet_soil    = 0.0_wp
      env%ground_emiss = 0.95_wp
      env%porous(POROUS_LITTER) = porous_env(o, POROUS_LITTER, 0.03_wp * o%porous(POROUS_LITTER)%bulk_density,      &
                                             1.0_wp, 280.0_wp)   ! litter alone: 1.35 kg/m2, 2.6 kJ/m2/K dry
      !--- Held at the water network's rate, evaporation or dew: the layer gives that rate, its energy closes, it takes its root.
      env%porous(POROUS_LITTER)%vapour = 1.0e-5_wp
      node = porous_node(env, 1.0_wp)
      call surface_node_balance(node, env, ofx)
      call check_true('a thin layer relaxes over an hour without passing its surroundings',                         &
                      ofx%t_surf > 268.0_wp .and. ofx%t_surf < 286.0_wp, ofx%t_surf)
      call check('held: it gives the network''s evaporation', ofx%w_flux, 1.0e-5_wp, 1.0e-18_wp)
      call check('held: its energy, what its surroundings give it', ofx%energy_in,                                 &
                 ofx%rad - ofx%h - ofx%le - ofx%g_below, 1.0e-9_wp)
      call check('held: the exchange it takes gives back its root', committed_temp(node, env%dt, ofx), ofx%t_surf,   &
                 1.0e-6_wp)
      env%porous(POROUS_LITTER)%vapour = -2.0e-6_wp
      node = porous_node(env, 1.0_wp)
      call surface_node_balance(node, env, ofx)
      call check('held: dew at the network''s rate', ofx%w_flux, -2.0e-6_wp, 1.0e-18_wp)
      call check_true('held: dew gives the layer its latent heat', ofx%le < 0.0_wp, ofx%le)
      call check('held: dew, the exchange it takes gives back its root', committed_temp(node, env%dt, ofx),         &
                 ofx%t_surf, 1.0e-6_wp)
      !--- Under a pack that covers all of it, the network's evaporation is still the layer's, per unit ground area.
      env%porous(POROUS_LITTER)%vapour = 1.0e-5_wp
      node = porous_node(env, 0.0_wp)
      call surface_node_balance(node, env, ofx)
      call check('held under a full cover: the network''s rate all the same', ofx%w_flux, 1.0e-5_wp, 1.0e-18_wp)
      call check_true('and its latent heat with it', ofx%le > 0.0_wp, ofx%le)
      call check('held under a full cover: its root', committed_temp(node, env%dt, ofx), ofx%t_surf, 1.0e-6_wp)
      node = porous_node(env, 0.5_wp)
      call surface_node_balance(node, env, ofx)
      call check('held under half a cover: the network''s rate', ofx%w_flux, 1.0e-5_wp, 1.0e-18_wp)
      !--- Its vapour the canopy air's: a trace of water evaporates no more than it holds, and the node takes its root.
      node = porous_node(env, 1.0_wp)
      node%vapour_held = .false.
      node%water  = 1.0e-6_wp
      node%energy = temp_to_internal_energy(node%dry_cap, node%water, 280.0_wp, 1.0_wp)
      call internal_energy_to_temp(node%energy, node%water, node%dry_cap, node%temp, node%fliq)
      call surface_node_balance(node, env, ofx)
      call check_true('a trace of water evaporates no more than it holds',                                         &
                      ofx%w_flux * DT <= node%water * (1.0_wp + 1.0e-12_wp), ofx%w_flux)
      call check('a trace of water: the node still takes its root', committed_temp(node, env%dt, ofx), ofx%t_surf,  &
                 1.0e-6_wp)
      call test_pack_closed_form(env, p)
      !--- The ground with the layer: the soil takes its conduction less the soil's own evaporation's heat.
      env%porous(POROUS_LITTER)%vapour = 0.0_wp
      env%soil_evap = 1.0e-6_wp
      env%snowfac   = 0.0_wp
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores, it_stack, fb_stack)
      s_in = stores%energy_in(STORE_WATER)
      o_in = stores%energy_in(STORE_WATER + POROUS_LITTER)
      node = porous_node(env, 1.0_wp)
      call surface_node_balance(node, env, ofx)
      call check('the soil takes what the layer conducts, less its evaporation''s heat', g_top,                    &
                 ofx%g_below - env%soil_evap * enthalpy_vapor(env%t_soil_top), 1.0e-9_wp)
      call check('the layer is the ground''s surface: its sensible heat, the ground''s', h_ground, ofx%h, 1.0e-12_wp)
      call check('no pack, no pack exchange', s_in, 0.0_wp, 0.0_wp)
      call check('the layer''s own energy', o_in, ofx%energy_in, 1.0e-12_wp)
      !--- A pack on the layer: it conducts into the layer, and what the ground absorbs is what both keep.
      env%params  = p
      env%pack    = cold_pack(30.0_wp, 265.0_wp, p)
      env%snowfac = 0.9_wp
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores, it_stack, fb_stack)
      call check('pack and layer: the ground''s energy closes', ground_rad - h_ground - le_ground - g_top,          &
                 sum(stores%energy_in(1:N_SURFACE_STORE)), 1.0e-9_wp)
   end subroutine test_porous_balance

   !---~---
   !   A cold pack on the ice line is its balance's closed form: its energy over the step is
   ! E^n + dt af (r0 + dr/dT (T - T^n)), each term per unit of the pack's own area, at its root
   ! T = (E^n + dt af (r0 - dr/dT T^n)) / (W c_ice - dt af dr/dT).
   !---~---
   subroutine test_pack_closed_form(env0, p)
      type(surface_env_t),    intent(in) :: env0
      type(water_material_t), intent(in) :: p
      type(surface_env_t)  :: env
      type(surface_node_t) :: node
      type(node_flux_t)    :: flux
      real(wp), parameter :: AF = 0.7_wp, G_BASE = 0.4_wp
      real(wp) :: dh_dt, w_n, dw_dt, l_vap, r0, drdt, t_root
      env = env0
      env%can_temp = 262.0_wp
      env%can_shv  = 1.0e-3_wp
      node%water   = 30.0_wp
      node%energy  = temp_to_internal_energy(0.0_wp, node%water, 258.0_wp, 0.0_wp)
      call internal_energy_to_temp(node%energy, node%water, 0.0_wp, node%temp, node%fliq)
      node%area    = AF
      node%emiss   = p%emissivity
      node%rad     = 40.0_wp + p%emissivity * (230.0_wp - stefan_boltzmann * node%temp ** 4)
      node%g_below = AF * G_BASE
      node%t_below = 271.0_wp
      call surface_node_balance(node, env, flux)
      dh_dt = env%ggnet * env%rho_air * env%cp_cas
      w_n   = env%ggnet * env%rho_air * (sat_specific_humidity(node%temp, env%press, 0.0_wp) - env%can_shv)
      dw_dt = env%ggnet * env%rho_air * sat_specific_humidity_temp_deriv(node%temp, env%press, 0.0_wp)
      l_vap = enthalpy_vapor(node%temp) - cp_ice * node%temp
      r0    = node%rad - dh_dt * (node%temp - env%can_temp) - l_vap * w_n - G_BASE * (node%temp - node%t_below)
      drdt  = -(4.0_wp * p%emissivity * stefan_boltzmann * node%temp ** 3 + dh_dt + l_vap * dw_dt + G_BASE)
      t_root = (node%energy + env%dt * AF * (r0 - drdt * node%temp)) / (node%water * cp_ice - env%dt * AF * drdt)
      call check_true('a cold pack stays on the ice line', flux%t_surf < t_3ple, flux%t_surf)
      call check('a pack on the ice line: its root is the closed form''s', flux%t_surf, t_root, 1.0e-9_wp)
   end subroutine test_pack_closed_form

   !--- The top porous layer as the ground surface's balance builds its node, at the pack's cover.
   function porous_node(env, open_share) result(node)
      type(surface_env_t), intent(in) :: env
      real(wp),            intent(in) :: open_share   ! [-] the share of the ground the pack leaves open
      type(surface_node_t) :: node
      node%dry_cap = env%porous(POROUS_LITTER)%dry_cap
      node%on_curve = .true.
      node%depth    = env%porous(POROUS_LITTER)%depth
      node%curve    = env%porous(POROUS_LITTER)%curve
      node%water   = env%porous(POROUS_LITTER)%water
      node%energy  = env%porous(POROUS_LITTER)%energy
      node%temp    = env%porous(POROUS_LITTER)%temp
      node%fliq    = env%porous(POROUS_LITTER)%fliq
      node%area        = open_share
      node%emiss       = env%ground_emiss
      node%rad         = env%sw_ground + node%emiss * (env%lw_down - stefan_boltzmann * node%temp ** 4)
      node%vapour_held = .true.
      node%vapour_rate = env%porous(POROUS_LITTER)%vapour
      node%g_below     = conductance_in_series(env%porous(POROUS_LITTER)%resistance, env%r_soil_top)
      node%t_below     = env%t_soil_top
   end function porous_node

   !---~---
   !   A node's temperature after it takes an exchange over the balance's step: a porous node's
   ! water is what it keeps, its held vapour apart, and its live vapour leaves from it.
   !---~---
   real(wp) function committed_temp(node, dt, flux) result(t)
      type(surface_node_t), intent(in) :: node
      real(wp),             intent(in) :: dt
      type(node_flux_t),    intent(in) :: flux
      real(wp) :: fl, held
      if (node%on_curve) then
         held = 0.0_wp
         if (node%vapour_held) held = node%vapour_rate
         call porous_layer_temp((node%energy + dt * flux%energy_in) / node%depth,                                   &
                                (node%water - dt * (flux%w_flux - held)) / (1000.0_wp * node%depth),                &
                                node%dry_cap / node%depth, node%curve, t, fl)
         return
      end if
      call internal_energy_to_temp(node%energy + dt * flux%energy_in, node%water - dt * flux%w_flux, node%dry_cap,  &
                                   t, fl)
   end function committed_temp

   !--- 12. The stack: a thin melting pack over a 0.3 mm litter over warm soil.
   subroutine test_stack()
      type(surface_env_t)  :: env
      type(water_material_t) :: p
      type(surface_node_t) :: node(2), one(1), base(2)
      type(node_flux_t)    :: flux(2), fl1(1), single
      integer(ik) :: it, fb
      real(wp) :: t_newton(2), gain, given
      print '(a)', 'test_stack:'
      call cold_env(env)
      env%dt      = 900.0_wp
      env%can_temp = 276.0_wp
      env%can_shv  = 4.0e-3_wp
      env%t_soil_top = 280.0_wp
      !--- The pack: 0.5 kg/m2 at the plateau, a tenth liquid, over half the ground.
      node(1)%water  = 0.5_wp
      node(1)%energy = temp_to_internal_energy(0.0_wp, 0.5_wp, t_3ple, 0.1_wp)
      node(1)%area   = 0.5_wp
      node(1)%emiss  = p%emissivity
      !--- The litter: 0.3 mm (13.5 g/m2), half its pores filled, at 278 K, its vapour the network's.
      node(2)%dry_cap     = cp_organic * 0.0135_wp
      node(2)%water       = 0.13_wp
      node(2)%energy      = temp_to_internal_energy(cp_organic * 0.0135_wp, 0.13_wp, 278.0_wp, 1.0_wp)
      node(2)%area        = 0.5_wp
      node(2)%emiss       = 0.95_wp
      node(2)%vapour_held = .true.
      node(2)%vapour_rate = 1.0e-6_wp
      call internal_energy_to_temp(node(1)%energy, node(1)%water, 0.0_wp, node(1)%temp, node(1)%fliq)
      call internal_energy_to_temp(node(2)%energy, node(2)%water, node(2)%dry_cap, node(2)%temp, node(2)%fliq)
      node(1)%rad     = 150.0_wp + node(1)%emiss * (env%lw_down - stefan_boltzmann * node(1)%temp ** 4)
      node(2)%rad     = 100.0_wp + node(2)%emiss * (env%lw_down - stefan_boltzmann * node(2)%temp ** 4)
      node(1)%g_below = 0.5_wp / 0.05_wp
      node(2)%g_above = node(1)%g_below
      node(2)%g_below = conductance_in_series(3.0e-3_wp, env%r_soil_top)
      node(2)%t_below = env%t_soil_top
      base = node
      call surface_stack_balance(node, 2_ik, env, flux, it, fb)
      t_newton = flux%t_surf
      call check_true('the stack: Newton converges', it > 0_ik .and. fb == 0_ik, real(it, wp))
      gain  = flux(1)%energy_in + flux(2)%energy_in
      given = flux(1)%rad + flux(2)%rad - flux(1)%h - flux(2)%h - flux(1)%le - flux(2)%le - flux(2)%g_below
      call check('the stack: what the nodes take is what they absorb less what they give', gain, given,            &
                 1.0e-9_wp * abs(given))
      call check('the stack: the pack gives the litter what the litter takes', flux(1)%g_below,                     &
                 base(2)%g_above * (flux(1)%t_surf - flux(2)%t_surf), 1.0e-9_wp)
      node = base
      call surface_stack_balance(node, 2_ik, env, flux, it, fb, newton_max=0_ik)
      call check_true('the stack: the Gauss-Seidel sweeps take over', fb == 1_ik)
      call check('the stack: their root is Newton''s, the pack', flux(1)%t_surf, t_newton(1), 1.0e-6_wp)
      call check('the stack: their root is Newton''s, the litter', flux(2)%t_surf, t_newton(2), 1.0e-6_wp)
      !--- One node: the stack's root is the node's own closed form.
      one(1) = base(2)
      one(1)%g_above = 0.0_wp
      call surface_node_balance(one(1), env, single)
      call surface_stack_balance(one, 1_ik, env, fl1, it, fb)
      call check('the stack of one: its closed form', fl1(1)%t_surf, single%t_surf, 1.0e-9_wp)
      call check('the stack of one: its gain', fl1(1)%energy_in, single%energy_in, 1.0e-7_wp)
   end subroutine test_stack

   !---~---
   !   12b. The stack with the water the network moves: a 0.3 mm litter over 5 cm of moss over
   ! cold soil, 2 kg/m2 of warm free water passing through the litter over the step, half of it
   ! into the moss, a little running off. What the stores and the soil take is what the ground
   ! absorbs less what it gives the air and the runoff carries off; the litter, of little heat
   ! capacity, ends between the water reaching it and the moss below; the Gauss-Seidel sweeps
   ! find Newton's root. A 3 cm litter over the moss, no water passing, ends at its own
   ! temperature, apart from the moss's.
   !---~---
   subroutine test_stack_water()
      real(wp), parameter :: DT = 900.0_wp, T_FREE = 283.0_wp, RUNOFF = 0.1_wp / 900.0_wp
      type(ground_surface_layer_params_t) :: o
      type(surface_env_t)    :: env
      type(store_exchange_t) :: stores
      type(surface_node_t)   :: node(2), base(2)
      type(node_flux_t)      :: flux(2)
      real(wp) :: h_ground, le_ground, g_top, ground_rad, t_newton(2)
      integer(ik) :: it, fb
      print '(a)', 'test_stack_water:'
      call cold_env(env)
      env%dt           = DT
      env%can_temp     = 276.0_wp
      env%can_shv      = 4.0e-3_wp
      env%t_soil_top   = 274.0_wp
      env%ggnet        = 0.02_wp
      env%sw_ground    = 80.0_wp
      env%ground_emiss = 0.95_wp
      env%porous(POROUS_LITTER) = porous_env(o, POROUS_LITTER, 0.0135_wp, 0.13_wp, 278.0_wp)
      env%porous(POROUS_MOSS)   = porous_env(o, POROUS_MOSS, 2.15_wp, 20.0_wp, 274.5_wp)
      env%free_temp = T_FREE
      env%free_down = 2.0_wp / DT                      ! into the litter
      env%runoff    = RUNOFF
      env%porous(POROUS_LITTER)%down = 2.0_wp / DT     ! all of it on into the moss
      env%porous(POROUS_MOSS)%down   = 1.0_wp / DT     ! half of it into the soil
      env%porous(POROUS_MOSS)%water  = env%porous(POROUS_MOSS)%water + 1.0_wp
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores, it, fb)
      call check('water passing: what the stores and the soil take is what the ground absorbs less what it gives', &
                 sum(stores%energy_in(1:N_SURFACE_STORE)) + g_top,                                                   &
                 ground_rad - h_ground - le_ground - RUNOFF * internal_energy_liquid(T_FREE), 1.0e-9_wp)
      call check('water passing: the free water gives what it sends at its temperature',                          &
                 stores%energy_in(STORE_WATER), -(2.0_wp / DT + RUNOFF) * internal_energy_liquid(T_FREE), 1.0e-9_wp)
      !--- The litter's temperature, back from its energy and the water it holds.
      call stack_nodes(env, node)
      base = node
      call surface_stack_balance(node, 2_ik, env, flux, it, fb)
      t_newton = flux%t_surf
      call check_true('a thin litter passing 2 kg/m2 ends between the water reaching it and the moss below',          &
                      flux(1)%t_surf < T_FREE .and. flux(1)%t_surf >= flux(2)%t_surf, flux(1)%t_surf)
      call check('the stack: the litter''s gain, its energy at its root', committed_temp(base(1), DT, flux(1)),       &
                 flux(1)%t_surf, 1.0e-6_wp)
      call check('the stack: the moss''s gain, its energy at its root', committed_temp(base(2), DT, flux(2)),         &
                 flux(2)%t_surf, 1.0e-6_wp)
      node = base
      call surface_stack_balance(node, 2_ik, env, flux, it, fb, newton_max=0_ik)
      call check_true('water passing: the Gauss-Seidel sweeps take over', fb == 1_ik)
      call check('water passing: their root is Newton''s, the litter', flux(1)%t_surf, t_newton(1), 1.0e-6_wp)
      call check('water passing: their root is Newton''s, the moss', flux(2)%t_surf, t_newton(2), 1.0e-6_wp)
      !--- A 3 cm litter, no water passing: its own temperature.
      env%porous(POROUS_LITTER) = porous_env(o, POROUS_LITTER, 1.35_wp, 5.0_wp, 278.0_wp)
      env%porous(POROUS_MOSS)   = porous_env(o, POROUS_MOSS, 2.15_wp, 20.0_wp, 274.5_wp)
      env%free_down = 0.0_wp
      call stack_nodes(env, node)
      call surface_stack_balance(node, 2_ik, env, flux, it, fb)
      call check_true('the litter and the moss each at its own temperature', flux(1)%t_surf - flux(2)%t_surf         &
                      > 1.0_wp, flux(1)%t_surf - flux(2)%t_surf)
   end subroutine test_stack_water

   !--- The litter and the moss as ground_surface_exchange builds them, with the water passing through.
   subroutine stack_nodes(env, node)
      type(surface_env_t),  intent(in)  :: env
      type(surface_node_t), intent(out) :: node(2)
      integer(ik) :: i, k
      do i = 1_ik, 2_ik
         k = merge(POROUS_LITTER, POROUS_MOSS, i == 1_ik)
         node(i)%dry_cap     = env%porous(k)%dry_cap
         node(i)%on_curve    = .true.
         node(i)%depth       = env%porous(k)%depth
         node(i)%curve       = env%porous(k)%curve
         node(i)%water       = env%porous(k)%water
         node(i)%energy      = env%porous(k)%energy
         node(i)%temp        = env%porous(k)%temp
         node(i)%fliq        = env%porous(k)%fliq
         node(i)%emiss       = env%ground_emiss
         node(i)%rad         = env%sw_ground + node(i)%emiss * (env%lw_down - stefan_boltzmann * node(i)%temp ** 4)
         node(i)%vapour_held = .true.
      end do
      node(1)%area    = 1.0_wp
      node(1)%g_below = 1.0_wp / env%porous(POROUS_LITTER)%resistance
      node(2)%g_above = node(1)%g_below
      node(2)%g_below = conductance_in_series(env%porous(POROUS_MOSS)%resistance, env%r_soil_top)
      node(2)%t_below = env%t_soil_top
      node(1)%q_in    = env%free_down * internal_energy_liquid(env%free_temp)
      node(1)%m_out   = env%porous(POROUS_LITTER)%down
      node(2)%m_above = env%porous(POROUS_LITTER)%down
      node(2)%m_out   = env%porous(POROUS_MOSS)%down
   end subroutine stack_nodes

   !---~---
   !   13. The ground's exchange across its litter's freezing onset: the state a tundra's late-summer
   ! snowfall left, a thin wet pack over 11 % of the ground on a saturated litter just above the
   ! triple point over wet moss, swept from 271 to 275 K in canopy-air temperatures 0.002 K apart,
   ! over which the ground's heat changes by about 0.1 W/m2 a step.
   !---~---
   subroutine test_freezing_onset()
      type(ground_surface_layer_params_t) :: o
      type(surface_env_t)    :: env
      type(store_exchange_t) :: stores
      real(wp) :: h_ground, le_ground, g_top, ground_rad, last, jump
      integer(ik) :: it, fb, i
      print '(a)', 'test_freezing_onset:'
      env%dt          = 900.0_wp
      env%lw_down     = 313.0_wp
      env%rnet_soil   = -3.3_wp
      env%ggnet       = 0.0423_wp
      env%rho_air     = 1.282_wp
      env%cp_cas      = 1008.0_wp
      env%press       = 101325.0_wp
      env%can_shv     = 4.04e-3_wp
      env%t_soil_top  = 276.53_wp
      env%r_soil_top  = 0.042_wp
      env%ground_emiss = 0.95_wp
      env%pack%water        = 1.05_wp
      env%pack%energy       = temp_to_internal_energy(0.0_wp, 1.05_wp, t_3ple, 0.1_wp)
      env%pack%matrix_depth = 3.9e-3_wp
      env%snowfac   = 0.113_wp
      env%free_temp = t_3ple
      env%free_down = 9.7e-6_wp
      env%porous(POROUS_LITTER) = porous_env(o, POROUS_LITTER, 292.3_wp / cp_organic, 3.314_wp, 273.33_wp)
      env%porous(POROUS_MOSS)   = porous_env(o, POROUS_MOSS, 4969.1_wp / cp_organic, 51.0_wp, 274.06_wp)
      env%porous(POROUS_LITTER)%vapour = 9.49e-6_wp
      env%porous(POROUS_LITTER)%down   = 2.4e-7_wp
      jump = 0.0_wp
      do i = 0, 2000
         env%can_temp = 271.0_wp + 0.002_wp * real(i, wp)
         call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores, it, fb)
         if (i > 0) jump = max(jump, abs(h_ground + le_ground - last))
         last = h_ground + le_ground
      end do
      call check_true('the freezing onset: the ground''s heat to the canopy air is continuous', jump < 0.5_wp, jump)
   end subroutine test_freezing_onset

end program test_surface
