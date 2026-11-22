! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_surface -- the ground's surface and its snow store (ed_surface;
! docs/science/snow_biophysics.md). The pack steps as the integrator steps it: the ground
! surface's exchange at the pack's state (ground_surface_exchange), then the pack takes the
! step's amounts (snow_take_exchange), with snow_accumulate and snow_drain_meltwater around
! them; its energy ledger is what the surroundings give up by the ground surface's balance.
!   1. COVER (snow_cover_fraction): 0 without snow, monotone in [0, 1], and Niu & Yang
!      (2007)'s fraction for 1 cm of snow.
!   2. ACCUMULATION: a cold snap's pack grows below t_3ple, all ice; its mass and energy
!      close against what fell, what it exchanged and what drained.
!   3. MELT: the pack pins at t_3ple, its liquid rises and drains to the pond, and it
!      melts out; both ledgers close through the plateau.
!   4. ISOTHERMAL: a pack at the temperature of its surroundings, under canopy air
!      saturated over ice, stays as it is with no heat into the soil; under canopy air
!      saturated over liquid it gains mass.
!   5. RAIN ON SNOW (snow_accumulate): warm rain on a cold pack adds its mass and exactly
!      its enthalpy.
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
!   9. THE POND (prepare_snow_stage): a pond below its freezing point starts a layer of its
!      ice and keeps its liquid; a cold pack takes a pond and freezes it, deeper by its water
!      as ice; a pack at the plateau takes none of a pond's liquid; a warm pond with no pack is left as it is; an empty
!      pond leaves a pack as it was; snow on ice mixes the densities by volume. The stores
!      the stage reports for the ledger, with the pond it leaves, are the stores it began from.
!  10. THE ORGANIC LAYER (ed_organic_lib and organic_surface_balance): a slot holds its
!      theta_sat over its depth; one temperature, the slots' energies summing to the total; each
!      slot gives the vapour the water network took from it; heat crosses each slot's dried top
!      at its dry conductivity, in series with its wet part. Each slot's curve is its material's
!      at its depth (the moss's at -10 kPa as Zhang et al. 2010 measured, the litter's near what
!      drained litter holds), its residual low enough that frozen or new matter starts 1 % of
!      the way to saturation; water above saturation (frost on a frozen slot) spills with its
!      energy, the layer keeping its temperature. Its balance closes on its energy, a thin layer of little heat
!      capacity relaxes over an hour without passing its surroundings, a dry layer does not
!      evaporate, and a wet one gives no more than it holds; held at the network's rate it gives
!      that rate, evaporation or dew, and takes its root; the soil takes what it conducts down,
!      less the soil's own evaporation's latent heat, and with a pack on it the ground's energy
!      closes on what the pack and the layer keep.
!  11. INSULATION: a deep pack colder than the soil, in balance with the air, passes the soil
!      the heat its whole depth and the top soil node conduct in series, since its one
!      temperature is its surface's; twice as deep, it passes less.
! The energy datum the pack is held in is test_therm_lib's.
!---~---
program test_surface
   use ed_test_assert,        only : check, check_true, test_report
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : t_3ple, stefan_boltzmann, rho_ice, cp_organic
   use ed_therm_lib,          only : temp_to_internal_energy, internal_energy_to_temp, internal_energy_ice,      &
                                     internal_energy_liquid, sat_specific_humidity, cas_heat_capacity,     &
                                     enthalpy_vapor
   use ed_patch_state_vars,   only : snow_state_t
   use ed_organic_lib,        only : organic_state_t, organic_params_t, SLOT_LITTER, SLOT_MOSS, organic_layer_sums, &
                                     organic_split_energy, organic_take_exchange, organic_slot_capacity,         &
                                     organic_slot_curve, organic_spill
   use ed_soil_hydr_lib,      only : soil_layer_curve_t, soil_head_to_theta
   use ed_biophysics_params, only : snow_params_t
   use ed_soil_water,         only : organic_cover_path, vapour_diffusivity
   use ed_surface,            only : surface_env_t, snow_out_t, organic_out_t, organic_surface_balance,          &
                                     snow_cover_fraction, snow_accumulate,                                      &
                                     snow_drain_meltwater, ground_surface_exchange, snow_take_exchange,         &
                                     snow_blend_ground_optics, ground_tile_radiation, ground_surface_fluxes,    &
                                     prepare_snow_stage
   implicit none

   call test_snow_cover()
   call test_accumulate()
   call test_melt_event()
   call test_isothermal()
   call test_rain_on_snow()
   call test_stiff_coupling()
   call test_tile_radiation()
   call test_ground_balance()
   call test_pond_stage()
   call test_organic_slots()
   call test_organic_balance()
   call test_insulation()

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
      env%k_soil_top  = 1.5_wp
      env%dz_soil_top = 0.05_wp
   end subroutine cold_env

   !--- A pack of swe [kg/m2], all ice at temp, at the default snow density.
   function cold_pack(swe, temp, p) result(snow)
      real(wp),            intent(in) :: swe, temp   ! [kg/m2], [K]
      type(snow_params_t), intent(in) :: p
      type(snow_state_t) :: snow
      snow = snow_state_t()
      snow%nlayer = 1_ik
      snow%swe(1) = swe
      snow%snow_energy(1) = temp_to_internal_energy(0.0_wp, swe, temp, 0.0_wp)
      snow%snow_depth(1)  = swe / p%rho_snow
      snow%density(1)     = p%rho_snow
   end function cold_pack

   !---~---
   !   One step of the pack as the integrator takes it: the ground surface's exchange at the
   ! pack's state and the given cover (ground_surface_exchange, the call the derivatives make),
   ! then the pack takes the step's amounts (snow_take_exchange, the update's). Returns the
   ! energy the surroundings gave the pack [J/m2], by the ground surface's own balance -- its
   ! net radiation less the sensible and latent heat to the canopy air and the heat into the
   ! soil -- the vapour the pack lost [kg/m2], and the heat flux into the soil top [W/m2].
   !---~---
   subroutine snow_step(snow, env, cover, dt, energy_gain, vapour_loss, g_top)
      type(snow_state_t),  intent(inout) :: snow          ! the pack
      type(surface_env_t), intent(in)    :: env           ! its surroundings; the pack and cover are set here
      real(wp),            intent(in)    :: cover, dt     ! [-], [s]
      real(wp),            intent(out)   :: energy_gain, vapour_loss, g_top   ! [J/m2], [kg/m2], [W/m2]
      type(surface_env_t) :: e
      real(wp) :: h_ground, le_ground, ground_rad, snow_energy_in, snow_vapour_out, org_in, org_out
      e = env
      e%pack    = snow
      e%snowfac = cover
      e%dt      = dt
      call ground_surface_exchange(e, h_ground, le_ground, g_top, ground_rad, snow_energy_in, snow_vapour_out,     &
                                   org_in, org_out)
      call snow_take_exchange(snow, dt * snow_energy_in, dt * snow_vapour_out, e%params)
      energy_gain = dt * (ground_rad - h_ground - le_ground - g_top)
      vapour_loss = dt * snow_vapour_out
   end subroutine snow_step

   subroutine test_snow_cover()
      real(wp), parameter :: GROUND_ROUGH = 0.01_wp   ! [m] the [aerodynamics] soil_rough default
      type(snow_params_t) :: p
      real(wp) :: f0, f1, f2
      print '(a)', 'test_snow_cover:'
      f0 = snow_cover_fraction(0.0_wp,  0.0_wp,  p%rho_snow, p, GROUND_ROUGH)
      f1 = snow_cover_fraction(1.0_wp,  1.0_wp/p%rho_snow,  p%rho_snow, p, GROUND_ROUGH)
      f2 = snow_cover_fraction(50.0_wp, 50.0_wp/p%rho_snow, p%rho_snow, p, GROUND_ROUGH)
      call check('snowfac(0)=0', f0, 0.0_wp, 1.0e-12_wp)
      call check_true('snowfac monotone 0<f1<f2<=1', f1 > 0.0_wp .and. f1 < f2 .and. f2 <= 1.0_wp, f2)
      call check_true('deep pack fully covered (~1)', f2 > 0.99_wp, f2)
      !--- 1 cm of snow on ground of 1 cm roughness covers about a sixth of it (Niu & Yang 2007).
      f1 = snow_cover_fraction(2.5_wp, 2.5_wp/p%rho_snow, p%rho_snow, p, GROUND_ROUGH)
      call check('1 cm of snow covers tanh(0.16) of the ground', f1, tanh(0.16_wp), 1.0e-12_wp)
   end subroutine test_snow_cover

   !--- Build a cold pack by accumulation; assert growth, sub-freezing, and closed ledgers.
   subroutine test_accumulate()
      type(snow_state_t)  :: snow
      type(surface_env_t) :: env
      type(snow_params_t) :: p
      type(snow_out_t)    :: melt
      real(wp) :: dt, snowfall, mass_in, mass_out, e_in, e_out, gain, vapour, g_top
      integer(ik) :: k
      print '(a)', 'test_accumulate:'
      call cold_env(env)
      snow = snow_state_t()
      dt = 900.0_wp
      snowfall = 1.0e-3_wp                                   ! [kg/m2/s] ~ 3.6 mm SWE/hr
      mass_in  = 0.0_wp
      mass_out = 0.0_wp
      e_in     = 0.0_wp
      e_out    = 0.0_wp
      do k = 1_ik, 40_ik
         call snow_accumulate(snow, snowfall, 0.0_wp, 264.0_wp, dt, p)
         mass_in = mass_in + snowfall * dt
         e_in    = e_in    + snowfall * dt * internal_energy_ice(min(t_3ple, 264.0_wp))
         call snow_step(snow, env, 1.0_wp, dt, gain, vapour, g_top)
         mass_out = mass_out + vapour
         e_in     = e_in     + gain
         call snow_drain_meltwater(snow, p, melt)
         mass_out = mass_out + melt%melt_mass + melt%dump_mass
         e_out    = e_out    + melt%melt_enth + melt%dump_enth
      end do
      call check_true('pack grew (swe > 30 kg/m2)', snow%swe(1) > 30.0_wp, snow%swe(1))
      call check_true('snow sub-freezing (temp < t_3ple)', snow%snow_temp(1) < t_3ple, snow%snow_temp(1))
      call check_true('fully frozen (fliq = 0)', snow%snow_fliq(1) < 1.0e-9_wp, snow%snow_fliq(1))
      call check('MASS closure: swe = in - out', snow%swe(1), mass_in - mass_out, 1.0e-9_wp)
      call check('ENERGY closure: energy = in - out', snow%snow_energy(1), e_in - e_out, 1.0e-4_wp)
   end subroutine test_accumulate

   !--- Warm a cold pack: temp pins at t_3ple, fliq rises, melt drains, swe -> melt-out.
   subroutine test_melt_event()
      type(snow_state_t)  :: snow
      type(surface_env_t) :: env
      type(snow_params_t) :: p
      type(snow_out_t)    :: melt
      real(wp) :: dt, swe0, mass_out, e_in, e_out, plateau_temp, tot_melt, gain, vapour, g_top
      integer(ik) :: k
      logical :: saw_plateau
      print '(a)', 'test_melt_event:'
      call cold_env(env)
      snow = cold_pack(20.0_wp, 270.0_wp, p)                 ! a 20 kg/m2 cold pack, all ice at 270 K
      swe0 = snow%swe(1)
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
      e_in = snow%snow_energy(1)                             ! initial store folded into the ledger baseline
      saw_plateau  = .false.
      plateau_temp = 0.0_wp
      do k = 1_ik, 120_ik
         call snow_step(snow, env, 1.0_wp, dt, gain, vapour, g_top)
         mass_out = mass_out + vapour
         e_in     = e_in     + gain
         if (snow%nlayer == 1_ik .and. snow%snow_fliq(1) > 0.05_wp .and. snow%snow_fliq(1) < 0.95_wp) then
            saw_plateau  = .true.
            plateau_temp = snow%snow_temp(1)
         end if
         call snow_drain_meltwater(snow, p, melt)
         mass_out = mass_out + melt%melt_mass + melt%dump_mass
         tot_melt = tot_melt + melt%melt_mass
         e_out    = e_out    + melt%melt_enth + melt%dump_enth
         if (snow%nlayer == 0_ik) exit
      end do
      call check_true('melted out (nlayer = 0)', snow%nlayer == 0_ik, real(snow%nlayer, wp))
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
      type(snow_state_t)  :: snow, probe
      type(surface_env_t) :: env
      type(snow_params_t) :: p
      type(snow_out_t)    :: melt
      real(wp) :: swe0, tsnow0, gain, vapour, g_top
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
      env%k_soil_top  = 1.5_wp
      env%dz_soil_top = 0.05_wp
      env%can_shv     = sat_specific_humidity(268.0_wp, env%press, 0.0_wp)
      swe0   = snow%swe(1)
      tsnow0 = 268.0_wp
      do k = 1_ik, 50_ik
         call snow_step(snow, env, 1.0_wp, 900.0_wp, gain, vapour, g_top)
         call snow_drain_meltwater(snow, p, melt)
      end do
      call check_true('no heat into the soil (no soil gradient)', abs(g_top) < 1.0e-6_wp, g_top)
      call check_true('swe steady (drift < 0.5 kg/m2)', abs(snow%swe(1) - swe0) < 0.5_wp, snow%swe(1) - swe0)
      call check_true('snow_temp steady (within 1 K)', abs(snow%snow_temp(1) - tsnow0) < 1.0_wp, snow%snow_temp(1))
      call check_true('stayed frozen (no spurious melt)', snow%snow_fliq(1) < 1.0e-6_wp, snow%snow_fliq(1))
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
      type(snow_state_t)  :: snow
      type(surface_env_t) :: env
      type(snow_params_t) :: p
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
      env%k_soil_top  = 1.5_wp
      env%dz_soil_top = 0.05_wp
      snow = cold_pack(SWE, T_PACK, p)
      g_expect = (T_PACK - T_SOIL) / (snow%snow_depth(1) / p%k_snow + env%dz_soil_top / env%k_soil_top)
      call snow_step(snow, env, 1.0_wp, 60.0_wp, gain, vapour, g_top)
      call check('the soil loses what the whole pack and its top node conduct', g_top, g_expect,           &
                 1.0e-3_wp * abs(g_expect))
      snow = cold_pack(2.0_wp * SWE, T_PACK, p)
      call snow_step(snow, env, 1.0_wp, 60.0_wp, gain, vapour, g_deep)
      call check_true('a pack twice as deep passes the soil less', abs(g_deep) < abs(g_top), g_deep)
   end subroutine test_insulation

   !--- Rain on a sub-freezing pack: adds mass + liquid enthalpy, refreezes (energy rises).
   subroutine test_rain_on_snow()
      type(snow_state_t)  :: snow
      type(snow_params_t) :: p
      real(wp) :: e0, swe0
      print '(a)', 'test_rain_on_snow:'
      snow = cold_pack(10.0_wp, 265.0_wp, p)
      e0   = snow%snow_energy(1)
      swe0 = snow%swe(1)
      call snow_accumulate(snow, 0.0_wp, 2.0e-3_wp, 278.0_wp, 900.0_wp, p)    ! warm rain (2 mm/900s)
      call check('rain adds mass', snow%swe(1), swe0 + 2.0e-3_wp * 900.0_wp, 1.0e-9_wp)
      call check_true('rain adds energy (warms/refreezes pack)', snow%snow_energy(1) > e0, snow%snow_energy(1) - e0)
      call check('rain enthalpy exact', snow%snow_energy(1) - e0,                                   &
                 2.0e-3_wp * 900.0_wp * internal_energy_liquid(278.0_wp), 1.0e-6_wp)
   end subroutine test_rain_on_snow

   !---~---
   !   A thin pack (0.3 kg/m2, a fifth of the ground) under the bare-ground conductance of an
   ! open patch: an explicit flux at t_3ple, or the nonlinear sublimation at a linearized
   ! temperature, would drive it hundreds of kelvin past its balance.
   !---~---
   subroutine test_stiff_coupling()
      type(snow_state_t)  :: snow
      type(surface_env_t) :: env
      type(snow_params_t) :: p
      type(snow_out_t)    :: melt
      real(wp), parameter :: SWE = 0.3_wp, DT = 900.0_wp, AF = 0.2_wp
      real(wp) :: e0, m_in, m_out, e_in, e_out, t_min, gain, vapour, g_top
      integer(ik) :: k
      print '(a)', 'test_stiff_coupling:'
      env%ggnet       = 0.065_wp
      env%rho_air     = 1.19_wp
      env%press       = 101325.0_wp
      env%t_soil_top  = t_3ple
      env%k_soil_top  = 1.5_wp
      env%dz_soil_top = 0.025_wp

      !--- (a) A cold pack (243 K) under warm, moist canopy air: it warms toward the air, not past it.
      snow = cold_pack(SWE, 243.0_wp, p)
      env%sw_snow  = 0.0_wp
      env%lw_down  = stefan_boltzmann * 243.0_wp ** 4 + 20.0_wp
      env%can_temp = 268.8_wp
      env%can_shv  = 1.79e-3_wp
      e0 = snow%snow_energy(1)
      call snow_step(snow, env, AF, DT, gain, vapour, g_top)
      call check_true('a cold pack under warm air and on warm soil warms, and stays frozen',             &
                      snow%snow_temp(1) > 243.0_wp .and. snow%snow_temp(1) < t_3ple .and.                &
                      snow%snow_fliq(1) == 0.0_wp, snow%snow_temp(1))
      call check('its energy change is what its surroundings gave up', snow%snow_energy(1) - e0, gain,       &
                 1.0e-9_wp * abs(e0))

      !---~---
      !   (b) A melting pack (10 % liquid) under cold, dry canopy air: it freezes through the plateau and
      ! cools toward the air; over a day it stays near the air, and both ledgers close.
      !---~---
      snow = snow_state_t()
      snow%nlayer = 1_ik
      snow%swe(1) = SWE
      snow%snow_energy(1) = temp_to_internal_energy(0.0_wp, SWE, t_3ple, 0.1_wp)
      snow%snow_depth(1)  = SWE / p%rho_snow
      snow%density(1)     = p%rho_snow
      env%lw_down  = stefan_boltzmann * t_3ple ** 4 - 60.0_wp
      env%can_temp = 255.0_wp
      env%can_shv  = 3.0e-4_wp
      m_in  = SWE
      m_out = 0.0_wp
      e_in  = snow%snow_energy(1)
      e_out = 0.0_wp
      t_min = t_3ple
      do k = 1_ik, 96_ik
         call snow_step(snow, env, AF, DT, gain, vapour, g_top)
         m_out = m_out + vapour
         e_in  = e_in  + gain
         if (k == 1_ik) call check_true('one step freezes the liquid without passing the air temperature', &
                                        snow%snow_fliq(1) == 0.0_wp .and. snow%snow_temp(1) > env%can_temp - 5.0_wp, &
                                        snow%snow_temp(1))
         call snow_drain_meltwater(snow, p, melt)
         m_out = m_out + melt%melt_mass + melt%dump_mass
         e_out = e_out + melt%melt_enth + melt%dump_enth
         if (snow%nlayer == 0_ik) exit
         t_min = min(t_min, snow%snow_temp(1))
      end do
      call check_true('over a day the pack never falls more than 5 K below the air', t_min > env%can_temp - 5.0_wp, &
                      t_min)
      call check('MASS closure', snow%swe(1), m_in - m_out, 1.0e-12_wp)
      call check('ENERGY closure', snow%snow_energy(1), e_in - e_out, 1.0e-6_wp * abs(e0))
   end subroutine test_stiff_coupling

   !--- 7. Each tile's radiation at its own optics; the cover-weighted sum is the blended ground's.
   subroutine test_tile_radiation()
      type(snow_state_t)  :: snow
      type(snow_params_t) :: p
      real(wp), parameter :: DN_VIS = 150.0_wp, DN_NIR = 150.0_wp, DN_LW = 280.0_wp   ! [W/m2] reaching the ground
      real(wp), parameter :: SOIL_VIS = 0.10_wp, SOIL_NIR = 0.20_wp, SOIL_EMISS = 0.95_wp, T_SOIL = 275.0_wp
      real(wp) :: fc, a_vis, a_nir, emiss, temp, rnet_soil, sw_snow, rnet_snow, blended
      print '(a)', 'test_tile_radiation:'
      snow = cold_pack(8.0_wp, 263.0_wp, p)                  ! 3.2 cm deep: cover ~0.47
      call internal_energy_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      fc = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), snow%density(1), p, 0.01_wp)
      call check_true('a patchy pack', fc > 0.05_wp .and. fc < 0.95_wp, fc)
      call ground_tile_radiation(DN_VIS, DN_NIR, DN_LW, SOIL_VIS, SOIL_NIR, SOIL_EMISS, T_SOIL, snow, p,          &
                                 rnet_soil, sw_snow)
      rnet_snow = sw_snow + p%snow_emiss * (DN_LW - stefan_boltzmann * snow%snow_temp(1) ** 4)
      a_vis = SOIL_VIS
      a_nir = SOIL_NIR
      emiss = SOIL_EMISS
      temp  = T_SOIL
      call snow_blend_ground_optics(snow, p, 0.01_wp, a_vis, a_nir, emiss, temp)
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
      real(wp) :: cp_cas, h_bare, le_soil, h_2k, h_ground, le_ground, g_top, ground_rad, energy_in, vapour_out
      real(wp) :: org_in, org_out
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
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, energy_in, vapour_out,           &
                                   org_in, org_out)
      call check('the ground surface''s exchange takes the canopy air''s specific heat', h_ground, h_bare,        &
                 1.0e-12_wp * h_bare)
   end subroutine test_ground_balance

   !--- 9. The pond and the water layer, before a step.
   subroutine test_pond_stage()
      real(wp), parameter :: DT = 900.0_wp
      type(surface_env_t) :: env
      type(snow_params_t) :: p
      type(snow_out_t)    :: st
      type(snow_state_t)  :: pack0
      real(wp) :: w, e, w0, e0, hold, mix
      print '(a)', 'test_pond_stage:'

      !--- (a) A pond half frozen at the triple point, alone: it starts a layer of ice.
      call cold_env(env)
      env%params = p
      env%dt     = DT
      w0 = 3.0_wp
      e0 = temp_to_internal_energy(0.0_wp, w0, t_3ple, 0.5_wp)
      w  = w0
      e  = e0
      call prepare_snow_stage(0.0_wp, 0.0_wp, 265.0_wp, 0.01_wp, env, st, w, e)
      hold = 0.5_wp * w0                                         ! the pond's liquid
      call check_true('a freezing pond starts a layer', st%exists .and. env%pack%nlayer == 1_ik)
      call check('the layer is ice', env%pack%density(1), rho_ice, 0.0_wp)
      call check('the pond keeps its liquid', w, hold, 1.0e-12_wp)
      call check('the layer is the pond''s ice', env%pack%swe(1), w0 - hold, 1.0e-12_wp)
      call check('MASS: layer + pond + drained = the pond', env%pack%swe(1) + w + st%melt_mass + st%dump_mass, w0, &
                 1.0e-12_wp)
      call check('ENERGY: layer + pond + drained = the pond', env%pack%snow_energy(1) + e + st%melt_enth          &
                 + st%dump_enth, e0, 1.0e-9_wp * abs(e0))
      call check('the ledger starts from the stores the stage began with', st%swe0 + w, w0, 1.0e-12_wp)
      call check('the layer is as deep as its water as ice', env%pack%snow_depth(1), env%pack%swe(1) / rho_ice,  &
                 1.0e-15_wp)

      !--- (b) A cold deep pack and a pond at the triple point: the pack freezes it.
      call cold_env(env)
      env%params = p
      env%dt     = DT
      env%pack   = cold_pack(40.0_wp, 255.0_wp, p)
      pack0      = env%pack
      w0 = 2.0_wp
      e0 = w0 * internal_energy_liquid(t_3ple)
      w  = w0
      e  = e0
      call prepare_snow_stage(0.0_wp, 0.0_wp, 255.0_wp, 0.01_wp, env, st, w, e)
      call check_true('the pond froze into the pack', w == 0.0_wp .and. e == 0.0_wp .and. st%melt_mass == 0.0_wp)
      call check('the pack holds both', env%pack%swe(1), 42.0_wp, 1.0e-12_wp)
      call check('the pack holds both energies', env%pack%snow_energy(1), pack0%snow_energy(1) + e0,              &
                 1.0e-9_wp * abs(e0))
      mix = 42.0_wp / (40.0_wp / p%rho_snow + 2.0_wp / rho_ice)
      call check('the pond adds its depth as ice', env%pack%density(1), mix, 1.0e-9_wp)
      call check('the ledger''s pack takes the pond', st%swe0, 42.0_wp, 1.0e-12_wp)
      call check('the ledger''s pack energy takes the pond''s', st%enth0, pack0%snow_energy(1) + e0,              &
                 1.0e-9_wp * abs(e0))
      call check_true('the frozen pond warmed the pack', env%pack%snow_temp(1) > 255.0_wp                         &
                      .and. env%pack%snow_fliq(1) == 0.0_wp, env%pack%snow_temp(1))

      !--- (b2) A pack at the plateau, holding little liquid, and a pond: it takes none of the pond.
      call cold_env(env)
      env%params = p
      env%dt     = DT
      env%pack   = cold_pack(40.0_wp, t_3ple, p)
      env%pack%snow_energy(1) = temp_to_internal_energy(0.0_wp, 40.0_wp, t_3ple, 0.02_wp)   ! 2 % liquid
      w0 = 2.0_wp
      e0 = w0 * internal_energy_liquid(t_3ple)
      w  = w0
      e  = e0
      call prepare_snow_stage(0.0_wp, 0.0_wp, t_3ple, 0.01_wp, env, st, w, e)
      call check('a pack at the plateau leaves the pond its water', w, w0, 1.0e-12_wp)
      call check('and keeps its own', env%pack%swe(1), 40.0_wp, 1.0e-12_wp)
      call check_true('and its density', env%pack%density(1) == p%rho_snow)

      !--- (c) A warm pond and no pack: left as it is.
      call cold_env(env)
      env%params = p
      env%dt     = DT
      w0 = 2.0_wp
      e0 = w0 * internal_energy_liquid(285.0_wp)
      w  = w0
      e  = e0
      call prepare_snow_stage(0.0_wp, 0.0_wp, 285.0_wp, 0.01_wp, env, st, w, e)
      call check_true('a warm pond with no pack is left as it is', .not. st%exists .and. w == w0 .and. e == e0     &
                      .and. env%pack%nlayer == 0_ik)

      !--- (d) A pack and an empty pond, snowing: the pack as before, fresh snow on fresh snow.
      call cold_env(env)
      env%params = p
      env%dt     = DT
      env%pack   = cold_pack(10.0_wp, 265.0_wp, p)
      w = 0.0_wp
      e = 0.0_wp
      call prepare_snow_stage(1.0e-3_wp, 0.0_wp, 265.0_wp, 0.01_wp, env, st, w, e)
      call check_true('fresh snow on fresh snow keeps its density', env%pack%density(1) == p%rho_snow)
      call check_true('an empty pond stays empty, and the ledger''s pack is the pack', w == 0.0_wp .and.          &
                      st%swe0 == 10.0_wp .and. st%pond_mass == 0.0_wp)

      !--- (e) Snow on a layer of ice: the densities mix by volume.
      call cold_env(env)
      env%params = p
      env%dt     = DT
      env%pack   = cold_pack(3.0_wp, 265.0_wp, p)
      env%pack%density(1)    = rho_ice
      env%pack%snow_depth(1) = 3.0_wp / rho_ice
      w = 0.0_wp
      e = 0.0_wp
      call prepare_snow_stage(1.0e-3_wp, 0.0_wp, 265.0_wp, 0.01_wp, env, st, w, e)
      mix = 3.9_wp / (3.0_wp / rho_ice + 0.9_wp / p%rho_snow)
      call check('snow on ice mixes the densities by volume', env%pack%density(1), mix, 1.0e-9_wp)
   end subroutine test_pond_stage

   !--- 10a. The organic layer's slots: 1.35 kg of litter (3 cm) over 2.15 kg of moss (5 cm), dry, at 280 K.
   function two_slots(o) result(org)
      type(organic_params_t), intent(in) :: o
      type(organic_state_t) :: org
      org = organic_state_t()
      org%density = [o%litter%bulk_density, o%moss%bulk_density]
      org%mass    = [0.03_wp * o%litter%bulk_density, 0.05_wp * o%moss%bulk_density]
      call organic_split_energy(org, temp_to_internal_energy(cp_organic * sum(org%mass), 0.0_wp, 280.0_wp, 1.0_wp))
   end function two_slots

   subroutine test_organic_slots()
      real(wp), parameter :: HEAD_10KPA = -1.0e4_wp / (1000.0_wp * 9.80665_wp)   ! [m]
      type(organic_params_t)   :: o
      type(organic_state_t)    :: org
      type(soil_layer_curve_t) :: c
      real(wp) :: cap_litter, cap_moss, dry_cap, w, e, w_max, resist, e_in, dry_depth, fill, w_out, e_out, t_in
      print '(a)', 'test_organic_slots:'
      org = two_slots(o)
      cap_litter = o%litter%theta_sat * 1000.0_wp * 0.03_wp            ! 29.4 kg/m2
      cap_moss   = o%moss%theta_sat   * 1000.0_wp * 0.05_wp            ! 42.5 kg/m2
      call check('a slot holds its theta_sat over its depth', organic_slot_capacity(o, org, SLOT_MOSS), cap_moss,      &
                 1.0e-12_wp)
      org%water = [cap_litter, cap_moss]
      call organic_split_energy(org, temp_to_internal_energy(cp_organic * sum(org%mass), sum(org%water), 283.0_wp, &
                                                             1.0_wp))
      call check('one temperature, the total''s', org%temp, 283.0_wp, 1.0e-9_wp)
      call check('the slots'' energies sum to the total', sum(org%energy),                                          &
                 temp_to_internal_energy(cp_organic * sum(org%mass), sum(org%water), 283.0_wp, 1.0_wp), 1.0e-6_wp)
      call organic_layer_sums(org, o, dry_cap, w, e, w_max, resist, dry_depth)
      call check('it holds the slots'' saturation', w_max, cap_litter + cap_moss, 1.0e-9_wp)
      call check('a saturated layer has no dried top', dry_depth, 0.0_wp, 1.0e-12_wp)
      call check('its resistance: each slot''s depth over its saturated conductivity', resist,                     &
                 0.03_wp / o%litter%k_wet + 0.05_wp / o%moss%k_wet, 1.0e-12_wp)
      e_in = sum(org%energy)
      call organic_take_exchange(org, -1.0e5_wp, [cap_litter, 2.0_wp])
      call check_true('each slot gives the vapour taken from it', org%water(SLOT_LITTER) == 0.0_wp                  &
                      .and. abs(org%water(SLOT_MOSS) - (cap_moss - 2.0_wp)) < 1.0e-9_wp)
      call check('the exchange''s energy reaches the layer', sum(org%energy), e_in - 1.0e5_wp, 1.0e-6_wp)
      call organic_layer_sums(org, o, dry_cap, w, e, w_max, resist, dry_depth)
      call check('its dried top: the empty litter, and the share of the moss its water does not fill', dry_depth, &
                 0.03_wp + 0.05_wp * 2.0_wp / cap_moss, 1.0e-12_wp)
      fill = 1.0_wp - 2.0_wp / cap_moss
      call check('its resistance: each dried top at the dry conductivity, the wet rest at the wet one', resist,    &
                 0.03_wp / o%litter%k_dry + 0.05_wp * ((1.0_wp - fill) / o%moss%k_dry + fill / o%moss%k_wet),        &
                 1.0e-12_wp)
      !--- Each slot's curve: its material's at its depth, the residual no more than the liquid it holds.
      c = organic_slot_curve(o, org, SLOT_MOSS)
      call check('the moss''s curve: its depth', c%dz, 0.05_wp, 1.0e-12_wp)
      call check('the moss at -10 kPa holds 0.29 (Zhang et al. 2010)',                                             &
                 soil_head_to_theta(c%retention, HEAD_10KPA, c%theta_sat, c%theta_res, c%par_a, c%par_n), 0.294_wp,  &
                 1.0e-3_wp)
      call check('a moss above its residual keeps it', c%theta_res, o%moss%theta_res, 0.0_wp)
      c = organic_slot_curve(o, org, SLOT_LITTER)
      call check('an empty litter starts 1 % of the way to saturation, where it can take water',                  &
                 (0.0_wp - c%theta_res) / (c%theta_sat - c%theta_res), 1.0e-2_wp, 1.0e-12_wp)
      org%water(SLOT_LITTER) = cap_litter
      c = organic_slot_curve(o, org, SLOT_LITTER)
      call check('the litter at -10 kPa holds 0.14, near drained litter''s 0.12 (Putuhena & Cordery 1996)',        &
                 soil_head_to_theta(c%retention, HEAD_10KPA, c%theta_sat, c%theta_res, c%par_a, c%par_n), 0.1415_wp, &
                 1.0e-3_wp)
      org%fliq = 0.0_wp
      c = organic_slot_curve(o, org, SLOT_MOSS)
      call check('a frozen slot''s liquid, none, starts 1 % of the way to saturation',                            &
                 (0.0_wp - c%theta_res) / (c%theta_sat - c%theta_res), 1.0e-2_wp, 1.0e-12_wp)
      !--- Frost on a frozen, saturated litter: what it holds above saturation spills, with its energy.
      org%water = [cap_litter + 0.5_wp, cap_moss]
      call organic_split_energy(org, temp_to_internal_energy(cp_organic * sum(org%mass), sum(org%water), 265.0_wp, &
                                                             0.0_wp))
      e_in = sum(org%energy)
      t_in = org%temp
      call organic_spill(org, o, w_out, e_out)
      call check('above saturation: the excess spills', w_out, 0.5_wp, 1.0e-12_wp)
      call check('above saturation: the slot keeps its saturation', org%water(SLOT_LITTER), cap_litter, 1.0e-9_wp)
      call check('above saturation: the spill takes its ice''s energy, the rest stays', sum(org%energy) + e_out,    &
                 e_in, 1.0e-6_wp)
      call check('above saturation: the layer keeps its temperature', org%temp, t_in, 1.0e-9_wp)
      call organic_spill(org, o, w_out, e_out)
      call check('a saturated layer spills nothing', w_out, 0.0_wp, 0.0_wp)
   end subroutine test_organic_slots

   !--- 10b. The organic layer's balance, and the ground surface with it.
   subroutine test_organic_balance()
      real(wp), parameter :: DT = 3600.0_wp
      type(organic_params_t) :: o
      type(organic_state_t)  :: org
      type(surface_env_t)    :: env
      type(organic_out_t)    :: ofx, ofx0
      type(snow_params_t)    :: p
      real(wp) :: w_max, h_ground, le_ground, g_top, ground_rad, s_in, s_out, o_in, o_out, w_save, e_save
      real(wp) :: dry_depth, r_save
      print '(a)', 'test_organic_balance:'
      org = organic_state_t()                                    ! litter alone: 1.35 kg/m2, 2.6 kJ/m2/K dry
      org%density(SLOT_LITTER) = o%litter%bulk_density
      org%mass(SLOT_LITTER)    = 0.03_wp * o%litter%bulk_density
      org%water(SLOT_LITTER)   = 1.0_wp
      call organic_split_energy(org, temp_to_internal_energy(cp_organic * org%mass(SLOT_LITTER), 1.0_wp, 280.0_wp, &
                                                             1.0_wp))
      call cold_env(env)
      env%dt           = DT
      env%can_temp     = 270.0_wp
      env%t_soil_top   = 285.0_wp
      env%ggnet        = 0.05_wp
      env%rnet_soil    = 0.0_wp
      env%has_org      = .true.
      env%ground_emiss = 0.95_wp
      call organic_layer_sums(org, o, env%org_dry_cap, env%org_water, env%org_energy, w_max, env%org_resistance, &
                              dry_depth)
      env%org_dry_resist = organic_cover_path(dry_depth, 0.0_wp, o%litter%theta_sat) / vapour_diffusivity(280.0_wp)
      call organic_surface_balance(env, 0.0_wp, ofx0)
      call check_true('a thin layer relaxes over an hour without passing its surroundings',                         &
                      ofx0%t_surf > 268.0_wp .and. ofx0%t_surf < 286.0_wp, ofx0%t_surf)
      call check('its energy: what its surroundings give it', ofx0%energy_in,                                       &
                 ofx0%rad - ofx0%h - ofx0%le - ofx0%g_soil, 1.0e-9_wp)
      call check_true('a wet layer evaporates, no more than it holds',                                              &
                      ofx0%w_flux > 0.0_wp .and. ofx0%w_flux * DT <= env%org_water, ofx0%w_flux)
      !---~---
      !   What the layer takes is its balance's root: its energy and water after the exchange give
      ! back the temperature the balance solved for, on a wet layer and on a trace of water whose
      ! evaporation the layer's water caps (a few joules per kelvin turn a mismatch into kelvins).
      !---~---
      call check('the exchange it takes gives back its root', committed_temp(env, ofx0), ofx0%t_surf, 1.0e-6_wp)
      w_save = env%org_water
      e_save = env%org_energy
      env%org_water  = 1.0e-6_wp
      env%org_energy = temp_to_internal_energy(env%org_dry_cap, env%org_water, 280.0_wp, 1.0_wp)
      call organic_surface_balance(env, 0.0_wp, ofx)
      call check('a trace of water: the layer still takes its root', committed_temp(env, ofx), ofx%t_surf, 1.0e-6_wp)
      call check_true('and evaporates no more than it holds', ofx%w_flux * DT <= env%org_water * (1.0_wp + 1.0e-12_wp), &
                      ofx%w_flux)
      env%org_water  = 0.0_wp
      env%org_energy = temp_to_internal_energy(env%org_dry_cap, 0.0_wp, 280.0_wp, 1.0_wp)
      call organic_surface_balance(env, 0.0_wp, ofx)
      call check('an empty layer does not evaporate', ofx%w_flux, 0.0_wp, 0.0_wp)
      !--- Its water's vapour crosses its dried top: the same layer full to its top exchanges as open water.
      env%org_water  = w_save
      env%org_energy = e_save
      r_save = env%org_dry_resist
      env%org_dry_resist = 0.0_wp
      call organic_surface_balance(env, 0.0_wp, ofx)
      call check_true('a dried top slows its evaporation', ofx0%w_flux > 0.0_wp .and. ofx0%w_flux < 0.5_wp * ofx%w_flux, &
                      ofx0%w_flux / ofx%w_flux)
      env%org_dry_resist = r_save
      !---~---
      !   Held at the water network's rate, evaporation or dew: the layer gives that rate, its energy
      ! closes, and it takes its root.
      !---~---
      env%org_held      = .true.
      env%org_evap_held = 1.0e-5_wp
      call organic_surface_balance(env, 0.0_wp, ofx)
      call check('held: it gives the network''s evaporation', ofx%w_flux, 1.0e-5_wp, 1.0e-18_wp)
      call check('held: its energy, what its surroundings give it', ofx%energy_in,                                 &
                 ofx%rad - ofx%h - ofx%le - ofx%g_soil, 1.0e-9_wp)
      call check('held: the exchange it takes gives back its root', committed_temp(env, ofx), ofx%t_surf, 1.0e-6_wp)
      env%org_evap_held = -2.0e-6_wp
      call organic_surface_balance(env, 0.0_wp, ofx)
      call check('held: dew at the network''s rate', ofx%w_flux, -2.0e-6_wp, 1.0e-18_wp)
      call check_true('held: dew gives the layer its latent heat', ofx%le < 0.0_wp, ofx%le)
      call check('held: dew, the exchange it takes gives back its root', committed_temp(env, ofx), ofx%t_surf,     &
                 1.0e-6_wp)
      env%org_held      = .false.
      env%org_evap_held = 0.0_wp
      !--- The ground with the layer: the soil takes its conduction less the soil's own evaporation's heat.
      env%soil_evap = 1.0e-6_wp
      env%snowfac   = 0.0_wp
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, s_in, s_out, o_in, o_out)
      call organic_surface_balance(env, 0.0_wp, ofx)
      call check('the soil takes what the layer conducts, less its evaporation''s heat', g_top,                    &
                 ofx%g_soil - env%soil_evap * enthalpy_vapor(env%t_soil_top), 1.0e-9_wp)
      call check('the layer is the ground''s surface: its sensible heat, the ground''s', h_ground, ofx%h, 1.0e-12_wp)
      call check('no pack, no pack exchange', s_in, 0.0_wp, 0.0_wp)
      call check('the layer''s own energy', o_in, ofx%energy_in, 1.0e-12_wp)
      !--- A pack on the layer: it conducts into the layer, and what the ground absorbs is what both keep.
      env%params  = p
      env%pack    = cold_pack(30.0_wp, 265.0_wp, p)
      env%snowfac = 0.9_wp
      call ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, s_in, s_out, o_in, o_out)
      call check('pack and layer: the ground''s energy closes', ground_rad - h_ground - le_ground - g_top,          &
                 s_in + o_in, 1.0e-9_wp)
   end subroutine test_organic_balance

   !--- The layer's temperature after it takes an exchange over the balance's step.
   real(wp) function committed_temp(env, ofx) result(t)
      type(surface_env_t), intent(in) :: env
      type(organic_out_t), intent(in) :: ofx
      real(wp) :: fl
      call internal_energy_to_temp(env%org_energy + env%dt * ofx%energy_in, env%org_water - env%dt * ofx%w_flux,       &
                                   env%org_dry_cap, t, fl)
   end function committed_temp

end program test_surface
