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
!      air, in proportion to the difference.
! The energy datum the pack is held in is test_therm_lib's.
!---~---
program test_surface
   use ed_test_assert,        only : check, check_true, test_report
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : t_3ple, stefan_boltzmann
   use ed_therm_lib,          only : temp_to_internal_energy, internal_energy_to_temp, internal_energy_ice,      &
                                     internal_energy_liquid, sat_specific_humidity
   use ed_patch_state_vars,   only : snow_state_t
   use ed_biophysics_params, only : snow_params_t
   use ed_surface,            only : surface_env_t, snow_out_t, snow_cover_fraction, snow_accumulate,            &
                                     snow_drain_meltwater, ground_surface_exchange, snow_take_exchange,         &
                                     snow_blend_ground_optics, ground_tile_radiation, ground_surface_fluxes
   implicit none

   call test_snow_cover()
   call test_accumulate()
   call test_melt_event()
   call test_isothermal()
   call test_rain_on_snow()
   call test_stiff_coupling()
   call test_tile_radiation()
   call test_ground_balance()

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
      real(wp),            intent(in) :: swe, temp   !--- [kg/m2], [K]
      type(snow_params_t), intent(in) :: p
      type(snow_state_t) :: snow
      snow = snow_state_t()
      snow%nlayer = 1_ik
      snow%swe(1) = swe
      snow%snow_energy(1) = temp_to_internal_energy(0.0_wp, swe, temp, 0.0_wp)
      snow%snow_depth(1)  = swe / p%rho_snow
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
      type(snow_state_t),  intent(inout) :: snow          !--- the pack
      type(surface_env_t), intent(in)    :: env           !--- its surroundings; the pack and cover are set here
      real(wp),            intent(in)    :: cover, dt     !--- [-], [s]
      real(wp),            intent(out)   :: energy_gain, vapour_loss, g_top   !--- [J/m2], [kg/m2], [W/m2]
      type(surface_env_t) :: e
      real(wp) :: h_ground, le_ground, ground_rad, snow_energy_in, snow_vapour_out
      e = env
      e%pack    = snow
      e%snowfac = cover
      e%dt      = dt
      call ground_surface_exchange(e, h_ground, le_ground, g_top, ground_rad, snow_energy_in, snow_vapour_out)
      call snow_take_exchange(snow, dt * snow_energy_in, dt * snow_vapour_out, e%params)
      energy_gain = dt * (ground_rad - h_ground - le_ground - g_top)
      vapour_loss = dt * snow_vapour_out
   end subroutine snow_step

   subroutine test_snow_cover()
      real(wp), parameter :: GROUND_ROUGH = 0.01_wp   !--- [m] the [aerodynamics] soil_rough default
      type(snow_params_t) :: p
      real(wp) :: f0, f1, f2
      print '(a)', 'test_snow_cover:'
      f0 = snow_cover_fraction(0.0_wp,  0.0_wp,  p, GROUND_ROUGH)
      f1 = snow_cover_fraction(1.0_wp,  1.0_wp/p%rho_snow,  p, GROUND_ROUGH)
      f2 = snow_cover_fraction(50.0_wp, 50.0_wp/p%rho_snow, p, GROUND_ROUGH)
      call check('snowfac(0)=0', f0, 0.0_wp, 1.0e-12_wp)
      call check_true('snowfac monotone 0<f1<f2<=1', f1 > 0.0_wp .and. f1 < f2 .and. f2 <= 1.0_wp, f2)
      call check_true('deep pack fully covered (~1)', f2 > 0.99_wp, f2)
      !--- 1 cm of snow on ground of 1 cm roughness covers about a sixth of it (Niu & Yang 2007).
      f1 = snow_cover_fraction(2.5_wp, 2.5_wp/p%rho_snow, p, GROUND_ROUGH)
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
      fc = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), p, 0.01_wp)
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

   !--- 8. The bare ground's skin fluxes to the canopy air.
   subroutine test_ground_balance()
      real(wp), parameter :: GGNET = 0.02_wp, RHO = 1.2_wp, T_GROUND = 300.0_wp
      real(wp) :: h_bare, le_soil, h_2k
      print '(a)', 'test_ground_balance:'
      call ground_surface_fluxes(T_GROUND, T_GROUND, GGNET, RHO, 0.0_wp, h_bare, le_soil)
      call check_true('H_bare = 0 (no gradient)', abs(h_bare) < 1.0e-9_wp, h_bare)
      call check_true('LE_soil = 0 (no evap)',    abs(le_soil) < 1.0e-30_wp, le_soil)
      call ground_surface_fluxes(T_GROUND, T_GROUND - 1.0_wp, GGNET, RHO, 0.0_wp, h_bare, le_soil)
      call ground_surface_fluxes(T_GROUND, T_GROUND - 2.0_wp, GGNET, RHO, 0.0_wp, h_2k, le_soil)
      call check_true('warm ground over cooler canopy air heats it', h_bare > 0.0_wp, h_bare)
      call check('twice the difference, twice the sensible flux', h_2k, 2.0_wp * h_bare, 1.0e-12_wp * h_2k)
   end subroutine test_ground_balance

end program test_surface
