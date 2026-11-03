! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_therm_lib -- the thermodynamics the stores share (ed_therm_lib).
!   1. TEMPERATURE RESPONSE (arrhenius_scale, peaked_arrhenius_scale,
!      kattge_knorr_entropy; docs/science/plant_respiration.md eq. (5),
!      docs/science/leaf_gas_exchange.md eq. (4)): each form is k25 at 25 C; Arrhenius
!      rises with temperature and the peaked form has an interior optimum; a warm-grown
!      plant's entropy shift keeps more capacity in heat and moves the optimum up.
!   2. DATUM (internal_energy_to_temp, temp_to_internal_energy, internal_energy_ice): the
!      phase-change inverter round-trips liquid and ice, is continuous at the ends of the
!      melt plateau and sits at t_3ple on it; ice holds internal_energy_ice per kilogram,
!      and melting it at t_3ple costs the latent heat of fusion
!      (docs/science/coupled_biophysics.md, principle 1).
!   3. CLAUSIUS-CLAPEYRON (sat_vapor_pressure_temp_deriv): the slope of e_sat against a
!      finite difference.
!   4. ICE SATURATION (sat_vapor_pressure; docs/science/snow_biophysics.md eq. (1)): an
!      absent liquid fraction is the liquid curve to the bit; the branches meet at the
!      triple point; over ice e_sat is Murphy & Koop (2005)'s; at 268 K the air saturated
!      over ice holds 3-8 % less vapour than over liquid.
!   5. SOIL CONDUCTIVITY (soil_thermal_conductivity; docs/science/soil_biophysics.md eq.
!      (11)): frozen soil conducts more than thawed; the Kersten number is Johansen's log
!      form thawed and the linear form frozen.
!   6. SOIL TEXTURE (soil_texture_to_thermal; soil_biophysics.md eq. (9)): CLM5's solids'
!      and dry soil's conductivities and the solids' heat capacity, computed apart.
!   7. MOIST AIR (sat_specific_humidity, rh_to_specific_humidity, specific_humidity_to_vpd,
!      precip_phase; docs/science/forcing.md): the humidity of a dewpoint is
!      eps e / (p - (1 - eps) e) and RH = 1 at the dewpoint gives it again; the VPD is read
!      at the air's pressure: air half saturated at 70 kPa has half e_sat of it there, and a
!      humidity at 70 kPa holds 0.7 of the vapour pressure it holds at 1 atm; the
!      rain/snow split keeps the total, all rain when warm, all snow when cold, part of each
!      at the triple point.
!   8. THE AIR AT ANOTHER HEIGHT (lapse_air_temperature, lapse_pressure,
!      lapse_specific_humidity; forcing.md §8): the temperature falls at the rate given, the
!      pressure drops, barometrically when isothermal, the humidity keeps its relative
!      humidity, and nothing changes at dz = 0.
!   9. A TISSUE'S ENERGY (tissue_temp_to_energy, tissue_energy_to_temp;
!      docs/science/vegetation_energy_dynamics.md eqs. (9)-(10)): it is
!      temp_to_internal_energy with liquid water, tissue_energy_to_temp inverts it, and
!      water that comes or goes at the tissue's temperature leaves the temperature there.
!---~---
program test_therm_lib
   use ed_test_assert,   only : check, check_true, test_report
   use ed_kinds,         only : wp, ik
   use ed_constants,     only : t_kelvin, t_3ple, cp_ice, latent_heat_fusion, k_water, k_ice, h2o_dry_ratio,        &
                                latent_heat_vap_3ple
   use ed_therm_lib,     only : arrhenius_scale, peaked_arrhenius_scale, kattge_knorr_entropy,                    &
                                temp_to_internal_energy, internal_energy_to_temp, internal_energy_ice,            &
                                sat_vapor_pressure, sat_vapor_pressure_temp_deriv, sat_specific_humidity,         &
                                soil_thermal_conductivity, soil_texture_to_thermal, rh_to_specific_humidity,      &
                                precip_phase, lapse_air_temperature, lapse_pressure, lapse_specific_humidity,     &
                                specific_humidity_to_vpd, tissue_temp_to_energy, tissue_energy_to_temp,          &
                                internal_energy_liquid, latent_heat_vap, enthalpy_vapor
   implicit none
   real(wp), parameter :: SCAN_LO = 5.0_wp, SCAN_HI = 55.0_wp   ! [C] the optimum scan's range

   call test_temperature_response()
   call test_datum()
   call test_latent_heat_vap()
   call test_clausius()
   call test_ice_saturation()
   call test_soil_conductivity()
   call test_soil_texture()
   call test_moist_air()
   call test_lapse()
   call test_tissue_energy()

   call test_report('test_therm_lib')

contains

   !---~---
   !   1. The temperature response of enzyme rates, at Vcmax's activation, deactivation and
   ! entropy (60 at 25 C; Ea 65330 J/mol, Hd 200000 J/mol, dS 650 J/mol/K). Kattge & Knorr
   ! (2007)'s acclimation shifts dS by the growth temperature; what it is for is moving where
   ! the response peaks, measured here from the function by a scan.
   !---~---
   subroutine test_temperature_response()
      real(wp), parameter :: VC25 = 60.0_wp, EA = 65330.0_wp, HD = 200000.0_wp, DS = 650.0_wp
      real(wp), parameter :: T_COLD = t_kelvin + 10.0_wp, T_WARM = t_kelvin + 30.0_wp
      real(wp) :: k25, k32, k45, ds_cold, ds_warm, k_cold, k_warm, topt_cold, topt_warm
      print '(a)', 'test_temperature_response:'
      call check('Arrhenius is k25 at 25 C', arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 25.0_wp), 40.49_wp,  &
                 1.0e-12_wp * 40.49_wp)
      call check_true('Arrhenius increases with temperature',                                                  &
                      arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 5.0_wp) <                              &
                      arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 35.0_wp))
      k25 = peaked_arrhenius_scale(VC25, EA, HD, DS, t_kelvin + 25.0_wp)
      k32 = peaked_arrhenius_scale(VC25, EA, HD, DS, t_kelvin + 32.0_wp)
      k45 = peaked_arrhenius_scale(VC25, EA, HD, DS, t_kelvin + 45.0_wp)
      call check('peaked is k25 at 25 C', k25, VC25, 1.0e-12_wp * VC25)
      call check_true('peaked has an interior thermal optimum', k32 > k25 .and. k32 > k45, k32)
      !---~---
      !   At a hot leaf temperature the warm-acclimated plant keeps more capacity than the
      ! cold-acclimated one, and its optimum sits higher.
      !---~---
      ds_cold = kattge_knorr_entropy(668.39_wp, 1.07_wp, T_COLD)
      ds_warm = kattge_knorr_entropy(668.39_wp, 1.07_wp, T_WARM)
      k_cold  = peaked_arrhenius_scale(VC25, EA, HD, ds_cold, t_kelvin + 38.0_wp)
      k_warm  = peaked_arrhenius_scale(VC25, EA, HD, ds_warm, t_kelvin + 38.0_wp)
      call check_true('warm acclimation retains more Vcmax at 38 C than cold acclimation', k_warm > k_cold,   &
                      k_warm - k_cold)
      topt_cold = peak_temperature(ds_cold)
      topt_warm = peak_temperature(ds_warm)
      call check_true('the thermal optimum moves up with growth temperature (> 1 K over 20 K)',                &
                      topt_warm > topt_cold + 1.0_wp, topt_warm - topt_cold)
      call check_true('both optima lie inside the scan, not at its ends',                                       &
                      min(topt_cold, topt_warm) > t_kelvin + SCAN_LO .and. max(topt_cold, topt_warm) < t_kelvin + SCAN_HI, &
                      topt_warm - t_kelvin)
   end subroutine test_temperature_response

   !---~---
   !   The peak of the peaked response [K] by a scan at 0.05 K over [SCAN_LO, SCAN_HI] C, so the
   ! optimum is measured from the function rather than derived from a formula the test would
   ! check against itself; the scan reaches well past a warm-grown stand's optimum (38 C).
   !---~---
   pure real(wp) function peak_temperature(ds) result(t_opt)
      real(wp), intent(in) :: ds   ! [J/mol/K] the entropy term
      real(wp)    :: t, k, k_best
      integer(ik) :: j
      k_best = -1.0_wp
      t_opt  = 0.0_wp
      do j = 0_ik, nint((SCAN_HI - SCAN_LO) / 0.05_wp, ik)
         t = t_kelvin + SCAN_LO + 0.05_wp * real(j, wp)
         k = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, ds, t)
         if (k > k_best) then
            k_best = k
            t_opt  = t
         end if
      end do
   end function peak_temperature

   !---~---
   !   2. The phase-change inverter, per m3 of soil with 1000 J/m3/K of dry matter and 0.3 kg of
   ! water: liquid and ice round-trip; at both ends of the melt plateau it reads t_3ple, and
   ! half way along it t_3ple and half liquid. And 5 kg of snow, with no dry matter: all ice,
   ! it holds internal_energy_ice per kilogram, and melting it at t_3ple costs 5 x 3.34e5 J.
   !---~---
   subroutine test_datum()
      real(wp), parameter :: DH = 1000.0_wp, WM = 0.3_wp, SWE = 5.0_wp
      real(wp) :: uext, t, fl, u_freeze, u_melt, u_ice, u_liq
      print '(a)', 'test_datum:'
      uext = temp_to_internal_energy(DH, WM, 290.0_wp, 1.0_wp)
      call internal_energy_to_temp(uext, WM, DH, t, fl)
      call check('liquid T round-trip', t, 290.0_wp, 1.0e-9_wp)
      call check_true('liquid fliq = 1', abs(fl - 1.0_wp) < 1.0e-12_wp, fl)
      uext = temp_to_internal_energy(DH, WM, 260.0_wp, 0.0_wp)
      call internal_energy_to_temp(uext, WM, DH, t, fl)
      call check('ice T round-trip', t, 260.0_wp, 1.0e-9_wp)
      call check_true('ice fliq = 0', abs(fl) < 1.0e-12_wp, fl)
      u_freeze = (DH + WM * cp_ice) * t_3ple
      u_melt   = u_freeze + WM * latent_heat_fusion
      call internal_energy_to_temp(u_freeze, WM, DH, t, fl)
      call check('continuity at u_freeze', t, t_3ple, 1.0e-9_wp)
      call internal_energy_to_temp(u_melt, WM, DH, t, fl)
      call check('continuity at u_melt', t, t_3ple, 1.0e-9_wp)
      call internal_energy_to_temp(0.5_wp * (u_freeze + u_melt), WM, DH, t, fl)
      call check('plateau temp = t_3ple', t, t_3ple, 1.0e-9_wp)
      call check('plateau fliq = 0.5', fl, 0.5_wp, 1.0e-9_wp)
      u_ice = temp_to_internal_energy(0.0_wp, SWE, 260.0_wp, 0.0_wp)
      call check('ice holds internal_energy_ice per kilogram', u_ice, SWE * internal_energy_ice(260.0_wp), 1.0e-6_wp)
      u_ice = temp_to_internal_energy(0.0_wp, SWE, t_3ple, 0.0_wp)
      u_liq = temp_to_internal_energy(0.0_wp, SWE, t_3ple, 1.0_wp)
      call check('melting at t_3ple costs the latent heat of fusion', u_liq - u_ice, SWE * 3.34e5_wp, 1.0e-3_wp)
   end subroutine test_datum

   !--- 3. The Clausius-Clapeyron slope of the liquid curve against a central difference.
   subroutine test_clausius()
      real(wp), parameter :: DT = 1.0e-3_wp
      real(wp) :: ana, fd
      print '(a)', 'test_clausius:'
      ana = sat_vapor_pressure_temp_deriv(295.0_wp)
      fd  = (sat_vapor_pressure(295.0_wp + DT) - sat_vapor_pressure(295.0_wp - DT)) / (2.0_wp * DT)
      call check('d(e_sat)/dT vs FD', ana, fd, 1.0e-3_wp * abs(fd) + 1.0e-6_wp)
   end subroutine test_clausius

   !---~---
   !   4. SATURATION OVER ICE (docs/science/snow_biophysics.md eq. (1)). The liquid curve over a
   ! frozen surface overstates e_sat by 10 % at -10 C and 34 % at -30 C, which would overdrive
   ! sublimation from the snow pack and from frozen soil. Three things must hold:
   !   (a) an ABSENT fliq gives the liquid curve (fliq = 1) BIT-IDENTICALLY;
   !   (b) the two branches meet EXACTLY at the triple point, so the blend puts no jump into the
   !       right-hand side an adaptive controller integrates;
   !   (c) the ice branch is actually accurate -- checked against Murphy & Koop (2005).
   !---~---
   subroutine test_ice_saturation()
      real(wp), parameter :: DT = 1.0e-3_wp
      real(wp) :: e_liq, e_ice, e_mix, fd, an
      print '(a)', 'test_ice_saturation:'

      !--- (a) absent fliq == fliq = 1 == the liquid curve, to the last bit.
      call check_true('absent fliq is bit-identical to fliq = 1 at 283 K',                        &
                      sat_vapor_pressure(283.15_wp) == sat_vapor_pressure(283.15_wp, 1.0_wp), 0.0_wp)
      call check_true('absent fliq is bit-identical to fliq = 1 at 253 K',                        &
                      sat_vapor_pressure(253.15_wp) == sat_vapor_pressure(253.15_wp, 1.0_wp), 0.0_wp)

      !---~---
      !   (b) the two branches share the 611.2 Pa constant, so they cross EXACTLY at tc = 0,
      ! i.e. at 273.15 K. (t_3ple is 273.16 K, one hundredth of a degree above, where they
      ! differ by 0.06 Pa = 0.01 % -- checked below.) The blend is a convex combination of
      ! two smooth curves, so it is continuous in BOTH temperature and fliq regardless; what
      ! this pins is that a surface crossing the phase boundary sees no step.
      !---~---
      call check('ice and liquid branches cross exactly at 273.15 K',                             &
                 sat_vapor_pressure(273.15_wp, 0.0_wp), sat_vapor_pressure(273.15_wp, 1.0_wp),    &
                 1.0e-12_wp)
      call check('half-frozen too, at 273.15 K', sat_vapor_pressure(273.15_wp, 0.5_wp),           &
                 sat_vapor_pressure(273.15_wp, 1.0_wp), 1.0e-12_wp)
      call check_true('branches agree to < 0.05 % at the triple point',                           &
                      abs(sat_vapor_pressure(t_3ple, 0.0_wp) - sat_vapor_pressure(t_3ple, 1.0_wp))&
                      < 5.0e-4_wp * sat_vapor_pressure(t_3ple, 1.0_wp),                           &
                      sat_vapor_pressure(t_3ple, 0.0_wp) / sat_vapor_pressure(t_3ple, 1.0_wp))

      !--- (c) accuracy against Murphy & Koop (2005) over ice, and the liquid curve's error.
      call check('e_sat over ice at -10 C (Murphy-Koop 259.892 Pa)',                              &
                 sat_vapor_pressure(t_3ple - 10.0_wp, 0.0_wp), 259.892_wp, 0.4_wp)
      call check('e_sat over ice at -20 C (Murphy-Koop 103.252 Pa)',                              &
                 sat_vapor_pressure(t_3ple - 20.0_wp, 0.0_wp), 103.252_wp, 0.5_wp)
      e_liq = sat_vapor_pressure(t_3ple - 10.0_wp, 1.0_wp)
      e_ice = sat_vapor_pressure(t_3ple - 10.0_wp, 0.0_wp)
      call check_true('liquid curve overstates a frozen surface by ~10 % at -10 C',               &
                      e_liq > 1.08_wp * e_ice .and. e_liq < 1.12_wp * e_ice, e_liq / e_ice)
      call check_true('at 268 K the air saturated over ice holds 3-8 % less vapour than over liquid',   &
                      sat_specific_humidity(268.0_wp, 101325.0_wp, 0.0_wp) <                              &
                      0.97_wp * sat_specific_humidity(268.0_wp, 101325.0_wp) .and.                        &
                      sat_specific_humidity(268.0_wp, 101325.0_wp, 0.0_wp) >                              &
                      0.92_wp * sat_specific_humidity(268.0_wp, 101325.0_wp),                             &
                      sat_specific_humidity(268.0_wp, 101325.0_wp, 0.0_wp) / sat_specific_humidity(268.0_wp, 101325.0_wp))

      !--- The blend is monotone in fliq and bracketed by the two pure branches.
      e_mix = sat_vapor_pressure(t_3ple - 15.0_wp, 0.4_wp)
      call check_true('blend lies between the ice and liquid curves',                             &
                      e_mix > sat_vapor_pressure(t_3ple - 15.0_wp, 0.0_wp) .and.                  &
                      e_mix < sat_vapor_pressure(t_3ple - 15.0_wp, 1.0_wp), e_mix)

      !--- The derivative follows the same branch as the value.
      fd = (sat_vapor_pressure(t_3ple - 20.0_wp + DT, 0.0_wp)                                     &
          - sat_vapor_pressure(t_3ple - 20.0_wp - DT, 0.0_wp)) / (2.0_wp * DT)
      an = sat_vapor_pressure_temp_deriv(t_3ple - 20.0_wp, 0.0_wp)
      call check('d(e_sat)/dT over ice vs FD', an, fd, 1.0e-4_wp * abs(fd) + 1.0e-9_wp)
      call check_true('absent fliq derivative is bit-identical to fliq = 1',                      &
                      sat_vapor_pressure_temp_deriv(253.15_wp) ==                                 &
                      sat_vapor_pressure_temp_deriv(253.15_wp, 1.0_wp), 0.0_wp)
   end subroutine test_ice_saturation

   !---~---
   !   5. Ice-aware conductivity: frozen saturated soil conducts more, as k_ice > k_water; and the
   ! Kersten number recovered from the conductivity is Johansen's log10 form thawed and
   ! Farouki's linear form in the relative saturation frozen.
   !---~---
   subroutine test_soil_conductivity()
      real(wp), parameter :: TS = 0.5_wp, TH = 0.15_wp, KSOL = 3.0_wp, KDRY = 0.15_wp
      real(wp) :: k_liq, k_frz, ksat_liq, ksat_frz, ke_liq, ke_frz, sr
      print '(a)', 'test_soil_conductivity:'
      k_liq = soil_thermal_conductivity(0.40_wp, 1.0_wp, 0.43_wp, 3.0_wp, 0.15_wp)   ! all-liquid saturation
      k_frz = soil_thermal_conductivity(0.40_wp, 0.0_wp, 0.43_wp, 3.0_wp, 0.15_wp)   ! all-ice saturation
      call check_true('frozen soil conducts more than liquid', k_frz > k_liq, k_frz - k_liq)
      sr       = TH / TS                                                     ! 0.30 (intermediate sat.)
      ksat_liq = KSOL ** (1.0_wp - TS) * k_water ** TS                       ! all-liquid geometric mean
      ksat_frz = KSOL ** (1.0_wp - TS) * k_ice   ** TS                       ! all-ice geometric mean
      ke_liq   = (soil_thermal_conductivity(TH, 1.0_wp, TS, KSOL, KDRY) - KDRY) / (ksat_liq - KDRY)
      ke_frz   = (soil_thermal_conductivity(TH, 0.0_wp, TS, KSOL, KDRY) - KDRY) / (ksat_frz - KDRY)
      call check('unfrozen Kersten is the log10 form',  ke_liq, log10(sr) + 1.0_wp, 1.0e-9_wp)
      call check('frozen Kersten is the linear s_r form', ke_frz, sr,               1.0e-9_wp)
   end subroutine test_soil_conductivity

   !---~---
   !   6. CLM5's mineral soil at Palo Verde's column (sand 33 %, clay 40 %, porosity 0.4456),
   ! computed apart from ED-MORE: solids (8.80 x 33 + 2.92 x 40) / 73 = 5.57808 W/m/K; dry bulk
   ! density 1496.88 kg/m3, so dry soil (0.135 x 1496.88 + 64.7) / (2700 - 0.947 x 1496.88)
   ! = 0.208022 W/m/K; and the solids' heat capacity (2.128 x 33 + 2.385 x 40) / 73 x 1e6 x
   ! 0.5544 = 1.257835e6 J/m3/K.
   !---~---
   subroutine test_soil_texture()
      real(wp) :: k_solid, k_dry, dry_cvol
      print '(a)', 'test_soil_texture:'
      call soil_texture_to_thermal(33.0_wp, 40.0_wp, 0.4456_wp, k_solid, k_dry, dry_cvol)
      call check('CLM5 solids conductivity', k_solid, 5.5780822_wp, 1.0e-6_wp)
      call check('CLM5 dry conductivity', k_dry, 0.20802202_wp, 1.0e-7_wp)
      call check('CLM5 solids heat capacity', dry_cvol, 1.2578349e6_wp, 1.0_wp)
   end subroutine test_soil_texture

   !---~---
   !   7. Moist air: the humidity of a dewpoint, the VPD at the air's pressure, and the rain/snow
   ! split. The VPD reads the vapour pressure a specific humidity holds at the air's pressure,
   ! which is proportional to it (Dalton): at 70 kPa, about 3,000 m, a humidity holds 0.7 of the
   ! vapour pressure it holds at 1 atm, so read at p_std its VPD would be low by 0.3 of that.
   !---~---
   subroutine test_moist_air()
      real(wp), parameter :: TD = 285.0_wp, P = 98000.0_wp, TOT = 1.0e-4_wp
      real(wp), parameter :: T_AIR = 293.0_wp, Q_AIR = 0.006_wp, P_LOW = 0.7_wp * 101325.0_wp
      real(wp) :: q, e, rain, snow, es, vpd_std, vpd_low
      print '(a)', '-- 7. moist air --'
      q = sat_specific_humidity(TD, P)
      e = sat_vapor_pressure(TD)
      call check('q from dewpoint == eps e/(p-(1-eps)e)', q, h2o_dry_ratio*e/(P-(1.0_wp-h2o_dry_ratio)*e), 1.0e-14_wp)
      call check_true('q in physical range', q > 0.001_wp .and. q < 0.02_wp, q)
      call check('rh=1 at T=Td matches dewpoint q', rh_to_specific_humidity(1.0_wp, TD, P), q, 1.0e-14_wp)
      es = sat_vapor_pressure(T_AIR)
      call check('air at half saturation at 70 kPa has a VPD of half e_sat there',                              &
                 specific_humidity_to_vpd(T_AIR, rh_to_specific_humidity(0.5_wp, T_AIR, P_LOW), P_LOW), 0.5_wp * es, &
                 1.0e-12_wp * es)
      vpd_std = specific_humidity_to_vpd(T_AIR, Q_AIR, 101325.0_wp)
      vpd_low = specific_humidity_to_vpd(T_AIR, Q_AIR, P_LOW)
      call check('at 70 kPa a humidity holds 0.7 of its vapour pressure at 1 atm (Dalton)', es - vpd_low,          &
                 0.7_wp * (es - vpd_std), 1.0e-12_wp * es)
      call precip_phase(TOT, t_3ple + 5.0_wp, rain, snow)
      call check('warm -> all rain', rain, TOT, 1.0e-16_wp)
      call check('warm -> no snow',  snow, 0.0_wp, 1.0e-16_wp)
      call precip_phase(TOT, t_3ple - 5.0_wp, rain, snow)
      call check('cold -> all snow', snow, TOT, 1.0e-16_wp)
      call precip_phase(TOT, t_3ple, rain, snow)
      call check('at t_3ple mass conserved', rain + snow, TOT, 1.0e-16_wp)
      call check_true('partial snow at t_3ple', snow > 0.0_wp .and. snow < TOT, snow)
   end subroutine test_moist_air

   !--- 8. The air at another height: a site 500 m above (dz > 0) or below (dz < 0) its cell.
   subroutine test_lapse()
      real(wp), parameter :: DZ = 500.0_wp, G = 0.0065_wp
      real(wp), parameter :: T0 = 285.0_wp, P0 = 95000.0_wp, Q0 = 0.007_wp
      real(wp) :: t, p, t1, p1, q1, e0, e1
      print '(a)', '-- 8. the air at another height --'
      t = lapse_air_temperature(290.0_wp, DZ, G)
      call check('lapse cools a higher site', t, 290.0_wp - G*DZ, 1.0e-9_wp)
      call check('lapse dz=0 -> T unchanged', lapse_air_temperature(290.0_wp, 0.0_wp, G), 290.0_wp, 1.0e-30_wp)
      p = lapse_pressure(101325.0_wp, 290.0_wp, DZ, G)
      call check_true('pressure drops at a higher site', p < 101325.0_wp, 101325.0_wp - p)
      call check('lapse dz=0 -> P unchanged', lapse_pressure(101325.0_wp, 290.0_wp, 0.0_wp, G), 101325.0_wp, 1.0e-6_wp)
      call check('isothermal (gamma=0) = barometric', lapse_pressure(101325.0_wp, 290.0_wp, DZ, 0.0_wp),        &
                 101325.0_wp*exp(-9.80665_wp*DZ/(287.04_wp*290.0_wp)), 1.0e-2_wp)
      !--- A site 500 m below its cell keeps the cell's relative humidity.
      t1 = lapse_air_temperature(T0, -DZ, G)
      p1 = lapse_pressure(P0, T0, -DZ, G)
      q1 = lapse_specific_humidity(Q0, T0, P0, t1, p1)
      e0 = Q0 * P0 / (h2o_dry_ratio + (1.0_wp - h2o_dry_ratio) * Q0)
      e1 = q1 * p1 / (h2o_dry_ratio + (1.0_wp - h2o_dry_ratio) * q1)
      call check('relative humidity is held', e1 / sat_vapor_pressure(t1), e0 / sat_vapor_pressure(T0), 1.0e-12_wp)
      call check_true('a site 500 m below its cell is warmer and holds more vapour', q1 > Q0, q1 - Q0)
      call check('dz = 0 leaves q unchanged', lapse_specific_humidity(Q0, T0, P0, T0, P0), Q0, 1.0e-15_wp)
   end subroutine test_lapse

   !---~---
   !   The latent heat of liquid equals enthalpy_vapor less internal_energy_liquid and that of ice
   ! enthalpy_vapor less internal_energy_ice, at every temperature, supercooled liquid included;
   ! a mixture is their average weighted by the liquid share.
   !---~---
   subroutine test_latent_heat_vap()
      real(wp) :: t
      integer  :: i
      print '(a)', 'test_latent_heat_vap:'
      do i = 0, 4
         t = 240.0_wp + 20.0_wp * real(i, wp)
         call check('liquid: enthalpy_vapor - internal_energy_liquid', latent_heat_vap(t, 1.0_wp),            &
                    enthalpy_vapor(t) - internal_energy_liquid(t), 1.0e-6_wp)
         call check('ice: enthalpy_vapor - internal_energy_ice', latent_heat_vap(t, 0.0_wp),                  &
                    enthalpy_vapor(t) - internal_energy_ice(t), 1.0e-6_wp)
         call check('mixture: the liquid-weighted average', latent_heat_vap(t, 0.3_wp),                       &
                    0.3_wp * latent_heat_vap(t, 1.0_wp) + 0.7_wp * latent_heat_vap(t, 0.0_wp), 1.0e-6_wp)
      end do
      call check('liquid at t_3ple', latent_heat_vap(t_3ple, 1.0_wp), latent_heat_vap_3ple, 1.0e-9_wp)
      call check('ice at t_3ple', latent_heat_vap(t_3ple, 0.0_wp), latent_heat_vap_3ple + latent_heat_fusion, 1.0e-9_wp)
      call check_true('liquid L_v falls as the air warms', latent_heat_vap(303.15_wp, 1.0_wp) < latent_heat_vap(t_3ple, 1.0_wp))
   end subroutine test_latent_heat_vap

   !---~---
   !   9. A tissue of 5e4 J/K of dry matter holding 12 kg of water at 297.3 K, and 0.7 kg of water
   ! arriving or leaving at that temperature, carrying internal_energy_liquid per kilogram.
   !---~---
   subroutine test_tissue_energy()
      real(wp), parameter :: DRY = 5.0e4_wp, WATER = 12.0_wp, T_TISSUE = 297.3_wp, DW = 0.7_wp
      real(wp) :: e
      print '(a)', 'test_tissue_energy:'
      e = tissue_temp_to_energy(DRY, WATER, T_TISSUE)
      call check_true('the tissue datum is the soil''s: temp_to_internal_energy with liquid water',             &
                      abs(e - temp_to_internal_energy(DRY, WATER, T_TISSUE, 1.0_wp)) <= 1.0e-12_wp * abs(e))
      call check_true('tissue_energy_to_temp inverts tissue_temp_to_energy',                                    &
                      abs(tissue_energy_to_temp(e, WATER, DRY) - T_TISSUE) < 1.0e-10_wp)
      call check_true('water that arrives at the tissue''s temperature leaves it there',                        &
                      abs(tissue_energy_to_temp(e + DW * internal_energy_liquid(T_TISSUE), WATER + DW, DRY)      &
                          - T_TISSUE) < 1.0e-10_wp)
      call check_true('and so does water that leaves at it',                                                    &
                      abs(tissue_energy_to_temp(e - DW * internal_energy_liquid(T_TISSUE), WATER - DW, DRY)      &
                          - T_TISSUE) < 1.0e-10_wp)
   end subroutine test_tissue_energy

end program test_therm_lib
