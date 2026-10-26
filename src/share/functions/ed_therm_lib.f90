! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_therm_lib -- the thermodynamics the model's stores share, so each uses one formula:   !
! the specific energies of liquid, ice, vapour and canopy air in one datum and the energy  !
! -> temperature inverter with phase change (docs/science/soil_biophysics.md); moist       !
! air: saturation vapour pressure and humidity, VPD and density                            !
! (docs/science/snow_biophysics.md eq. (1), docs/science/forcing.md eq. (11)); the         !
! rain/snow split; the air at another height (docs/science/forcing.md §8); plant           !
! tissue heat capacity and energy (docs/science/vegetation_energy_dynamics.md eqs.         !
! (8)-(10)); soil heat: the thermal properties of a texture and the soil's conductivity    !
! and heat capacity (docs/science/soil_biophysics.md eqs. (9) and (11)); and the           !
! temperature response of enzyme rates (docs/science/plant_respiration.md eq. (5),         !
! docs/science/leaf_gas_exchange.md §2). The Bolton (1980), Johansen and CLM5 coefficients !
! are fixed empirical constants, not model parameters.                                     !
!==========================================================================================!
module ed_therm_lib
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : mmdry, tiny_num, cp_air, cp_vap, cp_liq, cp_ice, latent_heat_fusion,  &
                              t_3ple, tsupercool_liq, tsupercool_vap, r_dry, rho_h2o, k_water,      &
                              k_ice, c2b_wood, r_gas, t_ref_photo, safe_exp, t_kelvin, grav,        &
                              h2o_dry_ratio, virtual_temp_factor
   implicit none
   private

   !----- Energy datum: specific energies and the phase-change inverter. ------------------!
   public :: internal_energy_to_temp, temp_to_internal_energy
   public :: enthalpy_vapor, internal_energy_liquid, internal_energy_ice, liquid_energy_to_temp
   public :: precip_enthalpy, precip_phase
   !----- Moist air. ----------------------------------------------------------------------!
   public :: sat_vapor_pressure, sat_specific_humidity, sat_vapor_pressure_temp_deriv
   public :: sat_specific_humidity_temp_deriv, air_vpd, specific_humidity_to_vpd
   public :: vapor_pressure_to_specific_humidity, specific_humidity_to_vapor_pressure, rh_to_specific_humidity
   public :: air_density, cas_temp_to_enthalpy, cas_enthalpy_to_temp, cas_molar_density
   !----- The air at another height. ------------------------------------------------------!
   public :: lapse_air_temperature, lapse_pressure, lapse_specific_humidity, dry_adiabat_temperature
   !----- Plant tissue heat. --------------------------------------------------------------!
   public :: tissue_heat_capacity, tissue_dry_heat_capacity, tissue_temp_to_energy, tissue_energy_to_temp
   !----- Soil heat. ----------------------------------------------------------------------!
   public :: soil_texture_to_thermal, soil_thermal_conductivity, soil_volumetric_heat_capacity
   !----- Temperature response of enzyme rates. -------------------------------------------!
   public :: arrhenius_scale, peaked_arrhenius_scale, temp_response
   public :: kattge_knorr_entropy, kattge_knorr_jv_ratio
   public :: TEMP_RESPONSE_ARRHENIUS, TEMP_RESPONSE_PEAKED

   !----- CLM5 mineral soil (Lawrence & Slater 2008; CLM5 technote 2.7.2): the solids'     !
   !      conductivity and heat capacity per sand and clay, and the dry soil's             !
   !      conductivity from its bulk density (Johansen 1975); grains are 2700 kg/m3. ------!
   real(wp), parameter :: K_QUARTZ          = 8.80_wp     !< [W/m/K]    sand grains
   real(wp), parameter :: K_CLAY_MINERAL    = 2.92_wp     !< [W/m/K]    clay grains
   real(wp), parameter :: CV_SAND           = 2.128e6_wp  !< [J/m3/K]   sand grains
   real(wp), parameter :: CV_CLAY           = 2.385e6_wp  !< [J/m3/K]   clay grains
   real(wp), parameter :: CLM_GRAIN_DENSITY = 2700.0_wp   !< [kg/m3]
   !----- [-] Floor on the relative saturation, where Johansen's log10 form holds. --------!
   real(wp), parameter :: SR_FLOOR = 0.05_wp

   !----- [K] Half-width of the rain/snow band about the triple point. -------------------!
   real(wp), parameter :: PHASE_BAND_K = 1.0_wp

   !----- The forms [leaf_physiology].temp_response_form selects. -------------------------!
   integer(ik), parameter :: TEMP_RESPONSE_ARRHENIUS = 1_ik  !< plain Arrhenius
   integer(ik), parameter :: TEMP_RESPONSE_PEAKED    = 2_ik  !< Arrhenius with high-temperature deactivation

contains

   !=======================================================================================!
   !     Energy datum: the specific energies of liquid, ice and vapour, and the            !
   !     phase-change inverter (docs/science/soil_biophysics.md).                          !
   !=======================================================================================!

   !----- Phase-change inverter: internal energy, water and dry heat capacity ->           !
   !      temperature and liquid fraction, per m3 of soil or per m2 of surface alike       !
   !      (docs/science/soil_biophysics.md). With no water the plateau is empty. ----------!
   elemental subroutine internal_energy_to_temp(uext, wmass, dry_hcap, temp, fliq)
      real(wp), intent(in)  :: uext, wmass, dry_hcap
      real(wp), intent(out) :: temp, fliq
      real(wp) :: u_freeze, u_melt
      u_freeze = (dry_hcap + wmass * cp_ice) * t_3ple
      u_melt   = u_freeze + wmass * latent_heat_fusion
      if (uext <= u_freeze) then                            ! all ice
         fliq = 0.0_wp
         temp = uext / (dry_hcap + wmass * cp_ice)
      else if (uext >= u_melt) then                         ! all liquid
         fliq = 1.0_wp
         temp = (uext + wmass * cp_liq * tsupercool_liq) / (dry_hcap + wmass * cp_liq)
      else                                                  ! mixed-phase plateau (wmass > 0 here)
         temp = t_3ple
         fliq = (uext - u_freeze) / (wmass * latent_heat_fusion)
      end if
   end subroutine internal_energy_to_temp

   !----- The inverse: temperature and liquid fraction -> internal energy. ----------------!
   elemental function temp_to_internal_energy(dry_hcap, wmass, temp, fliq) result(uext)
      real(wp), intent(in) :: dry_hcap, wmass, temp, fliq
      real(wp)             :: uext
      uext = dry_hcap * temp + wmass * (fliq * cp_liq * (temp - tsupercool_liq)                &
                                        + (1.0_wp - fliq) * cp_ice * temp)
   end function temp_to_internal_energy

   !----- Specific enthalpy of water vapour [J/kg] in the common datum, so a vapour flux   !
   !      carries its latent heat (docs/science/canopy_air_space_biophysics.md). ----------!
   elemental function enthalpy_vapor(t_k) result(h)
      real(wp), intent(in) :: t_k
      real(wp)             :: h
      h = cp_vap * (t_k - tsupercool_vap)
   end function enthalpy_vapor

   !----- Specific internal energy of liquid water [J/kg] (advected soil/xylem water). ----!
   elemental function internal_energy_liquid(t_k) result(u)
      real(wp), intent(in) :: t_k
      real(wp)             :: u
      u = cp_liq * (t_k - tsupercool_liq)
   end function internal_energy_liquid

   !----- The inverse of internal_energy_liquid, which turns a (mass, energy) transfer     !
   !      into the temperature a consumer takes (docs/science/soil_plant_water.md §5). ----!
   elemental function liquid_energy_to_temp(u) result(t_k)
      real(wp), intent(in) :: u
      real(wp)             :: t_k
      t_k = tsupercool_liq + u / cp_liq
   end function liquid_energy_to_temp

   !----- Specific internal energy of ice [J/kg], zero at 0 K: the datum of the inverter's !
   !      all-ice branch, so ice melting at t_3ple costs exactly latent_heat_fusion. ------!
   elemental function internal_energy_ice(t_k) result(u)
      real(wp), intent(in) :: t_k
      real(wp)             :: u
      u = cp_ice * t_k
   end function internal_energy_ice

   !----- The enthalpy precipitation brings, wherever it lands: rain as liquid at the air  !
   !      temperature above the canopy, snow as ice no warmer than freezing                !
   !      (docs/science/forcing.md §8). Rates give [W/m2], amounts [J/m2]. ----------------!
   elemental function precip_enthalpy(rain, snow, t_air) result(h)
      real(wp), intent(in) :: rain, snow   !< [kg/m2/s] or [kg/m2]
      real(wp), intent(in) :: t_air        !< [K] air temperature above the canopy
      real(wp)             :: h
      h = snow * internal_energy_ice(min(t_3ple, t_air)) + rain * internal_energy_liquid(t_air)
   end function precip_enthalpy

   !----- The rain/snow split (docs/science/forcing.md eq. (12)): a linear ramp across a 2 K   !
   !      band centred on the triple point; rainf + snowfall = the total. ---------------------!
   elemental subroutine precip_phase(precip_total, tair_k, rainf, snowfall)
      real(wp), intent(in)  :: precip_total, tair_k
      real(wp), intent(out) :: rainf, snowfall
      real(wp) :: frac_liq
      frac_liq = min(1.0_wp, max(0.0_wp, (tair_k - (t_3ple - PHASE_BAND_K)) / (2.0_wp * PHASE_BAND_K)))
      rainf = frac_liq * precip_total
      snowfall = precip_total - rainf
   end subroutine precip_phase

   !=======================================================================================!
   !     Moist air (docs/science/snow_biophysics.md eq. (1), docs/science/forcing.md eq.   !
   !     (11)).                                                                            !
   !=======================================================================================!

   !----- Saturation vapour pressure [Pa] (docs/science/snow_biophysics.md eq. (1)):       !
   !      Bolton (1980) over liquid, blended with the Magnus form over ice by the liquid   !
   !      fraction fliq. Without fliq it is over liquid, as dewpoint and VPD need          !
   !      (docs/science/forcing.md §7). ---------------------------------------------------!
   elemental function sat_vapor_pressure(t_k, fliq) result(esat)
      real(wp), intent(in)           :: t_k
      real(wp), intent(in), optional :: fliq    !< liquid fraction of the surface (absent = 1, all liquid)
      real(wp)                       :: esat, tc, fl
      tc   = t_k - t_kelvin
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
      if (present(fliq)) then
         fl = min(1.0_wp, max(0.0_wp, fliq))
         if (fl < 1.0_wp) esat = fl * esat                                                       &
                               + (1.0_wp - fl) * 611.2_wp * exp(21.8745584_wp * tc / (tc + 265.5_wp))
      end if
   end function sat_vapor_pressure

   !----- Saturation specific humidity [kg/kg]: q = eps e / (p - (1 - eps) e), eps =        !
   !      h2o_dry_ratio, the conversion of docs/science/forcing.md eq. (11). --------------!
   elemental function sat_specific_humidity(t_k, p_pa, fliq) result(qs)
      real(wp), intent(in)           :: t_k, p_pa
      real(wp), intent(in), optional :: fliq
      real(wp)                       :: qs, esat
      if (present(fliq)) then
         esat = sat_vapor_pressure(t_k, fliq)
      else
         esat = sat_vapor_pressure(t_k)
      end if
      qs   = vapor_pressure_to_specific_humidity(esat, p_pa)
   end function sat_specific_humidity

   !----- Specific humidity [kg/kg] of air at pressure p_pa holding vapour at pressure e_vap  !
   !      [Pa]: q = eps e / (p - (1 - eps) e), eps = h2o_dry_ratio (docs/science/forcing.md    !
   !      eq. (11)). -----------------------------------------------------------------------!
   elemental real(wp) function vapor_pressure_to_specific_humidity(e_vap, p_pa) result(q)
      real(wp), intent(in) :: e_vap, p_pa
      q = h2o_dry_ratio * e_vap / max(p_pa - (1.0_wp - h2o_dry_ratio) * e_vap, tiny_num)
   end function vapor_pressure_to_specific_humidity

   !----- Its inverse: the vapour pressure [Pa] of air at pressure p_pa holding specific      !
   !      humidity q, e = q p / (eps + (1 - eps) q). ----------------------------------------!
   elemental real(wp) function specific_humidity_to_vapor_pressure(q, p_pa) result(e_vap)
      real(wp), intent(in) :: q, p_pa
      e_vap = q * p_pa / max(h2o_dry_ratio + (1.0_wp - h2o_dry_ratio) * q, tiny_num)
   end function specific_humidity_to_vapor_pressure

   !----- Specific humidity [kg/kg] from relative humidity (a fraction, held to [0, 1]) at    !
   !      temperature t_k and pressure p_pa, over liquid water. -----------------------------!
   elemental real(wp) function rh_to_specific_humidity(rh, t_k, p_pa) result(q)
      real(wp), intent(in) :: rh, t_k, p_pa
      q = vapor_pressure_to_specific_humidity(min(1.0_wp, max(0.0_wp, rh)) * sat_vapor_pressure(t_k), p_pa)
   end function rh_to_specific_humidity

   !----- The slope d(e_sat)/dT [Pa/K], at fixed fliq: on the melt plateau, where fliq     !
   !      changes, the temperature is pinned at the triple point. -------------------------!
   elemental function sat_vapor_pressure_temp_deriv(t_k, fliq) result(desat)
      real(wp), intent(in)           :: t_k
      real(wp), intent(in), optional :: fliq
      real(wp)                       :: desat, tc, fl, e_liq, e_ice
      tc    = t_k - t_kelvin
      e_liq = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
      desat = e_liq * 17.67_wp * 243.5_wp / (tc + 243.5_wp) ** 2
      if (present(fliq)) then
         fl = min(1.0_wp, max(0.0_wp, fliq))
         if (fl < 1.0_wp) then
            e_ice = 611.2_wp * exp(21.8745584_wp * tc / (tc + 265.5_wp))
            desat = fl * desat                                                                   &
                  + (1.0_wp - fl) * e_ice * 21.8745584_wp * 265.5_wp / (tc + 265.5_wp) ** 2
         end if
      end if
   end function sat_vapor_pressure_temp_deriv

   !----- The slope d(q_sat)/dT [1/K] for the latent-flux linearizations                   !
   !      (docs/science/vegetation_energy_dynamics.md eq. (2)). ---------------------------!
   elemental function sat_specific_humidity_temp_deriv(t_k, p_pa, fliq) result(dqsdt)
      real(wp), intent(in)           :: t_k, p_pa
      real(wp), intent(in), optional :: fliq
      real(wp)                       :: dqsdt, esat, desat
      if (present(fliq)) then
         esat  = sat_vapor_pressure(t_k, fliq)
         desat = sat_vapor_pressure_temp_deriv(t_k, fliq)
      else
         esat  = sat_vapor_pressure(t_k)
         desat = sat_vapor_pressure_temp_deriv(t_k)
      end if
      dqsdt = h2o_dry_ratio * p_pa / max((p_pa - (1.0_wp - h2o_dry_ratio) * esat) ** 2, tiny_num) * desat
   end function sat_specific_humidity_temp_deriv

   !----- VPD [Pa] from temperature and the actual vapour pressure. -----------------------!
   pure elemental real(wp) function air_vpd(temp, e_vap) result(vpd)
      real(wp), intent(in) :: temp    !< [K]
      real(wp), intent(in) :: e_vap   !< [Pa] actual vapour pressure
      vpd = max(0.0_wp, sat_vapor_pressure(temp) - e_vap)
   end function air_vpd

   !----- VPD [Pa] from temperature, specific humidity and pressure. The vapour pressure e !
   !      = q p / (eps + (1 - eps) q) is the exact inverse of docs/science/forcing.md eq.  !
   !      (11), so a saturated forcing record reads back as VPD = 0. ----------------------!
   pure elemental real(wp) function specific_humidity_to_vpd(temp, shv, pressure) result(vpd)
      real(wp), intent(in) :: temp      !< [K]
      real(wp), intent(in) :: shv       !< [kg/kg] specific humidity
      real(wp), intent(in) :: pressure  !< [Pa]
      real(wp) :: e_vap
      e_vap = specific_humidity_to_vapor_pressure(shv, pressure)
      vpd   = air_vpd(temp, e_vap)
   end function specific_humidity_to_vpd

   !----- Moist-air density [kg/m3] via the virtual temperature. --------------------------!
   elemental function air_density(t_k, p_pa, shv) result(rho)
      real(wp), intent(in) :: t_k, p_pa, shv
      real(wp)             :: rho
      rho = p_pa / (r_dry * t_k * (1.0_wp + virtual_temp_factor * shv))
   end function air_density

   !----- Canopy-air specific enthalpy [J/kg] from temperature and specific humidity: dry  !
   !      air from 0 K, vapour in the common datum                                         !
   !      (docs/science/canopy_air_space_biophysics.md, the three variables). -------------!
   elemental function cas_temp_to_enthalpy(t_k, shv) result(enth)
      real(wp), intent(in) :: t_k, shv
      real(wp)             :: enth
      enth = (1.0_wp - shv) * cp_air * t_k + shv * cp_vap * (t_k - tsupercool_vap)
   end function cas_temp_to_enthalpy

   !----- The inverse: canopy-air temperature [K] from specific enthalpy and humidity. ----!
   elemental function cas_enthalpy_to_temp(enth, shv) result(t_k)
      real(wp), intent(in) :: enth, shv
      real(wp)             :: t_k
      t_k = (enth + shv * cp_vap * tsupercool_vap) / ((1.0_wp - shv) * cp_air + shv * cp_vap)
   end function cas_enthalpy_to_temp

   !----- Molar density of dry air in moist air [mol/m3], rho (1 - q)/mmdry: the canopy    !
   !      air's CO2 capacity per metre of depth and the scale of its exchange with the     !
   !      atmosphere (docs/science/canopy_air_space_biophysics.md, the CO2 variable). -----!
   elemental function cas_molar_density(rho, shv) result(dmol)
      real(wp), intent(in) :: rho, shv
      real(wp)             :: dmol
      dmol = rho * (1.0_wp - shv) / mmdry
   end function cas_molar_density

   !=======================================================================================!
   !     The air at another height (docs/science/forcing.md §8): the temperature, pressure   !
   !     and humidity a linear lapse gives a site dz above the air it is taken from, and   !
   !     the temperature dz above along the dry adiabat.                                   !
   !=======================================================================================!

   !----- Linear environmental lapse of air temperature (docs/science/forcing.md eq.     !
   !      (13)). dz = site - grid; gamma > 0 cools with height. -----------------------!
   elemental function lapse_air_temperature(tair_grid, dz, gamma) result(tair_site)
      real(wp), intent(in) :: tair_grid, dz, gamma
      real(wp) :: tair_site
      tair_site = tair_grid - gamma * dz
   end function lapse_air_temperature

   !----- Hypsometric pressure under the same linear T(z), so P and T stay ideal-gas     !
   !      consistent: P_site = P_grid (T_site/T_grid)^(g/(Rd gamma)), or                 !
   !      exp(-g dz/(Rd T)) as gamma -> 0 (docs/science/forcing.md eq. (13)). ---------!
   elemental function lapse_pressure(psurf_grid, tair_grid, dz, gamma) result(psurf_site)
      real(wp), intent(in) :: psurf_grid, tair_grid, dz, gamma
      real(wp) :: psurf_site, tair_site
      real(wp), parameter :: LAPSE_GAMMA_MIN = 1.0e-6_wp
      if (abs(gamma) > LAPSE_GAMMA_MIN) then
         tair_site  = tair_grid - gamma * dz
         psurf_site = psurf_grid * (tair_site / tair_grid) ** (grav / (r_dry * gamma))
      else
         psurf_site = psurf_grid * exp(-grav * dz / (r_dry * tair_grid))
      end if
   end function lapse_pressure

   !----- Specific humidity at the site, holding relative humidity across the lapse      !
   !      (docs/science/forcing.md eq. (14)); holding q would dry a site below its cell. !
   elemental function lapse_specific_humidity(q_grid, t_grid, p_grid, t_site, p_site) result(q_site)
      real(wp), intent(in) :: q_grid, t_grid, p_grid, t_site, p_site
      real(wp) :: q_site, e, rh
      e      = specific_humidity_to_vapor_pressure(max(q_grid, 0.0_wp), p_grid)
      rh     = e / sat_vapor_pressure(t_grid)
      q_site = rh_to_specific_humidity(rh, t_site, p_site)
   end function lapse_specific_humidity

   !----- Air temperature dz above, along the dry adiabat: T - (g/cp) dz, which conserves    !
   !      the potential temperature theta = T + (g/cp) z (docs/science/forcing.md eq. (15)). !
   elemental real(wp) function dry_adiabat_temperature(t_k, dz) result(t_dz)
      real(wp), intent(in) :: t_k, dz
      t_dz = t_k - (grav / cp_air) * dz
   end function dry_adiabat_temperature

   !=======================================================================================!
   !     Plant tissue heat (docs/science/vegetation_energy_dynamics.md eqs. (8)-(10)).     !
   !=======================================================================================!

   !----- A tissue's heat capacity [J/m2/K]: its dry matter plus the liquid water it holds !
   !      (docs/science/vegetation_energy_dynamics.md eq. (8)). ---------------------------!
   elemental pure real(wp) function tissue_heat_capacity(nplant, carbon, c_dry, water) result(cap)
      real(wp), intent(in) :: nplant   !< [plant/m2]
      real(wp), intent(in) :: carbon   !< [kgC/plant] the tissue's carbon
      real(wp), intent(in) :: c_dry    !< [J/kg/K]    specific heat of its dry matter
      real(wp), intent(in) :: water    !< [kg/plant]  the water it holds
      cap = tissue_dry_heat_capacity(nplant, carbon, c_dry) + water * nplant * cp_liq
   end function tissue_heat_capacity

   !----- The dry part of a tissue's heat capacity [J/m2/K]: its carbon as dry mass, at    !
   !      its dry specific heat (tissue_c_leaf or tissue_c_wood_dry). ---------------------!
   elemental pure real(wp) function tissue_dry_heat_capacity(nplant, carbon, c_dry) result(cap)
      real(wp), intent(in) :: nplant   !< [plant/m2]
      real(wp), intent(in) :: carbon   !< [kgC/plant] the tissue's carbon
      real(wp), intent(in) :: c_dry    !< [J/kg/K]    specific heat of its dry matter
      cap = carbon * nplant * c2b_wood * c_dry
   end function tissue_dry_heat_capacity

   !----- A tissue's internal energy [J per plant, or per m2 with per-m2 inputs], with its !
   !      water liquid at any temperature, since a live tissue is taken not to freeze      !
   !      (docs/science/vegetation_energy_dynamics.md eq. (9)). ---------------------------!
   elemental pure real(wp) function tissue_temp_to_energy(dry_hcap, wmass, temp) result(e)
      real(wp), intent(in) :: dry_hcap   !< [J/K]  dry heat capacity
      real(wp), intent(in) :: wmass      !< [kg]   the water it holds
      real(wp), intent(in) :: temp       !< [K]
      e = dry_hcap * temp + wmass * cp_liq * (temp - tsupercool_liq)
   end function tissue_temp_to_energy

   !----- The inverse (eq. (10)): the temperature of a tissue holding energy e. A tissue   !
   !      with neither dry matter nor water has no temperature of its own, and the caller  !
   !      gives it its parent's. ----------------------------------------------------------!
   elemental pure real(wp) function tissue_energy_to_temp(e, wmass, dry_hcap) result(temp)
      real(wp), intent(in) :: e          !< [J]
      real(wp), intent(in) :: wmass      !< [kg]
      real(wp), intent(in) :: dry_hcap   !< [J/K]
      temp = (e + wmass * cp_liq * tsupercool_liq) / (dry_hcap + wmass * cp_liq)
   end function tissue_energy_to_temp

   !=======================================================================================!
   !     Soil heat (docs/science/soil_biophysics.md §The soil thermal column).             !
   !=======================================================================================!

   !----- The thermal properties of one soil layer by CLM5's mineral-soil formulas         !
   !      (docs/science/soil_biophysics.md eq. (9)). The heat capacity holds the solids    !
   !      only, since the model counts the water and ice itself. The formulas weigh sand   !
   !      and clay alone, so a layer needs some of either: check_soil refuses one without. !
   elemental subroutine soil_texture_to_thermal(sand, clay, theta_sat, k_solid, k_dry, dry_cvol)
      real(wp), intent(in)  :: sand       !< [%]      sand share of the fine earth
      real(wp), intent(in)  :: clay       !< [%]      clay share of the fine earth
      real(wp), intent(in)  :: theta_sat  !< [m3/m3]  porosity of the layer
      real(wp), intent(out) :: k_solid    !< [W/m/K]  conductivity of the soil solids
      real(wp), intent(out) :: k_dry      !< [W/m/K]  conductivity of the dry soil
      real(wp), intent(out) :: dry_cvol   !< [J/m3/K] heat capacity of the solids per layer volume
      real(wp) :: bulk_density
      k_solid      = (K_QUARTZ*sand + K_CLAY_MINERAL*clay) / (sand + clay)
      bulk_density = CLM_GRAIN_DENSITY * (1.0_wp - theta_sat)
      k_dry        = (0.135_wp*bulk_density + 64.7_wp) / (CLM_GRAIN_DENSITY - 0.947_wp*bulk_density)
      dry_cvol     = (CV_SAND*sand + CV_CLAY*clay) / (sand + clay) * (1.0_wp - theta_sat)
   end subroutine soil_texture_to_thermal

   !----- Soil thermal conductivity [W/m/K], Johansen's made ice-aware                     !
   !      (docs/science/soil_biophysics.md eq. (11)). -------------------------------------!
   elemental function soil_thermal_conductivity(theta, fliq, theta_sat, k_solid, k_dry) result(kappa)
      real(wp), intent(in) :: theta, fliq, theta_sat, k_solid, k_dry
      real(wp)             :: kappa, s_r, k_e, kappa_sat
      s_r = min(max(theta / max(theta_sat, tiny_num), SR_FLOOR), 1.0_wp)      ! relative saturation
      !----- Kersten number: Johansen's log form unfrozen, Farouki's linear frozen. -------!
      k_e = fliq * (log10(s_r) + 1.0_wp) + (1.0_wp - fliq) * s_r
      k_e = min(max(k_e, 0.0_wp), 1.0_wp)
      !----- Saturated conductivity: the geometric mean of solids, liquid and ice. --------!
      kappa_sat = k_solid ** (1.0_wp - theta_sat)                                            &
                  * k_water ** (theta_sat * fliq) * k_ice ** (theta_sat * (1.0_wp - fliq))
      kappa = k_e * kappa_sat + (1.0_wp - k_e) * k_dry
   end function soil_thermal_conductivity

   !----- Soil volumetric heat capacity [J/m3/K]: the solids plus the water and ice. ------!
   elemental function soil_volumetric_heat_capacity(theta, fliq, dry_cvol) result(c_eff)
      real(wp), intent(in) :: theta, fliq, dry_cvol
      real(wp)             :: c_eff
      c_eff = dry_cvol + theta * rho_h2o * (fliq * cp_liq + (1.0_wp - fliq) * cp_ice)
   end function soil_volumetric_heat_capacity

   !=======================================================================================!
   !     Temperature response of enzyme rates: how a rate scales from its 25 degC value    !
   !     (t_ref_photo) to the tissue temperature, by the plain Arrhenius form or the       !
   !     peaked form, which adds a high-temperature deactivation (Medlyn et al. 2002;      !
   !     docs/science/plant_respiration.md eq. (5)). Leaf gas exchange and stem and root   !
   !     respiration all scale through these routines, so the tissues share one thermal    !
   !     shape; docs/science/leaf_gas_exchange.md §2 says which rate takes which form.     !
   !=======================================================================================!

   !----- Plain Arrhenius: k25 exp[Ea/(R Tref) (1 - Tref/T)], Tref = t_ref_photo. ---------!
   elemental pure function arrhenius_scale(k25, ea, t_leaf) result(k)
      real(wp), intent(in) :: k25       !< value at the reference temperature
      real(wp), intent(in) :: ea        !< [J/mol] activation energy
      real(wp), intent(in) :: t_leaf    !< [K] tissue temperature
      real(wp)             :: k
      k = k25 * safe_exp(ea / (r_gas * t_ref_photo) * (1.0_wp - t_ref_photo / t_leaf))
   end function arrhenius_scale

   !----- Peaked Arrhenius: the Arrhenius rise times a high-temperature deactivation,      !
   !      normalised so k(t_ref_photo) = k25. ---------------------------------------------!
   elemental pure function peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf) result(k)
      real(wp), intent(in) :: k25       !< value at the reference temperature
      real(wp), intent(in) :: ea        !< [J/mol]   activation energy
      real(wp), intent(in) :: hd        !< [J/mol]   deactivation energy
      real(wp), intent(in) :: ds        !< [J/mol/K] entropy term
      real(wp), intent(in) :: t_leaf    !< [K] tissue temperature
      real(wp)             :: k, fh_ref, fh_leaf
      fh_ref  = 1.0_wp + safe_exp((ds * t_ref_photo - hd) / (r_gas * t_ref_photo))
      fh_leaf = 1.0_wp + safe_exp((ds * t_leaf      - hd) / (r_gas * t_leaf))
      k = arrhenius_scale(k25, ea, t_leaf) * fh_ref / fh_leaf
   end function peaked_arrhenius_scale

   !----- Kattge & Knorr (2007) acclimation of the peaked form's entropy term, dS = a - b  !
   !      T_growth[degC] (docs/science/leaf_gas_exchange.md eq. (4)), with T_growth a      !
   !      running mean of air temperature, the quantity the fit was made on. --------------!
   elemental pure function kattge_knorr_entropy(a_coef, b_coef, t_growth) result(ds)
      real(wp), intent(in) :: a_coef, b_coef   !< [J/mol/K], [J/mol/K2] intercept and slope
      real(wp), intent(in) :: t_growth         !< [K] growth temperature
      real(wp)             :: ds
      ds = a_coef - b_coef * (t_growth - t_kelvin)
   end function kattge_knorr_entropy

   !----- Kattge & Knorr (2007) acclimation of Jmax25/Vcmax25 = a - b T_growth[degC]. -----!
   elemental pure function kattge_knorr_jv_ratio(a_coef, b_coef, t_growth) result(jv)
      real(wp), intent(in) :: a_coef, b_coef   !< [-], [1/K] intercept and slope
      real(wp), intent(in) :: t_growth         !< [K] growth temperature
      real(wp)             :: jv
      jv = max(0.1_wp, a_coef - b_coef * (t_growth - t_kelvin))   ! clamp-ok: positive far outside the 11-35 C fit
   end function kattge_knorr_jv_ratio

   !----- The configured form; the Arrhenius form ignores hd and ds. ----------------------!
   elemental pure function temp_response(form, k25, ea, hd, ds, t_leaf) result(k)
      integer(ik), intent(in) :: form   !< TEMP_RESPONSE_ARRHENIUS | TEMP_RESPONSE_PEAKED
      real(wp),    intent(in) :: k25, ea, hd, ds, t_leaf
      real(wp)                :: k
      select case (form)
      case (TEMP_RESPONSE_PEAKED)
         k = peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf)
      case default
         k = arrhenius_scale(k25, ea, t_leaf)
      end select
   end function temp_response

end module ed_therm_lib
