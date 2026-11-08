! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_constants -- the model's fixed numbers, in groups: mathematical constants, physical
! constants, reference conditions, conversions, array sizes and floating-point guards; and
! safe_exp, the exponential that stays inside those guards.
!
!   No mutable state lives here, and every name is lower case. Solver tolerances are not
! here: the ones a run may tune are config keys (docs/configuration.md), and a fixed one
! sits beside the routine that uses it.
!---~---
module ed_constants
   use ed_kinds, only : wp, ik
   implicit none
   public

   !--- Mathematical constants.
   real(wp), parameter :: pi   = 3.14159265358979323846_wp
   real(wp), parameter :: pio2 = 0.5_wp * pi                  ! [rad] pi/2
   real(wp), parameter :: pio4 = 0.25_wp * pi                 ! basal area = pio4 * dbh^2

   !--- Physical constants: universal.
   real(wp), parameter :: r_gas            = 8.314462618_wp      ! [J/mol/K] universal gas constant
   real(wp), parameter :: stefan_boltzmann = 5.670374419e-8_wp   ! [W/m2/K4] Stefan-Boltzmann constant
   real(wp), parameter :: von_karman       = 0.4_wp              ! [-] von Karman constant
   real(wp), parameter :: grav             = 9.80665_wp          ! [m/s2] standard gravity

   !--- Physical constants: dry air and water vapour.
   real(wp), parameter :: mmdry  = 0.0289655_wp               ! [kg/mol] dry-air molar mass
   real(wp), parameter :: mmh2o  = 0.01801528_wp              ! [kg/mol] water molar mass
   real(wp), parameter :: r_dry  = 287.04_wp                  ! [J/kg/K] dry-air gas constant
   real(wp), parameter :: r_wv   = 461.5_wp                   ! [J/kg/K] water-vapour gas constant
   real(wp), parameter :: cp_air = 1004.6_wp                  ! [J/kg/K] dry-air specific heat
   real(wp), parameter :: cp_vap = 1859.0_wp                  ! [J/kg/K] water-vapour specific heat
   !---~---
   !   The molar-mass ratios moist air needs: vapour pressure from specific humidity, and the
   ! virtual temperature T (1 + virtual_temp_factor q).
   !---~---
   real(wp), parameter :: h2o_dry_ratio       = mmh2o / mmdry           ! [-] 0.62196
   real(wp), parameter :: virtual_temp_factor = mmdry / mmh2o - 1.0_wp  ! [-] 0.60783

   !--- Physical constants: liquid water and ice.
   real(wp), parameter :: rho_h2o            = 1000.0_wp      ! [kg/m3] liquid-water density
   real(wp), parameter :: t_3ple             = 273.16_wp      ! [K] triple point of water (0.01 degC)
   real(wp), parameter :: cp_liq             = 4186.0_wp      ! [J/kg/K] liquid-water specific heat
   real(wp), parameter :: cp_ice             = 2093.0_wp      ! [J/kg/K] ice specific heat
   real(wp), parameter :: k_water            = 0.57_wp        ! [W/m/K] liquid-water thermal conductivity
   real(wp), parameter :: k_ice              = 2.29_wp        ! [W/m/K] ice thermal conductivity
   real(wp), parameter :: latent_heat_vap_3ple = 2.501e6_wp   ! [J/kg] latent heat of vaporization at t_3ple
   real(wp), parameter :: latent_heat_fusion = 3.34e5_wp      ! [J/kg] latent heat of fusion
   !--- [Pa] Saturation vapour pressure at t_3ple, Ambaum (2020, QJRMS 146:4252, doi:10.1002/qj.3899) eq. (14).
   real(wp), parameter :: es_3ple            = 611.655_wp
   real(wp), parameter :: grav_head          = rho_h2o * grav * 1.0e-6_wp  ! [MPa/m] water potential per metre of height
   !---~---
   !   One energy zero for ice, liquid and vapour: ice at 0 K. Liquid water then holds
   ! cp_liq * (T - tsupercool_liq) and vapour cp_vap * (T - tsupercool_vap); see
   ! ed_therm_lib and docs/science/canopy_air_space_biophysics.md.
   !---~---
   real(wp), parameter :: tsupercool_liq = t_3ple - (cp_ice * t_3ple + latent_heat_fusion) / cp_liq
   real(wp), parameter :: tsupercool_vap = t_3ple                                                     &
                          - (cp_ice * t_3ple + latent_heat_vap_3ple + latent_heat_fusion) / cp_vap

   !---~---
   !   Physical constants: water vapour and CO2 diffusing through a leaf, as the ratio of
   ! their conductances.
   !---~---
   real(wp), parameter :: gsw_2_gsc = 1.6_wp                  ! [-] stomata: gs(CO2) = gs(H2O) / gsw_2_gsc
   real(wp), parameter :: gbw_2_gbc = 1.4_wp                  ! [-] boundary layer: gb(CO2) = gb(H2O) / gbw_2_gbc

   !---~---
   !   Physical constants: plant tissue heat (docs/science/vegetation_energy_dynamics.md
   ! eq. (8)). The fast loop and the slow side share them, so both weigh a tissue's
   ! heat alike.
   !---~---
   real(wp), parameter :: c2b_wood          = 2.0_wp          ! [-] carbon -> dry biomass (carbon fraction 0.5)
   real(wp), parameter :: tissue_c_leaf     = 3200.0_wp       ! [J/kg/K] dry leaf specific heat
   real(wp), parameter :: tissue_c_wood_dry = 1217.0_wp       ! [J/kg/K] dry wood specific heat, 103.1 + 3.867 T at 15 degC

   !--- Reference conditions.
   real(wp), parameter :: p_std       = 101325.0_wp           ! [Pa] standard atmospheric pressure
   real(wp), parameter :: t_ref_photo = 298.15_wp             ! [K] 25 degC, the reference of the photosynthesis rates
   real(wp), parameter :: rho_air_ref = 1.2_wp                ! [kg/m3] air density the slow step values the canopy air at

   !--- Conversions: time, on the mean Gregorian year (365.2425 days).
   real(wp), parameter :: day_sec = 86400.0_wp                ! [s/day]
   real(wp), parameter :: yr_day  = 365.2425_wp               ! [day/yr]
   real(wp), parameter :: yr_sec  = yr_day * day_sec          ! [s/yr]

   !--- Conversions: units.
   real(wp), parameter :: t_kelvin       = 273.15_wp             ! [K] 0 degC in kelvin
   real(wp), parameter :: deg_to_rad     = pi / 180.0_wp         ! [rad/deg]
   real(wp), parameter :: cm2_to_m2      = 1.0e-4_wp             ! [m2/cm2] basal and sapwood areas are kept in cm2
   real(wp), parameter :: mol_2_umol     = 1.0e6_wp              ! [umol/mol] mole fraction -> ppm
   real(wp), parameter :: umol_2_kgc     = 1.20107e-8_wp         ! [kgC/umol] carbon in 1 umol CO2 (12.0107 g/mol)
   real(wp), parameter :: kgc_2_umol     = 1.0_wp / umol_2_kgc   ! [umol/kgC]
   real(wp), parameter :: kgcday_2_umols = kgc_2_umol / day_sec  ! [(umol/m2/s) per (kgC/m2/day)], about 963.6
   real(wp), parameter :: par_w_2_umol   = 4.6_wp                ! [umol/J] PAR energy -> photons (solar spectrum: ~4.57)

   !---~---
   !   Array sizes: the most soil and snow layers a column can hold. These are storage
   ! limits, not model parameters, and no result depends on them: the layers a run uses
   ! come from its config ([soil].depth and layer_bottoms), and the loader refuses
   ! a column with more layers than fit. The soil and snow stores are fixed-size arrays,
   ! so copying a patch's stores never allocates. The snow is one bulk layer
   ! (multi-layer snow is not modelled).
   !---~---
   integer(ik), parameter :: n_soil_layer_max = 20_ik
   integer(ik), parameter :: n_snow_layer_max = 1_ik

   !--- Floating-point guards: keep a division or an exp() inside its domain.
   real(wp), parameter :: tiny_num   = 1.0e-30_wp             ! keeps a divisor away from zero
   real(wp), parameter :: almost_one = 1.0_wp - 1.0e-7_wp     ! a fraction this close to one counts as one
   real(wp), parameter :: lnexp_min  = -38.0_wp               ! clamp for exp() arguments
   real(wp), parameter :: lnexp_max  =  38.0_wp

contains

   !--- Numerically safe exponential: clamps the argument to avoid under/overflow.
   elemental pure function safe_exp(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = exp(min(max(x, lnexp_min), lnexp_max))
   end function safe_exp

end module ed_constants
