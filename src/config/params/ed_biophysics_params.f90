! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_biophysics_params -- the [snow] and [aerodynamics] blocks: the snow store's and the
! canopy aerodynamics' parameters. A config leaf (ed_kinds only), so the device-eligible
! kernels import it without ed_config. The ground's are in ed_soil_params.
!---~---
module ed_biophysics_params
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: snow_params_t, aero_params_t

   !---~---
   !   The [snow] block (docs/science/snow_biophysics.md): one layer of fixed density
   ! and conductivity.
   !---~---
   type :: snow_params_t
      real(wp) :: rho_snow          = 250.0_wp    !--- [kg/m3] fixed seasonal-mean bulk density (depth = swe/rho)
      real(wp) :: k_snow            = 0.15_wp     !--- [W/m/K] fixed snow thermal conductivity
      real(wp) :: albedo_vis_fresh  = 0.85_wp     !--- [-] fresh-snow VIS albedo
      real(wp) :: albedo_nir_fresh  = 0.65_wp     !--- [-] fresh-snow NIR albedo
      real(wp) :: albedo_vis_aged   = 0.518_wp    !--- [-] aged VIS endpoint
      real(wp) :: albedo_nir_aged   = 0.435_wp    !--- [-] aged NIR endpoint
      real(wp) :: snow_emiss        = 0.97_wp     !--- [-] thermal (LW) emissivity
      real(wp) :: liquid_holding_frac = 0.10_wp   !--- [-] percolation holding capacity (LEAF-3 1:9)
      real(wp) :: min_new_snow_mass = 1.0e-3_wp   !--- [kg/m2] first-flake creation threshold
      real(wp) :: tiny_snow_mass    = 1.0e-6_wp   !--- [kg/m2] melt-out threshold (dump residual to soil)
      real(wp) :: ny07_a            = 2.5_wp      !--- [-] Niu-Yang07 snow-cover-fraction shape factor
      real(wp) :: ny07_m            = 1.0_wp      !--- [-] Niu-Yang07 snow-cover-fraction exponent
   end type snow_params_t

   !---~---
   !   The [aerodynamics] block (docs/science/canopy_aerodynamics.md); the von Karman
   ! constant is in ed_constants.
   !---~---
   type :: aero_params_t
      real(wp) :: z0m_ratio = 0.13_wp               !--- z0m / canopy height
      real(wp) :: d_ratio   = 0.63_wp               !--- displacement / canopy height
      real(wp) :: snow_rough = 0.001_wp             !--- [m] snow-surface roughness (blend floor)
      real(wp) :: soil_rough = 0.01_wp              !--- [m] bare soil's roughness length (CLM's zlnd)
      !--- Monin-Obukhov (CLM5).
      real(wp)    :: zeta_m = 1.574_wp, zeta_t = 0.465_wp   !--- momentum / scalar range transitions
      real(wp)    :: zeta_max_stable = 0.5_wp       !--- cap on the stable zeta, the first guess's and each sweep's
      integer(ik) :: n_iter_mo = 20_ik              !--- the most MO sweeps; they stop once the coupling settles
      real(wp)    :: wc = 0.5_wp                     !--- [m/s] convective velocity of the first stability guess
      !--- Floors.
      real(wp) :: ustmin = 0.10_wp, ubmin = 0.65_wp, ugbmin = 0.25_wp   !--- [m/s] u*, bulk-wind and in-canopy wind
      real(wp) :: gbhmos_min = 1.0e-9_wp            !--- [m/s] min boundary-layer conductance
      real(wp) :: min_canopy_depth = 5.0_wp         !--- [m] CAS depth floor
      real(wp) :: canopy_freeboard  = 5.0_wp        !--- [m] CAS depth above the tallest cohort
      !--- Ground conductance (CLM4 dense canopy).
      real(wp) :: cs_dense = 0.004_wp, gamma_g = 0.5_wp
      !--- Boundary-layer heat:vapour ratio + air properties (linear-in-T).
      real(wp) :: gbh_2_gbw = 1.075_wp
      real(wp) :: kin_visc0 = 1.33e-5_wp, dkin_visc = 0.007_wp   !--- kinematic viscosity [m2/s] @ t_ref_air
      real(wp) :: th_diff0  = 1.89e-5_wp, dth_diff  = 0.007_wp   !--- thermal diffusivity [m2/s] @ t_ref_air
      real(wp) :: t_ref_air = 273.15_wp             !--- [K] 0 C reference for the linear-in-T props
      !--- Nusselt coefficients -- flat plate (leaf).
      real(wp) :: aflat_lami = 0.60_wp, nflat_lami = 0.50_wp
      real(wp) :: aflat_turb = 0.032_wp, nflat_turb = 0.80_wp
      real(wp) :: bflat_lami = 0.50_wp, mflat_lami = 0.25_wp
      real(wp) :: bflat_turb = 0.19_wp, mflat_turb = 1.0_wp/3.0_wp
      !--- Nusselt coefficients -- cylinder (wood).
      real(wp) :: ocyli_lami = 0.32_wp, acyli_lami = 0.51_wp, ncyli_lami = 0.52_wp
      real(wp) :: ocyli_turb = 0.0_wp,  acyli_turb = 0.24_wp, ncyli_turb = 0.60_wp
      real(wp) :: bcyli_lami = 0.48_wp, mcyli_lami = 0.25_wp
      real(wp) :: bcyli_turb = 0.09_wp, mcyli_turb = 1.0_wp/3.0_wp
   end type aero_params_t

end module ed_biophysics_params
