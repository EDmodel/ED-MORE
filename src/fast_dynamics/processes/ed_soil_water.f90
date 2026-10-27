! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_soil_water -- how water moves through and out of the soil column, each flux with
! its derivatives (docs/science/soil_biophysics.md, "The soil water column"):
!   * soil_water_flux         -- Darcy flow between two layers, or the pond and the top
!                                layer;
!   * bottom_flux             -- drainage through the bottom face, or a two-way aquifer;
!   * ground_evaporation_rate -- evaporation from the top layer, on the boundary's terms
!                                that ground_evaporation_setup computes;
!   * pond_inflow             -- the water reaching the pond and its temperature.
! Layer water relations are in ed_soil_hydr_lib, root uptake in ed_plant_hydraulics.
!---~---
module ed_soil_water
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : rho_h2o, grav, grav_head, r_wv, tiny_num, p_std, safe_exp, t_kelvin
   use ed_soil_params,        only : soil_profile_t, soil_params_t, SOIL_BC_AQUIFER
   use ed_soil_hydr_lib,      only : soil_head_to_theta
   use ed_numerics,           only : smooth_max0
   use ed_therm_lib,          only : sat_specific_humidity, internal_energy_liquid, liquid_energy_to_temp,    &
                                     precip_enthalpy
   implicit none
   private

   public :: soil_water_flux, ground_evaporation_rate, bottom_flux
   public :: ground_evap_env_t, ground_evaporation_setup, pond_inflow

   !---~---
   !   Widths of the two rounded corners in ground_evaporation_rate: [kg/kg] where the surface's
   ! humidity meets the air's (evaporation stops; there is no dew), and [m3/m3] where the top
   ! layer dries below the moisture at which the dry surface layer starts to form.
   !---~---
   real(wp), parameter :: EVAP_HUMIDITY_WIDTH = 1.0e-6_wp
   real(wp), parameter :: DSL_THETA_WIDTH     = 1.0e-3_wp
   !--- CLM5 eq 5.78: the 'air dry' matric head [m], where the soil counts as dry (CLM5's 1e7 mm).
   real(wp), parameter :: PSI_AIR_DRY = -1.0e4_wp
   !---~---
   !   Water vapour's diffusivity in air, D = DVAP_AIR_0 (T/t_kelvin)^DVAP_AIR_TEXP, as CLM5
   ! takes it (Swenson & Lawrence 2014), and the Millington-Quirk (1961) tortuosity
   ! exponent.
   !---~---
   real(wp), parameter :: DVAP_AIR_0        = 2.12e-5_wp            ! [m2/s] at t_kelvin
   real(wp), parameter :: DVAP_AIR_TEXP     = 1.75_wp               ! [-]
   real(wp), parameter :: MQ_TORTUOSITY_EXP = 10.0_wp / 3.0_wp      ! [-]

   !---~---
   !   The boundary ground_evaporation_rate reads: the ground skin and the canopy air over it,
   ! and the terms that depend only on them and the top layer's curve.
   !---~---
   type :: ground_evap_env_t
      real(wp) :: t_ground = 298.15_wp                    ! [K] ground skin temperature
      real(wp) :: q_air    = 0.0_wp                       ! [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       ! [kg/m3] canopy-air density (soil evap)
      real(wp) :: press    = p_std                        ! [Pa] air pressure, for the saturation humidity
      real(wp) :: r_aero   = 100.0_wp                     ! [s/m] aerodynamic resistance of the bare-soil tile
                                                          ! (1/ggnet), not inflated for snow
      real(wp) :: snow_free_frac = 1.0_wp                 ! [-] snow-free area fraction (1-snowfac) that evaporates,
                                                          ! an area-weighted tile flux (snow_biophysics.md)
      !--- Set by ground_evaporation_setup from the fields above and the top layer's curve.
      real(wp) :: rv_t            = 0.0_wp   ! [J/kg] R_v T of the ground skin, the scale of Philip's humidity
      real(wp) :: q_sat           = 0.0_wp   ! [kg/kg] saturation specific humidity at the ground skin
      real(wp) :: theta_init      = 0.0_wp   ! [m3/m3] top-layer water content below which the dry layer forms
      real(wp) :: theta_span      = 0.0_wp   ! [m3/m3] theta_init less the air-dry water content
      real(wp) :: dsl_diffusivity = 0.0_wp   ! [m2/s] vapour diffusivity through the dry surface layer
   end type ground_evap_env_t

contains

   !---~---
   !   The Darcy flux of soil water [kg/m2/s] from node a to node b, two layers or the pond and
   ! the top layer, with its derivatives: the face conductivity k_a^(1-w_b) k_b^w_b times the
   ! potential difference plus the elevation drop z_drop of a above b, over the node distance
   ! dist (docs/science/soil_biophysics.md eq. (3)).
   !---~---
   elemental subroutine soil_water_flux(k_a, k_b, w_b, psi_a, psi_b, z_drop, dist, q,               &
                                        dq_dpsi_a, dq_dpsi_b, dq_dk_a, dq_dk_b)
      real(wp), intent(in)  :: k_a, k_b, w_b, psi_a, psi_b, z_drop, dist   ! [m/s], [m/s], [-], [MPa], [MPa], [m], [m]
      real(wp), intent(out) :: q, dq_dpsi_a, dq_dpsi_b, dq_dk_a, dq_dk_b   ! [kg/m2/s], [kg/m2/s/MPa] (twice), [kg/m3] (twice)
      real(wp) :: k_face, drive
      k_face = k_a ** (1.0_wp - w_b) * k_b ** w_b
      drive  = ((psi_a - psi_b) / grav_head + z_drop) / dist
      q      = rho_h2o * k_face * drive
      dq_dpsi_a =  rho_h2o * k_face / (grav_head * dist)
      dq_dpsi_b = -dq_dpsi_a
      dq_dk_a   = (1.0_wp - w_b) * q / max(k_a, tiny_num)
      dq_dk_b   = w_b * q / max(k_b, tiny_num)
   end subroutine soil_water_flux

   !---~---
   !   The bottom-face flux [m/s, down-positive] with its derivatives: drainage at
   ! drain_factor times K(theta_n), or the two-way aquifer boundary, the column's own face
   ! rule (soil_water_flux, eq. (3)) between layer n's node and the saturated table at the
   ! base, half a layer below it (docs/science/soil_biophysics.md, "Boundary conditions,
   ! and what the aquifer option actually means"). psi_n arrives in MPa. The path to the
   ! table lies in layer n, so the table's side takes the layer's ice impedance too.
   !---~---
   pure subroutine bottom_flux(profile, soil, n, psi_n, k_n, impedance, qbot, dq_dpsi, dq_dk)
      type(soil_profile_t), intent(in)  :: profile
      type(soil_params_t),  intent(in)  :: soil
      integer(ik),          intent(in)  :: n        ! the deepest active layer
      real(wp),             intent(in)  :: psi_n    ! [MPa] matric potential of the deepest layer
      real(wp),             intent(in)  :: k_n      ! [m/s] K(theta_n) of the deepest layer
      real(wp),             intent(in)  :: impedance   ! [-] its ice impedance (soil_ice_t)
      real(wp),             intent(out) :: qbot     ! [m/s] down-positive
      real(wp), optional,   intent(out) :: dq_dpsi  ! [m/s/MPa] d(qbot)/d(psi_n)
      real(wp), optional,   intent(out) :: dq_dk    ! [-]   d(qbot)/d(k_n)
      real(wp) :: delta, q, dq_dpsi_n, dq_dpsi_table, dq_dk_n, dq_dk_table
      qbot = 0.0_wp
      if (present(dq_dpsi)) dq_dpsi = 0.0_wp
      if (present(dq_dk))   dq_dk   = 0.0_wp
      if (soil%water_bottom_bc == SOIL_BC_AQUIFER) then
         !---~---
         !   The table is saturated (psi 0, K_sat), half a layer below the node; the face takes
         ! the geometric mean of the two conductivities, as the pond face does at the top.
         !---~---
         delta = 0.5_wp * profile%dz(n)
         call soil_water_flux(k_n, profile%ksat(n) * impedance, 0.5_wp, psi_n, 0.0_wp, delta, delta, q,  &
                              dq_dpsi_n, dq_dpsi_table, dq_dk_n, dq_dk_table)
         qbot = q / rho_h2o
         if (present(dq_dpsi)) dq_dpsi = dq_dpsi_n / rho_h2o
         if (present(dq_dk))   dq_dk   = dq_dk_n / rho_h2o
      else                                                              ! drainage, q = f K(theta_n)
         qbot = soil%drain_factor * k_n
         if (present(dq_dk)) dq_dk = soil%drain_factor
      end if
   end subroutine bottom_flux

   !---~---
   !   The terms of ground_evaporation_rate that depend only on its boundary and the top layer's
   ! curve, so they are computed once for a boundary: R_v T and the saturation humidity at the
   ! ground skin over liquid at the air's pressure, and the dry surface layer's water-content
   ! span and vapour diffusivity (docs/science/soil_plant_water.md §3). A frozen layer's ice
   ! enters through its liquid's potential, which the freezing curve puts at the ice's vapour
   ! pressure (docs/science/soil_biophysics.md, "Frozen soil").
   !---~---
   pure subroutine ground_evaporation_setup(profile, soil, forcing)
      type(soil_profile_t),    intent(in)    :: profile
      type(soil_params_t),     intent(in)    :: soil
      type(ground_evap_env_t), intent(inout) :: forcing   ! its boundary in, its derived terms set
      real(wp) :: theta_air, dvap, phi, phi_air, tau
      forcing%rv_t  = r_wv * forcing%t_ground
      forcing%q_sat = sat_specific_humidity(forcing%t_ground, forcing%press)
      !--- The dry surface layer, referenced to the air-dry water content.
      theta_air = soil_head_to_theta(profile%retention, PSI_AIR_DRY, profile%theta_sat(1), profile%theta_res(1),  &
                                     profile%curve_par_a(1), profile%curve_par_n(1))
      forcing%theta_init = soil%dsl_theta_init * profile%theta_sat(1)
      dvap    = DVAP_AIR_0 * (forcing%t_ground / t_kelvin) ** DVAP_AIR_TEXP   ! [m2/s] vapour diffusivity in air
      phi     = profile%theta_sat(1)
      phi_air = max(phi - theta_air, tiny_num)
      tau     = phi_air ** MQ_TORTUOSITY_EXP / max(phi * phi, tiny_num)   ! [-] Millington-Quirk tortuosity
      forcing%theta_span      = max(forcing%theta_init - theta_air, tiny_num)
      forcing%dsl_diffusivity = max(dvap * tau, tiny_num)
   end subroutine ground_evaporation_setup

   !---~---
   !   ground_evaporation_rate -- ground evaporation [kg/m2/s] from the top layer, with its
   ! derivatives: Philip's pore-space humidity through the Swenson-Lawrence
   ! dry-surface-layer resistance (CLM5 eqs. 5.77-5.80) in series with r_aero, from the
   ! snow-free share only (docs/science/soil_biophysics.md, "The soil water column"). Both
   ! corners are rounded, so the flux and its derivatives are continuous. The boundary's
   ! fixed terms are ground_evaporation_setup's.
   !---~---
   pure subroutine ground_evaporation_rate(theta1, psi1, forcing, soil, e_soil, de_dtheta, de_dpsi)
      real(wp),               intent(in)  :: theta1      ! [m3/m3] top-layer water content
      real(wp),               intent(in)  :: psi1        ! [MPa]   top-layer matric potential
      type(ground_evap_env_t), intent(in) :: forcing   ! the boundary, after ground_evaporation_setup
      type(soil_params_t),    intent(in)  :: soil
      real(wp),               intent(out) :: e_soil      ! [kg/m2/s] ground evaporation
      real(wp),               intent(out) :: de_dtheta   ! [kg/m2/s per m3/m3]
      real(wp),               intent(out) :: de_dpsi     ! [kg/m2/s per MPa]
      real(wp) :: alpha_soil, dalpha, drive, ddrive, dry, ddry, r_soil, dr_dtheta, resist
      alpha_soil = safe_exp(psi1 / grav_head * grav / forcing%rv_t)   ! Philip's humidity, on the head
      dalpha     = alpha_soil * grav / forcing%rv_t / grav_head          ! per MPa
      call smooth_max0(alpha_soil * forcing%q_sat - forcing%q_air, EVAP_HUMIDITY_WIDTH, drive, ddrive)
      !--- The dry surface layer, from theta_init down to the air-dry water content.
      call smooth_max0(forcing%theta_init - theta1, DSL_THETA_WIDTH, dry, ddry)
      r_soil    = soil%dsl_dmax * dry / forcing%theta_span / forcing%dsl_diffusivity
      dr_dtheta = -soil%dsl_dmax * ddry / forcing%theta_span / forcing%dsl_diffusivity
      resist    = forcing%r_aero + r_soil
      e_soil    = forcing%snow_free_frac * forcing%rho_air * drive / resist
      de_dpsi   = forcing%snow_free_frac * forcing%rho_air * ddrive * forcing%q_sat * dalpha / resist
      de_dtheta = -e_soil * dr_dtheta / resist
   end subroutine ground_evaporation_rate

   !---~---
   !   pond_inflow -- the water reaching the pond and the temperature that values it: the
   ! meltwater under a pack, the throughfall on bare ground, either joined by the water
   ! the tissues shed (docs/science/soil_plant_water.md §5, "The ground inflow"). When
   ! nothing falls, the throughfall is valued at the canopy-air temperature t_cas.
   ! boundary_enth is the part of the inflow's enthalpy that crosses the column's
   ! boundary: all of it on bare ground, the shed water's under a pack.
   !---~---
   pure subroutine pond_inflow(pack_exists, melt_rate, t_melt, throughfall, rainfall, snowfall, t_air, t_cas,     &
                               shed_rate, shed_enth, inflow, t_inflow, t_water_base, boundary_enth)
      logical,  intent(in)  :: pack_exists     !           a snow pack covers the ground
      real(wp), intent(in)  :: melt_rate       ! [kg/m2/s] the pack's meltwater
      real(wp), intent(in)  :: t_melt          ! [K]       its temperature
      real(wp), intent(in)  :: throughfall     ! [kg/m2/s] rain and snow through the canopy
      real(wp), intent(in)  :: rainfall, snowfall   ! [kg/m2/s] above the canopy
      real(wp), intent(in)  :: t_air           ! [K]       the air temperature precipitation is valued at
      real(wp), intent(in)  :: t_cas           ! [K]       canopy-air temperature (nothing falling)
      real(wp), intent(in)  :: shed_rate       ! [kg/m2/s] water the tissues shed to the ground
      real(wp), intent(in)  :: shed_enth       ! [W/m2]    its enthalpy
      real(wp), intent(out) :: inflow          ! [kg/m2/s] all the water reaching the pond
      real(wp), intent(out) :: t_inflow        ! [K]       the temperature that values it
      real(wp), intent(out) :: t_water_base    ! [K]       the meltwater's or the rain's, before the shed water
      real(wp), intent(out) :: boundary_enth   ! [W/m2]    what it brings across the column's boundary
      real(wp) :: water_base, inflow_enth
      if (pack_exists) then
         water_base = melt_rate
      else
         water_base = throughfall
      end if
      inflow = water_base + shed_rate
      t_water_base = t_cas
      if (pack_exists) then
         t_water_base = t_melt
      else if (rainfall + snowfall > tiny_num) then
         t_water_base = liquid_energy_to_temp(precip_enthalpy(rainfall, snowfall, t_air) / (rainfall + snowfall))
      end if
      inflow_enth = water_base * internal_energy_liquid(t_water_base) + shed_enth
      t_inflow = t_water_base
      if (shed_rate > 0.0_wp) t_inflow = liquid_energy_to_temp(inflow_enth / inflow)
      boundary_enth = inflow_enth
      if (pack_exists) boundary_enth = shed_enth
   end subroutine pond_inflow

end module ed_soil_water
