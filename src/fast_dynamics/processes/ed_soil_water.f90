! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_soil_water -- how water moves through and out of the soil column, each flux with
! its derivatives (docs/science/soil_biophysics.md, "The soil water column"):
!   * soil_water_flux         -- Darcy flow between two layers, or the pond and the top
!                                layer;
!   * bottom_flux             -- drainage through the bottom face, or a two-way aquifer;
!   * philip_evaporation      -- evaporation from a store of the ground's water through given
!                                resistances, the law the top soil layer and a porous layer share;
!   * ground_evaporation_rate -- the top layer's, through its dry surface layer, on the
!                                boundary's terms that ground_evaporation_setup computes;
!   * store_evaporation       -- a store's, supplied while it holds water above its residual.
! Layer water relations are in ed_soil_hydr_lib, root uptake in ed_plant_hydraulics.
!---~---
module ed_soil_water
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : rho_h2o, grav, grav_head, r_wv, tiny_num, p_std, safe_exp, t_kelvin
   use ed_soil_params,        only : soil_profile_t, soil_params_t, SOIL_BC_AQUIFER
   use ed_soil_hydr_lib,      only : soil_head_to_theta
   use ed_numerics,           only : smooth_max0, smooth_step
   use ed_therm_lib,          only : sat_specific_humidity
   implicit none
   private

   public :: soil_water_flux, ground_evaporation_rate, bottom_flux, evaporation_supply
   public :: philip_evaporation, dry_layer_resistance, store_evaporation
   public :: ground_evap_env_t, ground_evaporation_setup, porous_cover_path, vapour_diffusivity
   public :: slot_evaporation_setup

   !---~---
   !   Widths of the two rounded corners in ground_evaporation_rate: [kg/kg] where the surface's
   ! humidity meets the air's (evaporation stops; there is no dew), and [m3/m3] where the top
   ! layer dries below the moisture at which the dry surface layer starts to form.
   !---~---
   real(wp), parameter :: EVAP_HUMIDITY_WIDTH = 1.0e-6_wp
   real(wp), parameter :: DSL_THETA_WIDTH     = 1.0e-3_wp
   !---~---
   !   [-] The share of the top layer's span from residual to saturation over which its
   ! evaporation tapers to nothing as its water above residual, ice included, runs out
   ! (evaporation_supply).
   !---~---
   real(wp), parameter :: EVAP_SUPPLY_SPAN    = 1.0e-4_wp
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
      !---~---
      !   The porous layers above the soil: their depth as an equivalent path in free air, each
      ! layer's depth over its Millington-Quirk tortuosity (0 without one).
      !---~---
      real(wp) :: cover_path = 0.0_wp                     ! [m]
      !--- Set by ground_evaporation_setup from the fields above and the top layer's curve.
      real(wp) :: rv_t            = 0.0_wp   ! [J/kg] R_v T of the ground skin, the scale of Philip's humidity
      real(wp) :: q_sat           = 0.0_wp   ! [kg/kg] saturation specific humidity at the ground skin
      real(wp) :: theta_init      = 0.0_wp   ! [m3/m3] top-layer water content below which the dry layer forms
      real(wp) :: theta_span      = 0.0_wp   ! [m3/m3] theta_init less the air-dry water content
      real(wp) :: dsl_diffusivity = 0.0_wp   ! [m2/s] vapour diffusivity through the dry surface layer
      real(wp) :: r_cover         = 0.0_wp   ! [s/m]  the vapour's resistance through the porous layers
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
      dvap    = vapour_diffusivity(forcing%t_ground)
      phi     = profile%theta_sat(1)
      phi_air = max(phi - theta_air, tiny_num)
      tau     = phi_air ** MQ_TORTUOSITY_EXP / max(phi * phi, tiny_num)   ! [-] Millington-Quirk tortuosity
      forcing%theta_span      = max(forcing%theta_init - theta_air, tiny_num)
      forcing%dsl_diffusivity = max(dvap * tau, tiny_num)
      forcing%r_cover         = forcing%cover_path / dvap
   end subroutine ground_evaporation_setup

   !---~---
   !   The boundary a porous layer's evaporation reads, for philip_evaporation: the canopy air
   ! over the ground (the ground's boundary's), the layer's temperature, and its dried top's path
   ! for the vapour (porous_cover_path of the dried depth, in free air) as the cover in series. A
   ! porous layer has no dry surface layer of its own: its dried top is that path.
   !---~---
   pure function slot_evaporation_setup(ground, t_layer, dry_path) result(env)
      type(ground_evap_env_t), intent(in) :: ground     ! the ground's boundary (the canopy air, its conductance)
      real(wp),                intent(in) :: t_layer    ! [K] the porous layer's temperature
      real(wp),                intent(in) :: dry_path   ! [m] its dried top, as a depth of free air
      type(ground_evap_env_t) :: env
      env = ground
      env%t_ground        = t_layer
      env%rv_t            = r_wv * t_layer
      env%q_sat           = sat_specific_humidity(t_layer, ground%press)
      env%cover_path      = dry_path
      env%r_cover         = dry_path / vapour_diffusivity(t_layer)
   end function slot_evaporation_setup

   !--- [m2/s] Water vapour's diffusivity in air at `temp` [K].
   elemental function vapour_diffusivity(temp) result(dvap)
      real(wp), intent(in) :: temp
      real(wp) :: dvap
      dvap = DVAP_AIR_0 * (temp / t_kelvin) ** DVAP_AIR_TEXP
   end function vapour_diffusivity

   !---~---
   !   One porous layer's path for vapour [m], as a depth of free air: its depth over the
   ! Millington-Quirk tortuosity of its air-filled pores, its porosity less its water. The porosity
   ! is the layer's own, its retention curve's water at saturation, so a saturated layer has no
   ! open pores and passes no vapour.
   !---~---
   elemental function porous_cover_path(depth, theta, porosity) result(path)
      real(wp), intent(in) :: depth      ! [m]
      real(wp), intent(in) :: theta      ! [m3/m3] its water, per volume
      real(wp), intent(in) :: porosity   ! [m3/m3] its water at saturation
      real(wp) :: path, phi_air
      phi_air = max(porosity - theta, tiny_num)
      path    = depth * porosity ** 2 / phi_air ** MQ_TORTUOSITY_EXP
   end function porous_cover_path

   !---~---
   !   philip_evaporation -- evaporation [kg/m2/s] from a store of the ground's water, with its
   ! derivatives: Philip's pore-space humidity at the store's potential, smoothly stopping where it
   ! meets the canopy air's (there is no dew), through the aerodynamic resistance r_aero, the
   ! store's own dry-surface resistance r_dry and the cover above it in series, from the snow-free
   ! share only (docs/science/soil_biophysics.md, "The soil water column"). r_dry's slope in the
   ! store's water content, dr_dry, gives the flux's. The boundary's fixed terms are
   ! ground_evaporation_setup's, or slot_evaporation_setup's.
   !---~---
   pure subroutine philip_evaporation(psi, r_dry, dr_dry, forcing, e, de_dtheta, de_dpsi)
      real(wp),                intent(in)  :: psi         ! [MPa]   the store's matric potential
      real(wp),                intent(in)  :: r_dry       ! [s/m]   its dry-surface resistance
      real(wp),                intent(in)  :: dr_dry      ! [s/m per m3/m3] its slope in the store's water content
      type(ground_evap_env_t), intent(in)  :: forcing     ! the boundary
      real(wp),                intent(out) :: e           ! [kg/m2/s]
      real(wp),                intent(out) :: de_dtheta   ! [kg/m2/s per m3/m3]
      real(wp),                intent(out) :: de_dpsi     ! [kg/m2/s per MPa]
      real(wp) :: alpha, dalpha, drive, ddrive, resist
      alpha  = safe_exp(psi / grav_head * grav / forcing%rv_t)   ! Philip's humidity, on the head
      dalpha = alpha * grav / forcing%rv_t / grav_head           ! per MPa
      call smooth_max0(alpha * forcing%q_sat - forcing%q_air, EVAP_HUMIDITY_WIDTH, drive, ddrive)
      resist = forcing%r_aero + r_dry
      if (forcing%r_cover > 0.0_wp) resist = resist + forcing%r_cover      ! the porous layers above
      e         = forcing%snow_free_frac * forcing%rho_air * drive / resist
      de_dpsi   = forcing%snow_free_frac * forcing%rho_air * ddrive * forcing%q_sat * dalpha / resist
      de_dtheta = -e * dr_dry / resist
   end subroutine philip_evaporation

   !---~---
   !   The top layer's dry surface layer: its resistance to the layer's vapour [s/m] and the
   ! resistance's slope in the layer's water content, the Swenson-Lawrence form (CLM5 eqs.
   ! 5.79-5.80), from theta_init down to the air-dry water content, its corner rounded.
   !---~---
   pure subroutine dry_layer_resistance(theta1, forcing, soil, r_dry, dr_dry)
      real(wp),                intent(in)  :: theta1      ! [m3/m3] the water content the dry layer reads
      type(ground_evap_env_t), intent(in)  :: forcing     ! the boundary, after ground_evaporation_setup
      type(soil_params_t),     intent(in)  :: soil
      real(wp),                intent(out) :: r_dry       ! [s/m]
      real(wp),                intent(out) :: dr_dry      ! [s/m per m3/m3]
      real(wp) :: dry, ddry
      call smooth_max0(forcing%theta_init - theta1, DSL_THETA_WIDTH, dry, ddry)
      r_dry  = soil%dsl_dmax * dry / forcing%theta_span / forcing%dsl_diffusivity
      dr_dry = -soil%dsl_dmax * ddry / forcing%theta_span / forcing%dsl_diffusivity
   end subroutine dry_layer_resistance

   !--- ground_evaporation_rate -- the top layer's evaporation: philip_evaporation through its dry surface layer.
   pure subroutine ground_evaporation_rate(theta1, psi1, forcing, soil, e_soil, de_dtheta, de_dpsi)
      real(wp),               intent(in)  :: theta1      ! [m3/m3] top-layer water content
      real(wp),               intent(in)  :: psi1        ! [MPa]   top-layer matric potential
      type(ground_evap_env_t), intent(in) :: forcing   ! the boundary, after ground_evaporation_setup
      type(soil_params_t),    intent(in)  :: soil
      real(wp),               intent(out) :: e_soil      ! [kg/m2/s] ground evaporation
      real(wp),               intent(out) :: de_dtheta   ! [kg/m2/s per m3/m3]
      real(wp),               intent(out) :: de_dpsi     ! [kg/m2/s per MPa]
      real(wp) :: r_dry, dr_dry
      call dry_layer_resistance(theta1, forcing, soil, r_dry, dr_dry)
      call philip_evaporation(psi1, r_dry, dr_dry, forcing, e_soil, de_dtheta, de_dpsi)
   end subroutine ground_evaporation_rate

   !---~---
   !   store_evaporation -- the evaporation [kg/m2/s] of a store at the ground's top, the top soil
   ! layer or the top porous layer, depth deep and holding water w [kg/m2], and its slope in the
   ! store's unknown, given its water's and its potential's slopes: philip_evaporation through r_dry,
   ! supplied while the store holds water above its residual, ice included (evaporation_supply).
   !---~---
   pure subroutine store_evaporation(forcing, psi, dpsi, w, dw, depth, theta_res, theta_sat, r_dry, dr_dry, e, de)
      type(ground_evap_env_t), intent(in)  :: forcing
      real(wp),                intent(in)  :: psi, dpsi       ! [MPa] its potential and its slope
      real(wp),                intent(in)  :: w, dw           ! [kg/m2] its water and its slope
      real(wp),                intent(in)  :: depth           ! [m]
      real(wp),                intent(in)  :: theta_res, theta_sat   ! [m3/m3] its curve's ends
      real(wp),                intent(in)  :: r_dry, dr_dry   ! [s/m], [s/m per m3/m3] its dry-surface resistance
      real(wp),                intent(out) :: e, de           ! [kg/m2/s], its slope in the unknown
      real(wp) :: theta, de_dtheta, de_dpsi, span, supply, dsupply
      theta = w / (rho_h2o * depth)
      call philip_evaporation(psi, r_dry, dr_dry, forcing, e, de_dtheta, de_dpsi)
      span = theta_sat - theta_res
      call evaporation_supply((theta - theta_res) / span, supply, dsupply)
      de_dtheta = supply * de_dtheta + e * dsupply / span
      de_dpsi   = supply * de_dpsi
      e         = supply * e
      de        = de_dtheta * dw / (rho_h2o * depth) + de_dpsi * dpsi
   end subroutine store_evaporation

   !---~---
   !   The share of its demand a top layer's water can supply [-], and its slope in the layer's
   ! water above residual: one until that water, ice included, is EVAP_SUPPLY_SPAN of the span
   ! from residual to saturation, then a smooth step to nothing at residual. A steep retention
   ! curve keeps its pore air near saturation almost to residual, so Philip's humidity alone
   ! would stop the evaporation only a trillionth of the span above it; CLM5 too limits the
   ! soil's evaporation to the water its top layer holds (docs/science/soil_biophysics.md, "The
   ! soil water column").
   !---~---
   elemental subroutine evaporation_supply(avail, f, df, width)
      real(wp), intent(in)           :: avail   ! [-] the layer's water above residual over its span to saturation
      real(wp), intent(out)          :: f       ! [-]
      real(wp), intent(out)          :: df      ! [-] df / d avail
      real(wp), intent(in), optional :: width   ! [-] the step's width, EVAP_SUPPLY_SPAN unless given
      real(wp) :: w
      w = EVAP_SUPPLY_SPAN
      if (present(width)) w = width
      call smooth_step(avail / w, f, df)
      df = df / w
   end subroutine evaporation_supply

end module ed_soil_water
