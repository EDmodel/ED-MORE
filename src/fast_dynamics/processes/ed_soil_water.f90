! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_soil_water -- how water moves through and out of the soil column, each flux with
! its derivatives (docs/science/soil_biophysics.md, "The soil water column"):
!   * soil_water_flux         -- Darcy flow between two layers, or the pond and the top
!                                layer;
!   * bottom_flux             -- drainage through the bottom face, or a two-way aquifer;
!   * ground_vapour           -- the evaporation of the ground's stores on the snow-free
!                                share, each through its own path up to the ground's surface
!                                node and from there through the one aerodynamic leg;
!   * store_vapour_path       -- a store's path: its pore humidity and its conductance up to
!                                the node, through its dry surface layer and the porous layers
!                                above it, supplied while it holds water above its residual;
!   * ground_evaporation_rate -- the top soil layer's alone, on the boundary's terms that
!                                ground_evaporation_setup computes.
! Layer water relations are in ed_soil_hydr_lib, root uptake in ed_plant_hydraulics.
!---~---
module ed_soil_water
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : rho_h2o, grav, grav_head, r_wv, tiny_num, p_std, safe_exp, t_kelvin
   use ed_soil_params,        only : soil_profile_t, soil_params_t, SOIL_BC_AQUIFER
   use ed_soil_hydr_lib,      only : soil_head_to_theta, PSI_AIR_DRY
   use ed_numerics,           only : smooth_max0, smooth_step
   use ed_therm_lib,          only : sat_specific_humidity
   implicit none
   private

   public :: soil_water_flux, ground_evaporation_rate, bottom_flux, evaporation_supply
   public :: ground_vapour, store_vapour_path, vapour_path_t, dry_layer_resistance
   public :: ground_evap_env_t, ground_evaporation_setup, porous_cover_path, vapour_diffusivity
   public :: slot_evaporation_setup

   !---~---
   !   [m3/m3] The width of the rounded corner in dry_layer_resistance where a store dries below
   ! the moisture at which its dry surface layer starts to form.
   !---~---
   real(wp), parameter :: DSL_THETA_WIDTH     = 1.0e-3_wp
   !---~---
   !   [-] The share of the top layer's span from residual to saturation over which its
   ! evaporation tapers to nothing as its water above residual, ice included, runs out
   ! (evaporation_supply).
   !---~---
   real(wp), parameter :: EVAP_SUPPLY_SPAN    = 1.0e-4_wp
   !---~---
   !   Water vapour's diffusivity in air, D = DVAP_AIR_0 (T/t_kelvin)^DVAP_AIR_TEXP, as CLM5
   ! takes it (Swenson & Lawrence 2014), and the Millington-Quirk (1961) tortuosity
   ! exponent.
   !---~---
   real(wp), parameter :: DVAP_AIR_0        = 2.12e-5_wp            ! [m2/s] at t_kelvin
   real(wp), parameter :: DVAP_AIR_TEXP     = 1.75_wp               ! [-]
   real(wp), parameter :: MQ_TORTUOSITY_EXP = 10.0_wp / 3.0_wp      ! [-]

   !---~---
   !   The boundary a store's evaporation reads (store_vapour_path, ground_vapour): the store's
   ! temperature and the canopy air over the ground with the leg to it, and the terms that depend
   ! only on them and the store's curve.
   !---~---
   type :: ground_evap_env_t
      real(wp) :: t_ground = 298.15_wp                    ! [K] ground skin temperature
      real(wp) :: q_air    = 0.0_wp                       ! [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       ! [kg/m3] canopy-air density (soil evap)
      real(wp) :: press    = p_std                        ! [Pa] air pressure, for the saturation humidity
      real(wp) :: r_aero   = 100.0_wp                     ! [s/m] the ground's aerodynamic resistance to the canopy air
                                                          ! (1/ggnet), not inflated for snow: the one leg the
                                                          ! ground's stores share (ground_vapour)
      real(wp) :: snow_free_frac = 1.0_wp                 ! [-] snow-free area fraction (1-snowfac) that evaporates,
                                                          ! an area-weighted tile flux (ground_surface_layers.md)
      !---~---
      !   The porous layers above the store: their depth as an equivalent path in free air, each
      ! layer's depth over its Millington-Quirk tortuosity (0 without one).
      !---~---
      real(wp) :: cover_path = 0.0_wp                     ! [m]
      !--- Set by ground_evaporation_setup from the fields above and the top layer's curve.
      real(wp) :: rv_t            = 0.0_wp   ! [J/kg] R_v T of the ground skin, the scale of Philip's humidity
      real(wp) :: q_sat           = 0.0_wp   ! [kg/kg] saturation specific humidity at the ground skin
      real(wp) :: theta_init      = 0.0_wp   ! [m3/m3] the store's water content below which its dry layer forms
      real(wp) :: theta_span      = 0.0_wp   ! [m3/m3] theta_init less the water content at which it is whole
      real(wp) :: dsl_diffusivity = 0.0_wp   ! [m2/s] vapour diffusivity through the dry surface layer
      real(wp) :: r_cover         = 0.0_wp   ! [s/m]  the vapour's resistance through the porous layers
   end type ground_evap_env_t

   !---~---
   !   One store's path for its vapour as the shared leg reads it (ground_vapour): the humidity of
   ! its pore air and the conductance of its path up to the ground's surface node, each with its
   ! slope in the store's unknown.
   !---~---
   type :: vapour_path_t
      real(wp) :: q  = 0.0_wp, dq = 0.0_wp    ! [kg/kg], [kg/kg per unit of the unknown]
      real(wp) :: g  = 0.0_wp, dg = 0.0_wp    ! [m/s],   [m/s per unit of the unknown]
   end type vapour_path_t

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
   !   The terms of the top soil layer's evaporation that depend only on its boundary and the
   ! layer's curve, so they are computed once for a boundary: R_v T and the saturation humidity
   ! at the ground skin over liquid at the air's pressure, the dry surface layer's water contents and
   ! vapour diffusivity (dry_layer_setup: it forms below dsl_theta_init of saturation and is
   ! whole at the air-dry water content), and the porous layers' open pores above it as a
   ! resistance (docs/science/soil_plant_water.md §3). A frozen layer's ice enters through its
   ! liquid's potential, which the freezing curve puts at the ice's vapour pressure
   ! (docs/science/soil_biophysics.md, "Frozen soil").
   !---~---
   pure subroutine ground_evaporation_setup(profile, soil, forcing)
      type(soil_profile_t),    intent(in)    :: profile
      type(soil_params_t),     intent(in)    :: soil
      type(ground_evap_env_t), intent(inout) :: forcing   ! its boundary in, its derived terms set
      real(wp) :: theta_air
      forcing%rv_t  = r_wv * forcing%t_ground
      forcing%q_sat = sat_specific_humidity(forcing%t_ground, forcing%press)
      theta_air = soil_head_to_theta(profile%retention, PSI_AIR_DRY, profile%theta_sat(1), profile%theta_res(1),  &
                                     profile%curve_par_a(1), profile%curve_par_n(1))
      call dry_layer_setup(soil%dsl_theta_init * profile%theta_sat(1), theta_air, profile%theta_sat(1), forcing)
      forcing%r_cover = forcing%cover_path / vapour_diffusivity(forcing%t_ground)
   end subroutine ground_evaporation_setup

   !---~---
   !   The boundary a porous layer's evaporation reads (store_vapour_path through
   ! dry_layer_resistance): the canopy air over the ground (the ground's boundary's), the layer's
   ! temperature, the open pores of the porous layers above it as a path in free air
   ! (porous_cover_path), and its own dried top, by the soil's dry-surface-layer law with the
   ! layer's own terms: it forms as soon as the layer holds less than saturated and is the whole
   ! layer when it holds nothing, so the layer's water fills its lower part and the vapour crosses
   ! the dried rest's open pores (dry_layer_setup with theta_on its porosity, theta_off nothing).
   !---~---
   pure function slot_evaporation_setup(ground, t_layer, cover_path, porosity) result(env)
      type(ground_evap_env_t), intent(in) :: ground       ! the ground's boundary (the canopy air, its conductance)
      real(wp),                intent(in) :: t_layer      ! [K] the porous layer's temperature
      real(wp),                intent(in) :: cover_path   ! [m] the porous layers above it, as a depth of free air
      real(wp),                intent(in) :: porosity     ! [m3/m3] its water at saturation
      type(ground_evap_env_t) :: env
      env = ground
      env%t_ground   = t_layer
      env%rv_t       = r_wv * t_layer
      env%q_sat      = sat_specific_humidity(t_layer, ground%press)
      env%cover_path = cover_path
      env%r_cover    = cover_path / vapour_diffusivity(t_layer)
      call dry_layer_setup(porosity, 0.0_wp, porosity, env)
   end function slot_evaporation_setup

   !---~---
   !   A store's dry surface layer, the terms of dry_layer_resistance (the Swenson-Lawrence form,
   ! CLM5 eqs. 5.79-5.80): it forms below theta_on and is whole at theta_off, the vapour crossing
   ! its open pores at water's diffusivity in air at the store's temperature times their
   ! Millington-Quirk tortuosity at theta_off.
   !---~---
   pure subroutine dry_layer_setup(theta_on, theta_off, porosity, forcing)
      real(wp),                intent(in)    :: theta_on, theta_off, porosity   ! [m3/m3]
      type(ground_evap_env_t), intent(inout) :: forcing                         ! its t_ground in, the terms set
      real(wp) :: phi_air, tau
      phi_air = max(porosity - theta_off, tiny_num)
      tau     = phi_air ** MQ_TORTUOSITY_EXP / max(porosity * porosity, tiny_num)   ! [-] Millington-Quirk tortuosity
      forcing%theta_init      = theta_on
      forcing%theta_span      = max(theta_on - theta_off, tiny_num)
      forcing%dsl_diffusivity = max(vapour_diffusivity(forcing%t_ground) * tau, tiny_num)
   end subroutine dry_layer_setup

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
   !   A store's dry surface layer, the one law of the top soil layer and of every porous layer:
   ! its resistance to the store's vapour [s/m] and the resistance's slope in the store's water
   ! content, the Swenson-Lawrence form (CLM5 eqs. 5.79-5.80): as deep as depth_max over the span
   ! from theta_init down to its dry end, its corner rounded, crossed at dsl_diffusivity
   ! (dry_layer_setup). The top soil layer's depth_max is [soil].dsl_dmax, a porous layer's its own
   ! depth.
   !---~---
   pure subroutine dry_layer_resistance(theta, depth_max, forcing, r_dry, dr_dry)
      real(wp),                intent(in)  :: theta       ! [m3/m3] the water content the dry layer reads
      real(wp),                intent(in)  :: depth_max   ! [m] its depth where the store is dry
      type(ground_evap_env_t), intent(in)  :: forcing     ! the boundary, after dry_layer_setup
      real(wp),                intent(out) :: r_dry       ! [s/m]
      real(wp),                intent(out) :: dr_dry      ! [s/m per m3/m3]
      real(wp) :: dry, ddry
      call smooth_max0(forcing%theta_init - theta, DSL_THETA_WIDTH, dry, ddry)
      r_dry  = depth_max * dry / forcing%theta_span / forcing%dsl_diffusivity
      dr_dry = -depth_max * ddry / forcing%theta_span / forcing%dsl_diffusivity
   end subroutine dry_layer_resistance

   !---~---
   !   ground_vapour -- the evaporation [kg/m2/s] of the ground's stores on the snow-free share,
   ! the porous layers' and the top soil layer's, and its slopes (docs/science/soil_biophysics.md,
   ! "The soil water column"). Each store's vapour crosses its own path, its dried top and the
   ! open pores of the porous layers above it, up to the ground's surface node; the paths meet
   ! there in parallel, and the node passes their sum to the canopy air through the one
   ! aerodynamic leg the ground's sensible heat takes, r_aero = 1/g_a. The node's humidity q_s
   ! holds the balance in closed form,
   !     g_a (q_s - q_air) = sum_i g_i (q_i - q_s)
   ! over the stores more humid than the node: a store drier than it takes nothing, so there is
   ! no dew, from the air or between the stores. Those stores are the most humid ones, taken most
   ! humid first while each is more humid than the node of those before it. Store i evaporates
   ! rho g_i (q_i - q_s), written with q_s's closed form as rho g_i/G [g_a (q_i - q_air) +
   ! sum_k g_k (q_i - q_k)], G the sum of the conductances, which holds when a store's path is
   ! open (g_i unbounded, q_s = q_i). de(i, k) is store i's slope in store k's unknown, through
   ! store k's humidity and its path's conductance (vapour_path_t).
   !---~---
   pure subroutine ground_vapour(n, path, forcing, e, de)
      integer(ik),             intent(in)  :: n
      type(vapour_path_t),     intent(in)  :: path(n)
      type(ground_evap_env_t), intent(in)  :: forcing   ! the leg: r_aero, q_air, rho_air, snow_free_frac
      real(wp),                intent(out) :: e(n)      ! [kg/m2/s]
      real(wp),                intent(out) :: de(n, n)  ! [kg/m2/s per unit of the unknown]
      real(wp)    :: g_a, big_g, sum_gq, area_rho, b(n), rest(n)
      integer(ik) :: i, k, next
      logical     :: on(n)
      e  = 0.0_wp
      de = 0.0_wp
      g_a    = 1.0_wp / max(forcing%r_aero, tiny_num)
      big_g  = g_a
      sum_gq = g_a * forcing%q_air
      on = .false.
      do
         next = 0_ik
         do k = 1_ik, n
            if (on(k) .or. .not. (path(k)%g > 0.0_wp)) cycle
            if (next == 0_ik) then
               next = k
            else if (path(k)%q > path(next)%q) then
               next = k
            end if
         end do
         if (next == 0_ik) exit
         if (.not. (path(next)%q * big_g > sum_gq)) exit     ! no more humid than the node so far
         on(next) = .true.
         big_g  = big_g + path(next)%g
         sum_gq = sum_gq + path(next)%g * path(next)%q
      end do
      !--- b(i) = G (q_i - q_s), and rest(i) = G - g_i, both free of g_i.
      do i = 1_ik, n
         if (.not. on(i)) cycle
         b(i)    = g_a * (path(i)%q - forcing%q_air)
         rest(i) = g_a
         do k = 1_ik, n
            if (k == i .or. .not. on(k)) cycle
            b(i)    = b(i) + path(k)%g * (path(i)%q - path(k)%q)
            rest(i) = rest(i) + path(k)%g
         end do
      end do
      area_rho = forcing%snow_free_frac * forcing%rho_air
      do i = 1_ik, n
         if (.not. on(i)) cycle
         e(i)     = area_rho * path(i)%g / big_g * b(i)
         de(i, i) = area_rho * (path(i)%g / big_g * rest(i) * path(i)%dq + b(i) / big_g * rest(i) / big_g * path(i)%dg)
         do k = 1_ik, n
            if (k == i .or. .not. on(k)) cycle
            de(i, k) = -area_rho * path(i)%g / big_g * (path(k)%g * path(k)%dq + b(k) / big_g * path(k)%dg)
         end do
      end do
   end subroutine ground_vapour

   !---~---
   !   store_vapour_path -- the path of a store's vapour up to the ground's surface node
   ! (ground_vapour): the store depth deep, holding water w [kg/m2] at potential psi, with their
   ! slopes in its unknown. Its pore air holds Philip's humidity alpha(psi) q_sat at the store's
   ! temperature (alpha = exp(psi g / R_v T)); its path crosses its dry surface layer r_dry
   ! (dry_layer_resistance, whose slope in the water content is dr_dry) and the porous layers
   ! above it (r_cover), its conductance the share of its demand its water can supply
   ! (evaporation_supply, its water above residual, ice included) over their sum.
   !---~---
   pure function store_vapour_path(forcing, psi, dpsi, w, dw, depth, theta_res, theta_sat, r_dry, dr_dry) result(path)
      type(ground_evap_env_t), intent(in) :: forcing          ! the store's boundary
      real(wp),                intent(in) :: psi, dpsi        ! [MPa] its potential and its slope
      real(wp),                intent(in) :: w, dw            ! [kg/m2] its water and its slope
      real(wp),                intent(in) :: depth            ! [m]
      real(wp),                intent(in) :: theta_res, theta_sat   ! [m3/m3] its curve's ends
      real(wp),                intent(in) :: r_dry, dr_dry    ! [s/m], [s/m per m3/m3] its dry surface layer
      type(vapour_path_t) :: path
      real(wp) :: theta, dtheta, span, supply, dsupply, resist
      path%q  = safe_exp(psi / grav_head * grav / forcing%rv_t) * forcing%q_sat
      path%dq = path%q * grav / (forcing%rv_t * grav_head) * dpsi
      theta  = w / (rho_h2o * depth)
      dtheta = dw / (rho_h2o * depth)
      span   = theta_sat - theta_res
      call evaporation_supply((theta - theta_res) / span, supply, dsupply)
      resist  = max(r_dry + forcing%r_cover, tiny_num)
      path%g  = supply / resist
      path%dg = (dsupply / span - supply * dr_dry / resist) / resist * dtheta
   end function store_vapour_path

   !---~---
   !   ground_evaporation_rate -- the top soil layer's evaporation [kg/m2/s] alone through the
   ! leg (ground_vapour with that one store, its whole supply), and its slopes in the layer's water
   ! content and its potential.
   !---~---
   pure subroutine ground_evaporation_rate(theta1, psi1, forcing, soil, e_soil, de_dtheta, de_dpsi)
      real(wp),                intent(in)  :: theta1      ! [m3/m3] top-layer water content
      real(wp),                intent(in)  :: psi1        ! [MPa]   top-layer matric potential
      type(ground_evap_env_t), intent(in)  :: forcing     ! the boundary, after ground_evaporation_setup
      type(soil_params_t),     intent(in)  :: soil
      real(wp),                intent(out) :: e_soil      ! [kg/m2/s] ground evaporation
      real(wp),                intent(out) :: de_dtheta   ! [kg/m2/s per m3/m3]
      real(wp),                intent(out) :: de_dpsi     ! [kg/m2/s per MPa]
      type(vapour_path_t) :: path(1)
      real(wp) :: r_dry, dr_dry, resist, e(1), de(1, 1)
      call dry_layer_resistance(theta1, soil%dsl_dmax, forcing, r_dry, dr_dry)
      resist = max(r_dry + forcing%r_cover, tiny_num)
      path(1)%q  = safe_exp(psi1 / grav_head * grav / forcing%rv_t) * forcing%q_sat
      path(1)%g  = 1.0_wp / resist
      path(1)%dq = path(1)%q * grav / (forcing%rv_t * grav_head)      ! per MPa
      call ground_vapour(1_ik, path, forcing, e, de)
      e_soil  = e(1)
      de_dpsi = de(1, 1)
      path(1)%dq = 0.0_wp
      path(1)%dg = -dr_dry / resist ** 2                                ! per m3/m3
      call ground_vapour(1_ik, path, forcing, e, de)
      de_dtheta = de(1, 1)
   end subroutine ground_evaporation_rate

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
