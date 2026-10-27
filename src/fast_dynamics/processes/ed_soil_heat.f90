! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_soil_heat -- the soil thermal column, prognostic in each layer's internal energy,
! with temperature and liquid fraction diagnosed from it (docs/science/soil_biophysics.md,
! "The soil thermal column"). Two forms share the conductive and water-enthalpy faces:
!   * soil_energy_time_deriv    -- the explicit tendency dE_k/dt [W/m3], faces at T^n;
!   * soil_energy_step_implicit -- the backward-Euler advance over dt, faces at T^{n+1}.
! set_soil_heat_env assembles their forcing from the ground heat flux and the water
! network's flows; soil_layer_temp is one layer's temperature from its energy.
!---~---
module ed_soil_heat
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : rho_h2o, n_soil_layer_max, cp_liq, cp_ice, tsupercool_liq, t_3ple
   use ed_soil_params,        only : soil_profile_t, soil_params_t,                                    &
                                     ENERGY_BC_DIRICHLET
   use ed_patch_state_vars,     only : soil_energy_state_t
   use ed_therm_lib,          only : internal_energy_liquid, soil_thermal_conductivity,                         &
                                     soil_volumetric_heat_capacity
   use ed_soil_hydr_lib,      only : soil_layer_curve_t, soil_liquid_capacity
   use ed_linear_algebra,     only : thomas_solve
   implicit none
   private

   public :: soil_energy_step_implicit, soil_energy_time_deriv, set_soil_heat_env, soil_layer_temp
   public :: soil_layer_energy, soil_layer_energy_slope

   !---~---
   !   The inversion of a layer's energy on its freezing curve: the most Newton steps, and the
   ! temperature step [K] below which it has converged.
   !---~---
   integer(ik), parameter :: SOIL_TEMP_ITER_MAX = 60_ik
   real(wp),    parameter :: SOIL_TEMP_TOL      = 1.0e-10_wp
   !---~---
   !   The step's Newton on the layers' energies: the most iterations, the halvings of a step
   ! that does not lower the worst residual, and the residual [K, as the temperature change
   ! it would need] at which it has converged.
   !---~---
   integer(ik), parameter :: SOIL_STEP_ITER_MAX = 30_ik
   integer(ik), parameter :: SOIL_STEP_HALVINGS = 5_ik
   real(wp),    parameter :: SOIL_STEP_TOL      = 1.0e-8_wp
   public :: soil_heat_env_t, soil_heat_out_t

   !---~---
   !   The soil column's thermal forcing: the heat and water crossing its top, bottom and
   ! faces, and the soil moisture its thermal properties depend on.
   !---~---
   type :: soil_heat_env_t
      real(wp) :: g_top      = 0.0_wp                       !--- [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !--- [W/m2] prescribed bottom flux, positive up (held at 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !--- [m3/m3] soil moisture, for kappa and C_eff
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !--- [m/s]  inter-layer water flux, upward-positive like hf
                                                            !    (set_soil_heat_env flips the network's)
      real(wp) :: water_heat_out(n_soil_layer_max) = 0.0_wp !--- [W/m2] enthalpy leaving each layer with water other
                                                            !    than across a face (root uptake, exfiltration, drainage)
      !---~---
      !   The top face's water, applied by the same upwind rule and time level as the interior
      ! faces (docs/science/soil_biophysics.md, "Water-thermal coupling").
      !---~---
      real(wp) :: w_flux_top = 0.0_wp                       !--- [m/s] water across the top face (<0 = infiltration in)
      real(wp) :: t_water_top = 0.0_wp                      !--- [K] temperature of the infiltrating pond water
   end type soil_heat_env_t

   !--- The soil column's energy outputs and diagnostics.
   type :: soil_heat_out_t
      real(wp) :: bottom_heat  = 0.0_wp                     !--- [W/m2] heat conducted out through the bottom face
      real(wp) :: energy_resid = 0.0_wp                     !--- [J/m2] closed-budget residual (~0)
      real(wp) :: heat_capacity = 0.0_wp                    !--- [J/m2/K] the column's, at the step's start
      integer(ik) :: iterations = 0_ik                      !--- Newton iterations the step took
   end type soil_heat_out_t

   !---~---
   !   The soil column at one set of energies (soil_energy_step_implicit): each layer's
   ! temperature, apparent heat capacity and balance residual, the faces, and the worst
   ! residual as a temperature change.
   !---~---
   type :: soil_column_eval_t
      real(wp) :: temp(n_soil_layer_max) = 0.0_wp      !--- [K]
      real(wp) :: de(n_soil_layer_max)   = 0.0_wp      !--- [J/m3/K] dE/dT on the freezing curve
      real(wp) :: res(n_soil_layer_max)  = 0.0_wp      !--- [W/m2]
      real(wp) :: hf(0:n_soil_layer_max)  = 0.0_wp     !--- [W/m2] conduction, positive up
      real(wp) :: qwf(0:n_soil_layer_max) = 0.0_wp     !--- [W/m2] the water's enthalpy, positive up
      real(wp) :: norm = 0.0_wp                        !--- [K]
   end type soil_column_eval_t

contains

   !---~---
   !   set_soil_heat_env -- the soil column's thermal forcing from the ground heat flux and
   ! the water network's flows (docs/science/soil_biophysics.md, "Water-thermal
   ! coupling"): the enthalpy leaving each layer other than across a face (root uptake,
   ! the top layer's exfiltration, and the drainage through the bottom layer), the face
   ! flows turned from down-positive to the column's up-positive, and the infiltration
   ! across the top face. Runoff leaves the pond, not the column. water_heat and e_drain
   ! return the two parts of what leaves.
   !---~---
   pure subroutine set_soil_heat_env(nsl, g_top, theta, root_heat, exfil_enth, face_flux, infiltration,       &
                                     t_infil, drainage, t_drain, env, water_heat, e_drain)
      integer(ik),           intent(in)  :: nsl            !---           number of active layers
      real(wp),              intent(in)  :: g_top          !--- [W/m2]    net ground heat flux into the soil top
      real(wp),              intent(in)  :: theta(:)       !--- [m3/m3]   soil moisture (thermal properties)
      real(wp),              intent(in)  :: root_heat(:)   !--- [W/m2]    enthalpy leaving each layer with root uptake
      real(wp),              intent(in)  :: exfil_enth     !--- [W/m2]    enthalpy leaving the top layer by exfiltration
      real(wp),              intent(in)  :: face_flux(:)   !--- [m/s]     down-positive water flux at the layer faces
      real(wp),              intent(in)  :: infiltration   !--- [kg/m2/s] top-face infiltration (pond -> soil)
      real(wp),              intent(in)  :: t_infil        !--- [K]       temperature of the infiltrating water
      real(wp),              intent(in)  :: drainage       !--- [kg/m2/s] out of the bottom layer
      real(wp),              intent(in)  :: t_drain        !--- [K]       temperature it leaves at
      type(soil_heat_env_t), intent(out) :: env            !---           the column's thermal forcing
      real(wp),              intent(out) :: water_heat(:)  !--- [W/m2]    root uptake and exfiltration, by layer
      real(wp),              intent(out) :: e_drain        !--- [W/m2]    the drainage's
      integer(ik) :: k
      water_heat(1:nsl) = root_heat(1:nsl)
      water_heat(1)     = water_heat(1) + exfil_enth
      e_drain = drainage * internal_energy_liquid(t_drain)
      env%g_top = g_top
      do k = 1_ik, nsl
         env%soil_water(k)     = theta(k)
         env%water_heat_out(k) = water_heat(k)
         env%w_flux(k)         = -face_flux(k)
      end do
      env%w_flux_top  = -infiltration / rho_h2o
      env%t_water_top = t_infil
      env%water_heat_out(nsl) = env%water_heat_out(nsl) + e_drain
   end subroutine set_soil_heat_env

   !---~---
   !   One layer's temperature and liquid fraction from its internal energy and water content,
   ! on its freezing curve (soil_liquid_capacity; docs/science/soil_biophysics.md, "Frozen
   ! soil"): all liquid at and above t_3ple, or below it while the curve holds the layer's
   ! water; below that the liquid the curve holds at the temperature, the rest ice. Water a
   ! layer stores beyond saturation freezes at t_3ple itself, where the curve holds theta_sat,
   ! so an energy between the curve's there and the all-liquid one sits at t_3ple. The energy
   ! rises with the temperature, and a partly frozen layer is warmer than its all-liquid
   ! temperature and colder than its all-ice one and t_3ple, so Newton's iterates are kept
   ! inside that bracket.
   !---~---
   elemental pure subroutine soil_layer_temp(soil_energy, theta, dry_heat_capacity, curve, temp, fliq)
      real(wp),                 intent(in)  :: soil_energy         !--- [J/m3]   internal energy
      real(wp),                 intent(in)  :: theta               !--- [m3/m3]  water content
      real(wp),                 intent(in)  :: dry_heat_capacity   !--- [J/m3/K] heat capacity of the dry soil
      type(soil_layer_curve_t), intent(in)  :: curve               !--- the layer's retention curve
      real(wp),                 intent(out) :: temp                !--- [K]      temperature
      real(wp),                 intent(out) :: fliq                !--- [-]      liquid fraction of the water
      real(wp) :: de
      call soil_layer_state(soil_energy, theta, dry_heat_capacity, curve, temp, fliq, de)
   end subroutine soil_layer_temp

   !---~---
   !   soil_layer_temp, with the apparent heat capacity de [J/m3/K] at the temperature found
   ! (soil_layer_energy_slope's), which the step's Newton reads.
   !---~---
   elemental pure subroutine soil_layer_state(soil_energy, theta, dry_heat_capacity, curve, temp, fliq, de)
      real(wp),                 intent(in)  :: soil_energy, theta, dry_heat_capacity   !--- as soil_layer_temp's
      type(soil_layer_curve_t), intent(in)  :: curve
      real(wp),                 intent(out) :: temp, fliq, de      !--- [K], [-], [J/m3/K]
      real(wp)    :: wmass, theta_liq, dtheta, lo, hi, e, step
      integer(ik) :: it
      wmass = theta * rho_h2o
      temp  = (soil_energy + wmass * cp_liq * tsupercool_liq) / (dry_heat_capacity + wmass * cp_liq)
      fliq  = 1.0_wp
      de    = dry_heat_capacity + wmass * cp_liq
      if (wmass <= 0.0_wp .or. temp >= t_3ple) return        ! dry, or liquid at and above t_3ple
      call soil_liquid_capacity(curve, temp, theta_liq, dtheta)
      if (theta_liq >= theta) return                         ! the curve holds all of it: all liquid
      !--- At t_3ple, between the curve's energy there and the all-liquid one.
      call soil_layer_energy_slope(t_3ple, theta, dry_heat_capacity, curve, e, de, fliq)
      if (soil_energy >= e) then
         temp = t_3ple
         fliq = (soil_energy - (dry_heat_capacity + wmass * cp_ice) * t_3ple)                              &
                / (wmass * (internal_energy_liquid(t_3ple) - cp_ice * t_3ple))
         return
      end if
      lo = temp
      hi = min(t_3ple, soil_energy / (dry_heat_capacity + wmass * cp_ice))
      temp = hi
      do it = 1_ik, SOIL_TEMP_ITER_MAX
         call soil_layer_energy_slope(temp, theta, dry_heat_capacity, curve, e, de, fliq)
         if (e > soil_energy) then
            hi = temp
         else
            lo = temp
         end if
         step = (soil_energy - e) / de
         if (temp + step <= lo .or. temp + step >= hi) step = 0.5_wp * (lo + hi) - temp   ! bisect
         temp = temp + step
         if (abs(step) < SOIL_TEMP_TOL) exit
      end do
      call soil_layer_energy_slope(temp, theta, dry_heat_capacity, curve, e, de, fliq)
   end subroutine soil_layer_state

   !---~---
   !   The inverse: one layer's internal energy [J/m3] at a temperature, with the liquid its
   ! freezing curve holds there.
   !---~---
   elemental pure function soil_layer_energy(temp, theta, dry_heat_capacity, curve) result(soil_energy)
      real(wp),                 intent(in) :: temp, theta, dry_heat_capacity   !--- [K], [m3/m3], [J/m3/K]
      type(soil_layer_curve_t), intent(in) :: curve
      real(wp) :: soil_energy, de, fliq
      call soil_layer_energy_slope(temp, theta, dry_heat_capacity, curve, soil_energy, de, fliq)
   end function soil_layer_energy

   !---~---
   !   A layer's internal energy [J/m3] at a temperature, its slope, the apparent heat capacity
   ! [J/m3/K] (the sensible heat and the latent heat of the liquid the curve frees per
   ! kelvin), and its liquid fraction: all liquid above t_3ple, and at t_3ple the bottom of
   ! the plateau on which water stored beyond saturation freezes.
   !---~---
   elemental pure subroutine soil_layer_energy_slope(temp, theta, dry_heat_capacity, curve, e, de, fliq)
      real(wp),                 intent(in)  :: temp, theta, dry_heat_capacity   !--- [K], [m3/m3], [J/m3/K]
      type(soil_layer_curve_t), intent(in)  :: curve
      real(wp),                 intent(out) :: e, de, fliq                      !--- [J/m3], [J/m3/K], [-]
      real(wp) :: wmass, theta_liq, dtheta, dfliq
      wmass = theta * rho_h2o
      fliq  = 1.0_wp
      dfliq = 0.0_wp
      if (wmass > 0.0_wp .and. temp <= t_3ple) then
         call soil_liquid_capacity(curve, temp, theta_liq, dtheta)
         if (theta_liq < theta) then
            fliq  = theta_liq / theta
            dfliq = dtheta / theta
         end if
      end if
      e  = dry_heat_capacity * temp + wmass * (fliq * internal_energy_liquid(temp)                        &
                                             + (1.0_wp - fliq) * cp_ice * temp)
      de = dry_heat_capacity + wmass * (fliq * cp_liq + (1.0_wp - fliq) * cp_ice                          &
                                       + dfliq * (internal_energy_liquid(temp) - cp_ice * temp))
   end subroutine soil_layer_energy_slope

   !---~---
   !   Each layer's temperature, liquid fraction and ice-aware conductivity at its current
   ! energy, for the step and the tendency alike.
   !---~---
   pure subroutine layer_thermal_state(col, forcing, profile, n, temp, fliq, kappa)
      type(soil_energy_state_t), intent(in)  :: col
      type(soil_heat_env_t),     intent(in)  :: forcing
      type(soil_profile_t),      intent(in)  :: profile
      integer(ik),               intent(in)  :: n                  !--- number of active layers
      real(wp),                  intent(out) :: temp(:), fliq(:)   !--- [K], [-]
      real(wp),                  intent(out) :: kappa(:)           !--- [W/m/K]
      integer(ik) :: k
      do k = 1_ik, n
         call soil_layer_temp(col%soil_energy(k), forcing%soil_water(k), profile%dry_heat_capacity(k),          &
                              profile%layer_curve(k), temp(k), fliq(k))
         kappa(k) = soil_thermal_conductivity(forcing%soil_water(k), fliq(k), profile%theta_sat(k),            &
                                              profile%solid_conductivity(k), profile%dry_conductivity(k))
      end do
   end subroutine layer_thermal_state

   !---~---
   !   Each interior face's conductivity [W/m/K]: the two half-layers in series, the
   ! dz-weighted harmonic mean of the layers' conductivities.
   !---~---
   pure subroutine face_conductivity(dz, kappa, n, kf)
      real(wp),    intent(in)  :: dz(:), kappa(:)   !--- [m] layer thicknesses, [W/m/K] conductivities
      integer(ik), intent(in)  :: n                 !--- number of active layers
      real(wp),    intent(out) :: kf(:)   !--- [W/m/K] face k, between layers k and k+1 (k < n)
      integer(ik) :: k
      do k = 1_ik, n - 1_ik
         kf(k) = (dz(k) + dz(k+1)) / (dz(k) / kappa(k) + dz(k+1) / kappa(k+1))
      end do
   end subroutine face_conductivity

   !---~---
   !   water_enthalpy_faces -- the liquid enthalpy moving water carries across each soil
   ! face [W/m2, positive up], upwind on the given layer temperatures (the qwf terms of
   ! docs/science/soil_biophysics.md eq. (10)): one rule for the step and the tendency.
   !---~---
   pure subroutine water_enthalpy_faces(forcing, temp, n, qwf)
      type(soil_heat_env_t), intent(in)  :: forcing
      real(wp),               intent(in)  :: temp(:)   !--- [K] layer temperatures to upwind on (T^n or T^{n+1})
      integer(ik),            intent(in)  :: n         !--- number of active layers
      real(wp),               intent(out) :: qwf(0:)   !--- [W/m2] face 0 = top ... face n = bottom
      integer(ik) :: k
      !---~---
      !   Top face: infiltration carries the pond's water at t_water_top. Exfiltration leaves
      ! through water_heat_out, so the face never flows up: the second branch is unused.
      !---~---
      if (forcing%w_flux_top <= 0.0_wp) then
         qwf(0) = forcing%w_flux_top * rho_h2o * internal_energy_liquid(forcing%t_water_top)
      else
         qwf(0) = forcing%w_flux_top * rho_h2o * internal_energy_liquid(temp(1))
      end if
      do k = 1_ik, n - 1_ik
         if (forcing%w_flux(k) <= 0.0_wp) then
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(temp(k))
         else
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(temp(k+1))
         end if
      end do
      !---~---
      !   Bottom face: no advective term; drainage enthalpy leaves through the bottom layer's
      ! water_heat_out (set_soil_heat_env).
      !---~---
      qwf(n) = 0.0_wp
   end subroutine water_enthalpy_faces

   !---~---
   !   deep_conductance -- the deep-temperature anchor's conductance [W/m2/K]: the bottom
   ! layer's conductivity over the distance from the bottom node down to deep_depth; zero
   ! for the geothermal (Neumann) boundary (docs/science/soil_biophysics.md eq. (12)).
   !---~---
   pure function deep_conductance(soil, profile, kappa_bot, n) result(g_deep)
      type(soil_params_t),  intent(in) :: soil
      type(soil_profile_t), intent(in) :: profile
      real(wp),             intent(in) :: kappa_bot       !--- [W/m/K] bottom-layer conductivity
      integer(ik),          intent(in) :: n               !--- number of active layers
      real(wp) :: g_deep, l_deep
      if (soil%thermal_bottom_bc == ENERGY_BC_DIRICHLET) then
         l_deep = soil%deep_depth - abs(profile%z_node(n))
         g_deep = kappa_bot / l_deep
      else
         g_deep = 0.0_wp
      end if
   end function deep_conductance

   !---~---
   !   bottom_heat_face -- the bottom face's heat flux [W/m2, positive up], one expression
   ! for the step, the tendency and the solve's matrix: the prescribed geothermal flux
   ! plus conduction to the deep anchor, whose conductance is zero for the geothermal
   ! boundary (docs/science/soil_biophysics.md eq. (12)).
   !---~---
   pure function bottom_heat_face(soil, forcing, g_deep, t_bot) result(hf_bot)
      type(soil_params_t),    intent(in) :: soil
      type(soil_heat_env_t), intent(in) :: forcing
      real(wp),               intent(in) :: g_deep      !--- [W/m2/K] from deep_conductance
      real(wp),               intent(in) :: t_bot       !--- [K] bottom-layer temperature
      real(wp) :: hf_bot
      hf_bot = forcing%geothermal - g_deep * (t_bot - soil%deep_temp)
   end function bottom_heat_face

   !---~---
   !   soil_energy_step_implicit -- one patch's soil thermal column advanced over dt by
   ! backward Euler in the layers' energies (docs/science/soil_biophysics.md eq. (10)):
   ! Newton on each layer's balance dz (E - E^n)/dt = face divergence + source, with the
   ! faces at the temperatures the new energies give on the freezing curve and the
   ! conductivities held at the step's start. The energies are then committed from the
   ! faces at the solution, so the column's budget closes to round-off whatever Newton's
   ! residual, and the temperatures diagnosed from them.
   !---~---
   subroutine soil_energy_step_implicit(col, forcing, profile, soil, dt, flux)
      type(soil_energy_state_t),   intent(inout) :: col       !--- the layers' energy, temperature, liquid fraction
      type(soil_heat_env_t),      intent(in)    :: forcing
      type(soil_profile_t),        intent(in)    :: profile
      type(soil_params_t),         intent(in)    :: soil
      real(wp),                    intent(in)    :: dt        !--- [s]
      type(soil_heat_out_t),         intent(out)   :: flux    !--- bottom heat flux, budget residual, heat capacity

      integer(ik) :: n, k, it, ib
      real(wp), dimension(n_soil_layer_max) :: fl_n, kappa, kf, q_src, e_n, e_new, e_try, step
      real(wp), dimension(n_soil_layer_max) :: lo_t, di_t, up_t, a, b, c, rhs
      real(wp) :: g_deep, g_face, w_cp, lam, norm_try, e0, e1
      type(soil_column_eval_t) :: now, try

      n = profile%n_active

      !---~---
      !   Each layer's temperature, liquid fraction and apparent heat capacity at the step's
      ! start, and its ice-aware conductivity, held over the step.
      !---~---
      e_n(1:n) = col%soil_energy(1:n)
      call soil_layer_state(e_n(1:n), forcing%soil_water(1:n), profile%dry_heat_capacity(1:n),                &
                            profile%layer_curve(1:n), now%temp(1:n), fl_n(1:n), now%de(1:n))
      do k = 1_ik, n
         kappa(k) = soil_thermal_conductivity(forcing%soil_water(k), fl_n(k), profile%theta_sat(k),            &
                                              profile%solid_conductivity(k), profile%dry_conductivity(k))
      end do
      call face_conductivity(profile%dz, kappa, n, kf)
      g_deep = deep_conductance(soil, profile, kappa(n), n)
      do k = 1_ik, n
         q_src(k) = -forcing%water_heat_out(k) / profile%dz(k)               ! [W/m3] source (sink is negative)
      end do

      !---~---
      !   The balance's slopes in the temperatures [W/m2/K], fixed over the step: conduction
      ! across each face and to the deep anchor, and the water's enthalpy, whose upwind
      ! donor's temperature sets it. lo_t, di_t and up_t are row k's slopes in T(k-1), T(k)
      ! and T(k+1).
      !---~---
      lo_t(1:n) = 0.0_wp
      di_t(1:n) = 0.0_wp
      up_t(1:n) = 0.0_wp
      do k = 1_ik, n - 1_ik
         g_face = kf(k) / profile%dz_node(k)
         di_t(k)      = di_t(k) + g_face
         up_t(k)      = up_t(k) - g_face
         di_t(k+1)    = di_t(k+1) + g_face
         lo_t(k+1)    = lo_t(k+1) - g_face
         w_cp = forcing%w_flux(k) * rho_h2o * cp_liq
         if (forcing%w_flux(k) <= 0.0_wp) then                               ! down: layer k gives
            di_t(k)   = di_t(k) - w_cp
            lo_t(k+1) = lo_t(k+1) + w_cp
         else                                                                ! up: layer k+1 gives
            up_t(k)   = up_t(k) - w_cp
            di_t(k+1) = di_t(k+1) + w_cp
         end if
      end do
      if (forcing%w_flux_top > 0.0_wp) di_t(1) = di_t(1) + forcing%w_flux_top * rho_h2o * cp_liq
      di_t(n) = di_t(n) + g_deep

      !---~---
      !   Newton on the energies. Each step's Jacobian is the slopes above over each layer's
      ! apparent heat capacity, plus dz/dt on the diagonal; a step that does not lower the
      ! worst residual is halved, up to SOIL_STEP_HALVINGS times.
      !---~---
      e_new(1:n) = e_n(1:n)
      call column_balance(e_new, now)
      flux%iterations = 0_ik
      do it = 1_ik, SOIL_STEP_ITER_MAX
         if (now%norm <= SOIL_STEP_TOL) exit
         do k = 1_ik, n
            b(k) = profile%dz(k) / dt + di_t(k) / now%de(k)
            a(k) = 0.0_wp
            c(k) = 0.0_wp
            if (k > 1_ik) a(k) = lo_t(k) / now%de(k-1)
            if (k < n)    c(k) = up_t(k) / now%de(k+1)
            rhs(k) = -now%res(k)          ! a named array: an expression here is copied to a temporary
         end do
         call thomas_solve(a, b, c, rhs, step, n)
         lam = 1.0_wp
         do ib = 0_ik, SOIL_STEP_HALVINGS
            e_try(1:n) = e_new(1:n) + lam * step(1:n)
            call eval_column(e_try, try)
            norm_try = try%norm
            if (norm_try < now%norm) exit
            lam = 0.5_wp * lam
         end do
         e_new(1:n) = e_try(1:n)
         now = try
         flux%iterations = it
      end do

      !--- Conservative commit: each layer's energy from the faces at the solution.
      e0 = 0.0_wp
      e1 = 0.0_wp
      do k = 1_ik, n
         e0  = e0 + e_n(k) * profile%dz(k)
         col%soil_energy(k) = e_n(k) + dt / profile%dz(k) * ((now%hf(k) - now%hf(k-1)) + (now%qwf(k) - now%qwf(k-1))) &
                              + dt * q_src(k)
         e1  = e1 + col%soil_energy(k) * profile%dz(k)
      end do

      !--- Temperature and liquid fraction from the new energy.
      call soil_layer_temp(col%soil_energy(1:n), forcing%soil_water(1:n), profile%dry_heat_capacity(1:n),        &
                           profile%layer_curve(1:n), col%soil_temp(1:n), col%soil_fliq(1:n))

      !--- Diagnostics and the column's energy budget.
      flux%bottom_heat  = -now%hf(n)
      !---~---
      !   The column budget carries the boundary water enthalpy too; the interior faces
      ! telescope out of the sum exactly, as they must.
      !---~---
      flux%energy_resid = (e1 - e0) - dt * (forcing%g_top - flux%bottom_heat                   &
                          - sum(forcing%water_heat_out(1:n)) + (now%qwf(n) - now%qwf(0)))
      flux%heat_capacity = sum(soil_volumetric_heat_capacity(forcing%soil_water(1:n), fl_n(1:n),              &
                                                             profile%dry_heat_capacity(1:n)) * profile%dz(1:n))

   contains

      !---~---
      !   The column at energies e: each layer's temperature and apparent heat capacity on its
      ! freezing curve, the faces at those temperatures, each layer's balance residual
      ! [W/m2], and the worst residual as the temperature change [K] it would need.
      !---~---
      subroutine eval_column(e, ev)
         real(wp),                 intent(in)    :: e(:)
         type(soil_column_eval_t), intent(inout) :: ev
         real(wp)    :: fliq
         integer(ik) :: k
         do k = 1_ik, n
            call soil_layer_state(e(k), forcing%soil_water(k), profile%dry_heat_capacity(k), profile%layer_curve(k),&
                                  ev%temp(k), fliq, ev%de(k))
         end do
         call column_balance(e, ev)
      end subroutine eval_column

      !--- The faces and each layer's residual at ev's temperatures and apparent heat capacities.
      subroutine column_balance(e, ev)
         real(wp),                 intent(in)    :: e(:)
         type(soil_column_eval_t), intent(inout) :: ev
         integer(ik) :: k
         ev%hf(0) = -forcing%g_top
         do k = 1_ik, n - 1_ik
            ev%hf(k) = -kf(k) * (ev%temp(k) - ev%temp(k+1)) / profile%dz_node(k)
         end do
         ev%hf(n) = bottom_heat_face(soil, forcing, g_deep, ev%temp(n))
         call water_enthalpy_faces(forcing, ev%temp, n, ev%qwf)
         ev%norm = 0.0_wp
         do k = 1_ik, n
            ev%res(k) = profile%dz(k) * ((e(k) - e_n(k)) / dt - q_src(k))                                      &
                        - ((ev%hf(k) - ev%hf(k-1)) + (ev%qwf(k) - ev%qwf(k-1)))
            ev%norm = max(ev%norm, abs(ev%res(k)) * dt / (profile%dz(k) * ev%de(k)))
         end do
      end subroutine column_balance

   end subroutine soil_energy_step_implicit

   !---~---
   !   soil_energy_time_deriv -- the soil column's tendency dE_k/dt [W/m3]:
   ! soil_energy_step_implicit's flux divergence and sources with the faces at the current
   ! temperatures (docs/science/soil_biophysics.md eq. (10)). It changes no state;
   ! bottom_face is the bottom face's heat flux.
   !---~---
   pure subroutine soil_energy_time_deriv(col, forcing, profile, soil, dedt, bottom_face)
      type(soil_energy_state_t),   intent(in)  :: col
      type(soil_heat_env_t),      intent(in)  :: forcing
      type(soil_profile_t),        intent(in)  :: profile
      type(soil_params_t),         intent(in)  :: soil
      real(wp),                    intent(out) :: dedt(n_soil_layer_max)   !--- [W/m3] dE/dt per layer (0 for k>n)
      real(wp), optional,          intent(out) :: bottom_face   !--- [W/m2] bottom-face heat flux, positive up

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max)   :: t_n, fl_n, kappa, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: g_deep

      dedt = 0.0_wp
      n = profile%n_active

      !--- Each layer's temperature, liquid fraction and conductivity at its current energy.
      call layer_thermal_state(col, forcing, profile, n, t_n, fl_n, kappa)
      call face_conductivity(profile%dz, kappa, n, kf)

      !--- Conductive faces at the current temperatures (explicit).
      hf(0)  = -forcing%g_top
      do k = 1_ik, n - 1_ik
         hf(k) = -kf(k) * (t_n(k) - t_n(k+1)) / profile%dz_node(k)
      end do
      g_deep = deep_conductance(soil, profile, kappa(n), n)
      hf(n)  = bottom_heat_face(soil, forcing, g_deep, t_n(n))
      if (present(bottom_face)) bottom_face = hf(n)
      !--- Water-enthalpy faces at T^n, the same rule the implicit step applies at T^{n+1}.
      call water_enthalpy_faces(forcing, t_n, n, qwf)

      !--- dE_k/dt = flux divergence + source (q_src = -water_heat_out/dz).
      do k = 1_ik, n
         dedt(k) = ((hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))) / profile%dz(k)                        &
                   - forcing%water_heat_out(k) / profile%dz(k)
      end do
   end subroutine soil_energy_time_deriv

end module ed_soil_heat
