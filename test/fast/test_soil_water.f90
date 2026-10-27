! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_soil_water -- the soil's water flows (ed_soil_water;
! docs/science/soil_biophysics.md, docs/science/soil_plant_water.md §3).
!   1. DARCY (soil_water_flux, soil_biophysics.md eq. (3)): its four derivatives against
!      central differences; the pond face at zero head is the infiltration capacity.
!   2. EVAPORATION (ground_evaporation_rate, ground_evaporation_setup): its derivatives
!      against central differences, and the closed form where no dry surface layer forms.
!   3. SNOW-FREE SHARE (ground_evaporation_rate): evaporation is a tile flux, linear in
!      the snow-free share and zero under full cover.
!   4. MOISTURE (ground_evaporation_rate): a wetter soil evaporates more, approaching the
!      limit with no dry surface layer continuously.
!   5. BOTTOM FACE (bottom_flux; soil_biophysics.md §Boundary conditions): a leaky base
!      drains f K(theta_N), free drainage and sealed bedrock at its ends exactly; the
!      aquifer's face is smooth through zero gradient, with derivatives that match finite
!      differences, and carries the rise of a dry column across the geometric-mean face.
!   6. FROZEN GROUND (ground_evaporation_setup, ground_evaporation_rate;
!      soil_biophysics.md "Frozen soil"): a frozen layer's liquid, at the content its
!      freezing curve holds, evaporates at about the ice's saturation humidity, from -1 to
!      -20 C: saturation is taken over liquid, since the freezing curve already puts the
!      liquid at the ice's vapour pressure.
!   7. THE AIR'S PRESSURE (ground_evaporation_setup): at 70 kPa the ground's saturation
!      humidity holds the same saturation vapour pressure as at 1 atm, and the evaporation
!      rises by the saturation curve's ratio at that pressure.
!   8. ICE ON THE AQUIFER FACE (bottom_flux; soil_biophysics.md "Frozen soil"): a frozen
!      bottom layer's ice impedes the face to the table in full, both ways; a draining
!      base counts it once, in K(theta_N).
!---~---
program test_soil_water
   use ed_test_assert,     only : check_true, check_close, test_report
   use ed_kinds,           only : wp, ik
   use ed_constants,       only : grav, grav_head, rho_h2o, r_wv, p_std, n_soil_layer_max, t_3ple, h2o_dry_ratio
   use ed_soil_hydr_lib,   only : SOIL_RETENTION_VG, soil_theta_to_head, soil_liquid_capacity, soil_ice_state,  &
                                  soil_ice_t
   use ed_soil_water,      only : ground_evaporation_rate, ground_evaporation_setup, soil_water_flux, ground_evap_env_t, &
                                  bottom_flux
   use ed_therm_lib,       only : sat_specific_humidity, sat_vapor_pressure, specific_humidity_to_vapor_pressure
   use ed_soil_params,     only : soil_profile_t, soil_params_t, SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use ed_soil_derive,     only : build_profile_hydraulics
   use ed_test_soil_grid, only : geometric_layer_bottoms
   implicit none

   call test_soil_water_flux()
   call test_evaporation()
   call test_snow_free_evap()
   call test_evap_moisture_response()
   call test_bottom_face()
   call test_frozen_ground_evap()
   call test_evap_at_pressure()
   call test_aquifer_face_ice()

   call test_report('test_soil_water')

contains

   !--- Relative agreement of an analytic derivative with a central difference.
   logical function fd_agrees(analytic, numeric, rtol, atol)
      real(wp), intent(in) :: analytic, numeric, rtol, atol
      fd_agrees = abs(analytic - numeric) <= rtol * max(abs(analytic), abs(numeric)) + atol
   end function fd_agrees

   !---~---
   !   1. Darcy flow across the face between two layers 0.4 m apart (the geometric-mean face),
   ! and the pond's face at zero head.
   !---~---
   subroutine test_soil_water_flux()
      real(wp) :: q, qa, qb, dqa, dqb, dka, dkb, dum, h, ka, kb, pa, pb, psi1, k1, ks, z1, q_ref
      logical  :: ok
      ka = 3.0e-7_wp
      kb = 8.0e-9_wp
      pa = -0.02_wp
      pb = -0.6_wp
      call soil_water_flux(ka, kb, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, q, dqa, dqb, dka, dkb)
      ok = .true.
      h = 1.0e-7_wp
      call soil_water_flux(ka, kb, 0.4_wp, pa + h, pb, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka, kb, 0.4_wp, pa - h, pb, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dqa, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call soil_water_flux(ka, kb, 0.4_wp, pa, pb + h, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka, kb, 0.4_wp, pa, pb - h, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dqb, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      h = 1.0e-6_wp * ka
      call soil_water_flux(ka + h, kb, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka - h, kb, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dka, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      h = 1.0e-6_wp * kb
      call soil_water_flux(ka, kb + h, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka, kb - h, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dkb, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call check_true('darcy: derivatives match finite differences', ok)
      !---~---
      !   The pond face at zero head: the infiltration capacity,
      ! sqrt(K1 Ksat) (1 + (-psi1)/(-z1)), times the water density.
      !---~---
      psi1 = -3.0_wp
      k1 = 2.0e-9_wp
      ks = 4.5e-6_wp
      z1 = -0.015_wp
      call soil_water_flux(ks, k1, 0.5_wp, 0.0_wp, grav_head * psi1, -z1, -z1, q, dum, dum, dum, dum)
      q_ref = sqrt(k1 * ks) * (1.0_wp + (-psi1) / (-z1)) * rho_h2o
      call check_close('darcy: pond face at zero head = infiltration capacity', q, q_ref, 1.0e-12_wp * q_ref)
   end subroutine test_soil_water_flux

   !---~---
   !   2. Ground evaporation from BCI's clay top layer, over water contents from 0.11 to 0.43, at
   ! the potential its retention curve gives.
   !---~---
   subroutine test_evaporation()
      type(soil_profile_t)   :: sp
      type(ground_evap_env_t) :: f
      type(soil_params_t)    :: o
      real(wp) :: theta, head, psi, e, det, dep, ep, em, dum, ref, h
      integer  :: i
      logical  :: ok_d, ok_v
      sp%n_active = 3
      sp%retention = SOIL_RETENTION_VG
      sp%dz(1:3) = 0.03_wp
      sp%theta_sat(1:3) = 0.4314_wp
      sp%theta_res(1:3) = 0.0981_wp
      sp%curve_par_a(1:3) = 0.7607_wp
      sp%curve_par_n(1:3) = 1.3554_wp
      sp%ksat(1:3) = 4.54e-6_wp
      f%t_ground = 300.0_wp
      f%q_air = 0.012_wp
      f%rho_air = 1.15_wp
      f%r_aero = 60.0_wp
      call ground_evaporation_setup(sp, o, f)
      ok_d = .true.
      ok_v = .true.
      do i = 1, 30
         theta = 0.10_wp + 0.011_wp * real(i, wp)
         head  = -(1.0_wp / 0.7607_wp) * (((theta - 0.0981_wp) / (0.4314_wp - 0.0981_wp))              &
                  ** (-1.0_wp / (1.0_wp - 1.0_wp / 1.3554_wp)) - 1.0_wp) ** (1.0_wp / 1.3554_wp)
         psi   = head * grav_head                      ! [MPa], the potential the kernel takes
         call ground_evaporation_rate(theta, psi, f, o, e, det, dep)
         h = 1.0e-7_wp
         call ground_evaporation_rate(theta + h, psi, f, o, ep, dum, dum)
         call ground_evaporation_rate(theta - h, psi, f, o, em, dum, dum)
         if (.not. fd_agrees(det, (ep - em) / (2*h), 1.0e-4_wp, 1.0e-14_wp)) ok_d = .false.
         h = 1.0e-5_wp * max(abs(psi), grav_head)
         call ground_evaporation_rate(theta, psi + h, f, o, ep, dum, dum)
         call ground_evaporation_rate(theta, psi - h, f, o, em, dum, dum)
         if (.not. fd_agrees(dep, (ep - em) / (2*h), 1.0e-4_wp, 1.0e-14_wp)) ok_d = .false.
         !---~---
         !   Moist enough that no dry surface layer forms: the closed form rho (alpha q_sat - q_air)
         ! / r_aero, with Philip's alpha = exp(psi g / (R_v T)).
         !---~---
         if (theta > o%dsl_theta_init * 0.4314_wp + 0.01_wp) then
            ref = f%rho_air * (exp(head * grav / (r_wv * f%t_ground)) * sat_specific_humidity(f%t_ground, p_std, 1.0_wp) &
                               - f%q_air) / f%r_aero
            if (abs(e - ref) > 1.0e-10_wp * abs(ref)) ok_v = .false.
         end if
      end do
      call check_true('evaporation: derivatives match finite differences', ok_d)
      call check_true('evaporation: the closed form where no dry surface layer forms', ok_v)
   end subroutine test_evaporation

   !--- A 10-layer loam column, 2 m deep.
   subroutine loam_params(params)
      type(soil_profile_t), intent(out) :: params
      real(wp) :: bottoms(n_soil_layer_max)
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, bottoms, 0.43_wp, 0.078_wp,               &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, params)
   end subroutine loam_params

   !---~---
   !   3. Ground evaporation is an AREA-weighted tile flux, E = snow_free_frac rho dq / (r_aero +
   ! r_soil): linear in the snow-free share and exactly 0 under full cover, with a dry surface
   ! layer present (theta_1 = 0.30 is below dsl_theta_init * theta_sat = 0.344), where dividing
   ! r_aero by (1 - snowfac) instead would not be linear.
   !---~---
   subroutine test_snow_free_evap()
      type(soil_profile_t)    :: params
      type(soil_params_t)     :: opts
      type(ground_evap_env_t) :: f
      real(wp) :: e_full, e_half, e_none, de1, de2, psi
      print '(a)', 'test_snow_free_evap:'
      call loam_params(params)
      f%t_ground = 305.0_wp
      f%q_air    = 0.001_wp
      f%rho_air  = 1.2_wp
      f%r_aero   = 50.0_wp
      call ground_evaporation_setup(params, opts, f)
      psi = grav_head * soil_theta_to_head(SOIL_RETENTION_VG, 0.30_wp, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      f%snow_free_frac = 1.0_wp
      call ground_evaporation_rate(0.30_wp, psi, f, opts, e_full, de1, de2)
      f%snow_free_frac = 0.5_wp
      call ground_evaporation_rate(0.30_wp, psi, f, opts, e_half, de1, de2)
      f%snow_free_frac = 0.0_wp
      call ground_evaporation_rate(0.30_wp, psi, f, opts, e_none, de1, de2)
      call check_true('bare ground evaporates (non-trivial case)', e_full > 1.0e-8_wp, e_full)
      call check_close('half snow cover halves E exactly (area weighting)', e_half, 0.5_wp * e_full,         &
                       1.0e-14_wp * max(e_full, 1.0e-12_wp))
      call check_close('full snow cover gives exactly zero soil evaporation', e_none, 0.0_wp, 1.0e-30_wp)
   end subroutine test_snow_free_evap

   !---~---
   !   4. A wetter soil evaporates MORE. The dry-surface-layer resistance d_sl/(D_v tau) is built on
   ! the air-dry water content's tortuosity, so it falls as the layer wets and the response is
   ! monotone, approaching the no-DSL limit continuously at theta_init. Built from the bulk
   ! layer's air-filled porosity instead, the tortuosity would fall faster than the layer thins
   ! and a wetter soil would evaporate less.
   !---~---
   subroutine test_evap_moisture_response()
      type(soil_profile_t)    :: params
      type(soil_params_t)     :: opts
      type(ground_evap_env_t) :: f
      real(wp) :: theta, e_prev, e_now, e_sat, e_below, theta_init
      integer  :: i
      logical  :: monotone
      print '(a)', 'test_evap_moisture_response:'
      call loam_params(params)
      f%t_ground = 300.0_wp
      f%q_air    = 0.005_wp
      f%rho_air  = 1.2_wp
      f%r_aero   = 100.0_wp
      call ground_evaporation_setup(params, opts, f)
      theta_init = opts%dsl_theta_init * params%theta_sat(1)
      monotone = .true.
      e_prev   = -1.0_wp
      do i = 1, 40
         theta = params%theta_res(1) + real(i, wp) / 40.0_wp * (theta_init - params%theta_res(1))
         e_now = loam_evap(theta, f, opts)
         if (e_now < e_prev - 1.0e-14_wp) monotone = .false.
         e_prev = e_now
      end do
      call check_true('ground evaporation increases monotonically with soil moisture', monotone, 0.0_wp)
      e_sat   = loam_evap(theta_init, f, opts)
      e_below = loam_evap(theta_init - 1.0e-3_wp * (theta_init - params%theta_res(1)), f, opts)
      call check_true('approaches the DSL-free limit continuously at theta_init',                       &
                      e_below > 0.8_wp * e_sat, e_below / max(e_sat, 1.0e-30_wp))
   end subroutine test_evap_moisture_response

   !--- Ground evaporation of the loam's top layer at water content th, under f set up for it.
   real(wp) function loam_evap(th, f, opts) result(e)
      real(wp),                intent(in) :: th
      type(ground_evap_env_t), intent(in) :: f
      type(soil_params_t),     intent(in) :: opts
      real(wp) :: de1, de2
      call ground_evaporation_rate(th, grav_head * soil_theta_to_head(SOIL_RETENTION_VG, th, 0.43_wp, 0.078_wp, &
                                   3.6_wp, 1.56_wp), f, opts, e, de1, de2)
   end function loam_evap

   !---~---
   !   5. The bottom face of the loam's deepest layer. A leaky base drains f K(theta_N), the
   ! factor's ends being free drainage and sealed bedrock bit for bit. The aquifer's face is
   ! head-driven against a saturated table at the column's base: smooth through zero gradient
   ! (hydrostatic equilibrium, the bottom node's head minus half its thickness), with both
   ! derivatives against finite differences, and a dry column's rise at -2 m of head crossing
   ! the geometric-mean face.
   !---~---
   subroutine test_bottom_face()
      type(soil_profile_t) :: params
      type(soil_params_t)  :: opts
      real(wp), parameter :: F_LEAK = 0.05_wp          ! the ED soil maps' factor
      real(wp), parameter :: K_N = 1.234e-7_wp, PSI_N = -0.37_wp * grav_head   ! [m/s], [MPa] a deepest layer's state
      real(wp) :: q, dq_dpsi, dq_dk, delta, psi_eq, q_m, q_p, dqp_m, dqp_p, dqk, h_psi, q_k
      integer(ik) :: n
      print '(a)', 'test_bottom_face:'
      call loam_params(params)
      n = 10_ik
      opts%water_bottom_bc = SOIL_BC_DRAIN
      opts%drain_factor = 1.0_wp
      call bottom_flux(params, opts, n, PSI_N, K_N, 1.0_wp, q, dq_dpsi, dq_dk)
      call check_close('f = 1 is free drainage: q = K(theta_N), bit for bit', q, K_N, 0.0_wp)
      call check_close('f = 1: dq/dK = 1 and dq/dpsi = 0', abs(dq_dk - 1.0_wp) + abs(dq_dpsi), 0.0_wp, 0.0_wp)
      opts%drain_factor = 0.0_wp
      call bottom_flux(params, opts, n, PSI_N, K_N, 1.0_wp, q, dq_dpsi, dq_dk)
      call check_close('f = 0 is sealed bedrock: q = 0 and no derivative, bit for bit',                 &
                       abs(q) + abs(dq_dk) + abs(dq_dpsi), 0.0_wp, 0.0_wp)
      opts%drain_factor = F_LEAK
      call bottom_flux(params, opts, n, PSI_N, K_N, 1.0_wp, q, dq_dpsi, dq_dk)
      call check_close('a leaky base drains f K(theta_N)', q, F_LEAK * K_N, 0.0_wp)
      call check_close('and its flux moves with K(theta_N) by f', dq_dk, F_LEAK, 0.0_wp)
      !--- The aquifer's face.
      opts%water_bottom_bc = SOIL_BC_AQUIFER
      delta  = 0.5_wp * params%dz(n)
      psi_eq = -delta * grav_head
      h_psi  = 1.0e-6_wp
      call bottom_flux(params, opts, n, psi_eq - h_psi, 1.0e-7_wp, 1.0_wp, q_m, dqp_m, dqk)
      call bottom_flux(params, opts, n, psi_eq + h_psi, 1.0e-7_wp, 1.0_wp, q_p, dqp_p, dqk)
      call check_close('AQUIFER: dq/dpsi is continuous through zero gradient', dqp_m, dqp_p, 1.0e-12_wp * dqp_p)
      call check_close('AQUIFER: dq/dpsi matches its finite difference', (q_p - q_m) / (2.0_wp * h_psi), dqp_p,  &
                       1.0e-6_wp * dqp_p)
      call bottom_flux(params, opts, n, -2.0_wp * grav_head, 1.0e-7_wp, 1.0_wp, q_m, dqp_m, dqk)
      call bottom_flux(params, opts, n, -2.0_wp * grav_head, 1.0e-7_wp * (1.0_wp + 1.0e-6_wp), 1.0_wp, q_k)
      call check_close('AQUIFER: dq/dK matches its finite difference', (q_k - q_m) / 1.0e-13_wp, dqk,          &
                       1.0e-4_wp * abs(dqk))
      call check_close('AQUIFER: rise at -2 m head crosses the geometric-mean face',                          &
                       q_m, sqrt(1.0e-7_wp * params%ksat(n)) * (-2.0_wp / delta + 1.0_wp), 1.0e-12_wp * abs(q_m))
   end subroutine test_bottom_face

   !---~---
   !   6. A frozen loam top layer (docs/science/soil_biophysics.md, "Frozen soil"). Its liquid sits
   ! at the content its freezing curve holds at the ground's temperature (soil_liquid_capacity),
   ! at the head its retention curve gives that content; the ice fills the rest of its pores,
   ! so no dry surface layer forms. Into dry air its evaporation is rho alpha q_sat / r_aero.
   ! On the freezing curve the liquid's vapour pressure is the ice's (Clausius-Clapeyron), so
   ! alpha q_sat must be the ice's saturation humidity, to the error of the curve's constant
   ! latent heat and of the two saturation fits: within a fifth of the liquid-ice gap from -1 to
   ! -20 C (a twentieth at -10 C). Taking q_sat over ice as well would count the ice twice and
   ! land the whole gap below.
   !---~---
   subroutine test_frozen_ground_evap()
      type(soil_profile_t)    :: params
      type(soil_params_t)     :: opts
      type(ground_evap_env_t) :: f
      real(wp), parameter :: T_BELOW(4) = [1.0_wp, 5.0_wp, 10.0_wp, 20.0_wp]   ! [K] below the triple point
      real(wp) :: theta_liq, dtheta, psi, e, de1, de2, q_evap, q_liq, q_ice, worst
      integer(ik) :: k
      print '(a)', 'test_frozen_ground_evap:'
      call loam_params(params)
      f%q_air   = 0.0_wp
      f%rho_air = 1.3_wp
      f%r_aero  = 100.0_wp
      worst = 0.0_wp
      do k = 1_ik, size(T_BELOW, kind=ik)
         f%t_ground = t_3ple - T_BELOW(k)
         call ground_evaporation_setup(params, opts, f)
         call soil_liquid_capacity(params%layer_curve(1), f%t_ground, theta_liq, dtheta)
         psi = grav_head * soil_theta_to_head(SOIL_RETENTION_VG, theta_liq, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
         call ground_evaporation_rate(params%theta_sat(1), psi, f, opts, e, de1, de2)
         q_evap = e * f%r_aero / f%rho_air
         q_liq  = sat_specific_humidity(f%t_ground, p_std, 1.0_wp)
         q_ice  = sat_specific_humidity(f%t_ground, p_std, 0.0_wp)
         print '(a,f6.2,a,3es13.5,a,f7.4)', '   T = ', f%t_ground, ' K: evaporating, liquid, ice humidity',      &
               q_evap, q_liq, q_ice, '; off ice by (gap)', (q_evap - q_ice) / (q_liq - q_ice)
         worst = max(worst, abs(q_evap - q_ice) / (q_liq - q_ice))
      end do
      call check_true('frozen ground: its liquid evaporates at the ice''s saturation humidity, '//               &
                      'within a fifth of the liquid-ice gap (-1 to -20 C)', worst < 0.2_wp, worst)
   end subroutine test_frozen_ground_evap

   !---~---
   !   7. The ground and the air over it at 70 kPa, about 3,000 m, against the same at 1 atm, at
   ! one density and conductance, the loam moist enough that no dry surface layer forms, into
   ! dry air. The ground's saturation humidity at the air's pressure holds the saturation vapour
   ! pressure, a function of the temperature alone, so evaporation rises by q_sat(p)/q_sat(p_std)
   ! = (p_std - (1 - eps) e_s)/(p - (1 - eps) e_s); taken at p_std it would not move.
   !---~---
   subroutine test_evap_at_pressure()
      type(soil_profile_t)    :: params
      type(soil_params_t)     :: opts
      type(ground_evap_env_t) :: f_std, f_low
      real(wp), parameter :: P_LOW = 70000.0_wp, THETA = 0.40_wp   ! [Pa], [m3/m3]
      real(wp) :: psi, e_std, e_low, de1, de2, es, ratio
      print '(a)', 'test_evap_at_pressure:'
      call loam_params(params)
      f_std%t_ground = 290.0_wp
      f_std%q_air    = 0.0_wp
      f_std%rho_air  = 1.2_wp
      f_std%r_aero   = 80.0_wp
      f_low = f_std
      f_low%press = P_LOW
      call ground_evaporation_setup(params, opts, f_std)
      call ground_evaporation_setup(params, opts, f_low)
      es = sat_vapor_pressure(f_std%t_ground)
      call check_close('70 kPa: the ground''s saturation humidity holds the saturation vapour pressure',         &
                       specific_humidity_to_vapor_pressure(f_low%q_sat, P_LOW), es, 1.0e-12_wp * es)
      call check_close('and so does it at 1 atm', specific_humidity_to_vapor_pressure(f_std%q_sat, p_std), es,      &
                       1.0e-12_wp * es)
      psi = grav_head * soil_theta_to_head(SOIL_RETENTION_VG, THETA, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      call ground_evaporation_rate(THETA, psi, f_std, opts, e_std, de1, de2)
      call ground_evaporation_rate(THETA, psi, f_low, opts, e_low, de1, de2)
      ratio = (p_std - (1.0_wp - h2o_dry_ratio) * es) / (P_LOW - (1.0_wp - h2o_dry_ratio) * es)
      call check_close('70 kPa: the ground evaporates more, by the saturation curve''s ratio at that pressure',   &
                       e_low / e_std, ratio, 1.0e-6_wp * ratio)
   end subroutine test_evap_at_pressure

   !---~---
   !   8. The loam's deepest layer frozen, its ice 0.15 of 0.43, which impedes its conductivity by
   ! I = 10^(-6 theta_ice/theta_sat) (soil_ice_state). The face to the table lies inside the
   ! layer, so the table's side takes the impedance too and the face, their geometric mean,
   ! carries I times the thawed layer's flux at the same liquid potential, drawing up from a
   ! dry layer and draining a wet one. Impeding the node alone would leave sqrt(I). A draining
   ! base passes f K(theta_N) of the impeded conductivity: the ice counted once.
   !---~---
   subroutine test_aquifer_face_ice()
      type(soil_profile_t) :: params
      type(soil_params_t)  :: opts
      type(soil_ice_t)     :: ice
      real(wp), parameter :: K_LIQ = 1.0e-7_wp, THETA_ICE = 0.15_wp   ! [m/s] the liquid's K, [m3/m3]
      real(wp), parameter :: HEAD(2) = [-2.0_wp, -0.01_wp]            ! [m] a dry layer, a wet one
      character(len=*), parameter :: WHAT(2) = [character(len=32) :: 'the rise into a dry layer', 'the drainage of a wet layer']
      real(wp) :: q_thaw, q_froz(2), dqp_thaw, dqp_froz, dqk
      integer(ik) :: n, k
      print '(a)', 'test_aquifer_face_ice:'
      call loam_params(params)
      n = 10_ik
      ice = soil_ice_state(params%layer_curve(n), THETA_ICE)
      opts%water_bottom_bc = SOIL_BC_AQUIFER
      do k = 1_ik, 2_ik
         call bottom_flux(params, opts, n, HEAD(k) * grav_head, K_LIQ, 1.0_wp, q_thaw, dqp_thaw, dqk)
         call bottom_flux(params, opts, n, HEAD(k) * grav_head, ice%impedance * K_LIQ, ice%impedance, q_froz(k),    &
                          dqp_froz, dqk)
         call check_close('AQUIFER, frozen: the ice impedes '//trim(WHAT(k))//' in full', q_froz(k),                 &
                          ice%impedance * q_thaw, 1.0e-12_wp * abs(ice%impedance * q_thaw))
         call check_close('AQUIFER, frozen: and the flux''s slope in the potential with it', dqp_froz,               &
                          ice%impedance * dqp_thaw, 1.0e-12_wp * abs(ice%impedance * dqp_thaw))
      end do
      call check_true('AQUIFER, frozen: the dry layer draws up and the wet one drains',                            &
                      q_froz(1) < 0.0_wp .and. q_froz(2) > 0.0_wp, q_froz(1))
      opts%water_bottom_bc = SOIL_BC_DRAIN
      opts%drain_factor = 1.0_wp
      call bottom_flux(params, opts, n, -0.37_wp * grav_head, ice%impedance * K_LIQ, ice%impedance, q_froz(1))
      call check_close('a draining frozen base drains the impeded K(theta_N), the ice counted once', q_froz(1),     &
                       ice%impedance * K_LIQ, 0.0_wp)
   end subroutine test_aquifer_face_ice

end program test_soil_water
