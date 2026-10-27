! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_soil_hydr_lib -- the soil water laws (ed_soil_hydr_lib): a layer's relations in
! the water network's unknown (docs/science/soil_plant_water.md §2), its retention curve
! and conductivity (docs/science/soil_biophysics.md eq. (2)), the pond's depression store,
! and the pedotransfer function (soil_biophysics.md eqs. (6)-(8)).
!   1. LAYERS (soil_water_relations, soil_theta_to_unknown, build_soil_layer_curve): two
!      van Genuchten soils and a Campbell soil; every derivative against a central
!      difference; the two branches of the unknown meet; theta and K match the read-offs;
!      theta -> unknown -> potential round-trips; above saturation the layer stores
!      through the specific storage.
!   2. FROZEN SOIL (soil_liquid_capacity, soil_ice_state): the freezing curve against the
!      retention curve at the freezing head; a layer's ice adds its water, impedes K and
!      fills the pores smoothly; a zero-ice record is a layer with no ice.
!   3. POND (pond_water_relations, pond_water_to_head): smooth at empty and at the
!      depression depth; water -> head round-trips.
!   4. PEDOTRANSFER (soil_texture_to_hydraulics): Hodnett & Tomasella (2002)'s worked
!      example, Ottoni et al. (2019)'s Ks, coarse fragments, the residual-water floor, and
!      the elemental call over a profile.
!   5. RELATIVE SATURATION (soil_theta_to_saturation), the output's soil-water read-off:
!      (theta - theta_res) / (theta_sat - theta_res), clamped to [0, 1].
!---~---
program test_soil_hydr_lib
   use ed_test_assert,     only : check_true, check_close, test_report
   use ed_kinds,           only : wp, ik
   use ed_constants,       only : grav, grav_head, rho_h2o, t_3ple, latent_heat_fusion
   use ed_soil_hydr_lib,   only : soil_layer_curve_t, build_soil_layer_curve, soil_water_relations,           &
                                  soil_theta_to_unknown, pond_water_relations, pond_water_to_head,            &
                                  SOIL_SPECIFIC_STORAGE, POND_DEPRESSION_DEPTH, SOIL_RETENTION_VG,            &
                                  SOIL_RETENTION_CAMPBELL, soil_head_to_theta, VG_PORE_CONNECTIVITY,          &
                                  soil_ice_t, soil_ice_state, soil_liquid_capacity, soil_texture_t,           &
                                  soil_texture_to_hydraulics, soil_theta_to_saturation
   implicit none

   call test_soil_water_relations()
   call test_frozen_soil()
   call test_pond_water_relations()
   call test_pedotransfer()
   call test_relative_saturation()

   call test_report('test_soil_hydr_lib')

contains

   !---~---
   !   The conductivity's closed form, written independently of ed_soil_hydr_lib: Mualem-
   ! van Genuchten Ks Se^L [1 - (1 - Se^(1/m))^m]^2, or Campbell's Ks Se^(2b+3).
   !---~---
   pure real(wp) function conductivity_ref(c, se) result(k)
      type(soil_layer_curve_t), intent(in) :: c
      real(wp),                 intent(in) :: se
      real(wp) :: m
      if (c%retention == SOIL_RETENTION_CAMPBELL) then
         k = c%ksat * se ** (2.0_wp * c%par_n + 3.0_wp)
      else
         m = 1.0_wp - 1.0_wp / c%par_n
         k = c%ksat * se ** VG_PORE_CONNECTIVITY * (1.0_wp - (1.0_wp - se ** (1.0_wp / m)) ** m) ** 2
      end if
   end function conductivity_ref

   !--- Relative agreement of an analytic derivative with a central difference.
   logical function fd_agrees(analytic, numeric, rtol, atol)
      real(wp), intent(in) :: analytic, numeric, rtol, atol
      fd_agrees = abs(analytic - numeric) <= rtol * max(abs(analytic), abs(numeric)) + atol
   end function fd_agrees

   !---~---
   !   1. Each layer's relations in the solver's unknown tau for BCI's van Genuchten clay
   ! (example 04), a van Genuchten sand and a Campbell loam, from very dry, through the branch
   ! point, to 2 m of head above saturation: every derivative against a central difference;
   ! the branches' join; theta and K against the read-offs and the closed forms where their
   ! floors do not bind; theta -> tau -> potential round trips; and a supersaturated layer.
   !---~---
   subroutine test_soil_water_relations()
      type(soil_layer_curve_t) :: c(3)
      type(soil_ice_t)         :: no_ice
      character(len=24)  :: label(3)
      real(wp) :: tau, h, psi, dpsi, w, dw, k, dk, psip, psim, wp_, wm, kp, km, dum1, dum2
      real(wp) :: theta, theta_ref, k_ref, tau_back, head
      real(wp) :: psi_lo, psi_hi, dpsi_lo, dpsi_hi, w_lo, w_hi, dw_lo, dw_hi
      integer  :: ic, i
      logical  :: ok_psi, ok_w, ok_k, ok_val, ok_k_val, ok_round
      !--- BCI's clay (example 04), a sand, and a Campbell loam.
      call build_soil_layer_curve(SOIL_RETENTION_VG, 0.4314_wp, 0.0981_wp, 0.7607_wp, 1.3554_wp, 4.54e-6_wp, 0.05_wp, c(1))
      call build_soil_layer_curve(SOIL_RETENTION_VG, 0.38_wp, 0.05_wp, 3.5_wp, 2.7_wp, 8.0e-5_wp, 0.10_wp, c(2))
      call build_soil_layer_curve(SOIL_RETENTION_CAMPBELL, 0.45_wp, 0.0_wp, -0.25_wp, 6.0_wp, 5.0e-6_wp, 0.20_wp, c(3))
      label = [character(len=24) :: 'vG clay (BCI)', 'vG sand', 'Campbell loam']
      do ic = 1, 3
         ok_psi = .true.
         ok_w = .true.
         ok_k = .true.
         !--- tau from very dry, through the branch point, to 2 m of head above saturation.
         do i = 0, 400
            tau = 1.0e-4_wp * (c(ic)%se_switch / 1.0e-4_wp) ** (real(i, wp) / 300.0_wp)
            if (i > 300) tau = c(ic)%se_switch + real(i - 300, wp) / 100.0_wp                       &
                               * c(ic)%slope_switch * (2.0_wp - c(ic)%head_switch)
            if (abs(tau - c(ic)%se_switch) < 1.0e-6_wp) cycle            ! straddles the join
            if (abs(tau - (c(ic)%se_switch - c(ic)%slope_switch * c(ic)%head_switch)) < 1.0e-6_wp) cycle
            h = 1.0e-7_wp * max(abs(tau), 1.0e-3_wp)
            call soil_water_relations(c(ic), no_ice, tau, psi, dpsi, w, dw, k, dk)
            call soil_water_relations(c(ic), no_ice, tau + h, psip, dum1, wp_, dum2, kp, dum1)
            call soil_water_relations(c(ic), no_ice, tau - h, psim, dum1, wm, dum2, km, dum1)
            if (.not. fd_agrees(dpsi, (psip - psim) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok_psi = .false.
            if (.not. fd_agrees(dw,   (wp_ - wm) / (2*h),    1.0e-4_wp, 1.0e-9_wp))  ok_w = .false.
            !---~---
            !   Mualem's 1 - (1 - Se^(1/m))^m is the difference of two numbers within Se^(1/m) of
            ! one: where Se^(1/m) < 1e-4, K carries round-off of order 1e-16 / Se^(1/m), which a
            ! finite difference amplifies past the tolerance (with L = -1, K there is above the
            ! 1e-20 floor). The derivative is checked where a difference can resolve it.
            !---~---
            if (c(ic)%retention == SOIL_RETENTION_VG .and. tau < c(ic)%se_switch) then
               if (tau ** (c(ic)%par_n / (c(ic)%par_n - 1.0_wp)) < 1.0e-4_wp) cycle
            end if
            if (.not. fd_agrees(dk,   (kp - km) / (2*h),     1.0e-4_wp, 1.0e-20_wp)) ok_k = .false.
         end do
         call check_true(trim(label(ic))//': dpsi/dtau matches finite differences', ok_psi)
         call check_true(trim(label(ic))//': dW/dtau matches finite differences', ok_w)
         call check_true(trim(label(ic))//': dK/dtau matches finite differences', ok_k)
         !---~---
         !   The two branches meet with the same value and slope. Campbell's branch point is
         ! its air-entry head, where the curve reaches saturation with a finite slope and the
         ! layer then holds theta_sat: its water's slope steps there from the curve's to the
         ! storage term's, by construction, so that slope's continuity is van Genuchten's
         ! alone; the potential and its slope meet for all three soils.
         !---~---
         tau = c(ic)%se_switch
         call soil_water_relations(c(ic), no_ice, tau - 1.0e-9_wp, psi_lo, dpsi_lo, w_lo, dw_lo, dum1, dum2)
         call soil_water_relations(c(ic), no_ice, tau + 1.0e-9_wp, psi_hi, dpsi_hi, w_hi, dw_hi, dum1, dum2)
         call check_true(trim(label(ic))//': branches meet (potential)', abs(psi_hi - psi_lo) < 1.0e-8_wp)
         call check_true(trim(label(ic))//': branches meet (slope of potential)',                     &
                         abs(dpsi_hi - dpsi_lo) <= 1.0e-5_wp * abs(dpsi_lo))
         if (c(ic)%retention == SOIL_RETENTION_VG)                                                     &
            call check_true(trim(label(ic))//': branches meet (slope of water)',                       &
                            abs(dw_hi - dw_lo) <= 1.0e-5_wp * abs(dw_lo))
         !--- Values against ed_soil_hydr_lib, where its floors do not bind.
         ok_val = .true.
         ok_k_val = .true.
         ok_round = .true.
         do i = 1, 60
            head = -10.0_wp ** (-2.0_wp + 5.0_wp * real(i, wp) / 60.0_wp)          ! -1 cm .. -1 km
            theta_ref = soil_head_to_theta(c(ic)%retention, head, c(ic)%theta_sat, c(ic)%theta_res,   &
                                            c(ic)%par_a, c(ic)%par_n)
            k_ref = conductivity_ref(c(ic), (theta_ref - c(ic)%theta_res) / (c(ic)%theta_sat - c(ic)%theta_res))
            if ((theta_ref - c(ic)%theta_res) / (c(ic)%theta_sat - c(ic)%theta_res) < 1.0e-8_wp) cycle
            tau = soil_theta_to_unknown(c(ic), no_ice, theta_ref)
            call soil_water_relations(c(ic), no_ice, tau, psi, dpsi, w, dw, k, dk)
            theta = w / (rho_h2o * c(ic)%dz)
            if (abs(theta - theta_ref) > 1.0e-10_wp) ok_val = .false.
            if (k_ref > 1.0e-12_wp .and. abs(k - k_ref) > 1.0e-8_wp * k_ref) ok_k_val = .false.
            if (abs(psi / grav_head - head) > 1.0e-6_wp * abs(head) .and.                              &
                theta_ref < c(ic)%theta_sat - 1.0e-9_wp) ok_round = .false.
         end do
         call check_true(trim(label(ic))//': theta(psi) matches ed_soil_hydr_lib', ok_val)
         call check_true(trim(label(ic))//': K(theta) matches the closed form', ok_k_val)
         call check_true(trim(label(ic))//': theta -> unknown -> potential round-trips', ok_round)
         !--- Above saturation the layer holds theta_sat + S_s * head.
         tau = soil_theta_to_unknown(c(ic), no_ice, c(ic)%theta_sat + SOIL_SPECIFIC_STORAGE * 0.5_wp)
         call soil_water_relations(c(ic), no_ice, tau, psi, dpsi, w, dw, k, dk)
         call check_close(trim(label(ic))//': head above saturation', psi / grav_head, 0.5_wp, 1.0e-9_wp)
         call check_close(trim(label(ic))//': conductivity above saturation', k, c(ic)%ksat, 0.0_wp)
         tau_back = soil_theta_to_unknown(c(ic), no_ice, w / (rho_h2o * c(ic)%dz))
         call check_close(trim(label(ic))//': supersaturated round trip', tau_back, tau, 1.0e-9_wp)
      end do
   end subroutine test_soil_water_relations

   !---~---
   !   2. Frozen soil: the freezing curve against the retention curve at the freezing head; a
   ! layer's ice adds its water, impedes K by 10^(-6 theta_ice/theta_sat) and fills the pores
   ! smoothly; a zero-ice record is a layer with no ice; theta -> tau -> theta round-trips.
   !---~---
   subroutine test_frozen_soil()
      type(soil_layer_curve_t) :: c(3)
      type(soil_ice_t)         :: ice, no_ice, zero_ice
      character(len=24)        :: label(3)
      real(wp) :: theta_ice, tau, h, psi, dpsi, w, dw, k, dk, psip, psim, wp_, wm, kp, km, dum1, dum2
      real(wp) :: psi0, dpsi0, w0, dw0, k0, dk0, psi_lo, dpsi_lo, w_lo, psi_hi, dpsi_hi, w_hi
      real(wp) :: temp, theta_liq, dtheta, theta_p, theta_m, head, theta_ref, tau_back, theta, tau_zero
      integer  :: ic, i
      logical  :: ok_fd, ok_add, ok_same, ok_round
      call build_soil_layer_curve(SOIL_RETENTION_VG, 0.4314_wp, 0.0981_wp, 0.7607_wp, 1.3554_wp, 4.54e-6_wp, 0.05_wp, c(1))
      call build_soil_layer_curve(SOIL_RETENTION_VG, 0.38_wp, 0.05_wp, 3.5_wp, 2.7_wp, 8.0e-5_wp, 0.10_wp, c(2))
      call build_soil_layer_curve(SOIL_RETENTION_CAMPBELL, 0.45_wp, 0.0_wp, -0.25_wp, 6.0_wp, 5.0e-6_wp, 0.20_wp, c(3))
      label = [character(len=24) :: 'vG clay (BCI)', 'vG sand', 'Campbell loam']
      theta_ice = 0.15_wp
      do ic = 1, 3
         !---~---
         !   The freezing curve: saturation at and above t_3ple; below, the retention curve at
         ! the freezing head L_f (T - t_3ple) / (g T), with its slope.
         !---~---
         call soil_liquid_capacity(c(ic), t_3ple + 0.5_wp, theta_liq, dtheta)
         call check_true(trim(label(ic))//': thawed layer holds saturation', theta_liq == c(ic)%theta_sat      &
                         .and. dtheta == 0.0_wp)
         temp = 263.15_wp
         head = latent_heat_fusion * (temp - t_3ple) / (grav * temp)
         theta_ref = soil_head_to_theta(c(ic)%retention, head, c(ic)%theta_sat, c(ic)%theta_res, c(ic)%par_a,   &
                                        c(ic)%par_n)
         call soil_liquid_capacity(c(ic), temp, theta_liq, dtheta)
         call check_close(trim(label(ic))//': -10 C liquid = retention at the freezing head', theta_liq,        &
                          theta_ref, 1.0e-10_wp)
         call soil_liquid_capacity(c(ic), temp + 1.0e-4_wp, theta_p, dum1)
         call soil_liquid_capacity(c(ic), temp - 1.0e-4_wp, theta_m, dum1)
         call check_true(trim(label(ic))//': freezing curve slope matches finite differences',                  &
                         fd_agrees(dtheta, (theta_p - theta_m) / 2.0e-4_wp, 1.0e-5_wp, 1.0e-12_wp))

         !--- A zero-ice record is a layer with no ice, to the bit.
         zero_ice = soil_ice_state(c(ic), 0.0_wp)
         ok_same = .true.
         do i = 1, 50
            tau = c(ic)%se_switch * real(i, wp) / 40.0_wp
            call soil_water_relations(c(ic), no_ice, tau, psi0, dpsi0, w0, dw0, k0, dk0)
            call soil_water_relations(c(ic), zero_ice, tau, psi, dpsi, w, dw, k, dk)
            if (psi /= psi0 .or. dpsi /= dpsi0 .or. w /= w0 .or. dw /= dw0 .or. k /= k0 .or. dk /= dk0)     &
               ok_same = .false.
         end do
         call check_true(trim(label(ic))//': zero ice = no ice', ok_same)

         !---~---
         !   With ice: the impedance; below the pore-full point the liquid's relations plus the
         ! ice's water and the impedance on K; every slope against finite differences across
         ! that point and into the pore-full branch.
         !---~---
         ice = soil_ice_state(c(ic), theta_ice)
         call check_close(trim(label(ic))//': ice impedance 10^(-6 theta_ice/theta_sat)', ice%impedance,          &
                          10.0_wp ** (-6.0_wp * theta_ice / c(ic)%theta_sat), 1.0e-14_wp)
         ok_fd = .true.
         ok_add = .true.
         tau_zero = ice%tau_full - ice%head_full / ice%dhead_full      ! where the full layer's head passes zero
         do i = 1, 300
            tau = 0.2_wp * ice%tau_full + (2.0_wp * tau_zero - 1.2_wp * ice%tau_full) * real(i, wp) / 300.0_wp
            if (abs(tau - ice%tau_full) < 1.0e-6_wp * ice%tau_full) cycle
            if (abs(tau - tau_zero) < 1.0e-6_wp * tau_zero) cycle
            if (abs(tau - c(ic)%se_switch) < 1.0e-6_wp) cycle
            h = 1.0e-7_wp * max(abs(tau), 1.0e-3_wp)
            call soil_water_relations(c(ic), ice, tau, psi, dpsi, w, dw, k, dk)
            call soil_water_relations(c(ic), ice, tau + h, psip, dum1, wp_, dum2, kp, dum1)
            call soil_water_relations(c(ic), ice, tau - h, psim, dum1, wm, dum2, km, dum1)
            if (.not. fd_agrees(dpsi, (psip - psim) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok_fd = .false.
            if (.not. fd_agrees(dw,   (wp_ - wm) / (2*h),    1.0e-4_wp, 1.0e-9_wp))  ok_fd = .false.
            if (.not. fd_agrees(dk,   (kp - km) / (2*h),     1.0e-4_wp, 1.0e-20_wp)) ok_fd = .false.
            if (tau < ice%tau_full) then
               call soil_water_relations(c(ic), no_ice, tau, psi0, dpsi0, w0, dw0, k0, dk0)
               if (abs(w - w0 - rho_h2o * c(ic)%dz * theta_ice) > 1.0e-12_wp * w .or. psi /= psi0 .or.       &
                   abs(k - k0 * ice%impedance) > 1.0e-14_wp * k) ok_add = .false.
            end if
         end do
         call check_true(trim(label(ic))//': frozen slopes match finite differences', ok_fd)
         call check_true(trim(label(ic))//': frozen = liquid + ice, K x impedance', ok_add)
         call soil_water_relations(c(ic), ice, ice%tau_full * (1.0_wp - 1.0e-10_wp), psi_lo, dpsi_lo, w_lo, dum1, &
                                   dum2, dum2)
         call soil_water_relations(c(ic), ice, ice%tau_full * (1.0_wp + 1.0e-10_wp), psi_hi, dpsi_hi, w_hi, dum1, &
                                   dum2, dum2)
         call check_true(trim(label(ic))//': pores fill smoothly (potential, slope, water)',                    &
                         abs(psi_hi - psi_lo) <= 1.0e-6_wp * abs(psi_lo) + 1.0e-12_wp .and.                    &
                         abs(dpsi_hi - dpsi_lo) <= 1.0e-5_wp * abs(dpsi_lo) .and.                              &
                         abs(w_hi - w_lo) <= 1.0e-9_wp * w_lo)
         call check_close(trim(label(ic))//': pores full at theta_sat', w_lo / (rho_h2o * c(ic)%dz),            &
                          c(ic)%theta_sat, 1.0e-8_wp)
         call soil_water_relations(c(ic), ice, 0.5_wp * (ice%tau_full + tau_zero), psi, dpsi, w, dw, k, dk)
         call check_true(trim(label(ic))//': a full layer holds theta_sat until its head passes zero',          &
                         abs(w / (rho_h2o * c(ic)%dz) - c(ic)%theta_sat) < 1.0e-12_wp .and. dw == 0.0_wp)
         !--- theta -> unknown -> theta round-trips, below and above the pore-full point.
         ok_round = .true.
         do i = 1, 20
            theta = theta_ice + c(ic)%theta_res + (c(ic)%theta_sat - c(ic)%theta_res - theta_ice)             &
                    * real(i, wp) / 20.0_wp + merge(SOIL_SPECIFIC_STORAGE, 0.0_wp, i == 20)
            tau_back = soil_theta_to_unknown(c(ic), ice, theta)
            call soil_water_relations(c(ic), ice, tau_back, psi, dpsi, w, dw, k, dk)
            if (abs(w / (rho_h2o * c(ic)%dz) - theta) > 1.0e-10_wp) ok_round = .false.
         end do
         call check_true(trim(label(ic))//': frozen theta -> unknown -> theta round-trips', ok_round)
      end do
   end subroutine test_frozen_soil

   !---~---
   !   3. The pond's relations: smooth at empty and at the depression depth, water -> head round
   ! trips, and the depression store's half filled at 5 cm of head.
   !---~---
   subroutine test_pond_water_relations()
      real(wp) :: h, psi, dpsi, w, dw, wpl, wm, dum, hh
      integer  :: i
      logical  :: ok_d, ok_r
      ok_d = .true.
      ok_r = .true.
      do i = -20, 60
         h = 1.0e-4_wp * real(i, wp) + 1.0e-7_wp               ! -2 mm .. 6 mm, off the corners
         call pond_water_relations(h, psi, dpsi, w, dw)
         hh = 1.0e-9_wp
         call pond_water_relations(h + hh, dum, dum, wpl, dum)
         call pond_water_relations(h - hh, dum, dum, wm, dum)
         if (.not. fd_agrees(dw, (wpl - wm) / (2*hh), 1.0e-5_wp, 1.0e-6_wp)) ok_d = .false.
         if (h > 0.0_wp .and. abs(pond_water_to_head(w) - h) > 1.0e-12_wp) ok_r = .false.
      end do
      call check_true('pond: dW/dh matches finite differences, smooth at empty and full', ok_d)
      call check_true('pond: water -> head round-trips', ok_r)
      call pond_water_relations(0.0_wp, psi, dpsi, w, dw)
      call check_close('pond: empty pond holds nothing', w, 0.0_wp, 0.0_wp)
      call pond_water_relations(0.05_wp, psi, dpsi, w, dw)
      call check_close('pond: 5 cm of head holds 50 kg/m2 less half the depression depth', w,           &
                       rho_h2o * (0.05_wp - 0.5_wp * POND_DEPRESSION_DEPTH), 1.0e-10_wp)
   end subroutine test_pond_water_relations

   !---~---
   !   4. The pedotransfer function (soil_texture_to_hydraulics; docs/science/soil_biophysics.md
   ! eqs. (6)-(8)): Hodnett & Tomasella (2002)'s worked example, Ottoni et al. (2019)'s Ks of
   ! it, coarse fragments, the residual-water floor, and the elemental call over a profile
   ! against the call per layer.
   !---~---
   subroutine test_pedotransfer()
      type(soil_texture_t)   :: example, stony, bare_sand, profile(3)
      real(wp)               :: theta_sat, theta_res, vg_alpha, vg_n, ksat
      real(wp)               :: theta_sat_stony, theta_res_stony, alpha_stony, n_stony, ksat_stony
      real(wp), dimension(3) :: theta_sat_p, theta_res_p, alpha_p, n_p, ksat_p
      real(wp)               :: theta_sat_k, theta_res_k, alpha_k, n_k, ksat_k
      logical                :: same
      integer                :: k
      print '(a)', 'test_pedotransfer:'

      !---~---
      !   Hodnett & Tomasella (2002) Table 8's worked example, as Marthews et al. (2014, GMD 7:711)
      ! repeat it: alpha 0.2452 per kPa, n 1.3674, theta_sat 0.4993, theta_res 0.2134.
      !---~---
      example = soil_texture_t(sand=39.2_wp, silt=24.2_wp, clay=36.7_wp, org_carbon=14.0_wp, cec=15.9_wp,  &
                               ph=5.85_wp, bulk_density=1.2_wp, coarse_frag=0.0_wp)
      call soil_texture_to_hydraulics(example, theta_sat, theta_res, vg_alpha, vg_n, ksat)
      call check_close('worked example: theta_sat', theta_sat, 0.4993_wp, 1.0e-4_wp)
      call check_close('worked example: theta_res', theta_res, 0.2134_wp, 1.0e-4_wp)
      call check_close('worked example: alpha per kPa', vg_alpha / (grav_head*1000.0_wp), 0.2452_wp, 1.0e-4_wp)
      call check_close('worked example: n', vg_n, 1.3674_wp, 1.0e-4_wp)

      !---~---
      !   Ottoni et al. (2019) for the example: 19.31 m/day x (theta_sat - theta(-33 kPa))^1.948,
      ! with theta(-33 kPa) = 0.34409; 5.9343e-6 m/s (21.4 mm/h), computed apart from ED-MORE.
      !---~---
      call check_close('worked example: Ks [m/s]', ksat, 5.934313e-6_wp, 1.0e-11_wp)

      !---~---
      !   30 % coarse fragments by volume at a fine-earth bulk density of 1.2: the curve holds 0.7 of
      ! the fine earth's water; the fragments are 0.48624 of the mass, so Ks keeps 0.60870 of the
      ! fine earth's (Saxton & Rawls 2006 eqs. 19 and 22). The curve's shape does not change.
      !---~---
      stony = example
      stony%coarse_frag = 30.0_wp
      call soil_texture_to_hydraulics(stony, theta_sat_stony, theta_res_stony, alpha_stony, n_stony, ksat_stony)
      call check_close('coarse fragments: theta_sat', theta_sat_stony, 0.7_wp*theta_sat, 1.0e-12_wp)
      call check_close('coarse fragments: theta_res', theta_res_stony, 0.7_wp*theta_res, 1.0e-12_wp)
      call check_close('coarse fragments: Ks ratio', ksat_stony/ksat, 0.6086956522_wp, 1.0e-9_wp)
      call check_true('coarse fragments: alpha and n unchanged', alpha_stony == vg_alpha .and. n_stony == vg_n)

      !---~---
      !   A sandy, alkaline soil of little CEC: the residual-water regression is negative there
      ! (-0.0086), and the floor holds it at zero.
      !---~---
      bare_sand = soil_texture_t(sand=98.0_wp, silt=1.0_wp, clay=1.0_wp, org_carbon=1.0_wp, cec=0.5_wp,  &
                                 ph=9.5_wp, bulk_density=1.6_wp, coarse_frag=0.0_wp)
      call soil_texture_to_hydraulics(bare_sand, theta_sat_k, theta_res_k, alpha_k, n_k, ksat_k)
      call check_close('sandy alkaline soil: theta_res floored at 0', theta_res_k, 0.0_wp, 0.0_wp)
      call check_true('sandy alkaline soil: a curve that drains (theta_sat > theta_res, n > 1, Ks > 0)',  &
                      theta_sat_k > theta_res_k .and. n_k > 1.0_wp .and. ksat_k > 0.0_wp)

      !--- One elemental call over a profile gives exactly the per-layer results.
      profile = [example, stony, bare_sand]
      call soil_texture_to_hydraulics(profile, theta_sat_p, theta_res_p, alpha_p, n_p, ksat_p)
      same = .true.
      do k = 1, 3
         call soil_texture_to_hydraulics(profile(k), theta_sat_k, theta_res_k, alpha_k, n_k, ksat_k)
         same = same .and. theta_sat_p(k) == theta_sat_k .and. theta_res_p(k) == theta_res_k .and.     &
                alpha_p(k) == alpha_k .and. n_p(k) == n_k .and. ksat_p(k) == ksat_k
      end do
      call check_true('elemental call over a profile equals the per-layer calls', same)
   end subroutine test_pedotransfer

   !---~---
   !   5. The relative saturation the output reports (soil_theta_to_saturation): half way
   ! between theta_res 0.1 and theta_sat 0.5 is 0.5; above saturation it is 1, below the
   ! residual water 0.
   !---~---
   subroutine test_relative_saturation()
      print '(a)', 'test_relative_saturation:'
      call check_close('relative saturation', soil_theta_to_saturation(0.3_wp, 0.5_wp, 0.1_wp), 0.5_wp, 1.0e-12_wp)
      call check_close('oversaturated clamps to 1', soil_theta_to_saturation(0.6_wp, 0.5_wp, 0.1_wp), 1.0_wp,       &
                       1.0e-12_wp)
      call check_close('below residual clamps to 0', soil_theta_to_saturation(0.0_wp, 0.5_wp, 0.1_wp), 0.0_wp,      &
                       1.0e-12_wp)
   end subroutine test_relative_saturation

end program test_soil_hydr_lib
