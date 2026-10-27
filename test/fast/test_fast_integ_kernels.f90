! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_fast_integ_kernels -- the integrator's building blocks, below the schemes
! (src/fast_dynamics/integrator; docs/science/numerical_scheme.md §2-§3, eqs. (3)-(4)).
!   1. THE SURFACE RIGHT-HAND SIDE (surface_derivs): the diagnosed leaf temperature zeroes
!      the linearized leaf balance it was solved from; a hand-computable leaf with no
!      latent and no longwave forcing; a falling leaf that passes its film to its wood is
!      charged nothing for it.
!   2. THE ASSEMBLER (march_derivs): every tendency is finite; its canopy-air part is
!      surface_derivs at the state's soil-top temperature; the canopy hands the canopy air
!      what it absorbs; its soil-heat part is soil_energy_time_deriv on the water network's
!      flows, every term live.
!   3. DOMAIN SAFETY (march_derivs): off [theta_res, theta_sat] the tendencies stay finite,
!      the water tendency is the network's rate to the bit, and the soil-energy tendency
!      responds to the extra water only through its temperature.
!   4. THE BACKWARD-EULER STAGE (backward_euler_stage, advance_plant_water): stable and
!      physical at 900 s; close to the explicit RK45 reference at 4 s; its surface Newton
!      keeps the canopy air physical and below saturation under a harsh noon and from a
!      99 %-saturated start.
!   5. THE STATE VECTOR (state_init, state_sub, zero_like, state_err_diff, state_accum,
!      state_extrap, unpack_column_state, state_to_array, state_entry_rules): every field,
!      every layer and every cohort, through every operation; what each leaves out on
!      purpose is asserted left out.
!   6. THE ERROR CONTROL (build_tol_set, build_error_control, state_wrms_grouped,
!      step_control_factor): the config's dials reach every tolerance group; the grouped
!      norm reads 1 at one tolerance unit in every entry it counts, leaves out what the
!      march does not integrate, and is not diluted by the cohorts; the step factor takes
!      the order and the limits from the error control.
!   7. THE WATER NETWORK'S STORES (set_water_state, water_derivs): on a column of a PFT 1
!      and a PFT 2 tree, each cohort's wood and leaves sit at the potential their water
!      has on its own PFT's curves, and hold that water.
! The scheme-level properties (observed orders, the adaptive march) are
! test_fast_integ_schemes'; whole steps are test_fast_integ_step's.
!---~---
program test_fast_integ_kernels
   use ed_kinds,               only : wp, ik
   use ed_constants,           only : stefan_boltzmann, cp_air, tiny_num, rho_h2o, h2o_dry_ratio, n_soil_layer_max
   use ed_config,              only : ed_config_t
   use ed_therm_lib,           only : cas_temp_to_enthalpy, cas_enthalpy_to_temp, sat_specific_humidity,           &
                                      sat_vapor_pressure, sat_vapor_pressure_temp_deriv, enthalpy_vapor,          &
                                      temp_to_internal_energy, internal_energy_liquid
   use ed_soil_heat,           only : soil_heat_env_t, soil_energy_time_deriv, soil_layer_temp
   use ed_patch_state_vars,    only : soil_energy_state_t
   use ed_site_state_vars,     only : site_surface_t
   use ed_fast_integ_types,    only : surface_state_t, surface_tend_t, march_state_t, march_tend_t, frozen_t,      &
                                      error_control_t, tol_set_t, fast_dynamic_state_t, GRP_ENTH, GRP_SHV,         &
                                      GRP_CO2, GRP_SE, GRP_THETA, state_init, state_sub, state_err_diff,           &
                                      zero_like, state_accum, state_extrap, state_to_array, state_entry_rules,     &
                                      ensure_output_capacity, ensure_forcing_capacity, WOOD, LEAF
   use ed_fast_integ_derivs,   only : surface_derivs, march_derivs, set_water_state, water_derivs
   use ed_fast_integ_backward_euler, only : backward_euler_stage, advance_plant_water
   use ed_fast_integ_rk45,   only : rk45_step, RK45_P_ORDER
   use ed_fast_integ_control, only : build_tol_set, build_error_control, state_wrms_grouped, step_control_factor
   use ed_fast_integ_update,   only : unpack_column_state
   use ed_plant_hydr_lib,      only : tissue_psi_to_water, tissue_water_to_psi
   use ed_fast_state,          only : static_state_fixture
   use ed_test_fast_fixture,   only : fast_column_t, make_fast_column, kernel_record, reset_column, set_noon_forcing, &
                                      frozen_record
   use ed_test_assert,         only : check, check_true, test_report
   implicit none

   type(fast_column_t)   :: col          !--- the fixture column: its config and its site's soil
   type(error_control_t) :: ec_default   !--- the type's defaults: no budget stop

   call make_fast_column(col)

   call test_leaf_closure()
   call test_leaf_analytic()
   call test_falling_leaf_film()
   call test_assembler()
   call test_domain_safety()
   call test_be_stable_and_accurate()
   call test_be_coupled()
   call test_be_near_saturation()
   call test_state_vector()
   call test_tolerance_set()
   call test_error_norm()
   call test_step_factor()
   call test_water_curves_per_pft()

   call test_report('test_fast_integ_kernels')

contains

   !---~---
   !   1. THE SURFACE RIGHT-HAND SIDE (surface_derivs) on the kernel record, its tissues with
   !    no heat capacity, transpiring at the canopy air's humidity (not held).
   !---~---

   !---~---
   !   The diagnosed leaf temperature zeroes the linearized leaf balance it was solved from, for
   ! each of three cohorts: shortwave and longwave in, sensible heat, the full vapour
   ! enthalpy of the transpiration (docs/science/vegetation_energy_dynamics.md eq. (2)) and
   ! the emission about the leaf's start-of-step temperature out.
   !---~---
   subroutine test_leaf_closure()
      integer(ik), parameter :: N = 3_ik
      type(march_state_t)   :: y
      type(frozen_t)        :: frozen
      type(surface_state_t) :: ys
      type(surface_tend_t)  :: f
      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat, dtl, lw_slope, le_slope, le_ref, resid, worst
      integer(ik) :: i
      print '(a)', 'test_leaf_closure:'
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      frozen%tissue%transp_held = .false.
      ys%cas_enthalpy = cas_temp_to_enthalpy(298.0_wp, 0.012_wp)
      ys%cas_shv      = 0.012_wp
      ys%cas_co2      = 410.0_wp
      call surface_derivs(ys, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, N, f)
      tcas   = cas_enthalpy_to_temp(ys%cas_enthalpy, ys%cas_shv)
      qcas   = ys%cas_shv
      qsat_c = sat_specific_humidity(tcas, frozen%cas%press)
      esat   = sat_vapor_pressure(tcas)
      dqdt   = h2o_dry_ratio * frozen%cas%press / max((frozen%cas%press - (1.0_wp - h2o_dry_ratio) * esat) ** 2,  &
                                                      tiny_num) * sat_vapor_pressure_temp_deriv(tcas)
      worst = 0.0_wp
      do i = 1_ik, N
         dtl      = f%leaf_temp(i) - tcas
         lw_slope = 4.0_wp * frozen%tissue%leaf_emiss(i) * stefan_boltzmann * frozen%tissue%t_leaf0(i) ** 3      &
                    * frozen%tissue%lai(i)
         le_slope = enthalpy_vapor(tcas) * frozen%cas%rho * frozen%tissue%g_transp_leaf(i) * dqdt
         le_ref   = enthalpy_vapor(tcas) * frozen%cas%rho * frozen%tissue%g_transp_leaf(i) * (qsat_c - qcas)
         resid = frozen%tissue%abs_sw(i) + frozen%tissue%abs_lw(i) - frozen%tissue%h_coeff_leaf(i) * dtl          &
                 - (le_ref + le_slope * dtl) - lw_slope * (f%leaf_temp(i) - frozen%tissue%t_leaf0(i))
         worst = max(worst, abs(resid))
      end do
      call check_true('surface_derivs: the diagnosed leaf temperature zeroes the linearized balance',              &
                      worst < 1.0e-9_wp, worst)
      call check_true('surface_derivs: leaf temperatures physical (270-330 K)',                                    &
                      minval(f%leaf_temp) > 270.0_wp .and. maxval(f%leaf_temp) < 330.0_wp, f%leaf_temp(1))
   end subroutine test_leaf_closure

   !---~---
   !   No latent (g_transp_leaf = 0) and no longwave forcing (abs_lw = 0): the leaf balances its
   ! shortwave against sensible heat and its emission about its start temperature t_leaf0,
   ! dtl = (abs_sw + lw_slope*(t_leaf0 - tcas)) / (h + lw_slope).
   !---~---
   subroutine test_leaf_analytic()
      type(march_state_t)   :: y
      type(frozen_t)        :: frozen
      type(surface_state_t) :: ys
      type(surface_tend_t)  :: f
      real(wp) :: tcas, lw_slope, expect, t_leaf0
      print '(a)', 'test_leaf_analytic:'
      call kernel_record(col%site, col%cfg%pft%hydraulics, 1_ik, y, frozen)
      frozen%tissue%transp_held = .false.
      frozen%tissue%g_transp_leaf(1) = 0.0_wp
      frozen%tissue%abs_lw(1) = 0.0_wp
      frozen%tissue%abs_sw(1) = 200.0_wp
      frozen%tissue%lai(1) = 1.5_wp
      frozen%tissue%h_coeff_leaf(1) = 2.0_wp * frozen%tissue%lai(1) * 0.03_wp * 1.2_wp * cp_air
      ys%cas_enthalpy = cas_temp_to_enthalpy(299.0_wp, 0.010_wp)
      ys%cas_shv      = 0.010_wp
      call surface_derivs(ys, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, 1_ik, f)
      tcas     = cas_enthalpy_to_temp(ys%cas_enthalpy, ys%cas_shv)
      t_leaf0  = frozen%tissue%t_leaf0(1)
      lw_slope = 4.0_wp * frozen%tissue%leaf_emiss(1) * stefan_boltzmann * t_leaf0 ** 3 * frozen%tissue%lai(1)
      expect   = tcas + (frozen%tissue%abs_sw(1) + lw_slope * (t_leaf0 - tcas))                                    &
                        / (frozen%tissue%h_coeff_leaf(1) + lw_slope)
      call check('surface_derivs: leaf T = tcas + (Rn + lw_slope (t_leaf0 - tcas))/(h + lw_slope)', f%leaf_temp(1), &
                 expect, 1.0e-10_wp)
   end subroutine test_leaf_analytic

   !---~---
   !   A falling leaf's film passed to its wood (Palo Verde's dry season): the leaf, nearly gone,
   ! hands its whole film to the wood as the cohort's film is re-shared by capacity. The giver
   ! is charged nothing for the water it passes, so the leaf balances its sensible heat
   ! against its emission about t_leaf0 alone; charged the rain's valuation, it sat 4e10 K off.
   !---~---
   subroutine test_falling_leaf_film()
      type(march_state_t)   :: y
      type(frozen_t)        :: frozen
      type(surface_state_t) :: ys
      type(surface_tend_t)  :: f
      real(wp) :: tcas, t_leaf0, lw_slope, expect
      print '(a)', 'test_falling_leaf_film:'
      call kernel_record(col%site, col%cfg%pft%hydraulics, 1_ik, y, frozen)
      frozen%tissue%transp_held = .false.
      frozen%tissue%lai(1) = 1.0e-20_wp
      frozen%tissue%abs_sw(1) = 0.0_wp
      frozen%tissue%abs_lw(1) = 0.0_wp
      frozen%tissue%h_coeff_leaf(1)  = 2.0_wp * frozen%tissue%lai(1) * 0.03_wp * 1.2_wp * cp_air
      frozen%tissue%g_transp_leaf(1) = 0.004_wp * frozen%tissue%lai(1)
      frozen%tissue%wai(1) = 0.5_wp
      frozen%tissue%h_coeff_w(1) = 10.0_wp
      frozen%tissue%wood_hcap_per_dt(1) = 1.0e4_wp / 900.0_wp
      frozen%tissue%t_leaf0(1) = 303.9_wp
      frozen%tissue%t_wood0(1) = 304.5_wp
      frozen%film%f_wet_c(1) = 1.0_wp
      frozen%film%intercept_leaf(1) = -1.07e-12_wp            ! the leaf's whole film, passed on
      frozen%film%transfer_lw(1)    =  1.07e-12_wp
      frozen%film%intercept_wood(1) =  1.07e-12_wp
      frozen%film%rain_liquid_enthalpy = internal_energy_liquid(297.0_wp)
      ys%cas_enthalpy = cas_temp_to_enthalpy(299.0_wp, 0.010_wp)
      ys%cas_shv      = 0.010_wp
      call surface_derivs(ys, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, 1_ik, f)
      tcas     = cas_enthalpy_to_temp(ys%cas_enthalpy, ys%cas_shv)
      t_leaf0  = frozen%tissue%t_leaf0(1)
      lw_slope = 4.0_wp * frozen%tissue%leaf_emiss(1) * stefan_boltzmann * t_leaf0 ** 3 * frozen%tissue%lai(1)
      expect   = tcas + lw_slope * (t_leaf0 - tcas) / (frozen%tissue%h_coeff_leaf(1) + lw_slope)
      call check_true('surface_derivs: a falling leaf passing its film to its wood is charged nothing for it',      &
                      abs(f%leaf_temp(1) - expect) < 1.0e-8_wp, f%leaf_temp(1) - expect)
   end subroutine test_falling_leaf_film

   !---~---
   !   2. THE ASSEMBLER (march_derivs): two cohorts on the kernel record, the water network's
   !    flows set nonzero so that every term of the wiring check is live: a downward drift
   !    through every face, drainage, infiltration, the uptake's enthalpy leaving each layer,
   !    and water pushed up into the pond from the top layer. A check that omits a term is
   !    only green while that term is zero.
   !---~---
   subroutine test_assembler()
      integer(ik), parameter :: N = 2_ik
      type(march_state_t)       :: y
      type(frozen_t)            :: frozen
      type(march_tend_t)        :: f
      type(surface_state_t)     :: y_stage
      type(surface_tend_t)      :: surf_tend
      type(soil_heat_env_t)     :: eforc
      type(soil_energy_state_t) :: se
      real(wp)    :: dedt(n_soil_layer_max), canopy_to_cas, tcas, t_ground, fliq
      integer(ik) :: k, nsl
      logical     :: finite
      print '(a)', 'test_assembler:'
      nsl = col%site%soil_profile%n_active
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      frozen%hydrology%w_flux_frozen(1:nsl) = 2.0e-8_wp
      frozen%hydrology%w_flux_frozen(nsl)   = 0.0_wp
      frozen%hydrology%drainage     = 2.0e-8_wp * rho_h2o
      frozen%hydrology%infiltration = frozen%hydrology%q_top * rho_h2o
      frozen%hydrology%exfil_enth   = 1.0e-7_wp * internal_energy_liquid(296.0_wp)
      do k = 1_ik, nsl
         frozen%roots%root_heat(k) = frozen%roots%uptake * frozen%roots%root_share(k) * internal_energy_liquid(295.0_wp)
      end do
      call march_derivs(y, frozen, col%site, N, nsl, f)

      finite = ieee_ok(f%d_cas_enthalpy) .and. ieee_ok(f%d_cas_shv) .and. ieee_ok(f%d_cas_co2)
      do k = 1_ik, nsl
         finite = finite .and. ieee_ok(f%dedt(k)) .and. ieee_ok(f%dtheta_dt(k))
      end do
      finite = finite .and. all(ieee_ok(f%d_leaf_water_mass(1:N))) .and. all(ieee_ok(f%d_wood_water_mass(1:N)))
      call check_true('march_derivs: every tendency finite', finite)

      !--- The canopy-air part is surface_derivs at the state's diagnosed soil-top temperature.
      call soil_layer_temp(y%soil_energy(1), y%theta(1), col%site%soil_profile%dry_heat_capacity(1),             &
                           col%site%soil_profile%layer_curve(1), t_ground, fliq)
      y_stage%cas_enthalpy = y%cas_enthalpy
      y_stage%cas_shv      = y%cas_shv
      y_stage%cas_co2      = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, t_ground, N, surf_tend)
      call check('march_derivs: its canopy-air enthalpy tendency is surface_derivs''', f%d_cas_enthalpy,          &
                 surf_tend%d_cas_enthalpy, 1.0e-12_wp)

      !---~---
      !   The canopy is energy-neutral. With no advected enthalpy and no heat store in this record,
      ! what the canopy absorbs (coh_rnet) must equal what it hands the canopy air: sensible heat
      ! plus the full enthalpy of the vapour it sheds. Handing the canopy air transp*h_vapour
      ! while the leaf pays only latent_heat_vap, and charging the difference to the soil, breaks
      ! this by the proxy (~30 W/m2 at 3 mm/day) and makes the soil pay the transpired water's
      ! liquid enthalpy twice.
      !---~---
      tcas = cas_enthalpy_to_temp(y_stage%cas_enthalpy, y_stage%cas_shv)
      canopy_to_cas = surf_tend%src_enth - surf_tend%h_ground - surf_tend%le_ground                              &
                      + surf_tend%cond * internal_energy_liquid(tcas)
      call check_true('march_derivs: the canopy transpires in this record (the next check is live)',               &
                      sum(surf_tend%transp_c(1:N)) > 1.0e-7_wp, sum(surf_tend%transp_c(1:N)))
      call check('march_derivs: the canopy hands the canopy air what it absorbs (no soil proxy)', canopy_to_cas,   &
                 surf_tend%coh_rnet, 1.0e-9_wp * max(abs(surf_tend%coh_rnet), 1.0_wp))

      !---~---
      !   The soil-heat part is soil_energy_time_deriv on the same surface coupling (g_top) and the
      ! water network's flows: the enthalpy leaving each layer with the uptake, the time-mean
      ! faces, the infiltration at the pond's temperature and the drainage at the bottom's.
      !---~---
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top = surf_tend%g_top
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         eforc%water_heat_out(k) = frozen%roots%root_heat(k)
         eforc%w_flux(k)         = -frozen%hydrology%w_flux_frozen(k)
      end do
      eforc%water_heat_out(1)   = eforc%water_heat_out(1) + frozen%hydrology%exfil_enth
      eforc%water_heat_out(nsl) = eforc%water_heat_out(nsl)                                                      &
                                  + frozen%hydrology%drainage * internal_energy_liquid(frozen%hydrology%t_bot)
      eforc%w_flux_top  = -frozen%hydrology%infiltration / rho_h2o
      eforc%t_water_top = frozen%hydrology%t_infil
      call soil_energy_time_deriv(se, eforc, col%site%soil_profile, col%site%soil, dedt)
      call check('march_derivs: its soil-heat tendency is soil_energy_time_deriv on the network''s flows',        &
                 maxval(abs(f%dedt(1:nsl) - dedt(1:nsl))), 0.0_wp, 1.0e-9_wp)
   end subroutine test_assembler

   !---~---
   !   3. DOMAIN SAFETY (march_derivs). An explicit stage's theta can sit above theta_sat:
   ! theta is never clamped (theta_ood_max reports the excursion). Soil psi and K evaluated
   ! there must not steer the trajectory, and they cannot: a stage's water tendency is the
   ! water network's steady rate over the step, which does not read the stage's theta. So
   ! the water tendency at theta_sat + 0.05, far beyond any excursion the suite produces,
   ! equals the one at theta_sat to the bit. The soil-energy tendency is a flux divergence
   ! over internal energy, which theta reaches only through the temperature it diagnoses
   ! and the self-clamped Kersten number: it should respond to the extra water, whose heat
   ! capacity is real, so it is asserted finite and responsive.
   !---~---
   subroutine test_domain_safety()
      integer(ik), parameter :: N = 2_ik
      real(wp),    parameter :: EPS_OOD = 0.05_wp
      type(march_state_t) :: y, y_hi, y_lo, y_sat, y_res
      type(frozen_t)      :: frozen
      type(march_tend_t)  :: f_sat, f_hi, f_res, f_lo
      integer(ik) :: k, nsl
      logical     :: finite
      print '(a)', 'test_domain_safety:'
      nsl = col%site%soil_profile%n_active
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      call state_init(y, N, nsl, y_sat)
      call state_init(y, N, nsl, y_hi)
      call state_init(y, N, nsl, y_res)
      call state_init(y, N, nsl, y_lo)
      do k = 1_ik, nsl
         y_sat%theta(k) = col%site%soil_profile%theta_sat(k)
         y_hi%theta(k)  = col%site%soil_profile%theta_sat(k) + EPS_OOD
         y_res%theta(k) = col%site%soil_profile%theta_res(k)
         y_lo%theta(k)  = max(col%site%soil_profile%theta_res(k) - EPS_OOD, 0.0_wp)
      end do
      call march_derivs(y_sat, frozen, col%site, N, nsl, f_sat)
      call march_derivs(y_hi,  frozen, col%site, N, nsl, f_hi)
      call march_derivs(y_res, frozen, col%site, N, nsl, f_res)
      call march_derivs(y_lo,  frozen, col%site, N, nsl, f_lo)
      finite = all(ieee_ok(f_hi%dtheta_dt(1:nsl))) .and. all(ieee_ok(f_hi%dedt(1:nsl)))                          &
               .and. all(ieee_ok(f_lo%dtheta_dt(1:nsl))) .and. all(ieee_ok(f_lo%dedt(1:nsl)))
      call check_true('march_derivs: finite for theta far outside [theta_res, theta_sat]', finite)
      call check('march_derivs: the water tendency above saturation is the saturated one, to the bit',             &
                 maxval(abs(f_hi%dtheta_dt(1:nsl) - f_sat%dtheta_dt(1:nsl))), 0.0_wp, 0.0_wp)
      call check('march_derivs: the water tendency below residual is the residual one, to the bit',                &
                 maxval(abs(f_lo%dtheta_dt(1:nsl) - f_res%dtheta_dt(1:nsl))), 0.0_wp, 0.0_wp)

      !---~---
      !   Two all-liquid states at 290 K's internal energy, one saturated and one 0.05 above: the
      ! soil-energy tendency's response is asserted where the temperature depends on the water,
      ! not on a phase plateau, where it would be the same to the bit for the wrong reason.
      !---~---
      call state_init(y, N, nsl, y_sat)
      call state_init(y, N, nsl, y_hi)
      do k = 1_ik, nsl
         y_sat%theta(k) = col%site%soil_profile%theta_sat(k)
         y_hi%theta(k)  = col%site%soil_profile%theta_sat(k) + EPS_OOD
         y_sat%soil_energy(k) = temp_to_internal_energy(col%site%soil_profile%dry_heat_capacity(k),                &
                                                        y_sat%theta(k) * rho_h2o, 290.0_wp, 1.0_wp)
         y_hi%soil_energy(k)  = y_sat%soil_energy(k)
      end do
      call march_derivs(y_sat, frozen, col%site, N, nsl, f_sat)
      call march_derivs(y_hi,  frozen, col%site, N, nsl, f_hi)
      call check_true('march_derivs: the soil-energy tendency above saturation is finite and responds through T', &
                      maxval(abs(f_hi%dedt(1:nsl) - f_sat%dedt(1:nsl))) > 0.0_wp                                  &
                      .and. maxval(abs(f_hi%dedt(1:nsl) - f_sat%dedt(1:nsl))) < 1.0e6_wp,                          &
                      maxval(abs(f_hi%dedt(1:nsl) - f_sat%dedt(1:nsl))))
      call check('march_derivs: the water tendency is the network''s on an all-liquid state too',                  &
                 maxval(abs(f_hi%dtheta_dt(1:nsl) - f_sat%dtheta_dt(1:nsl))), 0.0_wp, 0.0_wp)
   end subroutine test_domain_safety

   !---~---
   !   4. THE BACKWARD-EULER STAGE. be_euler_step is the gamma = 1 degraded form of ESDIRK2,
   ! built from the two kernels esdirk2_step composes: one backward_euler_stage, then the
   ! closed-form plant-water split over the whole step on that stage's own transpiration
   ! (weight 1, the one stage). It is a probe, not an oracle: it shares backward_euler_stage
   ! with the scheme it degrades, so it is independent in the tableau only, and it belongs
   ! in the test that degrades it, not in the shipped source.
   !---~---

   !---~---
   !   Stable and physical at 900 s, where an explicit step blows up, over 6 h; and close to the
   ! explicit RK45 reference at 4 s over 4 min, which shares only the right-hand side with the
   ! stage, so agreement checks the implicit machinery against it. The canopy is ventilated
   ! four times as strongly as the kernel record's, so the source held over the full 900 s
   ! does not out-run the venting into supersaturation (test_be_coupled is the harsh case).
   !---~---
   subroutine test_be_stable_and_accurate()
      integer(ik), parameter :: N = 2_ik
      type(march_state_t) :: y, yi, yr, ytmp
      type(frozen_t)      :: frozen
      real(wp)    :: tcas, dcas, dmass
      integer(ik) :: step, nsl
      logical     :: physical
      print '(a)', 'test_be_stable_and_accurate:'
      nsl = col%site%soil_profile%n_active
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      call ventilate(frozen, 4.0_wp)
      call state_init(y, N, nsl, yi)
      physical = .true.
      do step = 1_ik, 24_ik                    ! 24 * 900 s = 6 h
         call be_euler_step(yi, frozen, N, nsl, 900.0_wp, ytmp)
         call state_init(ytmp, N, nsl, yi)
         tcas = cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. yi%cas_shv > 0.0_wp
         physical = physical .and. all(yi%theta(1:nsl) > 0.0_wp) .and. all(yi%theta(1:nsl) < 0.6_wp)
         physical = physical .and. all(yi%leaf_water_mass(1:N) > 0.0_wp)
      end do
      call check_true('backward_euler_stage: stable and physical at dt = 900 s (6 h)', physical, tcas)

      call state_init(y, N, nsl, yi)
      do step = 1_ik, 60_ik                    ! 60 * 4 s = 4 min
         call be_euler_step(yi, frozen, N, nsl, 4.0_wp, ytmp)
         call state_init(ytmp, N, nsl, yi)
      end do
      call march_rk45(y, frozen, N, nsl, 4.0_wp, 60_ik, yr)
      dcas  = abs(cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv) - cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv))
      dmass = max(maxval(abs(yi%leaf_water_mass(1:N) - yr%leaf_water_mass(1:N))),                                 &
                  maxval(abs(yi%wood_water_mass(1:N) - yr%wood_water_mass(1:N))))
      call check_true('backward_euler_stage ~ RK45 at dt = 4 s: canopy-air temperature within 5e-2 K', dcas < 5.0e-2_wp, &
                      dcas)
      call check_true('backward_euler_stage ~ RK45 at dt = 4 s: plant water within 1e-4 kg/plant', dmass < 1.0e-4_wp, &
                      dmass)
   end subroutine test_be_stable_and_accurate

   !---~---
   !   The leaf <-> canopy-air Newton keeps the canopy air physical and below saturation under a
   ! harsh constant noon at 900 s with modest venting, where an uncoupled pass over-humidifies.
   ! The leaves transpire at the canopy air's humidity, as the stages do whenever transpiration
   ! is not held at the network's: that leaf-vapour feedback is the coupling the Newton solves.
   !---~---
   subroutine test_be_coupled()
      integer(ik), parameter :: N = 2_ik
      type(march_state_t) :: y, yc, ytmp
      type(frozen_t)      :: frozen
      real(wp)    :: tcas, qsat
      integer(ik) :: step, nsl
      print '(a)', 'test_be_coupled:'
      nsl = col%site%soil_profile%n_active
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      frozen%tissue%transp_held = .false.
      call state_init(y, N, nsl, yc)
      do step = 1_ik, 12_ik                    ! 12 * 900 s = 3 h of constant noon
         call be_euler_step(yc, frozen, N, nsl, 900.0_wp, ytmp)
         call state_init(ytmp, N, nsl, yc)
      end do
      tcas = cas_enthalpy_to_temp(yc%cas_enthalpy, yc%cas_shv)
      qsat = sat_specific_humidity(tcas, frozen%cas%press)
      call check_true('backward_euler_stage: the coupled canopy air stays physical (270-325 K) at 900 s',          &
                      tcas > 270.0_wp .and. tcas < 325.0_wp, tcas)
      call check_true('backward_euler_stage: the coupled canopy air stays below saturation',                        &
                      yc%cas_shv <= qsat + 1.0e-4_wp, yc%cas_shv - qsat)
   end subroutine test_be_coupled

   !---~---
   !   The stage's surface Newton is robust near saturation (its line search): from canopy air at
   ! 99 % of saturation under the harsh noon it stays physical and below saturation at 900 s;
   ! and at 4 s it converges to the coupled solution the RK45 reference reaches.
   !---~---
   subroutine test_be_near_saturation()
      integer(ik), parameter :: N = 2_ik
      type(march_state_t) :: y, yi, yr, ytmp
      type(frozen_t)      :: frozen
      real(wp)    :: tcas, qsat, dcas, worst_super
      integer(ik) :: step, nsl
      logical     :: physical
      print '(a)', 'test_be_near_saturation:'
      nsl = col%site%soil_profile%n_active
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      tcas = 297.0_wp
      qsat = sat_specific_humidity(tcas, frozen%cas%press)
      y%cas_enthalpy = cas_temp_to_enthalpy(tcas, 0.99_wp * qsat)
      y%cas_shv      = 0.99_wp * qsat
      call state_init(y, N, nsl, yi)
      physical = .true.
      worst_super = -1.0_wp
      do step = 1_ik, 12_ik
         call be_euler_step(yi, frozen, N, nsl, 900.0_wp, ytmp)
         call state_init(ytmp, N, nsl, yi)
         tcas = cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv)
         qsat = sat_specific_humidity(tcas, frozen%cas%press)
         worst_super = max(worst_super, yi%cas_shv - qsat)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp
      end do
      call check_true('backward_euler_stage: below saturation from a 99 %-saturated start (900 s)',                &
                      worst_super <= 1.0e-4_wp, worst_super)
      call check_true('backward_euler_stage: physical (270-325 K) near saturation', physical, tcas)

      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      call state_init(y, N, nsl, yi)
      do step = 1_ik, 60_ik                    ! 60 * 4 s = 4 min
         call be_euler_step(yi, frozen, N, nsl, 4.0_wp, ytmp)
         call state_init(ytmp, N, nsl, yi)
      end do
      call march_rk45(y, frozen, N, nsl, 4.0_wp, 60_ik, yr)
      dcas = abs(cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv) - cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv))
      call check_true('backward_euler_stage''s Newton ~ RK45 at dt = 4 s: canopy-air temperature within 5e-2 K',   &
                      dcas < 5.0e-2_wp, dcas)
   end subroutine test_be_near_saturation

   !--- One first-order backward-Euler column step (see section 4's header).
   subroutine be_euler_step(y, frozen, n, nsl, dt, y_out)
      type(march_state_t), intent(in)    :: y
      type(frozen_t),      intent(in)    :: frozen
      integer(ik),         intent(in)    :: n, nsl
      real(wp),            intent(in)    :: dt
      type(march_state_t), intent(inout) :: y_out
      type(surface_tend_t) :: surf_tend
      call backward_euler_stage(y, frozen, col%site, ec_default, n, nsl, dt, y_out, surf_tend)
      call advance_plant_water(y, frozen, n, dt, surf_tend%transp_c(1:n), y_out)
   end subroutine be_euler_step

   !--- The explicit reference: nstep fixed Cash-Karp steps of dt, the 5th-order solution kept.
   subroutine march_rk45(y0, frozen, n, nsl, dt, nstep, y_out)
      type(march_state_t), intent(in)    :: y0
      type(frozen_t),      intent(in)    :: frozen
      integer(ik),         intent(in)    :: n, nsl, nstep
      real(wp),            intent(in)    :: dt
      type(march_state_t), intent(inout) :: y_out
      type(march_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call state_init(y0, n, nsl, y)
      do s = 1_ik, nstep
         call rk45_step(y, frozen, col%site, n, nsl, dt, ytmp, yerr)
         call state_init(ytmp, n, nsl, y)
      end do
      call state_init(y, n, nsl, y_out)
   end subroutine march_rk45

   !--- Scale the canopy air's three conductances to the atmosphere.
   subroutine ventilate(frozen, factor)
      type(frozen_t), intent(inout) :: frozen
      real(wp),       intent(in)    :: factor
      frozen%cas%box%atm_conductance_enthalpy = factor * frozen%cas%box%atm_conductance_enthalpy
      frozen%cas%box%atm_conductance_vapor    = factor * frozen%cas%box%atm_conductance_vapor
      frozen%cas%box%atm_conductance_co2      = factor * frozen%cas%box%atm_conductance_co2
   end subroutine ventilate

   !---~---
   !   5. THE STATE VECTOR. Its fields are enumerated in four routines (state_to_array,
   ! array_to_state, tend_to_array, state_entry_rules) and the combinators work on the flat
   ! array, so an omitted field is invisible three ways over: the compiler is silent, the
   ! ledgers still close (a field left at its default is a consistent zero), and the step
   ! controller reads a smaller error and takes a longer step. Fortran cannot enumerate a
   ! type's components, so each field, at every layer and every cohort, gets a distinct
   ! value and each operation's contract is asserted field by field. A new state field means
   ! new rows in same_state, fill_state and fill_tend.
   !---~---
   subroutine test_state_vector()
      integer(ik), parameter :: N = 2_ik, NL = 3_ik
      integer(ik), parameter :: NV = 5_ik + 2_ik*NL + 4_ik*N   !--- the flat layout's length
      real(wp),    parameter :: BB = 0.25_wp
      type(march_state_t) :: a, b, out, expect
      type(march_tend_t)  :: k
      type(fast_dynamic_state_t) :: dyn
      real(wp)    :: v(NV)
      integer(ik) :: j, grp, n_norm, n_err
      logical     :: counted
      print '(a)', 'test_state_vector:'
      call fill_state(a, 1.0_wp)
      call fill_state(b, 0.25_wp)

      !--- state_init copies every field.
      call state_init(a, N, NL, out)
      call same_state('state_init', out, a)
      !--- ... into an output already sized, keeping its arrays.
      call state_init(b, N, NL, out)
      call same_state('state_init into a reused output', out, b)

      !--- state_sub differences every field: a - b = fill(1 - 0.25).
      call state_sub(a, b, N, NL, out)
      call fill_state(expect, 0.75_wp)
      call same_state('state_sub', out, expect)

      !---~---
      !   zero_like allocates and zeroes every field: a film left unallocated is a field the norm
      ! cannot read, so the allocation is asserted, not just the values.
      !---~---
      out = zero_like(a, N, NL)
      call check_true('zero_like: the tissue-water and film fields allocated',                                    &
                      allocated(out%leaf_water_mass) .and. allocated(out%wood_water_mass)                          &
                      .and. allocated(out%leaf_surf_water) .and. allocated(out%wood_surf_water))
      call fill_state(expect, 0.0_wp)
      call same_state('zero_like', out, expect)

      !---~---
      !   state_err_diff: (a - b) - (b - a) = 2 (a - b) in the canopy air and the soil heat, which
      ! the march integrates; zero in the soil water, the tissue water, the films and the pond,
      ! which move by the water network's or the films' maps and are left out on purpose.
      !---~---
      call state_err_diff(a, b, b, a, N, NL, out)
      call fill_state(expect, 1.5_wp)
      expect%theta = 0.0_wp
      expect%leaf_water_mass = 0.0_wp
      expect%wood_water_mass = 0.0_wp
      expect%leaf_surf_water = 0.0_wp
      expect%wood_surf_water = 0.0_wp
      expect%w_surface = 0.0_wp
      expect%w_surface_energy = 0.0_wp
      call same_state('state_err_diff', out, expect)

      !---~---
      !   state_accum adds a*k to every field that has a stage tendency, and leaves the pond, which
      ! has none (it is committed from the water network), as it was rather than zeroed.
      !---~---
      call fill_state(out, 1.0_wp)
      call fill_tend(k)
      call state_accum(out, 2.0_wp, k, N, NL)
      call fill_state(expect, 1.0_wp)
      expect%cas_enthalpy = expect%cas_enthalpy + 2.0_wp * k%d_cas_enthalpy
      expect%cas_shv      = expect%cas_shv      + 2.0_wp * k%d_cas_shv
      expect%cas_co2      = expect%cas_co2      + 2.0_wp * k%d_cas_co2
      expect%soil_energy(1:NL) = expect%soil_energy(1:NL) + 2.0_wp * k%dedt(1:NL)
      expect%theta(1:NL)       = expect%theta(1:NL)       + 2.0_wp * k%dtheta_dt(1:NL)
      expect%leaf_water_mass = expect%leaf_water_mass + 2.0_wp * k%d_leaf_water_mass
      expect%wood_water_mass = expect%wood_water_mass + 2.0_wp * k%d_wood_water_mass
      expect%leaf_surf_water = expect%leaf_surf_water + 2.0_wp * k%d_leaf_surf_water
      expect%wood_surf_water = expect%wood_surf_water + 2.0_wp * k%d_wood_surf_water
      call same_state('state_accum (the pond passed through)', out, expect)

      !--- state_extrap blends every field, the pond included: (1 - b) a + b b = fill(0.8125).
      call state_extrap(a, BB, b, N, NL, out)
      call fill_state(expect, (1.0_wp - BB) * 1.0_wp + BB * 0.25_wp)
      call same_state('state_extrap', out, expect)

      !---~---
      !   unpack_column_state, the commit: a field it omits never reaches the patch state at all.
      ! It writes 9 of the 11 fields, and leaves the pond, which each scheme commits from the
      ! water network, as it found it. can_temp is re-diagnosed from the committed enthalpy and
      ! humidity, whatever that is: the synthetic fill is not physical.
      !---~---
      allocate(dyn%leaf_water_mass(N), dyn%wood_water_mass(N), dyn%leaf_surf_water(N), dyn%wood_surf_water(N))
      dyn%leaf_water_mass = 0.0_wp
      dyn%wood_water_mass = 0.0_wp
      dyn%leaf_surf_water = 0.0_wp
      dyn%wood_surf_water = 0.0_wp
      dyn%soil_water%w_surface = -1.0_wp
      dyn%soil_water%w_surface_energy = -2.0_wp
      call unpack_column_state(a, N, NL, dyn)
      call check('unpack: cas_enthalpy committed', dyn%cas%can_enthalpy, a%cas_enthalpy, 0.0_wp)
      call check('unpack: cas_shv committed', dyn%cas%can_shv, a%cas_shv, 0.0_wp)
      call check('unpack: cas_co2 committed', dyn%cas%can_co2, a%cas_co2, 0.0_wp)
      call check('unpack: soil_energy committed, every layer',                                                    &
                 maxval(abs(dyn%soil_energy%soil_energy(1:NL) - a%soil_energy(1:NL))), 0.0_wp, 0.0_wp)
      call check('unpack: theta committed, every layer', maxval(abs(dyn%soil_water%theta(1:NL) - a%theta(1:NL))),  &
                 0.0_wp, 0.0_wp)
      call check('unpack: leaf_water_mass committed, every cohort',                                                &
                 maxval(abs(dyn%leaf_water_mass - a%leaf_water_mass)), 0.0_wp, 0.0_wp)
      call check('unpack: wood_water_mass committed, every cohort',                                                &
                 maxval(abs(dyn%wood_water_mass - a%wood_water_mass)), 0.0_wp, 0.0_wp)
      call check('unpack: leaf_surf_water committed, every cohort',                                                &
                 maxval(abs(dyn%leaf_surf_water - a%leaf_surf_water)), 0.0_wp, 0.0_wp)
      call check('unpack: wood_surf_water committed, every cohort',                                                &
                 maxval(abs(dyn%wood_surf_water - a%wood_surf_water)), 0.0_wp, 0.0_wp)
      call check('unpack: can_temp re-diagnosed from the committed enthalpy', dyn%cas%can_temp,                   &
                 cas_enthalpy_to_temp(a%cas_enthalpy, a%cas_shv), 1.0e-12_wp)
      call check_true('unpack: the pond left as it was (the water network commits it)',                           &
                      dyn%soil_water%w_surface == -1.0_wp .and. dyn%soil_water%w_surface_energy == -2.0_wp)

      !---~---
      !   The layout: state_to_array places each field where state_entry_rules reads it, and the
      ! rules count the canopy air, the soil and the plant water in the step-size norm's groups
      ! (3 + 2*NL + 2*N entries) and the canopy air and the soil heat in the embedded error.
      !---~---
      call state_to_array(a, N, NL, v)
      call check('state_to_array: the canopy air first',                                                          &
                 maxval(abs(v(1:3) - [a%cas_enthalpy, a%cas_shv, a%cas_co2])), 0.0_wp, 0.0_wp)
      call check('state_to_array: then the soil heat and water, layer by layer',                                  &
                 maxval(abs(v(4:3+2*NL) - [a%soil_energy(1:NL), a%theta(1:NL)])), 0.0_wp, 0.0_wp)
      call check('state_to_array: then the tissue water and the films, cohort by cohort',                         &
                 maxval(abs(v(4+2*NL:3+2*NL+4*N) - [(a%leaf_water_mass(j), a%wood_water_mass(j), j = 1_ik, N),     &
                                                     (a%leaf_surf_water(j), a%wood_surf_water(j), j = 1_ik, N)])),  &
                 0.0_wp, 0.0_wp)
      call check('state_to_array: the pond last', abs(v(NV-1) - a%w_surface) + abs(v(NV) - a%w_surface_energy),   &
                 0.0_wp, 0.0_wp)
      n_norm = 0_ik
      n_err = 0_ik
      do j = 1_ik, NV
         call state_entry_rules(j, N, NL, grp, counted)
         if (grp /= 0_ik) n_norm = n_norm + 1_ik
         if (counted)     n_err  = n_err  + 1_ik
      end do
      call check('state_entry_rules: the norm groups the canopy air, the soil and the plant water', real(n_norm, wp), &
                 real(3_ik + 2_ik*NL + 2_ik*N, wp), 0.0_wp)
      call check('state_entry_rules: the embedded error counts the canopy air and the soil heat', real(n_err, wp),  &
                 real(3_ik + NL, wp), 0.0_wp)
   end subroutine test_state_vector

   !---~---
   !   Every field, at every layer and every cohort, a distinct non-zero multiple of s, so an
   ! operation that copies the wrong field, or leaves one at its default, cannot pass by
   ! coincidence. The values are quarters, so the combinations above are exact.
   !---~---
   subroutine fill_state(y, s)
      type(march_state_t), intent(out) :: y
      real(wp),            intent(in)  :: s
      integer(ik), parameter :: N = 2_ik, NL = 3_ik
      integer(ik) :: j
      y%cas_enthalpy     = 11.0_wp * s
      y%cas_shv          = 12.0_wp * s
      y%cas_co2          = 13.0_wp * s
      y%w_surface        = 14.0_wp * s
      y%w_surface_energy = 15.0_wp * s
      y%soil_energy = 0.0_wp
      y%theta = 0.0_wp
      do j = 1_ik, NL
         y%soil_energy(j) = (20.0_wp + real(j, wp)) * s
         y%theta(j)       = (30.0_wp + real(j, wp)) * s
      end do
      allocate(y%leaf_water_mass(N), y%wood_water_mass(N), y%leaf_surf_water(N), y%wood_surf_water(N))
      do j = 1_ik, N
         y%leaf_water_mass(j) = (40.0_wp + real(j, wp)) * s
         y%wood_water_mass(j) = (50.0_wp + real(j, wp)) * s
         y%leaf_surf_water(j) = (60.0_wp + real(j, wp)) * s
         y%wood_surf_water(j) = (70.0_wp + real(j, wp)) * s
      end do
   end subroutine fill_state

   subroutine fill_tend(t)
      type(march_tend_t), intent(out) :: t
      integer(ik), parameter :: N = 2_ik, NL = 3_ik
      integer(ik) :: j
      t%d_cas_enthalpy = 1.5_wp
      t%d_cas_shv = 2.5_wp
      t%d_cas_co2 = 3.5_wp
      t%dedt = 0.0_wp
      t%dtheta_dt = 0.0_wp
      do j = 1_ik, NL
         t%dedt(j) = 4.0_wp + real(j, wp)
         t%dtheta_dt(j) = 5.0_wp + real(j, wp)
      end do
      allocate(t%d_leaf_water_mass(N), t%d_wood_water_mass(N), t%d_leaf_surf_water(N), t%d_wood_surf_water(N))
      do j = 1_ik, N
         t%d_leaf_water_mass(j) = 6.0_wp + real(j, wp)
         t%d_wood_water_mass(j) = 7.0_wp + real(j, wp)
         t%d_leaf_surf_water(j) = 8.0_wp + real(j, wp)
         t%d_wood_surf_water(j) = 9.0_wp + real(j, wp)
      end do
   end subroutine fill_tend

   !---~---
   !   x equals y in every field, every layer and every cohort (to round-off on the fill's
   ! magnitudes). Layers past the state's own are zero in both.
   !---~---
   subroutine same_state(tag, x, y)
      character(len=*),    intent(in) :: tag
      type(march_state_t), intent(in) :: x, y
      real(wp), parameter :: TOL = 1.0e-12_wp
      call check(tag//': cas_enthalpy', x%cas_enthalpy, y%cas_enthalpy, TOL)
      call check(tag//': cas_shv', x%cas_shv, y%cas_shv, TOL)
      call check(tag//': cas_co2', x%cas_co2, y%cas_co2, TOL)
      call check(tag//': soil_energy, every layer', maxval(abs(x%soil_energy - y%soil_energy)), 0.0_wp, TOL)
      call check(tag//': theta, every layer', maxval(abs(x%theta - y%theta)), 0.0_wp, TOL)
      call check(tag//': leaf_water_mass, every cohort', maxval(abs(x%leaf_water_mass - y%leaf_water_mass)),       &
                 0.0_wp, TOL)
      call check(tag//': wood_water_mass, every cohort', maxval(abs(x%wood_water_mass - y%wood_water_mass)),       &
                 0.0_wp, TOL)
      call check(tag//': leaf_surf_water, every cohort', maxval(abs(x%leaf_surf_water - y%leaf_surf_water)),       &
                 0.0_wp, TOL)
      call check(tag//': wood_surf_water, every cohort', maxval(abs(x%wood_surf_water - y%wood_surf_water)),       &
                 0.0_wp, TOL)
      call check(tag//': w_surface (the pond)', x%w_surface, y%w_surface, TOL)
      call check(tag//': w_surface_energy (the pond)', x%w_surface_energy, y%w_surface_energy, TOL)
   end subroutine same_state

   !--- 6. THE ERROR CONTROL (ed_fast_integ_control).

   !---~---
   !   One tolerance set for every adaptive error estimate in the fast loop: build_tol_set is the
   ! single authority, the schemes take their error control from it (build_error_control) and
   ! the water network its GRP_THETA pair. A whole-model comparison cannot see this wiring (a
   ! tolerance change that flips no sub-step decision is byte-identical), so it is asserted
   ! directly, the atol_scale dial included: rtol alone saturates once atol dominates the
   ! atol + rtol|y| denominator.
   !---~---
   subroutine test_tolerance_set()
      type(ed_config_t)     :: c
      type(tol_set_t)       :: t, t1
      type(error_control_t) :: ec
      print '(a)', 'test_tolerance_set:'
      c = col%cfg
      t1 = build_tol_set(c)
      call check_true('build_tol_set: every group takes [fast].rtol', all(t1%rtol == c%rtol), maxval(abs(t1%rtol - c%rtol)))
      call check_true('build_tol_set: the theta group takes [fast].water_atol', t1%atol(GRP_THETA) == c%water_atol,  &
                      t1%atol(GRP_THETA))
      c%rtol       = 1.0e-7_wp
      c%atol_scale = 1.0e-2_wp
      t = build_tol_set(c)
      call check_true('build_tol_set: rtol sets every group', all(t%rtol == 1.0e-7_wp), maxval(abs(t%rtol - 1.0e-7_wp)))
      call check_true('build_tol_set: atol_scale multiplies every group''s atol', all(t%atol == 1.0e-2_wp * t1%atol), &
                      maxval(abs(t%atol - 1.0e-2_wp * t1%atol)))
      !---~---
      !   The soil heat has no integrator tolerance: its solve is one backward-Euler step and its
      ! accuracy is GRP_SE's. What the dial must still reach is the soil heat's budget stop, a
      ! closure threshold outside the groups, which atol_scale scales directly.
      !---~---
      ec = build_error_control(c)
      call check_true('build_error_control: atol_scale reaches the soil heat''s budget stop',                      &
                      ec%soil_energy_budget_atol == c%soil_energy_budget_atol * 1.0e-2_wp, ec%soil_energy_budget_atol)
   end subroutine test_tolerance_set

   !---~---
   !   The grouped error norm (state_wrms_grouped, docs/science/numerical_scheme.md eq. (3)): the
   ! RMS over the entries the march integrates, the canopy air and the soil heat, of each
   ! entry's difference in units of its group's tolerance, atol + rtol*|y_ref|. On the kernel
   ! record's state, whose magnitudes are physical, so a one-atol difference is well resolved.
   !---~---
   subroutine test_error_norm()
      integer(ik), parameter :: N = 2_ik, N_MANY = 6_ik
      type(march_state_t) :: y, a, y_neg, y_many, a_many
      type(frozen_t)      :: frozen
      type(tol_set_t)     :: t_abs, t_rel
      real(wp)            :: e1
      integer(ik)         :: j, nsl
      print '(a)', 'test_error_norm:'
      nsl = col%site%soil_profile%n_active
      call kernel_record(col%site, col%cfg%pft%hydraulics, N, y, frozen)
      y%leaf_surf_water = 0.1_wp
      y%wood_surf_water = 0.2_wp
      y%w_surface = 1.0_wp
      y%w_surface_energy = 1.0e5_wp
      t_abs%rtol = 0.0_wp
      t_rel%atol = 0.0_wp
      t_rel%rtol = 1.0e-3_wp

      call check('state_wrms_grouped: no difference, no error', state_wrms_grouped(y, y, y, N, nsl, t_abs), 0.0_wp, &
                 0.0_wp)
      !---~---
      !   What the march does not integrate is left out: the soil water and the plant water move at
      ! the water network's rate, the films and the pond by their maps.
      !---~---
      call state_init(y, N, nsl, a)
      a%theta = 0.5_wp * a%theta
      a%leaf_water_mass = 0.5_wp * a%leaf_water_mass
      a%wood_water_mass = 0.5_wp * a%wood_water_mass
      a%leaf_surf_water = 0.5_wp * a%leaf_surf_water
      a%wood_surf_water = 0.5_wp * a%wood_surf_water
      a%w_surface = 0.5_wp * a%w_surface
      a%w_surface_energy = 0.5_wp * a%w_surface_energy
      call check('state_wrms_grouped: the soil and plant water, the films and the pond are left out',              &
                 state_wrms_grouped(a, y, y, N, nsl, t_abs), 0.0_wp, 0.0_wp)
      !---~---
      !   One tolerance unit in every counted entry is an error of 1, with each group's own unit: the
      ! groups' atols span nine orders of magnitude, so a group read with another's is far off.
      !---~---
      call state_init(y, N, nsl, a)
      a%cas_enthalpy = y%cas_enthalpy + t_abs%atol(GRP_ENTH)
      a%cas_shv      = y%cas_shv      - t_abs%atol(GRP_SHV)
      a%cas_co2      = y%cas_co2      + t_abs%atol(GRP_CO2)
      a%soil_energy(1:nsl) = y%soil_energy(1:nsl) - t_abs%atol(GRP_SE)
      call check('state_wrms_grouped: one atol in every counted entry is 1', state_wrms_grouped(a, y, y, N, nsl, t_abs), &
                 1.0_wp, 1.0e-9_wp)
      call state_init(y, N, nsl, a)
      a%cas_enthalpy = y%cas_enthalpy * (1.0_wp + t_rel%rtol(GRP_ENTH))
      a%cas_shv      = y%cas_shv      * (1.0_wp + t_rel%rtol(GRP_SHV))
      a%cas_co2      = y%cas_co2      * (1.0_wp + t_rel%rtol(GRP_CO2))
      a%soil_energy(1:nsl) = y%soil_energy(1:nsl) * (1.0_wp + t_rel%rtol(GRP_SE))
      e1 = state_wrms_grouped(a, y, y, N, nsl, t_rel)
      call check('state_wrms_grouped: one rtol of |y_ref| in every counted entry is 1', e1, 1.0_wp, 1.0e-9_wp)
      !--- The relative part scales with y_ref's magnitude, whatever its sign.
      call state_sub(zero_like(y, N, nsl), y, N, nsl, y_neg)
      call check('state_wrms_grouped: the relative scale is |y_ref|', state_wrms_grouped(a, y, y_neg, N, nsl, t_rel), &
                 e1, 1.0e-12_wp)
      !---~---
      !   Six cohorts do not dilute the same canopy-air and soil-heat error: their entries are not
      ! counted.
      !---~---
      call kernel_record(col%site, col%cfg%pft%hydraulics, N_MANY, y_many, frozen)
      call state_init(y_many, N_MANY, nsl, a_many)
      a_many%cas_enthalpy = a%cas_enthalpy
      a_many%cas_shv      = a%cas_shv
      a_many%cas_co2      = a%cas_co2
      a_many%soil_energy  = a%soil_energy
      call check('state_wrms_grouped: six cohorts do not dilute the error',                                        &
                 state_wrms_grouped(a_many, y_many, y_many, N_MANY, nsl, t_rel), e1, 1.0e-12_wp)
      !--- Twice the difference, twice the error.
      call state_init(y, N, nsl, a)
      do j = 1_ik, nsl
         a%soil_energy(j) = y%soil_energy(j) + 2.0_wp * real(j, wp) * t_abs%atol(GRP_SE)
      end do
      e1 = state_wrms_grouped(a, y, y, N, nsl, t_abs)
      do j = 1_ik, nsl
         a%soil_energy(j) = y%soil_energy(j) + 4.0_wp * real(j, wp) * t_abs%atol(GRP_SE)
      end do
      call check('state_wrms_grouped: twice the difference, twice the error',                                      &
                 state_wrms_grouped(a, y, y, N, nsl, t_abs), 2.0_wp * e1, 1.0e-9_wp * e1)
   end subroutine test_error_norm

   !---~---
   !   The step factor (step_control_factor, docs/science/numerical_scheme.md eq. (4)) takes the
   ! embedded estimate's order and the growth limits from the error control: a step whose
   ! local error, of order p + 1 in the step, is 2^(p+1) tolerances is halved (times the
   ! safety factor), for ESDIRK2's first-order estimate and RK45's fourth-order one; an exact
   ! step grows by the control's fmax and a wild one shrinks by its fmin; and a larger error
   ! never gives a larger step. The I controller's own law is test_numerics'.
   !---~---
   subroutine test_step_factor()
      type(error_control_t) :: ec
      real(wp)    :: fac(40), err
      integer(ik) :: j
      print '(a)', 'test_step_factor:'
      ec%p_order = 1_ik
      call check('step_control_factor: a first-order estimate 4 tolerances off halves the step',                  &
                 step_control_factor(4.0_wp, ec), 0.5_wp * ec%safety, 1.0e-12_wp)
      ec%p_order = RK45_P_ORDER
      call check('step_control_factor: RK45''s fourth-order estimate 32 tolerances off halves the step',          &
                 step_control_factor(2.0_wp ** (RK45_P_ORDER + 1_ik), ec), 0.5_wp * ec%safety, 1.0e-12_wp)
      ec%fmin = 0.1_wp
      ec%fmax = 10.0_wp
      call check_true('step_control_factor: the control''s own limits bound it (an exact step, a wild one)',      &
                      step_control_factor(0.0_wp, ec) == ec%fmax                                          &
                      .and. step_control_factor(1.0e30_wp, ec) == ec%fmin,                                &
                      step_control_factor(0.0_wp, ec))
      do j = 1_ik, size(fac, kind=ik)
         err = 10.0_wp ** (-6.0_wp + 0.3_wp * real(j - 1_ik, wp))
         fac(j) = step_control_factor(err, ec)
      end do
      call check_true('step_control_factor: a larger error never gives a larger step', all(fac(2:) <= fac(:39)),   &
                      maxval(fac(2:) - fac(:39)))
   end subroutine test_step_factor

   !--- 7. THE WATER NETWORK'S STORES (ed_fast_integ_derivs: set_water_state, water_derivs)

   !---~---
   !   Each cohort's wood and leaves are read on its own PFT's pressure-volume curves (docs/
   ! science/soil_plant_water.md §1). The fixture's 20 cm tree, once of PFT 1 and once of PFT 2,
   ! in one column, each holding the water of PSI_WOOD in its wood and of PSI_LEAF in its leaves
   ! on its own curves: the network's state made from that water and evaluated there puts every
   ! store at its potential, holding its water. The fixture's PFTs differ in every curve
   ! parameter, so PFT 2's wood water read on PFT 1's curve sits at another potential. The step's
   ! drive (the soil's ice, the ground, the stomata's demand) is the frozen build's.
   !---~---
   subroutine test_water_curves_per_pft()
      integer(ik), parameter :: N = 2_ik
      integer(ik), parameter :: PFTS(N) = [1_ik, 2_ik]
      real(wp),    parameter :: PSI_WOOD = -0.6_wp, PSI_LEAF = -1.2_wp   ! [MPa]
      type(fast_column_t) :: two
      type(frozen_t)      :: fz
      type(march_state_t) :: ys
      real(wp)            :: leaf_w(N), wood_w(N), psi_off
      character(len=5)    :: tag
      integer(ik)         :: i, nsl
      print '(a)', 'test_water_curves_per_pft:'
      two = col
      call static_state_fixture(two%static, two%cfg%pft, PFTS, [20.0_wp, 20.0_wp], [0.0224_wp, 0.0224_wp],          &
                                two%site%soil_profile)
      two%static%geom%veg_height   = col%static%geom%veg_height
      two%static%geom%opencan_frac = col%static%geom%opencan_frac
      two%static%geom%snowfac      = col%static%geom%snowfac
      two%static%z_top             = col%static%z_top
      call ensure_output_capacity(two%out, N, .false.)
      call ensure_forcing_capacity(two%forcing, N)
      two%forcing%abs_sw_wood = 0.0_wp
      two%forcing%abs_lw_wood = 0.0_wp
      call reset_column(two)
      call set_noon_forcing(two)
      call frozen_record(two, fz, ys)
      nsl = two%site%soil_profile%n_active
      associate (s => two%static, hyd => two%cfg%pft%hydraulics, water => two%ws%water)
         do i = 1_ik, N
            wood_w(i) = tissue_psi_to_water(PSI_WOOD, hyd%pft(PFTS(i))%wood_curve, s%bwood_water(i))
            leaf_w(i) = tissue_psi_to_water(PSI_LEAF, hyd%pft(PFTS(i))%leaf_curve, s%bleaf(i))
         end do
         call set_water_state(s, two%site%soil_profile, hyd, two%dynamic%soil_water%theta(1:nsl),                    &
                              fz%water_drive%ice(1:nsl), 0.0_wp, 0.0_wp, leaf_w, wood_w, water)
         call water_derivs(s, two%site%soil_profile, two%site%soil, hyd, fz%water_drive, water%state%tau_soil,       &
                           water%state%tau_coh, water%ev(1), water%sys, .false.)
         do i = 1_ik, N
            write(tag, '(a,i0)') 'PFT ', s%pft(i)
            call check('water_derivs: '//tag//'''s wood at its water''s potential on its own curve',                &
                       water%ev(1)%psi_coh(WOOD, i), PSI_WOOD, 1.0e-9_wp)
            call check('water_derivs: '//tag//'''s wood holds its water', water%ev(1)%water_coh(WOOD, i),             &
                       wood_w(i) * s%nplant(i), 1.0e-12_wp * wood_w(i) * s%nplant(i))
            call check('water_derivs: '//tag//'''s leaves at their water''s potential on their own curve',           &
                       water%ev(1)%psi_coh(LEAF, i), PSI_LEAF, 1.0e-9_wp)
            call check('water_derivs: '//tag//'''s leaves hold their water', water%ev(1)%water_coh(LEAF, i),          &
                       leaf_w(i) * s%nplant(i), 1.0e-12_wp * leaf_w(i) * s%nplant(i))
         end do
         psi_off = tissue_water_to_psi(wood_w(2), hyd%pft(PFTS(1))%wood_curve, s%bwood_water(2)) - PSI_WOOD
         call check_true('water_derivs: PFT 2''s wood water read on PFT 1''s curve sits at another potential',        &
                         abs(psi_off) > 0.01_wp, psi_off)
      end associate
   end subroutine test_water_curves_per_pft

   elemental logical function ieee_ok(x)
      real(wp), intent(in) :: x
      ieee_ok = (x == x) .and. (abs(x) < huge(1.0_wp))
   end function ieee_ok

end program test_fast_integ_kernels
