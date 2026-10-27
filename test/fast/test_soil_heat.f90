! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_soil_heat -- the soil thermal column (ed_soil_heat;
! docs/science/soil_biophysics.md §The soil thermal column, §Frozen soil and §The bottom
! thermal boundary, eqs. (10)-(12)).
!   1. CONSERVATION (soil_energy_step_implicit): a column on an insulated base gains
!      exactly the heat let in at its top, with water moving between its layers or not,
!      and its top layer warms.
!   2. STEADY STATE (soil_energy_step_implicit): a sealed isothermal column does not
!      drift.
!   3. FREEZE (soil_energy_step_implicit on the freezing curve): a cooling wet layer cools
!      without a step, slowly just below t_3ple (the zero curtain), its liquid on the
!      curve and its temperature what its energy gives; its energy changes by exactly the
!      heat drawn out; Newton takes a few iterations; liquid remains at -5 C.
!   4. THAW: the liquid rises to all of the water below t_3ple, then the layer warms past
!      it; its energy changes by exactly the heat let in.
!   5. STEP INDEPENDENCE: a column frozen from the top for two days gives the same profile
!      with 30-minute and 1-minute steps.
!   6. ABOVE SATURATION (soil_layer_temp): a layer stored beyond theta_sat is liquid above
!      t_3ple, and its excess freezes at t_3ple.
!   7. MASS-CHANGE NEUTRALITY (soil_layer_temp): water that leaves or enters a layer with
!      its own enthalpy leaves the layer's temperature unchanged.
!   8. ANNUAL DAMPING: the amplitude of an annual surface wave decays with depth as
!      exp(-z/d), d = sqrt(2 alpha/omega), the semi-infinite analytic solution, under the
!      deep-temperature anchor and not under the insulated base, which reflects the wave.
!   9. THE TENDENCY (soil_energy_time_deriv): the explicit tendency the RK45 stages read
!      is the dt -> 0 limit of the backward-Euler step, with heat let in at the top, water
!      heat leaving every layer and a temperature gradient driving every face.
! The soil's thermal conductivity and heat capacity are test_therm_lib's.
!---~---
program test_soil_heat
   use ed_test_assert,      only : check, check_true, test_report
   use ed_kinds,            only : wp, ik
   use ed_constants,        only : rho_h2o, t_3ple, pi, n_soil_layer_max
   use ed_soil_heat,        only : soil_heat_env_t, soil_heat_out_t, soil_energy_step_implicit, soil_layer_energy,   &
                                   soil_layer_energy_slope, soil_layer_temp, soil_energy_time_deriv
   use ed_patch_state_vars,   only : soil_energy_state_t
   use ed_test_soil_grid,   only : geometric_layer_bottoms
   use ed_soil_params,      only : soil_profile_t, soil_params_t, ENERGY_BC_DIRICHLET
   use ed_soil_derive,      only : build_profile_hydraulics, build_profile_thermal
   use ed_soil_hydr_lib,    only : SOIL_RETENTION_VG, soil_liquid_capacity
   use ed_therm_lib,        only : soil_thermal_conductivity, soil_volumetric_heat_capacity,                         &
                                   temp_to_internal_energy, internal_energy_liquid
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks are computed on

   call test_soil_conserve()
   call test_soil_steady()
   call test_freeze_curve()
   call test_thaw_curve()
   call test_freeze_step_independence()
   call test_oversaturated_layer()
   call test_mass_change_neutrality()
   call test_soil_annual_damping()
   call test_soil_tendency()

   call test_report('test_soil_heat')

contains

   !--- A 2 m loam column of 10 layers at theta = 0.30, with no forcing.
   subroutine soil_setup(soil, forcing)
      type(soil_profile_t),  intent(out) :: soil
      type(soil_heat_env_t), intent(out) :: forcing
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,    &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soil)
      call build_profile_thermal(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, soil)
      forcing%soil_water(1:10) = 0.30_wp
      forcing%w_flux = 0.0_wp
      forcing%water_heat_out = 0.0_wp
      forcing%g_top = 0.0_wp
      forcing%geothermal = 0.0_wp
   end subroutine soil_setup

   subroutine init_col(col, soil, forcing, t_init)
      type(soil_energy_state_t), intent(out) :: col
      type(soil_profile_t),      intent(in)  :: soil
      type(soil_heat_env_t),     intent(in)  :: forcing
      real(wp),                  intent(in)  :: t_init
      integer(ik) :: k
      do k = 1_ik, 10_ik
         col%soil_energy(k) = temp_to_internal_energy(soil%dry_heat_capacity(k),                    &
                              forcing%soil_water(k) * rho_h2o, t_init, 1.0_wp)
      end do
   end subroutine init_col

   !--- The column's energy [J/m2], from its layers' state.
   real(wp) function column_energy(col, soil)
      type(soil_energy_state_t), intent(in) :: col
      type(soil_profile_t),      intent(in) :: soil
      column_energy = sum(col%soil_energy(1:10) * soil%dz(1:10))
   end function column_energy

   !---~---
   !   1. A 100 W/m2 warming on an insulated base: over 20 half-hour steps the column gains
   ! exactly 20 x 1800 x 100 J/m2, with no water moving and with a gentle downward flow
   ! between its layers, whose enthalpy only moves heat inside the column.
   !---~---
   subroutine test_soil_conserve()
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_energy_state_t) :: col
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: flux
      real(wp), parameter :: G_TOP = 100.0_wp, DT = 1800.0_wp
      integer(ik) :: step
      real(wp)    :: e0
      print '(a)', 'test_soil_conserve:'
      call soil_setup(soil, forcing)
      call init_col(col, soil, forcing, 285.0_wp)
      forcing%g_top = G_TOP
      e0 = column_energy(col, soil)
      do step = 1_ik, 20_ik
         call soil_energy_step_implicit(col, forcing, soil, opts, DT, flux)
      end do
      call check('the column gains the heat let in at its top', column_energy(col, soil) - e0, 20.0_wp * DT * G_TOP, &
                 1.0e-9_wp * abs(e0))
      call check_true('top layer warmed by G_top', col%soil_temp(1) > 285.0_wp, col%soil_temp(1))
      call init_col(col, soil, forcing, 285.0_wp)
      forcing%w_flux(1:9) = -1.0e-7_wp                            ! gentle downward flow
      do step = 1_ik, 20_ik
         call soil_energy_step_implicit(col, forcing, soil, opts, DT, flux)
      end do
      call check('with water moving between its layers, the same gain', column_energy(col, soil) - e0,          &
                 20.0_wp * DT * G_TOP, 1.0e-9_wp * abs(e0))
   end subroutine test_soil_conserve

   subroutine test_soil_steady()
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_energy_state_t) :: col
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: flux
      integer(ik) :: step, k
      real(wp)    :: drift
      print '(a)', 'test_soil_steady:'
      call soil_setup(soil, forcing)
      call init_col(col, soil, forcing, 285.0_wp)
      !--- g_top = geothermal = sink = w_flux = 0: an isothermal column must not drift.
      do step = 1_ik, 10_ik
         call soil_energy_step_implicit(col, forcing, soil, opts, 1800.0_wp, flux)
      end do
      drift = 0.0_wp
      do k = 1_ik, 10_ik
         drift = max(drift, abs(col%soil_temp(k) - 285.0_wp))
      end do
      call check_true('sealed isothermal column steady', drift < 1.0e-6_wp, drift)
   end subroutine test_soil_steady

   !---~---
   !   Single wet control volume, sealed bottom; the top Neumann G_top is the ONLY flux, so
   ! energy removed per step is exactly |G_top|*dt regardless of the BE predictor.
   !---~---
   subroutine setup_1layer(soil, forcing, theta, depth)
      type(soil_profile_t),  intent(out) :: soil
      type(soil_heat_env_t), intent(out) :: forcing
      real(wp),              intent(in)  :: theta, depth
      call geometric_layer_bottoms(depth, 1_ik, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(1_ik, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,      &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soil)
      call build_profile_thermal(1_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, soil)
      forcing%soil_water(1) = theta
      forcing%w_flux = 0.0_wp
      forcing%water_heat_out = 0.0_wp
      forcing%g_top = 0.0_wp
      forcing%geothermal = 0.0_wp
   end subroutine setup_1layer

   !---~---
   !   Cooling a wet layer through freezing on its freezing curve: it cools without a step, its
   ! liquid fraction on the curve all the way, slowly just below t_3ple (the zero curtain) as the
   ! latent heat comes out; its temperature is what its energy gives, which changes by exactly
   ! the heat drawn out; Newton takes a few iterations; and at -5 C the curve still holds some
   ! liquid.
   !---~---
   subroutine test_freeze_curve()
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_energy_state_t) :: col
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: flux
      real(wp), parameter :: THETA = 0.40_wp, DEPTH = 0.05_wp, DT = 1800.0_wp, G = -100.0_wp
      real(wp) :: e0, t_prev, worst_curve, worst_e, theta_liq, dtheta, t_curtain, e_back, de, fl_back
      integer(ik) :: step, iters_max, n_curtain
      logical :: monotone
      print '(a)', 'test_freeze_curve:'
      call setup_1layer(soil, forcing, THETA, DEPTH)
      col%soil_energy(1) = soil_layer_energy(t_3ple + 0.2_wp, THETA, soil%dry_heat_capacity(1), soil%layer_curve(1))
      e0 = col%soil_energy(1)
      forcing%g_top = G
      t_prev = t_3ple + 0.2_wp
      worst_curve = 0.0_wp
      worst_e = 0.0_wp
      iters_max = 0_ik
      n_curtain = 0_ik
      monotone = .true.
      do step = 1_ik, 200_ik
         call soil_energy_step_implicit(col, forcing, soil, opts, DT, flux)
         iters_max = max(iters_max, flux%iterations)
         if (col%soil_temp(1) > t_prev + 1.0e-12_wp) monotone = .false.
         t_prev = col%soil_temp(1)
         if (col%soil_temp(1) > t_3ple - 1.0_wp .and. col%soil_fliq(1) < 1.0_wp) n_curtain = n_curtain + 1_ik
         if (col%soil_fliq(1) < 1.0_wp) then
            call soil_liquid_capacity(soil%layer_curve(1), col%soil_temp(1), theta_liq, dtheta)
            worst_curve = max(worst_curve, abs(col%soil_fliq(1) * THETA - theta_liq))
         end if
         call soil_layer_energy_slope(col%soil_temp(1), THETA, soil%dry_heat_capacity(1), soil%layer_curve(1),   &
                                      e_back, de, fl_back)
         worst_e = max(worst_e, abs(e_back - col%soil_energy(1)) / de)
         if (col%soil_temp(1) < t_3ple - 5.0_wp) exit
      end do
      t_curtain = n_curtain * DT / 3600.0_wp
      call check('freeze: energy drawn out = heat flux x time', (col%soil_energy(1) - e0) * soil%dz(1),        &
                 G * DT * real(step, wp), 1.0e-9_wp * abs(e0 * soil%dz(1)))
      call check_true('freeze: cools monotonically, no step',             monotone, t_prev)
      call check_true('freeze: liquid on the freezing curve',             worst_curve < 1.0e-9_wp, worst_curve)
      call check_true('freeze: temperature is what the energy gives [K]', worst_e < 1.0e-9_wp, worst_e)
      call check_true('freeze: zero curtain, over 10 h within 1 K of t_3ple', t_curtain > 10.0_wp, t_curtain)
      call check_true('freeze: Newton converges in 4 iterations or fewer', iters_max <= 4_ik, real(iters_max, wp))
      call check_true('freeze: the curve still holds liquid at -5 C',     col%soil_fliq(1) > 0.05_wp .and.       &
                      col%soil_fliq(1) < 0.5_wp, col%soil_fliq(1))
   end subroutine test_freeze_curve

   !---~---
   !   Warming a frozen layer: its liquid fraction rises on the curve to 1 and the layer warms
   ! past t_3ple, all liquid; the energy closes.
   !---~---
   subroutine test_thaw_curve()
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_energy_state_t) :: col
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: flux
      real(wp), parameter :: THETA = 0.40_wp, DEPTH = 0.05_wp, DT = 1800.0_wp, G = 100.0_wp
      real(wp) :: e0, t_prev, worst_curve, theta_liq, dtheta
      integer(ik) :: step
      logical :: monotone
      print '(a)', 'test_thaw_curve:'
      call setup_1layer(soil, forcing, THETA, DEPTH)
      col%soil_energy(1) = soil_layer_energy(t_3ple - 5.0_wp, THETA, soil%dry_heat_capacity(1), soil%layer_curve(1))
      e0 = col%soil_energy(1)
      forcing%g_top = G
      t_prev = t_3ple - 5.0_wp
      worst_curve = 0.0_wp
      monotone = .true.
      do step = 1_ik, 200_ik
         call soil_energy_step_implicit(col, forcing, soil, opts, DT, flux)
         if (col%soil_temp(1) < t_prev - 1.0e-12_wp) monotone = .false.
         t_prev = col%soil_temp(1)
         if (col%soil_fliq(1) < 1.0_wp) then
            call soil_liquid_capacity(soil%layer_curve(1), col%soil_temp(1), theta_liq, dtheta)
            worst_curve = max(worst_curve, abs(col%soil_fliq(1) * THETA - theta_liq))
         end if
         if (col%soil_temp(1) > t_3ple + 1.0_wp) exit
      end do
      call check('thaw: energy added = heat flux x time', (col%soil_energy(1) - e0) * soil%dz(1),             &
                 G * DT * real(step, wp), 1.0e-9_wp * abs(e0 * soil%dz(1)))
      call check_true('thaw: warms monotonically, no step',          monotone, t_prev)
      call check_true('thaw: liquid on the freezing curve',          worst_curve < 1.0e-9_wp, worst_curve)
      call check_true('thaw: the layer warms past t_3ple, all liquid', col%soil_temp(1) > t_3ple + 1.0_wp .and.  &
                      col%soil_fliq(1) == 1.0_wp, t_prev)
   end subroutine test_thaw_curve

   !---~---
   !   A five-layer wet column frozen from the top for two days: the profile does not depend on
   ! the step, as Newton in the energies takes the latent heat into each step's balance (with
   ! only the sensible heat capacity in the step, a 30-minute step overshoots the freezing
   ! front).
   !---~---
   subroutine test_freeze_step_independence()
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_energy_state_t) :: col
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: flux
      real(wp) :: t_end(5, 2), dt(2), worst
      integer(ik) :: idt, step, nstep, k, iters_max
      print '(a)', 'test_freeze_step_independence:'
      dt = [1800.0_wp, 60.0_wp]
      iters_max = 0_ik
      do idt = 1, 2
         call geometric_layer_bottoms(1.0_wp, 5_ik, 3.0_wp, test_bottoms)
         call build_profile_hydraulics(5_ik, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,      &
              2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soil)
         call build_profile_thermal(5_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, soil)
         forcing%soil_water(1:5) = 0.38_wp
         forcing%w_flux = 0.0_wp
         forcing%water_heat_out = 0.0_wp
         forcing%g_top = -60.0_wp
         forcing%geothermal = 0.0_wp
         forcing%w_flux_top = 0.0_wp
         do k = 1_ik, 5_ik
            col%soil_energy(k) = soil_layer_energy(t_3ple + 1.0_wp, 0.38_wp, soil%dry_heat_capacity(k),   &
                                                   soil%layer_curve(k))
         end do
         nstep = nint(2.0_wp * 86400.0_wp / dt(idt), ik)
         do step = 1_ik, nstep
            call soil_energy_step_implicit(col, forcing, soil, opts, dt(idt), flux)
            iters_max = max(iters_max, flux%iterations)
         end do
         t_end(:, idt) = col%soil_temp(1:5)
      end do
      worst = maxval(abs(t_end(:, 1) - t_end(:, 2)))
      call check_true('freeze front: 30 min and 1 min steps agree within 0.1 K', worst < 0.1_wp, worst)
      call check_true('freeze front: the top layer froze', t_end(1, 2) < t_3ple - 1.0_wp, t_end(1, 2))
      call check_true('freeze front: Newton converges in 6 iterations or fewer', iters_max <= 6_ik,             &
                      real(iters_max, wp))
   end subroutine test_freeze_step_independence

   !---~---
   !   A layer stored beyond saturation (theta > theta_sat, through the specific storage): liquid
   ! at its temperature above t_3ple; below it on the curve; and between the curve's energy at
   ! t_3ple and the all-liquid one, at t_3ple with the excess partly frozen.
   !---~---
   subroutine test_oversaturated_layer()
      type(soil_profile_t)  :: soil
      type(soil_heat_env_t) :: forcing
      real(wp) :: theta, e, t, fl, e_lo, e_hi
      print '(a)', 'test_oversaturated_layer:'
      call setup_1layer(soil, forcing, 0.40_wp, 0.05_wp)
      theta = soil%theta_sat(1) + 5.0e-4_wp
      e = soil_layer_energy(300.0_wp, theta, soil%dry_heat_capacity(1), soil%layer_curve(1))
      call soil_layer_temp(e, theta, soil%dry_heat_capacity(1), soil%layer_curve(1), t, fl)
      call check('above saturation, warm: its temperature', t, 300.0_wp, 1.0e-9_wp)
      call check('above saturation, warm: all liquid', fl, 1.0_wp, 0.0_wp)
      e = soil_layer_energy(272.0_wp, theta, soil%dry_heat_capacity(1), soil%layer_curve(1))
      call soil_layer_temp(e, theta, soil%dry_heat_capacity(1), soil%layer_curve(1), t, fl)
      call check('above saturation, frozen: its temperature', t, 272.0_wp, 1.0e-9_wp)
      e_lo = soil_layer_energy(t_3ple, theta, soil%dry_heat_capacity(1), soil%layer_curve(1))
      e_hi = temp_to_internal_energy(soil%dry_heat_capacity(1), theta * rho_h2o, t_3ple, 1.0_wp)
      call soil_layer_temp(0.5_wp * (e_lo + e_hi), theta, soil%dry_heat_capacity(1), soil%layer_curve(1), t, fl)
      call check_true('the excess freezes at t_3ple', t == t_3ple .and. fl > soil%theta_sat(1) / theta .and.      &
                      fl < 1.0_wp, fl)
   end subroutine test_oversaturated_layer

   !---~---
   !   TEMPERATURE-NEUTRALITY of a mass change that carries its own enthalpy, on the soil's
   ! own inverter, soil_layer_temp.
   !
   !   The column inverts a layer's water content against its internal energy, so water
   ! that leaves or enters a layer without its matching enthalpy lands entirely in the
   ! diagnosed temperature. It is large: internal_energy_liquid carries the
   ! tsupercool_liq datum, so liquid water is ~1.0 MJ/kg in ABSOLUTE terms and only
   ! differences are physical. Moving the matching enthalpy with the mass, at the layer's
   ! OWN temperature, must leave the temperature exactly where it was.
   !
   !   Both branches are asserted: water moved without its enthalpy must MOVE the
   ! temperature (so the test fails if the datum is ever rebased to make u_liq ~ 0), and
   ! water moved with it must not.
   !---~---
   subroutine test_mass_change_neutrality()
      type(soil_profile_t)  :: soil
      type(soil_heat_env_t) :: forcing
      real(wp), parameter :: THETA0 = 0.40_wp, TEMP0 = 291.0_wp
      real(wp), parameter :: DTHETA = 0.01_wp     !--- [m3/m3] 1 kg/m2 out of or into a 0.1 m layer
      real(wp) :: e0, u_liq, temp_raw, temp_fix, fliq
      print '(a)', 'test_mass_change_neutrality:'
      call setup_1layer(soil, forcing, THETA0, 0.10_wp)
      e0 = soil_layer_energy(TEMP0, THETA0, soil%dry_heat_capacity(1), soil%layer_curve(1))

      !--- (a) mass removed, energy untouched: the layer's temperature jumps.
      call soil_layer_temp(e0, THETA0 - DTHETA, soil%dry_heat_capacity(1), soil%layer_curve(1), temp_raw, fliq)
      call check_true('water removed without its enthalpy visibly moves soil temperature',       &
                      abs(temp_raw - TEMP0) > 1.0_wp, temp_raw - TEMP0)

      !---~---
      !   (b) mass AND its enthalpy removed at the layer's own temperature: T is unchanged.
      ! internal_energy_liquid is what the model values every liquid-water flow with.
      !---~---
      u_liq = internal_energy_liquid(TEMP0)                      ! [J/kg] ~1.0e6, NOT ~0 -- that is the point
      call soil_layer_temp(e0 - DTHETA * rho_h2o * u_liq, THETA0 - DTHETA, soil%dry_heat_capacity(1),            &
                           soil%layer_curve(1), temp_fix, fliq)
      call check('water removed with its enthalpy leaves soil temperature unchanged', temp_fix, TEMP0, 1.0e-9_wp)
      call check_true('water removed with its enthalpy leaves the layer all-liquid',                &
                      abs(fliq - 1.0_wp) < 1.0e-12_wp, fliq)

      !--- Same for the OPPOSITE sign (water added to the layer).
      call soil_layer_temp(e0 + DTHETA * rho_h2o * u_liq, THETA0 + DTHETA, soil%dry_heat_capacity(1),            &
                           soil%layer_curve(1), temp_fix, fliq)
      call check('water added with its enthalpy leaves soil temperature unchanged', temp_fix, TEMP0, 1.0e-9_wp)
   end subroutine test_mass_change_neutrality

   !---~---
   !   8. ANNUAL-WAVE DAMPING vs the analytic semi-infinite solution.
   !
   !   A homogeneous column at constant theta has constant kappa and C, so the classical result
   ! applies exactly: a harmonic surface forcing of angular frequency omega decays with depth as
   ! exp(-z/d) with d = sqrt(2*alpha/omega), alpha = kappa/C. That makes the analytic profile an
   ! ORACLE for the bottom boundary condition, independent of the model's own machinery.
   !
   !   WHY THIS TEST EXISTS. A conservation check passes on BOTH bottom boundary conditions
   ! the insulated base conserves energy to round-off precisely BECAUSE it lets nothing out.
   ! What it does instead is REFLECT the annual wave back up the column, which no budget can
   ! see. Measured here: the insulated base holds 0.76 of the surface amplitude at the bottom
   ! node where the analytic answer is 0.42.
   !---~---
   subroutine test_soil_annual_damping()
      integer(ik), parameter :: NL = 10_ik, NYR = 4_ik
      real(wp),    parameter :: PERIOD = 365.0_wp * 86400.0_wp
      real(wp),    parameter :: OMEGA  = 2.0_wp * pi / PERIOD
      real(wp),    parameter :: DTS    = 3600.0_wp
      real(wp),    parameter :: Q0     = 8.0_wp       ! [W/m2] surface-wave amplitude
      real(wp),    parameter :: THETA  = 0.30_wp, TMEAN = 288.0_wp
      real(wp) :: amp_neu(NL), amp_dir(NL), ana(NL), z(NL), dd
      real(wp) :: err_neu, err_dir
      integer(ik) :: k
      print '(a)', 'test_soil_annual_damping:'

      call damping_profile(NL, NYR, PERIOD, OMEGA, DTS, Q0, THETA, TMEAN, -1.0_wp, amp_neu, z, dd)
      call damping_profile(NL, NYR, PERIOD, OMEGA, DTS, Q0, THETA, TMEAN, 3.12_wp, amp_dir, z, dd)
      do k = 1_ik, NL
         ana(k) = exp(-(abs(z(k)) - abs(z(1))) / dd)
      end do
      !---~---
      !   Normalise each profile to its own layer-1 amplitude: the analytic decay is a RATIO,
      ! so this is what makes the two boundary conditions comparable at all.
      !---~---
      amp_neu = amp_neu / amp_neu(1)
      amp_dir = amp_dir / amp_dir(1)
      err_neu = sqrt(sum((amp_neu - ana)**2) / real(NL, wp))
      err_dir = sqrt(sum((amp_dir - ana)**2) / real(NL, wp))

      !---~---
      !   The damping depth this test is built on, from the soil's properties by hand, as a guard
      ! on its premise: Johansen's conductivity at theta 0.30 of a 0.43 porosity, a Kersten
      ! number 1 + log10(0.30/0.43) = 0.8437 on the saturated 3.0^0.57 x 0.57^0.43 = 1.4689 W/m/K
      ! and the dry 0.15, gives kappa = 1.2627 W/m/K; with C = 2.0e6 + 0.30 x 4.186e6 = 3.2558e6
      ! J/m3/K, d = sqrt(2 kappa / (C omega)) = 1.9731 m.
      !---~---
      call check('analytic damping depth [m]', dd, 1.9731_wp, 1.0e-3_wp)

      !---~---
      !   The insulated base's reflection, stated as a number: it holds far too much amplitude
      ! at the bottom node. This pins the oracle's power to tell the two BCs apart.
      !---~---
      call check_true('insulated base over-amplifies (>1.6x analytic)',                            &
                      amp_neu(NL) > 1.6_wp * ana(NL), amp_neu(NL) / ana(NL))
      !--- The Dirichlet anchor tracks the analytic profile to a few percent.
      call check_true('dirichlet base within 10% of analytic',                                     &
                      abs(amp_dir(NL) - ana(NL)) < 0.10_wp * ana(NL), amp_dir(NL) / ana(NL))
      !---~---
      !   And it is better over the WHOLE profile, not just at the base -- by a wide enough
      ! margin that the assertion does not depend on the exact anchor depth.
      !---~---
      call check_true('dirichlet profile >5x closer to analytic (RMS)',                            &
                      err_dir < 0.2_wp * err_neu, err_neu / err_dir)
   end subroutine test_soil_annual_damping

   !---~---
   !   Drive one homogeneous column with a harmonic G_top to periodic steady state and return
   ! the peak-to-peak amplitude of each layer over the final year. `ddeep <= 0` selects the
   ! prescribed-flux (insulated) BC, otherwise the Dirichlet anchor at that depth.
   !---~---
   subroutine damping_profile(nl, nyr, period, omega, dts, q0, theta, tmean, ddeep, amp, z, dd)
      integer(ik), intent(in)  :: nl, nyr
      real(wp),    intent(in)  :: period, omega, dts, q0, theta, tmean, ddeep
      real(wp),    intent(out) :: amp(:), z(:), dd
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: eflux
      type(soil_energy_state_t) :: col
      real(wp)    :: tmin(n_soil_layer_max), tmax(n_soil_layer_max), tnow, kap, cef
      integer(ik) :: k, istep, nstep_yr, iyr

      call geometric_layer_bottoms(2.0_wp, nl, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(nl, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,     &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soil)
      call build_profile_thermal(nl, 3.0_wp, 0.15_wp, 2.0e6_wp, soil)
      forcing%soil_water(1:nl) = theta
      forcing%w_flux = 0.0_wp
      forcing%water_heat_out = 0.0_wp
      forcing%g_top = 0.0_wp
      forcing%geothermal = 0.0_wp
      forcing%w_flux_top = 0.0_wp
      if (ddeep > 0.0_wp) then
         opts%thermal_bottom_bc = ENERGY_BC_DIRICHLET
         opts%deep_temp = tmean
         opts%deep_depth = ddeep
      end if
      do k = 1_ik, nl
         col%soil_energy(k) = temp_to_internal_energy(soil%dry_heat_capacity(k),            &
                                 theta * rho_h2o, tmean, 1.0_wp)
      end do
      kap = soil_thermal_conductivity(theta, 1.0_wp, soil%theta_sat(1),                                   &
                              soil%solid_conductivity(1), soil%dry_conductivity(1))
      cef = soil_volumetric_heat_capacity(theta, 1.0_wp, soil%dry_heat_capacity(1))
      dd  = sqrt(2.0_wp * (kap / cef) / omega)

      nstep_yr = int(period / dts, ik)
      tmin = 1.0e30_wp
      tmax = -1.0e30_wp
      do iyr = 1_ik, nyr
         do istep = 1_ik, nstep_yr
            tnow = real(istep, wp) * dts
            forcing%g_top = q0 * sin(omega * tnow)
            call soil_energy_step_implicit(col, forcing, soil, opts, dts, eflux)
            if (iyr == nyr) then
               do k = 1_ik, nl
                  tmin(k) = min(tmin(k), col%soil_temp(k))
                  tmax(k) = max(tmax(k), col%soil_temp(k))
               end do
            end if
         end do
      end do
      do k = 1_ik, nl
         amp(k) = tmax(k) - tmin(k)
         z(k) = soil%z_node(k)
      end do
   end subroutine damping_profile

   !---~---
   !   9. The explicit tendency (soil_energy_time_deriv) is the dt -> 0 limit of the backward-Euler
   ! step (soil_energy_step_implicit): one 0.01 s step from a column cooling by 0.6 K a layer,
   ! 120 W/m2 let in at the top and 3 W/m3 of water heat leaving every layer, moves each layer's
   ! energy by its tendency to O(dt).
   !---~---
   subroutine test_soil_tendency()
      real(wp), parameter :: DT = 1.0e-2_wp
      type(soil_profile_t)      :: soil
      type(soil_heat_env_t)     :: forcing
      type(soil_energy_state_t) :: col, col0
      type(soil_params_t)       :: opts
      type(soil_heat_out_t)     :: flux
      real(wp)    :: dedt(n_soil_layer_max), worst
      integer(ik) :: k
      print '(a)', 'test_soil_tendency:'
      call soil_setup(soil, forcing)
      forcing%g_top = 120.0_wp
      forcing%water_heat_out(1:10) = 3.0_wp
      do k = 1_ik, 10_ik
         col%soil_energy(k) = temp_to_internal_energy(soil%dry_heat_capacity(k), forcing%soil_water(k) * rho_h2o,     &
                                                      291.0_wp - 0.6_wp * real(k - 1_ik, wp), 1.0_wp)
      end do
      col0 = col
      call soil_energy_time_deriv(col0, forcing, soil, opts, dedt)
      call soil_energy_step_implicit(col, forcing, soil, opts, DT, flux)
      worst = 0.0_wp
      do k = 1_ik, 10_ik
         worst = max(worst, abs((col%soil_energy(k) - col0%soil_energy(k)) / DT - dedt(k)) / max(abs(dedt(k)), 1.0_wp))
      end do
      call check_true('the tendency is the backward-Euler step''s dt -> 0 limit (O(dt))', worst < 1.0e-2_wp, worst)
   end subroutine test_soil_tendency

end program test_soil_heat
