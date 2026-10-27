! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_water_network -- the soil-plant water network: its state (set_water_state), flows
! (water_derivs) and step (water_step) (docs/science/soil_plant_water.md).
!
!   1. ARROWHEAD SOLVE : the Schur solve equals dense elimination on a random arrowhead system.
!   2. JACOBIAN        : the assembled Jacobian equals the residual's derivative (wet, dry,
!                        ponded, leafless).
!   3. CONSERVATION    : two days of rain, transpiration, evaporation and drainage close to
!                        round-off, step by step.
!   4. PLANT ONLY      : a soil too deep for the stand to draw down: under a steady demand the
!                        plant's steady state carries the transpiration through every path.
!   5. SOIL ONLY       : no cohorts, a sealed bottom: the column drains nothing and relaxes to
!                        hydrostatic equilibrium.
!   6. STIFF EXCHANGE  : a wet, densely rooted stand whose wood-soil exchange takes a fraction of
!                        a second costs a handful of sub-steps, not thousands.
!   7. INFEASIBLE DEMAND: a leaf asked for many times the water it and its xylem can supply stays
!                        above empty and transpires less than asked.
!   8. LEAFLESS        : a cohort with no leaves: its leaf sits at the wood less the crown head.
!   9. DRY-LAYER STORM : a thin top layer near residual, a stand refilling its wood, then a
!                        storm: converges, no store leaves its range.
!  10. SATURATED       : rain on a saturated column ponds, spills, and conserves.
!  11. ACCURACY        : the default tolerances track a tight-tolerance run.
!  12. VANISHING STORES: a near-leafless leaf and a trace pond stay in range, and water closes.
!  13. CAMPBELL        : on a Campbell clay the soil potentials are the Campbell curve's, finite
!                        and below zero; the stand takes up water, and every step closes.
! The soil column alone, on a 2 m loam (docs/science/soil_biophysics.md §The soil water column
! and §Boundary conditions):
!  14. DRAINAGE        : a wet column drains freely, and through a leaky base about f of that;
!                        its water closes against the drainage.
!  15. DOWNPOUR        : rain beyond a wet clay's infiltration ponds, and runs off past the pond.
!  16. REWETTING       : a dried top layer re-wets within 12 h of rain, with little runoff.
!  17. SUB-FREEZING    : rain below the triple point: the pond's enthalpy closes, and what
!                        infiltrates carries its own enthalpy.
!  18. AQUIFER         : the head-driven bottom is two-way and settles to hydrostatic equilibrium.
!  19. FROZEN TOP      : the ice impedes the rain's infiltration; the freezing curve's suction
!                        fills a frozen layer to saturation, not past it; the water closes.
!  20. FROZEN SURFACE  : the dry surface layer reads the liquid and dsl_theta_init of the ice
!                        (water_derivs), so it forms below dsl_theta_init of the pores the ice
!                        leaves (CLM5's effective porosity); above that the layer evaporates as
!                        one with no dry surface layer, at its liquid's potential.
!---~---
program test_water_network
   use ed_test_assert, only : check_true, check_close, check_int, test_report
   use, intrinsic :: iso_fortran_env, only : int64
   use ed_kinds,         only : wp, ik
   use ed_constants,     only : grav_head, rho_h2o, pi, tiny_num, n_soil_layer_max, grav, r_wv, p_std
   use ed_linear_algebra,  only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   use ed_fast_integ_types,  only : fast_static_state_t, water_drive_t, water_work_t, tol_set_t,                    &
                                ensure_static_capacity, POND, WOOD, LEAF, GRP_THETA, N_WORK, WORK_WATER_NSUB,   &
                                WORK_WATER_REJECT, WORK_WATER_NFAIL
   use ed_fast_state,    only : derive_cohort_terms
   use ed_fast_integ_derivs, only : set_water_state, water_derivs
   use ed_fast_integ_backward_euler, only : water_step, water_jacobian_error
   use ed_plant_hydr_lib,   only : water_curve_t, tissue_psi_to_rwc, tissue_psi_tlp
   use ed_pft_derive,       only : build_kirchhoff_table
   use ed_soil_hydr_lib,     only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, soil_theta_to_head, soil_ice_state,  &
                                    soil_liquid_capacity
   use ed_therm_lib,         only : internal_energy_liquid, liquid_energy_to_temp, sat_specific_humidity
   use ed_pft_params,    only : hydraulics_table_t
   use ed_soil_params, only : soil_profile_t, soil_params_t, SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use ed_soil_derive,        only : build_profile_hydraulics
   use ed_test_soil_grid, only : geometric_layer_bottoms
   use ed_soil_water,    only : ground_evap_env_t, ground_evaporation_setup
   implicit none

   !---~---
   !   A stand on the fixture soil: its cohorts (the static state the gather would build), what
   ! drives the network, the network itself, its tolerances and the work it last reported.
   !---~---
   type :: stand_t
      type(fast_static_state_t) :: static
      type(water_drive_t)       :: drive
      type(water_work_t)        :: water
      type(tol_set_t)           :: tols
      integer(ik)               :: work(N_WORK) = 0_ik
   end type stand_t

   type(hydraulics_table_t)   :: hyd
   type(soil_profile_t)       :: soil
   type(soil_params_t)        :: sopts
   type(ground_evap_env_t)    :: forc
   !--- The water network's tolerances: [fast].rtol and water_atol at their defaults.
   real(wp), parameter :: WATER_RTOL = 1.0e-3_wp, WATER_ATOL = 1.0e-4_wp

   call make_pft(hyd)
   call make_soil(soil, 10_ik, 1.44_wp, 3.0_wp)
   forc%t_ground = 299.0_wp
   forc%q_air = 0.016_wp
   forc%rho_air = 1.15_wp
   forc%r_aero = 80.0_wp

   call test_arrowhead_solve()
   call test_jacobian()
   call test_conservation()
   call test_plant_only()
   call test_soil_only()
   call test_stiff_exchange()
   call test_infeasible_demand()
   call test_leafless()
   call test_dry_layer_storm()
   call test_saturated()
   call test_accuracy()
   call test_vanishing_stores()
   call test_campbell()
   call test_column_scenarios()
   call test_report('test_water_network')

contains

   !--- Fixtures.

   !---~---
   !   One PFT with the shipped [hydraulics] values and a tropical wood conductivity; the
   ! stomata shut at twice the leaf's turgor-loss potential.
   !---~---
   subroutine make_pft(h)
      type(hydraulics_table_t),   intent(out) :: h
      h%n_pft = 1_ik
      allocate(h%pft(1))
      h%pft(1)%leaf_curve = water_curve_t(pi0=-1.5_wp, elastic_mod=12.0_wp, apoplast_frac=0.30_wp, water_sat=2.0_wp)
      h%pft(1)%wood_curve = water_curve_t(pi0=-1.0_wp, elastic_mod=8.0_wp,  apoplast_frac=0.20_wp, water_sat=1.0_wp)
      h%pft(1)%wood_psi50  = -2.0_wp
      h%pft(1)%wood_kexp = 2.0_wp
      h%pft(1)%wood_kmax   = 3.3_wp
      h%pft(1)%vessel_curl = 1.5_wp
      h%pft(1)%root_beta = 0.1_wp
      h%pft(1)%root_depth_b1 = 1.1140580_wp
      h%pft(1)%root_depth_b2 = 0.4223014_wp
      h%pft(1)%specific_root_length = 2.0e4_wp
      h%pft(1)%fine_root_radius = 2.5e-4_wp
      h%pft(1)%rwc_close = tissue_psi_to_rwc(2.0_wp * tissue_psi_tlp(h%pft(1)%leaf_curve%pi0,                    &
                                             h%pft(1)%leaf_curve%elastic_mod), h%pft(1)%leaf_curve%pi0,           &
                                             h%pft(1)%leaf_curve%elastic_mod)
      call build_kirchhoff_table(h%pft(1)%vuln_table, h%pft(1)%wood_kexp)
   end subroutine make_pft

   !--- BCI's clay (example 04), on a grid whose top layer is 2.6 cm.
   subroutine make_soil(s, nsl, depth, growth)
      type(soil_profile_t), intent(out) :: s
      integer(ik),          intent(in)  :: nsl
      real(wp),             intent(in)  :: depth, growth
      real(wp) :: bottoms(n_soil_layer_max)
      call geometric_layer_bottoms(depth, nsl, growth, bottoms)
      call build_profile_hydraulics(nsl, SOIL_RETENTION_VG, bottoms, 0.4314_wp, 0.0981_wp, 4.54e-6_wp,        &
                                  0.7607_wp, 1.3554_wp, -3.37_wp, s)
   end subroutine make_soil

   !---~---
   !   A stand of n cohorts on the fixture soil, a canopy tree first and smaller ones after,
   ! with the roots and water paths the gather derives for them, each cohort's demand its
   ! share of e_total by leaf, and the network started from the given water.
   !---~---
   subroutine make_stand(st, n, theta, w_pond, leaf_rwc, wood_rwc, e_total, leafless)
      type(stand_t),     intent(out) :: st
      integer(ik),       intent(in)  :: n
      real(wp),          intent(in)  :: theta(:), w_pond, leaf_rwc, wood_rwc, e_total
      logical, optional, intent(in)  :: leafless
      real(wp)    :: leaf_w(n), wood_w(n), lai_share(n)
      integer(ik) :: i
      call ensure_static_capacity(st%static, n)
      do i = 1_ik, n
         st%static%pft(i)      = 1_ik
         st%static%height(i)   = 28.0_wp / real(i, wp)
         st%static%nplant(i)   = 0.02_wp * real(i, wp) ** 1.5_wp
         st%static%bleaf(i)    = 4.0_wp / real(i, wp) ** 2
         st%static%bsap(i)     = 30.0_wp / real(i, wp) ** 2.5_wp
         st%static%broot(i)    = 2.0_wp / real(i, wp) ** 2
         st%static%sap_area(i) = 0.05_wp / real(i, wp) ** 2
      end do
      if (present(leafless)) then
         if (leafless) st%static%bleaf(n) = 0.0_wp
      end if
      call derive_cohort_terms(st%static, hyd, soil)
      associate (s => st%static, lc => hyd%pft(1)%leaf_curve, wc => hyd%pft(1)%wood_curve)
         lai_share = s%nplant * s%bleaf / max(sum(s%nplant * s%bleaf), tiny_num)
         leaf_w = s%bleaf * lc%water_sat * ((1.0_wp - lc%apoplast_frac) * leaf_rwc + lc%apoplast_frac)
         wood_w = s%bwood_water * wc%water_sat * ((1.0_wp - wc%apoplast_frac) * wood_rwc + wc%apoplast_frac)
      end associate
      allocate(st%drive%e_demand(n))
      st%drive%e_demand = e_total * lai_share
      st%drive%rain = 0.0_wp
      st%drive%t_soil = 299.0_wp
      st%drive%t_rain = 297.0_wp
      st%drive%ground = forc
      st%tols%rtol(GRP_THETA) = WATER_RTOL
      st%tols%atol(GRP_THETA) = WATER_ATOL
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, w_pond, 0.0_wp, leaf_w, wood_w, st%water)
   end subroutine make_stand

   !--- One fast step of the stand's network on the fixture soil, with the column's options.
   subroutine advance(st, opts, dt)
      type(stand_t),       intent(inout) :: st
      type(soil_params_t), intent(in)    :: opts
      real(wp),            intent(in)    :: dt
      call water_step(dt, st%static, soil, opts, hyd, st%drive, st%tols, st%water, st%work)
   end subroutine advance

   !--- All the water the network holds [kg/m2].
   real(wp) function total_water(st)
      type(stand_t), intent(in) :: st
      total_water = sum(st%water%state%water_soil) + sum(st%water%state%water_coh)
   end function total_water

   !--- The water that crossed the network's boundary over its last step [kg/m2].
   real(wp) function boundary_water(st, dt)
      type(stand_t), intent(in) :: st
      real(wp),      intent(in) :: dt
      associate (f => st%water%flows)
         boundary_water = st%drive%rain * dt - f%drain - f%evap - f%spill - sum(f%transp)
      end associate
   end function boundary_water

   subroutine test_arrowhead_solve()
      type(arrowhead_system_t) :: s
      integer(ik), parameter :: NC = 7, NG = 5, MG = 2, NT = NC + NG*MG
      real(wp) :: full(NT, NT), rhs(NT), x_ref(NT), r_c(NC), r_g(MG, NG), x_c(NC), x_g(MG, NG), v
      real(wp) :: add_c(NC), add_g(MG, NG), a_kept(NC, NC)
      integer(ik) :: i, j, g, l, seed
      logical :: ok
      call arrowhead_allocate(s, NC, NG, MG)
      seed = 12345_ik
      full = 0.0_wp
      do j = 1_ik, NC
         do i = 1_ik, NC
            v = pseudo(seed)
            s%a(i, j) = v
            full(i, j) = v
         end do
         s%a(j, j) = s%a(j, j) + 8.0_wp
         full(j, j) = full(j, j) + 8.0_wp
      end do
      do g = 1_ik, NG
         do i = 1_ik, NC            ! each group meets the core through its first unknown
            v = pseudo(seed)
            s%b(i, g) = v
            full(i, NC + (g-1)*MG + 1) = v
            v = pseudo(seed)
            s%c(i, g) = v
            full(NC + (g-1)*MG + 1, i) = v
         end do
         do l = 1_ik, MG
            do j = 1_ik, MG
               v = pseudo(seed) + merge(5.0_wp, 0.0_wp, j == l)
               s%d(l, j, g) = v
               full(NC + (g-1)*MG + l, NC + (g-1)*MG + j) = v
            end do
         end do
      end do
      !--- Additions to the blocks' diagonals, as a Newton matrix's storage terms are added.
      do i = 1_ik, NC
         add_c(i) = 1.5_wp + 0.1_wp * real(i, wp)
         full(i, i) = full(i, i) + add_c(i)
      end do
      do g = 1_ik, NG
         do l = 1_ik, MG
            add_g(l, g) = 0.7_wp
            full(NC + (g-1)*MG + l, NC + (g-1)*MG + l) = full(NC + (g-1)*MG + l, NC + (g-1)*MG + l) + 0.7_wp
         end do
      end do
      do i = 1_ik, NT
         rhs(i) = pseudo(seed)
      end do
      r_c = rhs(1:NC)
      r_g = reshape(rhs(NC+1:NT), [MG, NG])
      call dense_solve(full, rhs, x_ref, NT)
      a_kept = s%a
      call arrowhead_factor(s, add_c, add_g, ok)
      call check_true('arrowhead solve: factorisation succeeds', ok)
      call check_true('arrowhead solve: factoring leaves the blocks as they were', all(s%a == a_kept))
      call arrowhead_solve(s, r_c, r_g, x_c, x_g)
      call check_true('arrowhead solve: equals dense elimination',                                             &
                      maxval(abs([x_c, reshape(x_g, [NG*MG])] - x_ref)) < 1.0e-12_wp * maxval(abs(x_ref)))
   end subroutine test_arrowhead_solve

   subroutine test_jacobian()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), err
      theta = 0.36_wp
      call make_stand(st, 3_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 4.0e-5_wp)
      err = water_jacobian_error(st%static, soil, sopts, hyd, st%drive, st%water, 300.0_wp)
      call check_true('jacobian: wet soil, transpiring stand', err < 1.0e-5_wp, err)
      theta = 0.12_wp
      theta(1) = 0.1 + 0.003_wp
      call make_stand(st, 3_ik, theta, 0.0_wp, 0.4_wp, 0.35_wp, 2.0e-5_wp)
      err = water_jacobian_error(st%static, soil, sopts, hyd, st%drive, st%water, 300.0_wp)
      call check_true('jacobian: dry soil, stressed stand', err < 1.0e-5_wp, err)
      theta = 0.4314_wp - 1.0e-3_wp
      theta(1) = 0.4314_wp + 1.0e-5_wp
      call make_stand(st, 2_ik, theta, 12.0_wp, 0.98_wp, 0.97_wp, 0.0_wp)
      err = water_jacobian_error(st%static, soil, sopts, hyd, st%drive, st%water, 300.0_wp)
      call check_true('jacobian: ponded, supersaturated top layer', err < 1.0e-5_wp, err)
      theta = 0.3_wp
      call make_stand(st, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 1.0e-5_wp, leafless=.true.)
      err = water_jacobian_error(st%static, soil, sopts, hyd, st%drive, st%water, 300.0_wp)
      call check_true('jacobian: a leafless cohort', err < 1.0e-5_wp, err)
   end subroutine test_jacobian

   subroutine test_conservation()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, w1, dt, hour, resid, worst, worst_coh
      real(wp), allocatable :: wc0(:), wc1(:)
      integer(ik) :: step, nfail
      theta = 0.33_wp
      call make_stand(st, 4_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp)
      dt = 900.0_wp
      worst = 0.0_wp
      worst_coh = 0.0_wp
      nfail = 0_ik
      do step = 1_ik, 192
         hour = mod(real(step, wp) * dt / 3600.0_wp, 24.0_wp)
         st%drive%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                     &
                             * st%static%nplant * st%static%bleaf / sum(st%static%nplant * st%static%bleaf)
         st%drive%rain = 0.0_wp
         if (step >= 60 .and. step < 68) st%drive%rain = 2.0e-2_wp      ! 72 mm/h for two hours
         if (step >= 150 .and. step < 160) st%drive%rain = 1.0e-3_wp    ! a long drizzle
         w0 = total_water(st)
         wc0 = sum(st%water%state%water_coh, dim=1)
         call advance(st, sopts, dt)
         w1 = total_water(st)
         wc1 = sum(st%water%state%water_coh, dim=1)
         !--- Each cohort's leaf + wood water changes by its own uptake less its transpiration.
         worst_coh = max(worst_coh, maxval(abs((wc1 - wc0) - (sum(st%water%flows%uptake, dim=1)               &
                                                              - st%water%flows%transp))) / max(maxval(wc0), 1.0e-6_wp))
         resid = abs((w1 - w0) - boundary_water(st, dt))
         worst = max(worst, resid / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
      end do
      call check_true('conservation: every step closes to round-off', worst < 1.0e-13_wp, worst)
      call check_true('conservation: every cohort''s tissue water closes to round-off', worst_coh < 1.0e-12_wp, &
                      worst_coh)
      call check_int('conservation: no sub-step forced through', nfail, 0_ik)
   end subroutine test_conservation

   !---~---
   !   One cohort on a single soil layer 100 m deep, sealed and under still air, which a day's
   ! transpiration barely draws down. Under a steady demand the plant reaches its steady state,
   ! where every path carries the transpiration: root uptake = sapflow = transpiration, and
   ! the stores stop changing -- to within what the soil's slow drying still moves. The
   ! fixture's soil is swapped for the test and put back after it.
   !---~---
   subroutine test_plant_only()
      type(stand_t)          :: st
      type(soil_profile_t)   :: shallow
      type(soil_params_t)    :: sealed
      real(wp) :: theta(1), w_plant0, w_plant1, e, u, sap
      integer(ik) :: step
      shallow = soil
      call make_soil(soil, 1_ik, 100.0_wp, 3.0_wp)
      sealed = sopts
      sealed%water_bottom_bc = SOIL_BC_DRAIN
      sealed%drain_factor = 0.0_wp
      theta = 0.28_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.92_wp, 3.0e-5_wp)
      st%drive%ground%q_air = 0.05_wp                 ! air wetter than any surface: no evaporation
      st%tols%rtol(GRP_THETA) = 1.0e-8_wp
      st%tols%atol(GRP_THETA) = 1.0e-9_wp
      do step = 1_ik, 96_ik
         w_plant0 = sum(st%water%state%water_coh)
         call advance(st, sealed, 900.0_wp)
         w_plant1 = sum(st%water%state%water_coh)
      end do
      e = st%water%flows%transp(1)
      u = sum(st%water%flows%uptake(:, 1))
      sap = st%water%flows%sapflow(1)
      print '(a,3es11.3)', '  info : plant only: |u-e|/e, |sap-e|/e, |dW|/W: ', abs(u - e) / e, abs(sap - e) / e,   &
            abs(w_plant1 - w_plant0) / w_plant1
      soil = shallow
      call check_close('plant only: at steady state root uptake = transpiration', u, e, 1.0e-5_wp * e)
      call check_close('plant only: at steady state sapflow = transpiration', sap, e, 1.0e-5_wp * e)
      call check_close('plant only: the plant water stops changing', w_plant1 - w_plant0, 0.0_wp, 1.0e-7_wp * w_plant1)
      call check_true('plant only: the leaf sits below the wood by more than the crown head',                   &
                      st%water%state%psi_coh(LEAF, 1) < st%water%state%psi_coh(WOOD, 1) - st%static%crown_head(1), &
                      st%water%state%psi_coh(LEAF, 1))
   end subroutine test_plant_only

   !---~---
   !   No cohorts, a sealed bottom, no rain and no evaporation: the column relaxes to hydrostatic
   ! equilibrium, the total head (matric head plus elevation) the same in every layer.
   !---~---
   subroutine test_soil_only()
      type(stand_t)       :: st
      type(soil_params_t) :: sealed
      real(wp) :: theta(soil%n_active), head(soil%n_active), w0, drain_max
      integer(ik) :: step, nsl
      nsl = soil%n_active
      theta = 0.40_wp
      theta(1:3) = 0.425_wp                            ! wet: BCI's clay settles in days, not years
      sealed = sopts
      sealed%water_bottom_bc = SOIL_BC_DRAIN
      sealed%drain_factor = 0.0_wp
      call make_stand(st, 0_ik, theta, 0.0_wp, 0.9_wp, 0.9_wp, 0.0_wp)
      st%drive%ground%q_air = 0.05_wp                 ! air wetter than any surface: no evaporation
      w0 = total_water(st)
      drain_max = 0.0_wp
      do step = 1_ik, 96_ik * 30_ik
         call advance(st, sealed, 900.0_wp)
         drain_max = max(drain_max, abs(st%water%flows%drain))
      end do
      head = st%water%state%psi_soil(2:nsl+1) / grav_head + soil%z_node(1:nsl)
      call check_true('soil only: hydrostatic equilibrium (total head uniform to 1 cm)',                        &
                      maxval(head) - minval(head) < 1.0e-2_wp, maxval(head) - minval(head))
      call check_close('soil only: sealed column keeps its water', total_water(st), w0, 1.0e-11_wp * w0)
      call check_true('soil only: a sealed bottom drains nothing', drain_max < 1.0e-30_wp, drain_max)
   end subroutine test_soil_only

   subroutine test_stiff_exchange()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active)
      integer(ik) :: step, nsub
      theta = 0.42_wp
      call make_stand(st, 4_ik, theta, 0.0_wp, 0.97_wp, 0.95_wp, 0.0_wp)
      st%static%g_rhizo = 50.0_wp * st%static%g_rhizo      ! a fifty-fold root density: a sub-second exchange
      nsub = 0_ik
      do step = 1_ik, 8_ik
         st%drive%e_demand = merge(6.0e-5_wp, 0.0_wp, mod(step, 2_ik) == 0_ik) * st%static%nplant              &
                             / sum(st%static%nplant)
         call advance(st, sopts, 900.0_wp)
         nsub = nsub + st%work(WORK_WATER_NSUB) + st%work(WORK_WATER_REJECT)
      end do
      call check_true('stiff exchange: demand jumps cost at most 6 sub-steps a step', nsub <= 48_ik, real(nsub, wp))
   end subroutine test_stiff_exchange

   subroutine test_infeasible_demand()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), demand
      theta = 0.15_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.6_wp, 0.5_wp, 0.0_wp)
      st%static%bleaf(1) = 1.0e-3_wp                            ! a few grams of leaf
      st%static%k_xylem(1) = 1.0e-3_wp * st%static%k_xylem(1)   ! and a xylem that can hardly feed it
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice,                                                &
                           0.0_wp, 0.0_wp, [1.0e-3_wp * 2.0_wp * 0.75_wp],         &
                           [st%water%state%water_coh(WOOD, 1) / st%static%nplant(1)], st%water)
      demand = 1.0e-4_wp
      st%drive%e_demand(1) = demand
      call advance(st, sopts, 900.0_wp)
      associate (state => st%water%state)
         call check_int('infeasible demand: no sub-step forced through', st%work(WORK_WATER_NFAIL), 0_ik)
         call check_true('infeasible demand: the leaf stays above empty', state%tau_coh(LEAF, 1) > 0.0_wp,      &
                         state%tau_coh(LEAF, 1))
         call check_true('infeasible demand: transpiration falls short of the demand',                        &
                         st%water%flows%transp(1) < 0.5_wp * demand * 900.0_wp,                                &
                         st%water%flows%transp(1) / (demand * 900.0_wp))
         call check_true('infeasible demand: leaf potential is finite', abs(state%psi_coh(LEAF, 1)) < 1.0e6_wp, &
                         state%psi_coh(LEAF, 1))
      end associate
   end subroutine test_infeasible_demand

   subroutine test_leafless()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active)
      theta = 0.30_wp
      call make_stand(st, 3_ik, theta, 0.0_wp, 0.9_wp, 0.8_wp, 2.0e-5_wp, leafless=.true.)
      st%drive%e_demand(3) = 0.0_wp
      call advance(st, sopts, 900.0_wp)
      call check_close('leafless: leaf at the wood less the crown head', st%water%state%psi_coh(LEAF, 3),        &
                       st%water%state%psi_coh(WOOD, 3) - st%static%crown_head(3), 1.0e-6_wp)
      call check_close('leafless: its empty leaf store stays empty', st%water%state%water_coh(LEAF, 3), 0.0_wp,   &
                       1.0e-12_wp)
   end subroutine test_leafless

   !---~---
   !   A dry-layer storm, from a BCI regrowth stand: the top 2.6 cm a hair above residual, a
   ! stand whose wood is half empty and refilling, then 30 mm/h of rain for two steps.
   !---~---
   subroutine test_dry_layer_storm()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, flows, theta_min
      integer(ik) :: step, nsl
      logical :: finite
      nsl = soil%n_active
      theta = 0.16_wp
      theta(1) = 0.0981_wp + 2.0e-3_wp
      call make_stand(st, 6_ik, theta, 0.0_wp, 0.45_wp, 0.40_wp, 0.0_wp)
      w0 = total_water(st)
      flows = 0.0_wp
      finite = .true.
      theta_min = 1.0_wp
      do step = 1_ik, 6_ik
         st%drive%rain = merge(30.0_wp / 3600.0_wp, 0.0_wp, step == 3 .or. step == 4)
         call advance(st, sopts, 900.0_wp)
         flows = flows + boundary_water(st, 900.0_wp)
         associate (state => st%water%state)
            finite = finite .and. all(state%water_soil == state%water_soil) .and. all(state%psi_coh == state%psi_coh)
            theta_min = min(theta_min, minval(state%water_soil(2:nsl+1) / (rho_h2o * soil%dz(1:nsl))           &
                                              - soil%theta_res(1:nsl)))
         end associate
         call check_int('dry-layer storm: no sub-step forced through', st%work(WORK_WATER_NFAIL), 0_ik)
      end do
      call check_true('dry-layer storm: every state finite', finite)
      call check_true('dry-layer storm: no layer below residual', theta_min > 0.0_wp, theta_min)
      call check_close('dry-layer storm: water conserved', total_water(st) - w0, flows, 1.0e-10_wp * w0)
      call check_true('dry-layer storm: the wood refilled', st%water%state%tau_coh(WOOD, 1) > 0.40_wp,          &
                      st%water%state%tau_coh(WOOD, 1))
   end subroutine test_dry_layer_storm

   subroutine test_saturated()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, flows, spill, over
      integer(ik) :: step, nsl
      nsl = soil%n_active
      theta = 0.4314_wp - 1.0e-3_wp
      call make_stand(st, 2_ik, theta, 0.0_wp, 0.97_wp, 0.95_wp, 1.0e-5_wp)
      w0 = total_water(st)
      flows = 0.0_wp
      spill = 0.0_wp
      over = 0.0_wp
      do step = 1_ik, 8_ik
         st%drive%rain = merge(50.0_wp / 3600.0_wp, 0.0_wp, step <= 4)
         call advance(st, sopts, 900.0_wp)
         flows = flows + boundary_water(st, 900.0_wp)
         spill = spill + st%water%flows%spill
         over = max(over, maxval(st%water%state%water_soil(2:nsl+1) / (rho_h2o * soil%dz(1:nsl))              &
                                 - soil%theta_sat(1:nsl)))
         call check_int('saturated: no sub-step forced through', st%work(WORK_WATER_NFAIL), 0_ik)
      end do
      call check_true('saturated: the pond spilled', spill > 0.0_wp, spill)
      call check_true('saturated: no layer more than 1e-3 above saturation', over < 1.0e-3_wp, over)
      call check_true('saturated: pond at or below its capacity', st%water%state%water_soil(POND) <= sopts%w_pond_max, &
                      st%water%state%water_soil(POND))
      call check_close('saturated: water conserved', total_water(st) - w0, flows, 1.0e-10_wp * w0)
   end subroutine test_saturated

   !---~---
   !   A transpiring day with a midday storm at the default tolerances, against the same day at
   ! tolerances a thousand times tighter.
   !---~---
   subroutine test_accuracy()
      type(stand_t) :: st, ref
      real(wp) :: theta(soil%n_active), hour, worst_theta, worst_rwc
      integer(ik) :: step, nsl, isub
      nsl = soil%n_active
      theta = 0.30_wp
      call make_stand(st, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      call make_stand(ref, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      ref%tols%rtol(GRP_THETA) = 1.0e-3_wp * WATER_RTOL
      ref%tols%atol(GRP_THETA) = 1.0e-3_wp * WATER_ATOL
      isub = 0_ik
      do step = 1_ik, 96_ik
         hour = real(step, wp) * 0.25_wp
         st%drive%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                     &
                             * st%static%nplant * st%static%bleaf / sum(st%static%nplant * st%static%bleaf)
         ref%drive%e_demand = st%drive%e_demand
         st%drive%rain = merge(1.0e-2_wp, 0.0_wp, step >= 50 .and. step < 54)
         ref%drive%rain = st%drive%rain
         call advance(st, sopts, 900.0_wp)
         call advance(ref, sopts, 900.0_wp)
         isub = isub + st%work(WORK_WATER_NSUB)
      end do
      worst_theta = maxval(abs(st%water%state%water_soil(2:nsl+1) - ref%water%state%water_soil(2:nsl+1))         &
                           / (rho_h2o * soil%dz(1:nsl)))
      worst_rwc   = maxval(abs(st%water%state%tau_coh - ref%water%state%tau_coh))
      call check_true('accuracy: soil water within 10 atol of the tight run', worst_theta < 10.0_wp * WATER_ATOL,  &
                      worst_theta)
      call check_true('accuracy: tissue water within 1e-2 RWC of the tight run', worst_rwc < 1.0e-2_wp, worst_rwc)
      print '(a,i6)', '  info : sub-steps over the day at the default tolerances: ', isub
   end subroutine test_accuracy

   !---~---
   !   Two stores a BCI regrowth stand reaches: a cohort with 1e-22 kgC of leaves,
   ! whose leaf store holds 1e-26 kg/m2 against its plant's 1e-9 tolerance, and a pond holding a
   ! trace of water over a dry top layer in light rain, which drains tens of millimetres an hour.
   ! A day of them: no sub-step forced, the leaf never above saturation, the pond never below
   ! empty.
   !---~---
   subroutine test_vanishing_stores()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), leaf_w(3), wood_w(3), w0, flows, rwc_max, pond_min
      integer(ik) :: step, nfail
      theta = 0.16_wp
      call make_stand(st, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      st%static%bleaf(3) = 1.3e-22_wp
      leaf_w = st%water%state%water_coh(LEAF, :) / st%static%nplant
      wood_w = st%water%state%water_coh(WOOD, :) / st%static%nplant
      leaf_w(3) = 2.0_wp * st%static%bleaf(3) * 0.9_wp
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, 3.0e-8_wp, 0.0_wp, leaf_w, wood_w, st%water)
      w0 = total_water(st)
      flows = 0.0_wp
      nfail = 0_ik
      rwc_max = 0.0_wp
      pond_min = 0.0_wp
      do step = 1_ik, 96_ik
         st%drive%rain = merge(1.0e-6_wp, 1.2e-3_wp, step < 40)
         st%drive%e_demand(1:2) = 2.0e-5_wp
         st%drive%e_demand(3) = 1.0e-32_wp
         call advance(st, sopts, 900.0_wp)
         flows = flows + boundary_water(st, 900.0_wp)
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         rwc_max  = max(rwc_max, st%water%state%tau_coh(LEAF, 3))
         pond_min = min(pond_min, st%water%state%water_soil(POND))
      end do
      call check_int('vanishing stores: no sub-step forced through', nfail, 0_ik)
      call check_true('vanishing stores: the near-leafless leaf stays at or below saturation', rwc_max <= 1.0_wp, rwc_max)
      call check_true('vanishing stores: the pond never holds less than nothing', pond_min >= 0.0_wp, pond_min)
      call check_close('vanishing stores: water conserved', total_water(st) - w0, flows, 1.0e-10_wp * w0)
   end subroutine test_vanishing_stores

   !---~---
   !   Two transpiring days without rain on a Campbell clay (psi_sat = -0.26 m, b = 5.65), on
   ! the van Genuchten clay's layers. The column carries one parameter pair whose meaning the
   ! family sets, so the network must read (psi_sat, b) through it: each layer's potential is
   ! the Campbell curve's at its water content, finite and below zero, the stand draws water
   ! from the soil, and every step closes. No layer reaches saturation, where the curve reads
   ! zero and the network a head above air entry. The fixture's soil is swapped for the test
   ! and put back after it.
   !---~---
   subroutine test_campbell()
      type(stand_t)        :: st
      type(soil_profile_t) :: vg_soil
      real(wp), parameter  :: THETA_SAT = 0.44_wp, THETA_RES = 0.0_wp, PSI_SAT = -0.26_wp, B = 5.65_wp
      real(wp) :: bottoms(n_soil_layer_max), theta(soil%n_active), psi_curve(soil%n_active)
      real(wp) :: w0, w1, dt, hour, worst, worst_psi, uptake
      integer(ik) :: step, nfail, nsl
      logical :: below_zero
      vg_soil = soil
      nsl = vg_soil%n_active
      call geometric_layer_bottoms(1.44_wp, nsl, 3.0_wp, bottoms)
      call build_profile_hydraulics(nsl, SOIL_RETENTION_CAMPBELL, bottoms, THETA_SAT, THETA_RES, 4.54e-6_wp,    &
                                  PSI_SAT, B, -3.37_wp, soil)
      theta = 0.30_wp
      call make_stand(st, 4_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp)
      dt = 900.0_wp
      worst = 0.0_wp
      worst_psi = 0.0_wp
      uptake = 0.0_wp
      nfail = 0_ik
      below_zero = .true.
      do step = 1_ik, 192
         hour = mod(real(step, wp) * dt / 3600.0_wp, 24.0_wp)
         st%drive%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                     &
                             * st%static%nplant * st%static%bleaf / sum(st%static%nplant * st%static%bleaf)
         st%drive%rain = 0.0_wp
         w0 = total_water(st)
         call advance(st, sopts, dt)
         w1 = total_water(st)
         worst = max(worst, abs((w1 - w0) - boundary_water(st, dt)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         uptake = uptake + sum(st%water%flows%uptake)
         !---~---
         !   Each layer's potential against the Campbell curve at the water it holds, which the
         ! sub-step commits from its flows: the two agree to the Newton tolerance.
         !---~---
         associate (state => st%water%state)
            theta = state%water_soil(2:nsl+1) / (rho_h2o * soil%dz(1:nsl))
            psi_curve = grav_head * soil_theta_to_head(SOIL_RETENTION_CAMPBELL, theta, THETA_SAT, THETA_RES,  &
                                                        PSI_SAT, B)
            worst_psi = max(worst_psi, maxval(abs(state%psi_soil(2:nsl+1) - psi_curve) / abs(psi_curve)))
            below_zero = below_zero .and. all(state%psi_soil(2:nsl+1) < 0.0_wp)                              &
                                    .and. all(state%psi_soil(2:nsl+1) > -huge(1.0_wp))
         end associate
      end do
      soil = vg_soil
      call check_true('Campbell: every soil potential is finite and below zero', below_zero)
      call check_true('Campbell: every soil potential is the Campbell curve''s', worst_psi < 1.0e-5_wp, worst_psi)
      call check_true('Campbell: the stand takes up water from the soil', uptake > 0.0_wp, uptake)
      call check_true('Campbell: every step closes to round-off', worst < 1.0e-13_wp, worst)
      call check_int('Campbell: no sub-step forced through', nfail, 0_ik)
   end subroutine test_campbell

   !--- The soil column alone.

   !---~---
   !   The column's own scenarios, run with no stand on a 2 m loam of ten layers in place of the
   ! fixture's clay, which is put back after them.
   !---~---
   subroutine test_column_scenarios()
      type(soil_profile_t) :: clay
      real(wp) :: bottoms(n_soil_layer_max)
      clay = soil
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, bottoms, 0.43_wp, 0.078_wp,               &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soil)
      call test_column_drainage()
      call test_downpour()
      call test_dried_surface_rewets()
      call test_pond_subfreezing_inflow()
      call test_aquifer()
      call test_frozen_top()
      call test_frozen_surface_evaporation()
      soil = clay
   end subroutine test_column_scenarios

   !---~---
   !   No cohorts on the current soil at water content theta, an empty pond, the water it gives up
   ! and its rain at 290 K, under air wetter than any surface, so nothing evaporates.
   !---~---
   subroutine make_column(st, theta)
      type(stand_t), intent(out) :: st
      real(wp),      intent(in)  :: theta(:)
      call make_stand(st, 0_ik, theta, 0.0_wp, 0.9_wp, 0.9_wp, 0.0_wp)
      st%drive%t_soil = 290.0_wp
      st%drive%t_rain = 290.0_wp
      st%drive%ground%t_ground = 290.0_wp
      st%drive%ground%q_air    = 0.05_wp
      st%drive%ground%rho_air  = 1.2_wp
      st%drive%ground%r_aero   = 100.0_wp
   end subroutine make_column

   !--- make_column's, with the soil's evaporation suppressed by a huge aerodynamic resistance.
   subroutine aquifer_column(st, theta)
      type(stand_t), intent(out) :: st
      real(wp),      intent(in)  :: theta(:)
      call make_column(st, theta)
      st%drive%ground%r_aero = 1.0e6_wp
   end subroutine aquifer_column

   !---~---
   !   Freeze the column's top nfrozen layers at temperature temp: each keeps the liquid its
   ! freezing curve holds there, the rest ice, and the network restarts from that state.
   !---~---
   subroutine freeze_top(st, theta, temp, nfrozen)
      type(stand_t), intent(inout) :: st
      real(wp),      intent(in)    :: theta(:), temp
      integer(ik),   intent(in)    :: nfrozen
      real(wp)    :: theta_liq, dtheta, dum(0)
      integer(ik) :: k
      do k = 1_ik, nfrozen
         call soil_liquid_capacity(soil%layer_curve(k), temp, theta_liq, dtheta)
         st%drive%ice(k) = soil_ice_state(soil%layer_curve(k), max(theta(k) - theta_liq, 0.0_wp))
      end do
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, 0.0_wp, 0.0_wp, dum, dum, st%water)
   end subroutine freeze_top

   !---~---
   !   A wet column draining through its base (docs/science/soil_biophysics.md §Boundary
   ! conditions): freely, where its water closes against the drainage, and through a leaky
   ! base at f = 0.05, the ED soil maps' factor, where it drains about f of what the free one
   ! does and closes as well.
   !---~---
   subroutine test_column_drainage()
      type(soil_params_t) :: opts
      type(stand_t)       :: col
      real(wp), parameter :: F_LEAK = 0.05_wp
      real(wp) :: theta(10), w0, drained, drained_free
      integer(ik) :: step
      opts%water_bottom_bc = SOIL_BC_DRAIN
      theta = 0.40_wp
      !--- Free drainage.
      opts%drain_factor = 1.0_wp
      call make_column(col, theta)
      w0 = total_water(col)
      drained_free = 0.0_wp
      do step = 1_ik, 20_ik
         call advance(col, opts, 600.0_wp)
         drained_free = drained_free + col%water%flows%drain
      end do
      call check_true('free drain: a wet column drains', drained_free > 0.0_wp, drained_free)
      call check_close('free drain: its water closes against the drainage', w0 - total_water(col), drained_free,  &
                       1.0e-10_wp * w0)
      !--- A leaky base.
      opts%drain_factor = F_LEAK
      call make_column(col, theta)
      drained = 0.0_wp
      do step = 1_ik, 20_ik
         call advance(col, opts, 600.0_wp)
         drained = drained + col%water%flows%drain
      end do
      call check_true('leaky drain: a leaky column drains, by about f of what a free one does',                &
                      drained > F_LEAK * drained_free .and. drained < 2.0_wp * F_LEAK * drained_free,            &
                      drained / drained_free)
      call check_close('leaky drain: its water closes against the drainage', w0 - total_water(col), drained,      &
                       1.0e-10_wp * w0)
   end subroutine test_column_drainage

   !---~---
   !   WET clay (low Ksat, small suction) under a downpour: infiltration is limited by the
   ! clay's conductivity, so the excess ponds and, past the pond's capacity, runs off
   ! (Hortonian). A bone-dry clay would instead take it up by suction (Green-Ampt).
   !---~---
   subroutine test_downpour()
      type(soil_profile_t) :: loam
      type(soil_params_t)  :: opts
      type(stand_t)        :: col
      real(wp) :: theta(10), w0, rain, infil, runoff, drained, bottoms(n_soil_layer_max)
      integer(ik) :: step
      loam = soil
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, bottoms, 0.38_wp, 0.068_wp,          &
           5.6e-7_wp, 0.8_wp, 1.09_wp, -3.37_wp, soil)
      theta = 0.36_wp
      call make_column(col, theta)
      w0 = total_water(col)
      rain = 0.0_wp
      infil = 0.0_wp
      runoff = 0.0_wp
      drained = 0.0_wp
      do step = 1_ik, 4_ik
         col%drive%rain = 1.0e-2_wp                                     ! 36 mm/hr
         call advance(col, opts, 900.0_wp)
         rain    = rain + col%drive%rain * 900.0_wp
         infil   = infil + col%water%flows%infil - col%water%flows%exfil
         runoff  = runoff + col%water%flows%spill
         drained = drained + col%water%flows%drain
      end do
      call check_true('downpour: infiltration < rain (conductivity-limited)', infil < rain, infil / rain)
      call check_true('downpour: the pond filled to its capacity',                                  &
                      col%water%state%water_soil(POND) >= opts%w_pond_max - 1.0e-9_wp, col%water%state%water_soil(POND))
      call check_true('downpour: and the rest ran off', runoff > 0.0_wp, runoff)
      call check_close('downpour: rain = change in store + runoff + drainage', total_water(col) - w0,             &
                       rain - runoff - drained, 1.0e-10_wp * rain)
      soil = loam
   end subroutine test_downpour

   !---~---
   !   A top layer dried to residual takes up 2 mm/h of rain: within 12 h it has re-wetted
   ! and almost none of the rain has run off. The pond face's conductivity is the geometric
   ! mean of the saturated pond's and the top layer's, and the dry layer's suction drives it;
   ! the top layer's conductivity alone would seal a dried surface against rain.
   !---~---
   subroutine test_dried_surface_rewets()
      type(soil_params_t)    :: opts
      type(stand_t)          :: col
      real(wp) :: theta(10), rain, runoff
      integer(ik) :: step
      theta = 0.25_wp
      theta(1) = 0.08_wp
      call make_column(col, theta)
      col%drive%ground%t_ground = 298.0_wp
      col%drive%ground%q_air = 0.02_wp
      rain = 0.0_wp
      runoff = 0.0_wp
      do step = 1_ik, 48_ik
         col%drive%rain = 2.0_wp / 3600.0_wp                            ! 2 mm/h
         call advance(col, opts, 900.0_wp)
         rain = rain + col%drive%rain * 900.0_wp
         runoff = runoff + col%water%flows%spill
      end do
      theta(1) = col%water%state%water_soil(2) / (rho_h2o * soil%dz(1))
      call check_true('rewetting: the dried top layer re-wets within 12 h of rain', theta(1) > 0.20_wp, theta(1))
      call check_true('rewetting: and under 5 % of the rain runs off', runoff < 0.05_wp * rain, runoff / rain)
   end subroutine test_dried_surface_rewets

   !---~---
   !   Rain (or sub-threshold snow routed to the ground) below the triple point joins the pond:
   ! what the pond received in enthalpy is what it gave the soil, kept, and lost as runoff,
   ! and the water that infiltrated carried its own enthalpy -- an effective temperature below
   ! the melt plateau, not u_liq at the plateau, which would create energy all winter.
   !---~---
   subroutine test_pond_subfreezing_inflow()
      type(soil_params_t)    :: opts
      type(stand_t)          :: col
      real(wp) :: theta(10), dt, e_in, t_infil
      theta = 0.25_wp
      call make_column(col, theta)
      col%drive%rain = 5.0e-6_wp
      col%drive%t_rain = 270.0_wp                                             ! 18 mm/day at 3.16 K below the triple point
      col%drive%ground%t_ground = 271.0_wp
      col%drive%ground%q_air = 0.003_wp
      dt = 150.0_wp
      call advance(col, opts, dt)
      e_in = col%drive%rain * dt * internal_energy_liquid(270.0_wp)
      associate (fl => col%water%flows, st => col%water%state)
         call check_true('sub-freezing inflow: all of it infiltrated (pond empty after the step)',       &
                         st%water_soil(POND) <= 1.0e-9_wp .and. fl%infil > 0.0_wp, st%water_soil(POND))
         call check_close('sub-freezing inflow: enthalpy to soil + left in pond + runoff == enthalpy received', &
                          fl%infil_enth - fl%exfil_enth + st%e_pond + fl%spill_enth, e_in, 1.0e-9_wp * abs(e_in))
      end associate
      call check_close('sub-freezing inflow: the rain brought exactly that enthalpy', col%water%flows%rain_enth,     &
                       e_in, 0.0_wp)
      t_infil = liquid_energy_to_temp(col%water%flows%infil_enth / col%water%flows%infil)
      call check_true('sub-freezing inflow: infiltrating water carries its own (sub-plateau) enthalpy',    &
                      t_infil < 273.0_wp, t_infil)
   end subroutine test_pond_subfreezing_inflow

   !---~---
   !   The aquifer bottom is head-driven against a saturated zone at the column base: a dry
   ! column draws water UP (free drainage never can), a wet one still drains down, the water
   ! closes with the flux running upward, and the column relaxes to hydrostatic equilibrium
   ! with the base (the bottom node's head -> minus half its thickness).
   !---~---
   subroutine test_aquifer()
      type(soil_params_t)    :: opts
      type(stand_t)          :: col
      real(wp) :: theta(10), drain_dry, drain_wet, w0, drained, theta_bot0, head_n, delta
      integer(ik) :: step, n
      n = 10_ik
      opts%water_bottom_bc = SOIL_BC_AQUIFER
      theta = 0.12_wp
      call aquifer_column(col, theta)
      call advance(col, opts, 600.0_wp)
      drain_dry = col%water%flows%drain
      call check_true('AQUIFER: dry column draws water UP from the saturated base (drainage < 0)',        &
                      drain_dry < 0.0_wp, drain_dry)
      theta = 0.40_wp
      call aquifer_column(col, theta)
      call advance(col, opts, 600.0_wp)
      drain_wet = col%water%flows%drain
      call check_true('AQUIFER: wet column still drains DOWNWARD (drainage > 0)', drain_wet > 0.0_wp, drain_wet)
      theta = 0.15_wp
      call aquifer_column(col, theta)
      w0 = total_water(col)
      drained = 0.0_wp
      theta_bot0 = 0.15_wp
      do step = 1_ik, 40_ik
         call advance(col, opts, 600.0_wp)
         drained = drained + col%water%flows%drain
      end do
      call check_close('AQUIFER: water closes while drawing upward', total_water(col) - w0, -drained,              &
                       1.0e-10_wp * w0)
      theta(n) = col%water%state%water_soil(1+n) / (rho_h2o * soil%dz(n))
      call check_true('AQUIFER: the dry column WETS from below', theta(n) > theta_bot0, theta(n))
      do step = 1_ik, 1200_ik
         call advance(col, opts, 3600.0_wp)
      end do
      head_n = col%water%state%psi_soil(1+n) / grav_head
      delta  = 0.5_wp * soil%dz(n)
      call check_true('AQUIFER: bottom layer reaches hydrostatic equilibrium (head_n -> -Delta)',           &
                      abs(head_n + delta) < 0.2_wp, head_n + delta)
      call check_true('AQUIFER: and the residual bottom flux has died',                                         &
                      abs(col%water%flows%drain / 3600.0_wp) < 1.0e-5_wp, col%water%flows%drain / 3600.0_wp)
   end subroutine test_aquifer

   !---~---
   !   A loam column with its top two layers at -7 C (docs/science/soil_biophysics.md, "Frozen
   ! soil"). Rain of 10 mm/h for two hours on a moist column: the ice's impedance lets in a
   ! tenth of what the thawed column takes, or less. Then a dry frozen top over wet soil and
   ! no rain for four days: the freezing curve's suction draws water up into the frozen layer,
   ! which fills to saturation and no further. The water closes in both.
   !---~---
   subroutine test_frozen_top()
      type(soil_params_t)  :: opts
      type(stand_t)        :: thawed, frozen
      real(wp) :: theta(10), rain, infil_thawed, infil_frozen, w0, runoff, drained, theta1, theta1_max
      integer(ik) :: step
      theta = 0.36_wp
      call make_column(thawed, theta)
      call make_column(frozen, theta)
      call freeze_top(frozen, theta, 266.0_wp, 2_ik)
      w0 = total_water(frozen)
      rain = 0.0_wp
      infil_thawed = 0.0_wp
      infil_frozen = 0.0_wp
      runoff = 0.0_wp
      drained = 0.0_wp
      do step = 1_ik, 8_ik
         thawed%drive%rain = 10.0_wp / 3600.0_wp
         frozen%drive%rain = thawed%drive%rain
         call advance(thawed, opts, 900.0_wp)
         call advance(frozen, opts, 900.0_wp)
         rain = rain + frozen%drive%rain * 900.0_wp
         infil_thawed = infil_thawed + thawed%water%flows%infil - thawed%water%flows%exfil
         infil_frozen = infil_frozen + frozen%water%flows%infil - frozen%water%flows%exfil
         runoff  = runoff + frozen%water%flows%spill
         drained = drained + frozen%water%flows%drain
      end do
      call check_true('frozen top: the ice impedes infiltration (frozen < thawed / 10)', infil_frozen < 0.1_wp * infil_thawed,  &
                      infil_frozen / max(infil_thawed, tiny(1.0_wp)))
      call check_close('frozen top: rain = change in store + runoff + drainage', total_water(frozen) - w0,         &
                       rain - runoff - drained, 1.0e-10_wp * rain)

      !--- A dry frozen top over wet soil: cryosuction, up to saturation.
      theta = 0.40_wp
      theta(1) = 0.20_wp
      call make_column(frozen, theta)
      call freeze_top(frozen, theta, 266.0_wp, 1_ik)
      w0 = total_water(frozen)
      drained = 0.0_wp
      theta1_max = 0.0_wp
      do step = 1_ik, 384_ik
         call advance(frozen, opts, 900.0_wp)
         drained = drained + frozen%water%flows%drain
         theta1 = frozen%water%state%water_soil(2) / (rho_h2o * soil%dz(1))
         theta1_max = max(theta1_max, theta1)
      end do
      call check_true('frozen top: the frozen layer draws water up from the wet soil below', theta1 > 0.25_wp, theta1)
      call check_true('frozen top: and fills no further than saturation', theta1_max <= soil%theta_sat(1) + 1.0e-6_wp,   &
                      theta1_max)
      call check_close('frozen top: cryosuction, change in store = - drainage', total_water(frozen) - w0, -drained, &
                       1.0e-10_wp * w0)
   end subroutine test_frozen_top

   !---~---
   !   The loam's top layer frozen, holding 0.20 of ice, under dry air at 268 K (docs/science/
   ! soil_biophysics.md, "Frozen soil"). The dry surface layer reads the liquid and
   ! dsl_theta_init of the ice, so it forms once the liquid falls below dsl_theta_init of the
   ! pores the ice leaves, d (theta_sat - theta_ice). With the liquid 0.02 above that, the
   ! layer evaporates as one with no dry surface layer, rho (alpha q_sat - q_air) / r_aero,
   ! alpha being Philip's humidity at its liquid's potential on the full pore space and q_sat
   ! over liquid; 0.02 below it, a dry surface layer holds it back. Reading the liquid alone,
   ! a dry surface layer would already form at the first; reading all the water, none at the
   ! second.
   !---~---
   subroutine test_frozen_surface_evaporation()
      type(soil_params_t)  :: opts
      type(stand_t)        :: col
      real(wp), parameter  :: THETA_ICE = 0.20_wp, OFFSET = 0.02_wp, T_GROUND = 268.0_wp, Q_AIR = 1.0e-4_wp
      real(wp) :: theta(10), onset, theta_liq, head, e(2), e_free(2), dum(0)
      integer(ik) :: k
      onset = opts%dsl_theta_init * (soil%theta_sat(1) - THETA_ICE)
      do k = 1_ik, 2_ik
         theta_liq = onset + merge(OFFSET, -OFFSET, k == 1_ik)
         theta = 0.30_wp
         theta(1) = theta_liq + THETA_ICE
         call make_column(col, theta)
         col%drive%ground%t_ground = T_GROUND
         col%drive%ground%q_air    = Q_AIR
         col%drive%ice(1) = soil_ice_state(soil%layer_curve(1), THETA_ICE)
         call set_water_state(col%static, soil, hyd, theta, col%drive%ice, 0.0_wp, 0.0_wp, dum, dum, col%water)
         call ground_evaporation_setup(soil, opts, col%drive%ground)
         call water_derivs(col%static, soil, opts, hyd, col%drive, col%water%state%tau_soil, col%water%state%tau_coh, &
                           col%water%ev(1), col%water%sys, .false.)
         e(k) = col%water%ev(1)%evap
         head = soil_theta_to_head(SOIL_RETENTION_VG, theta_liq, soil%theta_sat(1), soil%theta_res(1),              &
                                   soil%curve_par_a(1), soil%curve_par_n(1))
         e_free(k) = col%drive%ground%rho_air * (exp(head * grav / (r_wv * T_GROUND))                               &
                     * sat_specific_humidity(T_GROUND, p_std) - Q_AIR) / col%drive%ground%r_aero
      end do
      call check_close('frozen surface: liquid above d of the pores the ice leaves evaporates with no dry layer',     &
                       e(1), e_free(1), 1.0e-9_wp * e_free(1))
      call check_true('frozen surface: liquid below it is held back by a dry surface layer', e(2) < 0.9_wp * e_free(2), &
                      e(2) / e_free(2))
   end subroutine test_frozen_surface_evaporation

   !--- Small helpers.
   !--- A linear congruential generator, in 64-bit arithmetic so the product cannot overflow.
   real(wp) function pseudo(seed)
      integer(ik), intent(inout) :: seed
      integer(int64) :: s
      s = mod(1103515245_int64 * int(seed, int64) + 12345_int64, 2147483648_int64)
      seed = int(s, ik)
      pseudo = real(s, wp) / 2147483648.0_wp - 0.5_wp
   end function pseudo

   subroutine dense_solve(a_in, b, x, n)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: a_in(n, n), b(n)
      real(wp),    intent(out) :: x(n)
      real(wp) :: a(n, n+1), row(n+1), f
      integer(ik) :: i, j, p
      a(:, 1:n) = a_in
      a(:, n+1) = b
      do j = 1_ik, n
         p = j
         do i = j + 1_ik, n
            if (abs(a(i, j)) > abs(a(p, j))) p = i
         end do
         row = a(j, :)
         a(j, :) = a(p, :)
         a(p, :) = row
         do i = j + 1_ik, n
            f = a(i, j) / a(j, j)
            a(i, :) = a(i, :) - f * a(j, :)
         end do
      end do
      do i = n, 1_ik, -1_ik
         x(i) = (a(i, n+1) - sum(a(i, i+1:n) * x(i+1:n))) / a(i, i)
      end do
   end subroutine dense_solve

end program test_water_network
