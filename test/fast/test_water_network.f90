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
!   8a. REFILL BURST   : dwarf shrubs whose wood is 7 % short over saturated soil refill within
!                        the step, through roots far quicker than any sub-step, forcing none.
!   9. DRY-LAYER STORM : a thin top layer near residual, a stand refilling its wood, then a
!                        storm: converges, no store leaves its range.
!  10. SATURATED       : rain on a saturated column ponds, spills, and conserves.
!  11. ACCURACY        : the default tolerances track a tight-tolerance run.
!  12. VANISHING STORES: a near-leafless leaf and a trace pond stay in range, and water closes.
!  13. CAMPBELL        : on a Campbell clay the soil potentials are the Campbell curve's, finite
!                        and below zero; the stand takes up water, and every step closes.
!  13a. FROZEN ROOTS   : a stand with dry wood over a wet soil takes up water; with every layer
!                        frozen (no root contact) it takes up none.
!  13b. GRASS          : a tree over a grass sward, whose leaves join its fine roots with no
!                        stem between: its path is its leaves' conductance alone; two days of
!                        transpiration and a shower, and an autumn storm on dry soil that
!                        refills its roots within the hour, close with no sub-step forced.
! The soil column alone, on a 2 m loam (docs/science/soil_biophysics.md §The soil water column
! and §Boundary conditions):
!  14. DRAINAGE        : a wet column drains freely, and through a leaky base about f of that;
!                        its water closes against the drainage.
!  15. DOWNPOUR        : rain beyond a wet clay's infiltration ponds, and runs off past the pond.
!  16. REWETTING       : a dried top layer re-wets within 12 h of rain, with little runoff.
!  17. SUB-FREEZING    : rain below the triple point infiltrates as any rain: the network moves
!                        water only.
!  18. AQUIFER         : the head-driven bottom is two-way and settles to hydrostatic equilibrium.
!  19. FROZEN TOP      : the ice impedes the rain's infiltration; the freezing curve's suction
!                        fills a frozen layer to saturation, not past it; the water closes.
!  20. FROZEN SURFACE  : the dry surface layer reads the liquid and dsl_theta_init of the ice
!                        (water_derivs), so it forms below dsl_theta_init of the pores the ice
!                        leaves (CLM5's effective porosity); above that the layer evaporates as
!                        one with no dry surface layer, at its liquid's potential.
!  21. FROZEN SUBLIMATION: a frozen top layer under dry air loses more water than its liquid holds
!                        above residual: the ice melts to supply it, the liquid stays at the
!                        freezing curve's potential, no sub-step is forced, and the water closes.
!  22. DRY PEAT        : a steep peat's top layer a hair above residual under dry, warm air
!                        evaporates what it can, stays at or above residual, forces no sub-step,
!                        and closes its water.
!  23. ICE-SEALED      : rain on a wet column frozen hard at the top: the empty pond, cut off
!                        from the network, holds the rain and spills past its capacity; next to
!                        nothing infiltrates, no sub-step is forced, and the water closes.
! The node classes (docs/science/soil_plant_water.md §3), on thin layers between the pond and
! the soil (a litter's Campbell curve and a moss's van Genuchten one, built here):
!  24. ZERO STORAGE    : a saturated moss, which holds next to no more water at zero head, between
!                        a litter and the soil under rain: it is eliminated, no sub-step is forced,
!                        and every step closes.
!  25. FROZEN GROUP    : an empty pond over an ice-filled litter and moss (their liquid on their
!                        freezing curve) over a layer sealed by its ice, under rain: the
!                        factorisation never fails, the rain ponds, and the water closes.
!  26. LITTER STORM    : 1 mm of dry litter under 50 mm/h: no sub-step forced; the litter passes the
!                        rain to the soil; its class changes as it wets and drains, and the water
!                        closes across every change.
!  27. DRYING LITTER   : a litter over a soil frozen hard at the top, under dry canopy air, dries by
!                        the soil's evaporation law to its equilibrium with the air's humidity,
!                        above its residual, with no sub-step forced.
!  28. FULL POND       : a 5 kg/m2 pond over 1 mm of saturated litter under rain: no sub-step
!                        forced, and the water closes.
!  29. FROZEN TO POND  : a frozen litter filled past its pores pushes water up into the pond; the
!                        water closes.
!  30. SMALL POND      : a small pond over a wet litter and moss empties within a step, never past
!                        empty.
!  31. ADJACENT ZERO STORAGE: an empty pond over a litter whose pores its ice and liquid fill, over
!                        a thawed moss, under light rain: the pond and the litter, neither of
!                        which holds water below zero head, are eliminated together; no sub-step
!                        is forced, and the water closes.
!  32. CRYOSUCTION     : a litter frozen at 268 K over a thawed, saturated moss: its freezing
!                        curve's suction draws the moss's water up, but no further than its pores
!                        hold; no sub-step is forced, and the water closes.
!  33. FORCED BOUNDS   : what a forced sub-step's unconverged iterate would commit is kept only if
!                        every store is finite, at or above empty and near its curve's water.
!  34. SMOOTH FACE     : the flow between a frozen litter and a thawed moss and its slope are
!                        continuous as the flow reverses: one conductance both ways.
!  35. FREEZE-UP       : a wet litter and moss freeze over four days, the litter first, the ice
!                        each holds rebuilt each step from its temperature on its freezing curve:
!                        no sub-step forced, every store in its range, and every step closes.
!---~---
program test_water_network
   use ed_test_assert, only : check_true, check_close, check_int, test_report
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ed_kinds,         only : wp, ik
   use ed_constants,     only : grav_head, rho_h2o, pi, tiny_num, n_soil_layer_max, grav, r_wv, p_std, mmh2o
   use ed_linear_algebra,  only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   use ed_fast_integ_types,  only : fast_static_state_t, water_drive_t, water_work_t, tol_set_t,                    &
                                ensure_static_capacity, POND, WOOD, LEAF, GRP_THETA, N_WORK, WORK_WATER_NSUB,   &
                                WORK_WATER_REJECT, WORK_WATER_NFAIL, water_porous_t,                           &
                                WORK_WATER_FAST_SLOT, WORK_WATER_CUT_POND, WORK_WATER_CUT_SLOT,                &
                                WORK_WATER_FAST_POND, water_state_t, water_tend_t
   use ed_fast_state,    only : derive_cohort_terms
   use ed_fast_integ_derivs, only : set_water_state, water_derivs
   use ed_fast_integ_backward_euler, only : water_step, water_jacobian_error, substep_in_bounds
   use ed_plant_hydr_lib,   only : water_curve_t, tissue_psi_to_rwc, tissue_psi_tlp
   use ed_pft_derive,       only : build_kirchhoff_table
   use ed_soil_hydr_lib,     only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, soil_theta_to_head, soil_ice_state,  &
                                    soil_liquid_capacity, soil_water_relations, build_soil_layer_curve, soil_ice_t, &
                                    soil_layer_curve_t
   use ed_therm_lib,         only : sat_specific_humidity
   use ed_pft_params,    only : hydraulics_table_t
   use ed_soil_params, only : soil_profile_t, soil_params_t, SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use ed_soil_derive,        only : build_profile_hydraulics
   use ed_test_soil_grid, only : geometric_layer_bottoms
   use ed_soil_water,    only : ground_evap_env_t, ground_evaporation_setup, slot_evaporation_setup
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
   call test_refill_burst()
   call test_dry_layer_storm()
   call test_saturated()
   call test_accuracy()
   call test_vanishing_stores()
   call test_campbell()
   call test_frozen_roots()
   call test_grass_sward()
   call test_column_scenarios()
   call test_zero_storage()
   call test_frozen_group()
   call test_litter_storm()
   call test_drying_litter()
   call test_full_pond()
   call test_frozen_slot_to_pond()
   call test_pond_drains_what_it_holds()
   call test_adjacent_zero_storage()
   call test_cryosuction()
   call test_forced_bounds()
   call test_smooth_face()
   call test_freeze_up()
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
      h%pft(1)%leaf_kmax = 10.0_wp
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
   ! share of e_total by leaf, and the network started from the given water. With grass, the
   ! last cohort is a 0.3 m grass sward of 2000 tillers per m2 at LAI 3.2: no stem, its fine
   ! roots its only store.
   !---~---
   subroutine make_stand(st, n, theta, w_pond, leaf_rwc, wood_rwc, e_total, leafless, grass, org)
      type(stand_t),     intent(out) :: st
      integer(ik),       intent(in)  :: n
      real(wp),          intent(in)  :: theta(:), w_pond, leaf_rwc, wood_rwc, e_total
      logical, optional, intent(in)  :: leafless, grass
      type(water_porous_t), optional, intent(in) :: org   ! thin layers between the pond and the soil
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
      if (present(grass)) then
         if (grass) then
            st%static%is_woody(n) = .false.
            st%static%height(n)   = 0.3_wp
            st%static%nplant(n)   = 2000.0_wp
            st%static%bleaf(n)    = 1.0e-4_wp
            st%static%broot(n)    = 1.0e-4_wp
            st%static%bsap(n)     = 0.0_wp
            st%static%sap_area(n) = 0.0_wp
         end if
      end if
      !--- Each cohort's leaf area at 16 m2/kgC, which the leaves' conductance scales with.
      st%static%leaf_area(1:n) = 16.0_wp * st%static%bleaf(1:n)
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
      st%drive%pond_capacity = 5.0_wp          ! [kg/m2] the shipped [ground_surface_layer.water].pond_capacity
      st%drive%ground = forc
      st%tols%rtol(GRP_THETA) = WATER_RTOL
      st%tols%atol(GRP_THETA) = WATER_ATOL
      if (present(org)) st%drive%porous = org
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, w_pond, leaf_w, wood_w, st%water,              &
                           st%drive%porous)
   end subroutine make_stand

   !---~---
   !   Thin layers over the soil, as the porous layers are built: a litter d_lit deep at
   ! theta_lit on Letts et al. (2000)'s Campbell curve over a moss d_moss deep at theta_moss on
   ! Zhang et al. (2010)'s van Genuchten one (a layer of no depth left out), frozen (all but the
   ! residual, which stays liquid) or not, with a held evaporation evap.
   !---~---
   function thin_layers(d_lit, theta_lit, d_moss, theta_moss, frozen) result(o)
      real(wp), intent(in) :: d_lit, theta_lit, d_moss, theta_moss
      logical,  intent(in) :: frozen
      type(water_porous_t) :: o
      real(wp) :: liq
      if (d_lit > 0.0_wp) then
         o%n = o%n + 1_ik
         call build_soil_layer_curve(SOIL_RETENTION_CAMPBELL, 0.98_wp, 0.01_wp, -0.0103_wp, 2.3_wp, 1.83e-3_wp,        &
                                     d_lit, o%curve(o%n))
         o%depth(o%n) = d_lit
         o%water(o%n) = theta_lit * rho_h2o * d_lit
         liq = theta_lit
         if (frozen) liq = min(theta_lit, frozen_liquid(o%curve(o%n)))
         o%ice(o%n) = soil_ice_state(o%curve(o%n), theta_lit - liq, liq)
      end if
      if (d_moss > 0.0_wp) then
         o%n = o%n + 1_ik
         call build_soil_layer_curve(SOIL_RETENTION_VG, 0.85_wp, 0.18_wp, 50.0_wp, 1.45_wp, 1.3e-3_wp, d_moss,        &
                                     o%curve(o%n))
         o%depth(o%n) = d_moss
         o%water(o%n) = theta_moss * rho_h2o * d_moss
         liq = theta_moss
         if (frozen) liq = min(theta_moss, frozen_liquid(o%curve(o%n)))
         o%ice(o%n) = soil_ice_state(o%curve(o%n), theta_moss - liq, liq)
      end if
      o%air   = slot_evaporation_setup(forc, 285.0_wp, 0.0_wp)
   end function thin_layers

   !--- The liquid a frozen layer holds at 268 K, on its freezing curve, as a frozen soil layer's.
   real(wp) function frozen_liquid(cv)
      type(soil_layer_curve_t), intent(in) :: cv
      real(wp) :: dtheta
      call soil_liquid_capacity(cv, 268.0_wp, frozen_liquid, dtheta)
   end function frozen_liquid

   !--- The water that crossed the network's boundary over its last step, the thin layers' evaporation too [kg/m2].
   real(wp) function boundary_water_thin(st, dt)
      type(stand_t), intent(in) :: st
      real(wp),      intent(in) :: dt
      boundary_water_thin = boundary_water(st, dt) - sum(st%water%flows%porous_evap)
   end function boundary_water_thin

   !---~---
   !   The stand's top nk soil layers frozen hard, their pores all but full of ice, their roots cut
   ! off, under a pond of w_pond, the network restarted.
   !---~---
   subroutine seal_top(st, nk, theta, w_pond)
      type(stand_t), intent(inout) :: st
      integer(ik),   intent(in)    :: nk
      real(wp),      intent(in)    :: theta(:), w_pond
      real(wp)    :: liq, th(size(theta))
      integer(ik) :: k
      th = theta
      do k = 1_ik, nk
         liq = soil%theta_res(k) + 1.0e-6_wp
         st%drive%ice(k) = soil_ice_state(soil%layer_curve(k), soil%theta_sat(k) - liq, liq)
         th(k) = soil%theta_sat(k)
      end do
      st%drive%root_contact(1:nk) = 0.0_wp
      call set_water_state(st%static, soil, hyd, th, st%drive%ice, w_pond,                                        &
                           st%water%state%water_coh(LEAF, :) / st%static%nplant(1:st%static%n),                    &
                           st%water%state%water_coh(WOOD, :) / st%static%nplant(1:st%static%n), st%water, st%drive%porous)
   end subroutine seal_top

   !--- The stand restarted from its network's end state, as the next step's start (the thin layers rebuilt).
   subroutine restart_stand(st, d_lit, d_moss, frozen)
      type(stand_t), intent(inout) :: st
      real(wp),      intent(in)    :: d_lit, d_moss
      logical,       intent(in)    :: frozen
      real(wp)    :: th_lit, th_moss
      integer(ik) :: n_porous, nsl, i
      n_porous = st%drive%porous%n
      nsl   = soil%n_active
      th_lit  = 0.0_wp
      th_moss = 0.0_wp
      i = 1_ik
      if (d_lit > 0.0_wp) then
         th_lit = st%water%state%water_soil(1_ik + i) / (rho_h2o * d_lit)
         i = i + 1_ik
      end if
      if (d_moss > 0.0_wp) th_moss = st%water%state%water_soil(1_ik + i) / (rho_h2o * d_moss)
      st%drive%porous = thin_layers(d_lit, th_lit, d_moss, th_moss, frozen)
      call set_water_state(st%static, soil, hyd,                                                                  &
                           st%water%state%water_soil(2+n_porous:1+n_porous+nsl) / (rho_h2o * soil%dz(1:nsl)), st%drive%ice,  &
                           st%water%state%water_soil(POND),                                                        &
                           st%water%state%water_coh(LEAF, :) / st%static%nplant(1:st%static%n),                    &
                           st%water%state%water_coh(WOOD, :) / st%static%nplant(1:st%static%n), st%water, st%drive%porous)
   end subroutine restart_stand

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
      st%static%k_plant(1) = 1.0e-3_wp * st%static%k_plant(1)   ! and a wood-to-leaf path that can hardly feed it
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice,                                                &
                           0.0_wp, [1.0e-3_wp * 2.0_wp * 0.75_wp],                 &
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
   !   A tundra stand of 10 cm dwarf shrubs, 4,500 stems per m2 with a trace of wood each, its
   ! wood 7 % short of full over a saturated soil and no transpiration. Its roots refill the wood
   ! in a fraction of a millisecond, and the one-way uptake shuts as the wood catches up, so at
   ! the step's end the Newton matrix no longer sees the exchange; the error filter's pass at the
   ! start still does. The step takes no sub-step forced through and ends with the wood full.
   !---~---
   subroutine test_refill_burst()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), wood_w(1), leaf_w(1)
      theta = soil%theta_sat(1:soil%n_active)
      call make_stand(st, 1_ik, theta, 0.0_wp, 1.0_wp, 0.93_wp, 0.0_wp)
      st%static%height(1)   = 0.1_wp
      st%static%nplant(1)   = 4500.0_wp
      st%static%bleaf(1)    = 2.0e-6_wp
      st%static%bsap(1)     = 3.0e-6_wp
      st%static%broot(1)    = 2.0e-5_wp
      st%static%sap_area(1) = 1.0e-6_wp
      st%static%leaf_area(1) = 16.0_wp * st%static%bleaf(1)
      call derive_cohort_terms(st%static, hyd, soil)
      associate (s => st%static, lc => hyd%pft(1)%leaf_curve, wc => hyd%pft(1)%wood_curve)
         leaf_w = s%bleaf(1:1) * lc%water_sat
         wood_w = s%bwood_water(1:1) * wc%water_sat * ((1.0_wp - wc%apoplast_frac) * 0.93_wp + wc%apoplast_frac)
      end associate
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, 0.0_wp, leaf_w, wood_w, st%water)
      call advance(st, sopts, 900.0_wp)
      call check_int('refill burst: no sub-step forced through', st%work(WORK_WATER_NFAIL), 0_ik)
      call check_true('refill burst: the wood refilled', st%water%state%tau_coh(WOOD, 1) > 0.99_wp,             &
                      st%water%state%tau_coh(WOOD, 1))
   end subroutine test_refill_burst

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
      call check_true('saturated: pond at or below its capacity', st%water%state%water_soil(POND) <= st%drive%pond_capacity, &
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
   subroutine test_grass_sward()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, flows, hour, k_leaf, psi_root_min
      integer(ik) :: step, nfail
      !---~---
      !   The stores start near the soil's potential, as a cohort is born. A sward's small root
      ! store behind its large root length relaxes to the soil within a step; from the fixture's
      ! default start, 0.9 MPa away, that first step forces sub-steps through.
      !---~---
      theta = 0.30_wp
      call make_stand(st, 2_ik, theta, 0.0_wp, 0.95_wp, 0.99_wp, 0.0_wp, grass=.true.)
      k_leaf = hyd%pft(1)%leaf_kmax * 1.0e-3_wp * mmh2o * st%static%leaf_area(2)
      call check_close('grass: its water path is its leaves'' conductance alone', st%static%k_plant(2),     &
                       st%static%nplant(2) * k_leaf, 1.0e-12_wp * st%static%k_plant(2))
      call check_close('grass: its only store is its fine roots', st%static%bwood_water(2), st%static%broot(2), 0.0_wp)
      w0 = total_water(st)
      flows = 0.0_wp
      nfail = 0_ik
      psi_root_min = 0.0_wp
      do step = 1_ik, 192_ik
         hour = mod(real(step, wp) * 900.0_wp / 3600.0_wp, 24.0_wp)
         st%drive%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                     &
                             * st%static%nplant * st%static%bleaf / sum(st%static%nplant * st%static%bleaf)
         st%drive%rain = merge(2.0e-3_wp, 0.0_wp, step >= 60 .and. step < 68)
         call advance(st, sopts, 900.0_wp)
         flows = flows + boundary_water(st, 900.0_wp)
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         psi_root_min = min(psi_root_min, st%water%state%psi_coh(WOOD, 2))
      end do
      call check_int('grass: no sub-step forced through', nfail, 0_ik)
      call check_close('grass: water conserved', total_water(st) - w0, flows, 1.0e-10_wp * w0)
      call check_true('grass: its root store stays on its curve, above -5 MPa', psi_root_min > -5.0_wp, psi_root_min)
      call check_true('grass: its leaves sit below its roots by the day''s draw',                              &
                      st%water%state%psi_coh(LEAF, 2) <= st%water%state%psi_coh(WOOD, 2), st%water%state%psi_coh(LEAF, 2))
      !---~---
      !   Autumn's first storm: a summer-dry soil (about -1.5 MPa), the stores at its potential,
      ! then 72 mm/h for two hours.
      !---~---
      theta = 0.16_wp
      call make_stand(st, 2_ik, theta, 0.0_wp, 0.893_wp, 0.685_wp, 0.0_wp, grass=.true.)
      w0 = total_water(st)
      flows = 0.0_wp
      nfail = 0_ik
      do step = 1_ik, 48_ik
         st%drive%e_demand = 0.0_wp
         st%drive%rain = merge(2.0e-2_wp, 0.0_wp, step >= 4 .and. step < 12)
         call advance(st, sopts, 900.0_wp)
         flows = flows + boundary_water(st, 900.0_wp)
         nfail = nfail + st%work(WORK_WATER_NFAIL)
      end do
      call check_int('grass storm: no sub-step forced through', nfail, 0_ik)
      call check_close('grass storm: water conserved', total_water(st) - w0, flows, 1.0e-10_wp * w0)
      call check_true('grass storm: the roots refill from the wetted soil', st%water%state%psi_coh(WOOD, 2) > -0.1_wp, &
                      st%water%state%psi_coh(WOOD, 2))
   end subroutine test_grass_sward

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
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, 3.0e-8_wp, leaf_w, wood_w, st%water)
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

   !---~---
   !   13a. A stand with dry wood over a wet clay, and no transpiration: its roots take up water.
   ! With every layer frozen its roots have no contact with the soil (root_soil_contact): no water
   ! passes, and the plants keep theirs to round-off.
   !---~---
   subroutine test_frozen_roots()
      type(stand_t) :: thawed, frozen
      real(wp) :: theta(soil%n_active), coh0, coh_t, coh_f
      integer(ik) :: step
      theta = 0.40_wp
      call make_stand(thawed, 2_ik, theta, 0.0_wp, 0.5_wp, 0.45_wp, 0.0_wp)
      call make_stand(frozen, 2_ik, theta, 0.0_wp, 0.5_wp, 0.45_wp, 0.0_wp)
      frozen%drive%root_contact = 0.0_wp
      coh0  = sum(frozen%water%state%water_coh)
      do step = 1_ik, 4_ik
         call advance(thawed, sopts, 900.0_wp)
         call advance(frozen, sopts, 900.0_wp)
      end do
      coh_t  = sum(thawed%water%state%water_coh)
      coh_f  = sum(frozen%water%state%water_coh)
      call check_true('frozen roots: in contact, the dry wood takes up water', coh_t > coh0 * (1.0_wp + 1.0e-6_wp),      &
                      coh_t / coh0)
      call check_close('frozen roots: frozen, the plants keep their water', coh_f, coh0, 1.0e-12_wp * coh0)
   end subroutine test_frozen_roots

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
      call test_ice_sealed()
      call test_frozen_surface_evaporation()
      call test_frozen_sublimation()
      call test_dry_peat()
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
         st%drive%ice(k) = soil_ice_state(soil%layer_curve(k), max(theta(k) - theta_liq, 0.0_wp),                 &
                                          min(theta(k), theta_liq))
      end do
      call set_water_state(st%static, soil, hyd, theta, st%drive%ice, 0.0_wp, dum, dum, st%water)
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
                      col%water%state%water_soil(POND) >= col%drive%pond_capacity - 1.0e-9_wp, col%water%state%water_soil(POND))
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
   !   Rain below the triple point reaches the pond and infiltrates as any rain: the network moves
   ! water only, and the heat the water carries is the energy solves' (the ground surface's stack
   ! and the soil's heat column, test_surface and test_soil_heat).
   !---~---
   subroutine test_pond_subfreezing_inflow()
      type(soil_params_t)    :: opts
      type(stand_t)          :: col
      real(wp) :: theta(10), dt
      theta = 0.25_wp
      call make_column(col, theta)
      col%drive%rain = 5.0e-6_wp                                              ! 18 mm/day
      col%drive%ground%t_ground = 271.0_wp
      col%drive%ground%q_air = 0.003_wp
      dt = 150.0_wp
      call advance(col, opts, dt)
      associate (fl => col%water%flows, st => col%water%state)
         call check_true('sub-freezing inflow: all of it infiltrated (pond empty after the step)',       &
                         st%water_soil(POND) <= 1.0e-9_wp .and. fl%infil > 0.0_wp, st%water_soil(POND))
         call check_close('sub-freezing inflow: the pond''s face carried what reached it', fl%face(POND),       &
                          col%drive%rain * dt - st%water_soil(POND), 1.0e-12_wp)
      end associate
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
   !   Rain of 5 mm/h for two hours on a loam column at saturation whose top two layers are
   ! frozen hard, at -20 C (docs/science/soil_plant_water.md §3, "Isolated stores"): the
   ! freezing curve leaves the top layer a trace of liquid and the ice's impedance seals it, so
   ! the empty pond over it is cut off from the network. Solved on its own, the pond holds the
   ! rain and spills what passes its capacity; next to nothing infiltrates, no sub-step is
   ! forced through, and the water closes.
   !---~---
   subroutine test_ice_sealed()
      type(soil_params_t) :: opts
      type(stand_t)       :: col
      real(wp) :: theta(10), w0, rain, runoff, drained, infil
      integer(ik) :: step, nfail
      theta = soil%theta_sat(1:10)
      call make_column(col, theta)
      call freeze_top(col, theta, 253.15_wp, 2_ik)
      w0 = total_water(col)
      rain = 0.0_wp
      runoff = 0.0_wp
      drained = 0.0_wp
      infil = 0.0_wp
      nfail = 0_ik
      do step = 1_ik, 8_ik
         col%drive%rain = 5.0_wp / 3600.0_wp
         call advance(col, opts, 900.0_wp)
         rain = rain + col%drive%rain * 900.0_wp
         runoff = runoff + col%water%flows%spill
         drained = drained + col%water%flows%drain
         infil = infil + col%water%flows%infil - col%water%flows%exfil
         nfail = nfail + col%work(WORK_WATER_NFAIL)
      end do
      call check_int('ice-sealed: no sub-step forced through', nfail, 0_ik)
      call check_true('ice-sealed: next to nothing infiltrates', abs(infil) < 1.0e-3_wp * rain, infil / rain)
      call check_true('ice-sealed: the pond spilled', runoff > 0.0_wp, runoff)
      call check_close('ice-sealed: what did not run off ponds', col%water%state%water_soil(POND), rain - runoff,     &
                       1.0e-3_wp * rain)
      call check_close('ice-sealed: rain = change in store + runoff + drainage', total_water(col) - w0,           &
                       rain - runoff - drained, 1.0e-10_wp * rain)
   end subroutine test_ice_sealed

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
         col%drive%ice(1) = soil_ice_state(soil%layer_curve(1), THETA_ICE, theta_liq)
         call set_water_state(col%static, soil, hyd, theta, col%drive%ice, 0.0_wp, dum, dum, col%water)
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

   !---~---
   !   A loam column whose top layer is frozen at -2 C, under dry air for two days: its
   ! sublimation outruns the liquid it holds above residual (docs/science/soil_biophysics.md,
   ! "Frozen soil"). The ice melts to supply the water: the layer's liquid stays at the
   ! freezing curve's potential, the evaporation goes on, no sub-step is forced through, and
   ! the column's water closes against what crosses its boundary (the base drains freely).
   ! Drawing on the liquid alone, the layer would empty of liquid within a day and its water
   ! step would fail.
   !---~---
   subroutine test_frozen_sublimation()
      type(soil_params_t)  :: opts
      type(stand_t)        :: col
      real(wp), parameter  :: T_FROZEN = 271.15_wp
      real(wp) :: theta(10), theta_liq, dtheta, w0, evap, evap_first, psi0, buffer, crossed, dum
      integer(ik) :: step, nfail
      theta = 0.30_wp
      call make_column(col, theta)
      col%drive%t_soil = T_FROZEN
      col%drive%ground%t_ground = T_FROZEN
      col%drive%ground%q_air    = 0.0_wp
      col%drive%ground%r_aero   = 10.0_wp
      call freeze_top(col, theta, T_FROZEN, 1_ik)
      call ground_evaporation_setup(soil, opts, col%drive%ground)
      call soil_liquid_capacity(soil%layer_curve(1), T_FROZEN, theta_liq, dtheta)
      buffer = rho_h2o * soil%dz(1) * (theta_liq - soil%theta_res(1))     ! [kg/m2] liquid above residual
      w0 = total_water(col)
      call soil_water_relations(soil%layer_curve(1), col%drive%ice(1), col%water%state%tau_soil(2), psi0, dum,    &
                                dum, dum, dum, dum)
      evap = 0.0_wp
      crossed = 0.0_wp
      nfail = 0_ik
      do step = 1_ik, 192_ik
         call advance(col, opts, 900.0_wp)
         evap = evap + col%water%flows%evap
         crossed = crossed + boundary_water(col, 900.0_wp)
         if (step == 1_ik) evap_first = col%water%flows%evap
         nfail = nfail + col%work(WORK_WATER_NFAIL)
      end do
      call check_true('frozen sublimation: the layer loses more than its liquid above residual', evap > buffer,      &
                      evap / buffer)
      call check_true('frozen sublimation: and less than its ice', evap < rho_h2o * soil%dz(1) * (theta(1) - theta_liq), &
                      evap)
      call check_int('frozen sublimation: no sub-step forced through', nfail, 0_ik)
      call check_close('frozen sublimation: the liquid stays at the freezing curve''s potential',                    &
                       col%water%state%psi_soil(2), psi0, 1.0e-9_wp * abs(psi0))
      call check_true('frozen sublimation: the evaporation goes on', col%water%flows%evap > 0.5_wp * evap_first,      &
                      col%water%flows%evap / evap_first)
      call check_close('frozen sublimation: change in store = the water crossing the boundary (evaporation, drainage)', &
                       total_water(col) - w0, crossed, 1.0e-10_wp * w0)
   end subroutine test_frozen_sublimation

   !---~---
   !   22. A 2 m column of steep peat (van Genuchten n 3.34, theta_res 0.141: the Imnavait fit to
   ! Toolik's organic horizon) whose top layer starts 0.001 above residual, under dry air at 20 C
   ! (docs/science/soil_biophysics.md, "The soil water column"). Its pore air stays near saturation
   ! almost to residual, so its evaporation runs at the demand until the layer's own water above
   ! residual runs out and then tapers off (evaporation_supply): no sub-step is forced through,
   ! the layer stays at or above residual, and the column's water closes against its boundary.
   !---~---
   subroutine test_dry_peat()
      type(soil_params_t) :: opts
      type(stand_t)       :: col
      real(wp) :: bottoms(n_soil_layer_max), theta(10), w0, crossed, theta_min
      integer(ik) :: step, nfail
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, bottoms, 0.895_wp, 0.141_wp, 2.0e-6_wp, 1.03_wp,   &
                                    3.34_wp, -3.37_wp, soil)
      theta = 0.30_wp
      theta(1) = 0.141_wp + 1.0e-3_wp
      call make_column(col, theta)
      col%drive%ground%t_ground = 293.15_wp
      col%drive%ground%q_air    = 0.0_wp
      col%drive%ground%r_aero   = 10.0_wp
      call ground_evaporation_setup(soil, opts, col%drive%ground)
      w0 = total_water(col)
      crossed = 0.0_wp
      nfail = 0_ik
      theta_min = huge(1.0_wp)
      do step = 1_ik, 96_ik
         call advance(col, opts, 900.0_wp)
         crossed = crossed + boundary_water(col, 900.0_wp)
         nfail = nfail + col%work(WORK_WATER_NFAIL)
         theta_min = min(theta_min, col%water%state%water_soil(2) / (rho_h2o * soil%dz(1)))
      end do
      call check_int('dry peat: no sub-step forced through', nfail, 0_ik)
      call check_true('dry peat: the top layer stays at or above residual', theta_min >= 0.141_wp - 1.0e-9_wp,       &
                      theta_min)
      call check_close('dry peat: change in store = the water crossing the boundary', total_water(col) - w0,       &
                       crossed, 1.0e-10_wp * w0)
   end subroutine test_dry_peat

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

   !--- 24. A saturated 1 mm litter between a pond and the soil holds no more water: it is eliminated.
   subroutine test_zero_storage()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, worst
      integer(ik) :: step, nfail, nfast
      print '(a)', 'test_zero_storage:'
      theta = 0.33_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=thin_layers(3.0e-3_wp, 0.90_wp, 0.06_wp, 0.85_wp, .false.))
      st%drive%rain = 1.0e-3_wp
      worst = 0.0_wp
      nfail = 0_ik
      nfast = 0_ik
      do step = 1_ik, 8_ik
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         nfast = nfast + st%work(WORK_WATER_FAST_SLOT)
         call restart_stand(st, 3.0e-3_wp, 0.06_wp, .false.)
      end do
      call check_true('zero storage: the saturated moss is eliminated', nfast > 0_ik, real(nfast, wp))
      call check_true('zero storage: every step closes to round-off', worst < 1.0e-13_wp, worst)
      call check_int('zero storage: no sub-step forced through', nfail, 0_ik)
   end subroutine test_zero_storage

   !--- 25. An empty pond over an ice-filled litter and moss over a layer sealed by its ice.
   subroutine test_frozen_group()
      type(stand_t)       :: st
      real(wp) :: theta(soil%n_active), w0, worst
      integer(ik) :: step, nfail, ncut
      print '(a)', 'test_frozen_group:'
      theta = 0.42_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=thin_layers(2.0e-2_wp, 0.98_wp, 0.06_wp, 0.85_wp, .true.))
      call seal_top(st, 2_ik, theta, 0.0_wp)
      st%drive%rain = 3.0e-4_wp
      worst = 0.0_wp
      nfail = 0_ik
      ncut  = 0_ik
      do step = 1_ik, 4_ik
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         ncut  = ncut + st%work(WORK_WATER_CUT_POND) + st%work(WORK_WATER_CUT_SLOT)
      end do
      if (ncut > 0_ik) print '(a,i8)', '  info : frozen group: detached node-sub-steps ', ncut
      call check_int('frozen group: no sub-step forced through', nfail, 0_ik)
      call check_true('frozen group: the rain ponds', st%water%state%water_soil(POND) + st%water%flows%spill         &
                      > 0.9_wp * 3.0e-4_wp * 900.0_wp, st%water%state%water_soil(POND))
      call check_true('frozen group: the water closes', worst < 1.0e-13_wp, worst)
   end subroutine test_frozen_group

   !--- 26. 1 mm of dry litter under 50 mm/h for an hour, then two dry hours.
   subroutine test_litter_storm()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, worst, passed
      integer(ik) :: step, nfail, nfast, nnormal
      print '(a)', 'test_litter_storm:'
      theta = 0.25_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=thin_layers(1.0e-3_wp, 0.02_wp, 0.0_wp, 0.0_wp, .false.))
      worst = 0.0_wp
      nfail = 0_ik
      nfast = 0_ik
      nnormal = 0_ik
      passed = 0.0_wp
      do step = 1_ik, 12_ik
         st%drive%rain = merge(50.0_wp / 3600.0_wp, 0.0_wp, step <= 4_ik)
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         nfast = nfast + st%work(WORK_WATER_FAST_SLOT)
         nnormal = nnormal + st%work(WORK_WATER_NSUB) - st%work(WORK_WATER_FAST_SLOT) - st%work(WORK_WATER_CUT_SLOT)
         passed = passed + st%water%flows%face(1_ik + st%drive%porous%n)
         call restart_stand(st, 1.0e-3_wp, 0.0_wp, .false.)
      end do
      call check_int('litter storm: no sub-step forced through', nfail, 0_ik)
      call check_true('litter storm: the litter passes the rain to the soil', passed > 0.5_wp * 50.0_wp, passed)
      call check_true('litter storm: the litter is eliminated and not, as it wets and drains',                      &
                      nfast > 0_ik .and. nnormal > 0_ik, real(nfast, wp))
      call check_true('litter storm: the water closes across every class change', worst < 1.0e-13_wp, worst)
   end subroutine test_litter_storm

   !---~---
   !   27. A litter over a soil frozen hard at the top, under canopy air at 24 % humidity, dries by
   ! Philip's pore humidity to its equilibrium with the air, near theta 0.012 (a head of about -2e4
   ! m), above its residual. (Over unfrozen soil it would draw the soil's water up through its
   ! upstream face, and dry only to the soil's potential.)
   !---~---
   subroutine test_drying_litter()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), th, th_min
      integer(ik) :: step, nfail
      print '(a)', 'test_drying_litter:'
      theta = 0.30_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=thin_layers(3.0e-3_wp, 0.03_wp, 0.0_wp, 0.0_wp, .false.))
      call seal_top(st, 2_ik, theta, 0.0_wp)
      st%drive%ground%q_air = 0.005_wp
      st%drive%porous%air = slot_evaporation_setup(st%drive%ground, 285.0_wp, 0.0_wp)
      nfail = 0_ik
      th_min = 1.0_wp
      do step = 1_ik, 48_ik
         call advance(st, sopts, 900.0_wp)
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         th = st%water%state%water_soil(2) / (rho_h2o * 3.0e-3_wp)
         th_min = min(th_min, th)
         st%drive%porous = thin_layers(3.0e-3_wp, th, 0.0_wp, 0.0_wp, .false.)
         st%drive%porous%air = slot_evaporation_setup(st%drive%ground, 285.0_wp, 0.0_wp)
         call set_water_state(st%static, soil, hyd,                                                                &
                              st%water%state%water_soil(3:2+soil%n_active) / (rho_h2o * soil%dz(1:soil%n_active)), &
                              st%drive%ice, st%water%state%water_soil(POND),                                       &
                              st%water%state%water_coh(LEAF, :) / st%static%nplant(1:1),                          &
                              st%water%state%water_coh(WOOD, :) / st%static%nplant(1:1), st%water, st%drive%porous)
      end do
      call check_true('drying litter: it dries to its equilibrium with the air, near 0.012',                         &
                      th_min > 0.0105_wp .and. th_min < 0.0140_wp, th_min)
      call check_int('drying litter: no sub-step forced through', nfail, 0_ik)
   end subroutine test_drying_litter

   !--- 28. A 5 kg/m2 pond over 1 mm of saturated litter under rain.
   subroutine test_full_pond()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, worst
      integer(ik) :: step, nfail
      print '(a)', 'test_full_pond:'
      theta = 0.4314_wp - 1.0e-3_wp
      call make_stand(st, 1_ik, theta, 5.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=thin_layers(1.0e-3_wp, 0.98_wp, 0.0_wp, 0.0_wp, .false.))
      st%drive%rain = 6.0e-4_wp
      worst = 0.0_wp
      nfail = 0_ik
      do step = 1_ik, 4_ik
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         call restart_stand(st, 1.0e-3_wp, 0.0_wp, .false.)
      end do
      call check_int('full pond: no sub-step forced through', nfail, 0_ik)
      call check_true('full pond: the water closes', worst < 1.0e-13_wp, worst)
   end subroutine test_full_pond

   !---~---
   !   29. A frozen litter filled past its pores over ground frozen hard pushes water up into a pond:
   ! the water reaches the pond, and the network's water closes.
   !---~---
   subroutine test_frozen_slot_to_pond()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0
      print '(a)', 'test_frozen_slot_to_pond:'
      theta = 0.42_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=thin_layers(3.0e-3_wp, 1.0_wp, 0.0_wp, 0.0_wp, .false.))
      call seal_top(st, 2_ik, theta, 0.5_wp)
      st%drive%rain = 0.0_wp
      w0 = total_water(st)
      call advance(st, sopts, 900.0_wp)
      associate (f => st%water%flows)
         call check_true('frozen slot to pond: the litter pushes water up to the pond', f%exfil > 1.0e-6_wp, f%exfil)
         call check_close('frozen slot to pond: the pond''s face carries it, net', f%face(POND), -f%exfil,             &
                          1.0e-12_wp * f%exfil)
      end associate
      call check_true('frozen slot to pond: the water closes',                                                     &
                      abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) < 1.0e-13_wp * max(w0, 1.0_wp))
   end subroutine test_frozen_slot_to_pond

   !---~---
   !   30. A small pond (0.022 kg/m2), as at Toolik's spring, over 11 mm of litter at theta 0.9 and 6
   ! cm of saturated moss: it empties into the litter within a step, and loses no more than it
   ! held, never past empty.
   !---~---
   subroutine test_pond_drains_what_it_holds()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w_min
      integer(ik) :: step
      print '(a)', 'test_pond_drains_what_it_holds:'
      theta = 0.42_wp
      call make_stand(st, 1_ik, theta, 0.022_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                            &
                      org=thin_layers(1.14e-2_wp, 0.90_wp, 0.06_wp, 0.85_wp, .false.))
      st%drive%rain = 0.0_wp
      w_min = huge(1.0_wp)
      do step = 1_ik, 4_ik
         call advance(st, sopts, 900.0_wp)
         w_min = min(w_min, st%water%state%water_soil(POND))
         call restart_stand(st, 1.14e-2_wp, 0.06_wp, .false.)
      end do
      call check_true('pond drains what it holds: never past empty', w_min >= 0.0_wp, w_min)
   end subroutine test_pond_drains_what_it_holds

   !---~---
   !   Thin layers as thin_layers builds them, each holding the ice its freezing curve leaves it at
   ! its own temperature (t_lit, t_moss), the rest liquid, as the porous layers' energy gives it.
   !---~---
   function layers_at(d_lit, theta_lit, t_lit, d_moss, theta_moss, t_moss) result(o)
      real(wp), intent(in) :: d_lit, theta_lit, t_lit, d_moss, theta_moss, t_moss
      type(water_porous_t) :: o
      real(wp)    :: th(2), t(2), liq, dliq
      integer(ik) :: j
      o = thin_layers(d_lit, theta_lit, d_moss, theta_moss, .false.)
      th = [theta_lit, theta_moss]
      t  = [t_lit, t_moss]
      if (d_lit <= 0.0_wp) then
         th(1) = theta_moss
         t(1)  = t_moss
      end if
      do j = 1_ik, o%n
         call soil_liquid_capacity(o%curve(j), t(j), liq, dliq)
         liq = min(th(j), liq)
         o%ice(j) = soil_ice_state(o%curve(j), th(j) - liq, liq)
      end do
   end function layers_at

   !--- The stand restarted from its network's end state, its litter and moss at their temperatures.
   subroutine restart_at(st, d_lit, t_lit, d_moss, t_moss)
      type(stand_t), intent(inout) :: st
      real(wp),      intent(in)    :: d_lit, t_lit, d_moss, t_moss
      integer(ik) :: nsl
      nsl = soil%n_active
      st%drive%porous = layers_at(d_lit, st%water%state%water_soil(2) / (rho_h2o * d_lit), t_lit,                &
                                  d_moss, st%water%state%water_soil(3) / (rho_h2o * d_moss), t_moss)
      call set_water_state(st%static, soil, hyd, st%water%state%water_soil(4:3+nsl) / (rho_h2o * soil%dz(1:nsl)),  &
                           st%drive%ice, st%water%state%water_soil(POND),                                          &
                           st%water%state%water_coh(LEAF, :) / st%static%nplant(1:st%static%n),                    &
                           st%water%state%water_coh(WOOD, :) / st%static%nplant(1:st%static%n), st%water, st%drive%porous)
   end subroutine restart_at

   !---~---
   !   31. An empty pond over a 3 mm litter on its melt plateau, a tenth of its pores ice and the
   ! rest liquid, its pores full (it holds no more water below zero head), over a thawed 6 cm
   ! moss, under rain.
   !---~---
   function plateau_litter(theta_moss) result(o)
      real(wp), intent(in) :: theta_moss
      type(water_porous_t) :: o
      o = thin_layers(3.0e-3_wp, 0.98_wp, 0.06_wp, theta_moss, .false.)
      o%ice(1) = soil_ice_state(o%curve(1), 0.10_wp, 0.88_wp)
   end function plateau_litter

   !--- The litter's unknown halfway up the head its full pores rise through to zero, where it holds no more water.
   subroutine litter_in_full_band(st)
      type(stand_t), intent(inout) :: st
      associate (ice => st%drive%porous%ice(1))
         st%water%state%tau_soil(2) = ice%tau_full - 0.5_wp * ice%head_full / ice%dhead_full
      end associate
   end subroutine litter_in_full_band

   subroutine test_adjacent_zero_storage()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, worst
      integer(ik) :: step, nfail, npond, nslot
      print '(a)', 'test_adjacent_zero_storage:'
      theta = 0.40_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp, org=plateau_litter(0.80_wp))
      call litter_in_full_band(st)
      st%drive%rain = 1.0e-4_wp
      worst = 0.0_wp
      nfail = 0_ik
      npond = 0_ik
      nslot = 0_ik
      do step = 1_ik, 4_ik
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         npond = npond + st%work(WORK_WATER_FAST_POND)
         nslot = nslot + st%work(WORK_WATER_FAST_SLOT)
         !--- The next step starts as the first: the litter's pores full, the moss as it ended.
         st%drive%porous = plateau_litter(st%water%state%water_soil(3) / (rho_h2o * 0.06_wp))
         st%water%state%water_soil(2) = st%drive%porous%water(1)
         call set_water_state(st%static, soil, hyd,                                                                &
                              st%water%state%water_soil(4:3+soil%n_active) / (rho_h2o * soil%dz(1:soil%n_active)), &
                              st%drive%ice, st%water%state%water_soil(POND),                                       &
                              st%water%state%water_coh(LEAF, :) / st%static%nplant(1:1),                          &
                              st%water%state%water_coh(WOOD, :) / st%static%nplant(1:1), st%water, st%drive%porous)
         call litter_in_full_band(st)
      end do
      print '(a,2i8)', '  info : adjacent zero storage: pond and litter eliminated (node-sub-steps) ', npond, nslot
      call check_true('adjacent zero storage: the pond and the litter are both eliminated',                          &
                      npond > 0_ik .and. nslot > 0_ik, real(min(npond, nslot), wp))
      call check_int('adjacent zero storage: no sub-step forced through', nfail, 0_ik)
      call check_true('adjacent zero storage: every step closes to round-off', worst < 1.0e-13_wp, worst)
   end subroutine test_adjacent_zero_storage

   !--- 32. A 1 cm litter frozen at 268 K (liquid on its curve, the rest ice) over a thawed, saturated moss, a day.
   subroutine test_cryosuction()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, worst, w_lit, cap_lit, w_start
      integer(ik) :: step, nfail
      print '(a)', 'test_cryosuction:'
      theta = 0.40_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=layers_at(1.0e-2_wp, 0.60_wp, 268.0_wp, 0.06_wp, 0.85_wp, 275.0_wp))
      cap_lit = 0.98_wp * rho_h2o * 1.0e-2_wp
      w_start = st%water%state%water_soil(2)
      worst = 0.0_wp
      nfail = 0_ik
      w_lit = 0.0_wp
      do step = 1_ik, 96_ik
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         w_lit = max(w_lit, st%water%state%water_soil(2))
         call restart_at(st, 1.0e-2_wp, 268.0_wp, 0.06_wp, 275.0_wp)
      end do
      print '(a,2es12.4)', '  info : cryosuction: litter water from, to most ', w_start, w_lit
      call check_true('cryosuction: the frozen litter holds no more than its pores (to 10 m of head)',              &
                      w_lit <= cap_lit + rho_h2o * 1.0e-2_wp * 1.0e-4_wp * 10.0_wp, w_lit - cap_lit)
      call check_int('cryosuction: no sub-step forced through', nfail, 0_ik)
      call check_true('cryosuction: every step closes to round-off', worst < 1.0e-13_wp, worst)
   end subroutine test_cryosuction

   !--- 33. The bounds a forced sub-step's iterate must meet to be committed.
   subroutine test_forced_bounds()
      type(water_state_t) :: state
      type(water_tend_t)  :: ev
      real(wp) :: scale_soil(3), scale_coh(2, 1)
      print '(a)', 'test_forced_bounds:'
      scale_soil = 26.0_wp
      scale_coh  = 50.0_wp
      allocate(state%water_soil(3), state%water_coh(2, 1), ev%water_soil(3), ev%water_coh(2, 1))
      ev%water_soil  = [0.0_wp, 9.8_wp, 11.0_wp]
      ev%water_coh   = reshape([20.0_wp, 4.0_wp], [2, 1])
      state%water_soil = ev%water_soil
      state%water_coh  = ev%water_coh
      call check_true('forced bounds: an iterate on its curves is kept',                                             &
                      substep_in_bounds(state, ev, scale_soil, scale_coh, WATER_RTOL, WATER_ATOL), 1.0_wp)
      state%water_soil(1) = -1.0e-3_wp
      call check_true('forced bounds: a pond overdrawn is refused',                                                 &
                      .not. substep_in_bounds(state, ev, scale_soil, scale_coh, WATER_RTOL, WATER_ATOL), 1.0_wp)
      state%water_soil = ev%water_soil
      state%water_soil(2) = ev%water_soil(2) + 7.0_wp
      call check_true('forced bounds: a layer 7 kg past its curve''s water is refused',                             &
                      .not. substep_in_bounds(state, ev, scale_soil, scale_coh, WATER_RTOL, WATER_ATOL), 1.0_wp)
      state%water_soil = ev%water_soil
      state%water_coh(LEAF, 1) = -0.1_wp
      call check_true('forced bounds: a leaf below empty is refused',                                               &
                      .not. substep_in_bounds(state, ev, scale_soil, scale_coh, WATER_RTOL, WATER_ATOL), 1.0_wp)
      state%water_coh = ev%water_coh
      state%water_soil(3) = ieee_value(1.0_wp, ieee_quiet_nan)
      call check_true('forced bounds: a store not finite is refused',                                               &
                      .not. substep_in_bounds(state, ev, scale_soil, scale_coh, WATER_RTOL, WATER_ATOL), 1.0_wp)
   end subroutine test_forced_bounds

   !---~---
   !   34. A litter at 272.5 K whose pores ice and liquid fill (its conductivity fixed at its
   ! ice-impeded value, many decades below the moss's) over a saturated moss: the flow down their
   ! face, and its slope in the litter's head, on either side of the head at which it reverses.
   !---~---
   subroutine test_smooth_face()
      type(stand_t) :: st
      type(water_tend_t) :: ev
      real(wp) :: theta(soil%n_active), tau(soil%n_active + 3), q(4), heads(4), slope_below, slope_above, h_rev
      integer(ik) :: k
      print '(a)', 'test_smooth_face:'
      theta = 0.40_wp
      call make_stand(st, 0_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=layers_at(1.0e-2_wp, 0.98_wp, 272.5_wp, 0.06_wp, 0.85_wp, 275.0_wp))
      ev = st%water%ev(1)
      tau = st%water%state%tau_soil
      !--- The moss saturated; the face reverses with the litter 3.5 cm (the nodes' distance) below the moss's head.
      associate (cm => st%drive%porous%curve(2), ice => st%drive%porous%ice(1))
         tau(3) = cm%se_switch + cm%slope_switch * (0.0_wp - cm%head_switch)
         call water_derivs(st%static, soil, sopts, hyd, st%drive, tau, st%water%state%tau_coh, ev, st%water%sys, .false.)
         h_rev = ev%psi_soil(3) / grav_head - 0.5_wp * (1.0e-2_wp + 0.06_wp)
         heads = h_rev + [-4.0e-3_wp, -2.0e-3_wp, 2.0e-3_wp, 4.0e-3_wp]
         call check_true('smooth face: the heads lie where the litter''s pores are full', ice%head_full < heads(1),  &
                         ice%head_full)
         do k = 1_ik, 4_ik
            tau(2) = ice%tau_full + (heads(k) - ice%head_full) / ice%dhead_full
            call water_derivs(st%static, soil, sopts, hyd, st%drive, tau, st%water%state%tau_coh, ev, st%water%sys,    &
                              .false.)
            q(k) = ev%face(2)
         end do
      end associate
      print '(a,5es12.4)', '  info : smooth face: reversal head, flows ', h_rev, q
      slope_below = (q(2) - q(1)) / (heads(2) - heads(1))
      slope_above = (q(4) - q(3)) / (heads(4) - heads(3))
      call check_true('smooth face: the flow reverses between the heads', q(1) < 0.0_wp .and. q(4) > 0.0_wp, q(1))
      call check_close('smooth face: its slope is the same on either side of the reversal', slope_above / slope_below, &
                       1.0_wp, 1.0e-6_wp)
   end subroutine test_smooth_face

   !---~---
   !   35. A wet 11 mm litter and 6 cm moss over a wet soil freeze over four days, as at Toolik in
   ! autumn: the litter from 273.5 to 266 K, the moss a day behind it and slower, the ice each
   ! holds rebuilt every step from its temperature on its freezing curve, as the porous layers'
   ! energy gives it.
   !---~---
   subroutine test_freeze_up()
      type(stand_t) :: st
      real(wp) :: theta(soil%n_active), w0, worst, t_lit, t_moss, w_min, over
      integer(ik) :: step, nfail
      print '(a)', 'test_freeze_up:'
      theta = 0.40_wp
      call make_stand(st, 1_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp,                                              &
                      org=layers_at(1.1e-2_wp, 0.90_wp, 273.5_wp, 0.06_wp, 0.85_wp, 274.0_wp))
      worst = 0.0_wp
      nfail = 0_ik
      w_min = huge(1.0_wp)
      over  = -huge(1.0_wp)
      do step = 1_ik, 384_ik
         w0 = total_water(st)
         call advance(st, sopts, 900.0_wp)
         worst = max(worst, abs((total_water(st) - w0) - boundary_water_thin(st, 900.0_wp)) / max(w0, 1.0_wp))
         nfail = nfail + st%work(WORK_WATER_NFAIL)
         w_min = min(w_min, minval(st%water%state%water_soil(2:)))
         over  = max(over, st%water%state%water_soil(2) - 0.98_wp * rho_h2o * 1.1e-2_wp,                             &
                     st%water%state%water_soil(3) - 0.85_wp * rho_h2o * 0.06_wp)
         t_lit  = 273.5_wp - 7.5_wp * min(1.0_wp, real(step, wp) / 384.0_wp)
         t_moss = 274.0_wp - 5.0_wp * min(1.0_wp, max(0.0_wp, real(step - 96_ik, wp) / 288.0_wp))
         call restart_at(st, 1.1e-2_wp, t_lit, 0.06_wp, t_moss)
      end do
      call check_int('freeze-up: no sub-step forced through', nfail, 0_ik)
      call check_true('freeze-up: no layer below empty', w_min >= 0.0_wp, w_min)
      call check_true('freeze-up: neither layer past its pores (to 10 m of head)',                                 &
                      over <= rho_h2o * 0.06_wp * 1.0e-4_wp * 10.0_wp, over)
      call check_true('freeze-up: every step closes to round-off', worst < 1.0e-13_wp, worst)
   end subroutine test_freeze_up

end program test_water_network
