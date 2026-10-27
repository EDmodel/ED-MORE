! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_slow_ledger -- the site conservation ledger's ARITHMETIC and its STORE COVERAGE.
!
!   The ledger reports numbers nobody can check by inspection, so what it must never do is be
! quietly BLIND. Two ways it can be, and one test for each:
!
!   1. THE STORE MISSES A RESERVOIR. If slow_site_store forgets a field, every operator that
!      moves that field reads as perfectly conservative -- the ledger prints zeros and the
!      defect it exists to find is the one thing it cannot see. Fortran cannot enumerate a
!      derived type's components, so this is done the way test_state_combinators does it:
!      perturb each reservoir in turn and assert the total MOVES, by the right amount and in
!      the right currency. Adding a state field means adding a row here.
!
!   2. THE DECLARATION ARITHMETIC IS WRONG. A declared term is subtracted from the residual, so
!      a sign error would let a real leak cancel against a real source. Asserted directly.
!
!   The third test binds the canopy-air open-volume terms to the store they are supposed to
! explain: resize the control volume, and the three declared terms must equal the three store
! changes EXACTLY. That is the one declaration in the skeleton that does real arithmetic rather
! than passing a quantity through, and a production run closes it to 1e-17.
!---~---
program test_slow_ledger
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_soil_derive,          only : build_profile_hydraulics, soil_layers, build_soil_profile,           &
                                       soil_start_theta, top_soil_psi
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area
   use ed_patch_state_vars, only : cas_set_depth, patch_physical_state, litter_input_t
   use ed_biogeochem_params,    only : N_SOIL_POOL, IP_SLOW, SOIL_POOL_NAME
   use ed_soil_params,          only : soil_profile_t
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_cohort_state_vars,  only : init_cohort, set_cohort_size_from_carbon, hold_cohort_tissue_temps
   use ed_slow_update, only : terminate_cohorts, apply_patch_disturbance, fuse_2_patches, refresh_canopy_depth
   use ed_slow_dynamics, only : vegetation_dynamics, restructure_stand, settle_tissue_water
   use ed_slow_state, only : kill_report_t, kill_litter, shed_turnover_water, slow_step_data_t
   use ed_slow_ledger,          only : slow_store_t, slow_ledger_t, slow_site_store,               &
                                       slow_ledger_open, slow_ledger_declare, slow_ledger_mark,    &
                                       slow_fast_carbon_handover, slow_tissue_energy, SLOW_PHASE_GROW, &
                                       SLOW_PHASE_ALLOCATE, SLOW_PHASE_RECRUIT, SLOW_PHASE_TISSUE_WATER, &
                                       litter_carbon_total
   use ed_plant_hydr_lib,       only : water_curve_t, tissue_water_to_psi
   use ed_constants,            only : grav_head, n_soil_layer_max, umol_2_kgc
   use ed_test_support, only : banner, build_test_config, check, check_close
   use ed_constants,          only : rho_air_ref
   implicit none
   type(kill_report_t)  :: kills   ! what the cull and the treefall killed
   type(slow_ledger_t) :: no_ledger    ! inactive: a test step books nothing
   integer(ik) :: which   ! the member of a group a perturbation moves (section 1)
   character(len=8), parameter :: TISSUE(2) = [character(len=8) :: 'leaf', 'wood']
   character(len=13), parameter :: LIVE_POOL(4) = [character(len=13) :: 'leaf', 'fine-root', 'wood', &
                                                   'nonstructural']
   real(wp)    :: soil_bottoms(n_soil_layer_max)   ! [m] the column's layer bottoms, as the ledger builds them
   integer(ik) :: n_soil_layer

   real(wp), parameter  :: RHO = rho_air_ref    ! the density the canopy-air store is valued at
   type(ed_config_t)    :: cfg
   type(polygon_t)      :: poly
   type(soil_profile_t) :: soil
   type(slow_store_t)   :: s0, s1
   type(slow_ledger_t)  :: led
   real(wp)             :: de, dw, dc, expect, hand

   call banner('slow ledger: store coverage and declaration arithmetic')

   cfg = build_test_config()
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 2_ik, 0.3_wp, 18.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.8_wp,  6.0_wp)
   call finalize_init(poly)
   call check(poly%cohort%n == 2_ik .and. poly%patch%n == 1_ik, 'fixture: 2 cohorts in 1 patch')

   associate (sc => cfg%soil)
      call soil_layers(sc, soil_bottoms, n_soil_layer)
      call build_profile_hydraulics(n_soil_layer, sc%retention, soil_bottoms, sc%theta_sat, sc%theta_res, &
                                  sc%ksat, sc%curve_par_a, sc%curve_par_n, sc%psi_fc, soil)
   end associate

   !---~---
   !   Give every reservoir a non-zero value, so a perturbation below is a change to a store
   ! that already exists rather than the creation of one.
   !---~---
   poly%patch%soil_water(1)%theta(1:soil%n_active)      = 0.30_wp
   poly%patch%soil_energy(1)%soil_energy(1:soil%n_active) = 5.0e7_wp
   poly%patch%soil_water(1)%w_surface      = 2.0_wp
   poly%patch%soil_water(1)%w_surface_energy = 8.0e5_wp
   poly%patch%snow(1)%swe(1)           = 5.0_wp
   poly%patch%snow(1)%snow_energy(1)   = -1.0e6_wp
   poly%patch%cas(1)%can_shv      = 0.008_wp
   poly%patch%cas(1)%can_enthalpy = 2.9e5_wp
   poly%patch%cas(1)%can_co2      = 400.0_wp
   poly%patch%cas(1)%can_depth    = 20.0_wp
   poly%patch%soil_carbon(1)%carbon(IP_SLOW) = 3.0_wp
   poly%patch%recruit_pool(:,1)          = 0.0_wp
   poly%cohort%leaf_water_mass(1:2) = 0.02_wp
   poly%cohort%wood_water_mass(1:2) = 0.05_wp
   poly%cohort%leaf_surf_water(1:2) = 0.01_wp
   poly%cohort%wood_surf_water(1:2) = 0.03_wp
   call hold_cohort_tissue_temps(poly%cohort, 1_ik)   ! the tissues' energy with that water and film
   call hold_cohort_tissue_temps(poly%cohort, 2_ik)

   !---~---
   !   1. STORE COVERAGE: every member of every reservoir must reach the total, in its own
   ! currency -- each tissue, each live and CENTURY pool, each PFT's recruit pool.
   !---~---
   call moves_water ('soil theta',        perturb_theta,        0.01_wp * soil%dz(1) * 1000.0_wp)
   call moves_water ('pond w_surface',    perturb_pond,         1.0_wp)
   call moves_water ('snow swe',          perturb_swe,          1.0_wp)
   call moves_water ('canopy vapour',     perturb_shv,          RHO * 20.0_wp * 0.001_wp)
   do which = 1_ik, 2_ik
      call moves_water (trim(TISSUE(which))//' water', perturb_tissue_water, poly%cohort%nplant(1) * 0.01_wp)
      call moves_water (trim(TISSUE(which))//' film',  perturb_film,         poly%cohort%nplant(1) * 0.01_wp)
   end do
   call moves_energy('soil energy',       perturb_soil_energy,  1.0e6_wp * soil%dz(1))
   call moves_energy('pond energy',       perturb_pond_energy,  1.0e5_wp)
   call moves_energy('snow energy',       perturb_snow_energy,  1.0e5_wp)
   call moves_energy('canopy enthalpy',   perturb_cas_enth,     RHO * 20.0_wp * 1.0e4_wp)
   do which = 1_ik, 4_ik
      call moves_carbon(trim(LIVE_POOL(which))//' carbon', perturb_live_c, poly%cohort%nplant(1) * 0.5_wp)
   end do
   do which = 1_ik, N_SOIL_POOL
      call moves_carbon('CENTURY '//trim(SOIL_POOL_NAME(which)), perturb_necromass, 1.0_wp)
   end do
   !---~---
   !   Tissue HEAT and canopy CO2 have no closed-form expectation worth restating here (they
   ! are the store's own formulae); assert only that they reach the total, which is the
   ! blindness this test is about.
   !---~---
   do which = 1_ik, 2_ik
      call moves(trim(TISSUE(which))//' energy', perturb_tissue_energy, 'energy')
   end do
   call moves('canopy CO2',    perturb_cas_co2,     'carbon')
   do which = 1_ik, poly%n_pft
      call moves('recruit pool',  perturb_recruit,     'carbon')
   end do

   !--- 2. DECLARATION ARITHMETIC: a declared term is subtracted from the residual.
   led%active = .true.
   call slow_ledger_open(led, poly, cfg)
   call slow_ledger_declare(led, carbon_in = 1.0_wp, water_out = 2.0_wp, energy_in = 3.0_wp)
   call slow_ledger_mark(led, poly, cfg, SLOW_PHASE_GROW)   ! nothing moved: residual = -(in - out)
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%carbon, -1.0_wp, 1.0e-12_wp, 'declared carbon_in is subtracted')
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%water,   2.0_wp, 1.0e-12_wp, 'declared water_out is added back')
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%energy, -3.0_wp, 1.0e-12_wp, 'declared energy_in is subtracted')

   !--- An UNDECLARED, UNCHANGED step must be exactly zero, or every phase carries a bias.
   led = slow_ledger_t()
   led%active = .true.
   call slow_ledger_open(led, poly, cfg)
   call slow_ledger_mark(led, poly, cfg, SLOW_PHASE_GROW)
   call check(led%resid_sum(SLOW_PHASE_GROW)%carbon == 0.0_wp .and.                                &
              led%resid_sum(SLOW_PHASE_GROW)%water  == 0.0_wp .and.                                &
              led%resid_sum(SLOW_PHASE_GROW)%energy == 0.0_wp, 'a step that changes nothing closes at exactly zero')

   !--- 3. The canopy-air OPEN-VOLUME terms must equal the store change they explain.
   s0 = slow_site_store(poly, cfg)
   call cas_set_depth(poly%patch%cas(1), 8.0_wp, rho_air=RHO, de_open=de, dw_open=dw, dc_open=dc)
   s1 = slow_site_store(poly, cfg)
   call check_close(s1%energy - s0%energy, de,             1.0e-6_wp, 'de_open explains the canopy energy change')
   call check_close(s1%water  - s0%water,  dw,             1.0e-12_wp, 'dw_open explains the canopy vapour change')
   call check_close(s1%carbon - s0%carbon, dc * umol_2_kgc, 1.0e-14_wp, 'dc_open explains the canopy CO2 change')
   call check(de < 0.0_wp, 'shrinking the canopy DETRAINS (negative)')

   !--- 4. The fast->slow handover is nplant- and area-weighted.
   cfg%fast_biophysics_on = .true.
   poly%cohort%gpp_accum(1:2)       = 1.0e-3_wp
   poly%cohort%leaf_resp_accum(1:2) = 2.0e-4_wp
   poly%cohort%stem_resp_accum(1:2) = 1.0e-4_wp
   poly%cohort%root_resp_accum(1:2) = 1.0e-4_wp
   expect = poly%patch%area(1) * (poly%cohort%nplant(1) + poly%cohort%nplant(2)) * 6.0e-4_wp
   hand   = slow_fast_carbon_handover(poly, cfg)
   call check_close(hand, expect, 1.0e-14_wp, 'handover = area * sum(nplant * (gpp - maintenance resp))')

   !---~---
   !   5. A CULL hands its water to the ground channel, and the store sees its heat go.
   ! The assertion is the TRANSFER, not the disappearance: what leaves the cohort must
   ! arrive in the patch's shed channel, and the tissue heat store must fall by what the
   ! cohort carried.
   !---~---
   call check_cull()

   !---~---
   !   5b. A TREEFALL reports the carbon it kills, which reaches the gap's litter.
   ! With soil carbon off this carbon once left undeclared; the ledger now declares the
   ! report, so the store's fall must be exactly it.
   !---~---
   call check_treefall_carbon()



   !--- 6. Tissue THERMAL MASS tracks biomass, and is the whole of the store's tissue term.
   call check_tissue_thermal_mass()



   !---~---
   !   7. The turnover shed takes its water's HEAT with it, and declares it.
   ! The leaf's heat capacity counts its water and the wood's its whole store, so
   ! shedding either moves the tissue heat store; the declared energy must
   ! explain that move exactly.
   !---~---
   call check_shed_heat()

   !---~---
   !   8. Tissue water is booked when it is made, shed or moved: a recruit's birth water, a
   ! canopy's snap to bare, and water above a store's capacity.
   !---~---
   call check_tissue_water_ledger()

   print '(a)', 'test_slow_ledger: ALL PASSED'

contains
   !---~---
   !   Every water amount a slow step moves into or out of the tissue is declared, so the
   ! ledger's water residual is zero across each: (a) a recruitment, whose recruits are born
   ! holding their patch's top-soil potential, drawn from outside like their carbon; (b) a
   ! canopy that snaps to bare, whose leaf water goes to the ground with the leaves; (c) a
   ! store left above its capacity, whose excess the settle before the fast loop sheds to the
   ! ground.
   !---~---
   subroutine check_tissue_water_ledger()
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      type(slow_ledger_t) :: lg
      type(slow_store_t)  :: s_now
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp) :: shed0, excess, n_new, tol, theta_start(n_soil_layer_max)
      integer(ik) :: p, ph(3)
      type(soil_profile_t) :: sp

      c = build_test_config()
      tol = 1.0e-12_wp

      !---~---
      !   (a) A recruitment: water in, declared, and the recruit born in its patch's state -- here
      ! a patch only the start has seeded, at the starting soil and the reference climate's
      ! canopy air, its depth set as the driver sets it before the first step (check_birth_temp
      ! covers a patch at another temperature).
      !---~---
      call init_bare_ground(st, c, 1_ik)
      call refresh_canopy_depth(st, c, no_ledger)
      st%patch%recruit_pool(:, 1) = 0.0_wp
      st%patch%recruit_pool(2, 1) = 5.0_wp * c%demography%min_recruit_size
      lg = slow_ledger_t()
      lg%active = .true.
      call slow_ledger_open(lg, st, c)
      call restructure_stand(st, c, .true., .false., lg)
      call check(st%cohort%n == 1_ik, 'tissue water (a): one recruit is born')
      call check(lg%decl_gross(SLOW_PHASE_RECRUIT)%water > 0.0_wp,                                 &
                 'tissue water (a): the recruit''s birth water is declared as drawn in')
      !---~---
      !   Zero to round-off of the store the ledger differences: the seeded soil holds hundreds
      ! of kg/m2 against the recruit's grams.
      !---~---
      s_now = slow_site_store(st, c)
      call check(abs(lg%resid_sum(SLOW_PHASE_RECRUIT)%water)                                        &
                 <= tol * (lg%decl_gross(SLOW_PHASE_RECRUIT)%water + s_now%water),                  &
                 'tissue water (a): the recruit phase''s water residual is zero')
      leaf_curve = c%pft%hydro(2_ik)%leaf_curve
      wood_curve = c%pft%hydro(2_ik)%wood_curve
      call build_soil_profile(c%soil, sp)
      call soil_start_theta(c%soil, sp, c%init%soil_theta, theta_start)
      call check_close(tissue_water_to_psi(st%cohort%leaf_water_mass(1), leaf_curve,            &
                       st%cohort%leaf_carbon(1)), top_soil_psi(sp, theta_start(1)), 1.0e-9_wp,     &
                       'tissue water (a): a recruit''s leaves are born at the starting top-soil potential')
      call check_close(st%cohort%leaf_temp(1), st%patch%cas(1)%can_temp, 1.0e-12_wp,               &
                       'tissue water (a): a recruit is born at its patch''s canopy-air temperature')

      !--- (b) A dormant canopy sheds every leaf in one step, its water with it.
      c%fast_biophysics_on = .true.
      c%demography%demography_on      = .false.
      c%pft%pheno(:)%shed_rate_max = 1.0_wp          ! [1/day]: the whole canopy in one day
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
      call finalize_init(st)
      st%cohort%leaf_flush_tendency(1) = 0.0_wp
      st%cohort%leaf_shed_tendency(1)  = 1.0_wp
      st%cohort%gpp_accum(1)       = 0.0_wp
      st%cohort%leaf_resp_accum(1) = 0.0_wp
      st%cohort%stem_resp_accum(1) = 0.0_wp
      st%cohort%root_resp_accum(1) = 0.0_wp
      lg = slow_ledger_t()
      lg%active = .true.
      call slow_ledger_open(lg, st, c)
      call vegetation_dynamics(st, c, lg)
      call settle_tissue_water(st, c, lg)
      call check(st%cohort%leaf_carbon(1) == 0.0_wp .and. st%cohort%leaf_water_mass(1) == 0.0_wp, &
                 'tissue water (b): the canopy went bare, and its leaves hold no water')
      ph = [SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW, SLOW_PHASE_TISSUE_WATER]
      do p = 1_ik, 3_ik
         call check(abs(lg%resid_sum(ph(p))%water) <= tol * max(lg%decl_gross(ph(p))%water, 1.0_wp), &
                    'tissue water (b): the snap to bare closes the water ledger, phase by phase')
      end do
      call check(lg%decl_gross(SLOW_PHASE_ALLOCATE)%water > 0.0_wp,                                &
                 'tissue water (b): the shed leaves'' water is declared as leaving')

      !--- (c) Leaf carbon taken without its water: the settle sheds the excess, declared.
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
      call finalize_init(st)
      leaf_curve = c%pft%hydro(1_ik)%leaf_curve
      wood_curve = c%pft%hydro(1_ik)%wood_curve
      st%cohort%leaf_carbon(1) = 0.5_wp * st%cohort%leaf_carbon(1)
      call set_cohort_size_from_carbon(st%cohort, 1_ik)
      excess = st%cohort%leaf_water_mass(1) - leaf_curve%water_sat * st%cohort%leaf_carbon(1)
      call check(excess > 0.0_wp, 'tissue water (c): the fixture leaves the leaves above their capacity')
      shed0 = shed_total(st, c)
      n_new = st%cohort%nplant(1)
      lg = slow_ledger_t()
      lg%active = .true.                                   ! no window yet: the settle takes its reference
      call settle_tissue_water(st, c, lg)
      call check_close(lg%decl_gross(SLOW_PHASE_TISSUE_WATER)%water, st%patch%area(1) * n_new * excess, &
                       1.0e-12_wp, 'tissue water (c): the excess is declared as leaving')
      call check_close(shed_total(st, c) - shed0, st%patch%area(1) * n_new * excess, 1.0e-12_wp,    &
                       'tissue water (c): and it arrives in the patch shed channel')
      call check(abs(lg%resid_sum(SLOW_PHASE_TISSUE_WATER)%water) <= tol * lg%decl_gross(SLOW_PHASE_TISSUE_WATER)%water, &
                 'tissue water (c): the settle phase''s water residual is zero')
      call check(abs(lg%resid_sum(SLOW_PHASE_TISSUE_WATER)%energy) <= 1.0e-12_wp                   &
                 * lg%decl_gross(SLOW_PHASE_TISSUE_WATER)%energy, 'tissue water (c): and its energy residual')

      !---~---
      !   (d) A bare cohort re-flushes from storage: its new leaves take their water from the
      ! wood, at the wood's potential less the crown's gravity head, rather than starting at
      ! the dry end of their curve, and the plant's water is unchanged.
      !---~---
      call check_reflush_water()
   end subroutine check_tissue_water_ledger
   subroutine check_reflush_water()
      type(ed_config_t)   :: c
      type(polygon_t)     :: st
      type(slow_ledger_t) :: lg
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp) :: w_plant0, w_plant1, psi_leaf, psi_wood
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography%demography_on      = .false.
      c%pft%fineroot_turnover_rate(:) = 0.0_wp      ! no fine-root shed: only the fill moves water
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
      call finalize_init(st)
      leaf_curve = c%pft%hydro(1_ik)%leaf_curve
      wood_curve = c%pft%hydro(1_ik)%wood_curve
      st%cohort%leaf_carbon(1)     = 0.0_wp           ! bare, with its storage, at a drier stem
      st%cohort%leaf_water_mass(1) = 0.0_wp
      call set_cohort_size_from_carbon(st%cohort, 1_ik)
      st%cohort%wood_water_mass(1) = 0.9_wp * st%cohort%wood_water_mass(1)
      st%cohort%leaf_flush_tendency(1) = 1.0_wp
      st%cohort%leaf_shed_tendency(1)  = 0.0_wp
      st%cohort%gpp_accum(1)       = 0.0_wp
      st%cohort%leaf_resp_accum(1) = 0.0_wp
      st%cohort%stem_resp_accum(1) = 0.0_wp
      st%cohort%root_resp_accum(1) = 0.0_wp
      w_plant0 = st%cohort%leaf_water_mass(1) + st%cohort%wood_water_mass(1)
      lg = slow_ledger_t()
      lg%active = .true.
      call slow_ledger_open(lg, st, c)
      call vegetation_dynamics(st, c, lg)
      call check(st%cohort%leaf_carbon(1) > 0.0_wp, 'reflush: the bare cohort grew leaves')
      w_plant1 = st%cohort%leaf_water_mass(1) + st%cohort%wood_water_mass(1)
      call check_close(w_plant1, w_plant0, 1.0e-14_wp, 'reflush: the plant''s water is unchanged')
      psi_leaf = tissue_water_to_psi(st%cohort%leaf_water_mass(1), leaf_curve, st%cohort%leaf_carbon(1))
      !--- The wood's potential before it gave the leaves their water.
      psi_wood = tissue_water_to_psi(st%cohort%wood_water_mass(1) + st%cohort%leaf_water_mass(1), wood_curve,  &
                                     st%cohort%sapwood_carbon(1) + st%cohort%fineroot_carbon(1))
      call check(psi_wood < -0.1_wp, 'reflush: the fixture''s stem is drier than the birth potential')
      call check(abs(psi_leaf - (psi_wood - grav_head * st%cohort%height(1))) <= 1.0e-9_wp,        &
                 'reflush: the new leaves start at the wood''s potential less the crown head')
      call check(abs(lg%resid_sum(SLOW_PHASE_GROW)%water) <= 1.0e-12_wp,                           &
                 'reflush: the fill moves water inside the plant; the ledger sees no residual')
   end subroutine check_reflush_water
   subroutine check_shed_heat()
      type(slow_step_data_t) :: step
      type(slow_ledger_t)     :: lg
      integer(ik) :: nc, i, ip
      real(wp)    :: fast_in
      nc = poly%cohort%n
      allocate(step%d_leaf(nc), step%d_fineroot(nc), step%d_wood(nc), step%d_nonstructural(nc))
      step%d_wood = 0.0_wp
      step%d_nonstructural = 0.0_wp
      step%d_leaf     = -0.25_wp * poly%cohort%leaf_carbon(1:nc)       ! a quarter of the leaves shed
      step%d_fineroot = -0.10_wp * poly%cohort%fineroot_carbon(1:nc)   ! a tenth of the fine roots
      poly%cohort%leaf_water_mass(1:nc) = 0.7_wp
      poly%cohort%wood_water_mass(1:nc) = 2.3_wp
      do i = 1_ik, nc
         call hold_cohort_tissue_temps(poly%cohort, i)
      end do
      lg%active = .true.
      call slow_ledger_open(lg, poly, cfg)
      call shed_turnover_water(poly, cfg, step, lg)
      call slow_ledger_mark(lg, poly, cfg, SLOW_PHASE_ALLOCATE)
      call check(lg%decl_gross(SLOW_PHASE_ALLOCATE)%energy > 0.0_wp, 'shed heat: the shed declares energy')
      call check(abs(lg%resid_sum(SLOW_PHASE_ALLOCATE)%energy)                                      &
                 <= 1.0e-12_wp * lg%decl_gross(SLOW_PHASE_ALLOCATE)%energy,                            &
                 'shed heat: the enthalpy handed to the ground is what the tissue store lost')
      call check(abs(lg%resid_sum(SLOW_PHASE_ALLOCATE)%water)                                       &
                 <= 1.0e-12_wp * lg%decl_gross(SLOW_PHASE_ALLOCATE)%water,                             &
                 'shed heat: the declared water is the tissue water store''s change')
      !--- What the fast tier takes in over the day it reads the rate, summed over the patches.
      fast_in = 0.0_wp
      do ip = 1_ik, poly%patch%n
         fast_in = fast_in + patch_polygon_area(poly, ip) * poly%patch%shed_water_enth(ip) * cfg%dt_slow
      end do
      call check_close(fast_in, lg%decl_gross(SLOW_PHASE_ALLOCATE)%energy, 1.0e-12_wp,                  &
                       'shed heat: the fast tier takes in the enthalpy the shed handed off')
   end subroutine check_shed_heat
   pure function shed_total(st, c) result(w)
      type(polygon_t),     intent(in) :: st
      type(ed_config_t),   intent(in) :: c
      real(wp)    :: w
      integer(ik) :: ip
      w = 0.0_wp
      do ip = 1_ik, st%patch%n
         w = w + st%patch%area(ip) * st%patch%shed_water_rate(ip) * c%dt_slow
      end do
   end function shed_total
   !---~---
   !   A cohort's heat capacity is a function of its biomass, so changing biomass alone moves
   ! the tissue store with no flux. That is the exchange the growth commit declares; these
   ! assert the two properties the declaration rests on -- that the term is real (biomass
   ! moves it) and that slow_tissue_energy is exactly the store's tissue slice (so declaring
   ! the one cancels the other).
   !---~---
   subroutine check_tissue_thermal_mass()
      type(slow_store_t) :: a, b
      real(wp) :: e0, e1, de_store, de_tissue
      e0 = slow_tissue_energy(poly)
      a  = slow_site_store(poly, cfg)
      !---~---
      !   Grow one cohort's wood at its temperature, as the growth commit does: ONLY the
      ! thermal mass changes.
      !---~---
      poly%cohort%wood_carbon(1) = poly%cohort%wood_carbon(1) * 1.5_wp
      call hold_cohort_tissue_temps(poly%cohort, 1_ik)
      e1 = slow_tissue_energy(poly)
      b  = slow_site_store(poly, cfg)
      de_tissue = e1 - e0
      de_store  = b%energy - a%energy
      call check(de_tissue > 0.0_wp, 'thermal mass: growing biomass raises the tissue heat store')
      call check_close(de_store, de_tissue, 1.0e-9_wp,                                             &
                       'thermal mass: slow_tissue_energy is exactly the store''s tissue slice')
      !---~---
      !   and it is the BIOMASS doing it, not the heat: take heat from the leaves and the store
      ! must fall, so a declaration written against biomass alone would not cancel that.
      !---~---
      poly%cohort%leaf_energy(1) = poly%cohort%leaf_energy(1) - 1.0e4_wp
      call check(slow_tissue_energy(poly) < e1, 'thermal mass: heat moves it too (so it is not frozen in)')
   end subroutine check_tissue_thermal_mass
   subroutine check_cull()
      real(wp) :: w_before, shed_before, w_shed, expect_w, heat_before, heat_after
      real(wp) :: lit_before, kill_c
      type(slow_store_t) :: before, after
      integer(ik) :: n0
      !--- Drive cohort 2 below the density floor so the cull takes it, and nothing else.
      poly%cohort%nplant(2)   = 0.5_wp * cfg%demography%negligible_nplant
      poly%cohort%leaf_temp(2) = 291.0_wp
      poly%cohort%wood_temp(2) = 289.0_wp
      n0          = poly%cohort%n
      shed_before = poly%patch%shed_water_rate(1)
      heat_before = slow_tissue_energy(poly)
      w_before    = poly%cohort%nplant(2) * (poly%cohort%leaf_water_mass(2) + poly%cohort%wood_water_mass(2) &
                                         + poly%cohort%leaf_surf_water(2) + poly%cohort%wood_surf_water(2))
      expect_w    = poly%patch%area(1) * w_before
      before      = slow_site_store(poly, cfg)
      lit_before  = litter_carbon_total(poly, poly%patch%litter_in)

      call terminate_cohorts(poly, cfg, w_shed, kills)
      call kill_litter(poly, kills, kill_c)
      after = slow_site_store(poly, cfg)

      !---~---
      !   Its carbon: what the stand lost is what the cull reports, for the ledger to
      ! declare, and it reaches the patch's litter.
      !---~---
      call check(kill_c > 0.0_wp, 'cull: the culled cohort carried carbon')
      call check_close(before%carbon - after%carbon, kill_c, 1.0e-12_wp,                          &
                       'cull: the carbon it reports is the carbon the stand lost')
      call check_close(litter_carbon_total(poly, poly%patch%litter_in) - lit_before, kill_c,      &
                       1.0e-12_wp, 'cull: its carbon reaches the patch litter')

      call check(poly%cohort%n == n0 - 1_ik, 'cull: the sub-floor cohort is removed')
      call check_close(w_shed, expect_w, 1.0e-12_wp, 'cull: reports the tissue + film water it carried')
      call check_close(poly%patch%shed_water_rate(1) - shed_before, w_before / cfg%dt_slow,       &
                       1.0e-12_wp, 'cull: the water reaches the patch shed channel, not the void')
      !---~---
      !   The cull's tissue HEAT is not reported by the operator: it is one part of the
      ! phase's thermal-mass change, which the driver brackets with slow_tissue_energy. Assert
      ! that the store SEES the loss, which is what makes that bracket able to declare it.
      !---~---
      heat_after = slow_tissue_energy(poly)
      call check(heat_after < heat_before, 'cull: the tissue heat store falls when a cohort goes')
   end subroutine check_cull
   !--- A treefall gap that kills every cohort's share, on a copy of the stand.
   subroutine check_treefall_carbon()
      type(polygon_t)      :: st
      logical, allocatable :: surv(:)
      real(wp)             :: lit_before, w_shed, kill_c
      type(slow_store_t)   :: before, after
      st = poly
      allocate(surv(st%cohort%n))
      surv = .false.
      before     = slow_site_store(st, cfg)
      lit_before = litter_carbon_total(st, st%patch%litter_in)
      call apply_patch_disturbance(st, cfg, 0.1_wp, surv, w_shed, kills)
      call kill_litter(st, kills, kill_c)
      after = slow_site_store(st, cfg)
      call check(kill_c > 0.0_wp, 'treefall: the gap killed carbon')
      call check_close(before%carbon - after%carbon, kill_c, 1.0e-12_wp,                          &
                       'treefall: the carbon it reports is the carbon the stand lost')
      call check_close(litter_carbon_total(st, st%patch%litter_in) - lit_before, kill_c,          &
                       1.0e-12_wp, 'treefall: its carbon reaches the litter')
   end subroutine check_treefall_carbon
   !--- Assert a perturbation moves the named currency by `amount` (and moves it at all).
   subroutine moves_water(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(poly, cfg)
      call perturb()
      b = slow_site_store(poly, cfg)
      call check_close(b%water - a%water, amount, 1.0e-9_wp * max(1.0_wp, abs(amount)),            &
                       'store: '//name//' reaches the water total')
   end subroutine moves_water
   subroutine moves_energy(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(poly, cfg)
      call perturb()
      b = slow_site_store(poly, cfg)
      call check_close(b%energy - a%energy, amount, 1.0e-6_wp * max(1.0_wp, abs(amount)),          &
                       'store: '//name//' reaches the energy total')
   end subroutine moves_energy
   subroutine moves_carbon(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(poly, cfg)
      call perturb()
      b = slow_site_store(poly, cfg)
      call check_close(b%carbon - a%carbon, amount, 1.0e-9_wp * max(1.0_wp, abs(amount)),          &
                       'store: '//name//' reaches the carbon total')
   end subroutine moves_carbon
   subroutine moves(name, perturb, currency)
      character(len=*), intent(in) :: name, currency
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      real(wp) :: d
      a = slow_site_store(poly, cfg)
      call perturb()
      b = slow_site_store(poly, cfg)
      select case (currency)
      case ('carbon')
         d = b%carbon - a%carbon
      case ('water')
         d = b%water  - a%water
      case default
         d = b%energy - a%energy
      end select
      call check(abs(d) > 0.0_wp, 'store: '//name//' reaches the '//currency//' total')
   end subroutine moves
   subroutine perturb_theta()
   poly%patch%soil_water(1)%theta(1) = poly%patch%soil_water(1)%theta(1) + 0.01_wp
   end subroutine perturb_theta
   subroutine perturb_pond()
   poly%patch%soil_water(1)%w_surface = poly%patch%soil_water(1)%w_surface + 1.0_wp
   end subroutine perturb_pond
   subroutine perturb_swe()
   poly%patch%snow(1)%swe(1) = poly%patch%snow(1)%swe(1) + 1.0_wp
   end subroutine perturb_swe
   subroutine perturb_shv()
   poly%patch%cas(1)%can_shv = poly%patch%cas(1)%can_shv + 0.001_wp
   end subroutine perturb_shv
   !--- The members of a group, one at a time: `which` picks the tissue, pool or PFT.
   subroutine perturb_tissue_water()
      select case (which)
      case (1_ik)
         poly%cohort%leaf_water_mass(1) = poly%cohort%leaf_water_mass(1) + 0.01_wp
      case default
         poly%cohort%wood_water_mass(1) = poly%cohort%wood_water_mass(1) + 0.01_wp
      end select
   end subroutine perturb_tissue_water
   subroutine perturb_film()
      select case (which)
      case (1_ik)
         poly%cohort%leaf_surf_water(1) = poly%cohort%leaf_surf_water(1) + 0.01_wp
      case default
         poly%cohort%wood_surf_water(1) = poly%cohort%wood_surf_water(1) + 0.01_wp
      end select
   end subroutine perturb_film
   subroutine perturb_soil_energy()
      poly%patch%soil_energy(1)%soil_energy(1) = poly%patch%soil_energy(1)%soil_energy(1) + 1.0e6_wp
   end subroutine perturb_soil_energy
   subroutine perturb_pond_energy()
      poly%patch%soil_water(1)%w_surface_energy = poly%patch%soil_water(1)%w_surface_energy + 1.0e5_wp
   end subroutine perturb_pond_energy
   subroutine perturb_snow_energy()
      poly%patch%snow(1)%snow_energy(1) = poly%patch%snow(1)%snow_energy(1) + 1.0e5_wp
   end subroutine perturb_snow_energy
   subroutine perturb_cas_enth()
   poly%patch%cas(1)%can_enthalpy = poly%patch%cas(1)%can_enthalpy + 1.0e4_wp
   end subroutine perturb_cas_enth
   subroutine perturb_cas_co2()
   poly%patch%cas(1)%can_co2 = poly%patch%cas(1)%can_co2 + 50.0_wp
   end subroutine perturb_cas_co2
   subroutine perturb_live_c()
      select case (which)
      case (1_ik)
         poly%cohort%leaf_carbon(1)          = poly%cohort%leaf_carbon(1)          + 0.5_wp
      case (2_ik)
         poly%cohort%fineroot_carbon(1)      = poly%cohort%fineroot_carbon(1)      + 0.5_wp
      case (3_ik)
         poly%cohort%wood_carbon(1)          = poly%cohort%wood_carbon(1)          + 0.5_wp
      case default
         poly%cohort%nonstructural_carbon(1) = poly%cohort%nonstructural_carbon(1) + 0.5_wp
      end select
   end subroutine perturb_live_c
   subroutine perturb_necromass()
      real(wp) :: x(N_SOIL_POOL)
      x = poly%patch%soil_carbon(1)%carbon
      x(which) = x(which) + 1.0_wp
      poly%patch%soil_carbon(1)%carbon = x
   end subroutine perturb_necromass
   subroutine perturb_tissue_energy()
      select case (which)
      case (1_ik)
         poly%cohort%leaf_energy(1) = poly%cohort%leaf_energy(1) + 1.0e3_wp
      case default
         poly%cohort%wood_energy(1) = poly%cohort%wood_energy(1) + 1.0e3_wp
      end select
   end subroutine perturb_tissue_energy
   subroutine perturb_recruit()
   poly%patch%recruit_pool(which,1) = poly%patch%recruit_pool(which,1) + 0.5_wp
   end subroutine perturb_recruit
end program test_slow_ledger
