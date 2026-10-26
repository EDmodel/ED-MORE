! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! The commit: the applier advances only the carbon and derives the geometry from it and   !
! holds each tissue's temperature through its change of mass, the driver's tendency       !
! follows the applier's rule, a decline aimed at the floor ends on it, and the cull        !
! removes the cohort there.                                                                !
!==========================================================================================!
module test_slow_update_commit
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : tissue_c_wood_dry
   use ed_allometry_lib,        only : carbon_to_structure
   use ed_therm_lib,            only : tissue_dry_heat_capacity
   use ed_config,               only : ed_config_t, growth_window_steps
   use ed_cohort_state_vars,    only : cohort_deriv_alloc, cohort_bad_field, cohort_field_name, cohort_tissue_temps
   use ed_polygon_state_vars,   only : polygon_t
   use ed_slow_update,          only : update_cohort_states, terminate_cohorts,                          &
                                       update_cohort_derivatives
   use ed_slow_state,           only : kill_report_t, slow_step_data_t
   use ed_init_bare_ground,     only : init_bare_ground, add_cohort, finalize_init
   use ed_test_support,         only : banner, build_test_config, check, check_close
   use, intrinsic :: ieee_arithmetic,         only : ieee_value, ieee_quiet_nan
   implicit none
   private

   public :: run_commit
   type(kill_report_t)  :: kills   !< what the cull and the treefall killed
   real(wp)            :: shed_water   !< an operator's shed water, unread here

   type(ed_config_t)   :: cfg
   type(polygon_t)     :: poly

contains

   subroutine run_commit()
      call banner('the commit: the applier, the held temperatures, the tendency rule, the floor and the cull')
      cfg = build_test_config()

      !=== The applier advances only the carbon and derives the geometry from it: a cohort         !
      !    that sheds every leaf ends with leaf area exactly 0, never a round-off negative, and its     !
      !    dbh and AGB are the allometry of its committed wood carbon.  ===============================!
      block
         real(wp) :: dt, dbh_c, h_c, ba_c, agb_c, la_c
         dt = 1.0_wp / 365.0_wp
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
         call finalize_init(poly)
         call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
         !----- Shed a hair more leaf carbon than there is, as a rounded loss rate can. ------------!
         poly%deriv%d_leaf_carbon_dt(1)          = -poly%cohort%leaf_carbon(1) / dt * (1.0_wp + 1.0e-15_wp)
         poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
         poly%deriv%d_wood_carbon_dt(1)          = 0.01_wp
         poly%deriv%d_nonstructural_carbon_dt(1) = 0.0_wp
         poly%deriv%dln_nplant_dt(1)             = 0.0_wp
         call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%demography%negligible_nplant)
         call check(poly%cohort%leaf_carbon(1) == 0.0_wp, 'applier: a full shed leaves no leaf carbon')
         call check(poly%cohort%leaf_area(1) == 0.0_wp, 'applier: and exactly no leaf area')
         call carbon_to_structure(poly%cohort%pft_table%allom(poly%cohort%pft(1)), poly%cohort%wood_carbon(1),        &
                                  poly%cohort%leaf_carbon(1), cfg%pft%wood_density(poly%cohort%pft(1)),     &
                                  poly%cohort%sla(1), dbh_c, h_c, ba_c, agb_c, la_c)
         call check(poly%cohort%dbh(1) == dbh_c .and. poly%cohort%height(1) == h_c                     &
                    .and. poly%cohort%agb(1) == agb_c, 'applier: geometry is the allometry of the carbon')

         !----- An overdraw beyond round-off is an operator bug: the applier does not floor it (a   !
         !      floor would create the carbon it took), and the bounds guard names the pool. --------!
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
         call finalize_init(poly)
         call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
         poly%deriv%d_leaf_carbon_dt(1)          = 0.0_wp
         poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
         poly%deriv%d_wood_carbon_dt(1)          = 0.0_wp
         poly%deriv%d_nonstructural_carbon_dt(1) = -poly%cohort%nonstructural_carbon(1) / dt * 1.01_wp
         poly%deriv%dln_nplant_dt(1)             = 0.0_wp
         call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'applier: the cohort starts in bounds')
         call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%demography%negligible_nplant)
         call check(poly%cohort%nonstructural_carbon(1) < 0.0_wp, 'applier: a 1 % overdraw is not floored')
         call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'nonstructural_carbon', &
                    'applier: the bounds guard names the overdrawn pool')
      end block

      !=== The applier holds each tissue's temperature through its change of mass                   !
      !    (update_cohort_states, through hold_cohort_tissue_temps; docs/science/                      !
      !    vegetation_energy_dynamics.md eqs. (9)-(10)): the tissue made or lost came or went at the    !
      !    tissue's temperature. Wood that grows keeps its temperature, its energy rising by the new    !
      !    dry wood's heat at that temperature; a crown that sheds its last leaf, its water gone,       !
      !    holds no energy and takes its wood's temperature. =========================================!
      block
         real(wp) :: dt, tl, tw, t_wood, e_wood, c_wood
         dt = 1.0_wp / 365.0_wp
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
         call finalize_init(poly)
         call cohort_tissue_temps(poly%cohort, 1_ik, tl, tw)
         call check(tl == poly%cohort%leaf_temp(1) .and. tw == poly%cohort%wood_temp(1),               &
                    'held temperatures: a cohort is born with the temperatures its energies give')
         t_wood = poly%cohort%wood_temp(1)
         e_wood = poly%cohort%wood_energy(1)
         c_wood = poly%cohort%wood_carbon(1)
         call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
         poly%deriv%d_leaf_carbon_dt(1)          = 0.0_wp
         poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
         poly%deriv%d_wood_carbon_dt(1)          = 0.2_wp * c_wood / dt
         poly%deriv%d_nonstructural_carbon_dt(1) = 0.0_wp
         poly%deriv%dln_nplant_dt(1)             = 0.0_wp
         call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%demography%negligible_nplant)
         call check(abs(poly%cohort%wood_temp(1) - t_wood) < 1.0e-10_wp,                              &
                    'held temperatures: wood that grows at its temperature keeps it')
         call check(abs(poly%cohort%wood_energy(1) - e_wood                                            &
                        - tissue_dry_heat_capacity(1.0_wp, poly%cohort%wood_carbon(1) - c_wood, tissue_c_wood_dry) &
                          * t_wood) <= 1.0e-9_wp * e_wood,                                             &
                    'held temperatures: and its energy rises by the new dry wood''s heat at that temperature')
         !----- The last leaf shed, as a rounded full shed is (it lands on 0), its water gone. -----!
         poly%cohort%leaf_water_mass(1)          = 0.0_wp
         poly%deriv%d_leaf_carbon_dt(1)          = -poly%cohort%leaf_carbon(1) / dt * (1.0_wp + 1.0e-15_wp)
         poly%deriv%d_wood_carbon_dt(1)          = 0.0_wp
         call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%demography%negligible_nplant)
         call check(poly%cohort%leaf_energy(1) == 0.0_wp, 'held temperatures: a leafless crown holds no energy')
         call check(poly%cohort%leaf_temp(1) == poly%cohort%wood_temp(1),                              &
                    'held temperatures: and takes its wood''s temperature')
      end block

      !=== The driver's tendency follows the applier's rule: a step that empties a pool lands it on  !
      !    0, and an overdraw reaches the applier as an overdraw, which the bounds guard names. A     !
      !    floor on the driver's tentative pool would create the carbon before the guard could see it. =!
      block
         type(slow_step_data_t) :: step
         real(wp) :: dt
         dt = 1.0_wp / 365.0_wp
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
         call finalize_init(poly)
         allocate(step%d_leaf(1), step%d_fineroot(1), step%d_wood(1), step%d_nonstructural(1))
         step%d_leaf(1)          = -poly%cohort%leaf_carbon(1) * (1.0_wp + 1.0e-15_wp)   ! a rounded full shed
         step%d_fineroot(1)      = 0.0_wp
         step%d_wood(1)          = 0.0_wp
         step%d_nonstructural(1) = -poly%cohort%nonstructural_carbon(1) * 1.01_wp        ! a 1 % overdraw
         call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
         step%mortality = [0.0_wp]
         call update_cohort_derivatives(poly, step, dt, growth_window_steps(cfg), 1_ik)
         call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%demography%negligible_nplant)
         call check(poly%cohort%leaf_carbon(1) == 0.0_wp, 'driver + applier: a rounded full shed lands on 0')
         call check(poly%cohort%nonstructural_carbon(1) < 0.0_wp,                                     &
                    'driver + applier: a 1 % overdraw reaches the applier, not floored on the way')
         call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'nonstructural_carbon', &
                    'driver + applier: the bounds guard names the overdrawn pool')
      end block

      !=== A density that declines to the floor ends ON it, and the cull removes a cohort there. ====!
      !    nplant * exp(log(floor/nplant)) is the floor only to an epsilon either side; a cohort an   !
      !    epsilon above the floor would never be culled, and would decline no further. ===========!
      block
         real(wp)    :: dt, n0
         integer(ik) :: k, n_on, n_cut
         dt = 1.0_wp / 365.0_wp
         n_on = 0_ik
         do k = 1_ik, 200_ik
            !----- 200 starting densities across six decades, each aimed exactly at the floor. -------!
            n0 = cfg%demography%negligible_nplant * 10.0_wp ** (6.0_wp * real(k, wp) / 200.0_wp)
            call init_bare_ground(poly, cfg, 1_ik)
            call add_cohort(poly, cfg, 1_ik, 2_ik, n0, 3.0_wp)
            call finalize_init(poly)
            call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
            poly%deriv%d_leaf_carbon_dt(1)          = 0.0_wp
            poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
            poly%deriv%d_wood_carbon_dt(1)          = 0.0_wp
            poly%deriv%d_nonstructural_carbon_dt(1) = 0.0_wp
            poly%deriv%dln_nplant_dt(1)             = -(log(n0 / cfg%demography%negligible_nplant) / dt)
            call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%demography%negligible_nplant)
            if (poly%cohort%nplant(1) == cfg%demography%negligible_nplant) n_on = n_on + 1_ik
         end do
         call check(n_on == 200_ik, 'applier: a decline aimed at the floor ends exactly on it')

         !----- Two cohorts of a tracked size: one at the floor (culled), one a hair above (kept). ---!
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 1.0_wp, 60.0_wp)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 1.0_wp, 61.0_wp)
         call finalize_init(poly)
         poly%cohort%nplant(1) = cfg%demography%negligible_nplant
         poly%cohort%nplant(2) = nearest(cfg%demography%negligible_nplant, 1.0_wp)
         call check(all(poly%cohort%nplant(1:2) * poly%cohort%agb(1:2) >= cfg%demography%min_cohort_agb),        &
                    'cull fixture: both cohorts clear the AGB floor, so only the density floor acts')
         n_cut = poly%cohort%n
         call terminate_cohorts(poly, cfg, shed_water, kills)
         call check(poly%cohort%n == n_cut - 1_ik, 'cull: exactly one cohort removed')
         call check(poly%cohort%nplant(1) > cfg%demography%negligible_nplant, 'cull: the cohort AT the floor is the one removed')
      end block
   end subroutine run_commit


end module test_slow_update_commit

!==========================================================================================!
! Cohort fusion and fission: what each conserves, the tolerance, the lockstep sort, carbon !
! and leaf area through a fusion, geometry's one writer, the bounds guard and the growth   !
! window.                                                                                  !
!==========================================================================================!
module test_slow_update_cohorts
   use ed_kinds,             only : wp, ik
   use ed_constants,         only : pio4
   use ed_config,            only : ed_config_t
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ed_cohort_state_vars,            only : set_cohort_size, cohort_tissue_heat_capacity, set_cohort_size_from_carbon, &
                                               cohort_bad_field, cohort_field_name, cohort_deriv_block, cohort_deriv_alloc, &
                                               init_cohort, hold_cohort_tissue_temps
   use ed_polygon_state_vars,           only : polygon_t
   use ed_slow_update,   only : fill_cohort_deriv, fuse_2_cohorts, fuse_cohorts, split_cohorts,      &
                                max_cohort_count
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_output_reduce, only : total_nplant, total_agb
   use ed_diag_vars,        only : alloc_cohort_diag, cohort_diag_value, CO_LEAF_TEMP, CO_WOOD_TEMP
   use ed_allometry_lib,   only : storage_carbon_target
   use ed_test_support, only : banner, build_test_config, check, check_close
   implicit none
   private

   public :: run_cohort_operators

   type(ed_config_t) :: cfg
   type(polygon_t)  :: poly
   real(wp)            :: agb_tot, n0, agb0, dbh_avg
   real(wp)            :: wr, wd, leafmass_exp, woodmass_exp, ltemp_exp, wtemp_exp
   real(wp)            :: cl1, cw1, cl2, cw2, clf, cwf, e_leaf0, e_wood0
   real(wp) :: np1, np2, xdiag(4), film0, e_sum0
   integer(ik) :: ndiag
   integer(ik)         :: j, pf

contains

   subroutine run_cohort_operators()
      call banner('cohort fusion/fission conservation')
      cfg = build_test_config()

      !=== 1. fuse_2_cohorts conserves N and AGB; DBH is re-derived, not averaged. ============!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
      call finalize_init(poly)                          ! sorted: index 1 = taller (dbh 12)
      agb_tot = poly%cohort%nplant(1)*poly%cohort%agb(1) + poly%cohort%nplant(2)*poly%cohort%agb(2)
      dbh_avg = 0.5_wp*(poly%cohort%dbh(1) + poly%cohort%dbh(2))
      !----- Seed distinct fast state + predict the merge. The tissue TEMPERATURES are intensive per   !
      !      unit HEAT CAPACITY, so they fuse on the capacity that carries them -- what that conserves  !
      !      is the tissue ENERGY, which is asserted directly below as well as through the formula.     !
      !      A leaf-area weighting is roughly right for leaves and wrong for wood (whose capacity       !
      !      follows wood carbon and the sapwood ring, not leaf area).                                  !
      !      leaf_water_mass/wood_water_mass are EXTENSIVE like AGB (nplant-weighted, so the fused      !
      !      site-TOTAL water is conserved -- leaf-area-weighting them would not conserve it). ---------!
      poly%cohort%leaf_water_mass(1) = 2.0_wp ; poly%cohort%leaf_water_mass(2) = 5.0_wp
      poly%cohort%wood_water_mass(1) = 3.0_wp ; poly%cohort%wood_water_mass(2) = 7.0_wp
      poly%cohort%leaf_temp(1) = 300.0_wp ; poly%cohort%leaf_temp(2) = 305.0_wp
      poly%cohort%wood_temp(1) = 298.0_wp ; poly%cohort%wood_temp(2) = 306.0_wp
      !----- And the canopy films [kg/plant]: nplant-weighted like the tissue water, and part of the  !
      !      tissues' heat capacity, so set before the energies. ------------------------------------!
      poly%cohort%leaf_surf_water(1) = 0.02_wp
      poly%cohort%leaf_surf_water(2) = 0.06_wp
      poly%cohort%wood_surf_water(1) = 0.01_wp
      poly%cohort%wood_surf_water(2) = 0.00_wp
      !----- ...through their energies, the state: each tissue holds what it does at those temperatures. !
      call hold_cohort_tissue_temps(poly%cohort, 1_ik)
      call hold_cohort_tissue_temps(poly%cohort, 2_ik)
      e_sum0 = poly%cohort%nplant(1) * (poly%cohort%leaf_energy(1) + poly%cohort%wood_energy(1))         &
             + poly%cohort%nplant(2) * (poly%cohort%leaf_energy(2) + poly%cohort%wood_energy(2))
      !----- The four per-plant flux accumulators are EXTENSIVE too. Assert them, so the declared     !
      !      policy in fuse_cohort_fast_state has a witness for every kind it names.                  !
      poly%cohort%gpp_accum(1)       = 1.0_wp ; poly%cohort%gpp_accum(2)       = 4.0_wp
      poly%cohort%leaf_resp_accum(1) = 0.2_wp ; poly%cohort%leaf_resp_accum(2) = 0.9_wp
      poly%cohort%stem_resp_accum(1) = 0.1_wp ; poly%cohort%stem_resp_accum(2) = 0.6_wp
      poly%cohort%root_resp_accum(1) = 0.3_wp ; poly%cohort%root_resp_accum(2) = 0.8_wp
      wr = poly%cohort%nplant(1) * poly%cohort%leaf_area(1)
      wd = poly%cohort%nplant(2) * poly%cohort%leaf_area(2)
      np1 = poly%cohort%nplant(1) ; np2 = poly%cohort%nplant(2)   ! captured BEFORE the fuse
      leafmass_exp = (poly%cohort%nplant(1)*2.0_wp + poly%cohort%nplant(2)*5.0_wp)                     &
                   / (poly%cohort%nplant(1) + poly%cohort%nplant(2))
      woodmass_exp = (poly%cohort%nplant(1)*3.0_wp + poly%cohort%nplant(2)*7.0_wp)                     &
                   / (poly%cohort%nplant(1) + poly%cohort%nplant(2))
      call cohort_tissue_heat_capacity(poly%cohort, 1_ik, cl1, cw1)
      call cohort_tissue_heat_capacity(poly%cohort, 2_ik, cl2, cw2)
      ltemp_exp = (cl1*300.0_wp + cl2*305.0_wp) / (cl1 + cl2)
      wtemp_exp = (cw1*298.0_wp + cw2*306.0_wp) / (cw1 + cw2)
      e_leaf0   = cl1*300.0_wp + cl2*305.0_wp
      e_wood0   = cw1*298.0_wp + cw2*306.0_wp
      !----- The cohort DIAGNOSTIC block, holding one sample of the same tissue temperatures, must fuse  !
      !      them on the same weights as the prognostic state, so the two agree after the merge. --------!
      call alloc_cohort_diag(poly%cohort%diag, poly%cohort%cap, .true.)
      poly%cohort%diag%n = 2_ik
      poly%cohort%diag%v(CO_LEAF_TEMP, 1:2) = [300.0_wp, 305.0_wp]
      poly%cohort%diag%v(CO_WOOD_TEMP, 1:2) = [298.0_wp, 306.0_wp]
      poly%cohort%diag%w(1:2) = 1.0_wp
      call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%demography%conservation_tol)
      call check_close(poly%cohort%nplant(1), 0.8_wp, 1.0e-12_wp, 'fused nplant must be summed')
      call check_close(poly%cohort%nplant(1)*poly%cohort%agb(1), agb_tot, 1.0e-12_wp,                  &
                       'fused total AGB not conserved')
      call check_close(poly%cohort%basal_area(1), pio4*poly%cohort%dbh(1)**2, 1.0e-12_wp,              &
                       'basal area inconsistent with re-derived DBH')
      call check(abs(poly%cohort%dbh(1) - dbh_avg) > 1.0e-6_wp, 'DBH must NOT be a plain average')
      call check_close(poly%cohort%leaf_water_mass(1), leafmass_exp, 1.0e-12_wp,                       &
                       'leaf_water_mass not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%wood_water_mass(1), woodmass_exp, 1.0e-12_wp,                       &
                       'wood_water_mass not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%leaf_surf_water(1), (np1*0.02_wp + np2*0.06_wp) / (np1 + np2), 1.0e-12_wp, &
                       'leaf_surf_water not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%wood_surf_water(1), np1*0.01_wp / (np1 + np2), 1.0e-12_wp,          &
                       'wood_surf_water not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%leaf_temp(1), ltemp_exp, 1.0e-9_wp,                                 &
                       'leaf_temp not heat-capacity-weighted on cohort fusion')
      call check_close(poly%cohort%wood_temp(1), wtemp_exp, 1.0e-9_wp,                                 &
                       'wood_temp not heat-capacity-weighted on cohort fusion')
      call cohort_diag_value(poly%cohort%diag, CO_LEAF_TEMP, xdiag, ndiag, 0.0_wp)
      call check_close(xdiag(1), poly%cohort%leaf_temp(1), 1.0e-9_wp,                                  &
                       'diagnostic leaf temperature fuses as the prognostic one (leaf heat capacity)')
      call cohort_diag_value(poly%cohort%diag, CO_WOOD_TEMP, xdiag, ndiag, 0.0_wp)
      call check_close(xdiag(1), poly%cohort%wood_temp(1), 1.0e-9_wp,                                  &
                       'diagnostic wood temperature fuses as the prognostic one (wood heat capacity)')
      !----- The formulae above are the MEANS to this end: the merged tissue must hold the energy the  !
      !      two cohorts held. Asserted separately so a future change of weighting is judged on what   !
      !      it conserves, not on whether it matches today's expression. Both capacities are additive  !
      !      through the merge (the carbon pools and the tissue waters are all nplant-weighted), so    !
      !      both energies close exactly. -----------------------------------------------------------!
      call cohort_tissue_heat_capacity(poly%cohort, 1_ik, clf, cwf)
      call check_close(clf*poly%cohort%leaf_temp(1), e_leaf0, 1.0e-12_wp,                              &
                       'cohort fusion must conserve LEAF tissue energy exactly')
      call check_close(cwf*poly%cohort%wood_temp(1), e_wood0, 1.0e-12_wp,                              &
                       'cohort fusion must conserve WOOD tissue energy exactly')
      call check_close(poly%cohort%nplant(1) * (poly%cohort%leaf_energy(1) + poly%cohort%wood_energy(1)),   &
                       e_sum0, 1.0e-12_wp, 'cohort fusion sums the tissue energies, the state')
      block
         real(wp) :: nt
         nt = np1 + np2
         call check_close(poly%cohort%gpp_accum(1),       (np1*1.0_wp + np2*4.0_wp)/nt, 1.0e-12_wp,      &
                          'gpp_accum not nplant-weighted (extensive) on cohort fusion')
         call check_close(poly%cohort%leaf_resp_accum(1), (np1*0.2_wp + np2*0.9_wp)/nt, 1.0e-12_wp,      &
                          'leaf_resp_accum not nplant-weighted (extensive) on cohort fusion')
         call check_close(poly%cohort%stem_resp_accum(1), (np1*0.1_wp + np2*0.6_wp)/nt, 1.0e-12_wp,      &
                          'stem_resp_accum not nplant-weighted (extensive) on cohort fusion')
         call check_close(poly%cohort%root_resp_accum(1), (np1*0.3_wp + np2*0.8_wp)/nt, 1.0e-12_wp,      &
                          'root_resp_accum not nplant-weighted (extensive) on cohort fusion')
      end block

      !=== 2. fuse_cohorts reduces count to <= max_cohort, conserving N and AGB. ==========!
      call init_bare_ground(poly, cfg, 1_ik)
      do j = 1_ik, 200_ik
         !----- Low density so per-cohort LAI stays well under cohort_lai_cap (realistic). ---!
         call add_cohort(poly, cfg, 1_ik, 1_ik, 2.0e-4_wp, 2.0_wp + 0.55_wp*real(j-1_ik, wp))  ! 2..111 cm
      end do
      call finalize_init(poly)
      n0   = total_nplant(poly)
      agb0 = total_agb(poly)
      call fuse_cohorts(poly, cfg)
      call check(max_cohort_count(poly) <= abs(cfg%demography%max_cohort), 'fusion did not reach max_cohort')
      call check_close(total_nplant(poly), n0,   1.0e-10_wp, 'fusion broke nplant conservation')
      call check_close(total_agb(poly),    agb0, cfg%demography%conservation_tol, 'fusion broke AGB conservation')

      !=== 2b. The fusion tolerance is a share of the patch's top height. Two saplings 1.0 and 1.1 cm  !
      !        wide (about 0.15 m apart) stay apart in a stand of saplings, where 0.02 of the top is   !
      !        about 5 cm, and merge once the patch holds a 60 cm tree, where it is about 0.6 m. -----!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.1_wp)
      call finalize_init(poly)
      call fuse_cohorts(poly, cfg)
      call check(poly%cohort%n == 2_ik, 'saplings 0.15 m apart merged in a sapling stand')
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.1_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 1.0e-3_wp, 60.0_wp)
      call finalize_init(poly)
      call fuse_cohorts(poly, cfg)
      call check(poly%cohort%n == 2_ik, 'saplings 0.15 m apart not merged under a 60 cm tree')

      !=== 3. split_cohorts conserves N and AGB and increases the count. =====================!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 3_ik, 2.0_wp, 50.0_wp)   ! nplant*leaf_area >> cohort_lai_cap
      call finalize_init(poly)
      n0   = total_nplant(poly)
      agb0 = total_agb(poly)
      !----- The interception films are per plant: each daughter keeps its parent's, so the patch's     !
      !      film water is what it was. -------------------------------------------------------------!
      poly%cohort%leaf_surf_water(1) = 0.30_wp ; poly%cohort%wood_surf_water(1) = 0.05_wp
      film0 = poly%cohort%nplant(1) * 0.35_wp
      !----- The tissues away from their birth temperature, held to their energies, films included:   !
      !      each daughter keeps its parent's temperatures (the wood carbon it gained or lost held at  !
      !      them), and together they hold the parent's tissue energy. ------------------------------!
      poly%cohort%leaf_temp(1) = 301.0_wp
      poly%cohort%wood_temp(1) = 296.0_wp
      call hold_cohort_tissue_temps(poly%cohort, 1_ik)
      ltemp_exp = poly%cohort%leaf_temp(1)
      wtemp_exp = poly%cohort%wood_temp(1)
      e_sum0    = poly%cohort%nplant(1) * (poly%cohort%leaf_energy(1) + poly%cohort%wood_energy(1))
      call split_cohorts(poly, cfg)
      call check(poly%cohort%n >= 2_ik, 'split did not create a second cohort')
      call check_close(total_nplant(poly), n0,   1.0e-12_wp, 'split broke nplant conservation')
      call check_close(total_agb(poly),    agb0, cfg%demography%conservation_tol, 'split broke AGB conservation')
      call check_close(sum(poly%cohort%nplant(1:poly%cohort%n) * (poly%cohort%leaf_surf_water(1:poly%cohort%n)        &
                                                                + poly%cohort%wood_surf_water(1:poly%cohort%n))),   &
                       film0, 1.0e-12_wp, 'split broke canopy film water conservation')
      do j = 1_ik, poly%cohort%n
         call check(abs(poly%cohort%leaf_temp(j) - ltemp_exp) < 1.0e-9_wp .and.                            &
                    abs(poly%cohort%wood_temp(j) - wtemp_exp) < 1.0e-9_wp,                                 &
                    'split: a daughter lost its parent''s tissue temperatures')
      end do
      call check_close(sum(poly%cohort%nplant(1:poly%cohort%n) * (poly%cohort%leaf_energy(1:poly%cohort%n)            &
                                                                + poly%cohort%wood_energy(1:poly%cohort%n))),       &
                       e_sum0, 1.0e-12_wp, 'split broke tissue energy conservation')

      !=== 4. Carbon pools/traits thread correctly through the lockstep sort reorder. ==========!
      !     3 PFTs (sla 16/13/10) are inserted then height-sorted (reordered); re-deriving with   !
      !     each cohort's THREADED p_sla / p_aboveground_frac must reproduce its own PFT's         !
      !     on-allometry pools -- a mis-threaded per-cohort trait would misalign them.             !
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 1.0e-3_wp,  8.0_wp)   ! PFT 1 (shortest)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 1.0e-3_wp, 40.0_wp)   ! PFT 2 (tallest -> sorts first)
      call add_cohort(poly, cfg, 1_ik, 3_ik, 1.0e-3_wp, 20.0_wp)   ! PFT 3
      call finalize_init(poly)                                     ! sort_cohorts reorders EVERY field
      do j = 1_ik, poly%cohort%n
         call set_cohort_size(poly%cohort, j)                      ! re-derive with the threaded traits
         pf = poly%cohort%pft(j)
         call check_close(poly%cohort%leaf_carbon(j) * cfg%pft%sla(pf), poly%cohort%leaf_area(j),        &
                          1.0e-9_wp, 'leaf_carbon*sla(pft) /= leaf_area after reorder (trait mis-threaded?)')
         call check_close(poly%cohort%wood_carbon(j) * cfg%pft%allom(pf)%aboveground_frac, poly%cohort%agb(j), &
                          1.0e-9_wp, 'wood_carbon*aboveground_frac(pft) /= agb after reorder')
      end do

      !=== 5. Fusion conserves the FOUR prognostic carbon pools: the storage buffer is            !
      !     nplant-weighted and CARRIED, not snapped back onto allometry (set_cohort_size_from_carbon !
      !     path); the empirical set_cohort_size would reset the pools from size. ====================!
      block
         real(wp) :: lc_t, fc_t, wc_t, nc_t, nr5, nd5, alloc_store
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
         call finalize_init(poly)
         !----- Push one cohort's storage OFF allometry (the starvation buffer a reset would lose). -!
         poly%cohort%nonstructural_carbon(1) = 3.0_wp * poly%cohort%nonstructural_carbon(1)
         nr5 = poly%cohort%nplant(1) ; nd5 = poly%cohort%nplant(2)
         lc_t = nr5*poly%cohort%leaf_carbon(1)          + nd5*poly%cohort%leaf_carbon(2)
         fc_t = nr5*poly%cohort%fineroot_carbon(1)      + nd5*poly%cohort%fineroot_carbon(2)
         wc_t = nr5*poly%cohort%wood_carbon(1)          + nd5*poly%cohort%wood_carbon(2)
         nc_t = nr5*poly%cohort%nonstructural_carbon(1) + nd5*poly%cohort%nonstructural_carbon(2)
         call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%demography%conservation_tol)
         call check_close(poly%cohort%nplant(1)*poly%cohort%leaf_carbon(1),          lc_t, 1.0e-9_wp, &
                          'carbon fuse: leaf_carbon not conserved')
         call check_close(poly%cohort%nplant(1)*poly%cohort%fineroot_carbon(1),      fc_t, 1.0e-9_wp, &
                          'carbon fuse: fineroot_carbon not conserved')
         call check_close(poly%cohort%nplant(1)*poly%cohort%wood_carbon(1),          wc_t, 1.0e-9_wp, &
                          'carbon fuse: wood_carbon not conserved')
         call check_close(poly%cohort%nplant(1)*poly%cohort%nonstructural_carbon(1), nc_t, 1.0e-9_wp, &
                          'carbon fuse: nonstructural (storage) not conserved')
         !----- Prove storage was CARRIED, not reset to the allometric storage target. -------------!
         alloc_store = storage_carbon_target(cfg%pft%storage_cushion(1), poly%cohort%leaf_carbon(1),   &
                                             poly%cohort%fineroot_carbon(1))
         call check(abs(poly%cohort%nonstructural_carbon(1) - alloc_store) > 1.0e-6_wp,               &
                    'carbon fuse must CARRY off-allometry storage, not snap it to allometry')
      end block

      !=== 6. Fusion of cohorts with different SLA conserves leaf AREA: sla blends on leaf          !
      !     carbon, the quantity it converts, so sum(nplant*leaf_carbon*sla) is what the fused        !
      !     cohort shows. A leaf-area-weighted sla would give sum(LA*sla)/sum(LA) instead. The leaf    !
      !     turnover rate is the share of the leaf carbon lost a year, so it too blends on leaf        !
      !     carbon, and the fused cohort loses the leaf carbon its two parents did; with the two      !
      !     SLAs differing, a leaf-area weight would not. =============================================!
      block
         real(wp) :: lai0, loss0
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
         call finalize_init(poly)
         poly%cohort%sla(1) = 12.0_wp ; poly%cohort%sla(2) = 24.0_wp      ! shade-acclimated donor
         poly%cohort%leaf_turnover_rate(1) = 1.0_wp ; poly%cohort%leaf_turnover_rate(2) = 0.4_wp
         call set_cohort_size_from_carbon(poly%cohort, 1_ik)
         call set_cohort_size_from_carbon(poly%cohort, 2_ik)
         lai0 = poly%cohort%nplant(1)*poly%cohort%leaf_area(1) + poly%cohort%nplant(2)*poly%cohort%leaf_area(2)
         loss0 = poly%cohort%nplant(1)*poly%cohort%leaf_carbon(1)*poly%cohort%leaf_turnover_rate(1)       &
               + poly%cohort%nplant(2)*poly%cohort%leaf_carbon(2)*poly%cohort%leaf_turnover_rate(2)
         call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%demography%conservation_tol)
         call check_close(poly%cohort%nplant(1)*poly%cohort%leaf_area(1), lai0, 1.0e-12_wp,              &
                          'fusion of cohorts with different SLA must conserve leaf area')
         call check_close(poly%cohort%nplant(1)*poly%cohort%leaf_carbon(1)*poly%cohort%leaf_turnover_rate(1), &
                          loss0, 1.0e-12_wp, 'fusion must conserve the leaf carbon turned over a year')
         !----- The same at any size: two LEAFLESS cohorts, whose leaf carbon has decayed to ~1e-193   !
         !      (as in a BCI run). A blend that skipped weights below 1e-30 would keep one SLA here. -!
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
         call finalize_init(poly)
         poly%cohort%sla(1) = 12.0_wp ; poly%cohort%sla(2) = 24.0_wp
         poly%cohort%leaf_carbon(1) = 1.0e-193_wp ; poly%cohort%leaf_carbon(2) = 3.0e-193_wp
         call set_cohort_size_from_carbon(poly%cohort, 1_ik)
         call set_cohort_size_from_carbon(poly%cohort, 2_ik)
         lai0 = poly%cohort%nplant(1)*poly%cohort%leaf_area(1) + poly%cohort%nplant(2)*poly%cohort%leaf_area(2)
         call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%demography%conservation_tol)
         call check(abs(poly%cohort%nplant(1)*poly%cohort%leaf_area(1) - lai0) <= 1.0e-12_wp * lai0,      &
                    'fusion of leafless cohorts with different SLA must conserve their leaf area too')
      end block

      !=== 7. Geometry has ONE writer: a cohort placed on the allometry from its dbh carries        !
      !     exactly the geometry its carbon gives, so a restart that re-derives it gets the same bits. !
      block
         real(wp) :: g0(5)
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, 0.4_wp, 17.3_wp)
         call finalize_init(poly)
         g0 = [poly%cohort%dbh(1), poly%cohort%height(1), poly%cohort%basal_area(1),                    &
               poly%cohort%agb(1), poly%cohort%leaf_area(1)]
         call set_cohort_size_from_carbon(poly%cohort, 1_ik)
         call check(all([poly%cohort%dbh(1), poly%cohort%height(1), poly%cohort%basal_area(1),          &
                         poly%cohort%agb(1), poly%cohort%leaf_area(1)] == g0),                          &
                    'birth geometry must be the geometry derived from the birth carbon, bit for bit')

         !=== 8. The bounds guard names the first field out of bounds, and passes a sound cohort. ==!
         call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'a fresh cohort is in bounds')
         poly%cohort%leaf_surf_water(1) = -1.0e-30_wp
         call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'leaf_surf_water',  &
                    'a negative leaf film is out of bounds, however small')
         poly%cohort%nplant(1) = ieee_value(1.0_wp, ieee_quiet_nan)
         call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'nplant',           &
                    'a NaN density is out of bounds, and checked first')
      end block

      !=== 9. The growth window's sum is the buffer's: once a large sample has been evicted, a       !
      !     window of zeros averages exactly zero; and a recruit in a reused slot starts empty. ======!
      block
         type(cohort_deriv_block) :: deriv
         integer(ik) :: k, pos, nwin
         real(wp)    :: sample
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.4_wp, 17.3_wp)
         call finalize_init(poly)
         nwin = int(size(poly%cohort%growth_hist, 1), ik)
         call cohort_deriv_alloc(deriv, poly%cohort%n)
         pos = 0_ik
         do k = 1_ik, 3_ik * nwin
            pos = mod(pos, nwin) + 1_ik
            sample = 0.0_wp
            if (k <= 2_ik * nwin) sample = 0.1_wp
            if (k == 1_ik)        sample = 1.0e8_wp
            call fill_cohort_deriv(poly%cohort, 1_ik, deriv, cfg%dt_years, 0.0_wp, sample,             &
                                   poly%cohort%leaf_carbon(1), poly%cohort%fineroot_carbon(1),         &
                                   poly%cohort%wood_carbon(1), poly%cohort%nonstructural_carbon(1),    &
                                   nwin, pos)
         end do
         call check(poly%cohort%growth_avg(1) == 0.0_wp,                                                &
                    'a window of zero growth averages exactly zero after a large sample left it')
         call check(poly%cohort%growth_accum(1) == sum(poly%cohort%growth_hist(1:nwin, 1)),             &
                    'the window''s sum is the buffer''s')
         !----- Any birth state: the window is the point. ---------------------------------------!
         call init_cohort(poly%cohort, 1_ik, cfg%pft, 1_ik, 1_ik, 0.1_wp, 1.0_wp, 288.15_wp, -0.1_wp)
         call check(all(poly%cohort%growth_hist(1:nwin, 1) == 0.0_wp) .and.                             &
                    poly%cohort%growth_count(1) == 0_ik, 'a recruit in a reused slot starts with an empty window')
      end block
   end subroutine run_cohort_operators


end module test_slow_update_cohorts

!==========================================================================================!
! Treefall gaps: one gap per site, the survivors' density, the donors' remainder, the area !
! conserved and the diagnostic rows.                                                       !
!==========================================================================================!
module test_slow_update_gaps
   use ed_kinds,              only : wp, ik
   use ed_config,             only : ed_config_t
   use ed_patch_state_vars,   only : DIST_TREEFALL
   use ed_polygon_state_vars, only : polygon_t
   use ed_slow_update, only : apply_patch_disturbance
   use ed_diag_vars,                   only : PA_DISTURB_AREA, PA_LE, PA_MORT_C_CULL, alloc_patch_diag, patch_diag_value
   use ed_constants,                   only : yr_sec
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_output_reduce, only : total_area, total_nplant
   use ed_test_support, only : banner, build_test_config, check, check_close
   use ed_slow_state,  only : kill_report_t, decide_disturbance
   implicit none
   private

   public :: run_treefall_gaps
   type(kill_report_t)  :: kills   !< what the cull and the treefall killed
   real(wp)             :: dist_frac          !< the treefall decision's share
   logical, allocatable :: dist_survives(:)   !< and its survivors
   real(wp)            :: shed_water   !< an operator's shed water, unread here

   type(ed_config_t) :: cfg
   type(polygon_t)  :: poly
   integer(ik)         :: ig, ip, i0, i1, i, n_gap_cohorts, n_slot
   real(wp)            :: disturb_rate(8)
   real(wp)            :: n_before, h_tall, h_short, film_before, film_after, frac_expect
   real(wp)            :: le(8), le_site, cull_before, cull_after, w_canopy

contains

   subroutine run_treefall_gaps()
      call banner('treefall patch disturbance')
      cfg = build_test_config()

      !----- One patch with a tall canopy cohort and a short understorey cohort. -------------!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall  (height >> threshold)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short (height <  threshold)
      call finalize_init(poly)
      h_tall  = poly%cohort%height(1)                                ! sorted tallest-first
      h_short = poly%cohort%height(poly%cohort%n)
      call check(h_tall  >= cfg%demography%disturbance_survive_height, 'tall cohort should exceed threshold')
      call check(h_short <  cfg%demography%disturbance_survive_height, 'short cohort should be below threshold')
      n_before = total_nplant(poly)
      !----- The survivor's film [kg/plant] goes into the gap with it; the site total of the        !
      !      SURVIVOR's film is conserved across the split. -----------------------------------------!
      poly%cohort%leaf_surf_water(poly%cohort%n) = 0.25_wp
      film_before = poly%patch%area(1) * poly%cohort%nplant(poly%cohort%n) * 0.25_wp

      !----- Turn the per-patch diagnostic block ON. It defaults OFF (a run with no diagnostics    !
      !      allocates nothing and every entry point no-ops), so a test that wants to read a slot     !
      !      must ask for it: a slot nothing reads can sit at a silent zero unseen. -------------------!
      call alloc_patch_diag(poly%patch%diag, poly%patch%n + 1_ik, .true.)
      !----- The donor's own shed rate, which the gap must not inherit, and the tissue water the dying  !
      !      canopy holds, which it must receive. ------------------------------------------------------!
      poly%patch%shed_water_rate(1) = 7.0_wp
      w_canopy = poly%cohort%nplant(1) * (poly%cohort%leaf_water_mass(1) + poly%cohort%wood_water_mass(1))
      call check(w_canopy > 0.0_wp, 'the dying canopy holds tissue water')

      call decide_disturbance(poly, cfg, 1.0_wp, dist_frac, dist_survives)
      call apply_patch_disturbance(poly, cfg, dist_frac, dist_survives, shed_water, kills)

      !----- A new age-0 treefall gap patch was opened; area is conserved. --------------------!
      call check(poly%patch%n == 2_ik, 'disturbance should add exactly one gap patch')
      call check(poly%site%n == 1_ik .and. poly%site%patch_count(1) == 2_ik .and.                             &
                 all(poly%patch%owner_site(1:2) == 1_ik), 'the gap patch is in the polygon''s one site')
      call check_close(total_area(poly), 1.0_wp, 1.0e-9_wp, 'disturbance broke area conservation')

      !----- The disturbed-area flux is REPORTED. A slot declared and never written reads as a       !
      !      silent zero -- which no conservation check can see, because zero disturbed area is a    !
      !      perfectly conservative answer. Assert the VALUE, not merely that it is finite.          !
      !                                                                                              !
      !      Asserted through patch_diag_value -- the SAME reader the output layer uses -- and not   !
      !      off the raw accumulator. The block is (value, weight) and the reader returns            !
      !      value/weight, so reading `v` directly tests the contribution rather than the number the !
      !      variable emits, and can pass a contribution that is a factor yr_sec wrong. The weight   !
      !      is what a year of fast steps accumulates, because that is the window a once-a-year      !
      !      disturbance is normalized over; the answer is then the declared [1/yr] rate.            !
      poly%patch%diag%n      = poly%patch%n
      poly%patch%diag%w(1:poly%patch%n) = yr_sec
      frac_expect = 1.0_wp - exp(-cfg%demography%patch_disturbance_rate * 1.0_wp)
      call check(frac_expect > 1.0e-6_wp, 'the fixture must actually disturb something')
      call patch_diag_value(poly%patch%diag, PA_DISTURB_AREA, disturb_rate, n_slot)
      call check_close(disturb_rate(1), frac_expect, 1.0e-12_wp,                                       &
                       'PA_DISTURB_AREA must report the fraction the donor lost, in 1/yr')

      ig = 0_ik
      do ip = 1_ik, poly%patch%n
         if (poly%patch%dist_type(ip) == DIST_TREEFALL) ig = ip
      end do
      call check(ig > 0_ik, 'no treefall gap patch found')
      call check_close(poly%patch%age(ig), 0.0_wp, 1.0e-12_wp, 'gap patch should have age 0')
      call check_close(poly%patch%shed_water_rate(ig), w_canopy / cfg%dt_slow, 1.0e-12_wp,             &
                       'gap patch starts a fresh shed_water_rate, like age, and receives the dead canopy''s water')

      !----- The gap holds ONLY survivors (height < threshold): the short cohort. -------------!
      i0 = poly%patch%cohort_offset(ig) ; i1 = i0 + poly%patch%cohort_count(ig) - 1_ik
      n_gap_cohorts = poly%patch%cohort_count(ig)
      call check(n_gap_cohorts == 1_ik, 'gap should contain exactly the one survivor cohort')
      film_after = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip) ; i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            if (poly%cohort%height(i) < cfg%demography%disturbance_survive_height)                                &
               film_after = film_after + poly%patch%area(ip) * poly%cohort%nplant(i) * poly%cohort%leaf_surf_water(i)
         end do
      end do
      call check_close(film_after, film_before, 1.0e-12_wp, 'disturbance broke the survivor film-water conservation')
      i0 = poly%patch%cohort_offset(ig) ; i1 = i0 + poly%patch%cohort_count(ig) - 1_ik
      do i = i0, i1
         call check(poly%cohort%height(i) < cfg%demography%disturbance_survive_height,                        &
                    'tall canopy cohort must not survive in the gap')
      end do

      !----- Net plant number drops (the disturbed canopy fraction is killed). ----------------!
      call check(total_nplant(poly) < n_before, 'canopy mortality should reduce plant number')

      !----- TWO donors of unequal area and unequal film: the gap's survivor copies take their donors'  !
      !      densities diluted by frac*area(d)/new_area and keep their per-plant films, so the site's   !
      !      survivor film water is what it was. -------------------------------------------------------!
      call init_bare_ground(poly, cfg, 2_ik)
      poly%patch%area(1) = 0.6_wp ; poly%patch%area(2) = 0.4_wp
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall, dies
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short survivor, WET
      call add_cohort(poly, cfg, 2_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall, dies
      call add_cohort(poly, cfg, 2_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short survivor, dry
      call finalize_init(poly)
      film_before = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip) ; i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            if (poly%cohort%height(i) < cfg%demography%disturbance_survive_height .and. ip == 1_ik) then
               poly%cohort%leaf_surf_water(i) = 0.25_wp
               film_before = film_before + poly%patch%area(ip) * poly%cohort%nplant(i) * 0.25_wp
            end if
         end do
      end do
      !----- Inside a step (the C API's patch operators), the step's patch diagnostics must survive the !
      !      carving: the gap's ground was donor ground for the step so far, so the site mean          !
      !      Sum(area * value) is unchanged and the gap reads the donors' area-weighted value; a        !
      !      cleared gap slot would read 0 there and pull the site mean low by the disturbed fraction.  !
      !      Donor latent heat 100 and 300 W/m2 on areas 0.6 and 0.4: 180 W/m2 before and after. -------!
      call alloc_patch_diag(poly%patch%diag, poly%patch%n + 1_ik, .true.)
      poly%patch%diag%n = poly%patch%n
      poly%patch%diag%w(1:2) = cfg%dt_slow
      poly%patch%diag%v(PA_LE, 1) = 100.0_wp * cfg%dt_slow
      poly%patch%diag%v(PA_LE, 2) = 300.0_wp * cfg%dt_slow
      call decide_disturbance(poly, cfg, 0.5_wp, dist_frac, dist_survives)
      call apply_patch_disturbance(poly, cfg, dist_frac, dist_survives, shed_water, kills)
      call patch_diag_value(poly%patch%diag, PA_LE, le, n_slot)
      call check(n_slot == poly%patch%n, 'the gap''s diagnostic slot is readable right after the disturbance')
      le_site = sum(poly%patch%area(1:n_slot) * le(1:n_slot))
      call check_close(le_site, 180.0_wp, 1.0e-9_wp,                                                   &
                       'disturbance leaves the area-weighted site mean of a patch diagnostic unchanged')
      ig = 0_ik
      do ip = 1_ik, poly%patch%n
         if (poly%patch%dist_type(ip) == DIST_TREEFALL) ig = ip
      end do
      call check_close(le(ig), 180.0_wp, 1.0e-9_wp, 'the gap inherits the donors'' area-weighted diagnostics')
      film_after = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip) ; i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            film_after = film_after + poly%patch%area(ip) * poly%cohort%nplant(i) * poly%cohort%leaf_surf_water(i)
         end do
      end do
      call check_close(total_area(poly), 1.0_wp, 1.0e-9_wp, 'two-donor disturbance broke area conservation')
      call check_close(film_after, film_before, 1.0e-12_wp,                                            &
                       'two-donor disturbance: survivor films go into the gap with their plants')

      !----- At a calendar boundary, as a run does it: the step's diagnostics were read and reset, so  !
      !      every weight is 0 and the slots hold only this boundary's events, here a cull. The gap     !
      !      takes its area-weighted share, so the stand total Sum(area * v) is unchanged. -------------!
      poly%patch%diag%w(1:poly%patch%n) = 0.0_wp
      poly%patch%diag%v(:, 1:poly%patch%n) = 0.0_wp
      do ip = 1_ik, poly%patch%n
         poly%patch%diag%v(PA_MORT_C_CULL, ip) = 1.0e-3_wp * real(ip, wp)
      end do
      cull_before = sum(poly%patch%area(1:poly%patch%n) * poly%patch%diag%v(PA_MORT_C_CULL, 1:poly%patch%n))
      call decide_disturbance(poly, cfg, 0.5_wp, dist_frac, dist_survives)
      call apply_patch_disturbance(poly, cfg, dist_frac, dist_survives, shed_water, kills)
      cull_after = sum(poly%patch%area(1:poly%patch%n) * poly%patch%diag%v(PA_MORT_C_CULL, 1:poly%patch%n))
      call check_close(cull_after, cull_before, 1.0e-15_wp,                                            &
                       'a boundary disturbance keeps the stand total of the events it carves')
      call check(all(poly%patch%diag%w(1:poly%patch%n) == 0.0_wp), 'and leaves every weight at 0')
   end subroutine run_treefall_gaps


end module test_slow_update_gaps

!==========================================================================================!
! Patch fusion and termination: plant number and area conserved, a sliver merged rather    !
! than dropped, the stores riding the patch lockstep, the tolerance ceiling.               !
!==========================================================================================!
module test_slow_update_patches
   use ed_kinds,            only : wp, ik
   use ed_config,           only : ed_config_t, derive_parameters
   use ed_polygon_state_vars,        only : polygon_t
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_slow_update, only : fuse_patches, terminate_patches, sort_patches, patch_light_profile
   use ed_output_reduce, only : total_nplant, total_area
   use ed_cohort_state_vars, only : cohort_film_water
   use ed_test_support, only : banner, build_test_config, check, check_close
   implicit none
   private

   public :: run_patch_operators

   type(ed_config_t) :: cfg
   type(polygon_t)  :: poly
   real(wp)            :: n0, w_before, w_after, cas_before, cas_after, davg, film0
   real(wp), allocatable :: light1(:), light2(:)
   integer(ik)         :: ip

contains

   subroutine run_patch_operators()
      call banner('patch fusion & termination conservation')
      cfg = build_test_config()

      !=== Two identical patches must fuse, conserving the polygon's N and total area. =========!
      call init_bare_ground(poly, cfg, 2_ik)   ! area 0.5 each
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.4_wp, 18.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.7_wp,  6.0_wp)
      call add_cohort(poly, cfg, 2_ik, 2_ik, 0.4_wp, 18.0_wp)
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.7_wp,  6.0_wp)
      call finalize_init(poly)
      n0 = total_nplant(poly)
      !----- Films are per plant: fusion rescales the densities and leaves the films alone, and the   !
      !      site's film water, the area- and nplant-weighted sum, is what it was. -------------------!
      poly%cohort%leaf_surf_water(1) = 0.20_wp ; poly%cohort%wood_surf_water(2) = 0.10_wp   ! patch 1
      poly%cohort%leaf_surf_water(3) = 0.40_wp                                               ! patch 2
      film0 = site_film(poly)
      call fuse_patches(poly, cfg)
      call check(poly%patch%n == 1_ik, 'identical patches should fuse to one')
      call check_close(total_area(poly),   1.0_wp, 1.0e-9_wp, 'total area not conserved by fusion')
      call check_close(total_nplant(poly), n0,     1.0e-9_wp, 'the polygon''s plant number not conserved by fusion')
      call check(poly%site%patch_offset(1) == 1_ik .and. poly%site%patch_count(1) == poly%patch%n,            &
                 'the site map covers the fused patch block')
      call check_close(site_film(poly), film0, 1.0e-12_wp, 'patch fusion broke canopy film water conservation')

      !=== A negligible-area patch is removed and the remaining areas renormalize to 1. ======!
      !----- TWO doomed patches, not one: the merge loop reads the survivor's CSR slice, so a       !
      !      second sliver is the case a stale map would corrupt.  ----------------------------------!
      call init_bare_ground(poly, cfg, 4_ik)
      poly%patch%area(1) = 0.6_wp
      poly%patch%area(2) = 0.4_wp - 2.0e-6_wp
      poly%patch%area(3) = 1.0e-6_wp                  ! below min_patch_area
      poly%patch%area(4) = 1.0e-6_wp                  ! and another
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(poly, cfg, 3_ik, 1_ik, 0.4_wp,  8.0_wp)
      call add_cohort(poly, cfg, 4_ik, 1_ik, 0.3_wp, 15.0_wp)   ! taller than the survivor's own cohort
      call finalize_init(poly)
      !----- Give the doomed patch content UNLIKE its neighbours', which is the realistic case (a    !
      !      sliver is usually a fresh gap). Dropping it and renormalizing the survivors would       !
      !      destroy that content outright; it is merged into the largest survivor, so the site      !
      !      total must come through untouched. ------------------------------------------------------!
      poly%patch%soil_water(1)%theta(1) = 0.20_wp
      poly%patch%soil_water(2)%theta(1) = 0.20_wp
      poly%patch%soil_water(3)%theta(1) = 0.45_wp        ! the slivers are much wetter
      poly%patch%soil_water(4)%theta(1) = 0.50_wp
      n0 = total_nplant(poly)
      w_before = sum(poly%patch%area(1:4) * [(poly%patch%soil_water(ip)%theta(1), ip = 1, 4)])
      call terminate_patches(poly, cfg)
      call check(poly%patch%n == 2_ik, 'sub-threshold patch not removed')
      call check_close(total_area(poly), 1.0_wp, 1.0e-9_wp, 'areas not renormalized to 1')
      w_after = sum(poly%patch%area(1:poly%patch%n)                                                   &
                    * [(poly%patch%soil_water(ip)%theta(1), ip = 1, poly%patch%n)])
      call check_close(w_after, w_before, 1.0e-12_wp,                                                 &
                       'terminate_patches must MERGE the sliver, not discard its content')
      call check_close(total_nplant(poly), n0, 1.0e-12_wp,                                            &
                       'terminate_patches must carry the slivers'' plants into the survivor')
      call check(tallest_first(poly), 'terminate_patches leaves each patch''s cohorts tallest first, '  &
                 // 'the order the cohort fusion after it scans in')

      !=== Fast-biophysics reservoirs ride the patch lockstep. ================================!
      !    (a) area-weighted MERGE on fusion conserves the area-extensive store; (b) SORT keeps  !
      !    each reservoir bound to its patch. Both are silent-wrong if the lockstep misses them. !
      call init_bare_ground(poly, cfg, 2_ik)
      poly%patch%area(1) = 0.75_wp ; poly%patch%area(2) = 0.25_wp
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 12.0_wp)   ! identical cohort in each patch
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.5_wp, 12.0_wp)   ! -> identical light profiles -> fuse
      call finalize_init(poly)
      poly%patch%soil_water(1)%theta(1)       = 0.20_wp   ; poly%patch%soil_water(2)%theta(1)       = 0.40_wp
      poly%patch%cas(1)%can_enthalpy      = 1000.0_wp ; poly%patch%cas(2)%can_enthalpy      = 2000.0_wp
      poly%patch%soil_energy(1)%soil_energy(3) = 1.0e6_wp  ; poly%patch%soil_energy(2)%soil_energy(3) = 3.0e6_wp
      !----- DIFFERENT canopy depths, or this fixture cannot tell the two weightings apart. The      !
      !      canopy air is intensive per kg of AIR, and the air mass is area x depth, so the merged  !
      !      value is mass-weighted: (0.75*30*1000 + 0.25*10*2000)/(0.75*30 + 0.25*10) = 1100, not   !
      !      the 1250 an area weighting gives. Equal depths are exactly the case where the wrong     !
      !      answer and the right one agree. ---------------------------------------------------------!
      poly%patch%cas(1)%can_depth = 30.0_wp ; poly%patch%cas(2)%can_depth = 10.0_wp
      cas_before = poly%patch%area(1) * poly%patch%cas(1)%can_depth * poly%patch%cas(1)%can_enthalpy  &
                 + poly%patch%area(2) * poly%patch%cas(2)%can_depth * poly%patch%cas(2)%can_enthalpy
      call fuse_patches(poly, cfg)
      call check(poly%patch%n == 1_ik, 'reservoir patches should fuse to one')
      call check_close(poly%patch%soil_water(1)%theta(1),       0.25_wp,   1.0e-9_wp, 'soil water not area-weighted on fusion')
      call check_close(poly%patch%soil_energy(1)%soil_energy(3), 1.5e6_wp,  1.0e-3_wp, 'soil energy not area-weighted on fusion')
      call check_close(poly%patch%cas(1)%can_depth,         25.0_wp,   1.0e-9_wp, 'CAS depth not area-weighted on fusion')
      call check_close(poly%patch%cas(1)%can_enthalpy,      1100.0_wp, 1.0e-9_wp, 'CAS enthalpy not AIR-MASS-weighted on fusion')
      !----- and the point of that weighting: the extensive content survives the merge. -------------!
      cas_after = poly%patch%area(1) * poly%patch%cas(1)%can_depth * poly%patch%cas(1)%can_enthalpy
      call check_close(cas_after, cas_before, 1.0e-12_wp, 'patch fusion must conserve canopy-air energy')

      call init_bare_ground(poly, cfg, 2_ik)
      poly%patch%age(1) = 1.0_wp ; poly%patch%age(2) = 5.0_wp             ! patch 2 is older
      poly%patch%soil_water(1)%theta(1) = 0.11_wp ; poly%patch%soil_water(2)%theta(1) = 0.55_wp
      call finalize_init(poly)
      call sort_patches(poly)                                             ! sorts age-descending: older first
      call check_close(poly%patch%age(1), 5.0_wp, 1.0e-12_wp, 'sort did not put the older patch first')
      call check_close(poly%patch%soil_water(1)%theta(1), 0.55_wp, 1.0e-12_wp, 'reservoir did not follow its patch through sort')

      !=== The tolerance ceiling: dissimilar patches stay apart even above max_patch. ===========!
      !    A closed stand and a sparse understory, and max_patch = 1 asking for one patch. Under the !
      !    shipped ceiling (0.15) they stay two; under an uncapped schedule, 0.10 x 1.5 per pass     !
      !    (0.76 by the sixth pass), they fuse.                                                      !
      cfg = build_test_config()
      cfg%demography%max_patch = 1_ik
      call check_close(cfg%demography%patch_light_tol                                                    &
                       * cfg%demography%patch_light_tol_mult**(cfg%demography%n_patch_fusion_iter - 1_ik),  &
                       cfg%demography%patch_light_tol_max, 1.0e-12_wp, 'the last pass must use exactly the ceiling')
      call two_unlike_patches()
      allocate(light1(cfg%demography%n_height_layers), light2(cfg%demography%n_height_layers))
      call patch_light_profile(poly, cfg, 1_ik, light1)
      call patch_light_profile(poly, cfg, 2_ik, light2)
      davg = sum(abs(light1 - light2), mask=(light1 < 1.0_wp - 1.0e-7_wp .or. light2 < 1.0_wp - 1.0e-7_wp)) &
           / real(count(light1 < 1.0_wp - 1.0e-7_wp .or. light2 < 1.0_wp - 1.0e-7_wp), wp)
      call check(davg > 0.3_wp .and. davg < 0.7_wp,                                                  &
                 'fixture: the two light profiles must differ by between 0.3 and 0.7 on average')
      call fuse_patches(poly, cfg)
      call check(poly%patch%n == 2_ik, 'patches more different than the ceiling must stay apart above max_patch')
      !----- The uncapped schedule, spelled as a ceiling: 0.10 * 1.5**5 over six passes is a step of 1.5. -!
      cfg%demography%patch_light_tol_max = 0.10_wp * 1.5_wp**5
      call derive_parameters(cfg)
      call check_close(cfg%demography%patch_light_tol_mult, 1.5_wp, 1.0e-12_wp, 'ceiling 0.759375 must give the old 1.5 step')
      call two_unlike_patches()
      call fuse_patches(poly, cfg)
      call check(poly%patch%n == 1_ik, 'without the ceiling the same two patches fuse down to max_patch')
   end subroutine run_patch_operators

   !----- Each patch's cohorts in height order, tallest first. -------------------------------!
   logical function tallest_first(p)
      type(polygon_t), intent(in) :: p
      integer(ik) :: q, i
      tallest_first = .true.
      do q = 1_ik, p%patch%n
         do i = p%patch%cohort_offset(q), p%patch%cohort_offset(q) + p%patch%cohort_count(q) - 2_ik
            if (p%cohort%height(i) < p%cohort%height(i + 1_ik)) tallest_first = .false.
         end do
      end do
   end function tallest_first
   !----- Two equal-area patches: a closed stand of 30 cm trees, and 5 cm saplings alone. ------!
   subroutine two_unlike_patches()
      call init_bare_ground(poly, cfg, 2_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.012_wp, 30.0_wp)
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.050_wp,  5.0_wp)
      call finalize_init(poly)
   end subroutine two_unlike_patches
   !----- The site's canopy film water [kg/m2 of site]: each cohort's, area-weighted by its patch. -!
   real(wp) function site_film(p) result(w)
      type(polygon_t), intent(in) :: p
      integer(ik) :: i
      w = 0.0_wp
      do i = 1_ik, p%cohort%n
         w = w + p%patch%area(p%cohort%owner_patch(i)) * cohort_film_water(p%cohort, i)
      end do
   end function site_film

end module test_slow_update_patches

!==========================================================================================!
! A recruit is born in its patch's state: the canopy air's temperature and the top soil    !
! layer's potential, not a global constant.                                                !
!==========================================================================================!
module test_slow_update_birth
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_polygon_state_vars,   only : polygon_t
   use ed_patch_state_vars,     only : patch_physical_state
   use ed_soil_derive,          only : top_soil_psi
   use ed_cohort_state_vars,    only : init_cohort, expand_cohort_block
   use ed_plant_hydr_lib,       only : water_curve_t, tissue_water_to_psi
   use ed_init_bare_ground,     only : init_bare_ground, add_cohort, finalize_init
   use ed_test_support,         only : banner, build_test_config, check, check_close
   implicit none
   private

   public :: run_birth_state
   type(ed_config_t)   :: cfg
   type(polygon_t)     :: poly

contains

   subroutine run_birth_state()
      call banner('a recruit is born in its patch''s state')
      cfg = build_test_config()
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.3_wp, 18.0_wp)
      call finalize_init(poly)
      call expand_cohort_block(poly%cohort, poly%cohort%n + 1_ik)   ! room for the recruit
      call check_birth_temp()
   end subroutine run_birth_state

   !----- A cohort is born in its patch's state (patch_physical_state): at the canopy-air temperature, !
   !      its tissues at the top soil layer's potential. -------------------------------------------!
   subroutine check_birth_temp()
      integer(ik) :: m
      real(wp)    :: t_birth, psi_birth, t_save
      type(water_curve_t) :: leaf_curve, wood_curve
      m = poly%cohort%n + 1_ik
      t_save = poly%patch%cas(1)%can_temp
      poly%patch%cas(1)%can_temp = 268.0_wp        ! a cold, seeded patch
      call patch_physical_state(poly%patch%cas(1), poly%patch%soil_water(1), poly%site%surface(1),       &
                                t_birth, psi_birth)
      call check_close(t_birth, 268.0_wp, 1.0e-12_wp, 'birth: a seeded patch gives its canopy-air temperature')
      call check_close(psi_birth, top_soil_psi(poly%site%surface(1)%soil_profile, poly%patch%soil_water(1)%theta(1)), &
                       1.0e-12_wp,                                                                        &
                       'birth: a seeded patch gives its top soil layer''s potential')
      call init_cohort(poly%cohort, m, cfg%pft, 1_ik, 1_ik, 0.01_wp, 0.45_wp, t_birth, psi_birth)
      call check_close(poly%cohort%leaf_temp(m), 268.0_wp, 1.0e-12_wp, 'birth: leaf starts at the patch temperature')
      call check_close(poly%cohort%wood_temp(m), 268.0_wp, 1.0e-12_wp, 'birth: wood starts at the patch temperature')
      leaf_curve = cfg%pft%hydro(1_ik)%leaf_curve
      wood_curve = cfg%pft%hydro(1_ik)%wood_curve
      call check_close(tissue_water_to_psi(poly%cohort%leaf_water_mass(m), leaf_curve,                &
                       poly%cohort%leaf_carbon(m)), psi_birth, 1.0e-9_wp,                               &
                       'birth: the leaves hold the water of that potential')
      poly%patch%cas(1)%can_temp = t_save
   end subroutine check_birth_temp

end module test_slow_update_birth

!==========================================================================================!
! The restructuring through the slow ledger, with soil carbon on: one month-and-year       !
! boundary and the step after it. Every phase runs and closes on all three currencies by   !
! the report's own test (slow_ledger_failures), each operator declares the carbon it       !
! moves, ids stay distinct and a new id is never an old one, and the boundary's dead       !
! carbon reaches the soil with the next step (docs/science/conservation.md §3, §5).        !
!==========================================================================================!
module test_slow_update_ledger
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area
   use ed_biogeochem_params,    only : N_SOIL_POOL
   use ed_init_bare_ground,     only : init_bare_ground, add_cohort, finalize_init
   use ed_slow_dynamics,        only : advance_boundary_dynamics, advance_slow_dynamics
   use ed_slow_ledger,          only : slow_ledger_t, slow_ledger_open, slow_ledger_failures,             &
                                       litter_carbon_total, SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW,         &
                                       SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT, SLOW_PHASE_DISTURB,         &
                                       SLOW_PHASE_PATCH, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC
   use ed_test_support,         only : banner, build_test_config, check
   implicit none
   private

   public :: run_restructure_ledger

contains

   subroutine run_restructure_ledger()
      integer(ik), parameter :: AT_BOUNDARY(5) = [SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT, SLOW_PHASE_DISTURB,   &
                                                  SLOW_PHASE_PATCH, SLOW_PHASE_CANOPY]
      integer(ik), parameter :: IN_STEP(4) = [SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW, SLOW_PHASE_CANOPY,        &
                                              SLOW_PHASE_SOILC]
      character(len=8), parameter :: BOUNDARY_NAME(5) = [character(len=8) :: 'recruit', 'cohort', 'disturb',  &
                                                         'patch', 'canopy']
      character(len=8), parameter :: STEP_NAME(4) = [character(len=8) :: 'allocate', 'grow', 'canopy', 'soilc']
      type(ed_config_t)   :: cfg
      type(polygon_t)     :: poly
      type(slow_ledger_t) :: led
      integer(ik), allocatable :: patch_id0(:), cohort_id0(:)
      integer(ik) :: k, culled_id
      real(wp)    :: soil0, litter0

      call banner('the restructuring through the slow ledger, soil carbon on')
      cfg = build_test_config()
      cfg%soil_carbon%soil_carbon_on = .true.
      cfg%demography%patch_disturbance_rate = 0.2_wp   ! a gap large enough to kill canopy trees

      !----- Canopy trees the treefall kills, understory it spares, a cohort at the density   !
      !      floor the cull removes, and a recruit pool that founds a cohort. ----------------!
      call init_bare_ground(poly, cfg, 3_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.05_wp, 40.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.50_wp,  3.0_wp)
      call add_cohort(poly, cfg, 2_ik, 2_ik, 0.05_wp, 35.0_wp)
      call add_cohort(poly, cfg, 3_ik, 1_ik, 0.80_wp,  2.0_wp)
      call finalize_init(poly)
      poly%cohort%nplant(poly%patch%cohort_offset(3)) = cfg%demography%negligible_nplant
      culled_id = poly%cohort%global_id(poly%patch%cohort_offset(3))
      poly%patch%recruit_pool(:, 1) = 2.0_wp * cfg%demography%min_recruit_size
      patch_id0  = poly%patch%global_id(1:poly%patch%n)
      cohort_id0 = poly%cohort%global_id(1:poly%cohort%n)

      !----- The boundary: every phase runs and closes, and each operator declared what it moved. !
      led%active = .true.
      call slow_ledger_open(led, poly, cfg)
      call advance_boundary_dynamics(poly, cfg, .true., .true., led)
      do k = 1_ik, size(AT_BOUNDARY, kind=ik)
         call check(led%n_mark(AT_BOUNDARY(k)) > 0_ik, 'the boundary''s '//trim(BOUNDARY_NAME(k))//' phase ran')
         call check(slow_ledger_failures(led, AT_BOUNDARY(k)) == 0_ik,                                     &
                    'the boundary''s '//trim(BOUNDARY_NAME(k))//' phase closes on carbon, water and energy')
      end do
      call check(led%decl_gross(SLOW_PHASE_RECRUIT)%carbon > 0.0_wp, 'a recruit''s carbon is declared')
      call check(led%decl_gross(SLOW_PHASE_COHORT)%carbon > 0.0_wp, 'the cull''s dead carbon is declared')
      call check(led%decl_gross(SLOW_PHASE_DISTURB)%carbon > 0.0_wp, 'the treefall''s dead carbon is declared')

      !----- Persistent ids: distinct, and a new one is never one used before. ----------------!
      call check(all_distinct(poly%patch%global_id(1:poly%patch%n)), 'every patch has a distinct id')
      call check(all_distinct(poly%cohort%global_id(1:poly%cohort%n)), 'every cohort has a distinct id')
      call check(new_ids_are_new(poly%patch%global_id(1:poly%patch%n), patch_id0),                      &
                 'a gap takes an id never used before')
      call check(new_ids_are_new(poly%cohort%global_id(1:poly%cohort%n), cohort_id0),                   &
                 'a recruit or a split takes an id never used before')
      call check(.not. any(poly%cohort%global_id(1:poly%cohort%n) == culled_id), 'the culled cohort''s id is gone')

      !----- The step after: the boundary's dead carbon reaches the soil, and its phases close. !
      litter0 = litter_carbon_total(poly, poly%patch%litter_in)
      call check(litter0 > 0.0_wp, 'the boundary''s dead carbon waits in the litter')
      soil0 = soil_carbon_total(poly)
      call advance_slow_dynamics(poly, cfg, led)
      call check(soil_carbon_total(poly) - soil0 >= litter0 * (1.0_wp - 1.0e-12_wp),                    &
                 'the next soil step takes the boundary''s dead carbon up')
      call check(litter_carbon_total(poly, poly%patch%litter_in) == 0.0_wp,                             &
                 'the litter accumulator ends the step empty')
      do k = 1_ik, size(IN_STEP, kind=ik)
         call check(slow_ledger_failures(led, IN_STEP(k)) == 0_ik,                                         &
                    'the step''s '//trim(STEP_NAME(k))//' phase closes on carbon, water and energy')
      end do
   end subroutine run_restructure_ledger

   logical function all_distinct(ids)
      integer(ik), intent(in) :: ids(:)
      integer(ik) :: i
      all_distinct = .true.
      do i = 2_ik, size(ids, kind=ik)
         if (any(ids(1:i - 1) == ids(i))) all_distinct = .false.
      end do
   end function all_distinct

   !----- Every id not in `old` exceeds all of them: ids are never reused. -----------------!
   logical function new_ids_are_new(ids, old)
      integer(ik), intent(in) :: ids(:), old(:)
      integer(ik) :: i
      new_ids_are_new = .true.
      do i = 1_ik, size(ids, kind=ik)
         if (.not. any(old == ids(i)) .and. ids(i) <= maxval(old)) new_ids_are_new = .false.
      end do
   end function new_ids_are_new

   !----- The CENTURY pools over the polygon [kgC/m2]. --------------------------------------!
   real(wp) function soil_carbon_total(poly)
      type(polygon_t), intent(in) :: poly
      real(wp)    :: x(N_SOIL_POOL)
      integer(ik) :: ip
      soil_carbon_total = 0.0_wp
      do ip = 1_ik, poly%patch%n
         x = poly%patch%soil_carbon(ip)%carbon
         soil_carbon_total = soil_carbon_total + patch_polygon_area(poly, ip) * sum(x)
      end do
   end function soil_carbon_total

end module test_slow_update_ledger

!==========================================================================================!
! test_slow_update -- the operators of ed_slow_update, the one writer of the stand's       !
! carbon and structure (docs/science/demography.md §2-§4): the commit, cohort fusion and   !
! fission, the treefall gaps, patch fusion and termination, a recruit's birth state, and   !
! the restructuring through the slow ledger with soil carbon on. One module per section.   !
!==========================================================================================!
program test_slow_update
   use test_slow_update_commit, only : run_commit
   use test_slow_update_cohorts, only : run_cohort_operators
   use test_slow_update_gaps, only : run_treefall_gaps
   use test_slow_update_patches, only : run_patch_operators
   use test_slow_update_birth, only : run_birth_state
   use test_slow_update_ledger, only : run_restructure_ledger
   use ed_test_support, only : test_report
   implicit none

   call run_commit()
   call run_cohort_operators()
   call run_treefall_gaps()
   call run_patch_operators()
   call run_birth_state()
   call run_restructure_ledger()

   call test_report('test_slow_update')
end program test_slow_update
