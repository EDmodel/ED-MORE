! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_leaf_gas_exchange -- one leaf's gas exchange (meds_leaf_gas_exchange;               !
! docs/science/leaf_gas_exchange.md), on the in-code config fixture's leaf table.          !
!   1. C3 DEMAND (assimilation_demand_c3, eq. (2)): a hand-computed Ac, Aj and Ap, and    !
!      the co-limitation curvatures' cost against the sharp minimum.                       !
!   2. REGIMES (assimilation_demand_c3): Rubisco, RuBP and product limitation where each   !
!      should bind.                                                                        !
!   3. THE SOLVE (solve_leaf_gas_exchange, §5): a midday C3 leaf converges inside its      !
!      bounds, and its gs is the Medlyn law's at the solution.                             !
!   4. STOMATAL MODELS (§3): Leuning, Medlyn and Katul converge, and gs falls with VPD.    !
!   5. C4 (Collatz, §2): it draws Ci lower than C3; its PEP limb rises with temperature;   !
!      transpiration crosses gs and gb in series.                                          !
!   6. WATER STRESS (§4, eqs. (9)-(10)): the fast limb on the leaf potential and the slow  !
!      limb on the predawn potential, for every stomatal model; under Katul, every solve,  !
!      open or pinned at g0, lies on the CO2 diffusion line (eq. (1)); Katul behind a      !
!      thin boundary layer and at the dew point.                                           !
!   7. LIGHT: the night branch at PAR = 0, and a monotone light response.                  !
!   8. O2 (eq. (3)): oxygen reaches Gamma* as well as the Michaelis term.                  !
!   9. PRESSURE (rubisco_kinetics, eq. (3)): Gamma* as a mole fraction does not depend on  !
!      pressure; Kc and Ko are partial pressures.                                          !
!  10. A COHORT'S LEAF PAR (cohort_leaf_par): its absorbed PAR per unit of its own leaf    !
!      area, the same for a thin cohort as a thick one, none for a leafless one.           !
! The temperature responses are test_therm_lib's.                                          !
!==========================================================================================!
program test_leaf_gas_exchange
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_kelvin
   use meds_config,             only : meds_config_t, derive_parameters
   use meds_pft_derive,         only : pft_leaf_psi_tlp, pft_stomata_psi_onset, pft_wstress_psi_open,  &
                                       pft_wstress_psi_close
   use meds_pft_params,         only : HYD_UNSET, leaf_photo_pft_t
   use meds_leaf_params,        only : SM_LEUNING, SM_MEDLYN, SM_KATUL, COLIM_MIN, COLIM_QUADRATIC
   use meds_leaf_gas_exchange,  only : assimilation_demand_c3, rubisco_kinetics, stomata_gs_medlyn,          &
                                       solve_leaf_gas_exchange, leaf_env_t, leaf_out_t, LIM_NONE,         &
                                       LIM_RUBISCO, LIM_RUBP, LIM_C4_PEP, cohort_leaf_par
   use meds_test_support,       only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t) :: cfg
   type(leaf_env_t)    :: env
   type(leaf_out_t)   :: flux, flux2
   real(wp) :: a_gross, ac, aj, ap, an0, an1, an2, prev
   real(wp) :: psi_tlp, jump, gs_top
   integer(ik) :: i, j, n_pinned
   integer(ik), dimension(3) :: sms = [ SM_LEUNING, SM_MEDLYN, SM_KATUL ]

   call banner('leaf gas exchange (photosynthesis + stomata)')
   cfg = build_test_config()

   !=== 1. C3 demand: hand-computed Ac/Aj known-answer (vcmax=100, j=120, ci=300, gstar=40, ==!
   !       kc=400, ko=275000, o2=210000, tpu=15). Ac=100*260/(300+400*(1+210000/275000))      !
   !       = 25.8593; Aj=120*260/(4*300+8*40)=20.5263; Ap=3*15=45. =============================!
   call assimilation_demand_c3(300.0_wp, 100.0_wp, 120.0_wp, 15.0_wp, 40.0_wp, 400.0_wp, 275000.0_wp, &
                        210000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check_close(ac, 25.85930_wp, 1.0e-4_wp, 'C3 Ac known-answer')
   call check_close(aj, 20.52632_wp, 1.0e-4_wp, 'C3 Aj known-answer')
   call check_close(ap, 45.0_wp,     1.0e-6_wp, 'C3 Ap known-answer')
   call check_close(a_gross, 20.52632_wp, 1.0e-4_wp, 'C3 min co-limitation picks Aj')

   !=== 1b. CO-LIMITATION strength. The two C3 smoothings must use the CO-LIMITATION ========!
   !        curvatures, not theta_j (docs/science/leaf_gas_exchange.md §2). This is a REGIME     !
   !        assertion: at ambient CO2 Ac and Aj sit close together, so the smoothing penalty     !
   !        lands exactly where the model spends most of its time, and borrowing theta_j = 0.85  !
   !        costs ~29 % of assimilation against the sharp min(Ac,Aj,Ap). Measured here at PFT-1  !
   !        kinetics, Ci = 280, saturating light.                                                !
   block
      real(wp) :: a_min, a_new, a_old
      real(wp), parameter :: CI = 280.0_wp, VC = 60.0_wp, JR = 108.0_wp, TP = 15.0_wp
      real(wp), parameter :: GS = 42.75_wp, KC = 404.9_wp, KO = 278400.0_wp, O2 = 209000.0_wp
      call assimilation_demand_c3(CI, VC, JR, TP, GS, KC, KO, O2, COLIM_MIN,                     &
                           0.98_wp, 0.95_wp, a_min, ac, aj, ap)
      !----- The shipped curvatures. --------------------------------------------------------!
      call assimilation_demand_c3(CI, VC, JR, TP, GS, KC, KO, O2, COLIM_QUADRATIC,               &
                           0.98_wp, 0.95_wp, a_new, ac, aj, ap)
      !----- What borrowing theta_j would do: the size of the error it makes. ----------------!
      call assimilation_demand_c3(CI, VC, JR, TP, GS, KC, KO, O2, COLIM_QUADRATIC,               &
                           0.85_wp, 0.85_wp, a_old, ac, aj, ap)
      call check(a_old < 0.80_wp * a_min,                                                        &
                 'theta_j = 0.85 as a co-limitation curvature must cost > 20 % vs min()')
      call check(a_new > 0.88_wp * a_min,                                                        &
                 'co-limitation curvatures must keep the smoothed rate within 12 % of min()')
      call check(a_new < a_min,                                                                  &
                 'smoothing must still sit below the sharp min(), not above it')
      print '(a,3(f8.4,a))', '   [co-limitation] min(Ac,Aj,Ap) = ', a_min, ' | theta_*_c3 = ', a_new,     &
            ' | theta_j borrowed = ', a_old, ''
   end block

   !=== 2. C3 limitation regimes (raw demand, COLIM_MIN). ===================================!
   !----- High light, low Ci -> Rubisco-limited. ---------------------------------------------!
   call assimilation_demand_c3(120.0_wp, 60.0_wp, 104.0_wp, 10.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,  &
                        209000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check(ac < aj .and. ac < ap, 'low Ci, high light should be Rubisco-limited')
   !----- High Ci, low light (small J) -> RuBP/light-limited. ---------------------------------!
   call assimilation_demand_c3(300.0_wp, 60.0_wp, 20.0_wp, 10.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,   &
                        209000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check(aj < ac .and. aj < ap, 'high Ci, low light should be RuBP-limited')
   !----- High Ci, high light, low TPU -> product-limited. ------------------------------------!
   call assimilation_demand_c3(600.0_wp, 60.0_wp, 104.0_wp, 2.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,   &
                        209000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check(ap < ac .and. ap < aj, 'high Ci, low TPU should be product-limited')

   !=== 3. Full solve (PFT 1, C3, Medlyn default): convergence, bounds, the Medlyn law. =====!
   env = std_env()
   call leaf_solve(env, 1_ik, flux)
   call check(flux%converged, 'C3 solve must converge')
   call check(flux%A_net > 0.0_wp, 'C3 net assimilation must be positive at midday')
   call check(flux%gs >= cfg%pft%stomatal_g0(1), 'gs must be at least the cuticular g0')
   call check(flux%ci > 0.0_wp .and. flux%ci < flux%cs, 'Ci must lie between 0 and Cs')
   call check(flux%limitation == LIM_RUBISCO .or. flux%limitation == LIM_RUBP,                 &
              'midday C3 leaf should be Rubisco- or RuBP-limited')
   !----- The converged gs reproduces the Medlyn law (root really solved the coupled system). -!
   call check_close(flux%gs, stomata_gs_medlyn(flux%A_net, flux%cs, env%vpd,                   &
                    cfg%pft%stomatal_g0(1), cfg%pft%stomatal_g1(1), cfg%leaf_physiology%vpd_min), 1.0e-4_wp, &
                    'back-computed gs must match the Medlyn model at the solution')

   !=== 4. All three stomatal models converge; gs decreases as VPD rises. ===================!
   do i = 1_ik, 3_ik
      cfg%leaf_physiology%stomatal_model = sms(i)
      call derive_parameters(cfg)
      env = std_env()
      env%vpd = 800.0_wp
      call leaf_solve(env, 1_ik, flux)
      env%vpd = 2500.0_wp
      call leaf_solve(env, 1_ik, flux2)
      call check(flux%converged .and. flux2%converged, 'each stomatal model must converge')
      call check(flux%gs >= cfg%pft%stomatal_g0(1) .and. flux%A_net > 0.0_wp,                  &
                 'each stomatal model must give a positive, open-stomata solution')
      call check(flux2%gs < flux%gs, 'gs must decrease as VPD increases')
   end do
   cfg%leaf_physiology%stomatal_model = SM_MEDLYN
   call derive_parameters(cfg)

   !=== 5. C3 vs C4 contrast (PFT 1 vs PFT 3): C4 runs and draws Ci lower (concentrating). ==!
   env = std_env()
   call leaf_solve(env, 1_ik, flux)    ! C3
   call leaf_solve(env, 3_ik, flux2)   ! C4
   call check(flux2%converged .and. flux2%A_net > 0.0_wp, 'C4 leaf must converge and assimilate')
   call check(flux2%limitation == LIM_RUBISCO .or. flux2%limitation == LIM_RUBP .or.           &
              flux2%limitation == LIM_C4_PEP, 'C4 limitation flag must be a C4-valid regime')
   call check(flux2%ci / flux2%cs < flux%ci / flux%cs, 'C4 should operate at a lower Ci/Cs than C3')

   !=== 5b. C4 kp temperature response: at a CO2/PEP-limiting Ci the PEP-limited gross rate     !
   !     must RISE with leaf temperature, which it cannot with kp held at its 25 degC value. =====!
   block
      real(wp) :: a_cold, a_warm
      env = std_env()
      env%ca = 40.0_wp                                       ! very low CO2 -> C4 PEP-limited
      env%leaf_temp = t_kelvin + 18.0_wp
      call leaf_solve(env, 3_ik, flux)           ! C4 (PFT 3), cold
      a_cold = flux%A_gross
      env%leaf_temp = t_kelvin + 34.0_wp
      call leaf_solve(env, 3_ik, flux2)          ! C4 (PFT 3), warm
      a_warm = flux2%A_gross
      call check(flux%converged .and. flux2%converged, 'C4 low-CO2 solve must converge at both temperatures')
      call check(flux%limitation == LIM_C4_PEP, 'low-CO2 C4 leaf should be PEP-limited (isolates kp)')
      call check(a_warm > a_cold, 'C4 PEP-limited gross assimilation must rise with temperature (kp now temp-scaled)')
   end block

   !=== 5c. Transpiration puts stomata gs and the boundary layer gb in SERIES when use_boundary_layer is on   !
   !     (gs alone would overestimate the water flux). ===========================================!
   block
      real(wp) :: e_series
      cfg%leaf_physiology%use_boundary_layer = .true.
      call derive_parameters(cfg)
      env = std_env()
      env%gb = 0.6_wp                                        ! finite water-vapour boundary conductance [mol/m2/s]
      call leaf_solve(env, 1_ik, flux)
      e_series = flux%gs * env%gb / (flux%gs + env%gb) * env%vpd / env%pressure
      call check_close(flux%transpiration, e_series, 1.0e-9_wp, 'transpiration uses the gs-gb series conductance')
      call check(flux%transpiration < flux%gs * env%vpd / env%pressure,                            &
                 'series transpiration is below the gs-only value (gb resistance applied)')
      cfg%leaf_physiology%use_boundary_layer = .false.
      call derive_parameters(cfg)
   end block

   !=== 6-. The fast limb is idle while the leaf holds turgor: psi_leaf above the turgor-loss   !
   !     point is inert. ====================================================================!
   psi_tlp = pft_leaf_psi_tlp(cfg%pft, 1_ik)
   env = std_env()
   env%psi_leaf =  0.0_wp
   call leaf_solve(env, 1_ik, flux)
   an0 = flux%A_net
   env%psi_leaf = 0.9_wp * psi_tlp
   call leaf_solve(env, 1_ik, flux)
   an1 = flux%A_net
   call check_close(an1, an0, 1.0e-12_wp, 'psi_leaf above the turgor-loss point is inert')

   !=== 6. Water stress (PFT 1, Medlyn): A_net falls monotonically as the leaf wilts past its      !
   !     turgor-loss point, to -Rd once the fast limb has shut the biochemistry. ===================!
   env = std_env()
   env%psi_leaf =  0.0_wp
   call leaf_solve(env, 1_ik, flux)
   an0 = flux%A_net
   env%psi_leaf = 1.5_wp * psi_tlp
   call leaf_solve(env, 1_ik, flux)
   an1 = flux%A_net
   env%psi_leaf = 2.5_wp * psi_tlp
   call leaf_solve(env, 1_ik, flux)
   an2 = flux%A_net
   call check(an0 > an1 .and. an1 > an2, 'A_net must decrease as the leaf dries')
   call check_close(an2, -flux%rd, 1.0e-6_wp, 'fully stressed C3 leaf nets -Rd (beta_nonstomata = 0)')

   !=== 6a. The fast limb on a C4 leaf (PFT 3): it scales Vcmax, the light slope and kp, so a      !
   !     light-limited leaf at beta_nonstomata = 0.5 keeps about half its gross rate. =============!
   block
      real(wp) :: psi_tlp_c4, ag_wet
      psi_tlp_c4 = pft_leaf_psi_tlp(cfg%pft, 3_ik)
      env = std_env()
      env%par = 200.0_wp
      call leaf_solve(env, 3_ik, flux)
      ag_wet = flux%A_gross
      env%psi_leaf = 1.5_wp * psi_tlp_c4
      call leaf_solve(env, 3_ik, flux)
      call check_close(flux%beta_nonstomata, 0.5_wp, 1.0e-12_wp, 'C4 fast limb: 0.5 half way')
      call check(abs(flux%A_gross / ag_wet - 0.5_wp) < 0.05_wp,                                        &
                 'a light-limited C4 leaf at beta_nonstomata 0.5 keeps about half its gross rate')
   end block

   !=== 6b. Katul under full water stress must CLOSE (g0 fallback), not return an open flux. =!
   cfg%leaf_physiology%stomatal_model = SM_KATUL
   call derive_parameters(cfg)
   env = std_env()
   env%psi_leaf =  0.0_wp
   call leaf_solve(env, 1_ik, flux)
   an0 = flux%A_net
   env%psi_leaf = -5.0_wp
   call leaf_solve(env, 1_ik, flux2)
   call check(flux2%converged, 'Katul full-stress solve must converge via the g0 fallback')
   call check(flux2%A_net < an0, 'Katul full water stress must reduce assimilation')
   call check_close(flux2%gs, cfg%pft%stomatal_g0(1), 1.0e-2_wp, 'Katul full-stress leaf must close to ~g0')

   !=== 6c. Under water stress the Katul optimum can fall below the cuticular floor g0; the solve   !
   !     then finds Ci again with gs pinned at g0, so A, gs and Ci stay on the CO2 diffusion line   !
   !     A_net = gs (Cs - Ci) / 1.6 (docs/science/leaf_gas_exchange.md eq. (1)). An open solve is on !
   !     it by construction, its gs being back-computed from it; a pinned solve with the leaf still  !
   !     assimilating is on it only through that second solve, to the Ci solve's tolerance. The    !
   !     sweep, fine enough to reach the pinned case, must find every solve within 0.1 umol/mol of !
   !     Ci of the line. ==============================================================================!
   n_pinned = 0_ik
   do i = 0_ik, 2000_ik
      env = std_env()
      env%psi_leaf = -0.002_wp * real(i, wp)                                ! 0 .. -4 MPa
      call leaf_solve(env, 1_ik, flux)
      call check(flux%converged, 'Katul stress sweep must converge')
      if (flux%gs <= cfg%pft%stomatal_g0(1) * (1.0_wp + 1.0e-6_wp) .and. flux%A_net > 0.0_wp)               &
         n_pinned = n_pinned + 1_ik
      call check(abs(flux%gs * (flux%cs - flux%ci) / 1.6_wp - flux%A_net) <= 0.1_wp * flux%gs / 1.6_wp,    &
                 'Katul: every solve, open or pinned at g0, lies on the CO2 diffusion line')
   end do
   print '(a,i0)', '   [Katul] solves pinned at g0 with the leaf assimilating: ', n_pinned
   call check(n_pinned > 0_ik, 'the Katul stress sweep reaches the g0-pinned solve with the leaf assimilating')

   !=== 6c2. Katul behind a thin boundary layer and at the dew point. With gb = 0.2 mol/m2/s the    !
   !     bracket ends where Ci = Cs, so a lit leaf keeps an open solve instead                        !
   !     of the cuticular fallback; at VPD = 0 the conductance is the floored VPD's, as at 1 Pa. ======!
   block
      logical  :: bl_saved
      real(wp) :: gs_dew
      bl_saved = cfg%leaf_physiology%use_boundary_layer
      cfg%leaf_physiology%use_boundary_layer = .true.
      call derive_parameters(cfg)
      env = std_env()
      env%gb = 0.2_wp
      call leaf_solve(env, 1_ik, flux)
      call check(flux%converged .and. flux%gs > 2.0_wp * cfg%pft%stomatal_g0(1),                        &
                 'Katul behind a thin boundary layer keeps its stomata open')
      cfg%leaf_physiology%use_boundary_layer = bl_saved
      call derive_parameters(cfg)
      env = std_env()
      env%vpd = 0.0_wp
      call leaf_solve(env, 1_ik, flux)
      gs_dew = flux%gs
      env%vpd = 1.0_wp
      call leaf_solve(env, 1_ik, flux)
      call check(flux%A_net > 0.0_wp .and. gs_dew > 2.0_wp * cfg%pft%stomatal_g0(1),                    &
                 'Katul at the dew point keeps its stomata open')
      call check_close(gs_dew, flux%gs, 1.0e-12_wp, 'Katul gs is continuous at the dew point')
   end block
   cfg%leaf_physiology%stomatal_model = SM_MEDLYN
   call derive_parameters(cfg)

   !=== 6d. Stomatal-limb stress (Sabot beta_stomata via psi): with psi_leaf = 0 the fast limb  !
   !     is idle, so a drop in the predawn potential must close stomata (gs falls) through g1. ===!
   env = std_env()
   env%psi_leaf = 0.0_wp
   env%psi =  0.0_wp
   call leaf_solve(env, 1_ik, flux)
   an0 = flux%gs
   env%psi = -1.0_wp
   call leaf_solve(env, 1_ik, flux)
   an1 = flux%gs
   env%psi = -3.0_wp
   call leaf_solve(env, 1_ik, flux2)
   an2 = flux2%gs
   call check(flux2%converged, 'stomatal-limb (psi) solve must converge')
   call check(an0 > an1 .and. an1 > an2, 'gs must fall as soil water potential drops (beta_stomata on g1)')

   !=== 6e. The fast limb (beta_nonstomata) on the leaf's current potential: 1 down to the     !
   !     turgor-loss point, falling linearly to 0 at twice it (the defaults; the PFT file may set  !
   !     wstress_psi_open/_close). Every stomatal model ties gs to A, so as the leaf wilts gs falls !
   !     continuously to the cuticular g0, and a wilted leaf nets -Rd by day. ======================!
   call check_close(pft_wstress_psi_open(cfg%pft, 1_ik),  psi_tlp,          1.0e-14_wp,                 &
                    'fast limb starts at the turgor-loss point')
   call check_close(pft_wstress_psi_close(cfg%pft, 1_ik), 2.0_wp * psi_tlp, 1.0e-14_wp,                 &
                    'fast limb completes at twice the turgor-loss point')
   env = std_env()
   env%psi_leaf = 0.9_wp * psi_tlp
   call leaf_solve(env, 1_ik, flux)
   call check_close(flux%beta_nonstomata, 1.0_wp, 0.0_wp, 'fast limb: 1 while the leaf holds turgor')
   env%psi_leaf = 1.5_wp * psi_tlp
   call leaf_solve(env, 1_ik, flux)
   call check_close(flux%beta_nonstomata, 0.5_wp, 1.0e-12_wp, 'fast limb: 0.5 half way')
   env%psi_leaf = 2.0_wp * psi_tlp
   call leaf_solve(env, 1_ik, flux)
   call check_close(flux%beta_nonstomata, 0.0_wp, 0.0_wp, 'fast limb: 0 at twice the turgor-loss point')
   do j = 1_ik, size(sms)
      cfg%leaf_physiology%stomatal_model = sms(j)
      call derive_parameters(cfg)
      env = std_env()
      env%psi_leaf = psi_tlp
      call leaf_solve(env, 1_ik, flux)
      gs_top = flux%gs
      jump = 0.0_wp
      prev = huge(1.0_wp)
      do i = 0_ik, 280_ik
         env%psi_leaf = psi_tlp * (0.8_wp + 1.4_wp * real(i, wp) / 280.0_wp)
         call leaf_solve(env, 1_ik, flux)
         call check(flux%converged .and. flux%gs >= cfg%pft%stomatal_g0(1) * (1.0_wp - 1.0e-12_wp),   &
                    'fast-limb sweep: converged, gs >= g0')
         !----- to the Ci solve's tolerance (1e-3 ppm): near the g0 floor Cs - Ci is ~30 ppm, so the  !
         !      back-computed gs carries ~1e-5 relative noise (Katul switches to its g0 pass there). !
         if (i > 0_ik) call check(flux%gs <= prev * (1.0_wp + 1.0e-4_wp) + 1.0e-12_wp,            &
                                  'fast-limb sweep: gs never rises as the leaf dries')
         if (i > 0_ik) jump = max(jump, prev - flux%gs)
         prev = flux%gs
      end do
      call check(jump < 0.05_wp * gs_top, 'fast-limb sweep: no step larger than 5% of gs at the turgor-loss point')
      env%psi_leaf = 2.1_wp * psi_tlp
      call leaf_solve(env, 1_ik, flux)
      call check_close(flux%gs, cfg%pft%stomatal_g0(1), 1.0e-9_wp, 'wilted leaf by day: gs is the cuticular g0')
      call check_close(flux%A_gross, 0.0_wp, 1.0e-12_wp, 'wilted leaf by day: no photosynthesis')
      call check_close(flux%A_net, -flux%rd, 1.0e-12_wp, 'wilted leaf by day: the leaf respires')
   end do
   cfg%leaf_physiology%stomatal_model = SM_MEDLYN
   call derive_parameters(cfg)

   !=== 6f. The stomatal stress begins at an onset: half the turgor-loss point unless the PFT     !
   !     file sets stomata_psi_onset. Above it beta_stomata is 1, so the gravity head a tall tree   !
   !     has in wet soil does not close its stomata; below it beta = exp(sref*(psi - onset)). =====!
   psi_tlp = pft_leaf_psi_tlp(cfg%pft, 1_ik)
   call check_close(pft_stomata_psi_onset(cfg%pft, 1_ik), 0.5_wp * psi_tlp, 1.0e-14_wp,                 &
                    'onset defaults to half the turgor-loss point')
   env = std_env()
   env%psi_leaf = 0.0_wp
   env%psi = 0.0_wp
   call leaf_solve(env, 1_ik, flux)
   gs_top = flux%gs
   env%psi = 0.99_wp * pft_stomata_psi_onset(cfg%pft, 1_ik)
   call leaf_solve(env, 1_ik, flux)
   call check(flux%beta_stomata == 1.0_wp .and. flux%gs == gs_top, 'onset: no stress just above it')
   env%psi = pft_stomata_psi_onset(cfg%pft, 1_ik) - 0.5_wp
   call leaf_solve(env, 1_ik, flux)
   call check_close(flux%beta_stomata, exp(-0.5_wp * cfg%pft%wstress_sref_stomata(1)), 1.0e-14_wp,     &
                    'onset: beta = exp(sref*(psi - onset)) below it')
   cfg%pft%stomata_psi_onset(1) = -0.3_wp
   call derive_parameters(cfg)
   call check_close(pft_stomata_psi_onset(cfg%pft, 1_ik), -0.3_wp, 0.0_wp, 'onset: a set value is used')
   env%psi = -0.5_wp
   call leaf_solve(env, 1_ik, flux)
   call check_close(flux%beta_stomata, exp(-0.2_wp * cfg%pft%wstress_sref_stomata(1)), 1.0e-14_wp,     &
                    'onset: a set value moves the threshold')
   cfg%pft%stomata_psi_onset(1) = HYD_UNSET
   call derive_parameters(cfg)

   !=== 7. PAR sweep: night/closed branch at PAR=0, monotone rise, no NaNs, all converge. ===!
   env = std_env()
   env%par = 0.0_wp
   call leaf_solve(env, 1_ik, flux)
   call check(flux%limitation == LIM_NONE, 'PAR=0 must hit the closed/night branch')
   call check_close(flux%A_net, -flux%rd, 1.0e-6_wp, 'night A_net must equal -Rd')
   call check_close(flux%gs, cfg%pft%stomatal_g0(1), 1.0e-12_wp, 'night gs must be g0')
   prev = -1.0e30_wp
   do i = 0_ik, 20_ik
      env%par = 100.0_wp * real(i, wp)
      call leaf_solve(env, 1_ik, flux)
      call check(flux%converged, 'every PAR level must converge')
      call check(flux%A_net == flux%A_net, 'A_net must never be NaN')   ! NaN /= itself
      call check(flux%A_net >= prev - 1.0e-9_wp, 'A_net must rise monotonically with PAR')
      prev = flux%A_net
   end do

   !=== 8. O2 propagates to BOTH places oxygen enters the C3 demand. =======================!
   !      Gamma* = 0.5*O/S_(c/o) is proportional to the O2 partial pressure
   !      (docs/science/leaf_gas_exchange.md eq. (3)). If o2_mol_frac reached only the Michaelis
   !      term Kc(1+O/Ko), raising O2 would inhibit carboxylation while leaving the
   !      photorespiratory penalty on Aj untouched.
   env = std_env()
   cfg%leaf_physiology%o2_mol_frac = 0.209_wp
   call derive_parameters(cfg)
   call leaf_solve(env, 1_ik, flux)
   an0 = flux%A_net                                           ! the calibration point
   cfg%leaf_physiology%o2_mol_frac = 0.105_wp
   call derive_parameters(cfg)
   call leaf_solve(env, 1_ik, flux)
   an1 = flux%A_net                                           ! half O2
   cfg%leaf_physiology%o2_mol_frac = 0.350_wp
   call derive_parameters(cfg)
   call leaf_solve(env, 1_ik, flux2)
   an2 = flux2%A_net                                          ! paleo-high O2
   call check(flux%converged .and. flux2%converged, 'the O2 sweep must converge')
   !----- Less O2 => less photorespiration => MORE assimilation, and vice versa. Monotone in O2. -!
   call check(an1 > an0, 'halving O2 must RAISE A_net (less photorespiration)')
   call check(an2 < an0, 'raising O2 to 35% must LOWER A_net (more photorespiration)')
   !----- The response must run through Gamma*, NOT just the Michaelis term. This is the assertion -!
   !      a frozen Gamma* fails, and the thresholds are set from measuring BOTH on this exact       !
   !      fixture rather than guessed:                                                              !
   !                                                                                          !
   !                        O2 = 10.5%      O2 = 35%                                                !
   !        with Gamma*(O2)   +22.0 %       -23.7 %                                                  !
   !        frozen Gamma*      +8.2 %       -10.5 %   <- Kc(1+O/Ko) acting alone                      !
   !                                                                                          !
   !      A >5% gate would pass on BOTH; 15% separates them with room on either side.                 !
   call check((an1 - an0) / abs(an0) > 0.15_wp,                                                  &
              'halving O2 must raise A_net > 15% -- the response runs through Gamma*, not just Kc')
   call check((an2 - an0) / abs(an0) < -0.15_wp,                                                 &
              'raising O2 to 35% must cut A_net > 15% -- likewise')
   cfg%leaf_physiology%o2_mol_frac = 0.209_wp                 ! back to the shipped default
   call derive_parameters(cfg)

   !=== 9. Gamma* as a mole fraction does not depend on pressure (docs/science/leaf_gas_exchange.md!
   !       eq. (3)): it is proportional to the O2 partial pressure, which scales with pressure as   !
   !       the CO2 partial pressure does. At 25 C and 21 % O2 it is Bernacchi et al. (2001)'s       !
   !       42.75 umol/mol at sea level and at 3000 m (70 kPa) alike. Kc and Ko are partial          !
   !       pressures, so as mole fractions they scale as 1/P. ======================================!
   block
      type(leaf_photo_pft_t)   :: pp
      real(wp) :: kc_low, ko_low, gstar_low, o2_low, kc_high, ko_high, gstar_high, o2_high
      pp = cfg%pft%leaf_photo%pft(1)
      call rubisco_kinetics(pp, t_kelvin + 25.0_wp, 101325.0_wp, kc_low, ko_low, gstar_low, o2_low)
      call rubisco_kinetics(pp, t_kelvin + 25.0_wp, 70000.0_wp, kc_high, ko_high, gstar_high, o2_high)
      call check(gstar_high == gstar_low, 'Gamma* in umol/mol does not depend on pressure')
      call check_close(gstar_low, 42.75_wp, 1.0e-12_wp, 'Gamma* at 25 C and 21 % O2 is 42.75 umol/mol')
      call check(o2_high == o2_low, 'O2 in umol/mol does not depend on pressure')
      call check_close(kc_high * 70000.0_wp, kc_low * 101325.0_wp, 1.0e-12_wp, 'Kc is a partial pressure')
      call check_close(ko_high * 70000.0_wp, ko_low * 101325.0_wp, 1.0e-12_wp, 'Ko is a partial pressure')
      !----- The same Gamma* at 35 C, where the Arrhenius factor is not 1. -----------------------!
      call rubisco_kinetics(pp, t_kelvin + 35.0_wp, 101325.0_wp, kc_low, ko_low, gstar_low, o2_low)
      call rubisco_kinetics(pp, t_kelvin + 35.0_wp, 70000.0_wp, kc_high, ko_high, gstar_high, o2_high)
      call check(gstar_high == gstar_low, 'Gamma* in umol/mol does not depend on pressure at 35 C')
   end block

   !=== 10. A cohort's leaf PAR is its absorbed PAR per unit of its own leaf area: a thin cohort's  !
   !        leaves see the same light as a thick one's when each absorbs in proportion to its leaf  !
   !        area, and a leafless cohort sees none. 5 W/m2 absorbed by LAI 0.05 at 4.6 umol/J is     !
   !        460 umol/m2 leaf/s. =================================================================!
   call check_close(cohort_leaf_par(5.0_wp, 0.05_wp, 4.6_wp), 460.0_wp, 1.0e-12_wp,                    &
                    'the leaf PAR of a thin cohort (LAI 0.05) is its own')
   call check_close(cohort_leaf_par(5.0_wp, 0.05_wp, 4.6_wp), cohort_leaf_par(50.0_wp, 0.5_wp, 4.6_wp), 1.0e-12_wp, &
                    'the leaf PAR does not depend on the cohort''s thickness')
   call check(cohort_leaf_par(0.0_wp, 0.0_wp, 4.6_wp) == 0.0_wp, 'a leafless cohort absorbs nothing and sees no PAR')

   write(*,'(a)') '   PASS'

contains

   !----- A representative midday tropical leaf environment. ------------------------------!
   function std_env() result(e)
      type(leaf_env_t) :: e
      e%par = 1500.0_wp
      e%leaf_temp = t_kelvin + 25.0_wp
      e%vpd = 1500.0_wp
      e%ca = 400.0_wp
      e%pressure = 101325.0_wp
      e%psi_leaf = 0.0_wp
      e%gb = 0.0_wp
      e%psi = 0.0_wp
   end function std_env

   !----- One leaf's gas exchange for PFT ipft, as the fast loop solves it: the configuration's  !
   !      leaf table, which the test derives again after changing a key it reads. --------------!
   subroutine leaf_solve(e, ipft, f)
      type(leaf_env_t),  intent(in)  :: e
      integer(ik),       intent(in)  :: ipft
      type(leaf_out_t), intent(out) :: f
      associate (t => cfg%pft%leaf_photo)
         call solve_leaf_gas_exchange(e, t%pft(ipft), t%stomatal_model, t%temp_response_form, t%colimitation,  &
                                      t%use_boundary_layer, f)
      end associate
   end subroutine leaf_solve

end program test_leaf_gas_exchange
