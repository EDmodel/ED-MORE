! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_plant_energy -- a tissue's energy balance (ed_plant_energy: tissue_energy_balance;
! docs/science/vegetation_energy_dynamics.md eqs. (1)-(8)).
!   1. EXACT RELAXATION (eqs. (3)-(7)): the closed-form tissue step, its diagnostic and
!      frozen limits, its energy identity and small-x branch, marched exactly against the
!      analytic relaxation; a store with no coupling warms by its forcing.
!   2. LATENT PATHWAYS: the stomata's and the wetted film's; f_wet = 0 leaves no trace of
!      the film, f_wet = 1 is the film alone, and partial wetting closes the balance.
!   3. THIN TISSUE (eqs. (1), (8)): energy closes at every area index, halving the area
!      halves every flux and keeps the temperature, and a tissue with no area takes the
!      temperature of what it grows from.
!   4. THE CANOPY'S FILMS (canopy_tissue_energy, eq. (2)): energy is conserved across the
!      canopy, its films capturing, passing between leaf and wood and evaporating, with
!      heat stores and the plant water's heat; rain captured at a tissue's own temperature
!      brings it no heat and warmer rain warms it; a film passed from leaf to wood charges
!      the giver nothing and brings the receiver the giver's enthalpy; an evaporating film
!      costs its tissue the vapour's enthalpy at the canopy air less the liquid's at the
!      tissue's start-of-step temperature.
!   5. THE TISSUE'S LONGWAVE (canopy_tissue_energy; "The emission base"): linearized about
!      the tissue's own start-of-step temperature, it is the exact Stefan-Boltzmann
!      balance when the tissue ends there, and second order in the distance from it.
!   6. THE COMMIT (commit_tissue_energy; eqs. (9)-(10)): the energies gain exactly the
!      heat each tissue relaxed through and its water's enthalpy, the gained water arrives
!      at the tissue's start-of-step temperature, and the committed temperatures give back
!      the energies, so the next step starts where this one ended; a film far heavier than
!      its tissue that evaporates within the step leaves the tissue at the step's
!      temperature, its heat capacity over the step being the tissue's own.
!---~---
program test_plant_energy
   use ed_test_assert,    only : check, check_true, test_report
   use ed_kinds,          only : wp, ik
   use ed_constants,      only : cp_air, cp_liq, stefan_boltzmann, pi
   use ed_therm_lib,      only : enthalpy_vapor, internal_energy_liquid, latent_heat_vap, sat_specific_humidity,   &
                                 sat_specific_humidity_temp_deriv, tissue_temp_to_energy
   use ed_canopy_water,   only : canopy_water_out_t
   use ed_plant_energy,   only : tissue_energy_balance, latent_pathway_t, plant_energy_env_t, canopy_tissue_energy,  &
                                 commit_tissue_energy, tissue_hcap_per_dt
   implicit none

   !--- The canopy air the canopy tests run under.
   real(wp), parameter :: T_CAS = 298.0_wp, Q_CAS = 0.012_wp, P_AIR = 101325.0_wp, RHO_AIR = 1.2_wp
   !---~---
   !   The longwave test's leaf: its shortwave, its area, emissivity and sensible coupling, and the
   ! longwave it absorbs, two sides under 285 K surroundings.
   !---~---
   real(wp), parameter :: LEAF_SW = 300.0_wp, LEAF_LAI = 2.0_wp, LEAF_EMISS = 0.97_wp, LEAF_H = 40.0_wp
   real(wp), parameter :: LEAF_LW_IN = 2.0_wp * 0.97_wp * stefan_boltzmann * 285.0_wp ** 4

   call test_veg_exponential()
   call test_veg_energy_diagnostic_wetted()
   call test_veg_energy_thin_tissue()
   call test_canopy_energy_closes()
   call test_film_capture()
   call test_film_transfer()
   call test_film_evaporation()
   call test_tissue_longwave()
   call test_commit()
   call test_commit_heavy_film()

   call test_report('test_plant_energy')

contains

   !---~---
   !   EXACT EXPONENTIAL RELAXATION in tissue_energy_balance.
   !
   !   Under the Category-0 freeze the tissue equation is linear with tau = cap/denom, so the
   ! step has a closed form and the kernel uses it: endpoint weight exp(-x), flux weight
   ! (1-exp(-x))/x, x = denom/store_hcap_per_dt (docs/science/vegetation_energy_dynamics.md
   ! eqs. (3)-(7)). These checks pin the properties that make it usable -- the two limits, the
   ! conservation identity, the small-x branch, and (the one that matters) that marching it is
   ! EXACT, not merely convergent.
   !---~---
   subroutine test_veg_exponential()
      real(wp) :: dtt, ts, tr, dh, drn, dt_diag_ref, ts_ref
      real(wp) :: t_cas, t0, sw, lw, h, les, lws, ler, cap, dt, denom, tau, a_st
      real(wp) :: t_exact, t_march, worst, fe
      type(latent_pathway_t) :: dry, none   ! the stomata's pathway; none: no film
      integer(ik) :: k, nstep
      print '(a)', 'test_veg_exponential:'
      t_cas = 295.0_wp
      sw = 300.0_wp
      lw = -40.0_wp
      h = 120.0_wp
      les = 30.0_wp
      lws = 5.0_wp
      ler = 25.0_wp
      dry = latent_pathway_t(slope=les, ref=ler)
      denom = h + les + lws
      dt    = 1800.0_wp

      !---~---
      !   (a) store_hcap_per_dt -> 0 is the DIAGNOSTIC limit, exactly. This is why there is no
      ! diagnostic/prognostic selector: diagnostic is a limit of one formula, not a second
      ! branch.
      !---~---
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, 0.0_wp, t_cas, &
                                 dtt, ts, tr, fe, dh, drn)
      dt_diag_ref = (sw + lw - ler) / denom
      call check('store_hcap_per_dt = 0 gives the exact diagnostic offset', dtt, dt_diag_ref, 1.0e-12_wp)
      ts_ref = t_cas + dt_diag_ref

      !--- (b) store_hcap_per_dt -> huge FREEZES the store at its entry temperature.
      t0 = 288.0_wp
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, 1.0e12_wp, t0, &
                                 dtt, ts, tr, fe, dh, drn)
      call check('store_hcap_per_dt -> infinity holds the store at t_store0', ts, t0, 1.0e-6_wp)

      !---~---
      !   (c) THE CONSERVATION IDENTITY. The endpoint and step-average weights are different
      ! numbers, and pairing them correctly is exactly what makes the balance close:
      !   store_hcap_per_dt*(dt_end - dt_prev) + denom*dt_avg == numer.
      ! dt_avg is not returned, so recover it from dh = h*dt_avg.
      !---~---
      cap  = 4.0e4_wp                       ! a big cohort: tau = cap/denom is order dt
      a_st = cap / dt
      tau  = cap / denom
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, a_st, t0, &
                                 dtt, ts, tr, fe, dh, drn)
      call check('energy balance closes with the endpoint/average pair',                            &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-9_wp)

      !---~---
      !   (d) THE POINT: marching the kernel reproduces the analytic relaxation EXACTLY, not
      ! approximately. A backward-Euler storage term would drift here; this does not.
      ! T(t) = T_eq + (T0 - T_eq)*exp(-t/tau).
      !---~---
      nstep = 20_ik
      t_march = t0
      worst = 0.0_wp
      do k = 1_ik, nstep
         call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, a_st, t_march, &
                                    dtt, ts, tr, fe, dh, drn)
         t_march = ts
         t_exact = ts_ref + (t0 - ts_ref) * exp(-real(k, wp) * dt / tau)
         worst   = max(worst, abs(t_march - t_exact))
      end do
      call check_true('marching 20 steps is EXACT against the analytic relaxation', worst < 1.0e-10_wp, worst)

      !---~---
      !   (e) The small-x series branch (a very large capacity) must not lose precision where
      ! 1 - exp(-x) would cancel. tau = 1e6 s against dt = 1800 s puts x ~ 2.7e-4.
      !---~---
      a_st = 1.0e6_wp * denom / dt
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, a_st, t0, &
                                 dtt, ts, tr, fe, dh, drn)
      call check('large-capacity limit stays accurate (series branch)',                             &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-6_wp)

      !---~---
      !   (f) A store with NO coupling (x = 0) warms by its forcing over the step, numer*dt/cap,
      ! and passes nothing on: the relaxation needs no division by the coupling.
      !---~---
      a_st = cap / dt
      call tissue_energy_balance(sw, lw, 0.0_wp, 0.0_wp, 0.0_wp, latent_pathway_t(ref=ler), none, 0.0_wp, t_cas, &
                                 t_cas, a_st, t0, dtt, ts, tr, fe, dh, drn)
      call check('no coupling: the store warms by numer*dt/cap', ts, t0 + (sw + lw - ler) / a_st, 1.0e-12_wp)
      call check_true('no coupling: no sensible heat leaves', dh == 0.0_wp, dh)
   end subroutine test_veg_exponential

   subroutine test_veg_energy_diagnostic_wetted()
      !---~---
      !   Direct unit test of tissue_energy_balance's two latent pathways (the column tests exercise
      ! them only through the full integration). f_wet splits the latent flux into a (1-f_wet)
      ! DRY (stomatal) share and an f_wet WET (boundary-layer film) share; DRY and WET carry two
      ! DIFFERENT conductances, so the three properties below are non-trivial (a bug that swapped
      ! or mixed them would still "look" plausible).
      !---~---
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_temp2, t_store2, transp2, dh2, drnet2, film_evap2
      real(wp), parameter :: ABS_SW = 400.0_wp, ABS_LW = -50.0_wp, H_COEFF = 12.0_wp
      real(wp), parameter :: LW_SLOPE = 6.0_wp, T_CAS = 298.0_wp, T_PARENT = 298.0_wp
      type(latent_pathway_t) :: DRY, WET, none
      real(wp) :: h_vap
      print '(a)', 'test_veg_energy_diagnostic_wetted:'
      !--- Both pathways evaporate liquid at the canopy air's temperature.
      h_vap = latent_heat_vap(T_CAS, 1.0_wp)
      DRY   = latent_pathway_t(slope=8.0_wp,  ref=15.0_wp, h_vap=h_vap)   ! stomatal
      WET   = latent_pathway_t(slope=20.0_wp, ref=40.0_wp, h_vap=h_vap)   ! film

      !---~---
      !   1. f_wet = 0 leaves no trace of the film's pathway: the same result, bit for bit, as no
      ! film pathway at all.
      !---~---
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, none, 0.0_wp,           &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp, t_store, transp, film_evap, dh, drnet)
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, WET, 0.0_wp,            &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp2, t_store2, transp2, film_evap2, dh2, &
                                 drnet2)
      call check('f_wet=0 == no film pathway (dt_temp)', dt_temp2, dt_temp, 1.0e-14_wp)
      call check('f_wet=0 == no film pathway (transp)',  transp2,  transp,  1.0e-14_wp)
      call check_true('f_wet=0 -> film_evap = 0 exactly', film_evap2 == 0.0_wp, film_evap2)

      !---~---
      !   2. FULLY WET (f_wet=1): the dry (stomatal) pathway must vanish -- a saturated film
      ! transpires ~0 and evaporates its whole latent capacity through the wet pathway instead.
      ! SWAP-EQUIVALENCE: this must be algebraically IDENTICAL to calling the kernel with the film's
      ! pathway as the DRY one and no film -- not just approximately equal (both reduce to the same
      ! 2-unknown linear solve).
      !---~---
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, WET, 1.0_wp,            &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp, t_store, transp, film_evap, dh, drnet)
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, WET, none, 0.0_wp,           &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp2, t_store2, transp2, film_evap2, dh2, &
                                 drnet2)
      call check_true('f_wet=1: dry pathway (transp) vanishes', abs(transp) < 1.0e-14_wp, transp)
      call check('f_wet=1: film_evap == swap-equivalent transp', film_evap, transp2, 1.0e-12_wp)
      call check('f_wet=1: dt_temp == swap-equivalent dt_temp',  dt_temp,   dt_temp2, 1.0e-12_wp)

      !---~---
      !   3. PARTIAL wetting (f_wet=0.4): the energy balance must still close -- absorbed radiation
      ! net of LW emission (drnet) splits EXACTLY into sensible + BOTH latent pathways, each at the
      ! pathway's h_vap (here latent_heat_vap at the canopy-air temperature) -- the identity the surface path relies
      ! on to make the CAS's temperature-dependent vapour-enthalpy credit balance against a store's
      ! liquid-enthalpy debit (for transpiration the leaf pays the full vapour enthalpy itself; the
      ! surface-water accounting does the same for film_evap).
      !---~---
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, WET, 0.4_wp,            &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp, t_store, transp, film_evap, dh, drnet)
      call check('partial f_wet: energy balance closes (drnet = dh + latent)', drnet,                  &
                 dh + (transp + film_evap) * h_vap, 1.0e-9_wp)
      call check_true('partial f_wet: film_evap > 0 (wet pathway active)', film_evap > 0.0_wp, film_evap)
      call check_true('partial f_wet: transp > 0 (dry pathway still active)', transp > 0.0_wp, transp)
   end subroutine test_veg_energy_diagnostic_wetted


   subroutine test_veg_energy_thin_tissue()
      !---~---
      !   Every term of a tissue's balance -- absorbed radiation, the sapflow's enthalpy, the
      ! conductances -- and its heat capacity scale with the tissue's own area
      ! (docs/science/vegetation_energy_dynamics.md eqs. (1), (8)), so a thin tissue is simply a
      ! small one. Three properties, swept over an area index from 1e-6 to 10 for leaf and wood,
      ! with and without a heat store:
      !
      !   * CONSERVATION: within the step everything that enters leaves or is stored,
      !       drnet + q_extra == dh + Lv*(transp + film_evap) + store_hcap_per_dt*(t_store - t0),
      !     to round-off. surface_derivs relies on it: the whole-column ledger books coh_rnet on
      !     the input side and coh_h/coh_qw as they reach the CAS, and nothing re-checks the
      !     tissue in between.
      !   * SPLIT INVARIANCE: halving the area (one cohort of two at half density) keeps the
      !     temperature and halves every flux. A floor in ground units on the coupling or the
      !     capacity breaks this for a thin tissue, pulling it toward the canopy air.
      !   * NO AREA: a tissue with no coupling and no heat capacity takes the temperature of what
      !     it grows from (t_parent) and carries no flux.
      !---~---
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_h, t_h, transp_h, dh_h, drnet_h, film_h
      real(wp) :: ai, ast, resid, worst, worst_split, h_vap
      integer(ik) :: ia, il, is
      real(wp), parameter :: T_CAS = 298.0_wp, T0 = T_CAS - 0.5_wp, T_PARENT = 301.0_wp
      real(wp), parameter :: CAP_PER_AREA = 1500.0_wp   ! [J/m2 tissue/K] a leaf's dry tissue and water
      logical :: is_leaf
      print '(a)', 'test_veg_energy_thin_tissue:'
      h_vap = latent_heat_vap(T_CAS, 1.0_wp)   ! liquid at the canopy air's temperature
      worst = 0.0_wp
      worst_split = 0.0_wp
      do is = 0_ik, 1_ik
         do il = 1_ik, 2_ik
            is_leaf = (il == 1_ik)
            do ia = 1_ik, 15_ik
               ai = 1.0e-6_wp * (10.0_wp ** (0.5_wp * real(ia - 1_ik, wp)))
               ast = 0.0_wp
               if (is == 1_ik) ast = CAP_PER_AREA * ai / 900.0_wp
               call thin_tissue(ai, ast, is_leaf, dt_temp, t_store, transp, dh, drnet, film_evap)
               resid = drnet + 25.0_wp*ai - dh - h_vap*(transp + film_evap) - ast * (t_store - T0)
               worst = max(worst, abs(resid) / (abs(drnet) + 25.0_wp*ai))
               call thin_tissue(0.5_wp*ai, 0.5_wp*ast, is_leaf, dt_h, t_h, transp_h, dh_h, drnet_h, film_h)
               worst_split = max(worst_split, abs(t_h - t_store), abs(2.0_wp*dh_h - dh) / abs(dh),          &
                                 abs(2.0_wp*drnet_h - drnet) / abs(drnet))
            end do
         end do
      end do
      call check_true('a thin tissue conserves energy at every area index (60 cases)', worst < 1.0e-12_wp, worst)
      call check_true('half the area: the same temperature, half of every flux', worst_split < 1.0e-12_wp,       &
                      worst_split)
      !--- No area: no coupling and no heat capacity.
      call thin_tissue(0.0_wp, 0.0_wp, .true., dt_temp, t_store, transp, dh, drnet, film_evap)
      call check('a tissue with no area takes the temperature of what it grows from', t_store, T_PARENT, 0.0_wp)
      call check_true('a tissue with no area carries no flux',                                                  &
                      dh == 0.0_wp .and. drnet == 0.0_wp .and. transp == 0.0_wp .and. film_evap == 0.0_wp, dh)
   end subroutine test_veg_energy_thin_tissue

   !---~---
   !   One tissue of area index a at 298 K canopy air, starting 0.5 K below it, whose parent sits at
   ! 301 K: every coupling and forcing term proportional to a.
   !---~---
   subroutine thin_tissue(a, a_store, leaf, dtt, ts, tr, h_out, rn, fe)
      real(wp), intent(in)  :: a, a_store
      logical,  intent(in)  :: leaf
      real(wp), intent(out) :: dtt, ts, tr, h_out, rn, fe
      real(wp), parameter   :: GBH = 0.03_wp, RHO = 1.2_wp, EMISS = 0.97_wp
      real(wp), parameter   :: T_CAS = 298.0_wp, T0 = T_CAS - 0.5_wp, T_PARENT = 301.0_wp
      real(wp) :: hc, lws, les, ler, h_vap
      hc  = merge(2.0_wp, pi, leaf) * a * GBH * RHO * cp_air
      lws = 4.0_wp * EMISS * stefan_boltzmann * T_CAS**3 * a
      les = merge( 8.0_wp * a, 0.0_wp, leaf)    ! stomatal (wood does not transpire)
      ler = merge(15.0_wp * a, 0.0_wp, leaf)
      h_vap = latent_heat_vap(T_CAS, 1.0_wp)   ! liquid at the canopy air's temperature
      call tissue_energy_balance(400.0_wp*a, -50.0_wp*a, 25.0_wp*a, hc, lws, latent_pathway_t(slope=les, ref=ler, h_vap=h_vap), &
                                 latent_pathway_t(slope=20.0_wp*a, ref=40.0_wp*a, h_vap=h_vap), 0.4_wp, T_CAS, T_PARENT,      &
                                 a_store, T0, dtt, ts, tr, fe, h_out, rn)
   end subroutine thin_tissue

   !---~---
   !   4. THE CANOPY'S FILMS (canopy_tissue_energy). A film is part of its tissue's energy, at
   ! the tissue's start-of-step liquid enthalpy u0 (as the commit values its change); its
   ! capture brings the rain's liquid enthalpy, its evaporation takes the vapour's at the
   ! canopy air's temperature, and water passed between a cohort's two films leaves the giver
   ! at the giver's enthalpy. Then film, tissue and canopy air close with no slack
   ! (docs/science/vegetation_energy_dynamics.md eq. (2)).
   !---~---

   !---~---
   !   Two cohorts with everything on: radiation, sensible heat, heat stores, the plant water's
   ! heat, partly wet films capturing rain and passing water leaf to wood in one cohort and
   ! wood to leaf in the other, transpiration held at the network's or free. What enters (the
   ! net radiation, the plant water's heat, the captured rain's liquid enthalpy) less what
   ! leaves into the canopy air (sensible heat, the vapour's enthalpy) is what the tissues
   ! store: their heat, cap*(T - T0), and their films' change at u0.
   !---~---
   subroutine test_canopy_energy_closes()
      type(plant_energy_env_t) :: env
      type(canopy_water_out_t) :: film
      real(wp) :: worst, resid, scale
      integer(ik) :: mode
      print '(a)', 'test_canopy_energy_closes:'
      worst = 0.0_wp
      do mode = 1_ik, 2_ik
         call canopy_inputs(2_ik, env, film)
         env%h_coeff_leaf  = [30.0_wp, 20.0_wp]
         env%g_transp_leaf = [0.006_wp, 0.004_wp]
         env%abs_sw        = [300.0_wp, 150.0_wp]
         env%abs_lw        = [-40.0_wp, -20.0_wp]
         env%lai           = [3.0_wp, 1.5_wp]
         env%h_coeff_w     = [8.0_wp, 5.0_wp]
         env%abs_sw_wood   = [20.0_wp, 10.0_wp]
         env%abs_lw_wood   = [-5.0_wp, -3.0_wp]
         env%wai           = [0.6_wp, 0.3_wp]
         env%leaf_hcap_per_dt = [5.0_wp, 3.0_wp]
         env%wood_hcap_per_dt = [60.0_wp, 40.0_wp]
         env%t_leaf0 = [301.0_wp, 299.0_wp]
         env%t_wood0 = [299.5_wp, 300.5_wp]
         env%qwflux_wl  = [2.0_wp, 1.0_wp]
         env%q_wood_net = [-1.5_wp, 0.5_wp]
         env%transp_held = mode == 1_ik
         env%transp_net  = [3.0e-5_wp, 1.5e-5_wp]
         film%f_wet_c     = [0.3_wp, 0.6_wp]
         film%g_film_leaf = [0.02_wp, 0.01_wp]
         film%g_film_w    = [0.004_wp, 0.002_wp]
         film%rain_liquid_enthalpy = internal_energy_liquid(293.0_wp)
         film%transfer_lw    = [5.0e-6_wp, -4.0e-6_wp]               ! leaf -> wood in one, wood -> leaf in the other
         film%intercept_leaf = [2.0e-5_wp, 1.0e-5_wp] - film%transfer_lw
         film%intercept_wood = [1.0e-5_wp, 6.0e-6_wp] + film%transfer_lw
         call canopy_closure(env, film, resid, scale)
         worst = max(worst, abs(resid) / scale)
      end do
      call check_true('canopy_tissue_energy: energy is conserved across the canopy, its films included '//        &
                      '(transpiration held and free)', worst < 1.0e-12_wp, worst)
   end subroutine test_canopy_energy_closes

   !--- The canopy's energy residual: what entered less what left less what was stored.
   subroutine canopy_closure(env, film, resid, scale)
      type(plant_energy_env_t), intent(in)  :: env
      type(canopy_water_out_t), intent(in)  :: film
      real(wp),                 intent(out) :: resid, scale
      real(wp) :: tl(size(env%lai)), tw(size(env%lai)), tr(size(env%lai)), fl(size(env%lai)), fw(size(env%lai))
      real(wp) :: h, qw, transp, fe, rnet, u_l0(size(env%lai)), u_w0(size(env%lai)), gain, loss, stored
      call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      u_l0 = internal_energy_liquid(env%t_leaf0)
      u_w0 = internal_energy_liquid(env%t_wood0)
      gain   = rnet + sum(env%qwflux_wl + env%q_wood_net)                                                       &
               + sum(film%intercept_leaf + film%intercept_wood) * film%rain_liquid_enthalpy
      loss   = h + qw
      stored = sum(env%leaf_hcap_per_dt * (tl - env%t_leaf0) + env%wood_hcap_per_dt * (tw - env%t_wood0))       &
               + sum((film%intercept_leaf - fl) * u_l0 + (film%intercept_wood - fw) * u_w0)
      resid = gain - loss - stored
      scale = abs(rnet) + sum(abs(env%qwflux_wl) + abs(env%q_wood_net)) + abs(h) + abs(qw)                       &
              + sum(abs(film%intercept_leaf) + abs(film%intercept_wood)) * film%rain_liquid_enthalpy
   end subroutine canopy_closure

   !---~---
   !   Rain captured at a tissue's own temperature brings it no heat: capture on both tissues of a
   ! cohort at 300 K, the rain at 300 K, leaves both temperatures as they are with no capture,
   ! only the films growing; rain at 310 K warms both.
   !---~---
   subroutine test_film_capture()
      type(plant_energy_env_t) :: env
      type(canopy_water_out_t) :: film
      real(wp) :: tl0(1), tw0(1), tl(1), tw(1), tr(1), fl(1), fw(1), h, qw, transp, fe, rnet
      print '(a)', 'test_film_capture:'
      call one_wet_cohort(env, film, 300.0_wp, 300.0_wp)
      call run_canopy(env, film, tl0, tw0, tr, fl, fw, h, qw, transp, fe, rnet)
      film%intercept_leaf = 1.0e-5_wp
      film%intercept_wood = 1.0e-5_wp
      film%rain_liquid_enthalpy = internal_energy_liquid(300.0_wp)
      call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      call check('canopy_tissue_energy: rain captured at the tissues'' own temperature brings them no heat',       &
                 max(abs(tl(1) - tl0(1)), abs(tw(1) - tw0(1))), 0.0_wp, 1.0e-12_wp)
      film%rain_liquid_enthalpy = internal_energy_liquid(310.0_wp)
      call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      call check_true('canopy_tissue_energy: warmer rain warms the leaf and the wood that capture it',              &
                      tl(1) > tl0(1) .and. tw(1) > tw0(1), min(tl(1) - tl0(1), tw(1) - tw0(1)))
   end subroutine test_film_capture

   !---~---
   !   A film passed from leaf to wood (no net capture) leaves the leaf at the leaf's enthalpy and
   ! reaches the wood at the leaf's: with leaf and wood at one temperature it moves neither, and
   ! from a warmer leaf it leaves the leaf as it was and warms the wood.
   !---~---
   subroutine test_film_transfer()
      type(plant_energy_env_t) :: env
      type(canopy_water_out_t) :: film
      real(wp) :: tl0(1), tw0(1), tl(1), tw(1), tr(1), fl(1), fw(1), h, qw, transp, fe, rnet
      print '(a)', 'test_film_transfer:'
      call one_wet_cohort(env, film, 300.0_wp, 300.0_wp)
      call run_canopy(env, film, tl0, tw0, tr, fl, fw, h, qw, transp, fe, rnet)
      call pass_leaf_to_wood(film, 5.0e-6_wp)
      call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      call check('canopy_tissue_energy: a film passed between tissues at one temperature moves neither',           &
                 max(abs(tl(1) - tl0(1)), abs(tw(1) - tw0(1))), 0.0_wp, 1.0e-12_wp)
      call one_wet_cohort(env, film, 305.0_wp, 300.0_wp)
      call run_canopy(env, film, tl0, tw0, tr, fl, fw, h, qw, transp, fe, rnet)
      call pass_leaf_to_wood(film, 5.0e-6_wp)
      call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      call check('canopy_tissue_energy: the warmer leaf giving its film is charged nothing for it', tl(1), tl0(1),  &
                 1.0e-12_wp)
      call check_true('canopy_tissue_energy: the wood receiving it gains the leaf''s warmer water', tw(1) > tw0(1), &
                      tw(1) - tw0(1))
   end subroutine test_film_transfer

   subroutine pass_leaf_to_wood(film, rate)
      type(canopy_water_out_t), intent(inout) :: film
      real(wp),                 intent(in)    :: rate   ! [kg/m2/s]
      film%transfer_lw    = rate
      film%intercept_leaf = -rate
      film%intercept_wood = rate
   end subroutine pass_leaf_to_wood

   !---~---
   !   An evaporating film costs its tissue the vapour's enthalpy at the canopy air less the liquid's
   ! at the tissue's start-of-step temperature, the rest being the film's own: a wood with only its
   ! film's pathway (no radiation, no sensible heat), fully wet, loses heat cap*(T - T0) per
   ! kilogram evaporated of exactly h_vapour(T_cas) - u_liquid(T0), from 2 K above the canopy air
   ! and 3 K below it.
   !---~---
   subroutine test_film_evaporation()
      type(plant_energy_env_t) :: env
      type(canopy_water_out_t) :: film
      real(wp) :: tl(1), tw(1), tr(1), fl(1), fw(1), h, qw, transp, fe, rnet, t0, per_kg, worst
      integer(ik) :: k
      print '(a)', 'test_film_evaporation:'
      worst = 0.0_wp
      do k = 1_ik, 2_ik
         t0 = merge(T_CAS + 2.0_wp, T_CAS - 3.0_wp, k == 1_ik)
         call canopy_inputs(1_ik, env, film)
         env%t_leaf0 = t0
         env%t_wood0 = t0
         env%wood_hcap_per_dt = 50.0_wp
         film%f_wet_c = 1.0_wp
         film%g_film_w = 0.004_wp
         call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
         per_kg = -env%wood_hcap_per_dt(1) * (tw(1) - t0) / fw(1)
         worst = max(worst, abs(per_kg - (enthalpy_vapor(T_CAS) - internal_energy_liquid(t0)))                     &
                            / (enthalpy_vapor(T_CAS) - internal_energy_liquid(t0)))
      end do
      call check_true('canopy_tissue_energy: an evaporating film costs its tissue h_vapour(T_cas) - u_liquid(T0)',  &
                      worst < 1.0e-10_wp, worst)
   end subroutine test_film_evaporation

   !---~---
   !   5. THE TISSUE'S LONGWAVE. The two-stream's abs_lw is the tissue's net longwave with it
   ! emitting at its start-of-step temperature t0, abs_lw = LW_in - eps sigma A t0^4, and the
   ! balance carries the emission's change linearized about t0. One leaf (LAI 2, no latent,
   ! no heat store) in sun and a 298 K canopy air: when the tissue ends at t0, the
   ! linearization is exact and its temperature solves the Stefan-Boltzmann balance
   !   abs_sw + LW_in - eps sigma A T^4 - h (T - T_cas) = 0,
   ! found here by bisection, independently; t0 a distance d from there gives a temperature
   ! off by O(d^2), the linearization's error. A linearization about the canopy air's
   ! temperature instead is off at first order, and not exact anywhere but at T_cas.
   !---~---
   subroutine test_tissue_longwave()
      real(wp) :: t_star, t_fix, t, e(3), p1, p2, rnet, lo, hi, mid
      integer(ik) :: k
      print '(a)', 'test_tissue_longwave:'
      lo = 250.0_wp
      hi = 350.0_wp
      do k = 1_ik, 200_ik
         mid = 0.5_wp * (lo + hi)
         if (exact_balance(mid) > 0.0_wp) then
            lo = mid
         else
            hi = mid
         end if
      end do
      t_star = 0.5_wp * (lo + hi)
      !--- Start the emission at the canopy air and let it follow the tissue to where they meet.
      t_fix = T_CAS
      do k = 1_ik, 20_ik
         call leaf_in_sun(t_fix, t, rnet)
         t_fix = t
      end do
      call leaf_in_sun(t_fix, t, rnet)
      call check('canopy_tissue_energy: at its emission temperature the tissue solves the exact balance', t, t_star, &
                 1.0e-8_wp)
      call check('canopy_tissue_energy: there its net radiation is the exact emission''s', rnet,                    &
                 LEAF_SW + LEAF_LW_IN - LEAF_EMISS * stefan_boltzmann * LEAF_LAI * t ** 4, 1.0e-9_wp * LEAF_SW)
      do k = 1_ik, 3_ik
         call leaf_in_sun(t_star + 2.0_wp / 2.0_wp ** (k - 1_ik), t, rnet)
         e(k) = abs(t - t_star)
      end do
      p1 = log(e(1) / e(2)) / log(2.0_wp)
      p2 = log(e(2) / e(3)) / log(2.0_wp)
      print '(a,3es10.3,a,2f6.2)', '   temperature error at t0 - T* = 2, 1, 0.5 K:', e, ' ; order', p1, p2
      call check_true('canopy_tissue_energy: the linearized longwave is second order in t0 - T',                    &
                      abs(p1 - 2.0_wp) < 0.2_wp .and. abs(p2 - 2.0_wp) < 0.2_wp, p2)
   end subroutine test_tissue_longwave

   !---~---
   !   6. THE COMMIT (commit_tissue_energy). One cohort, one plant per m2, whose leaf and wood
   ! gain water over a 900 s step and relax from 300 K and 296 K towards 301.5 K and 296.8 K,
   ! each at its heat capacity over the step (tissue_hcap_per_dt: its dry matter and its
   ! start-of-step water).
   !---~---
   subroutine test_commit()
      real(wp), parameter :: DT = 900.0_wp
      real(wp), parameter :: DL = 2.0e3_wp, DW = 4.0e4_wp          ! [J/K] dry heat capacities
      real(wp), parameter :: WL0 = 0.8_wp, WW0 = 9.0_wp            ! [kg] start-of-step water
      real(wp), parameter :: DWL = 0.05_wp, DWW = 0.4_wp           ! [kg] the water each gains
      real(wp), parameter :: TL0 = 300.0_wp, TW0 = 296.0_wp, TL_AVG = 301.5_wp, TW_AVG = 296.8_wp
      type(plant_energy_env_t) :: env
      real(wp) :: e_leaf(1), e_wood(1), t_leaf(1), t_wood(1)
      real(wp) :: store0, store1, heat, enth, tl_expect, tw_expect
      print '(a)', 'test_commit:'
      call commit_env(env, DL, WL0, DW, WW0, TL0, TW0, DT)
      env%leaf_water_enth = DWL / DT * internal_energy_liquid(TL0)
      env%wood_water_enth = DWW / DT * internal_energy_liquid(TW0)
      e_leaf = tissue_temp_to_energy(DL, WL0, TL0)
      e_wood = tissue_temp_to_energy(DW, WW0, TW0)
      t_leaf = TL0
      t_wood = TW0
      call commit_tissue_energy(env, [1.0_wp], DT, [TL_AVG], [TW_AVG], 1_ik, [WL0 + DWL], [WW0 + DWW],          &
                                [0.0_wp], [0.0_wp], [0.0_wp], [0.0_wp], e_leaf, e_wood, t_leaf, t_wood,          &
                                store0, store1)
      heat = (DL + WL0 * cp_liq) * (TL_AVG - TL0) + (DW + WW0 * cp_liq) * (TW_AVG - TW0)
      enth = DWL * internal_energy_liquid(TL0) + DWW * internal_energy_liquid(TW0)
      tl_expect = ((DL + WL0 * cp_liq) * TL_AVG + DWL * cp_liq * TL0) / (DL + (WL0 + DWL) * cp_liq)
      tw_expect = ((DW + WW0 * cp_liq) * TW_AVG + DWW * cp_liq * TW0) / (DW + (WW0 + DWW) * cp_liq)
      call check_true('commit_tissue_energy: the tissues gain the heat they relaxed through and their '//      &
                      'water''s enthalpy', abs((store1 - store0) - (heat + enth)) <= 1.0e-12_wp * abs(heat + enth))
      call check_true('commit_tissue_energy: the gained water arrived at the tissues'' start-of-step '//       &
                      'temperatures', abs(t_leaf(1) - tl_expect) < 1.0e-9_wp .and. abs(t_wood(1) - tw_expect) < 1.0e-9_wp)
      call check_true('commit_tissue_energy: the committed temperatures give back the energies, so the next '// &
                      'step starts where this one ended',                                                      &
                      abs(tissue_temp_to_energy(DL, WL0 + DWL, t_leaf(1)) - e_leaf(1)) <= 1.0e-12_wp * abs(e_leaf(1)) &
                      .and. abs(tissue_temp_to_energy(DW, WW0 + DWW, t_wood(1)) - e_wood(1))                   &
                            <= 1.0e-12_wp * abs(e_wood(1)))
   end subroutine test_commit

   !---~---
   !   A nearly leafless cohort at Palo Verde's first rains (one plant per m2): its leaf holds a
   ! film of a thousand times its own heat capacity, which evaporates within the step while the
   ! leaf cools from 298 K to 297.5 K. The leaf's heat capacity over the step is its own, the
   ! film left out, and the commit reads the temperature back with the end-of-step film, so the
   ! leaf reads back at the step's temperature.
   !---~---
   subroutine test_commit_heavy_film()
      real(wp), parameter :: DT = 900.0_wp
      real(wp), parameter :: DL = 3.4e-8_wp, WL = 1.0e-11_wp, F0 = 8.4e-9_wp   ! [J/K], [kg], [kg]
      real(wp), parameter :: DW = 4.0e4_wp, WW = 9.0_wp                       ! [J/K], [kg]
      real(wp), parameter :: TL0 = 298.0_wp, TL_AVG = 297.5_wp, TW0 = 298.0_wp
      type(plant_energy_env_t) :: env
      real(wp) :: e_leaf(1), e_wood(1), t_leaf(1), t_wood(1), store0, store1
      print '(a)', 'test_commit_heavy_film:'
      call commit_env(env, DL, WL, DW, WW, TL0, TW0, DT)
      e_leaf = tissue_temp_to_energy(DL, WL + F0, TL0)
      e_wood = tissue_temp_to_energy(DW, WW, TW0)
      t_leaf = TL0
      t_wood = TW0
      call commit_tissue_energy(env, [1.0_wp], DT, [TL_AVG], [TW0], 1_ik, [WL], [WW], [F0], [0.0_wp], [0.0_wp],    &
                                [0.0_wp], e_leaf, e_wood, t_leaf, t_wood, store0, store1)
      call check('commit_tissue_energy: a film heavier than its tissue, evaporated within the step, leaves '//  &
                 'the tissue at the step''s temperature', t_leaf(1), TL_AVG, 1.0e-6_wp)
   end subroutine test_commit_heavy_film

   !---~---
   !   The commit's inputs for one cohort: each tissue's dry heat capacity, its heat capacity over
   ! the step (tissue_hcap_per_dt on its start-of-step water) and its start-of-step temperature,
   ! with no water arriving.
   !---~---
   subroutine commit_env(env, dry_leaf, w_leaf, dry_wood, w_wood, t_leaf0, t_wood0, dt)
      type(plant_energy_env_t), intent(out) :: env
      real(wp),                 intent(in)  :: dry_leaf, w_leaf, dry_wood, w_wood, t_leaf0, t_wood0, dt
      allocate(env%leaf_hcap_per_dt(1), env%wood_hcap_per_dt(1), env%t_leaf0(1), env%t_wood0(1),                 &
               env%leaf_dry_hcap(1), env%wood_dry_hcap(1), env%leaf_water_enth(1), env%wood_water_enth(1))
      env%leaf_dry_hcap    = dry_leaf
      env%wood_dry_hcap    = dry_wood
      env%leaf_hcap_per_dt = tissue_hcap_per_dt(dry_leaf, w_leaf, dt)
      env%wood_hcap_per_dt = tissue_hcap_per_dt(dry_wood, w_wood, dt)
      env%t_leaf0          = t_leaf0
      env%t_wood0          = t_wood0
      env%leaf_water_enth  = 0.0_wp
      env%wood_water_enth  = 0.0_wp
   end subroutine commit_env

   !--- The exact balance's residual at T.
   pure real(wp) function exact_balance(t) result(g)
      real(wp), intent(in) :: t
      g = LEAF_SW + LEAF_LW_IN - LEAF_EMISS * stefan_boltzmann * LEAF_LAI * t ** 4 - LEAF_H * (t - T_CAS)
   end function exact_balance

   !--- The leaf emitting at t0: its temperature and net radiation from canopy_tissue_energy.
   subroutine leaf_in_sun(t0, t_leaf, rnet)
      real(wp), intent(in)  :: t0
      real(wp), intent(out) :: t_leaf, rnet
      type(plant_energy_env_t) :: env
      type(canopy_water_out_t) :: film
      real(wp) :: tl(1), tw(1), tr(1), fl(1), fw(1), h, qw, transp, fe
      call canopy_inputs(1_ik, env, film)
      env%lai = LEAF_LAI
      env%leaf_emiss = LEAF_EMISS
      env%abs_sw = LEAF_SW
      env%h_coeff_leaf = LEAF_H
      env%t_leaf0 = t0
      env%abs_lw = LEAF_LW_IN - LEAF_EMISS * stefan_boltzmann * LEAF_LAI * t0 ** 4
      call run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      t_leaf = tl(1)
   end subroutine leaf_in_sun

   !---~---
   !   One cohort, partly wet, with leaf and wood exchanging sensible heat and longwave and holding
   ! heat, at the given start-of-step temperatures; no capture.
   !---~---
   subroutine one_wet_cohort(env, film, t_leaf0, t_wood0)
      type(plant_energy_env_t), intent(out) :: env
      type(canopy_water_out_t), intent(out) :: film
      real(wp),                 intent(in)  :: t_leaf0, t_wood0
      call canopy_inputs(1_ik, env, film)
      env%h_coeff_leaf = 30.0_wp
      env%g_transp_leaf = 0.005_wp
      env%abs_sw = 250.0_wp
      env%abs_lw = -30.0_wp
      env%lai = 2.5_wp
      env%h_coeff_w = 8.0_wp
      env%abs_sw_wood = 15.0_wp
      env%abs_lw_wood = -4.0_wp
      env%wai = 0.5_wp
      env%leaf_hcap_per_dt = 4.0_wp
      env%wood_hcap_per_dt = 50.0_wp
      env%t_leaf0 = t_leaf0
      env%t_wood0 = t_wood0
      film%f_wet_c = 0.5_wp
      film%g_film_leaf = 0.02_wp
      film%g_film_w = 0.004_wp
      film%rain_liquid_enthalpy = internal_energy_liquid(290.0_wp)
   end subroutine one_wet_cohort

   !---~---
   !   n cohorts with nothing on: every coefficient, flux, film and heat store zero, the tissues at
   ! the canopy air's temperature with the shipped emissivities, transpiration free.
   !---~---
   subroutine canopy_inputs(n, env, film)
      integer(ik),              intent(in)  :: n
      type(plant_energy_env_t), intent(out) :: env
      type(canopy_water_out_t), intent(out) :: film
      allocate(env%h_coeff_leaf(n), env%g_transp_leaf(n), env%abs_sw(n), env%abs_lw(n), env%lai(n),            &
               env%h_coeff_w(n), env%abs_sw_wood(n), env%abs_lw_wood(n), env%wai(n), env%leaf_hcap_per_dt(n),  &
               env%wood_hcap_per_dt(n), env%t_leaf0(n), env%t_wood0(n), env%qwflux_wl(n), env%q_wood_net(n),   &
               env%transp_net(n), env%leaf_emiss(n), env%wood_emiss(n), source=0.0_wp)
      env%t_leaf0 = T_CAS
      env%t_wood0 = T_CAS
      env%leaf_emiss = 0.97_wp
      env%wood_emiss = 0.90_wp
      env%transp_held = .false.
      allocate(film%f_wet_c(n), film%g_film_leaf(n), film%g_film_w(n), film%intercept_leaf(n),                   &
               film%intercept_wood(n), film%transfer_lw(n), source=0.0_wp)
      film%rain_liquid_enthalpy = internal_energy_liquid(T_CAS)
   end subroutine canopy_inputs

   !--- canopy_tissue_energy under the test canopy air.
   subroutine run_canopy(env, film, tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
      type(plant_energy_env_t), intent(in)  :: env
      type(canopy_water_out_t), intent(in)  :: film
      real(wp),                 intent(out) :: tl(:), tw(:), tr(:), fl(:), fw(:), h, qw, transp, fe, rnet
      tl = 0.0_wp
      tw = 0.0_wp
      tr = 0.0_wp
      fl = 0.0_wp
      fw = 0.0_wp
      call canopy_tissue_energy(env, film, T_CAS, Q_CAS, sat_specific_humidity(T_CAS, P_AIR),                    &
                                sat_specific_humidity_temp_deriv(T_CAS, P_AIR), RHO_AIR, size(env%lai, kind=ik), &
                                tl, tw, tr, fl, fw, h, qw, transp, fe, rnet)
   end subroutine run_canopy

end program test_plant_energy
