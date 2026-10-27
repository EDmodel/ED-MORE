! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_plant_hydraulics -- how water moves through a cohort's plants
! (ed_plant_hydraulics; docs/science/plant_hydraulics.md §2-§4,
! docs/science/soil_plant_water.md §1 and §4). The tissue curves and conductances it is
! built from are test_plant_hydr_lib's; the network that solves it is test_water_network's.
!   1. KIRCHHOFF (xylem_sapflow, kirchhoff_potential; eqs. (8)-(9)): the mean retained
!      fraction lies between the ends' and is the potential's difference quotient; equal
!      ends give the pointwise value; with gravity the law overstates the exact steady
!      flux by the margin §3 states (a = 2, closed form).
!   2. TABLE (kirchhoff_potential on build_kirchhoff_table): a tabled exponent against a
!      Simpson reference, and one table for any psi50.
!   3. FAR TAIL: far beyond the table the potential reaches its asymptote and a desiccated
!      leaf keeps a positive sapflow.
!   4. XYLEM (xylem_sapflow): the sapflow is the Kirchhoff mean conductance times the
!      drive, exponents 2 and 3, with derivatives that match finite differences.
!   5. ROOT UPTAKE (root_water_uptake, eq. (15)): derivatives against finite differences
!      on and off the ramp; linear in the drive past it; nothing flows backwards.
!   6. TRANSPIRATION (supply_limited_transpiration; soil_plant_water.md eq. (3)): zero at
!      an empty leaf, the demand above the closing water content, half of it at half, and
!      a derivative that matches finite differences.
!   7. COHORT ROOTS (cohort_root_profile, eqs. (12)-(13)): a 35 m tree roots to ED2's 5 m,
!      a sapling shallower than 2 m, and the depth stops at the soil column; the profile
!      sums to one within the rooting depth.
!   8. LEAFLESS (leaf_node_psi): a leafless cohort's leaf node sits at the wood less the
!      crown head; with leaves it reads the leaf curve.
!---~---
program test_plant_hydraulics
   use ed_test_assert,        only : check, check_true, test_report
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : grav_head, pi
   use ed_plant_hydr_lib,     only : conductance_retained_fraction, xylem_max_conductance, rhizosphere_conductance,  &
                                     tissue_psi_to_water, tissue_water_to_psi, water_curve_t
   use ed_pft_params,         only : hydraulics_pft_t, kirchhoff_table_t
   use ed_pft_derive,         only : build_kirchhoff_table
   use ed_plant_hydraulics,   only : leaf_node_psi, kirchhoff_potential, xylem_sapflow, cohort_root_profile,       &
                                     root_water_uptake, supply_limited_transpiration
   implicit none

   call test_xylem_kirchhoff_mean()
   call test_hydro_table()
   call test_kirchhoff_far_tail()
   call test_xylem()
   call test_root_uptake()
   call test_transpiration()
   call test_cohort_roots()
   call test_leafless_leaf_psi()

   call test_report('test_plant_hydraulics')

contains

   !--- A pioneer-like default trait set + cohort environment.
   subroutine defaults(p)
      type(hydraulics_pft_t), intent(out) :: p
      p%leaf_curve = water_curve_t(pi0 = -1.5_wp, elastic_mod = 12.0_wp, apoplast_frac = 0.30_wp, water_sat = 2.0_wp)
      p%wood_curve = water_curve_t(pi0 = -1.0_wp, elastic_mod =  8.0_wp, apoplast_frac = 0.20_wp, water_sat = 1.0_wp)
      p%wood_psi50 = -2.0_wp
      p%wood_kexp  = 2.0_wp
      p%wood_kmax  = 9.0_wp
      p%vessel_curl = 1.5_wp      ! 3.0e-3 kg/s/MPa over the 20 m path
      p%root_beta = 0.1_wp
      p%root_depth_b1 = 1.1140580_wp
      p%root_depth_b2 = 0.4223014_wp
      p%specific_root_length = 2.0e4_wp
      p%fine_root_radius = 2.5e-4_wp
   end subroutine defaults

   !--- Relative agreement of an analytic derivative with a central difference.
   logical function fd_agrees(analytic, numeric, rtol, atol)
      real(wp), intent(in) :: analytic, numeric, rtol, atol
      fd_agrees = abs(analytic - numeric) <= rtol * max(abs(analytic), abs(numeric)) + atol
   end function fd_agrees

   subroutine test_xylem_kirchhoff_mean()
      real(wp), parameter :: PSI50 = -2.0_wp, KEXP = 2.0_wp, KCOND = 3.0e-3_wp
      type(kirchhoff_table_t) :: tab
      real(wp) :: psi_l, psi_w, avg_plc, s, s0, dum1, dum2, pmid
      print '(a)', '-- Kirchhoff mean conductance of the xylem --'
      call build_kirchhoff_table(tab, KEXP)
      psi_l = -1.2_wp
      psi_w = -0.4_wp
      !--- With no gravity head the sapflow over the drop is the Kirchhoff mean conductance.
      call xylem_sapflow(KCOND, psi_w, psi_l, 0.0_wp, PSI50, KEXP, tab, s, dum1, dum2)
      avg_plc = s / ((psi_w - psi_l) * KCOND)
      !--- <plc> lies between the two endpoint retained fractions.
      call check_true('<plc> in (plc(leaf), plc(wood))',                                       &
           avg_plc > conductance_retained_fraction(psi_l,PSI50,KEXP)                                &
           .and. avg_plc < conductance_retained_fraction(psi_w,PSI50,KEXP))
      !--- <plc> = (Phi(up)-Phi(down))/(up-down).
      call check('<plc> = DeltaPhi/Deltapsi', avg_plc,                                          &
           (kirchhoff_potential(psi_w,PSI50,KEXP,tab)-kirchhoff_potential(psi_l,PSI50,KEXP,tab))/(psi_w-psi_l),  &
           1.0e-11_wp)
      !--- |dpsi|->0 limit is the pointwise value at the mean.
      pmid  = -0.9_wp
      call xylem_sapflow(KCOND, pmid, pmid - 1.0e-9_wp, 0.0_wp, PSI50, KEXP, tab, s0, dum1, dum2)
      call check('Kirchhoff dpsi->0 limit', s0 / (1.0e-9_wp * KCOND), conductance_retained_fraction(pmid,PSI50,KEXP), &
                 1.0e-6_wp)
      !---~---
      !   With gravity the law is not the steady flux. For a = 2 the exact one solves
      ! psi50 / sqrt(c Q) [atan(r_w lam) - atan(r_l lam)] = 1, c = Q + G, lam = sqrt(Q/c),
      ! r = psi/psi50 (Q = flux/k_cond). At psi_w -1, psi_l -3, G 0.5 and psi50 -1.5 MPa the law
      ! overstates it by 3.4 % (docs/science/plant_hydraulics.md §3).
      !---~---
      block
         real(wp), parameter :: P50G = -1.5_wp, PW = -1.0_wp, PL = -3.0_wp, G = 0.5_wp
         type(kirchhoff_table_t) :: tab_g
         real(wp) :: q_law, q_lo, q_hi, q_mid
         integer(ik) :: it
         call build_kirchhoff_table(tab_g, 2.0_wp)
         call xylem_sapflow(1.0_wp, PW, PL, G, P50G, 2.0_wp, tab_g, q_law, dum1, dum2)
         q_lo = 1.0e-9_wp
         q_hi = 2.0_wp * q_law
         do it = 1_ik, 200_ik
            q_mid = 0.5_wp * (q_lo + q_hi)
            if (steady_residual(q_mid) > 0.0_wp) then
               q_lo = q_mid
            else
               q_hi = q_mid
            end if
         end do
         call check('with gravity the law overstates the steady flux by 3.4 %', q_law / q_lo - 1.0_wp,    &
                    0.0336_wp, 5.0e-4_wp)
      end block
   end subroutine test_xylem_kirchhoff_mean

   !---~---
   !   The exact a = 2 steady-flux condition with gravity at psi_w -1, psi_l -3, G 0.5 and psi50 -1.5
   ! MPa: positive below the root Q, negative above it.
   !---~---
   pure real(wp) function steady_residual(q) result(f)
      real(wp), intent(in) :: q   ! [MPa] the flux per k_cond
      real(wp), parameter :: P50 = -1.5_wp, PW = -1.0_wp, PL = -3.0_wp, G = 0.5_wp
      real(wp) :: c, lam
      c   = q + G
      lam = sqrt(q / c)
      f   = P50 / sqrt(c * q) * (atan(lam * PW / P50) - atan(lam * PL / P50)) - 1.0_wp
   end function steady_residual

   !--- integral_0^r du/(1+u^kexp) by composite Simpson on 2000 intervals: the reference.
   pure real(wp) function simpson_reference(r, kexp) result(g)
      real(wp), intent(in) :: r, kexp
      integer(ik), parameter :: NINT = 2000_ik
      real(wp)    :: h, u
      integer(ik) :: k
      h = r / real(NINT, wp)
      g = 1.0_wp + 1.0_wp / (1.0_wp + r**kexp)
      do k = 1_ik, NINT - 1_ik
         u = real(k, wp) * h
         g = g + merge(4.0_wp, 2.0_wp, mod(k, 2_ik) == 1_ik) / (1.0_wp + u**kexp)
      end do
      g = g * h / 3.0_wp
   end function simpson_reference

   !--- The Kirchhoff potential [MPa] from the Simpson reference, for r = psi/psi50 <= 8.
   pure real(wp) function kirchhoff_reference(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      phi = psi50 * simpson_reference(psi / psi50, kexp)
   end function kirchhoff_reference

   !---~---
   !   The precomputed lookup table (build_kirchhoff_table, read through kirchhoff_potential)
   ! reproduces a Simpson reference of the Kirchhoff potential for a general (non-{1,2})
   ! exponent, and the table-backed xylem conductance matches the reference's. kexp is stamped for
   ! the staleness guard, and psi50 is a runtime scale (one table serves any psi50).
   !---~---
   subroutine test_hydro_table()
      real(wp), parameter :: KEXP = 3.0_wp, PSI50 = -2.0_wp, KCOND = 3.0e-3_wp
      type(kirchhoff_table_t) :: tab
      real(wp) :: psi, err, emax, psi_l, psi_w, s, dum1, dum2
      integer(ik) :: i
      print '(a)', '-- Kirchhoff lookup table --'
      call build_kirchhoff_table(tab, KEXP)
      call check('table stamps kexp (staleness guard)', tab%kexp, KEXP, 0.0_wp)
      !--- The table matches the reference across the operating range.
      emax = 0.0_wp
      do i = 0_ik, 80_ik
         psi = PSI50 * (0.05_wp + real(i,wp)*0.05_wp)          ! r in [0.05, 4.05]
         err = abs(kirchhoff_potential(psi, PSI50, KEXP, tab) - kirchhoff_reference(psi, PSI50, KEXP))
         emax = max(emax, err)
      end do
      call check_true('linear table within 1e-4 MPa of the reference (kexp=3)', emax < 1.0e-4_wp)
      !--- The table-backed xylem conductance matches the reference's Kirchhoff mean.
      psi_l = -1.7_wp
      psi_w = -0.6_wp
      call xylem_sapflow(KCOND, psi_w, psi_l, 0.0_wp, PSI50, KEXP, tab, s, dum1, dum2)
      call check('table-backed xylem conductance ~ reference (kexp=3)', s / (psi_w - psi_l),      &
           KCOND * (kirchhoff_reference(psi_w, PSI50, KEXP) - kirchhoff_reference(psi_l, PSI50, KEXP)) &
           / (psi_w - psi_l), 1.0e-6_wp)
      !--- psi50 is a runtime scale: the SAME table serves a different psi50.
      call check('table serves a different psi50 (r-normalized)',                                &
           kirchhoff_potential(-3.0_wp, -3.0_wp, KEXP, tab), kirchhoff_reference(-3.0_wp, -3.0_wp, KEXP), &
           1.0e-4_wp)
   end subroutine test_hydro_table

   !---~---
   !   Far beyond the table (a leaf dried to the apoplastic floor reads about pi0 * 1e4 MPa), the
   ! Kirchhoff potential approaches its asymptote psi50 * (pi/a)/sin(pi/a), and the wood-to-leaf
   ! edge keeps a positive conductance, so the leaf can refill. At its edge the table equals the
   ! Simpson reference.
   !---~---
   subroutine test_kirchhoff_far_tail()
      real(wp), parameter :: KEXP = 4.0_wp, PSI50 = -3.209_wp, KCOND = 5.0e-3_wp
      type(kirchhoff_table_t) :: tab
      real(wp) :: g_inf, psi_far, s, dum1, dum2
      print '(a)', '-- Kirchhoff potential far beyond the table --'
      call build_kirchhoff_table(tab, KEXP)
      g_inf   = (pi/KEXP) / sin(pi/KEXP)
      psi_far = -25040.0_wp
      call check('table reaches the asymptote', kirchhoff_potential(psi_far, PSI50, KEXP, tab), PSI50*g_inf,   &
                 1.0e-4_wp)
      call check('table at its edge = reference', kirchhoff_potential(8.0_wp*PSI50, PSI50, KEXP, tab),      &
                 kirchhoff_reference(8.0_wp*PSI50, PSI50, KEXP), 1.0e-9_wp)
      call xylem_sapflow(KCOND, -1.33_wp, psi_far, 0.0_wp, PSI50, KEXP, tab, s, dum1, dum2)
      call check_true('desiccated-leaf sapflow is positive', s > 0.0_wp)
      call check('desiccated-leaf sapflow = k (Phi_wood - Phi_inf)', s,                                   &
                 KCOND * (kirchhoff_reference(-1.33_wp, PSI50, KEXP) - PSI50*g_inf), 1.0e-6_wp)
   end subroutine test_kirchhoff_far_tail

   !---~---
   !   The sapflow of an 18 m cohort over drops from -0.6 to +0.6 MPa (through equal ends and the
   ! midpoint expansion), against the Kirchhoff mean conductance times the drive, and its two
   ! derivatives against central differences: exactly for exponent 2, and to the table's
   ! resolution for a tabled 3, whose slope is piecewise.
   !---~---
   subroutine test_xylem()
      type(kirchhoff_table_t) :: tab3, tab2
      real(wp) :: kc, s, dsw, dsl, sp, sm, dum, h, pw, pl, ref, grav
      integer  :: i
      logical  :: ok2, ok3, okv
      print '(a)', '-- Xylem flow --'
      call build_kirchhoff_table(tab3, 3.0_wp)
      call build_kirchhoff_table(tab2, 2.0_wp)
      kc = xylem_max_conductance(3.3_wp, 2.0e-3_wp, 18.0_wp, 1.0_wp) * 0.05_wp
      grav = grav_head * 18.0_wp
      ok2 = .true.
      ok3 = .true.
      okv = .true.
      do i = 0, 40
         pw = -0.1_wp - 0.15_wp * real(i, wp)
         pl = pw - 0.6_wp + 0.03_wp * real(i, wp)          ! crosses pw - pl = 0 at i = 20
         if (i == 20) pl = pw - 3.0e-7_wp                  ! inside the midpoint expansion
         h = 1.0e-7_wp
         call xylem_sapflow(kc, pw, pl, grav, -2.5_wp, 2.0_wp, tab2, s, dsw, dsl)
         ref = kirchhoff_mean(pw, pl, kc, -2.5_wp, 2.0_wp, tab2) * (pw - pl - grav)
         if (abs(s - ref) > 1.0e-9_wp * abs(ref) + 1.0e-11_wp * kc) okv = .false.
         call xylem_sapflow(kc, pw + h, pl, grav, -2.5_wp, 2.0_wp, tab2, sp, dum, dum)
         call xylem_sapflow(kc, pw - h, pl, grav, -2.5_wp, 2.0_wp, tab2, sm, dum, dum)
         if (.not. fd_agrees(dsw, (sp - sm) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok2 = .false.
         call xylem_sapflow(kc, pw, pl + h, grav, -2.5_wp, 2.0_wp, tab2, sp, dum, dum)
         call xylem_sapflow(kc, pw, pl - h, grav, -2.5_wp, 2.0_wp, tab2, sm, dum, dum)
         if (.not. fd_agrees(dsl, (sp - sm) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok2 = .false.
         !---~---
         !   A tabled vulnerability exponent: the table's slope is piecewise, so the analytic
         ! slope (the exact retained fraction) agrees only to the table's resolution.
         !---~---
         call xylem_sapflow(kc, pw, pl, grav, -2.5_wp, 3.0_wp, tab3, s, dsw, dsl)
         !---~---
         !   At near-equal ends the reference leaves its table for the exact retained fraction;
         ! the solver stays on the table, so the two differ there by its resolution.
         !---~---
         ref = kirchhoff_mean(pw, pl, kc, -2.5_wp, 3.0_wp, tab3) * (pw - pl - grav)
         if (i /= 20 .and. abs(s - ref) > 1.0e-9_wp * abs(ref) + 1.0e-11_wp * kc) okv = .false.
         h = 1.0e-4_wp
         call xylem_sapflow(kc, pw + h, pl, grav, -2.5_wp, 3.0_wp, tab3, sp, dum, dum)
         call xylem_sapflow(kc, pw - h, pl, grav, -2.5_wp, 3.0_wp, tab3, sm, dum, dum)
         if (.not. fd_agrees(dsw, (sp - sm) / (2*h), 2.0e-2_wp, 1.0e-8_wp * kc)) then
            ok3 = .false.
            print '(a,i3,4es14.6)', '    table case i, pw, pl, analytic, fd: ', i, pw, pl, dsw, (sp - sm) / (2*h)
         end if
      end do
      call check_true('xylem: sapflow = Kirchhoff mean conductance x driving force (exponents 2 and 3)', okv)
      call check_true('xylem: derivatives match finite differences (exponent 2, incl. equal ends)', ok2)
      call check_true('xylem: derivatives match the table to its resolution (exponent 3)', ok3)
   end subroutine test_xylem

   !---~---
   !   The reference Kirchhoff mean conductance, written independently of xylem_sapflow:
   ! k (Phi(up) - Phi(down)) / (up - down) from kirchhoff_potential (the closed form for
   ! exponent 2, the table for 3), and the retained fraction at the midpoint for near-equal
   ! ends.
   !---~---
   pure real(wp) function kirchhoff_mean(psi_up, psi_down, k, psi50, kexp, tab) result(keff)
      real(wp),            intent(in) :: psi_up, psi_down, k, psi50, kexp
      type(kirchhoff_table_t), intent(in) :: tab
      if (abs(psi_up - psi_down) <= 1.0e-6_wp) then
         keff = k * conductance_retained_fraction(0.5_wp * (psi_up + psi_down), psi50, kexp)
      else
         keff = k * (kirchhoff_potential(psi_up, psi50, kexp, tab) - kirchhoff_potential(psi_down, psi50, kexp, tab)) &
                / (psi_up - psi_down)
      end if
   end function kirchhoff_mean

   !---~---
   !   Root uptake from a layer 0.4 m down, at Gardner's conductance for 0.02 plants/m2 with 35 m of
   ! root each: its three derivatives against central differences, off the ramp and half a width
   ! onto it; past the ramp, twice the drive takes up twice the water; and a soil drier than the
   ! wood gives nothing.
   !---~---
   subroutine test_root_uptake()
      real(wp) :: g, u, dus, duw, duk, up, um, dum, h, k, ps, pw, z, u2
      logical  :: ok
      print '(a)', '-- Root uptake --'
      g = rhizosphere_conductance(0.02_wp, 35.0_wp, 0.01_wp, 2.5e-4_wp)
      k  = 1.0e-8_wp
      ps = -0.3_wp
      pw = -0.9_wp
      z  = -0.4_wp
      call root_water_uptake(g, k, ps, z, pw, u, dus, duw, duk)
      ok = .true.
      h = 1.0e-7_wp
      call root_water_uptake(g, k, ps + h, z, pw, up, dum, dum, dum)
      call root_water_uptake(g, k, ps - h, z, pw, um, dum, dum, dum)
      if (.not. fd_agrees(dus, (up - um) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call root_water_uptake(g, k, ps, z, pw + h, up, dum, dum, dum)
      call root_water_uptake(g, k, ps, z, pw - h, um, dum, dum, dum)
      if (.not. fd_agrees(duw, (up - um) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call root_water_uptake(g, k * (1.0_wp + 1.0e-6_wp), ps, z, pw, up, dum, dum, dum)
      call root_water_uptake(g, k * (1.0_wp - 1.0e-6_wp), ps, z, pw, um, dum, dum, dum)
      if (.not. fd_agrees(duk, (up - um) / (2.0e-6_wp * k), 1.0e-6_wp, 0.0_wp)) ok = .false.
      !--- On the ramp: half a width of driving force.
      call root_water_uptake(g, k, pw - grav_head * z + 5.0e-3_wp, z, pw, u2, dus, duw, duk)
      call root_water_uptake(g, k, pw - grav_head * z + 5.0e-3_wp + h, z, pw, up, dum, dum, dum)
      call root_water_uptake(g, k, pw - grav_head * z + 5.0e-3_wp - h, z, pw, um, dum, dum, dum)
      if (.not. fd_agrees(dus, (up - um) / (2*h), 1.0e-5_wp, 0.0_wp)) ok = .false.
      call check_true('root: derivatives match finite differences, on and off the ramp', ok)
      !---~---
      !   Past the ramp the uptake is linear in the drive: the soil 0.6 MPa wetter than the wood
      ! less its depth's head, and then 1.2 MPa.
      !---~---
      call root_water_uptake(g, k, pw - grav_head * z + 0.6_wp, z, pw, u, dus, duw, duk)
      call root_water_uptake(g, k, pw - grav_head * z + 1.2_wp, z, pw, u2, dus, duw, duk)
      call check('root: past the ramp, twice the drive takes up twice the water', u2, 2.0_wp * u, 1.0e-12_wp * u2)
      call root_water_uptake(g, k, -2.0_wp, z, pw, u, dus, duw, duk)
      call check('root: no flow from soil drier than the wood', u, 0.0_wp, 0.0_wp)
   end subroutine test_root_uptake

   !---~---
   !   A leaf transpiring 2e-5 kg/m2/s on demand, whose stomata close at a relative water content of
   ! 0.42.
   !---~---
   subroutine test_transpiration()
      real(wp) :: e, de, ep, em, dum, r, h
      integer  :: i
      logical  :: ok
      print '(a)', '-- Transpiration --'
      ok = .true.
      do i = 1, 99
         r = 0.01_wp * real(i, wp) + 1.0e-6_wp
         h = 1.0e-7_wp
         call supply_limited_transpiration(2.0e-5_wp, r, 0.42_wp, e, de)
         call supply_limited_transpiration(2.0e-5_wp, r + h, 0.42_wp, ep, dum)
         call supply_limited_transpiration(2.0e-5_wp, r - h, 0.42_wp, em, dum)
         if (.not. fd_agrees(de, (ep - em) / (2*h), 1.0e-5_wp, 1.0e-15_wp)) ok = .false.
      end do
      call check_true('transpiration: derivative matches finite differences', ok)
      call supply_limited_transpiration(2.0e-5_wp, 0.0_wp, 0.42_wp, e, de)
      call check('transpiration: an empty leaf transpires nothing', e, 0.0_wp, 0.0_wp)
      call supply_limited_transpiration(2.0e-5_wp, 0.5_wp, 0.42_wp, e, de)
      call check('transpiration: the demand stands above the closing water content', e, 2.0e-5_wp, 0.0_wp)
      call supply_limited_transpiration(2.0e-5_wp, 0.21_wp, 0.42_wp, e, de)
      call check('transpiration: half the closing water content, half the demand', e, 1.0e-5_wp, 1.0e-15_wp)
   end subroutine test_transpiration

   !---~---
   !   Cohort roots through the profile the network reads: the rooting depth follows ED2's height
   ! allometry (IALLOM 1, 5 m for a 35 m tree), so the share above 1 m is ED2's
   ! (1 - beta^(1/5)) / (1 - beta) and nothing lies below 5 m; a 3 m sapling roots above 2 m; and a
   ! column shallower than the tree stops the depth at its bottom, its share above 1 m
   ! (1 - beta^(1/2)) / (1 - beta). The profile sums to one.
   !---~---
   subroutine test_cohort_roots()
      type(hydraulics_pft_t) :: hp
      real(wp), parameter :: ED2_DEPTH_35M = 5.0_wp   ! [m] ED2's rooting depth of a 35 m tree
      real(wp)    :: edges(13), frac(12), frac2(4)
      integer(ik) :: k
      print '(a)', '-- Cohort roots --'
      call defaults(hp)                             ! ED2 IALLOM 1, beta 0.1
      edges = [(0.5_wp * k, k = 0, 12)]                ! 0 to 6 m in 0.5 m layers
      call cohort_root_profile(35.0_wp, hp, edges(1:12), edges(2:13), frac)
      call check('the profile sums to one', sum(frac), 1.0_wp, 1.0e-12_wp)
      call check('a 35 m tree roots to 5 m: ED2''s share above 1 m', sum(frac(1:2)),                       &
                 (1.0_wp - hp%root_beta**(1.0_wp / ED2_DEPTH_35M)) / (1.0_wp - hp%root_beta), 1.0e-4_wp)
      call check_true('a 35 m tree has no roots below 5 m', all(frac(11:12) == 0.0_wp) .and. frac(10) > 0.0_wp,  &
                      frac(10))
      call cohort_root_profile(3.0_wp, hp, edges(1:12), edges(2:13), frac)
      call check_true('a 3 m sapling roots shallower than 2 m', all(frac(5:12) == 0.0_wp) .and. frac(4) > 0.0_wp, &
                      frac(4))
      call cohort_root_profile(35.0_wp, hp, edges(1:4), edges(2:5), frac2)
      call check('a 2 m column stops the rooting depth at its bottom', sum(frac2(1:2)),                     &
                 (1.0_wp - hp%root_beta**(1.0_wp / 2.0_wp)) / (1.0_wp - hp%root_beta), 1.0e-12_wp)
   end subroutine test_cohort_roots

   !---~---
   !   A leafless cohort (leaf carbon 0, an empty leaf store) has no leaf curve to read: its leaf
   ! node takes the wood potential less the gravity head to the crown, not the floored reading
   ! (about pi0 * 1e4 MPa). A cohort with leaves still reads its own leaf water.
   !---~---
   subroutine test_leafless_leaf_psi()
      type(hydraulics_pft_t) :: p
      real(wp) :: bwood, w_wood, w_leaf, psi_w, psi_l
      real(wp), parameter :: BLEAF = 0.5_wp, HEIGHT = 20.0_wp
      print '(a)', '-- Leafless cohort: leaf node at the wood potential less gravity --'
      call defaults(p)
      bwood  = 7.0_wp
      w_wood = tissue_psi_to_water(-1.2_wp, p%wood_curve, bwood)
      psi_w  = tissue_water_to_psi(w_wood, p%wood_curve, bwood)
      psi_l  = leaf_node_psi(0.0_wp, 0.0_wp, w_wood, bwood, HEIGHT, p)
      call check('leafless: wood psi less gravity', psi_l, psi_w - grav_head*HEIGHT, 1.0e-12_wp)
      w_leaf = tissue_psi_to_water(-0.8_wp, p%leaf_curve, BLEAF)
      psi_l  = leaf_node_psi(w_leaf, BLEAF, w_wood, bwood, HEIGHT, p)
      call check('with leaves: the leaf curve', psi_l, tissue_water_to_psi(w_leaf, p%leaf_curve, BLEAF), 0.0_wp)
      call check('with leaves: -0.8 MPa recovered', psi_l, -0.8_wp, 1.0e-9_wp)
   end subroutine test_leafless_leaf_psi

end program test_plant_hydraulics
