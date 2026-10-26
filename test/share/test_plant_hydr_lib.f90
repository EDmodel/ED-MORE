! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_hydr_lib -- the plant water laws (ed_plant_hydr_lib;                          !
! docs/science/plant_hydraulics.md §2-§4).                                                 !
!   1. PRESSURE-VOLUME (tissue_rwc_to_psi, tissue_psi_to_rwc, tissue_psi_to_water,         !
!      tissue_capacitance; eqs. (3)-(6)): psi and RWC round-trip in the turgid and flaccid !
!      regions; a tissue at full turgor holds its saturated water; the capacitance is      !
!      dW/dpsi.                                                                            !
!   2. TISSUE RELATIONS (tissue_water_relations, tissue_water_to_rwc;                      !
!      soil_plant_water.md §2): the water network's tissue unknown gives the               !
!      pressure-volume curve's potential and water, with derivatives that match finite     !
!      differences.                                                                        !
!   3. SAPWOOD (xylem_max_conductance, eq. (10)): the intact conductance scales with the   !
!      conductivity and the sapwood area and inversely with the path, and stays positive   !
!      with no sapwood.                                                                    !
!   4. ROOT PROFILE (root_fraction_profile, eq. (13)): shallower layers hold more roots,   !
!      the column to the rooting depth holds 1 - beta of them, and none lie below it.      !
!   5. GARDNER (rhizosphere_conductance, eq. (14)): a layer split in two keeps its         !
!      conductance, and twice the plants at one spacing double it.                         !
!   6. GROWTH SEAM (tissue_water_to_psi): growth at conserved water lowers the potential.  !
!   7. CAPACITY CLAMP (tissue_capped_water): a no-op on growth, the new ceiling on a       !
!      shrink, zero for a bare tissue.                                                     !
!==========================================================================================!
program test_plant_hydr_lib
   use ed_test_assert,      only : check, check_true, test_report
   use ed_kinds,            only : wp, ik
   use ed_constants,        only : pi
   use ed_plant_hydr_lib,   only : water_curve_t, tissue_psi_to_rwc, tissue_rwc_to_psi, tissue_psi_to_water,      &
                                   tissue_capacitance, tissue_water_to_psi, tissue_capped_water,                  &
                                   tissue_water_relations, tissue_water_to_rwc, xylem_max_conductance,            &
                                   root_fraction_profile, rhizosphere_conductance
   implicit none

   call test_pv_curve()
   call test_tissue_water_relations()
   call test_sapwood_conductance()
   call test_root_profile()
   call test_gardner()
   call test_growth_seam()
   call test_seam_capacity_clamp()

   call test_report('test_plant_hydr_lib')

contains

   !----- Relative agreement of an analytic derivative with a central difference. -----------!
   logical function fd_agrees(analytic, numeric, rtol, atol)
      real(wp), intent(in) :: analytic, numeric, rtol, atol
      fd_agrees = abs(analytic - numeric) <= rtol * max(abs(analytic), abs(numeric)) + atol
   end function fd_agrees

   !----- 1. A leaf of osmotic potential -1.5 MPa, elastic modulus 12 MPa, 30 % apoplast and 2 kg   !
   !      of water per kgC at saturation, 0.5 kgC of it. -----------------------------------------!
   subroutine test_pv_curve()
      real(wp), parameter :: PI0 = -1.5_wp, EPS = 12.0_wp, AF = 0.3_wp, WS = 2.0_wp, BM = 0.5_wp
      type(water_curve_t), parameter :: CURVE = water_curve_t(PI0, EPS, AF, WS)
      real(wp) :: r, psi, cnum, cana, dp, w1, w2
      print '(a)', '-- PV curve --'
      !----- psi<->RWC round trip in the turgid region. -------------------------------!
      r   = 0.9_wp
      psi = tissue_rwc_to_psi(r, PI0, EPS)
      call check('rwc round-trip (turgid)', tissue_psi_to_rwc(psi, PI0, EPS), r, 1.0e-9_wp)
      !----- ... and in the flaccid region. -------------------------------------------!
      r   = 0.5_wp*((PI0+EPS)/EPS)
      psi = tissue_rwc_to_psi(r, PI0, EPS)
      call check('rwc round-trip (flaccid)', tissue_psi_to_rwc(psi, PI0, EPS), r, 1.0e-9_wp)
      !----- Full-turgor water = saturated water. -------------------------------------!
      call check('W(psi=0) = water_sat*biomass', tissue_psi_to_water(0.0_wp, CURVE, BM),  &
                 WS*BM, 1.0e-9_wp)
      !----- Capacitance = finite-difference dW/dpsi. ---------------------------------!
      psi  = -0.8_wp
      dp   = 1.0e-5_wp
      w1   = tissue_psi_to_water(psi - dp, CURVE, BM)
      w2   = tissue_psi_to_water(psi + dp, CURVE, BM)
      cnum = (w2 - w1)/(2.0_wp*dp)
      cana = tissue_capacitance(psi, CURVE, BM)
      call check('C(psi) ~ dW/dpsi (turgid)', cana, cnum, 1.0e-6_wp)
   end subroutine test_pv_curve

   !----- 2. The network's tissue relations over relative water contents 0.01 to 1.2, on a leaf of  !
   !      0.3 kgC per plant at 0.05 plants/m2: derivatives against central differences (off the     !
   !      turgor-loss point, where the slope of psi has a corner), and the potential and water      !
   !      against the pressure-volume read-offs. -----------------------------------------------!
   subroutine test_tissue_water_relations()
      type(water_curve_t) :: leaf
      real(wp) :: r, h, psi, dpsi, w, dw, psip, psim, wpl, wm, dum
      integer  :: i
      logical  :: ok_d, ok_v
      print '(a)', '-- tissue relations --'
      leaf = water_curve_t(pi0=-1.8_wp, elastic_mod=12.0_wp, apoplast_frac=0.2_wp, water_sat=2.5_wp)
      ok_d = .true.
      ok_v = .true.
      do i = 1, 120
         r = 0.01_wp * real(i, wp)                                  ! 0.01 .. 1.2
         call tissue_water_relations(leaf, 0.3_wp, 0.05_wp, r, psi, dpsi, w, dw)
         h = 1.0e-7_wp
         call tissue_water_relations(leaf, 0.3_wp, 0.05_wp, r + h, psip, dum, wpl, dum)
         call tissue_water_relations(leaf, 0.3_wp, 0.05_wp, r - h, psim, dum, wm, dum)
         if (abs(r - (leaf%pi0 + leaf%elastic_mod) / leaf%elastic_mod) > 1.0e-6_wp) then
            if (.not. fd_agrees(dpsi, (psip - psim) / (2*h), 1.0e-5_wp, 1.0e-10_wp)) ok_d = .false.
         end if
         if (.not. fd_agrees(dw, (wpl - wm) / (2*h), 1.0e-6_wp, 1.0e-12_wp)) ok_d = .false.
         if (r >= 1.0e-4_wp .and. abs(psi - tissue_rwc_to_psi(r, leaf%pi0, leaf%elastic_mod)) > 1.0e-10_wp)    &
            ok_v = .false.
         if (abs(w - 0.05_wp * tissue_psi_to_water(psi, leaf, 0.3_wp)) > 1.0e-10_wp) ok_v = .false.
         if (abs(tissue_water_to_rwc(leaf, 0.3_wp, w / 0.05_wp) - r) > 1.0e-10_wp) ok_v = .false.
      end do
      call check_true('tissue: derivatives match finite differences', ok_d)
      call check_true('tissue: potential and water match the pressure-volume curve', ok_v)
   end subroutine test_tissue_water_relations

   !----- 3. The maximum conductance scales with wood_kmax and the sapwood area and inversely with   !
   !      the path factor: doubling the conductivity, doubling the sapwood area and halving the path  !
   !      factor give the same step, and a faster one than the base. -----------------------------!
   subroutine test_sapwood_conductance()
      real(wp) :: k0
      print '(a)', '-- Sapwood conductance --'
      k0 = xylem_max_conductance(9.0_wp, 0.01_wp, 20.0_wp, 1.5_wp)
      call check('2 x sapwood area = 2 x wood_kmax', xylem_max_conductance(9.0_wp, 0.02_wp, 20.0_wp, 1.5_wp), &
                 xylem_max_conductance(18.0_wp, 0.01_wp, 20.0_wp, 1.5_wp), 1.0e-15_wp)
      call check('vessel_curl / 2 = 2 x the conductance', xylem_max_conductance(9.0_wp, 0.01_wp, 20.0_wp, 0.75_wp), &
                 2.0_wp * k0, 1.0e-15_wp)
      call check_true('a plant with no sapwood keeps a positive conductance',                         &
                      xylem_max_conductance(9.0_wp, 0.0_wp, 20.0_wp, 1.5_wp) > 0.0_wp)
   end subroutine test_sapwood_conductance

   !----- 4. ED2's root-fraction profile, beta = 0.96 to a 2 m rooting depth. -----------------------!
   subroutine test_root_profile()
      real(wp) :: f1, f2, f3
      print '(a)', '-- Root profile --'
      f1 = root_fraction_profile(0.96_wp, 2.0_wp, 0.0_wp, 0.5_wp)   ! equal-thickness (0.5 m) layers
      f2 = root_fraction_profile(0.96_wp, 2.0_wp, 0.5_wp, 1.0_wp)
      f3 = root_fraction_profile(0.96_wp, 2.0_wp, 1.0_wp, 1.5_wp)
      call check_true('shallower (equal-thickness) layer holds more roots', f1 > f2 .and. f2 > f3)
      call check('full column [0,D] telescopes to 1-beta',                                       &
           root_fraction_profile(0.96_wp,2.0_wp,0.0_wp,1.0_wp) +                                  &
           root_fraction_profile(0.96_wp,2.0_wp,1.0_wp,2.0_wp), 1.0_wp-0.96_wp, 1.0e-12_wp)
      call check_true('below rooting depth => 0 roots',                                          &
           root_fraction_profile(0.96_wp,2.0_wp,2.5_wp,3.0_wp) == 0.0_wp)
   end subroutine test_root_profile

   !----- 5. Gardner's single-root conductance is per root length: 0.01 plants/m2 with 1000 m of root  !
   !      each in a 0.4 m layer, against the same density in two 0.2 m halves, each plant with 500 m  !
   !      per half at the same half-spacing; and twice the plants at that spacing. -----------------!
   subroutine test_gardner()
      real(wp) :: lv, half, g, g_half
      print '(a)', '-- Gardner conductance --'
      lv   = 0.01_wp * 1000.0_wp / 0.4_wp
      half = 1.0_wp / sqrt(pi * lv)
      g      = rhizosphere_conductance(1.0_wp, 1000.0_wp, half, 1.5e-4_wp)
      g_half = rhizosphere_conductance(1.0_wp,  500.0_wp, half, 1.5e-4_wp)
      call check('a layer split in two keeps its conductance', 2.0_wp * g_half, g, 1.0e-14_wp * g)
      call check('twice the plants at the same spacing, twice the conductance',                    &
                 rhizosphere_conductance(2.0_wp, 1000.0_wp, half, 1.5e-4_wp), 2.0_wp * g, 1.0e-14_wp * g)
   end subroutine test_gardner

   !=======================================================================================!
   ! Biomass fast/slow SEAM: mass, not psi, is conserved across a slow-loop biomass update, so a !
   ! day of wood growth at FIXED water mass reads as a drier wood, which then draws more water   !
   ! from the soil (the network's uptake follows its potential). --------------------------------!
   !=======================================================================================!
   subroutine test_growth_seam()
      type(water_curve_t), parameter :: WOOD = water_curve_t(pi0 = -1.0_wp, elastic_mod = 8.0_wp,              &
                                                             apoplast_frac = 0.20_wp, water_sat = 1.0_wp)
      real(wp) :: bwood0, bwood1, psi_wood0, psi_wood1, w0
      print '(a)', '-- Biomass fast/slow seam: mass-conserving growth --'
      bwood0    = 7.0_wp
      psi_wood0 = tissue_rwc_to_psi(0.90_wp, WOOD%pi0, WOOD%elastic_mod)
      w0        = tissue_psi_to_water(psi_wood0, WOOD, bwood0)
      bwood1    = bwood0 * 1.05_wp
      psi_wood1 = tissue_water_to_psi(w0, WOOD, bwood1)
      call check_true('growth at conserved mass lowers psi_wood (drier)', psi_wood1 < psi_wood0)
   end subroutine test_growth_seam

   !=======================================================================================!
   ! Saturation-CEILING clamp (the guard mass conservation needs): a persisted mass can only       !
   ! exceed the current capacity after biomass SHRINKS (e.g. the phenology dormant-canopy leaf       !
   ! snap-to-bare), never after it grows. tissue_capped_water must be an exact no-op on growth,       !
   ! must cap exactly at the new ceiling on a shrink, and must reduce a fully-bare tissue to exactly  !
   ! 0 (no NaN/negative) -- see ed_plant_hydr_lib.f90.                                                 !
   !=======================================================================================!
   subroutine test_seam_capacity_clamp()
      real(wp), parameter :: PI0 = -1.0_wp, EPS = 8.0_wp, AF = 0.20_wp, WS = 1.0_wp
      type(water_curve_t), parameter :: CURVE = water_curve_t(PI0, EPS, AF, WS)
      real(wp) :: biomass0, biomass_shrunk, biomass_grown, w0, ceiling_shrunk, w_capped
      print '(a)', '-- Biomass fast/slow seam: saturation-ceiling clamp --'
      biomass0 = 7.0_wp
      w0 = tissue_psi_to_water(tissue_rwc_to_psi(0.97_wp, PI0, EPS), CURVE, biomass0)   ! near-saturated

      !----- (a) GROWING biomass: the clamp must be a strict no-op. ---------------------------!
      biomass_grown = biomass0 * 1.05_wp
      call check('clamp is a no-op when biomass grows', tissue_capped_water(w0, CURVE, biomass_grown), &
                 w0, 1.0e-12_wp)

      !----- (b) SHRINKING biomass (a snap-to-bare-like 10x drop): the carried-forward mass would   !
      !          now be wildly supersaturated (w0 above the new ceiling) without the clamp. --------!
      biomass_shrunk = biomass0 * 0.10_wp
      ceiling_shrunk = WS * biomass_shrunk
      call check_true('unclamped carried-forward mass would exceed the new ceiling', w0 > ceiling_shrunk)
      w_capped = tissue_capped_water(w0, CURVE, biomass_shrunk)
      call check('clamp caps mass exactly at the new ceiling', w_capped, ceiling_shrunk, 1.0e-12_wp)
      call check_true('released excess is non-negative', (w0 - w_capped) >= 0.0_wp)
      !----- The clamped mass reads back at exactly rwc=1 (sane, bounded), not the pathological     !
      !      super-turgid extrapolation the unclamped mass would imply. ----------------------------!
      call check('clamped mass diagnoses psi at exactly rwc=1',                                      &
                 tissue_water_to_psi(w_capped, CURVE, biomass_shrunk),                  &
                 tissue_rwc_to_psi(1.0_wp, PI0, EPS), 1.0e-9_wp)

      !----- (c) Full snap-to-bare (biomass -> 0): clamps to exactly 0, not NaN/negative. ----------!
      w_capped = tissue_capped_water(w0, CURVE, 0.0_wp)
      call check('fully-bare tissue clamps to exactly 0', w_capped, 0.0_wp, 0.0_wp)
   end subroutine test_seam_capacity_clamp

end program test_plant_hydr_lib
