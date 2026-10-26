! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_canopy_radiation -- the canopy two-stream (meds_canopy_radiation;                   !
! docs/science/canopy_radiation_transfer.md).                                              !
!                                                                                          !
! Checks the physical invariants that a correct two-stream must satisfy (canopy_radiation,   !
! ground_optics, and layer_beam_sources for a thin layer's scattering):                     !
!   1. ENERGY CONSERVATION  (per shortwave band): absorbed + net-ground + reflected = incident.!
!   2. BEER'S LAW           : a single black cohort under a direct beam attenuates as exp(-k L).!
!   3. ZERO-LAI ALBEDO      : an empty canopy has albedo = ground albedo.                       !
!   4. MULTI-PFT PATCH       : energy conserves, and each cohort's diffuse light follows its OWN    !
!                              PFT's optics, not the top cohort's (the guard against ED2 bug B1 of  !
!                              docs/ed2_comparison.md).                                              !
!   5. LONGWAVE EQUILIBRIUM  : an isothermal canopy+ground+sky has ~zero net thermal absorption.    !
!   6. COHORT ORDER          : the upper of two equal cohorts absorbs more and sees the top light.  !
!   7. TIED COHORTS          : equal-height cohorts share one layer, whatever their order.          !
!   8. CUTTING              : one LAI-5 cohort and five LAI-1 cohorts of a PFT give the same light; !
!                             a thin layer's scattered beam is the single-scatter split.            !
!   9. TISSUE LONGWAVE      : leaves and wood, or two tied cohorts, exchange longwave reciprocally,  !
!                             each charged its own emission.                                        !
!  10. NO COHORTS           : a patch with no cohort (ncoh = 0): the ground takes what its albedo   !
!                             leaves of the shortwave and its emissivity's share of the sky less its !
!                             emission, the canopy nothing, every band closes, and the light is a   !
!                             cohort of no leaf or wood area's.                                     !
! The leaf-angle distribution and its moments are test_optics_lib's.                       !
!==========================================================================================!
program test_canopy_radiation
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : stefan_boltzmann
   use meds_canopy_radiation, only : rad_env_t, rad_out_t, alloc_rad_env, ground_optics_env_t
   use meds_pft_params,   only : rad_pft_optics_t, RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   use meds_optics_lib,     only : gfun_direct
   use meds_pft_derive,     only : derive_rad_optics, moments_to_beta_params
   use meds_canopy_radiation, only : canopy_radiation, ground_optics, layer_beam_sources
   implicit none

   real(wp), parameter :: TOL = 1.0e-9_wp

   call test_energy_conservation()
   call test_beers_law()
   call test_zero_lai_albedo()
   call test_multi_pft()
   call test_longwave_equilibrium()
   call test_cohort_order()
   call test_tied_cohorts_codominant()
   call test_cutting_invariance()
   call test_tissue_longwave()
   call test_no_cohorts()

   call test_report('test_canopy_radiation')

contains

   !----- Build a 2-PFT, 3-band optics table (VIS/NIR/LW) with a chosen leaf-angle mean. ------!
   subroutine build_optics(mean_deg, opt)
      real(wp),               intent(in)  :: mean_deg
      type(rad_pft_optics_t), intent(out) :: opt
      integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT, NP = 2_ik
      real(wp) :: rl(NB,NP), tl(NB,NP), rw(NB,NP), tw(NB,NP)
      real(wp) :: cl(NP), cw(NP), bp(NP), bq(NP)
      logical  :: hb(NB), he(NB)
      integer(ik) :: p
      real(wp) :: emiss_leaf, emiss_wood
      !----- Leaf spectral: PFT 1 darker (low VIS reflectance), PFT 2 brighter. -------------!
      rl(RAD_VIS,:) = [0.09_wp, 0.14_wp] ; tl(RAD_VIS,:) = [0.05_wp, 0.08_wp]
      rl(RAD_NIR,:) = [0.45_wp, 0.50_wp] ; tl(RAD_NIR,:) = [0.25_wp, 0.30_wp]
      emiss_leaf = 0.97_wp
      rl(RAD_LW,:)  = 1.0_wp - emiss_leaf ; tl(RAD_LW,:) = 0.0_wp
      !----- Wood spectral (near-opaque). ---------------------------------------------------!
      rw(RAD_VIS,:) = 0.11_wp ; tw(RAD_VIS,:) = 0.001_wp
      rw(RAD_NIR,:) = 0.25_wp ; tw(RAD_NIR,:) = 0.001_wp
      emiss_wood = 0.90_wp
      rw(RAD_LW,:)  = 1.0_wp - emiss_wood ; tw(RAD_LW,:) = 0.0_wp
      cl = [0.80_wp, 0.85_wp] ; cw = [0.50_wp, 0.50_wp]
      do p = 1_ik, NP
         call moments_to_beta_params(mean_deg, 20.0_wp, bp(p), bq(p))
      end do
      hb = [.true.,  .true.,  .false.]        ! VIS/NIR have a beam; LW does not
      he = [.false., .false., .true. ]        ! only LW emits
      call derive_rad_optics(NB, NP, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, opt)
   end subroutine build_optics

   !----- Fill shortwave + longwave forcing from a surface state. ----------------------------!
   subroutine build_forcing(cosz, t_sky, t_soil, f)
      real(wp),            intent(in)  :: cosz, t_sky, t_soil
      type(rad_env_t), intent(out) :: f
      type(ground_optics_env_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)
      call alloc_rad_env(f, N_RAD_BAND_DEFAULT)
      f%cosz = cosz
      f%incid_beam(RAD_VIS) = 300.0_wp ; f%incid_diff(RAD_VIS) = 100.0_wp
      f%incid_beam(RAD_NIR) = 250.0_wp ; f%incid_diff(RAD_NIR) =  80.0_wp
      f%incid_beam(RAD_LW)  = 0.0_wp   ; f%incid_diff(RAD_LW)  = stefan_boltzmann * t_sky ** 4
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = [0.15_wp, 0.30_wp, 0.0_wp]
      surf%soil_emiss  = 0.95_wp ; surf%soil_temp = t_soil
      he = [.false., .false., .true.]
      call ground_optics(surf, N_RAD_BAND_DEFAULT, he, f%grnd_refl, f%grnd_emiss)
   end subroutine build_forcing

   !=======================================================================================!
   subroutine test_energy_conservation()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 4_ik
      integer(ik) :: pft(NC), b
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC), absorbed, net_ground, reflected, incid, resid
      integer :: ih
      print '(a)', '[1] shortwave energy conservation (4 cohorts)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      pft  = [1_ik, 2_ik, 1_ik, 2_ik]
      lai  = [2.0_wp, 1.5_wp, 1.0_wp, 0.5_wp]   ! bottom -> top
      wai  = 0.1_wp * lai
      tcan = 298.0_wp
      hgt = [(2.0_wp*real(ih,wp), ih = 1, NC)]      ! distinct heights: no ties, per-cohort layers
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      do b = RAD_VIS, RAD_NIR
         absorbed   = sum(flux%abs_leaf(b,:)) + sum(flux%abs_wood(b,:))
         net_ground = flux%dn_ground(b) - flux%up_ground(b)
         incid      = f%incid_beam(b) + f%incid_diff(b)
         reflected  = flux%albedo(b) * incid
         resid      = incid - (absorbed + net_ground + reflected)
         call check('energy residual band', resid, 0.0_wp, 1.0e-8_wp * incid)
         if (flux%albedo(b) < 0.0_wp .or. flux%albedo(b) > 1.0_wp) then
            call check_true('albedo within [0,1]', .false., flux%albedo(b))
         end if
      end do
   end subroutine test_energy_conservation

   !=======================================================================================!
   !  Cohort ORDER: two EQUAL-LAI cohorts stacked BOTTOM(1)->TOP(2). The top cohort intercepts !
   !  the incoming beam+diffuse first, so per band it must absorb strictly MORE than the shaded  !
   !  bottom cohort. This vertical-gradient property is exactly what the fast loop's              !
   !  apply_rt_forcing `perm` relies on to route the higher absorbed PAR to the taller cohort.     !
   subroutine test_cohort_order()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 2_ik
      integer(ik) :: pft(NC)
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC)
      integer :: ih
      print '(a)', '[6] cohort order: top of two equal-LAI cohorts absorbs more than the bottom'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      pft = [1_ik, 1_ik] ; lai = [2.0_wp, 2.0_wp] ; wai = [0.0_wp, 0.0_wp] ; tcan = 298.0_wp  ! idx1=bottom, idx2=top
      hgt = [(2.0_wp*real(ih,wp), ih = 1, NC)]      ! distinct heights: no ties, per-cohort layers
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      call check_true('VIS: top cohort (idx2) absorbs more than bottom (idx1)',                 &
                      flux%abs_leaf(RAD_VIS,2) > flux%abs_leaf(RAD_VIS,1),                       &
                      flux%abs_leaf(RAD_VIS,2) - flux%abs_leaf(RAD_VIS,1))
      call check_true('NIR: top cohort (idx2) absorbs more than bottom (idx1)',                 &
                      flux%abs_leaf(RAD_NIR,2) > flux%abs_leaf(RAD_NIR,1),                       &
                      flux%abs_leaf(RAD_NIR,2) - flux%abs_leaf(RAD_NIR,1))
      !----- The light REACHING each cohort's top (the phenology PAR cue's driver): the top     !
      !      cohort sees exactly the canopy-top incident, the one below it less. ---------------!
      call check('VIS reaching the top cohort = canopy-top incident', flux%incid_top(RAD_VIS,2), &
                 f%incid_beam(RAD_VIS) + f%incid_diff(RAD_VIS), 1.0e-12_wp)
      call check_true('VIS reaching the lower cohort is less than at the top',                   &
                      flux%incid_top(RAD_VIS,1) < flux%incid_top(RAD_VIS,2),                     &
                      flux%incid_top(RAD_VIS,1) / flux%incid_top(RAD_VIS,2))
   end subroutine test_cohort_order

   !=======================================================================================!
   !  EXACTLY-TIED cohorts are CO-DOMINANT, not stacked: cohorts of equal height share one RT      !
   !  layer (docs/science/canopy_radiation_transfer.md §4, eq. (14)), as `update_overtopping_lai`  !
   !  treats them. Recruits tie BIT-EXACTLY -- one scalar `recruit_dbh` for every PFT -- and the   !
   !  sort is stable, so stacking them would put the lower PFT INDEX on top: the lower slot then   !
   !  absorbs 12.6 % less shortwave at zero overstory, rising to 15.1 % under LAI 6.               !
   !                                                                                               !
   !  Asserted two ways, because either alone is weak. Equal absorption for optically identical    !
   !  twins could be met by a broken solver returning zeros; permutation invariance could be met   !
   !  by ignoring PFT. Together they pin it.                                                       !
   !=======================================================================================!
   subroutine test_tied_cohorts_codominant()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 3_ik
      integer(ik) :: pft(NC)
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC)
      real(wp) :: a_pft1, a_pft2, b_pft1, b_pft2, total_a, total_b
      print '(a)', '[7] exactly-tied cohorts are co-dominant, not stacked'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      !----- slots 1,2 tie exactly; slot 3 is a real overstory, the regime where the bias is worst. -!
      hgt  = [2.0_wp, 2.0_wp, 20.0_wp]
      lai  = [0.25_wp, 0.25_wp, 4.0_wp]
      wai  = [0.02_wp, 0.02_wp, 0.40_wp]
      tcan = 298.0_wp

      !----- (a) two OPTICALLY IDENTICAL tied cohorts must absorb equally. Same PFT, so any         !
      !      difference is positional and nothing else. -----------------------------------------!
      pft = [1_ik, 1_ik, 1_ik]
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      call check('tied twins absorb equally (VIS)', flux%abs_leaf(RAD_VIS,1),                    &
                 flux%abs_leaf(RAD_VIS,2), 1.0e-12_wp * max(flux%abs_leaf(RAD_VIS,2), 1.0e-12_wp))
      call check('tied twins absorb equally (NIR)', flux%abs_leaf(RAD_NIR,1),                    &
                 flux%abs_leaf(RAD_NIR,2), 1.0e-12_wp * max(flux%abs_leaf(RAD_NIR,2), 1.0e-12_wp))
      call check_true('tied twins absorb a non-trivial amount', flux%abs_leaf(RAD_NIR,1) > 1.0_wp, &
                      flux%abs_leaf(RAD_NIR,1))

      !----- (b) PERMUTING the PFT order of the tied pair must not move either one's absorption.    !
      !      This is the defect it guards: creation order is `do pf = 1, n_pft`, so a stacking sort  !
      !      hands PFT 1 the upper slot at every recruitment event. -------------------------------!
      pft = [1_ik, 2_ik, 1_ik]
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      a_pft1 = flux%abs_leaf(RAD_NIR,1) ; a_pft2 = flux%abs_leaf(RAD_NIR,2)
      total_a = a_pft1 + a_pft2
      pft = [2_ik, 1_ik, 1_ik]
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      b_pft2 = flux%abs_leaf(RAD_NIR,1) ; b_pft1 = flux%abs_leaf(RAD_NIR,2)
      total_b = b_pft2 + b_pft1
      call check('permuting tied PFTs leaves PFT 1 unchanged', a_pft1, b_pft1,                   &
                 1.0e-12_wp * max(a_pft1, 1.0e-12_wp))
      call check('permuting tied PFTs leaves PFT 2 unchanged', a_pft2, b_pft2,                   &
                 1.0e-12_wp * max(a_pft2, 1.0e-12_wp))
      call check('permuting tied PFTs conserves the pair total', total_a, total_b,               &
                 1.0e-12_wp * max(total_a, 1.0e-12_wp))
   end subroutine test_tied_cohorts_codominant

   !=======================================================================================!
   !  CUTTING: a uniform canopy's light does not depend on how it is cut into cohorts, since     !
   !  each layer's scattered beam is solved through its depth. And a thin layer's scattered beam !
   !  is the single-scatter split of what it intercepts.                                        !
   !=======================================================================================!
   subroutine test_cutting_invariance()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)        :: f
      type(rad_out_t)        :: one, five
      integer(ik) :: b, ih
      integer(ik) :: pft1(1), pft5(5)
      real(wp)    :: hgt1(1), hgt5(5), lai1(1), lai5(5), wai1(1), wai5(5), t1(1), t5(5)
      real(wp)    :: up, dn, intercept
      print '(a)', '[8] cutting invariance (one LAI-5 cohort against five LAI-1 cohorts)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.87_wp, 290.0_wp, 298.0_wp, f)
      pft1 = 1_ik ; hgt1 = 20.0_wp ; lai1 = 5.0_wp ; wai1 = 0.5_wp ; t1 = 298.0_wp
      pft5 = 1_ik ; lai5 = 1.0_wp  ; wai5 = 0.1_wp ; t5 = 298.0_wp
      hgt5 = [(10.0_wp + 2.0_wp * real(ih, wp), ih = 1, 5)]     ! distinct heights: five layers
      call canopy_radiation(opt, f, 1_ik, pft1, hgt1, lai1, wai1, t1, t1, one)
      call canopy_radiation(opt, f, 5_ik, pft5, hgt5, lai5, wai5, t5, t5, five)
      do b = RAD_VIS, RAD_LW
         call check('albedo, one cohort or five', five%albedo(b), one%albedo(b), TOL)
         call check('light reaching the ground, one cohort or five', five%dn_ground(b), one%dn_ground(b),  &
                    TOL * 1.0e3_wp)
         call check('absorbed, one cohort or five', sum(five%abs_leaf(b,:)) + sum(five%abs_wood(b,:)),     &
                    one%abs_leaf(b,1) + one%abs_wood(b,1), TOL * 1.0e3_wp)
      end do
      call layer_beam_sources(0.15_wp, 0.5_wp, 0.4_wp, 0.6_wp, 1.0e-6_wp, 100.0_wp, up, dn)
      intercept = 100.0_wp * (1.0_wp - exp(-0.6_wp * 1.0e-6_wp))
      call check('thin layer: up = omega beta0 x intercepted', up, 0.15_wp * 0.4_wp * intercept,           &
                 1.0e-4_wp * up)
      call check('thin layer: down = omega (1 - beta0) x intercepted', dn, 0.15_wp * 0.6_wp * intercept,  &
                 1.0e-4_wp * dn)
   end subroutine test_cutting_invariance

   !=======================================================================================!
   subroutine test_beers_law()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 1_ik
      integer(ik) :: pft(NC)
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC), kdir, gee, absorbed_vis, expect
      integer :: ih
      real(wp) :: elai
      print '(a)', '[2] Beers law (single black cohort, direct beam)'
      call build_optics(45.0_wp, opt)
      call build_forcing(1.0_wp, 290.0_wp, 298.0_wp, f)
      !----- Black leaves: force VIS single-scatter albedo to ~0 by zeroing spectra. --------!
      opt%omega_leaf(RAD_VIS,:) = 0.0_wp ; opt%g_leaf(RAD_VIS,:) = 0.0_wp
      opt%omega_wood(RAD_VIS,:) = 0.0_wp ; opt%g_wood(RAD_VIS,:) = 0.0_wp
      f%incid_diff(RAD_VIS) = 0.0_wp                       ! beam only
      f%grnd_refl(RAD_VIS)  = 0.0_wp                       ! black soil: no reflection to re-absorb
      pft = [1_ik] ; lai = [3.0_wp] ; wai = [0.0_wp] ; tcan = 298.0_wp
      hgt = [(2.0_wp*real(ih,wp), ih = 1, NC)]      ! distinct heights: no ties, per-cohort layers
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      elai         = opt%clumping_leaf(1) * lai(1)
      gee          = gfun_direct(opt%lidf(:,1), 1.0_wp)
      kdir         = gee / 1.0_wp
      expect       = f%incid_beam(RAD_VIS) * (1.0_wp - exp(-kdir * elai))
      absorbed_vis = flux%abs_leaf(RAD_VIS,1) + flux%abs_wood(RAD_VIS,1)
      call check('absorbed = incid*(1-exp(-kL))', absorbed_vis, expect, 1.0e-7_wp)
      call check('below-canopy = incid*exp(-kL)', flux%dn_ground(RAD_VIS),                     &
                 f%incid_beam(RAD_VIS) * exp(-kdir * elai), 1.0e-7_wp)
   end subroutine test_beers_law

   !=======================================================================================!
   subroutine test_zero_lai_albedo()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 1_ik
      integer(ik) :: pft(NC)
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC)
      integer :: ih
      print '(a)', '[3] zero-LAI canopy has albedo = ground albedo'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      pft = [1_ik] ; lai = [0.0_wp] ; wai = [0.0_wp] ; tcan = 298.0_wp
      hgt = [(2.0_wp*real(ih,wp), ih = 1, NC)]      ! distinct heights: no ties, per-cohort layers
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      call check('albedo(VIS) = soil albedo', flux%albedo(RAD_VIS), f%grnd_refl(RAD_VIS), 1.0e-9_wp)
      call check('albedo(NIR) = soil albedo', flux%albedo(RAD_NIR), f%grnd_refl(RAD_NIR), 1.0e-9_wp)
   end subroutine test_zero_lai_albedo

   !=======================================================================================!
   subroutine test_multi_pft()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 2_ik
      integer(ik) :: pft(NC)
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC), absorbed, net_ground, reflected, incid, resid
      real(wp) :: a_own1, a_own2, a_top1, a_top2
      integer :: ih
      print '(a)', '[4] multi-PFT patch: conservation, and each cohort on its own optics (ED2 bug B1)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.6_wp, 290.0_wp, 298.0_wp, f)
      pft = [1_ik, 2_ik] ; lai = [1.5_wp, 1.5_wp] ; wai = 0.15_wp ; tcan = 298.0_wp
      call check_true('the two PFTs differ in VIS omega',                                              &
                      abs(opt%omega_leaf(RAD_VIS,1) - opt%omega_leaf(RAD_VIS,2)) > 1.0e-3_wp,           &
                      opt%omega_leaf(RAD_VIS,2) - opt%omega_leaf(RAD_VIS,1))
      hgt = [(2.0_wp*real(ih,wp), ih = 1, NC)]      ! distinct heights: no ties, per-cohort layers
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      incid      = f%incid_beam(RAD_VIS) + f%incid_diff(RAD_VIS)
      absorbed   = sum(flux%abs_leaf(RAD_VIS,:)) + sum(flux%abs_wood(RAD_VIS,:))
      net_ground = flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS)
      reflected  = flux%albedo(RAD_VIS) * incid
      resid      = incid - (absorbed + net_ground + reflected)
      call check('energy residual (2 PFT)', resid, 0.0_wp, 1.0e-8_wp * incid)

      !----- ED2 bug B1 gave every cohort's DIFFUSE light the top cohort's PFT optics. Under      !
      !      diffuse light only, cohort 1 (2 m) sits below cohort 2 (4 m). (a) Changing only the   !
      !      lower cohort's PFT must change its absorption: under B1 it would not. (b) Changing   !
      !      only the PFT of a top cohort too thin to shade must leave the lower cohort's          !
      !      absorption alone: under B1 it would move with the top's optics. Conservation alone   !
      !      cannot see this, since B1 conserves energy. ---------------------------------------!
      f%incid_beam = 0.0_wp
      lai = [1.5_wp, 0.5_wp] ; wai = [0.15_wp, 0.05_wp]
      call canopy_radiation(opt, f, NC, [1_ik, 1_ik], hgt, lai, wai, tcan, tcan, flux)
      a_own1 = flux%abs_leaf(RAD_VIS,1)
      call canopy_radiation(opt, f, NC, [2_ik, 1_ik], hgt, lai, wai, tcan, tcan, flux)
      a_own2 = flux%abs_leaf(RAD_VIS,1)
      call check_true('B1: the lower cohort''s diffuse absorption follows its own PFT',                   &
                      abs(a_own2 - a_own1) > 1.0e-3_wp * a_own1, (a_own2 - a_own1) / a_own1)
      lai = [1.5_wp, 1.0e-6_wp] ; wai = [0.15_wp, 0.0_wp]
      call canopy_radiation(opt, f, NC, [1_ik, 1_ik], hgt, lai, wai, tcan, tcan, flux)
      a_top1 = flux%abs_leaf(RAD_VIS,1)
      call canopy_radiation(opt, f, NC, [1_ik, 2_ik], hgt, lai, wai, tcan, tcan, flux)
      a_top2 = flux%abs_leaf(RAD_VIS,1)
      call check_true('B1: a thin top cohort''s PFT does not set the lower cohort''s optics',              &
                      abs(a_top2 - a_top1) < 1.0e-4_wp * a_top1, (a_top2 - a_top1) / a_top1)
   end subroutine test_multi_pft

   !=======================================================================================!
   subroutine test_longwave_equilibrium()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)    :: f
      type(rad_out_t)       :: flux
      integer(ik), parameter :: NC = 3_ik
      integer(ik) :: pft(NC), i
      real(wp) :: hgt(NC), lai(NC), wai(NC), tcan(NC), net_lw, maxnet, max_tissue
      integer :: ih
      real(wp), parameter :: TEQ = 300.0_wp
      print '(a)', '[5] longwave radiative equilibrium (isothermal -> zero net absorption)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.5_wp, TEQ, TEQ, f)          ! sky and soil both at Teq
      pft  = [1_ik, 2_ik, 1_ik] ; lai = [2.0_wp, 1.0_wp, 0.5_wp] ; wai = 0.1_wp * lai
      tcan = TEQ                                        ! canopy also at Teq
      hgt = [(2.0_wp*real(ih,wp), ih = 1, NC)]      ! distinct heights: no ties, per-cohort layers
      call canopy_radiation(opt, f, NC, pft, hgt, lai, wai, tcan, tcan, flux)
      maxnet = 0.0_wp
      max_tissue = 0.0_wp
      do i = 1_ik, NC
         net_lw = flux%abs_leaf(RAD_LW,i) + flux%abs_wood(RAD_LW,i)
         maxnet = max(maxnet, abs(net_lw))
         max_tissue = max(max_tissue, abs(flux%abs_leaf(RAD_LW,i)), abs(flux%abs_wood(RAD_LW,i)))
      end do
      call check('max |net LW| at equilibrium', maxnet, 0.0_wp, 1.0e-7_wp * stefan_boltzmann * TEQ ** 4)
      call check('max |net LW| of a tissue at equilibrium', max_tissue, 0.0_wp, 1.0e-7_wp * stefan_boltzmann * TEQ ** 4)
   end subroutine test_longwave_equilibrium

   !=======================================================================================!
   !  TISSUE LONGWAVE: a cohort's leaves and wood meet the same radiation, share what it takes  !
   !  in by leaf_frac and are each charged their own emission, as tied cohorts share their    !
   !  layer's (docs/science/canopy_radiation_transfer.md §4, eqs. (14)-(15)). With sky and     !
   !  ground at the leaves' temperature, wood 10 K warmer warms the leaves, and by reciprocity !
   !  it gives them exactly what leaves 10 K warmer give the wood; likewise between two tied   !
   !  cohorts. Charged a share of the cohort's net instead, the cool leaves lost to the warm   !
   !  wood.                                                                                    !
   !=======================================================================================!
   subroutine test_tissue_longwave()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)        :: f
      type(rad_out_t)        :: hot_wood, hot_leaf, hot_a, hot_b
      integer(ik) :: pft1(1), pft2(2)
      real(wp)    :: hgt1(1), lai1(1), wai1(1), hgt2(2), lai2(2), wai2(2)
      real(wp)    :: net, gain_a, gain_b
      real(wp), parameter :: T_COOL = 295.0_wp, T_WARM = 305.0_wp
      print '(a)', '[9] tissue longwave: each tissue charged its own emission'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.5_wp, T_COOL, T_COOL, f)          ! sky and soil at the cool temperature
      pft1 = 1_ik ; hgt1 = 10.0_wp ; lai1 = 3.0_wp ; wai1 = 0.6_wp
      call canopy_radiation(opt, f, 1_ik, pft1, hgt1, lai1, wai1, [T_COOL], [T_WARM], hot_wood)
      call canopy_radiation(opt, f, 1_ik, pft1, hgt1, lai1, wai1, [T_WARM], [T_COOL], hot_leaf)
      call check_true('warm wood warms its cool leaves', hot_wood%abs_leaf(RAD_LW,1) > 0.0_wp,            &
                      hot_wood%abs_leaf(RAD_LW,1))
      call check_true('warm wood loses longwave', hot_wood%abs_wood(RAD_LW,1) < 0.0_wp,                   &
                      hot_wood%abs_wood(RAD_LW,1))
      call check('leaves gain from warm wood what wood gains from warm leaves', hot_wood%abs_leaf(RAD_LW,1), &
                 hot_leaf%abs_wood(RAD_LW,1), 1.0e-9_wp * abs(hot_leaf%abs_wood(RAD_LW,1)))
      !----- The tissues' nets add up to the band's: in at the top less out, less what the ground   !
      !      takes. --------------------------------------------------------------------------------!
      net = f%incid_diff(RAD_LW) * (1.0_wp - hot_wood%albedo(RAD_LW))                                 &
            - (hot_wood%dn_ground(RAD_LW) - hot_wood%up_ground(RAD_LW))
      call check('leaf + wood longwave = the band''s in - out', hot_wood%abs_leaf(RAD_LW,1)                &
                 + hot_wood%abs_wood(RAD_LW,1), net, 1.0e-9_wp * stefan_boltzmann * T_WARM ** 4)

      !----- Two tied cohorts of different PFT and area: the cool one gains from the warm one   !
      !      exactly what it would give were the temperatures swapped. -------------------------!
      pft2 = [1_ik, 2_ik] ; hgt2 = 4.0_wp ; lai2 = [1.0_wp, 2.5_wp] ; wai2 = [0.3_wp, 0.2_wp]
      call canopy_radiation(opt, f, 2_ik, pft2, hgt2, lai2, wai2, [T_WARM, T_COOL], [T_WARM, T_COOL], hot_a)
      call canopy_radiation(opt, f, 2_ik, pft2, hgt2, lai2, wai2, [T_COOL, T_WARM], [T_COOL, T_WARM], hot_b)
      gain_b = hot_a%abs_leaf(RAD_LW,2) + hot_a%abs_wood(RAD_LW,2)
      gain_a = hot_b%abs_leaf(RAD_LW,1) + hot_b%abs_wood(RAD_LW,1)
      call check_true('a cool tied cohort gains from a warm one', gain_b > 0.0_wp, gain_b)
      call check('tied cohorts exchange reciprocally', gain_b, gain_a, 1.0e-9_wp * abs(gain_a))
   end subroutine test_tissue_longwave

   !=======================================================================================!
   !  NO COHORTS (canopy_radiation, ncoh = 0): a bare patch under the sun at cos(z) = 0.7, a 280 K  !
   !  sky and 295 K soil of albedo 0.15 / 0.30 and emissivity 0.95. The ground absorbs (1 - albedo) !
   !  of each shortwave band and, by Kirchhoff, eps (L_sky - sigma T^4) of the longwave; the canopy  !
   !  has nothing to absorb; what reaches the ground and what it returns close each band. A single !
   !  cohort of no leaf or wood area must give the same light, so the branch is the canopy's limit.  !
   !=======================================================================================!
   subroutine test_no_cohorts()
      type(rad_pft_optics_t) :: opt
      type(rad_env_t)        :: f
      type(rad_out_t)        :: bare, empty
      integer(ik) :: pft0(0), b
      real(wp)    :: r0(0), incid, l_sky
      real(wp), parameter :: T_SKY = 280.0_wp, T_SOIL = 295.0_wp, EMISS = 0.95_wp
      real(wp), parameter :: ALBEDO(2) = [0.15_wp, 0.30_wp]   ! VIS, NIR, as build_forcing sets them
      print '(a)', '[10] no cohorts: the ground takes it all'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, T_SKY, T_SOIL, f)
      call canopy_radiation(opt, f, 0_ik, pft0, r0, r0, r0, r0, r0, bare)
      call check_true('no cohorts: the canopy has no tissue to absorb anything',                                   &
                      bare%n_coh == 0_ik .and. size(bare%abs_leaf, 2) == 0 .and. size(bare%abs_wood, 2) == 0)
      do b = RAD_VIS, RAD_NIR
         incid = f%incid_beam(b) + f%incid_diff(b)
         call check('no cohorts: all the shortwave reaches the ground', bare%dn_ground(b), incid, 1.0e-12_wp * incid)
         call check('no cohorts: the ground absorbs 1 - albedo of it', bare%dn_ground(b) - bare%up_ground(b),       &
                    (1.0_wp - ALBEDO(b)) * incid, 1.0e-12_wp * incid)
         call check('no cohorts: the patch''s albedo is the ground''s', bare%albedo(b), ALBEDO(b), 1.0e-15_wp)
         call check('no cohorts: incident = absorbed by the ground + reflected', incid,                            &
                    (bare%dn_ground(b) - bare%up_ground(b)) + bare%albedo(b) * incid, 1.0e-12_wp * incid)
      end do
      l_sky = stefan_boltzmann * T_SKY ** 4
      call check('no cohorts: the ground nets eps (L_sky - sigma T^4) of the longwave (Kirchhoff)',                  &
                 bare%dn_ground(RAD_LW) - bare%up_ground(RAD_LW), EMISS * (l_sky - stefan_boltzmann * T_SOIL ** 4),  &
                 1.0e-12_wp * l_sky)
      call check('no cohorts: the sky''s longwave in = the ground''s net + what leaves the top',                    &
                 l_sky, (bare%dn_ground(RAD_LW) - bare%up_ground(RAD_LW)) + bare%albedo(RAD_LW) * l_sky,          &
                 1.0e-12_wp * l_sky)
      !----- The canopy's limit: one cohort with no leaf and no wood. ---------------------------------!
      call canopy_radiation(opt, f, 1_ik, [1_ik], [5.0_wp], [0.0_wp], [0.0_wp], [T_SOIL], [T_SOIL], empty)
      do b = RAD_VIS, RAD_LW
         call check('no cohorts = a cohort with no leaf or wood area: the light reaching the ground',               &
                    bare%dn_ground(b), empty%dn_ground(b), 1.0e-12_wp * max(bare%dn_ground(b), 1.0_wp))
         call check('no cohorts = a cohort with no leaf or wood area: the light leaving the ground',                &
                    bare%up_ground(b), empty%up_ground(b), 1.0e-12_wp * max(bare%up_ground(b), 1.0_wp))
         call check('no cohorts = a cohort with no leaf or wood area: the albedo', bare%albedo(b), empty%albedo(b),  &
                    1.0e-12_wp)
      end do
   end subroutine test_no_cohorts

end program test_canopy_radiation
