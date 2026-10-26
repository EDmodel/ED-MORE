! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_respiration -- the respiration of ed_respiration: the plants' maintenance           !
! respiration (docs/science/plant_respiration.md eqs. (3)-(5)) and the soil's              !
! heterotrophic respiration from its carbon pools (docs/science/soil_carbon.md §3-§4).     !
!   1. 25 C (root_zone_temp_scale, fine_root_maintenance_respiration,                      !
!      stem_maintenance_respiration): the temperature response is 1 at 25 C, so the roots  !
!      respire their rate times their carbon and the sapwood its rate times its volume.    !
!   2. LAYERS (root_zone_temp_scale, eq. (4)): the response summed over the root layers is !
!      the mean temperature's on a uniform profile, and higher on a profile spread about  !
!      the same mean, the response being convex below its optimum.                         !
!   3. SAPWOOD (stem_maintenance_respiration, eq. (3)): the sapwood term adds to the       !
!      surface term.                                                                       !
!   4. GRASS (stem_maintenance_respiration): a cohort that is not woody respires no stem.  !
!   5. ABOVEGROUND FRACTION (stem_maintenance_respiration): the stem area is the plant's   !
!      above ground over its aboveground fraction, so halving the fraction doubles it.     !
!   6. SCHEME TOPOLOGY (decomposition_transfer_matrix, soil_carbon.md §3): the default     !
!      scheme sends only the structural pools to slow and holds microbial and passive      !
!      inert; scheme 5 routes by soil texture, at the hand values for sand 0.35, clay 0.25. !
!   7. Rh COMPLEMENT (heterotrophic_respiration, decomposition_transfer_matrix, eq. (2)):   !
!      Rh = sum_j er_j xi_j K_j X_j, and every column of A sums to -er_j with a_jj = -1.    !
!   8. SCALAR PLACEMENT (heterotrophic_respiration, §4): the network's loss with each      !
!      pool's scalar on its donor, A (xi K X), is Rh; on the receiver, xi (A K X), it is   !
!      not.                                                                                !
!   9. TEMPERATURE AND MOISTURE (decomposition_env_scalar, eqs. (3)-(4)): scheme 5's Q10   !
!      doubles xi per 10 K at Q10 = 2; the default scheme caps its temperature factor at 1 !
!      (47 C); xi peaks at resp_opt_water; with no lignin a structural pool's xi is its    !
!      layer's.                                                                            !
!  10. FAST/SLOW SEAM (decomposition_env_scalar, heterotrophic_respiration, §6-§7): Rh is  !
!      linear in xi, so Rh at the day's integral of xi is the integral of the sub-step Rh; !
!      under a diurnal cycle xi at the daily-mean temperature falls short of the day's     !
!      mean xi (Jensen: the Q10 response is convex), so the day must integrate xi.         !
! Growth respiration is the daily allocator's (test_plant_carbon_allocation), and the      !
! pools' daily step ed_soil_carbon's (test_soil_carbon).                                   !
!==========================================================================================!
program test_respiration
   use ed_test_assert,         only : check_close, check_true, test_report
   use ed_kinds,               only : wp, ik
   use ed_constants,           only : t_ref_photo, pi
   use ed_patch_state_vars,    only : soil_carbon_state_t
   use ed_biogeochem_params,   only : soil_carbon_params_t, DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5, N_SOIL_POOL,  &
                                      IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW,  &
                                      IP_PASSIVE
   use ed_respiration,         only : stem_maintenance_respiration, fine_root_maintenance_respiration,              &
                                      root_zone_temp_scale, decomposition_env_scalar, decomposition_transfer_matrix, &
                                      heterotrophic_respiration
   implicit none
   !----- Leaf Rd's peaked response, as the reference configuration sets it ([leaf_physiology]). -!
   real(wp), parameter :: EA_RD = 46390.0_wp, HD_RD = 200000.0_wp, DS_RD = 490.0_wp
   real(wp), parameter :: AGF  = 0.7_wp    !< [-] the PFT aboveground fraction, passed per cohort
   real(wp), parameter :: SRF  = 0.06_wp   !< [umol/m2 stem/s @25C] the PFT stem baseline
   real(wp), parameter :: RRF  = 0.30_wp   !< [umol/kgC root/s @25C] the PFT root baseline
   real(wp), parameter :: RSAP = 40.0_wp   !< [umol/m3 sapwood/s @25C] the PFT sapwood rate
   real(wp), parameter :: VSAP = 0.05_wp   !< [m3/plant] a cohort's sapwood volume
   logical,  parameter :: WOODY = .true.

   call test_reference_temperature()
   call test_root_layers()
   call test_sapwood_adds()
   call test_grass()
   call test_aboveground_frac()
   call test_scheme_topology()
   call test_rh_complement()
   call test_scalar_placement()
   call test_temperature_moisture()
   call test_fast_slow_seam()

   call test_report('test_respiration')

contains

   !----- 1. At 25 C the peaked response is 1: the roots respire their 25 C rate times their      !
   !      carbon, and the sapwood, with no surface rate, its rate times its volume. --------------!
   subroutine test_reference_temperature()
      real(wp) :: r_root, r_sap
      call fine_root_maintenance_respiration(root_zone_temp_scale([t_ref_photo], [1.0_wp], 1_ik, EA_RD, HD_RD, DS_RD), &
                                             2.0_wp, RRF, r_root)
      call check_close('root @25C = factor25 * broot', r_root, RRF * 2.0_wp)
      call stem_maintenance_respiration(t_ref_photo, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,   &
                                        0.0_wp, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_sap)
      call check_close('sapwood @25C = rate x volume', r_sap, RSAP * VSAP)
   end subroutine test_reference_temperature

   !----- 2. The response is summed OVER LAYERS (eq. (4)), not taken at a mean temperature. A      !
   !      UNIFORM profile must give exactly the mean-temperature answer (the layer sum is inert     !
   !      where there is no gradient), and a SPLIT profile with the same root-weighted mean (285 K, !
   !      over a 40 K range) must NOT: below the optimum the response is convex, so the layered sum !
   !      is higher, which a mean-temperature formulation cannot express. --------------------------!
   subroutine test_root_layers()
      real(wp), parameter :: FRAC(4) = [0.4_wp, 0.3_wp, 0.2_wp, 0.1_wp]
      real(wp), parameter :: T_UNIFORM(4) = [285.0_wp, 285.0_wp, 285.0_wp, 285.0_wp]
      real(wp), parameter :: T_SPLIT(4) = [275.0_wp, 295.0_wp, 275.0_wp, 315.0_wp]
      real(wp) :: f_mean, f_layer
      f_mean  = root_zone_temp_scale([285.0_wp], [1.0_wp], 1_ik, EA_RD, HD_RD, DS_RD)
      f_layer = root_zone_temp_scale(T_UNIFORM, FRAC, 4_ik, EA_RD, HD_RD, DS_RD)
      call check_close('uniform profile: layered == mean-temperature', f_layer, f_mean)
      f_layer = root_zone_temp_scale(T_SPLIT, FRAC, 4_ik, EA_RD, HD_RD, DS_RD)
      call check_true('a temperature GRADIENT changes the response at fixed mean',              &
                      abs(f_layer - f_mean) > 1.0e-3_wp * f_mean, f_layer / f_mean)
      call check_true('convex below the optimum => layered exceeds mean-temperature',           &
                      f_layer > f_mean, f_layer / f_mean)
   end subroutine test_root_layers

   !----- 3. The surface and sapwood terms share the temperature response and add. ----------------!
   subroutine test_sapwood_adds()
      real(wp) :: r_sap, r_surf, r_both
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        0.0_wp, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_sap)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        SRF, 0.0_wp, 0.0_wp, EA_RD, HD_RD, DS_RD, r_surf)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        SRF, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_both)
      call check_close('surface and sapwood terms add', r_both, r_surf + r_sap)
   end subroutine test_sapwood_adds

   !----- 4. is_woody is a per-PFT trait: a grass cohort respires no stem, surface or sapwood. ----!
   subroutine test_grass()
      real(wp) :: r_grass
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, .false.,    &
                                        SRF, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_grass)
      call check_true('a grass cohort respires no stem', r_grass == 0.0_wp, r_grass)
   end subroutine test_grass

   !----- 5. The aboveground fraction is a per-cohort trait and must scale the stem area. A        !
   !      run-uniform value would make two PFTs with different allocation respire identically      !
   !      while demography used their real values. The kernel divides the stem area by it, so      !
   !      halving the fraction must double the respiration exactly. -------------------------------!
   subroutine test_aboveground_frac()
      real(wp) :: r_high, r_low
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, 0.7_wp,  WOODY, SRF, 0.0_wp,      &
                                        0.0_wp, EA_RD, HD_RD, DS_RD, r_high)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, 0.35_wp, WOODY, SRF, 0.0_wp,      &
                                        0.0_wp, EA_RD, HD_RD, DS_RD, r_low)
      call check_close('aboveground_frac scales stem respiration inversely', r_low, 2.0_wp * r_high)
   end subroutine test_aboveground_frac

   !----- 6. The routing of decayed carbon, against hand values. The default scheme's only        !
   !      transfers are the structural pools' 0.7 to slow (er_struct 0.3), the rest is respired,   !
   !      and microbial and passive are inert (K = 0). Scheme 5 routes by texture: microbial's     !
   !      respired share rises with sand, and microbial and slow send clay's share to passive. ----!
   subroutine test_scheme_topology()
      type(soil_carbon_state_t)  :: pools
      type(soil_carbon_params_t) :: opts
      real(wp) :: a_mat(N_SOIL_POOL,N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL), offsum
      real(wp) :: er_micr, ex_mp, ex_sp
      integer(ik) :: i, j
      opts = soil_carbon_params_t()                            ! the default scheme
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      call check_close('scheme0 K_micr = 0', k_diag(IP_MICR),    0.0_wp, 0.0_wp)
      call check_close('scheme0 K_pass = 0', k_diag(IP_PASSIVE), 0.0_wp, 0.0_wp)
      call check_close('scheme0 a(slow,struct_grnd)=0.7', a_mat(IP_SLOW,IP_STRUCT_GRND), 0.7_wp, 1.0e-14_wp)
      call check_close('scheme0 a(slow,struct_soil)=0.7', a_mat(IP_SLOW,IP_STRUCT_SOIL), 0.7_wp, 1.0e-14_wp)
      call check_close('scheme0 er(struct_grnd)=0.3', er(IP_STRUCT_GRND), 0.3_wp, 1.0e-14_wp)
      call check_close('scheme0 er(fast_soil)=1.0',   er(IP_FAST_SOIL),   1.0_wp, 1.0e-14_wp)
      !----- Everything but those two transfers and the -1 diagonal is 0. ----------------------!
      offsum = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         do i = 1_ik, N_SOIL_POOL
            if (i /= j .and. .not. (i == IP_SLOW .and. (j == IP_STRUCT_GRND .or. j == IP_STRUCT_SOIL))) then
               offsum = offsum + abs(a_mat(i,j))
            end if
         end do
      end do
      call check_true('scheme0 no other transfers', offsum < 1.0e-14_wp, offsum)

      !----- Scheme 5 at sand 0.35 and clay 0.25, no lignin. ------------------------------------!
      opts = century5_params()
      opts%xsand = 0.35_wp
      opts%xclay = 0.25_wp
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      er_micr = 0.60_wp + 0.17_wp * 0.35_wp
      ex_mp   = min(1.0_wp - er_micr, 0.003_wp + 0.032_wp * 0.25_wp)
      ex_sp   = min(1.0_wp - 0.55_wp, 0.003_wp + 0.009_wp * 0.25_wp)
      call check_close('scheme5 a(micr,fast_grnd)=0.45', a_mat(IP_MICR,IP_FAST_GRND), 0.45_wp, 1.0e-14_wp)
      call check_close('scheme5 a(micr,struct_grnd)=0.50', a_mat(IP_MICR,IP_STRUCT_GRND), 0.50_wp, 1.0e-14_wp)
      call check_close('scheme5 a(pass,micr)=ex_mp', a_mat(IP_PASSIVE,IP_MICR), ex_mp, 1.0e-14_wp)
      call check_close('scheme5 a(slow,micr)=(1-er_micr)-ex_mp', a_mat(IP_SLOW,IP_MICR), (1.0_wp - er_micr) - ex_mp, &
                       1.0e-14_wp)
      call check_close('scheme5 a(pass,slow)=ex_sp', a_mat(IP_PASSIVE,IP_SLOW), ex_sp, 1.0e-14_wp)
      call check_close('scheme5 a(slow,pass)=0.45', a_mat(IP_SLOW,IP_PASSIVE), 0.45_wp, 1.0e-14_wp)
      call check_close('scheme5 er(micr)=er_micr', er(IP_MICR), er_micr, 1.0e-13_wp)
   end subroutine test_scheme_topology

   !----- 7. Rh is the carbon the pools lose that A passes to no other pool: the sum over pools of  !
   !      er_j xi_j K_j X_j, at er from the routine that builds A. And A's columns say the same     !
   !      thing: each sums to -er_j, its diagonal -1 removing the donor's whole loss. --------------!
   subroutine test_rh_complement()
      real(wp), parameter :: XI(N_SOIL_POOL) = [0.585_wp, 0.471_wp, 0.656_wp, 0.542_wp, 0.427_wp, 0.612_wp, 0.498_wp]
      type(soil_carbon_state_t)  :: pools
      type(soil_carbon_params_t) :: opts
      real(wp) :: a_mat(N_SOIL_POOL,N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: rh_fun, rh_hand, worst
      integer(ik) :: j
      opts  = century5_params()
      pools = century5_pools()
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      rh_fun  = heterotrophic_respiration(a_mat, k_diag, XI, pools)
      rh_hand = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         rh_hand = rh_hand + er(j) * XI(j) * k_diag(j) * pools%carbon(j)
      end do
      call check_close('Rh matrix == sum er*xi*K*X', rh_fun, rh_hand, 1.0e-14_wp * max(abs(rh_hand), 1.0_wp))
      worst = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         worst = max(worst, abs(sum(a_mat(:, j)) + er(j)))
         worst = max(worst, abs(a_mat(j,j) + 1.0_wp))
      end do
      call check_true('column identity sum_i a_ij = -er_j, a_jj=-1', worst < 1.0e-14_wp, worst)
   end subroutine test_rh_complement

   !----- 8. Each pool's scalar belongs to the DONOR: dX_i = sum_j a_ij (xi_j K_j X_j). With a      !
   !      surface scalar unlike the soil's, the network then loses exactly Rh; putting the scalar  !
   !      on the receiver, xi_i (sum_j a_ij K_j X_j), does not conserve carbon to it. --------------!
   subroutine test_scalar_placement()
      real(wp), parameter :: XI(N_SOIL_POOL) = [0.9_wp, 0.3_wp, 0.9_wp, 0.3_wp, 0.3_wp, 0.3_wp, 0.3_wp]
      type(soil_carbon_state_t)  :: pools
      type(soil_carbon_params_t) :: opts
      real(wp) :: a_mat(N_SOIL_POOL,N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL), x(N_SOIL_POOL)
      real(wp) :: d_donor(N_SOIL_POOL), d_receiver(N_SOIL_POOL), rh
      integer(ik) :: i, j
      opts = soil_carbon_params_t()                            ! the default scheme
      pools%carbon(IP_FAST_GRND)   = 1.0_wp
      pools%carbon(IP_FAST_SOIL)   = 1.5_wp
      pools%carbon(IP_STRUCT_GRND) = 2.0_wp
      pools%carbon(IP_STRUCT_SOIL) = 3.0_wp
      pools%carbon(IP_SLOW)        = 5.0_wp
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      x = pools%carbon
      d_donor    = 0.0_wp
      d_receiver = 0.0_wp
      do i = 1_ik, N_SOIL_POOL
         do j = 1_ik, N_SOIL_POOL
            d_donor(i)    = d_donor(i)    + a_mat(i,j) * (XI(j) * k_diag(j) * x(j))
            d_receiver(i) = d_receiver(i) + a_mat(i,j) * (XI(i) * k_diag(j) * x(j))
         end do
      end do
      rh = heterotrophic_respiration(a_mat, k_diag, XI, pools)
      call check_close('donor scalar: the network loses Rh', -sum(d_donor), rh, 1.0e-14_wp * max(rh, 1.0_wp))
      call check_true('receiver scalar: it does not', abs(-sum(d_receiver) - rh) > 1.0e-6_wp, -sum(d_receiver) - rh)
   end subroutine test_scalar_placement

   !----- 9. The environmental scalar's laws. Scheme 5's Q10 is uncapped, so at Q10 = 2 xi doubles  !
   !      over 10 K. The default scheme's exponential is capped at 1 from 45 C, so at 47 C and the   !
   !      moisture optimum xi is 1. The moisture factor peaks at resp_opt_water. And with no lignin  !
   !      there is no brake, so a structural pool's xi is the fast pool's in its layer. -------------!
   subroutine test_temperature_moisture()
      type(soil_carbon_state_t)  :: pools
      type(soil_carbon_params_t) :: opts
      real(wp) :: xi_a(N_SOIL_POOL), xi_b(N_SOIL_POOL), xi_opt(N_SOIL_POOL), xi_dry(N_SOIL_POOL), theta_opt
      opts = soil_carbon_params_t()
      opts%decomp_scheme = DECOMP_SCHEME_CENTURY5
      opts%rh_q10   = 2.0_wp
      opts%rh_t_ref = 288.15_wp
      call decomposition_env_scalar(298.15_wp, 298.15_wp, 0.5_wp, 0.1_wp, 0.6_wp, pools, opts, xi_a)
      call decomposition_env_scalar(308.15_wp, 308.15_wp, 0.5_wp, 0.1_wp, 0.6_wp, pools, opts, xi_b)
      call check_close('Q10=2 doubling per 10K', xi_b(IP_FAST_SOIL) / xi_a(IP_FAST_SOIL), 2.0_wp, 1.0e-12_wp)
      !----- The default scheme, at the water content whose relative saturation is the optimum. --!
      opts = soil_carbon_params_t()
      theta_opt = 0.1_wp + opts%resp_opt_water * (0.6_wp - 0.1_wp)
      call decomposition_env_scalar(320.15_wp, 320.15_wp, theta_opt, 0.1_wp, 0.6_wp, pools, opts, xi_opt)
      call check_close('ED2 cap: f_temp<=1 at 47C, rel=opt => xi=1', xi_opt(IP_FAST_SOIL), 1.0_wp, 1.0e-12_wp)
      call decomposition_env_scalar(298.15_wp, 298.15_wp, theta_opt, 0.1_wp, 0.6_wp, pools, opts, xi_opt)
      call decomposition_env_scalar(298.15_wp, 298.15_wp, 0.2_wp, 0.1_wp, 0.6_wp, pools, opts, xi_dry)
      call check_true('moisture peak at opt (dry side lower)', xi_dry(IP_FAST_SOIL) < xi_opt(IP_FAST_SOIL),     &
                      xi_opt(IP_FAST_SOIL) - xi_dry(IP_FAST_SOIL))
      call check_close('no lignin: struct xi == fast xi', xi_opt(IP_STRUCT_SOIL), xi_opt(IP_FAST_SOIL), 1.0e-14_wp)
   end subroutine test_temperature_moisture

   !----- 10. The fast/slow seam (soil_carbon.md §6-§7). The fast loop respires each sub-step at its !
   !       own xi; the daily step debits the pools with the day's integral of xi. The two agree only !
   !       because Rh is linear in xi, so Rh at that integral must be the integral of the sub-step   !
   !       Rh. Over a diurnal cycle of 16 K the day's mean xi exceeds xi at the daily-mean           !
   !       temperature, the Q10 response being convex (Jensen), so a step that took xi from the      !
   !       daily mean would debit the pools by less than the air received. ---------------------------!
   subroutine test_fast_slow_seam()
      integer(ik), parameter :: NSUB   = 48_ik                     ! sub-steps in the day
      real(wp),    parameter :: DT_SUB = 1.0_wp / real(NSUB, wp)   ! [day]
      real(wp),    parameter :: THETA  = 0.45_wp                   ! [m3/m3] the day's water content
      type(soil_carbon_state_t)  :: pools
      type(soil_carbon_params_t) :: opts
      real(wp) :: a_mat(N_SOIL_POOL,N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: xi(N_SOIL_POOL), xi_int(N_SOIL_POOL), xi_at_mean(N_SOIL_POOL)
      real(wp) :: t_surf, t_mean, rh_fast, rh_at_mean
      integer(ik) :: is
      opts  = century5_params()
      pools = century5_pools()
      pools%lignin = [0.6_wp, 0.9_wp]
      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      xi_int  = 0.0_wp
      rh_fast = 0.0_wp
      t_mean  = 0.0_wp
      do is = 1_ik, NSUB
         t_surf = 295.15_wp + 8.0_wp * sin(2.0_wp * pi * (real(is, wp) - 0.5_wp) / real(NSUB, wp))
         call decomposition_env_scalar(t_surf, t_surf - 2.0_wp, THETA, 0.1_wp, 0.6_wp, pools, opts, xi)
         xi_int  = xi_int  + xi * DT_SUB
         rh_fast = rh_fast + heterotrophic_respiration(a_mat, k_diag, xi, pools) * DT_SUB
         t_mean  = t_mean  + t_surf * DT_SUB
      end do
      call check_close('Rh at the day''s integral of xi == the integral of the sub-step Rh',                     &
                       heterotrophic_respiration(a_mat, k_diag, xi_int, pools), rh_fast, 1.0e-13_wp * rh_fast)
      !----- Over one day the integral of xi [day] is the day's mean xi. -----------------------!
      call decomposition_env_scalar(t_mean, t_mean - 2.0_wp, THETA, 0.1_wp, 0.6_wp, pools, opts, xi_at_mean)
      call check_true('every pool''s xi at the daily-mean temperature is below its mean xi (Jensen)',            &
                      all(xi_at_mean < xi_int), minval(xi_int - xi_at_mean))
      rh_at_mean = heterotrophic_respiration(a_mat, k_diag, xi_at_mean, pools)
      call check_true('so Rh from the daily mean falls short of the fast loop''s', rh_at_mean < rh_fast,          &
                      rh_fast - rh_at_mean)
   end subroutine test_fast_slow_seam

   !----- The CENTURY scheme with five active pools, at the rates [1/yr] and respired fractions    !
   !      the soil sections share. ---------------------------------------------------------------!
   function century5_params() result(opts)
      type(soil_carbon_params_t) :: opts
      opts = soil_carbon_params_t()
      opts%decomp_scheme    = DECOMP_SCHEME_CENTURY5
      opts%k_fast           = 12.0_wp
      opts%k_struct         = 1.5_wp
      opts%k_micr           = 6.0_wp
      opts%k_slow           = 0.15_wp
      opts%k_pass           = 0.012_wp
      opts%er_fast          = 0.55_wp
      opts%er_struct_nonlig = 0.50_wp
      opts%er_struct_lig    = 0.20_wp
      opts%er_slow          = 0.55_wp
      opts%er_pass          = 0.55_wp
      opts%er_micr_int      = 0.60_wp
      opts%er_micr_slp      = 0.17_wp
   end function century5_params

   !----- Carbon [kgC/m2] in all seven pools, no lignin. --------------------------------------!
   function century5_pools() result(pools)
      type(soil_carbon_state_t) :: pools
      pools = soil_carbon_state_t()
      pools%carbon(IP_FAST_GRND)   = 0.4_wp
      pools%carbon(IP_FAST_SOIL)   = 0.6_wp
      pools%carbon(IP_STRUCT_GRND) = 2.0_wp
      pools%carbon(IP_STRUCT_SOIL) = 3.0_wp
      pools%carbon(IP_MICR)        = 0.5_wp
      pools%carbon(IP_SLOW)        = 8.0_wp
      pools%carbon(IP_PASSIVE)     = 20.0_wp
   end function century5_pools

end program test_respiration
