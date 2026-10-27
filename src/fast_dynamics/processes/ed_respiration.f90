! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_respiration -- the patch's respiration besides the leaf's: the plants' maintenance
! respiration of stems and fine roots, per plant, at rates referenced to 25 degC
! (docs/science/plant_respiration.md §2-§3), and the soil's heterotrophic respiration
! from the decomposing carbon pools (docs/science/soil_carbon.md §3-§4):
!   * stem_maintenance_respiration      -- stem surface and sapwood terms, eq. (3);
!   * fine_root_maintenance_respiration -- per unit fine-root carbon, eq. (4), with the
!                                          layer-summed scale root_zone_temp_scale;
!   * decomposition_env_scalar          -- each pool's environmental scalar xi;
!   * decomposition_transfer_matrix     -- the pools' decay rates K, the routing A of the
!                                          decayed carbon and the respired fractions er;
!   * heterotrophic_respiration         -- Rh = sum_j er_j xi_j K_j X_j.
! Leaf Rd is the leaf solver's (ed_leaf_gas_exchange), growth respiration the daily
! allocator's, and the pools' daily step ed_soil_carbon's.
!---~---
module ed_respiration
   use ed_kinds,               only : wp, ik
   use ed_constants,           only : pi, tiny_num, yr_day
   use ed_therm_lib,           only : peaked_arrhenius_scale
   use ed_patch_state_vars,    only : soil_carbon_state_t
   use ed_biogeochem_params,   only : soil_carbon_params_t, DECOMP_SCHEME_CENTURY5, N_SOIL_POOL, IP_FAST_GRND,      &
                                      IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration, root_zone_temp_scale
   public :: decomposition_env_scalar, decomposition_transfer_matrix, heterotrophic_respiration

contains

   !--- Plant (autotrophic) maintenance respiration (plant_respiration.md)

   !---~---
   !   Stem maintenance respiration [umol CO2/plant/s] (docs/science/plant_respiration.md
   ! eq. (3)). Elemental: a scalar call does one cohort, an array call a patch, with the
   ! parameters broadcast. A grass (is_woody = .false.) respires nothing here.
   !---~---
   elemental pure subroutine stem_maintenance_respiration(wood_temp, dbh, height, wai, nplant,          &
                                                         aboveground_frac, is_woody, resp_factor25,     &
                                                         sap_volume, resp_sapwood25, ea, hd, ds, stem_resp)
      real(wp),            intent(in)  :: wood_temp   ! [K]  woody-tissue temperature
      real(wp),            intent(in)  :: dbh         ! [cm] stem diameter at breast height
      real(wp),            intent(in)  :: height      ! [m]  cohort height
      real(wp),            intent(in)  :: wai         ! [m2/m2 ground] wood area index
      real(wp),            intent(in)  :: nplant      ! [plant/m2] stem density
      real(wp),            intent(in)  :: aboveground_frac  ! [-] cohort PFT's aboveground fraction
      logical,             intent(in)  :: is_woody          ! cohort PFT is woody (grass => 0)
      real(wp),            intent(in)  :: resp_factor25     ! [umol CO2/m2 stem/s @25C] cohort PFT's surface rate
      !--- Sapwood term: the living wood respires by volume, coarse roots included (§2).
      real(wp),            intent(in)  :: sap_volume        ! [m3/plant] sapwood volume
      real(wp),            intent(in)  :: resp_sapwood25    ! [umol CO2/m3 sapwood/s @25C] cohort PFT's sapwood rate
      !--- Leaf Rd's temperature parameters, always peaked (plant_respiration.md §4).
      real(wp),            intent(in)  :: ea          ! [J/mol]   activation energy
      real(wp),            intent(in)  :: hd          ! [J/mol]   deactivation energy
      real(wp),            intent(in)  :: ds          ! [J/mol/K] entropy term
      real(wp),            intent(out) :: stem_resp   ! [umol CO2/plant/s]
      real(wp) :: tscale, stem_area

      if (.not. is_woody) then
         stem_resp = 0.0_wp
         return
      end if

      !--- Peaked temperature response (= 1 at 25 degC) on leaf Rd's parameters.
      tscale = peaked_arrhenius_scale(1.0_wp, ea, hd, ds, wood_temp)
      !---~---
      !   Per-plant stem surface [m2/plant]: the cylinder plus the WAI branch term (WAI is per
      ! ground, hence /nplant), over the aboveground fraction to cover the whole woody plant;
      ! both divisors are kept above zero.
      !---~---
      stem_area = ( pi * (dbh * 1.0e-2_wp) * height                                             &
                  + pi * wai / max(nplant, tiny_num) ) / max(aboveground_frac, tiny_num)
      !---~---
      !   The surface and sapwood terms share the temperature response; a PFT with no sapwood
      ! rate respires by its surface alone.
      !---~---
      stem_resp = resp_factor25 * tscale * stem_area + resp_sapwood25 * tscale * sap_volume
   end subroutine stem_maintenance_respiration

   !---~---
   !   Fine-root maintenance respiration [umol CO2/plant/s]
   ! (docs/science/plant_respiration.md eq. (4)). Elemental; the patch-uniform tscale_root
   ! broadcasts over the cohorts.
   !---~---
   elemental pure subroutine fine_root_maintenance_respiration(tscale_root, broot, resp_factor25,      &
                                                              root_resp)
      real(wp),            intent(in)  :: tscale_root ! [-] root-weighted temperature scale (root_zone_temp_scale)
      real(wp),            intent(in)  :: broot       ! [kgC/plant] fine-root biomass
      real(wp),            intent(in)  :: resp_factor25 ! [umol CO2/kgC root/s @25C] cohort PFT's baseline
      real(wp),            intent(out) :: root_resp   ! [umol CO2/plant/s]
      root_resp = resp_factor25 * tscale_root * broot
   end subroutine fine_root_maintenance_respiration

   !---~---
   !   The root-weighted temperature scale sum_k root_frac_k * f(T_k), with f the peaked
   ! response (= 1 at 25 degC): summed over layers, not f at the mean temperature, whose
   ! error flips sign with the season (docs/science/plant_respiration.md §3). It is the
   ! same for every cohort of a patch.
   !---~---
   pure function root_zone_temp_scale(soil_temp, root_frac, n, ea, hd, ds) result(tscale)
      real(wp),            intent(in) :: soil_temp(:)  ! [K] per-layer soil temperature
      real(wp),            intent(in) :: root_frac(:)  ! [-] per-layer root fraction (sums to 1)
      integer(ik),         intent(in) :: n             ! [-] active layers
      real(wp),            intent(in) :: ea, hd, ds    ! [J/mol], [J/mol], [J/mol/K] leaf Rd's response
      real(wp)    :: tscale                            ! [-]
      integer(ik) :: k
      tscale = 0.0_wp
      do k = 1_ik, n
         tscale = tscale + root_frac(k)                                                          &
                * peaked_arrhenius_scale(1.0_wp, ea, hd, ds, soil_temp(k))
      end do
   end function root_zone_temp_scale

   !--- Soil (heterotrophic) respiration (soil_carbon.md)

   !---~---
   !   decomposition_env_scalar -- each pool's environmental scalar xi, the product of a
   ! temperature factor, a moisture-and-oxygen factor and, for the structural pools
   ! (3, 4), the lignin brake exp(-e_lignin*f_lignin) (docs/science/soil_carbon.md
   ! eqs. (3)-(4)). The above-ground pools (1, 3) take the surface temperature and the
   ! rest the soil's (ED2's A_decomp / B_decomp). xi > 0; under the default scheme it is
   ! also <= 1, while scheme 5's Q10 is uncapped.
   !---~---
   pure subroutine decomposition_env_scalar(soil_temp_surf, soil_temp_active, theta_active, theta_dry, &
                                       theta_sat, pools, params, xi)
      real(wp),                   intent(in)  :: soil_temp_surf, soil_temp_active     ! [K] ground surface, soil
      real(wp),                   intent(in)  :: theta_active, theta_dry, theta_sat   ! [m3/m3]
      type(soil_carbon_state_t),  intent(in)  :: pools                               ! [kgC/m2] for f_lignin
      type(soil_carbon_params_t), intent(in)  :: params                              ! the [soil_carbon] parameters
      real(wp),                   intent(out) :: xi(N_SOIL_POOL)                     ! [-] diagonal scalar, > 0
      real(wp) :: f_temp_surf, f_temp_active, f_water, rel, f_lig_g, f_lig_s, brake_g, brake_s

      f_temp_surf   = decomp_temperature_factor(soil_temp_surf,   params)
      f_temp_active = decomp_temperature_factor(soil_temp_active, params)
      rel     = min(1.0_wp, max(0.0_wp, (theta_active - theta_dry) / max(theta_sat - theta_dry, tiny_num)))
      f_water = decomp_water_factor(rel, params)

      f_lig_g = lignin_fraction(pools%lignin(1), pools%carbon(IP_STRUCT_GRND))
      f_lig_s = lignin_fraction(pools%lignin(2), pools%carbon(IP_STRUCT_SOIL))
      brake_g = exp(-params%e_lignin * f_lig_g)
      brake_s = exp(-params%e_lignin * f_lig_s)

      xi(IP_FAST_GRND)   = f_temp_surf   * f_water
      xi(IP_FAST_SOIL)   = f_temp_active * f_water
      xi(IP_STRUCT_GRND) = f_temp_surf   * f_water * brake_g
      xi(IP_STRUCT_SOIL) = f_temp_active * f_water * brake_s
      xi(IP_MICR)        = f_temp_active * f_water
      xi(IP_SLOW)        = f_temp_active * f_water
      xi(IP_PASSIVE)     = f_temp_active * f_water
   end subroutine decomposition_env_scalar

   !--- f_T of eq. (3): scheme 5's Q10, uncapped; otherwise ED2's exponential, capped at 1.
   pure function decomp_temperature_factor(soil_temp, params) result(f_temp)
      real(wp),                   intent(in) :: soil_temp   ! [K]
      type(soil_carbon_params_t), intent(in) :: params
      real(wp) :: f_temp                                    ! [-]
      if (params%decomp_scheme == DECOMP_SCHEME_CENTURY5) then         ! Collatz/K13 Q10 (uncapped, > 1 above ref)
         f_temp = params%rh_q10 ** ((soil_temp - params%rh_t_ref) * 0.1_wp)
      else                                                           ! ED2 scheme-0 capped exp, min(1, .)
         f_temp = min(1.0_wp, exp(params%resp_temp_increase * (soil_temp - params%resp_temp_ref)))
      end if
   end function decomp_temperature_factor

   !---~---
   !   f_theta of eq. (3): one exponential on each side of resp_opt_water, substrate diffusion
   ! on the dry side and anoxia on the wet side.
   !---~---
   pure function decomp_water_factor(rel, params) result(f_water)
      real(wp),                   intent(in) :: rel         ! [-] relative saturation m, in [0, 1]
      type(soil_carbon_params_t), intent(in) :: params
      real(wp) :: f_water                                   ! [-]
      if (rel <= params%resp_opt_water) then
         f_water = exp((rel - params%resp_opt_water) * params%resp_water_below_opt)      ! dry side
      else
         f_water = exp((params%resp_opt_water - rel) * params%resp_water_above_opt)      ! wet side (anoxia)
      end if
   end function decomp_water_factor

   !--- A structural pool's lignin share of its carbon, f_lignin = L/C in [0, 1].
   pure function lignin_fraction(lignin, carbon) result(f_lignin)
      real(wp), intent(in) :: lignin, carbon                          ! [kgC/m2]
      real(wp) :: f_lignin                                            ! [-]
      if (carbon > tiny_num) then
         f_lignin = min(1.0_wp, max(0.0_wp, lignin / carbon))
      else
         f_lignin = 0.0_wp                                            ! C=0 => f_lignin undefined, read as 0
      end if
   end function lignin_fraction

   !---~---
   !   decomposition_transfer_matrix -- the transfer matrix A, the decay diagonal K and the
   ! respired fractions er (docs/science/soil_carbon.md eq. (2)). A(i,j) is the share of
   ! pool j's decayed carbon entering pool i, with a_jj = -1 for every pool, so
   ! er_j = 1 - sum_{i!=j} a_ij. K is in [1/day]; params carries it in [1/yr]. The default
   ! scheme (ED2's 0-4) sends structural carbon to slow and respires all the rest; scheme
   ! 5 (Bolker/Koven) routes by soil texture, with back-transfers.
   !---~---
   pure subroutine decomposition_transfer_matrix(pools, params, a_mat, k_diag, er)
      type(soil_carbon_state_t),  intent(in)  :: pools                             ! [kgC/m2] for the lignin fractions
      type(soil_carbon_params_t), intent(in)  :: params                            ! the [soil_carbon] parameters
      real(wp),                   intent(out) :: a_mat(N_SOIL_POOL, N_SOIL_POOL)   ! [-] a_ii=-1, a_ij transfer frac
      real(wp),                   intent(out) :: k_diag(N_SOIL_POOL)               ! [1/day] baseline decay rates
      real(wp),                   intent(out) :: er(N_SOIL_POOL)                   ! [-] respired fraction 1 - sum_i a_ij
      real(wp)    :: f_lig_g, f_lig_s, er_stg, er_sts, er_micr, ex_mp, ex_sp
      integer(ik) :: j

      !--- The decay diagonal K [1/day] (params in /yr).
      k_diag(IP_FAST_GRND)   = params%k_fast   / yr_day
      k_diag(IP_FAST_SOIL)   = params%k_fast   / yr_day
      k_diag(IP_STRUCT_GRND) = params%k_struct / yr_day
      k_diag(IP_STRUCT_SOIL) = params%k_struct / yr_day
      k_diag(IP_MICR)        = params%k_micr   / yr_day
      k_diag(IP_SLOW)        = params%k_slow   / yr_day
      k_diag(IP_PASSIVE)     = params%k_pass   / yr_day

      !--- The structural pools' respired fractions, lignin-weighted.
      f_lig_g = lignin_fraction(pools%lignin(1), pools%carbon(IP_STRUCT_GRND))
      f_lig_s = lignin_fraction(pools%lignin(2), pools%carbon(IP_STRUCT_SOIL))
      er_stg  = (1.0_wp - f_lig_g) * params%er_struct_nonlig + f_lig_g * params%er_struct_lig
      er_sts  = (1.0_wp - f_lig_s) * params%er_struct_nonlig + f_lig_s * params%er_struct_lig

      !--- a_jj = -1 for every pool (all decayed carbon leaves the donor), off-diagonals 0.
      a_mat = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         a_mat(j, j) = -1.0_wp
      end do

      if (params%decomp_scheme == DECOMP_SCHEME_CENTURY5) then
         !--- Scheme 5, CENTURY with five active pools (Bolker/Koven): routing set by texture.
         er_micr = params%er_micr_int + params%er_micr_slp * params%xsand             ! sand -> microbial respired frac
         ex_mp   = min(1.0_wp - er_micr, params%fx_micr_pass_int + params%fx_micr_pass_slp * params%xclay)
         ex_sp   = min(1.0_wp - params%er_slow, params%fx_slow_pass_int + params%fx_slow_pass_slp * params%xclay)
         !--- fast litter (1,2) -> microbial.
         a_mat(IP_MICR, IP_FAST_GRND) = 1.0_wp - params%er_fast
         a_mat(IP_MICR, IP_FAST_SOIL) = 1.0_wp - params%er_fast
         !--- structural (3,4): non-lignin -> microbial, lignin -> slow.
         a_mat(IP_MICR, IP_STRUCT_GRND) = (1.0_wp - er_stg) * (1.0_wp - f_lig_g)
         a_mat(IP_MICR, IP_STRUCT_SOIL) = (1.0_wp - er_sts) * (1.0_wp - f_lig_s)
         a_mat(IP_SLOW, IP_STRUCT_GRND) = (1.0_wp - er_stg) * f_lig_g
         a_mat(IP_SLOW, IP_STRUCT_SOIL) = (1.0_wp - er_sts) * f_lig_s
         !--- microbial (5): -> passive (clay), remainder -> slow.
         a_mat(IP_PASSIVE, IP_MICR) = ex_mp
         a_mat(IP_SLOW,    IP_MICR) = (1.0_wp - er_micr) - ex_mp
         !--- slow (6): -> passive (clay), remainder -> microbial.
         a_mat(IP_PASSIVE, IP_SLOW) = ex_sp
         a_mat(IP_MICR,    IP_SLOW) = (1.0_wp - params%er_slow) - ex_sp
         !--- passive (7): -> slow.
         a_mat(IP_SLOW, IP_PASSIVE) = 1.0_wp - params%er_pass
      else
         !--- Schemes 0-4 (three active pools): structural -> slow, the rest respired.
         a_mat(IP_SLOW, IP_STRUCT_GRND) = 1.0_wp - er_stg
         a_mat(IP_SLOW, IP_STRUCT_SOIL) = 1.0_wp - er_sts
      end if

      !--- Respired fraction er_j = 1 - sum_{i!=j} a_ij (= -sum_i a_ij, since a_jj = -1).
      do j = 1_ik, N_SOIL_POOL
         er(j) = 1.0_wp - (sum(a_mat(:, j)) + 1.0_wp)      ! sum_i a_ij = a_jj(-1) + sum_{i!=j}; er = 1 - sum_{i!=j}
      end do
   end subroutine decomposition_transfer_matrix

   !---~---
   !   heterotrophic_respiration -- the soil's heterotrophic respiration at the pools'
   ! current state, Rh = sum_j er_j*xi_j*K_j*X_j = -1^T*A*xi*K*X [kgC/m2/day]
   ! (docs/science/soil_carbon.md eq. (2)): the carbon the pools lose that A does not pass
   ! to another pool. Its A, K and xi are those the pools' daily step debits with, so the
   ! flux and the pools' loss agree (§7). kgcday_2_umols converts it to umol/m2/s.
   !---~---
   pure function heterotrophic_respiration(a_mat, k_diag, xi, pools) result(rh)
      real(wp),            intent(in) :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)   ! [-], [1/day]
      real(wp),            intent(in) :: xi(N_SOIL_POOL)                                       ! [-]
      type(soil_carbon_state_t), intent(in) :: pools                                                 ! [kgC/m2]
      real(wp) :: rh, x(N_SOIL_POOL), er_j                                                     ! [kgC/m2/day] Rh
      integer(ik) :: j
      x = pools%carbon
      rh = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         er_j = -sum(a_mat(:, j))                    ! column sum of A is -er_j (a_jj=-1)
         rh   = rh + er_j * xi(j) * k_diag(j) * x(j)
      end do
   end function heterotrophic_respiration

end module ed_respiration
