! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biogeochem_params -- the [soil_carbon] block (docs/science/soil_carbon.md): the     !
! decomposition selectors and their codes, the decay rates, the respired fractions and the !
! environmental-scalar parameters; and the pool layout they index, N_SOIL_POOL and IP_*.   !
! A config leaf (meds_kinds only): the soil-carbon kernels and the patch state import it  !
! without meds_config.                                                                     !
!==========================================================================================!
module meds_biogeochem_params
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: N_SOIL_POOL
   public :: IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   public :: SOIL_POOL_NAME
   public :: DECOMP_STEP_EULER, DECOMP_STEP_EXPM
   public :: DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   public :: soil_carbon_params_t

   !==========================================================================================!
   !  The soil-carbon matrix's pools: the count and the index of each (the matrix state vector X). !
   !  Ordered litter -> SOM -> passive so the default-scheme ACTIVE block of A is lower-triangular  !
   !  (back-substitution / GE solvable). The 7-pool ceiling matches ED2's CENTURY set exactly.      !
   !==========================================================================================!
   integer(ik), parameter :: N_SOIL_POOL = 7_ik   !< length of the carbon matrix state vector X

   integer(ik), parameter :: IP_FAST_GRND   = 1_ik   !< metabolic litter, above-ground (flammable)
   integer(ik), parameter :: IP_FAST_SOIL   = 2_ik   !< metabolic litter, below-ground
   integer(ik), parameter :: IP_STRUCT_GRND = 3_ik   !< structural litter + CWD, above-ground
   integer(ik), parameter :: IP_STRUCT_SOIL = 4_ik   !< structural litter + CWD, below-ground
   integer(ik), parameter :: IP_MICR        = 5_ik   !< microbial SOM        (active only in scheme 5)
   integer(ik), parameter :: IP_SLOW        = 6_ik   !< slow / humified SOM
   integer(ik), parameter :: IP_PASSIVE     = 7_ik   !< passive SOM          (active only in scheme 5)
   !----- Each pool's name, by index: the restart file's soilc_<name>_pa and the reports.  !
   !      Index 0 is no pool, for a check that found none. --------------------------------!
   character(len=11), parameter :: SOIL_POOL_NAME(0:N_SOIL_POOL) = [character(len=11) ::          &
      '(none)', 'fast_grnd', 'fast_soil', 'struct_grnd', 'struct_soil', 'microbial', 'slow', 'passive']

   !----- Daily-step solver (docs/science/soil_carbon.md §6). -----------------------------!
   integer(ik), parameter :: DECOMP_STEP_EULER = 1_ik  !< forward Euler (the coupled daily step)
   integer(ik), parameter :: DECOMP_STEP_EXPM  = 2_ik  !< exact matrix exponential (large accelerated or spin-up steps)

   !----- Decomposition scheme: the pool topology and the temperature form, not the rates  !
   !      (docs/science/soil_carbon.md §2). -----------------------------------------------!
   integer(ik), parameter :: DECOMP_SCHEME_ED2      = 0_ik  !< 3-active, capped-exponential temperature
   integer(ik), parameter :: DECOMP_SCHEME_CENTURY5 = 5_ik  !< 5-active, Collatz Q10 temperature

   !=======================================================================================!
   !  The [soil_carbon] block; every key is optional, at the default below. Decay rates are!
   !  in [1/yr], and the kernel converts them to [1/day]. One set of environmental-scalar  !
   !  parameters serves the sub-daily Rh and the daily step, so the two reconcile          !
   !  (docs/science/soil_carbon.md §4, §7).                                                !
   !=======================================================================================!
   type :: soil_carbon_params_t
      logical     :: soil_carbon_on = .true.              !< .false.: no soil carbon at all
      integer(ik) :: decomp_scheme = DECOMP_SCHEME_ED2    !< DECOMP_SCHEME_ED2 (3-active) or DECOMP_SCHEME_CENTURY5 (5-active)
      integer(ik) :: step_solver   = DECOMP_STEP_EULER    !< {EULER, EXPM}
      !----- Baseline decay rates [1/yr]. -------------------------------------------------!
      real(wp)    :: k_fast   = 11.0_wp    !< metabolic litter
      real(wp)    :: k_struct = 4.5_wp     !< structural litter and CWD
      real(wp)    :: k_micr   = 0.0_wp     !< microbial SOM; 0: inert
      real(wp)    :: k_slow   = 0.2_wp     !< slow SOM
      real(wp)    :: k_pass   = 0.0_wp     !< passive SOM; 0: inert
      !----- Respired fractions [-]: the share of a pool's decay that leaves as CO2. ------!
      real(wp)    :: er_fast         = 1.0_wp   !< metabolic litter
      real(wp)    :: er_struct_nonlig= 0.3_wp   !< structural litter, its non-lignin part
      real(wp)    :: er_struct_lig   = 0.3_wp   !< structural litter, its lignin part
      real(wp)    :: er_slow         = 1.0_wp   !< slow SOM
      real(wp)    :: er_pass         = 1.0_wp   !< passive SOM
      real(wp)    :: er_micr_int     = 0.0_wp   !< microbial SOM (century5): er_micr = int + slp*xsand
      real(wp)    :: er_micr_slp     = 0.0_wp   !< its sand slope
      real(wp)    :: e_lignin        = 3.0_wp   !< lignin brake exponent  L = exp(-e_lignin*f_lignin)
      !----- Texture (century5): sand sets er_micr, clay the transfers to passive. --------!
      real(wp)    :: xsand = 0.35_wp, xclay = 0.25_wp
      real(wp)    :: fx_micr_pass_int = 0.003_wp, fx_micr_pass_slp = 0.032_wp   !< ex_micr->pass = int + slp*xclay
      real(wp)    :: fx_slow_pass_int = 0.003_wp, fx_slow_pass_slp = 0.009_wp   !< ex_slow->pass = int + slp*xclay
      !----- Environmental scalar xi (docs/science/soil_carbon.md eqs. (3)-(4)). ----------!
      real(wp)    :: rh_q10               = 1.5_wp        !< [-] Q10 (century5)
      real(wp)    :: rh_t_ref             = 288.15_wp     !< [K] 15 C Q10 reference (century5)
      real(wp)    :: resp_temp_increase   = 0.0757_wp     !< [1/K] capped-exponential slope (ed2)
      real(wp)    :: resp_temp_ref        = 318.15_wp     !< [K] 45 C saturation (ed2)
      real(wp)    :: resp_opt_water       = 0.8938_wp     !< [-] optimum relative moisture
      real(wp)    :: resp_water_below_opt = 5.0786_wp     !< [-] moisture shape below the optimum
      real(wp)    :: resp_water_above_opt = 4.5139_wp     !< [-] and above it
      !----- Cold start: zero, or the steady state (docs/science/soil_carbon.md §8). ------!
      logical     :: spinup_steady      = .false.         !< .true.: start from the steady state
      real(wp)    :: spinup_xi          = 1.0_wp          !< [-] its environmental scalar, every pool; also a
                                                          !<     run's without the fast loop (soil_carbon.md §7)
      real(wp)    :: spinup_labile_grnd = 0.0_wp          !< [kgC/m2/day] its litter input, per pool
      real(wp)    :: spinup_labile_soil = 0.0_wp
      real(wp)    :: spinup_struct_grnd = 0.0_wp
      real(wp)    :: spinup_struct_soil = 0.0_wp
   end type soil_carbon_params_t

end module meds_biogeochem_params
