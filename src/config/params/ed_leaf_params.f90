! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_leaf_params -- the [leaf_physiology] block's record and its selector codes, as a
! config leaf: ed_config reads and validates them, and the leaf kernel
! (ed_leaf_gas_exchange) switches on the codes without importing ed_config.
!---~---
module ed_leaf_params
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: leaf_physiology_params_t
   public :: SM_LEUNING, SM_MEDLYN, SM_KATUL
   public :: COLIM_MIN, COLIM_QUADRATIC

   !--- [leaf_physiology].stomatal_model (docs/science/leaf_gas_exchange.md §3).
   integer(ik), parameter :: SM_LEUNING = 1_ik      ! Leuning (1995) BWB-VPD semi-empirical
   integer(ik), parameter :: SM_MEDLYN  = 2_ik      ! Medlyn et al. (2011) unified optimization (USO)
   integer(ik), parameter :: SM_KATUL   = 3_ik      ! Katul et al. (2010) analytical optimization

   !--- [leaf_physiology].colimitation (docs/science/leaf_gas_exchange.md §2).
   integer(ik), parameter :: COLIM_MIN       = 1_ik ! sharp minimum
   integer(ik), parameter :: COLIM_QUADRATIC = 2_ik ! smoothed co-limitation quadratics

   !---~---
   !   The temp_response_form codes, TEMP_RESPONSE_*, are in ed_therm_lib with the
   ! responses they select.
   !---~---

   !---~---
   !   [leaf_physiology]: the leaf model's selectors and the biochemistry every PFT
   ! shares (docs/science/leaf_gas_exchange.md), one field per key. A field without an
   ! initializer is a required key.
   !---~---
   type :: leaf_physiology_params_t
      integer(ik) :: stomatal_model                       ! SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form                   ! TEMP_RESPONSE_ARRHENIUS | TEMP_RESPONSE_PEAKED
      integer(ik) :: colimitation                         ! COLIM_MIN | COLIM_QUADRATIC
      logical     :: use_boundary_layer                   ! couple through the leaf boundary layer
      real(wp)    :: kc25, ko25, gstar25                  ! [Pa]
      real(wp)    :: ea_kc, ea_ko, ea_gstar               ! [J/mol]
      real(wp)    :: ea_vcmax, ea_jmax, ea_rd             ! [J/mol]
      real(wp)    :: hd_vcmax, hd_jmax, hd_rd             ! [J/mol]
      real(wp)    :: ds_vcmax, ds_jmax, ds_rd             ! [J/mol/K]
      real(wp)    :: o2_mol_frac                          ! [mol/mol]
      real(wp)    :: leaf_absorptance                     ! [-]
      real(wp)    :: phi_psii                             ! [-]
      real(wp)    :: vpd_min = 0.05_wp                    ! [kPa]
      logical     :: thermal_acclimation = .false.        ! Kattge & Knorr (2007), peaked form only
      real(wp)    :: acclim_ds_vcmax_a = 668.39_wp, acclim_ds_vcmax_b = 1.07_wp   ! [J/mol/K], [J/mol/K2]
      real(wp)    :: acclim_ds_jmax_a  = 659.70_wp, acclim_ds_jmax_b  = 0.75_wp   ! [J/mol/K], [J/mol/K2]
      real(wp)    :: acclim_jv_a       = 2.59_wp,   acclim_jv_b       = 0.035_wp  ! [-], [1/K]
      real(wp)    :: acclim_window_days = 30.0_wp      ! [day]
   end type leaf_physiology_params_t
end module ed_leaf_params
