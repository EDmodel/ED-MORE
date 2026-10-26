! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_cohort_state_vars -- the state of every cohort of a site, and every routine that      !
! lists the cohort fields.                                                                 !
!                                                                                          !
! All cohorts of the site live in one Structure-of-Arrays, `cohort_block`: one contiguous  !
! 1-D array per field, so the daily kernels sweep a field with unit stride. The patch      !
! block groups them with a CSR (compressed sparse row) map, the layout of a sparse         !
! matrix's rows: the cohorts are sorted by patch, cohort_offset(ip) is the slot of patch   !
! ip's first cohort and cohort_count(ip) how many it has, so a patch's cohorts are one     !
! slice, offset .. offset+count-1 (rebuild_cohort_map builds the map).                     !
!                                                                                          !
! The type is laid out in sections, and every routine that lists all its fields uses the   !
! same order. Every structural change -- a sort, a fusion, a compaction, a grow -- goes    !
! through those routines, so a field added to each of them moves with its cohort. The      !
! values a field holds before it has a real one are in ed_state_defaults.                  !
!==========================================================================================!
module ed_cohort_state_vars
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ed_kinds,             only : wp, ik
   use ed_constants,         only : tissue_c_leaf, tissue_c_wood_dry, cp_liq
   use ed_therm_lib,         only : tissue_heat_capacity, tissue_dry_heat_capacity, tissue_temp_to_energy,  &
                                    tissue_energy_to_temp, internal_energy_liquid
   use ed_plant_hydr_lib,    only : water_curve_t, tissue_psi_to_water
   use ed_numerics,          only : pair_mean
   use ed_allometry_lib,     only : dbh_to_height, size_to_leaf_carbon, size_to_wood_carbon,  &
                                    carbon_to_structure, dbh_to_wai, sapwood_fraction, storage_carbon_target
   use ed_pft_params,        only : pft_table_t
   use ed_diag_vars,         only : cohort_diag_block, expand_cohort_diag, reorder_cohort_diag,              &
                                    copy_cohort_diag, clear_cohort_diag
   use ed_state_defaults,    only : UNSET_POSITIVE, UNSET_NEGATIVE, DMAX_PSI_LEAF_ACCUM_RESET,            &
                                    PHENO_FLUSH_INIT, PHENO_SHED_INIT
   implicit none
   private

   public :: cohort_block, alloc_cohort_block, free_cohort_block, expand_cohort_block, reorder_cohort_block
   public :: compact_cohort_block, copy_cohort, set_cohort_size, set_cohort_size_from_carbon, init_cohort
   public :: fuse_cohort_fast_state
   public :: cohort_deriv_block, cohort_deriv_alloc
   public :: cohort_tissue_heat_capacity, cohort_tissue_water, cohort_film_water, cohort_water_enthalpy
   public :: cohort_wood_store_carbon, cohort_carbon_supply
   public :: set_cohort_tissue_water, cohort_tissue_temps, set_cohort_tissue_temps, hold_cohort_tissue_temps
   public :: cohort_bad_field, cohort_field_finite, cohort_field_name

   !----- All cohorts of a site. The model integrates nplant and the four carbon pools; the  !
   !      size fields are a function of the pools, written only by set_cohort_size_from_carbon. !
   type :: cohort_block
      !==========================================================================================!
      ! 0. BLOCK RECORD: one per block, not per cohort                                           !
      !==========================================================================================!
      integer(ik) :: n   = 0_ik            !< cohorts in use
      integer(ik) :: cap = 0_ik            !< allocated slots
      type(pft_table_t) :: pft_table       !< the site's PFT traits (cfg%pft), read as pft_table%x(pft(i));
                                           !< its allometry records carry the polygon's chave_e
      !==========================================================================================!
      ! 1. IDENTITY AND DENSITY                                                                  !
      !==========================================================================================!
      integer(ik), allocatable :: pft(:)          !< PFT index
      integer(ik), allocatable :: global_id(:)    !< persistent id, never reused
      integer(ik), allocatable :: owner_patch(:)  !< the patch it belongs to
      real(wp),    allocatable :: nplant(:)       !< [plant/m2] plant density
      !==========================================================================================!
      ! 2. CARBON POOLS                                                                          !
      !==========================================================================================!
      real(wp),    allocatable :: leaf_carbon(:)          !< [kgC/plant]
      real(wp),    allocatable :: fineroot_carbon(:)      !< [kgC/plant]
      real(wp),    allocatable :: wood_carbon(:)          !< [kgC/plant] the size anchor
      real(wp),    allocatable :: nonstructural_carbon(:) !< [kgC/plant] storage
      !==========================================================================================!
      ! 3. SIZE, DERIVED FROM THE CARBON POOLS (never integrated)                                !
      !==========================================================================================!
      real(wp),    allocatable :: dbh(:)              !< [cm]
      real(wp),    allocatable :: height(:)           !< [m]
      real(wp),    allocatable :: basal_area(:)       !< [cm2/plant]
      real(wp),    allocatable :: agb(:)              !< [kgC/plant] aboveground biomass
      real(wp),    allocatable :: leaf_area(:)        !< [m2/plant]  LAI = nplant * leaf_area
      real(wp),    allocatable :: wood_area(:)        !< [m2/plant]  WAI = nplant * wood_area
      real(wp),    allocatable :: sapwood_carbon(:)   !< [kgC/plant]
      real(wp),    allocatable :: sapwood_area(:)     !< [cm2/plant]
      real(wp),    allocatable :: overtopping_lai(:)  !< [m2/m2] LAI of the taller cohorts in its patch
      !==========================================================================================!
      ! 4. WATER                                                                                 !
      !==========================================================================================!
      real(wp),    allocatable :: leaf_water_mass(:)  !< [kg/plant] water inside the leaves
      real(wp),    allocatable :: wood_water_mass(:)  !< [kg/plant] water inside the sapwood and fine roots
      real(wp),    allocatable :: leaf_surf_water(:)  !< [kg/plant] intercepted water on the leaves
      real(wp),    allocatable :: wood_surf_water(:)  !< [kg/plant] intercepted water on the wood
      real(wp),    allocatable :: dmax_psi_leaf(:)    !< [MPa] yesterday's daily-max (predawn) leaf potential
      !==========================================================================================!
      ! 5. ENERGY                                                                                !
      !==========================================================================================!
      real(wp),    allocatable :: leaf_energy(:)      !< [J/plant] the leaves' internal energy
      real(wp),    allocatable :: wood_energy(:)      !< [J/plant] the wood's internal energy
      real(wp),    allocatable :: leaf_temp(:)        !< [K] derived from leaf_energy
      real(wp),    allocatable :: wood_temp(:)        !< [K] derived from wood_energy
      !==========================================================================================!
      ! 6. PLASTIC LEAF TRAITS (acclimate to light)                                              !
      !==========================================================================================!
      real(wp),    allocatable :: sla(:)                !< [m2/kgC]    specific leaf area
      real(wp),    allocatable :: vcmax25(:)            !< [umol/m2/s] at 25 degC
      real(wp),    allocatable :: rd25(:)               !< [umol/m2/s] leaf dark respiration at 25 degC
      real(wp),    allocatable :: leaf_turnover_rate(:) !< [1/yr]
      !==========================================================================================!
      ! 7. LEAF PHENOLOGY                                                                        !
      !==========================================================================================!
      real(wp),    allocatable :: leaf_flush_tendency(:)  !< [-]       smoothed flush signal, 0-1
      real(wp),    allocatable :: leaf_shed_tendency(:)   !< [-]       smoothed shed signal, 0-1
      real(wp),    allocatable :: growing_degree_days(:)  !< [K day]   warmth since midwinter
      real(wp),    allocatable :: cold_degree_days(:)     !< [K day]   cold since midsummer
      real(wp),    allocatable :: dry_psi_sum(:)          !< [MPa day] predawn potential below the turgor-loss point
      real(wp),    allocatable :: wet_psi_sum(:)          !< [MPa day] predawn potential above it
      real(wp),    allocatable :: light_hours_mean(:)     !< [h/day]   running mean of the hours of light
      !==========================================================================================!
      ! 8. GROWTH HISTORY (the mortality law reads growth_avg)                                   !
      !==========================================================================================!
      real(wp),    allocatable :: growth_avg(:)       !< [cm/yr] moving-average dbh growth, read once growth_count > 0
      real(wp),    allocatable :: growth_accum(:)     !< [cm/yr] sum of the samples in the window
      integer(ik), allocatable :: growth_count(:)     !< samples in the window
      real(wp),    allocatable :: growth_hist(:,:)    !< [cm/yr] (window, cohort) ring buffer of samples
      !==========================================================================================!
      ! 9. FAST-TO-SLOW ACCUMULATORS (summed by the fast loop, then read and reset)              !
      !==========================================================================================!
      real(wp),    allocatable :: gpp_accum(:)            !< [kgC/plant] gross GPP
      real(wp),    allocatable :: leaf_resp_accum(:)      !< [kgC/plant] leaf dark respiration
      real(wp),    allocatable :: stem_resp_accum(:)      !< [kgC/plant] stem maintenance respiration
      real(wp),    allocatable :: root_resp_accum(:)      !< [kgC/plant] fine-root maintenance respiration
      real(wp),    allocatable :: light_hours_accum(:)    !< [h] time the PAR at its top exceeded par_min
      real(wp),    allocatable :: dmax_psi_leaf_accum(:)  !< [MPa] today's running max leaf potential
      !==========================================================================================!
      ! 10. OUTPUT ACCUMULATORS (ed_diag_vars; never restarted)                                  !
      !==========================================================================================!
      type(cohort_diag_block) :: diag    !< fast-loop diagnostics
      type(cohort_diag_block) :: sdiag   !< slow-step diagnostics: growth, mortality, NPP allocation
   end type cohort_block

   !----- The tendencies the slow step computes and update_cohort_states applies [unit/yr], for the !
   !      five integrated fields (dln_nplant_dt in log space). Transient -- filled and applied in    !
   !      one step -- so it is not part of cohort_block and is never reordered.                       !
   type :: cohort_deriv_block
      real(wp), allocatable :: d_leaf_carbon_dt(:), d_fineroot_carbon_dt(:)
      real(wp), allocatable :: d_wood_carbon_dt(:), d_nonstructural_carbon_dt(:)
      real(wp), allocatable :: dln_nplant_dt(:)
   end type cohort_deriv_block

   !----- The cohort fields cohort_bad_field checks, in its order: the first N_COHORT_AMOUNTS  !
   !      are amounts, finite and never below zero; the tissue energies are finite. ---------!
   integer(ik), parameter :: N_COHORT_AMOUNTS = 10_ik, N_COHORT_CHECKED = 12_ik
   character(len=20), parameter :: COHORT_CHECKED_NAME(N_COHORT_CHECKED) = [character(len=20) ::   &
      'nplant', 'leaf_carbon', 'fineroot_carbon', 'wood_carbon', 'nonstructural_carbon', 'sla',     &
      'leaf_water_mass', 'wood_water_mass', 'leaf_surf_water', 'wood_surf_water', 'leaf_energy',    &
      'wood_energy']

contains

   !=======================================================================================!
   !  THE COHORT LOCKSTEP: the routines that list every cohort field, in section order     !
   !=======================================================================================!

   !----- Allocate a block for `cap` cohorts with a growth window of `nwin` samples, section by    !
   !      section, every slot at its unset value. ------------------------------------------------!
   subroutine alloc_cohort_block(cohort, cap, nwin)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: cap, nwin
      cohort%cap = cap
      cohort%n = 0_ik
      !----- 1. Identity and density ------------------------------------------------------------!
      allocate(cohort%pft(cap), cohort%global_id(cap), cohort%owner_patch(cap), cohort%nplant(cap))
      cohort%pft = 0_ik
      cohort%global_id = 0_ik
      cohort%owner_patch = 0_ik
      cohort%nplant = 0.0_wp
      !----- 2. Carbon pools --------------------------------------------------------------------!
      allocate(cohort%leaf_carbon(cap), cohort%fineroot_carbon(cap), cohort%wood_carbon(cap),                  &
               cohort%nonstructural_carbon(cap))
      cohort%leaf_carbon = 0.0_wp
      cohort%fineroot_carbon = 0.0_wp
      cohort%wood_carbon = 0.0_wp
      cohort%nonstructural_carbon = 0.0_wp
      !----- 3. Size, derived from the carbon pools ---------------------------------------------!
      allocate(cohort%dbh(cap), cohort%height(cap), cohort%basal_area(cap), cohort%agb(cap),                   &
               cohort%leaf_area(cap), cohort%wood_area(cap), cohort%sapwood_carbon(cap),                       &
               cohort%sapwood_area(cap), cohort%overtopping_lai(cap))
      cohort%dbh = 0.0_wp
      cohort%height = 0.0_wp
      cohort%basal_area = 0.0_wp
      cohort%agb = 0.0_wp
      cohort%leaf_area = 0.0_wp
      cohort%wood_area = 0.0_wp
      cohort%sapwood_carbon = 0.0_wp
      cohort%sapwood_area = 0.0_wp
      cohort%overtopping_lai = 0.0_wp
      !----- 4. Water ---------------------------------------------------------------------------!
      allocate(cohort%leaf_water_mass(cap), cohort%wood_water_mass(cap), cohort%leaf_surf_water(cap),          &
               cohort%wood_surf_water(cap), cohort%dmax_psi_leaf(cap))
      cohort%leaf_water_mass = 0.0_wp
      cohort%wood_water_mass = 0.0_wp
      cohort%leaf_surf_water = 0.0_wp
      cohort%wood_surf_water = 0.0_wp
      cohort%dmax_psi_leaf = UNSET_NEGATIVE
      !----- 5. Energy --------------------------------------------------------------------------!
      allocate(cohort%leaf_energy(cap), cohort%wood_energy(cap), cohort%leaf_temp(cap), cohort%wood_temp(cap))
      cohort%leaf_energy = 0.0_wp
      cohort%wood_energy = 0.0_wp
      cohort%leaf_temp = 0.0_wp
      cohort%wood_temp = 0.0_wp      ! a slot's tissue state is set at birth
      !----- 6. Plastic leaf traits -------------------------------------------------------------!
      allocate(cohort%sla(cap), cohort%vcmax25(cap), cohort%rd25(cap), cohort%leaf_turnover_rate(cap))
      cohort%sla = 0.0_wp
      cohort%vcmax25 = 0.0_wp
      cohort%rd25 = 0.0_wp
      cohort%leaf_turnover_rate = 0.0_wp
      !----- 7. Leaf phenology ------------------------------------------------------------------!
      allocate(cohort%leaf_flush_tendency(cap), cohort%leaf_shed_tendency(cap),                                &
               cohort%growing_degree_days(cap), cohort%cold_degree_days(cap), cohort%dry_psi_sum(cap),         &
               cohort%wet_psi_sum(cap), cohort%light_hours_mean(cap))
      cohort%leaf_flush_tendency = PHENO_FLUSH_INIT
      cohort%leaf_shed_tendency = PHENO_SHED_INIT
      cohort%growing_degree_days = 0.0_wp
      cohort%cold_degree_days = 0.0_wp
      cohort%dry_psi_sum = 0.0_wp
      cohort%wet_psi_sum = 0.0_wp
      cohort%light_hours_mean = UNSET_POSITIVE
      !----- 8. Growth history ------------------------------------------------------------------!
      allocate(cohort%growth_avg(cap), cohort%growth_accum(cap), cohort%growth_count(cap),                     &
               cohort%growth_hist(nwin, cap))
      cohort%growth_avg = UNSET_POSITIVE
      cohort%growth_accum = 0.0_wp
      cohort%growth_count = 0_ik
      cohort%growth_hist = 0.0_wp
      !----- 9. Fast-to-slow accumulators -------------------------------------------------------!
      allocate(cohort%gpp_accum(cap), cohort%leaf_resp_accum(cap), cohort%stem_resp_accum(cap),                &
               cohort%root_resp_accum(cap), cohort%light_hours_accum(cap), cohort%dmax_psi_leaf_accum(cap))
      cohort%gpp_accum = 0.0_wp
      cohort%leaf_resp_accum = 0.0_wp
      cohort%stem_resp_accum = 0.0_wp
      cohort%root_resp_accum = 0.0_wp
      cohort%light_hours_accum = 0.0_wp
      cohort%dmax_psi_leaf_accum = DMAX_PSI_LEAF_ACCUM_RESET
   end subroutine alloc_cohort_block

   !----- Free every per-cohort array. ------------------------------------------------------!
   subroutine free_cohort_block(cohort)
      type(cohort_block), intent(inout) :: cohort
      cohort%n = 0_ik
      cohort%cap = 0_ik
      if (allocated(cohort%nplant)) deallocate(cohort%pft, cohort%global_id, cohort%owner_patch,               &
                                               cohort%nplant, cohort%leaf_carbon, cohort%fineroot_carbon,      &
                                               cohort%wood_carbon, cohort%nonstructural_carbon, cohort%dbh,    &
                                               cohort%height, cohort%basal_area, cohort%agb, cohort%leaf_area, &
                                               cohort%wood_area, cohort%sapwood_carbon, cohort%sapwood_area,   &
                                               cohort%overtopping_lai, cohort%leaf_water_mass,                 &
                                               cohort%wood_water_mass, cohort%leaf_surf_water,                 &
                                               cohort%wood_surf_water, cohort%dmax_psi_leaf,                   &
                                               cohort%leaf_energy, cohort%wood_energy, cohort%leaf_temp,       &
                                               cohort%wood_temp, cohort%sla, cohort%vcmax25, cohort%rd25,      &
                                               cohort%leaf_turnover_rate, cohort%leaf_flush_tendency,          &
                                               cohort%leaf_shed_tendency, cohort%growing_degree_days,          &
                                               cohort%cold_degree_days, cohort%dry_psi_sum,                    &
                                               cohort%wet_psi_sum, cohort%light_hours_mean, cohort%growth_avg, &
                                               cohort%growth_accum, cohort%growth_count, cohort%growth_hist,   &
                                               cohort%gpp_accum, cohort%leaf_resp_accum,                       &
                                               cohort%stem_resp_accum, cohort%root_resp_accum,                 &
                                               cohort%light_hours_accum, cohort%dmax_psi_leaf_accum)
   end subroutine free_cohort_block

   !=======================================================================================!
   !  expand_cohort_block -- ensure the cohort SoA can hold at least `need` cohorts,       !
   !  growing it only if it must. This is the growable-array primitive (like C++           !
   !  vector::reserve) for the flat cohort block, where `n` is the LIVE count and          !
   !  `cap` the ALLOCATED count (slots n+1..cap are dead reserve). Callers about to        !
   !  APPEND cohorts call it first; it NO-OPs when room already exists (need <= cap),      !
   !  so it is safe to call unconditionally; otherwise it grows GEOMETRICALLY (1.5x)       !
   !  for amortized-O(1) appends, copying the live 1:n prefix and move_alloc-ing in.       !
   !  LOCKSTEP SITE: a new per-cohort field must be added to the 1:n copy list AND to      !
   !  move_cohort_block; a field in move_cohort_block alone would be reset by every grow.   !
   !=======================================================================================!
   subroutine expand_cohort_block(cohort, need)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: need
      type(cohort_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= cohort%cap) return
      newcap = max(need, (cohort%cap * 3_ik) / 2_ik + 1_ik)
      !----- Grow the diagnostic block IN PLACE first, so move_cohort_block below can hand its     !
      !      (already correctly sized, prefix-preserving) arrays straight across.  ---------------!
      call expand_cohort_diag(cohort%diag,  newcap)
      call expand_cohort_diag(cohort%sdiag, newcap)
      call alloc_cohort_block(tmp, newcap, int(size(cohort%growth_hist, 1), ik))  ! keep the window size
      m = cohort%n
      tmp%n = m
      !----- 1. Identity and density ------------------------------------------------------------!
      tmp%pft(1:m)         = cohort%pft(1:m)
      tmp%global_id(1:m)   = cohort%global_id(1:m)
      tmp%owner_patch(1:m) = cohort%owner_patch(1:m)
      tmp%nplant(1:m)      = cohort%nplant(1:m)
      !----- 2. Carbon pools --------------------------------------------------------------------!
      tmp%leaf_carbon(1:m)          = cohort%leaf_carbon(1:m)
      tmp%fineroot_carbon(1:m)      = cohort%fineroot_carbon(1:m)
      tmp%wood_carbon(1:m)          = cohort%wood_carbon(1:m)
      tmp%nonstructural_carbon(1:m) = cohort%nonstructural_carbon(1:m)
      !----- 3. Size, derived from the carbon pools ---------------------------------------------!
      tmp%dbh(1:m)             = cohort%dbh(1:m)
      tmp%height(1:m)          = cohort%height(1:m)
      tmp%basal_area(1:m)      = cohort%basal_area(1:m)
      tmp%agb(1:m)             = cohort%agb(1:m)
      tmp%leaf_area(1:m)       = cohort%leaf_area(1:m)
      tmp%wood_area(1:m)       = cohort%wood_area(1:m)
      tmp%sapwood_carbon(1:m)  = cohort%sapwood_carbon(1:m)
      tmp%sapwood_area(1:m)    = cohort%sapwood_area(1:m)
      tmp%overtopping_lai(1:m) = cohort%overtopping_lai(1:m)
      !----- 4. Water ---------------------------------------------------------------------------!
      tmp%leaf_water_mass(1:m) = cohort%leaf_water_mass(1:m)
      tmp%wood_water_mass(1:m) = cohort%wood_water_mass(1:m)
      tmp%leaf_surf_water(1:m) = cohort%leaf_surf_water(1:m)
      tmp%wood_surf_water(1:m) = cohort%wood_surf_water(1:m)
      tmp%dmax_psi_leaf(1:m)   = cohort%dmax_psi_leaf(1:m)
      !----- 5. Energy --------------------------------------------------------------------------!
      tmp%leaf_energy(1:m) = cohort%leaf_energy(1:m)
      tmp%wood_energy(1:m) = cohort%wood_energy(1:m)
      tmp%leaf_temp(1:m)   = cohort%leaf_temp(1:m)
      tmp%wood_temp(1:m)   = cohort%wood_temp(1:m)
      !----- 6. Plastic leaf traits -------------------------------------------------------------!
      tmp%sla(1:m)                = cohort%sla(1:m)
      tmp%vcmax25(1:m)            = cohort%vcmax25(1:m)
      tmp%rd25(1:m)               = cohort%rd25(1:m)
      tmp%leaf_turnover_rate(1:m) = cohort%leaf_turnover_rate(1:m)
      !----- 7. Leaf phenology ------------------------------------------------------------------!
      tmp%leaf_flush_tendency(1:m) = cohort%leaf_flush_tendency(1:m)
      tmp%leaf_shed_tendency(1:m)  = cohort%leaf_shed_tendency(1:m)
      tmp%growing_degree_days(1:m) = cohort%growing_degree_days(1:m)
      tmp%cold_degree_days(1:m)    = cohort%cold_degree_days(1:m)
      tmp%dry_psi_sum(1:m)         = cohort%dry_psi_sum(1:m)
      tmp%wet_psi_sum(1:m)         = cohort%wet_psi_sum(1:m)
      tmp%light_hours_mean(1:m)    = cohort%light_hours_mean(1:m)
      !----- 8. Growth history ------------------------------------------------------------------!
      tmp%growth_avg(1:m)    = cohort%growth_avg(1:m)
      tmp%growth_accum(1:m)  = cohort%growth_accum(1:m)
      tmp%growth_count(1:m)  = cohort%growth_count(1:m)
      tmp%growth_hist(:,1:m) = cohort%growth_hist(:,1:m)
      !----- 9. Fast-to-slow accumulators -------------------------------------------------------!
      tmp%gpp_accum(1:m)           = cohort%gpp_accum(1:m)
      tmp%leaf_resp_accum(1:m)     = cohort%leaf_resp_accum(1:m)
      tmp%stem_resp_accum(1:m)     = cohort%stem_resp_accum(1:m)
      tmp%root_resp_accum(1:m)     = cohort%root_resp_accum(1:m)
      tmp%light_hours_accum(1:m)   = cohort%light_hours_accum(1:m)
      tmp%dmax_psi_leaf_accum(1:m) = cohort%dmax_psi_leaf_accum(1:m)
      !----- 10. Output accumulators ------------------------------------------------------------!
      !----- Carry the (already grown) diagnostic block over the fresh tmp, whose own diag is      !
      !      inactive -- move_cohort_block copies tmp INTO cohort, so it must be the one holding it. !
      tmp%diag = cohort%diag
      tmp%sdiag = cohort%sdiag
      call move_cohort_block(tmp, cohort)
   end subroutine expand_cohort_block

   !----- Move every allocatable from src into dst (src is reset). ------------------------!
   subroutine move_cohort_block(src, dst)
      type(cohort_block), intent(inout) :: src
      type(cohort_block), intent(inout) :: dst
      dst%n = src%n
      dst%cap = src%cap
      !----- 1. Identity and density ------------------------------------------------------------!
      call move_alloc(src%pft, dst%pft)
      call move_alloc(src%global_id, dst%global_id)
      call move_alloc(src%owner_patch, dst%owner_patch)
      call move_alloc(src%nplant, dst%nplant)
      !----- 2. Carbon pools --------------------------------------------------------------------!
      call move_alloc(src%leaf_carbon, dst%leaf_carbon)
      call move_alloc(src%fineroot_carbon, dst%fineroot_carbon)
      call move_alloc(src%wood_carbon, dst%wood_carbon)
      call move_alloc(src%nonstructural_carbon, dst%nonstructural_carbon)
      !----- 3. Size, derived from the carbon pools ---------------------------------------------!
      call move_alloc(src%dbh, dst%dbh)
      call move_alloc(src%height, dst%height)
      call move_alloc(src%basal_area, dst%basal_area)
      call move_alloc(src%agb, dst%agb)
      call move_alloc(src%leaf_area, dst%leaf_area)
      call move_alloc(src%wood_area, dst%wood_area)
      call move_alloc(src%sapwood_carbon, dst%sapwood_carbon)
      call move_alloc(src%sapwood_area, dst%sapwood_area)
      call move_alloc(src%overtopping_lai, dst%overtopping_lai)
      !----- 4. Water ---------------------------------------------------------------------------!
      call move_alloc(src%leaf_water_mass, dst%leaf_water_mass)
      call move_alloc(src%wood_water_mass, dst%wood_water_mass)
      call move_alloc(src%leaf_surf_water, dst%leaf_surf_water)
      call move_alloc(src%wood_surf_water, dst%wood_surf_water)
      call move_alloc(src%dmax_psi_leaf, dst%dmax_psi_leaf)
      !----- 5. Energy --------------------------------------------------------------------------!
      call move_alloc(src%leaf_energy, dst%leaf_energy)
      call move_alloc(src%wood_energy, dst%wood_energy)
      call move_alloc(src%leaf_temp, dst%leaf_temp)
      call move_alloc(src%wood_temp, dst%wood_temp)
      !----- 6. Plastic leaf traits -------------------------------------------------------------!
      call move_alloc(src%sla, dst%sla)
      call move_alloc(src%vcmax25, dst%vcmax25)
      call move_alloc(src%rd25, dst%rd25)
      call move_alloc(src%leaf_turnover_rate, dst%leaf_turnover_rate)
      !----- 7. Leaf phenology ------------------------------------------------------------------!
      call move_alloc(src%leaf_flush_tendency, dst%leaf_flush_tendency)
      call move_alloc(src%leaf_shed_tendency, dst%leaf_shed_tendency)
      call move_alloc(src%growing_degree_days, dst%growing_degree_days)
      call move_alloc(src%cold_degree_days, dst%cold_degree_days)
      call move_alloc(src%dry_psi_sum, dst%dry_psi_sum)
      call move_alloc(src%wet_psi_sum, dst%wet_psi_sum)
      call move_alloc(src%light_hours_mean, dst%light_hours_mean)
      !----- 8. Growth history ------------------------------------------------------------------!
      call move_alloc(src%growth_avg, dst%growth_avg)
      call move_alloc(src%growth_accum, dst%growth_accum)
      call move_alloc(src%growth_count, dst%growth_count)
      call move_alloc(src%growth_hist, dst%growth_hist)
      !----- 9. Fast-to-slow accumulators -------------------------------------------------------!
      call move_alloc(src%gpp_accum, dst%gpp_accum)
      call move_alloc(src%leaf_resp_accum, dst%leaf_resp_accum)
      call move_alloc(src%stem_resp_accum, dst%stem_resp_accum)
      call move_alloc(src%root_resp_accum, dst%root_resp_accum)
      call move_alloc(src%light_hours_accum, dst%light_hours_accum)
      call move_alloc(src%dmax_psi_leaf_accum, dst%dmax_psi_leaf_accum)
   end subroutine move_cohort_block

   !=======================================================================================!
   !  Centralized lockstep reorder: place old index perm(k) at new position k, for k=1..m.  !
   !  Array assignment with a vector subscript on the RHS is evaluated to a temporary       !
   !  first, so self-permutation is safe. THIS is the one place that knows every field.     !
   !=======================================================================================!
   subroutine reorder_cohort_block(cohort, perm, m)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: perm(:)
      integer(ik),        intent(in)    :: m
      !----- 1. Identity and density ------------------------------------------------------------!
      cohort%pft(1:m)         = cohort%pft(perm(1:m))
      cohort%global_id(1:m)   = cohort%global_id(perm(1:m))
      cohort%owner_patch(1:m) = cohort%owner_patch(perm(1:m))
      cohort%nplant(1:m)      = cohort%nplant(perm(1:m))
      !----- 2. Carbon pools --------------------------------------------------------------------!
      cohort%leaf_carbon(1:m)          = cohort%leaf_carbon(perm(1:m))
      cohort%fineroot_carbon(1:m)      = cohort%fineroot_carbon(perm(1:m))
      cohort%wood_carbon(1:m)          = cohort%wood_carbon(perm(1:m))
      cohort%nonstructural_carbon(1:m) = cohort%nonstructural_carbon(perm(1:m))
      !----- 3. Size, derived from the carbon pools ---------------------------------------------!
      cohort%dbh(1:m)             = cohort%dbh(perm(1:m))
      cohort%height(1:m)          = cohort%height(perm(1:m))
      cohort%basal_area(1:m)      = cohort%basal_area(perm(1:m))
      cohort%agb(1:m)             = cohort%agb(perm(1:m))
      cohort%leaf_area(1:m)       = cohort%leaf_area(perm(1:m))
      cohort%wood_area(1:m)       = cohort%wood_area(perm(1:m))
      cohort%sapwood_carbon(1:m)  = cohort%sapwood_carbon(perm(1:m))
      cohort%sapwood_area(1:m)    = cohort%sapwood_area(perm(1:m))
      cohort%overtopping_lai(1:m) = cohort%overtopping_lai(perm(1:m))
      !----- 4. Water ---------------------------------------------------------------------------!
      cohort%leaf_water_mass(1:m) = cohort%leaf_water_mass(perm(1:m))
      cohort%wood_water_mass(1:m) = cohort%wood_water_mass(perm(1:m))
      cohort%leaf_surf_water(1:m) = cohort%leaf_surf_water(perm(1:m))
      cohort%wood_surf_water(1:m) = cohort%wood_surf_water(perm(1:m))
      cohort%dmax_psi_leaf(1:m)   = cohort%dmax_psi_leaf(perm(1:m))
      !----- 5. Energy --------------------------------------------------------------------------!
      cohort%leaf_energy(1:m) = cohort%leaf_energy(perm(1:m))
      cohort%wood_energy(1:m) = cohort%wood_energy(perm(1:m))
      cohort%leaf_temp(1:m)   = cohort%leaf_temp(perm(1:m))
      cohort%wood_temp(1:m)   = cohort%wood_temp(perm(1:m))
      !----- 6. Plastic leaf traits -------------------------------------------------------------!
      cohort%sla(1:m)                = cohort%sla(perm(1:m))
      cohort%vcmax25(1:m)            = cohort%vcmax25(perm(1:m))
      cohort%rd25(1:m)               = cohort%rd25(perm(1:m))
      cohort%leaf_turnover_rate(1:m) = cohort%leaf_turnover_rate(perm(1:m))
      !----- 7. Leaf phenology ------------------------------------------------------------------!
      cohort%leaf_flush_tendency(1:m) = cohort%leaf_flush_tendency(perm(1:m))
      cohort%leaf_shed_tendency(1:m)  = cohort%leaf_shed_tendency(perm(1:m))
      cohort%growing_degree_days(1:m) = cohort%growing_degree_days(perm(1:m))
      cohort%cold_degree_days(1:m)    = cohort%cold_degree_days(perm(1:m))
      cohort%dry_psi_sum(1:m)         = cohort%dry_psi_sum(perm(1:m))
      cohort%wet_psi_sum(1:m)         = cohort%wet_psi_sum(perm(1:m))
      cohort%light_hours_mean(1:m)    = cohort%light_hours_mean(perm(1:m))
      !----- 8. Growth history ------------------------------------------------------------------!
      cohort%growth_avg(1:m)    = cohort%growth_avg(perm(1:m))
      cohort%growth_accum(1:m)  = cohort%growth_accum(perm(1:m))
      cohort%growth_count(1:m)  = cohort%growth_count(perm(1:m))
      cohort%growth_hist(:,1:m) = cohort%growth_hist(:,perm(1:m))
      !----- 9. Fast-to-slow accumulators -------------------------------------------------------!
      cohort%gpp_accum(1:m)           = cohort%gpp_accum(perm(1:m))
      cohort%leaf_resp_accum(1:m)     = cohort%leaf_resp_accum(perm(1:m))
      cohort%stem_resp_accum(1:m)     = cohort%stem_resp_accum(perm(1:m))
      cohort%root_resp_accum(1:m)     = cohort%root_resp_accum(perm(1:m))
      cohort%light_hours_accum(1:m)   = cohort%light_hours_accum(perm(1:m))
      cohort%dmax_psi_leaf_accum(1:m) = cohort%dmax_psi_leaf_accum(perm(1:m))
      !----- 10. Output accumulators ------------------------------------------------------------!
      !----- The diagnostic accumulators ride the SAME permutation. One call, and it cannot omit  !
      !      a field: they are rows of one 2-D array (ed_diag_vars, principle 4).  --------        !
      call reorder_cohort_diag(cohort%diag,  perm, m)
      call reorder_cohort_diag(cohort%sdiag, perm, m)
      cohort%n = m
   end subroutine reorder_cohort_block

   !----- Compact by keep-mask (length n): drop entries where keep is .false. -------------!
   subroutine compact_cohort_block(cohort, keep)
      type(cohort_block), intent(inout) :: cohort
      logical,            intent(in)    :: keep(:)
      integer(ik), allocatable :: perm(:)
      integer(ik)              :: i, m
      allocate(perm(cohort%n))
      m = 0_ik
      do i = 1_ik, cohort%n
         if (keep(i)) then
            m = m + 1_ik
            perm(m) = i
         end if
      end do
      call reorder_cohort_block(cohort, perm, m)
   end subroutine compact_cohort_block

   !----- Copy every per-cohort field from slot src to slot dst (centralized). ------------!
   subroutine copy_cohort(cohort, dst, src)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: dst, src
      !----- 1. Identity and density ------------------------------------------------------------!
      cohort%pft(dst)         = cohort%pft(src)
      cohort%global_id(dst)   = cohort%global_id(src)
      cohort%owner_patch(dst) = cohort%owner_patch(src)
      cohort%nplant(dst)      = cohort%nplant(src)
      !----- 2. Carbon pools --------------------------------------------------------------------!
      cohort%leaf_carbon(dst)          = cohort%leaf_carbon(src)
      cohort%fineroot_carbon(dst)      = cohort%fineroot_carbon(src)
      cohort%wood_carbon(dst)          = cohort%wood_carbon(src)
      cohort%nonstructural_carbon(dst) = cohort%nonstructural_carbon(src)
      !----- 3. Size, derived from the carbon pools ---------------------------------------------!
      cohort%dbh(dst)             = cohort%dbh(src)
      cohort%height(dst)          = cohort%height(src)
      cohort%basal_area(dst)      = cohort%basal_area(src)
      cohort%agb(dst)             = cohort%agb(src)
      cohort%leaf_area(dst)       = cohort%leaf_area(src)
      cohort%wood_area(dst)       = cohort%wood_area(src)
      cohort%sapwood_carbon(dst)  = cohort%sapwood_carbon(src)
      cohort%sapwood_area(dst)    = cohort%sapwood_area(src)
      cohort%overtopping_lai(dst) = cohort%overtopping_lai(src)
      !----- 4. Water ---------------------------------------------------------------------------!
      cohort%leaf_water_mass(dst) = cohort%leaf_water_mass(src)
      cohort%wood_water_mass(dst) = cohort%wood_water_mass(src)
      cohort%leaf_surf_water(dst) = cohort%leaf_surf_water(src)
      cohort%wood_surf_water(dst) = cohort%wood_surf_water(src)
      cohort%dmax_psi_leaf(dst)   = cohort%dmax_psi_leaf(src)
      !----- 5. Energy --------------------------------------------------------------------------!
      cohort%leaf_energy(dst) = cohort%leaf_energy(src)
      cohort%wood_energy(dst) = cohort%wood_energy(src)
      cohort%leaf_temp(dst)   = cohort%leaf_temp(src)
      cohort%wood_temp(dst)   = cohort%wood_temp(src)
      !----- 6. Plastic leaf traits -------------------------------------------------------------!
      cohort%sla(dst)                = cohort%sla(src)
      cohort%vcmax25(dst)            = cohort%vcmax25(src)
      cohort%rd25(dst)               = cohort%rd25(src)
      cohort%leaf_turnover_rate(dst) = cohort%leaf_turnover_rate(src)
      !----- 7. Leaf phenology ------------------------------------------------------------------!
      cohort%leaf_flush_tendency(dst) = cohort%leaf_flush_tendency(src)
      cohort%leaf_shed_tendency(dst)  = cohort%leaf_shed_tendency(src)
      cohort%growing_degree_days(dst) = cohort%growing_degree_days(src)
      cohort%cold_degree_days(dst)    = cohort%cold_degree_days(src)
      cohort%dry_psi_sum(dst)         = cohort%dry_psi_sum(src)
      cohort%wet_psi_sum(dst)         = cohort%wet_psi_sum(src)
      cohort%light_hours_mean(dst)    = cohort%light_hours_mean(src)
      !----- 8. Growth history ------------------------------------------------------------------!
      cohort%growth_avg(dst)    = cohort%growth_avg(src)
      cohort%growth_accum(dst)  = cohort%growth_accum(src)
      cohort%growth_count(dst)  = cohort%growth_count(src)
      cohort%growth_hist(:,dst) = cohort%growth_hist(:,src)
      !----- 9. Fast-to-slow accumulators -------------------------------------------------------!
      cohort%gpp_accum(dst)           = cohort%gpp_accum(src)
      cohort%leaf_resp_accum(dst)     = cohort%leaf_resp_accum(src)
      cohort%stem_resp_accum(dst)     = cohort%stem_resp_accum(src)
      cohort%root_resp_accum(dst)     = cohort%root_resp_accum(src)
      cohort%light_hours_accum(dst)   = cohort%light_hours_accum(src)
      cohort%dmax_psi_leaf_accum(dst) = cohort%dmax_psi_leaf_accum(src)
      !----- 10. Output accumulators ------------------------------------------------------------!
      call copy_cohort_diag(cohort%diag,  dst, src)
      call copy_cohort_diag(cohort%sdiag, dst, src)
   end subroutine copy_cohort

   !=======================================================================================!
   !  SIZE: the geometry derived from the carbon pools                                     !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Place ONE cohort slot on the allometry from its dbh (setup, census, recruitment, a state  !
   ! file without pools): the four carbon pools a plant of that dbh holds, then the geometry   !
   ! derived from those pools by set_cohort_size_from_carbon, exactly as after every other     !
   ! carbon change. Geometry has that one writer, so a cohort's geometry is always the same    !
   ! function of its carbon and a restart that re-derives it gets the same bits.               !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp)    :: height
      integer(ik) :: pf
      pf = cohort%pft(i)
      associate (pft => cohort%pft_table)
         height = dbh_to_height(pft%allom(pf), cohort%dbh(i))
         cohort%leaf_carbon(i)          = size_to_leaf_carbon(pft%allom(pf), cohort%dbh(i), height, cohort%sla(i))
         cohort%wood_carbon(i)          = size_to_wood_carbon(pft%allom(pf), cohort%dbh(i), height,             &
                                                           pft%wood_density(pf))
         cohort%fineroot_carbon(i)      = pft%root_to_leaf_ratio(pf) * cohort%leaf_carbon(i)
         cohort%nonstructural_carbon(i) = storage_carbon_target(pft%storage_cushion(pf), cohort%leaf_carbon(i), &
                                                                cohort%fineroot_carbon(i))
      end associate
      call set_cohort_size_from_carbon(cohort, i)
   end subroutine set_cohort_size

   !---------------------------------------------------------------------------------------!
   ! The ONE writer of a cohort's geometry: with wood_carbon the size anchor, derive dbh from   !
   ! it (wood_to_dbh) and the rest of the geometry, and take leaf_area straight from leaf_carbon. !
   ! Called after every operator that moves carbon or SLA (birth, the daily update, trait       !
   ! acclimation, fusion, fission, a restart read). The pools are the inputs and are not        !
   ! changed.                                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size_from_carbon(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      !----- The carbon->geometry flip is the single shared kernel carbon_to_structure             !
      !       (ed_allometry_lib); this per-slot wrapper just packs/unpacks the SoA fields.         !
      integer(ik) :: pf
      pf = cohort%pft(i)
      call carbon_to_structure(cohort%pft_table%allom(pf), cohort%wood_carbon(i), cohort%leaf_carbon(i), &
                               cohort%pft_table%wood_density(pf), cohort%sla(i),                 &
                               cohort%dbh(i), cohort%height(i), cohort%basal_area(i),             &
                               cohort%agb(i), cohort%leaf_area(i))
      call set_cohort_wood_allometry(cohort, i)
   end subroutine set_cohort_size_from_carbon


   !---------------------------------------------------------------------------------------!
   ! The cached WOOD geometry of one slot, from the just-set dbh / basal_area / wood_carbon.  !
   ! Part of set_cohort_size_from_carbon, so the leaf and wood caches can never disagree about !
   ! which dbh they were derived from. Everything here is PER PLANT (see the declarations):   !
   ! `wood_area` is the per-plant twin of `leaf_area`, so WAI = nplant*wood_area at the point  !
   ! of use and no cached field carries a stale plant density.                                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine set_cohort_wood_allometry(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp)    :: f_sap
      integer(ik) :: pf
      pf = cohort%pft(i)
      associate (pft => cohort%pft_table)
         !----- dbh_to_wai is exactly linear in nplant, so passing 1 gives the per-plant area. ---!
         cohort%wood_area(i)      = dbh_to_wai(pft%allom(pf), cohort%dbh(i), 1.0_wp)
         f_sap                    = sapwood_fraction(pft%allom(pf), cohort%dbh(i))
         cohort%sapwood_carbon(i) = f_sap * cohort%wood_carbon(i)
         cohort%sapwood_area(i)   = f_sap * cohort%basal_area(i)
      end associate
   end subroutine set_cohort_wood_allometry

   !=======================================================================================!
   !  TISSUE WATER AND HEAT                                                                !
   !=======================================================================================!

   !----- The carbon whose water a cohort's wood store holds [kgC/plant]: the sapwood ring and the   !
   !      fine roots, one store in the plant-water network (docs/science/soil_plant_water.md §1). ---!
   elemental pure real(wp) function cohort_wood_store_carbon(sapwood_carbon, fineroot_carbon) result(c)
      real(wp), intent(in) :: sapwood_carbon, fineroot_carbon   !< [kgC/plant]
      c = sapwood_carbon + fineroot_carbon
   end function cohort_wood_store_carbon

   !----- Fill a cohort's leaf and wood stores with the water they hold at potential `psi` on their  !
   !      own curves (docs/science/plant_hydraulics.md eq. (5)): a cohort's water at birth. ---------!
   pure subroutine set_cohort_tissue_water(cohort, i, psi, leaf_curve, wood_curve)
      type(cohort_block),  intent(inout) :: cohort
      integer(ik),         intent(in)    :: i
      real(wp),            intent(in)    :: psi                    !< [MPa]
      type(water_curve_t), intent(in)    :: leaf_curve, wood_curve
      cohort%leaf_water_mass(i) = tissue_psi_to_water(psi, leaf_curve, cohort%leaf_carbon(i))
      cohort%wood_water_mass(i) = tissue_psi_to_water(psi, wood_curve,                                     &
                                     cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
   end subroutine set_cohort_tissue_water

   !----- One cohort's INTERNAL (hydraulic) tissue water [kg/m2 ground]: the stores the plant    !
   !       hydraulics carries, both of which its heat capacity counts. --------------------------!
   pure function cohort_tissue_water(cohort, i) result(w)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      real(wp)                       :: w
      w = cohort%nplant(i) * (max(cohort%leaf_water_mass(i), 0.0_wp)                               &
                            + max(cohort%wood_water_mass(i), 0.0_wp))
   end function cohort_tissue_water

   !----- One cohort's canopy films [kg/m2 ground]: the water intercepted on its leaves and wood. --!
   pure function cohort_film_water(cohort, i) result(w)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      real(wp)                       :: w
      w = cohort%nplant(i) * (cohort%leaf_surf_water(i) + cohort%wood_surf_water(i))
   end function cohort_film_water

   !----- The enthalpy of the water a cohort holds [J/plant]: its tissue water and its canopy films,  !
   !      each at the temperature of the tissue it is in or on. What a shed hands the ground with  !
   !      the water. ------------------------------------------------------------------------------!
   pure function cohort_water_enthalpy(cohort, i) result(h)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      real(wp)                       :: h
      h = (max(cohort%leaf_water_mass(i), 0.0_wp) + cohort%leaf_surf_water(i))                    &
          * internal_energy_liquid(cohort%leaf_temp(i))                                           &
        + (max(cohort%wood_water_mass(i), 0.0_wp) + cohort%wood_surf_water(i))                    &
          * internal_energy_liquid(cohort%wood_temp(i))
   end function cohort_water_enthalpy

   !----- One cohort's leaf and wood heat capacity [J/m2 ground/K], from tissue_heat_capacity. The !
   !      wood holds the whole wood store's water, its fine roots' share included, and each tissue !
   !      the film it carries. The fusion weights, the demography operators and the slow ledger    !
   !      all take it from here, so they agree on a cohort's thermal mass. ------------------------!
   pure subroutine cohort_tissue_heat_capacity(cohort, i, cap_leaf, cap_wood)
      type(cohort_block), intent(in)  :: cohort
      integer(ik),        intent(in)  :: i
      real(wp),           intent(out) :: cap_leaf, cap_wood             !< [J/m2 ground/K]
      cap_leaf = tissue_heat_capacity(cohort%nplant(i), cohort%leaf_carbon(i), tissue_c_leaf,           &
                                      max(cohort%leaf_water_mass(i), 0.0_wp) + cohort%leaf_surf_water(i))
      cap_wood = tissue_heat_capacity(cohort%nplant(i), cohort%wood_carbon(i), tissue_c_wood_dry,       &
                                      max(cohort%wood_water_mass(i), 0.0_wp) + cohort%wood_surf_water(i))
   end subroutine cohort_tissue_heat_capacity

   !----- A cohort's per-plant tissue masses, on which its energies are counted: each tissue's dry   !
   !      heat capacity [J/K/plant] and its water [kg/plant], its own and the film it carries         !
   !      (docs/science/vegetation_energy_dynamics.md eq. (8)). -------------------------------------!
   pure subroutine cohort_tissue_masses(cohort, i, dry_leaf, w_leaf, dry_wood, w_wood)
      type(cohort_block), intent(in)  :: cohort
      integer(ik),        intent(in)  :: i
      real(wp),           intent(out) :: dry_leaf, w_leaf, dry_wood, w_wood
      dry_leaf = tissue_dry_heat_capacity(1.0_wp, cohort%leaf_carbon(i), tissue_c_leaf)
      dry_wood = tissue_dry_heat_capacity(1.0_wp, cohort%wood_carbon(i), tissue_c_wood_dry)
      w_leaf   = max(cohort%leaf_water_mass(i), 0.0_wp) + cohort%leaf_surf_water(i)
      w_wood   = max(cohort%wood_water_mass(i), 0.0_wp) + cohort%wood_surf_water(i)
   end subroutine cohort_tissue_masses

   !----- The tissue temperatures a cohort's energies give with the masses it holds (eq. (10)). A    !
   !      leaf with neither dry matter nor water -- a leafless crown -- takes its wood's temperature,  !
   !      as a bud takes its twig's; wood with no heat capacity keeps the temperature it has. -------!
   pure subroutine cohort_tissue_temps(cohort, i, t_leaf, t_wood)
      type(cohort_block), intent(in)  :: cohort
      integer(ik),        intent(in)  :: i
      real(wp),           intent(out) :: t_leaf, t_wood   !< [K]
      real(wp) :: dry_leaf, w_leaf, dry_wood, w_wood
      call cohort_tissue_masses(cohort, i, dry_leaf, w_leaf, dry_wood, w_wood)
      t_wood = cohort%wood_temp(i)
      if (dry_wood + w_wood * cp_liq > 0.0_wp) t_wood = tissue_energy_to_temp(cohort%wood_energy(i), w_wood, dry_wood)
      t_leaf = t_wood
      if (dry_leaf + w_leaf * cp_liq > 0.0_wp) t_leaf = tissue_energy_to_temp(cohort%leaf_energy(i), w_leaf, dry_leaf)
   end subroutine cohort_tissue_temps

   !----- Set a cohort's tissue temperatures from its energies. The temperatures are derived state:   !
   !      after birth only this and hold_cohort_tissue_temps write them. ---------------------------!
   pure subroutine set_cohort_tissue_temps(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp) :: t_leaf, t_wood
      call cohort_tissue_temps(cohort, i, t_leaf, t_wood)
      cohort%leaf_temp(i) = t_leaf
      cohort%wood_temp(i) = t_wood
   end subroutine set_cohort_tissue_temps

   !----- Hold a cohort's tissue temperatures through a change of its masses, for matter that came or !
   !      went at the tissue's own temperature: each energy becomes what its tissue holds at its       !
   !      temperature with the masses it now has, and the temperatures are then taken back from the   !
   !      energies, so they stay exactly what the energies give. -------------------------------------!
   pure subroutine hold_cohort_tissue_temps(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp) :: dry_leaf, w_leaf, dry_wood, w_wood
      call cohort_tissue_masses(cohort, i, dry_leaf, w_leaf, dry_wood, w_wood)
      cohort%leaf_energy(i) = tissue_temp_to_energy(dry_leaf, w_leaf, cohort%leaf_temp(i))
      cohort%wood_energy(i) = tissue_temp_to_energy(dry_wood, w_wood, cohort%wood_temp(i))
      call set_cohort_tissue_temps(cohort, i)
   end subroutine hold_cohort_tissue_temps

   !=======================================================================================!
   !  BIRTH                                                                                !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Initialize ONE freshly-created cohort slot at BIRTH from (pft, patch, dbh): set the      !
   ! per-slot identity + prognostic fields to their birth values, then derive                 !
   ! the on-allometry geometry via set_cohort_size. The SINGLE canonical birth-field policy    !
   ! shared by add_cohort (setup / census) and apply_recruitment -- add any new per-cohort     !
   ! field that needs a birth value HERE. overtopping_lai is zeroed (a reused, culled slot may  !
   ! hold a stale value; it is recomputed by update_overtopping_lai every slow step anyway).    !
   ! The caller ensures capacity, stamps the global id, and bumps cohort%n.                    !
   !---------------------------------------------------------------------------------------!
   subroutine init_cohort(cohort, m, pft, ipft, owner_patch, nplant, dbh, birth_temp, birth_psi)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: m, ipft, owner_patch
      type(pft_table_t),  intent(in)    :: pft
      real(wp),           intent(in)    :: nplant, dbh
      !----- The state the cohort is born in, its patch's (patch_physical_state): its leaves and wood !
      !      start at birth_temp, holding the water their tissues hold at birth_psi, and with the    !
      !      energy those tissues then hold. ------------------------------------------------------!
      real(wp),           intent(in)    :: birth_temp   !< [K]
      real(wp),           intent(in)    :: birth_psi    !< [MPa]
      cohort%pft(m)              = ipft
      cohort%owner_patch(m)      = owner_patch
      cohort%nplant(m)           = nplant
      cohort%dbh(m)              = dbh
      cohort%growth_avg(m)       = UNSET_POSITIVE     ! read once growth_count > 0
      !----- UNSET: the fast loop seeds it from the surface-layer soil potential on the cohort's  !
      !      first step, so a recruit starts at its patch's actual water status, not at 0 (turgid). !
      cohort%dmax_psi_leaf(m)   = UNSET_NEGATIVE
      cohort%dmax_psi_leaf_accum(m) = DMAX_PSI_LEAF_ACCUM_RESET
      cohort%growth_accum(m)     = 0.0_wp
      cohort%growth_count(m)     = 0_ik
      cohort%growth_hist(:, m)   = 0.0_wp             ! a reused slot keeps a dead cohort's samples otherwise
      cohort%overtopping_lai(m)  = 0.0_wp             ! fresh competition context (recomputed each slow step)
      cohort%light_hours_accum(m)        = 0.0_wp             ! no light seen yet this step
      call clear_cohort_diag(cohort%diag,  m)    ! fresh diagnostics (slot may be a reused, stale cull)
      call clear_cohort_diag(cohort%sdiag, m)
      cohort%leaf_flush_tendency(m) = PHENO_FLUSH_INIT  ! born flushing
      cohort%leaf_shed_tendency(m)  = PHENO_SHED_INIT   ! not senescing at birth
      cohort%growing_degree_days(m) = 0.0_wp            ! fresh phenology memory
      cohort%cold_degree_days(m)    = 0.0_wp
      cohort%dry_psi_sum(m)         = 0.0_wp
      cohort%wet_psi_sum(m)         = 0.0_wp
      cohort%light_hours_mean(m) = UNSET_POSITIVE  ! no light memory: its first day sets it
      cohort%sla(m)                = pft%sla(ipft)
      cohort%vcmax25(m)              = pft%vcmax25(ipft)          ! dynamic leaf traits born at
      cohort%rd25(m)                 = pft%rd25(ipft)             ! top-of-canopy; plasticity relaxes
      cohort%leaf_turnover_rate(m)   = pft%leaf_turnover_rate_toc(ipft) ! them toward the shaded target
      cohort%leaf_temp(m)        = birth_temp         ! fresh fast state (slot may be a reused, stale cull)
      cohort%wood_temp(m)        = birth_temp
      cohort%leaf_surf_water(m)  = 0.0_wp             ! a new cohort's canopy starts dry
      cohort%wood_surf_water(m)  = 0.0_wp
      call set_cohort_size(cohort, m)                 ! carbon pools from dbh, then the geometry from them
      call set_cohort_tissue_water(cohort, m, birth_psi, pft%hydro(ipft)%leaf_curve, pft%hydro(ipft)%wood_curve)
      call hold_cohort_tissue_temps(cohort, m)
   end subroutine init_cohort

   !=======================================================================================!
   !  FUSION AND FISSION: the per-field policy for the fast-loop-owned state               !
   !=======================================================================================!

   !=========================================================================================!
   !  THE PER-FIELD FUSION POLICY OF THE FAST-LOOP-OWNED COHORT STATE -- declared once,      !
   !  here.                                                                                  !
   !                                                                                         !
   !  These fields are written by the fast loop and read back by it; the slow loop only has  !
   !  to carry them correctly through a fusion. Each is a per-PLANT amount (FSK_EXTENSIVE):  !
   !  the tissue energies and waters, the canopy films and the per-plant fluxes, blended by  !
   !  nplant, so the site total (nplant * value, summed) is conserved. The tissue            !
   !  temperatures are not fused: they are derived from the fused energies, by the caller,   !
   !  once the carbon pools are merged. Getting a field's kind wrong is invisible: the AGB   !
   !  assert still passes, the conservation ledgers still close, and the answer is quietly   !
   !  wrong.                                                                                 !
   !                                                                                         !
   !  The diagnostic twins are fused from their own declared table (CO_FUSION_KIND,          !
   !  ed_diag_vars) on weights the caller computes: the tissue heat capacities for           !
   !  the leaf and wood temperatures, and leaf area for a per-LEAF-AREA diagnostic           !
   !  (transpiration per m2 of leaf, say), which is not a temperature. Adding a fast-loop    !
   !  field is one line here; the compiler will not remind you, so the line is the reminder. !
   !=========================================================================================!
   !----- THE per-field policy for the fast-loop cohort state, declared once: every field below is!
   !      a per-plant amount, blended by nplant.                                                  !
   !                                                                                              !
   !      It is centralised but not enforced: a per-cohort fast field left out of this routine    !
   !      compiles clean and is silently wrong at the next fusion. A new field needs a line here  !
   !      and a witness in test_slow_update.   A table the blend iterated would be no safer (the  !
   !      compiler cannot check that a table lists every field either) and would be longer than   !
   !      these lines. ---------------------------------------------------------------------------!
   pure subroutine fuse_cohort_fast_state(cohort, recc, donc, np_r, np_d)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: recc, donc   !< survivor, donor
      real(wp),           intent(in)    :: np_r, np_d   !< nplant weights
      !----- The tissue energies [J/plant]: nplant-weighted, so the merged tissue holds the      !
      !      energy the two cohorts held, exactly. Its temperatures follow from it once its      !
      !      carbon pools are merged (fuse_2_cohorts). ------------------------------------------!
      call blend(cohort%leaf_energy,      np_r, np_d)
      call blend(cohort%wood_energy,      np_r, np_d)
      !----- EXTENSIVE: tissue water and the canopy films [kg/plant], and the per-plant flux     !
      !      accumulators. ----------------------------------------------------------------------!
      call blend(cohort%leaf_water_mass,  np_r, np_d)
      call blend(cohort%wood_water_mass,  np_r, np_d)
      call blend(cohort%gpp_accum,        np_r, np_d)
      call blend(cohort%leaf_resp_accum,  np_r, np_d)
      call blend(cohort%stem_resp_accum,  np_r, np_d)
      call blend(cohort%root_resp_accum,  np_r, np_d)
      call blend(cohort%leaf_surf_water,  np_r, np_d)
      call blend(cohort%wood_surf_water,  np_r, np_d)
      !----- SURVIVOR-KEEPS, declared rather than left implicit: the phenology memories            !
      !      (the two tendencies, the degree-day and psi sums, the light mean), light_hours_accum    !
      !      (already consumed this step) and dmax_psi_leaf are absent from this routine ON PURPOSE,  !
      !      so the survivor keeps its own. Two cohorts only fuse when they already match in size,   !
      !      height and PFT, and every one of these is driven by a signal both of them saw (the PAR at a !
      !      shared height included), so the two memories are near-identical at fusion and a blend would !
      !      move nothing. Blending is also not obviously right for a sum that resets on a trigger:  !
      !      an nplant-weighted average of a sum that has just reset and one that has not is a       !
      !      history neither cohort had. If a future cue is driven by something genuinely           !
      !      per-cohort and slow, revisit this -- and add it here rather than leaving it implicit.   !
   contains
      pure subroutine blend(a, wr, wd)
         real(wp), intent(inout) :: a(:)
         real(wp), intent(in)    :: wr, wd
         a(recc) = pair_mean(a(recc), wr, a(donc), wd)
      end subroutine blend
   end subroutine fuse_cohort_fast_state

   !=======================================================================================!
   !  BOUNDS: the check the driver runs every step                                         !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! The bounds of a cohort's state: its plant density, carbon pools, SLA, and the water it    !
   ! holds (internal stores and the films on its leaves and wood) are amounts, never below     !
   ! zero and never infinite; its tissue energies are finite. Every operator keeps them in     !
   ! bounds by construction -- the applier zeroes only a pool's round-off, mortality never     !
   ! takes nplant below `negligible`, the fast loop clamps the films into [0, capacity] -- so  !
   ! a violation is an operator bug, and a non-finite energy is a fast-loop integrator that    !
   ! could not take even its smallest step. The driver checks this every step and names the   !
   ! field, so the bug is caught the step it happens. The geometry and the temperatures are    !
   ! not checked: they are functions of these. Returns the first bad field's index, 0 if none. !
   !---------------------------------------------------------------------------------------!
   pure function cohort_bad_field(cohort, i) result(k)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      integer(ik)                    :: k
      real(wp)    :: x(N_COHORT_CHECKED)
      integer(ik) :: j
      x = checked_fields(cohort, i)
      k = 0_ik
      do j = 1_ik, N_COHORT_CHECKED
         !----- Not finite first, by classification, as soil_carbon_bad_pool does: a comparison   !
         !      with a NaN raises IEEE invalid, which a trapping build stops on. ------------------!
         if (.not. ieee_is_finite(x(j))) then
            k = j
            return
         end if
         if (j <= N_COHORT_AMOUNTS .and. x(j) < 0.0_wp) then
            k = j
            return
         end if
      end do
   end function cohort_bad_field

   !----- Whether field k (cohort_bad_field's index) of cohort i is finite: the driver reports a  !
   !      non-finite field as a NaN and a finite one out of bounds as a bound. --------------------!
   pure logical function cohort_field_finite(cohort, i, k) result(yes)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i, k
      real(wp) :: x(N_COHORT_CHECKED)
      x = checked_fields(cohort, i)
      yes = ieee_is_finite(x(k))
   end function cohort_field_finite

   !----- Name for the index cohort_bad_field returns, so a failure message is readable. -------!
   pure function cohort_field_name(k) result(nm)
      integer(ik), intent(in) :: k
      character(len=20)       :: nm
      if (k >= 1_ik .and. k <= N_COHORT_CHECKED) then
         nm = COHORT_CHECKED_NAME(k)
      else
         nm = '(none)'
      end if
   end function cohort_field_name

   !----- The fields cohort_bad_field checks, in COHORT_CHECKED_NAME's order. --------------------!
   pure function checked_fields(cohort, i) result(x)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      real(wp)                       :: x(N_COHORT_CHECKED)
      x = [cohort%nplant(i), cohort%leaf_carbon(i), cohort%fineroot_carbon(i),                      &
           cohort%wood_carbon(i), cohort%nonstructural_carbon(i), cohort%sla(i),                    &
           cohort%leaf_water_mass(i), cohort%wood_water_mass(i),                                    &
           cohort%leaf_surf_water(i), cohort%wood_surf_water(i),                                    &
           cohort%leaf_energy(i), cohort%wood_energy(i)]
   end function checked_fields

   !=======================================================================================!
   !  TENDENCIES: the bundle the slow step fills and update_cohort_states applies          !
   !=======================================================================================!


   !----- Allocate (or resize-on-grow) the tendency bundle to hold n cohorts. Resize-on-grow so   !
   !      the site-carried scratch buffer allocates once and then reuses. The arrays are ALWAYS    !
   !      allocated/freed together, so a partially-allocated bundle (e.g. after a piecemeal edit)   !
   !      forces a full re-allocation rather than a stale reuse. -------------------------------!
   subroutine cohort_deriv_alloc(deriv, n)
      type(cohort_deriv_block), intent(inout) :: deriv
      integer(ik),              intent(in)    :: n
      integer(ik) :: m
      logical     :: ok
      m  = max(n, 1_ik)
      ok = allocated(deriv%d_leaf_carbon_dt)       .and. allocated(deriv%d_fineroot_carbon_dt) .and. &
           allocated(deriv%d_wood_carbon_dt)       .and. allocated(deriv%d_nonstructural_carbon_dt) .and. &
           allocated(deriv%dln_nplant_dt)
      if (ok) then
         if (size(deriv%d_leaf_carbon_dt) >= m) return   ! big enough (all co-sized) -> reuse
         deallocate(deriv%d_leaf_carbon_dt, deriv%d_fineroot_carbon_dt, deriv%d_wood_carbon_dt,      &
                    deriv%d_nonstructural_carbon_dt, deriv%dln_nplant_dt)
      else                                       ! none or partial -> drop whatever is there
         if (allocated(deriv%d_leaf_carbon_dt))          deallocate(deriv%d_leaf_carbon_dt)
         if (allocated(deriv%d_fineroot_carbon_dt))      deallocate(deriv%d_fineroot_carbon_dt)
         if (allocated(deriv%d_wood_carbon_dt))          deallocate(deriv%d_wood_carbon_dt)
         if (allocated(deriv%d_nonstructural_carbon_dt)) deallocate(deriv%d_nonstructural_carbon_dt)
         if (allocated(deriv%dln_nplant_dt))             deallocate(deriv%dln_nplant_dt)
      end if
      allocate(deriv%d_leaf_carbon_dt(m), deriv%d_fineroot_carbon_dt(m), deriv%d_wood_carbon_dt(m),  &
               deriv%d_nonstructural_carbon_dt(m), deriv%dln_nplant_dt(m))
   end subroutine cohort_deriv_alloc

   !----- Cohort i's carbon for the slow step [kgC/plant]: the fast loop's GPP and maintenance  !
   !      respiration, or, with the fast loop off, gpp_ref per unit leaf area and no maintenance !
   !      (docs/science/plant_carbon_allocation.md §1). The allocator and the slow ledger both   !
   !      read it, so the ledger declares what the allocator consumes. -----------------------!
   pure subroutine cohort_carbon_supply(cohort, i, fast_biophysics_on, gpp_ref, dt_yr, gross_gpp, resp_maint)
      type(cohort_block), intent(in)  :: cohort
      integer(ik),        intent(in)  :: i
      logical,            intent(in)  :: fast_biophysics_on
      real(wp),           intent(in)  :: gpp_ref      !< [kgC/m2 leaf/yr] the fast-off stub's GPP
      real(wp),           intent(in)  :: dt_yr        !< [yr] the slow step
      real(wp),           intent(out) :: gross_gpp, resp_maint
      if (fast_biophysics_on) then
         gross_gpp  = cohort%gpp_accum(i)
         resp_maint = cohort%leaf_resp_accum(i) + cohort%stem_resp_accum(i) + cohort%root_resp_accum(i)
      else
         gross_gpp  = gpp_ref * cohort%leaf_area(i) * dt_yr
         resp_maint = 0.0_wp
      end if
   end subroutine cohort_carbon_supply

end module ed_cohort_state_vars
