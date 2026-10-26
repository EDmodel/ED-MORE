! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_patch_state_vars -- what one patch holds, and the block of all the site's patches.    !
!                                                                                          !
! A patch holds five stores -- canopy air, snow, soil water, soil energy and soil carbon   !
! -- and two daily accumulators. Each is a derived type, and `patch_block` keeps an array  !
! of each, one element per patch, rather than one array per field as the cohort block      !
! does. A store is used whole and one patch at a time: a kernel takes one patch's store as !
! its argument and sees nothing else of the state, the fast loop copies a patch's stores   !
! into its working column once a day and back, and a fusion or a gap mixes two stores with !
! one blend. No hot loop reads one field across all patches, so an array per field would   !
! gain no speed, and the block's field lists would name about 46 fields instead of 7       !
! stores. The patch's own numbers (area, age, the daily rates) are plain arrays, and its   !
! cohort map is the CSR map explained in ed_cohort_state_vars. The stores' parameters      !
! are config (ed_soil_params), not state.                                                  !
!==========================================================================================!
module ed_patch_state_vars
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : n_soil_layer_max, n_snow_layer_max, tiny_num, t_3ple, rho_h2o
   use ed_therm_lib,          only : cas_molar_density
   use ed_diag_vars,          only : patch_diag_block, expand_patch_diag, reorder_patch_diag
   use ed_state_defaults,     only : UNSET_POSITIVE, CAN_CO2_UNSEEDED, CAN_DEPTH_UNSET
   use ed_soil_derive,        only : top_soil_psi
   use ed_site_state_vars,      only : site_surface_t
   use ed_biogeochem_params,    only : N_SOIL_POOL
   implicit none
   private

   public :: cas_state_t, snow_state_t, soil_water_state_t, soil_energy_state_t, soil_carbon_state_t
   public :: decomp_xi_t, litter_input_t
   public :: blend_cas, cas_set_depth, cas_capacities, blend_snow, blend_soil_water, blend_soil_energy, blend_soil_carbon
   public :: soil_water_store, soil_energy_store
   public :: blend_decomp_xi, blend_litter_input, blend_patch_stores, add_shed_water
   public :: patch_block, alloc_patch_block, expand_patch_block, reorder_patch_block, free_patch_block
   public :: patch_physical_state
   public :: canopy_air_bad_field, canopy_air_field_name
   public :: DIST_PRIMARY, DIST_TREEFALL

   !----- A patch's disturbance class (patch_block%dist_type). ------------------------------------!
   integer(ik), parameter :: DIST_PRIMARY  = 1_ik   !< undisturbed / primary stand
   integer(ik), parameter :: DIST_TREEFALL = 2_ik   !< treefall-gap (age-0) patch

   !----- The canopy air: three prognostic twins, and the depth that sizes the box         !
   !      (docs/science/canopy_air_space_biophysics.md). ----------------------------------!
   type :: cas_state_t
      real(wp) :: can_enthalpy = 0.0_wp                     !< [J/kg]     specific enthalpy (prognostic)
      real(wp) :: can_shv      = 0.0_wp                     !< [kg/kg]    specific humidity (prognostic)
      real(wp) :: can_co2      = CAN_CO2_UNSEEDED           !< [umol/mol] CO2 mixing ratio in dry air (prognostic)
      real(wp) :: can_temp     = UNSET_POSITIVE             !< [K]        diagnosed; unset until a start seeds it
      real(wp) :: can_depth    = CAN_DEPTH_UNSET            !< [m]        tallest cohort + freeboard, set by the slow step
   end type cas_state_t

   !----- Snow, or temporary surface water: mass and energy per ground area, with          !
   !      temperature and liquid fraction read off them; nlayer = 0 is no snow             !
   !      (docs/science/snow_biophysics.md). ----------------------------------------------!
   type :: snow_state_t
      real(wp) :: swe(n_snow_layer_max)         = 0.0_wp   !< [kg/m2] water-equivalent mass (prognostic)
      real(wp) :: snow_energy(n_snow_layer_max) = 0.0_wp   !< [J/m2]  internal energy (prognostic)
      real(wp) :: snow_depth(n_snow_layer_max)  = 0.0_wp   !< [m]     swe over the snow density
      real(wp) :: snow_temp(n_snow_layer_max)   = t_3ple   !< [K]     diagnosed
      real(wp) :: snow_fliq(n_snow_layer_max)   = 0.0_wp   !< [-]     liquid fraction, diagnosed
      integer(ik) :: nlayer                     = 0_ik     !< layers in use (0 = no snow)
   end type snow_state_t

   !----- Soil water, and the pond on top of it. The pond holds its own energy, so water   !
   !      that crosses into or out of it carries heat (docs/science/numerical_scheme.md,   !
   !      "The pond holds heat"). ---------------------------------------------------------!
   type :: soil_water_state_t
      real(wp) :: theta(n_soil_layer_max) = 0.0_wp   !< [m3/m3] volumetric water content (prognostic)
      real(wp) :: w_surface = 0.0_wp                 !< [kg/m2] ponded water (prognostic)
      real(wp) :: w_surface_energy = 0.0_wp          !< [J/m2]  ponded water's internal energy (prognostic)
   end type soil_water_state_t

   !----- Soil heat: internal energy per volume of each layer, with temperature and liquid !
   !      fraction read off it. -----------------------------------------------------------!
   type :: soil_energy_state_t
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp    !< [J/m3] internal energy (prognostic)
      real(wp) :: soil_temp(n_soil_layer_max)   = 0.0_wp    !< [K]    diagnosed
      real(wp) :: soil_fliq(n_soil_layer_max)   = 1.0_wp    !< [-]    liquid fraction, diagnosed
   end type soil_energy_state_t

   !----- The soil-carbon pools, the matrix state vector X (docs/science/soil_carbon.md    !
   !      §2) by IP_* index: written daily by soil_carbon_step, held fixed for the fast    !
   !      loop's heterotrophic respiration. -----------------------------------------------!
   type :: soil_carbon_state_t
      real(wp) :: carbon(N_SOIL_POOL) = 0.0_wp   !< [kgC/m2] by IP_* index, litter -> SOM -> passive
      ! lignin of the structural pools [kgC/m2], above (1) and below (2) ground, 0 <= L <= C:
      ! its fraction slows their decay
      real(wp) :: lignin(2) = 0.0_wp
   end type soil_carbon_state_t

   !----- The day's decomposition scalars for soil_carbon_step: each pool's xi integrated  !
   !      over the day (docs/science/soil_carbon.md §7, eq. (5)), summed by the fast loop, !
   !      by pool (IP_*). -----------------------------------------------------------------!
   type :: decomp_xi_t
      real(wp) :: xi_days(N_SOIL_POOL) = 0.0_wp   !< [day] INT_day xi_j dt
      real(wp) :: rh_fast = 0.0_wp   !< [kgC/m2] today's fast-loop Rh, a check on soil_carbon_step's own
   end type decomp_xi_t

   !----- The day's litter [kgC/m2/day], already split among the pools it enters           !
   !      (docs/science/soil_carbon.md §5). -----------------------------------------------!
   type :: litter_input_t
      real(wp) :: labile_grnd  = 0.0_wp      !< -> X(1) fast_grnd  (labile leaf, fruit; above-ground storage + wood)
      real(wp) :: labile_soil  = 0.0_wp      !< -> X(2) fast_soil  (labile fine root; below-ground storage + wood)
      real(wp) :: struct_grnd  = 0.0_wp      !< -> X(3) struct_grnd (structural leaf + CWD, above)
      real(wp) :: struct_soil  = 0.0_wp      !< -> X(4) struct_soil (structural fineroot + belowground CWD)
      real(wp) :: lignin_grnd  = 0.0_wp      !< lignin flux to struct_grnd (sets f_lignin of incoming litter)
      real(wp) :: lignin_soil  = 0.0_wp      !< lignin flux to struct_soil
      ! the same carbon by the tissue it came from, for output; they sum to the four pool streams.
      ! Leaf litter is counted after resorption, as a litter trap collects it.
      real(wp) :: leaf_c     = 0.0_wp      !< leaves (turnover, drop, dead plants)
      real(wp) :: fineroot_c = 0.0_wp      !< fine roots
      real(wp) :: wood_c     = 0.0_wp      !< wood of dead plants (stem + coarse roots)
      real(wp) :: storage_c  = 0.0_wp      !< non-structural carbon of dead plants
      real(wp) :: repro_c    = 0.0_wp      !< reproduction carbon that made no recruit (fruit, seed)
      real(wp) :: exudate_c  = 0.0_wp      !< root exudate (carbon the growth sink could not use)
   end type litter_input_t

   !----- All patches of a polygon, sorted by site: an array per number, and an array of     !
   !      each store.                                                                       !
   type :: patch_block
      !==========================================================================================!
      ! 0. BLOCK RECORD                                                                          !
      !==========================================================================================!
      integer(ik) :: n   = 0_ik       !< patches in use
      integer(ik) :: cap = 0_ik       !< allocated capacity
      !==========================================================================================!
      ! 1. IDENTITY AND AREA                                                                     !
      !==========================================================================================!
      integer(ik), allocatable :: global_id(:)      !< persistent id, never reused
      integer(ik), allocatable :: dist_type(:)      !< disturbance/land-use class
      integer(ik), allocatable :: owner_site(:)     !< its site in the polygon's site block
      real(wp),    allocatable :: area(:)           !< fraction of its site, sum = 1 (kept real64)
      real(wp),    allocatable :: age(:)            !< [yr]
      !==========================================================================================!
      ! 2. COHORT MAP: a CSR map into the cohort block (rebuild_cohort_map derives it)           !
      !==========================================================================================!
      integer(ik), allocatable :: cohort_offset(:)  !< first cohort of the patch (1-based)
      integer(ik), allocatable :: cohort_count(:)   !< cohorts in the patch
      !==========================================================================================!
      ! 3. RECRUIT BANK                                                                          !
      !==========================================================================================!
      real(wp),    allocatable :: recruit_pool(:,:)  !< [plant/m2] (pft, patch) recruits waiting to reach min_recruit_size
      !==========================================================================================!
      ! 4. CANOPY AIR                                                                            !
      !==========================================================================================!
      type(cas_state_t),          allocatable :: cas(:)          !< canopy air
      !==========================================================================================!
      ! 5. SNOW                                                                                  !
      !==========================================================================================!
      type(snow_state_t),         allocatable :: snow(:)         !< snow, or temporary surface water
      !==========================================================================================!
      ! 6. SOIL WATER AND HEAT                                                                   !
      !==========================================================================================!
      type(soil_water_state_t),   allocatable :: soil_water(:)   !< soil water and the pond
      type(soil_energy_state_t),  allocatable :: soil_energy(:)  !< soil heat
      !==========================================================================================!
      ! 7. SOIL CARBON                                                                           !
      !==========================================================================================!
      type(soil_carbon_state_t),  allocatable :: soil_carbon(:)  !< soil-carbon pools
      !==========================================================================================!
      ! 8. DAILY SLOW-TO-FAST RATES (set once per slow step, held fixed all day)                 !
      !==========================================================================================!
      real(wp),    allocatable :: shed_water_rate(:)   !< [kg/m2/s] water of shed tissue, to the ground (soil_plant_water.md)
      real(wp),    allocatable :: shed_water_enth(:)   !< [W/m2] its enthalpy, at the temperatures of the tissues it left
      real(wp),    allocatable :: slow_co2_rate(:)     !< [umol/m2/s] growth respiration less starvation deficit, to the canopy air
      !==========================================================================================!
      ! 9. WITHIN-STEP ACCUMULATORS (reset each step, never restarted)                           !
      !==========================================================================================!
      !----- On the patch, not local to the slow driver, so a patch made or fused partway        !
      !      through the step keeps its own sums. -----------------------------------------------!
      type(decomp_xi_t),          allocatable :: decomp_xi(:)  !< the day's decomposition scalars
      type(litter_input_t),       allocatable :: litter_in(:)  !< [kgC/m2/day] this step's litter
      !==========================================================================================!
      ! 10. SOLVER MEMORY                                                                        !
      !==========================================================================================!
      !----- Per patch, so a patch's first sub-step does not depend on which patch ran before    !
      !      it, or on the thread schedule (docs/science/numerical_scheme.md §3). ---------------!
      real(wp),    allocatable :: adapt_dt_last(:)     !< [s] last accepted fast-loop sub-step
      !----- The patch's wall time in the last fast window, which orders the next window's patch   !
      !      loop costliest first. Scheduling only: no result reads it, and it is never restarted. --!
      real(wp),    allocatable :: fast_cost(:)         !< [s]
      !==========================================================================================!
      ! 11. OUTPUT ACCUMULATORS (ed_diag_vars; never restarted)                                  !
      !==========================================================================================!
      type(patch_diag_block) :: diag   !< the fast loop's rows, then the slow step's
   end type patch_block

contains

   !=======================================================================================!
   !  BLENDS: the store of a patch made from two, in a fusion or a disturbance gap. w1     !
   !  and w2 are area fractions, w1 = a1/(a1+a2), so a quantity per ground area keeps its  !
   !  site total; the canopy air is weighted by air mass instead. Diagnosed fields are     !
   !  mixed only for a provisional value: the caller re-diagnoses them                     !
   !  (docs/science/coupled_biophysics.md, "Mixing two patches").                          !
   !=======================================================================================!

   !----- Mix two canopy airs on their air mass, area x depth: enthalpy, humidity and CO2  !
   !      are per kg of air. The depth is area-weighted, the merged volume over the merged !
   !      ground, and that is also the total mass weight. ---------------------------------!
   pure function blend_cas(w1, a, w2, b) result(c)
      real(wp),          intent(in) :: w1, w2
      type(cas_state_t), intent(in) :: a, b
      type(cas_state_t)             :: c
      real(wp) :: m1, m2, mt
      m1 = w1 * a%can_depth
      m2 = w2 * b%can_depth
      mt = m1 + m2
      c%can_depth = mt                                   ! = w1*d1 + w2*d2, the area-weighted depth
      if (mt > tiny_num) then        ! below it, a blend has no air mass to conserve
         c%can_enthalpy = (m1 * a%can_enthalpy + m2 * b%can_enthalpy) / mt
         c%can_shv      = (m1 * a%can_shv      + m2 * b%can_shv     ) / mt
         c%can_co2      = (m1 * a%can_co2      + m2 * b%can_co2     ) / mt
         c%can_temp     = (m1 * a%can_temp     + m2 * b%can_temp    ) / mt
      else
         ! both boxes are empty: nothing to conserve, so area weights
         c%can_enthalpy = w1 * a%can_enthalpy + w2 * b%can_enthalpy
         c%can_shv      = w1 * a%can_shv      + w2 * b%can_shv
         c%can_co2      = w1 * a%can_co2      + w2 * b%can_co2
         c%can_temp     = w1 * a%can_temp     + w2 * b%can_temp
      end if
   end function blend_cas

   !----- Resize the canopy-air box when the canopy height changes. The box is open: a     !
   !      deeper canopy entrains air from above at the canopy-air state, so the intensive  !
   !      state stays and the content changes. The optional *_open outputs, which need     !
   !      rho_air, report that exchange for a caller to book                               !
   !      (docs/science/canopy_air_space_biophysics.md, the canopy air's three variables). !
   pure subroutine cas_set_depth(cas, depth_new, rho_air, de_open, dw_open, dc_open)
      type(cas_state_t),  intent(inout) :: cas
      real(wp),           intent(in)    :: depth_new  !< [m] new CAS depth
      real(wp), optional, intent(in)    :: rho_air    !< [kg/m3]
      real(wp), optional, intent(out)   :: de_open    !< [J/m2]     energy    entrained (+) / detrained (-)
      real(wp), optional, intent(out)   :: dw_open    !< [kg/m2]    vapour    entrained (+) / detrained (-)
      real(wp), optional, intent(out)   :: dc_open    !< [umol/m2]  CO2       entrained (+) / detrained (-)
      real(wp) :: dmass, dmol
      dmass = 0.0_wp
      dmol = 0.0_wp
      if (present(rho_air)) then
         dmass = rho_air * (depth_new - cas%can_depth)
         dmol  = cas_molar_density(rho_air, cas%can_shv) * (depth_new - cas%can_depth)
      end if
      if (present(de_open)) de_open = dmass * cas%can_enthalpy
      if (present(dw_open)) dw_open = dmass * cas%can_shv
      if (present(dc_open)) dc_open = dmol  * cas%can_co2
      cas%can_depth = depth_new
   end subroutine cas_set_depth

   !----- The canopy air's capacities for its three variables: its air mass [kg/m2] and its dry-air  !
   !      moles [mol/m2], the depth times the density and the dry-air molar density at its humidity  !
   !      (docs/science/canopy_air_space_biophysics.md eq. (1)). The fast step and the slow ledger   !
   !      both take them from here. -------------------------------------------------------------!
   elemental subroutine cas_capacities(rho_air, shv, depth, molar_density, mass_capacity, molar_capacity)
      real(wp), intent(in)  :: rho_air          !< [kg/m3]  air density
      real(wp), intent(in)  :: shv              !< [kg/kg]  specific humidity
      real(wp), intent(in)  :: depth            !< [m]      canopy-air depth
      real(wp), intent(out) :: molar_density    !< [mol/m3] dry-air molar density
      real(wp), intent(out) :: mass_capacity    !< [kg/m2]
      real(wp), intent(out) :: molar_capacity   !< [mol/m2]
      molar_density  = cas_molar_density(rho_air, shv)
      mass_capacity  = rho_air * depth
      molar_capacity = molar_density * depth
   end subroutine cas_capacities

   !----- The soil column's water [kg/m2] and internal energy [J/m2], from theta [m3/m3] and the energy  !
   !      per volume [J/m3] over the layer thicknesses: plain ascending loops, which the fast step's and  !
   !      the slow step's ledgers both take from here. --------------------------------------------!
   pure function soil_water_store(theta, dz, nsl) result(w)
      real(wp),    intent(in) :: theta(:), dz(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: w
      integer(ik) :: k
      w = 0.0_wp
      do k = 1_ik, nsl
         w = w + theta(k) * dz(k) * rho_h2o
      end do
   end function soil_water_store

   pure function soil_energy_store(soil_energy, dz, nsl) result(e)
      real(wp),    intent(in) :: soil_energy(:), dz(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: e
      integer(ik) :: k
      e = 0.0_wp
      do k = 1_ik, nsl
         e = e + soil_energy(k) * dz(k)
      end do
   end function soil_energy_store

   !----- nlayer is the larger; the caller sets it to 0 when no snow is left. -------------!
   pure function blend_snow(w1, a, w2, b) result(c)
      real(wp),            intent(in) :: w1, w2
      type(snow_state_t),  intent(in) :: a, b
      type(snow_state_t)              :: c
      c%swe         = w1 * a%swe         + w2 * b%swe
      c%snow_energy = w1 * a%snow_energy + w2 * b%snow_energy
      c%snow_depth  = w1 * a%snow_depth  + w2 * b%snow_depth
      c%snow_temp   = w1 * a%snow_temp   + w2 * b%snow_temp   ! provisional
      c%snow_fliq   = w1 * a%snow_fliq   + w2 * b%snow_fliq   ! provisional
      c%nlayer      = max(a%nlayer, b%nlayer)
   end function blend_snow

   pure function blend_soil_water(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_water_state_t),  intent(in) :: a, b
      type(soil_water_state_t)              :: c
      c%theta     = w1 * a%theta     + w2 * b%theta
      c%w_surface = w1 * a%w_surface + w2 * b%w_surface
      c%w_surface_energy = w1 * a%w_surface_energy + w2 * b%w_surface_energy
   end function blend_soil_water

   pure function blend_soil_energy(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_energy_state_t),  intent(in) :: a, b
      type(soil_energy_state_t)              :: c
      c%soil_energy = w1 * a%soil_energy + w2 * b%soil_energy
      c%soil_temp   = w1 * a%soil_temp   + w2 * b%soil_temp
      c%soil_fliq   = w1 * a%soil_fliq   + w2 * b%soil_fliq
   end function blend_soil_energy

   pure function blend_soil_carbon(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_carbon_state_t), intent(in) :: a, b
      type(soil_carbon_state_t)             :: c
      c%carbon = w1 * a%carbon + w2 * b%carbon
      c%lignin = w1 * a%lignin + w2 * b%lignin
   end function blend_soil_carbon

   pure function blend_decomp_xi(w1, a, w2, b) result(c)
      real(wp),          intent(in) :: w1, w2
      type(decomp_xi_t), intent(in) :: a, b
      type(decomp_xi_t)             :: c
      c%xi_days = w1 * a%xi_days + w2 * b%xi_days
      c%rh_fast = w1 * a%rh_fast + w2 * b%rh_fast
   end function blend_decomp_xi

   pure function blend_litter_input(w1, a, w2, b) result(c)
      real(wp),             intent(in) :: w1, w2
      type(litter_input_t), intent(in) :: a, b
      type(litter_input_t)             :: c
      c%labile_grnd   = w1 * a%labile_grnd   + w2 * b%labile_grnd
      c%labile_soil   = w1 * a%labile_soil   + w2 * b%labile_soil
      c%struct_grnd   = w1 * a%struct_grnd   + w2 * b%struct_grnd
      c%struct_soil   = w1 * a%struct_soil   + w2 * b%struct_soil
      c%lignin_grnd   = w1 * a%lignin_grnd   + w2 * b%lignin_grnd
      c%lignin_soil   = w1 * a%lignin_soil   + w2 * b%lignin_soil
      c%leaf_c        = w1 * a%leaf_c        + w2 * b%leaf_c
      c%fineroot_c    = w1 * a%fineroot_c    + w2 * b%fineroot_c
      c%wood_c        = w1 * a%wood_c        + w2 * b%wood_c
      c%storage_c     = w1 * a%storage_c     + w2 * b%storage_c
      c%repro_c       = w1 * a%repro_c       + w2 * b%repro_c
      c%exudate_c     = w1 * a%exudate_c     + w2 * b%exudate_c
   end function blend_litter_input


   !=======================================================================================!
   !  THE PATCH LOCKSTEP: the routines that list every patch field, in section order       !
   !=======================================================================================!

   subroutine alloc_patch_block(patch, cap, n_pft)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: cap, n_pft
      patch%cap = cap
      patch%n = 0_ik
      !----- 1. Identity and area ---------------------------------------------------------------!
      allocate(patch%global_id(cap), patch%dist_type(cap), patch%owner_site(cap), patch%area(cap),         &
               patch%age(cap))
      patch%global_id = 0_ik
      patch%dist_type = DIST_PRIMARY
      patch%owner_site = 1_ik                                    ! the first site, until placed in another
      patch%area = 0.0_wp
      patch%age = 0.0_wp
      !----- 2. Cohort map ----------------------------------------------------------------------!
      allocate(patch%cohort_offset(cap), patch%cohort_count(cap))
      patch%cohort_offset = 0_ik
      patch%cohort_count = 0_ik
      !----- 3. Recruit bank --------------------------------------------------------------------!
      allocate(patch%recruit_pool(n_pft, cap))
      patch%recruit_pool = 0.0_wp
      !----- 4. Canopy air ----------------------------------------------------------------------!
      allocate(patch%cas(cap))                                   ! default-initialised by its type
      !----- 5. Snow ----------------------------------------------------------------------------!
      allocate(patch%snow(cap))                                  ! default-initialised by its type
      !----- 6. Soil water and heat -------------------------------------------------------------!
      allocate(patch%soil_water(cap), patch%soil_energy(cap))    ! default-initialised by its type
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      allocate(patch%soil_carbon(cap))                           ! default-initialised by its type
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      allocate(patch%shed_water_rate(cap), patch%shed_water_enth(cap), patch%slow_co2_rate(cap))
      patch%shed_water_rate = 0.0_wp
      patch%shed_water_enth = 0.0_wp
      patch%slow_co2_rate = 0.0_wp
      !----- 9. Within-step accumulators --------------------------------------------------------!
      allocate(patch%decomp_xi(cap), patch%litter_in(cap))       ! default-initialised by its type
      !----- 10. Solver memory ------------------------------------------------------------------!
      allocate(patch%adapt_dt_last(cap), patch%fast_cost(cap))
      patch%adapt_dt_last = 0.0_wp
      patch%fast_cost     = 0.0_wp
   end subroutine alloc_patch_block

   !----- Make room for at least `need` patches, growing 1.5x; a no-op when there is room  !
   !      (as expand_cohort_block). A new patch field is listed here, in alloc_patch_block,!
   !      reorder_patch_block and free_patch_block. ---------------------------------------!
   subroutine expand_patch_block(patch, need, n_pft)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: need, n_pft
      type(patch_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= patch%cap) return
      newcap = max(need, (patch%cap * 3_ik) / 2_ik + 1_ik)
      !----- Grow the diagnostic block in place; the move_alloc list below leaves it alone.  -----!
      call expand_patch_diag(patch%diag, newcap)
      call alloc_patch_block(tmp, newcap, n_pft)
      m = patch%n
      tmp%n = m
      !----- 1. Identity and area ---------------------------------------------------------------!
      tmp%global_id(1:m) = patch%global_id(1:m)
      tmp%dist_type(1:m)  = patch%dist_type(1:m)
      tmp%owner_site(1:m) = patch%owner_site(1:m)
      tmp%area(1:m)       = patch%area(1:m)
      tmp%age(1:m)       = patch%age(1:m)
      !----- 2. Cohort map ----------------------------------------------------------------------!
      tmp%cohort_offset(1:m) = patch%cohort_offset(1:m)
      tmp%cohort_count(1:m)  = patch%cohort_count(1:m)
      !----- 3. Recruit bank --------------------------------------------------------------------!
      tmp%recruit_pool(:,1:m) = patch%recruit_pool(:,1:m)
      !----- 4. Canopy air ----------------------------------------------------------------------!
      tmp%cas(1:m) = patch%cas(1:m)
      !----- 5. Snow ----------------------------------------------------------------------------!
      tmp%snow(1:m) = patch%snow(1:m)
      !----- 6. Soil water and heat -------------------------------------------------------------!
      tmp%soil_water(1:m)  = patch%soil_water(1:m)
      tmp%soil_energy(1:m) = patch%soil_energy(1:m)
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      tmp%soil_carbon(1:m) = patch%soil_carbon(1:m)
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      tmp%shed_water_rate(1:m) = patch%shed_water_rate(1:m)
      tmp%shed_water_enth(1:m) = patch%shed_water_enth(1:m)
      tmp%slow_co2_rate(1:m)   = patch%slow_co2_rate(1:m)
      !----- 9. Within-step accumulators --------------------------------------------------------!
      tmp%decomp_xi(1:m) = patch%decomp_xi(1:m)
      tmp%litter_in(1:m) = patch%litter_in(1:m)
      !----- 10. Solver memory ------------------------------------------------------------------!
      tmp%adapt_dt_last(1:m) = patch%adapt_dt_last(1:m)
      tmp%fast_cost(1:m)     = patch%fast_cost(1:m)
      patch%n = tmp%n
      patch%cap = tmp%cap
      !----- 1. Identity and area ---------------------------------------------------------------!
      call move_alloc(tmp%global_id, patch%global_id)
      call move_alloc(tmp%dist_type, patch%dist_type)
      call move_alloc(tmp%owner_site, patch%owner_site)
      call move_alloc(tmp%area, patch%area)
      call move_alloc(tmp%age, patch%age)
      !----- 2. Cohort map ----------------------------------------------------------------------!
      call move_alloc(tmp%cohort_offset, patch%cohort_offset)
      call move_alloc(tmp%cohort_count, patch%cohort_count)
      !----- 3. Recruit bank --------------------------------------------------------------------!
      call move_alloc(tmp%recruit_pool, patch%recruit_pool)
      !----- 4. Canopy air ----------------------------------------------------------------------!
      call move_alloc(tmp%cas, patch%cas)
      !----- 5. Snow ----------------------------------------------------------------------------!
      call move_alloc(tmp%snow, patch%snow)
      !----- 6. Soil water and heat -------------------------------------------------------------!
      call move_alloc(tmp%soil_water, patch%soil_water)
      call move_alloc(tmp%soil_energy, patch%soil_energy)
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      call move_alloc(tmp%soil_carbon, patch%soil_carbon)
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      call move_alloc(tmp%shed_water_rate, patch%shed_water_rate)
      call move_alloc(tmp%shed_water_enth, patch%shed_water_enth)
      call move_alloc(tmp%slow_co2_rate, patch%slow_co2_rate)
      !----- 9. Within-step accumulators --------------------------------------------------------!
      call move_alloc(tmp%decomp_xi, patch%decomp_xi)
      call move_alloc(tmp%litter_in, patch%litter_in)
      !----- 10. Solver memory ------------------------------------------------------------------!
      call move_alloc(tmp%adapt_dt_last, patch%adapt_dt_last)
      call move_alloc(tmp%fast_cost, patch%fast_cost)
   end subroutine expand_patch_block

   !----- Put old patch perm(k) at position k, k = 1..m, in every per-patch array and the  !
   !      diagnostics (as reorder_cohort_block). A vector subscript on the right is        !
   !      evaluated first, so the permutation is safe in place. ---------------------------!
   subroutine reorder_patch_block(patch, perm, m)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: perm(:)
      integer(ik),       intent(in)    :: m
      !----- 1. Identity and area ---------------------------------------------------------------!
      patch%global_id(1:m) = patch%global_id(perm(1:m))
      patch%dist_type(1:m)  = patch%dist_type(perm(1:m))
      patch%owner_site(1:m) = patch%owner_site(perm(1:m))
      patch%area(1:m)       = patch%area(perm(1:m))
      patch%age(1:m)       = patch%age(perm(1:m))
      !----- 2. Cohort map ----------------------------------------------------------------------!
      !      (not permuted: rebuild_cohort_map derives the map from the cohorts afterwards)
      !----- 3. Recruit bank --------------------------------------------------------------------!
      patch%recruit_pool(:,1:m) = patch%recruit_pool(:,perm(1:m))
      !----- 4. Canopy air ----------------------------------------------------------------------!
      patch%cas(1:m) = patch%cas(perm(1:m))
      !----- 5. Snow ----------------------------------------------------------------------------!
      patch%snow(1:m) = patch%snow(perm(1:m))
      !----- 6. Soil water and heat -------------------------------------------------------------!
      patch%soil_water(1:m)  = patch%soil_water(perm(1:m))
      patch%soil_energy(1:m) = patch%soil_energy(perm(1:m))
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      patch%soil_carbon(1:m) = patch%soil_carbon(perm(1:m))
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      patch%shed_water_rate(1:m) = patch%shed_water_rate(perm(1:m))
      patch%shed_water_enth(1:m) = patch%shed_water_enth(perm(1:m))
      patch%slow_co2_rate(1:m)   = patch%slow_co2_rate(perm(1:m))
      !----- 9. Within-step accumulators --------------------------------------------------------!
      patch%decomp_xi(1:m) = patch%decomp_xi(perm(1:m))
      patch%litter_in(1:m) = patch%litter_in(perm(1:m))
      !----- 10. Solver memory ------------------------------------------------------------------!
      patch%adapt_dt_last(1:m) = patch%adapt_dt_last(perm(1:m))
      patch%fast_cost(1:m)     = patch%fast_cost(perm(1:m))
      !----- 11. Output accumulators ------------------------------------------------------------!
      call reorder_patch_diag(patch%diag, perm, m)
      patch%n = m
   end subroutine reorder_patch_block

   !----- Free every per-patch array. -----------------------------------------------------!
   subroutine free_patch_block(patch)
      type(patch_block), intent(inout) :: patch
      patch%n = 0_ik
      patch%cap = 0_ik
      if (allocated(patch%area)) deallocate(patch%global_id, patch%dist_type, patch%owner_site, patch%area,    &
                                            patch%age,                                                         &
                                            patch%cohort_offset, patch%cohort_count, patch%recruit_pool,       &
                                            patch%cas, patch%snow, patch%soil_water, patch%soil_energy,        &
                                            patch%soil_carbon, patch%shed_water_rate, patch%shed_water_enth,   &
                                            patch%slow_co2_rate,                                               &
                                            patch%decomp_xi, patch%litter_in, patch%adapt_dt_last,      &
                                            patch%fast_cost)
   end subroutine free_patch_block

   !=======================================================================================!
   !  THE PHYSICAL STATE A NEW COHORT IS BORN INTO                                          !
   !=======================================================================================!

   !----- A patch's canopy-air temperature [K] and the water potential of its top soil     !
   !      layer [MPa] (top_soil_psi), which a new cohort takes so it starts in equilibrium !
   !      with where it stands. Every patch is seeded from its start (init_bare_sites, a   !
   !      restart's file, a gap's donors). `surf` is the patch's site's surface. ----------!
   pure subroutine patch_physical_state(cas, soil_water, surf, temp, psi)
      type(cas_state_t),         intent(in)  :: cas          !< the patch's canopy air
      type(soil_water_state_t), intent(in)  :: soil_water   !< the patch's soil water
      type(site_surface_t),      intent(in)  :: surf
      real(wp),                  intent(out) :: temp, psi
      temp = cas%can_temp
      psi  = top_soil_psi(surf%soil_profile, soil_water%theta(1))
   end subroutine patch_physical_state

   !=======================================================================================!
   !  BOUNDS: the check the driver runs every step                                         !
   !=======================================================================================!

   !----- The first of a canopy air's prognostic numbers -- enthalpy, humidity, CO2 -- that is   !
   !      not finite, 0 if none: what a fast-loop integrator that could not take even its smallest !
   !      step leaves behind. The driver checks it every step. -------------------------------------!
   pure integer(ik) function canopy_air_bad_field(cas) result(k)
      type(cas_state_t), intent(in) :: cas
      real(wp)    :: x(3)
      integer(ik) :: j
      x = [cas%can_enthalpy, cas%can_shv, cas%can_co2]
      k = 0_ik
      do j = 1_ik, size(x, kind=ik)
         if (.not. ieee_is_finite(x(j))) then
            k = j
            return
         end if
      end do
   end function canopy_air_bad_field

   !----- Name for the index canopy_air_bad_field returns. ----------------------------------------!
   pure function canopy_air_field_name(k) result(nm)
      integer(ik), intent(in) :: k
      character(len=12)       :: nm
      character(len=12), parameter :: NAMES(3) = [character(len=12) :: 'can_enthalpy', 'can_shv', 'can_co2']
      nm = '(none)'
      if (k >= 1_ik .and. k <= size(NAMES, kind=ik)) nm = NAMES(k)
   end function canopy_air_field_name

   !=======================================================================================!
   !  One patch's stores from two others', and its shed channel.                           !
   !=======================================================================================!

   !----- Patch dst's stores as w1 x src1's plus w2 x src2's, in the file's section order: the  !
   !      recruit pool, the canopy air, snow, soil water and heat, soil carbon, the day's      !
   !      decomposition integral and litter. Patch fusion and a new treefall gap both blend     !
   !      through it, so a store added to the patch is blended in both. ----------------------!
   subroutine blend_patch_stores(patch, dst, w1, src1, w2, src2)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: dst, src1, src2
      real(wp),          intent(in)    :: w1, w2
      patch%recruit_pool(:, dst) = w1 * patch%recruit_pool(:, src1) + w2 * patch%recruit_pool(:, src2)
      patch%cas(dst)         = blend_cas(w1, patch%cas(src1), w2, patch%cas(src2))
      patch%snow(dst)        = blend_snow(w1, patch%snow(src1), w2, patch%snow(src2))
      patch%soil_water(dst)  = blend_soil_water(w1, patch%soil_water(src1), w2, patch%soil_water(src2))
      patch%soil_energy(dst) = blend_soil_energy(w1, patch%soil_energy(src1), w2, patch%soil_energy(src2))
      patch%soil_carbon(dst) = blend_soil_carbon(w1, patch%soil_carbon(src1), w2, patch%soil_carbon(src2))
      patch%decomp_xi(dst)   = blend_decomp_xi(w1, patch%decomp_xi(src1), w2, patch%decomp_xi(src2))
      patch%litter_in(dst)   = blend_litter_input(w1, patch%litter_in(src1), w2, patch%litter_in(src2))
   end subroutine blend_patch_stores

   !----- Add water [kg/m2 of patch] and its enthalpy [J/m2] to patch ip's shed channel, as rates  !
   !      over the slow step; the next fast window lands them on the ground (conservation.md §3). !
   pure subroutine add_shed_water(patch, ip, water, enthalpy, dt_slow)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: ip
      real(wp),          intent(in)    :: water, enthalpy, dt_slow
      patch%shed_water_rate(ip) = patch%shed_water_rate(ip) + water / dt_slow
      patch%shed_water_enth(ip) = patch%shed_water_enth(ip) + enthalpy / dt_slow
   end subroutine add_shed_water

end module ed_patch_state_vars
