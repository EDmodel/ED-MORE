! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_types -- the fast loop's records, and the march state's arithmetic. What
! fast_step takes besides the configuration and the site: the static state (the patch's
! cohorts and what they fix for the slow step), the dynamic state it advances, the
! forcing, the output (the step report, the aerodynamics, the per-cohort fluxes) and the
! workspace. What the integrator works with inside: the error control's tolerances, the
! frozen record (frozen_t), the march state and its tendency (march_state_t,
! march_tend_t), the boundary-flux rates and amounts (stage_bflux_t, boundary_flux_t), and
! the soil-plant water network's records (water_drive_t, water_tend_t, water_state_t,
! water_flows_t, water_work_t). The march state's fields are listed once, in
! state_to_array, array_to_state, tend_to_array and state_entry_rules; the combinators
! work on the flat array and never name a field. A new state field goes in those four, the
! type, pack_march_start, unpack_column_state and test_fast_integ_kernels. The module sits
! below every other integrator module (docs/science/numerical_scheme.md §2-§3).
!---~---
module ed_fast_integ_types
   use ed_kinds, only : wp, ik
   use ed_constants, only : n_soil_layer_max
   use ed_pft_params, only : leaf_photo_table_t
   use ed_patch_state_vars,   only : cas_state_t, soil_water_state_t, soil_energy_state_t, snow_state_t, soil_carbon_state_t
   use ed_therm_lib, only : cas_temp_to_enthalpy, cas_enthalpy_to_temp
   use ed_biophysics_params,   only : aero_params_t
   use ed_biogeochem_params,   only : N_SOIL_POOL
   use ed_budget_check, only : budget_t
   use ed_linear_algebra, only : arrowhead_system_t
   use ed_soil_water, only : ground_evap_env_t
   use ed_soil_hydr_lib, only : soil_ice_t
   use ed_canopy_aerodynamics, only : aero_env_t, aero_args_t, aero_out_t, ensure_aero_out_capacity
   use ed_canopy_air_space, only : cas_env_t
   use ed_plant_energy,   only : plant_energy_env_t
   use ed_canopy_water,   only : canopy_water_out_t
   use ed_surface, only : surface_env_t, snow_out_t
   use ed_soil_heat, only : soil_layer_temp
   use ed_soil_params, only : soil_profile_t
   use ed_diag_vars, only : N_CO_FAST_ROWS
   use ed_state_defaults,    only : UNSET_NEGATIVE
   implicit none
   private

   public :: fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t, fast_workspace_t
   public :: step_report_t, step_ledgers_t, step_exports_t
   public :: WORK_INTEG_STEPS, WORK_INTEG_REJECT, WORK_WATER_NSUB, WORK_WATER_NEWTON, WORK_WATER_REJECT,         &
             WORK_WATER_NFAIL, WORK_RK45_RESCUE, WORK_CLAMP_STAGE, N_WORK
   public :: GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W, GRP_THETA, N_TOL_GROUP
   public :: tol_set_t, error_control_t
   public :: ensure_static_capacity, ensure_dynamic_capacity, ensure_forcing_capacity, ensure_output_capacity
   public :: surface_state_t, surface_tend_t
   public :: cas_boundary_t
   public :: soil_hydrology_t, root_zone_t, plant_water_t, esdirk2_stages_t
   public :: march_state_t, frozen_t, march_tend_t
   public :: stage_bflux_t, boundary_flux_t
   public :: POND, WOOD, LEAF, COHORT_STORES
   public :: water_drive_t, water_tend_t, water_state_t, water_flows_t, water_work_t, ensure_water_capacity
   public :: alloc_static_state, alloc_dynamic_state   ! used only by the tests
   !--- The march state's arithmetic.
   public :: state_to_array, state_entry_rules
   public :: state_init, state_accum, state_extrap, state_sub, state_err_diff, zero_like
   public :: bflux_zero, bflux_add, bflux_bweight
   public :: clamp_cas, diagnose_soil_temps, at_temp_bound, STATE_TEMP_MIN, STATE_TEMP_MAX

   !---~---
   !   The temperature bounds of a march state [K]: clamp_cas holds a stage's canopy air inside
   ! them, and a march that ends at one has railed (at_temp_bound).
   !---~---
   real(wp), parameter :: STATE_TEMP_MIN = 180.0_wp
   real(wp), parameter :: STATE_TEMP_MAX = 350.0_wp

   !---~---
   !   The tolerance groups, one per kind of state, each scaling its entries in the error norm
   ! (docs/science/numerical_scheme.md eq. (3)). GRP_THETA also sets the soil-plant water
   ! solver's own sub-step control, so one tolerance set governs both.
   !---~---
   integer(ik), parameter :: GRP_ENTH    = 1_ik   ! CAS specific enthalpy   [J/kg]
   integer(ik), parameter :: GRP_SHV     = 2_ik   ! CAS specific humidity   [kg/kg]
   integer(ik), parameter :: GRP_CO2     = 3_ik   ! CAS CO2 mole fraction   [umol/mol]
   integer(ik), parameter :: GRP_SE      = 4_ik   ! soil internal energy    [J/m3]
   integer(ik), parameter :: GRP_LEAF_W  = 5_ik   ! leaf internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_WOOD_W  = 6_ik   ! wood internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_THETA   = 7_ik   ! soil moisture           [m3/m3] (soil-water sub-solver)
   !--- No soil-temperature group: the state carries soil energy (GRP_SE).
   integer(ik), parameter :: N_TOL_GROUP = 7_ik

   !--- Per-field default absolute tolerances: the group defaults unless config overrides them.
   real(wp), parameter :: ATOL_ENTH_DEF   = 5.0e1_wp    ! [J/kg]      (~0.05 K in enthalpy)
   real(wp), parameter :: ATOL_SHV_DEF    = 1.0e-6_wp   ! [kg/kg]
   real(wp), parameter :: ATOL_CO2_DEF    = 1.0e-1_wp   ! [umol/mol]
   real(wp), parameter :: ATOL_SE_DEF     = 1.0e3_wp    ! [J/m3]
   real(wp), parameter :: ATOL_LEAF_W_DEF = 1.0e-4_wp   ! [kg/plant]
   real(wp), parameter :: ATOL_WOOD_W_DEF = 1.0e-4_wp   ! [kg/plant]
   real(wp), parameter :: ATOL_THETA_DEF  = 1.0e-4_wp   ! [m3/m3] (== [fast].water_atol's default)

   !---~---
   !   The soil-plant water network's nodes (docs/science/soil_plant_water.md §1): the pond is
   ! soil node POND and soil layer k is node 1 + k; each cohort has two stores, its wood and its
   ! leaves.
   !---~---
   integer(ik), parameter :: POND = 1_ik
   integer(ik), parameter :: WOOD = 1_ik
   integer(ik), parameter :: LEAF = 2_ik
   integer(ik), parameter :: COHORT_STORES = 2_ik

   !--- Per-group (rtol, atol). The WRMS normalizes state group g by atol(g) + rtol(g)*|y|.
   type :: tol_set_t
      real(wp) :: rtol(N_TOL_GROUP) = 1.0e-3_wp
      real(wp) :: atol(N_TOL_GROUP) = [ATOL_ENTH_DEF, ATOL_SHV_DEF, ATOL_CO2_DEF, ATOL_SE_DEF,   &
                                       ATOL_LEAF_W_DEF, ATOL_WOOD_W_DEF, ATOL_THETA_DEF]
   end type tol_set_t

   !--- What an adaptive march runs under: the step-factor limits and the tolerances.
   type :: error_control_t
      real(wp)         :: safety     = 0.9_wp
      real(wp)         :: fmin       = 0.2_wp
      real(wp)         :: fmax       = 5.0_wp
      !---~---
      !   p of docs/science/numerical_scheme.md eq. (4), the embedded estimate's order: 1 ESDIRK2, 4
      ! RK45.
      !---~---
      integer(ik)      :: p_order    = 1_ik
      type(tol_set_t)  :: tols
      !---~---
      !   [fast].debug_error's budget stops: a breached ledger stops the run, and so does a soil
      ! heat step whose budget misses by more than soil_energy_budget_atol [K] of the column's
      ! heat capacity (atol_scale applied).
      !---~---
      logical          :: debug_error = .false.
      real(wp)         :: soil_energy_budget_atol = 1.0e-2_wp
   end type error_control_t

   !---~---
   !   The step's state of the patch's cohorts and everything else a fast step reads that is fixed
   ! for the slow step, gathered once per patch (gather_patch): each cohort's structure and slow-
   ! varying traits, the canopy geometry and the canopy-air top, the slow loop's daily inputs,
   ! and the site's leaf table.
   !---~---
   type :: fast_static_state_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       ! PFT index (into the PFT table)
      !---~---
      !   The cohorts from the ground up: order(k) is the k-th lowest (ascending_order of height). The
      ! cohorts are gathered tallest first; the radiation and the aerodynamics take them bottom to
      ! top, so the order is formed once, with the heights, and both read it.
      !---~---
      integer(ik), allocatable :: order(:)
      real(wp),    allocatable :: lai(:), wai(:), height(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: aboveground_frac(:)          ! [--] gathered per-PFT (stem respiration)
      logical,     allocatable :: is_woody(:)                   ! gathered per-PFT (stem respiration off for grass)
      real(wp),    allocatable :: stem_resp_factor25(:)         ! [umol CO2/m2 stem/s @25C] gathered per-PFT
      real(wp),    allocatable :: stem_resp_sapwood25(:)        ! [umol CO2/m3 sapwood/s @25C] gathered per-PFT
      real(wp),    allocatable :: sap_volume(:)                 ! [m3/plant] sapwood volume (stem respiration)
      real(wp),    allocatable :: root_resp_factor25(:)         ! [umol CO2/kgC root/s @25C] gathered per-PFT
      real(wp),    allocatable :: leaf_emissivity(:), wood_emissivity(:)   ! [-] gathered per-PFT (tissue longwave)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   ! [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              ! [kgC/plant],[kgC/plant],[m2] (hydraulics)
      !---~---
      !   Total wood carbon, not interchangeable with bsap: the sapwood sets the hydraulic
      ! capacitance, while all the wood, thin branches included, is the one thermal store
      ! (docs/science/vegetation_energy_dynamics.md §The store is ON).
      !---~---
      real(wp),    allocatable :: bwood(:)                                    ! [kgC/plant] total wood (thermal store)
      real(wp),    allocatable :: vcmax25(:), rd25(:)                         ! [umol/m2/s] per-cohort (plastic) capacities
      !---~---
      !   Yesterday's daily-maximum leaf water potential, the predawn potential behind the slow
      ! stomatal control (docs/science/leaf_gas_exchange.md §4.2). The frozen build always
      ! passes it, since leaf_gas_exchange_batch reads an absent `psi` as no stress.
      !---~---
      real(wp),    allocatable :: dmax_psi_leaf(:)                        ! [MPa] <= 0; UNSET_NEGATIVE => seed from soil
      !---~---
      !   What the cohorts' structure fixes for the slow step, derived once at the gather
      ! (ed_fast_state::derive_cohort_terms): each cohort's fine roots by soil layer and the
      ! patch's, weighted by fine-root carbon; the water network's paths, the intact xylem,
      ! each layer's rhizosphere and the head from the wood to the crown; and the dry heat
      ! capacity of the leaves and of all the wood.
      !---~---
      real(wp),    allocatable :: root_frac(:,:)          ! (n_soil_layer_max, n) [-] each cohort's fine roots
      real(wp)                 :: root_frac_patch(n_soil_layer_max) = 0.0_wp   ! [-] the patch's, summing to 1
      real(wp),    allocatable :: k_xylem(:)              ! [kg/s/MPa per m2 ground] intact xylem conductance
      real(wp),    allocatable :: g_rhizo(:,:)            ! (n_soil_layer_max, n) rhizosphere conductance per
                                                          ! unit soil conductivity [kg/s/MPa per m2, per m/s]
      real(wp),    allocatable :: crown_head(:)           ! [MPa] hydrostatic head from the wood to the leaves
      real(wp),    allocatable :: bwood_water(:)          ! [kgC/plant] sapwood and fine roots, the carbon of the
                                                          ! network's wood store
      real(wp),    allocatable :: leaf_dry_hcap(:), wood_dry_hcap(:)   ! [J/K per m2 ground]
      !---~---
      !   The canopy geometry, and the canopy-air top the forcing is moved to with the patch's
      ! roughness and displacement (docs/science/forcing.md §8).
      !---~---
      type(aero_args_t)        :: geom
      real(wp)                 :: z_top    = 30.0_wp   ! [m] the canopy-air top, the forcing's reference height
      real(wp)                 :: rough    = 0.0_wp    ! [m] roughness length
      real(wp)                 :: displace = 0.0_wp    ! [m] zero-plane displacement height
      !---~---
      !   The patch's soil-carbon pools (docs/science/soil_carbon.md §7), held over the day's fast
      ! sub-steps: the fast Rh respires against them, and the daily soil_carbon_step alone writes
      ! the real pools. 0 when off.
      !---~---
      type(soil_carbon_state_t) :: soil_carbon
      !--- The day's leaf/root-turnover shed water, added to the ground's water input every sub-step.
      real(wp)                 :: shed_water_rate = 0.0_wp  ! [kg/m2 ground/s]
      real(wp)                 :: shed_water_enth = 0.0_wp  ! [W/m2 ground] its enthalpy, from the tissues it left
      !---~---
      !   The slow tier's CO2 owing (growth respiration and its corrections,
      ! docs/science/plant_respiration.md §6).
      !---~---
      real(wp)                 :: slow_co2_rate   = 0.0_wp  ! [umol/m2 ground/s]
      type(leaf_photo_table_t) :: leaf_photo   ! the site's leaf table, acclimated to its growth temperature
   end type fast_static_state_t

   !---~---
   !   The patch's stores and its cohorts' tissue state, the working copy a fast step advances:
   ! gathered from the persisted patch and cohort state and written back after (gather_patch,
   ! scatter_patch).
   !---~---
   type :: fast_dynamic_state_t
      type(cas_state_t)          :: cas               ! canopy-air-space twins (enthalpy/shv/co2)
      type(soil_energy_state_t)  :: soil_energy       ! soil thermal column (internal energy; temp diagnosed)
      type(soil_water_state_t)   :: soil_water        ! soil water column (theta; psi_soil diagnosed)
      type(snow_state_t)         :: snow              ! temporary-surface-water / snow store (swe + energy)
      real(wp), allocatable      :: leaf_energy(:)    ! [J/m2 ground] per-cohort leaf internal energy
      real(wp), allocatable      :: wood_energy(:)    ! [J/m2 ground] per-cohort wood internal energy
      real(wp), allocatable      :: leaf_temp(:)      ! [K] derived from leaf_energy
      real(wp), allocatable      :: wood_temp(:)      ! [K] derived from wood_energy
      !---~---
      !   Internal water mass, the prognostic hydraulic state; psi is diagnosed from it
      ! (tissue_water_to_psi) and never persisted (docs/science/soil_plant_water.md §2).
      !---~---
      real(wp), allocatable      :: leaf_water_mass(:) ! [kg/plant] internal leaf water
      real(wp), allocatable      :: wood_water_mass(:) ! [kg/plant] internal wood water
      !--- The interception films, per m2 of ground here; the cohort holds them per plant.
      real(wp), allocatable      :: leaf_surf_water(:) ! [kg/m2 ground] leaf interception film
      real(wp), allocatable      :: wood_surf_water(:) ! [kg/m2 ground] wood interception film
      !---~---
      !   Warm start for the adaptive march (docs/science/numerical_scheme.md §3, `esdirk2`): the step
      ! that seeds the next dt_fast's first sub-step; 0 = no history, start from dt_fast.
      !---~---
      real(wp)                   :: adapt_dt_last = 0.0_wp   ! [s] last accepted sub-step not cut by the step's end
   end type fast_dynamic_state_t

   !---~---
   !   The drivers of one fast step: the air at the patch's canopy-air top and the radiation the
   ! two-stream gives each cohort and the ground (patch_forcing).
   !---~---
   type :: fast_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   ! [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   ! [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp ! [umol/mol] free-atmosphere CO2
      real(wp)              :: air_temp      = 288.0_wp ! [K]        reference-level air temperature
      real(wp)              :: wind          = 2.0_wp   ! [m/s]      wind at the reference level
      real(wp)              :: press         = 101325.0_wp  ! [Pa]   air pressure
      real(wp)              :: rho_air       = 1.2_wp   ! [kg/m3]    air density
      real(wp)              :: abs_sw_ground = 0.0_wp   ! [W/m2] shortwave the ground absorbs, snow and soil
      real(wp)              :: abs_lw_ground = 0.0_wp   ! [W/m2] net longwave at the ground, snow and soil
      !--- Each ground tile's radiation per unit of its own area (ground_tile_radiation).
      real(wp)              :: rnet_soil     = 0.0_wp   ! [W/m2] the bare soil's net all-wave radiation
      real(wp)              :: sw_snow       = 0.0_wp   ! [W/m2] the shortwave the snow absorbs
      real(wp)              :: lw_down_ground = 0.0_wp  ! [W/m2] the longwave reaching the ground
      real(wp)              :: rainfall      = 0.0_wp   ! [kg/m2/s] reference-level rainfall (intercepted downstream)
      real(wp)              :: snowfall      = 0.0_wp   ! [kg/m2/s] snowfall (drives snow accumulation)
      real(wp)              :: par_per_w     = 2.1_wp   ! [umol photon / (W absorbed)] absorbed->PAR-photon factor
      !---~---
      !   Top-of-canopy radiative fluxes per band, kept from the two-stream solve for comparison with a
      ! radiometer or a satellite product (docs/science/canopy_radiation_transfer.md §Surface radiative
      ! fluxes on the diagnostic path). Upwelling longwave includes the surface's own emission.
      !---~---
      real(wp) :: sw_in_vis = 0.0_wp, sw_in_nir = 0.0_wp   ! [W/m2] incident SW at canopy top, per band
      real(wp) :: sw_up_vis = 0.0_wp, sw_up_nir = 0.0_wp   ! [W/m2] upwelling SW leaving the canopy top
      real(wp) :: lw_up     = 0.0_wp                        ! [W/m2] upwelling LW leaving the canopy top
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     ! [W/m2] absorbed SW (VIS+NIR) / net LW per cohort (leaf energy)
      real(wp), allocatable :: abs_par(:)               ! [W/m2] incident-equivalent PAR (VIS) per cohort; the leaf
                                                        ! model applies leaf_absorptance itself (photosynthesis)
      real(wp)              :: apar          = 0.0_wp   ! [W/m2] PAR (VIS) the leaves absorbed, all cohorts
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) ! [W/m2] absorbed SW / net LW per cohort (wood energy)
      real(wp), allocatable :: par_top(:)               ! [W/m2] PAR (VIS) reaching the cohort's top: the
                                                        ! phenology PAR cue's driver (canopy light gradient)
   end type fast_forcing_t

   !---~---
   !   The work counters of one fast step, the cost axis of a numerics comparison: indices of
   ! step_report_t%work, set by the step and accumulated site-wide by the fast driver.
   !---~---
   integer(ik), parameter :: WORK_INTEG_STEPS  = 1_ik   ! accepted integrator sub-steps
   integer(ik), parameter :: WORK_INTEG_REJECT = 2_ik   ! rejected integrator steps
   integer(ik), parameter :: WORK_WATER_NSUB   = 3_ik   ! the water network's accepted sub-steps
   integer(ik), parameter :: WORK_WATER_NEWTON = 4_ik   ! its Newton iterations
   integer(ik), parameter :: WORK_WATER_REJECT = 5_ik   ! its rejected sub-steps (too large an error, or Newton failed)
   integer(ik), parameter :: WORK_WATER_NFAIL  = 6_ik   ! its sub-steps taken at the shortest length although unconverged
   !---~---
   !   dt_fast steps whose RK45 step bailed or committed a railed state and were redone on
   ! ESDIRK2 (docs/science/numerical_scheme.md §3, `rk45`); a persistently high count flags a
   ! stiff regime.
   !---~---
   integer(ik), parameter :: WORK_RK45_RESCUE  = 7_ik
   !---~---
   !   Stage-clamp activations: how often clamp_cas bounded a throwaway stage input's canopy air
   ! so the right-hand side stays evaluable. A rising count means the trial step is too big;
   ! neither scheme clamps the committed state.
   !---~---
   integer(ik), parameter :: WORK_CLAMP_STAGE  = 8_ik
   integer(ik), parameter :: N_WORK            = 8_ik

   !---~---
   !   The per-patch conservation ledgers. The per-kernel ledgers close by construction;
   ! whole_energy and whole_water set every store against the column's true boundary fluxes,
   ! which is what catches a leak across a seam (docs/science/numerical_scheme.md §5b).
   !---~---
   type :: step_ledgers_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
      !---~---
      !   Per-layer face mass closure, summed over layers: |the water a layer's theta moved - the
      ! water its faces were charged for|. It sees the vertical misplacement the whole-column
      ! ledgers cannot (docs/science/numerical_scheme.md §5b). A pre-formed residual, so only
      ! resid/worst/abs_sum/n_check/n_fail are set; n_fail counts unconverged water solves.
      !---~---
      type(budget_t) :: soil_face_mass
   end type step_ledgers_t

   !--- What the step exported, and the slow tier's inputs it produced.
   type :: step_exports_t
      real(wp)       :: gpp_last = 0.0_wp, nee_last = 0.0_wp   ! [umol/m2/s] this step's patch GPP and NEE
      !--- Net CAS -> atmosphere export over the step, b-weighted like the ledgers: the reported H, LE.
      real(wp)       :: atm_heat_export = 0.0_wp   ! [J/m2]  sensible heat over this dt_fast
      real(wp)       :: atm_vap_export  = 0.0_wp   ! [kg/m2] over this dt_fast
      !---~---
      !   Water leaving the soil column over the step, as the water network moved it: the pond's
      ! spill (surface runoff) and the bottom flux (drainage; negative when an aquifer feeds it).
      !---~---
      real(wp)       :: runoff_surf_export = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: drainage_export    = 0.0_wp   ! [kg/m2] over this dt_fast
      !---~---
      !   The canopy air's vapour sources over the step, as its ledger booked them: the cohorts'
      ! transpiration, the canopy films' evaporation (net of dew onto them), the ground's (bare
      ! soil and snow pack), and the condensate (dew, fog) it lost to soil layer 1. ET
      ! (atm_vap_export) is the first three less the condensate, less the canopy air's own
      ! vapour gain.
      !---~---
      real(wp)       :: transp_export      = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: veg_surface_evap_export   = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: ground_evap_export = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: cond_export        = 0.0_wp   ! [kg/m2] over this dt_fast
      !---~---
      !   The whole-column water budget's other terms over the step (docs/science/diagnostics.md
      ! §4, The budgets): the films' net capture of the rain and snow, the pond's infiltration
      ! into the soil, the pack's meltwater and sublimation, the tissue water the slow step
      ! shed, the films' water above capacity that the update removed, and the change of every
      ! store; and three stores at the step's end.
      !---~---
      real(wp)       :: intercept_export    = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: infiltration_export = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: melt_export         = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: sublimation_export  = 0.0_wp   ! [kg/m2] over this dt_fast (deposition < 0)
      real(wp)       :: shed_water_export   = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: veg_surface_overflow_export  = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: water_store_gain    = 0.0_wp   ! [kg/m2] over this dt_fast
      real(wp)       :: veg_surface_water_store          = 0.0_wp   ! [kg/m2] at the step's end
      real(wp)       :: soil_water_store    = 0.0_wp   ! [kg/m2] at the step's end
      real(wp)       :: plant_water_store   = 0.0_wp   ! [kg/m2] at the step's end
      !---~---
      !   The whole-column energy budget's other terms over the step: the heat through the soil
      ! column's bottom (+ up, into it), the enthalpy the rain and snow carried in and the
      ! enthalpy the outflows carried out (runoff, drainage, the surface water's overflow and the
      ! exported vapour's beyond LE), and the change of the soil's, the canopy's (its air and tissues) and the surface's
      ! (snow pack and pond) stores; and the sensible heat from the canopy and the ground to
      ! the canopy air.
      !---~---
      real(wp)       :: heat_bottom_export  = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: precip_enth_export  = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: outflow_enth_export = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: heat_gain_soil      = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: heat_gain_canopy    = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: heat_gain_surface   = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: h_canopy_export     = 0.0_wp   ! [J/m2] over this dt_fast
      real(wp)       :: h_ground_export     = 0.0_wp   ! [J/m2] over this dt_fast
      !---~---
      !   Soil-carbon Rh diagnostics (docs/science/soil_carbon.md §7), set by the frozen build when
      ! soil carbon is on (else 0). The fast driver accumulates them into the patch's day
      ! integrals, decomp_xi.
      !---~---
      real(wp)       :: xi_step(N_SOIL_POOL) = 0.0_wp   ! [-] this sub-step's per-pool env scalar
      real(wp)       :: rh_matrix_step       = 0.0_wp   ! [kgC/m2/day] this sub-step's matrix Rh
   end type step_exports_t

   !--- The step's report: its ledgers, its exports and its work.
   type :: step_report_t
      type(step_ledgers_t) :: ledgers
      type(step_exports_t) :: exports
      integer(ik)          :: work(N_WORK) = 0_ik   ! WORK_*
      !---~---
      !   How far theta lay outside its constitutive domain at an RK45 stage-1 input, which no
      ! clamp sees; telemetry, not a correction (docs/science/numerical_scheme.md §5b). 0 on
      ! ESDIRK2.
      !---~---
      real(wp)             :: theta_ood_max = 0.0_wp   ! [m3/m3] max excursion outside [theta_res, theta_sat]
   end type step_report_t

   !---~---
   !   What one fast step reports: the report, the aerodynamics, each cohort's carbon fluxes and
   ! leaf potential, the canopy air's exchange with the atmosphere, and the per-cohort
   ! diagnostic capture (docs/science/diagnostics.md §1, stage [2]) when the run asks for it.
   !---~---
   type :: fast_output_t
      type(step_report_t)   :: report
      type(aero_out_t)      :: aero
      real(wp), allocatable :: gpp_coh(:)         ! [umol CO2/plant/s] gross GPP (fast->slow)
      real(wp), allocatable :: leaf_resp_coh(:)   ! [umol CO2/plant/s] leaf dark respiration
      real(wp), allocatable :: stem_resp_coh(:)   ! [umol CO2/plant/s] stem maintenance respiration
      real(wp), allocatable :: root_resp_coh(:)   ! [umol CO2/plant/s] fine-root maintenance respiration
      real(wp), allocatable :: psi_leaf_coh(:)    ! [MPa] the leaf potential at state^n (daily-max accumulator)
      real(wp)              :: le_flux = 0.0_wp   ! [W/m2] CAS -> atmosphere latent heat (ET), the step's mean
      real(wp)              :: h_flux  = 0.0_wp   ! [W/m2] CAS -> atmosphere sensible heat, the step's mean
      real(wp)              :: rnet    = 0.0_wp   ! [W/m2] net radiation the canopy and ground booked, the step's mean
      logical               :: capture = .false.  ! fill cdiag; off, nothing extra is computed
      real(wp), allocatable :: cdiag(:,:)         ! (N_CO_FAST_ROWS, cohort) per-cohort diagnostics
   end type fast_output_t

   !--- The prognostic CAS surface state advanced by the fast loop.
   type :: surface_state_t
      real(wp) :: cas_enthalpy = 0.0_wp    ! [J/kg]      canopy-air specific enthalpy
      real(wp) :: cas_shv      = 0.0_wp    ! [kg/kg]     canopy-air specific humidity
      real(wp) :: cas_co2      = 0.0_wp    ! [umol/mol]  canopy-air CO2 mixing ratio
   end type surface_state_t

   !---~---
   !   The frozen record, by physical content (docs/science/numerical_scheme.md §2): what the
   ! pre-pass (build_frozen) evaluates once at state^n and every stage of either
   ! integrator reads. A kernel takes only the pieces it reads (surface_derivs: cas, tissue,
   ! film, ground, snow). The CAS <-> atmosphere conductances are re-solved per stage.
   !---~---

   !---~---
   !   CAS <-> atmosphere boundary: the canopy-air box (its capacities, its conductances to the
   ! atmosphere and the reference-level state, the box kernel's own record), the surface-layer
   ! inputs a per-stage re-solve needs, and the biotic CO2 source. Plain scalars, so a stage
   ! can copy it and override the conductances.
   !---~---
   type :: cas_boundary_t
      type(cas_env_t) :: box   ! its sources are each evaluation's, set where the box is advanced
      !---~---
      !   The box's conductances are the state^n values, which seed the march and the diagnostics; a stage
      ! re-solves them at its own canopy-air state (docs/science/numerical_scheme.md §2, item 3;
      ! docs/science/canopy_aerodynamics.md §2). mo_live says env and geom below are set: it
      ! is false for a hand-built record, which supplies its own conductances (a re-solve from
      ! unset roughness and wind would not converge), true from build_frozen, and false
      ! again in ESDIRK2's stage copy once refreshed, so the Newton's repeated evaluations reuse
      ! it.
      !---~---
      logical :: mo_live = .false.            ! env and geom are set => re-solve per evaluation
      type(aero_params_t) :: aero_cfg         ! MO constants (copy; POD)
      type(aero_env_t)    :: env              ! the step's aerodynamic environment: the reference level, the canopy air
      type(aero_args_t)   :: geom             ! the canopy geometry the roughness and displacement come from
      real(wp) :: nee_biotic    = 0.0_wp      ! [umol/m2/s] frozen biotic CO2 source (Ra+Rh-GPP)
      real(wp) :: rho           = 0.0_wp      ! [kg/m3]   canopy-air density
      real(wp) :: press         = 0.0_wp      ! [Pa]      canopy-air pressure
   end type cas_boundary_t

   !---~---
   !   Surface-block tendencies and the diagnostics the ledgers and the plant-water update use
   ! (transp_c is each cohort's leaf water loss).
   !---~---
   type :: surface_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp     ! [J/kg/s]     dH/dt
      real(wp) :: d_cas_shv      = 0.0_wp     ! [kg/kg/s]    dq/dt
      real(wp) :: d_cas_co2      = 0.0_wp     ! [umol/mol/s] dC/dt
      real(wp) :: src_enth       = 0.0_wp     ! [W/m2]     summed surface enthalpy source into the CAS
      real(wp) :: src_vap        = 0.0_wp     ! [kg/m2/s]  summed surface vapour source into the CAS
      real(wp) :: g_top          = 0.0_wp     ! [W/m2]     net energy into the soil-top store
      real(wp) :: h_ground       = 0.0_wp     ! [W/m2]     ground sensible flux to the CAS
      real(wp) :: le_ground      = 0.0_wp     ! [W/m2]     ground latent flux to the CAS
      real(wp) :: h_canopy       = 0.0_wp     ! [W/m2]     the leaves' and wood's sensible flux to the CAS
      real(wp) :: coh_rnet       = 0.0_wp     ! [W/m2]     net radiation absorbed by the canopy
      real(wp) :: cond           = 0.0_wp     ! [kg/m2/s]  smooth condensation sink (dew) draining CAS supersat
      real(wp) :: cond_enth      = 0.0_wp     ! [W/m2]     liquid enthalpy the sink debited from the CAS
      real(wp) :: ground_rad     = 0.0_wp     ! [W/m2]     the ground's net radiation, snow's and bare soil's
      real(wp) :: snow_enth_in   = 0.0_wp     ! [W/m2]     energy into the snow pack (0 without one)
      real(wp) :: snow_vap_out   = 0.0_wp     ! [kg/m2/s]  vapour from the snow pack to the CAS
      real(wp), allocatable :: leaf_temp(:)   ! [K]        diagnosed per-cohort leaf temperature
      real(wp), allocatable :: wood_temp(:)   ! [K]        diagnosed per-cohort wood temperature
      real(wp), allocatable :: transp_c(:)    ! [kg/m2/s]  per-cohort transpiration demand
      !---~---
      !   Per-cohort film evaporation (dew if negative), the films' analogue of transp_c. Zero
      ! when canopy_water_on is off.
      !---~---
      real(wp), allocatable :: leaf_surface_evap(:), wood_surface_evap(:)   ! [kg/m2 ground/s]
   end type surface_tend_t

   !---~---
   !   The prognostic column state the integrators advance: a flat copy of the patch's stores
   ! that the schemes combine and measure as one vector. Plant water is carried as mass; psi
   ! is diagnosed from it once per dt_fast for the pre-pass (docs/science/soil_plant_water.md).
   !---~---
   type :: march_state_t
      real(wp) :: cas_enthalpy = 0.0_wp                    ! [J/kg]
      real(wp) :: cas_shv      = 0.0_wp                    ! [kg/kg]
      real(wp) :: cas_co2      = 0.0_wp                    ! [umol/mol]
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp   ! [J/m3]   per soil layer
      real(wp) :: theta(n_soil_layer_max)       = 0.0_wp   ! [m3/m3]  per soil layer
      real(wp), allocatable :: leaf_water_mass(:) ! [kg/plant] internal leaf water (ncoh)
      real(wp), allocatable :: wood_water_mass(:) ! [kg/plant] internal wood water (ncoh)
      !---~---
      !   The pond: carried through the stages unchanged and committed from the water network.
      ! Like the films it carries no error group, so the controller scores no frozen quantity.
      !---~---
      real(wp) :: w_surface      = 0.0_wp         ! [kg/m2] ponded surface water
      real(wp) :: w_surface_energy = 0.0_wp         ! [J/m2]  enthalpy of the ponded water
      !---~---
      !   The interception films, distinct from the internal water above: a film evaporates with
      ! no stomatal resistance. Per m2 of ground, unlike the per-plant fields above.
      !---~---
      real(wp), allocatable :: leaf_surf_water(:) ! [kg/m2 ground] canopy interception film on leaf
      real(wp), allocatable :: wood_surf_water(:) ! [kg/m2 ground] canopy interception film on wood
   end type march_state_t

   !---~---
   !   Soil water over this dt_fast, from the soil-plant water network's step (water_step, run by
   ! fast_step after build_frozen): its boundary flows, end-of-step stores, the temperatures that value
   ! the water crossing each boundary, and its interior faces. Every scheme commits its theta and
   ! advects soil enthalpy on these same faces (docs/science/numerical_scheme.md §5b).
   !---~---
   type :: soil_hydrology_t
      real(wp) :: q_top         = 0.0_wp          ! [m/s]     net top-face water flux (infiltration - exfiltration - evaporation)
      real(wp) :: infiltration  = 0.0_wp          ! [kg/m2/s] pond -> top layer
      real(wp) :: drainage      = 0.0_wp          ! [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf   = 0.0_wp          ! [kg/m2/s] surface runoff
      real(wp) :: ground_inflow_enth = 0.0_wp     ! [W/m2]    the enthalpy it brings across the column's boundary
      real(wp) :: t_infil       = 0.0_wp          ! [K]       temperature of the infiltrating (pond) water
      real(wp) :: w_surface_energy1 = 0.0_wp        ! [J/m2]    end-of-step pond ENTHALPY
      real(wp) :: runoff_enth   = 0.0_wp          ! [W/m2]    enthalpy leaving with surface runoff
      !---~---
      !   Water a supersaturated top layer pushed up into the pond, and the enthalpy it took, valued at
      ! the layer's start-of-step temperature (the pond received the same number).
      !---~---
      real(wp) :: exfiltration  = 0.0_wp          ! [kg/m2/s] top layer -> pond
      real(wp) :: exfil_enth    = 0.0_wp          ! [W/m2]    enthalpy leaving the top layer with it
      real(wp) :: t_bot         = 0.0_wp          ! [K]       bottom-layer soil temperature @ state^n
      real(wp) :: w_surface1    = 0.0_wp          ! [kg/m2]   end-of-step ponded surface water
      real(wp) :: w_flux_frozen(n_soil_layer_max) = 0.0_wp  ! [m/s]   DOWNWARD interior face flux, k=1..nsl-1
      !--- The network's end-of-step soil moisture, which every scheme commits.
      real(wp), allocatable :: theta1(:)          ! [m3/m3]   committed post-step soil moisture (per layer)
      !---~---
      !   (theta1 - theta^n) / dt_fast: the steady rate at which the stages move soil water, so the
      ! water a layer holds keeps pace with the enthalpy its faces carry in (see
      ! backward_euler_stage).
      !---~---
      real(wp) :: theta_rate(n_soil_layer_max) = 0.0_wp   ! [1/s]
   end type soil_hydrology_t

   !---~---
   !   Root zone over this dt_fast, from the water network: the patch's root uptake, its shares by
   ! layer (sum 1; the root profile when the roots took nothing), and the liquid enthalpy it
   ! carries out of the soil, by cohort (what each wood is credited) and by layer, at each
   ! layer's start-of-step temperature.
   !---~---
   type :: root_zone_t
      real(wp) :: uptake        = 0.0_wp          ! [kg/m2/s] the patch's root uptake
      real(wp), allocatable :: root_share(:)      ! [-]       per-layer root-sink shares (sum = 1)
      real(wp), allocatable :: qloss_frozen(:)    ! [W/m2 ground] (ncoh) advected enthalpy leaving the soil with uptake
      !--- By layer [W/m2], so the soil loses heat where it loses water (docs/science/soil_plant_water.md §5).
      real(wp), allocatable :: root_heat(:)       ! (nsl)
   end type root_zone_t

   !---~---
   !   Plant water over this dt_fast, from the water network, per plant and constant across the
   ! stages (sapflow wood -> leaf, uptake soil -> wood). With the held transpiration (transp_net)
   ! they move each cohort's water exactly to the network's end-of-step values.
   !---~---
   type :: plant_water_t
      real(wp), allocatable :: nplant(:)                            ! [plants/m2] (ncoh)
      real(wp), allocatable :: sapflow_frozen(:), uptake_frozen(:)   ! [kg/plant/s] (ncoh)
   end type plant_water_t

   !---~---
   !   What drives the soil-plant water network over the step, held fixed across it: the water
   ! reaching the pond and the temperature that values it, the temperature of the water each
   ! layer gives up, each cohort's transpiration demand and the ground evaporation's boundary
   ! (docs/science/soil_plant_water.md §4-§5).
   !---~---
   type :: water_drive_t
      real(wp) :: rain   = 0.0_wp                          ! [kg/m2/s] water reaching the pond
      real(wp) :: t_rain = 298.15_wp                       ! [K] the temperature that values it
      real(wp) :: t_soil(n_soil_layer_max) = 298.15_wp     ! [K] temperature of the water a layer gives up
      type(soil_ice_t) :: ice(n_soil_layer_max)              ! each layer's ice, held over the step
      real(wp), allocatable :: e_demand(:)                 ! (ncoh) [kg/m2/s] the stomata's transpiration demand
      type(ground_evap_env_t) :: ground                    ! the ground evaporation's boundary
   end type water_drive_t

   !---~---
   !   The water network at one state (water_derivs): each store's potential [MPa] and water
   ! [kg/m2] with their slopes in its unknown, each soil node's conductivity [m/s], the net
   ! inflow into each store [kg/m2/s], and the flow along each path [kg/m2/s].
   !---~---
   type :: water_tend_t
      real(wp), allocatable :: psi_soil(:), dpsi_soil(:), water_soil(:), dwater_soil(:)   ! (nsoil)
      real(wp), allocatable :: k_soil(:), dk_soil(:), inflow_soil(:)                       ! (nsoil)
      real(wp), allocatable :: psi_coh(:,:), dpsi_coh(:,:), water_coh(:,:), dwater_coh(:,:) ! (2, ncoh)
      real(wp), allocatable :: inflow_coh(:,:)                                             ! (2, ncoh)
      real(wp), allocatable :: face_flux(:)        ! (nsl) down through the face below layer k
      real(wp), allocatable :: uptake(:,:)         ! (nsl, ncoh) from each layer into each cohort
      real(wp), allocatable :: sapflow(:), transp(:)   ! (ncoh)
      real(wp) :: pond_flux = 0.0_wp               ! pond -> top layer
      real(wp) :: drain = 0.0_wp, evap = 0.0_wp    ! out of the bottom; out of the top layer to the air
   end type water_tend_t

   !---~---
   !   The water network's state: each store's unknown and the water it holds [kg/m2], the pond's
   ! enthalpy, and each store's potential at the end of the step.
   !---~---
   type :: water_state_t
      real(wp), allocatable :: tau_soil(:), tau_coh(:,:)       ! (nsoil), (2, ncoh)
      real(wp), allocatable :: water_soil(:), water_coh(:,:)   ! (nsoil), (2, ncoh)
      real(wp), allocatable :: psi_soil(:), psi_coh(:,:)       ! (nsoil), (2, ncoh) [MPa]
      real(wp) :: e_pond = 0.0_wp                              ! [J/m2] the pond's enthalpy
   end type water_state_t

   !--- What each of the network's paths moved over the step [kg/m2, J/m2].
   type :: water_flows_t
      real(wp), allocatable :: face(:)          ! (nsl) down through the face below layer k
      real(wp), allocatable :: uptake(:,:)      ! (nsl, ncoh) from each layer into each cohort
      real(wp), allocatable :: sapflow(:), transp(:)   ! (ncoh)
      real(wp) :: infil = 0.0_wp, infil_enth = 0.0_wp  ! pond -> top layer, and its enthalpy
      real(wp) :: exfil = 0.0_wp, exfil_enth = 0.0_wp  ! top layer -> pond
      real(wp) :: drain = 0.0_wp, evap = 0.0_wp        ! out of the bottom; out of the top layer to the air
      real(wp) :: spill = 0.0_wp, spill_enth = 0.0_wp  ! pond runoff over its capacity
      real(wp) :: rain_enth = 0.0_wp                   ! enthalpy the rain brought
   end type water_flows_t

   !---~---
   !   One thread's water network, kept from step to step so its arrays are allocated again only
   ! when the shape changes: its size, its state and flows, and the Newton scratch. ev(i_start)
   ! is the network at a sub-step's start; Newton evaluates its iterates into the other. sys
   ! holds the flows' Jacobian of the latest evaluation (its factorisation is kept apart). A
   ! sub-step starts where the last accepted one ended, so while start_ready holds, its start's
   ! evaluation and Jacobian are already there.
   !---~---
   type :: water_work_t
      integer(ik)              :: nsl   = 0_ik      ! soil layers
      integer(ik)              :: ncoh  = 0_ik      ! cohorts
      integer(ik)              :: nsoil = 0_ik      ! soil nodes: the pond and the nsl layers
      type(water_state_t)      :: state
      type(water_flows_t)      :: flows
      type(arrowhead_system_t) :: sys
      type(water_tend_t)       :: ev(2)
      integer(ik)              :: i_start = 1_ik
      logical                  :: start_ready = .false.
   end type water_work_t

   !--- The container: everything held constant over one dt_fast, by physical content.
   type :: frozen_t
      type(cas_boundary_t)         :: cas          ! CAS <-> atmosphere boundary
      type(plant_energy_env_t)     :: tissue       ! per-cohort leaf/wood energy coefficients + heat store
      type(canopy_water_out_t)     :: veg_surface         ! canopy interception film
      type(surface_env_t)          :: ground       ! what the ground surface reads, the snow pack among it
      type(snow_out_t)             :: snow         ! the snow stage's hand-off (all zero without a pack)
      type(soil_hydrology_t)       :: hydrology    ! the water network's soil flows over the step
      type(root_zone_t)            :: roots        ! the network's root uptake and its placement by layer
      type(plant_water_t)          :: plant        ! the network's sapflow and uptake per plant
      type(water_drive_t)          :: water_drive  ! what drives the water network over the step
   end type frozen_t

   !---~---
   !   Storage one thread's fast step reuses from step to step: the aerodynamic environment, the
   ! soil-plant water network, the frozen record, and the ESDIRK2 march's states and stages,
   ! filled in place and reallocated only when the cohort count changes. The stages are their own record so the
   ! march passes them beside its states.
   !---~---
   type :: esdirk2_stages_t
      type(march_state_t) :: y2, base3, y3            ! the ARS(2,2,2) stages
      type(surface_tend_t) :: sf2, sf3                 ! the stages' surface tendencies
   end type esdirk2_stages_t

   type :: fast_workspace_t
      type(water_work_t)    :: water                         ! the soil-plant water network
      type(frozen_t)        :: frozen                        ! the frozen record (prepare_frozen resets it)
      type(march_state_t)  :: y_march, y_new, y_err, y_lo   ! the ESDIRK2 march's state and trial step
      type(esdirk2_stages_t)    :: stages                        ! the ESDIRK2 step's stages
   end type fast_workspace_t

   !--- The whole-column tendency vector + diagnostics.
   type :: march_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp
      real(wp) :: d_cas_shv      = 0.0_wp
      real(wp) :: d_cas_co2      = 0.0_wp
      real(wp) :: dedt(n_soil_layer_max)   = 0.0_wp   ! [W/m3] dsoil_energy/dt
      real(wp) :: dtheta_dt(n_soil_layer_max) = 0.0_wp! [1/s]  dtheta/dt
      real(wp), allocatable :: d_leaf_water_mass(:)   ! [kg/plant/s] frozen_sapflow - transp(stage)
      real(wp), allocatable :: d_wood_water_mass(:)   ! [kg/plant/s] frozen_uptake  - frozen_sapflow
      !--- The films' analogue of d_leaf_water_mass/d_wood_water_mass.
      real(wp), allocatable :: d_leaf_surf_water(:)   ! [kg/m2 ground/s] frozen_intercept_leaf - film_evap
      real(wp), allocatable :: d_wood_surf_water(:)   ! [kg/m2 ground/s] frozen_intercept_wood - film_evap
      real(wp) :: g_top = 0.0_wp, drainage_rate = 0.0_wp, uptake_rate = 0.0_wp
      real(wp) :: g_bottom = 0.0_wp   ! [W/m2] heat across the soil column's bottom face, positive UP
      real(wp), allocatable :: leaf_temp(:)
   end type march_tend_t

   !---~---
   !   ESDIRK2 conservation ledger: per-stage boundary-flux rates (from backward_euler_stage), and
   ! their b-weighted amounts summed over the accepted sub-steps (boundary_flux_t). Since Y3 - y
   ! = (1-gamma)*h*K2 + gamma*h*K3 (docs/science/numerical_scheme.md eq. (1)), the amounts
   ! telescope against the committed store change to machine precision.
   !---~---
   type :: stage_bflux_t                                    ! per-stage RATES
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp    ! [W/m2]
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp    ! [kg/m2/s]
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp    ! [umol/m2/s]
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp  ! [W/m2]
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp  ! [kg/m2/s]
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp! [W/m2]
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp! [kg/m2/s]
      real(wp) :: whole_cond    = 0.0_wp                         ! [kg/m2/s] condensate
      real(wp) :: whole_cond_enth = 0.0_wp                       ! [W/m2] its liquid enthalpy at the stage CAS temperature
      real(wp) :: snow_enth_in = 0.0_wp                          ! [W/m2]    energy into the snow pack
      real(wp) :: snow_vap_out = 0.0_wp                          ! [kg/m2/s] vapour from it to the CAS
      real(wp) :: veg_surface_evap = 0.0_wp                      ! [kg/m2/s] leaf and wood surface water evaporated (dew < 0)
      real(wp) :: h_canopy     = 0.0_wp                          ! [W/m2]    the leaves' and wood's sensible heat to the CAS
      real(wp) :: h_ground     = 0.0_wp                          ! [W/m2]    the ground's (soil, snow) sensible heat to the CAS
      !---~---
      !   Net CAS -> atmosphere export at the ledger's own stage states and conductances: the
      ! reported H and LE (ed_fast_integ_stepper%atm_fluxes; docs/science/numerical_scheme.md §5b).
      !---~---
      real(wp) :: atm_heat_out = 0.0_wp                          ! [W/m2]    g_atm_heat*cp_air*(T_cas - theta_atm)
      real(wp) :: atm_vap_out  = 0.0_wp                          ! [kg/m2/s] g_atm_vapour*(q_cas - q_atm)
      !---~---
      !   The net radiation the canopy and the ground absorbed, as their balances booked it: the
      ! reported net radiation (docs/science/vegetation_energy_dynamics.md, "The emission base").
      !---~---
      real(wp) :: rad_net = 0.0_wp                               ! [W/m2]
   end type stage_bflux_t

   type :: boundary_flux_t                                  ! accumulated AMOUNTS (J/m2, kg/m2, umol/m2)
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp
      !---~---
      !   Condensate, apart from whole_wat_out: dew and fog land inside the column, and the caller
      ! deposits them into soil layer 1 (docs/science/numerical_scheme.md §7, item 3).
      !---~---
      real(wp) :: whole_cond    = 0.0_wp   ! [kg/m2] condensed vapour over the step (>= 0)
      !---~---
      !   ...and the liquid enthalpy it left the CAS with, b-weighted at each stage's own CAS
      ! temperature, as surface_derivs debited it; the deposit carries this, not u_liq(T_end).
      !---~---
      real(wp) :: whole_cond_enth = 0.0_wp ! [J/m2]
      !--- The snow pack's exchange over the step, which the update commits to the pack.
      real(wp) :: snow_enth_in = 0.0_wp    ! [J/m2]  energy into the pack
      real(wp) :: snow_vap_out = 0.0_wp    ! [kg/m2] vapour from it to the CAS
      real(wp) :: veg_surface_evap    = 0.0_wp    ! [kg/m2] leaf and wood surface water evaporated (dew < 0)
      real(wp) :: h_canopy     = 0.0_wp    ! [J/m2]  the leaves' and wood's sensible heat to the CAS
      real(wp) :: h_ground     = 0.0_wp    ! [J/m2]  the ground's (soil, snow) sensible heat to the CAS
      real(wp) :: atm_heat_out = 0.0_wp    ! [J/m2]  net CAS -> atmosphere sensible-heat export (reported H)
      real(wp) :: atm_vap_out  = 0.0_wp    ! [kg/m2] net CAS -> atmosphere vapour export (reported LE / L_v)
      real(wp) :: rad_net      = 0.0_wp    ! [J/m2]  net radiation the canopy and the ground booked (reported)
      !---~---
      !   Per-cohort tissue-temperature time integrals [K*s], b-weighted over the stages and summed
      ! over the accepted sub-steps. A tissue's energy is committed from this time average, the
      ! complement of what the canopy air received (docs/science/numerical_scheme.md §5b).
      !---~---
      real(wp), allocatable :: tissue_leaf_int(:), tissue_wood_int(:)   ! [K*s]
      !---~---
      !   Each cohort's transpiration over the step [kg/m2 ground], b-weighted the same way: the
      ! water its leaves were debited and the canopy air credited, which the diagnostics report.
      !---~---
      real(wp), allocatable :: transp_coh(:)
   end type boundary_flux_t

contains

   !--- Allocate a fast_static_state_t's cohort arrays and set their defaults.
   subroutine alloc_static_state(static, n)
      type(fast_static_state_t), intent(out) :: static
      integer(ik),               intent(in)  :: n
      integer(ik) :: i
      static%n = n
      allocate(static%pft(n), static%order(n), static%lai(n), static%wai(n), static%height(n),                    &
               static%leaf_width(n), static%branch_diam(n), static%aboveground_frac(n),                   &
               static%is_woody(n), static%stem_resp_factor25(n), static%root_resp_factor25(n),            &
               static%stem_resp_sapwood25(n), static%sap_volume(n),                                         &
               static%leaf_emissivity(n), static%wood_emissivity(n),                                         &
               static%leaf_area(n), static%nplant(n),                                                          &
               static%dbh(n), static%broot(n), static%bleaf(n), static%bsap(n), static%sap_area(n),           &
               static%bwood(n), static%vcmax25(n), static%rd25(n), static%dmax_psi_leaf(n))
      allocate(static%root_frac(n_soil_layer_max, n), static%k_xylem(n), static%g_rhizo(n_soil_layer_max, n),      &
               static%crown_head(n), static%bwood_water(n), static%leaf_dry_hcap(n), static%wood_dry_hcap(n))
      static%pft = 1_ik
      static%order = [(n - i + 1_ik, i = 1_ik, n)]   ! tallest first, as the gather leaves them
      static%lai = 0.0_wp
      static%wai = 0.0_wp
      static%height = 0.0_wp
      static%leaf_width = 0.04_wp
      static%branch_diam = 0.02_wp
      static%aboveground_frac = 0.7_wp
      static%is_woody = .true.
      static%stem_resp_factor25 = 0.0_wp
      static%root_resp_factor25 = 0.0_wp
      static%stem_resp_sapwood25 = 0.0_wp
      static%sap_volume = 0.0_wp
      static%leaf_emissivity = 1.0_wp
      static%wood_emissivity = 1.0_wp
      static%leaf_area = 0.0_wp
      static%nplant = 0.0_wp
      static%dbh = 0.0_wp
      static%broot = 0.0_wp
      static%bleaf = 0.0_wp
      static%bsap = 0.0_wp
      static%sap_area = 0.0_wp
      static%bwood = 0.0_wp
      static%vcmax25 = 0.0_wp
      static%rd25 = 0.0_wp
      !---~---
      !   Unset, not 0: a 0 would read as fully turgid and disable the stomatal stress limb for a
      ! caller that forgets to fill it; the sentinel makes the frozen build seed from the soil.
      !---~---
      static%dmax_psi_leaf = UNSET_NEGATIVE
      static%root_frac = 0.0_wp
      static%k_xylem = 0.0_wp
      static%g_rhizo = 0.0_wp
      static%crown_head = 0.0_wp
      static%bwood_water = 0.0_wp
      static%leaf_dry_hcap = 0.0_wp
      static%wood_dry_hcap = 0.0_wp
   end subroutine alloc_static_state

   !---~---
   !   Grow-only capacity: reallocate only when the arrays hold fewer than n cohorts, else just
   ! set the active count. Every reader loops to static%n and the gather overwrites 1..n, so
   ! the caller sizes it once per fast_dynamics call, to the largest patch.
   !---~---
   subroutine ensure_static_capacity(static, n)
      type(fast_static_state_t), intent(inout) :: static
      integer(ik),               intent(in)    :: n
      if (.not. allocated(static%pft)) then
         call alloc_static_state(static, n)
      else if (size(static%pft) < n) then
         call alloc_static_state(static, n)
      else
         static%n = n
      end if
   end subroutine ensure_static_capacity

   !--- Allocate a fast_dynamic_state_t and seed it from an initial canopy-air state.
   subroutine alloc_dynamic_state(dynamic, n_coh, can_temp0, can_shv0, can_co2, leaf_temp0)
      type(fast_dynamic_state_t), intent(out) :: dynamic
      integer(ik),                intent(in)  :: n_coh
      real(wp),                   intent(in)  :: can_temp0, can_shv0, can_co2, leaf_temp0
      allocate(dynamic%leaf_energy(n_coh), dynamic%wood_energy(n_coh))
      allocate(dynamic%leaf_temp(n_coh), dynamic%wood_temp(n_coh))
      allocate(dynamic%leaf_water_mass(n_coh), dynamic%wood_water_mass(n_coh))
      allocate(dynamic%leaf_surf_water(n_coh), dynamic%wood_surf_water(n_coh))
      dynamic%leaf_energy      = 0.0_wp    ! scratch seed: the gather sets it, a column fixture with tissue_temp_to_energy
      dynamic%wood_energy      = 0.0_wp
      dynamic%leaf_temp        = leaf_temp0
      dynamic%wood_temp        = leaf_temp0
      dynamic%leaf_water_mass  = 0.0_wp    ! scratch seeds: the next gather overwrites them
      dynamic%wood_water_mass  = 0.0_wp
      dynamic%leaf_surf_water  = 0.0_wp
      dynamic%wood_surf_water  = 0.0_wp
      dynamic%cas%can_temp     = can_temp0
      dynamic%cas%can_shv      = can_shv0
      dynamic%cas%can_co2      = can_co2
      dynamic%cas%can_enthalpy = cas_temp_to_enthalpy(can_temp0, can_shv0)
   end subroutine alloc_dynamic_state

   !---~---
   !   Grow-only capacity for the per-cohort arrays, as ensure_static_capacity. Nothing is seeded:
   ! every caller overwrites the patch stores and cohorts 1..n with the patch's persisted state
   ! right after.
   !---~---
   !---~---
   !   Size the water network for nsl soil layers and ncoh cohorts, keeping the storage when the
   ! shape is unchanged.
   !---~---
   subroutine ensure_water_capacity(water, nsl, ncoh)
      type(water_work_t), intent(inout) :: water
      integer(ik),        intent(in)    :: nsl, ncoh
      integer(ik) :: ns
      if (allocated(water%state%tau_soil) .and. water%nsl == nsl .and. water%ncoh == ncoh) return
      water%nsl   = nsl
      water%ncoh  = ncoh
      water%nsoil = 1_ik + nsl
      ns = water%nsoil
      associate (st => water%state, fl => water%flows)
         if (allocated(st%tau_soil)) deallocate(st%tau_soil, st%tau_coh, st%water_soil, st%water_coh,          &
                                                st%psi_soil, st%psi_coh)
         allocate(st%tau_soil(ns), st%tau_coh(COHORT_STORES, ncoh), st%water_soil(ns),                      &
                  st%water_coh(COHORT_STORES, ncoh), st%psi_soil(ns), st%psi_coh(COHORT_STORES, ncoh))
         if (allocated(fl%face)) deallocate(fl%face, fl%uptake, fl%sapflow, fl%transp)
         allocate(fl%face(nsl), fl%uptake(nsl, ncoh), fl%sapflow(ncoh), fl%transp(ncoh))
      end associate
      call tend_allocate(water%ev(1))
      call tend_allocate(water%ev(2))
   contains
      subroutine tend_allocate(ev)
         type(water_tend_t), intent(inout) :: ev
         if (allocated(ev%psi_soil)) deallocate(ev%psi_soil, ev%dpsi_soil, ev%water_soil, ev%dwater_soil,       &
                                                ev%k_soil, ev%dk_soil, ev%inflow_soil, ev%psi_coh, ev%dpsi_coh, &
                                                ev%water_coh, ev%dwater_coh, ev%inflow_coh, ev%face_flux,       &
                                                ev%uptake, ev%sapflow, ev%transp)
         allocate(ev%psi_soil(ns), ev%dpsi_soil(ns), ev%water_soil(ns), ev%dwater_soil(ns), ev%k_soil(ns),     &
                  ev%dk_soil(ns), ev%inflow_soil(ns))
         allocate(ev%psi_coh(COHORT_STORES, ncoh), ev%dpsi_coh(COHORT_STORES, ncoh),                           &
                  ev%water_coh(COHORT_STORES, ncoh), ev%dwater_coh(COHORT_STORES, ncoh),                       &
                  ev%inflow_coh(COHORT_STORES, ncoh))
         allocate(ev%face_flux(nsl), ev%uptake(nsl, ncoh), ev%sapflow(ncoh), ev%transp(ncoh))
         ev%face_flux = 0.0_wp
      end subroutine tend_allocate
   end subroutine ensure_water_capacity

   subroutine ensure_dynamic_capacity(dynamic, n_coh)
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      integer(ik),                intent(in)    :: n_coh
      if (allocated(dynamic%leaf_temp)) then
         if (size(dynamic%leaf_temp) >= n_coh) return
         deallocate(dynamic%leaf_energy, dynamic%wood_energy, dynamic%leaf_temp, dynamic%wood_temp,             &
                    dynamic%leaf_water_mass, dynamic%wood_water_mass, dynamic%leaf_surf_water, dynamic%wood_surf_water)
      end if
      allocate(dynamic%leaf_energy(n_coh), dynamic%wood_energy(n_coh), dynamic%leaf_temp(n_coh),               &
               dynamic%wood_temp(n_coh), dynamic%leaf_water_mass(n_coh), dynamic%wood_water_mass(n_coh),        &
               dynamic%leaf_surf_water(n_coh), dynamic%wood_surf_water(n_coh), source=0.0_wp)
   end subroutine ensure_dynamic_capacity

   !---~---
   !   Grow-only capacity of the per-cohort forcing; its writers loop to n, so a larger capacity
   ! is exact.
   !---~---
   subroutine ensure_forcing_capacity(forcing, n)
      type(fast_forcing_t), intent(inout) :: forcing
      integer(ik),          intent(in)    :: n
      if (allocated(forcing%abs_sw)) then
         if (size(forcing%abs_sw) >= n) return
         deallocate(forcing%abs_sw, forcing%abs_lw, forcing%abs_par, forcing%abs_sw_wood, forcing%abs_lw_wood,  &
                    forcing%par_top)
      end if
      allocate(forcing%abs_sw(n), forcing%abs_lw(n), forcing%abs_par(n), forcing%abs_sw_wood(n),               &
               forcing%abs_lw_wood(n), forcing%par_top(n))
   end subroutine ensure_forcing_capacity

   !---~---
   !   Grow-only capacity of the per-cohort outputs, the aerodynamics' and, with `capture`, the
   ! diagnostic rows; every writer loops to n.
   !---~---
   subroutine ensure_output_capacity(out, n, capture)
      type(fast_output_t), intent(inout) :: out
      integer(ik),         intent(in)    :: n
      logical,             intent(in)    :: capture
      call ensure_aero_out_capacity(out%aero, n)
      if (allocated(out%gpp_coh)) then
         if (size(out%gpp_coh) < n) deallocate(out%gpp_coh, out%leaf_resp_coh, out%stem_resp_coh, out%root_resp_coh, &
                                               out%psi_leaf_coh)
      end if
      if (.not. allocated(out%gpp_coh)) allocate(out%gpp_coh(n), out%leaf_resp_coh(n), out%stem_resp_coh(n),       &
                                                 out%root_resp_coh(n), out%psi_leaf_coh(n), source=0.0_wp)
      out%capture = capture
      if (.not. capture) return
      if (allocated(out%cdiag)) then
         if (size(out%cdiag, 2) < n) deallocate(out%cdiag)
      end if
      if (.not. allocated(out%cdiag)) allocate(out%cdiag(N_CO_FAST_ROWS, n), source=0.0_wp)
   end subroutine ensure_output_capacity


   !---~---
   !   The state's fields, listed once. The integrator's bookkeeping works on one flat array per
   ! state; the three routines below are the only ones that name march_state_t's fields, and
   ! state_entry_rules the only one that says how each is treated. Fortran cannot check that a
   ! field is listed, so test_fast_integ_kernels runs every field through every combinator.
   ! Layout: CAS enthalpy, humidity, CO2 | soil energy (1:nsl) | soil water (1:nsl) |
   !   per cohort: leaf, wood water | per cohort: leaf, wood film | pond: water, enthalpy
   !---~---
   pure subroutine state_to_array(y, n, nsl, v)
      type(march_state_t), intent(in)  :: y
      integer(ik),          intent(in)  :: n, nsl
      real(wp),             intent(out) :: v(:)
      integer(ik) :: i, o
      v(1) = y%cas_enthalpy
      v(2) = y%cas_shv
      v(3) = y%cas_co2
      v(4:3+nsl)       = y%soil_energy(1:nsl)
      v(4+nsl:3+2*nsl) = y%theta(1:nsl)
      o = 3_ik + 2_ik*nsl
      do i = 1_ik, n
         v(o+2*i-1) = y%leaf_water_mass(i)
         v(o+2*i) = y%wood_water_mass(i)
      end do
      o = o + 2_ik*n
      do i = 1_ik, n
         v(o+2*i-1) = y%leaf_surf_water(i)
         v(o+2*i) = y%wood_surf_water(i)
      end do
      o = o + 2_ik*n
      v(o+1) = y%w_surface
      v(o+2) = y%w_surface_energy
   end subroutine state_to_array

   !---~---
   !   The inverse. The per-cohort arrays are allocated only when they are missing or of another
   ! length, so a state reused from call to call is filled in place. Soil layers past nsl are 0.
   !---~---
   pure subroutine array_to_state(v, n, nsl, y)
      real(wp),             intent(in)    :: v(:)
      integer(ik),          intent(in)    :: n, nsl
      type(march_state_t), intent(inout) :: y
      integer(ik) :: i, o
      call ensure_cohort_arrays(y, n)
      y%cas_enthalpy = v(1)
      y%cas_shv = v(2)
      y%cas_co2 = v(3)
      y%soil_energy = 0.0_wp
      y%theta = 0.0_wp
      y%soil_energy(1:nsl) = v(4:3+nsl)
      y%theta(1:nsl)       = v(4+nsl:3+2*nsl)
      o = 3_ik + 2_ik*nsl
      do i = 1_ik, n
         y%leaf_water_mass(i) = v(o+2*i-1)
         y%wood_water_mass(i) = v(o+2*i)
      end do
      o = o + 2_ik*n
      do i = 1_ik, n
         y%leaf_surf_water(i) = v(o+2*i-1)
         y%wood_surf_water(i) = v(o+2*i)
      end do
      o = o + 2_ik*n
      y%w_surface = v(o+1)
      y%w_surface_energy = v(o+2)
   end subroutine array_to_state

   !---~---
   !   A tendency in the same layout. The pond has none (it is committed from the water network),
   ! so its entries are 0 and y + a*k leaves it as it was.
   !---~---
   pure subroutine tend_to_array(k, n, nsl, v)
      type(march_tend_t), intent(in)  :: k
      integer(ik),         intent(in)  :: n, nsl
      real(wp),            intent(out) :: v(:)
      integer(ik) :: i, o
      v(1) = k%d_cas_enthalpy
      v(2) = k%d_cas_shv
      v(3) = k%d_cas_co2
      v(4:3+nsl)       = k%dedt(1:nsl)
      v(4+nsl:3+2*nsl) = k%dtheta_dt(1:nsl)
      o = 3_ik + 2_ik*nsl
      do i = 1_ik, n
         v(o+2*i-1) = k%d_leaf_water_mass(i)
         v(o+2*i) = k%d_wood_water_mass(i)
      end do
      o = o + 2_ik*n
      do i = 1_ik, n
         v(o+2*i-1) = k%d_leaf_surf_water(i)
         v(o+2*i) = k%d_wood_surf_water(i)
      end do
      o = o + 2_ik*n
      v(o+1) = 0.0_wp
      v(o+2) = 0.0_wp
   end subroutine tend_to_array

   !---~---
   !   How entry j of the flat array is treated: norm_group, its group's tolerances (0: none), and
   ! in_embedded_error, whether the march integrates it and so has an error to measure: the
   ! canopy air and the soil heat. The soil and plant water move at the water network's fixed
   ! rate on both schemes, and the films and the pond by operator-split maps, so their embedded
   ! difference is structurally zero: it is set to 0, not left as rounding noise, and the norm
   ! leaves them out.
   !---~---
   pure subroutine state_entry_rules(j, n, nsl, norm_group, in_embedded_error)
      integer(ik), intent(in)  :: j, n, nsl
      integer(ik), intent(out) :: norm_group
      logical,     intent(out) :: in_embedded_error
      integer(ik) :: o
      o = 3_ik + 2_ik*nsl
      if (j == 1_ik) then
         norm_group = GRP_ENTH
         in_embedded_error = .true.
      else if (j == 2_ik) then
         norm_group = GRP_SHV
         in_embedded_error = .true.
      else if (j == 3_ik) then
         norm_group = GRP_CO2
         in_embedded_error = .true.
      else if (j <= 3_ik + nsl) then
         norm_group = GRP_SE
         in_embedded_error = .true.
      else if (j <= o) then
         norm_group = GRP_THETA
         in_embedded_error = .false.
      else if (j <= o + 2_ik*n) then
         norm_group = merge(GRP_LEAF_W, GRP_WOOD_W, mod(j - o, 2_ik) == 1_ik)
         in_embedded_error = .false.
      else
         norm_group = 0_ik
         in_embedded_error = .false.         ! films and pond
      end if
   end subroutine state_entry_rules

   !--- Allocate a state's per-cohort arrays for n cohorts, unless they already are.
   pure subroutine ensure_cohort_arrays(y, n)
      type(march_state_t), intent(inout) :: y
      integer(ik),          intent(in)    :: n
      if (allocated(y%leaf_water_mass) .and. allocated(y%wood_water_mass) .and.                   &
          allocated(y%leaf_surf_water) .and. allocated(y%wood_surf_water)) then
         if (size(y%leaf_water_mass) == n .and. size(y%wood_water_mass) == n .and.               &
             size(y%leaf_surf_water) == n .and. size(y%wood_surf_water) == n) return
      end if
      if (allocated(y%leaf_water_mass)) deallocate(y%leaf_water_mass)
      if (allocated(y%wood_water_mass)) deallocate(y%wood_water_mass)
      if (allocated(y%leaf_surf_water)) deallocate(y%leaf_surf_water)
      if (allocated(y%wood_surf_water)) deallocate(y%wood_surf_water)
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n), y%leaf_surf_water(n), y%wood_surf_water(n))
   end subroutine ensure_cohort_arrays

   !--- The combinators, on the flat array. Each fills its output in place.
   !--- copy: y_stage = y.
   pure subroutine state_init(y, n, nsl, y_stage)
      type(march_state_t), intent(in)    :: y
      integer(ik),          intent(in)    :: n, nsl
      type(march_state_t), intent(inout) :: y_stage
      real(wp) :: v(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y, n, nsl, v)
      call array_to_state(v, n, nsl, y_stage)
   end subroutine state_init

   !---~---
   !   y_stage = y_stage + a*k, in place; a stage of several k's is state_init and repeated
   ! state_accum.
   !---~---
   pure subroutine state_accum(y_stage, a, k, n, nsl)
      type(march_state_t), intent(inout) :: y_stage
      real(wp),             intent(in)    :: a
      type(march_tend_t),  intent(in)    :: k
      integer(ik),          intent(in)    :: n, nsl
      real(wp) :: vy(5_ik + 2_ik*nsl + 4_ik*n), vk(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y_stage, n, nsl, vy)
      call tend_to_array(k, n, nsl, vk)
      vy = vy + a * vk
      call array_to_state(vy, n, nsl, y_stage)
   end subroutine state_accum

   !---~---
   !   Ledger helpers: b-weight two stage rates into amounts over dt (weights (1-gamma, gamma)
   ! times dt); zero an accumulator; add one sub-step's amounts.
   !---~---
   pure subroutine bflux_bweight(acc, s2, s3, dt, gam)
      type(boundary_flux_t), intent(out) :: acc
      type(stage_bflux_t),  intent(in)  :: s2, s3
      real(wp),             intent(in)  :: dt, gam
      real(wp) :: b2, b3
      b2 = (1.0_wp - gam) * dt
      b3 = gam * dt
      acc%cas_enth_in   = b2*s2%cas_enth_in   + b3*s3%cas_enth_in
      acc%cas_enth_out  = b2*s2%cas_enth_out  + b3*s3%cas_enth_out
      acc%cas_vap_in    = b2*s2%cas_vap_in    + b3*s3%cas_vap_in
      acc%cas_vap_out   = b2*s2%cas_vap_out   + b3*s3%cas_vap_out
      acc%cas_co2_in    = b2*s2%cas_co2_in    + b3*s3%cas_co2_in
      acc%cas_co2_out   = b2*s2%cas_co2_out   + b3*s3%cas_co2_out
      acc%soil_enth_in  = b2*s2%soil_enth_in  + b3*s3%soil_enth_in
      acc%soil_enth_out = b2*s2%soil_enth_out + b3*s3%soil_enth_out
      acc%soil_wat_in   = b2*s2%soil_wat_in   + b3*s3%soil_wat_in
      acc%soil_wat_out  = b2*s2%soil_wat_out  + b3*s3%soil_wat_out
      acc%whole_enth_in = b2*s2%whole_enth_in + b3*s3%whole_enth_in
      acc%whole_enth_out= b2*s2%whole_enth_out+ b3*s3%whole_enth_out
      acc%whole_wat_in  = b2*s2%whole_wat_in  + b3*s3%whole_wat_in
      acc%whole_wat_out = b2*s2%whole_wat_out + b3*s3%whole_wat_out
      acc%whole_cond    = b2*s2%whole_cond    + b3*s3%whole_cond
      acc%whole_cond_enth = b2*s2%whole_cond_enth + b3*s3%whole_cond_enth
      acc%snow_enth_in  = b2*s2%snow_enth_in  + b3*s3%snow_enth_in
      acc%snow_vap_out  = b2*s2%snow_vap_out  + b3*s3%snow_vap_out
      acc%veg_surface_evap = b2*s2%veg_surface_evap + b3*s3%veg_surface_evap
      acc%h_canopy      = b2*s2%h_canopy      + b3*s3%h_canopy
      acc%h_ground      = b2*s2%h_ground      + b3*s3%h_ground
      acc%atm_heat_out  = b2*s2%atm_heat_out  + b3*s3%atm_heat_out
      acc%atm_vap_out   = b2*s2%atm_vap_out   + b3*s3%atm_vap_out
      acc%rad_net       = b2*s2%rad_net       + b3*s3%rad_net
   end subroutine bflux_bweight

   pure subroutine bflux_zero(acc, n)
      type(boundary_flux_t), intent(out) :: acc
      integer(ik), optional, intent(in) :: n   ! allocate + zero the per-cohort tissue integrals
      !--- intent(out) zeroes the scalars; no empty structure constructor (docs/building.md traps).
      if (present(n)) then
         allocate(acc%tissue_leaf_int(n), acc%tissue_wood_int(n), acc%transp_coh(n))
         acc%tissue_leaf_int = 0.0_wp
         acc%tissue_wood_int = 0.0_wp
         acc%transp_coh = 0.0_wp
      end if
   end subroutine bflux_zero

   pure subroutine bflux_add(acc, s)
      type(boundary_flux_t), intent(inout) :: acc
      type(boundary_flux_t), intent(in)    :: s
      acc%cas_enth_in   = acc%cas_enth_in   + s%cas_enth_in
      acc%cas_enth_out  = acc%cas_enth_out  + s%cas_enth_out
      acc%cas_vap_in    = acc%cas_vap_in    + s%cas_vap_in
      acc%cas_vap_out   = acc%cas_vap_out   + s%cas_vap_out
      acc%cas_co2_in    = acc%cas_co2_in    + s%cas_co2_in
      acc%cas_co2_out   = acc%cas_co2_out   + s%cas_co2_out
      acc%soil_enth_in  = acc%soil_enth_in  + s%soil_enth_in
      acc%soil_enth_out = acc%soil_enth_out + s%soil_enth_out
      acc%soil_wat_in   = acc%soil_wat_in   + s%soil_wat_in
      acc%soil_wat_out  = acc%soil_wat_out  + s%soil_wat_out
      acc%whole_enth_in = acc%whole_enth_in + s%whole_enth_in
      acc%whole_enth_out= acc%whole_enth_out+ s%whole_enth_out
      acc%whole_wat_in  = acc%whole_wat_in  + s%whole_wat_in
      acc%whole_wat_out = acc%whole_wat_out + s%whole_wat_out
      acc%whole_cond    = acc%whole_cond    + s%whole_cond
      acc%whole_cond_enth = acc%whole_cond_enth + s%whole_cond_enth
      acc%snow_enth_in  = acc%snow_enth_in  + s%snow_enth_in
      acc%snow_vap_out  = acc%snow_vap_out  + s%snow_vap_out
      acc%veg_surface_evap = acc%veg_surface_evap + s%veg_surface_evap
      acc%h_canopy      = acc%h_canopy      + s%h_canopy
      acc%h_ground      = acc%h_ground      + s%h_ground
      acc%atm_heat_out  = acc%atm_heat_out  + s%atm_heat_out
      acc%atm_vap_out   = acc%atm_vap_out   + s%atm_vap_out
      acc%rad_net       = acc%rad_net       + s%rad_net
      !---~---
      !   Only accepted sub-steps reach here, so the tissue integrals cover the same sub-steps as
      ! every other amount.
      !---~---
      if (allocated(acc%tissue_leaf_int) .and. allocated(s%tissue_leaf_int)) then
         acc%tissue_leaf_int = acc%tissue_leaf_int + s%tissue_leaf_int
         acc%tissue_wood_int = acc%tissue_wood_int + s%tissue_wood_int
      end if
      if (allocated(acc%transp_coh) .and. allocated(s%transp_coh)) acc%transp_coh = acc%transp_coh + s%transp_coh
   end subroutine bflux_add

   !---~---
   !   out = (1-b)*y + b*Y2, the ARS stage-3 extrapolation base (b = (1-gamma)/gamma,
   ! docs/science/numerical_scheme.md eq. (1)). The plant water, the films and the pond are
   ! frozen through the stages, so for them this is y again (up to rounding).
   !---~---
   pure subroutine state_extrap(y, b, Y2, n, nsl, out)
      type(march_state_t), intent(in)    :: y, Y2
      real(wp),             intent(in)    :: b
      integer(ik),          intent(in)    :: n, nsl
      type(march_state_t), intent(inout) :: out
      real(wp) :: vy(5_ik + 2_ik*nsl + 4_ik*n), v2(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y, n, nsl, vy)
      call state_to_array(Y2, n, nsl, v2)
      vy = (1.0_wp - b)*vy + b*v2
      call array_to_state(vy, n, nsl, out)
   end subroutine state_extrap

   !---~---
   !   Clamp the extrapolated CAS enthalpy and humidity into a wide physical range, so an overshoot
   ! of the ARS base (b = 1 + sqrt(2)) cannot reach a temperature where qsat overflows. A state
   ! in range is left as it is, bit for bit; one out of range keeps its temperature (clamped)
   ! at its clamped humidity. nfire counts the calls that moved the state (WORK_CLAMP_STAGE).
   !---~---
   pure subroutine clamp_cas(s, nfire)
      type(march_state_t), intent(inout) :: s
      integer(ik), optional, intent(inout) :: nfire
      real(wp) :: t, shv_c
      real(wp), parameter :: SHV_LO = 1.0e-8_wp, SHV_HI = 0.06_wp
      shv_c = min(max(s%cas_shv, SHV_LO), SHV_HI)
      t     = cas_enthalpy_to_temp(s%cas_enthalpy, shv_c)
      if (shv_c == s%cas_shv .and. t >= STATE_TEMP_MIN .and. t <= STATE_TEMP_MAX) return
      t     = min(max(t, STATE_TEMP_MIN), STATE_TEMP_MAX)
      s%cas_shv      = shv_c
      s%cas_enthalpy = cas_temp_to_enthalpy(t, shv_c)
      if (present(nfire)) nfire = nfire + 1_ik
   end subroutine clamp_cas

   !---~---
   !   Is a temperature at or beyond a state bound? Within a microkelvin counts as at it, since a
   ! clamped canopy air's temperature comes back from its enthalpy a few bits off the bound.
   !---~---
   elemental logical function at_temp_bound(t)
      real(wp), intent(in) :: t
      real(wp), parameter :: AT_BOUND = 1.0e-6_wp   ! [K]
      at_temp_bound = t <= STATE_TEMP_MIN + AT_BOUND .or. t >= STATE_TEMP_MAX - AT_BOUND
   end function at_temp_bound

   !---~---
   !   The embedded error estimate of the ARS(2,2,2) step, err = (Y3 - base3) - (Y2 - y)
   ! (docs/science/numerical_scheme.md eq. (2)), for the entries state_entry_rules counts, and
   ! 0 for the rest.
   !---~---
   pure subroutine state_err_diff(Y3, base3, Y2, y, n, nsl, err)
      type(march_state_t), intent(in)    :: Y3, base3, Y2, y
      integer(ik),          intent(in)    :: n, nsl
      type(march_state_t), intent(inout) :: err
      real(wp)    :: v3(5_ik + 2_ik*nsl + 4_ik*n), vb(5_ik + 2_ik*nsl + 4_ik*n)
      real(wp)    :: v2(5_ik + 2_ik*nsl + 4_ik*n), vy(5_ik + 2_ik*nsl + 4_ik*n)
      integer(ik) :: j, grp
      logical     :: counted
      call state_to_array(Y3, n, nsl, v3)
      call state_to_array(base3, n, nsl, vb)
      call state_to_array(Y2, n, nsl, v2)
      call state_to_array(y, n, nsl, vy)
      do j = 1_ik, size(v3, kind=ik)
         call state_entry_rules(j, n, nsl, grp, counted)
         if (counted) then
            v3(j) = (v3(j) - vb(j)) - (v2(j) - vy(j))
         else
            v3(j) = 0.0_wp
         end if
      end do
      call array_to_state(v3, n, nsl, err)
   end subroutine state_err_diff

   !--- out = a - b, the state difference that forms the low-order embedded solution.
   pure subroutine state_sub(a, b, n, nsl, out)
      type(march_state_t), intent(in)    :: a, b
      integer(ik),          intent(in)    :: n, nsl
      type(march_state_t), intent(inout) :: out
      real(wp) :: va(5_ik + 2_ik*nsl + 4_ik*n), vb(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(a, n, nsl, va)
      call state_to_array(b, n, nsl, vb)
      va = va - vb
      call array_to_state(va, n, nsl, out)
   end subroutine state_sub

   !---~---
   !   A state for n cohorts with every field zero, so the step-size norm can read an error that
   ! is already a difference (state_wrms_grouped takes two states).
   !---~---
   pure function zero_like(ref, n, nsl) result(z)
      type(march_state_t), intent(in) :: ref
      integer(ik),          intent(in) :: n, nsl
      type(march_state_t) :: z
      real(wp) :: v(5_ik + 2_ik*nsl + 4_ik*n)
      v = 0.0_wp
      call array_to_state(v, n, nsl, z)
   end function zero_like

   !--- Re-diagnose every layer's temperature and liquid fraction from the committed energy + water.
   pure subroutine diagnose_soil_temps(y_out, profile, nsl, soil_temp, soil_fliq)
      type(march_state_t), intent(in)    :: y_out
      type(soil_profile_t), intent(in)    :: profile
      integer(ik),          intent(in)    :: nsl
      real(wp),             intent(inout) :: soil_temp(:), soil_fliq(:)
      integer(ik) :: k
      do k = 1_ik, nsl
         call soil_layer_temp(y_out%soil_energy(k), y_out%theta(k), profile%dry_heat_capacity(k),             &
                              profile%layer_curve(k), soil_temp(k), soil_fliq(k))
      end do
   end subroutine diagnose_soil_temps

end module ed_fast_integ_types
