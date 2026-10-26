! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_diag_vars -- the DIAGNOSTIC accumulators, per cohort, per patch and per polygon.      !
!                                                                                          !
! Stage [2] of the diagnostic wall, CAPTURE (docs/science/diagnostics.md §1). The problem  !
! it solves: the fast loop computes ~40 per-cohort and per-patch quantities every dt_fast  !
! -- leaf assimilation, stomatal conductance, tissue water potentials, sapflow, absorbed   !
! radiation, turbulent fluxes, runoff, budget residuals -- and sub-daily resolution exists !
! only inside the fast loop's sub-step, so without a home they would be recomputed and     !
! discarded ~48 times a day. These blocks are where they land instead.                     !
!                                                                                          !
! SIX DESIGN PRINCIPLES:                                                                    !
!                                                                                          !
!  1. dt-WEIGHTED SUMS, normalized on read. Every field accumulates x*dt and `w`           !
!     accumulates dt; the reader divides. That is correct under the adaptive integrator    !
!     (unequal sub-steps) and under AGG_TMEAN chaining for free, where a plain running     !
!     mean would not be.                                                                   !
!                                                                                          !
!  2. RESET ONCE PER SLOW STEP -- the same lifecycle as gpp_accum / et_accum / decomp_xi.  !
!     No second reset cadence is invented.                                                 !
!                                                                                          !
!  3. TRANSIENT: never written to the restart file. A checkpoint is prognostic state at an !
!     instant; a time-averaged diagnostic is not state.                                    !
!                                                                                          !
!  4. THEY STILL RIDE THE LOCKSTEP. The reset is per slow step, but the slow step re-sorts !
!     the cohorts after the fast loop has filled these and before the output tick reads    !
!     them, and the calendar restructuring (fuse / split / cull / recruit / disturb)       !
!     permutes the axis at a boundary. So slot i must keep meaning cohort i across every   !
!     permutation. That obligation is REAL and is the main cost of this design.            !
!                                                                                          !
!     It is discharged by STORAGE LAYOUT rather than by discipline: the fields are rows of !
!     ONE 2-D array v(N_CO_FAST_ROWS, cap), so every lockstep operation is a single whole-array   !
!     statement that cannot omit a field. Adding a diagnostic is a new index parameter     !
!     plus a fill line -- ZERO edits to the reorder machinery. The alternative, one        !
!     allocatable per field, would need a new line per field in each of alloc_cohort_block /     !
!     expand_cohort_block / move_cohort_block / reorder_cohort_block / copy_cohort --               !
!     precisely the "forgot to reallocate" class CLAUDE.md names as the reason the         !
!     centralized reorder exists.                                                          !
!                                                                                          !
!  5. ALLOCATED ONLY ON DEMAND. A run that reports no per-cohort ecophysiology pays        !
!     nothing: main sets `active` from the registry, and every entry point is a no-op      !
!     while it is .false.                                                                  !
!                                                                                          !
!  6. FUSION BLENDS BY THE FIELD'S OWN KIND. Each field declares the weight that carries   !
!     it: leaf area (INTENSIVE -- a conductance, a potential, a per-leaf-area flux), plant !
!     number (EXTENSIVE -- a per-plant flux), the tissue's heat capacity (a leaf or wood   !
!     temperature, as its prognostic twin fuses), or none (GROUND-REFERENCED: summed,      !
!     already per m2 of patch). Blending a temperature by nplant, or an extensive flux by  !
!     leaf area, is the classic silent fusion bug; the kind is DATA here, so the blend     !
!     cannot disagree with the quantity.                                                   !
!==========================================================================================!
module ed_diag_vars
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : tiny_num
   use ed_numerics,    only : pair_mean
   implicit none
   private

   public :: cohort_diag_block, patch_diag_block, polygon_diag_block, site_diag_block
   public :: reset_polygon_diag, polygon_diag_value, alloc_site_diag, reset_site_diag, site_diag_value
   public :: alloc_cohort_diag, free_cohort_diag, expand_cohort_diag, reset_cohort_diag
   public :: reorder_cohort_diag, copy_cohort_diag, fuse_cohort_diag, clear_cohort_diag
   public :: alloc_patch_diag, free_patch_diag, expand_patch_diag, reset_patch_diag
   public :: reorder_patch_diag, copy_patch_diag, blend_patch_diag, clear_patch_diag
   public :: inherit_patch_diag
   public :: cohort_diag_value, patch_diag_value, patch_diag_slow_row

   !----- Fusion kinds (principle 6 above). ---------------------------------------------------!
   integer(ik), parameter, public :: FSK_INTENSIVE = 1_ik  !< leaf-area-weighted mean (gs, psi, per-leaf-area fluxes)
   integer(ik), parameter, public :: FSK_EXTENSIVE = 2_ik  !< nplant-weighted mean (per-plant fluxes)
   integer(ik), parameter, public :: FSK_GROUND    = 3_ik  !< summed (already per m2 of patch ground)
   integer(ik), parameter, public :: FSK_LEAF_HEAT = 4_ik  !< leaf-heat-capacity-weighted mean (leaf temperature)
   integer(ik), parameter, public :: FSK_WOOD_HEAT = 5_ik  !< wood-heat-capacity-weighted mean (wood temperature)

   !==========================================================================================!
   !  PER-COHORT diagnostic rows, numbered once: the fast loop's rows are 1..N_CO_FAST_ROWS    !
   !  (cohort%diag), the slow step's follow up to N_CO_ROWS (cohort%sdiag). Each block holds   !
   !  only its own rows, so a fast row used on the slow block is out of bounds, which a Debug  !
   !  build stops on. Order is free; the parameters are what bind a row to a meaning.          !
   !==========================================================================================!
   !----- 1. Carbon -----------------------------------------------------------------------------!
   integer(ik), parameter, public :: CO_ANET         =  1_ik  !< [umol CO2/m2 leaf/s] net assimilation
   integer(ik), parameter, public :: CO_AGROSS       =  2_ik  !< [umol CO2/m2 leaf/s] gross assimilation
   integer(ik), parameter, public :: CO_RD           =  3_ik  !< [umol/m2 leaf/s] leaf dark respiration
   integer(ik), parameter, public :: CO_GPP_RATE     =  4_ik  !< [umol CO2/plant/s]
   !----- 2. Stomata and leaf exchange ----------------------------------------------------------!
   integer(ik), parameter, public :: CO_GSW          =  5_ik  !< [mol H2O/m2 leaf/s]  stomatal conductance
   integer(ik), parameter, public :: CO_GBW          =  6_ik  !< [m/s] leaf boundary-layer conductance
   integer(ik), parameter, public :: CO_CI           =  7_ik  !< [umol/mol] intercellular CO2
   integer(ik), parameter, public :: CO_CS           =  8_ik  !< [umol/mol] leaf-surface CO2
   integer(ik), parameter, public :: CO_TRANSP       =  9_ik  !< [mol H2O/m2 leaf/s] leaf transpiration
   integer(ik), parameter, public :: CO_BETA_STOM    = 10_ik  !< [-] stomatal water-stress limb
   integer(ik), parameter, public :: CO_BETA_NONSTOM = 11_ik  !< [-] non-stomatal (capacity) limb
   !----- 3. Plant water ------------------------------------------------------------------------!
   integer(ik), parameter, public :: CO_PSI_LEAF     = 12_ik  !< [MPa]
   integer(ik), parameter, public :: CO_PSI_WOOD     = 13_ik  !< [MPa]
   integer(ik), parameter, public :: CO_PLC          = 14_ik  !< [-] percent loss of conductance
   integer(ik), parameter, public :: CO_SAPFLOW      = 15_ik  !< [kg/plant/s] wood -> leaf
   integer(ik), parameter, public :: CO_ROOT_UPTAKE  = 16_ik  !< [kg/plant/s] soil -> root
   integer(ik), parameter, public :: CO_LEAF_WATER   = 17_ik  !< [kg/plant] internal leaf water
   integer(ik), parameter, public :: CO_WOOD_WATER   = 18_ik  !< [kg/plant] internal wood water
   !----- 4. Energy and radiation ---------------------------------------------------------------!
   integer(ik), parameter, public :: CO_LEAF_TEMP    = 19_ik  !< [K]
   integer(ik), parameter, public :: CO_WOOD_TEMP    = 20_ik  !< [K]
   integer(ik), parameter, public :: CO_ABS_PAR      = 21_ik  !< [W/m2 ground] absorbed PAR
   integer(ik), parameter, public :: CO_ABS_SW       = 22_ik  !< [W/m2 ground] absorbed shortwave
   integer(ik), parameter, public :: CO_ABS_LW       = 23_ik  !< [W/m2 ground] net longwave
   !----- 5. Canopy environment -----------------------------------------------------------------!
   integer(ik), parameter, public :: CO_LEAF_VPD     = 24_ik  !< [Pa] leaf-to-air VPD
   integer(ik), parameter, public :: CO_WIND         = 25_ik  !< [m/s] in-canopy wind
   !----- 6. The leaf solve's other inputs ------------------------------------------------------!
   !----- The rest of the leaf solve's own inputs, as it saw them, so the solve can be repeated     !
   !      outside the model (the calibration's canopy of leaf kernels): with CO_LEAF_TEMP,           !
   !      CO_LEAF_VPD and CO_PSI_LEAF they are every driver of leaf_gas_exchange_batch. ------------!
   integer(ik), parameter, public :: CO_LEAF_PAR     = 26_ik  !< [umol photon/m2 leaf/s] incident-equivalent PAR per leaf area
   integer(ik), parameter, public :: CO_GB_MOL       = 27_ik  !< [mol H2O/m2 leaf/s] boundary-layer conductance
   integer(ik), parameter, public :: CO_PSI_PREDAWN  = 28_ik  !< [MPa] the stomatal limb's predawn leaf potential
   integer(ik), parameter, public :: CO_CA           = 29_ik  !< [umol/mol] canopy-air CO2 the leaf sees
   integer(ik), parameter, public :: CO_PRESSURE     = 30_ik  !< [Pa] air pressure the leaf sees
   integer(ik), parameter, public :: N_CO_FAST_ROWS  = 30_ik

   !==========================================================================================!
   !  PER-COHORT SLOW-STEP diagnostic rows, N_CO_FAST_ROWS+1..N_CO_ROWS -- a SECOND block of the   !
   !  same shape, written once per slow step by the vegetation-dynamics driver.                     !
   !                                                                                          !
   !  WHY A SEPARATE BLOCK AND NOT MORE ROWS ABOVE. The two have different sample semantics: the     !
   !  fast rows accumulate ~48 sub-step samples per slow step, the slow rows get exactly one. One     !
   !  shared dt weight cannot normalize both.                                                         !
   !                                                                                          !
   !  WHY NOT READ poly%deriv DIRECTLY. `cohort_deriv_block` is TRANSIENT and deliberately NOT        !
   !  lockstep-reordered -- which is fine for its own consumer, since update_cohort_states applies it  !
   !  immediately. But the step re-sorts the cohorts after that, so at the output tick `deriv(i)` and  !
   !  `cohort(i)` need not be the same plant, and an error of that kind shows only when the ordering   !
   !  changes, as thread count can make it. Storing the tendencies in a block that DOES ride the       !
   !  lockstep removes the failure mode rather than timing around it.                                  !
   !==========================================================================================!
   !----- 7. Growth and mortality ---------------------------------------------------------------!
   integer(ik), parameter, public :: CO_DDBH_DT      = 31_ik  !< [cm/yr]         diameter growth
   integer(ik), parameter, public :: CO_DAGB_DT      = 32_ik  !< [kgC/plant/yr]  AGB growth
   integer(ik), parameter, public :: CO_MORT_RATE    = 33_ik  !< [1/yr]          mortality hazard (positive)
   !----- 8. Carbon allocation ------------------------------------------------------------------!
   integer(ik), parameter, public :: CO_NPP_LEAF     = 34_ik  !< [kgC/plant/yr]  NPP to leaf
   integer(ik), parameter, public :: CO_NPP_FINEROOT = 35_ik  !< [kgC/plant/yr]  NPP to fine root
   integer(ik), parameter, public :: CO_NPP_WOOD     = 36_ik  !< [kgC/plant/yr]  NPP to wood
   integer(ik), parameter, public :: CO_NPP_STORAGE  = 37_ik  !< [kgC/plant/yr]  NPP to non-structural storage
   integer(ik), parameter, public :: CO_NPP_REPRO    = 38_ik  !< [kgC/plant/yr]  NPP to reproduction
   !----- 9. Carbon losses ----------------------------------------------------------------------!
   integer(ik), parameter, public :: CO_GROWTH_RESP  = 39_ik  !< [kgC/plant/yr]  growth respiration
   integer(ik), parameter, public :: CO_EXUDATE      = 40_ik  !< [kgC/plant/yr]  root exudate (carbon the growth sink could not use)
   integer(ik), parameter, public :: N_CO_ROWS       = 40_ik

   !----- Fusion kind of each cohort row. The leaf and wood temperatures fuse on the heat capacity   !
   !      of their tissue, as their prognostic twins do (fuse_cohort_fast_state), so a fused        !
   !      cohort's mean temperature over the output window is the one its tissue energy implies.    !
   !      Conductances, potentials and the leaf-area-normalized fluxes are INTENSIVE (leaf-area-    !
   !      weighted). The per-PLANT quantities (GPP rate, sapflow, uptake, tissue water mass, and    !
   !      every slow row) are EXTENSIVE (nplant-weighted), matching how their prognostic twins      !
   !      fuse. The absorbed-radiation terms are per m2 of GROUND, so two cohorts' contributions    !
   !      to the same patch simply add -- the same convention leaf_surf_water follows. -------------!
   integer(ik), parameter, public :: CO_FUSION_KIND(N_CO_ROWS) = [                                               &
        ! 1. Carbon
        FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_EXTENSIVE,                                              &
        ! 2. Stomata and leaf exchange
        FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, &
        ! 3. Plant water
        FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE, &
        ! 4. Energy and radiation
        FSK_LEAF_HEAT, FSK_WOOD_HEAT, FSK_GROUND, FSK_GROUND, FSK_GROUND,                                        &
        ! 5. Canopy environment
        FSK_INTENSIVE, FSK_INTENSIVE,                                                                            &
        ! 6. The leaf solve's other inputs
        FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE, FSK_INTENSIVE,                               &
        ! 7. Growth and mortality
        FSK_EXTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE,                                                             &
        ! 8. Carbon allocation
        FSK_EXTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE, FSK_EXTENSIVE,                               &
        ! 9. Carbon losses
        FSK_EXTENSIVE, FSK_EXTENSIVE                                                                             ]

   !==========================================================================================!
   !  PER-PATCH diagnostic fields. The fast loop fills the fast rows once per (patch, sub-step)  !
   !  (patch_diag_row): the coarse tiers accumulate that row, and the FAST tier stages it as it  !
   !  is, so the two read one list.                                                             !
   !==========================================================================================!
   !----- 1. Forcing at the canopy-air top ------------------------------------------------------!
   !----- The patch's own forcing: the sample moved to its canopy-air top (docs/science/forcing.md  !
   !      §8). Only the wind and the air temperature differ between patches.                       !
   integer(ik), parameter, public :: PA_WIND_CAS_TOP      =  1_ik  !< [m/s] wind speed at the canopy-air top
   integer(ik), parameter, public :: PA_TAIR_CAS_TOP      =  2_ik  !< [K]   air temperature at the canopy-air top
   integer(ik), parameter, public :: PA_Z_CAS_TOP         =  3_ik  !< [m]   that top: the canopy-air depth the fast loop used
   !----- 2. Top-of-canopy radiation ------------------------------------------------------------!
   !----- TOP-OF-CANOPY radiative fluxes. Kept as FLUXES, not as a time-averaged albedo: a        !
   !      period-mean albedo is the mean of a RATIO, which is not the ratio of the means, and     !
   !      at night the shortwave ratio is 0/0. The albedo a reader wants is sum(up)/sum(down)     !
   !      over the period (docs/science/canopy_radiation_transfer.md eq. (16)), which these       !
   !      make computable -- and which is what a satellite product is. ---------------------------!
   integer(ik), parameter, public :: PA_SW_IN_VIS         =  4_ik  !< [W/m2] incident VIS at canopy top
   integer(ik), parameter, public :: PA_SW_IN_NIR         =  5_ik  !< [W/m2] incident NIR at canopy top
   integer(ik), parameter, public :: PA_SW_UP_VIS         =  6_ik  !< [W/m2] upwelling VIS leaving the canopy top
   integer(ik), parameter, public :: PA_SW_UP_NIR         =  7_ik  !< [W/m2] upwelling NIR leaving the canopy top
   integer(ik), parameter, public :: PA_SW_UP             =  8_ik  !< [W/m2] upwelling shortwave, both bands (VIS + NIR)
   integer(ik), parameter, public :: PA_LW_UP             =  9_ik  !< [W/m2] upwelling LW (surface emission included)
   !----- The SKIN temperature, in the land-model sense (CLM's TSKIN): the black-body             !
   !      temperature of the longwave leaving the canopy top, what an infrared thermometer or     !
   !      a satellite land-surface temperature sees (docs/science/diagnostics.md §4, "The         !
   !      skin temperature"). ED-MORE has no skin layer of its own -- snow-free, the ground       !
   !      surface is the top soil layer (PA_SOIL_TEMP_TOP) -- and over a forest this is           !
   !      mostly the canopy. It is a fourth root, so it is formed per sub-step, like the VPD. ----!
   integer(ik), parameter, public :: PA_SKIN_TEMP         = 10_ik  !< [K] (longwave up / sigma)^(1/4)
   !----- 3. Turbulence -------------------------------------------------------------------------!
   integer(ik), parameter, public :: PA_USTAR             = 11_ik  !< [m/s]
   integer(ik), parameter, public :: PA_GGNET             = 12_ik  !< [m/s] ground conductance
   integer(ik), parameter, public :: PA_ROUGH             = 13_ik  !< [m] roughness length
   integer(ik), parameter, public :: PA_DISPLACE          = 14_ik  !< [m] displacement height
   !----- 4. Surface energy balance -------------------------------------------------------------!
   integer(ik), parameter, public :: PA_RNET              = 15_ik  !< [W/m2] net all-wave radiation
   integer(ik), parameter, public :: PA_H                 = 16_ik  !< [W/m2] sensible heat (CAS -> atm)
   integer(ik), parameter, public :: PA_LE                = 17_ik  !< [W/m2] latent heat (CAS -> atm)
   integer(ik), parameter, public :: PA_SW_GROUND         = 18_ik  !< [W/m2] shortwave reaching the ground
   integer(ik), parameter, public :: PA_LW_GROUND         = 19_ik  !< [W/m2] net longwave at the ground
   !----- 5. Carbon and water fluxes ------------------------------------------------------------!
   integer(ik), parameter, public :: PA_GPP               = 20_ik  !< [umol CO2/m2/s] patch GPP rate
   integer(ik), parameter, public :: PA_NPP               = 21_ik  !< [umol CO2/m2/s] GPP net of maintenance respiration
   integer(ik), parameter, public :: PA_RECO              = 22_ik  !< [umol CO2/m2/s] ecosystem respiration, NEE + GPP
   integer(ik), parameter, public :: PA_NEE               = 23_ik  !< [umol CO2/m2/s] net ecosystem exchange
   integer(ik), parameter, public :: PA_TRANSP            = 24_ik  !< [kg/m2/s] canopy transpiration
   !----- 6. Canopy air and ground state --------------------------------------------------------!
   !----- Fast-varying CAS / surface state, accumulated dt-weighted rather than read at the      !
   !      output tick. A state read at the tick is ONE instantaneous sample per output           !
   !      window, taken at whatever local time the boundary falls on, so its bias is a           !
   !      function of the site's longitude -- measured at +1.4 K for the top soil layer at       !
   !      Ithaca. -------------------------------------------------------------------------------!
   integer(ik), parameter, public :: PA_CAS_TEMP          = 25_ik  !< [K]
   integer(ik), parameter, public :: PA_CAS_SHV           = 26_ik  !< [kg/kg]
   integer(ik), parameter, public :: PA_CAS_CO2           = 27_ik  !< [umol/mol]
   integer(ik), parameter, public :: PA_CAS_VPD           = 28_ik  !< [Pa]
   integer(ik), parameter, public :: PA_SOIL_TEMP_TOP     = 29_ik  !< [K] top soil layer (the ground surface when snow-free)
   integer(ik), parameter, public :: PA_W_SURFACE         = 30_ik  !< [kg/m2] ponded surface water
   !----- 7. Budget residuals -------------------------------------------------------------------!
   integer(ik), parameter, public :: PA_RESID_ENERGY      = 31_ik  !< [W/m2] whole-column energy residual
   integer(ik), parameter, public :: PA_RESID_WATER       = 32_ik  !< [kg/m2/s] whole-column water residual
   integer(ik), parameter, public :: N_PA_FAST_ROWS         = 32_ik  !< the fast loop's last row; the slow rows follow
   !----- 8. SLOW rows: litter ------------------------------------------------------------------!
   !----- SLOW-loop patch diagnostics. They share this block (and its weight) with the fast rows    !
   !      above because the reader normalizes by w either way, and the slow driver writes them as   !
   !      a dt-weighted contribution over the same slow step the fast rows span. -------------------!
   integer(ik), parameter, public :: PA_LITTER_LEAF       = 33_ik  !< [kgC/m2/yr] leaf litterfall (after resorption)
   integer(ik), parameter, public :: PA_LITTER_FINEROOT   = 34_ik  !< [kgC/m2/yr] fine-root litter (after resorption)
   !----- Litter by the tissue it came from, completing PA_LITTER_LEAF / PA_LITTER_FINEROOT:     !
   !      with mean_root_exudate_py they sum to the carbon entering the soil pools                  !
   !      (docs/science/soil_carbon.md §5). -----------------------------------------------------!
   integer(ik), parameter, public :: PA_LITTER_WOOD       = 35_ik  !< [kgC/m2/yr] wood + storage of dead plants
   integer(ik), parameter, public :: PA_LITTER_REPRO      = 36_ik  !< [kgC/m2/yr] fruit and seed that made no recruit
   integer(ik), parameter, public :: PA_LITTER_STRUCT     = 37_ik  !< [kgC/m2/yr] structural litter + CWD
   !----- 9. SLOW rows: demography --------------------------------------------------------------!
   integer(ik), parameter, public :: PA_RECRUIT_NPLANT    = 38_ik  !< [plant/m2/yr] recruitment density flux
   integer(ik), parameter, public :: PA_DISTURB_AREA      = 39_ik  !< [1/yr] area fraction disturbed
   !----- MORTALITY CARBON BY PATHWAY: the three ways an ED-MORE plant can die, kept apart so a   !
   !      user can tell a stand thinning continuously from one being knocked over. The three      !
   !      are mutually exclusive and jointly exhaustive -- every kilogram of carbon that          !
   !      leaves the live pools as a whole-individual death passes through exactly one of         !
   !      them -- so they SUM to the total, which is the property test_slow_step                  !
   !      asserts.                                                                                !
   !                                                                                              !
   !      Written UNCONDITIONALLY, unlike the litter each pathway also generates: two of the      !
   !      three litter sites sit behind [soil_carbon].soil_carbon_on, and how much biomass        !
   !      died is a DEMOGRAPHIC question that does not stop being asked when the soil pools       !
   !      are switched off. Whole-individual death carries every pool (leaf + fine root +         !
   !      wood + storage), as the litter partitioning does. --------------------------------------!
   integer(ik), parameter, public :: PA_MORT_C_BACKGROUND = 40_ik  !< [kgC/m2/yr] continuous hazard mortality
   integer(ik), parameter, public :: PA_MORT_C_CULL       = 41_ik  !< [kgC/m2/yr] cohorts culled below the tracking floor
   integer(ik), parameter, public :: PA_MORT_C_DISTURB    = 42_ik  !< [kgC/m2/yr] canopy killed by patch disturbance
   integer(ik), parameter, public :: N_PA_ROWS              = 42_ik

   !==========================================================================================!
   !  PER-POLYGON diagnostic fields: quantities that are the same everywhere in the polygon,  !
   !  so they are accumulated ONCE, not per patch. The prefixes: CD_ cohort (fast loop), CS_  !
   !  cohort (slow loop), PD_ patch, PY_ polygon. The forcing echo (the `forcing` group,      !
   !  docs/science/diagnostics.md §5) is the forcing the fast loop used, after the reader's   !
   !  shortwave partition, rain/snow split and optional height and lapse corrections, dt-     !
   !  weighted over the sub-steps.                                                            !
   !==========================================================================================!
   !----- 1. Radiation --------------------------------------------------------------------------!
   integer(ik), parameter, public :: PY_SW_IN       =  1_ik  !< [W/m2] incident shortwave at canopy top
   integer(ik), parameter, public :: PY_PAR_BEAM    =  2_ik  !< [W/m2] direct-beam PAR
   integer(ik), parameter, public :: PY_PAR_DIFFUSE =  3_ik  !< [W/m2] diffuse PAR
   integer(ik), parameter, public :: PY_NIR_BEAM    =  4_ik  !< [W/m2] direct-beam NIR
   integer(ik), parameter, public :: PY_NIR_DIFFUSE =  5_ik  !< [W/m2] diffuse NIR
   integer(ik), parameter, public :: PY_LWDOWN      =  6_ik  !< [W/m2] downward longwave
   integer(ik), parameter, public :: PY_COSZ        =  7_ik  !< [-] cosine of the solar zenith angle (0 with the sun down)
   !----- 2. Precipitation ----------------------------------------------------------------------!
   integer(ik), parameter, public :: PY_PRECIP      =  8_ik  !< [kg/m2/s] total precipitation (rain + snow)
   integer(ik), parameter, public :: PY_RAINF       =  9_ik  !< [kg/m2/s] liquid precipitation
   integer(ik), parameter, public :: PY_SNOWFALL    = 10_ik  !< [kg/m2/s] frozen precipitation
   !----- 3. Air --------------------------------------------------------------------------------!
   integer(ik), parameter, public :: PY_TAIR        = 11_ik  !< [K] reference-height air temperature
   integer(ik), parameter, public :: PY_QAIR        = 12_ik  !< [kg/kg] specific humidity
   integer(ik), parameter, public :: PY_PSURF       = 13_ik  !< [Pa] surface pressure
   integer(ik), parameter, public :: PY_WIND        = 14_ik  !< [m/s] wind speed at the reference height
   integer(ik), parameter, public :: PY_CO2         = 15_ik  !< [umol/mol] free-atmosphere CO2
   integer(ik), parameter, public :: PY_RHO_AIR     = 16_ik  !< [kg/m3] air density
   integer(ik), parameter, public :: N_PY_ROWS       = 16_ik
   !----- SI_: the forcing each site receives, the polygon's after the site's terrain              !
   !      (docs/science/forcing.md §8). -----------------------------------------------------------!
   integer(ik), parameter, public :: SI_SW_IN    = 1_ik  !< [W/m2] incident shortwave at canopy top
   integer(ik), parameter, public :: SI_PAR_BEAM = 2_ik  !< [W/m2] direct-beam PAR on the slope
   integer(ik), parameter, public :: SI_NIR_BEAM = 3_ik  !< [W/m2] direct-beam NIR on the slope
   integer(ik), parameter, public :: SI_LWDOWN   = 4_ik  !< [W/m2] downward longwave
   integer(ik), parameter, public :: SI_TAIR     = 5_ik  !< [K] reference-height air temperature
   integer(ik), parameter, public :: SI_QAIR     = 6_ik  !< [kg/kg] specific humidity
   integer(ik), parameter, public :: SI_PSURF    = 7_ik  !< [Pa] surface pressure
   integer(ik), parameter, public :: SI_RHO_AIR  = 8_ik  !< [kg/m3] air density
   integer(ik), parameter, public :: N_SI_ROWS   = 8_ik

   !==========================================================================================!
   !  The blocks themselves. `v` is (field, slot): field-major so a lockstep permutation of the  !
   !  SLOT axis is one contiguous whole-array statement (principle 4). Extraction reads one field  !
   !  across slots, which is strided -- acceptable, because this is the once-per-step diagnostic   !
   !  path, not the inner integrator.                                                              !
   !==========================================================================================!
   type :: cohort_diag_block
      integer(ik) :: n = 0_ik, cap = 0_ik
      logical     :: active = .false.        !< .false. => every entry point is a no-op, nothing allocated
      real(wp), allocatable :: v(:,:)        !< (rows, cap) running Sum(x*dt): the fast rows, or the slow rows
      real(wp), allocatable :: w(:)          !< (cap)         running Sum(dt)
   end type cohort_diag_block

   type :: patch_diag_block
      integer(ik) :: n = 0_ik, cap = 0_ik
      logical     :: active = .false.
      !----- Whether the fast loop fills the block. A slow-only run fills only the slow operators'  !
      !      rows (patch_diag_slow_row); the rest are not simulated and read as missing. -----------!
      logical     :: fast_rows = .true.
      real(wp), allocatable :: v(:,:)        !< (N_PA_ROWS, cap) running Sum(x*dt)
      real(wp), allocatable :: w(:)          !< (cap)          running Sum(dt)
   end type patch_diag_block

   !----- One slot: the polygon. Fixed size, so nothing to allocate or blend. --------------------!
   type :: polygon_diag_block
      logical  :: active = .false.
      real(wp) :: v(N_PY_ROWS) = 0.0_wp      !< running Sum(x*dt)
      real(wp) :: w = 0.0_wp                !< running Sum(dt)
   end type polygon_diag_block

   !----- One slot per site. The sites never change, so nothing to reorder or blend. ------------!
   type :: site_diag_block
      logical  :: active = .false.
      real(wp), allocatable :: v(:,:)        !< (N_SI_ROWS, site) running Sum(x*dt)
      real(wp), allocatable :: w(:)          !< (site)            running Sum(dt)
   end type site_diag_block

contains

   !=======================================================================================!
   !  COHORT block lifecycle. Every routine no-ops while `active` is .false., so a run that   !
   !  reports no per-cohort ecophysiology never allocates and never pays.                     !
   !=======================================================================================!

   subroutine alloc_cohort_diag(d, cap, active, slow_rows)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: cap
      logical,                 intent(in)    :: active
      logical, optional,       intent(in)    :: slow_rows   !< the slow step's rows (default: the fast loop's)
      integer(ik) :: first, last
      first = 1_ik
      last  = N_CO_FAST_ROWS
      if (present(slow_rows)) then
         if (slow_rows) then
            first = N_CO_FAST_ROWS + 1_ik
            last  = N_CO_ROWS
         end if
      end if
      d%active = active
      if (.not. active) return
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      allocate(d%v(first:last, max(cap,1_ik)), d%w(max(cap,1_ik)))
      d%v = 0.0_wp
      d%w = 0.0_wp          ! the first window starts empty: nothing zeroes it later
      d%cap = max(cap, 1_ik)
      d%n = 0_ik
      call reset_cohort_diag(d)
   end subroutine alloc_cohort_diag

   subroutine free_cohort_diag(d)
      type(cohort_diag_block), intent(inout) :: d
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      d%cap = 0_ik
      d%n = 0_ik
      d%active = .false.
   end subroutine free_cohort_diag

   !----- Grow to at least `need` slots, preserving the live prefix (the ensure_capacity twin). !
   subroutine expand_cohort_diag(d, need)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: need
      real(wp), allocatable :: v2(:,:), w2(:)
      integer(ik) :: m
      if (.not. d%active) return
      if (need <= d%cap) return
      !----- Keep every existing slot, not just the first n: a slow operator may have written rows    !
      !      ahead of the count (a boundary's restructuring before any fast loop has set n). ---------!
      m = d%cap
      allocate(v2(lbound(d%v,1):ubound(d%v,1), need), w2(need))
      v2 = 0.0_wp
      w2 = 0.0_wp
      if (m > 0_ik) then
         v2(:, 1:m) = d%v(:, 1:m)
         w2(1:m) = d%w(1:m)
      end if
      call move_alloc(v2, d%v)
      call move_alloc(w2, d%w)
      d%cap = need
   end subroutine expand_cohort_diag

   !----- Zero the whole block for a fresh slow step. ---------------------------------------!
   subroutine reset_cohort_diag(d)
      type(cohort_diag_block), intent(inout) :: d
      if (.not. d%active) return
      if (allocated(d%v)) d%v = 0.0_wp
      if (allocated(d%w)) d%w = 0.0_wp
   end subroutine reset_cohort_diag

   !----- Permute the live prefix in lockstep with the cohort SoA. ONE statement -- it cannot   !
   !      omit a field, which is the whole point of the 2-D layout (principle 4).                 !
   !                                                                                          !
   !      THE COUNT RIDES THE LOCKSTEP TOO. `m` is the post-operator cohort count, and         !
   !      reorder_cohort_block sets cohort%n = m right after calling this; a block that kept its     !
   !      old d%n would claim more live slots than the cohort array has. Everything            !
   !      downstream believes d%n: cohort_diag_value fills 1..d%n into a caller buffer         !
   !      that extract_variable sizes from cohort%n, so after any operator that SHRINKS        !
   !      the array -- a cull is the common one -- that write would run off the end: an        !
   !      out-of-bounds WRITE, silent in Release. Assigned BEFORE the m <= 0 bail, so          !
   !      emptying the block records the emptying. --------------------------------------------!
   subroutine reorder_cohort_diag(d, perm, m)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: perm(:), m
      if (.not. d%active) return
      d%n = max(m, 0_ik)
      if (m <= 0_ik) return
      d%v(:, 1:m) = d%v(:, perm(1:m))
      d%w(1:m)    = d%w(perm(1:m))
   end subroutine reorder_cohort_diag

   subroutine copy_cohort_diag(d, dst, src)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: dst, src
      if (.not. d%active) return
      d%v(:, dst) = d%v(:, src)
      d%w(dst)    = d%w(src)
   end subroutine copy_cohort_diag

   !----- Zero one slot (a fresh birth, or a reused stale cull slot). -----------------------!
   subroutine clear_cohort_diag(d, i)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: i
      if (.not. d%active) return
      if (i < 1_ik .or. i > d%cap) return
      d%v(:, i) = 0.0_wp
      d%w(i) = 0.0_wp
   end subroutine clear_cohort_diag

   !----- Blend the donor into the survivor on cohort fusion, each field by its own kind.      !
   !                                                                                           !
   !      `li_*` (leaf area) and `np_*` (nplant) are the weights the caller,                   !
   !      fuse_2_cohorts, computes. The extensive rows share their weight with their           !
   !      prognostic twins; the intensive rows, temperatures included, are leaf-area           !
   !      means, while the prognostic temperatures fuse on heat capacity                       !
   !      (fuse_cohort_fast_state). The dt weight `w` is a property of the window, not of      !
   !      either cohort (see below). ----------------------------------------------------------!
   subroutine fuse_cohort_diag(d, recc, donc, li_r, li_d, np_r, np_d, hl_r, hl_d, hw_r, hw_d)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: recc, donc
      real(wp),                intent(in)    :: li_r, li_d   !< leaf-area weights (nplant*leaf_area)
      real(wp),                intent(in)    :: np_r, np_d   !< nplant weights
      real(wp),                intent(in)    :: hl_r, hl_d   !< [J/m2/K] leaf heat capacities, before the merge
      real(wp),                intent(in)    :: hw_r, hw_d   !< [J/m2/K] wood heat capacities, before the merge
      integer(ik) :: f
      if (.not. d%active) return
      do f = int(lbound(d%v, 1), ik), int(ubound(d%v, 1), ik)
         select case (CO_FUSION_KIND(f))
         case (FSK_INTENSIVE)
            d%v(f, recc) = pair_mean(d%v(f, recc), li_r, d%v(f, donc), li_d)
         case (FSK_EXTENSIVE)
            d%v(f, recc) = pair_mean(d%v(f, recc), np_r, d%v(f, donc), np_d)
         case (FSK_LEAF_HEAT)
            d%v(f, recc) = pair_mean(d%v(f, recc), hl_r, d%v(f, donc), hl_d)
         case (FSK_WOOD_HEAT)
            d%v(f, recc) = pair_mean(d%v(f, recc), hw_r, d%v(f, donc), hw_d)
         case default   ! FSK_GROUND -- already per m2 of the SAME patch ground, so they add
            d%v(f, recc) = d%v(f, recc) + d%v(f, donc)
         end select
      end do
      !----- The dt weight is a property of the WINDOW, not of either cohort, and both slots saw  !
      !      the same sub-steps -- so the survivor keeps its own rather than summing (summing would !
      !      halve every normalized mean).  -------------------------------------------------!
   end subroutine fuse_cohort_diag

   !----- Extract one field, normalized by its own weight. `n` is what the caller may          !
   !      READ, and it is bounded by the caller's array as well as by the block: a kernel      !
   !      that writes past a caller's buffer on a count mismatch should not be possible,       !
   !      and this costs nothing on the once-per-step diagnostic path. With the count          !
   !      riding the lockstep (reorder_cohort_diag) the clamp never bites; it is what          !
   !      turns the next instance of that class into a short read instead of memory            !
   !      corruption. -------------------------------------------------------------------------!
   pure subroutine cohort_diag_value(d, field, x, n, empty)
      type(cohort_diag_block), intent(in)  :: d
      integer(ik),             intent(in)  :: field
      real(wp),                intent(out) :: x(:)
      integer(ik),             intent(out) :: n
      real(wp),                intent(in)  :: empty   !< what a slot with no samples reads
      integer(ik) :: i
      n = 0_ik
      if (.not. d%active) return
      !----- A cohort with no samples reads `empty`: one recruited in the slow step, after the      !
      !      fast loop, or a slot past the block's count. Reading 0 instead would make a            !
      !      recruit's leaf temperature 0 K and pull every canopy mean towards 0. ------------------!
      x = empty
      n = min(d%n, int(size(x), ik))
      do i = 1_ik, n
         if (d%w(i) > tiny_num) x(i) = d%v(field, i) / d%w(i)
      end do
   end subroutine cohort_diag_value

   !=======================================================================================!
   !  PATCH block lifecycle -- the same shape, one axis smaller.                             !
   !=======================================================================================!

   subroutine alloc_patch_diag(d, cap, active)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: cap
      logical,                intent(in)    :: active
      d%active = active
      if (.not. active) return
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      allocate(d%v(N_PA_ROWS, max(cap,1_ik)), d%w(max(cap,1_ik)))
      d%v = 0.0_wp
      d%w = 0.0_wp          ! as alloc_cohort_diag
      d%cap = max(cap, 1_ik)
      d%n = 0_ik
      call reset_patch_diag(d)
   end subroutine alloc_patch_diag

   subroutine free_patch_diag(d)
      type(patch_diag_block), intent(inout) :: d
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      d%cap = 0_ik
      d%n = 0_ik
      d%active = .false.
   end subroutine free_patch_diag

   subroutine expand_patch_diag(d, need)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: need
      real(wp), allocatable :: v2(:,:), w2(:)
      integer(ik) :: m
      if (.not. d%active) return
      if (need <= d%cap) return
      m = d%cap                              ! every existing slot (see expand_cohort_diag)
      allocate(v2(N_PA_ROWS, need), w2(need))
      v2 = 0.0_wp
      w2 = 0.0_wp
      if (m > 0_ik) then
         v2(:, 1:m) = d%v(:, 1:m)
         w2(1:m) = d%w(1:m)
      end if
      call move_alloc(v2, d%v)
      call move_alloc(w2, d%w)
      d%cap = need
   end subroutine expand_patch_diag

   subroutine reset_patch_diag(d)
      type(patch_diag_block), intent(inout) :: d
      if (.not. d%active) return
      if (allocated(d%v)) d%v = 0.0_wp
      if (allocated(d%w)) d%w = 0.0_wp
   end subroutine reset_patch_diag

   !----- The patch analogue, and the count rides it for the same reason. ------------------!
   subroutine reorder_patch_diag(d, perm, m)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: perm(:), m
      if (.not. d%active) return
      d%n = max(m, 0_ik)
      if (m <= 0_ik) return
      d%v(:, 1:m) = d%v(:, perm(1:m))
      d%w(1:m)    = d%w(perm(1:m))
   end subroutine reorder_patch_diag

   subroutine copy_patch_diag(d, dst, src)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: dst, src
      if (.not. d%active) return
      d%v(:, dst) = d%v(:, src)
      d%w(dst)    = d%w(src)
   end subroutine copy_patch_diag

   subroutine clear_patch_diag(d, i)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: i
      if (.not. d%active) return
      if (i < 1_ik .or. i > d%cap) return
      d%v(:, i) = 0.0_wp
      d%w(i) = 0.0_wp
   end subroutine clear_patch_diag

   !----- Area-weighted blend on patch fusion. Every patch diagnostic is per m2 of ITS OWN     !
   !      ground, so all of them blend the same way -- there is no extensive/intensive split at  !
   !      this level, which is why there is no kind array here.                                  !
   subroutine blend_patch_diag(d, recp, donp, area_r, area_d)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: recp, donp
      real(wp),               intent(in)    :: area_r, area_d
      real(wp) :: wtot
      if (.not. d%active) return
      wtot = area_r + area_d
      if (wtot <= tiny_num) return
      d%v(:, recp) = (area_r * d%v(:, recp) + area_d * d%v(:, donp)) / wtot
   end subroutine blend_patch_diag

   !----- A patch carved out of others inherits their slots, sums and weight, blended by the area  !
   !      each donor gave. Every donor loses the same fraction, so blending by the donors' areas is  !
   !      exact and the stand total Sum(area * v) is unchanged by the carving; a cleared slot would  !
   !      lose the carved area's share. At a calendar boundary the donors hold only that boundary's !
   !      events (weight 0: the step's sums were read and reset); inside a step (the C API) they    !
   !      hold the step so far, and the site mean Sum(area * v/w) is unchanged as well. `n` is      !
   !      raised to cover the new slot, so a reader sees it before the next reorder.                !
   subroutine inherit_patch_diag(d, dst, donor_area)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: dst
      real(wp),               intent(in)    :: donor_area(:)   !< donors are slots 1..size(donor_area)
      real(wp)    :: atot, wd
      integer(ik) :: i
      if (.not. d%active) return
      if (dst < 1_ik .or. dst > d%cap) return
      d%v(:, dst) = 0.0_wp
      d%w(dst) = 0.0_wp
      d%n = max(d%n, dst)
      atot = sum(donor_area)
      if (atot <= tiny_num) return
      do i = 1_ik, size(donor_area, kind=ik)
         wd = donor_area(i) / atot
         d%v(:, dst) = d%v(:, dst) + wd * d%v(:, i)
         d%w(dst)    = d%w(dst)    + wd * d%w(i)
      end do
   end subroutine inherit_patch_diag

   !----- The rows the slow operators write (litterfall, recruitment, disturbance, mortality      !
   !      carbon): the ones a slow-only run fills. Every other row is the fast loop's. ------------!
   pure logical function patch_diag_slow_row(field) result(slow)
      integer(ik), intent(in) :: field
      !----- The slow driver's rows follow the fast loop's (see N_PA_FAST_ROWS). ------------------!
      slow = field > N_PA_FAST_ROWS
   end function patch_diag_slow_row

   pure subroutine patch_diag_value(d, field, x, n)
      type(patch_diag_block), intent(in)  :: d
      integer(ik),            intent(in)  :: field
      real(wp),               intent(out) :: x(:)
      integer(ik),            intent(out) :: n
      integer(ik) :: i
      n = 0_ik
      if (.not. d%active) return
      x = 0.0_wp                         ! slots with no entry read 0
      n = min(d%n, int(size(x), ik))     ! bounded by the caller's array too (see cohort_diag_value)
      do i = 1_ik, n
         if (d%w(i) > tiny_num) x(i) = d%v(field, i) / d%w(i)
      end do
   end subroutine patch_diag_value

   !=======================================================================================!
   !  POLYGON block. Reset at the start of each slow step's fast loop, read as the step mean.   !
   !=======================================================================================!
   subroutine reset_polygon_diag(d)
      type(polygon_diag_block), intent(inout) :: d
      if (.not. d%active) return
      d%v = 0.0_wp
      d%w = 0.0_wp
   end subroutine reset_polygon_diag

   !----- The step mean of one field; 0 if nothing was accumulated, as the patch block reads. --!
   pure real(wp) function polygon_diag_value(d, field) result(x)
      type(polygon_diag_block), intent(in) :: d
      integer(ik),              intent(in) :: field
      x = 0.0_wp
      if (d%active .and. d%w > tiny_num) x = d%v(field) / d%w
   end function polygon_diag_value

   !=======================================================================================!
   !  SITE block. Reset at the start of each slow step's fast loop, read as the step mean. !
   !=======================================================================================!
   subroutine alloc_site_diag(d, n_site, on)
      type(site_diag_block), intent(inout) :: d
      integer(ik),           intent(in)    :: n_site
      logical,               intent(in)    :: on
      d%active = on
      if (.not. on .or. allocated(d%v)) return
      allocate(d%v(N_SI_ROWS, n_site), d%w(n_site))
      d%v = 0.0_wp
      d%w = 0.0_wp
   end subroutine alloc_site_diag

   subroutine reset_site_diag(d)
      type(site_diag_block), intent(inout) :: d
      if (.not. d%active) return
      d%v = 0.0_wp
      d%w = 0.0_wp
   end subroutine reset_site_diag

   !----- The step mean of one field at each site; `valid` false where nothing was accumulated. -!
   pure subroutine site_diag_value(d, field, n_site, x, valid)
      type(site_diag_block), intent(in)  :: d
      integer(ik),           intent(in)  :: field, n_site
      real(wp),              intent(out) :: x(:)
      logical,               intent(out) :: valid(:)
      integer(ik) :: is
      do is = 1_ik, n_site
         valid(is) = d%active .and. d%w(is) > tiny_num
         x(is) = 0.0_wp
         if (valid(is)) x(is) = d%v(field, is) / d%w(is)
      end do
   end subroutine site_diag_value

end module ed_diag_vars
