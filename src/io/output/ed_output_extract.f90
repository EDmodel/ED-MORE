! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_extract -- what a diagnostic variable reads: the source ids (FLD_*, SRC_*),
! each naming a field and, by its range, the entity it lives on; the field accessors,
! which copy a field out of live state and derive the few closed forms no physics library
! owns (stage [1]); and extract_variable, which reduces it to the scale a variable is
! written at (stage [3]). extract_fast_scalar reads the FAST tier's staging instead. A new
! field is a source id here and an add_variable line in ed_output_registry
! (docs/science/diagnostics.md §1, §8). netCDF-free.
!---~---
module ed_output_extract
   use ed_kinds,          only : wp, ik
   use ed_biogeochem_params,   only : IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, &
                                      IP_MICR, IP_SLOW, IP_PASSIVE
   use ed_output_types,     only : var_desc_t, diag_setup_t, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT,            &
                                   DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE, MISSING_VALUE
   use ed_polygon_state_vars,   only : polygon_t
   use ed_diag_vars,            only : N_CO_FAST_ROWS, N_CO_ROWS, N_PA_ROWS, N_PY_ROWS, cohort_diag_value, patch_diag_value, &
                                       patch_diag_slow_row, polygon_diag_value, N_SI_ROWS, site_diag_value
   use ed_state_defaults,       only : is_unset_positive
   use ed_constants,     only : n_soil_layer_max
   use ed_soil_hydr_lib,      only : soil_matric_potential, soil_theta_to_saturation
   use ed_output_reduce,  only : reduce_cohort_to_polygon, reduce_cohort_to_patch,             &
                                       reduce_cohort_to_pft, reduce_cohort_to_size,               &
                                       reduce_patch_to_polygon, reduce_patch_columns_to_sites,       &
                                       gather_patch_columns, W_NPLANT
   implicit none
   private

   public :: extract_variable, extract_fast_scalar
   public :: src_class, SRCK_COHORT, SRCK_PATCH, SRCK_LAYER, SRCK_POLYGON, SRCK_FAST, SRCK_SITE, SRC_NONE
   !--- Cohort fields (1000-1999).
   public :: FLD_C_NPLANT, FLD_C_DBH, FLD_C_HEIGHT, FLD_C_BASAL_AREA, FLD_C_AGB, FLD_C_LEAF_AREA
   public :: FLD_C_GROWTH_AVG, FLD_C_PFT, FLD_C_OWNER_PATCH, FLD_C_GLOBAL_ID
   public :: FLD_C_LAI, FLD_C_WAI, FLD_C_LEAF_CARBON, FLD_C_FINEROOT_CARBON, FLD_C_WOOD_CARBON
   public :: FLD_C_NONSTRUCTURAL_CARBON, FLD_C_BGB, FLD_C_VEG_CARBON, FLD_C_ONE
   public :: FLD_C_SLA, FLD_C_VCMAX25, FLD_C_RD25, FLD_C_LEAF_TURNOVER_RATE, FLD_C_OVERTOPPING_LAI
   public :: FLD_C_GPP_ACCUM, FLD_C_NPP_ACCUM, FLD_C_LEAF_RESP, FLD_C_STEM_RESP, FLD_C_ROOT_RESP
   public :: FLD_C_DMAX_PSI_LEAF, FLD_C_LEAF_FLUSH_TENDENCY, FLD_C_LEAF_SHED_TENDENCY
   public :: FLD_C_DIAG0, FLD_P_DIAG0, FLD_PY_DIAG0
   !--- Patch fields (2000-2999).
   public :: FLD_P_AREA, FLD_P_AGE, FLD_P_DIST_TYPE, FLD_P_COHORT_OFFSET, FLD_P_COHORT_COUNT
   public :: FLD_P_GLOBAL_ID, FLD_P_SITE, FLD_P_CAS_DEPTH, FLD_P_SWE, FLD_P_SNOW_DEPTH
   public :: FLD_P_SOILC_FAST_GRND, FLD_P_SOILC_FAST_SOIL, FLD_P_SOILC_STRUCT_GRND
   public :: FLD_P_SOILC_STRUCT_SOIL, FLD_P_SOILC_MICROBIAL, FLD_P_SOILC_SLOW
   public :: FLD_P_SOILC_PASSIVE, FLD_P_SOILC_TOTAL, FLD_P_RHET
   !--- Soil layer fields (3000-3999).
   public :: FLD_L_SOIL_TEMP, FLD_L_SOIL_WATER, FLD_L_SOIL_PSI, FLD_L_SOIL_WETNESS, FLD_L_SOIL_FLIQ
   !--- The polygon's own scalars (4000-4999).
   public :: SRC_PY_N_COHORT, SRC_PY_N_PATCH, SRC_PY_CANOPY_HEIGHT, SRC_PY_ECO_CARBON
   public :: SRC_PY_WORK_STEPS, SRC_PY_WORK_REJ, SRC_PY_WORK_WATER_NSUB,                              &
             SRC_PY_WORK_WATER_NEWTON, SRC_PY_WORK_WATER_NFAIL, SRC_PY_WORK_WATER_REJECT
   public :: SRC_PY_WORK_RK45_RESCUE, SRC_PY_WORK_CLAMP_STAGE
   !--- FAST-tier sources (5000-5999): read from the fast loop's staging.
   public :: SRC_F_SOIL_TEMP, SRC_F_SOIL_WATER, SRC_F_SOIL_LAYER
   public :: SRC_F_COH_GPP, SRC_F_COH_HEIGHT, SRC_F_COH_NPP, SRC_F_CO0
   public :: SRC_F_PY0, SRC_F_PA0
   !--- Per-site fields (6000-6999).
   public :: FLD_SI_GROWTH_TEMP, FLD_SI_DIAG0

   !---~---
   !   Source ids: an id names a field, and its numeric range the entity the field lives on, so
   ! the variable's dim alone picks the reduction. The ranges are disjoint so the coarse and
   ! FAST readers cannot resolve each other's ids (docs/science/diagnostics.md §8, Source ids).
   !---~---
   integer(ik), parameter :: SRCK_NONE   = 0_ik
   integer(ik), parameter :: SRCK_COHORT = 1_ik   ! 1000-1999 : per-cohort field
   integer(ik), parameter :: SRCK_PATCH  = 2_ik   ! 2000-2999 : per-patch field
   integer(ik), parameter :: SRCK_LAYER  = 3_ik   ! 3000-3999 : per-soil-layer field (per patch)
   integer(ik), parameter :: SRCK_POLYGON   = 4_ik   ! 4000-4999 : the polygon's own scalar (no reduction)
   integer(ik), parameter :: SRCK_FAST   = 5_ik   ! 5000-5999 : FAST-tier staged sample
   integer(ik), parameter :: SRC_NONE    = 0_ik   ! no source: a variable without a value on those tiers
   integer(ik), parameter :: SRCK_SITE  = 6_ik   ! 6000-6999 : one value per site of the polygon

   !--- Cohort fields (SRCK_COHORT). Raw SoA columns first, then derived ones.
   integer(ik), parameter :: FLD_C_NPLANT          = 1001_ik
   integer(ik), parameter :: FLD_C_DBH             = 1002_ik
   integer(ik), parameter :: FLD_C_HEIGHT          = 1003_ik
   integer(ik), parameter :: FLD_C_BASAL_AREA      = 1004_ik
   integer(ik), parameter :: FLD_C_AGB             = 1005_ik
   integer(ik), parameter :: FLD_C_LEAF_AREA       = 1006_ik
   integer(ik), parameter :: FLD_C_GROWTH_AVG      = 1007_ik
   integer(ik), parameter :: FLD_C_PFT             = 1008_ik
   integer(ik), parameter :: FLD_C_OWNER_PATCH     = 1009_ik
   integer(ik), parameter :: FLD_C_GLOBAL_ID       = 1010_ik
   integer(ik), parameter :: FLD_C_LEAF_CARBON     = 1011_ik
   integer(ik), parameter :: FLD_C_FINEROOT_CARBON = 1012_ik
   integer(ik), parameter :: FLD_C_WOOD_CARBON     = 1013_ik
   integer(ik), parameter :: FLD_C_NONSTRUCTURAL_CARBON  = 1014_ik
   integer(ik), parameter :: FLD_C_SLA             = 1015_ik
   integer(ik), parameter :: FLD_C_VCMAX25         = 1016_ik
   integer(ik), parameter :: FLD_C_RD25            = 1017_ik
   integer(ik), parameter :: FLD_C_LEAF_TURNOVER_RATE   = 1018_ik
   integer(ik), parameter :: FLD_C_OVERTOPPING_LAI     = 1019_ik
   integer(ik), parameter :: FLD_C_GPP_ACCUM       = 1020_ik
   integer(ik), parameter :: FLD_C_LEAF_RESP       = 1021_ik
   integer(ik), parameter :: FLD_C_STEM_RESP       = 1022_ik
   integer(ik), parameter :: FLD_C_ROOT_RESP       = 1023_ik
   integer(ik), parameter :: FLD_C_DMAX_PSI_LEAF   = 1024_ik
   integer(ik), parameter :: FLD_C_LEAF_FLUSH_TENDENCY     = 1025_ik
   integer(ik), parameter :: FLD_C_LEAF_SHED_TENDENCY      = 1026_ik
   !---~---
   !   Leaf and wood temperature are not read here at the output tick: one sample per window,
   ! at the boundary's local hour, biases them. Their time means are rows of the
   ! cohort block (FLD_C_DIAG0 + CO_LEAF_TEMP, CO_WOOD_TEMP), sampled every sub-step.
   !---~---
   !--- Derived cohort fields: closed forms of the cohort state.
   integer(ik), parameter :: FLD_C_ONE             = 1050_ik ! constant 1 (stem counts on any axis)
   integer(ik), parameter :: FLD_C_LAI             = 1051_ik ! nplant*leaf_area [m2/m2]
   integer(ik), parameter :: FLD_C_NPP_ACCUM       = 1052_ik ! gpp - (leaf+stem+root) maintenance resp
   integer(ik), parameter :: FLD_C_BGB             = 1053_ik ! belowground biomass per plant
   integer(ik), parameter :: FLD_C_VEG_CARBON      = 1054_ik ! leaf+fineroot+wood+storage per plant
   integer(ik), parameter :: FLD_C_WAI             = 1055_ik ! nplant*wood_area [m2/m2]
   !---~---
   !   Per-cohort diagnostics (1100+): the id is FLD_C_DIAG0 + the row (CO_*). The fast loop's
   ! rows are read from poly%cohort%diag, the slow step's (above N_CO_FAST_ROWS) from
   ! poly%cohort%sdiag; never from poly%deriv, which is not lockstep-reordered, so at the
   ! output tick its index need not match the cohort axis (docs/science/diagnostics.md §8).
   !---~---
   integer(ik), parameter :: FLD_C_DIAG0           = 1100_ik
   !--- Per-patch diagnostics (2100+), read from poly%patch%diag (FLD_P_DIAG0 + PA_*).
   integer(ik), parameter :: FLD_P_DIAG0           = 2100_ik
   !---~---
   !   The polygon block's rows (4100+, among the polygon's own scalars), read from poly%diag
   ! (FLD_PY_DIAG0 + PY_*).
   !---~---
   integer(ik), parameter :: FLD_PY_DIAG0          = 4100_ik
   !---~---
   !   Per-site fields (6000+, SRCK_SITE): the growth temperature, and the site block's
   ! forcing rows (SI_* + FLD_SI_DIAG0).
   !---~---
   integer(ik), parameter :: FLD_SI_GROWTH_TEMP    = 6001_ik
   integer(ik), parameter :: FLD_SI_DIAG0          = 6100_ik

   !--- Patch fields (SRCK_PATCH).
   integer(ik), parameter :: FLD_P_AREA            = 2001_ik
   integer(ik), parameter :: FLD_P_AGE             = 2002_ik
   integer(ik), parameter :: FLD_P_DIST_TYPE       = 2003_ik
   integer(ik), parameter :: FLD_P_COHORT_OFFSET   = 2004_ik
   integer(ik), parameter :: FLD_P_COHORT_COUNT    = 2005_ik
   integer(ik), parameter :: FLD_P_GLOBAL_ID       = 2006_ik
   integer(ik), parameter :: FLD_P_SITE            = 2007_ik   ! the patch's site (owner_site)
   !---~---
   !   The canopy-air, soil-top and surface-water states are not read here at the output tick:
   ! one sample per window, at the boundary's local hour, biases them. Their time means
   ! are rows of the patch block (FLD_P_DIAG0 + PA_*), accumulated through the step.
   !---~---
   integer(ik), parameter :: FLD_P_CAS_DEPTH       = 2014_ik
   integer(ik), parameter :: FLD_P_SWE             = 2016_ik
   integer(ik), parameter :: FLD_P_SNOW_DEPTH      = 2017_ik
   integer(ik), parameter :: FLD_P_SOILC_FAST_GRND   = 2020_ik
   integer(ik), parameter :: FLD_P_SOILC_FAST_SOIL   = 2021_ik
   integer(ik), parameter :: FLD_P_SOILC_STRUCT_GRND = 2022_ik
   integer(ik), parameter :: FLD_P_SOILC_STRUCT_SOIL = 2023_ik
   integer(ik), parameter :: FLD_P_SOILC_MICROBIAL   = 2024_ik
   integer(ik), parameter :: FLD_P_SOILC_SLOW        = 2025_ik
   integer(ik), parameter :: FLD_P_SOILC_PASSIVE     = 2026_ik
   integer(ik), parameter :: FLD_P_SOILC_TOTAL       = 2027_ik
   integer(ik), parameter :: FLD_P_RHET                = 2028_ik

   !---~---
   !   Soil layer fields (SRCK_LAYER): per (layer, patch); reduced to each site's column
   ! for DIM_SOIL, or flattened for DIM_SOIL_PATCH.
   !---~---
   integer(ik), parameter :: FLD_L_SOIL_TEMP    = 3001_ik
   integer(ik), parameter :: FLD_L_SOIL_WATER   = 3002_ik
   integer(ik), parameter :: FLD_L_SOIL_PSI     = 3003_ik
   integer(ik), parameter :: FLD_L_SOIL_WETNESS = 3004_ik
   integer(ik), parameter :: FLD_L_SOIL_FLIQ    = 3005_ik

   !--- The polygon's own scalars (SRCK_POLYGON): already polygon quantities, no reduction.
   integer(ik), parameter :: SRC_PY_N_COHORT         = 4002_ik
   integer(ik), parameter :: SRC_PY_N_PATCH          = 4003_ik
   integer(ik), parameter :: SRC_PY_CANOPY_HEIGHT    = 4004_ik ! [m] tallest live cohort
   integer(ik), parameter :: SRC_PY_ECO_CARBON       = 4008_ik ! [kgC/m2] live vegetation + soil carbon
   integer(ik), parameter :: SRC_PY_WORK_STEPS       = 4030_ik
   integer(ik), parameter :: SRC_PY_WORK_REJ         = 4031_ik
   integer(ik), parameter :: SRC_PY_WORK_WATER_NSUB   = 4032_ik
   integer(ik), parameter :: SRC_PY_WORK_WATER_NEWTON = 4033_ik
   integer(ik), parameter :: SRC_PY_WORK_WATER_NFAIL  = 4034_ik
   integer(ik), parameter :: SRC_PY_WORK_RK45_RESCUE = 4035_ik
   integer(ik), parameter :: SRC_PY_WORK_CLAMP_STAGE = 4036_ik
   integer(ik), parameter :: SRC_PY_WORK_WATER_REJECT = 4040_ik

   !---~---
   !   FAST-tier staged sources (SRCK_FAST): read from the fast loop's staging, never from the
   ! polygon's live state, the end of the slow step by replay time. Two are rows of the
   ! blocks' own tables: the forcing (docs/science/diagnostics.md §5), SRC_F_PY0 + PY_*, and
   ! the patch block's fast rows, SRC_F_PA0 + PA_*, the polygon mean or per patch by the dim.
   !---~---
   integer(ik), parameter :: SRC_F_PY0           = 5100_ik
   integer(ik), parameter :: SRC_F_PA0           = 5200_ik
   integer(ik), parameter :: SRC_F_SOIL_TEMP     = 5010_ik  ! DIM_SOIL or DIM_SOIL_PATCH, the fast soil slabs
   integer(ik), parameter :: SRC_F_SOIL_WATER    = 5011_ik  ! DIM_SOIL or DIM_SOIL_PATCH, the fast soil slabs
   !--- A soil-layer field formed from the staged water content and liquid fraction (layer_source_field).
   integer(ik), parameter :: SRC_F_SOIL_LAYER    = 5012_ik  ! DIM_SOIL or DIM_SOIL_PATCH
   integer(ik), parameter :: SRC_F_COH_GPP       = 5021_ik  ! DIM_COHORT
   integer(ik), parameter :: SRC_F_COH_HEIGHT    = 5022_ik  ! DIM_COHORT
   integer(ik), parameter :: SRC_F_COH_NPP       = 5023_ik  ! DIM_COHORT
   !---~---
   !   The per-cohort diagnostic block's rows at the FAST tier, SRC_F_CO0 + CO_*, DIM_COHORT: the
   ! fast loop stages every row of every sub-step when the block is on (fast_cohort_diag).
   !---~---
   integer(ik), parameter :: SRC_F_CO0           = 5300_ik

contains

   !---~---
   !   Extract: copy one variable's current value out of live state, with no pointers.
   ! nvfortran -O2 can merge two extracts across an inline store to poly: extract each variable
   ! once per step.
   !---~---
   !---~---
   !   A diagnostic-block row whose block this run does not fill: a block that is off, or a fast
   ! row of the patch block in a slow-only run. Its variables read as missing, not as a
   ! 0 nothing computed.
   !---~---
   pure logical function not_simulated(poly, src) result(ns)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in) :: src
      ns = .false.
      if (src > FLD_C_DIAG0 .and. src <= FLD_C_DIAG0 + N_CO_FAST_ROWS) then
         ns = .not. poly%cohort%diag%active
      else if (src > FLD_C_DIAG0 + N_CO_FAST_ROWS .and. src <= FLD_C_DIAG0 + N_CO_ROWS) then
         ns = .not. poly%cohort%sdiag%active
      else if (src > FLD_P_DIAG0 .and. src <= FLD_P_DIAG0 + N_PA_ROWS) then
         ns = .not. poly%patch%diag%active .or.                                                    &
              (.not. poly%patch%diag%fast_rows .and. .not. patch_diag_slow_row(src - FLD_P_DIAG0))
      else if (src > FLD_PY_DIAG0 .and. src <= FLD_PY_DIAG0 + N_PY_ROWS) then
         ns = .not. poly%diag%active
      else if (src > FLD_SI_DIAG0 .and. src <= FLD_SI_DIAG0 + N_SI_ROWS) then
         ns = .not. poly%site%diag%active
      end if
   end function not_simulated

   !---~---
   !   Which entity a source id names, from its range (the Source ids block above). Every
   ! dispatch asks it, so an id outside the ranges is SRCK_NONE, never a neighbour's field.
   !---~---
   pure integer(ik) function src_class(src) result(k)
      integer(ik), intent(in) :: src
      if      (src >= 1000_ik .and. src < 2000_ik) then
         k = SRCK_COHORT
      else if (src >= 2000_ik .and. src < 3000_ik) then
         k = SRCK_PATCH
      else if (src >= 3000_ik .and. src < 4000_ik) then
         k = SRCK_LAYER
      else if (src >= 4000_ik .and. src < 5000_ik) then
         k = SRCK_POLYGON
      else if (src >= 5000_ik .and. src < 6000_ik) then
         k = SRCK_FAST
      else if (src >= 6000_ik .and. src < 7000_ik) then
         k = SRCK_SITE
      else
         k = SRCK_NONE
      end if
   end function src_class

   !---~---
   !   Field accessors: copy one entity's column out of live state, deriving the closed
   ! forms of stage [1], with no reduction; the descriptor's dim chooses that.
   !---~---

   !--- Per-cohort field over the polygon's whole cohort array [1:poly%cohort%n].
   pure subroutine cohort_source_field(poly, src, x, n)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in)  :: src
      real(wp),     intent(out) :: x(:)
      integer(ik),  intent(out) :: n
      n = poly%cohort%n
      if (n <= 0_ik) return
      select case (src)
      case (FLD_C_NPLANT)
         x(1:n) = poly%cohort%nplant(1:n)
      case (FLD_C_DBH)
         x(1:n) = poly%cohort%dbh(1:n)
      case (FLD_C_HEIGHT)
         x(1:n) = poly%cohort%height(1:n)
      case (FLD_C_BASAL_AREA)
         x(1:n) = poly%cohort%basal_area(1:n)
      case (FLD_C_AGB)
         x(1:n) = poly%cohort%agb(1:n)
      case (FLD_C_LEAF_AREA)
         x(1:n) = poly%cohort%leaf_area(1:n)
      case (FLD_C_GROWTH_AVG)
         x(1:n) = poly%cohort%growth_avg(1:n)
      case (FLD_C_PFT)
         x(1:n) = real(poly%cohort%pft(1:n), wp)
      case (FLD_C_OWNER_PATCH)
         x(1:n) = real(poly%cohort%owner_patch(1:n), wp)
      case (FLD_C_GLOBAL_ID)
         x(1:n) = real(poly%cohort%global_id(1:n), wp)
      case (FLD_C_LEAF_CARBON)
         x(1:n) = poly%cohort%leaf_carbon(1:n)
      case (FLD_C_FINEROOT_CARBON)
         x(1:n) = poly%cohort%fineroot_carbon(1:n)
      case (FLD_C_WOOD_CARBON)
         x(1:n) = poly%cohort%wood_carbon(1:n)
      case (FLD_C_NONSTRUCTURAL_CARBON)
         x(1:n) = poly%cohort%nonstructural_carbon(1:n)
      case (FLD_C_SLA)
         x(1:n) = poly%cohort%sla(1:n)
      case (FLD_C_VCMAX25)
         x(1:n) = poly%cohort%vcmax25(1:n)
      case (FLD_C_RD25)
         x(1:n) = poly%cohort%rd25(1:n)
      case (FLD_C_LEAF_TURNOVER_RATE)
         x(1:n) = poly%cohort%leaf_turnover_rate(1:n)
      case (FLD_C_OVERTOPPING_LAI)
         x(1:n) = poly%cohort%overtopping_lai(1:n)
      case (FLD_C_GPP_ACCUM)
         x(1:n) = poly%cohort%gpp_accum(1:n)
      case (FLD_C_LEAF_RESP)
         x(1:n) = poly%cohort%leaf_resp_accum(1:n)
      case (FLD_C_STEM_RESP)
         x(1:n) = poly%cohort%stem_resp_accum(1:n)
      case (FLD_C_ROOT_RESP)
         x(1:n) = poly%cohort%root_resp_accum(1:n)
      case (FLD_C_DMAX_PSI_LEAF)
         x(1:n) = poly%cohort%dmax_psi_leaf(1:n)
      case (FLD_C_LEAF_FLUSH_TENDENCY)
         x(1:n) = poly%cohort%leaf_flush_tendency(1:n)
      case (FLD_C_LEAF_SHED_TENDENCY)
         x(1:n) = poly%cohort%leaf_shed_tendency(1:n)
      !--- derived
      case (FLD_C_ONE)
         x(1:n) = 1.0_wp
      case (FLD_C_LAI)
         x(1:n) = poly%cohort%nplant(1:n) * poly%cohort%leaf_area(1:n)
      case (FLD_C_WAI)
         x(1:n) = poly%cohort%nplant(1:n) * poly%cohort%wood_area(1:n)
      case (FLD_C_NPP_ACCUM)
         !---~---
         !   Gross GPP less the three maintenance respirations; growth respiration is charged in
         ! the slow step's allocation (mean_growth_resp_py), so subtracting it here would count
         ! it twice.
         !---~---
         x(1:n) = poly%cohort%gpp_accum(1:n) - poly%cohort%leaf_resp_accum(1:n)                     &
                  - poly%cohort%stem_resp_accum(1:n) - poly%cohort%root_resp_accum(1:n)
      case (FLD_C_BGB)
         x(1:n) = poly%cohort%fineroot_carbon(1:n)                  &
                                    + poly%cohort%wood_carbon(1:n)                             &
                                      * (1.0_wp - poly%cohort%pft_table%allom(poly%cohort%pft(1:n))%aboveground_frac)
      case (FLD_C_VEG_CARBON)
         x(1:n) = poly%cohort%leaf_carbon(1:n)                      &
                                    + poly%cohort%fineroot_carbon(1:n)                         &
                                    + poly%cohort%wood_carbon(1:n)                             &
                                    + poly%cohort%nonstructural_carbon(1:n)
      case default
         !---~---
         !   Diagnostic-block rows (FLD_C_DIAG0 + CO_*): read the dt-weighted
         ! accumulator and normalize. A block that is off never reaches here (not_simulated).
         !---~---
         if (src > FLD_C_DIAG0 .and. src <= FLD_C_DIAG0 + N_CO_FAST_ROWS) then
            call cohort_diag_value(poly%cohort%diag, src - FLD_C_DIAG0, x, n, MISSING_VALUE)
         else if (src > FLD_C_DIAG0 + N_CO_FAST_ROWS .and. src <= FLD_C_DIAG0 + N_CO_ROWS) then
            call cohort_diag_value(poly%cohort%sdiag, src - FLD_C_DIAG0, x, n, MISSING_VALUE)
         else
            x(1:n) = MISSING_VALUE
         end if
      end select
   end subroutine cohort_source_field

   !--- Per-patch field [1:poly%patch%n].
   pure subroutine patch_source_field(poly, src, x, n)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in)  :: src
      real(wp),     intent(out) :: x(:)
      integer(ik),  intent(out) :: n
      integer(ik) :: ip
      n = poly%patch%n
      if (n <= 0_ik) return
      select case (src)
      case (FLD_P_AREA)
         x(1:n) = poly%patch%area(1:n)
      case (FLD_P_AGE)
         x(1:n) = poly%patch%age(1:n)
      case (FLD_P_DIST_TYPE)
         x(1:n) = real(poly%patch%dist_type(1:n), wp)
      case (FLD_P_COHORT_OFFSET)
         x(1:n) = real(poly%patch%cohort_offset(1:n), wp)
      case (FLD_P_COHORT_COUNT)
         x(1:n) = real(poly%patch%cohort_count(1:n), wp)
      case (FLD_P_GLOBAL_ID)
         x(1:n) = real(poly%patch%global_id(1:n), wp)
      case (FLD_P_SITE)
         x(1:n) = real(poly%patch%owner_site(1:n), wp)
      case (FLD_P_CAS_DEPTH)
         do ip = 1_ik, n
            x(ip) = poly%patch%cas(ip)%can_depth
         end do
      case (FLD_P_SWE)
         do ip = 1_ik, n
            x(ip) = sum(poly%patch%snow(ip)%swe(:))
         end do
      case (FLD_P_SNOW_DEPTH)
         do ip = 1_ik, n
            x(ip) = sum(poly%patch%snow(ip)%snow_depth(:))
         end do
      case (FLD_P_SOILC_FAST_GRND)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_FAST_GRND)
         end do
      case (FLD_P_SOILC_FAST_SOIL)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_FAST_SOIL)
         end do
      case (FLD_P_SOILC_STRUCT_GRND)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_STRUCT_GRND)
         end do
      case (FLD_P_SOILC_STRUCT_SOIL)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_STRUCT_SOIL)
         end do
      case (FLD_P_SOILC_MICROBIAL)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_MICR)
         end do
      case (FLD_P_SOILC_SLOW)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_SLOW)
         end do
      case (FLD_P_SOILC_PASSIVE)
         do ip = 1_ik, n
            x(ip) = poly%patch%soil_carbon(ip)%carbon(IP_PASSIVE)
         end do
      case (FLD_P_SOILC_TOTAL)
         do ip = 1_ik, n
            x(ip) = sum(poly%patch%soil_carbon(ip)%carbon)
         end do
      case (FLD_P_RHET)
         do ip = 1_ik, n
            x(ip) = poly%patch%decomp_xi(ip)%rh_fast
         end do
      case default
         if (src > FLD_P_DIAG0 .and. src <= FLD_P_DIAG0 + N_PA_ROWS) then
            call patch_diag_value(poly%patch%diag, src - FLD_P_DIAG0, x, n)
         else
            x(1:n) = MISSING_VALUE
         end if
      end select
   end subroutine patch_source_field

   !---~---
   !   Per-(soil layer, patch) field: col is (n_soil_layer_max, n_patch), each patch's column
   ! on its site's active layers. psi and wetness take the site's retention curve, the one
   ! the fast loop integrates, and are missing without a fast loop (setup%soil_ready).
   ! Only the active layers are data and marked valid: the padding would average zeros
   ! into a column, and soil_matric_potential on it (theta_sat 0) divides by zero.
   !---~---
   pure subroutine layer_source_field(poly, setup, src, col, ok, theta_x, fliq_x)
      type(polygon_t),     intent(in)  :: poly
      type(diag_setup_t),  intent(in)  :: setup
      integer(ik),         intent(in)  :: src
      real(wp),            intent(out) :: col(:,:)
      logical,             intent(out) :: ok
      !--- In place of the live columns: a staged sub-step's water content and liquid fraction.
      real(wp), optional,  intent(in)  :: theta_x(:,:), fliq_x(:,:)   ! (n_soil_layer_max, patch)
      real(wp)    :: theta(n_soil_layer_max), fliq(n_soil_layer_max)
      integer(ik) :: ip, np, nl
      np = poly%patch%n
      ok = np > 0_ik
      if (.not. ok) return
      if ((src == FLD_L_SOIL_PSI .or. src == FLD_L_SOIL_WETNESS) .and. .not. setup%soil_ready) then
         ok = .false.
         return
      end if
      do ip = 1_ik, np
         associate (soil => poly%site%surface(poly%patch%owner_site(ip))%soil_profile)
            nl = soil%n_active
            if (present(theta_x)) then
               theta(1:nl) = theta_x(1:nl, ip)
               fliq(1:nl)  = fliq_x(1:nl, ip)
            else
               theta(1:nl) = poly%patch%soil_water(ip)%theta(1:nl)
               fliq(1:nl)  = poly%patch%soil_energy(ip)%soil_fliq(1:nl)
            end if
            select case (src)
            case (FLD_L_SOIL_TEMP)
               col(1:nl, ip) = poly%patch%soil_energy(ip)%soil_temp(1:nl)
            case (FLD_L_SOIL_WATER)
               col(1:nl, ip) = theta(1:nl)
            case (FLD_L_SOIL_FLIQ)
               col(1:nl, ip) = fliq(1:nl)
            case (FLD_L_SOIL_PSI)
               !---~---
               !   The liquid's potential: in a frozen layer, the freezing curve's
               ! (docs/science/soil_biophysics.md, "Frozen soil").
               !---~---
               col(1:nl, ip) = soil_matric_potential(soil%retention, theta(1:nl) * fliq(1:nl),         &
                                                     soil%theta_sat(1:nl), soil%theta_res(1:nl),          &
                                                     soil%curve_par_a(1:nl), soil%curve_par_n(1:nl))
            case (FLD_L_SOIL_WETNESS)
               col(1:nl, ip) = soil_theta_to_saturation(theta(1:nl), soil%theta_sat(1:nl), soil%theta_res(1:nl))
            case default
               ok = .false.
            end select
         end associate
      end do
   end subroutine layer_source_field

   !---~---
   !   Per-site field [1:poly%site%n]: a row of the site block's forcing (SI_*), or the growth
   ! temperature; a site with nothing to report this step is not valid.
   !---~---
   pure subroutine site_source_field(poly, src, x, valid, n)
      type(polygon_t), intent(in)  :: poly
      integer(ik),     intent(in)  :: src
      real(wp),        intent(out) :: x(:)
      logical,         intent(out) :: valid(:)
      integer(ik),     intent(out) :: n
      integer(ik) :: is
      n = poly%site%n
      if (src > FLD_SI_DIAG0 .and. src <= FLD_SI_DIAG0 + N_SI_ROWS) then
         call site_diag_value(poly%site%diag, src - FLD_SI_DIAG0, n, x, valid)
      else if (src == FLD_SI_GROWTH_TEMP) then
         do is = 1_ik, n
            valid(is) = .not. is_unset_positive(poly%site%t_growth_avg(is))
            x(is) = poly%site%t_growth_avg(is)
         end do
      else
         valid(1:n) = .false.
         x(1:n) = MISSING_VALUE
      end if
   end subroutine site_source_field

   !---~---
   !   The dispatcher: one variable -> one scalar or slab, its field reduced to the scale its
   ! dim names with its descriptor's weight, mean and scale. One field serves every scale:
   ! mean_agb_co, mean_agb_pa, mean_agb_py, mean_agb_pft_py and mean_agb_size_py all read FLD_C_AGB. Extract each
   ! variable once per step (the nvfortran trap above).
   !---~---
   subroutine extract_variable(poly, setup, v, scalar_out, slab_out, valid_out, n_out, cohort_x, cohort_w, soil_theta, &
                               soil_fliq)
      type(polygon_t),     intent(in)  :: poly
      type(diag_setup_t),  intent(in)  :: setup
      type(var_desc_t),    intent(in)  :: v
      real(wp),            intent(out) :: scalar_out
      real(wp),            intent(out) :: slab_out(:)
      logical,             intent(out) :: valid_out(:)
      integer(ik),         intent(out) :: n_out
      real(wp), optional,  intent(in)  :: cohort_x(:)   ! a cohort field to reduce in place of the live one
                                                        ! (a staged sub-step's row, for a snapshot)
      real(wp), optional,  intent(in)  :: cohort_w(:,:) ! with cohort_x: the weights the fast loop used
                                                        ! (cohort, W_NPLANT:W_AGB), in the live order
      real(wp), optional,  intent(in)  :: soil_theta(:,:), soil_fliq(:,:)   ! a staged sub-step's soil columns
                                                        ! (n_soil_layer_max, patch), in place of the live ones
      real(wp)    :: x(max(poly%cohort%n, poly%patch%n, 1_ik))
      real(wp)    :: col(n_soil_layer_max, max(poly%patch%n, 1_ik))
      integer(ik) :: nf, kcls
      logical     :: ok
      scalar_out = 0.0_wp
      n_out = 0_ik
      kcls = src_class(v%source_id)
      !---~---
      !   A source this run does not simulate reads as missing: no slab entries, and a scalar
      ! that integrate_scalar skips, so the whole window closes as _FillValue.
      !---~---
      if (not_simulated(poly, v%source_id)) then
         scalar_out = MISSING_VALUE
         return
      end if
      !--- A cohort field reduced to a coarser axis: the live field, or the one given.
      if (kcls == SRCK_COHORT .and. v%dim /= DIM_COHORT) then
         if (present(cohort_x)) then
            nf = size(cohort_x, kind=ik)
            x(1:nf) = cohort_x
         else
            call cohort_source_field(poly, v%source_id, x, nf)
         end if
      end if

      select case (v%dim)

      !--- RAW cohort slab: the field itself, no reduction.
      case (DIM_COHORT)
         call cohort_source_field(poly, v%source_id, slab_out, n_out)
         valid_out(1:max(n_out,1_ik)) = .false.
         valid_out(1:n_out) = slab_out(1:n_out) /= MISSING_VALUE   ! a cohort with no samples: fill
         call scale_slab(slab_out, n_out, v%scale)

      !--- Patch axis: either a per-patch field verbatim, or a cohort field reduced to it.
      case (DIM_PATCH)
         if (kcls == SRCK_COHORT) then
            call reduce_cohort_to_patch(poly, x, v%weight, v%mean, slab_out, valid_out, n_out, &
                                        scale=v%scale, w_given=cohort_w)
         else
            call patch_source_field(poly, v%source_id, slab_out, n_out)
            valid_out(1:max(n_out,1_ik)) = .true.
            call scale_slab(slab_out, n_out, v%scale)
         end if

      !--- PFT axis (run-time length; empty PFTs are a true 0 for sums, see the reducer).
      case (DIM_PFT)
         call reduce_cohort_to_pft(poly, x, v%weight, v%mean, setup%n_pft, slab_out, valid_out,   &
                                   n_out, scale=v%scale, w_given=cohort_w)

      !--- DBH size class: each cohort whole in the class of its mean dbh.
      case (DIM_SIZE)
         call reduce_cohort_to_size(poly, x, v%weight, v%mean, setup%dbh_edges, setup%n_dbh_class,   &
                                    slab_out, valid_out, n_out, scale=v%scale, w_given=cohort_w)

      !--- Each site's column, its patches area-weighted, at (is-1)*n_soil_layer_max + k.
      case (DIM_SOIL)
         call layer_source_field(poly, setup, v%source_id, col, ok, soil_theta, soil_fliq)
         if (ok) then
            call reduce_patch_columns_to_sites(poly, col, slab_out, valid_out, n_out)
            call scale_slab(slab_out, n_out, v%scale)
         else
            n_out = 0_ik
         end if

      !--- One value per site.
      case (DIM_SITE)
         call site_source_field(poly, v%source_id, slab_out, valid_out, n_out)
         call scale_slab(slab_out, n_out, v%scale)

      !--- 2-D (soil layer x patch), flattened at (ip-1)*n_soil_layer_max + k.
      case (DIM_SOIL_PATCH)
         call layer_source_field(poly, setup, v%source_id, col, ok, soil_theta, soil_fliq)
         if (ok) then
            call gather_patch_columns(poly, col, slab_out, valid_out, n_out)
            call scale_slab(slab_out, n_out, v%scale)
         else
            n_out = 0_ik
         end if

      !--- Polygon scalar: one of its own, or a cohort/patch field reduced to the polygon.
      case default
         select case (kcls)
         case (SRCK_COHORT)
            call reduce_cohort_to_polygon(poly, x, v%weight, v%mean, scalar_out, ok, scale=v%scale,     &
                                          w_given=cohort_w)
         case (SRCK_PATCH)
            call patch_source_field(poly, v%source_id, x, nf)
            call reduce_patch_to_polygon(poly, x, v%mean, scalar_out, ok)
            scalar_out = scalar_out * v%scale
         case default
            scalar_out = extract_scalar_source(poly, v%source_id)
            if (scalar_out /= MISSING_VALUE) scalar_out = scalar_out * v%scale
            ok = .true.
         end select
         !---~---
         !   A mean over nothing (bare ground, or only cohorts with no samples) has nothing to
         ! report this step: missing, which integrate_scalar skips, not a 0 it would average.
         !---~---
         if (.not. ok) scalar_out = MISSING_VALUE
      end select
   end subroutine extract_variable

   !---~---
   !   The polygon's ecosystem carbon [kgC/m2]: its live vegetation (leaf, fine root, wood and
   ! storage of every plant) and its soil pools, litter included, each patch weighted by its
   ! share of the polygon, as mean_veg_carbon_py and mean_soilc_total_py reduce them.
   !---~---
   pure real(wp) function ecosystem_carbon(poly) result(c)
      type(polygon_t), intent(in) :: poly
      real(wp)    :: x(max(poly%cohort%n, poly%patch%n, 1_ik)), veg, soil
      integer(ik) :: nf
      logical     :: ok
      call cohort_source_field(poly, FLD_C_VEG_CARBON, x, nf)
      call reduce_cohort_to_polygon(poly, x, W_NPLANT, .false., veg, ok)
      call patch_source_field(poly, FLD_P_SOILC_TOTAL, x, nf)
      call reduce_patch_to_polygon(poly, x, .false., soil, ok)
      c = veg + soil
   end function ecosystem_carbon

   !--- The polygon's own scalars: already polygon quantities, no reduction.
   pure real(wp) function extract_scalar_source(poly, source_id) result(val)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in) :: source_id
      select case (source_id)
      case (SRC_PY_ECO_CARBON)
         val = ecosystem_carbon(poly)
      case (SRC_PY_N_COHORT)
         val = real(poly%cohort%n, wp)
      case (SRC_PY_N_PATCH)
         val = real(poly%patch%n, wp)
      case (SRC_PY_CANOPY_HEIGHT)
         if (poly%cohort%n > 0_ik) then
            val = maxval(poly%cohort%height(1:poly%cohort%n))
         else
            val = 0.0_wp
         end if
      case (SRC_PY_WORK_STEPS)
         val = poly%work_integ_steps
      case (SRC_PY_WORK_REJ)
         val = poly%work_integ_rej
      case (SRC_PY_WORK_WATER_NSUB)
         val = poly%work_water_nsub
      case (SRC_PY_WORK_WATER_NEWTON)
         val = poly%work_water_newton
      case (SRC_PY_WORK_WATER_NFAIL)
         val = poly%work_water_nfail
      case (SRC_PY_WORK_WATER_REJECT)
         val = poly%work_water_reject
      case (SRC_PY_WORK_RK45_RESCUE)
         val = poly%work_rk45_rescue
      case (SRC_PY_WORK_CLAMP_STAGE)
         val = poly%work_clamp_stage
      case default
         if (source_id > FLD_PY_DIAG0 .and. source_id <= FLD_PY_DIAG0 + N_PY_ROWS) then
            val = MISSING_VALUE
            if (poly%diag%active) val = polygon_diag_value(poly%diag, source_id - FLD_PY_DIAG0)
         else
            val = MISSING_VALUE
         end if
      end select
   end function extract_scalar_source

   !---~---
   !   FAST tier: one scalar variable's value in a staged sub-step, a patch-block row's
   ! polygon mean or a forcing row; slab sources are read in output_integrate_fast. The
   ! FAST tier shares everything else with the coarse tiers; only this source differs, the
   ! staging, never live state (docs/science/diagnostics.md §4, How the FAST tier is sampled).
   !---~---
   !--- A variable's unit factor on a slab its reducer did not scale: every value but the fill.
   pure subroutine scale_slab(x, n, scale)
      real(wp),    intent(inout) :: x(:)
      integer(ik), intent(in)    :: n
      real(wp),    intent(in)    :: scale
      where (x(1:n) /= MISSING_VALUE) x(1:n) = x(1:n) * scale
   end subroutine scale_slab

   pure real(wp) function extract_fast_scalar(source_id, site_row, forcing) result(val)
      integer(ik), intent(in) :: source_id
      real(wp),    intent(in) :: site_row(:)   ! the sub-step's polygon mean of the patch rows (PA_*)
      real(wp),    intent(in) :: forcing(:)    ! the sub-step's forcing (PY_*)
      if (source_id > SRC_F_PA0 .and. source_id <= SRC_F_PA0 + N_PA_ROWS) then
         val = site_row(source_id - SRC_F_PA0)
      else if (source_id > SRC_F_PY0 .and. source_id <= SRC_F_PY0 + N_PY_ROWS) then
         val = forcing(source_id - SRC_F_PY0)
      else
         val = MISSING_VALUE
      end if
   end function extract_fast_scalar

end module ed_output_extract
