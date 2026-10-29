! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_registry -- the one registration list of diagnostic variables (register_*) and
! the selection that turns it into each tier's live variables: the defaults, checked by the
! annual and fast guards, the io_config lines (apply_selections), and the tiers switched off
! and the axes' bookkeeping (resolve_streams) (docs/science/diagnostics.md §5); and
! dump_io_config, the list of what a selection can name. netCDF-free. Adding a variable:
! docs/science/diagnostics.md §8.
!---~---
module ed_output_registry
   use ed_kinds,            only : wp, ik
   use ed_config,           only : ed_config_t
   use ed_diag_vars,            only : SI_SW_IN, SI_PAR_BEAM, SI_NIR_BEAM, SI_LWDOWN, SI_TAIR, SI_QAIR,           &
                                       SI_PSURF, SI_RHO_AIR
   use ed_output_params,    only : N_FREQ, N_GRP, TIER_FAST, FREQ_FAST, FREQ_DAILY,                               &
                                   FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE, GRP_STRUCTURE, GRP_CARBON, GRP_WATER,    &
                                   GRP_ENERGY, GRP_RADIATION, GRP_ECOPHYS, GRP_BIOGEOCHEM, GRP_NUMERICS,          &
                                   GRP_FORCING
   use ed_output_types,     only : var_desc_t, output_registry_t, MAX_OUTPUT_VARS, AGG_MEAN, AGG_LAST,            &
                                   AGG_TMEAN, AGG_RATE, DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT,      &
                                   DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE, XTYPE_DOUBLE, XTYPE_INT, ragged_dim
   use ed_output_reduce, only : W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA, W_AGB
   use ed_constants,         only : cm2_to_m2, kgc_2_umol
   use ed_output_extract,     only : src_class, SRCK_FAST, SRC_NONE, FLD_C_NPLANT, FLD_C_DBH, FLD_C_HEIGHT,                 &
                                     FLD_C_BASAL_AREA, FLD_C_AGB, FLD_C_LEAF_AREA, FLD_C_GROWTH_AVG, FLD_C_PFT,   &
                                     FLD_C_OWNER_PATCH, FLD_C_GLOBAL_ID, FLD_C_LAI, FLD_C_WAI,                    &
                                     FLD_C_LEAF_CARBON, FLD_C_FINEROOT_CARBON, FLD_C_WOOD_CARBON,                 &
                                     FLD_C_NONSTRUCTURAL_CARBON, FLD_C_BGB, FLD_C_VEG_CARBON, FLD_C_SLA,          &
                                     FLD_C_VCMAX25, FLD_C_RD25, FLD_C_LEAF_TURNOVER_RATE,                         &
                                     FLD_C_OVERTOPPING_LAI, FLD_C_GPP_ACCUM, FLD_C_NPP_ACCUM, FLD_C_LEAF_RESP,    &
                                     FLD_C_STEM_RESP, FLD_C_ROOT_RESP, FLD_C_DMAX_PSI_LEAF,                       &
                                     FLD_C_LEAF_FLUSH_TENDENCY, FLD_C_LEAF_SHED_TENDENCY, FLD_P_AREA,             &
                                     FLD_P_AGE, FLD_P_DIST_TYPE, FLD_P_COHORT_OFFSET, FLD_P_COHORT_COUNT,         &
                                     FLD_P_GLOBAL_ID, FLD_P_SITE, FLD_P_CAS_DEPTH, FLD_P_SWE, FLD_P_SNOW_DEPTH,   &
                                     FLD_P_SOILC_FAST_GRND, FLD_P_SOILC_FAST_SOIL, FLD_P_SOILC_STRUCT_GRND,       &
                                     FLD_P_SOILC_STRUCT_SOIL, FLD_P_SOILC_MICROBIAL, FLD_P_SOILC_SLOW,            &
                                     FLD_P_SOILC_PASSIVE, FLD_P_SOILC_TOTAL, FLD_P_RH, FLD_L_SOIL_TEMP,           &
                                     FLD_L_SOIL_WATER, FLD_L_SOIL_PSI, FLD_L_SOIL_WETNESS, FLD_L_SOIL_FLIQ,       &
                                     SRC_PY_N_COHORT, SRC_PY_N_PATCH, SRC_PY_CANOPY_HEIGHT, SRC_PY_ECO_CARBON,     &
                                     SRC_PY_WORK_STEPS, SRC_PY_WORK_REJ, SRC_PY_WORK_WATER_NSUB,                   &
                                     SRC_PY_WORK_WATER_NEWTON, SRC_PY_WORK_WATER_NFAIL, SRC_PY_WORK_WATER_REJECT,  &
                                     SRC_PY_WORK_RK45_RESCUE, SRC_PY_WORK_CLAMP_STAGE, SRC_F_SOIL_TEMP,            &
                                     SRC_F_SOIL_WATER, SRC_F_SOIL_LAYER, SRC_F_PY0, SRC_F_PA0, SRC_F_COH_GPP,      &
                                     SRC_F_COH_HEIGHT, SRC_F_COH_NPP, SRC_F_CO0, FLD_C_DIAG0, FLD_P_DIAG0,         &
                                     FLD_PY_DIAG0, FLD_SI_DIAG0, FLD_SI_GROWTH_TEMP
   use ed_diag_vars,         only : N_CO_FAST_ROWS, N_PA_FAST_ROWS, N_PY_ROWS, N_SI_ROWS,                    &
                                    CO_ANET, CO_AGROSS, CO_GSW, CO_GBW, CO_CI, CO_CS, CO_RD, CO_TRANSP,           &
                                    CO_BETA_STOM, CO_BETA_NONSTOM, CO_LEAF_VPD, CO_PSI_LEAF, CO_PSI_WOOD,         &
                                    CO_PLC, CO_SAPFLOW, CO_ROOT_UPTAKE, CO_ABS_PAR, CO_ABS_SW, CO_ABS_LW,         &
                                    CO_WIND, CO_LEAF_WATER, CO_WOOD_WATER, CO_LEAF_TEMP, CO_WOOD_TEMP,            &
                                    CO_LEAF_PAR, CO_GB_MOL, CO_PSI_PREDAWN, CO_CA, CO_PRESSURE, PA_LE, PA_H,      &
                                    PA_RNET, PA_SW_GROUND, PA_LW_GROUND, PA_SW_IN_VIS, PA_SW_IN_NIR,              &
                                    PA_SW_UP_VIS, PA_SW_UP_NIR, PA_LW_UP, PA_USTAR, PA_GGNET, PA_ROUGH,           &
                                    PA_DISPLACE, PA_GPP, PA_NEE, PA_TRANSP, PA_SOIL_TEMP_TOP, PA_RESID_ENERGY,    &
                                    PA_EVAP_CANOPY, PA_EVAP_GROUND, PA_DEW,                                        &
                                    PA_CAS_VPD, PA_W_SURFACE, PA_CAS_TEMP, PA_CAS_SHV, PA_CAS_CO2,                &
                                    PA_WIND_CAS_TOP, PA_TAIR_CAS_TOP, PA_Z_CAS_TOP, PA_RESID_WATER, PA_NPP,       &
                                    PA_RECO, PA_SW_UP, PA_SKIN_TEMP, PY_SW_IN, PY_PRECIP, PY_TAIR, PY_QAIR,       &
                                    PY_PSURF, PY_WIND, PY_LWDOWN, PY_PAR_BEAM, PY_PAR_DIFFUSE, PY_NIR_BEAM,       &
                                    PY_NIR_DIFFUSE, PY_SNOWFALL, PY_CO2, PY_COSZ, PY_RHO_AIR, PY_RAINF,           &
                                    CO_DDBH_DT, CO_DAGB_DT, CO_MORT_RATE, CO_NPP_LEAF, CO_NPP_FINEROOT,           &
                                    CO_NPP_WOOD, CO_NPP_STORAGE, CO_NPP_REPRO, CO_GROWTH_RESP,                    &
                                    CO_EXUDATE, PA_LITTER_LEAF, PA_LITTER_FINEROOT, PA_MORT_C_BACKGROUND,         &
                                    PA_MORT_C_CULL, PA_MORT_C_DISTURB, PA_LITTER_STRUCT, PA_LITTER_WOOD,          &
                                    PA_LITTER_REPRO, PA_RECRUIT_NPLANT, PA_DISTURB_AREA, PA_APAR, PA_H_CANOPY,    &
                                    PA_H_GROUND, PA_HEAT_BOTTOM, PA_ENTH_ADVECT, PA_HEAT_GAIN_SOIL,               &
                                    PA_HEAT_GAIN_CANOPY, PA_HEAT_GAIN_SURFACE, PA_RA, PA_LEAF_RESP,               &
                                    PA_STEM_RESP, PA_ROOT_RESP, PA_RH, PA_ET, PA_RUNOFF_SURF, PA_DRAINAGE,        &
                                    PA_RUNOFF, PA_INTERCEPT, PA_THROUGHFALL, PA_INFILTRATION, PA_SNOWMELT,        &
                                    PA_SUBLIMATION, PA_SHED_WATER, PA_FILM_EXCESS, PA_WATER_STORE_GAIN,           &
                                    PA_CANOPY_WATER, PA_SOIL_WATER_TOTAL, PA_PLANT_WATER, CO_WET_FRACTION
   implicit none
   private

   public :: build_output_registry, build_freq_index, find_var_index, parse_stream_mask
   public :: selection_t, apply_selections, apply_variable_override, resolve_streams, usable_tiers
   public :: glob_match, freq_bit, OVR_TRUE, OVR_FALSE, OVR_MASK
   public :: dump_io_config

   !---~---
   !   Named default stream masks. A FAST bit is written only with [output.fast].enabled, which
   ! is off by default.
   !---~---
   integer(ik), parameter :: DAY_MON    = ior(FREQ_DAILY, FREQ_MONTHLY)
   integer(ik), parameter :: MON        = FREQ_MONTHLY
   integer(ik), parameter :: MON_YR     = ior(FREQ_MONTHLY, FREQ_ANNUAL)
   integer(ik), parameter :: DAY_MON_YR = ior(ior(FREQ_DAILY, FREQ_MONTHLY), FREQ_ANNUAL)
   integer(ik), parameter :: FAST       = FREQ_FAST

   !---~---
   !   Per-variable override kinds (ed_io_config.toml value grammar,
   ! docs/science/diagnostics.md §5).
   !---~---
   integer(ik), parameter :: OVR_TRUE  = 1_ik  ! bool true  -> restore registry-default streams
   integer(ik), parameter :: OVR_FALSE = 2_ik  ! bool false -> disable everywhere
   integer(ik), parameter :: OVR_MASK  = 3_ik  ! stream-string -> replace the stream mask

   !---~---
   !   One io_config line: a variable's name or a pattern of names (`*` any run of characters),
   ! and what it selects: OVR_TRUE, OVR_FALSE or OVR_MASK with its tiers.
   !---~---
   type :: selection_t
      character(len=64) :: key  = ''
      integer(ik)       :: kind = OVR_FALSE
      integer(ik)       :: mask = FREQ_NONE
   end type selection_t

contains

   !--- FREQ_* bit for a 1-based tier index (TIER_FAST .. TIER_ANNUAL).
   pure integer(ik) function freq_bit(tier) result(b)
      integer(ik), intent(in) :: tier
      b = ishft(1_ik, tier - 1_ik)
   end function freq_bit

   !--- Build the registry: the source list at its default tiers, less the tiers switched off.
   subroutine build_output_registry(reg, cfg)
      type(output_registry_t), intent(out) :: reg
      type(ed_config_t),     intent(in)    :: cfg

      allocate(reg%var(MAX_OUTPUT_VARS))
      reg%nvar = 0_ik
      allocate(reg%idx_freq(MAX_OUTPUT_VARS, N_FREQ))
      reg%idx_freq = 0_ik
      reg%nidx = 0_ik

      call register_structure_cohort(reg)
      call register_structure_patch(reg)
      call register_structure_site(reg)
      call register_structure_pft(reg)
      call register_structure_size(reg)
      call register_carbon(reg)
      call register_water(reg)
      call register_energy(reg)
      call register_biogeochem(reg)
      call register_numerics(reg, cfg%dt_fast)
      call register_ecophys(reg)
      call register_allocation(reg)
      call register_patch_fluxes(reg)
      call register_budgets(reg)
      call register_fast(reg)
      call register_forcing(reg)
      call register_snapshots(reg)                ! each mean_<x>'s end-of-period <x>, off by default

      call enforce_annual_guard(reg)              ! no cohort or patch variable on FREQ_ANNUAL (docs/science/diagnostics.md §3)
      call enforce_fast_guard(reg)                ! FAST tier <-> FAST-staged source, both ways
      reg%tier_on = cfg%output%freq_on            ! [output.<tier>].enabled: a tier off writes nothing
      call resolve_streams(reg)                   ! the axes' bookkeeping, each tier's live list
   end subroutine build_output_registry

   !---~---
   !   The registration list, one routine per (group, axis). A variable is one add_variable
   ! line, plus a field case in ed_output_integrate for a new field
   ! (docs/science/diagnostics.md §8):
   !
   !   add_variable(reg, name, long_name, units, dim, agg, group, streams, source
   !              [, xtype][, weight][, mean][, scale][, fast_source][, flux][, coord])
   !
   !   streams are the default tiers: FREQ_NONE for a variable written only when a selection
   ! names it, as the radiation and ecophys groups, the soil-by-patch columns and the
   ! snapshots are.
   !
   !   A variable whose name has no mean_ prefix is a snapshot (AGG_LAST); a state's snapshot
   ! is registered here with its mean, and register_snapshots adds the others. A flux or a
   ! rate (flux=, implied by AGG_RATE) is a mean only.
   !---~---

   !--- Per-cohort structure: the demographic atom's own state, written raw.
   subroutine register_structure_cohort(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'nplant_co', 'plant number density', 'plant/m2', DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON,     &
                        FLD_C_NPLANT)
      call add_variable(reg, 'mean_nplant_co', 'plant number density', 'plant/m2', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE,         &
                        FREQ_NONE, FLD_C_NPLANT)
      call add_variable(reg, 'mean_dbh_co', 'diameter at breast height', 'cm', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON,    &
                        FLD_C_DBH)
      call add_variable(reg, 'mean_height_co', 'height', 'm', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, ior(FAST, DAY_MON),          &
                        FLD_C_HEIGHT, fast_source=SRC_F_COH_HEIGHT)
      call add_variable(reg, 'mean_basal_area_co', 'basal area per plant', 'cm2/plant', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE,    &
                        MON, FLD_C_BASAL_AREA)
      call add_variable(reg, 'mean_agb_co', 'aboveground biomass per plant', 'kgC/plant', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE,  &
                        DAY_MON, FLD_C_AGB)
      call add_variable(reg, 'mean_leaf_area_co', 'leaf area per plant', 'm2/plant', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON,  &
                        FLD_C_LEAF_AREA)
      call add_variable(reg, 'mean_lai_co', 'cohort leaf area index', 'm2/m2', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON,        &
                        FLD_C_LAI)
      call add_variable(reg, 'mean_wai_co', 'cohort wood (stem and branch) area index', 'm2/m2', DIM_COHORT, AGG_MEAN,          &
                        GRP_STRUCTURE, MON, FLD_C_WAI)
      call add_variable(reg, 'mean_growth_avg_co', 'moving-average growth', 'cm/yr', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON,  &
                        FLD_C_GROWTH_AVG)
      call add_variable(reg, 'mean_overtopping_lai_co', 'cumulative LAI of taller cohorts', 'm2/m2', DIM_COHORT, AGG_MEAN,      &
                        GRP_STRUCTURE, MON, FLD_C_OVERTOPPING_LAI)
      call add_variable(reg, 'pft_co', 'plant functional type index', '-', DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON,        &
                        FLD_C_PFT, xtype=XTYPE_INT)
      call add_variable(reg, 'owner_patch', 'owning patch (1-based)', '-', DIM_COHORT, AGG_LAST, GRP_STRUCTURE, FREQ_NONE,      &
                        FLD_C_OWNER_PATCH, xtype=XTYPE_INT, coord=.true.)
      call add_variable(reg, 'global_cohort_id', 'persistent cohort id', '-', DIM_COHORT, AGG_LAST, GRP_STRUCTURE, FREQ_NONE,   &
                        FLD_C_GLOBAL_ID, xtype=XTYPE_INT, coord=.true.)
      !--- Per-plant carbon pools (the four pools the size anchor rides on).
      call add_variable(reg, 'mean_leaf_carbon_co', 'leaf carbon per plant', 'kgC/plant', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE,  &
                        MON, FLD_C_LEAF_CARBON)
      call add_variable(reg, 'mean_fineroot_carbon_co', 'fine-root carbon per plant', 'kgC/plant', DIM_COHORT, AGG_MEAN,        &
                        GRP_STRUCTURE, MON, FLD_C_FINEROOT_CARBON)
      call add_variable(reg, 'mean_wood_carbon_co', 'wood carbon per plant', 'kgC/plant', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE,  &
                        MON, FLD_C_WOOD_CARBON)
      call add_variable(reg, 'mean_storage_carbon_co', 'non-structural carbon per plant', 'kgC/plant', DIM_COHORT, AGG_MEAN,    &
                        GRP_STRUCTURE, MON, FLD_C_NONSTRUCTURAL_CARBON)
      !--- Dynamic (light-acclimating) leaf traits.
      call add_variable(reg, 'mean_sla_co', 'specific leaf area', 'm2/kgC', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_SLA)
      call add_variable(reg, 'mean_vcmax25_co', 'max carboxylation rate at 25 degC', 'umol/m2/s', DIM_COHORT, AGG_MEAN,         &
                        GRP_STRUCTURE, MON, FLD_C_VCMAX25)
      call add_variable(reg, 'mean_rd25_co', 'leaf dark respiration at 25 degC', 'umol/m2/s', DIM_COHORT, AGG_MEAN,             &
                        GRP_STRUCTURE, MON, FLD_C_RD25)
      call add_variable(reg, 'mean_leaf_turnover_rate_co', 'baseline leaf turnover rate', '1/yr', DIM_COHORT, AGG_MEAN,         &
                        GRP_STRUCTURE, MON, FLD_C_LEAF_TURNOVER_RATE)
      call add_variable(reg, 'mean_leaf_flush_tendency_co', 'smoothed leaf flush tendency', '-', DIM_COHORT, AGG_MEAN,          &
                        GRP_STRUCTURE, MON, FLD_C_LEAF_FLUSH_TENDENCY)
      call add_variable(reg, 'mean_leaf_shed_tendency_co', 'smoothed leaf senescence tendency', '-', DIM_COHORT, AGG_MEAN,      &
                        GRP_STRUCTURE, MON, FLD_C_LEAF_SHED_TENDENCY)
      call add_variable(reg, 'mean_dmax_psi_leaf_co', 'predawn (daily-max) leaf water potential', 'MPa', DIM_COHORT, AGG_MEAN,  &
                        GRP_STRUCTURE, DAY_MON, FLD_C_DMAX_PSI_LEAF)
      !--- Slow-loop rates off the per-step diagnostic block (poly%cohort%sdiag).
      call add_variable(reg, 'mean_dbh_growth_co', 'diameter growth rate', 'cm/yr', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON,   &
                        FLD_C_DIAG0 + CO_DDBH_DT, flux=.true.)
      call add_variable(reg, 'mean_agb_growth_co', 'per-plant AGB growth rate', 'kgC/plant/yr', DIM_COHORT, AGG_MEAN,           &
                        GRP_STRUCTURE, MON, FLD_C_DIAG0 + CO_DAGB_DT, flux=.true.)
      call add_variable(reg, 'mean_mort_rate_co', 'mortality hazard rate', '1/yr', DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON,    &
                        FLD_C_DIAG0 + CO_MORT_RATE, flux=.true.)
   end subroutine register_structure_cohort

   !--- Per-patch structure: geometry + identity, plus cohort fields reduced to the patch.
   subroutine register_structure_patch(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'area_pa', 'patch area fraction', '-', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_P_AREA)
      call add_variable(reg, 'mean_area_pa', 'patch area fraction', '-', DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, FREQ_NONE,         &
                        FLD_P_AREA)
      call add_variable(reg, 'age_pa', 'time since disturbance', 'yr', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, MON, FLD_P_AGE)
      call add_variable(reg, 'mean_age_pa', 'time since disturbance', 'yr', DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, FREQ_NONE,      &
                        FLD_P_AGE)
      call add_variable(reg, 'dist_type_pa', 'disturbance type', '-', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, MON, FLD_P_DIST_TYPE, &
                        xtype=XTYPE_INT)
      call add_variable(reg, 'cohort_offset', 'first cohort of patch (CSR)', '-', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, FREQ_NONE,&
                        FLD_P_COHORT_OFFSET, xtype=XTYPE_INT, coord=.true.)
      call add_variable(reg, 'cohort_count', 'cohorts in patch (CSR)', '-', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, FREQ_NONE,      &
                        FLD_P_COHORT_COUNT, xtype=XTYPE_INT, coord=.true.)
      call add_variable(reg, 'global_patch_id', 'persistent patch id', '-', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, FREQ_NONE,      &
                        FLD_P_GLOBAL_ID, xtype=XTYPE_INT, coord=.true.)
      call add_variable(reg, 'patch_site', 'the patch''s site (1-based)', '-', DIM_PATCH, AGG_LAST, GRP_STRUCTURE, FREQ_NONE,   &
                        FLD_P_SITE, xtype=XTYPE_INT, coord=.true.)
      !---~---
      !   Cohort fields reduced onto the patch axis: weighted sums per m2 of that patch's own
      ! ground, the gap-versus-closed-canopy contrast the polygon scalars average away.
      !---~---
      call add_variable(reg, 'mean_nplant_pa', 'patch plant number density', 'plant/m2', DIM_PATCH, AGG_MEAN, GRP_STRUCTURE,    &
                        DAY_MON, FLD_C_NPLANT, weight=W_NONE)
      call add_variable(reg, 'mean_lai_pa', 'patch leaf area index', 'm2/m2', DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, DAY_MON,      &
                        FLD_C_LEAF_AREA, weight=W_NPLANT)
      call add_variable(reg, 'mean_agb_pa', 'patch aboveground biomass', 'kgC/m2', DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, DAY_MON, &
                        FLD_C_AGB, weight=W_NPLANT)
      call add_variable(reg, 'mean_basal_area_pa', 'patch basal area', 'm2/m2', DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, MON,        &
                        FLD_C_BASAL_AREA, weight=W_NPLANT, scale=cm2_to_m2)
   end subroutine register_structure_patch

   !--- The polygon's structure scalars.
   subroutine register_structure_site(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_nplant_py', 'polygon total plant number', 'plant/m2', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE,      &
                        DAY_MON_YR, FLD_C_NPLANT, weight=W_NONE)
      call add_variable(reg, 'mean_basal_area_py', 'polygon total basal area', 'm2/m2', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE,       &
                        DAY_MON_YR, FLD_C_BASAL_AREA, weight=W_NPLANT, scale=cm2_to_m2)
      call add_variable(reg, 'mean_agb_py', 'polygon aboveground biomass', 'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE,          &
                        DAY_MON_YR, FLD_C_AGB, weight=W_NPLANT)
      call add_variable(reg, 'mean_lai_py', 'polygon leaf area index', 'm2/m2', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR,   &
                        FLD_C_LEAF_AREA, weight=W_NPLANT)
      !---~---
      !   Intensive structural means, basal-area-weighted: dominated by the trees that hold
      ! the stand, not by a regenerating understory (docs/science/diagnostics.md §2).
      !---~---
      call add_variable(reg, 'mean_dbh_py', 'basal-area-weighted mean DBH', 'cm', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE,          &
                        DAY_MON_YR, FLD_C_DBH, weight=W_BASAL_AREA, mean=.true.)
      call add_variable(reg, 'mean_height_py', 'basal-area-weighted mean height', 'm', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE,     &
                        DAY_MON_YR, FLD_C_HEIGHT, weight=W_BASAL_AREA, mean=.true.)
      call add_variable(reg, 'mean_canopy_height_py', 'height of the tallest cohort', 'm', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, &
                        DAY_MON_YR, SRC_PY_CANOPY_HEIGHT)
      call add_variable(reg, 'n_cohort_py', 'live cohort count', '-', DIM_SCALAR, AGG_LAST, GRP_STRUCTURE, DAY_MON_YR,          &
                        SRC_PY_N_COHORT, xtype=XTYPE_INT)
      call add_variable(reg, 'n_patch_py', 'live patch count', '-', DIM_SCALAR, AGG_LAST, GRP_STRUCTURE, DAY_MON_YR,            &
                        SRC_PY_N_PATCH, xtype=XTYPE_INT)
      call add_variable(reg, 'mean_n_cohort_py', 'live cohort count', '-', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, FREQ_NONE,      &
                        SRC_PY_N_COHORT)
      call add_variable(reg, 'mean_n_patch_py', 'live patch count', '-', DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, FREQ_NONE,        &
                        SRC_PY_N_PATCH)
   end subroutine register_structure_site

   !--- PFT axis. The composition view: is any PFT actually competing, or has one won?
   subroutine register_structure_pft(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_nplant_pft_py', 'plant number density by PFT', 'plant/m2', DIM_PFT, AGG_MEAN, GRP_STRUCTURE, &
                        MON_YR, FLD_C_NPLANT, weight=W_NONE)
      call add_variable(reg, 'mean_agb_pft_py', 'aboveground biomass by PFT', 'kgC/m2', DIM_PFT, AGG_MEAN, GRP_STRUCTURE,       &
                        MON_YR, FLD_C_AGB, weight=W_NPLANT)
      call add_variable(reg, 'mean_lai_pft_py', 'leaf area index by PFT', 'm2/m2', DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR,    &
                        FLD_C_LEAF_AREA, weight=W_NPLANT)
      call add_variable(reg, 'mean_basal_area_pft_py', 'basal area by PFT', 'm2/m2', DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR,  &
                        FLD_C_BASAL_AREA, weight=W_NPLANT, scale=cm2_to_m2)
      call add_variable(reg, 'mean_dbh_pft_py', 'basal-area-weighted mean DBH by PFT', 'cm', DIM_PFT, AGG_MEAN, GRP_STRUCTURE,  &
                        MON_YR, FLD_C_DBH, weight=W_BASAL_AREA, mean=.true.)
   end subroutine register_structure_pft

   !---~---
   !   DBH size-class axis: a stem-density distribution comparable to a forest inventory
   ! (docs/science/diagnostics.md §3).
   !---~---
   subroutine register_structure_size(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_nplant_size_py', 'plant number density by DBH class', 'plant/m2', DIM_SIZE, AGG_MEAN,        &
                        GRP_STRUCTURE, MON_YR, FLD_C_NPLANT, weight=W_NONE)
      call add_variable(reg, 'mean_agb_size_py', 'aboveground biomass by DBH class', 'kgC/m2', DIM_SIZE, AGG_MEAN,              &
                        GRP_STRUCTURE, MON_YR, FLD_C_AGB, weight=W_NPLANT)
      call add_variable(reg, 'mean_lai_size_py', 'leaf area index by DBH class', 'm2/m2', DIM_SIZE, AGG_MEAN, GRP_STRUCTURE,    &
                        MON_YR, FLD_C_LEAF_AREA, weight=W_NPLANT)
      call add_variable(reg, 'mean_basal_area_size_py', 'basal area by DBH class', 'm2/m2', DIM_SIZE, AGG_MEAN, GRP_STRUCTURE,  &
                        MON_YR, FLD_C_BASAL_AREA, weight=W_NPLANT, scale=cm2_to_m2)
      call add_variable(reg, 'mean_agb_growth_size_py', 'AGB growth rate by DBH class', 'kgC/m2/yr', DIM_SIZE, AGG_MEAN,        &
                        GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_DAGB_DT, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_mort_rate_size_py', 'nplant-weighted mortality rate by DBH class', '1/yr', DIM_SIZE,         &
                        AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, weight=W_NPLANT, mean=.true., flux=.true.)
   end subroutine register_structure_size

   !---~---
   !   Carbon. gpp_accum and the *_resp_accum are per-slow-step amounts the fast loop fills
   ! [kgC]; AGG_RATE reports their mean rate over the period, in umol CO2 (all 0 when the
   ! fast loop is off).
   !---~---
   subroutine register_carbon(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_gpp_py', 'polygon gross primary productivity', 'umol/m2/s', DIM_SCALAR, AGG_RATE, GRP_CARBON,   &
                        ior(FAST, DAY_MON_YR), FLD_C_GPP_ACCUM, weight=W_NPLANT, scale=kgc_2_umol,                              &
                        fast_source=SRC_F_PA0 + PA_GPP)
      call add_variable(reg, 'mean_npp_py', 'polygon net primary productivity', 'umol/m2/s', DIM_SCALAR, AGG_RATE, GRP_CARBON,     &
                        ior(FAST, DAY_MON_YR), FLD_C_NPP_ACCUM, weight=W_NPLANT, scale=kgc_2_umol,                              &
                        fast_source=SRC_F_PA0 + PA_NPP)
      call add_variable(reg, 'mean_leaf_resp_py', 'polygon leaf dark respiration', 'umol/m2/s', DIM_SCALAR, AGG_RATE, GRP_CARBON,  &
                        DAY_MON_YR, FLD_C_LEAF_RESP, weight=W_NPLANT, scale=kgc_2_umol,      &
                        fast_source=SRC_F_PA0 + PA_LEAF_RESP)
      call add_variable(reg, 'mean_stem_resp_py', 'polygon stem maintenance respiration', 'umol/m2/s', DIM_SCALAR, AGG_RATE,       &
                        GRP_CARBON, DAY_MON_YR, FLD_C_STEM_RESP, weight=W_NPLANT, scale=kgc_2_umol,      &
                        fast_source=SRC_F_PA0 + PA_STEM_RESP)
      call add_variable(reg, 'mean_root_resp_py', 'polygon fine-root maintenance respiration', 'umol/m2/s', DIM_SCALAR, AGG_RATE,  &
                        GRP_CARBON, DAY_MON_YR, FLD_C_ROOT_RESP, weight=W_NPLANT, scale=kgc_2_umol,      &
                        fast_source=SRC_F_PA0 + PA_ROOT_RESP)
      !--- Vegetation carbon stocks (AGG_MEAN: a stock, not a flux).
      call add_variable(reg, 'mean_leaf_carbon_py', 'polygon leaf carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, &
                        FLD_C_LEAF_CARBON, weight=W_NPLANT)
      call add_variable(reg, 'mean_fineroot_carbon_py', 'polygon fine-root carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_CARBON,    &
                        DAY_MON_YR, FLD_C_FINEROOT_CARBON, weight=W_NPLANT)
      call add_variable(reg, 'mean_wood_carbon_py', 'polygon wood carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, &
                        FLD_C_WOOD_CARBON, weight=W_NPLANT)
      call add_variable(reg, 'mean_storage_carbon_py', 'polygon non-structural carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN,         &
                        GRP_CARBON, DAY_MON_YR, FLD_C_NONSTRUCTURAL_CARBON, weight=W_NPLANT)
      call add_variable(reg, 'mean_bgb_py', 'polygon belowground biomass', 'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, &
                        FLD_C_BGB, weight=W_NPLANT)
      call add_variable(reg, 'mean_veg_carbon_py', 'polygon total live vegetation carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN,         &
                        GRP_CARBON, DAY_MON_YR, FLD_C_VEG_CARBON, weight=W_NPLANT)
      !--- PFT-resolved carbon fluxes.
      call add_variable(reg, 'mean_gpp_pft_py', 'gross primary productivity by PFT', 'umol/m2/s', DIM_PFT, AGG_RATE, GRP_CARBON, &
                        MON_YR, FLD_C_GPP_ACCUM, weight=W_NPLANT, scale=kgc_2_umol)
      call add_variable(reg, 'mean_npp_pft_py', 'net primary productivity by PFT', 'umol/m2/s', DIM_PFT, AGG_RATE, GRP_CARBON,  &
                        MON_YR, FLD_C_NPP_ACCUM, weight=W_NPLANT, scale=kgc_2_umol)
      !--- Per-patch carbon fluxes.
      call add_variable(reg, 'mean_gpp_pa', 'patch gross primary productivity', 'umol/m2/s', DIM_PATCH, AGG_RATE, GRP_CARBON,   &
                        ior(FAST, DAY_MON), FLD_C_GPP_ACCUM, weight=W_NPLANT, scale=kgc_2_umol, fast_source=SRC_F_PA0 + PA_GPP)
      call add_variable(reg, 'mean_npp_pa', 'patch net primary productivity', 'umol/m2/s', DIM_PATCH, AGG_RATE, GRP_CARBON,     &
                        ior(FAST, DAY_MON), FLD_C_NPP_ACCUM, weight=W_NPLANT, scale=kgc_2_umol, fast_source=SRC_F_PA0 + PA_NPP)
      !--- Per-cohort carbon fluxes (per plant).
      call add_variable(reg, 'mean_gpp_co', 'per-plant GPP', 'umol/plant/s', DIM_COHORT, AGG_RATE, GRP_CARBON,                  &
                        ior(FAST, DAY_MON), FLD_C_GPP_ACCUM, scale=kgc_2_umol, fast_source=SRC_F_COH_GPP)
      call add_variable(reg, 'mean_npp_co', 'per-plant NPP', 'umol/plant/s', DIM_COHORT, AGG_RATE, GRP_CARBON, DAY_MON,         &
                        FLD_C_NPP_ACCUM, scale=kgc_2_umol, fast_source=SRC_F_COH_NPP)
   end subroutine register_carbon

   !--- Water.
   subroutine register_water(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_et_py', 'polygon evapotranspiration', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_WATER,          &
                        DAY_MON_YR, FLD_P_DIAG0 + PA_ET, flux=.true.)
      call add_variable(reg, 'mean_et_pa', 'patch evapotranspiration', 'kg/m2/s', DIM_PATCH, AGG_TMEAN, GRP_WATER,             &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_ET, flux=.true.)
      !---~---
      !   The water leaving the soil column: surface runoff (the pond's spill above its capacity)
      ! and drainage (the bottom flux: free drainage, 0 on bedrock, net and signed with an
      ! aquifer), and their sum.
      !---~---
      call add_variable(reg, 'mean_runoff_py', 'polygon runoff: surface runoff + drainage', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN,      &
                        GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_RUNOFF, flux=.true.)
      call add_variable(reg, 'mean_runoff_surf_py', 'polygon surface runoff', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_WATER,         &
                        DAY_MON_YR, FLD_P_DIAG0 + PA_RUNOFF_SURF, flux=.true.)
      call add_variable(reg, 'mean_drainage_py', 'polygon drainage at the soil column''s bottom', 'kg/m2/s', DIM_SCALAR,          &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_DRAINAGE, flux=.true.)
      !--- Each site's soil column, its patches area-weighted (DIM_SOIL).
      call add_variable(reg, 'mean_soil_water_layer_si', 'each site''s area-weighted volumetric soil moisture', 'm3/m3',        &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, ior(FAST, DAY_MON_YR), FLD_L_SOIL_WATER, fast_source=SRC_F_SOIL_WATER)
      call add_variable(reg, 'mean_soil_psi_layer_si', 'each site''s area-weighted soil matric potential', 'MPa', DIM_SOIL,     &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_PSI,                    &
                        fast_source=SRC_F_SOIL_LAYER)
      call add_variable(reg, 'mean_soil_wetness_layer_si', 'each site''s area-weighted relative saturation', '-', DIM_SOIL,     &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_WETNESS,                    &
                        fast_source=SRC_F_SOIL_LAYER)
      !--- Surface water stores.
      call add_variable(reg, 'mean_w_surface_py', 'ponded surface water', 'kg/m2', DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, &
                        FLD_P_DIAG0 + PA_W_SURFACE)
      call add_variable(reg, 'mean_swe_py', 'snow water equivalent', 'kg/m2', DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR,     &
                        FLD_P_SWE)
      call add_variable(reg, 'mean_snow_depth_py', 'snow depth', 'm', DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR,             &
                        FLD_P_SNOW_DEPTH)
      call add_variable(reg, 'mean_w_surface_pa', 'patch ponded surface water', 'kg/m2', DIM_PATCH, AGG_TMEAN, GRP_WATER, MON,  &
                        FLD_P_DIAG0 + PA_W_SURFACE)
      call add_variable(reg, 'mean_swe_pa', 'patch snow water equivalent', 'kg/m2', DIM_PATCH, AGG_TMEAN, GRP_WATER, MON,       &
                        FLD_P_SWE)
      !---~---
      !   2-D (soil layer x patch) columns, off by default (the largest after the cohorts): the
      ! gap-vs-closed-canopy drydown contrast the area-weighted site column hides.
      !---~---
      call add_variable(reg, 'mean_soil_water_layer_pa', 'volumetric soil moisture by layer and patch', 'm3/m3', DIM_SOIL_PATCH, &
                        AGG_TMEAN, GRP_WATER, FREQ_NONE, FLD_L_SOIL_WATER, fast_source=SRC_F_SOIL_WATER)
      call add_variable(reg, 'mean_soil_psi_layer_pa', 'soil matric potential by layer and patch', 'MPa', DIM_SOIL_PATCH,       &
                        AGG_TMEAN, GRP_WATER, FREQ_NONE, FLD_L_SOIL_PSI,                    &
                        fast_source=SRC_F_SOIL_LAYER)
      call add_variable(reg, 'mean_soil_wetness_layer_pa', 'relative saturation by layer and patch', '-', DIM_SOIL_PATCH,       &
                        AGG_TMEAN, GRP_WATER, FREQ_NONE, FLD_L_SOIL_WETNESS,                    &
                        fast_source=SRC_F_SOIL_LAYER)
   end subroutine register_water

   !--- Energy.
   subroutine register_energy(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_cas_temp_py', 'polygon canopy-air-space temperature', 'K', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY,   &
                        ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_CAS_TEMP)
      call add_variable(reg, 'mean_cas_shv_py', 'polygon canopy-air specific humidity', 'kg/kg', DIM_SCALAR, AGG_TMEAN,         &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_CAS_SHV)
      call add_variable(reg, 'mean_cas_co2_py', 'polygon canopy-air CO2 mixing ratio', 'umol/mol', DIM_SCALAR, AGG_TMEAN,          &
                        GRP_ENERGY, ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_CAS_CO2)
      call add_variable(reg, 'mean_cas_vpd_py', 'polygon canopy-air vapour-pressure deficit', 'Pa', DIM_SCALAR, AGG_TMEAN,         &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_CAS_VPD)
      call add_variable(reg, 'mean_cas_depth_py', 'polygon canopy-air-space depth', 'm', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY,     &
                        MON_YR, FLD_P_CAS_DEPTH)
      call add_variable(reg, 'mean_soil_temp_top_py', 'polygon soil-top (layer 1) temperature', 'K', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_ENERGY, ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_SOIL_TEMP_TOP)
      !---~---
      !   The skin temperature (docs/science/diagnostics.md §4): what an infrared thermometer
      ! or a satellite land-surface temperature sees, not the soil.
      !---~---
      call add_variable(reg, 'mean_skin_temp_py', 'skin temperature: black-body temperature of the longwave up', 'K',           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_SKIN_TEMP)
      call add_variable(reg, 'mean_soil_temp_layer_si', 'each site''s area-weighted soil temperature', 'K', DIM_SOIL, AGG_TMEAN, &
                        GRP_ENERGY, ior(FAST, DAY_MON_YR), FLD_L_SOIL_TEMP, fast_source=SRC_F_SOIL_TEMP)
      call add_variable(reg, 'mean_soil_fliq_layer_si', 'each site''s area-weighted soil liquid fraction', '-', DIM_SOIL,       &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_L_SOIL_FLIQ,                    &
                        fast_source=SRC_F_SOIL_LAYER)

      !---~---
      !   Canopy temperatures, leaf-area-weighted: a bare sapling does not pull the canopy
      ! mean as hard as a closed overstory.
      !---~---
      call add_variable(reg, 'mean_leaf_temp_py', 'leaf-area-weighted canopy leaf temperature', 'K', DIM_SCALAR, AGG_TMEAN,     &
                        GRP_ENERGY, DAY_MON_YR, FLD_C_DIAG0 + CO_LEAF_TEMP, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_wood_temp_py', 'leaf-area-weighted wood temperature', 'K', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, &
                        DAY_MON_YR, FLD_C_DIAG0 + CO_WOOD_TEMP, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_leaf_temp_co', 'per-cohort leaf temperature', 'K', DIM_COHORT, AGG_TMEAN, GRP_ENERGY,        &
                        ior(FAST, DAY_MON), FLD_C_DIAG0 + CO_LEAF_TEMP)
      call add_variable(reg, 'mean_wood_temp_co', 'per-cohort wood temperature', 'K', DIM_COHORT, AGG_TMEAN, GRP_ENERGY,        &
                        DAY_MON, FLD_C_DIAG0 + CO_WOOD_TEMP)
      call add_variable(reg, 'mean_cas_temp_pa', 'patch canopy-air-space temperature', 'K', DIM_PATCH, AGG_TMEAN, GRP_ENERGY,   &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_CAS_TEMP)
      call add_variable(reg, 'mean_soil_temp_top_pa', 'patch soil-top temperature', 'K', DIM_PATCH, AGG_TMEAN, GRP_ENERGY,      &
                        ior(FAST, MON), FLD_P_DIAG0 + PA_SOIL_TEMP_TOP)
      call add_variable(reg, 'mean_skin_temp_pa', 'patch skin temperature (longwave up)', 'K', DIM_PATCH, AGG_TMEAN, GRP_ENERGY, &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_SKIN_TEMP)
      call add_variable(reg, 'mean_soil_temp_layer_pa', 'soil temperature by layer and patch', 'K', DIM_SOIL_PATCH, AGG_TMEAN,  &
                        GRP_ENERGY, FREQ_NONE, FLD_L_SOIL_TEMP, fast_source=SRC_F_SOIL_TEMP)
      call add_variable(reg, 'mean_soil_fliq_layer_pa', 'soil liquid fraction by layer and patch', '-', DIM_SOIL_PATCH,         &
                        AGG_TMEAN, GRP_ENERGY, FREQ_NONE, FLD_L_SOIL_FLIQ,                    &
                        fast_source=SRC_F_SOIL_LAYER)
   end subroutine register_energy

   !--- Biogeochemistry: the CENTURY matrix pools and Rh, all 0 when soil_carbon_on is off.
   subroutine register_biogeochem(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_soilc_fast_grnd_py', 'polygon fast/metabolic litter carbon, above-ground', 'kgC/m2',         &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_FAST_GRND)
      call add_variable(reg, 'mean_soilc_fast_soil_py', 'polygon fast/metabolic litter carbon, below-ground', 'kgC/m2',         &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_FAST_SOIL)
      call add_variable(reg, 'mean_soilc_struct_grnd_py', 'polygon structural litter + CWD carbon, above-ground', 'kgC/m2',        &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_STRUCT_GRND)
      call add_variable(reg, 'mean_soilc_struct_soil_py', 'polygon structural litter + CWD carbon, below-ground', 'kgC/m2',        &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_STRUCT_SOIL)
      call add_variable(reg, 'mean_soilc_microbial_py', 'polygon microbial SOM carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN,            &
                        GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_MICROBIAL)
      call add_variable(reg, 'mean_soilc_slow_py', 'polygon slow/humified SOM carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN,             &
                        GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_SLOW)
      call add_variable(reg, 'mean_soilc_passive_py', 'polygon passive SOM carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN,             &
                        GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_PASSIVE)
      call add_variable(reg, 'mean_soilc_total_py', 'polygon total soil organic carbon', 'kgC/m2', DIM_SCALAR, AGG_MEAN,           &
                        GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_TOTAL)
      call add_variable(reg, 'mean_rh_py', 'polygon heterotrophic respiration', 'umol/m2/s', DIM_SCALAR, AGG_RATE, GRP_BIOGEOCHEM, &
                        DAY_MON_YR, FLD_P_RH, scale=kgc_2_umol, fast_source=SRC_F_PA0 + PA_RH)
      call add_variable(reg, 'mean_soilc_total_pa', 'patch total soil organic carbon', 'kgC/m2', DIM_PATCH, AGG_MEAN,           &
                        GRP_BIOGEOCHEM, MON, FLD_P_SOILC_TOTAL)
      call add_variable(reg, 'mean_rh_pa', 'patch heterotrophic respiration', 'umol/m2/s', DIM_PATCH, AGG_RATE, GRP_BIOGEOCHEM, &
                        MON, FLD_P_RH, scale=kgc_2_umol, fast_source=SRC_F_PA0 + PA_RH)
   end subroutine register_biogeochem

   !---~---
   !   Numerics: the integrators' work, each count's mean per fast step over the period (the
   ! daily tier shows the one day that went wrong).
   !---~---
   subroutine register_numerics(reg, dt_fast)
      type(output_registry_t), intent(inout) :: reg
      real(wp),                intent(in)    :: dt_fast   ! [s] a count per second times this: per fast step
      call add_variable(reg, 'mean_work_integ_steps_py', 'accepted integrator sub-steps per fast step, over the patches', '1',  &
                        DIM_SCALAR, AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_STEPS, scale=dt_fast)
      call add_variable(reg, 'mean_work_integ_rej_py', 'rejected integrator steps per fast step, over the patches', '1',        &
                        DIM_SCALAR, AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_REJ, scale=dt_fast)
      call add_variable(reg, 'mean_work_water_nsub_py', 'soil-plant water network sub-steps per fast step', '1', DIM_SCALAR,    &
                        AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_WATER_NSUB, scale=dt_fast)
      call add_variable(reg, 'mean_work_water_newton_py', 'soil-plant water network Newton iterations per fast step', '1',      &
                        DIM_SCALAR, AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_WATER_NEWTON, scale=dt_fast)
      call add_variable(reg, 'mean_work_water_reject_py',                                                                       &
                        'soil-plant water network sub-steps rejected and retried shorter, per fast step', '1', DIM_SCALAR,      &
                        AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_WATER_REJECT, scale=dt_fast)
      call add_variable(reg, 'mean_work_water_nfail_py',                                                                        &
                        'soil-plant water network sub-steps taken although unconverged, per fast step', '1', DIM_SCALAR,        &
                        AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_WATER_NFAIL, scale=dt_fast)
      call add_variable(reg, 'mean_work_rk45_rescue_py', 'RK45 steps rescued to the implicit path, per fast step', '1',         &
                        DIM_SCALAR, AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_RK45_RESCUE, scale=dt_fast)
      call add_variable(reg, 'mean_work_clamp_stage_py', 'stage-input clamp activations per fast step', '1', DIM_SCALAR,        &
                        AGG_RATE, GRP_NUMERICS, DAY_MON_YR, SRC_PY_WORK_CLAMP_STAGE, scale=dt_fast)
   end subroutine register_numerics


   !---~---
   !   Ecophysiology (GRP_ECOPHYS): the per-cohort leaf gas exchange and hydraulics, from the
   ! fast loop's diagnostic block; off by default, the highest-volume group
   ! (docs/science/diagnostics.md §5). With them, their leaf-area-weighted canopy means.
   !---~---
   subroutine register_ecophys(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_anet_co', 'net leaf assimilation', 'umol/m2 leaf/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,     &
                        FREQ_NONE, FLD_C_DIAG0 + CO_ANET, flux=.true.)
      call add_variable(reg, 'mean_agross_co', 'gross leaf assimilation', 'umol/m2 leaf/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, &
                        FREQ_NONE, FLD_C_DIAG0 + CO_AGROSS, flux=.true.)
      call add_variable(reg, 'mean_gsw_co', 'stomatal conductance to water vapour', 'mol/m2 leaf/s', DIM_COHORT, AGG_TMEAN,     &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_GSW)
      call add_variable(reg, 'mean_gbw_co', 'leaf boundary-layer conductance', 'm/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,        &
                        FREQ_NONE, FLD_C_DIAG0 + CO_GBW)
      call add_variable(reg, 'mean_ci_co', 'intercellular CO2 mole fraction', 'umol/mol', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,   &
                        FREQ_NONE, FLD_C_DIAG0 + CO_CI)
      call add_variable(reg, 'mean_cs_co', 'leaf-surface CO2 mole fraction', 'umol/mol', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,    &
                        FREQ_NONE, FLD_C_DIAG0 + CO_CS)
      call add_variable(reg, 'mean_rd_co', 'leaf dark respiration', 'umol/m2 leaf/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,        &
                        FREQ_NONE, FLD_C_DIAG0 + CO_RD, flux=.true.)
      call add_variable(reg, 'mean_transp_co', 'leaf transpiration (the water the leaves lose)', 'mol/m2 leaf/s', DIM_COHORT,   &
                        AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_TRANSP, flux=.true.)
      call add_variable(reg, 'mean_beta_stomata_co', 'stomatal water-stress factor', '-', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,   &
                        FREQ_NONE, FLD_C_DIAG0 + CO_BETA_STOM)
      call add_variable(reg, 'mean_beta_nonstomata_co', 'non-stomatal (capacity) water-stress factor', '-', DIM_COHORT,         &
                        AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_BETA_NONSTOM)
      call add_variable(reg, 'mean_leaf_vpd_co', 'leaf-to-air vapour-pressure deficit', 'Pa', DIM_COHORT, AGG_TMEAN,            &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_LEAF_VPD)
      call add_variable(reg, 'mean_psi_leaf_co', 'leaf water potential', 'MPa', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,             &
                        FREQ_NONE, FLD_C_DIAG0 + CO_PSI_LEAF)
      call add_variable(reg, 'mean_psi_wood_co', 'wood water potential', 'MPa', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE,  &
                        FLD_C_DIAG0 + CO_PSI_WOOD)
      call add_variable(reg, 'mean_plc_co', 'plant loss of hydraulic conductance', '-', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,     &
                        FREQ_NONE, FLD_C_DIAG0 + CO_PLC)
      call add_variable(reg, 'mean_sapflow_co', 'wood-to-leaf sapflow', 'kg/plant/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,       &
                        FREQ_NONE, FLD_C_DIAG0 + CO_SAPFLOW, flux=.true.)
      call add_variable(reg, 'mean_root_uptake_co', 'soil-to-root water uptake', 'kg/plant/s', DIM_COHORT, AGG_TMEAN,           &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_ROOT_UPTAKE, flux=.true.)
      call add_variable(reg, 'mean_leaf_water_co', 'internal leaf water', 'kg/plant', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,        &
                        FREQ_NONE, FLD_C_DIAG0 + CO_LEAF_WATER)
      call add_variable(reg, 'mean_wood_water_co', 'internal wood water', 'kg/plant', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,        &
                        FREQ_NONE, FLD_C_DIAG0 + CO_WOOD_WATER)
      call add_variable(reg, 'mean_abs_par_co', 'absorbed photosynthetically active radiation', 'W/m2', DIM_COHORT, AGG_TMEAN,  &
                        GRP_RADIATION, FREQ_NONE, FLD_C_DIAG0 + CO_ABS_PAR, flux=.true.)
      call add_variable(reg, 'mean_abs_sw_co', 'absorbed shortwave radiation', 'W/m2', DIM_COHORT, AGG_TMEAN, GRP_RADIATION,    &
                        FREQ_NONE, FLD_C_DIAG0 + CO_ABS_SW, flux=.true.)
      call add_variable(reg, 'mean_abs_lw_co', 'net longwave radiation', 'W/m2', DIM_COHORT, AGG_TMEAN, GRP_RADIATION, FREQ_NONE, &
                        FLD_C_DIAG0 + CO_ABS_LW, flux=.true.)
      call add_variable(reg, 'mean_wind_co', 'in-canopy wind speed', 'm/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE,      &
                        FLD_C_DIAG0 + CO_WIND)
      !--- Canopy twins: leaf-area-weighted means, the standard flux-tower comparison.
      call add_variable(reg, 'mean_anet_py', 'leaf-area-weighted canopy net assimilation', 'umol/m2 leaf/s', DIM_SCALAR,        &
                        AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_ANET, weight=W_LEAF_AREA, mean=.true., flux=.true.)
      call add_variable(reg, 'mean_gsw_py', 'leaf-area-weighted canopy stomatal conductance', 'mol/m2 leaf/s', DIM_SCALAR,      &
                        AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_GSW, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_psi_leaf_py', 'leaf-area-weighted canopy leaf water potential', 'MPa', DIM_SCALAR, AGG_TMEAN, &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_PSI_LEAF, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_plc_py', 'leaf-area-weighted canopy loss of conductance', '-', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_PLC, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_beta_stomata_py', 'leaf-area-weighted stomatal water-stress factor', '-', DIM_SCALAR,        &
                        AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_BETA_STOM, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_root_uptake_py', 'polygon root water uptake', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_WATER,      &
                        DAY_MON_YR, FLD_C_DIAG0 + CO_ROOT_UPTAKE, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_transp_leaf_py', 'leaf-area-weighted canopy transpiration', 'mol/m2 leaf/s', DIM_SCALAR,     &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_C_DIAG0 + CO_TRANSP, weight=W_LEAF_AREA, mean=.true., flux=.true.)
   end subroutine register_ecophys

   !---~---
   !   Slow-loop carbon allocation: where compute_carbon_allocation sends each cohort's carbon
   ! every step, so the carbon budget closes from an output file: npp_site reconciles
   ! against its five destinations plus growth respiration.
   !---~---
   subroutine register_allocation(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_npp_leaf_py', 'NPP allocated to leaf', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_CARBON,       &
                        MON_YR, FLD_C_DIAG0 + CO_NPP_LEAF, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_npp_fineroot_py', 'NPP allocated to fine root', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN,          &
                        GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_FINEROOT, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_npp_wood_py', 'NPP allocated to wood', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_CARBON,       &
                        MON_YR, FLD_C_DIAG0 + CO_NPP_WOOD, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_npp_storage_py', 'NPP to non-structural storage', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_STORAGE, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_npp_repro_py', 'NPP allocated to reproduction', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN,          &
                        GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_REPRO, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_growth_resp_py', 'growth respiration', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_CARBON,       &
                        MON_YR, FLD_C_DIAG0 + CO_GROWTH_RESP, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_root_exudate_py', 'root exudate: carbon the growth sink could not use', 'kgC/m2/yr',         &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_EXUDATE, weight=W_NPLANT, flux=.true.)
      !---~---
      !   The polygon's demographic rates, the terms that make the AGB trajectory: d(agb)/dt
      ! against growth - mortality from the file (docs/science/diagnostics.md eq. (4)).
      !---~---
      call add_variable(reg, 'mean_agb_growth_py', 'polygon AGB growth rate', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE,   &
                        MON_YR, FLD_C_DIAG0 + CO_DAGB_DT, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_agb_mort_py', 'polygon AGB mortality loss rate', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN,            &
                        GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, weight=W_AGB, flux=.true.)
      call add_variable(reg, 'mean_mort_rate_py', 'nplant-weighted mortality rate', '1/yr', DIM_SCALAR, AGG_TMEAN,              &
                        GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, weight=W_NPLANT, mean=.true., flux=.true.)
      !---~---
      !   Litterfall, the vegetation -> soil carbon seam, by the tissue the carbon came from,
      ! read from the accumulator the biogeochemistry consumes (docs/science/diagnostics.md §6,
      ! The carbon seam). litter_struct_site is a destination, so it overlaps the others.
      !---~---
      call add_variable(reg, 'mean_litter_leaf_py', 'leaf litterfall carbon (after resorption)', 'kgC/m2/yr', DIM_SCALAR,       &
                        AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_LEAF, flux=.true.)
      call add_variable(reg, 'mean_litter_fineroot_py', 'fine-root litter carbon (after resorption)', 'kgC/m2/yr', DIM_SCALAR,  &
                        AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_FINEROOT, flux=.true.)
      call add_variable(reg, 'mean_litter_wood_py', 'wood and storage litter carbon of dead plants', 'kgC/m2/yr', DIM_SCALAR,   &
                        AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_WOOD, flux=.true.)
      call add_variable(reg, 'mean_litter_repro_py', 'fruit and seed litter carbon', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_REPRO, flux=.true.)
      call add_variable(reg, 'mean_litter_struct_py', 'structural litter + CWD carbon', 'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN,     &
                        GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_STRUCT, flux=.true.)
      call add_variable(reg, 'mean_nplant_recruit_py', 'recruitment density flux', 'plant/m2/yr', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_RECRUIT_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_disturb_area_py', 'patch area fraction disturbed', '1/yr', DIM_SCALAR, AGG_TMEAN,            &
                        GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_DISTURB_AREA, flux=.true.)
      !---~---
      !   Mortality carbon by pathway: three exclusive parts of the whole-individual mortality
      ! carbon, every pool included (docs/science/diagnostics.md §6, The carbon seam). Written
      ! whatever soil_carbon_on says: this is demography.
      !---~---
      call add_variable(reg, 'mean_mort_carbon_background_py', 'mortality carbon, background hazard', 'kgC/m2/yr', DIM_SCALAR,  &
                        AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_MORT_C_BACKGROUND, flux=.true.)
      call add_variable(reg, 'mean_mort_carbon_cull_py', 'mortality carbon, cohorts culled below the tracking floor',           &
                        'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_MORT_C_CULL, flux=.true.)
      call add_variable(reg, 'mean_mort_carbon_disturb_py', 'mortality carbon, canopy killed by disturbance', 'kgC/m2/yr',      &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_MORT_C_DISTURB, flux=.true.)
      call add_variable(reg, 'mean_mort_carbon_background_pa', 'patch mortality carbon, background hazard', 'kgC/m2/yr',        &
                        DIM_PATCH, AGG_TMEAN, GRP_STRUCTURE, MON, FLD_P_DIAG0 + PA_MORT_C_BACKGROUND, flux=.true.)
      call add_variable(reg, 'mean_litter_leaf_pa', 'patch leaf litterfall carbon', 'kgC/m2/yr', DIM_PATCH, AGG_TMEAN,          &
                        GRP_BIOGEOCHEM, MON, FLD_P_DIAG0 + PA_LITTER_LEAF, flux=.true.)
      call add_variable(reg, 'mean_agb_growth_pft_py', 'AGB growth rate by PFT', 'kgC/m2/yr', DIM_PFT, AGG_TMEAN, GRP_STRUCTURE, &
                        MON_YR, FLD_C_DIAG0 + CO_DAGB_DT, weight=W_NPLANT, flux=.true.)
      call add_variable(reg, 'mean_agb_mort_pft_py', 'AGB mortality loss rate by PFT', 'kgC/m2/yr', DIM_PFT, AGG_TMEAN,         &
                        GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, weight=W_AGB, flux=.true.)
   end subroutine register_allocation

   !---~---
   !   The per-patch fast-loop fluxes and states, and their polygon (_site) twins: the energy,
   ! water and carbon exchange a run is judged on, at every coarse tier.
   !---~---
   subroutine register_patch_fluxes(reg)
      type(output_registry_t), intent(inout) :: reg
      !--- energy
      call add_variable(reg, 'mean_le_py', 'latent heat flux', 'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, ior(FAST, DAY_MON_YR), &
                        FLD_P_DIAG0 + PA_LE, flux=.true.)
      call add_variable(reg, 'mean_h_py', 'sensible heat flux', 'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY,                      &
                        ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_H, flux=.true.)
      call add_variable(reg, 'mean_rnet_py', 'net all-wave radiation', 'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY,               &
                        ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_RNET, flux=.true.)
      !---~---
      !   Top-of-canopy radiative fluxes, for a radiometer or satellite comparison: fluxes, not
      ! a mean albedo, which is a mean of ratios (docs/science/canopy_radiation_transfer.md
      ! eq. (16)).
      !---~---
      call add_variable(reg, 'mean_sw_in_vis_py', 'incident VIS shortwave at canopy top', 'W/m2', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_IN_VIS, flux=.true.)
      call add_variable(reg, 'mean_sw_in_nir_py', 'incident NIR shortwave at canopy top', 'W/m2', DIM_SCALAR, AGG_TMEAN,        &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_IN_NIR, flux=.true.)
      call add_variable(reg, 'mean_sw_up_vis_py', 'upwelling VIS shortwave at canopy top', 'W/m2', DIM_SCALAR, AGG_TMEAN,       &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_UP_VIS, flux=.true.)
      call add_variable(reg, 'mean_sw_up_nir_py', 'upwelling NIR shortwave at canopy top', 'W/m2', DIM_SCALAR, AGG_TMEAN,       &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_UP_NIR, flux=.true.)
      call add_variable(reg, 'mean_lw_up_py', 'upwelling longwave at canopy top (emission included)', 'W/m2', DIM_SCALAR,       &
                        AGG_TMEAN, GRP_ENERGY, ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_LW_UP, flux=.true.)
      call add_variable(reg, 'mean_sw_in_py', 'incident shortwave at canopy top', 'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_FORCING,   &
                        ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_SW_IN, flux=.true.)
      call add_variable(reg, 'mean_sw_ground_py', 'shortwave absorbed at the ground', 'W/m2', DIM_SCALAR, AGG_TMEAN,            &
                        GRP_RADIATION, FREQ_NONE, FLD_P_DIAG0 + PA_SW_GROUND, flux=.true.)
      call add_variable(reg, 'mean_lw_ground_py', 'net longwave at the ground', 'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_RADIATION,   &
                        FREQ_NONE, FLD_P_DIAG0 + PA_LW_GROUND, flux=.true.)
      call add_variable(reg, 'mean_ustar_py', 'friction velocity', 'm/s', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY,                    &
                        ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_USTAR)
      call add_variable(reg, 'mean_ggnet_py', 'ground conductance', 'm/s', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON,              &
                        FLD_P_DIAG0 + PA_GGNET)
      call add_variable(reg, 'mean_rough_py', 'aerodynamic roughness length', 'm', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR,   &
                        FLD_P_DIAG0 + PA_ROUGH)
      call add_variable(reg, 'mean_displace_py', 'zero-plane displacement height', 'm', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY,      &
                        MON_YR, FLD_P_DIAG0 + PA_DISPLACE)
      !---~---
      !   Per patch too: with its canopy-air depth these are what move the forcing to the patch's
      ! canopy-air top (docs/science/forcing.md §8), so each patch's forcing can be rebuilt.
      !---~---
      call add_variable(reg, 'mean_rough_pa', 'patch aerodynamic roughness length', 'm', DIM_PATCH, AGG_TMEAN, GRP_ENERGY,      &
                        DAY_MON, FLD_P_DIAG0 + PA_ROUGH)
      call add_variable(reg, 'mean_displace_pa', 'patch zero-plane displacement height', 'm', DIM_PATCH, AGG_TMEAN, GRP_ENERGY, &
                        DAY_MON, FLD_P_DIAG0 + PA_DISPLACE)
      !--- carbon
      call add_variable(reg, 'mean_nee_py', 'net ecosystem exchange', 'umol/m2/s', DIM_SCALAR, AGG_TMEAN,                        &
                        GRP_CARBON, ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_NEE, flux=.true.)
      !--- water
      call add_variable(reg, 'mean_precip_py', 'total precipitation', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_FORCING,             &
                        ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_PRECIP, flux=.true.)
      !---~---
      !   ET's sources, as the canopy air's vapour ledger booked them: mean_et_py is the first three
      ! less the dew, less the canopy air's own vapour gain (docs/science/diagnostics.md §4).
      !---~---
      call add_variable(reg, 'mean_transp_py', 'canopy transpiration', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_WATER,             &
                        DAY_MON_YR, FLD_P_DIAG0 + PA_TRANSP, flux=.true.)
      call add_variable(reg, 'mean_evap_canopy_py', 'evaporation of the canopy''s intercepted water (dew onto it < 0)',          &
                        'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_EVAP_CANOPY, flux=.true.)
      call add_variable(reg, 'mean_evap_ground_py', 'ground evaporation: bare soil and snow pack', 'kg/m2/s', DIM_SCALAR,       &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_EVAP_GROUND, flux=.true.)
      call add_variable(reg, 'mean_dew_py', 'condensation from the canopy air (dew, fog), to the soil', 'kg/m2/s', DIM_SCALAR,  &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_DEW, flux=.true.)
      !--- budget health (GRP_NUMERICS): the numbers that say whether anything above is real.
      call add_variable(reg, 'mean_resid_energy_py', 'signed whole-column energy-budget residual (+ = appearing)', 'W/m2',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_NUMERICS, DAY_MON_YR, FLD_P_DIAG0 + PA_RESID_ENERGY, flux=.true.)
      call add_variable(reg, 'mean_resid_water_py', 'signed whole-column water-budget residual (+ = appearing)', 'kg/m2/s',      &
                        DIM_SCALAR, AGG_TMEAN, GRP_NUMERICS, DAY_MON_YR, FLD_P_DIAG0 + PA_RESID_WATER, flux=.true.)
      !--- per-patch twins of the four that vary most between a gap and a closed canopy.
      call add_variable(reg, 'mean_le_pa', 'patch latent heat flux', 'W/m2', DIM_PATCH, AGG_TMEAN, GRP_ENERGY,                  &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_LE, flux=.true.)
      call add_variable(reg, 'mean_h_pa', 'patch sensible heat flux', 'W/m2', DIM_PATCH, AGG_TMEAN, GRP_ENERGY,                 &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_H, flux=.true.)
      call add_variable(reg, 'mean_rnet_pa', 'patch net all-wave radiation', 'W/m2', DIM_PATCH, AGG_TMEAN, GRP_ENERGY,          &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_RNET, flux=.true.)
      call add_variable(reg, 'mean_nee_pa', 'patch net ecosystem exchange', 'umol/m2/s', DIM_PATCH, AGG_TMEAN, GRP_CARBON,      &
                        ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_NEE, flux=.true.)
   end subroutine register_patch_fluxes

   !---~---
   !   The budgets' terms (docs/science/diagnostics.md §4, The budgets): every flux the
   ! whole-column water and energy ledgers book, so each budget closes from the files, its
   ! residual (mean_resid_*_py) included; the canopy's carbon split; and the stocks.
   !---~---
   subroutine register_budgets(reg)
      type(output_registry_t), intent(inout) :: reg
      !--- water: the fluxes
      call add_variable(reg, 'mean_intercept_py', 'precipitation the canopy films captured, net', 'kg/m2/s', DIM_SCALAR,     &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_INTERCEPT, flux=.true.)
      call add_variable(reg, 'mean_throughfall_py', 'precipitation reaching the ground', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN,   &
                        GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_THROUGHFALL, flux=.true.)
      call add_variable(reg, 'mean_infiltration_py', 'infiltration: pond to soil layer 1', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN, &
                        GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_INFILTRATION, flux=.true.)
      call add_variable(reg, 'mean_snowmelt_py', 'meltwater the snow pack drains to the pond', 'kg/m2/s', DIM_SCALAR,        &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_SNOWMELT, flux=.true.)
      call add_variable(reg, 'mean_sublimation_py', 'sublimation from the snow pack (deposition < 0)', 'kg/m2/s',            &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_SUBLIMATION, flux=.true.)
      call add_variable(reg, 'mean_shed_water_py', 'tissue water the slow step shed to the pond', 'kg/m2/s', DIM_SCALAR,     &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_SHED_WATER, flux=.true.)
      call add_variable(reg, 'mean_film_excess_py', 'canopy film water above capacity, removed', 'kg/m2/s', DIM_SCALAR,      &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_FILM_EXCESS, flux=.true.)
      call add_variable(reg, 'mean_water_store_gain_py', 'gain of every water store', 'kg/m2/s', DIM_SCALAR, AGG_TMEAN,     &
                        GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_WATER_STORE_GAIN, flux=.true.)
      !--- water: the stores
      call add_variable(reg, 'mean_canopy_water_py', 'water on the leaves'' and wood''s surfaces', 'kg/m2', DIM_SCALAR,      &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_CANOPY_WATER)
      call add_variable(reg, 'mean_soil_water_total_py', 'the soil column''s water, liquid and ice', 'kg/m2', DIM_SCALAR,    &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_SOIL_WATER_TOTAL)
      call add_variable(reg, 'mean_plant_water_py', 'the plants'' tissue water', 'kg/m2', DIM_SCALAR, AGG_TMEAN, GRP_WATER,  &
                        DAY_MON_YR, FLD_P_DIAG0 + PA_PLANT_WATER)
      !--- the wet leaf: each cohort's film-covered fraction, and the canopy's leaf-area-weighted
      call add_variable(reg, 'mean_wet_fraction_py', 'leaf-area-weighted wet fraction of the leaves', '-', DIM_SCALAR,       &
                        AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_C_DIAG0 + CO_WET_FRACTION, weight=W_LEAF_AREA, mean=.true.)
      call add_variable(reg, 'mean_wet_fraction_co', 'per-cohort wet fraction of the leaves', '-', DIM_COHORT, AGG_TMEAN,    &
                        GRP_WATER, FREQ_NONE, FLD_C_DIAG0 + CO_WET_FRACTION)
      !--- energy: the sensible heat's sources, and the terms beyond Rn, H and LE
      call add_variable(reg, 'mean_h_canopy_py', 'sensible heat, leaves and wood to the canopy air', 'W/m2', DIM_SCALAR,     &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_H_CANOPY, flux=.true.)
      call add_variable(reg, 'mean_h_ground_py', 'sensible heat, ground (soil, snow) to the canopy air', 'W/m2', DIM_SCALAR, &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_H_GROUND, flux=.true.)
      call add_variable(reg, 'mean_heat_bottom_py', 'heat into the soil column through its bottom', 'W/m2', DIM_SCALAR,      &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_HEAT_BOTTOM, flux=.true.)
      call add_variable(reg, 'mean_enth_advect_py', 'enthalpy the water carries in, beyond H and LE', 'W/m2', DIM_SCALAR,    &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_ENTH_ADVECT, flux=.true.)
      call add_variable(reg, 'mean_heat_gain_soil_py', 'the soil column''s energy gain', 'W/m2', DIM_SCALAR, AGG_TMEAN,      &
                        GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_HEAT_GAIN_SOIL, flux=.true.)
      call add_variable(reg, 'mean_heat_gain_canopy_py', 'the canopy air''s and tissues'' energy gain', 'W/m2', DIM_SCALAR,  &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_HEAT_GAIN_CANOPY, flux=.true.)
      call add_variable(reg, 'mean_heat_gain_surface_py', 'the snow pack''s and pond''s energy gain', 'W/m2', DIM_SCALAR,    &
                        AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_HEAT_GAIN_SURFACE, flux=.true.)
      call add_variable(reg, 'mean_apar_py', 'PAR (VIS) the leaves absorbed', 'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_RADIATION, &
                        DAY_MON_YR, FLD_P_DIAG0 + PA_APAR, flux=.true.)
      !--- carbon: autotrophic respiration (Reco - Rh) and the ecosystem's carbon
      call add_variable(reg, 'mean_ra_py', 'autotrophic respiration (ecosystem less heterotrophic)', 'umol/m2/s',           &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, DAY_MON_YR, FLD_P_DIAG0 + PA_RA, flux=.true.)
      call add_variable(reg, 'mean_eco_carbon_py', 'ecosystem carbon: live vegetation + soil organic carbon', 'kgC/m2',      &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, SRC_PY_ECO_CARBON)
   end subroutine register_budgets

   !---~---
   !   FAST-tier (sub-daily) diagnostics. They read the staged sub-step samples, not live
   ! state, so they carry the FAST bit only (docs/science/diagnostics.md §4, How the FAST
   ! tier is sampled).
   !---~---
   subroutine register_fast(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'mean_sw_up_py', 'shortwave leaving the canopy top (VIS + NIR)', 'W/m2', DIM_SCALAR, AGG_TMEAN,    &
                        GRP_ENERGY, FAST, FLD_P_DIAG0 + PA_SW_UP, flux=.true.)
      call add_variable(reg, 'mean_reco_py', 'ecosystem respiration (autotrophic + heterotrophic)', 'umol/m2/s', DIM_SCALAR,    &
                        AGG_TMEAN, GRP_CARBON, ior(FAST, DAY_MON_YR), FLD_P_DIAG0 + PA_RECO, flux=.true.)
      !---~---
      !   The leaf solve's inputs and output per cohort, sub-daily, as the solve saw them (rows of
      ! the per-cohort diagnostic block, SRC_F_CO0 + CO_*): enough to repeat the solve outside
      ! the model, which the fast calibration's canopy of leaf kernels does (calibrate_fast).
      !---~---
      call add_variable(reg, 'mean_leaf_par_co', 'leaf solve: incident-equivalent PAR per leaf area', 'umol/m2/s', DIM_COHORT,  &
                        AGG_TMEAN, GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_LEAF_PAR, flux=.true.)
      call add_variable(reg, 'mean_ca_co', 'leaf solve: canopy-air CO2', 'umol/mol', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,        &
                        FREQ_NONE, FLD_C_DIAG0 + CO_CA)
      call add_variable(reg, 'mean_pressure_co', 'leaf solve: air pressure', 'Pa', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS,          &
                        FREQ_NONE, FLD_C_DIAG0 + CO_PRESSURE)
      call add_variable(reg, 'mean_psi_predawn_co', 'leaf solve: predawn leaf water potential', 'MPa', DIM_COHORT, AGG_TMEAN,   &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_PSI_PREDAWN)
      call add_variable(reg, 'mean_gb_mol_co', 'leaf solve: boundary-layer conductance', 'mol/m2/s', DIM_COHORT, AGG_TMEAN,     &
                        GRP_ECOPHYS, FREQ_NONE, FLD_C_DIAG0 + CO_GB_MOL)
      !---~---
      !   Per-patch twins: the same samples before the sum over patches, named and switched as
      ! docs/science/diagnostics.md §4 (The FAST tier by patch) says; sub-daily by default, and
      ! daily or monthly on request.
      !---~---
      call add_variable(reg, 'mean_sw_up_pa', 'patch shortwave leaving the canopy top (VIS + NIR)', 'W/m2', DIM_PATCH,          &
                        AGG_TMEAN, GRP_ENERGY, FAST, FLD_P_DIAG0 + PA_SW_UP, flux=.true.)
      call add_variable(reg, 'mean_lw_up_pa', 'patch longwave leaving the canopy top (emission included)', 'W/m2', DIM_PATCH,   &
                        AGG_TMEAN, GRP_ENERGY, FAST, FLD_P_DIAG0 + PA_LW_UP, flux=.true.)
      call add_variable(reg, 'mean_ustar_pa', 'patch friction velocity', 'm/s', DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST,    &
                        FLD_P_DIAG0 + PA_USTAR)
      call add_variable(reg, 'mean_reco_pa', 'patch ecosystem respiration (autotrophic + heterotrophic)', 'umol/m2/s',          &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, ior(FAST, DAY_MON), FLD_P_DIAG0 + PA_RECO, flux=.true.)
      call add_variable(reg, 'mean_cas_co2_pa', 'patch canopy-air CO2 mixing ratio', 'umol/mol', DIM_PATCH, AGG_TMEAN,          &
                        GRP_CARBON, FAST, FLD_P_DIAG0 + PA_CAS_CO2)
   end subroutine register_fast

   !---~---
   !   The forcing the run used (docs/science/diagnostics.md §5), after the reader's partitions
   ! and corrections: the coarse tiers' dt-weighted means, kept once per polygon (poly%diag,
   ! PY_*; the totals are sw_in_site and precip_site), and the FAST tier's sub-step samples.
   !---~---
   subroutine register_forcing(reg)
      type(output_registry_t), intent(inout) :: reg
      !---~---
      !   Each site's forcing: the polygon's after the site's terrain (docs/science/forcing.md
      ! §8), and the growth temperature its leaves acclimate to.
      !---~---
      call add_variable(reg, 'mean_air_temp_si', 'each site''s air temperature at the forcing''s height', 'K', DIM_SITE,        &
                        AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_TAIR)
      call add_variable(reg, 'mean_qair_si', 'each site''s specific humidity', 'kg/kg', DIM_SITE, AGG_TMEAN, GRP_FORCING,       &
                        DAY_MON_YR, FLD_SI_DIAG0 + SI_QAIR)
      call add_variable(reg, 'mean_psurf_si', 'each site''s surface pressure', 'Pa', DIM_SITE, AGG_TMEAN, GRP_FORCING,          &
                        DAY_MON_YR, FLD_SI_DIAG0 + SI_PSURF)
      call add_variable(reg, 'mean_rho_air_si', 'each site''s air density', 'kg/m3', DIM_SITE, AGG_TMEAN, GRP_FORCING,          &
                        DAY_MON_YR, FLD_SI_DIAG0 + SI_RHO_AIR)
      call add_variable(reg, 'mean_lwdown_si', 'each site''s downward longwave', 'W/m2', DIM_SITE, AGG_TMEAN, GRP_FORCING,      &
                        DAY_MON_YR, FLD_SI_DIAG0 + SI_LWDOWN, flux=.true.)
      call add_variable(reg, 'mean_sw_in_si', 'each site''s incident shortwave, its beam on its slope', 'W/m2', DIM_SITE,       &
                        AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_SW_IN, flux=.true.)
      call add_variable(reg, 'mean_par_beam_si', 'each site''s direct-beam PAR on its slope', 'W/m2', DIM_SITE, AGG_TMEAN,      &
                        GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_PAR_BEAM, flux=.true.)
      call add_variable(reg, 'mean_nir_beam_si', 'each site''s direct-beam NIR on its slope', 'W/m2', DIM_SITE, AGG_TMEAN,      &
                        GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_NIR_BEAM, flux=.true.)
      call add_variable(reg, 'mean_growth_temp_si', 'each site''s running-mean growth temperature', 'K', DIM_SITE, AGG_TMEAN,   &
                        GRP_FORCING, DAY_MON_YR, FLD_SI_GROWTH_TEMP)
      call add_variable(reg, 'mean_air_temp_py', 'air temperature at the forcing''s own height (forcing)', 'K', DIM_SCALAR,     &
                        AGG_TMEAN, GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_TAIR)
      call add_variable(reg, 'mean_qair_py', 'specific humidity at the forcing''s own height (forcing)', 'kg/kg', DIM_SCALAR,   &
                        AGG_TMEAN, GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_QAIR)
      call add_variable(reg, 'mean_psurf_py', 'surface pressure (forcing)', 'Pa', DIM_SCALAR, AGG_TMEAN, GRP_FORCING,           &
                        ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_PSURF)
      call add_variable(reg, 'mean_wind_py', 'wind speed at the forcing''s own height (forcing)', 'm/s', DIM_SCALAR, AGG_TMEAN, &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_WIND)
      call add_variable(reg, 'mean_lwdown_py', 'downward longwave at canopy top (forcing)', 'W/m2', DIM_SCALAR, AGG_TMEAN,      &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_LWDOWN, flux=.true.)
      call add_variable(reg, 'mean_par_beam_py', 'direct-beam PAR at canopy top (forcing)', 'W/m2', DIM_SCALAR, AGG_TMEAN,      &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_PAR_BEAM, flux=.true.)
      call add_variable(reg, 'mean_par_diffuse_py', 'diffuse PAR at canopy top (forcing)', 'W/m2', DIM_SCALAR, AGG_TMEAN,       &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_PAR_DIFFUSE, flux=.true.)
      call add_variable(reg, 'mean_nir_beam_py', 'direct-beam NIR at canopy top (forcing)', 'W/m2', DIM_SCALAR, AGG_TMEAN,      &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_NIR_BEAM, flux=.true.)
      call add_variable(reg, 'mean_nir_diffuse_py', 'diffuse NIR at canopy top (forcing)', 'W/m2', DIM_SCALAR, AGG_TMEAN,       &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_NIR_DIFFUSE, flux=.true.)
      call add_variable(reg, 'mean_snowfall_py', 'frozen precipitation (split at each patch''s canopy-air top)',                 &
                        'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_SNOWFALL,        &
                        flux=.true.)
      call add_variable(reg, 'mean_atm_co2_py', 'free-atmosphere CO2 (forcing)', 'umol/mol', DIM_SCALAR, AGG_TMEAN, GRP_FORCING, &
                        ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_CO2)
      call add_variable(reg, 'mean_cosz_py', 'cosine of the solar zenith angle (0 with the sun down)', '1', DIM_SCALAR,         &
                        AGG_TMEAN, GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_COSZ)
      call add_variable(reg, 'mean_rho_air_py', 'air density at the forcing''s own height', 'kg/m3', DIM_SCALAR, AGG_TMEAN,     &
                        GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_RHO_AIR)
      !---~---
      !   The patch's own forcing: the sample moved to its canopy-air top
      ! (docs/science/forcing.md §8). Only the wind and the air temperature differ between
      ! patches; the rest is the polygon's, above.
      !---~---
      call add_variable(reg, 'mean_wind_cas_top_pa', 'wind speed at the patch canopy-air top', 'm/s', DIM_PATCH, AGG_TMEAN,     &
                        GRP_FORCING, DAY_MON, FLD_P_DIAG0 + PA_WIND_CAS_TOP)
      call add_variable(reg, 'mean_air_temp_cas_top_pa', 'air temperature at the patch canopy-air top', 'K', DIM_PATCH,         &
                        AGG_TMEAN, GRP_FORCING, DAY_MON, FLD_P_DIAG0 + PA_TAIR_CAS_TOP)
      call add_variable(reg, 'mean_cas_depth_pa', 'patch canopy-air depth (the top its forcing is moved to)', 'm', DIM_PATCH,   &
                        AGG_TMEAN, GRP_FORCING, DAY_MON, FLD_P_DIAG0 + PA_Z_CAS_TOP)
      !--- The rain: the precipitation less mean_snowfall_py, as the patches split it.
      call add_variable(reg, 'mean_rainf_py', 'liquid precipitation (split at each patch''s canopy-air top)',                    &
                        'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_FORCING, ior(FAST, DAY_MON_YR), FLD_PY_DIAG0 + PY_RAINF,           &
                        flux=.true.)
   end subroutine register_forcing

   !---~---
   !   Append one descriptor. weight, mean and scale say how the field reduces (var_desc_t);
   ! their defaults (W_NONE, sum, 1.0) suit a raw slab, so the cohort and patch rows mostly
   ! omit them. A fast-loop diagnostic row is staged every sub-step, so its F source follows
   ! from its coarse one (staged_row); fast_source names any other.
   !---~---
   subroutine add_variable(reg, name, long_name, units, dim, agg, group, streams, source, xtype, weight, mean, scale, &
                           fast_source, flux, coord)
      type(output_registry_t), intent(inout) :: reg
      character(len=*),        intent(in)    :: name, long_name, units
      integer(ik),             intent(in)    :: dim, agg, group, streams
      integer(ik),             intent(in)    :: source        ! the D, M and Y tiers' field (SRC_NONE: none)
      integer(ik), optional,   intent(in)    :: fast_source   ! the F tier's staged value
      integer(ik), optional,   intent(in)    :: xtype, weight
      logical,     optional,   intent(in)    :: mean
      real(wp),    optional,   intent(in)    :: scale
      logical,     optional,   intent(in)    :: flux          ! a flux or a rate (AGG_RATE is one)
      logical,     optional,   intent(in)    :: coord         ! the axis's bookkeeping, written with it
      type(var_desc_t) :: v
      v%name            = name
      v%long_name       = long_name
      v%units           = units
      v%dim             = dim
      v%agg             = agg
      v%group           = group
      v%streams         = streams
      v%streams_default = streams
      v%enabled         = .true.
      v%source_id       = source
      v%source_fast     = staged_row(dim, source)
      if (present(fast_source)) v%source_fast = fast_source
      v%xtype           = XTYPE_DOUBLE
      if (present(xtype))  v%xtype  = xtype
      v%weight          = W_NONE
      if (present(weight)) v%weight = weight
      v%mean            = .false.
      if (present(mean))   v%mean   = mean
      v%scale           = 1.0_wp
      if (present(scale))  v%scale  = scale
      v%flux            = agg == AGG_RATE
      if (present(flux))   v%flux   = flux
      if (present(coord))  v%coord  = coord
      call append_variable(reg, v)
   end subroutine add_variable

   subroutine append_variable(reg, v)
      type(output_registry_t), intent(inout) :: reg
      type(var_desc_t),        intent(in)    :: v
      if (reg%nvar >= MAX_OUTPUT_VARS) error stop 'ed_output_registry: MAX_OUTPUT_VARS exceeded'
      if (find_var_index(reg, v%name) > 0_ik)                                                    &
         error stop 'ed_output_registry: duplicate variable name ('//trim(v%name)//')'
      reg%nvar = reg%nvar + 1_ik
      reg%var(reg%nvar) = v
   end subroutine append_variable

   !---~---
   !   The staged twin of a fast-loop diagnostic row. Each sub-step the fast loop stages its
   ! per-cohort rows, each patch's row and their polygon mean (the area-weighted sum the
   ! coarse tiers reduce a patch row to), and the forcing, so a variable on one of these has
   ! a sub-step value on its own axis; a cohort row reduced to the patch or the polygon is
   ! reduced from the staged row, with the weights the fast loop used. SRC_NONE for any
   ! other source.
   !---~---
   pure integer(ik) function staged_row(dm, src) result(f)
      integer(ik), intent(in) :: dm, src
      f = SRC_NONE
      if ((dm == DIM_COHORT .or. dm == DIM_PATCH .or. dm == DIM_SCALAR) .and. src > FLD_C_DIAG0 .and.   &
          src <= FLD_C_DIAG0 + N_CO_FAST_ROWS) then
         f = SRC_F_CO0 + (src - FLD_C_DIAG0)
      else if ((dm == DIM_PATCH .or. dm == DIM_SCALAR) .and. src > FLD_P_DIAG0 .and.            &
               src <= FLD_P_DIAG0 + N_PA_FAST_ROWS) then
         f = SRC_F_PA0 + (src - FLD_P_DIAG0)
      else if (dm == DIM_SCALAR .and. src > FLD_PY_DIAG0 .and. src <= FLD_PY_DIAG0 + N_PY_ROWS) then
         f = SRC_F_PY0 + (src - FLD_PY_DIAG0)
      end if
   end function staged_row

   !---~---
   !   The snapshots: for each state's mean_<x>, <x>, its value at the end of the period, off
   ! unless asked for (docs/science/diagnostics.md §4). A state of the slow step is read as
   ! the period's last step left it; one the fast loop records each sub-step (the canopy
   ! air, a leaf's temperature or potential) is its value at the period's last fast step,
   ! read from the staging (snap_staged), so on an F tier one sub-step long a snapshot
   ! equals its mean. A flux or a rate is a mean only, and so is each site's forcing, which
   ! the fast loop keeps no sub-step of.
   !---~---
   subroutine register_snapshots(reg)
      type(output_registry_t), intent(inout) :: reg
      type(var_desc_t) :: base, snap
      integer(ik) :: k, n_base
      n_base = reg%nvar
      do k = 1_ik, n_base
         base = reg%var(k)
         if (base%name(1:5) /= 'mean_') cycle
         if (find_var_index(reg, base%name(6:)) > 0_ik) cycle      ! a state's, registered with its mean
         if (base%flux) cycle
         if (base%source_id > FLD_SI_DIAG0 .and. base%source_id <= FLD_SI_DIAG0 + N_SI_ROWS) cycle
         snap = base
         snap%name            = base%name(6:)
         snap%agg             = AGG_LAST
         snap%streams         = FREQ_NONE
         snap%streams_default = FREQ_NONE
         snap%snap_staged     = base%source_id == SRC_NONE .or. fast_row(base%source_id)
         !---~---
         !   The staging holds a sub-step of each source but a cohort row reduced to a coarser
         ! axis, which the snapshot reduces itself from the staged row (staged_snapshot).
         !---~---
         if (snap%snap_staged .and. snap%source_fast == SRC_NONE .and.                             &
             .not. (base%source_id > FLD_C_DIAG0 .and. base%source_id <= FLD_C_DIAG0 + N_CO_FAST_ROWS)) &
            error stop 'ed_output_registry: no staged sub-step for the snapshot '//trim(snap%name)
         call append_variable(reg, snap)
      end do
   end subroutine register_snapshots

   !---~---
   !   A row the fast loop fills each sub-step, which the coarse tiers read as its mean over
   ! the slow step: the per-cohort and per-patch blocks' fast rows and the forcing echo.
   !---~---
   pure logical function fast_row(src)
      integer(ik), intent(in) :: src
      fast_row = (src > FLD_C_DIAG0  .and. src <= FLD_C_DIAG0 + N_CO_FAST_ROWS) .or.                &
                 (src > FLD_P_DIAG0  .and. src <= FLD_P_DIAG0 + N_PA_FAST_ROWS) .or.                &
                 (src > FLD_PY_DIAG0 .and. src <= FLD_PY_DIAG0 + N_PY_ROWS)
   end function fast_row



   !---~---
   !   Guard: no cohort- or patch-dimensioned variable on the annual stream
   ! (docs/science/diagnostics.md §3).
   !---~---
   subroutine enforce_annual_guard(reg)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: k
      do k = 1_ik, reg%nvar
         if (ragged_dim(reg%var(k)%dim) .and. iand(reg%var(k)%streams_default, FREQ_ANNUAL) /= 0_ik) &
            error stop 'ed_output_registry: cohort/patch variable on the annual stream ('//trim(reg%var(k)%name)//')'
      end do
   end subroutine enforce_annual_guard

   !---~---
   !   Guard: the FAST tier and the FAST-staged sources go together, in both directions, or
   ! the variable writes _FillValue for the whole run (docs/science/diagnostics.md §4, How
   ! the FAST tier is sampled). apply_variable_override runs the same check on each override.
   !---~---
   subroutine enforce_fast_guard(reg)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: k
      do k = 1_ik, reg%nvar
         call check_fast_tier(reg%var(k), reg%var(k)%streams)
      end do
   end subroutine enforce_fast_guard

   !---~---
   !   One variable's streams against its sources: the F tier reads a staged sub-step value,
   ! the D, M and Y tiers a field of the state at the step's end, or, for a snapshot that
   ! reads the staging, the last sub-step.
   !---~---
   subroutine check_fast_tier(v, streams)
      type(var_desc_t), intent(in) :: v
      integer(ik),      intent(in) :: streams
      if (iand(streams, FREQ_FAST) /= 0_ik .and. src_class(v%source_fast) /= SRCK_FAST)             &
         error stop 'ed_output_registry: "'//trim(v%name)//'" has no sub-daily value and cannot be '// &
                    'put on the F stream'
      if (iand(streams, not(FREQ_FAST)) /= 0_ik .and. v%source_id == SRC_NONE .and. .not. v%snap_staged) &
         error stop 'ed_output_registry: "'//trim(v%name)//'" exists only sub-daily and can only be '// &
                    'on the F stream'
   end subroutine check_fast_tier



   !---~---
   !   Selection (docs/science/diagnostics.md §5). An io_config line names a variable or a
   ! pattern of names, `*` matching any run of characters, and gives it tiers, true (its
   ! default) or false. Each variable takes the most specific line that matches it: its own
   ! name over any pattern, a pattern with more literal characters over one with fewer. A
   ! name asks for exactly its tiers; a pattern gives each of its variables the tiers it
   ! asks for that the variable can be on (usable_tiers). Two lines as specific as each
   ! other that disagree, a line that matches no variable and a tier none of a line's
   ! variables can be on stop the run. The axes' bookkeeping is not selected: it follows
   ! its axis (resolve_streams).
   !---~---
   subroutine apply_selections(reg, sel)
      type(output_registry_t), intent(inout) :: reg
      type(selection_t),       intent(in)    :: sel(:)
      integer(ik) :: k, e, b, s, s_best
      integer(ik) :: best(reg%nvar), granted(size(sel)), n_won(size(sel))
      logical     :: hit(size(sel)), hit_coord(size(sel))
      character(len=8) :: missing
      hit = .false.
      hit_coord = .false.
      granted = FREQ_NONE
      n_won = 0_ik
      do k = 1_ik, reg%nvar
         best(k) = 0_ik
         s_best = -1_ik
         do e = 1_ik, size(sel, kind=ik)
            if (.not. glob_match(trim(sel(e)%key), trim(reg%var(k)%name))) cycle
            if (reg%var(k)%coord) then
               if (index(sel(e)%key, '*') == 0)                                                     &
                  error stop 'ed_output_registry: io_config "'//trim(sel(e)%key)//'" is written with '// &
                             'its axis, wherever a variable on the axis is, and is not selected'
               hit_coord(e) = .true.
               cycle
            end if
            hit(e) = .true.
            s = specificity(sel(e)%key)
            if (s > s_best) then
               best(k) = e
               s_best = s
            else if (s == s_best) then
               b = best(k)
               if (sel(e)%kind /= sel(b)%kind .or. sel(e)%mask /= sel(b)%mask)                        &
                  error stop 'ed_output_registry: io_config "'//trim(sel(b)%key)//'" and "'//            &
                             trim(sel(e)%key)//'" select '//trim(reg%var(k)%name)//' alike and disagree; '// &
                             'make one more specific'
            end if
         end do
      end do
      do e = 1_ik, size(sel, kind=ik)
         if (hit(e)) cycle
         if (hit_coord(e))                                                                        &
            error stop 'ed_output_registry: io_config "'//trim(sel(e)%key)//'" matches only the '//     &
                       'axes'' bookkeeping, which is written with its axis'
         error stop 'ed_output_registry: io_config "'//trim(sel(e)%key)//'" matches no variable (typo?)'
      end do
      do k = 1_ik, reg%nvar
         e = best(k)
         if (e == 0_ik) cycle
         call select_variable(reg%var(k), sel(e)%kind, sel(e)%mask, index(sel(e)%key, '*') == 0)
         n_won(e) = n_won(e) + 1_ik
         granted(e) = ior(granted(e), reg%var(k)%streams)
      end do
      do e = 1_ik, size(sel, kind=ik)
         if (sel(e)%kind /= OVR_MASK .or. n_won(e) == 0_ik) cycle
         missing = mask_to_string(iand(sel(e)%mask, not(granted(e))))
         if (len_trim(missing) > 0)                                                                &
            error stop 'ed_output_registry: io_config "'//trim(sel(e)%key)//'": none of its '//         &
                       'variables can be on '//trim(missing)
      end do
   end subroutine apply_selections

   !---~---
   !   One line's choice for one variable. A name asks for exactly its tiers, held to the
   ! variable's sources and the annual rule; a pattern gets the tiers it can be on. `true`
   ! on a name off by default is refused: there are no default tiers to restore.
   !---~---
   subroutine select_variable(v, kind, mask, named)
      type(var_desc_t), intent(inout) :: v
      integer(ik),      intent(in)    :: kind   ! OVR_TRUE | OVR_FALSE | OVR_MASK
      integer(ik),      intent(in)    :: mask   ! used for OVR_MASK
      logical,          intent(in)    :: named  ! the line is the variable's own name, not a pattern
      select case (kind)
      case (OVR_TRUE)
         if (named .and. v%streams_default == FREQ_NONE)                                        &
            error stop 'ed_output_registry: "'//trim(v%name)//'" is off by default, so true has no '// &
                       'tiers to restore; name them, e.g. "D M"'
         v%streams = v%streams_default
      case (OVR_FALSE)
         v%streams = FREQ_NONE
      case (OVR_MASK)
         if (named) then
            if (ragged_dim(v%dim) .and. iand(mask, FREQ_ANNUAL) /= 0_ik)                          &
               error stop 'ed_output_registry: annual (Y) stream forbidden on cohort/patch variable ('// &
                          trim(v%name)//')'
            call check_fast_tier(v, mask)
            v%streams = mask
         else
            v%streams = iand(mask, usable_tiers(v))
         end if
      end select
   end subroutine select_variable

   !---~---
   !   One named variable's override, found = .false. for a name no variable has (the
   ! caller decides what that means); the axes' bookkeeping is not selected.
   !---~---
   subroutine apply_variable_override(reg, name, kind, mask, found)
      type(output_registry_t), intent(inout) :: reg
      character(len=*),        intent(in)    :: name
      integer(ik),             intent(in)    :: kind   ! OVR_TRUE | OVR_FALSE | OVR_MASK
      integer(ik),             intent(in)    :: mask   ! used for OVR_MASK
      logical,                 intent(out)   :: found
      integer(ik) :: k
      k = find_var_index(reg, name)
      found = (k > 0_ik)
      if (.not. found) return
      if (reg%var(k)%coord)                                                                       &
         error stop 'ed_output_registry: "'//trim(name)//'" is written with its axis and is not selected'
      call select_variable(reg%var(k), kind, mask, .true.)
   end subroutine apply_variable_override

   !---~---
   !   The tiers a variable can be on: F with a staged sub-step value, D and M with a value
   ! at the step's end or a staged snapshot, and Y too off the cohort and patch axes.
   !---~---
   pure integer(ik) function usable_tiers(v) result(m)
      type(var_desc_t), intent(in) :: v
      m = FREQ_NONE
      !--- The axes' bookkeeping is fixed within a slow step: the F tier reads it from the stand.
      if (src_class(v%source_fast) == SRCK_FAST .or. v%coord) m = FREQ_FAST
      if (v%source_id /= SRC_NONE .or. v%snap_staged) then
         m = ior(m, ior(FREQ_DAILY, FREQ_MONTHLY))
         if (.not. ragged_dim(v%dim)) m = ior(m, FREQ_ANNUAL)
      end if
   end function usable_tiers

   !---~---
   !   How specific an io_config line is: a name beats every pattern, and a pattern counts
   ! its literal characters.
   !---~---
   pure integer(ik) function specificity(key) result(s)
      character(len=*), intent(in) :: key
      integer :: i
      if (index(key, '*') == 0) then
         s = huge(1_ik)
         return
      end if
      s = 0_ik
      do i = 1, len_trim(key)
         if (key(i:i) /= '*') s = s + 1_ik
      end do
   end function specificity

   !--- Does name match pat, where `*` matches any run of characters, none included?
   pure logical function glob_match(pat, name) result(ok)
      character(len=*), intent(in) :: pat, name
      integer :: p, i, star, mark
      p = 1
      i = 1
      star = 0
      mark = 0
      do while (i <= len(name))
         if (p <= len(pat)) then
            if (pat(p:p) == '*') then
               star = p
               mark = i
               p = p + 1
               cycle
            else if (pat(p:p) == name(i:i)) then
               p = p + 1
               i = i + 1
               cycle
            end if
         end if
         if (star == 0) then
            ok = .false.
            return
         end if
         p = star + 1                             ! let the last * take one more character
         mark = mark + 1
         i = mark
      end do
      do while (p <= len(pat))
         if (pat(p:p) /= '*') exit
         p = p + 1
      end do
      ok = p > len(pat)
   end function glob_match

   !---~---
   !   The tiers each variable is written at, from what selected it: a tier switched off
   ! writes nothing, and the axes' bookkeeping (coord) is written at each tier a variable on
   ! its axis is, the soil-by-patch columns counting as the patch axis. Then each tier's
   ! live list. Run after any change of selection.
   !---~---
   subroutine resolve_streams(reg)
      type(output_registry_t), intent(inout) :: reg
      integer(ik) :: k, t, on, cohort_tiers, patch_tiers
      on = FREQ_NONE
      do t = 1_ik, N_FREQ
         if (reg%tier_on(t)) on = ior(on, freq_bit(t))
      end do
      cohort_tiers = FREQ_NONE
      patch_tiers = FREQ_NONE
      do k = 1_ik, reg%nvar
         reg%var(k)%streams = iand(reg%var(k)%streams, on)
         if (reg%var(k)%coord .or. .not. reg%var(k)%enabled) cycle
         select case (reg%var(k)%dim)
         case (DIM_COHORT)
            cohort_tiers = ior(cohort_tiers, reg%var(k)%streams)
         case (DIM_PATCH, DIM_SOIL_PATCH)
            patch_tiers = ior(patch_tiers, reg%var(k)%streams)
         end select
      end do
      do k = 1_ik, reg%nvar
         if (.not. reg%var(k)%coord) cycle
         if (reg%var(k)%dim == DIM_COHORT) then
            reg%var(k)%streams = iand(cohort_tiers, usable_tiers(reg%var(k)))
         else
            reg%var(k)%streams = iand(patch_tiers, usable_tiers(reg%var(k)))
         end if
      end do
      call build_freq_index(reg)
   end subroutine resolve_streams

   !--- 1-based registry index of a variable by name; 0 if absent.
   pure integer(ik) function find_var_index(reg, name) result(idx)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik) :: k
      idx = 0_ik
      do k = 1_ik, reg%nvar
         if (trim(reg%var(k)%name) == trim(name)) then
            idx = k
            return
         end if
      end do
   end function find_var_index

   !---~---
   !   Parse a "F D M Y" stream-string into an ior(FREQ_*) mask (the value grammar of
   ! docs/science/diagnostics.md §5).
   ! status = 0 ok; status = 1 unrecognized token (the offender returned in bad_token).
   !---~---
   subroutine parse_stream_mask(s, mask, status, bad_token)
      character(len=*), intent(in)  :: s
      integer(ik),      intent(out) :: mask, status
      character(len=*), intent(out) :: bad_token
      integer(ik) :: i
      character(len=1) :: c
      mask = FREQ_NONE
      status = 0_ik
      bad_token = ''
      do i = 1_ik, int(len_trim(s), ik)
         c = s(i:i)
         select case (c)
         case ('F','f')
            mask = ior(mask, FREQ_FAST)
         case ('D','d')
            mask = ior(mask, FREQ_DAILY)
         case ('M','m')
            mask = ior(mask, FREQ_MONTHLY)
         case ('Y','y')
            mask = ior(mask, FREQ_ANNUAL)
         case (' ')     ! separator
         case default
            status = 1_ik
            bad_token = c
            return
         end select
      end do
   end subroutine parse_stream_mask

   !--- Last: precompute, per tier, the list of live (enabled + in-tier) variable indices.
   subroutine build_freq_index(reg)
      type(output_registry_t), intent(inout) :: reg
      integer(ik) :: t, k, bit
      reg%idx_freq = 0_ik
      reg%nidx = 0_ik
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         do k = 1_ik, reg%nvar
            if (reg%var(k)%enabled .and. iand(reg%var(k)%streams, bit) /= 0_ik) then
               reg%nidx(t) = reg%nidx(t) + 1_ik
               reg%idx_freq(reg%nidx(t), t) = k
            end if
         end do
      end do
   end subroutine build_freq_index

   !---~---
   !   Write a commented-out ed_io_config.toml listing every registered variable with its
   ! default streams, group and axis: the names the overrides accept, generated from the
   ! registry so it cannot drift from it.
   !---~---
   !---~---
   !   Write a commented-out ed_io_config.toml listing every variable a selection can name,
   ! by group: its default tiers, the tiers it can be on, its axis, long name and units;
   ! then the axes' bookkeeping. Generated from the registry, so it cannot drift from it.
   !---~---
   subroutine dump_io_config(reg, path)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: path
      integer     :: u
      integer(ik) :: k, g
      character(len=52) :: lhs      ! left-hand side, blank-padded so the comments line up
      character(len=8)  :: can
      open(newunit=u, file=path, status='replace', action='write')
      write(u,'(a)') '# ED-MORE output selection -- GENERATED by `edmore --dump-io-config`.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Point the main config at a file like this with [output].io_config = "ed_io_config.toml".'
      write(u,'(a)') '# A line names a variable, or a pattern of names with * (quoted: "mean_*_py"), and what is'
      write(u,'(a)') '# written:'
      write(u,'(a)') '#   name = "F D M Y" -> at these tiers: F sub-daily, D daily, M monthly, Y yearly'
      write(u,'(a)') '#   name = true      -> at its default tiers, below'
      write(u,'(a)') '#   name = false     -> nowhere'
      write(u,'(a)') '# A variable takes the most specific line that matches it: its own name over any pattern, a'
      write(u,'(a)') '# pattern with more literal characters over one with fewer; with none, its default. A name'
      write(u,'(a)') '# may ask only for the tiers it can be on (below; a cohort or patch variable is never on Y);'
      write(u,'(a)') '# a pattern gives each of its variables the tiers it asks that the variable can be on. Two'
      write(u,'(a)') '# lines as specific as each other that disagree, a line that matches nothing, a tier none'
      write(u,'(a)') '# of a line''s variables can be on, and `true` for a name off by default stop the run. A tier'
      write(u,'(a)') '# whose [output.<tier>].enabled is false writes nothing, whatever this file says.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Each line below: the default (false: off), the tiers it can be on, [axis], long name (units).'
      write(u,'(a)') '# Every line is commented out, so this file as generated changes nothing.'
      write(u,'(a)') ''
      write(u,'(a)') '[variables]'
      do g = 1_ik, N_GRP
         write(u,'(a)') ''
         write(u,'(a)') '# ---- '//trim(group_name(g))//' '//repeat('-', 90 - len_trim(group_name(g)))
         do k = 1_ik, reg%nvar
            if (reg%var(k)%group /= g .or. reg%var(k)%coord) cycle
            if (reg%var(k)%streams_default == FREQ_NONE) then
               lhs = '# '//trim(reg%var(k)%name)//' = false'
            else
               lhs = '# '//trim(reg%var(k)%name)//' = "'//trim(mask_to_string(reg%var(k)%streams_default))//'"'
            end if
            can = mask_to_string(usable_tiers(reg%var(k)))
            write(u,'(a)') lhs//'  # '//can//'  ['//trim(dim_name(reg%var(k)%dim))//'] '//                &
                           trim(reg%var(k)%long_name)//' ('//trim(reg%var(k)%units)//')'
         end do
      end do
      write(u,'(a)') ''
      write(u,'(a)') '# ---- the axes'' bookkeeping, written with its axis at each tier a variable on it is ---------'
      do k = 1_ik, reg%nvar
         if (.not. reg%var(k)%coord) cycle
         write(u,'(a)') '#   '//reg%var(k)%name(1:20)//'  ['//trim(dim_name(reg%var(k)%dim))//'] '//            &
                        trim(reg%var(k)%long_name)
      end do
      close(u)
      write(*,'(2a)') ' wrote ', trim(path)
   end subroutine dump_io_config

   pure function mask_to_string(mask) result(s)
      integer(ik), intent(in) :: mask
      character(len=8) :: s
      integer :: p
      s = ''
      p = 0
      if (iand(mask, FREQ_FAST)    /= 0_ik) then
         p = p+1
         s(p:p) = 'F'
         p = p+1
         s(p:p) = ' '
      end if
      if (iand(mask, FREQ_DAILY)   /= 0_ik) then
         p = p+1
         s(p:p) = 'D'
         p = p+1
         s(p:p) = ' '
      end if
      if (iand(mask, FREQ_MONTHLY) /= 0_ik) then
         p = p+1
         s(p:p) = 'M'
         p = p+1
         s(p:p) = ' '
      end if
      if (iand(mask, FREQ_ANNUAL)  /= 0_ik) then
         p = p+1
         s(p:p) = 'Y'
         p = p+1
         s(p:p) = ' '
      end if
   end function mask_to_string

   pure function group_name(g) result(s)
      integer(ik), intent(in) :: g
      character(len=10) :: s
      select case (g)
      case (GRP_STRUCTURE)
         s = 'structure'
      case (GRP_CARBON)
         s = 'carbon'
      case (GRP_WATER)
         s = 'water'
      case (GRP_ENERGY)
         s = 'energy'
      case (GRP_NUMERICS)
         s = 'numerics'
      case (GRP_RADIATION)
         s = 'radiation'
      case (GRP_ECOPHYS)
         s = 'ecophys'
      case (GRP_BIOGEOCHEM)
         s = 'biogeochem'
      case (GRP_FORCING)
         s = 'forcing'
      case default
         s = '?'
      end select
   end function group_name

   pure function dim_name(d) result(s)
      integer(ik), intent(in) :: d
      character(len=10) :: s
      select case (d)
      case (DIM_SCALAR)
         s = 'polygon'
      case (DIM_COHORT)
         s = 'cohort'
      case (DIM_PATCH)
         s = 'patch'
      case (DIM_SOIL)
         s = 'soil'
      case (DIM_PFT)
         s = 'pft'
      case (DIM_SIZE)
         s = 'dbh_class'
      case (DIM_SOIL_PATCH)
         s = 'soil,patch'
      case (DIM_SITE)
         s = 'by site'
      case default
         s = '?'
      end select
   end function dim_name

end module ed_output_registry
