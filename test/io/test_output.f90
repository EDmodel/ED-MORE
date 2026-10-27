! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_output -- the diagnostic output's netCDF-free half (src/io/output;
! docs/science/diagnostics.md), one section per module:
!   1. THE REGISTRY AND THE SELECTION (ed_output_registry: build_output_registry,
!      apply_selections, resolve_streams, glob_match, apply_variable_override,
!      parse_stream_mask; §5): the registration list and its defaults, the io_config lines
!      (names and patterns, the most specific winning), the axes' bookkeeping, the tier
!      switches, the unknown-key trap, and the snapshot each state has.
!   2. THE REDUCTION (ed_output_reduce: reduce_cohort_to_polygon, reduce_cohort_to_patch,
!      reduce_cohort_to_pft, reduce_cohort_to_size, reduce_patch_columns_to_sites,
!      gather_patch_columns, dbh_class_index; §2 eqs. (1)-(3)): the weighted aggregation,
!      extensive against intensive, the PFT and size-class closure, empty sets, unsampled
!      cohorts, the class edges and the soil columns.
!   3. THE EXTRACT (ed_output_extract: extract_variable, extract_fast_scalar): a
!      variable's live slab or site value from the state, and a sub-step's row from the
!      staging.
!   4. THE INTEGRATION (ed_output_integrate: integrate_scalar, integrate_slab,
!      normalize_scalar, normalize_slab, output_integrate_fast, close_tier,
!      output_integrate; with ed_output_setup's sizing; §4): the operators and the
!      zero-sample guard, a slab over a fixed-slot window, the FAST tier, a slab sized
!      after the overrides, two polygons' buffers, the month boundary and the snapshots.
!
!   What section 2 is FOR. The reducer's whole job is to apply the right weight to the right
! quantity, and the two ways that goes wrong are both silent: an EXTENSIVE quantity
! aggregated as a mean (or an INTENSIVE one as a sum) produces a plausible number with no
! error anywhere, and an empty bin reported as 0 is indistinguishable in the output file from
! a real measurement of zero. So its assertions are chosen so that a wrong weight cannot pass:
! the fixture uses DIFFERENT nplant per cohort and different values per PFT and size class,
! which makes the sum, the stem-weighted mean and the basal-area-weighted mean three distinct
! numbers.
!
!   nvfortran note: the slab tests mutate the site's agb through the OPAQUE setter below (a
! separate module scope) between steps -- mirroring production, where state changes happen
! deep inside commit_vegetation_step. Mutating an allocatable component INLINE between two
! extract_variable() calls lets nvfortran -O2 wrongly CSE the reads (an allocatable-component
! aliasing blind spot); an opaque mutation defeats it, as production does.
!---~---
module test_output_support
   use ed_kinds,                only : wp, ik
   use ed_polygon_state_vars,   only : polygon_t
   implicit none
contains
   !--- Opaque per-cohort agb setter (separate compilation scope; nvfortran cannot CSE across it).
   subroutine set_cohort_agb(poly, vals, n)
      type(polygon_t), intent(inout) :: poly
      real(wp),        intent(in)    :: vals(:)
      integer(ik),     intent(in)    :: n
      poly%cohort%agb(1:n) = vals(1:n)
   end subroutine set_cohort_agb
end module test_output_support

program test_output
   use test_output_support,     only : set_cohort_agb
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_constants,            only : n_soil_layer_max
   use ed_polygon_state_vars,   only : polygon_t, alloc_polygon, free_polygon, rebuild_patch_map
   use ed_output_params,        only : FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE, N_FREQ
   use ed_output_types,         only : output_registry_t, var_desc_t, integ_buffer_t, output_files_t,          &
                                       output_buffers_t, diag_setup_t, slab_col, MAX_OUTPUT_VARS,               &
                                       MISSING_VALUE, AGG_MEAN, AGG_RATE, AGG_LAST, AGG_TMEAN, DIM_SCALAR,      &
                                       DIM_COHORT, DIM_PATCH
   use ed_output_registry,      only : build_output_registry, build_freq_index, find_var_index,                &
                                       apply_variable_override, parse_stream_mask, selection_t,                 &
                                       apply_selections, resolve_streams, glob_match, OVR_TRUE, OVR_FALSE,      &
                                       OVR_MASK
   use ed_output_reduce,        only : reduce_cohort_to_polygon, reduce_cohort_to_patch, reduce_cohort_to_pft,  &
                                       reduce_cohort_to_size, reduce_patch_columns_to_sites,                    &
                                       gather_patch_columns, dbh_class_index, W_NONE, W_NPLANT, W_LEAF_AREA,    &
                                       W_BASAL_AREA, total_agb, total_lai, total_nplant
   use ed_output_extract,       only : extract_variable, extract_fast_scalar, FLD_C_AGB, FLD_C_DIAG0,           &
                                       SRC_F_PA0, SRC_F_PY0, SRC_F_CO0, SRC_NONE
   use ed_output_integrate,     only : alloc_integ_buffer, reset_buffer, integrate_scalar, integrate_slab,      &
                                       normalize_scalar, normalize_slab, output_integrate_fast, close_tier,     &
                                       output_integrate
   use ed_output_setup,         only : alloc_output_set, alloc_output_buffers, setup_output_files,              &
                                       finalize_output_files
   use ed_diag_vars,            only : alloc_cohort_diag, CO_LEAF_TEMP, CO_GSW, N_PY_ROWS, PY_TAIR, PY_CO2,      &
                                       N_PA_ROWS, N_CO_FAST_ROWS, PA_CAS_TEMP, PA_LE, PA_H, PA_GPP, PA_SW_UP,   &
                                       PA_LW_UP, PA_CAS_CO2
   use ed_time,                 only : ed_time_t, time_advance_days
   use ed_test_support,         only : banner, build_test_config, test_pft_table, check, check_close
   implicit none

   call banner('output')
   !--- 1. The registry and the selection.
   call test_registry()
   call test_snapshot_registry()
   !--- 2. The reduction.
   call test_extensive_vs_intensive()
   call test_pft_and_size_closure()
   call test_empty_sets()
   call test_unsampled_cohorts()
   call test_class_edges()
   call test_soil_columns()
   !--- 3. The extract.
   call test_extract()
   call test_fast_extract()
   !--- 4. The integration.
   call test_scalar_operators()
   call test_zero_sample_guard()
   call test_slab_window()
   call test_fast_tier()
   call test_slab_sized_after_overrides()
   call test_soil_patch_record()
   call test_two_buffers()
   call test_boundary_step()
   call test_snapshots()
   write(*,'(a)') 'test_output: ALL PASSED'

contains

   !---~---
   !   1. THE REGISTRY AND THE SELECTION (ed_output_registry; docs/science/diagnostics.md
   ! §5): build_output_registry's list and defaults, then the io_config lines through
   ! apply_selections and resolve_streams, the tier switches, apply_variable_override and
   ! parse_stream_mask.
   !---~---
   subroutine test_registry()
      type(ed_config_t)       :: cfg
      type(output_registry_t) :: reg
      integer(ik)             :: mask, status
      character(len=8)        :: bad
      logical                 :: found

      cfg = build_test_config(86400.0_wp)     ! 1-day slow step; cfg%output holds the [output] defaults

      !--- Build with the defaults (freq_on=[F,T,T,T]).
      call build_output_registry(reg, cfg)
      !---~---
      !   Size sanity, not an exact count. A hard-coded total has to be edited every time a
      ! variable is added, which trains the reader to update the number without checking what
      ! changed -- the assertion then guards nothing. Bound it instead, and let the NAMED
      ! presence checks below carry the real meaning.
      !---~---
      call check(reg%nvar > 80_ik .and. reg%nvar <= MAX_OUTPUT_VARS,                                 &
                 'registry populated and within MAX_OUTPUT_VARS')
      call check(find_var_index(reg, 'mean_agb_co') > 0_ik, 'mean_agb_co present')
      call check(find_var_index(reg, 'mean_agb_py') > 0_ik, 'mean_agb_py present')
      call check(find_var_index(reg, 'mean_gpp_py') > 0_ik, 'mean_gpp_py (carbon) present')
      call check(find_var_index(reg, 'mean_soilc_slow_py') > 0_ik, 'mean_soilc_slow_py (carbon) present')
      call check(find_var_index(reg, 'mean_rh_py') > 0_ik, 'mean_rh_py (carbon) present')
      call check(find_var_index(reg, 'mean_soil_temp_layer_si') > 0_ik, 'mean_soil_temp_layer_si (energy) present')
      call check(find_var_index(reg, 'nope') == 0_ik, 'unknown name -> 0')
      !---~---
      !   The defaults (docs/science/diagnostics.md §5): the radiation and ecophys groups, the
      ! soil-by-patch columns and the snapshots are off until a selection names them; the
      ! axes' bookkeeping is written with the default cohort and patch variables.
      !---~---
      call check(in_tier(reg, 'mean_soil_temp_layer_si', 3_ik),                                     &
                 'energy on by default -> mean_soil_temp_layer_si in MONTHLY')
      call check(.not. live(reg, 'mean_gsw_co'), 'the ecophys group off by default')
      call check(.not. live(reg, 'mean_abs_sw_co'), 'the radiation group off by default')
      call check(.not. live(reg, 'mean_soil_water_layer_pa'), 'the soil by patch off by default')
      call check(.not. live(reg, 'agb_py'), 'a snapshot off by default')
      call check(in_tier(reg, 'global_cohort_id', 2_ik) .and. in_tier(reg, 'global_patch_id', 3_ik),         &
                 'the bookkeeping with the default cohort and patch variables')

      !---~---
      !   FAST tier off by default. The ANNUAL tier must contain NO cohort/patch-dimensioned
      ! variable (docs/science/diagnostics.md §3) -- that invariant is what the check below
      ! states, rather than a count that has to be re-derived on every registry edit.
      !---~---
      call check(reg%nidx(1) == 0_ik, 'FAST tier empty by default (freq_on(FAST)=false)')
      call check(reg%nidx(2) > 0_ik,  'DAILY tier populated')
      call check(reg%nidx(4) > 0_ik,  'ANNUAL tier populated')
      call check(annual_is_site_level(reg), 'ANNUAL tier carries no cohort/patch variable')
      !---~---
      !   A pattern: "*_co" = false drops the cohort axis, its bookkeeping included, and nothing
      ! else.
      !---~---
      call choose(reg, cfg, [selection_t('*_co', OVR_FALSE, FREQ_NONE)])
      call check(.not. live(reg, 'mean_agb_co'), '"*_co" = false -> mean_agb_co in no tier')
      call check(.not. live(reg, 'global_cohort_id'), '"*_co" = false -> no cohort bookkeeping')
      call check(in_tier(reg, 'mean_agb_py', 2_ik), '"*_co" = false leaves the polygon variables alone')
      call check(in_tier(reg, 'global_patch_id', 2_ik), '"*_co" = false leaves the patch axis alone')
      !--- The most specific line wins: a name over a pattern, more literal characters over fewer.
      call choose(reg, cfg, [selection_t('*', OVR_FALSE, FREQ_NONE),                                  &
                             selection_t('mean_*_py', OVR_MASK, FREQ_MONTHLY),                        &
                             selection_t('mean_gpp_py', OVR_MASK, ior(FREQ_DAILY, FREQ_ANNUAL))])
      call check(in_tier(reg, 'mean_agb_py', 3_ik) .and. .not. in_tier(reg, 'mean_agb_py', 2_ik),              &
                 '"mean_*_py" wins over "*"')
      call check(in_tier(reg, 'mean_gpp_py', 2_ik) .and. in_tier(reg, 'mean_gpp_py', 4_ik) .and.               &
                 .not. in_tier(reg, 'mean_gpp_py', 3_ik), 'a name wins over a pattern')
      call check(.not. live(reg, 'nplant_co') .and. .not. live(reg, 'owner_patch'),                             &
                 '"*" = false: no other variable, so no bookkeeping')
      !--- A pattern gives each variable the tiers it asks that the variable can be on.
      call choose(reg, cfg, [selection_t('mean_agb_*', OVR_MASK, ior(FREQ_MONTHLY, FREQ_ANNUAL))])
      call check(in_tier(reg, 'mean_agb_py', 4_ik) .and. in_tier(reg, 'mean_agb_co', 3_ik) .and.               &
                 .not. in_tier(reg, 'mean_agb_co', 4_ik), 'a pattern gives Y only off the cohort axis')
      call choose(reg, cfg, [selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('mean_*_co', OVR_MASK, FREQ_FAST)])
      call check(.not. live(reg, 'mean_gsw_co'), 'the F tier switched off writes nothing')
      !--- The bookkeeping follows its axis, tier by tier.
      call choose(reg, cfg, [selection_t('*', OVR_FALSE, FREQ_NONE),                                  &
                             selection_t('mean_agb_co', OVR_MASK, FREQ_MONTHLY)])
      call check(in_tier(reg, 'owner_patch', 3_ik) .and. .not. in_tier(reg, 'owner_patch', 2_ik),              &
                 'the cohort bookkeeping at the cohort variables'' tiers only')
      call check(.not. live(reg, 'cohort_offset'), 'no patch variable, no patch bookkeeping')
      call choose(reg, cfg, [selection_t('*', OVR_FALSE, FREQ_NONE),                                  &
                             selection_t('mean_soil_water_layer_pa', OVR_MASK, FREQ_DAILY)])
      call check(in_tier(reg, 'cohort_offset', 2_ik), 'the soil by patch is on the patch axis')
      !--- `true`: a variable's default; a pattern's true leaves one off by default off.
      call choose(reg, cfg, [selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('*_py', OVR_TRUE, FREQ_NONE)])
      call check(in_tier(reg, 'mean_agb_py', 4_ik) .and. .not. live(reg, 'agb_py'),                            &
                 'a pattern''s true restores each default, off included')
      !--- The patterns (glob_match).
      call check(glob_match('*_py', 'mean_gpp_py') .and. glob_match('mean_*_py', 'mean_gpp_pft_py') .and.      &
                 glob_match('*', 'x') .and. glob_match('a*b*c', 'axxbyyc') .and. glob_match('a*b', 'ab'),       &
                 'glob_match: a * matches any run of characters, none included')
      call check(.not. glob_match('mean_*_co', 'mean_gpp_py') .and. .not. glob_match('a*b', 'abc') .and.       &
                 .not. glob_match('gpp_py', 'mean_gpp_py'), 'glob_match: the literal parts must match')
      call build_output_registry(reg, cfg)

      !--- PFT + DBH-size axes are registered and live at the monthly cadence.
      call check(find_var_index(reg, 'mean_agb_pft_py')  > 0_ik, 'mean_agb_pft_py present')
      call check(find_var_index(reg, 'mean_agb_size_py') > 0_ik, 'mean_agb_size_py present')
      call check(in_tier(reg, 'mean_agb_pft_py',  3_ik), 'mean_agb_pft_py in MONTHLY')
      call check(in_tier(reg, 'mean_agb_size_py', 3_ik), 'mean_agb_size_py in MONTHLY')
      call check(in_tier(reg, 'mean_agb_pft_py',  4_ik), 'mean_agb_pft_py in ANNUAL (site-level axis, allowed)')

      call check(in_tier(reg, 'mean_agb_co', 2_ik),       'mean_agb_co in DAILY')
      call check(.not. in_tier(reg, 'mean_agb_co', 4_ik), 'mean_agb_co NOT in ANNUAL (no cohort dim annually)')
      call check(in_tier(reg, 'mean_agb_py', 4_ik),       'mean_agb_py in ANNUAL')
      call check(in_tier(reg, 'mean_gpp_py', 4_ik),       'mean_gpp_py in ANNUAL')

      !---~---
      !   Per-tier enable: turn DAILY off -> DAILY empties, MONTHLY unaffected, and a selection
      ! cannot put a variable back on it.
      !---~---
      cfg%output%freq_on(2) = .false.
      call build_output_registry(reg, cfg)
      call check(reg%nidx(2) == 0_ik, 'DAILY off -> empty')
      call check(reg%nidx(3) > 0_ik,  'MONTHLY still populated')
      call choose(reg, cfg, [selection_t('mean_agb_py', OVR_MASK, ior(FREQ_DAILY, FREQ_MONTHLY))])
      call check(reg%nidx(2) == 0_ik, 'DAILY off -> a selection writes nothing there')
      cfg%output%freq_on(2) = .true.

      !--- Per-variable override (apply_variable_override): OFF removes it from every tier.
      call build_output_registry(reg, cfg)
      call apply_variable_override(reg, 'mean_growth_avg_co', OVR_FALSE, 0_ik, found)
      call check(found, 'override found mean_growth_avg_co')
      call build_freq_index(reg)
      call check(.not. in_tier(reg, 'mean_growth_avg_co', 3_ik), 'mean_growth_avg_co OFF -> absent from MONTHLY')

      !--- Per-variable override: MASK re-targets streams (mean_basal_area_co M -> D M).
      call build_output_registry(reg, cfg)
      call check(.not. in_tier(reg, 'mean_basal_area_co', 2_ik), 'mean_basal_area_co default NOT in DAILY')
      call parse_stream_mask('D M', mask, status, bad)
      call check(status == 0_ik, 'parse_stream_mask("D M") ok')
      call apply_variable_override(reg, 'mean_basal_area_co', OVR_MASK, mask, found)
      call build_freq_index(reg)
      call check(in_tier(reg, 'mean_basal_area_co', 2_ik), 'mean_basal_area_co override -> now in DAILY')

      !--- Unknown-key trap: a typo matches no registry variable.
      call apply_variable_override(reg, 'grwoth_avg_cohort', OVR_FALSE, 0_ik, found)
      call check(.not. found, 'typo key -> not found (unknown-key trap)')

      !--- parse_stream_mask flags an unrecognized token.
      call parse_stream_mask('D X', mask, status, bad)
      call check(status == 1_ik .and. trim(bad) == 'X', 'bad stream token flagged')
   end subroutine test_registry

   !---~---
   !   The snapshots (docs/science/diagnostics.md §4): each state's mean_<x> has its <x>, off by
   ! default, and a flux or a rate none; a slow-step state's snapshot reads the state, a fast-loop
   ! state's the staging; and a fast-loop row is staged, so it may be on the F tier.
   !---~---
   subroutine test_snapshot_registry()
      type(ed_config_t)    :: cfg
      type(output_files_t) :: files
      integer(ik) :: k
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled = .true.
      call setup_output_files(files, cfg)
      k = find_var_index(files%reg, 'cas_temp_py')
      call check(k > 0_ik, 'mean_cas_temp_py has its snapshot cas_temp_py')
      if (k > 0_ik) then
         call check(files%reg%var(k)%agg == AGG_LAST .and. files%reg%var(k)%streams == FREQ_NONE        &
                    .and. files%reg%var(k)%snap_staged, 'cas_temp_py: a snapshot, off, read from the staging')
      end if
      k = find_var_index(files%reg, 'dbh_co')
      call check(k > 0_ik, 'mean_dbh_co has its snapshot dbh_co')
      if (k > 0_ik) call check(.not. files%reg%var(k)%snap_staged, 'dbh_co, a state: read from the state')
      k = find_var_index(files%reg, 'height_co')
      if (k > 0_ik) call check(.not. files%reg%var(k)%snap_staged,                                    &
                               'height_co: the state at the step''s end, though its F tier is staged')
      call check(find_var_index(files%reg, 'gpp_py') == 0_ik .and. find_var_index(files%reg, 'le_py') == 0_ik &
                 .and. find_var_index(files%reg, 'mort_rate_co') == 0_ik .and.                            &
                 find_var_index(files%reg, 'air_temp_si') == 0_ik,                                     &
                 'a flux, a rate and a site''s forcing are means only')
      call check(find_var_index(files%reg, 'mean_nplant_co') > 0_ik .and.                              &
                 find_var_index(files%reg, 'mean_area_pa') > 0_ik, 'a state snapshot has its mean')
      k = find_var_index(files%reg, 'mean_gsw_co')
      call check(files%reg%var(k)%source_fast == SRC_F_CO0 + CO_GSW,                                  &
                 'a fast-loop cohort row is staged, so the F tier reads it')
      k = find_var_index(files%reg, 'cas_co2_pa')
      call check(k > 0_ik, 'an F-only state has a snapshot')
      if (k > 0_ik) call check(files%reg%var(k)%source_id == SRC_NONE .and. files%reg%var(k)%snap_staged, &
                               'cas_co2_pa reads the staging on every tier')
   end subroutine test_snapshot_registry

   !--- The registry from the config, with these io_config lines applied.
   subroutine choose(reg, cfg, sel)
      type(output_registry_t), intent(inout) :: reg
      type(ed_config_t),       intent(in)    :: cfg
      type(selection_t),       intent(in)    :: sel(:)
      call build_output_registry(reg, cfg)
      call apply_selections(reg, sel)
      call resolve_streams(reg)
   end subroutine choose

   !--- .true. if the variable is live in some tier.
   logical function live(reg, name) result(yes)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik) :: t
      yes = .false.
      do t = 1_ik, N_FREQ
         if (in_tier(reg, name, t)) yes = .true.
      end do
   end function live

   !---~---
   !   The annual stream is SITE-LEVEL ONLY (docs/science/diagnostics.md §3). A cohort/patch
   ! window longer than a month would straddle the annual disturbance restructuring, so the
   ! slot set it averaged would not be the slot set present at flush.
   !---~---
   logical function annual_is_site_level(reg) result(ok)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: j, k
      ok = .true.
      do j = 1_ik, reg%nidx(4_ik)
         k = reg%idx_freq(j, 4_ik)
         if (reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) ok = .false.
      end do
   end function annual_is_site_level

   !--- .true. if the variable is on the tier.
   logical function in_tier(reg, name, tier) result(yes)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik),             intent(in) :: tier
      integer(ik) :: idx, j
      idx = find_var_index(reg, name)
      yes = .false.
      do j = 1_ik, reg%nidx(tier)
         if (reg%idx_freq(j, tier) == idx) then
            yes = .true.
            return
         end if
      end do
   end function in_tier

   !--- 2. THE REDUCTION (ed_output_reduce; docs/science/diagnostics.md §2, eqs. (1)-(3)).

   !---~---
   !   Two patches, three cohorts, deliberately unequal nplant / basal_area so that a SUM, a
   ! stem-weighted MEAN and a basal-area-weighted MEAN are three different numbers.
   !---~---
   subroutine build_fixture(poly)
      type(polygon_t), intent(out) :: poly
      call alloc_polygon(poly, 3_ik, 16_ik, 4_ik, 4_ik, test_pft_table())
      poly%n_pft = 3_ik
      !--- patch 1 (area 0.75): cohorts 1,2 ; patch 2 (area 0.25): cohort 3.
      poly%patch%n = 2_ik
      poly%patch%area(1) = 0.75_wp
      poly%patch%area(2) = 0.25_wp
      poly%patch%cohort_offset(1) = 1_ik
      poly%patch%cohort_count(1)  = 2_ik
      poly%patch%cohort_offset(2) = 3_ik
      poly%patch%cohort_count(2)  = 1_ik
      poly%cohort%n = 3_ik
      poly%cohort%owner_patch(1:3) = [1_ik, 1_ik, 2_ik]
      poly%cohort%pft(1:3)         = [1_ik, 2_ik, 3_ik]
      poly%cohort%nplant(1:3)      = [2.0_wp,  1.0_wp,  4.0_wp]
      poly%cohort%agb(1:3)         = [10.0_wp, 20.0_wp, 30.0_wp]
      poly%cohort%leaf_area(1:3)   = [ 1.0_wp,  3.0_wp,  5.0_wp]
      poly%cohort%basal_area(1:3)  = [ 5.0_wp,  1.0_wp,  2.0_wp]
      poly%cohort%dbh(1:3)         = [ 5.0_wp, 15.0_wp, 45.0_wp]   ! -> classes 1, 2, 4
      call rebuild_patch_map(poly)                                 ! both patches in the one site
   end subroutine build_fixture

   !---~---
   !   EXTENSIVE vs INTENSIVE (reduce_cohort_to_polygon, reduce_cohort_to_patch). The same
   ! field, reduced with the two different contracts, must give the two different right
   ! answers -- and swapping the weight must CHANGE the answer, which is what makes the
   ! assertion able to fail.
   !---~---
   subroutine test_extensive_vs_intensive()
      type(polygon_t) :: poly
      real(wp) :: v, v2, out(8)
      logical  :: ok, vld(8)
      integer(ik) :: n
      call build_fixture(poly)

      !---~---
      !   EXTENSIVE: mean_agb_py = SUM area*nplant*agb
      ! = 0.75*(2*10 + 1*20) + 0.25*(4*30) = 0.75*40 + 0.25*120 = 30 + 30 = 60.
      !---~---
      call reduce_cohort_to_polygon(poly, poly%cohort%agb, W_NPLANT, .false., v, ok)
      call check_close(v, 60.0_wp, 1.0e-12_wp, 'mean_agb_py = weighted SUM over nplant')
      call check(ok, 'extensive reduction is always valid')
      call check_close(total_agb(poly), 60.0_wp, 1.0e-12_wp, 'total_agb wrapper agrees')

      !---~---
      !   The same call with mean=.true. is a DIFFERENT number (mean per-plant agb):
      ! 60 / (0.75*3 + 0.25*4) = 60/3.25.
      !---~---
      call reduce_cohort_to_polygon(poly, poly%cohort%agb, W_NPLANT, .true., v2, ok)
      call check_close(v2, 60.0_wp/3.25_wp, 1.0e-12_wp, 'mean=.true. gives the per-plant mean')
      call check(abs(v - v2) > 1.0_wp, 'sum and mean are genuinely different (a wrong flag cannot pass)')

      !---~---
      !   INTENSIVE with a WEIGHT that matters: basal-area-weighted mean dbh
      ! w = nplant*basal_area = [10, 1, 8]; area-weighted:
      ! num = 0.75*(10*5 + 1*15) + 0.25*(8*45) = 0.75*65 + 0.25*360 = 48.75 + 90 = 138.75
      ! den = 0.75*(10 + 1)      + 0.25*8      = 8.25 + 2 = 10.25
      !---~---
      call reduce_cohort_to_polygon(poly, poly%cohort%dbh, W_BASAL_AREA, .true., v, ok)
      call check_close(v, 138.75_wp/10.25_wp, 1.0e-12_wp, 'BA-weighted mean dbh')
      !---~---
      !   The STEM-weighted mean of the same field is a different number, so the weight kind
      ! is doing real work rather than being decorative.
      !---~---
      call reduce_cohort_to_polygon(poly, poly%cohort%dbh, W_NPLANT, .true., v2, ok)
      call check(abs(v - v2) > 1.0_wp, 'BA-weighted /= stem-weighted mean dbh')

      !--- LEAF-AREA weighting, the contract canopy temperatures and conductances use.
      call reduce_cohort_to_polygon(poly, poly%cohort%dbh, W_LEAF_AREA, .true., v, ok)
      call check(abs(v - v2) > 0.1_wp, 'leaf-area-weighted /= stem-weighted mean dbh')

      !--- W_NONE + SUM is the stem-density contract: 0.75*(2+1) + 0.25*4 = 2.25 + 1 = 3.25.
      call reduce_cohort_to_polygon(poly, poly%cohort%nplant, W_NONE, .false., v, ok)
      call check_close(v, 3.25_wp, 1.0e-12_wp, 'mean_nplant_py = plain area-weighted sum')
      call check_close(total_nplant(poly), 3.25_wp, 1.0e-12_wp, 'total_nplant wrapper agrees')

      !---~---
      !   PATCH axis carries NO area factor (a patch value is per m2 of its OWN ground):
      ! patch 1 agb = 2*10 + 1*20 = 40 ; patch 2 = 4*30 = 120.
      !---~---
      call reduce_cohort_to_patch(poly, poly%cohort%agb, W_NPLANT, .false., out, vld, n)
      call check(n == 2_ik, 'patch reduction length')
      call check_close(out(1),  40.0_wp, 1.0e-12_wp, 'mean_agb_pa(1) has no area factor')
      call check_close(out(2), 120.0_wp, 1.0e-12_wp, 'mean_agb_pa(2) has no area factor')

      call free_polygon(poly)
   end subroutine test_extensive_vs_intensive

   !---~---
   !   CLOSURE (reduce_cohort_to_pft, reduce_cohort_to_size). sum over PFTs and sum over size
   ! classes must both reproduce the site total EXACTLY -- that identity is the entire reason
   ! those axes are trustworthy, and it is the one property a reader will assume without
   ! checking.
   !---~---
   subroutine test_pft_and_size_closure()
      type(polygon_t) :: poly
      real(wp) :: site_val, pft_out(8), size_out(8), edges(8)
      logical  :: ok, vld(8)
      integer(ik) :: n
      call build_fixture(poly)
      edges(1:5) = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp, 100.0_wp]     ! 4 classes

      call reduce_cohort_to_polygon(poly, poly%cohort%agb, W_NPLANT, .false., site_val, ok)
      call reduce_cohort_to_pft (poly, poly%cohort%agb, W_NPLANT, .false., 3_ik, pft_out, vld, n)
      call check(n == 3_ik, 'pft axis length = n_pft')
      call check_close(sum(pft_out(1:3)), site_val, 1.0e-12_wp, 'sum over PFTs == mean_agb_py')

      call reduce_cohort_to_size(poly, poly%cohort%agb, W_NPLANT, .false., edges, 4_ik,          &
                                 size_out, vld, n)
      call check(n == 4_ik, 'size axis length = n_class')
      call check_close(sum(size_out(1:4)), site_val, 1.0e-12_wp, 'sum over DBH classes == mean_agb_py')

      !--- The stem count must be conserved by the binning too: no cohort split, none dropped.
      call reduce_cohort_to_polygon(poly, poly%cohort%nplant, W_NONE, .false., site_val, ok)
      call reduce_cohort_to_size(poly, poly%cohort%nplant, W_NONE, .false., edges, 4_ik,         &
                                 size_out, vld, n)
      call check_close(sum(size_out(1:4)), site_val, 1.0e-12_wp, 'sum over classes == mean_nplant_py')
      !--- Each cohort landed in exactly ONE bin (dbh 5, 15, 45 -> classes 1, 2, 4).
      call check(size_out(3) == 0.0_wp, 'the empty class-3 bin holds nothing')

      !--- LAI closes on both axes as well (a second field, so the test is not agb-specific).
      call reduce_cohort_to_polygon(poly, poly%cohort%leaf_area, W_NPLANT, .false., site_val, ok)
      call reduce_cohort_to_pft (poly, poly%cohort%leaf_area, W_NPLANT, .false., 3_ik, pft_out, vld, n)
      call check_close(sum(pft_out(1:3)), site_val, 1.0e-12_wp, 'sum over PFTs == mean_lai_py')
      call check_close(total_lai(poly), site_val, 1.0e-12_wp, 'total_lai wrapper agrees')

      call free_polygon(poly)
   end subroutine test_pft_and_size_closure

   !---~---
   !   EMPTY SETS. A patch with no cohorts, and a MEAN over a PFT with no members, must report
   ! invalid (-> _FillValue), never 0/0 and never a bare 0 that a reader would take for data.
   ! But an empty-PFT SUM is a true 0, because reporting fill there would break the closure
   ! asserted above the moment a PFT went locally extinct. The two conventions are different
   ! on purpose and both are pinned here.
   !---~---
   subroutine test_empty_sets()
      type(polygon_t) :: poly
      real(wp) :: out(8), v
      logical  :: vld(8), ok
      integer(ik) :: n
      call build_fixture(poly)
      !--- Make patch 2 empty (its cohort moves to patch 1).
      poly%patch%cohort_count(1) = 3_ik
      poly%patch%cohort_count(2) = 0_ik
      poly%patch%cohort_offset(2) = 4_ik
      poly%cohort%owner_patch(3) = 1_ik

      call reduce_cohort_to_patch(poly, poly%cohort%agb, W_NPLANT, .false., out, vld, n)
      call check(n == 2_ik,      'empty patch still occupies a slot')
      call check(vld(1),         'occupied patch is valid')
      call check(.not. vld(2),   'EMPTY patch is invalid -> _FillValue, not 0')

      !--- A MEAN over a PFT with no members is invalid; a SUM over it is a true 0.
      poly%n_pft = 4_ik                                   ! PFT 4 exists but has no cohorts
      call reduce_cohort_to_pft(poly, poly%cohort%dbh, W_NPLANT, .true., 4_ik, out, vld, n)
      call check(.not. vld(4), 'MEAN over an empty PFT is invalid')
      call reduce_cohort_to_pft(poly, poly%cohort%agb, W_NPLANT, .false., 4_ik, out, vld, n)
      call check(vld(4),            'SUM over an empty PFT is VALID (closure would break otherwise)')
      call check(out(4) == 0.0_wp,  'SUM over an empty PFT is exactly 0')

      !--- A site with no cohorts at all: the extensive sum is 0, the mean is invalid.
      poly%cohort%n = 0_ik
      poly%patch%cohort_count(1:2) = 0_ik
      call reduce_cohort_to_polygon(poly, poly%cohort%agb, W_NPLANT, .false., v, ok)
      call check(v == 0.0_wp, 'bare ground: extensive site total is 0')
      call reduce_cohort_to_polygon(poly, poly%cohort%dbh, W_BASAL_AREA, .true., v, ok)
      call check(.not. ok,    'bare ground: weighted MEAN is invalid, not 0/0')

      call free_polygon(poly)
   end subroutine test_empty_sets

   !---~---
   !   A cohort with NO SAMPLES (missing x: recruited after the fast loop) is left out: of a
   ! mean's numerator AND denominator, so the mean is that of the sampled cohorts, and of a
   ! sum, which it would only have added 0 to. A patch whose cohorts are all unsampled has
   ! nothing to report. Counting it as 0 would pull a canopy temperature towards 0 K
   ! (docs/science/diagnostics.md §2).
   !---~---
   subroutine test_unsampled_cohorts()
      type(polygon_t) :: poly
      real(wp) :: x(3), out(8), v
      logical  :: vld(8), ok
      integer(ik) :: n
      call build_fixture(poly)
      x = [300.0_wp, 310.0_wp, MISSING_VALUE]          ! cohort 3, alone in patch 2, has no samples
      !--- Leaf-area weights nplant*leaf_area = 2, 3 (patch 1, area 0.75): (600 + 930) / 5.
      call reduce_cohort_to_polygon(poly, x, W_LEAF_AREA, .true., v, ok)
      call check(ok, 'mean over the sampled cohorts is valid')
      call check_close(v, 306.0_wp, 1.0e-12_wp, 'an unsampled cohort is out of the mean, not a 0')
      call reduce_cohort_to_patch(poly, x, W_LEAF_AREA, .true., out, vld, n)
      call check(vld(1) .and. .not. vld(2), 'a patch of only unsampled cohorts reads fill')
      x = [1.0_wp, 2.0_wp, MISSING_VALUE]
      call reduce_cohort_to_polygon(poly, x, W_NONE, .false., v, ok)
      call check_close(v, 2.25_wp, 1.0e-12_wp, 'an unsampled cohort adds nothing to a sum')
      x = MISSING_VALUE
      call reduce_cohort_to_polygon(poly, x, W_LEAF_AREA, .true., v, ok)
      call check(.not. ok, 'a mean over only unsampled cohorts is invalid')
      call free_polygon(poly)
   end subroutine test_unsampled_cohorts

   !---~---
   !   SIZE-CLASS EDGES (dbh_class_index). Half-open [lo, hi) except the last class, which is
   ! CLOSED at the top so the largest tree in the stand is never dropped; below the first edge
   ! clamps into bin 1 for the same reason. Both are closure-preserving choices, so both are
   ! pinned.
   !---~---
   subroutine test_class_edges()
      real(wp) :: edges(5)
      edges = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp, 100.0_wp]
      call check(dbh_class_index( 0.0_wp, edges, 4_ik) == 1_ik, 'exactly the first edge -> bin 1')
      call check(dbh_class_index( 9.9_wp, edges, 4_ik) == 1_ik, 'just below an edge -> lower bin')
      call check(dbh_class_index(10.0_wp, edges, 4_ik) == 2_ik, 'exactly on an edge -> UPPER bin')
      call check(dbh_class_index(29.9_wp, edges, 4_ik) == 3_ik, 'interior bin')
      call check(dbh_class_index(30.0_wp, edges, 4_ik) == 4_ik, 'exactly the last lower edge -> last bin')
      call check(dbh_class_index(1.0e4_wp, edges, 4_ik) == 4_ik, 'above the top edge -> last bin (never dropped)')
      call check(dbh_class_index(-1.0_wp, edges, 4_ik) == 1_ik, 'below the bottom edge -> bin 1 (never dropped)')
   end subroutine test_class_edges

   !---~---
   !   SOIL COLUMNS (reduce_patch_columns_to_sites, gather_patch_columns). The area-weighted
   ! site column, and the FLATTENED 2-D (layer, patch) slab whose stride is the COMPILE-TIME
   ! layer ceiling. The closure below -- area-weighted mean over the patch axis of the 2-D
   ! slab == the site column -- is the acceptance test, and it is what proves the flattening
   ! index is right. Varying the ACTIVE layer count while the stride stays fixed is the
   ! specific mis-mapping it has to rule out. Two sites on columns of different depths each
   ! keep their own layers, at their own stride.
   !---~---
   subroutine test_soil_columns()
      type(polygon_t) :: poly, two
      real(wp) :: col(n_soil_layer_max, 2), sitecol(2*n_soil_layer_max)
      real(wp) :: flat(2*n_soil_layer_max), wmean
      logical  :: vld(2*n_soil_layer_max), svld(2*n_soil_layer_max)
      integer(ik) :: n, k, nlayer
      call build_fixture(poly)
      nlayer = 5_ik                                     ! ACTIVE layers < the ceiling, on purpose
      poly%site%surface(1)%soil_profile%n_active = nlayer
      do k = 1_ik, nlayer
         col(k, 1) = 0.10_wp * real(k, wp)              ! patch 1 profile
         col(k, 2) = 0.50_wp * real(k, wp)              ! patch 2 profile (deliberately different)
      end do

      call reduce_patch_columns_to_sites(poly, col, sitecol, svld, n)
      call check(n == n_soil_layer_max, 'one site column = the layer ceiling')
      call check(all(svld(1:nlayer)) .and. .not. any(svld(nlayer + 1_ik:n)), 'only the active layers are valid')
      !--- area-weighted: 0.75*0.1k + 0.25*0.5k = 0.2k
      do k = 1_ik, nlayer
         call check_close(sitecol(k), 0.2_wp*real(k,wp), 1.0e-12_wp, 'area-weighted site soil column')
      end do

      call gather_patch_columns(poly, col, flat, vld, n)
      call check(n == 2_ik*n_soil_layer_max, 'flat slab length = n_patch * n_soil_layer_max')
      !---~---
      !   The STRIDE is the compile-time ceiling, NOT nlayer: patch 2 layer 1 sits at
      ! n_soil_layer_max + 1, not nlayer + 1. Striding by the live count is the mis-mapping
      ! that would produce plausible shifted profiles instead of an error.
      !---~---
      call check_close(flat(1),                       col(1, 1), 1.0e-12_wp, 'patch 1 layer 1')
      call check_close(flat(nlayer),                  col(nlayer, 1), 1.0e-12_wp, 'patch 1 last active layer')
      call check_close(flat(n_soil_layer_max + 1_ik), col(1, 2), 1.0e-12_wp, 'patch 2 starts at the CEILING stride')
      call check(.not. vld(nlayer + 1_ik), 'inactive layers of patch 1 are invalid -> _FillValue')
      call check(vld(n_soil_layer_max + 1_ik), 'patch 2 active layers are valid')

      !--- The acceptance identity, computed the way a reader would.
      do k = 1_ik, nlayer
         wmean = poly%patch%area(1) * flat(k)                                                    &
               + poly%patch%area(2) * flat(n_soil_layer_max + k)
         call check_close(wmean, sitecol(k), 1.0e-12_wp, 'area-weighted 2-D slab == site column')
      end do

      !---~---
      !   Two sites, a patch each: site 1 on 5 layers, site 2 on 3. Each site's column is its
      ! patch's, at the site's own stride, whatever the sites' areas.
      !---~---
      call alloc_polygon(two, 3_ik, 4_ik, 4_ik, 4_ik, test_pft_table(), n_site=2_ik)
      two%site%area(1:2) = [0.75_wp, 0.25_wp]
      two%site%surface(1)%soil_profile%n_active = 5_ik
      two%site%surface(2)%soil_profile%n_active = 3_ik
      two%patch%n = 2_ik
      two%patch%area(1:2) = 1.0_wp
      two%patch%owner_site(1:2) = [1_ik, 2_ik]
      call rebuild_patch_map(two)
      call reduce_patch_columns_to_sites(two, col, sitecol, svld, n)
      call check(n == 2_ik*n_soil_layer_max, 'two site columns, each the layer ceiling long')
      call check_close(sitecol(3), col(3, 1), 1.0e-12_wp, 'site 1 is its patch''s column')
      call check_close(sitecol(n_soil_layer_max + 3_ik), col(3, 2), 1.0e-12_wp,                    &
                       'site 2 is its patch''s column, at the ceiling stride')
      call check(svld(5) .and. .not. svld(n_soil_layer_max + 4_ik), 'each site has its own layers')
      call gather_patch_columns(two, col, flat, vld, n)
      call check(vld(5) .and. vld(n_soil_layer_max + 3_ik) .and. .not. vld(n_soil_layer_max + 4_ik),    &
                 'each patch has its site''s layers')

      call free_polygon(two)
      call free_polygon(poly)
   end subroutine test_soil_columns

   !--- 3. THE EXTRACT (ed_output_extract).

   !--- A fixed 3-cohort / 1-patch site (slot set constant across a window), agb 1, 2, 3.
   subroutine extract_fixture(poly)
      type(polygon_t), intent(out) :: poly
      call alloc_polygon(poly, 2_ik, 64_ik, 8_ik, 4_ik, test_pft_table())
      poly%cohort%n = 3_ik
      poly%cohort%nplant(1:3)      = [1.0_wp, 1.0_wp, 1.0_wp]
      call set_cohort_agb(poly, [1.0_wp, 2.0_wp, 3.0_wp], 3_ik)   ! opaque (see header note)
      poly%cohort%pft(1:3)         = [1_ik, 1_ik, 2_ik]
      poly%cohort%owner_patch(1:3) = [1_ik, 1_ik, 1_ik]
      poly%patch%n = 1_ik
      poly%patch%area(1)          = 1.0_wp
      poly%patch%cohort_offset(1) = 1_ik
      poly%patch%cohort_count(1)  = 3_ik
   end subroutine extract_fixture

   !---~---
   !   extract_variable: the live cohort slab and its count; the site twin of the same field,
   ! reduced by its weight; and a cohort-block mean, where a cohort with no samples is fill.
   !---~---
   subroutine test_extract()
      type(polygon_t)    :: poly
      type(var_desc_t)   :: v_agb_c, v_agb_s, v_lt
      real(wp)           :: slab(64), scal
      logical            :: vslab(64)
      type(diag_setup_t) :: dp
      integer(ik)        :: n_out
      call extract_fixture(poly)

      !--- extract_variable returns the live cohort slab + count.
      v_agb_c%name = 'mean_agb_co'
      v_agb_c%dim = DIM_COHORT
      v_agb_c%agg = AGG_TMEAN
      v_agb_c%source_id = FLD_C_AGB
      call extract_variable(poly, dp, v_agb_c, scal, slab, vslab, n_out)
      call check(n_out == 3_ik, 'extract cohort n')
      call check_close(slab(2), 2.0_wp, 1.0e-12_wp, 'extract cohort agb(2)')

      !---~---
      !   Site scalar reduction via extract (total_agb = sum area*nplant*agb = 1*(1+2+3)). The
      ! site twin is the SAME field with a different dim + weight: agb is EXTENSIVE, so it
      ! reduces as a weighted SUM over nplant (W_NPLANT), landing in kgC/m2 ground.
      !---~---
      v_agb_s%name = 'mean_agb_py'
      v_agb_s%dim = DIM_SCALAR
      v_agb_s%agg = AGG_MEAN
      v_agb_s%source_id = FLD_C_AGB
      v_agb_s%weight = W_NPLANT
      v_agb_s%mean = .false.
      call extract_variable(poly, dp, v_agb_s, scal, slab, vslab, n_out)
      call check(n_out == 0_ik, 'scalar extract n_out=0')
      call check_close(scal, 6.0_wp, 1.0e-12_wp, 'extract site agb')

      !---~---
      !   Leaf temperature is the cohort block's within-step mean. A cohort with no
      ! samples -- recruited after the fast loop -- is fill on the cohort axis and out of the
      ! leaf-area-weighted site mean; a mean over nothing is missing, not 0 K.
      !---~---
      call alloc_cohort_diag(poly%cohort%diag, 64_ik, .true.)
      poly%cohort%diag%n = 3_ik
      poly%cohort%leaf_area(1:3) = 1.0_wp
      poly%cohort%diag%v(CO_LEAF_TEMP, 1:2) = [300.0_wp, 310.0_wp] * 900.0_wp
      poly%cohort%diag%w(1:2) = 900.0_wp
      poly%cohort%diag%w(3) = 0.0_wp
      v_lt%name = 'mean_leaf_temp_co'
      v_lt%dim = DIM_COHORT
      v_lt%agg = AGG_TMEAN
      v_lt%source_id = FLD_C_DIAG0 + CO_LEAF_TEMP
      call extract_variable(poly, dp, v_lt, scal, slab, vslab, n_out)
      call check(n_out == 3_ik .and. vslab(1) .and. vslab(2) .and. .not. vslab(3),               &
                 'a cohort with no samples is fill on the cohort axis')
      v_lt%name = 'mean_leaf_temp_py'
      v_lt%dim = DIM_SCALAR
      v_lt%weight = W_LEAF_AREA
      v_lt%mean = .true.
      call extract_variable(poly, dp, v_lt, scal, slab, vslab, n_out)
      call check_close(scal, 305.0_wp, 1.0e-12_wp, 'the site leaf temperature leaves it out')
      poly%cohort%diag%w(1:2) = 0.0_wp
      call extract_variable(poly, dp, v_lt, scal, slab, vslab, n_out)
      call check(scal == MISSING_VALUE, 'a mean over no sampled cohort is missing, not 0 K')

      call free_polygon(poly)
   end subroutine test_extract

   !---~---
   !   extract_fast_scalar: a patch-block id reads its row of the sub-step's site mean, and a
   ! forcing id its row of the sub-step's forcing table.
   !---~---
   subroutine test_fast_extract()
      real(wp) :: f0(N_PY_ROWS), s(N_PA_ROWS)
      f0 = 0.0_wp
      f0(PY_CO2) = 415.0_wp
      s = 0.0_wp
      s(PA_CAS_TEMP) = 290.0_wp
      s(PA_LE)       = 100.0_wp
      s(PA_H)        = 50.0_wp
      s(PA_GPP)      = 12.0_wp
      s(PA_SW_UP)    = 83.0_wp
      s(PA_LW_UP)    = 455.0_wp
      call check_close(extract_fast_scalar(SRC_F_PA0 + PA_CAS_TEMP, s, f0), 290.0_wp, 1.0e-12_wp, 'fast extract cas_temp')
      call check_close(extract_fast_scalar(SRC_F_PA0 + PA_LE,       s, f0), 100.0_wp, 1.0e-12_wp, 'fast extract le')
      call check_close(extract_fast_scalar(SRC_F_PA0 + PA_H,        s, f0),  50.0_wp, 1.0e-12_wp, 'fast extract h')
      call check_close(extract_fast_scalar(SRC_F_PA0 + PA_GPP,      s, f0),  12.0_wp, 1.0e-12_wp, 'fast extract gpp_rate')
      call check_close(extract_fast_scalar(SRC_F_PA0 + PA_SW_UP,    s, f0),  83.0_wp, 1.0e-12_wp, 'fast extract sw_up')
      call check_close(extract_fast_scalar(SRC_F_PA0 + PA_LW_UP,    s, f0), 455.0_wp, 1.0e-12_wp, 'fast extract lw_up')
      call check_close(extract_fast_scalar(SRC_F_PY0 + PY_CO2, s, f0), 415.0_wp, 1.0e-12_wp, 'fast extract forcing CO2')
   end subroutine test_fast_extract

   !---~---
   !   4. THE INTEGRATION (ed_output_integrate, with ed_output_setup's sizing;
   ! docs/science/diagnostics.md §4).
   !---~---

   !--- Build a scalar buffer for one operator, fold a sequence, normalize.
   subroutine run_scalar(agg, x, dt, out, valid)
      integer(ik), intent(in)  :: agg
      real(wp),    intent(in)  :: x(:), dt(:)
      real(wp),    intent(out) :: out
      logical,     intent(out) :: valid
      type(var_desc_t)     :: v
      type(integ_buffer_t) :: buf
      integer(ik) :: i
      v%dim = DIM_SCALAR
      v%agg = agg
      call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      do i = 1_ik, int(size(x), ik)
         call integrate_scalar(buf, x(i), dt(i))
      end do
      call normalize_scalar(buf, out, valid)
   end subroutine run_scalar

   !--- integrate_scalar + normalize_scalar, one operator at a time.
   subroutine test_scalar_operators()
      real(wp) :: out
      logical  :: valid
      !--- MEAN (equal weight).
      call run_scalar(AGG_MEAN, [2.0_wp,4.0_wp,6.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check(valid, 'MEAN valid')
      call check_close(out, 4.0_wp, 1.0e-12_wp, 'AGG_MEAN')
      !--- TMEAN on NON-uniform dt (17.5, not the plain mean 15).
      call run_scalar(AGG_TMEAN, [10.0_wp,20.0_wp], [1.0_wp,3.0_wp], out, valid)
      call check_close(out, 17.5_wp, 1.0e-12_wp, 'AGG_TMEAN dt-weighted')
      !--- RATE: amounts 1, 2 and 3 over steps of 1, 1 and 2 s -> 6 per 4 s, 1.5 per s.
      call run_scalar(AGG_RATE, [1.0_wp,2.0_wp,3.0_wp], [1.0_wp,1.0_wp,2.0_wp], out, valid)
      call check_close(out, 1.5_wp, 1.0e-12_wp, 'AGG_RATE: the amounts over the time they took')
      !--- LAST.
      call run_scalar(AGG_LAST, [7.0_wp,8.0_wp,9.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check_close(out, 9.0_wp, 1.0e-12_wp, 'AGG_LAST')
   end subroutine test_scalar_operators

   !--- normalize_scalar on a buffer that saw nothing.
   subroutine test_zero_sample_guard()
      type(var_desc_t)     :: v
      type(integ_buffer_t) :: buf
      real(wp) :: out
      logical  :: valid
      !--- Never-touched MEAN buffer -> invalid, MISSING, no NaN / no huge leak.
      v%dim = DIM_SCALAR
      v%agg = AGG_MEAN
      call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      call normalize_scalar(buf, out, valid)
      call check(.not. valid, 'zero-sample -> invalid')
      call check(out == MISSING_VALUE, 'zero-sample -> MISSING')
      !--- A RATE buffer that saw no step reports nothing, not a 0/0.
      v%agg = AGG_RATE
      call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      call reset_buffer(buf)
      call normalize_scalar(buf, out, valid)
      call check(.not. valid .and. out == MISSING_VALUE, 'a RATE with no step -> MISSING')
   end subroutine test_zero_sample_guard

   !---~---
   !   integrate_slab + normalize_slab over a fixed-slot window: TMEAN per slot, dt-weighted,
   ! of the slabs extract_variable reads. step1 agb=[1,2,3] dt=1 ; step2 agb=[3,4,5] dt=3 ->
   ! slot1 = (1+9)/4 = 2.5. The agb mutation goes through the OPAQUE set_cohort_agb
   ! (production mutates via veg dynamics).
   !---~---
   subroutine test_slab_window()
      type(polygon_t)      :: poly
      type(var_desc_t)     :: v_agb_c
      type(integ_buffer_t) :: buf
      real(wp)             :: slab(64), out(64), scal
      logical              :: valid(64), vslab(64)
      type(diag_setup_t)   :: dp
      integer(ik)          :: n_out
      call extract_fixture(poly)
      v_agb_c%name = 'mean_agb_co'
      v_agb_c%dim = DIM_COHORT
      v_agb_c%agg = AGG_TMEAN
      v_agb_c%source_id = FLD_C_AGB
      call alloc_integ_buffer(buf, v_agb_c, FREQ_MONTHLY, 64_ik)
      call set_cohort_agb(poly, [1.0_wp, 2.0_wp, 3.0_wp], 3_ik)
      call extract_variable(poly, dp, v_agb_c, scal, slab, vslab, n_out)
      call integrate_slab(buf, slab, n_out, 1.0_wp)
      call set_cohort_agb(poly, [3.0_wp, 4.0_wp, 5.0_wp], 3_ik)
      call extract_variable(poly, dp, v_agb_c, scal, slab, vslab, n_out)
      call integrate_slab(buf, slab, n_out, 3.0_wp)
      call normalize_slab(buf, out, valid, n_out)
      call check(n_out == 3_ik, 'slab n_out')
      call check(all(valid(1:3)), 'slab slots valid')
      call check_close(out(1), 2.5_wp, 1.0e-12_wp, 'slab TMEAN slot 1')
      call check_close(out(2), 3.5_wp, 1.0e-12_wp, 'slab TMEAN slot 2')
      call check_close(out(3), 4.5_wp, 1.0e-12_wp, 'slab TMEAN slot 3')
      !--- Slots past the live count are invalid (fill tail).
      call check(.not. valid(4), 'fill-tail slot invalid')
      call free_polygon(poly)
   end subroutine test_slab_window

   !---~---
   !   The FAST (sub-daily) tier: a 2-sub-step window folded by output_integrate_fast +
   ! close_tier(1) -> TMEAN in pending(1), across the scalar, soil-column, and per-cohort
   ! slab paths. netCDF-free; mirrors what main's replay loop does.
   !---~---
   subroutine test_fast_tier()
      type(ed_config_t)      :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      integer(ik) :: k_cas, k_soil, k_air, k_hgt, k_cas_p, k_soil_p, nl, c_soil, c_hgt, c_cas_p, c_soil_p
      real(wp), parameter :: DT = 900.0_wp    ! uniform sub-step -> TMEAN == plain mean

      !---~---
      !   Build a manager with the FAST tier on (the defaults' F variables) and the soil by patch,
      ! then fold a 2-sub-step window and close it.
      !---~---
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.               ! FAST tier on
      cfg%output%cohort_max = 8_ik
      call alloc_output_set(files, bufs, cfg, [selection_t('mean_soil_temp_layer_pa', OVR_MASK, FREQ_FAST)])
      call check(files%reg%nidx(1) > 0_ik, 'FAST tier has live variables')

      !--- Stage 2 sub-steps of two patches with known values (as fast_dynamics would).
      nl = n_soil_layer_max
      call alloc_fast_staging(bufs, 2_ik, 2_ik, 8_ik, 2_ik)
      bufs%fast_n_cohort = 2_ik
      bufs%fast_site(PA_CAS_TEMP, :) = [290.0_wp, 294.0_wp]                 ! mean 292
      bufs%fast_site(PA_LE, :)       = [100.0_wp, 200.0_wp]
      bufs%fast_patch(PA_CAS_TEMP, 1, :) = [290.0_wp, 294.0_wp]             ! patch 1 mean 292
      bufs%fast_patch(PA_CAS_TEMP, 2, :) = [300.0_wp, 304.0_wp]             ! patch 2 mean 302
      bufs%fast_forcing(PY_TAIR, 1) = 300.0_wp                              ! mean 301
      bufs%fast_forcing(PY_TAIR, 2) = 302.0_wp
      bufs%fast_soil_temp(1:2,1) = [280.0_wp, 281.0_wp]                     ! slot1 mean 281
      bufs%fast_soil_temp(1:2,2) = [282.0_wp, 283.0_wp]                     ! slot2 mean 282
      bufs%fast_soil_temp_patch(1:2,1)       = [280.0_wp, 281.0_wp]         ! patch 1: 281, 282
      bufs%fast_soil_temp_patch(1:2,2)       = [282.0_wp, 283.0_wp]
      bufs%fast_soil_temp_patch(nl+1:nl+2,1) = [290.0_wp, 291.0_wp]         ! patch 2: 291, 292
      bufs%fast_soil_temp_patch(nl+1:nl+2,2) = [292.0_wp, 293.0_wp]
      bufs%fast_cohort_height(1:2,1) = [10.0_wp, 12.0_wp]                   ! slot1=10, slot2=12 (constant)
      bufs%fast_cohort_height(1:2,2) = [10.0_wp, 12.0_wp]

      call output_integrate_fast(files, bufs, 1_ik, DT)
      call output_integrate_fast(files, bufs, 2_ik, DT)
      call close_tier(files, bufs, 1_ik)

      !---~---
      !   Scalar path: mean_cas_temp_py = TMEAN(290,294) = 292. On the F tier a variable reads its
      ! fast source, the staged sub-step samples, where the coarse tiers read live state.
      !---~---
      k_cas = find_var_index(files%reg, 'mean_cas_temp_py')
      call check(k_cas > 0_ik, 'mean_cas_temp_py registered')
      call check(bufs%pending(1)%svalid(k_cas), 'FAST cas_temp valid after close')
      call check_close(bufs%pending(1)%sval(k_cas), 292.0_wp, 1.0e-10_wp, 'FAST cas_temp TMEAN')
      !--- Soil-column slab path: mean_soil_temp_layer_si slot means.
      k_soil = find_var_index(files%reg, 'mean_soil_temp_layer_si')
      c_soil = slab_col(bufs%pending(1), k_soil)
      call check(bufs%pending(1)%nslab(k_soil) == nl .and. bufs%pending(1)%slabvalid(2,c_soil)          &
                 .and. .not. bufs%pending(1)%slabvalid(3,c_soil), 'FAST soil slab: the site''s 2 layers of nl')
      call check_close(bufs%pending(1)%slab(1,c_soil), 281.0_wp, 1.0e-10_wp, 'FAST soil_temp slot 1')
      call check_close(bufs%pending(1)%slab(2,c_soil), 282.0_wp, 1.0e-10_wp, 'FAST soil_temp slot 2')
      !--- Per-cohort slab path: n_cohort (the slots are checked on mean_height_co below).
      call check(bufs%pending(1)%n_cohort == 2_ik, 'FAST n_cohort = 2')
      !--- Forcing air-temp scalar path + per-cohort HEIGHT slab path (tallest-cohort post-proc).
      k_air = find_var_index(files%reg, 'mean_air_temp_py')
      call check(k_air > 0_ik, 'mean_air_temp_py registered')
      call check_close(bufs%pending(1)%sval(k_air), 301.0_wp, 1.0e-10_wp, 'FAST air_temp TMEAN')
      k_hgt = find_var_index(files%reg, 'mean_height_co')
      c_hgt = slab_col(bufs%pending(1), k_hgt)
      call check(k_hgt > 0_ik, 'mean_height_co registered')
      call check_close(bufs%pending(1)%slab(1,c_hgt), 10.0_wp, 1.0e-10_wp, 'FAST height cohort 1')
      call check_close(bufs%pending(1)%slab(2,c_hgt), 12.0_wp, 1.0e-10_wp, 'FAST height cohort 2')
      !---~---
      !   The patch axis: each patch's own mean, and each patch's soil column in the
      ! soil-by-patch layout, with the inactive layers as fill.
      !---~---
      k_cas_p = find_var_index(files%reg, 'mean_cas_temp_pa')
      c_cas_p = slab_col(bufs%pending(1), k_cas_p)
      call check(k_cas_p > 0_ik, 'mean_cas_temp_pa registered')
      call check(bufs%pending(1)%n_patch == 2_ik, 'FAST n_patch = 2')
      call check_close(bufs%pending(1)%slab(1,c_cas_p), 292.0_wp, 1.0e-10_wp, 'FAST cas_temp patch 1')
      call check_close(bufs%pending(1)%slab(2,c_cas_p), 302.0_wp, 1.0e-10_wp, 'FAST cas_temp patch 2')
      k_soil_p = find_var_index(files%reg, 'mean_soil_temp_layer_pa')
      c_soil_p = slab_col(bufs%pending(1), k_soil_p)
      call check(k_soil_p > 0_ik, 'mean_soil_temp_layer_pa registered')
      call check_close(bufs%pending(1)%slab(1,c_soil_p),    281.0_wp, 1.0e-10_wp, 'FAST soil patch 1 layer 1')
      call check_close(bufs%pending(1)%slab(2,c_soil_p),    282.0_wp, 1.0e-10_wp, 'FAST soil patch 1 layer 2')
      call check_close(bufs%pending(1)%slab(nl+1,c_soil_p), 291.0_wp, 1.0e-10_wp, 'FAST soil patch 2 layer 1')
      call check_close(bufs%pending(1)%slab(nl+2,c_soil_p), 292.0_wp, 1.0e-10_wp, 'FAST soil patch 2 layer 2')
      call check(.not. bufs%pending(1)%slabvalid(3,c_soil_p), 'FAST soil patch 1 layer 3 is fill (inactive)')
   end subroutine test_fast_tier

   !---~---
   !   The driver's order: setup_output_files, then an [output].io_config override that switches on a
   ! slab variable setup left off, then finalize_output_files and alloc_output_buffers. The slab must
   ! hold every soil layer, read from the QUEUED record through slab_col, as the writer reads it.
   !---~---
   subroutine test_slab_sized_after_overrides()
      type(ed_config_t)      :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      integer(ik) :: k, c, i, nl, n_live_slab
      logical     :: found, every_layer
      nl = n_soil_layer_max
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.               ! FAST tier on
      call setup_output_files(files, cfg)
      call select_nothing(files)                   ! no slab variable is live yet
      n_live_slab = 0_ik
      do k = 1_ik, files%reg%nvar
         if (files%reg%var(k)%dim /= DIM_SCALAR .and. files%reg%var(k)%streams /= FREQ_NONE)        &
            n_live_slab = n_live_slab + 1_ik
      end do
      call check(n_live_slab == 0_ik, 'override test premise: setup leaves no slab variable live')

      call apply_variable_override(files%reg, 'mean_soil_temp_layer_si', OVR_MASK, FREQ_FAST, found)
      call check(found, 'mean_soil_temp_layer_si is a registry variable')
      call resolve_streams(files%reg)              ! what apply_io_overrides does after its lines
      call finalize_output_files(files)
      call alloc_output_buffers(files, bufs)
      call check(files%max_slab >= nl, 'max_slab covers the soil axis an override switched on')

      call alloc_fast_staging(bufs, 1_ik, 0_ik, 1_ik, nl)
      bufs%fast_n_cohort = 0_ik
      bufs%fast_soil_temp(:,1) = [(270.0_wp + real(i, wp), i = 1_ik, nl)]
      call output_integrate_fast(files, bufs, 1_ik, 900.0_wp)
      call close_tier(files, bufs, 1_ik)

      k = find_var_index(files%reg, 'mean_soil_temp_layer_si')
      associate (r => bufs%queue(1)%rec(bufs%queue(1)%n))
         c = slab_col(r, k)
         call check(r%nslab(k) == nl, 'the queued soil slab has every layer')
         call check(size(r%slab, 1) >= nl, 'the queued record has a row for every layer')
         every_layer = size(r%slab, 1) >= nl
         if (every_layer) every_layer = all(r%slab(1:nl, c) == [(270.0_wp + real(i, wp), i = 1_ik, nl)])
         call check(every_layer, 'every soil layer reads back through slab_col, as the writer reads it')
      end associate
   end subroutine test_slab_sized_after_overrides

   !---~---
   !   A tier whose only patch output is soil by patch still records its patch count: the slab
   ! is patch-major, so it holds nslab / n_soil_layer_max patches. If only a DIM_PATCH
   ! variable set the count, such a record would write one patch and every value as fill.
   !---~---
   subroutine test_soil_patch_record()
      type(ed_config_t)      :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      integer(ik) :: k, nl
      logical     :: found
      nl = n_soil_layer_max
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.               ! FAST tier on
      call setup_output_files(files, cfg)
      call select_nothing(files)                   ! no slab variable is live yet
      call apply_variable_override(files%reg, 'mean_soil_temp_layer_pa', OVR_MASK, FREQ_FAST, found)
      call check(found, 'mean_soil_temp_layer_pa is a registry variable')
      call resolve_streams(files%reg)
      call finalize_output_files(files)
      call alloc_output_buffers(files, bufs)
      call alloc_fast_staging(bufs, 1_ik, 3_ik, 1_ik, 2_ik)
      bufs%fast_soil_temp_patch(:,1) = 285.0_wp
      call output_integrate_fast(files, bufs, 1_ik, 900.0_wp)
      call close_tier(files, bufs, 1_ik)
      k = find_var_index(files%reg, 'mean_soil_temp_layer_pa')
      associate (r => bufs%queue(1)%rec(bufs%queue(1)%n))
         call check(r%n_patch == 3_ik, 'a soil-by-patch-only record counts its patches')
         call check(r%nslab(k) == 3_ik * nl, 'the soil-by-patch slab holds every patch''s column')
      end associate
   end subroutine test_soil_patch_record

   !---~---
   !   Two polygons' buffers for one file set: each reduces only its own samples, and closing
   ! one leaves the other's window open. Folds are interleaved the way a region's month
   ! loop would interleave its polygons.
   !---~---
   subroutine test_two_buffers()
      type(ed_config_t)      :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: a, b
      integer(ik) :: k_cas
      real(wp), parameter :: DT = 900.0_wp
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.
      cfg%output%cohort_max = 8_ik
      call alloc_output_set(files, a, cfg)
      call alloc_output_buffers(files, b)
      call check(a%fast_on .and. b%fast_on, 'both polygons'' buffers stage the FAST tier')
      call stage_cas(a, [290.0_wp, 294.0_wp])                              ! mean 292
      call stage_cas(b, [270.0_wp, 280.0_wp])                              ! mean 275
      call output_integrate_fast(files, a, 1_ik, DT)
      call output_integrate_fast(files, b, 1_ik, DT)
      call output_integrate_fast(files, a, 2_ik, DT)
      call close_tier(files, a, 1_ik)
      k_cas = find_var_index(files%reg, 'mean_cas_temp_py')
      call check_close(a%pending(1)%sval(k_cas), 292.0_wp, 1.0e-10_wp, 'polygon a: its own mean')
      call check(a%queue(1)%n == 1_ik .and. b%queue(1)%n == 0_ik, 'closing polygon a''s window queues only a''s record')
      call check(b%has_data(1) .and. .not. a%has_data(1), 'polygon b''s window stays open')
      call output_integrate_fast(files, b, 2_ik, DT)
      call close_tier(files, b, 1_ik)
      call check_close(b%pending(1)%sval(k_cas), 275.0_wp, 1.0e-10_wp, 'polygon b: its own mean')
      call check_close(a%queue(1)%rec(1)%sval(k_cas), 292.0_wp, 1.0e-10_wp, 'polygon a''s queued record is untouched')
   end subroutine test_two_buffers

   !---~---
   !   The slow tick across a month boundary (docs/science/diagnostics.md §4, Which steps a
   ! record holds). The step from 31 January to 1 February is January's: it folds into
   ! January, which then closes. The stand is restructured only after that tick
   ! (the next step's start), so January's cohort slab holds one slot set, and February opens on
   ! 1 February with the new one.
   !---~---
   subroutine test_boundary_step()
      type(ed_config_t)      :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      type(polygon_t)        :: poly
      type(ed_time_t)        :: jan30, jan31, feb1
      integer(ik) :: k_s, k_c, k_n, c_c
      real(wp), parameter :: DT = 86400.0_wp
      cfg = build_test_config(DT)
      cfg%output%enabled    = .true.
      cfg%output%freq_on    = [.false., .true., .true., .false.]    ! daily + monthly
      cfg%output%cohort_max = 8_ik
      call alloc_output_set(files, bufs, cfg, [selection_t('*', OVR_FALSE, FREQ_NONE),              &
                            selection_t('mean_agb_*', OVR_MASK, ior(FREQ_DAILY, FREQ_MONTHLY)),          &
                            selection_t('n_cohort_py', OVR_MASK, ior(FREQ_DAILY, FREQ_MONTHLY))])
      k_s = find_var_index(files%reg, 'mean_agb_py')
      k_c = find_var_index(files%reg, 'mean_agb_co')
      c_c = slab_col(bufs%pending(3), k_c)
      k_n = find_var_index(files%reg, 'n_cohort_py')
      call check(k_s > 0_ik .and. k_c > 0_ik .and. k_n > 0_ik, 'boundary fixture registers its variables')

      call alloc_polygon(poly, 2_ik, 8_ik, 8_ik, 4_ik, cfg%pft)
      poly%patch%n = 1_ik
      poly%patch%area(1) = 1.0_wp
      poly%patch%cohort_offset(1) = 1_ik
      jan30 = ed_time_t(year=2000_ik, month=1_ik, day=30_ik)
      jan31 = time_advance_days(jan30, 1_ik)
      feb1  = time_advance_days(jan31, 1_ik)

      !---~---
      !   30 and 31 January: one cohort, agb 10 then 12 at the ends of the two steps. The
      ! second step ends on the month boundary.
      !---~---
      call set_cohorts(poly, [10.0_wp])
      call output_integrate(files, bufs, poly, jan30, DT, .true., .false., .false.)
      call set_cohorts(poly, [12.0_wp])
      call output_integrate(files, bufs, poly, jan31, DT, .true., .true., .false.)

      call check(bufs%pending(2)%t_open%month == 1_ik .and. bufs%pending(2)%t_open%day == 31_ik,     &
                 'the step from 31 January to 1 February is the 31 January daily record')
      call check(bufs%queue(3)%n == 1_ik, 'the month turning closes January')
      call check(bufs%pending(3)%t_open%month == 1_ik .and. bufs%pending(3)%t_open%day == 30_ik,     &
                 'January opens at the start of its first step')
      call check_close(bufs%pending(3)%sval(k_s), 11.0_wp, 1.0e-12_wp,                             &
                       'January''s site mean includes its last step: (10 + 12) / 2')
      call check(bufs%pending(3)%nslab(k_c) == 1_ik, 'January''s cohort slab: its one slot set')
      call check_close(bufs%pending(3)%slab(1,c_c), 11.0_wp, 1.0e-12_wp,                           &
                       'January''s cohort slab includes its last step too')
      call check(.not. bufs%has_data(3), 'nothing of the boundary step is left for February')

      !---~---
      !   The boundary: a recruit joins (restructure_stand runs here in the driver). Then 1
      ! February's own step; close the month by hand.
      !---~---
      call set_cohorts(poly, [20.0_wp, 5.0_wp])
      call output_integrate(files, bufs, poly, feb1, DT, .true., .false., .false.)
      call check(bufs%t_open(3)%month == 2_ik .and. bufs%t_open(3)%day == 1_ik, 'February opens on 1 February')
      call close_tier(files, bufs, 3_ik)
      call check_close(bufs%pending(3)%sval(k_s), 25.0_wp, 1.0e-12_wp, 'February''s site mean: its own step')
      call check(bufs%pending(3)%nslab(k_c) == 2_ik, 'February''s cohort slab: the new slot set')
      call check_close(bufs%pending(3)%slab(2,c_c), 5.0_wp, 1.0e-12_wp, 'February''s recruit')
      call check_close(bufs%pending(3)%sval(k_n), 2.0_wp, 1.0e-12_wp, 'February''s cohort count')
      call free_polygon(poly)
   end subroutine test_boundary_step

   !---~---
   !   Snapshots on the tiers. Daily: the last staged sub-step of the step, for a polygon value
   ! (cas_temp_py), a patch slab (cas_temp_pa), an F-only state (cas_co2_pa) and a cohort row
   ! reduced as its mean is (leaf_temp_py, leaf-area-weighted); a step with nothing staged
   ! reads as fill. F: a snapshot keeps the last sub-step where its mean averages them.
   !---~---
   subroutine test_snapshots()
      type(ed_config_t)      :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      type(polygon_t)        :: poly
      type(ed_time_t)        :: day1
      integer(ik) :: k_t, k_tp, k_r, k_lt, k_m, c_tp
      logical     :: found
      real(wp), parameter :: DT = 86400.0_wp
      cfg = build_test_config(DT)
      cfg%output%enabled    = .true.
      cfg%output%freq_on    = [.true., .true., .false., .false.]     ! F + daily
      cfg%output%cohort_max = 8_ik
      call setup_output_files(files, cfg)
      call select_nothing(files)
      call apply_variable_override(files%reg, 'cas_temp_py', OVR_MASK, ior(FREQ_FAST, FREQ_DAILY), found)
      call apply_variable_override(files%reg, 'mean_cas_temp_py', OVR_MASK, FREQ_FAST, found)
      call apply_variable_override(files%reg, 'cas_temp_pa', OVR_MASK, FREQ_DAILY, found)
      call apply_variable_override(files%reg, 'cas_co2_pa', OVR_MASK, FREQ_DAILY, found)
      call apply_variable_override(files%reg, 'leaf_temp_py', OVR_MASK, FREQ_DAILY, found)
      call check(found, 'the snapshots are registry variables')
      call resolve_streams(files%reg)
      call finalize_output_files(files)
      call alloc_output_buffers(files, bufs)
      call check(bufs%fast_on, 'a snapshot read from the staging turns the staging on')
      k_t  = find_var_index(files%reg, 'cas_temp_py')
      k_m  = find_var_index(files%reg, 'mean_cas_temp_py')
      k_tp = find_var_index(files%reg, 'cas_temp_pa')
      k_r  = find_var_index(files%reg, 'cas_co2_pa')
      k_lt = find_var_index(files%reg, 'leaf_temp_py')

      !--- One patch, two cohorts of leaf area 1 and 3, the cohort block on; two staged sub-steps.
      call alloc_polygon(poly, 2_ik, 8_ik, 8_ik, 4_ik, cfg%pft)
      poly%patch%n = 1_ik
      poly%patch%area(1) = 1.0_wp
      poly%patch%cohort_offset(1) = 1_ik
      call set_cohorts(poly, [10.0_wp, 20.0_wp])
      poly%cohort%leaf_area(1:2) = [1.0_wp, 3.0_wp]
      call alloc_cohort_diag(poly%cohort%diag, 8_ik, .true.)
      poly%cohort%diag%n = 2_ik
      call stage_cas(bufs, [290.0_wp, 294.0_wp])
      bufs%fast_n_cohort = 2_ik
      bufs%fast_patch(PA_CAS_CO2, 1, :) = [400.0_wp, 410.0_wp]
      allocate(bufs%fast_cohort_diag(8, N_CO_FAST_ROWS, 2), source=0.0_wp)
      bufs%fast_cohort_diag(1:2, CO_LEAF_TEMP, 2) = [290.0_wp, 298.0_wp]
      bufs%fast_ready = .true.

      call output_integrate_fast(files, bufs, 1_ik, 900.0_wp)
      call output_integrate_fast(files, bufs, 2_ik, 900.0_wp)
      call close_tier(files, bufs, 1_ik)
      call check_close(bufs%pending(1)%sval(k_t), 294.0_wp, 1.0e-12_wp, 'F: the snapshot keeps the last sub-step')
      call check_close(bufs%pending(1)%sval(k_m), 292.0_wp, 1.0e-12_wp, 'F: the mean averages them')

      day1 = ed_time_t(year=2000_ik, month=1_ik, day=1_ik)
      call output_integrate(files, bufs, poly, day1, DT, .true., .false., .false.)
      call check_close(bufs%pending(2)%sval(k_t), 294.0_wp, 1.0e-12_wp, 'D: a polygon value''s last sub-step')
      c_tp = slab_col(bufs%pending(2), k_tp)
      call check(bufs%pending(2)%nslab(k_tp) == 1_ik, 'D: the patch snapshot has the staged patch')
      call check_close(bufs%pending(2)%slab(1,c_tp), 294.0_wp, 1.0e-12_wp, 'D: a patch slab''s last sub-step')
      call check_close(bufs%pending(2)%slab(1, slab_col(bufs%pending(2), k_r)), 410.0_wp, 1.0e-12_wp,     &
                       'D: an F-only state''s last sub-step')
      call check_close(bufs%pending(2)%sval(k_lt), 296.0_wp, 1.0e-12_wp,                             &
                       'D: a cohort row''s last sub-step, leaf-area-weighted as its mean is')

      bufs%fast_ready = .false.
      call output_integrate(files, bufs, poly, time_advance_days(day1, 1_ik), DT, .true., .false., .false.)
      call check(.not. bufs%pending(2)%svalid(k_t), 'D: a step with nothing staged has no snapshot')
      call free_polygon(poly)
   end subroutine test_snapshots

   !--- `"*" = false`: nothing selected, the axes' bookkeeping included.
   subroutine select_nothing(files)
      type(output_files_t), intent(inout) :: files
      call apply_selections(files%reg, [selection_t('*', OVR_FALSE, FREQ_NONE)])
      call resolve_streams(files%reg)
   end subroutine select_nothing

   !--- One patch holding `agb`'s cohorts, one plant each.
   subroutine set_cohorts(poly, agb)
      type(polygon_t), intent(inout) :: poly
      real(wp),        intent(in)    :: agb(:)
      integer(ik) :: n
      n = size(agb, kind=ik)
      poly%cohort%n = n
      poly%patch%cohort_count(1) = n
      poly%cohort%nplant(1:n) = 1.0_wp
      poly%cohort%pft(1:n) = 1_ik
      poly%cohort%owner_patch(1:n) = 1_ik
      call set_cohort_agb(poly, agb, n)
   end subroutine set_cohorts

   !--- Stage two FAST sub-steps with the given CAS temperatures and zero everything else.
   subroutine stage_cas(bufs, cas)
      type(output_buffers_t), intent(inout) :: bufs
      real(wp),               intent(in)    :: cas(2)
      call alloc_fast_staging(bufs, 2_ik, 1_ik, 8_ik, 2_ik)
      bufs%fast_n_cohort = 1_ik
      bufs%fast_site(PA_CAS_TEMP, :) = cas
      bufs%fast_patch(PA_CAS_TEMP, 1, :) = cas
      bufs%fast_soil_temp = 280.0_wp
      bufs%fast_soil_temp_patch = 280.0_wp
      bufs%fast_cohort_height = 10.0_wp
   end subroutine stage_cas

   !---~---
   !   One polygon's FAST staging as fast_dynamics sizes it (size_fast_staging): nsub sub-steps,
   ! np patches and nc cohort slots, all zero, in one site of n_soil layers.
   !---~---
   subroutine alloc_fast_staging(bufs, nsub, np, nc, n_soil)
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: nsub, np, nc, n_soil
      integer(ik) :: nl, ip
      nl = n_soil_layer_max
      allocate(bufs%fast_time(nsub), bufs%fast_forcing(N_PY_ROWS, nsub), bufs%fast_site(N_PA_ROWS, nsub),  &
               bufs%fast_patch(N_PA_ROWS, max(np, 1_ik), nsub))
      allocate(bufs%fast_soil_temp(nl, nsub), bufs%fast_soil_water(nl, nsub),                        &
               bufs%fast_soil_temp_patch(nl * max(np, 1_ik), nsub),                                  &
               bufs%fast_soil_water_patch(nl * max(np, 1_ik), nsub))
      allocate(bufs%fast_soil_live(nl), bufs%fast_soil_patch_live(nl * max(np, 1_ik)), source=.false.)
      bufs%fast_n_site = 1_ik
      bufs%fast_soil_live(1:n_soil) = .true.
      do ip = 1_ik, np
         bufs%fast_soil_patch_live((ip - 1_ik) * nl + 1_ik:(ip - 1_ik) * nl + n_soil) = .true.
      end do
      allocate(bufs%fast_cohort_gpp(nc, nsub), bufs%fast_cohort_height(nc, nsub))
      bufs%fast_forcing          = 0.0_wp
      bufs%fast_site             = 0.0_wp
      bufs%fast_patch            = 0.0_wp
      bufs%fast_soil_temp        = 0.0_wp
      bufs%fast_soil_water       = 0.0_wp
      bufs%fast_soil_temp_patch  = 0.0_wp
      bufs%fast_soil_water_patch = 0.0_wp
      bufs%fast_cohort_gpp       = 0.0_wp
      bufs%fast_cohort_height    = 0.0_wp
      bufs%n_fast_sub  = nsub
      bufs%fast_n_patch = np
   end subroutine alloc_fast_staging

end program test_output
