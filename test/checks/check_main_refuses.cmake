# SPDX-License-Identifier: Apache-2.0
#----- edmore refuses a run it cannot start, and names the cause. ------------------------------   #
#
# Each case edits the reference config (ed_config_main.toml), or a copy of the PFT file it names
# (ed_config_pft.toml), into a run that must stop before its first step, runs edmore on it, and
# asserts that it stopped and that its report names the cause:
#   config_unknown_keys        a misspelt key, a key that moved block, and a [region] block in a site run
#   config_bad_choice          a key set to none of its choices (run.mode = "globe")
#   config_inapplicable_keys   keys a setting leaves nothing to do (check_key_use): a lapse rate with the
#                              lapse off, and the [soil] keys the [maps] soil maps provide (the
#                              reference's own, and deep_temp)
#   config_rk45_loose_rtol     RK45 at a tolerance looser than it can integrate (fast.rtol > 1e-4)
#   config_malformed_value     a value the reader cannot parse as the number its key takes
#   output_fast_window         a sub-daily record that spans more fast steps than divide a slow step, so
#                              its window would straddle two slow steps
#   io_config_*                an [output].io_config line the selection refuses (ed_output_registry):
#                              a variable with no sub-daily value on F (coarse_on_f), one that exists
#                              only sub-daily on D (fast_on_coarse), `true` on a variable off by default
#                              (true_off), a pattern that matches nothing (no_match), two lines as
#                              specific as each other that disagree (tie), a tier none of a pattern's
#                              variables can be on (tier_none), and the axes' bookkeeping (axis)
#   maps_chave_e               a Chave E map with no PFT on the Chave height allometry
#   config_fineroot_plastic    fine roots turning over with the leaves' plastic rate
#                              (fineroot_turnover_plastic) while trait plasticity is off
#   hydraulics_pv_traits       a leaf elastic modulus below |pi0|, whose turgor-loss point is not
#                              negative, in a PFT file that sets its stress points by hand
#   hydraulics_vulnerability   a xylem psi50 above zero
#   soil_nan_layer             a NaN soil conductivity, as a map layer with no sand or clay gives; NaN
#                              fails every "<= 0" test, so the check must state what holds (Release
#                              builds only: a Debug build traps the NaN at that comparison)
#   toml_unclosed_string       a string value with no closing quote, which would otherwise read as ""
#   toml_duplicate_key         a key set twice ([soil] depth), which TOML forbids
#   init_census_missing        a census file that is not there (init mode 1)
#   init_census_empty          a census file with no row for the site
#   init_restart_missing       a state file that is not there (init mode 2)
# Every case runs in its own directory of the build tree, on a one-day run whose input paths are
# absolute: a refusal that stops working then costs a day's run there, never a long run in the
# source tree.
#
# Usage: cmake -DEDMORE_MAIN=<edmore> -DSOURCE_DIR=<repo> -DWORK_DIR=<build-tree dir> -DCASE=<case>
#              -P check_main_refuses.cmake
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.20)

#----- The edits the reference config must still accept; a key that moved stops the test here. ---#
function(replace_or_fail var pattern replacement what)
   string(REGEX MATCH "${pattern}" found "${${var}}")
   if(NOT found)
      message(FATAL_ERROR "the reference configuration no longer holds ${what}: the test cannot build its run")
   endif()
   string(REGEX REPLACE "${pattern}" "${replacement}" edited "${${var}}")
   set(${var} "${edited}" PARENT_SCOPE)
endfunction()

#----- The forcing on, from the ED_ERA5land archive at a site, the [site] keys given ending the   #
#      block: the keys load and check, and the archive is never opened, for each case using it     #
#      stops before the run starts. ------------------------------------------------------------#
function(archive_forcing var site_keys)
   set(${var} "
[site]
latitude = 42.44
longitude = -76.50
elevation = 320.0
apply_elevation_lapse = false
${site_keys}[forcing]
forcing_on = true
format = \"ERA5land\"
data_path = \"${WORK_DIR}/no_archive\"
max_distance_km = 15.0
tq_height = 2.0
wind_height = 10.0
height_above = \"zero_plane\"
wind_exposure = \"open_terrain\"
wind_exposure_z0 = 0.03
wind_blending_height = 40.0
timestep = \"3600s\"
avg_convention = \"end\"
sw_partition = \"clearidx\"
lwdown_source = \"file\"
co2_const = 420.0
recycle = false
start_clamp = \"hold\"
" PARENT_SCOPE)
endfunction()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(READ "${SOURCE_DIR}/ed_config_main.toml" text)

#----- A one-day run, its input files named by absolute path. ------------------------------------#
replace_or_fail(text "\nstart_time *= *\"[^\"]*\"" "\nstart_time = \"2000-01-01\"" "start_time")
replace_or_fail(text "\nend_time *= *\"[^\"]*\"" "\nend_time = \"2000-01-02\"" "end_time")
replace_or_fail(text "\npft_config *= *\"([^/\"][^\"]*)\"" "\npft_config = \"${SOURCE_DIR}/\\1\"" "pft_config")
replace_or_fail(text "\ncensus_file *= *\"([^/\"][^\"]*)\"" "\ncensus_file = \"${SOURCE_DIR}/\\1\"" "census_file")

if(CASE STREQUAL "config_unknown_keys")
   string(APPEND text "\n[forcing]\ntimestpe = \"3600s\"\n[state]\ncohort_max = 2048\n[region]\nland_fraction_min = 0.5\n")
   set(expects "forcing.timestpe .*did you mean forcing.timestep\\?"
               "state.cohort_max .*did you mean output.cohort_max\\?"
               "\\[region\\] .*needs \\[run\\]\\.mode = \"region\""
               "the configuration has missing, unknown or inapplicable keys")
elseif(CASE STREQUAL "config_bad_choice")
   replace_or_fail(text "\n\\[run\\][^\n]*\n" "\n[run]\nmode = \"globe\"\n" "a [run] block")
   set(expects "run\\.mode = \"globe\" is not one of site \\| region")
elseif(CASE STREQUAL "config_inapplicable_keys")
   #----- The lapse off with a lapse rate, and the soil maps beside the reference's [soil] keys. ---#
   archive_forcing(forcing "lapse_rate_tair = 0.0065\n")
   replace_or_fail(text "\n\\[soil\\][^\n]*\n" "\n[soil]\ndeep_temp = 290.0\n" "a [soil] block")
   string(APPEND text "${forcing}\n[maps]\nsoil_depth_map = \"${WORK_DIR}/soil_depth.nc\"\n"
                      "soil_texture_map = \"${WORK_DIR}/soil_texture.nc\"\n"
                      "soil_chemistry_map = \"${WORK_DIR}/soil_chemistry.nc\"\n")
   set(expects "site\\.lapse_rate_tair [^\n]*applies only with site\\.apply_elevation_lapse = true"
               "soil\\.theta_sat [^\n]*comes from the \\[maps\\] soil maps"
               "soil\\.deep_temp [^\n]*comes from the \\[maps\\] soil maps")
elseif(CASE STREQUAL "maps_chave_e")
   archive_forcing(forcing "")
   string(APPEND text "${forcing}\n[maps]\nchave_e_map = \"${WORK_DIR}/chave_e.nc\"\n")
   set(expects "maps\\.chave_e_map needs a PFT with allometry\\.height_allometry = \"chave2014\"")
elseif(CASE STREQUAL "output_fast_window")
   #----- The fast loop and the diagnostics on: 900 s steps, 96 to the day, and a 5-step record. ---#
   replace_or_fail(text "\nfast_biophysics_on *= *false" "\nfast_biophysics_on = true" "fast_biophysics_on = false")
   replace_or_fail(text "\nenabled *= *false" "\nenabled = true" "[output] enabled")
   replace_or_fail(text "\nfast_interval_steps *= *[0-9]+" "\nfast_interval_steps = 5" "[output] fast_interval_steps")
   set(expects "output\\.fast_interval_steps must divide the fast steps in a slow step")
elseif(CASE MATCHES "^io_config_")
   #----- The fast loop and the diagnostics on, and a file of the case's [variables] lines. --------#
   if(CASE STREQUAL "io_config_coarse_on_f")
      set(lines "mean_dbh_co = \"F D\"")
      set(expects "\"mean_dbh_co\" has no sub-daily value and cannot be put on the F stream")
   elseif(CASE STREQUAL "io_config_fast_on_coarse")
      set(lines "mean_reco_pa = \"D\"")
      set(expects "\"mean_reco_pa\" exists only sub-daily and can only be on the F stream")
   elseif(CASE STREQUAL "io_config_true_off")
      set(lines "cas_temp_py = true")
      set(expects "\"cas_temp_py\" is off by default, so true has no tiers to restore")
   elseif(CASE STREQUAL "io_config_no_match")
      set(lines "\"*_xyz\" = false")
      set(expects "io_config \"\\*_xyz\" matches no variable")
   elseif(CASE STREQUAL "io_config_tie")
      set(lines "\"mean_a*\" = \"D\"\n\"*agb_py\" = \"M\"")
      set(expects "select [a-z_]+ alike and disagree; make one more specific")
   elseif(CASE STREQUAL "io_config_tier_none")
      set(lines "\"*_co\" = \"Y\"")
      set(expects "io_config \"\\*_co\": none of its variables can be on Y")
   elseif(CASE STREQUAL "io_config_axis")
      set(lines "owner_patch = \"D\"")
      set(expects "io_config \"owner_patch\" is written with its axis")
   else()
      message(FATAL_ERROR "check_main_refuses: unknown CASE '${CASE}'")
   endif()
   file(WRITE "${WORK_DIR}/io_config.toml" "[variables]\n${lines}\n")
   replace_or_fail(text "\nfast_biophysics_on *= *false" "\nfast_biophysics_on = true" "fast_biophysics_on = false")
   replace_or_fail(text "\nenabled *= *false" "\nenabled = true" "[output] enabled")
   replace_or_fail(text "\nio_config *= *\"\"[^\n]*" "\nio_config = \"${WORK_DIR}/io_config.toml\"" "[output] io_config")
elseif(CASE STREQUAL "config_rk45_loose_rtol")
   #----- The fast loop on, RK45, rtol 1e-3: any line naming either key goes, then both are set. -#
   replace_or_fail(text "\nfast_biophysics_on *= *false" "\nfast_biophysics_on = true" "fast_biophysics_on = false")
   string(REGEX REPLACE "\n#? *time_integrator *=[^\n]*" "" text "${text}")
   string(REGEX REPLACE "\n#? *rtol *=[^\n]*" "" text "${text}")
   replace_or_fail(text "\n\\[fast\\][^\n]*\n" "\n[fast]\ntime_integrator = \"rk45\"\nrtol = 1.0e-3\n" "a [fast] block")
   set(expects "fast\\.rtol > 1e-4 with time_integrator = \"rk45\"")
elseif(CASE STREQUAL "config_malformed_value")
   replace_or_fail(text "\ndepth *= *[0-9.]+" "\ndepth = two" "[soil] depth")
   set(expects "cannot parse real for key \"soil\\.depth\"" "raw value: \"two\"" "malformed configuration value")
elseif(CASE STREQUAL "config_fineroot_plastic")
   #----- Trait plasticity off, whatever the reference says, and the fine roots plastic. ---------#
   string(REGEX REPLACE "\n#? *trait_plasticity_on *=[^\n]*" "" text "${text}")
   replace_or_fail(text "\n# *fineroot_turnover_plastic *= *false" "\nfineroot_turnover_plastic = true"
                   "[trait_dynamics] fineroot_turnover_plastic = false, commented out")
   set(expects "trait_dynamics\\.fineroot_turnover_plastic needs trait_plasticity_on = true")
elseif(CASE STREQUAL "soil_nan_layer")
   replace_or_fail(text "\nsolid_conductivity *= *[0-9.eE+-]+" "\nsolid_conductivity = nan" "[soil] solid_conductivity")
   set(expects "soil conductivities must be positive \\(a map layer needs sand or clay\\)")
elseif(CASE STREQUAL "toml_unclosed_string")
   #----- The diagnostics on, and the file of their variables named with no closing quote. --------#
   replace_or_fail(text "\nenabled *= *false" "\nenabled = true" "[output] enabled")
   replace_or_fail(text "\nio_config *= *\"\"[^\n]*" "\nio_config = \"my_variables.txt" "[output] io_config")
   set(expects "cannot parse string[^\n]*for key \"output\\.io_config\"")
elseif(CASE STREQUAL "toml_duplicate_key")
   replace_or_fail(text "\ndepth *= *([0-9.]+)" "\ndepth = \\1\ndepth = 3.0" "[soil] depth")
   set(expects "soil\\.depth[^\n]*twice")
elseif(CASE MATCHES "^hydraulics_")
   #----- The PFT file's copy beside the run, with one trait made impossible for every PFT. --------#
   file(READ "${SOURCE_DIR}/ed_config_pft.toml" pft)
   if(CASE STREQUAL "hydraulics_pv_traits")
      #----- The stress points set by hand, as a PFT file may: left to their default, the positive    #
      #      turgor-loss point would already be refused through them. -------------------------------#
      foreach(key wstress_psi_open wstress_psi_close stomata_psi_onset)
         replace_or_fail(pft "\n# *${key}( *=)" "\n${key}\\1" "[pft] ${key}, commented out")
      endforeach()
      replace_or_fail(pft "\nleaf_elastic_mod *= *[^ #\n]+" "\nleaf_elastic_mod = 1.0" "[hydraulics] leaf_elastic_mod")
      set(expects "hydraulics\\.leaf_elastic_mod and wood_elastic_mod must exceed \\|pi0\\|")
   else()
      replace_or_fail(pft "\nwood_psi50 *= *[^ #\n]+" "\nwood_psi50 = 0.5" "[hydraulics] wood_psi50")
      set(expects "hydraulics\\.wood_psi50 must be < 0")
   endif()
   file(WRITE "${WORK_DIR}/ed_config_pft.toml" "${pft}")
   replace_or_fail(text "\npft_config *= *\"[^\"]*\"" "\npft_config = \"${WORK_DIR}/ed_config_pft.toml\"" "pft_config")
elseif(CASE MATCHES "^init_")
   if(CASE STREQUAL "init_restart_missing")
      set(mode 2)
      set(key "restart_file")
      set(file "${WORK_DIR}/no_such_state.nc")
   else()
      set(mode 1)
      set(key "census_file")
      if(CASE STREQUAL "init_census_empty")
         set(file "${WORK_DIR}/empty_census.csv")
         file(WRITE "${file}" "site_id,patch_id,cohort_id,dbh,height,pft,nplant\n")
      else()
         set(file "${WORK_DIR}/no_such_census.csv")
      endif()
   endif()
   string(REGEX REPLACE "([][+.*()^$?|\\\\])" "\\\\\\1" file_re "${file}")   # the path, matched literally
   if(CASE STREQUAL "init_restart_missing")
      set(expects "restart \\(mode 2\\): the state file ${file_re} cannot be opened")
   elseif(CASE STREQUAL "init_census_empty")
      set(expects "census \\(mode 1\\): the census file ${file_re} holds no cohort for the site")
   else()
      set(expects "census \\(mode 1\\): the census file ${file_re} cannot be opened")
   endif()
   string(REGEX REPLACE "\n${key} *= *\"[^\"]*\"[^\n]*" "" text "${text}")
   replace_or_fail(text "\ninit_mode *= *[0-9]" "\ninit_mode = ${mode}\n${key} = \"${file}\"" "init_mode")
else()
   message(FATAL_ERROR "check_main_refuses: unknown CASE '${CASE}'")
endif()

file(WRITE "${WORK_DIR}/ed_config_main.toml" "${text}")
execute_process(COMMAND "${EDMORE_MAIN}" "${WORK_DIR}/ed_config_main.toml" WORKING_DIRECTORY "${WORK_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
set(log "${out}${err}")
set(problems "")
if(rc EQUAL 0)
   string(APPEND problems "\n  the run did not stop")
endif()
if(CASE MATCHES "^init_" AND log MATCHES "bare ground")
   string(APPEND problems "\n  the run started from bare ground")
endif()
foreach(expect IN LISTS expects)
   if(NOT log MATCHES "${expect}")
      string(APPEND problems "\n  the report lacks: ${expect}")
   endif()
endforeach()
if(problems)
   message(FATAL_ERROR "edmore did not refuse case ${CASE} as expected:${problems}\n--- output:\n${log}")
endif()
message(STATUS "edmore refused case ${CASE}, naming the cause")
