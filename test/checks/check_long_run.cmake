# SPDX-License-Identifier: Apache-2.0
#----- A long run completes: example 05's Palo Verde inventory or its 150-year regrowth. ----------#
#
# The two runs of examples/example05_forest_regeneration/palo_verde (run_palo_verde.py), through
# edmore, each from a copy of the example's config in its own directory of the build tree:
#   inventory  the 2009 inventory, July 2008 - December 2010 (ed_config_inventory.toml), which reads
#              the ED_ERA5land archive (ARCHIVE) itself
#   regrowth   150 years from bare ground, June 1871 - May 2021 (ed_config_regrowth.toml), on Palo
#              Verde's cell of the archive, 2003-2022, repeated: the file FORCING names, or, without one,
#              the cut this script makes into the work directory once, as run_palo_verde.py makes it
# The copy names the example's input files by absolute path and writes its output, its log
# (<case>.log) and nothing else into the work directory. A case passes when edmore completes (it
# reports "OK: simulation completed") and the run's slow ledger closes: its verdict is that every phase
# closes within tolerance, and it reports no impossible store.
#
# Usage: cmake -DEDMORE_MAIN=<edmore> -DSOURCE_DIR=<repo> -DWORK_DIR=<build-tree dir> -DCASE=<case>
#              -DARCHIVE=<ED_ERA5land folder> [-DFORCING=<file>] [-DPYTHON=<python>] -P check_long_run.cmake
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.20)

set(EXAMPLE "${SOURCE_DIR}/examples/example05_forest_regeneration/palo_verde")
#----- Palo Verde's plots [deg N, deg E], where run_palo_verde.py cuts the forcing. ----------------#
set(PV_LAT 10.35)
set(PV_LON -85.35)

#----- The edits the example's config must still accept; a key that moved stops the test here. ---#
function(replace_or_fail var pattern replacement what)
   string(REGEX MATCH "${pattern}" found "${${var}}")
   if(NOT found)
      message(FATAL_ERROR "ed_config_${CASE}.toml no longer holds ${what}: the test cannot build its run")
   endif()
   string(REGEX REPLACE "${pattern}" "${replacement}" edited "${${var}}")
   set(${var} "${edited}" PARENT_SCOPE)
endfunction()

if(NOT CASE MATCHES "^(inventory|regrowth)$")
   message(FATAL_ERROR "check_long_run: unknown CASE '${CASE}'")
endif()
if(NOT IS_DIRECTORY "${ARCHIVE}")
   message(FATAL_ERROR "the ${CASE} run needs the ED_ERA5land archive (\"${ARCHIVE}\" is not a folder): configure "
                       "with -DED_ERA5LAND_ARCHIVE=<folder>, or set ERA5LAND_ARCHIVE before configuring")
endif()
file(MAKE_DIRECTORY "${WORK_DIR}")
file(REMOVE_RECURSE "${WORK_DIR}/output")
file(READ "${EXAMPLE}/ed_config_${CASE}.toml" text)

#----- Every input by absolute path: the census, the PFT file, the output selection and the CO2      #
#      history are the example's; the output stays relative, in the work directory. ----------------#
foreach(key census_file pft_config io_config co2_file)
   replace_or_fail(text "\n${key}( *= *)\"([^/\"][^\"]*)\"" "\n${key}\\1\"${EXAMPLE}/\\2\"" "${key}")
endforeach()

if(CASE STREQUAL "inventory")
   replace_or_fail(text "\ndata_path( *= *)\"[^\"]*\"" "\ndata_path\\1\"${ARCHIVE}\"" "[forcing].data_path")
else()
   #----- The forcing: the file named, or Palo Verde's cell cut from the archive, 2003-2022 -- once,   #
   #      to a part file renamed when the cut completes, so a cut that stopped is never read. --------#
   if(NOT FORCING)
      set(FORCING "${WORK_DIR}/data/pv_era5land_2003-2022.nc")
      if(NOT EXISTS "${FORCING}")
         if(NOT PYTHON)
            message(FATAL_ERROR "the regrowth's forcing must be cut from the archive, and no Python was found: "
                                "configure with -DED_TEST_PV_FORCING=<its file>, cut by run_palo_verde.py")
         endif()
         file(MAKE_DIRECTORY "${WORK_DIR}/data")
         string(TIMESTAMP t_cut "%s")
         execute_process(COMMAND "${PYTHON}" "${SOURCE_DIR}/scripts/prepare_input/era5/make_forcing_file.py"
                                 --data-path "${ARCHIVE}" --start 2003-01-01 --end 2022-12-31
                                 --lat ${PV_LAT} --lon ${PV_LON} --out "${FORCING}.part"
                         OUTPUT_FILE "${WORK_DIR}/forcing_cut.log" ERROR_FILE "${WORK_DIR}/forcing_cut.log"
                         RESULT_VARIABLE rc)
         if(NOT rc EQUAL 0)
            message(FATAL_ERROR "cutting Palo Verde's forcing from ${ARCHIVE} failed; see ${WORK_DIR}/forcing_cut.log")
         endif()
         file(RENAME "${FORCING}.part" "${FORCING}")
         string(TIMESTAMP t_now "%s")
         math(EXPR dt_cut "${t_now} - ${t_cut}")
         message(STATUS "cut Palo Verde's forcing, 2003-2022, in ${dt_cut} s: ${FORCING}")
      endif()
   elseif(NOT EXISTS "${FORCING}")
      message(FATAL_ERROR "the regrowth's forcing file ${FORCING} (ED_TEST_PV_FORCING) is not there")
   endif()
   replace_or_fail(text "\npath( *= *)\"[^\"]*\"" "\npath\\1\"${FORCING}\"" "[forcing].path")
endif()

file(WRITE "${WORK_DIR}/ed_config_${CASE}.toml" "${text}")
set(log_file "${WORK_DIR}/${CASE}.log")
string(TIMESTAMP t_start "%s")
execute_process(COMMAND "${EDMORE_MAIN}" "${WORK_DIR}/ed_config_${CASE}.toml" WORKING_DIRECTORY "${WORK_DIR}"
                OUTPUT_FILE "${log_file}" ERROR_FILE "${log_file}" RESULT_VARIABLE rc)
string(TIMESTAMP t_end "%s")
math(EXPR dt_run "${t_end} - ${t_start}")
file(READ "${log_file}" log)

set(problems "")
if(NOT rc EQUAL 0)
   string(APPEND problems "\n  edmore stopped (exit status ${rc})")
endif()
if(NOT log MATCHES "OK: simulation completed")
   string(APPEND problems "\n  the run did not report completing")
endif()
if(NOT log MATCHES "VERDICT: every phase closes within tolerance on all three currencies")
   string(APPEND problems "\n  the slow ledger did not report every phase closing")
endif()
if(log MATCHES "IMPOSSIBLE STORE")
   string(APPEND problems "\n  the slow ledger reports an impossible store")
endif()
if(problems)
   message(FATAL_ERROR "the Palo Verde ${CASE} run failed after ${dt_run} s:${problems}\n  its log: ${log_file}")
endif()
message(STATUS "the Palo Verde ${CASE} run completed in ${dt_run} s, its slow ledger closed; its log: ${log_file}")
