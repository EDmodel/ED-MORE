# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# The MEDS test suite (docs/building.md, "Running the tests"): one test per src module, in   #
# folders that mirror src/, its support libraries, the source-tree checks and the Python     #
# tools' tests.                                                                              #
#                                                                                            #
# Every test is registered by meds_add_test (below) with exactly one label:                  #
#   unit         a module's test, and meds_main's refusals: seconds each                     #
#   integration  whole runs: several sites, a region, the fast loop's driver, the C-API run  #
#                and a restart                                                               #
#   lint         the source-tree rules and the tracked copies                                #
#   python       the Python tools' and the package's own tests (pytest)                      #
#   long         opt-in, -DMEDS_TEST_LONG=ON: example 05's Palo Verde inventory and its       #
#                150-year regrowth                                                           #
#                                                                                            #
# Paths: CMAKE_CURRENT_SOURCE_DIR is this folder; PROJECT_SOURCE_DIR the repository and      #
# PROJECT_BINARY_DIR the build tree's top, where every test runs unless it names another      #
# directory and where the work directories are.                                              #
#==========================================================================================#

#----- meds_add_test(<name> [SOURCES <file>...] [LINK <library>...] [COMMAND <command>...]      #
#                    [ARGS <argument>...] [WORKDIR <dir>] [ENVIRONMENT <VAR=value>...]          #
#                    LABELS <label>)                                                            #
#      The test <name>. With SOURCES it builds the program test_<name> from them, linked to     #
#      LINK, and runs that; otherwise it runs COMMAND. ARGS follow on the command line. It runs  #
#      in WORKDIR, made at configure time (the build tree's top when none is given), with        #
#      ENVIRONMENT. It takes exactly one of the labels above. The TIMEOUT every test gets is set #
#      at the end of this file.                                                                  #
function(meds_add_test name)
   cmake_parse_arguments(PARSE_ARGV 1 arg "" "WORKDIR" "SOURCES;LINK;COMMAND;ARGS;ENVIRONMENT;LABELS")
   if(arg_UNPARSED_ARGUMENTS)
      message(FATAL_ERROR "meds_add_test(${name}): unknown arguments: ${arg_UNPARSED_ARGUMENTS}")
   endif()
   list(LENGTH arg_LABELS n_labels)
   if(NOT n_labels EQUAL 1 OR NOT arg_LABELS MATCHES "^(unit|integration|lint|python|long)$")
      message(FATAL_ERROR "meds_add_test(${name}): give exactly one label of unit, integration, lint, python "
                          "and long")
   endif()
   if(arg_SOURCES)
      add_executable(test_${name} ${arg_SOURCES})
      target_link_libraries(test_${name} PRIVATE ${arg_LINK})
      meds_fortran_flags(test_${name})
      set(command test_${name})
   else()
      set(command ${arg_COMMAND})
   endif()
   if(NOT arg_WORKDIR)
      set(arg_WORKDIR ${PROJECT_BINARY_DIR})
   endif()
   file(MAKE_DIRECTORY ${arg_WORKDIR})
   add_test(NAME ${name} COMMAND ${command} ${arg_ARGS} WORKING_DIRECTORY ${arg_WORKDIR})
   set_tests_properties(${name} PROPERTIES LABELS ${arg_LABELS})
   if(arg_ENVIRONMENT)
      set_tests_properties(${name} PROPERTIES ENVIRONMENT "${arg_ENVIRONMENT}")
   endif()
endfunction()

#==========================================================================================#
#  support/: what the tests share.                                                            #
#==========================================================================================#
#----- The assertions, the soil grid and the derived configs (meds_test_run), depending on nothing #
#      but a kind. Kept apart from meds_testsupport because many tests deliberately link ONE      #
#      narrow library each to keep that layer standalone-buildable, and they must be able to      #
#      assert without linking the whole model.                                                    #
add_library(meds_testassert STATIC support/meds_test_assert.f90 support/meds_test_soil_grid.f90
                                   support/meds_test_run.f90)
target_link_libraries(meds_testassert PUBLIC meds_share)
meds_fortran_flags(meds_testassert)
#----- The in-code config fixture; re-exports the assertions so most tests need one `use`. ------#
add_library(meds_testsupport STATIC support/meds_test_support.f90)
target_link_libraries(meds_testsupport PUBLIC meds_config meds_testassert)
meds_fortran_flags(meds_testsupport)
#----- The synthetic ED_ERA5land archive, shared by the reader, the maps and the region tests. ---#
add_library(meds_testarchive STATIC support/meds_test_archive.f90)
target_link_libraries(meds_testarchive PUBLIC meds_met)
meds_fortran_flags(meds_testarchive)
#----- The fast column the integrator tests run on: assembled as a run assembles one, its forcing, #
#      its frozen record and the ledger check.                                                      #
add_library(meds_testfastfixture STATIC support/meds_test_fast_fixture.f90)
target_link_libraries(meds_testfastfixture PUBLIC meds_fast meds_testsupport)
meds_fortran_flags(meds_testfastfixture)

#==========================================================================================#
#  share/: the libraries of src/share, one test each, named after its module. Each links      #
#  meds_share alone; test_optics_lib also the config layer, whose meds_pft_derive turns a     #
#  PFT's leaf-angle moments into the LIDF's Beta parameters.                                  #
#==========================================================================================#
foreach(t numerics therm_lib plant_hydr_lib soil_hydr_lib optics_lib radiation_lib)
   meds_add_test(${t} SOURCES share/test_${t}.f90 LINK meds_share meds_testassert LABELS unit)
endforeach()
target_link_libraries(test_optics_lib PRIVATE meds_config)
#----- meds_allometry_lib: the size laws, and the shipped [allometry] block equal to the library's  #
#      initializers (the Python library's config-free values), which are Chave 2014 in carbon and   #
#      ED2's leaf fit over C2B. meds_config only for the TOML reader.                                #
meds_add_test(allometry_lib SOURCES share/test_allometry_lib.f90 LINK meds_config meds_testassert
              ARGS ${PROJECT_SOURCE_DIR}/meds_config_pft.toml LABELS unit)

#==========================================================================================#
#  config/: the configuration layer (src/config).                                             #
#==========================================================================================#
#----- nearest_key, and the misspelt, moved and unknown keys a run reports. ----------------------#
meds_add_test(config_keys SOURCES config/test_config_keys.f90 LINK meds_config meds_testsupport LABELS unit)
#----- The config blocks reach their records: the [soil]/[fast]/[snow]/[aerodynamics] loaders, and #
#      the [soil] column through build_soil_profile.                                                #
meds_add_test(params_config SOURCES config/test_params_config.f90 LINK meds_config meds_testassert LABELS unit)
#----- The PFT file's loader and the per-PFT tables derived from it: the shipped file through the  #
#      loader, and the optics and hydraulics that differentiate PFTs. It runs from the SOURCE dir,  #
#      where meds_config_main.toml names the PFT file, and writes its CSV in a directory of its own #
#      in the build tree.                                                                            #
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/work/pft_config)
meds_add_test(pft_config SOURCES config/test_pft_config.f90 LINK meds_model meds_testsupport
              ARGS ${PROJECT_BINARY_DIR}/work/pft_config WORKDIR ${PROJECT_SOURCE_DIR} LABELS unit)

#==========================================================================================#
#  state/: the state layer's containers (src/state).                                          #
#==========================================================================================#
meds_add_test(containers SOURCES state/test_containers.f90 LINK meds_testsupport meds_model LABELS unit)

#==========================================================================================#
#  fast/: the fast loop (src/fast_dynamics).                                                  #
#==========================================================================================#
#----- The process modules (processes/), one test each, named after its module. Each links the    #
#      kernel library alone, so the kernels build and test standalone; the leaf test also the      #
#      in-code config fixture, whose leaf table it solves on.                                      #
foreach(t canopy_aerodynamics canopy_air_space canopy_radiation canopy_water leaf_gas_exchange plant_energy
          plant_hydraulics respiration soil_heat soil_water surface)
   meds_add_test(${t} SOURCES fast/test_${t}.f90 LINK meds_fast_kernels meds_testassert LABELS unit)
endforeach()
target_link_libraries(test_leaf_gas_exchange PRIVATE meds_testsupport)
#----- The integrator (integrator/), one test per concern, named after it: its building blocks,   #
#      its two schemes and the march that drives them, and one fast step end to end. Each links    #
#      the fast loop and the shared fast fixture (support/).                                       #
foreach(t fast_integ_kernels fast_integ_schemes fast_integ_step)
   meds_add_test(${t} SOURCES fast/test_${t}.f90 LINK meds_testfastfixture LABELS unit)
endforeach()
#----- The soil-plant water network on synthetic stands and soils, and the soil column alone. -----#
meds_add_test(water_network SOURCES fast/test_water_network.f90 LINK meds_fast meds_testassert LABELS unit)
#----- The fast loop's driver (driver/) and the stepper that runs it, on a forcing file the test   #
#      writes in a directory of its own in the build tree.                                         #
meds_add_test(fast_driver SOURCES fast/test_fast_driver.f90 LINK meds_testsupport meds_model
              WORKDIR ${PROJECT_BINARY_DIR}/work/fast_driver LABELS integration)

#==========================================================================================#
#  slow/: the slow loop (src/slow_dynamics).                                                  #
#==========================================================================================#
#----- The daily plant kernels (allocation, traits, demography) and the soil carbon pools' daily   #
#      step, which debits them through the decomposition matrix of meds_respiration (in            #
#      meds_fast_kernels): the kernel libraries, state-free.                                       #
foreach(t plant_carbon_allocation plant_traits demography soil_carbon)
   meds_add_test(${t} SOURCES slow/test_${t}.f90 LINK meds_fast_kernels meds_slow_kernels meds_testassert
                 LABELS unit)
endforeach()
#----- The phenology kernel and its C-API shim, so the shim compiles under ctest. ----------------#
meds_add_test(plant_phenology SOURCES slow/test_plant_phenology.f90
              LINK meds_slow_kernels meds_c_api meds_testassert LABELS unit)
#----- The driver: the stand's operators, the slow-loop site conservation ledger and the step. ----#
foreach(t slow_update slow_step)
   meds_add_test(${t} SOURCES slow/test_${t}.f90 LINK meds_testsupport meds_model LABELS unit)
endforeach()
meds_add_test(slow_ledger SOURCES slow/test_slow_ledger.f90 LINK meds_model meds_testsupport LABELS unit)

#==========================================================================================#
#  io/: the initial state, the forcing, the maps, the output and the state file (src/io).     #
#==========================================================================================#
#----- How a run's polygon starts: the census reader and every reason it refuses a file, and the  #
#      initial soil state.                                                                         #
meds_add_test(init SOURCES io/test_init.f90 LINK meds_testsupport meds_model LABELS unit)
#----- A MEDS forcing file: its contract, the reader, and the forcing at a site. The met tests     #
#      write their files in a directory of their own in the build tree.                            #
meds_add_test(met_file SOURCES io/test_met_file.f90 LINK meds_met meds_testassert
              ARGS ${PROJECT_SOURCE_DIR}/data/co2/co2_cmip7_global_annual_1000-2022.txt
              WORKDIR ${PROJECT_BINARY_DIR}/work/met_file LABELS unit)
#----- The ED_ERA5land reader, on the synthetic archive. -----------------------------------------#
meds_add_test(met_era5land SOURCES io/test_met_era5land.f90 LINK meds_met meds_testassert meds_testarchive
              WORKDIR ${PROJECT_BINARY_DIR}/work/met_era5land LABELS unit)
#----- The static-map reader and the cell mask, on the synthetic archive's grid. ------------------#
meds_add_test(maps SOURCES io/test_maps.f90 LINK meds_maps meds_testassert meds_testarchive
              ARGS ${PROJECT_BINARY_DIR}/work/maps LABELS unit)
#----- The diagnostic output's netCDF-free half (src/io/output): the registry and the selection,   #
#      the reduction, the extract and the per-tier integration, one section per module. It links   #
#      meds_output rather than the serializer.                                                      #
meds_add_test(output SOURCES io/test_output.f90 LINK meds_output meds_testsupport LABELS unit)
#----- The serializer's round trip: writes and reopens netCDF, so it links the io layer. ---------#
meds_add_test(output_roundtrip SOURCES io/test_output_roundtrip.f90 LINK meds_output_stream meds_testsupport
              LABELS unit)
#----- The state file: written, read back and written again, every variable (links the IO layer). #
meds_add_test(state_file SOURCES io/test_state_file.f90 LINK meds_model meds_output_stream meds_testsupport
              WORKDIR ${PROJECT_BINARY_DIR}/work/state_file LABELS unit)

#==========================================================================================#
#  driver/: the coupled model's runs (src/main).                                              #
#==========================================================================================#
#----- A polygon of several sites: its sites from [[sites]], their files, a census, their stands,   #
#      threads and state file, and what meds_main refuses. Its configs derive from example 04's,    #
#      so it runs from that example's directory; it writes only under the build tree.               #
meds_add_test(sites SOURCES driver/test_sites.f90 LINK meds_testsupport meds_run
              ARGS ${PROJECT_BINARY_DIR}/work/sites $<TARGET_FILE:meds_main>
              WORKDIR ${PROJECT_SOURCE_DIR}/examples/example04_column_biophysics LABELS integration)
#----- A region against its site runs, on the synthetic archive. It runs from driver/region, whose  #
#      site config it derives its configs from; it writes everything under the work directory in   #
#      the build tree.                                                                               #
meds_add_test(region SOURCES driver/test_region.f90 LINK meds_run meds_testassert meds_testarchive
              ARGS ${PROJECT_BINARY_DIR}/region_test_tmp $<TARGET_FILE:meds_main>
              WORKDIR ${CMAKE_CURRENT_SOURCE_DIR}/driver/region LABELS integration)

#==========================================================================================#
#  c_api/: the C-API shims (src/c_api). Linking every shim INTO a ctest target is the point: a #
#  shim built only by the optional -DMEDS_BUILD_PYLIB=ON shared library breaks unseen -- a     #
#  component inserted mid-type in leaf_photo_args_t, or a signature changed under the          #
#  demography shim, leaves the whole suite green on both back ends while the C API is broken.  #
#  Linking it here makes that class of break a normal build failure. Needs no .so and no       #
#  Python: each test calls the bind(c) procedures directly, and links meds_c_api, the shims'   #
#  one compile, so a default build compiles and links every shim in src/c_api/.                #
#==========================================================================================#
#----- Runs from the SOURCE dir: the canopy shim opens the reference configuration and the PFT     #
#      file it names, relative to it, and writes nothing.                                           #
meds_add_test(c_api_leaf SOURCES c_api/test_c_api_leaf.f90 LINK meds_testsupport meds_c_api
              WORKDIR ${PROJECT_SOURCE_DIR} LABELS unit)
#----- Runs from the SOURCE dir: it loads examples/example03_demography's shipped config, which    #
#      gives that config a build-time consumer. A config nothing loads drifts behind the schema     #
#      exactly like a shim nothing compiles.                                                         #
meds_add_test(c_api_demography SOURCES c_api/test_c_api_demography.f90 LINK meds_testsupport meds_c_api
              WORKDIR ${PROJECT_SOURCE_DIR} LABELS unit)
#----- Runs from examples/example04_column_biophysics: it loads both of that example's shipped     #
#      configs (the default and the calibrated run), so the build notices when they fall behind    #
#      the config reader. It writes only under the build tree.                                      #
meds_add_test(c_api_run SOURCES c_api/test_c_api_run.f90 LINK meds_testsupport meds_c_api
              ARGS ${PROJECT_BINARY_DIR}/c_api_run
              WORKDIR ${PROJECT_SOURCE_DIR}/examples/example04_column_biophysics LABELS integration)

#==========================================================================================#
#  checks/: the CMake scripts that check the source tree and meds_main as a program.          #
#==========================================================================================#
#----- No procedure is passed as an argument anywhere in src/: the ifx cost of doing so serialises #
#      threads. The script explains the rule it checks.                                            #
meds_add_test(no_procedure_arguments
              COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${PROJECT_SOURCE_DIR}/src
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_no_procedure_arguments.cmake
              LABELS lint)
#----- The two reference configs list exactly the keys the loader reads: a run stops on any key its #
#      reference does not list. The script explains the rule it checks.                              #
meds_add_test(config_keys_listed
              COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${PROJECT_SOURCE_DIR}/src
                      -DMAIN_TOML=${PROJECT_SOURCE_DIR}/meds_config_main.toml
                      -DPFT_TOML=${PROJECT_SOURCE_DIR}/meds_config_pft.toml
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_config_keys.cmake
              LABELS lint)
#----- The shipped per-variable override EXAMPLE is a generated file that is also tracked, so it can #
#      drift silently from the registry: it must be what meds_main --dump-io-config writes from the  #
#      reference configs. The script explains the rule.                                              #
meds_add_test(io_config_example
              COMMAND ${CMAKE_COMMAND} -DMEDS_MAIN=$<TARGET_FILE:meds_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                      -DWORK_DIR=${PROJECT_BINARY_DIR}/work/io_config_example
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_io_config_example.cmake
              LABELS lint)
#----- The Python wheel ships COPIES of LICENSE and NOTICE, because python/pyproject.toml cannot    #
#      list license files above its own directory. A tracked copy drifts silently, for the same     #
#      reason as the io_config_example above, so the build compares them.                            #
foreach(f LICENSE NOTICE)
   string(TOLOWER ${f} f_lower)
   meds_add_test(${f_lower}_copy COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
                 ${PROJECT_SOURCE_DIR}/${f} ${PROJECT_SOURCE_DIR}/python/${f} LABELS lint)
endforeach()
#----- meds_main refuses a run it cannot start and names the cause: unknown, moved or misplaced     #
#      keys, a choice none of its key's, keys a setting makes inapplicable; RK45 at a loose          #
#      tolerance; a value that is not a number, a string with no closing quote or a key set twice;  #
#      a sub-daily record across two slow steps; io_config lines the selection refuses; a Chave E   #
#      map no PFT uses; plastic fine roots without plastic traits; impossible hydraulic traits; a   #
#      NaN soil layer; a census or state file that is missing or empty. Each case is a one-day run  #
#      in its own build-tree directory; the script explains the cases. A Debug build traps the NaN  #
#      at the comparison that refuses it, so that case runs in the other builds.                     #
set(refused_cases config_unknown_keys config_bad_choice config_inapplicable_keys config_rk45_loose_rtol
                  config_malformed_value output_fast_window io_config_coarse_on_f io_config_fast_on_coarse
                  io_config_true_off io_config_no_match io_config_tie io_config_tier_none io_config_axis
                  maps_chave_e hydraulics_pv_traits hydraulics_vulnerability init_census_missing
                  init_census_empty init_restart_missing)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
   list(APPEND refused_cases soil_nan_layer)
endif()
list(APPEND refused_cases toml_unclosed_string toml_duplicate_key config_fineroot_plastic)
foreach(case IN LISTS refused_cases)
   meds_add_test(${case}_refused
                 COMMAND ${CMAKE_COMMAND} -DMEDS_MAIN=$<TARGET_FILE:meds_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                         -DWORK_DIR=${PROJECT_BINARY_DIR}/refused/${case} -DCASE=${case}
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_main_refuses.cmake
                 LABELS unit)
endforeach()
#----- A step that fails stops meds_main after the run is finalized, which closes the files. The    #
#      script explains the case.                                                                      #
meds_add_test(main_failure_finalizes
              COMMAND ${CMAKE_COMMAND} -DMEDS_MAIN=$<TARGET_FILE:meds_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                      -DWORK_DIR=${PROJECT_BINARY_DIR}/main_failure
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_main_failure.cmake
              LABELS unit)

#==========================================================================================#
#  The rules in scripts/lint, and the Python tests, when the Python CMake finds has what each   #
#  needs; otherwise they are left out, and configure says so.                                   #
#==========================================================================================#
find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
   #----- No bare literal clamps in src/: a guard uses tiny_num, a threshold is a named setting, or  #
   #      the line says why. scripts/lint/check_clamps.py explains the rule.                          #
   meds_add_test(clamps_named COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/lint/check_clamps.py
                 ${PROJECT_SOURCE_DIR}/src LABELS lint)
   #----- No issue or pull-request numbers outside CHANGELOG.md, docs/ROADMAP.md and docs/dev_plans/:  #
   #      a comment says why the code is as it is; the history is in the CHANGELOG.                     #
   meds_add_test(no_issue_refs COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/lint/check_issue_refs.py
                 ${PROJECT_SOURCE_DIR} LABELS lint)
   #----- The kernels take only the per-patch store types from the state layer, never polygon_t: the  #
   #      shared module directory cannot enforce it, so this does, on every run.                      #
   meds_add_test(kernels_patch_only
                 COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/lint/check_kernel_state.py
                 ${PROJECT_SOURCE_DIR}/src LABELS lint)
   #----- The integrator steps one patch: it takes the site's surface and the per-patch store types,  #
   #      never the polygon, the forcing reader, the output or the driver.                            #
   meds_add_test(integrator_boundary
                 COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/lint/check_integrator_boundary.py
                 ${PROJECT_SOURCE_DIR}/src LABELS lint)
   #----- A named constant is spelled by its kind: lower case in meds_constants and meds_kinds, UPPER  #
   #      everywhere else. scripts/lint/check_constant_case.py explains the rule.                     #
   meds_add_test(constants_named
                 COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/lint/check_constant_case.py
                 ${PROJECT_SOURCE_DIR} LABELS lint)

   #----- The flux-tower forcing tool's own tests (scripts/prepare_flux_tower/tests). ----------------#
   execute_process(COMMAND ${Python3_EXECUTABLE} -c
                   "import sys, numpy, pandas, netCDF4, pytest; sys.version_info >= (3, 11) or __import__('tomli')"
                   RESULT_VARIABLE meds_tower_tool_deps OUTPUT_QUIET ERROR_QUIET)
   if (meds_tower_tool_deps EQUAL 0)
      meds_add_test(prepare_flux_tower
                    COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                            ${PROJECT_SOURCE_DIR}/scripts/prepare_flux_tower/tests
                    ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 LABELS python)
      #----- The fast-parameter calibration tool (scripts/calibrate_fast/tests): the unit tests, a   #
      #      smoke test that runs a chain, a trial and a Jacobian column through meds_main (and,    #
      #      when libmeds.so is built, through the Python API), and a short end-to-end fit.         #
      set(meds_calib_env "PYTHONDONTWRITEBYTECODE=1" "MEDS_MAIN=$<TARGET_FILE:meds_main>")
      if (MEDS_BUILD_PYLIB)
         list(APPEND meds_calib_env "MEDS_LIB=$<TARGET_FILE:meds_py>")
      endif()
      meds_add_test(calibrate_fast
                    COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                            ${PROJECT_SOURCE_DIR}/scripts/calibrate_fast/tests
                    ENVIRONMENT ${meds_calib_env} LABELS python)
      #----- The Python package's own tests (python/tests): the leaf and phenology kernels,        #
      #      meds.config, and meds.model against meds_main, bit for bit.                            #
      if (MEDS_BUILD_PYLIB)
         meds_add_test(python_api
                       COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider ${PROJECT_SOURCE_DIR}/python/tests
                       ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=${PROJECT_SOURCE_DIR}/python
                                   MEDS_LIB=$<TARGET_FILE:meds_py> MEDS_MAIN=$<TARGET_FILE:meds_main>
                       LABELS python)
      endif()
   else()
      message(STATUS "MEDS: prepare_flux_tower tests skipped (${Python3_EXECUTABLE} lacks numpy, pandas, "
                     "netCDF4, pytest or tomli)")
   endif()
   #----- The census preparation tool's tests (scripts/prepare_census/tests), on the same terms. ----#
   execute_process(COMMAND ${Python3_EXECUTABLE} -c
                   "import sys, numpy, pandas, pytest; sys.version_info >= (3, 11) or __import__('tomli')"
                   RESULT_VARIABLE meds_census_tool_deps OUTPUT_QUIET ERROR_QUIET)
   if (meds_census_tool_deps EQUAL 0)
      meds_add_test(prepare_census
                    COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                            ${PROJECT_SOURCE_DIR}/scripts/prepare_census/tests
                    ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 LABELS python)
   else()
      message(STATUS "MEDS: prepare_census tests skipped (${Python3_EXECUTABLE} lacks numpy, pandas, "
                     "pytest or tomli)")
   endif()
   #----- A restart continues exactly, and re-acclimates its traits on request                        #
   #      (python/test_restart_exact.py; meds_main on the constant reference climate). A whole run,    #
   #      so it is labelled integration, though pytest drives it. It needs only numpy, netCDF4 and     #
   #      pytest, so it does not wait on the tools' pandas.                                              #
   execute_process(COMMAND ${Python3_EXECUTABLE} -c
                   "import sys, numpy, netCDF4, pytest; sys.version_info >= (3, 11) or __import__('tomli')"
                   RESULT_VARIABLE meds_restart_test_deps OUTPUT_QUIET ERROR_QUIET)
   if (meds_restart_test_deps EQUAL 0)
      meds_add_test(restart_exact
                    COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_restart_exact.py
                    ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 MEDS_MAIN=$<TARGET_FILE:meds_main> LABELS integration)
   else()
      message(STATUS "MEDS: restart_exact skipped (${Python3_EXECUTABLE} lacks numpy, netCDF4, pytest or tomli)")
   endif()
endif()

#----- Every test stops at 10 minutes. The longest takes about 1.5 in a Debug build, so a test that #
#      hangs, or slows by an order of magnitude, fails instead of holding the suite. The long runs, #
#      registered after this, set their own.                                                         #
get_property(meds_all_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${meds_all_tests} PROPERTIES TIMEOUT 600)

#==========================================================================================#
#  long (opt-in, -DMEDS_TEST_LONG=ON): example 05's two Palo Verde runs through meds_main, the  #
#  runs to make before a physics or tissue-energy change merges. Each runs a copy of the        #
#  example's config in its own build-tree directory and passes when the run completes and its  #
#  slow ledger closes; checks/check_long_run.cmake explains them. The inventory (85 s on the     #
#  example's 12 threads, ifx dev) reads the ED_ERA5land archive; the regrowth (31 min) reads    #
#  Palo Verde's cell cut from it, 2003-2022: MEDS_TEST_PV_FORCING names a cut file, and         #
#  without one the test cuts it into its directory once (under 2 min).                          #
#==========================================================================================#
if (MEDS_TEST_LONG)
   set(MEDS_ERA5LAND_ARCHIVE "$ENV{ERA5LAND_ARCHIVE}" CACHE PATH
       "The processed ED_ERA5land archive the long tests read (default: the environment's ERA5LAND_ARCHIVE)")
   set(MEDS_TEST_PV_FORCING "" CACHE FILEPATH
       "Palo Verde's forcing file, 2003-2022, as run_palo_verde.py cuts it; empty: the test cuts it from the archive")
   if (NOT IS_DIRECTORY "${MEDS_ERA5LAND_ARCHIVE}")
      message(WARNING "MEDS: MEDS_TEST_LONG is on, but MEDS_ERA5LAND_ARCHIVE (\"${MEDS_ERA5LAND_ARCHIVE}\") is not a "
                      "folder: the long tests will stop at their start. Configure with "
                      "-DMEDS_ERA5LAND_ARCHIVE=<the ED_ERA5land folder>, or set ERA5LAND_ARCHIVE.")
   endif()
   foreach(case inventory regrowth)
      meds_add_test(palo_verde_${case}
                    COMMAND ${CMAKE_COMMAND} -DMEDS_MAIN=$<TARGET_FILE:meds_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                            -DWORK_DIR=${PROJECT_BINARY_DIR}/long/palo_verde_${case} -DCASE=${case}
                            -DARCHIVE=${MEDS_ERA5LAND_ARCHIVE} -DFORCING=${MEDS_TEST_PV_FORCING}
                            -DPYTHON=${Python3_EXECUTABLE}
                            -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_long_run.cmake
                    LABELS long)
   endforeach()
   #----- The runs' own thread count ([run].n_threads), so a parallel ctest books it; a generous    #
   #      limit, a Debug build being many times slower than the times above.                       #
   set_tests_properties(palo_verde_inventory PROPERTIES TIMEOUT 3600 PROCESSORS 12)
   set_tests_properties(palo_verde_regrowth PROPERTIES TIMEOUT 21600 PROCESSORS 12)
endif()
