# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# The ED-MORE test suite (docs/building.md, "Running the tests"): one test per src module,   #
# in folders that mirror src/, its support libraries and the repository's checks.            #
#                                                                                            #
# Every test is registered by ed_add_test (below) with exactly one label:                    #
#   unit         a module's test, and edmore's refusals: seconds each                        #
#   integration  whole runs: a region and the fast loop's driver                             #
#   lint         the repository's rules                                                      #
#   python       pytest suites                                                               #
#   long         opt-in runs too long for every build                                        #
#                                                                                            #
# Paths: CMAKE_CURRENT_SOURCE_DIR is this folder; PROJECT_SOURCE_DIR the repository and      #
# PROJECT_BINARY_DIR the build tree's top, where every test runs unless it names another      #
# directory and where the work directories are.                                              #
#==========================================================================================#

#----- ed_add_test(<name> [SOURCES <file>...] [LINK <library>...] [COMMAND <command>...]        #
#                    [ARGS <argument>...] [WORKDIR <dir>] [ENVIRONMENT <VAR=value>...]          #
#                    LABELS <label>)                                                            #
#      The test <name>. With SOURCES it builds the program test_<name> from them, linked to     #
#      LINK, and runs that; otherwise it runs COMMAND. ARGS follow on the command line. It runs  #
#      in WORKDIR, made at configure time (the build tree's top when none is given), with        #
#      ENVIRONMENT. It takes exactly one of the labels above. The TIMEOUT every test gets is set #
#      at the end of this file.                                                                  #
function(ed_add_test name)
   cmake_parse_arguments(PARSE_ARGV 1 arg "" "WORKDIR" "SOURCES;LINK;COMMAND;ARGS;ENVIRONMENT;LABELS")
   if(arg_UNPARSED_ARGUMENTS)
      message(FATAL_ERROR "ed_add_test(${name}): unknown arguments: ${arg_UNPARSED_ARGUMENTS}")
   endif()
   list(LENGTH arg_LABELS n_labels)
   if(NOT n_labels EQUAL 1 OR NOT arg_LABELS MATCHES "^(unit|integration|lint|python|long)$")
      message(FATAL_ERROR "ed_add_test(${name}): give exactly one label of unit, integration, lint, python "
                          "and long")
   endif()
   if(arg_SOURCES)
      add_executable(test_${name} ${arg_SOURCES})
      target_link_libraries(test_${name} PRIVATE ${arg_LINK})
      ed_fortran_flags(test_${name})
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
#----- The assertions, the soil grid and the derived configs (ed_test_run), depending on nothing   #
#      but a kind. Kept apart from ed_testsupport because many tests deliberately link ONE        #
#      narrow library each to keep that layer standalone-buildable, and they must be able to      #
#      assert without linking the whole model.                                                    #
add_library(ed_testassert   STATIC support/ed_test_assert.f90 support/ed_test_soil_grid.f90
                                   support/ed_test_run.f90)
target_link_libraries(ed_testassert PUBLIC ed_share)
ed_fortran_flags(ed_testassert)
#----- The in-code config fixture; re-exports the assertions so most tests need one `use`. ------#
add_library(ed_testsupport STATIC support/ed_test_support.f90)
target_link_libraries(ed_testsupport PUBLIC ed_config ed_testassert)
ed_fortran_flags(ed_testsupport)
#----- The synthetic ED_ERA5land archive, shared by the reader, the maps and the region tests. ---#
add_library(ed_testarchive STATIC support/ed_test_archive.f90)
target_link_libraries(ed_testarchive PUBLIC ed_met)
ed_fortran_flags(ed_testarchive)
#----- The fast column the integrator tests run on: assembled as a run assembles one, its forcing, #
#      its frozen record and the ledger check.                                                      #
add_library(ed_testfastfixture STATIC support/ed_test_fast_fixture.f90)
target_link_libraries(ed_testfastfixture PUBLIC ed_fast ed_testsupport)
ed_fortran_flags(ed_testfastfixture)

#==========================================================================================#
#  share/: the libraries of src/share, one test each, named after its module. Each links      #
#  ed_share alone; test_optics_lib also the config layer, whose ed_pft_derive turns a         #
#  PFT's leaf-angle moments into the LIDF's Beta parameters.                                  #
#==========================================================================================#
foreach(t numerics therm_lib plant_hydr_lib soil_hydr_lib optics_lib radiation_lib)
   ed_add_test(${t} SOURCES share/test_${t}.f90 LINK ed_share ed_testassert LABELS unit)
endforeach()
target_link_libraries(test_optics_lib PRIVATE ed_config)
#----- ed_allometry_lib: the size laws, and the shipped [allometry] block equal to the library's    #
#      initializers (the Python library's config-free values), which are Chave 2014 in carbon and   #
#      ED2's leaf fit over C2B. ed_config only for the TOML reader.                                  #
ed_add_test(allometry_lib SOURCES share/test_allometry_lib.f90 LINK ed_config ed_testassert
            ARGS ${PROJECT_SOURCE_DIR}/ed_config_pft.toml LABELS unit)

#==========================================================================================#
#  config/: the configuration layer (src/config).                                             #
#==========================================================================================#
#----- nearest_key, and the misspelt, moved and unknown keys a run reports. ----------------------#
ed_add_test(config_keys SOURCES config/test_config_keys.f90 LINK ed_config ed_testsupport LABELS unit)
#----- The config blocks reach their records: the [soil]/[fast]/[snow]/[aerodynamics] loaders, and #
#      the [soil] column through build_soil_profile.                                                #
ed_add_test(params_config SOURCES config/test_params_config.f90 LINK ed_config ed_testassert LABELS unit)
#----- The PFT file's loader and the per-PFT tables derived from it: the shipped file through the  #
#      loader, and the optics and hydraulics that differentiate PFTs. It runs from the SOURCE dir,  #
#      where ed_config_main.toml names the PFT file, and writes its CSV in a directory of its own   #
#      in the build tree.                                                                            #
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/work/pft_config)
ed_add_test(pft_config SOURCES config/test_pft_config.f90 LINK ed_model ed_testsupport
            ARGS ${PROJECT_BINARY_DIR}/work/pft_config WORKDIR ${PROJECT_SOURCE_DIR} LABELS unit)

#==========================================================================================#
#  state/: the state layer's containers (src/state).                                          #
#==========================================================================================#
ed_add_test(containers SOURCES state/test_containers.f90 LINK ed_testsupport ed_model LABELS unit)

#==========================================================================================#
#  fast/: the fast loop (src/fast_dynamics).                                                  #
#==========================================================================================#
#----- The process modules (processes/), one test each, named after its module. Each links the    #
#      kernel library alone, so the kernels build and test standalone; the leaf test also the      #
#      in-code config fixture, whose leaf table it solves on.                                      #
foreach(t canopy_aerodynamics canopy_air_space canopy_radiation canopy_water leaf_gas_exchange plant_energy
          plant_hydraulics respiration soil_heat soil_water surface)
   ed_add_test(${t} SOURCES fast/test_${t}.f90 LINK ed_fast_kernels ed_testassert LABELS unit)
endforeach()
target_link_libraries(test_leaf_gas_exchange PRIVATE ed_testsupport)
#----- The integrator (integrator/), one test per concern, named after it: its building blocks,   #
#      its two schemes and the march that drives them, and one fast step end to end. Each links    #
#      the fast loop and the shared fast fixture (support/).                                       #
foreach(t fast_integ_kernels fast_integ_schemes fast_integ_step)
   ed_add_test(${t} SOURCES fast/test_${t}.f90 LINK ed_testfastfixture LABELS unit)
endforeach()
#----- The soil-plant water network on synthetic stands and soils, and the soil column alone. -----#
ed_add_test(water_network SOURCES fast/test_water_network.f90 LINK ed_fast ed_testassert LABELS unit)
#----- The fast loop's driver (driver/) and the stepper that runs it, on a forcing file the test   #
#      writes in a directory of its own in the build tree.                                         #
ed_add_test(fast_driver SOURCES fast/test_fast_driver.f90 LINK ed_testsupport ed_model
            WORKDIR ${PROJECT_BINARY_DIR}/work/fast_driver LABELS integration)

#==========================================================================================#
#  slow/: the slow loop (src/slow_dynamics).                                                  #
#==========================================================================================#
#----- The daily plant kernels (allocation, traits, demography) and the soil carbon pools' daily   #
#      step, which debits them through the decomposition matrix of ed_respiration (in              #
#      ed_fast_kernels): the kernel libraries, state-free.                                         #
foreach(t plant_carbon_allocation plant_traits demography soil_carbon)
   ed_add_test(${t} SOURCES slow/test_${t}.f90 LINK ed_fast_kernels ed_slow_kernels ed_testassert
               LABELS unit)
endforeach()
#----- The phenology kernel and its C-API shim, so the shim compiles under ctest. ----------------#
ed_add_test(plant_phenology SOURCES slow/test_plant_phenology.f90
            LINK ed_slow_kernels ed_c_api ed_testassert LABELS unit)
#----- The driver: the stand's operators, the slow-loop site conservation ledger and the step. ----#
foreach(t slow_update slow_step)
   ed_add_test(${t} SOURCES slow/test_${t}.f90 LINK ed_testsupport ed_model LABELS unit)
endforeach()
ed_add_test(slow_ledger SOURCES slow/test_slow_ledger.f90 LINK ed_model ed_testsupport LABELS unit)

#==========================================================================================#
#  io/: the initial state, the forcing, the maps, the output and the state file (src/io).     #
#==========================================================================================#
#----- How a run's polygon starts: the census reader and every reason it refuses a file, and the  #
#      initial soil state.                                                                         #
ed_add_test(init SOURCES io/test_init.f90 LINK ed_testsupport ed_model LABELS unit)
#----- An ED-MORE forcing file: its contract, the reader, and the forcing at a site. The met tests #
#      write their files in a directory of their own in the build tree.                            #
ed_add_test(met_file SOURCES io/test_met_file.f90 LINK ed_met ed_testassert
            ARGS ${PROJECT_SOURCE_DIR}/data/co2/co2_cmip7_global_annual_1000-2022.txt
            WORKDIR ${PROJECT_BINARY_DIR}/work/met_file LABELS unit)
#----- The ED_ERA5land reader, on the synthetic archive. -----------------------------------------#
ed_add_test(met_era5land SOURCES io/test_met_era5land.f90 LINK ed_met ed_testassert ed_testarchive
            WORKDIR ${PROJECT_BINARY_DIR}/work/met_era5land LABELS unit)
#----- The static-map reader and the cell mask, on the synthetic archive's grid. ------------------#
ed_add_test(maps SOURCES io/test_maps.f90 LINK ed_maps ed_testassert ed_testarchive
            ARGS ${PROJECT_BINARY_DIR}/work/maps LABELS unit)
#----- The diagnostic output's netCDF-free half (src/io/output): the registry and the selection,   #
#      the reduction, the extract and the per-tier integration, one section per module. It links   #
#      ed_output rather than the serializer.                                                        #
ed_add_test(output SOURCES io/test_output.f90 LINK ed_output ed_testsupport LABELS unit)
#----- The serializer's round trip: writes and reopens netCDF, so it links the io layer. ---------#
ed_add_test(output_roundtrip SOURCES io/test_output_roundtrip.f90 LINK ed_output_stream ed_testsupport
            LABELS unit)
#----- The state file: written, read back and written again, every variable (links the IO layer). #
ed_add_test(state_file SOURCES io/test_state_file.f90 LINK ed_model ed_output_stream ed_testsupport
            WORKDIR ${PROJECT_BINARY_DIR}/work/state_file LABELS unit)

#==========================================================================================#
#  driver/: the coupled model's runs (src/main).                                              #
#==========================================================================================#
#----- A region against its site runs, on the synthetic archive. It runs from driver/region, whose  #
#      site config it derives its configs from; it writes everything under the work directory in   #
#      the build tree.                                                                               #
ed_add_test(region SOURCES driver/test_region.f90 LINK ed_run ed_testassert ed_testarchive
            ARGS ${PROJECT_BINARY_DIR}/region_test_tmp $<TARGET_FILE:ed_main>
            WORKDIR ${CMAKE_CURRENT_SOURCE_DIR}/driver/region LABELS integration)

#==========================================================================================#
#  c_api/: the C-API shims (src/c_api). Linking every shim INTO a ctest target is the point: a #
#  shim built only by the optional -DED_BUILD_PYLIB=ON shared library breaks unseen -- a       #
#  component inserted mid-type in leaf_photo_args_t, or a signature changed under the          #
#  demography shim, leaves the whole suite green on both back ends while the C API is broken.  #
#  Linking it here makes that class of break a normal build failure. Needs no .so and no       #
#  Python: each test calls the bind(c) procedures directly, and links ed_c_api, the shims'     #
#  one compile, so a default build compiles and links every shim in src/c_api/.                #
#==========================================================================================#
#----- Runs from the SOURCE dir: the canopy shim opens the reference configuration and the PFT     #
#      file it names, relative to it, and writes nothing.                                           #
ed_add_test(c_api_leaf SOURCES c_api/test_c_api_leaf.f90 LINK ed_testsupport ed_c_api
            WORKDIR ${PROJECT_SOURCE_DIR} LABELS unit)

#==========================================================================================#
#  checks/: the CMake scripts that check the source tree and edmore as a program.             #
#==========================================================================================#
#----- No procedure is passed as an argument anywhere in src/: the ifx cost of doing so serialises #
#      threads. The script explains the rule it checks.                                            #
ed_add_test(no_procedure_arguments
            COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${PROJECT_SOURCE_DIR}/src
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_no_procedure_arguments.cmake
              LABELS lint)
#----- The two reference configs list exactly the keys the loader reads: a run stops on any key its #
#      reference does not list. The script explains the rule it checks.                              #
ed_add_test(config_keys_listed
            COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${PROJECT_SOURCE_DIR}/src
                      -DMAIN_TOML=${PROJECT_SOURCE_DIR}/ed_config_main.toml
                      -DPFT_TOML=${PROJECT_SOURCE_DIR}/ed_config_pft.toml
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_config_keys.cmake
              LABELS lint)
#----- The shipped per-variable override EXAMPLE is a generated file that is also tracked, so it can #
#      drift silently from the registry: it must be what edmore --dump-io-config writes from the     #
#      reference configs. The script explains the rule.                                              #
ed_add_test(io_config_example
            COMMAND ${CMAKE_COMMAND} -DEDMORE_MAIN=$<TARGET_FILE:ed_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                      -DWORK_DIR=${PROJECT_BINARY_DIR}/work/io_config_example
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_io_config_example.cmake
              LABELS lint)
#----- edmore refuses a run it cannot start and names the cause: unknown, moved or misplaced        #
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
   ed_add_test(${case}_refused
               COMMAND ${CMAKE_COMMAND} -DEDMORE_MAIN=$<TARGET_FILE:ed_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                         -DWORK_DIR=${PROJECT_BINARY_DIR}/refused/${case} -DCASE=${case}
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_main_refuses.cmake
               LABELS unit)
endforeach()
#----- A step that fails stops edmore after the run is finalized, which closes the files. The       #
#      script explains the case.                                                                      #
ed_add_test(main_failure_finalizes
            COMMAND ${CMAKE_COMMAND} -DEDMORE_MAIN=$<TARGET_FILE:ed_main> -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
                      -DWORK_DIR=${PROJECT_BINARY_DIR}/main_failure
                      -P ${CMAKE_CURRENT_SOURCE_DIR}/checks/check_main_failure.cmake
            LABELS unit)

#==========================================================================================#
#  The repository's checks (scripts/lint), one test each and named after the check, when   #
#  CMake finds a Python: `ctest -L lint` runs every rule that a check enforces. A check    #
#  added to scripts/lint becomes a test at the next build.                                 #
#==========================================================================================#
find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
   file(GLOB ed_lint_checks CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/scripts/lint/check_*.py)
   foreach(script IN LISTS ed_lint_checks)
      get_filename_component(check ${script} NAME_WE)
      string(REGEX REPLACE "^check_" "" check ${check})
      ed_add_test(${check} COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/lint/run_checks.py
                  --root ${PROJECT_SOURCE_DIR} ${check} LABELS lint)
   endforeach()
else()
   message(STATUS "ED-MORE: the repository checks are not registered as tests (CMake found no Python)")
endif()

#==========================================================================================#
#  python/: the package's own tests, the restart test that pytest drives, and the wheel's  #
#  copies of LICENSE and NOTICE. The pytest suites need a Python with numpy, netCDF4 and   #
#  pytest (and tomli below 3.11); without one they are left out, and configure says so.    #
#==========================================================================================#
#----- The wheel ships COPIES of LICENSE and NOTICE, because python/pyproject.toml cannot           #
#      list license files above its own directory. A tracked copy drifts silently, so the           #
#      build compares them.                                                                         #
foreach(f LICENSE NOTICE)
   string(TOLOWER ${f} f_lower)
   ed_add_test(${f_lower}_copy COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
               ${PROJECT_SOURCE_DIR}/${f} ${PROJECT_SOURCE_DIR}/python/${f} LABELS lint)
endforeach()
if (Python3_Interpreter_FOUND)
   execute_process(COMMAND ${Python3_EXECUTABLE} -c
                   "import sys, numpy, netCDF4, pytest; sys.version_info >= (3, 11) or __import__('tomli')"
                   RESULT_VARIABLE ed_python_test_deps OUTPUT_QUIET ERROR_QUIET)
endif()
if (Python3_Interpreter_FOUND AND ed_python_test_deps EQUAL 0)
   #----- The package's tests (python/tests): the leaf, phenology and canopy kernels,               #
   #      edmore.config, and edmore.model against edmore, bit for bit. They load libedmore.so.      #
   if (ED_BUILD_PYLIB)
      ed_add_test(python_api
                  COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider ${PROJECT_SOURCE_DIR}/python/tests
                  ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=${PROJECT_SOURCE_DIR}/python
                              EDMORE_LIB=$<TARGET_FILE:ed_py> EDMORE_MAIN=$<TARGET_FILE:ed_main>
                  LABELS python)
   endif()
   #----- A restart continues exactly, and re-acclimates its traits on request                      #
   #      (python/test_restart_exact.py: edmore on the constant reference climate, and on a         #
   #      forcing file edmore.forcing writes). A whole run, so it is labelled integration,          #
   #      though pytest drives it.                                                                  #
   ed_add_test(restart_exact
               COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                       ${CMAKE_CURRENT_SOURCE_DIR}/python/test_restart_exact.py
               ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=${PROJECT_SOURCE_DIR}/python
                           EDMORE_MAIN=$<TARGET_FILE:ed_main>
               LABELS integration)
else()
   message(STATUS "ED-MORE: python_api and restart_exact are not registered (they need a Python with "
                  "numpy, netCDF4 and pytest)")
endif()

#----- Every test stops at 10 minutes. The longest takes about 1.5 in a Debug build, so a test that #
#      hangs, or slows by an order of magnitude, fails instead of holding the suite.                #
get_property(ed_all_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${ed_all_tests} PROPERTIES TIMEOUT 600)
