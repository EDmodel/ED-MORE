# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# MEDS -- Modular Ecosystem Dynamics Simulator / Modern Ecosystem Demography Simulator        #
#         (standalone demographic core).                                                     #
#                                                                                          #
# CMake handles Fortran module-dependency ordering automatically (no hand-maintained object #
# lists, no repeated-build hack).                                                            #
#                                                                                          #
# THE PARALLEL PATH IS THE CPU, not the GPU. Patch-axis OpenMP threading on the host is what  #
# ships and what is byte-identical at any thread count. GPU offload was evaluated and is not   #
# viable as scoped: the offload build ran 1.4x slower than the CPU (49.8 s against 36.2 s),    #
# one kernel at 0.4% occupancy, so its build option is retired                                 #
# (docs/dev_plans/archive/MEDS_GPU_EVALUATION.md).                                             #
#==========================================================================================#
cmake_minimum_required(VERSION 3.20)
project(MEDS VERSION 0.4.0 LANGUAGES Fortran
        DESCRIPTION "Modular Ecosystem Dynamics Simulator (Modern Ecosystem Demography Simulator)")

#      A configure that names no build type gets the dev kind (docs/building.md, "Build kinds"):
#      Release, the project's -O2, code for any x86-64 CPU. The debug kind is asked for by name.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Release (dev, profile, production) or Debug" FORCE)
endif()

option(MEDS_BUILD_TESTS "Build the CTest test suite" ON)
# The long runs (ctest label `long`): example 05's Palo Verde inventory and its 150-year regrowth,
# which read the ED_ERA5land archive (MEDS_ERA5LAND_ARCHIVE). OFF by default: they take over half an hour.
option(MEDS_TEST_LONG "Register the long tests (label long): the Palo Verde inventory and 150-year regrowth" OFF)
# netCDF output + I/O are ALWAYS compiled: the standard (and only) build needs the netCDF C library.
# Point CMake at it with -DCMAKE_PREFIX_PATH=<prefix> (find_package(netCDF CONFIG REQUIRED), below).

# Optional Python C-API shared libraries (ctypes). OFF by default -> normal builds are foreign-
# call-free. When ON, ALL static libs are built position-independent so they link into the .so.
option(MEDS_BUILD_PYLIB "Build the C-API shared libraries for Python (ctypes)" OFF)
if(MEDS_BUILD_PYLIB)
   set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# All Fortran .mod files in one directory, visible to every target.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

#----- Per-compiler, per-config flags. ----------------------------------------------------#
#
#  THE AUTOMATIC-ARRAY STACK CEILING (opt-in flag below).
#
#  The fast loop carries ~740 B of automatic (stack) arrays per COHORT, live simultaneously across
#  four nested frames (fast_step -> build_frozen -> canopy_leaf_gas_exchange /
#  aero_bottom_to_top). On an 8 MB stack that is a hard SIGSEGV ceiling near 10-12k cohorts, and in
#  a threaded run the ceiling follows each thread's stack, OMP_STACKSIZE, instead.
#
#  PREFERRED FIX IS A RUNTIME ONE -- raise the stack rather than change codegen:
#      ulimit -s unlimited          (or e.g. `ulimit -s 65536` for ~90k cohorts)
#      export OMP_STACKSIZE=512M    (a threaded run)
#
#  -DMEDS_HEAP_ARRAYS=ON moves large automatic arrays to the heap instead. It is OFF by default
#  because ifx's -heap-arrays SEGFAULTS test_output (its registry section) at the exit of
#  apply_variable_override in a Debug build -- a compiler codegen interaction, not a MEDS bug (the
#  routine only assigns plain integer components, and the same code is clean without the flag).
#  Turn it on only if the runtime fix is unavailable, and re-run the Debug suite when you do.
#-----------------------------------------------------------------------------------------#
#
#  AUTOMATIC (STACK) LOCALS ARE MANDATORY WHEN THE PATCH LOOP IS THREADED.
#
#  It is the hardest trap here to find, so it is written down. Intel Fortran defaults to
#  `-auto-scalar`: local SCALARS go on the stack, but local ARRAYS and DERIVED TYPES go in STATIC
#  storage. Static locals are shared by every thread. The fast loop calls code in meds_fast_kernels,
#  meds_share and meds_state, and every routine with a local array or a local
#  derived type (e.g. `rad_forcing_t rf` in apply_rt_forcing -> alloc_rad_forcing) then has all
#  threads writing the same object. It does NOT crash: it produces plausible, silently
#  thread-count-dependent answers -- precisely what the byte-identity-at-any-thread-count requirement
#  exists to catch. ThreadSanitizer finds it; a green test suite does not.
#
#  `-qopenmp`/`-fopenmp` implies automatic locals, but only for the target it is applied to, so
#  flagging just the target that owns meds_fast_dynamics is NOT enough -- the whole call tree needs
#  it. MEDS_OPENMP therefore adds the per-compiler "all locals on the stack" flag to EVERY target.
#  Corollary: more goes on the stack, so the OMP_STACKSIZE note above applies with more force.
#-----------------------------------------------------------------------------------------#
option(MEDS_HEAP_ARRAYS "Put large automatic arrays on the heap (see the note above)" OFF)
option(MEDS_OPENMP "Compile the patch (one polygon) and polygon (several) parallelism; serial if the compiler has none" ON)
#      ON BY DEFAULT. The thread count is a RUNTIME setting ([run].n_threads, default 1), so an OpenMP
#      build still runs serially until a config asks for threads, and its output is byte-identical at
#      any thread count. Only an OpenMP build lets test_fast_driver's 4-thread check test anything: in a
#      serial build the directives are comments and it passes trivially. The price is the automatic-
#      locals flag on every target, measured in docs/building.md "Parallel builds". -DMEDS_OPENMP=OFF
#      builds serial; a compiler with no Fortran OpenMP falls back to serial with a warning.
set(MEDS_OPENMP_ACTIVE OFF)
if(MEDS_OPENMP)
   find_package(OpenMP COMPONENTS Fortran)
   if(OpenMP_Fortran_FOUND)
      set(MEDS_OPENMP_ACTIVE ON)
   else()
      message(WARNING "MEDS: ${CMAKE_Fortran_COMPILER_ID} offers no Fortran OpenMP, so this build is "
                      "SERIAL and [run].n_threads has no effect. Configure with -DMEDS_OPENMP=OFF "
                      "to ask for that and silence this warning.")
   endif()
endif()
#      Any build that can run host threads needs automatic locals.
set(MEDS_AUTO_LOCALS ${MEDS_OPENMP_ACTIVE})

#----- THE BUILD KINDS (docs/building.md, "Build kinds"): debug and dev are the Debug and Release   #
#      build types; production is Release plus the two options below; CMakePresets.json names all   #
#      four. Both options are offered for ifx and gfortran only.                                    #
#                                                                                                   #
#      MEDS_IPO optimises across files. Each file compiles to an intermediate form at the project's #
#      -O2 and the link generates the code for the whole program at once, so the link names its     #
#      own level, -O3, rather than inheriting whatever CMAKE_Fortran_FLAGS_RELEASE holds. The level #
#      is ifx's: gfortran keeps each file's -O2 through the link, and its binary is the same with   #
#      the link at -O2 or -O3. Every link redoes the whole program: with ifx an edit to one file    #
#      takes ~35 s to rebuild instead of ~4 s, which is why everyday work builds without it.        #
#      MEDS_CPU_TARGET names the CPU family to generate code for, so nobody writes a compiler's     #
#      flag by hand: empty is any x86-64 CPU, `cascadelake` the AVX-512 Xeons, `host` the machine   #
#      that builds. Such a binary may not start on an older CPU.                                    #
option(MEDS_IPO "Optimise across files, generating the code at the link at -O3 (the production kind)" OFF)
set(MEDS_CPU_TARGET "" CACHE STRING "CPU family to generate code for: empty (any x86-64), cascadelake or host")
set_property(CACHE MEDS_CPU_TARGET PROPERTY STRINGS "" cascadelake host)
if(MEDS_IPO OR NOT MEDS_CPU_TARGET STREQUAL "")
   if(NOT (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
      message(FATAL_ERROR "MEDS: MEDS_IPO and MEDS_CPU_TARGET are offered for ifx and gfortran, not "
                          "${CMAKE_Fortran_COMPILER_ID}. Configure without them.")
   endif()
endif()
#      Each CPU target's flag: ifx's, then gfortran's.
set(meds_cpu_cascadelake -xCORE-AVX512 -march=cascadelake)
set(meds_cpu_host        -xHost        -march=native)
set(MEDS_CPU_FLAG "")
if(NOT MEDS_CPU_TARGET STREQUAL "")
   if(NOT MEDS_CPU_TARGET MATCHES "^(cascadelake|host)$")
      message(FATAL_ERROR "MEDS: MEDS_CPU_TARGET=${MEDS_CPU_TARGET} is not a CPU target MEDS knows. Use "
                          "cascadelake or host, or leave it empty for a binary that runs on any x86-64 CPU.")
   endif()
   if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      list(GET meds_cpu_${MEDS_CPU_TARGET} 1 MEDS_CPU_FLAG)
   else()
      list(GET meds_cpu_${MEDS_CPU_TARGET} 0 MEDS_CPU_FLAG)
   endif()
   message(STATUS "MEDS: code for the ${MEDS_CPU_TARGET} CPU family (${MEDS_CPU_FLAG})")
endif()
if(MEDS_IPO)
   include(CheckIPOSupported)
   check_ipo_supported(RESULT meds_ipo_supported OUTPUT meds_ipo_output LANGUAGES Fortran)
   if(NOT meds_ipo_supported)
      message(FATAL_ERROR "MEDS: MEDS_IPO=ON, but ${CMAKE_Fortran_COMPILER_ID} "
                          "${CMAKE_Fortran_COMPILER_VERSION} cannot optimise across files here:\n"
                          "${meds_ipo_output}")
   endif()
   #      Set before any target is made: each target takes its INTERPROCEDURAL_OPTIMIZATION from it.
   set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
   #      A link line carries the build type's flags but none of a target's compile options, so the
   #      project's -O2 never reaches it. The level goes last in the Release flags, where it wins over
   #      whatever they held; files still compile at -O2, which follows them on every compile line.
   string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O3")
   message(STATUS "MEDS: interprocedural optimisation, the code generated at the link at -O3")
endif()

function(meds_fortran_flags tgt)
   set(id ${CMAKE_Fortran_COMPILER_ID})
   if(id STREQUAL "IntelLLVM" OR id STREQUAL "Intel")
      target_compile_options(${tgt} PRIVATE -stand f18 -warn all -diag-disable=7712
         $<$<BOOL:${MEDS_HEAP_ARRAYS}>:-heap-arrays;256>
         $<$<BOOL:${MEDS_AUTO_LOCALS}>:-auto>
         $<$<CONFIG:Debug>:-O0;-g;-check;all;-fpe0;-traceback>
         $<$<CONFIG:Release>:-O2>)
   elseif(id STREQUAL "NVHPC" OR id STREQUAL "PGI")
      target_compile_options(${tgt} PRIVATE -Minform=warn
         $<$<BOOL:${MEDS_HEAP_ARRAYS}>:-Mnostack_arrays>
         $<$<BOOL:${MEDS_AUTO_LOCALS}>:-Mrecursive>
         $<$<CONFIG:Debug>:-O0;-g;-Mbounds;-Ktrap=fp>
         $<$<CONFIG:Release>:-O2>)
   elseif(id STREQUAL "GNU")
      target_compile_options(${tgt} PRIVATE -std=f2018 -Wall -Wno-unused-dummy-argument
         $<$<BOOL:${MEDS_HEAP_ARRAYS}>:-fmax-stack-var-size=262144>
         $<$<BOOL:${MEDS_AUTO_LOCALS}>:-frecursive>
         $<$<CONFIG:Debug>:-O0;-g;-fcheck=all;-ffpe-trap=invalid,zero,overflow;-fbacktrace>
         $<$<CONFIG:Release>:-O2>)
   endif()
   #  The build kind's CPU target goes on the link too: under IPO the link generates the code.
   if(NOT MEDS_CPU_FLAG STREQUAL "")
      target_compile_options(${tgt} PRIVATE ${MEDS_CPU_FLAG})
      target_link_options(${tgt} PRIVATE ${MEDS_CPU_FLAG})
   endif()
endfunction()

#----- FOUNDATION (src/share/): the kinds, the constants, the numerics, the calendar and the sun's #
#      position, and the budget checks at the top, and the stateless constitutive laws in         #
#      functions/, one library per domain: allometry, plant hydraulics, soil hydraulics,           #
#      thermodynamics (moist air and its lapse, with the temperature response of enzyme rates),    #
#      canopy optics and the radiation above the canopy.                                           #
#                                                                                                  #
#      One folder, one library, and it links NOTHING. A module belongs in src/share/ only if it    #
#      `use`s nothing outside src/share/ -- no model state, no configuration, no external          #
#      library. That is what lets every layer above, and the tests that link one narrow library,   #
#      reach it freely, and it is the test that keeps this folder from becoming the place a file   #
#      goes when no other place fits. The temperature response sits here so leaf, respiration and  #
#      any other tissue share one code path without a plant-domain -> plant-domain library edge;   #
#      meds_allometry_lib because it is a geometric CONSTRAINT that state, demography, plant and   #
#      config all read.                                                                            #
file(GLOB SHARE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/share/*.f90
                                          ${CMAKE_CURRENT_SOURCE_DIR}/src/share/functions/*.f90)
add_library(meds_share STATIC ${SHARE_SOURCES})
meds_fortran_flags(meds_share)

#----- CONFIGURATION: the meds_config_t aggregator, the loader that fills it from two TOML files, #
#      the checks it passes and the minimal TOML reader; in params/, one leaf per block (the PFT   #
#      table, the soil, forcing, output, the fast loop, snow and aerodynamics, leaf,               #
#      biogeochemistry, maps, region, sites). `toml` and `config_io` live here,                    #
#      not in src/io/:                                                                             #
#      they are configuration INPUT, not I/O of model results, and they carry no external          #
#      dependency, unlike the netCDF output layer. Links meds_share only; nothing here imports     #
#      model state or a kernel.                                                                    #
file(GLOB CONFIG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/config/*.f90
                                           ${CMAKE_CURRENT_SOURCE_DIR}/src/config/params/*.f90)
add_library(meds_config STATIC ${CONFIG_SOURCES})
#      The keys a config may hold are the ones the two reference configs list (meds_config_check): #
#      written into an include at configure time, and re-written when either reference changes.   #
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/meds_config_keys.cmake)
meds_write_config_keys(${CMAKE_CURRENT_SOURCE_DIR}/meds_config_main.toml
                       ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_pft.toml
                       ${CMAKE_CURRENT_BINARY_DIR}/generated/meds_config_keys.inc)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
             ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_main.toml ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_pft.toml)
target_include_directories(meds_config PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_link_libraries(meds_config PUBLIC meds_share)
meds_fortran_flags(meds_config)

#----- STATE, one LAYER in six files by model hierarchy (src/README.md, rule 3). --------------#
#                                                                                            #
#  meds_cohort_state_vars  -- the flat cohort SoA and the cohort lockstep.                    #
#  meds_patch_state_vars   -- ONE patch's stores (canopy air, snow, soil water and heat,      #
#                             soil carbon), their blends, and the patch SoA with its CSR map. #
#  meds_polygon_state_vars -- polygon_t: its site, patch and cohort blocks, and its numbers.  #
#  meds_site_state_vars    -- the site block: each site's area and its slice of the patches.  #
#  meds_diag_vars          -- the diagnostic accumulators, per cohort, patch and polygon.     #
#                                                                                            #
#  The kernels take only the per-patch store types (from meds_patch_state_vars), never        #
#  polygon_t, which keeps them OpenMP-`target` device-eligible; ctest kernels_patch_only       #
#  checks it, because the shared module directory cannot. The parameters that describe the    #
#  stores are config (meds_soil_params): a parameter is not state.                            #
file(GLOB STATE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/state/*.f90)
add_library(meds_state STATIC ${STATE_SOURCES})
target_link_libraries(meds_state PUBLIC meds_share meds_config)
meds_fortran_flags(meds_state)



#==========================================================================================#
#  PROCESSES, TIMESCALE-FIRST (src/README.md, rule 2).                                        #
#                                                                                            #
#  MEDS is an operator-split fast/slow model -- that is its dominant organizing axis -- and    #
#  the tree says so. Each half holds its process kernels (processes/) and the loop that drives   #
#  them (driver/); the fast half also holds the integrator machinery (integrator/).               #
#==========================================================================================#

#----- FAST kernels (sub-daily): the process modules in processes/, one per process: canopy      #
#      radiation, aerodynamics and the canopy air space; leaf gas exchange, respiration (the       #
#      plants' maintenance and the soil's heterotrophic), the plant's energy, the canopy's water   #
#      and plant hydraulics; soil water, soil heat and the surface (snow and bare ground).         #
#      State-free apart from the patch store types; builds standalone.                             #
#                                                                                                  #
#      It links meds_state and meds_config only, and imports nothing from slow_dynamics/.         #
file(GLOB FAST_KERNEL_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/processes/*.f90)
add_library(meds_fast_kernels STATIC ${FAST_KERNEL_SOURCES})
target_link_libraries(meds_fast_kernels PUBLIC meds_state meds_config)
meds_fortran_flags(meds_fast_kernels)

#----- SLOW kernels (daily / annual): the process modules in processes/, one per science page:  #
#      leaf phenology, carbon allocation, the plastic traits, the soil carbon pools' daily step   #
#      and the demographic laws. Stateless and `pure`-where-possible; links the patch store      #
#      types it reads/writes, the config leaves, and the fast kernels for the soil's respiration  #
#      chemistry (meds_respiration), which the fast loop respires and the daily step integrates,  #
#      so the two debit the same CENTURY matrix (docs/science/soil_carbon.md §7). Builds          #
#      standalone.                                                                                 #
file(GLOB SLOW_KERNEL_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/processes/*.f90)
add_library(meds_slow_kernels STATIC ${SLOW_KERNEL_SOURCES})
target_link_libraries(meds_slow_kernels PUBLIC meds_state meds_config meds_fast_kernels)
meds_fortran_flags(meds_slow_kernels)

#----- Diagnostic output (src/io/output), its netCDF-free half: every module of the folder but the  #
#      writer -- the types, the registry, the file set's setup, the extract, the reductions and the  #
#      per-tier integrators. Links meds_state only, so the stepper edge that ticks the               #
#      integrators pulls no C dependency (src/README.md, "The library graph"). The writer is         #
#      meds_output_stream (below). The stages: docs/science/diagnostics.md §1.                      #
file(GLOB OUTPUT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/output/*.f90)
list(REMOVE_ITEM OUTPUT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/io/output/meds_output_stream.f90)
add_library(meds_output STATIC ${OUTPUT_SOURCES})
target_link_libraries(meds_output PUBLIC meds_state)
meds_fortran_flags(meds_output)





#----- netCDF C bindings, in their own minimal target so BOTH the output layer                  #
#      (meds_init_restart, state serialization) AND the forcing reader (meds_met) link the netCDF C     #
#      API without the rest of the output layer. netCDF is a hard dependency;                          #
#      its CMake config pulls in HDF5 -> FindThreads, which needs the C language enabled.          #
enable_language(C)
find_package(netCDF CONFIG REQUIRED)
add_library(meds_netcdf_c STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_c.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_write.f90)
target_link_libraries(meds_netcdf_c PUBLIC meds_share netCDF::netcdf)
meds_fortran_flags(meds_netcdf_c)

#----- Meteorological FORCING (src/io/met): the reader, its samples at a site and at a patch's   #
#      canopy-air top, the grid and the runtime types; the laws are share/'s. Links meds_config   #
#      + the netCDF bindings ONLY (no model state): a prescribed external driver, low in the DAG. #
#      The science: docs/science/forcing.md.                                                     #
file(GLOB FORCING_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/met/*.f90)
add_library(meds_met STATIC ${FORCING_SOURCES})
target_link_libraries(meds_met PUBLIC meds_config meds_netcdf_c)
meds_fortran_flags(meds_met)
#----- Static per-cell MAPS (src/io/maps): the one module that opens a map -- the cell       #
#      mask, and each cell's soil and height-allometry inputs.                                     #
#      Links meds_met for the archive grid it checks maps against and the cell list it reads   #
#      them at; the driver links it.                                                                 #
file(GLOB MAPS_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/maps/*.f90)
add_library(meds_maps STATIC ${MAPS_SOURCES})
target_link_libraries(meds_maps PUBLIC meds_met meds_config meds_netcdf_c)
meds_fortran_flags(meds_maps)


#----- The FAST loop: the integrator machinery (integrator/ -- the records, the frozen build, the #
#      water network's step, the ESDIRK2 and RK45 steps, the march and the update) and the        #
#      driver that walks a slow step in dt_fast sub-steps over the patch axis (driver/).          #
#                                                                                                  #
#      The frozen build's process rates sit in integrator/, not driver/, because the step calls   #
#      them once per dt_fast, so filing them under driver/ would make a cycle.                    #
file(GLOB FAST_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/integrator/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/driver/*.f90)
add_library(meds_fast STATIC ${FAST_SOURCES})
target_link_libraries(meds_fast PUBLIC meds_fast_kernels meds_state meds_met meds_output)
meds_fortran_flags(meds_fast)

#----- The SLOW loop (driver/): the orchestrator, the per-process reads and writes of the stand, #
#      the update that writes it, the slow rows and the site ledger. It is the one place a        #
#      demographic rate meets its application (src/README.md rule 4).                            #
file(GLOB SLOW_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/driver/*.f90)
add_library(meds_slow STATIC ${SLOW_SOURCES})
target_link_libraries(meds_slow PUBLIC meds_slow_kernels)
meds_fortran_flags(meds_slow)

#----- The initial state (src/io/init): bare ground, a census or a restart file, and           #
#      init_polygon, which starts a polygon from the one [init].init_mode names. It links       #
#      meds_slow for a census stand's restructuring and trait acclimation, and the netCDF        #
#      bindings for the restart file, which the checkpoints write.                               #
file(GLOB INIT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/init/*.f90)
add_library(meds_init STATIC ${INIT_SOURCES})
target_link_libraries(meds_init PUBLIC meds_slow meds_config meds_netcdf_c)
meds_fortran_flags(meds_init)

#----- src/main: the master stepper (the cadence owner -- one call advances the coupled model by  #
#      one slow step, sub-stepping the fast loop inside it) and the program that drives it. The    #
#      stepper is a library because the integration tests drive it directly; meds_main.f90 is the   #
#      PROGRAM and is the only thing in the executable.                                             #
add_library(meds_stepper STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_stepper.f90)
target_link_libraries(meds_stepper PUBLIC meds_fast meds_slow)
meds_fortran_flags(meds_stepper)

#----- The whole coupled model, as one name for everything above it to link. -----------------#
add_library(meds_model INTERFACE)
target_link_libraries(meds_model INTERFACE meds_stepper meds_fast meds_slow meds_init)

#----- HOST OpenMP over the fast-loop PATCH axis (src/fast_dynamics/driver/meds_fast_dynamics).  #
#      ON by default (the option, beside meds_fortran_flags, says why). Applied to meds_fast -- the #
#      target that owns the fast loop -- and PUBLIC so meds_main and the tests inherit the link     #
#      flags. find_package, run with the option, resolves the per-compiler flag (-qopenmp         #
#      IntelLLVM, -fopenmp GNU, -mp NVHPC). The thread count is a RUNTIME knob ([run].n_threads,     #
#      default 1), not a build-time one.                                                            #
#      The option itself is declared beside meds_fortran_flags(), which is where its OTHER (and more #
#      important) half lives: the automatic-locals flag on every target.                             #
if(MEDS_OPENMP_ACTIVE)
   target_link_libraries(meds_fast PUBLIC OpenMP::OpenMP_Fortran)
   message(STATUS "MEDS: fast-loop OpenMP -> ${OpenMP_Fortran_FLAGS} (patch axis; [run].n_threads)")
else()
   message(STATUS "MEDS: serial build; [run].n_threads has no effect")
endif()

#----- Diagnostic output's writer (src/io/output/meds_output_stream.f90): the netCDF files. Links  #
#      the netCDF bindings (meds_netcdf_c, above). netCDF is always compiled.                        #
#      (The C-binding layer passes array temporaries to the netCDF C `(*)` dummies --                #
#      correct for interop, but Intel Debug `-check all` emits a harmless arg_temp_created warning   #
#      per call; a Release build, the default, is clean.)                                            #
add_library(meds_output_stream STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/output/meds_output_stream.f90)
target_link_libraries(meds_output_stream PUBLIC meds_netcdf_c meds_output)
meds_fortran_flags(meds_output_stream)
message(STATUS "MEDS: netCDF output (always compiled)")

#----- The coupled model as OPEN / STEP / FINALIZE objects: one polygon's state and step           #
#      (meds_polygon) and the run of one polygon or many (meds_region). A LIBRARY rather than part  #
#      of the program because it has two callers: the executable below, and the C-API shim         #
#      src/c_api/meds_c_api_run.f90 that drives the full coupled model from Python. It sits above   #
#      meds_model + the IO layer -- the one place that owns configuration, the calendar, the output #
#      files and the closing conservation reports together.                                         #
add_library(meds_run STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_polygon.f90
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_region.f90)
target_link_libraries(meds_run PUBLIC meds_model meds_maps meds_output_stream)
meds_fortran_flags(meds_run)

#----- The main entry point: parse the command line, drive the model, exit. ----------------#
add_executable(meds_main ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_main.f90)
target_link_libraries(meds_main PRIVATE meds_run)

#----- The C-API shims (src/c_api): the bind(c) procedures Python calls, compiled ONCE, into an #
#      object library that both of their users link: libmeds.so (below) and the C-API tests.    #
#      Every source is compiled by exactly one target because a Fortran compile also writes the  #
#      source's .mod files into the one shared module directory: two targets compiling the same  #
#      source would write the same .mod at once, and a parallel build would fail at random. It   #
#      links meds_run, the whole model, because the run shim holds a whole run.               #
#                                                                                                #
#      Plain GLOB, no exclusions: every src/c_api/*.f90 belongs in the one library.              #
if(MEDS_BUILD_PYLIB OR MEDS_BUILD_TESTS)
   file(GLOB C_API_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/c_api/*.f90)
   add_library(meds_c_api OBJECT ${C_API_SOURCES})
   target_link_libraries(meds_c_api PUBLIC meds_run)
   meds_fortran_flags(meds_c_api)
endif()

#----- THE Python C-API shared library (opt-in): ONE libmeds.so for the whole package.        #
#                                                                                              #
#      Per-subsystem independence lives at the CALL level (verbs, not libraries), not the LINK #
#      level: `meds.plant`, `meds.canopy`, `meds.demography` and `meds.model` are sub-packages #
#      over one backend. A kernels-only library beside a site-driver one would force every     #
#      domain kernel library to link meds_share ONLY, and the column state types into          #
#      src/share/. The accepted cost: the library links meds_model + IO, so even the leaf      #
#      shim carries the whole model, and the wheel declares netCDF4 a runtime dependency       #
#      rather than vendoring it.                                                               #
#                                                                                              #
#      It compiles nothing of its own: it links the shims' objects (meds_c_api, above), and    #
#      meds_fortran_flags gives its link the build kind's CPU target.                          #
if(MEDS_BUILD_PYLIB)
   add_library(meds_py SHARED $<TARGET_OBJECTS:meds_c_api>)
   set_target_properties(meds_py PROPERTIES OUTPUT_NAME meds LINKER_LANGUAGE Fortran)
   target_link_libraries(meds_py PRIVATE meds_run)
   meds_fortran_flags(meds_py)
   #----- The same numbers as meds_main, in any process. An Intel build calls `exp`, `sin`, `pow`  #
   #      and the rest by their C names, and meds_main finds them in Intel's math library. A       #
   #      library loaded into Python would find glibc's first, because Python has already loaded  #
   #      it, so without this a run through `meds.model` differs from the same run of meds_main    #
   #      in the last digits (BCI: up to 7e-4 W m-2 in LE over ten days). Linking Intel's          #
   #      runtime into the library and keeping its names inside it makes every call land where     #
   #      meds_main's does.                                                                        #
   #      gfortran builds use glibc's math in both, so they need nothing here.                    #
   if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
      target_link_options(meds_py PRIVATE -static-intel LINKER:--exclude-libs,ALL)
   endif()
   #----- RPATH so the INSTALLED library resolves its own dependencies.                         #
   #                                                                                            #
   #      libmeds.so links netCDF, HDF5 and the Fortran runtime. Declaring netCDF4 rather than     #
   #      vendoring the ~67 MB closure is necessary, but on its own NOT enough: the pip netCDF4    #
   #      wheel bundles its libnetcdf under a mangled soname inside its own package, so it cannot  #
   #      satisfy a `libnetcdf.so.19` NEEDED by a different .so. Freshly installed into a bare     #
   #      venv, the library without this RPATH fails to load on ten unresolved sonames.            #
   #                                                                                            #
   #      USE_LINK_PATH writes the directories of the actually-linked externals into the install  #
   #      RPATH, so the wheel finds them with no LD_LIBRARY_PATH and no `source setvars.sh`.      #
   #      $ORIGIN comes first so a co-bundled copy would still win.                                #
   #                                                                                            #
   #      THIS MAKES THE WHEEL MACHINE-LOCAL: those paths are the build machine's conda prefix    #
   #      and oneAPI install. That is the right trade for `pip install python/` from a source      #
   #      checkout, which is how MEDS is used; it is NOT a redistributable manylinux wheel, and    #
   #      building one would mean auditwheel plus the vendoring declined above.                    #
   #      USE_LINK_PATH alone is NOT enough: it captures the directories of libraries CMake was  #
   #      told to link (netCDF, HDF5), but the compiler's own shared runtime (the OpenMP runtime  #
   #      in an Intel build, which stays shared; libgfortran in a GNU one) is added implicitly by #
   #      the compiler driver, so it is absent from the RPATH and the library still fails to       #
   #      dlopen outside a `source setvars.sh` shell. CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES is   #
   #      exactly that set. (`ldd` will NOT show this: run from a sourced shell it inherits        #
   #      LD_LIBRARY_PATH and reports everything resolved. Check with `env -u LD_LIBRARY_PATH      #
   #      ldd`, or by importing.)                                                                  #
   set_target_properties(meds_py PROPERTIES
                         INSTALL_RPATH "$ORIGIN;${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
                         INSTALL_RPATH_USE_LINK_PATH TRUE)
   #----- INSTALL it INTO the package directory, so a wheel carries meds/libmeds.so and         #
   #      meds._libmeds finds it as a sibling of __init__.py -- the first entry in its search    #
   #      order, needing no env var. The `pylib` component is what python/pyproject.toml         #
   #      installs; nothing else in this project is part of the wheel.                            #
   install(TARGETS meds_py LIBRARY DESTINATION meds COMPONENT pylib)
   message(STATUS "MEDS: Python C-API ENABLED -> libmeds.so (meds.plant, meds.canopy, meds.demography, meds.model)")
endif()
meds_fortran_flags(meds_main)

#----- Tests: the suite, its support libraries and the checks, in test/ (test/CMakeLists.txt).  #
#      enable_testing() stays here, so ctest runs the suite from the build tree's top.             #
if(MEDS_BUILD_TESTS)
   enable_testing()
   add_subdirectory(test)
endif()
