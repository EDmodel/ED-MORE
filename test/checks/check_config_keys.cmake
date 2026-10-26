# SPDX-License-Identifier: Apache-2.0
# Fails unless the two reference configs list exactly the keys the loader reads.
#
# Why: a run stops on any key its reference does not list (meds_config_check), so a key the loader
# reads but the references omit could never be set, and a key the references list but nothing
# reads would be accepted and do nothing. Both directions are checked, and the two files may not
# share a section, so a key cannot be listed in the file the loader does not read it from.
#
# Usage: cmake -DSRC_DIR=<src> -DMAIN_TOML=<main> -DPFT_TOML=<pft> -P check_config_keys.cmake
cmake_minimum_required(VERSION 3.20)
include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/meds_config_keys.cmake)

meds_toml_keys("${MAIN_TOML}" main_keys)
meds_toml_keys("${PFT_TOML}" pft_keys)

# The keys the loader reads: the literal second argument of every toml_* / req_* / opt_* reader
# call outside a comment (a key read through toml_has alone only refuses or gates on it).
set(read_keys "")
file(GLOB sources "${SRC_DIR}/config/*.f90" "${SRC_DIR}/main/*.f90")
foreach(source IN LISTS sources)
   file(STRINGS "${source}" lines REGEX "(toml_|req_|opt_)[a-z_]+ *\\( *[a-z]+ *, *'")
   foreach(line IN LISTS lines)
      if(line MATCHES "^[ \t]*!")
         continue()
      endif()
      string(REGEX MATCHALL "(toml_|req_|opt_)[a-z_]+ *\\( *[a-z]+ *, *'[^']+'" calls "${line}")
      foreach(call IN LISTS calls)
         if(call MATCHES "^toml_has")
            continue()
         endif()
         string(REGEX REPLACE ".*'([^']+)'$" "\\1" key "${call}")
         if(key MATCHES "\\.")
            list(APPEND read_keys "${key}")
         endif()
      endforeach()
   endforeach()
endforeach()
list(REMOVE_DUPLICATES read_keys)

set(listed ${main_keys} ${pft_keys})
set(problems "")
foreach(key IN LISTS read_keys)
   if(NOT key IN_LIST listed)
      string(APPEND problems "\n  read but not listed in either reference: ${key}")
   endif()
endforeach()
foreach(key IN LISTS listed)
   if(NOT key IN_LIST read_keys)
      string(APPEND problems "\n  listed in a reference but never read: ${key}")
   endif()
endforeach()

# The two files may not share a section.
foreach(key IN LISTS main_keys)
   string(REGEX REPLACE "\\..*" "" section "${key}")
   list(APPEND main_sections "${section}")
endforeach()
foreach(key IN LISTS pft_keys)
   string(REGEX REPLACE "\\..*" "" section "${key}")
   if(section IN_LIST main_sections)
      string(APPEND problems "\n  section [${section}] is in both meds_config_main.toml and meds_config_pft.toml")
   endif()
endforeach()

if(problems)
   message(FATAL_ERROR "The reference configs and the loader disagree (see the comment at the top "
                       "of test/checks/check_config_keys.cmake):${problems}")
endif()
list(LENGTH read_keys n)
message(STATUS "the references list all ${n} keys the loader reads, and only those")
