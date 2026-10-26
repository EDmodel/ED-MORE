# SPDX-License-Identifier: Apache-2.0
#----- The shipped ed_io_config.toml is what edmore --dump-io-config writes. ------------------     #
#
# The example lists every output variable an io_config may switch, with its default files, and
# docs/science/diagnostics.md points users at it. It is generated (edmore --dump-io-config) and
# also tracked, so it can fall behind the output registry without anything noticing: a variable
# it omits keeps its default files when a user copies the example and switches everything off.
# This regenerates it from the two reference configs in a directory of the build tree and asserts
# that the tracked file is the same, byte for byte. When the registry changes, regenerate the
# tracked file as the failure message says.
#
# Usage: cmake -DEDMORE_MAIN=<edmore> -DSOURCE_DIR=<repo> -DWORK_DIR=<build-tree dir>
#              -P check_io_config_example.cmake
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.20)

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(COPY "${SOURCE_DIR}/ed_config_main.toml" "${SOURCE_DIR}/ed_config_pft.toml" DESTINATION "${WORK_DIR}")
execute_process(COMMAND "${EDMORE_MAIN}" --dump-io-config ed_config_main.toml WORKING_DIRECTORY "${WORK_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
if(NOT rc EQUAL 0 OR NOT EXISTS "${WORK_DIR}/ed_io_config.toml")
   message(FATAL_ERROR "edmore --dump-io-config did not write ed_io_config.toml\n--- output:\n${out}${err}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${WORK_DIR}/ed_io_config.toml"
                        "${SOURCE_DIR}/ed_io_config.toml" RESULT_VARIABLE differ)
if(NOT differ EQUAL 0)
   message(FATAL_ERROR "ed_io_config.toml is not what edmore --dump-io-config writes from the reference "
                       "configs. Regenerate it:\n  cp ${WORK_DIR}/ed_io_config.toml ${SOURCE_DIR}/ed_io_config.toml")
endif()
message(STATUS "ed_io_config.toml is what edmore --dump-io-config writes")
