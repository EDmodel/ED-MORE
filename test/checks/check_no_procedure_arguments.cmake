# SPDX-License-Identifier: Apache-2.0
# Fails if any Fortran source under SRC_DIR declares a procedure argument or a procedure pointer,
# that is, a "procedure(...)" declaration outside a comment.
#
# Why: under ifx, a routine that hands one of its contained functions to another routine (a
# residual to a root finder, an integrand to a quadrature) allocates a lock-guarded record on every
# call. Called per cohort and per step, that lock makes more than four threads slower than four.
# Solvers in MEDS take their problem as an explicit record instead.
#
# Usage: cmake -DSRC_DIR=<src> -P check_no_procedure_arguments.cmake
file(GLOB_RECURSE sources "${SRC_DIR}/*.f90")
set(found "")
foreach(source IN LISTS sources)
   file(STRINGS "${source}" lines REGEX "[Pp][Rr][Oo][Cc][Ee][Dd][Uu][Rr][Ee] *\\(")
   foreach(line IN LISTS lines)
      string(REGEX REPLACE "!.*" "" code "${line}")
      if (code MATCHES "[Pp][Rr][Oo][Cc][Ee][Dd][Uu][Rr][Ee] *\\(")
         string(APPEND found "\n  ${source}: ${line}")
      endif()
   endforeach()
endforeach()
if (found)
   message(FATAL_ERROR "A procedure argument or pointer is declared (see the comment at the top of "
                       "test/checks/check_no_procedure_arguments.cmake):${found}")
endif()
message(STATUS "no procedure arguments under ${SRC_DIR}")
