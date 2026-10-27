! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_kinds -- the numeric kinds every module uses. The working precision wp is set here,
! on one line.
!---~---
module ed_kinds
   use, intrinsic :: iso_fortran_env, only : real64, real32, int32, int64
   implicit none
   public

   integer, parameter :: wp = real64   !--- working precision of the model's state and arithmetic
   integer, parameter :: sp = real32   !--- single precision: the forcing archive's buffers, as it stores them
   integer, parameter :: ik = int32    !--- indices and counts
   integer, parameter :: lk = int64    !--- wide integer: the calendar's day numbers and spans in seconds
end module ed_kinds
