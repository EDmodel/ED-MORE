! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_state_defaults -- the values a state field holds before it has a real one, in one
! place: the two "no value yet" markers and their tests, a new cohort's phenology at
! birth, and the canopy-air placeholders a patch holds until the loops first set them. An
! empty amount (0 kg, a zeroed sum) or a physical starting value (a liquid fraction of 1)
! is the real value, not a marker, and is written where it is used.
!---~---
module ed_state_defaults
   use ed_kinds, only : wp
   implicit none
   private

   public :: is_unset_positive, is_unset_negative

   !---~---
   !   "No value yet", by the sign of the field marked: a field that is never negative
   ! once set (a temperature, an hour count) holds UNSET_POSITIVE until then, and one
   ! never positive (a water potential) UNSET_NEGATIVE. The tests are is_unset_positive
   ! and is_unset_negative: a value on the wrong side of zero has not been set.
   !---~---
   real(wp), parameter, public :: UNSET_POSITIVE = -1.0_wp
   real(wp), parameter, public :: UNSET_NEGATIVE =  1.0_wp

   !--- The daily maximum's start: below any leaf potential, so max() replaces it.
   real(wp), parameter, public :: DMAX_PSI_LEAF_ACCUM_RESET = -1.0e30_wp  !--- [MPa]

   !--- Cohort: a newborn's phenology -- flushing, not shedding.
   real(wp), parameter, public :: PHENO_FLUSH_INIT = 1.0_wp
   real(wp), parameter, public :: PHENO_SHED_INIT  = 0.0_wp

   !--- Patch: the canopy air until the fast loop seeds it and the slow loop sizes it.
   real(wp), parameter, public :: CAN_CO2_UNSEEDED = 400.0_wp   !--- [umol/mol]
   real(wp), parameter, public :: CAN_DEPTH_UNSET  =  20.0_wp   !--- [m]

contains

   elemental logical function is_unset_positive(x)
      real(wp), intent(in) :: x
      is_unset_positive = x < 0.0_wp
   end function is_unset_positive

   elemental logical function is_unset_negative(x)
      real(wp), intent(in) :: x
      is_unset_negative = x > 0.0_wp
   end function is_unset_negative

end module ed_state_defaults
