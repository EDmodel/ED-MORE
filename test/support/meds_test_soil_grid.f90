! SPDX-License-Identifier: Apache-2.0
!----- Test fixture: geometrically thickening soil layers for the column tests                  !
!      (depth D in n layers, layer k's bottom at D (e^(g k / n) - 1) / (e^g - 1)). The model's   !
!      own layers are SOIL_LAYER_TEMPLATE cut at the depth; these tests pass this grid as an     !
!      explicit list of bottoms, because the numbers they check are computed on this grid. -----!
module meds_test_soil_grid
   use meds_kinds, only : wp, ik
   implicit none
   private
   public :: geometric_layer_bottoms

contains

   pure subroutine geometric_layer_bottoms(depth, n_layer, growth, bottoms)
      real(wp),    intent(in)  :: depth, growth
      integer(ik), intent(in)  :: n_layer
      real(wp),    intent(out) :: bottoms(:)
      real(wp)    :: denom
      integer(ik) :: k
      bottoms = 0.0_wp
      denom   = exp(growth) - 1.0_wp
      do k = 1_ik, n_layer
         bottoms(k) = depth * (exp(growth * real(k, wp) / real(n_layer, wp)) - 1.0_wp) / denom
      end do
   end subroutine geometric_layer_bottoms

end module meds_test_soil_grid
