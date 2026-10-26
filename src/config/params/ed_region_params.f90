! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_region_params -- [run].mode and the [region] block, as a config leaf                  !
! (docs/configuration.md, "Regional runs"). A site run simulates one location from [site]; !
! a region run simulates every selected forcing cell of a box as its own polygon, at its   !
! cell's centre.                                                                           !
!==========================================================================================!
module ed_region_params
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: region_params_t, RUN_MODE_SITE, RUN_MODE_REGION

   !----- [run].mode. ---------------------------------------------------------------------!
   integer(ik), parameter :: RUN_MODE_SITE   = 1_ik   !< one location, the [site] block
   integer(ik), parameter :: RUN_MODE_REGION = 2_ik   !< every selected cell of [region].box_nwse

   type :: region_params_t
      real(wp)    :: box_nwse(4) = 0.0_wp           !< [deg] N, W, S, E; W > E crosses 180 degrees
      real(wp)    :: land_fraction_min = 0.5_wp     !< [-] cells of less static land are skipped
      !----- Cell mask: a cell runs where cell_mask_var >= cell_mask_min ('' = all). ------!
      character(len=512) :: cell_mask     = ''
      character(len=64)  :: cell_mask_var = 'mask'
      real(wp)           :: cell_mask_min = 0.5_wp
   end type region_params_t

end module ed_region_params
