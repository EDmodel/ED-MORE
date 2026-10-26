! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_trait_dynamics_params -- the [trait_dynamics] block: light acclimation of the leaf  !
! traits (docs/science/plant_traits.md).                                                   !
!==========================================================================================!
module meds_trait_dynamics_params
   implicit none
   private

   public :: trait_dynamics_params_t

   type :: trait_dynamics_params_t
      logical :: trait_plasticity_on       = .false.   !< light acclimation of the leaf traits
      logical :: fineroot_turnover_plastic = .false.   !< fine roots turn over with their leaves
   end type trait_dynamics_params_t

end module meds_trait_dynamics_params
