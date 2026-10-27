! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_demography_params -- the [demography] block: one switch per operator of the stand,
! and their settings: recruitment, cohort fusion and fission, the cull, the treefall gaps,
! patch fusion and termination (docs/science/demography.md; docs/configuration.md,
! "Cohort and patch fusion").
!---~---
module ed_demography_params
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: demography_params_t

   type :: demography_params_t
      !---~---
      !   The switches. Each operator runs only with demography_on; the cull and patch
      ! termination, which keep the blocks bounded, have no switch of their own.
      !---~---
      logical     :: demography_on        = .true.   !--- .false.: the stand's structure is frozen
      logical     :: recruitment_on       = .true.   !--- monthly founding from the recruit pools
      logical     :: cohort_fusion_on     = .true.   !--- monthly, and yearly after patch fusion
      logical     :: cohort_fission_on    = .true.   !--- monthly
      logical     :: patch_disturbance_on = .true.   !--- yearly treefall gaps
      logical     :: patch_fusion_on      = .true.   !--- yearly
      !--- Recruitment.
      real(wp)    :: min_recruit_size                !--- [plant/m2] the pool that founds a cohort
      !--- Cohort fusion and fission, and the cull.
      integer(ik) :: max_cohort                      !--- cohorts in the most crowded patch, a target
      integer(ik) :: n_cohort_fusion_iter
      real(wp)    :: cohort_size_tol_min, cohort_size_tol_max
      real(wp)    :: cohort_lai_cap                  !--- [m2/m2]
      real(wp)    :: split_eps
      real(wp)    :: min_cohort_agb                  !--- [kgC/m2]
      real(wp)    :: negligible_nplant               !--- [plant/m2]
      !--- Treefall gaps.
      real(wp)    :: patch_disturbance_rate          !--- [1/yr]
      real(wp)    :: disturbance_survive_height      !--- [m]
      !--- Patch fusion and termination.
      integer(ik) :: max_patch                       !--- patches per site, a target
      integer(ik) :: n_patch_fusion_iter, n_height_layers
      real(wp)    :: patch_light_tol, patch_light_tol_max
      real(wp)    :: light_ext_coeff_fusefiss = 0.5_wp   !--- [-]
      real(wp)    :: patch_light_maxdev_factor
      real(wp)    :: min_patch_area
      !--- The operators' conservation tolerance, and the growth memory mortality reads.
      real(wp)    :: conservation_tol
      real(wp)    :: growth_memory_days              !--- [day]
      !--- DERIVED (derive_config).
      real(wp)    :: cohort_size_tol_mult            !--- geometric step of the cohort tolerance
      real(wp)    :: patch_light_tol_mult            !--- geometric step of the patch tolerance
      real(wp), allocatable :: height_edges(:)       !--- [m] interior edges of the height layers
   end type demography_params_t

end module ed_demography_params
