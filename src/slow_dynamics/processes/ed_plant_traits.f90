! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_plant_traits -- changes to a cohort's leaf traits within its lifetime.
!   * light_plastic_traits       -- the light-acclimated targets of SLA, Vcmax25, Rd25 and
!       the leaf turnover rate (docs/science/plant_traits.md eq. (1)).
!   * update_plastic_trait       -- one step of a live trait toward its target, as fast as
!       leaves are replaced, or a jump to it on a census start (eq. (2)).
!   * acclimate_leaf_photo_table -- a site's leaf table at its growth temperature (Kattge &
!       Knorr 2007; docs/science/leaf_gas_exchange.md §2, "Thermal acclimation").
! The light routines are elemental; the loop over the cohorts and their trait state is
! ed_slow_state's.
!---~---
module ed_plant_traits
   use ed_kinds,        only : wp, ik
   use ed_constants,    only : safe_exp
   use ed_pft_params,   only : leaf_photo_table_t
   use ed_therm_lib,    only : kattge_knorr_entropy, kattge_knorr_jv_ratio
   implicit none
   private

   public :: light_plastic_traits, update_plastic_trait
   public :: acclimate_leaf_photo_table

contains

   !---~---
   !   Light-acclimated targets for one cohort: each top-of-canopy trait times
   ! exp(kplastic * cum_lai_above), cum_lai_above floored at 0. kplastic < 0 (Vcmax, Rd,
   ! leaf turnover) lowers a trait in shade, kplastic > 0 (SLA) raises it (eq. (1)).
   !---~---
   elemental pure subroutine light_plastic_traits(cum_lai_above, sla_toc, vcmax25_toc, rd25_toc,  &
            turnover_toc, kplastic_sla, kplastic_vcmax, kplastic_rd, kplastic_turnover,              &
            sla_target, vcmax25_target, rd25_target, turnover_target)
      real(wp), intent(in)  :: cum_lai_above                                 ! [m2/m2] LAI above the cohort
      real(wp), intent(in)  :: sla_toc, vcmax25_toc, rd25_toc, turnover_toc  ! top-of-canopy trait values
      real(wp), intent(in)  :: kplastic_sla, kplastic_vcmax, kplastic_rd, kplastic_turnover  ! light-response slopes
      real(wp), intent(out) :: sla_target, vcmax25_target, rd25_target, turnover_target
      sla_target      = sla_toc      * light_gradient(kplastic_sla,      cum_lai_above)
      vcmax25_target  = vcmax25_toc  * light_gradient(kplastic_vcmax,      cum_lai_above)
      rd25_target     = rd25_toc     * light_gradient(kplastic_rd,       cum_lai_above)
      turnover_target = turnover_toc * light_gradient(kplastic_turnover, cum_lai_above)
   end subroutine light_plastic_traits

   !--- The exponential canopy-light trait gradient, guarded against exp overflow.
   elemental pure function light_gradient(kplastic, cum_lai_above) result(f)
      real(wp), intent(in) :: kplastic, cum_lai_above
      real(wp)             :: f
      f = safe_exp(kplastic * max(cum_lai_above, 0.0_wp))
   end function light_gradient

   !---~---
   !   One update of a live trait toward its target by the fraction of leaves replaced this
   ! step, f = 1 - exp(-turnover_rate dt_yr), at the cohort's current turnover rate [1/yr];
   ! with `instant` (a census start, no trait history) it jumps to the target (eq. (2)).
   !---~---
   elemental pure function update_plastic_trait(current, target, turnover_rate, dt_yr, instant) result(updated)
      real(wp), intent(in) :: current, target, turnover_rate, dt_yr
      logical,  intent(in) :: instant
      real(wp)             :: updated, f
      if (instant) then
         updated = target
      else
         f       = 1.0_wp - safe_exp(-turnover_rate * dt_yr)
         updated = current + f * (target - current)
      end if
   end function update_plastic_trait

   !---~---
   !   A leaf table's entropy terms and Jmax:Vcmax ratio at a growth temperature, each
   ! a - b T with T in degC (leaf_gas_exchange.md §2, "Thermal acclimation", eq. (4)).
   ! Every other parameter stays the PFT table's, so a repeat call changes nothing.
   !---~---
   pure subroutine acclimate_leaf_photo_table(t_growth, ds_vcmax_a, ds_vcmax_b, ds_jmax_a, ds_jmax_b,   &
                                              jv_a, jv_b, table)
      real(wp),                 intent(in)    :: t_growth     ! [K] running-mean growth temperature
      real(wp),                 intent(in)    :: ds_vcmax_a, ds_vcmax_b, ds_jmax_a, ds_jmax_b, jv_a, jv_b
      type(leaf_photo_table_t), intent(inout) :: table
      integer(ik) :: ipft
      do ipft = 1_ik, table%n_pft
         table%pft(ipft)%ds_vcmax = kattge_knorr_entropy(ds_vcmax_a, ds_vcmax_b, t_growth)
         table%pft(ipft)%ds_jmax  = kattge_knorr_entropy(ds_jmax_a,  ds_jmax_b,  t_growth)
         table%jmax_vcmax_ratio(ipft) = kattge_knorr_jv_ratio(jv_a, jv_b, t_growth)
      end do
   end subroutine acclimate_leaf_photo_table

end module ed_plant_traits
