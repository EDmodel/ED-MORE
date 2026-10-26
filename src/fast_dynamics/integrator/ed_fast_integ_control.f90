! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_fast_integ_control -- the fast loop's error control, shared by both schemes so a      !
! target accuracy means the same on each: the grouped error norm state_wrms_grouped and    !
! the step factor step_control_factor (the integral controller). The march state mixes     !
! units, so each entry's error is scaled by its group's (rtol, atol) (the seven groups of  !
! ed_fast_integ_types; the films and the pond carry none). The run's error control is      !
! built from [fast] (build_error_control), with its budget stops; default_error_control    !
! serves the tests (docs/science/numerical_scheme.md eqs. (3)-(4)). water_scales sets the  !
! scale of each of the soil-plant water network's stores, against which its sub-step       !
! tolerances are read (docs/science/soil_plant_water.md §3).                               !
!==========================================================================================!
module ed_fast_integ_control
   use ed_kinds,         only : wp, ik
   use ed_constants,     only : tiny_num, rho_h2o
   use ed_numerics,      only : adaptive_step_update
   use ed_config,        only : ed_config_t
   use ed_fast_integ_types,  only : march_state_t, tol_set_t, error_control_t, GRP_THETA, fast_static_state_t, POND, &
                                state_to_array, state_entry_rules
   use ed_soil_params,   only : soil_profile_t
   use ed_pft_params,    only : hydraulics_table_t
   implicit none
   private

   public :: state_wrms_grouped, step_control_factor
   public :: build_error_control, build_tol_set
   public :: water_scales
   public :: default_error_control   ! used only by the tests

contains

   !---------------------------------------------------------------------------------------!
   ! The single tolerance set of the fast loop, the per-group tolerances of its error norm !
   ! (docs/science/numerical_scheme.md eq. (3)): [fast].rtol for every group, each group's !
   ! default atol, and [fast].water_atol for GRP_THETA, which the water network's own       !
   ! sub-step control takes too (docs/science/soil_plant_water.md §3). ESDIRK2 advances    !
   ! plant water outside its stages, so GRP_LEAF_W and GRP_WOOD_W act only under RK45.     !
   !---------------------------------------------------------------------------------------!
   pure function build_tol_set(cfg) result(tols)
      type(ed_config_t),   intent(in) :: cfg
      type(tol_set_t)                 :: tols
      tols%rtol             = cfg%rtol
      tols%atol(GRP_THETA)  = cfg%water_atol
      !----- Every atol scales with atol_scale: rtol alone stops tightening once atol dominates the   !
      !      WRMS denominator atol + rtol*|y| (docs/science/numerical_scheme.md §6). ----------------!
      tols%atol = tols%atol * cfg%atol_scale
   end function build_tol_set

   !----- The run's error control from [fast]: the tolerance set, the step controller and the budget !
   !      stops (the soil heat's threshold, a budget closure rather than a tolerance, scaled by       !
   !      atol_scale directly). ------------------------------------------------------------------!
   pure function build_error_control(cfg) result(ec)
      type(ed_config_t),   intent(in) :: cfg
      type(error_control_t)           :: ec
      ec%tols       = build_tol_set(cfg)
      ec%debug_error = cfg%debug_error
      ec%soil_energy_budget_atol = cfg%soil_energy_budget_atol * cfg%atol_scale
   end function build_error_control

   !----- The default error control: one rtol for every group, the default atols, the I controller, !
   !      level L1. -------------------------------------------------------------------------------!
   pure function default_error_control(rtol) result(ec)
      real(wp), intent(in)  :: rtol
      type(error_control_t) :: ec
      ec%tols = default_tol_set(rtol)
   end function default_error_control

   !----- The default tolerance set: one `rtol` for every group and the per-field default atols    !
   !      (tol_set_t's initial values, ed_fast_integ_types). ------------------------------------------  !
   pure function default_tol_set(rtol) result(tols)
      real(wp), intent(in) :: rtol
      type(tol_set_t)      :: tols
      tols%rtol = rtol              ! the same rtol for every group
      ! tols%atol keeps its default (the per-field constants)
   end function default_tol_set

   !---------------------------------------------------------------------------------------!
   ! Grouped WRMS norm of (a - b), each entry scaled by its group's atol + rtol*|y_ref|:        !
   ! docs/science/numerical_scheme.md eq. (3), over the entries the march integrates            !
   ! (state_entry_rules): the canopy air and the soil heat, on both schemes. The soil and plant  !
   ! water move at the water network's fixed rate, so their terms would be exact zeros that      !
   ! dilute the norm by a factor that grows with the cohort count; the network's own sub-step    !
   ! control holds them to the same GRP_THETA tolerances.                                       !
   !---------------------------------------------------------------------------------------!
   pure function state_wrms_grouped(a, b, y_ref, n, nsl, tols) result(err)
      type(march_state_t), intent(in) :: a, b, y_ref
      integer(ik),          intent(in) :: n, nsl
      type(tol_set_t),      intent(in) :: tols
      real(wp)    :: err, s
      real(wp)    :: va(5_ik + 2_ik*nsl + 4_ik*n), vb(5_ik + 2_ik*nsl + 4_ik*n), vr(5_ik + 2_ik*nsl + 4_ik*n)
      integer(ik) :: j, grp, cnt
      logical     :: counted
      call state_to_array(a, n, nsl, va)
      call state_to_array(b, n, nsl, vb)
      call state_to_array(y_ref, n, nsl, vr)
      s = 0.0_wp
      cnt = 0_ik
      do j = 1_ik, size(va, kind=ik)            ! the layout's order is the order the terms are added
         call state_entry_rules(j, n, nsl, grp, counted)
         if (grp == 0_ik .or. .not. counted) cycle
         s = s + ((va(j) - vb(j)) / (tols%atol(grp) + tols%rtol(grp)*abs(vr(j))))**2
         cnt = cnt + 1_ik
      end do
      err = sqrt(s / real(cnt, wp))
   end function state_wrms_grouped

   !---------------------------------------------------------------------------------------!
   ! Step-size factor from the normalised error err (<= 1 is within tolerance),                 !
   ! docs/science/numerical_scheme.md eq. (4): the integral controller adaptive_step_update,   !
   ! exponent -1/(p+1). The error is floored at tiny_num, so an exact step gives fmax rather    !
   ! than Inf.                                                                                  !
   !---------------------------------------------------------------------------------------!
   pure function step_control_factor(err, ec) result(fac)
      real(wp),              intent(in) :: err
      type(error_control_t), intent(in) :: ec
      real(wp) :: fac
      fac = adaptive_step_update(max(err, tiny_num), ec%safety, ec%fmin, ec%fmax, ec%p_order)
   end function step_control_factor

   !----- Each water store's scale [kg/m2], against which its tolerances are set: a soil layer's or  !
   !      the pond's is a metre of water per metre of the (top) layer, so its tolerances read as     !
   !      water content; a cohort's is the water its wood and leaf cells hold at full turgor, shared  !
   !      by both stores so a leafless cohort's empty leaf is weighed against its plant. ---------!
   pure subroutine water_scales(static, profile, hyd, scale_soil, scale_coh)
      type(fast_static_state_t), intent(in)  :: static
      type(soil_profile_t),      intent(in)  :: profile
      type(hydraulics_table_t),  intent(in)  :: hyd
      real(wp),                  intent(out) :: scale_soil(:), scale_coh(:,:)
      integer(ik) :: k, i
      real(wp)    :: s
      scale_soil(POND) = rho_h2o * profile%dz(1)
      do k = 1_ik, profile%n_active
         scale_soil(1_ik + k) = rho_h2o * profile%dz(k)
      end do
      do i = 1_ik, static%n
         associate (lc => hyd%pft(static%pft(i))%leaf_curve, wc => hyd%pft(static%pft(i))%wood_curve)
            s = static%nplant(i) * ((1.0_wp - lc%apoplast_frac) * lc%water_sat * static%bleaf(i)                &
                                  + (1.0_wp - wc%apoplast_frac) * wc%water_sat * static%bwood_water(i))
         end associate
         scale_coh(:, i) = max(s, tiny_num)
      end do
   end subroutine water_scales

end module ed_fast_integ_control
