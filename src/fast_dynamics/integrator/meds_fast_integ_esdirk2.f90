! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_integ_esdirk2 -- the default fast-loop scheme, an L-stable two-stage ESDIRK2   !
! with the ARS(2,2,2) tableau (config name "esdirk2"; its explicit part is empty, so it is !
! not an IMEX scheme): esdirk2_step is one step with its embedded error and its b-weighted !
! boundary-flux amounts. The stepper marches it (docs/science/numerical_scheme.md §3       !
! `esdirk2`, eqs. (1)-(2), §5b).                                                           !
!==========================================================================================!
module meds_fast_integ_esdirk2
   use meds_kinds, only : wp, ik
   use meds_site_state_vars, only : site_surface_t
   use meds_fast_integ_types, only : march_state_t, frozen_t, stage_bflux_t, boundary_flux_t, error_control_t,   &
                               esdirk2_stages_t, state_init, state_extrap, state_err_diff, bflux_bweight, clamp_cas
   use meds_fast_integ_backward_euler, only : backward_euler_stage, advance_plant_water, advance_surf_water_full
   implicit none
   private

   public :: esdirk2_step

contains

   !---------------------------------------------------------------------------------------!
   ! esdirk2_step -- one ESDIRK2 step, docs/science/numerical_scheme.md eqs. (1)-(2): the  !
   ! ARS(2,2,2) tableau (Ascher, Ruuth & Spiteri 1997), two backward-Euler stages over gamma*dt !
   ! (backward_euler_stage), stiffly accurate, so Y3 is the step. Soil water moves at the  !
   ! network's rate in the stages; plant water and the films are advanced once over dt on  !
   ! the stages' b-weighted fluxes, so neither adds to the embedded error y_err. ----------!
   subroutine esdirk2_step(y, frozen, site, ec, n, nsl, dt, y_out, y_err, stages, bf, clamp_n)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      type(error_control_t), intent(in)  :: ec     !< its budget stop (backward_euler_stage)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(march_state_t),  intent(inout) :: y_out, y_err   !< filled in place
      type(esdirk2_stages_t),    intent(inout) :: stages    !< this thread's stage storage, filled in place
      type(boundary_flux_t), optional, intent(out) :: bf   !< b-weighted boundary-flux AMOUNTS over dt (ledger)
      !----- Stage-clamp count (the report's WORK_CLAMP_STAGE): only the stage-3 base is clamped. ---!
      integer(ik), optional, intent(inout) :: clamp_n
      real(wp), parameter :: GAMMA = 0.2928932188134524_wp   ! 1 - 1/sqrt(2)
      real(wp), parameter :: BETA  = 2.4142135623730951_wp   ! (1-gamma)/gamma = 1 + sqrt(2)
      type(stage_bflux_t)   :: bf2, bf3
      real(wp)              :: transp_bw(n)
      real(wp)              :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      associate (Y2 => stages%y2, base3 => stages%base3, Y3 => stages%y3, sf2 => stages%sf2, sf3 => stages%sf3)
         !----- Stage 2: a gamma*dt backward-Euler stage from y_n. --------------------------------------!
         call backward_euler_stage(y, frozen, site, ec, n, nsl, GAMMA*dt, Y2, sf2, bf=bf2)
         !----- Stage 3 from the extrapolated base. Extrapolating by BETA can carry the canopy-air      !
         !      enthalpy to a temperature where qsat overflows, so the base's canopy air is clamped to a  !
         !      physical range: an identity in range, and an oversized step is rejected anyway. Soil     !
         !      water needs no clamp: it moves at a steady rate, so the base stays between its ends. ----!
         call state_extrap(y, BETA, Y2, n, nsl, base3)
         call clamp_cas(base3, nfire=clamp_n)
         call backward_euler_stage(base3, frozen, site, ec, n, nsl, GAMMA*dt, Y3, sf3, bf=bf3)
         call state_init(Y3, n, nsl, y_out)
         !----- Plant water over the full dt at the network's rates, debited the same b-weighted         !
         !      transpiration the canopy air's vapour balance was credited. -----------------------------!
         transp_bw(1:n) = (1.0_wp - GAMMA)*sf2%transp_c(1:n) + GAMMA*sf3%transp_c(1:n)
         call advance_plant_water(y, frozen, n, dt, transp_bw, y_out)
         !----- The films the same way (a no-op with canopy_water_on off). ------------------------------!
         film_evap_leaf_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_leaf(1:n) + GAMMA*sf3%film_evap_leaf(1:n)
         film_evap_wood_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_wood(1:n) + GAMMA*sf3%film_evap_wood(1:n)
         call advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
         !----- Embedded first-order error, eq. (2). ------------------------------------------------------!
         call state_err_diff(Y3, base3, Y2, y, n, nsl, y_err)
         !----- b-weighted boundary-flux amounts, b = (0, 1-gamma, gamma), which telescope exactly. -----!
         if (present(bf)) then
            call bflux_bweight(bf, bf2, bf3, dt, GAMMA)
            !----- ...the tissue-temperature time integrals the tissues' energy is committed from. -----!
            allocate(bf%tissue_leaf_int(n), bf%tissue_wood_int(n))
            bf%tissue_leaf_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%leaf_temp(1:n) + GAMMA*sf3%leaf_temp(1:n))
            bf%tissue_wood_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%wood_temp(1:n) + GAMMA*sf3%wood_temp(1:n))
            !----- ...and each cohort's transpiration, the transp_bw the leaf water was debited. -------!
            allocate(bf%transp_coh(n))
            bf%transp_coh(1:n) = dt * transp_bw(1:n)
         end if
      end associate
   end subroutine esdirk2_step

end module meds_fast_integ_esdirk2
