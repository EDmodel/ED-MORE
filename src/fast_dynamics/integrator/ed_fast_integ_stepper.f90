! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_stepper -- one dt_fast of one patch. fast_step takes the run's
! configuration, the patch's site, its static state (what is fixed for the slow step), the
! forcing, the dynamic state it advances, the output it fills and this thread's workspace,
! and runs the step in order: the frozen record, the water network's step between its two
! parts, the march, and the update that commits it (docs/science/numerical_scheme.md
! §2-§3). adaptive_march is one controller for both schemes, ESDIRK2 (the default) or RK45
! (the explicit Cash-Karp integrator and the tests' explicit reference); a step RK45
! cannot resolve is marched again with ESDIRK2 on the same frozen record. Also each
! cohort's psi_leaf, and H and LE (atm_fluxes). There is no operator-split scheme.
!---~---
module ed_fast_integ_stepper
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : latent_heat_vap, tiny_num
   use ed_config,             only : ed_config_t
   use ed_fast_params,        only : INTEG_ESDIRK2, INTEG_RK45
   use ed_site_state_vars,    only : site_surface_t
   use ed_fast_integ_types,       only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,     &
                                     fast_workspace_t, step_exports_t, march_state_t, frozen_t,           &
                                     boundary_flux_t, error_control_t, WORK_INTEG_STEPS, WORK_INTEG_REJECT,        &
                                     WORK_RK45_RESCUE, WORK_CLAMP_STAGE, state_init, state_sub, zero_like,         &
                                     bflux_zero, bflux_add
   use ed_fast_integ_control,       only : build_error_control, state_wrms_grouped, step_control_factor
   use ed_fast_integ_frozen,        only : build_frozen, record_water_flows
   use ed_fast_integ_backward_euler,    only : water_step
   use ed_fast_integ_esdirk2,         only : esdirk2_step
   use ed_fast_integ_rk45,          only : rk45_step, rk45_state_railed, RK45_P_ORDER, RK45_WORK_CAP
   use ed_fast_integ_update,        only : commit_amounts_t, update_dynamic_state, close_ledgers
   use ed_plant_hydraulics, only : leaf_node_psi
   implicit none
   private

   public :: fast_step
   public :: adaptive_march   ! used only by the tests

   !--- The ESDIRK2 march's trial cap, a backstop that should never trigger.
   integer(ik), parameter :: ESDIRK2_TRIAL_CAP = 4096_ik

contains

   !---~---
   !   Advance one dt_fast with the configured scheme ([fast].time_integrator):
   !   "esdirk2" (default) -- the 2-solve ESDIRK2, gamma = 1 - 1/sqrt(2); not IMEX, since the
   !                          biotic CO2 source is implicit and the explicit tableau is empty
   !                          (docs/science/numerical_scheme.md eq. (1)).
   !   "rk45"              -- the adaptive explicit Cash-Karp march, a full alternative a run
   !                          can choose and the explicit reference the tests compare against.
   !---~---
   subroutine fast_step(dt_fast, cfg, site, static, forcing, dynamic, out, work)
      real(wp),                   intent(in)    :: dt_fast
      type(ed_config_t),          intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site      ! the patch's site: its soil and its ground
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      type(fast_output_t),        intent(inout) :: out       ! sized by the caller (ensure_output_capacity)
      type(fast_workspace_t),     intent(inout) :: work      ! this thread's step storage, reused step to step
      type(march_state_t)   :: y, y_out
      type(boundary_flux_t)   :: acc
      type(commit_amounts_t) :: moved
      type(error_control_t)  :: ec, ec_rk45
      real(wp)    :: dt0, dt_warm
      integer(ik) :: jcoh, n, nsl, nsteps, nrej
      logical     :: kept

      n   = static%n
      nsl = site%soil_profile%n_active
      out%report%work(WORK_RK45_RESCUE) = 0_ik
      if (out%capture) out%cdiag(:, 1:n) = 0.0_wp

      !---~---
      !   Each cohort's psi_leaf for the daily-maximum accumulator behind the next day's stomatal
      ! control (docs/science/leaf_gas_exchange.md §4.2), at state^n, the value the leaf kernel
      ! is handed.
      !---~---
      do jcoh = 1_ik, n
         out%psi_leaf_coh(jcoh) = leaf_node_psi(dynamic%leaf_water_mass(jcoh), static%bleaf(jcoh),              &
               dynamic%wood_water_mass(jcoh), static%bwood_water(jcoh),                                            &
               static%height(jcoh), cfg%pft%hydraulics%pft(static%pft(jcoh)))
      end do

      !--- The run's error control: tolerances and controller, eqs. (3)-(4), and the budget stops.
      ec = build_error_control(cfg)

      !---~---
      !   The frozen record, with the soil-plant water network stepped between its two parts: its
      ! time-mean flows are constants of the march (operator splitting; docs/science/
      ! numerical_scheme.md §2).
      !---~---
      associate (frozen => work%frozen)
      call build_frozen(dt_fast, cfg, site, static, forcing, dynamic, out, work%water, frozen, y)
      call water_step(dt_fast, static, site%soil_profile, site%soil, cfg%pft%hydraulics, frozen%water_drive,     &
                      ec%tols, work%water, out%report%work)
      call record_water_flows(dt_fast, cfg, site, static, dynamic, work%water, out, frozen)

      !---~---
      !   Warm start (docs/science/numerical_scheme.md §3 `esdirk2`): from the step this patch settled
      ! on last dt_fast; a first call starts from dt_fast. A fixed-sub-step march keeps it.
      !---~---
      dt0 = dt_fast
      if (dynamic%adapt_dt_last > tiny_num) dt0 = min(dynamic%adapt_dt_last, dt_fast)
      dt_warm = dynamic%adapt_dt_last

      kept = .false.
      if (cfg%time_integrator == INTEG_RK45) then
         !---~---
         !   RK45 is explicit over the whole column. On a dense, cold canopy the leaf <-> canopy-air
         ! coupling can be too stiff even at the sub-step floor: the march hits its work cap, or
         ! its 5th- and 4th-order solutions rail together and the controller accepts the railed
         ! state. Then the step is marched again with ESDIRK2, whose Newton solves that coupled
         ! pair, from the same frozen record (docs/science/numerical_scheme.md §3, `rk45`). The
         ! dynamic state is written only by the update, so nothing needs rolling back.
         !---~---
         ec_rk45 = ec
         ec_rk45%p_order = RK45_P_ORDER
         out%report%work(WORK_CLAMP_STAGE) = 0_ik
         out%report%theta_ood_max = 0.0_wp
         call adaptive_march(INTEG_RK45, y, frozen, site, n, nsl, dt_fast, ec_rk45, dt0, y_out, work, nsteps, nrej, &
                             acc=acc, dt_warm_out=dt_warm, clamp_n=out%report%work(WORK_CLAMP_STAGE),               &
                             ood_max=out%report%theta_ood_max)
         kept = nsteps + nrej < RK45_WORK_CAP .and. .not. rk45_state_railed(y_out, site, nsl)
         if (.not. kept) then
            out%report%work(WORK_RK45_RESCUE) = out%report%work(WORK_RK45_RESCUE) + 1_ik
            dt_warm = dynamic%adapt_dt_last
         end if
      end if

      if (.not. kept) then
         !---~---
         !   ESDIRK2: the default, and the RK45 rescue. It has no explicit stage, so theta_ood_max
         ! stays 0.
         !---~---
         out%report%work(WORK_CLAMP_STAGE) = 0_ik
         out%report%theta_ood_max = 0.0_wp
         call adaptive_march(INTEG_ESDIRK2, y, frozen, site, n, nsl, dt_fast, ec, dt0, y_out, work, nsteps, nrej,       &
                             acc=acc, dt_warm_out=dt_warm, clamp_n=out%report%work(WORK_CLAMP_STAGE))
      end if
      !--- Work counters; the water network's own work is counted in water_step.
      out%report%work(WORK_INTEG_STEPS)  = nsteps
      out%report%work(WORK_INTEG_REJECT) = nrej

      call update_dynamic_state(dt_fast, cfg, site, static, frozen, acc, dt_warm, y_out, dynamic, out, moved)
      call close_ledgers(dt_fast, site, forcing, static, frozen, y, y_out, acc, moved, ec%debug_error, out%report)
      call atm_fluxes(out%report%exports, dt_fast, out%le_flux, out%h_flux)
      out%rnet = acc%rad_net / dt_fast
      end associate
   end subroutine fast_step

   !---~---
   !   adaptive_march -- march the column to t_end with one scheme (INTEG_ESDIRK2, ESDIRK2; or
   ! INTEG_RK45, Cash-Karp), accepting a step when the grouped norm of its embedded error is
   ! within tolerance (docs/science/numerical_scheme.md eq. (3)) and sizing the next by
   ! step_control_factor (eq. (4), p from ec). ESDIRK2's error is measured against its
   ! first-order embedded solution, relative to the step's start state; RK45's is the 5th - 4th
   ! difference, relative to the 5th-order solution. The march stops after trial_cap trials:
   ! RK45's work cap (the stepper then marches ESDIRK2) or ESDIRK2's backstop. Returns the step
   ! and reject counts; acc sums the accepted steps' boundary-flux amounts.
   !---~---
   subroutine adaptive_march(scheme, y0, frozen, site, n, nsl, t_end, ec, dt_init, y_out, work, nsteps, nrej,       &
                             acc, dt_warm_out, clamp_n, ood_max)
      integer(ik),           intent(in)  :: scheme
      type(march_state_t),  intent(in)  :: y0
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec       ! tolerances and controller (ed_fast_integ_control)
      type(march_state_t),  intent(inout) :: y_out  ! filled in place
      type(fast_workspace_t), intent(inout) :: work  ! this thread's march and stage storage
      integer(ik),           intent(out) :: nsteps, nrej
      type(boundary_flux_t), optional, intent(out) :: acc  ! accumulated boundary-flux amounts (ledger)
      !---~---
      !   Warm start (docs/science/numerical_scheme.md §3 `esdirk2`): the last accepted step that the
      ! interval's end did not cut short, or larger, what an end-cut step's controller proposes,
      ! to seed the next call.
      !---~---
      real(wp),    optional, intent(out) :: dt_warm_out
      !---~---
      !   Stage-clamp count over the march, rejected trials included: a clamp still says the
      ! controller probed too far. Accumulated, so the caller zeroes it.
      !---~---
      integer(ik), optional, intent(inout) :: clamp_n
      !---~---
      !   [m3/m3] running max of theta's excursion outside its domain at each trial's start, each a
      ! place an explicit right-hand side was evaluated (theta_ood_max). Accumulated.
      !---~---
      real(wp),    optional, intent(inout) :: ood_max

      type(boundary_flux_t) :: bfsub
      real(wp)             :: t, dt, err, fac, dt_floor
      real(wp)             :: dt_try, dt_warm
      integer(ik)          :: trial_cap, k
      logical              :: clamped
      trial_cap = ESDIRK2_TRIAL_CAP
      if (scheme == INTEG_RK45) trial_cap = RK45_WORK_CAP
      associate (y => work%y_march, y_new => work%y_new, y_err => work%y_err, y_lo => work%y_lo)
         if (present(acc)) call bflux_zero(acc, n)
         !---~---
         !   Sub-step floor t_end/64, at least 1e-2 s. ESDIRK2's backward-Euler stages are L-stable,
         ! so a floor step stays bounded even when its error is above tolerance; the floor caps a
         ! stiff transient at 64 sub-steps and surfaces a non-finite state promptly. For RK45 a
         ! floor step is not guaranteed stable, so reaching the floor marks a pathological step.
         !---~---
         dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)   ! clamp-ok: the smallest sub-step [s] (above)

         call state_init(y0, n, nsl, y)
         t = 0.0_wp
         dt = min(dt_init, t_end)
         nsteps = 0_ik
         nrej = 0_ik
         dt_warm = dt                                                ! fallback if every step is end-clamped
         do
            if (t >= t_end - tiny_num) exit
            dt_try = dt
            dt = min(dt, t_end - t)
            clamped = dt < dt_try - tiny_num
            if (present(ood_max)) then
               do k = 1_ik, nsl
                  ood_max = max(ood_max, y%theta(k) - site%soil_profile%theta_sat(k),                             &
                                site%soil_profile%theta_res(k) - y%theta(k))
               end do
            end if
            !---~---
            !   The step and its error under the grouped norm, eq. (3), over the states the march
            ! integrates (docs/science/numerical_scheme.md §3 `esdirk2`).
            !---~---
            select case (scheme)
            case (INTEG_RK45)
               call rk45_step(y, frozen, site, n, nsl, dt, y_new, y_err, bf=bfsub, clamp_stage_n=clamp_n)
               y_lo = zero_like(y_err, n, nsl)
               err = state_wrms_grouped(y_err, y_lo, y_new, n, nsl, ec%tols)
            case default
               call esdirk2_step(y, frozen, site, ec, n, nsl, dt, y_new, y_err, work%stages, bf=bfsub,      &
                                     clamp_n=clamp_n)
               call state_sub(y_new, y_err, n, nsl, y_lo)               ! the 1st-order embedded solution
               err = state_wrms_grouped(y_new, y_lo, y, n, nsl, ec%tols)
            end select
            !---~---
            !   A non-finite error means the step was too big (typically the canopy air overshot to
            ! where qsat overflows): reject it and shrink by fmin, since the NaN poisons the factor.
            ! A non-finite floor step is committed, so the driver's state guard reports it.
            !---~---
            if (err /= err .or. dt /= dt) then
               if (dt <= dt_floor) then
                  call state_init(y_new, n, nsl, y)
                  t = t + dt_floor
                  nsteps = nsteps + 1_ik
                  exit
               end if
               nrej = nrej + 1_ik
               dt = max(dt * ec%fmin, dt_floor)
               cycle
            end if
            fac = step_control_factor(err, ec)
            if (err <= 1.0_wp .or. dt <= dt_floor) then
               call state_init(y_new, n, nsl, y)
               if (present(acc)) call bflux_add(acc, bfsub)          ! accumulate only accepted sub-steps
               t = t + dt
               nsteps = nsteps + 1_ik
               !---~---
               !   The warm-start seed is the step just accepted, before the growth factor: a grown
               ! proposal over-estimates (docs/science/numerical_scheme.md §3 `esdirk2`). A step the
               ! interval's end cut short can only raise it, to the step the controller asked for
               ! before the cut, capped at the interval: otherwise a seed below the interval stays
               ! there, every later interval taking a step at the seed and an end-cut one.
               !---~---
               if (clamped) then
                  dt_warm = max(dt_warm, min(dt_try, t_end))
               else
                  dt_warm = dt
               end if
               dt = dt * fac
            else
               nrej = nrej + 1_ik
               dt = dt * fac
            end if
            if (nsteps + nrej >= trial_cap) exit
         end do
         call state_init(y, n, nsl, y_out)
         if (present(dt_warm_out)) dt_warm_out = dt_warm
      end associate
   end subroutine adaptive_march

   !---~---
   !   atm_fluxes -- the canopy-air -> atmosphere fluxes the run reports, from the ledger's own
   ! boundary export over the step, so the reported ET is the conserved vapour export:
   ! LE = L_v * vapour export / dt and H = sensible export / dt, the b-weighted
   ! g_atm_heat*cp_air*(T_cas - theta_atm), not the enthalpy export minus LE
   ! (docs/science/numerical_scheme.md §5b).
   !---~---
   pure subroutine atm_fluxes(exports, dt_fast, le_flux, h_flux)
      type(step_exports_t), intent(in)  :: exports
      real(wp),             intent(in)  :: dt_fast
      real(wp),             intent(out) :: le_flux, h_flux   ! [W/m2] step-mean CAS -> atmosphere
      le_flux = exports%atm_vap_export  / dt_fast * latent_heat_vap
      h_flux  = exports%atm_heat_export / dt_fast
   end subroutine atm_fluxes

end module ed_fast_integ_stepper
