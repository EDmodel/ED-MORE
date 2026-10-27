! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_plant_phenology -- leaf phenology (docs/science/plant_phenology.md). Daily cues and
! a PFT's traits become two smoothed tendencies in [0,1], leaf_flush_tendency and
! leaf_shed_tendency; the carbon allocation turns them into leaf growth and leaf loss.
!   * phenology_kernel   -- one daily update: accumulate the cue memory, switch each cue,
!                           combine the switches, smooth (§1-§4, eqs. (1)-(8)).
!   * leaf_turnover_step -- one step of leaf loss and the flush cap (§5, eq. (9)), also
!                           exposed through the C API; turnover_gate gives the share of
!                           its background turnover that acts.
! It reads the PFT's own record, phenology_traits_t (ed_pft_params, with the cue bits
! CUE_*), and its own: pheno_env_t (the day's cues), pheno_state_t (the cohort's memory,
! read and written back) and pheno_out_t (the potential rates).
!---~---
module ed_plant_phenology
   use ed_kinds,            only : wp, ik
   use ed_constants,        only : yr_day
   use ed_numerics,         only : logistic, clamp01
   use ed_time,             only : doy_effective
   use ed_state_defaults,   only : is_unset_positive
   use ed_pft_params,       only : phenology_traits_t, CUE_TEMP, CUE_LIGHT, CUE_WATER
   implicit none
   private

   public :: pheno_env_t, pheno_state_t, pheno_out_t
   public :: phenology_kernel, leaf_turnover_step, turnover_gate

   !--- Daily environmental drivers (read-only).
   type :: pheno_env_t
      real(wp)    :: temp_day         = 0.0_wp    ! [K]   daily-mean air temperature
      real(wp)    :: par_hours        = 12.0_wp   ! [h/day] hours with PAR at the cohort's top above par_min
      real(wp)    :: predawn_leaf_psi = 0.0_wp    ! [MPa, <=0] predawn (daily-max) leaf water potential
      integer(ik) :: doy              = 1_ik      ! [-]   day of year
      logical     :: hemis_north      = .true.    ! northern hemisphere (season windows)
   end type pheno_env_t

   !--- The prognostic memory: two smoothed tendencies + the cue accumulators.
   type :: pheno_state_t
      real(wp) :: leaf_flush_tendency = 1.0_wp   ! [-] smoothed flush signal; a cohort is born flushing
      real(wp) :: leaf_shed_tendency  = 0.0_wp   ! [-] smoothed senescence signal
      real(wp) :: growing_degree_days = 0.0_wp   ! [K day] warmth sum since midwinter     (TEMP flush)
      real(wp) :: cold_degree_days    = 0.0_wp   ! [K day] cold sum since midsummer       (TEMP shed)
      real(wp) :: wet_psi_sum         = 0.0_wp   ! [MPa day] sum of psi above the TLP     (WATER flush)
      real(wp) :: dry_psi_sum         = 0.0_wp   ! [MPa day] sum of psi below the TLP     (WATER shed)
      real(wp) :: light_hours_mean    = -1.0_wp  ! [h/day] running-mean par_hours (LIGHT); < 0: unset, the first day sets it
   end type pheno_state_t

   !--- Outputs: the potential relative rates the tendencies allow.
   type :: pheno_out_t
      real(wp) :: leaf_flush_potential = 0.0_wp   ! [1/day] flush_rate_max * leaf_flush_tendency
      real(wp) :: leaf_shed_potential  = 0.0_wp   ! [1/day] shed_rate_max  * leaf_shed_tendency
   end type pheno_out_t

   !--- Northern-equivalent day of the summer solstice, where the cold sum starts (§2).
   integer(ik), parameter :: MIDSUMMER_DOY     = 172_ik
   !---~---
   !   At or below this flush rate a canopy is dormant and may snap to bare (§5): a numerical
   ! threshold, not a trait.
   !---~---
   real(wp),    parameter :: DORMANT_FLUSH_EPS = 1.0e-6_wp    ! [1/day]
   !--- At or below this shed tendency no leaves drop, so background turnover runs (§5).
   real(wp),    parameter :: NO_DROP_SHED_TENDENCY = 0.01_wp   ! [-]

contains

   !--- Advance one cohort's phenology over one step of dt days.
   pure subroutine phenology_kernel(env, params, dt, state, out)
      type(pheno_env_t),    intent(in)    :: env
      type(phenology_traits_t), intent(in) :: params
      real(wp),             intent(in)    :: dt
      type(pheno_state_t),  intent(inout) :: state
      type(pheno_out_t),    intent(out)   :: out
      integer(ik) :: both
      real(wp)    :: wet_before, dry_before, wet_switch, dry_switch, s_flush, seasonal, s_shed, w

      both = ior(params%flush_cue_mask, params%shed_cue_mask)
      wet_before = wet_water_switch(params, state)
      dry_before = dry_water_switch(params, state)
      call accumulate(env, params, dt, both, state)

      !---~---
      !   Water has no calendar: each sum restarts when the other side's switch crosses 0.5
      ! upward today, the crossing and not the level (§2).
      !---~---
      if (iand(both, CUE_WATER) /= 0_ik) then
         if (wet_before <= 0.5_wp .and. wet_water_switch(params, state) > 0.5_wp) state%dry_psi_sum = 0.0_wp
         if (dry_before <= 0.5_wp .and. dry_water_switch(params, state) > 0.5_wp) state%wet_psi_sum = 0.0_wp
      end if
      wet_switch = wet_water_switch(params, state)
      dry_switch = dry_water_switch(params, state)

      !--- (3a) Flush: every enabled cue must permit it (an empty mask always flushes).
      s_flush = 1.0_wp
      if (iand(params%flush_cue_mask, CUE_TEMP) /= 0_ik) s_flush = s_flush *                     &
         logistic(params%flush_temp_sharpness * (state%growing_degree_days - params%flush_degree_days))
      if (iand(params%flush_cue_mask, CUE_LIGHT) /= 0_ik) s_flush = s_flush *                    &
         logistic(params%flush_light_sharpness * (state%light_hours_mean - params%flush_light_hours))
      if (iand(params%flush_cue_mask, CUE_WATER) /= 0_ik) s_flush = s_flush * wet_switch

      !---~---
      !   (3b) Shed: the larger of the seasonal trigger (cold and light, each if enabled) and
      ! the water trigger; an empty mask never senesces (eq. (7)).
      !---~---
      seasonal = 0.0_wp
      if (iand(params%shed_cue_mask, CUE_TEMP + CUE_LIGHT) /= 0_ik) then
         seasonal = 1.0_wp
         if (iand(params%shed_cue_mask, CUE_TEMP) /= 0_ik) seasonal = seasonal *                 &
            logistic(params%shed_temp_sharpness * (state%cold_degree_days - params%shed_degree_days))
         if (iand(params%shed_cue_mask, CUE_LIGHT) /= 0_ik) seasonal = seasonal *                &
            logistic(params%shed_light_sharpness * (state%light_hours_mean - params%shed_light_hours))
      end if
      s_shed = seasonal
      if (iand(params%shed_cue_mask, CUE_WATER) /= 0_ik) s_shed = max(s_shed, dry_switch)

      !--- (4) Smooth into the tendencies (eq. (8)) and map to the potential rates (eq. (1)).
      w = dt / max(params%flush_cue_timescale, dt)
      state%leaf_flush_tendency = clamp01(state%leaf_flush_tendency + w * (s_flush - state%leaf_flush_tendency))
      w = dt / max(params%shed_cue_timescale, dt)
      state%leaf_shed_tendency  = clamp01(state%leaf_shed_tendency  + w * (s_shed  - state%leaf_shed_tendency))
      out%leaf_flush_potential = params%flush_rate_max * state%leaf_flush_tendency
      out%leaf_shed_potential  = params%shed_rate_max  * state%leaf_shed_tendency
   end subroutine phenology_kernel

   !--- The two water switches, from the sums of predawn leaf psi above and below the TLP.
   pure real(wp) function wet_water_switch(params, state) result(sw)
      type(phenology_traits_t), intent(in) :: params
      type(pheno_state_t),  intent(in) :: state
      sw = logistic(params%flush_water_sharpness * (state%wet_psi_sum - params%flush_water_sum))
   end function wet_water_switch

   pure real(wp) function dry_water_switch(params, state) result(sw)
      type(phenology_traits_t), intent(in) :: params
      type(pheno_state_t),  intent(in) :: state
      sw = logistic(params%shed_water_sharpness * (state%dry_psi_sum - params%shed_water_sum))
   end function dry_water_switch

   !--- Accumulate the cue memory from today's drivers, for the cues either side uses (§2).
   pure subroutine accumulate(env, params, dt, both, state)
      type(pheno_env_t),    intent(in)    :: env
      type(phenology_traits_t), intent(in) :: params
      real(wp),             intent(in)    :: dt
      integer(ik),          intent(in)    :: both
      type(pheno_state_t),  intent(inout) :: state
      integer(ik) :: de
      real(wp)    :: w

      !---~---
      !   Temperature: both sums restart at midwinter (the first step of the northern-
      ! equivalent year); warmth counts all year, cold only once days shorten (eq. (3)).
      !---~---
      if (iand(both, CUE_TEMP) /= 0_ik) then
         de = doy_effective(env%doy, env%hemis_north)
         if (real(de, wp) < 1.0_wp + dt) then
            state%growing_degree_days = 0.0_wp
            state%cold_degree_days    = 0.0_wp
         end if
         state%growing_degree_days = state%growing_degree_days + max(0.0_wp, env%temp_day - params%flush_base_temp) * dt
         if (de >= MIDSUMMER_DOY)                                                                &
            state%cold_degree_days = state%cold_degree_days + max(0.0_wp, params%shed_base_temp - env%temp_day) * dt
      end if

      !--- Water: predawn leaf psi above and below the turgor-loss point (eq. (5)).
      if (iand(both, CUE_WATER) /= 0_ik) then
         state%wet_psi_sum = state%wet_psi_sum + max(0.0_wp, env%predawn_leaf_psi - params%leaf_psi_tlp) * dt
         state%dry_psi_sum = state%dry_psi_sum + max(0.0_wp, params%leaf_psi_tlp - env%predawn_leaf_psi) * dt
      end if

      !---~---
      !   Light: a running mean of the daily hours of light (eq. (4)); a cohort with no light
      ! memory yet (a cold start or a recruit) starts from today's hours.
      !---~---
      if (iand(both, CUE_LIGHT) /= 0_ik) then
         if (is_unset_positive(state%light_hours_mean)) then
            state%light_hours_mean = env%par_hours
         else
            w = dt / max(params%light_window, dt)
            state%light_hours_mean = state%light_hours_mean + w * (env%par_hours - state%light_hours_mean)
         end if
      end if
   end subroutine accumulate

   !---~---
   !   One step of leaf loss and the flush cap, in whatever carbon unit `leaf` is in: the
   ! background turnover (gated by turnover_gate), the senescence that stops at
   ! min_leaf_cover, and the flush cap. A dormant canopy that senescence would leave
   ! below bare_leaf_cover snaps to bare, when its floor lies below that cover (§5, eq. (9)).
   !---~---
   pure subroutine leaf_turnover_step(leaf, leaf_full, leaf_flush_tendency, leaf_shed_tendency,   &
                                      flush_rate_max, shed_rate_max, leaf_turnover_rate,          &
                                      min_leaf_cover, bare_leaf_cover, dt,                        &
                                      senescence, background, flush_cap)
      real(wp), intent(in)  :: leaf, leaf_full                          ! [carbon] now and full canopy
      real(wp), intent(in)  :: leaf_flush_tendency, leaf_shed_tendency  ! [-]
      real(wp), intent(in)  :: flush_rate_max, shed_rate_max            ! [1/day]
      real(wp), intent(in)  :: leaf_turnover_rate                       ! [1/yr] background loss
      real(wp), intent(in)  :: min_leaf_cover, bare_leaf_cover          ! [-] fractions of leaf_full
      real(wp), intent(in)  :: dt                                       ! [day]
      real(wp), intent(out) :: senescence, background, flush_cap        ! [carbon] this step
      real(wp) :: pool

      pool       = max(leaf, 0.0_wp)
      background = min(max(leaf_turnover_rate, 0.0_wp) / yr_day                                    &
                       * turnover_gate(leaf_flush_tendency, leaf_shed_tendency) * pool * dt, pool)
      senescence = min(max(shed_rate_max, 0.0_wp) * leaf_shed_tendency * pool * dt,                  &
                       max(0.0_wp, pool - background - min_leaf_cover * leaf_full))
      if (flush_rate_max * leaf_flush_tendency <= DORMANT_FLUSH_EPS .and. senescence > 0.0_wp      &
          .and. min_leaf_cover < bare_leaf_cover                                                  &
          .and. pool - background - senescence < bare_leaf_cover * leaf_full)                     &
         senescence = pool - background
      flush_cap  = max(flush_rate_max, 0.0_wp) * leaf_flush_tendency * leaf_full * dt
   end subroutine leaf_turnover_step

   !---~---
   !   The share of the background turnover rate that acts today: the flush tendency while no
   ! leaves drop, nothing while they do; the fine roots turn over on the same gate (§5).
   !---~---
   pure real(wp) function turnover_gate(leaf_flush_tendency, leaf_shed_tendency) result(gate)
      real(wp), intent(in) :: leaf_flush_tendency, leaf_shed_tendency   ! [-]
      gate = 0.0_wp
      if (leaf_shed_tendency <= NO_DROP_SHED_TENDENCY) gate = leaf_flush_tendency
   end function turnover_gate

end module ed_plant_phenology
