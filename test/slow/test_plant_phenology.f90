! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_plant_phenology -- unit tests for the phenology kernel, leaf_turnover_step and gate.
!
!   The kernel turns daily cues into two smoothed tendencies (leaf_flush_tendency,
! leaf_shed_tendency) and the potential rates rate_max * tendency
! (docs/science/plant_phenology.md eqs. (1)-(9)). The tests:
!
!   1. NO CUES             : both masks empty => flush potential = flush_rate_max, shed 0.
!   2. TEMPERATE DECIDUOUS : warmth x light flush, cold x short-day shed (the hours of light
!                            at a low par_min are the day length); the sums reset
!                            at midwinter and the cold sum waits for midsummer.
!   3. SOUTHERN HEMISPHERE : the same PFT half a year later.
!   4. DROUGHT DECIDUOUS   : wet/dry psi sums; a long wet season does not block a drought, a
!                            brief rain does not wipe one.
!   5. LIGHT EXCHANGING    : many bright hours trigger senescence, flush stays on; a cohort
!                            with no light memory starts from its first day's hours.
!   6. RATE MAPPING        : potential == rate_max * tendency, tendencies in [0,1].
!   7. DEGENERATE / FPE    : extreme drivers keep the potentials finite and >= 0.
!   8. LEAF TURNOVER STEP  : background, senescence, the leaf-cover floor (an emergent
!                            evergreen), the dormant snap to bare, the pool clamp.
!   9. THE C-API SHIM      : the bind(c) mirrors carry every value to the kernel and back.
!---~---
program test_plant_phenology
   use ed_test_assert, only : check, check_close, check_true, banner, test_report
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use ed_c_api_phenology,   only : pheno_env_c, pheno_args_c, pheno_state_c, pheno_out_c,        &
                                    ed_phenology_step, ed_leaf_turnover_step, ed_daylength
   use ed_kinds,             only : wp, ik
   use ed_constants,         only : yr_day
   use ed_time,              only : daylength
   use ed_plant_phenology,   only : pheno_env_t, pheno_state_t, pheno_out_t,               &
                                    phenology_kernel, leaf_turnover_step, turnover_gate
   use ed_pft_params, only : phenology_traits_t, CUE_NONE, CUE_TEMP, CUE_LIGHT, CUE_WATER, CUE_ALL
   implicit none

   real(wp), parameter :: TWOPI = 6.283185307179586_wp
   real(wp), parameter :: LAT   = 42.5_wp                 ! Harvard Forest

   call test_no_cues()
   call test_temperate_deciduous()
   call test_southern_hemisphere()
   call test_drought_deciduous()
   call test_light_exchanging()
   call test_rate_mapping()
   call test_degenerate()
   call test_leaf_turnover_step()
   call test_c_api_shim()

   call test_report('test_plant_phenology')

contains
   !--- Annual air temperature (northern hemisphere; warmest ~ doy 201).
   pure real(wp) function annual_temp(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 12.0_wp * cos(TWOPI * real(doy - 201_ik, wp) / 365.0_wp)
   end function annual_temp
   !---~---
   !   A temperate deciduous PFT (warmth x light, cold x short days). The tests feed the
   ! day length as the hours of light, which is what a low par_min counts.
   !---~---
   pure function temperate_deciduous() result(p)
      type(phenology_traits_t)   :: p
      p%flush_cue_mask        = CUE_TEMP + CUE_LIGHT
      p%shed_cue_mask         = CUE_TEMP + CUE_LIGHT
      p%flush_degree_days     = 92.0_wp
      p%shed_base_temp        = 290.37_wp
      p%shed_degree_days      = 48.0_wp
      p%flush_light_hours     = 10.35_wp
      p%shed_light_hours      = 9.83_wp
   end function temperate_deciduous
   !--- 1. No cues: always flushing, never senescing, under ANY drivers.
   subroutine test_no_cues()
      type(pheno_env_t)    :: env
      type(phenology_traits_t)   :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: flush_ok, shed_ok
      print '(a)', '-- 1. no cues (flush={}, shed={}) --'
      params%flush_cue_mask = CUE_NONE
      params%shed_cue_mask = CUE_NONE
      flush_ok = .true.
      shed_ok = .true.
      do d = 1_ik, 400_ik
         env%doy              = modulo(d - 1_ik, 365_ik) + 1_ik
         env%temp_day         = annual_temp(env%doy)
         env%predawn_leaf_psi = -8.0_wp                  ! cavitated ...
         env%par_hours        = 24.0_wp                  ! ... and never dark -- none of it is a cue
         call phenology_kernel(env, params, 1.0_wp, state, out)
         if (abs(out%leaf_flush_potential - params%flush_rate_max) > 1.0e-12_wp) flush_ok = .false.
         if (out%leaf_shed_potential /= 0.0_wp) shed_ok = .false.
      end do
      call check_true('no cues: flush potential = flush_rate_max regardless of drivers', flush_ok)
      call check_true('no cues: shed potential = 0 regardless of drivers',               shed_ok)
   end subroutine test_no_cues
   !--- Run a temperate deciduous PFT for two years and sample the second.
   subroutine run_two_years(north, offset, fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171)
      logical,     intent(in)  :: north
      integer(ik), intent(in)  :: offset          !--- calendar lag of the climate [day] (183: doy_effective)
      real(wp),    intent(out) :: fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171
      type(pheno_env_t)    :: env
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      type(phenology_traits_t)   :: params
      integer(ik) :: d, de
      real(wp)    :: lat_signed
      params = temperate_deciduous()
      lat_signed = merge(LAT, -LAT, north)
      env%hemis_north = north
      do d = 1_ik, 730_ik
         env%doy       = modulo(d - 1_ik + offset, 365_ik) + 1_ik
         de            = modulo(d - 1_ik, 365_ik) + 1_ik          ! the season, northern-equivalent
         env%temp_day  = annual_temp(de)
         env%par_hours = daylength(lat_signed, env%doy)
         call phenology_kernel(env, params, 1.0_wp, state, out)
         if (d == 366_ik) gdd_1   = state%growing_degree_days
         if (d == 536_ik) cdd_171 = state%cold_degree_days
         if (d == 565_ik) then
            fl_200 = state%leaf_flush_tendency
            sh_200 = state%leaf_shed_tendency
            gdd_200 = state%growing_degree_days
         end if
         if (d == 705_ik) then
            fl_340 = state%leaf_flush_tendency
            sh_340 = state%leaf_shed_tendency
         end if
      end do
   end subroutine run_two_years
   !--- 2. Temperate deciduous: flush in summer, senescence in late autumn.
   subroutine test_temperate_deciduous()
      real(wp) :: fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171
      print '(a)', '-- 2. temperate deciduous (flush={TEMP,LIGHT}, shed={TEMP,LIGHT}) --'
      call run_two_years(.true., 0_ik, fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171)
      call check_true('deciduous: flush tendency HIGH in mid-summer', fl_200 > 0.9_wp, fl_200)
      call check_true('deciduous: shed tendency  LOW  in mid-summer', sh_200 < 0.05_wp, sh_200)
      call check_true('deciduous: shed tendency  HIGH in late autumn', sh_340 > 0.5_wp, sh_340)
      !--- The warmth sum is still high in autumn: it is the light gate that closes.
      call check_true('deciduous: flush tendency LOW in late autumn (light gate)', fl_340 < 0.5_wp, fl_340)
      call check_true('deciduous: warmth sum restarts at midwinter',   gdd_1 < 1.0_wp, gdd_1)
      call check_true('deciduous: warmth sum accumulated by summer',   gdd_200 > 92.0_wp, gdd_200)
      call check_close('deciduous: no cold sum before midsummer',      cdd_171, 0.0_wp, 0.0_wp)
   end subroutine test_temperate_deciduous
   !--- 3. Southern hemisphere: the same seasons half a year later on the calendar.
   subroutine test_southern_hemisphere()
      real(wp) :: fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171
      print '(a)', '-- 3. southern hemisphere --'
      call run_two_years(.false., 183_ik, fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171)
      call check_true('south: flush tendency HIGH in its mid-summer', fl_200 > 0.9_wp, fl_200)
      call check_true('south: shed tendency  LOW  in its mid-summer', sh_200 < 0.05_wp, sh_200)
      call check_true('south: shed tendency  HIGH in its late autumn', sh_340 > 0.5_wp, sh_340)
      call check_true('south: warmth sum restarts at its midwinter',  gdd_1 < 1.0_wp, gdd_1)
   end subroutine test_southern_hemisphere
   !--- 4. Drought deciduous: the wet and dry sums of predawn psi against the TLP.
   subroutine test_drought_deciduous()
      type(pheno_env_t)    :: env
      type(phenology_traits_t)   :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      real(wp)    :: dry_before_rain
      print '(a)', '-- 4. drought deciduous (flush={WATER}, shed={WATER}) --'
      params%flush_cue_mask = CUE_WATER
      params%shed_cue_mask = CUE_WATER                                       ! tlp -2, sums 10
      !--- (a) A long wet season: the wet sum grows far past its centre.
      do d = 1_ik, 60_ik
         env%predawn_leaf_psi = -0.5_wp                     ! 1.5 MPa above the TLP each day
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_true('drought-decid: watered => no senescence', out%leaf_shed_potential < 0.05_wp * params%shed_rate_max)
      call check_true('drought-decid: watered => flushing',      out%leaf_flush_potential > 0.9_wp * params%flush_rate_max)
      !---~---
      !   (b) Drought: the dry sum crosses its centre DESPITE the large wet sum (the
      ! crossing resets the wet sum; a level test would wipe the dry sum every day).
      !---~---
      do d = 1_ik, 30_ik
         env%predawn_leaf_psi = -3.0_wp                     ! 1 MPa below the TLP each day
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_true('drought-decid: drought after a wet season => senescence', &
                      out%leaf_shed_potential > 0.9_wp * params%shed_rate_max, out%leaf_shed_potential)
      call check_true('drought-decid: drought => flush falls',   out%leaf_flush_potential < 0.1_wp * params%flush_rate_max)
      !--- (c) A brief rain adds wet credit but does not wipe the drought.
      dry_before_rain = state%dry_psi_sum
      do d = 1_ik, 2_ik
         env%predawn_leaf_psi = -0.5_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('drought-decid: brief rain keeps the dry sum', state%dry_psi_sum, dry_before_rain, 0.0_wp)
      call check_true('drought-decid: brief rain keeps senescing', out%leaf_shed_potential > 0.5_wp * params%shed_rate_max)
      !--- (d) Rewetting: the wet sum crosses, the dry sum resets, flushing resumes.
      do d = 1_ik, 30_ik
         env%predawn_leaf_psi = -0.5_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('drought-decid: rewet resets the dry sum', state%dry_psi_sum, 0.0_wp, 0.0_wp)
      call check_true('drought-decid: rewet => flush recovers',  out%leaf_flush_potential > 0.9_wp * params%flush_rate_max)
      call check_true('drought-decid: rewet => senescence stops', out%leaf_shed_potential < 0.05_wp * params%shed_rate_max)
   end subroutine test_drought_deciduous
   !--- 5. Light-driven leaf exchange: flush stays on, senescence rises with light.
   subroutine test_light_exchanging()
      type(pheno_env_t)    :: env
      type(phenology_traits_t)   :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 5. light-driven leaf exchange (flush={}, shed={LIGHT}) --'
      params%flush_cue_mask       = CUE_NONE
      params%shed_cue_mask = CUE_LIGHT
      params%shed_light_hours     = 6.0_wp
      params%shed_light_sharpness = 2.0_wp               ! > 0: many bright hours trigger senescence
      call check_true('leaf-exch: a new cohort has no light memory', state%light_hours_mean < 0.0_wp, &
                      state%light_hours_mean)
      env%par_hours = 2.0_wp
      call phenology_kernel(env, params, 1.0_wp, state, out)
      call check_close('leaf-exch: its first day sets the light memory', state%light_hours_mean, 2.0_wp, 0.0_wp)
      do d = 1_ik, 40_ik
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('leaf-exch: running mean tracks the hours of light', state%light_hours_mean, 2.0_wp, 1.0e-9_wp)
      call check_true('leaf-exch: few bright hours => little senescence', out%leaf_shed_potential < 0.01_wp * params%shed_rate_max)
      call check_close('leaf-exch: flush potential = flush_rate_max', out%leaf_flush_potential, params%flush_rate_max, 1.0e-12_wp)
      do d = 1_ik, 40_ik
         env%par_hours = 10.0_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_true('leaf-exch: bright => senescence', out%leaf_shed_potential > 0.9_wp * params%shed_rate_max)
      call check_close('leaf-exch: flush unchanged by light', out%leaf_flush_potential, params%flush_rate_max, 1.0e-12_wp)
   end subroutine test_light_exchanging
   !--- 6. Rate mapping: the potentials are exactly rate_max times the tendencies.
   subroutine test_rate_mapping()
      type(pheno_env_t)    :: env
      type(phenology_traits_t)   :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 6. rate mapping (potential = rate_max * tendency) --'
      params%flush_cue_mask = CUE_WATER
      params%shed_cue_mask = CUE_LIGHT
      params%shed_light_hours = 6.0_wp
      params%shed_light_sharpness = 0.5_wp
      do d = 1_ik, 25_ik
         env%predawn_leaf_psi = -1.7_wp
         env%par_hours = 6.2_wp                                    ! partial flush + partial shed
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('flush potential == flush_rate_max * flush tendency', out%leaf_flush_potential, &
                       params%flush_rate_max * state%leaf_flush_tendency, 1.0e-12_wp)
      call check_close('shed potential  == shed_rate_max  * shed tendency',  out%leaf_shed_potential,  &
                       params%shed_rate_max  * state%leaf_shed_tendency,  1.0e-12_wp)
      call check_true('tendencies are partial and in [0,1]',                                           &
                      state%leaf_flush_tendency > 0.0_wp .and. state%leaf_flush_tendency < 1.0_wp .and.  &
                      state%leaf_shed_tendency  > 0.0_wp .and. state%leaf_shed_tendency  < 1.0_wp)
   end subroutine test_rate_mapping
   !--- 7. Degenerate drivers: potentials stay finite and non-negative, no FP trap.
   subroutine test_degenerate()
      type(pheno_env_t)    :: env
      type(phenology_traits_t)   :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: ok
      print '(a)', '-- 7. degenerate drivers (every cue on both sides) --'
      params%flush_cue_mask = CUE_ALL
      params%shed_cue_mask = CUE_ALL
      ok = .true.
      do d = 1_ik, 60_ik
         env%doy = modulo(d - 1_ik, 365_ik) + 1_ik
         if (mod(d, 2_ik) == 0_ik) then
            env%temp_day = 350.0_wp
            env%predawn_leaf_psi = 5.0_wp
            env%par_hours = 30.0_wp
         else
            env%temp_day = 200.0_wp
            env%predawn_leaf_psi = -100.0_wp
            env%par_hours = -5.0_wp
         end if
         call phenology_kernel(env, params, 1.0_wp, state, out)
         if (out%leaf_flush_potential < 0.0_wp .or. out%leaf_shed_potential < 0.0_wp) ok = .false.
         if (out%leaf_flush_potential /= out%leaf_flush_potential .or.                            &
             out%leaf_shed_potential  /= out%leaf_shed_potential) ok = .false.        ! NaN
      end do
      call check_true('degenerate: potentials finite and >= 0, no trap', ok)
   end subroutine test_degenerate
   !--- 9. leaf_turnover_step: the one place the tendencies become leaf loss.
   subroutine test_leaf_turnover_step()
      real(wp) :: sen, bg, cap
      real(wp), parameter :: TURN = 0.5_wp            ! [1/yr]
      print '(a)', '-- 8. leaf_turnover_step --'
      !--- Flushing, not senescing: background turnover only, and the flush cap.
      call leaf_turnover_step(1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.06_wp, 0.3_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('flushing: no senescence',             sen, 0.0_wp, 0.0_wp)
      call check_close('flushing: background = turnover/yr',  bg,  TURN / yr_day, 1.0e-15_wp)
      call check_close('flushing: flush cap = rate*full*dt',  cap, 0.06_wp, 1.0e-15_wp)
      !--- Background turnover follows the flush tendency: none in dormancy.
      call leaf_turnover_step(1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.06_wp, 0.3_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('dormant: no background turnover', bg,  0.0_wp, 0.0_wp)
      call check_close('dormant: no flush cap',           cap, 0.0_wp, 0.0_wp)
      !---~---
      !   Senescence = shed_rate_max * tendency * pool * dt; background stops while leaves
      ! drop, so a drop is not also counted as ageing (turnover_gate).
      !---~---
      call leaf_turnover_step(0.5_wp, 1.0_wp, 0.5_wp, 1.0_wp, 0.06_wp, 0.1_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('senescing: senescence = rate*pool*dt', sen, 0.05_wp, 1.0e-15_wp)
      call check_close('senescing: no background turnover',   bg,  0.0_wp, 0.0_wp)
      !--- The gate: the flush tendency while no leaves drop, 0 while they do or in dormancy.
      call check_close('gate: no drop -> flush tendency',      turnover_gate(0.7_wp, 0.0_wp),   0.7_wp, 0.0_wp)
      call check_close('gate: a trace of drop still counts as none', turnover_gate(0.7_wp, 0.005_wp), 0.7_wp, 0.0_wp)
      call check_close('gate: dropping -> 0',                  turnover_gate(1.0_wp, 0.2_wp),   0.0_wp, 0.0_wp)
      call check_close('gate: dormant -> 0',                   turnover_gate(0.0_wp, 0.0_wp),   0.0_wp, 0.0_wp)
      !---~---
      !   The leaf-cover floor: senescence stops at min_leaf_cover of the full canopy,
      ! which is what makes a PFT evergreen (an emergent habit, not a flag).
      !---~---
      call leaf_turnover_step(0.82_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 0.3_wp, 0.0_wp, 0.8_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('floor: senescence stops at min_leaf_cover', 0.82_wp - sen, 0.8_wp, 1.0e-15_wp)
      call leaf_turnover_step(0.7_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 0.3_wp, 0.0_wp, 0.8_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('floor: none below min_leaf_cover', sen, 0.0_wp, 0.0_wp)
      !--- Dormant snap: a dormant canopy left below bare_leaf_cover goes bare.
      call leaf_turnover_step(0.021_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.06_wp, 0.1_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('snap: a dormant canopy goes bare', sen + bg, 0.021_wp, 1.0e-15_wp)
      call leaf_turnover_step(0.021_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('snap: none while flushing', sen, 0.0021_wp, 1.0e-15_wp)
      call leaf_turnover_step(0.021_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.06_wp, 0.1_wp, 0.0_wp, 0.02_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('snap: none when the floor is not below bare', sen, 0.001_wp, 1.0e-15_wp)
      !--- A rate * dt above 1 removes at most the pool.
      call leaf_turnover_step(0.4_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 5.0_wp, 0.0_wp, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('clamp: senescence removes at most the pool', sen, 0.4_wp, 1.0e-15_wp)
   end subroutine test_leaf_turnover_step

   !---~---
   !   9. The C-API shim (ed_c_api_phenology). Only the optional Python library compiles
   ! it, so a field inserted into pheno_args_c or a renamed component could leave the C
   ! API unable to compile while the suite stays green; it is compiled here, so an
   ! interface mismatch fails the build (docs/building.md, Running the tests). The checks
   ! cover what a compile cannot: that the bind(c) mirrors carry values to the kernel and
   ! back rather than reading the wrong member -- the cue sums and the leaf-turnover terms
   ! of docs/science/plant_phenology.md eqs. (3)-(5) and (9). It calls the bind(c)
   ! procedures directly as Fortran, with no .so and no Python.
   !---~---
   subroutine test_c_api_shim()
      type(pheno_env_c)    :: env
      type(pheno_args_c)   :: p
      type(pheno_state_c)  :: st
      type(pheno_out_c)    :: out
      real(c_double), parameter :: DT = 1.0_c_double
      real(c_double) :: sen, bg, cap

      call banner('phenology C-API shim')

      !---~---
      !   1. A no-cue canopy flushes at flush_rate_max and does not senesce.
      ! This pins the plainest path through the shim: the masks and the rate maxima reach
      ! the kernel and the two potentials come back through pheno_out_c.
      !---~---
      env = warm_day()
      p   = base_params(CUE_NONE, CUE_NONE)
      st  = birth_state()
      call ed_phenology_step(env, p, DT, st, out)
      call check_close(out%leaf_flush_potential, p%flush_rate_max, 1.0e-12_wp,                         &
                       'no-cue canopy flushes at flush_rate_max (reaches pheno_out_c)')
      call check_close(out%leaf_shed_potential, 0.0_wp, 1.0e-12_wp, 'no-cue canopy does not senesce')

      !---~---
      !   2. The STATE struct is in/out: the warmth sum must come back changed by exactly
      ! (temp_day - flush_base_temp), so BOTH fields mapped. A shim that passed the state
      ! by value, or unpacked it into the wrong member, would look fine on a single call.
      !---~---
      p  = base_params(CUE_TEMP, CUE_TEMP)
      st = birth_state()
      call ed_phenology_step(env, p, DT, st, out)
      call check_close(st%growing_degree_days, env%temp_day - p%flush_base_temp, 1.0e-9_wp,            &
                       'the warmth sum grows by (temp_day - flush_base_temp)')

      !---~---
      !   3. A COLD day must not add warmth -- the value reaches the kernel, not just the
      ! struct. A mirror that read a neighbouring member would still "work" in test 2.
      !---~---
      env%temp_day = p%flush_base_temp - 5.0_c_double
      st = birth_state()
      call ed_phenology_step(env, p, DT, st, out)
      call check_close(st%growing_degree_days, 0.0_wp, 1.0e-12_wp, 'a sub-base day adds no warmth')

      !--- 4. The water fields: the dry sum grows by (leaf_psi_tlp - predawn_leaf_psi).
      env = warm_day()
      env%predawn_leaf_psi = -3.5_c_double
      p   = base_params(CUE_WATER, CUE_WATER)
      st  = birth_state()
      call ed_phenology_step(env, p, DT, st, out)
      call check_close(st%dry_psi_sum, p%leaf_psi_tlp - env%predawn_leaf_psi, 1.0e-12_wp,               &
                       'the dry sum grows by (leaf_psi_tlp - predawn_leaf_psi)')
      call check_close(st%wet_psi_sum, 0.0_wp, 0.0_wp, 'a dry day adds no wet credit')

      !---~---
      !   4b. The light field and its window: from 0, one step moves the running mean by
      ! par_hours / light_window; with no memory yet (< 0), the first day sets it.
      !---~---
      env = warm_day()
      p   = base_params(CUE_NONE, CUE_LIGHT)
      st  = birth_state()
      call ed_phenology_step(env, p, DT, st, out)
      call check_close(st%light_hours_mean, env%par_hours / p%light_window, 1.0e-12_wp,                &
                       'the light running mean moves by par_hours / light_window in one day')
      st  = birth_state()
      st%light_hours_mean = -1.0_c_double
      call ed_phenology_step(env, p, DT, st, out)
      call check_close(st%light_hours_mean, env%par_hours, 0.0_wp,                                     &
                       'with no light memory, the first day sets it to par_hours')

      !---~---
      !   5. ed_leaf_turnover_step carries every argument through (background turnover,
      ! senescence, the flush cap), each a distinct product of its inputs. Background
      ! turnover acts only while no leaves drop, so it is checked with no shed tendency.
      !---~---
      call ed_leaf_turnover_step(0.5_c_double, 1.0_c_double, 0.5_c_double, 1.0_c_double,              &
                                 0.06_c_double, 0.1_c_double, 0.5_c_double, 0.0_c_double,                &
                                 0.02_c_double, 1.0_c_double, sen, bg, cap)
      call check_close(sen, 0.05_wp, 1.0e-15_wp, 'leaf turnover: senescence = shed_rate_max*tendency*leaf*dt')
      call check_close(bg, 0.0_wp, 0.0_wp, 'leaf turnover: no background while leaves drop')
      call check_close(cap, 0.03_wp, 1.0e-15_wp, 'leaf turnover: flush cap = flush_rate_max*tendency*full*dt')
      call ed_leaf_turnover_step(0.5_c_double, 1.0_c_double, 0.5_c_double, 0.0_c_double,              &
                                 0.06_c_double, 0.1_c_double, 0.5_c_double, 0.0_c_double,                &
                                 0.02_c_double, 1.0_c_double, sen, bg, cap)
      call check_close(bg, 0.5_wp / yr_day * 0.5_wp * 0.5_wp, 1.0e-15_wp,                              &
                       'leaf turnover: background = rate/yr * flush tendency * leaf * dt')

      !---~---
      !   6. ed_daylength passes latitude and day of year through: 12 h at the equator,
      ! a long day at 60 N in June, polar night at 80 N in December.
      !---~---
      call check_close(ed_daylength(0.0_c_double, 172_c_int), 12.0_wp, 0.5_wp, 'daylength: equator ~12 h')
      call check(ed_daylength(60.0_c_double, 172_c_int) > 18.0_wp, 'daylength: 60 N in June > 18 h')
      call check(ed_daylength(80.0_c_double, 355_c_int) < 0.5_wp, 'daylength: 80 N in December ~0 h')
   end subroutine test_c_api_shim


   type(pheno_env_c) function warm_day() result(e)
      e%temp_day         = 293.15_c_double
      e%par_hours        = 14.0_c_double
      e%predawn_leaf_psi = -0.2_c_double
      e%doy              = 150_c_int
      e%hemis_north      = 1_c_int
   end function warm_day

   !--- Parameters with the two cue masks under test and everything else benign.
   type(pheno_args_c)   function base_params(flush_mask, shed_mask) result(q)
      integer(ik), intent(in) :: flush_mask, shed_mask
      q%flush_cue_mask        = int(flush_mask, c_int)
      q%shed_cue_mask         = int(shed_mask,  c_int)
      q%flush_cue_timescale   = 5.0_c_double
      q%shed_cue_timescale    = 5.0_c_double
      q%flush_rate_max        = 0.05_c_double
      q%shed_rate_max         = 0.04_c_double
      q%flush_base_temp       = 278.15_c_double
      q%flush_degree_days     = 100.0_c_double
      q%flush_temp_sharpness  = 0.04_c_double
      q%shed_base_temp        = 290.15_c_double
      q%shed_degree_days      = 50.0_c_double
      q%shed_temp_sharpness   = 0.1_c_double
      q%par_min               = 5.0_c_double
      q%flush_light_hours     = 12.0_c_double
      q%flush_light_sharpness = 1.0_c_double
      q%shed_light_hours      = 11.0_c_double
      q%shed_light_sharpness  = -1.0_c_double
      q%light_window          = 10.0_c_double
      q%leaf_psi_tlp          = -2.0_c_double
      q%flush_water_sum       = 10.0_c_double
      q%flush_water_sharpness = 0.5_c_double
      q%shed_water_sum        = 10.0_c_double
      q%shed_water_sharpness  = 0.5_c_double
   end function base_params

   type(pheno_state_c) function birth_state() result(s)
      s%leaf_flush_tendency = 1.0_c_double
      s%leaf_shed_tendency = 0.0_c_double
      s%growing_degree_days = 0.0_c_double
      s%cold_degree_days   = 0.0_c_double
      s%wet_psi_sum         = 0.0_c_double
      s%dry_psi_sum        = 0.0_c_double
      s%light_hours_mean    = 0.0_c_double
   end function birth_state

end program test_plant_phenology
