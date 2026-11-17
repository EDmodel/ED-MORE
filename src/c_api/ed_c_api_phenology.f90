! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_c_api_phenology -- the C-API shim for the PHENOLOGY kernel (`edmore.plant.pheno`).
!
!   One shim per subsystem, mirroring the Fortran tree (src/README.md). The bind(c) structs are
! an ABI contract with `python/edmore/plant/pheno.py` and are compiled by a mandatory ctest
! target (`test_plant_phenology`) so a field-order change breaks the BUILD -- see the header
! of ed_c_api_leaf for why.
!---~---
module ed_c_api_phenology
   use iso_c_binding,        only : c_double, c_int
   use ed_kinds,             only : wp, ik
   use ed_plant_phenology,   only : pheno_env_t, pheno_state_t, pheno_out_t, phenology_kernel,  &
                                    leaf_turnover_step
   use ed_pft_params,        only : phenology_traits_t
   use ed_time,              only : daylength
   implicit none
   private

   public :: pheno_env_c, pheno_args_c, pheno_state_c, pheno_out_c
   public :: ed_phenology_step, ed_phenology_run, ed_leaf_turnover_step, ed_daylength

   !--- C mirror of pheno_env_t (3 doubles + 2 ints; hemis_north 0/1).
   type, bind(c) :: pheno_env_c
      real(c_double) :: temp_day, par_hours, predawn_leaf_psi
      integer(c_int) :: doy, hemis_north
   end type pheno_env_c

   !--- C mirror of the kernel's PFT record, phenology_traits_t, less its leaf-cover pair.
   type, bind(c) :: pheno_args_c
      integer(c_int) :: flush_cue_mask, shed_cue_mask
      real(c_double) :: flush_cue_timescale, shed_cue_timescale, flush_rate_max, shed_rate_max
      real(c_double) :: flush_base_temp, flush_degree_days, flush_temp_sharpness
      real(c_double) :: shed_base_temp, shed_degree_days, shed_temp_sharpness
      real(c_double) :: par_min, flush_light_hours, flush_light_sharpness
      real(c_double) :: shed_light_hours, shed_light_sharpness, light_window
      real(c_double) :: leaf_psi_tlp, flush_water_sum, flush_water_sharpness
      real(c_double) :: shed_water_sum, shed_water_sharpness
   end type pheno_args_c

   !--- C mirror of pheno_state_t (7 doubles; the prognostic memory, in/out).
   type, bind(c) :: pheno_state_c
      real(c_double) :: leaf_flush_tendency, leaf_shed_tendency, growing_degree_days,     &
                        cold_degree_days, wet_psi_sum, dry_psi_sum, light_hours_mean
   end type pheno_state_c

   !--- C mirror of pheno_out_t (2 doubles).
   type, bind(c) :: pheno_out_c
      real(c_double) :: leaf_flush_potential, leaf_shed_potential
   end type pheno_out_c

contains

   !---~---
   !   Advance one cohort's phenology ONE step: unpack the C structs, call phenology_kernel
   ! (docs/science/plant_phenology.md §1-§4; state advanced in place), pack the state and the
   ! potential rates back out. `state_c` is intent(inout): the caller keeps it across days
   ! (the phenological memory).
   !---~---
   subroutine ed_phenology_step(env_c, p_c, dt, state_c, out_c) bind(c, name="ed_phenology_step")
      type(pheno_env_c),    intent(in)    :: env_c
      type(pheno_args_c),   intent(in)    :: p_c
      real(c_double), value, intent(in)   :: dt
      type(pheno_state_c),  intent(inout) :: state_c
      type(pheno_out_c),    intent(out)   :: out_c
      type(pheno_env_t)    :: env
      type(phenology_traits_t) :: p
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out

      env%temp_day         = env_c%temp_day
      env%par_hours        = env_c%par_hours
      env%predawn_leaf_psi = env_c%predawn_leaf_psi
      env%doy              = int(env_c%doy, ik)
      env%hemis_north      = env_c%hemis_north /= 0_c_int
      p = traits_from_c(p_c)
      call state_from_c(state_c, state)

      call phenology_kernel(env, p, real(dt, wp), state, out)

      call state_to_c(state, state_c)
      out_c%leaf_flush_potential  = out%leaf_flush_potential
      out_c%leaf_shed_potential   = out%leaf_shed_potential
   end subroutine ed_phenology_step

   !---~---
   !   The kernel and the leaf rule over many columns and days in one call: each column a cohort
   ! (or a site) with its own daily drivers, all of one PFT's traits. Each day is what
   ! edmore.plant.pheno.Phenology.step does: the kernel (docs/science/plant_phenology.md §1-§4),
   ! then the leaf rule on the leaf cover (§5, eq. (9)) with carbon never limiting the flush. The
   ! drivers, the leaf turnover rate among them (a cohort's, which follows its light), are
   ! (n_day, n_col) in Fortran order, a C (n_col, n_day) array; `state_c` and `cover`
   ! hold each column's memory and leaf cover in and out, so a run can continue in a later call.
   ! `cover_out` is the leaf cover after each day, `loss_out` the cover lost that day (senescence
   ! plus background turnover, before resorption).
   !---~---
   subroutine ed_phenology_run(n_col, n_day, p_c, min_leaf_cover, bare_leaf_cover, dt, hemis_north, doy,  &
                               temp, hours, psi, turnover, state_c, cover, cover_out, loss_out)           &
                               bind(c, name="ed_phenology_run")
      integer(c_int), value, intent(in)  :: n_col, n_day, hemis_north
      type(pheno_args_c),    intent(in)  :: p_c
      real(c_double), value, intent(in)  :: min_leaf_cover, bare_leaf_cover, dt
      integer(c_int),        intent(in)  :: doy(n_day)
      real(c_double),        intent(in)  :: temp(n_day, n_col), hours(n_day, n_col), psi(n_day, n_col)
      real(c_double),        intent(in)  :: turnover(n_day, n_col)   ! [1/yr] leaf turnover rate
      type(pheno_state_c),   intent(inout) :: state_c(n_col)
      real(c_double),        intent(inout) :: cover(n_col)
      real(c_double),        intent(out) :: cover_out(n_day, n_col), loss_out(n_day, n_col)
      type(phenology_traits_t) :: p
      type(pheno_env_t)    :: env
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      real(wp)    :: c, sen, bg, cap, post
      integer(ik) :: ic, id

      p = traits_from_c(p_c)
      env%hemis_north = hemis_north /= 0_c_int
      do ic = 1_ik, int(n_col, ik)
         call state_from_c(state_c(ic), state)
         c = cover(ic)
         do id = 1_ik, int(n_day, ik)
            env%temp_day         = temp(id, ic)
            env%par_hours        = hours(id, ic)
            env%predawn_leaf_psi = psi(id, ic)
            env%doy              = int(doy(id), ik)
            call phenology_kernel(env, p, real(dt, wp), state, out)
            call leaf_turnover_step(c, 1.0_wp, state%leaf_flush_tendency, state%leaf_shed_tendency,   &
                                    p%flush_rate_max, p%shed_rate_max, real(turnover(id, ic), wp),    &
                                    real(min_leaf_cover, wp), real(bare_leaf_cover, wp), real(dt, wp),  &
                                    sen, bg, cap)
            post = c - sen - bg
            c = post + min(max(0.0_wp, 1.0_wp - post), cap)
            cover_out(id, ic)  = c
            loss_out(id, ic) = sen + bg
         end do
         call state_to_c(state, state_c(ic))
         cover(ic) = c
      end do
   end subroutine ed_phenology_run

   !--- The kernel's PFT record from its C mirror.
   function traits_from_c(p_c) result(p)
      type(pheno_args_c), intent(in) :: p_c
      type(phenology_traits_t)       :: p
      p%flush_cue_mask        = int(p_c%flush_cue_mask, ik)
      p%shed_cue_mask         = int(p_c%shed_cue_mask, ik)
      p%flush_cue_timescale   = p_c%flush_cue_timescale
      p%shed_cue_timescale    = p_c%shed_cue_timescale
      p%flush_rate_max        = p_c%flush_rate_max
      p%shed_rate_max         = p_c%shed_rate_max
      p%flush_base_temp       = p_c%flush_base_temp
      p%flush_degree_days     = p_c%flush_degree_days
      p%flush_temp_sharpness  = p_c%flush_temp_sharpness
      p%shed_base_temp        = p_c%shed_base_temp
      p%shed_degree_days      = p_c%shed_degree_days
      p%shed_temp_sharpness   = p_c%shed_temp_sharpness
      p%par_min               = p_c%par_min
      p%flush_light_hours     = p_c%flush_light_hours
      p%flush_light_sharpness = p_c%flush_light_sharpness
      p%shed_light_hours      = p_c%shed_light_hours
      p%shed_light_sharpness  = p_c%shed_light_sharpness
      p%light_window          = p_c%light_window
      p%leaf_psi_tlp          = p_c%leaf_psi_tlp
      p%flush_water_sum       = p_c%flush_water_sum
      p%flush_water_sharpness = p_c%flush_water_sharpness
      p%shed_water_sum        = p_c%shed_water_sum
      p%shed_water_sharpness  = p_c%shed_water_sharpness
   end function traits_from_c

   !--- The prognostic memory from and to its C mirror.
   subroutine state_from_c(state_c, state)
      type(pheno_state_c), intent(in)  :: state_c
      type(pheno_state_t), intent(out) :: state
      state%leaf_flush_tendency = state_c%leaf_flush_tendency
      state%leaf_shed_tendency  = state_c%leaf_shed_tendency
      state%growing_degree_days = state_c%growing_degree_days
      state%cold_degree_days    = state_c%cold_degree_days
      state%wet_psi_sum         = state_c%wet_psi_sum
      state%dry_psi_sum         = state_c%dry_psi_sum
      state%light_hours_mean    = state_c%light_hours_mean
   end subroutine state_from_c

   subroutine state_to_c(state, state_c)
      type(pheno_state_t), intent(in)  :: state
      type(pheno_state_c), intent(out) :: state_c
      state_c%leaf_flush_tendency = state%leaf_flush_tendency
      state_c%leaf_shed_tendency  = state%leaf_shed_tendency
      state_c%growing_degree_days = state%growing_degree_days
      state_c%cold_degree_days    = state%cold_degree_days
      state_c%wet_psi_sum         = state%wet_psi_sum
      state_c%dry_psi_sum         = state%dry_psi_sum
      state_c%light_hours_mean    = state%light_hours_mean
   end subroutine state_to_c

   !---~---
   !   One step of leaf loss and the flush cap (ed_plant_phenology::leaf_turnover_step,
   ! docs/science/plant_phenology.md §5, eq. (9)), so Python applies the carbon layer's own
   ! rule rather than a copy of it.
   !---~---
   subroutine ed_leaf_turnover_step(leaf, leaf_full, leaf_flush_tendency, leaf_shed_tendency,    &
                                    flush_rate_max, shed_rate_max, leaf_turnover_rate,           &
                                    min_leaf_cover, bare_leaf_cover, dt,                         &
                                    senescence, background, flush_cap)                           &
                                    bind(c, name="ed_leaf_turnover_step")
      real(c_double), value, intent(in) :: leaf, leaf_full, leaf_flush_tendency, leaf_shed_tendency
      real(c_double), value, intent(in) :: flush_rate_max, shed_rate_max, leaf_turnover_rate
      real(c_double), value, intent(in) :: min_leaf_cover, bare_leaf_cover, dt
      real(c_double), intent(out)       :: senescence, background, flush_cap
      real(wp) :: sen, bg, cap
      call leaf_turnover_step(real(leaf, wp), real(leaf_full, wp), real(leaf_flush_tendency, wp),     &
                              real(leaf_shed_tendency, wp), real(flush_rate_max, wp),                 &
                              real(shed_rate_max, wp), real(leaf_turnover_rate, wp),                  &
                              real(min_leaf_cover, wp), real(bare_leaf_cover, wp), real(dt, wp),      &
                              sen, bg, cap)
      senescence = sen
      background = bg
      flush_cap = cap
   end subroutine ed_leaf_turnover_step

   !---~---
   !   Day length [h] at a latitude and day of year (ed_time::daylength), the light cue's
   ! driver, so Python builds its forcing with the model's own formula.
   !---~---
   real(c_double) function ed_daylength(latitude_deg, doy) bind(c, name="ed_daylength")
      real(c_double), value, intent(in) :: latitude_deg
      integer(c_int), value, intent(in) :: doy
      ed_daylength = daylength(real(latitude_deg, wp), int(doy, ik))
   end function ed_daylength

end module ed_c_api_phenology
