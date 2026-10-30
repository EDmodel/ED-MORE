! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_config_check -- the checks a configuration passes before a run starts
! (docs/configuration.md), in the order a load meets them:
!   1. key names (check_config_keys): each key is one its reference lists, else it is
!      reported with the listed key it most likely meant;
!   2. keys a setting makes inapplicable (check_key_use), refused as unknown keys are;
!   3. values (validate_config): at load, after the soil maps fill a polygon's ground,
!      and for each region polygon.
! 1 and 2 are reported with the missing required keys; 3 stops at the first bad value.
!---~---
module ed_config_check
   use ed_kinds,              only : wp, ik
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ed_constants,            only : n_soil_layer_max, day_sec
   use ed_output_params,        only : STATE_DAYS, FC_DAY, FC_MONTH, TIER_MONTHLY, TIER_ANNUAL
   use ed_time,                 only : time_lt, time_valid, time_to_string, whole_years_between
   use ed_toml,                 only : toml_table_t, toml_has, toml_has_section
   use ed_allometry_lib,        only : HEIGHT_CHAVE
   use ed_therm_lib,            only : TEMP_RESPONSE_PEAKED
   use ed_soil_hydr_lib,        only : SOIL_RETENTION_VG
   use ed_pft_params,           only : PATH_C3, PATH_C4
   use ed_soil_params,          only : soil_texture_layers_t, soil_profile_t, ENERGY_BC_DIRICHLET
   use ed_forcing_params,       only : LW_SYNTHESIZE, MET_BACKEND_ERA5LAND, WIND_EXPOSURE_OPEN_TERRAIN,           &
                                       WIND_EXPOSURE_LOCAL, ARCHIVE_DT_SEC, CO2_SOURCE_CONST, CO2_SOURCE_FILE,    &
                                       METAVG_END, MAX_RECYCLE_YEARS
   use ed_radiation_lib,        only : SWPART_PASSTHROUGH
   use ed_region_params,        only : RUN_MODE_REGION
   use ed_maps_params,          only : soil_maps_on
   use ed_config,               only : ed_config_t, site_config, INIT_CENSUS, INIT_RESTART
   use ed_pft_derive,           only : pft_wstress_psi_open, pft_wstress_psi_close
   use ed_soil_derive,          only : soil_layers, soil_layer_textures, build_soil_profile,          &
                                       SOIL_LAYER_TEMPLATE
   use ed_fast_params,          only : INTEG_ESDIRK2, INTEG_RK45, RK45_RTOL
   implicit none
   private

   public :: key_report_t, check_config_keys, check_key_use, nearest_key
   public :: MAIN_KEYS, PFT_KEYS, SITE_GROUND_KEYS
   public :: validate_config

   !--- 1. Key names
   integer, parameter :: KEY_LEN = 64            ! ed_toml's key length
   include 'ed_config_keys.inc'                  ! MAIN_KEYS, PFT_KEYS, from the two references

   !---~---
   !   The [soil] keys a [[sites]] table may set, under the same names: the site's soil column
   ! and its bottom. The ground's surface and its thermal boundary are the run's.
   !---~---
   character(len=*), parameter :: SITE_GROUND_KEYS(15) = [character(len=24) ::                     &
      'depth', 'layer_bottoms', 'retention', 'theta_sat', 'theta_res', 'ksat', 'curve_par_a',        &
      'curve_par_n', 'psi_fc', 'solid_conductivity', 'dry_conductivity', 'dry_heat_capacity',         &
      'water_bottom_bc', 'deep_temp', 'deep_depth']

   !--- The keys the soil maps provide, refused while they are set.
   character(len=*), parameter :: SOIL_MAP_KEYS(15) = [character(len=32) ::                          &
      'soil.depth', 'soil.layer_bottoms', 'soil.theta_sat', 'soil.theta_res', 'soil.ksat',             &
      'soil.curve_par_a', 'soil.curve_par_n', 'soil.solid_conductivity', 'soil.dry_conductivity',       &
      'soil.dry_heat_capacity', 'soil.retention', 'soil.water_bottom_bc', 'soil.deep_temp',             &
      'soil.deep_depth', 'init.soil_theta']

   !--- What the key checks found, one line per key, reported with the missing required keys.
   integer, parameter :: MAX_REPORT = 256
   type :: key_report_t
      integer(ik)        :: n = 0_ik
      character(len=480) :: line(MAX_REPORT) = ''
   end type key_report_t

   character(len=*), parameter :: TAG = 'ed_config: '

contains

   !---~---
   !   Check every key of `t` against `known`, the keys its reference lists (MAIN_KEYS or
   ! PFT_KEYS), adding a line to `report` for each unknown one. A key of an array's table,
   ! `sites.2.depth`, is checked as `sites.depth`; a [[sites]] table also takes
   ! SITE_GROUND_KEYS.
   !---~---
   subroutine check_config_keys(t, known, report)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: known(:)
      type(key_report_t), intent(inout) :: report
      character(len=KEY_LEN) :: hint, key
      character(len=480)     :: line
      integer(ik) :: i
      do i = 1_ik, t%n
         key = table_key(t%key(i))
         if (any(known == key)) cycle
         if (key(1:6) == 'sites.' .and. any(SITE_GROUND_KEYS == key(7:))) cycle
         line = trim(t%key(i))//' (in '//trim(t%source)//') is not a key ED-MORE reads'
         hint = nearest_key(key, known)
         if (len_trim(hint) > 0) line = trim(line)//'; did you mean '//trim(hint)//'?'
         call add_line(report, line)
      end do
   end subroutine check_config_keys

   subroutine add_line(report, line)
      type(key_report_t), intent(inout) :: report
      character(len=*),   intent(in)    :: line
      if (report%n >= MAX_REPORT) return
      report%n = report%n + 1_ik
      report%line(report%n) = line
   end subroutine add_line

   !--- A key without its table's number: `sites.2.depth` is `sites.depth`.
   pure function table_key(key) result(k)
      character(len=*), intent(in) :: key
      character(len=KEY_LEN) :: k
      integer :: d1, d2
      k = key
      d1 = index(key, '.')
      if (d1 == 0) return
      d2 = index(key(d1 + 1:), '.')
      if (d2 <= 1) return
      if (verify(key(d1 + 1:d1 + d2 - 1), '0123456789') == 0) k = key(1:d1)//key(d1 + d2 + 1:)
   end function table_key

   !---~---
   !   The listed key an unknown one most likely meant. First the same name under another
   ! section (a key in the wrong block, like state.cohort_max for output.cohort_max),
   ! the closest in spelling if there are several; else the closest spelling within a
   ! third of the key's length (a typo). Blank when nothing is that close.
   !---~---
   pure function nearest_key(key, known) result(best)
      character(len=*), intent(in) :: key, known(:)
      character(len=KEY_LEN) :: best
      integer :: i, d, d_best
      best = ''
      d_best = huge(1)
      do i = 1, size(known)
         if (leaf(known(i)) /= leaf(key)) cycle
         d = edit_distance(trim(key), trim(known(i)))
         if (d < d_best) then
            d_best = d
            best = known(i)
         end if
      end do
      if (len_trim(best) > 0) return
      d_best = max(2, len_trim(key) / 3) + 1
      do i = 1, size(known)
         d = edit_distance(trim(key), trim(known(i)))
         if (d < d_best) then
            d_best = d
            best = known(i)
         end if
      end do
   end function nearest_key

   !--- A key's name without its section: what follows the last dot.
   pure function leaf(key) result(name)
      character(len=*), intent(in) :: key
      character(len=len(key)) :: name
      name = key(index(trim(key), '.', back=.true.) + 1:)
   end function leaf

   !---~---
   !   The number of single-character insertions, deletions and substitutions that turn a
   ! into b (Levenshtein), by the two-row recurrence.
   !---~---
   pure integer function edit_distance(a, b) result(d)
      character(len=*), intent(in) :: a, b
      integer :: prev(0:len(b)), cur(0:len(b)), i, j
      prev = [(j, j = 0, len(b))]
      do i = 1, len(a)
         cur(0) = i
         do j = 1, len(b)
            cur(j) = min(prev(j) + 1, cur(j-1) + 1, prev(j-1) + merge(0, 1, a(i:i) == b(j:j)))
         end do
         prev = cur
      end do
      d = prev(len(b))
   end function edit_distance

   !--- 2. Keys a setting makes inapplicable

   !---~---
   !   Report every key of the main file `tm` or the PFT file `tp` that the loaded `cfg`
   ! gives nothing to do: it would parse and change nothing.
   !---~---
   subroutine check_key_use(tm, tp, cfg, report)
      type(toml_table_t),  intent(in)    :: tm, tp
      type(ed_config_t),   intent(in)    :: cfg
      type(key_report_t),  intent(inout) :: report
      character(len=*), parameter :: ARCHIVE = 'does not apply to forcing.format = "ERA5land": the '//           &
                                               'ED_ERA5land archive finds the site''s cell and its elevation itself'
      character(len=*), parameter :: CELL    = 'does not apply to [run].mode = "region": each polygon sits at '// &
                                               'its cell''s centre, at the cell''s orography'
      character(len=*), parameter :: CHAVE   = 'applies only to a PFT with height_allometry = "chave2014"'
      integer(ik) :: i
      if (cfg%mode /= RUN_MODE_REGION .and. toml_has_section(tm, 'region'))                       &
         call add_line(report, '[region] (in '//trim(tm%source)//') needs [run].mode = "region"')
      !---~---
      !   [forcing]/[site]: the keys of the other source, location or mode. The CO2 keys apply
      ! with forcing on or off.
      !---~---
      if (cfg%forcing%co2_source == CO2_SOURCE_CONST)                                                 &
         call refuse(report, tm, 'forcing.co2_file', 'needs forcing.co2_source = "file"')
      if (cfg%forcing%co2_source == CO2_SOURCE_FILE)                                                  &
         call refuse(report, tm, 'forcing.co2_const', 'does not apply to forcing.co2_source = "file"')
      if (cfg%forcing%forcing_on) then
         if (cfg%forcing%format == MET_BACKEND_ERA5LAND) then
            call refuse(report, tm, 'forcing.path', ARCHIVE)
            call refuse(report, tm, 'forcing.grid_index', ARCHIVE)
            call refuse(report, tm, 'forcing.grid_match', ARCHIVE)
            call refuse(report, tm, 'site.grid_elevation', ARCHIVE)
         end if
         if (cfg%mode == RUN_MODE_REGION) then
            call refuse(report, tm, 'site.latitude', CELL)
            call refuse(report, tm, 'site.longitude', CELL)
            call refuse(report, tm, 'site.elevation', CELL)
            call refuse(report, tm, 'forcing.max_distance_km', CELL)
         end if
         if (cfg%forcing%wind_exposure == WIND_EXPOSURE_LOCAL) then
            call refuse(report, tm, 'forcing.wind_exposure_z0', 'applies only to wind_exposure = "open_terrain"')
            call refuse(report, tm, 'forcing.wind_blending_height',                                   &
                        'applies only to wind_exposure = "open_terrain"')
         end if
         if (.not. cfg%forcing%apply_elevation_lapse) then
            call refuse(report, tm, 'site.lapse_rate_tair', 'applies only with site.apply_elevation_lapse = true')
            if (cfg%forcing%format /= MET_BACKEND_ERA5LAND)                                          &
               call refuse(report, tm, 'site.grid_elevation',                                         &
                           'applies only with site.apply_elevation_lapse = true')
         end if
      end if
      !--- [maps]: a key the maps provide; a site's own soil.
      if (soil_maps_on(cfg%maps)) then
         do i = 1_ik, size(SOIL_MAP_KEYS, kind=ik)
            call refuse(report, tm, trim(SOIL_MAP_KEYS(i)), 'comes from the [maps] soil maps')
         end do
         if (allocated(cfg%sites)) then
            if (any(cfg%sites%own_ground))                                                            &
               call add_line(report, '[[sites]] [soil] keys (in '//trim(tm%source)//') with the soil '// &
                                     'maps: every site takes its cell''s soil')
         end if
      end if
      !--- [allometry]: the log-quadratic's keys need a PFT on it; its E comes from a map when set.
      if (.not. any(cfg%pft%allom(1:cfg%pft%n)%height_allometry == HEIGHT_CHAVE)) then
         call refuse(report, tp, 'allometry.chave_a', CHAVE)
         call refuse(report, tp, 'allometry.chave_b', CHAVE)
         call refuse(report, tp, 'allometry.chave_c', CHAVE)
         call refuse(report, tp, 'allometry.chave_rse', CHAVE)
         call refuse(report, tp, 'allometry.chave_e', CHAVE)
      else if (len_trim(cfg%maps%chave_e_map) > 0) then
         call refuse(report, tp, 'allometry.chave_e', 'comes from maps.chave_e_map')
      end if
   end subroutine check_key_use

   subroutine refuse(report, t, key, why)
      type(key_report_t), intent(inout) :: report
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key, why
      if (toml_has(t, key)) call add_line(report, key//' (in '//trim(t%source)//') '//why)
   end subroutine refuse

   !--- 3. Values

   !---~---
   !   Validate a configuration; halt on a setting that would corrupt the run. Each site
   ! with its own ground is validated as the run with that ground.
   !---~---
   recursive subroutine validate_config(cfg)
      type(ed_config_t), intent(in) :: cfg
      call check_sites(cfg)
      call check_run(cfg)
      call check_fast(cfg)
      call check_soil(cfg)
      call check_demography(cfg)
      call check_init(cfg)
      call check_output(cfg)
      call check_trait_dynamics(cfg)
      call check_leaf_physiology(cfg)
      call check_forcing(cfg)
      call check_region(cfg)
      call check_maps(cfg)
      call check_pft(cfg)
      call check_phenology(cfg)
      call check_hydraulics(cfg)
   end subroutine validate_config

   !--- [run]: a span, a step and a thread count.
   subroutine check_run(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (.not. time_lt(cfg%start_time, cfg%end_time)) error stop TAG//'end_time must be after start_time'
      if (cfg%dt_slow <= 0.0_wp)                       error stop TAG//'dt_slow <= 0'
      if (cfg%n_threads < 1_ik)                        error stop TAG//'n_threads < 1'
   end subroutine check_run

   !---~---
   !   [fast]: dt_fast is an accuracy setting, measured to 900 s (docs/science/numerical_scheme.md
   ! §5a, §6), and divides dt_slow. The sub-daily probe holds one shared, order-significant CSV,
   ! so it runs single-threaded (§6a).
   !---~---
   subroutine check_fast(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (.not. cfg%fast_biophysics_on) return
      if (cfg%dt_fast <= 0.0_wp) error stop TAG//'dt_fast <= 0'
      if (cfg%dt_fast > 900.0_wp) then
         print '(a)', 'WARNING [ed_config]: dt_fast > 900 s is beyond the measured range. At'
         print '(a)', '  900 s the canopy air is stable and carbon and daily leaf water potential'
         print '(a)', '  converge; above it neither has been measured.'
      end if
      if (cfg%dt_fast > 1800.0_wp) then
         print '(a)', 'WARNING [ed_config]: dt_fast > 1800 s is outside the measured range entirely.'
      end if
      if (cfg%dt_fast > cfg%dt_slow) error stop TAG//'dt_fast > dt_slow'
      if (abs(cfg%dt_slow / cfg%dt_fast - real(nint(cfg%dt_slow / cfg%dt_fast, ik), wp)) > 1.0e-6_wp)   &
         error stop TAG//'dt_slow must be an integer multiple of dt_fast'
      if (cfg%time_integrator /= INTEG_ESDIRK2 .and. cfg%time_integrator /= INTEG_RK45)                    &
         error stop TAG//'time_integrator out of range'
      if (cfg%rtol <= 0.0_wp)       error stop TAG//'fast.rtol <= 0'
      if (cfg%time_integrator == INTEG_RK45 .and. cfg%rtol > RK45_RTOL)                                &
         error stop TAG//'fast.rtol > 1e-4 with time_integrator = "rk45", whose fluxes drift at a looser '// &
                    'rtol (docs/science/numerical_scheme.md §6)'
      if (cfg%atol_scale <= 0.0_wp) error stop TAG//'atol_scale <= 0 (1 = unset)'
      if (cfg%forcing_sample_frac < 0.0_wp .or. cfg%forcing_sample_frac > 1.0_wp)                       &
         error stop TAG//'forcing_sample_frac outside [0,1]'
   end subroutine check_fast

   !---~---
   !   [soil]: the ground's optics, its column and its thermal bottom. Each column fault
   ! gives a silently wrong column rather than a crash (docs/science/soil_biophysics.md).
   !---~---
   subroutine check_soil(cfg)
      type(ed_config_t), intent(in) :: cfg
      type(soil_texture_layers_t) :: tex
      type(soil_profile_t)        :: sp
      real(wp)    :: bottoms(n_soil_layer_max)
      integer(ik) :: n
      if (cfg%soil%ground_albedo_vis < 0.0_wp .or. cfg%soil%ground_albedo_vis >= 1.0_wp .or.          &
          cfg%soil%ground_albedo_nir < 0.0_wp .or. cfg%soil%ground_albedo_nir >= 1.0_wp)              &
         error stop TAG//'soil.ground_albedo_vis and ground_albedo_nir must lie in [0, 1)'
      if (cfg%soil%ground_emissivity <= 0.0_wp .or. cfg%soil%ground_emissivity > 1.0_wp)              &
         error stop TAG//'soil.ground_emissivity must lie in (0, 1]'
      associate (sc => cfg%soil)
         if (sc%depth <= 0.0_wp) error stop TAG//'soil.depth <= 0'
         if (sc%n_layer_bottoms > 0_ik) then
            if (sc%layer_bottoms(1) <= 0.0_wp .or. any(sc%layer_bottoms(2:sc%n_layer_bottoms)          &
                <= sc%layer_bottoms(1:sc%n_layer_bottoms - 1)))                                         &
               error stop TAG//'soil.layer_bottoms must be positive depths [m] that increase'
            if (sc%layer_bottoms(sc%n_layer_bottoms) < sc%depth)                                        &
               error stop TAG//'soil.layer_bottoms must reach soil.depth'
         else if (sc%depth > SOIL_LAYER_TEMPLATE(size(SOIL_LAYER_TEMPLATE))) then
            error stop TAG//'soil.depth is below the deepest default layer (27 m): set soil.layer_bottoms'
         end if
         call soil_layers(sc, bottoms, n)
         if (sc%n_layer_values > 0_ik .and. sc%n_layer_values /= n)                                     &
            error stop TAG//'soil: the per-layer texture does not match the column''s layers'
         call soil_layer_textures(sc, tex)
         !---~---
         !   Each layer check states what must hold, so a NaN fails it: a map layer with no sand
         ! or clay (pure silt, or zero-filled) gives NaN thermal properties.
         !---~---
         if (.not. all(tex%theta_sat(1:n) > tex%theta_res(1:n))) error stop TAG//'soil.theta_sat <= theta_res'
         if (.not. all(tex%theta_res(1:n) >= 0.0_wp))            error stop TAG//'soil.theta_res < 0'
         !---~---
         !   With a texture per layer each layer starts at its field capacity, so init.soil_theta
         ! does not apply.
         !---~---
         if (sc%n_layer_values == 0_ik) then
            if (cfg%init%soil_theta <= sc%theta_res .or. cfg%init%soil_theta > sc%theta_sat)            &
               error stop TAG//'init.soil_theta must lie in (soil.theta_res, theta_sat]'
         end if
         if (.not. all(tex%ksat(1:n) > 0.0_wp)) error stop TAG//'soil.ksat <= 0'
         !---~---
         !   curve_par_a and curve_par_n are van Genuchten's alpha [1/m] > 0 and n > 1, or
         ! Campbell's air-entry head psi_sat [m] < 0 and exponent b > 0.
         !---~---
         if (sc%retention == SOIL_RETENTION_VG) then
            if (.not. all(tex%par_a(1:n) > 0.0_wp))                                                     &
               error stop TAG//'soil.curve_par_a (van Genuchten alpha, 1/m) must be > 0'
            if (.not. all(tex%par_n(1:n) > 1.0_wp))                                                     &
               error stop TAG//'soil.curve_par_n must exceed 1 for van Genuchten'
         else
            if (.not. all(tex%par_a(1:n) < 0.0_wp))                                                     &
               error stop TAG//'soil.curve_par_a (Campbell psi_sat, m) must be < 0'
            if (.not. all(tex%par_n(1:n) > 0.0_wp))                                                     &
               error stop TAG//'soil.curve_par_n (Campbell b) must be > 0'
         end if
         if (sc%psi_fc >= 0.0_wp) error stop TAG//'soil.psi_fc must be negative (a suction head)'
         if (.not. (all(tex%k_solid(1:n) > 0.0_wp) .and. all(tex%k_dry(1:n) > 0.0_wp)))                &
            error stop TAG//'soil conductivities must be positive (a map layer needs sand or clay)'
         if (.not. all(tex%dry_cvol(1:n) > 0.0_wp)) error stop TAG//'soil.dry_heat_capacity <= 0'
      end associate
      !---~---
      !   The Dirichlet thermal bottom needs the site's deep temperature and an anchor below the
      ! bottom node (docs/science/soil_biophysics.md, "The bottom thermal boundary"). With the
      ! soil maps a polygon's anchor comes from its cell, so a configuration the maps have not
      ! filled yet has none to check.
      !---~---
      if (cfg%soil%thermal_bottom_bc == ENERGY_BC_DIRICHLET .and.                                              &
          .not. (soil_maps_on(cfg%maps) .and. cfg%soil%n_layer_values == 0_ik)) then
         if (cfg%soil%deep_temp <= 0.0_wp)                                                             &
            error stop TAG//'soil.thermal_bottom_bc = "dirichlet" requires soil.deep_temp [K] -- '//     &
                       'the mean annual soil temperature below the damping depth. It is a site '//     &
                       'property, like latitude; there is no safe default.'
         if (cfg%soil%deep_temp < 200.0_wp .or. cfg%soil%deep_temp > 330.0_wp)                         &
            error stop TAG//'soil.deep_temp outside 200..330 K -- it is an absolute temperature'
         call build_soil_profile(cfg%soil, sp)
         if (cfg%soil%deep_depth <= abs(sp%z_node(sp%n_active)))                                       &
            error stop TAG//'soil.deep_depth must lie BELOW the bottom soil node '//                    &
                       '(it is measured down from the surface, and the anchor conducts to it)'
      end if
   end subroutine check_soil

   !---~---
   !   [demography]. A recruit must survive its own birth: the pool threshold exceeds the
   ! cull. An operator is turned off by its switch, not by a zero setting.
   !---~---
   subroutine check_demography(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (cfg%demography%max_cohort < 1_ik)                  error stop TAG//'max_cohort < 1'
      if (cfg%demography%max_patch < 1_ik)                   error stop TAG//'max_patch < 1'
      if (cfg%demography%cohort_size_tol_min <= 0.0_wp)                 error stop TAG//'cohort_size_tol_min <= 0'
      if (cfg%demography%cohort_size_tol_max < cfg%demography%cohort_size_tol_min) error stop TAG//'cohort_size_tol_max < min'
      if (cfg%demography%n_cohort_fusion_iter < 1_ik)                   error stop TAG//'n_cohort_fusion_iter < 1'
      if (cfg%demography%n_patch_fusion_iter < 1_ik)                    error stop TAG//'n_patch_fusion_iter < 1'
      if (cfg%demography%patch_light_tol <= 0.0_wp)                     error stop TAG//'patch_light_tol <= 0'
      if (cfg%demography%patch_light_tol_max < cfg%demography%patch_light_tol)                                              &
         error stop TAG//'patch_light_tol_max < patch_light_tol'
      if (cfg%demography%n_height_layers < 2_ik)               error stop TAG//'n_height_layers < 2'
      if (cfg%demography%min_patch_area <= 0.0_wp)             error stop TAG//'min_patch_area <= 0'
      if (cfg%demography%cohort_lai_cap <= 0.0_wp)             error stop TAG//'cohort_lai_cap <= 0'
      if (cfg%demography%growth_memory_days <= 0.0_wp)         error stop TAG//'growth_memory_days <= 0'
      if (cfg%demography%patch_disturbance_on .and. cfg%demography%patch_disturbance_rate <= 0.0_wp)    &
         error stop TAG//'patch_disturbance_rate <= 0: turn the gaps off with patch_disturbance_on'
      if (cfg%demography%disturbance_survive_height <= 0.0_wp) error stop TAG//'disturbance_survive_height <= 0'
      if (cfg%demography%min_recruit_size <= cfg%demography%negligible_nplant)                                              &
         error stop TAG//'min_recruit_size must exceed negligible_nplant'
   end subroutine check_demography

   !--- [init].
   subroutine check_init(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (cfg%init%soil_temp < 233.0_wp .or. cfg%init%soil_temp > 333.0_wp)                           &
         error stop TAG//'init.soil_temp outside 233-333 K'
      if (cfg%init%reacclimate_traits .and. cfg%init%init_mode /= INIT_RESTART)                            &
         error stop TAG//'init.reacclimate_traits applies to a restart (init.init_mode = 2) only'
   end subroutine check_init

   !---~---
   !   [output]: a state file every N days falls on a step's end; a file holds at least one of
   ! its tier's records, so a monthly file spans a month or more and an annual one a year or
   ! more; and a fast-tier window does not straddle two slow steps, whose restructuring would
   ! move the cohort and patch slots it opened with.
   !---~---
   subroutine check_output(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (cfg%output%state_write .and. cfg%output%state_unit == STATE_DAYS) then
         if (mod(real(cfg%output%state_every, wp) * day_sec, cfg%dt_slow) /= 0.0_wp)                 &
            error stop TAG//'output.state_interval in days must be whole slow steps ([run].dt_slow)'
      end if
      if (cfg%output%enabled) then
         if (cfg%output%file_chunk(TIER_MONTHLY) == FC_DAY)                                           &
            error stop TAG//'output.monthly.file_chunk = "day" is shorter than a monthly record: use '//   &
                            '"month", "year" or "all"'
         if (cfg%output%file_chunk(TIER_ANNUAL) == FC_DAY .or. cfg%output%file_chunk(TIER_ANNUAL) == FC_MONTH) &
            error stop TAG//'output.annual.file_chunk is shorter than an annual record: use "year" or "all"'
      end if
      if (.not. (cfg%fast_biophysics_on .and. cfg%output%enabled)) return
      if (cfg%output%fast_interval_steps < 1_ik) error stop TAG//'output.fast_interval_steps < 1'
      if (mod(cfg%n_fast_per_slow, cfg%output%fast_interval_steps) /= 0_ik)                           &
         error stop TAG//'output.fast_interval_steps must divide the fast steps in a slow step '//      &
                         '(dt_slow / dt_fast), so no fast-tier window straddles two slow steps'
   end subroutine check_output

   !--- [trait_dynamics].
   subroutine check_trait_dynamics(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (cfg%trait_dynamics%fineroot_turnover_plastic .and. .not. cfg%trait_dynamics%trait_plasticity_on)    &
         error stop TAG//'trait_dynamics.fineroot_turnover_plastic needs trait_plasticity_on = true'
   end subroutine check_trait_dynamics

   !---~---
   !   [leaf_physiology]: Kc, Ko and Gamma* are denominators; thermal acclimation shifts the
   ! peaked form's entropy term, which the Arrhenius form does not have
   ! (docs/science/leaf_gas_exchange.md eq. (4)).
   !---~---
   subroutine check_leaf_physiology(cfg)
      type(ed_config_t), intent(in) :: cfg
      associate (lp => cfg%leaf_physiology)
         if (lp%kc25 <= 0.0_wp)             error stop TAG//'kc25 <= 0'
         if (lp%ko25 <= 0.0_wp)             error stop TAG//'ko25 <= 0'
         if (lp%gstar25 <= 0.0_wp)          error stop TAG//'gstar25 <= 0'
         if (lp%o2_mol_frac <= 0.0_wp)      error stop TAG//'o2_mol_frac <= 0'
         if (lp%leaf_absorptance <= 0.0_wp) error stop TAG//'leaf_absorptance <= 0'
         if (lp%phi_psii <= 0.0_wp)         error stop TAG//'phi_psii <= 0'
         if (lp%vpd_min <= 0.0_wp)          error stop TAG//'vpd_min <= 0 [kPa]'
         if (lp%thermal_acclimation .and. lp%temp_response_form /= TEMP_RESPONSE_PEAKED)          &
            error stop TAG//'leaf_physiology.thermal_acclimation requires temp_response_form = '//      &
                       '"peaked" (the Arrhenius form has no entropy term to shift)'
         if (lp%thermal_acclimation .and. lp%acclim_window_days <= 0.0_wp)                        &
            error stop TAG//'leaf_physiology.acclim_window_days must be > 0'
      end associate
   end subroutine check_leaf_physiology

   !---~---
   !   [forcing]/[site] (docs/science/forcing.md). The archive is hourly, end-stamped and holds
   ! total shortwave, and its reader loads a month at a time (§1); a recycle window is declared
   ! and spans whole calendar years, or every wrap drifts the hour and the season (§9).
   !---~---
   subroutine check_forcing(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (.not. cfg%forcing%forcing_on) return
      if (cfg%forcing%tq_height <= 0.0_wp)   error stop TAG//'forcing.tq_height must be > 0'
      if (cfg%forcing%wind_height <= 0.0_wp) error stop TAG//'forcing.wind_height must be > 0'
      if (cfg%forcing%wind_exposure == WIND_EXPOSURE_OPEN_TERRAIN) then
         if (cfg%forcing%wind_exposure_z0 <= 0.0_wp .or.                                               &
             cfg%forcing%wind_exposure_z0 >= cfg%forcing%wind_height)                                   &
            error stop TAG//'forcing.wind_exposure_z0 must lie in (0, wind_height)'
         if (cfg%forcing%wind_blending_height < cfg%forcing%wind_height)                                &
            error stop TAG//'forcing.wind_blending_height must be >= wind_height'
      end if
      if (any(cfg%forcing%lapse_rate_tair < 0.0_wp) .or. any(cfg%forcing%lapse_rate_tair > 0.01_wp))     &
         error stop TAG//'site.lapse_rate_tair must lie in [0, 0.01] K/m'
      if (cfg%forcing%lwdown_source == LW_SYNTHESIZE .and. cfg%forcing%lw_cloud_a < 0.0_wp)             &
         error stop TAG//'forcing.lw_cloud_a must be >= 0 (a cloudy sky emits MORE, not less)'
      if (cfg%forcing%format == MET_BACKEND_ERA5LAND) then
         if (len_trim(cfg%forcing%data_path) == 0) error stop TAG//'forcing.data_path is empty'
         if (cfg%mode /= RUN_MODE_REGION .and. cfg%forcing%max_distance_km <= 0.0_wp)               &
            error stop TAG//'forcing.max_distance_km must be > 0'
         if (abs(cfg%forcing%timestep - ARCHIVE_DT_SEC) > 0.5_wp)                                      &
            error stop TAG//'forcing.timestep must be 1 hour for format = "ERA5land"'
         if (cfg%forcing%avg_convention /= METAVG_END)                                                   &
            error stop TAG//'forcing.avg_convention must be "end" for format = "ERA5land"'
         if (cfg%forcing%sw_partition == SWPART_PASSTHROUGH)                                             &
            error stop TAG//'forcing.sw_partition cannot be "passthrough" for format = "ERA5land" '//     &
                            '(the archive stores total shortwave)'
         if (abs(cfg%dt_slow - 86400.0_wp) > 0.5_wp)                                                     &
            error stop TAG//'format = "ERA5land" needs [run].dt_slow = "1d" (the reader loads a month at a time)'
         if (cfg%start_time%hour /= 0_ik .or. cfg%start_time%minute /= 0_ik .or.                         &
             cfg%start_time%second /= 0_ik)                                                             &
            error stop TAG//'format = "ERA5land" needs [run].start_time at 00:00:00'
      end if
      if (cfg%forcing%recycle) then
         if (.not. time_valid(cfg%forcing%recycle_start))                                              &
            error stop TAG//'forcing.recycle_start is not a valid date/time'
         if (.not. time_valid(cfg%forcing%recycle_end))                                                &
            error stop TAG//'forcing.recycle_end is not a valid date/time'
         if (.not. time_lt(cfg%forcing%recycle_start, cfg%forcing%recycle_end))                        &
            error stop TAG//'forcing.recycle_end must be after forcing.recycle_start'
         if (whole_years_between(cfg%forcing%recycle_start, cfg%forcing%recycle_end,                   &
                                 MAX_RECYCLE_YEARS) < 1_ik) then
            write(*,'(4a)') ' ed_config: forcing recycle window ',                                    &
               time_to_string(cfg%forcing%recycle_start), ' .. ', time_to_string(cfg%forcing%recycle_end)
            write(*,'(a)')  '   is not an exact whole number of calendar years. Recycling a'
            write(*,'(a)')  '   partial-year window drifts hour-of-day and day-of-year on every'
            write(*,'(a)')  '   wrap. Set recycle_end to the SAME month/day/time as recycle_start,'
            write(*,'(a)')  '   N years later (it is the exclusive end of the window).'
            error stop TAG//'forcing recycle window is not a whole number of years'
         end if
      end if
   end subroutine check_forcing

   !---~---
   !   [region] (docs/configuration.md, "Regional runs"): every polygon reads its own cell of the
   ! archive and the maps, starts from bare ground or a region state file, and loads a month's
   ! forcing once, so a recycle seam falls on a month boundary.
   !---~---
   subroutine check_region(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (cfg%mode /= RUN_MODE_REGION) return
      if (.not. (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on))                                 &
         error stop TAG//'[run].mode = "region" needs fast.fast_biophysics_on and forcing.forcing_on'
      if (cfg%forcing%format /= MET_BACKEND_ERA5LAND)                                                  &
         error stop TAG//'[run].mode = "region" needs forcing.format = "ERA5land"'
      if (cfg%init%init_mode == INIT_CENSUS)                                                                 &
         error stop TAG//'[run].mode = "region" starts from bare ground or a region state file '//      &
                         '(init.init_mode = 0 or 2), not a census'
      associate (b => cfg%region%box_nwse)
         if (.not. (b(1) > b(3) .and. b(1) <= 90.0_wp .and. b(3) >= -90.0_wp))                         &
            error stop TAG//'region.box_nwse needs -90 <= S < N <= 90 ([N, W, S, E])'
         if (any(abs(b([2, 4])) > 360.0_wp))                                                           &
            error stop TAG//'region.box_nwse longitudes must lie within [-360, 360]'
      end associate
      if (cfg%region%land_fraction_min < 0.0_wp .or. cfg%region%land_fraction_min > 1.0_wp)             &
         error stop TAG//'region.land_fraction_min must lie in [0, 1]'
      if (len_trim(cfg%region%cell_mask) > 0 .and. len_trim(cfg%region%cell_mask_var) == 0)             &
         error stop TAG//'region.cell_mask_var must name the mask''s variable'
      if (.not. ieee_is_finite(cfg%region%cell_mask_min))                                               &
         error stop TAG//'region.cell_mask_min must be a finite number'
      if (.not. soil_maps_on(cfg%maps))                                                                 &
         error stop TAG//'[run].mode = "region" needs the soil maps: maps.soil_depth_map, '//          &
                         'soil_texture_map and soil_chemistry_map'
      if (any(cfg%pft%allom(1:cfg%pft%n)%height_allometry == HEIGHT_CHAVE)                                   &
          .and. len_trim(cfg%maps%chave_e_map) == 0)                                                    &
         error stop TAG//'[run].mode = "region" with allometry.height_allometry = "chave2014" '//       &
                         'needs maps.chave_e_map'
      if (cfg%forcing%recycle .and. time_valid(cfg%forcing%recycle_start)) then
         associate (a => cfg%forcing%recycle_start)
            if (a%day /= 1_ik .or. a%hour > 1_ik .or. a%minute /= 0_ik .or. a%second /= 0_ik)            &
               error stop TAG//'[run].mode = "region" needs forcing.recycle_start at 00:00 or '//        &
                               '01:00 on the 1st of a month'
         end associate
      end if
   end subroutine check_region

   !---~---
   !   [maps]: on the ED_ERA5land grid, read at the forcing cell, so they need the archive; a
   ! chave_e map serves the Chave height allometry only.
   !---~---
   subroutine check_maps(cfg)
      type(ed_config_t), intent(in) :: cfg
      if (soil_maps_on(cfg%maps) .or. len_trim(cfg%maps%chave_e_map) > 0) then
         if (.not. (cfg%forcing%forcing_on .and. cfg%forcing%format == MET_BACKEND_ERA5LAND))           &
            error stop TAG//'the [maps] files need forcing.format = "ERA5land" (they are on its grid)'
      end if
      if (len_trim(cfg%maps%chave_e_map) > 0 .and.                                                      &
          .not. any(cfg%pft%allom(1:cfg%pft%n)%height_allometry == HEIGHT_CHAVE))                            &
         error stop TAG//'maps.chave_e_map needs a PFT with allometry.height_allometry = "chave2014"'
   end subroutine check_maps

   !---~---
   !   [[sites]]: shares of the polygon that sum to 1, a terrain that is a surface, the census
   ! rows of each site in a census run, and a valid ground for each site that has its own. A
   ! region takes one site per cell until a site map gives its sites.
   !---~---
   recursive subroutine check_sites(cfg)
      type(ed_config_t), intent(in) :: cfg
      type(ed_config_t) :: scfg
      integer(ik) :: is, ns
      if (.not. allocated(cfg%sites)) return
      ns = size(cfg%sites, kind=ik)
      if (cfg%mode == RUN_MODE_REGION)                                                              &
         error stop TAG//'a region refuses [[sites]]: each cell is one site until a site map gives its sites'
      if (any(cfg%sites%area_fraction <= 0.0_wp)) error stop TAG//'[[sites]].area_fraction must be > 0'
      if (abs(sum(cfg%sites%area_fraction) - 1.0_wp) > 1.0e-6_wp)                                       &
         error stop TAG//'the [[sites]] area_fraction values must sum to 1'
      if (any(cfg%sites%slope < 0.0_wp .or. cfg%sites%slope >= 90.0_wp))                        &
         error stop TAG//'[[sites]].slope must lie in [0, 90) degrees'
      if (any(cfg%sites%aspect < 0.0_wp .or. cfg%sites%aspect > 360.0_wp))                      &
         error stop TAG//'[[sites]].aspect must lie in [0, 360] degrees'
      if (cfg%init%init_mode == INIT_CENSUS) then
         if (any(cfg%sites%census_site_id <= 0_ik))                                                     &
            error stop TAG//'a census run needs [[sites]].census_site_id > 0 for every site'
         do is = 2_ik, ns
            if (any(cfg%sites(1:is - 1_ik)%census_site_id == cfg%sites(is)%census_site_id))             &
               error stop TAG//'two [[sites]] read the same census_site_id'
         end do
      end if
      do is = 1_ik, ns
         if (.not. cfg%sites(is)%own_ground) cycle
         write(*,'(a,i0,a)') ' config: checking the ground of [[sites]] ', is, ''
         scfg = site_config(cfg, is)
         deallocate(scfg%sites)
         call validate_config(scfg)
      end do
   end subroutine check_sites

   !---~---
   !   [pft]: the photosynthesis curvatures lie in (0, 1], or the smoothing quadratic has no
   ! root in the physical branch; reflectance + transmittance is the single-scatter albedo,
   ! below 1 for a two-stream that conserves energy
   ! (docs/science/canopy_radiation_transfer.md); a leaf turnover rate of 0 would keep
   ! every leaf, and the plasticity slope divides by it.
   !---~---
   subroutine check_pft(cfg)
      type(ed_config_t), intent(in) :: cfg
      integer(ik) :: pf
      if (cfg%pft%n < 1_ik) error stop TAG//'empty PFT table'
      associate (p => cfg%pft)
         do pf = 1_ik, p%n
            if (p%theta_j(pf)     <= 0.0_wp .or. p%theta_j(pf)     > 1.0_wp) error stop TAG//'pft.theta_j outside (0, 1]'
            if (p%theta_cj_c3(pf) <= 0.0_wp .or. p%theta_cj_c3(pf) > 1.0_wp)                            &
               error stop TAG//'pft.theta_cj_c3 outside (0, 1]'
            if (p%theta_ip_c3(pf) <= 0.0_wp .or. p%theta_ip_c3(pf) > 1.0_wp)                            &
               error stop TAG//'pft.theta_ip_c3 outside (0, 1]'
            if (p%theta_cj_c4(pf) <= 0.0_wp .or. p%theta_cj_c4(pf) > 1.0_wp)                            &
               error stop TAG//'pft.theta_cj_c4 outside (0, 1]'
            if (p%theta_ic_c4(pf) <= 0.0_wp .or. p%theta_ic_c4(pf) > 1.0_wp)                            &
               error stop TAG//'pft.theta_ic_c4 outside (0, 1]'
            if (p%leaf_reflect_vis(pf) + p%leaf_transmit_vis(pf) >= 1.0_wp)                             &
               error stop TAG//'pft.leaf_reflect_vis + leaf_transmit_vis >= 1'
            if (p%leaf_reflect_nir(pf) + p%leaf_transmit_nir(pf) >= 1.0_wp)                             &
               error stop TAG//'pft.leaf_reflect_nir + leaf_transmit_nir >= 1'
            if (p%wood_reflect_vis(pf) + p%wood_transmit_vis(pf) >= 1.0_wp)                             &
               error stop TAG//'pft.wood_reflect_vis + wood_transmit_vis >= 1'
            if (p%wood_reflect_nir(pf) + p%wood_transmit_nir(pf) >= 1.0_wp)                             &
               error stop TAG//'pft.wood_reflect_nir + wood_transmit_nir >= 1'
            if (p%leaf_emissivity(pf) <= 0.0_wp .or. p%leaf_emissivity(pf) > 1.0_wp)                    &
               error stop TAG//'pft.leaf_emissivity outside (0,1]'
            if (p%wood_emissivity(pf) <= 0.0_wp .or. p%wood_emissivity(pf) > 1.0_wp)                    &
               error stop TAG//'pft.wood_emissivity outside (0,1]'
            if (p%leaf_clumping(pf) <= 0.0_wp .or. p%leaf_clumping(pf) > 1.0_wp)                        &
               error stop TAG//'pft.leaf_clumping outside (0,1]'
            if (p%wood_clumping(pf) <= 0.0_wp .or. p%wood_clumping(pf) > 1.0_wp)                        &
               error stop TAG//'pft.wood_clumping outside (0,1]'
            if (p%leaf_angle_mean(pf) <= 0.0_wp .or. p%leaf_angle_mean(pf) >= 90.0_wp)                  &
               error stop TAG//'pft.leaf_angle_mean outside (0,90) degrees'
            if (p%leaf_angle_std(pf) < 0.0_wp) error stop TAG//'pft.leaf_angle_std < 0'
            if (pft_wstress_psi_open(p, pf) > 0.0_wp) error stop TAG//'wstress_psi_open must be <= 0'
            if (pft_wstress_psi_close(p, pf) >= pft_wstress_psi_open(p, pf))                            &
               error stop TAG//'wstress_psi_close must be below wstress_psi_open'
         end do
         if (any(p%nbg_init_density < 0.0_wp))       error stop TAG//'pft.nbg_init_density must be >= 0'
         if (any(p%wood_density <= 0.0_wp))          error stop TAG//'wood_density <= 0'
         if (any(p%leaf_turnover_rate_toc <= 0.0_wp)) error stop TAG//'leaf_turnover_rate_toc <= 0'
         if (any(p%stem_resp_sapwood25 < 0.0_wp))    error stop TAG//'stem_resp_sapwood25 < 0'
         if (any(p%max_relative_growth_rate < 0.0_wp)) error stop TAG//'max_relative_growth_rate < 0'
         if (any(p%max_absolute_growth_rate < 0.0_wp)) error stop TAG//'max_absolute_growth_rate < 0'
         if (any(p%photosynthetic_pathway /= PATH_C3 .and. p%photosynthetic_pathway /= PATH_C4))         &
            error stop TAG//'photosynthetic_pathway not in {1,2}'
         if (any(p%vcmax25 <= 0.0_wp))          error stop TAG//'vcmax25 <= 0'
         if (any(p%jmax_vcmax_ratio <= 0.0_wp)) error stop TAG//'jmax_vcmax_ratio <= 0'
         if (any(p%tpu_vcmax_ratio <= 0.0_wp))  error stop TAG//'tpu_vcmax_ratio <= 0'
         if (any(p%rd_vcmax_ratio < 0.0_wp))    error stop TAG//'rd_vcmax_ratio < 0'
         if (any(p%stomatal_g0 < 0.0_wp))       error stop TAG//'stomatal_g0 < 0'
         if (any(p%stomatal_g1 < 0.0_wp))       error stop TAG//'stomatal_g1 < 0'
         if (any(p%stomatal_d0 <= 0.0_wp))      error stop TAG//'stomatal_d0 <= 0 (Leuning divides by it)'
         if (any(p%katul_lambda25 <= 0.0_wp))   error stop TAG//'katul_lambda25 <= 0'
         if (any(p%wstress_lambda_exp < 0.0_wp .or. p%wstress_lambda_exp > 8.0_wp))                      &
            error stop TAG//'wstress_lambda_exp must be in [0,8]'
         if (any(p%wstress_sref_stomata <= 0.0_wp))                                                     &
            error stop TAG//'wstress_sref_stomata must be > 0 (beta_stomata = exp(sref*(psi - psi_onset)))'
         if (any(p%stomata_psi_onset > 0.0_wp))                                                         &
            error stop TAG//'stomata_psi_onset must be <= 0 (a predawn leaf water potential)'
         if (any(p%leaf_surf_water_max < 0.0_wp) .or. any(p%wood_surf_water_max < 0.0_wp))               &
            error stop TAG//'leaf_surf_water_max and wood_surf_water_max must be >= 0 [kg/m2 of surface]'
         if (any(p%photosynthetic_pathway == PATH_C3 .and. (p%theta_j <= 0.0_wp .or. p%theta_j >= 1.0_wp))) &
            error stop TAG//'C3 theta_j must be in (0,1)'
         if (any(p%photosynthetic_pathway == PATH_C4 .and. p%kp25 <= 0.0_wp))                            &
            error stop TAG//'C4 PFT needs kp25 > 0'
         if (any(p%photosynthetic_pathway == PATH_C4 .and. p%quantum_yield_c4 <= 0.0_wp))                &
            error stop TAG//'C4 PFT needs quantum_yield_c4 > 0'
         if (any(p%photosynthetic_pathway == PATH_C4 .and.                                               &
                 (p%theta_cj_c4 <= 0.0_wp .or. p%theta_cj_c4 >= 1.0_wp)))                                &
            error stop TAG//'C4 theta_cj_c4 must be in (0,1)'
         if (any(p%photosynthetic_pathway == PATH_C4 .and.                                               &
                 (p%theta_ic_c4 <= 0.0_wp .or. p%theta_ic_c4 >= 1.0_wp)))                                &
            error stop TAG//'C4 theta_ic_c4 must be in (0,1)'
      end associate
   end subroutine check_pft

   !---~---
   !   [phenology]: each cue mask is a set of TEMP(1), LIGHT(2), WATER(4); a sharpness of 0
   ! would make its switch a constant 0.5 (docs/science/plant_phenology.md eq. (2)).
   !---~---
   subroutine check_phenology(cfg)
      type(ed_config_t), intent(in) :: cfg
      associate (p => cfg%pft, n => cfg%pft%n)
         if (any(p%pheno(1:n)%flush_cue_mask < 0_ik .or. p%pheno(1:n)%flush_cue_mask > 7_ik) .or.        &
             any(p%pheno(1:n)%shed_cue_mask  < 0_ik .or. p%pheno(1:n)%shed_cue_mask  > 7_ik))            &
            error stop TAG//'phenology flush_cue_mask / shed_cue_mask out of range [0,7]'
         if (any(p%pheno(1:n)%flush_cue_timescale <= 0.0_wp) .or.                                       &
             any(p%pheno(1:n)%shed_cue_timescale  <= 0.0_wp) .or.                                       &
             any(p%pheno(1:n)%light_window        <= 0.0_wp))                                           &
            error stop TAG//'phenology cue timescales and light_window must be > 0'
         if (any(p%pheno(1:n)%par_min < 0.0_wp))        error stop TAG//'phenology.par_min must be >= 0'
         if (any(p%pheno(1:n)%flush_rate_max <= 0.0_wp)) error stop TAG//'phenology.flush_rate_max must be > 0'
         if (any(p%pheno(1:n)%shed_rate_max < 0.0_wp))   error stop TAG//'phenology.shed_rate_max must be >= 0'
         if (any(p%pheno(1:n)%min_leaf_cover < 0.0_wp .or. p%pheno(1:n)%min_leaf_cover >= 1.0_wp) .or.   &
             any(p%pheno(1:n)%bare_leaf_cover < 0.0_wp .or. p%pheno(1:n)%bare_leaf_cover >= 1.0_wp))     &
            error stop TAG//'phenology min_leaf_cover / bare_leaf_cover must be in [0,1)'
         if (any(p%pheno(1:n)%flush_temp_sharpness == 0.0_wp) .or.                                     &
             any(p%pheno(1:n)%shed_temp_sharpness  == 0.0_wp) .or.                                     &
             any(p%pheno(1:n)%flush_light_sharpness == 0.0_wp) .or.                                    &
             any(p%pheno(1:n)%shed_light_sharpness  == 0.0_wp) .or.                                    &
             any(p%pheno(1:n)%flush_water_sharpness == 0.0_wp) .or.                                    &
             any(p%pheno(1:n)%shed_water_sharpness  == 0.0_wp))                                        &
            error stop TAG//'a phenology sharpness is 0: its switch would be a constant 0.5'
      end associate
   end subroutine check_phenology

   !---~---
   !   [hydraulics]: the rooting traits set the root profile over the soil column, the sapwood
   ! conductance needs a positive conductivity and path factor, and each tissue's
   ! pressure-volume curve and the xylem's vulnerability curve need their signs: with
   ! elastic_mod <= |pi0| the turgor-loss point pi0 eps / (pi0 + eps) is not negative
   ! (docs/science/plant_hydraulics.md).
   !---~---
   subroutine check_hydraulics(cfg)
      type(ed_config_t), intent(in) :: cfg
      integer(ik) :: ipft
      associate (hydro => cfg%pft%hydro(1:cfg%pft%n))
         if (any(hydro%root_beta <= 0.0_wp .or. hydro%root_beta >= 1.0_wp))                             &
            error stop TAG//'hydraulics.root_beta must lie in (0, 1)'
         if (any(hydro%root_depth_b1 <= 0.0_wp)) error stop TAG//'hydraulics.root_depth_b1 must be > 0'
         if (any(hydro%root_depth_b2 < 0.0_wp))  error stop TAG//'hydraulics.root_depth_b2 must be >= 0'
         if (any(hydro%specific_root_length <= 0.0_wp))                                                 &
            error stop TAG//'hydraulics.specific_root_length must be > 0'
         if (any(hydro%fine_root_radius <= 0.0_wp)) error stop TAG//'hydraulics.fine_root_radius must be > 0'
         if (any(hydro%vessel_curl <= 0.0_wp))      error stop TAG//'hydraulics.vessel_curl must be > 0'
         if (any(hydro%wood_kmax <= 0.0_wp))        error stop TAG//'hydraulics.wood_kmax must be > 0'
         !---~---
         !   The curves one PFT at a time: an array expression over hydro%leaf_curve%<key>,
         ! a key of a record inside each element of the PFT array, reads other PFTs' keys
         ! in ifx 2026.1's AVX-512 code (the production build) and refused every multi-PFT
         ! file.
         !---~---
         do ipft = 1_ik, cfg%pft%n
            associate (leaf => cfg%pft%hydro(ipft)%leaf_curve, wood => cfg%pft%hydro(ipft)%wood_curve)
               if (.not. (leaf%pi0 < 0.0_wp .and. wood%pi0 < 0.0_wp))                                   &
                  error stop TAG//'hydraulics.leaf_pi0 and wood_pi0 must be < 0'
               if (.not. (leaf%elastic_mod > -leaf%pi0 .and. wood%elastic_mod > -wood%pi0))             &
                  error stop TAG//'hydraulics.leaf_elastic_mod and wood_elastic_mod must exceed |pi0|'
               if (.not. (leaf%apoplast_frac >= 0.0_wp .and. leaf%apoplast_frac < 1.0_wp .and.          &
                          wood%apoplast_frac >= 0.0_wp .and. wood%apoplast_frac < 1.0_wp))              &
                  error stop TAG//'hydraulics.leaf_apoplast_frac and wood_apoplast_frac must lie in [0, 1)'
               if (.not. (leaf%water_sat > 0.0_wp .and. wood%water_sat > 0.0_wp))                       &
                  error stop TAG//'hydraulics.leaf_water_sat and wood_water_sat must be > 0'
            end associate
         end do
         if (.not. all(hydro%wood_psi50 < 0.0_wp)) error stop TAG//'hydraulics.wood_psi50 must be < 0'
         if (.not. all(hydro%wood_kexp > 0.0_wp))  error stop TAG//'hydraulics.wood_kexp must be > 0'
      end associate
   end subroutine check_hydraulics

end module ed_config_check
