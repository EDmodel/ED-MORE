! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_output_roundtrip -- the netCDF serializer end to end: build a manager, walk a 3-day        !
! calendar folding site state, flush per-stream files, then RE-OPEN them and assert record counts,  !
! averaged values, calendar, and CF layout (docs/science/diagnostics.md §3 and §4).                !
!==========================================================================================!
module test_ro_support
   use ed_kinds,            only : wp, ik
   use ed_polygon_state_vars, only : polygon_t
   implicit none
contains
   !----- Opaque per-cohort agb mutation (production mutates through veg dynamics; see the        !
   !      test_output header note on the nvfortran -O2 aliasing hazard).                          !
   subroutine set_site_agb(poly, a)
      type(polygon_t), intent(inout) :: poly
      real(wp),     intent(in)    :: a
      poly%cohort%agb(1) = a
   end subroutine set_site_agb

   !----- Opaque setters for the flux/state sources (fast-loop accumulators + soil column). -------!
   subroutine set_site_gpp(poly, g)
      type(polygon_t), intent(inout) :: poly
      real(wp),     intent(in)    :: g
      poly%cohort%gpp_accum(1) = g
   end subroutine set_site_gpp

   subroutine set_site_soil_temp(poly, t)
      type(polygon_t), intent(inout) :: poly
      real(wp),     intent(in)    :: t
      poly%patch%soil_energy(1)%soil_temp(:) = t
   end subroutine set_site_soil_temp
end module test_ro_support

program test_output_roundtrip
   use iso_c_binding,         only : c_int, c_size_t, c_double
   use ed_kinds,              only : wp, ik
   use ed_config,             only : ed_config_t
   use ed_soil_derive,        only : build_profile_hydraulics
   use ed_time,               only : ed_time_t, time_advance_days
   use ed_polygon_state_vars, only : polygon_t, alloc_polygon, free_polygon, rebuild_patch_map
   use ed_output_params,      only : FC_DAY, FC_MONTH, FC_YEAR, FC_RUN, SYNC_FLUSH
   use ed_output_types,       only : output_files_t, output_buffers_t
   use ed_output_setup,       only : alloc_output_set
   use ed_output_registry,    only : selection_t, OVR_FALSE, OVR_MASK
   use ed_output_params,      only : FREQ_NONE, FREQ_DAILY, FREQ_ANNUAL
   use ed_output_integrate,   only : output_integrate
   use ed_output_stream,   only : output_serialize, output_close
   use ed_netcdf_c
   use ed_test_soil_grid, only : geometric_layer_bottoms
   use ed_constants,          only : n_soil_layer_max, kgc_2_umol, day_sec
   use ed_soil_params,        only : soil_profile_t
   use ed_soil_hydr_lib,      only : SOIL_RETENTION_VG
   use ed_output_types,       only : MISSING_VALUE
   use ed_test_support, only : banner, check, check_close, test_pft_table
   use test_ro_support,       only : set_site_agb, set_site_gpp, set_site_soil_temp
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks are computed on

   type(ed_config_t)      :: cfg
   type(polygon_t)        :: poly
   type(output_files_t)   :: files
   type(output_buffers_t) :: bufs(1)     !< one site's buffers: an array of one, as the driver holds them
   type(ed_time_t)        :: prev, now
   integer(ik)            :: iday
   real(wp)               :: dt
   type(soil_profile_t)   :: soilp
   integer(ik), parameter :: N_ACTIVE = 6_ik

   call banner('output_roundtrip')

   !----- A fixed 1-cohort / 1-patch site (constant slot set). -----!
   call alloc_polygon(poly, 1_ik, 16_ik, 4_ik, 4_ik, test_pft_table())
   poly%cohort%n = 1_ik
   poly%patch%n = 1_ik
   poly%cohort%nplant(1) = 1.0_wp
   poly%cohort%dbh(1) = 5.0_wp
   poly%cohort%height(1) = 4.0_wp
   poly%cohort%basal_area(1) = 19.6_wp
   poly%cohort%leaf_area(1) = 2.0_wp
   poly%cohort%growth_avg(1) = 0.1_wp
   poly%cohort%pft(1) = 1_ik
   poly%cohort%owner_patch(1) = 1_ik
   poly%cohort%global_id(1) = 1_ik
   poly%patch%area(1) = 1.0_wp
   poly%patch%cohort_offset(1) = 1_ik
   poly%patch%cohort_count(1) = 1_ik
   poly%patch%age(1) = 0.0_wp
   poly%patch%dist_type(1) = 1_ik
   poly%patch%global_id(1) = 1_ik
   call rebuild_patch_map(poly)

   !----- Constant carbon + soil sources (fast loop is off in this unit test, so set them here    !
   !      via opaque setters, as the fast loop would). gpp_accum=2 kgC/m2 a day -> daily mean_gpp_py,  !
   !      its mean rate, 2 kgC/m2/day in umol/m2/s;                                                 !
   !      soil_temp=290 -> daily AGG_TMEAN mean_soil_temp_layer_si=290 per layer.                            !
   call set_site_gpp(poly, 2.0_wp)
   call set_site_soil_temp(poly, 290.0_wp)

   !----- Output config: daily + annual streams; the variables the checks read. -----!
   cfg = build_cfg()
   dt  = 86400.0_wp
   call alloc_output_set(files, bufs(1), cfg, [selection_t('*', OVR_FALSE, FREQ_NONE),                      &
                         selection_t('mean_agb_*', OVR_MASK, ior(FREQ_DAILY, FREQ_ANNUAL)),                    &
                         selection_t('mean_gpp_py', OVR_MASK, ior(FREQ_DAILY, FREQ_ANNUAL)),                   &
                         selection_t('mean_soil_temp_layer_si', OVR_MASK, ior(FREQ_DAILY, FREQ_ANNUAL))])

   !----- A soil column with FEWER active layers than the compile-time ceiling. That gap is the   !
   !      whole point: a soil-dimensioned variable emitted over all n_soil_layer_max slots would  !
   !      write the inactive tail as DATA (0 for a state, NaN for a derived quantity) instead of  !
   !      leaving it for netCDF to fill. Reducing over the axis then gives nonsense -- a 136 K    !
   !      soil column on a real run -- and soil_matric_potential evaluated on padding whose       !
   !      theta_sat is 0 aborts a -fpe0 build. ---------------------------------------------------!
   call geometric_layer_bottoms(2.0_wp, N_ACTIVE, 3.0_wp, test_bottoms)
   call build_profile_hydraulics(N_ACTIVE, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,   &
                               2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soilp)
   !----- The site's soil column gives the extraction its layers and curve; the polygon's row of  !
   !      the file set gives the files their layer depths, as attach_output does in a run. -------!
   poly%site%surface(1)%soil_profile = soilp
   bufs(1)%diag%soil_ready = .true.
   files%site_n_soil(1, 1)       = soilp%n_active
   files%site_soil_z(:, 1, 1)    = soilp%z_node
   files%site_soil_face(:, 1, 1) = soilp%soil_layer_z
   call check(N_ACTIVE < n_soil_layer_max, 'fixture must have an inactive tail to be worth anything')

   !----- Walk 3 daily steps of 2000-01, as polygon_step calls the tick: the step from `prev` to  !
   !      `now` folds into the day that holds `prev`, then closes it, because `now` left that day.  !
   !      agb 10 / 20 / 30 are the states at the end of the steps of 1, 2 and 3 January. ---------!
   now = ed_time_t(year=2000_ik, month=1_ik, day=1_ik)
   do iday = 1_ik, 3_ik
      prev = now
      now = time_advance_days(prev, 1_ik)
      call set_site_agb(poly, 10.0_wp * real(iday, wp))
      call output_integrate(files, bufs(1), poly, prev, dt, .true., .false., .false.)
      call output_serialize(files, bufs)
   end do

   call output_close(files, bufs)              ! flush the 2000 annual partial

   !----- Re-open and assert. -----!
   call check_daily('test_ro-D-200001.nc')
   call check_annual('test_ro-Y.nc')

   call free_polygon(poly)
   write(*,'(a)') 'test_output_roundtrip: ALL PASSED'

contains

   function build_cfg() result(c)
      type(ed_config_t) :: c
      c%output%enabled    = .true.
      c%output%dir        = '.'
      c%output%prefix     = 'test_ro'
      c%output%cohort_max = 16_ik
      c%output%patch_max  = 4_ik
      c%output%sync_every = SYNC_FLUSH
      c%output%freq_on    = [.false., .true., .false., .true.]     ! daily + annual
      c%output%file_chunk = [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN]
   end function build_cfg

   !----- Daily file: 3 records, mean_agb_py = [10,20,30], n_cohort=1, mean_agb_co slab = [10,20,30], !
   !      each behind the file's one polygon. ---------------------------------------------------!
   subroutine check_daily(path)
      character(len=*), intent(in) :: path
      integer(c_int)    :: ncid, vid
      integer(c_size_t) :: nt
      real(c_double)    :: agbs(3), agbc(3), gpp(3), soilt(3)
      integer(c_int)    :: ncoh(3), days(3)
      call nc_check(nc_open_f(trim(path), NC_NOWRITE, ncid), 'open daily')
      call nc_check(nc_inq_dimlen_f(ncid, 'time', nt), 'daily time len')
      call check(int(nt, ik) == 3_ik, 'daily has 3 records')
      call nc_check(nc_inq_varid_f(ncid, 'mean_agb_py', vid), 'daily mean_agb_py id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [3_c_size_t, 1_c_size_t], agbs), &
                    'get mean_agb_py')
      call check_close(real(agbs(1), wp), 10.0_wp, 1.0e-9_wp, 'daily mean_agb_py day1')
      call check_close(real(agbs(2), wp), 20.0_wp, 1.0e-9_wp, 'daily mean_agb_py day2')
      call check_close(real(agbs(3), wp), 30.0_wp, 1.0e-9_wp, 'daily mean_agb_py day3')
      !----- Each record is stamped with the day its step started: the step from 1 to 2            !
      !      January is the 1 January record, not the 2 January one. ------------------------------!
      call nc_check(nc_inq_varid_f(ncid, 'day', vid), 'daily day id')
      call nc_check(nc_get_vara_int(ncid, vid, [0_c_size_t], [3_c_size_t], days), 'get day')
      call check(all(days == [1_c_int, 2_c_int, 3_c_int]), 'daily records are stamped 1, 2, 3 January')
      call nc_check(nc_inq_varid_f(ncid, 'n_cohort', vid), 'daily n_cohort id')
      call nc_check(nc_get_vara_int(ncid, vid, [0_c_size_t, 0_c_size_t], [3_c_size_t, 1_c_size_t], ncoh), &
                    'get n_cohort')
      call check(all(ncoh == 1_c_int), 'daily n_cohort == 1')
      call nc_check(nc_inq_varid_f(ncid, 'mean_agb_co', vid), 'daily mean_agb_co id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t, 0_c_size_t],               &
                    [3_c_size_t, 1_c_size_t, 1_c_size_t], agbc), 'get mean_agb_co')
      call check_close(real(agbc(1), wp), 10.0_wp, 1.0e-9_wp, 'daily mean_agb_co day1 slot1')
      call check_close(real(agbc(3), wp), 30.0_wp, 1.0e-9_wp, 'daily mean_agb_co day3 slot1')
      !----- Carbon: mean_gpp_py, 2 kgC/m2 a day as its mean rate [umol/m2/s]. -----!
      call nc_check(nc_inq_varid_f(ncid, 'mean_gpp_py', vid), 'daily mean_gpp_py id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [3_c_size_t, 1_c_size_t], gpp), &
                    'get mean_gpp_py')
      call check_close(real(gpp(2), wp), 2.0_wp * kgc_2_umol / day_sec, 1.0e-9_wp,                   &
                       'daily mean_gpp_py = 2 kgC/m2/day as umol/m2/s')
      !----- Soil (DIM_SOIL): mean_soil_temp_layer_si(time, polygon, site, soil) AGG_TMEAN = 290 K per layer. -!
      call nc_check(nc_inq_varid_f(ncid, 'mean_soil_temp_layer_si', vid), 'daily mean_soil_temp_layer_si id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t, 0_c_size_t, 0_c_size_t],   &
                    [3_c_size_t, 1_c_size_t, 1_c_size_t, 1_c_size_t], soilt), 'get soil_temp')
      call check_close(real(soilt(1), wp), 290.0_wp, 1.0e-9_wp, 'daily mean_soil_temp_layer_si layer1 = 290')

      !----- ...and the INACTIVE tail is fill, not data. The source array is 290 K in every       !
      !      slot including the padding, so if the writer emitted the full ceiling this would read  !
      !      290 and pass for the wrong reason. It has to read _FillValue, as the cohort and patch  !
      !      axes do with their own unused capacity. -----------------------------------------------!
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t, 0_c_size_t, int(N_ACTIVE, c_size_t)], &
                    [1_c_size_t, 1_c_size_t, 1_c_size_t, 1_c_size_t], soilt), 'get soil_temp padding')
      call check_close(real(soilt(1), wp), MISSING_VALUE, 1.0e-9_wp,                              &
                       'mean_soil_temp_layer_si padding slot is _FillValue, not data')

      !----- The coordinate must describe the same number of layers. A soil_z running on past the  !
      !      active column with 0.0 m node depths reads as extra layers at the surface.  -----------!
      call nc_check(nc_inq_varid_f(ncid, 'soil_z', vid), 'daily soil_z id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t, int(N_ACTIVE, c_size_t)], &
                    [1_c_size_t, 1_c_size_t, 1_c_size_t], soilt), 'get soil_z padding')
      call check_close(real(soilt(1), wp), MISSING_VALUE, 1.0e-9_wp,                              &
                       'soil_z padding slot is _FillValue, not 0 m')
      call nc_check(nc_close(ncid), 'close daily')
   end subroutine check_daily

   !----- Annual file: 1 record, site-level ONLY (no cohort dim), mean_agb_py = mean(10,20,30)=20. --!
   subroutine check_annual(path)
      character(len=*), intent(in) :: path
      integer(c_int)    :: ncid, vid, st
      integer(c_size_t) :: nt
      real(c_double)    :: agbs(1)
      call nc_check(nc_open_f(trim(path), NC_NOWRITE, ncid), 'open annual')
      call nc_check(nc_inq_dimlen_f(ncid, 'time', nt), 'annual time len')
      call check(int(nt, ik) == 1_ik, 'annual has 1 record')
      call nc_check(nc_inq_varid_f(ncid, 'mean_agb_py', vid), 'annual mean_agb_py id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [1_c_size_t, 1_c_size_t], agbs), &
                    'get annual mean_agb_py')
      call check_close(real(agbs(1), wp), 20.0_wp, 1.0e-9_wp, 'annual mean_agb_py = mean(10,20,30)')
      !----- annual is site-level only: mean_agb_co must NOT exist. -----!
      st = nc_inq_varid_f(ncid, 'mean_agb_co', vid)
      call check(st /= NC_NOERR, 'annual has NO mean_agb_co (site-level only)')
      call nc_check(nc_close(ncid), 'close annual')
   end subroutine check_annual

end program test_output_roundtrip
