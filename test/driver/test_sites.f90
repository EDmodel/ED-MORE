! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_sites -- a polygon of several sites (docs/configuration.md, Several sites).
!
!   1. THE SITES FROM [[sites]] (ed_config_io: load_sites; ed_region: region_open): two
!    sites, 0.25 of the polygon on its own 0.6 m soil and terrain, 0.75 on [soil]. The
!    config holds them, the polygon is built with them, and they run three days.
! 2. THE SITE AXIS IN THE FILES (ed_output_setup: attach_output): the daily file holds
!    the two sites, each one's share, and the shallow site's soil coordinate filled below
!    its bottom layer.
! 3. A CENSUS BY census_site_id (ed_init_census: init_from_census): each site is built
!    from its own rows.
! 4. THE STANDS EVOLVE AS ONE-SITE POLYGONS (ed_stepper: advance_one_step). Two
!    identical sites of areas 0.3 and 0.7 each evolve as a polygon of one site does, over
!    31 days of fast and slow steps on a short calendar -- a month boundary every 15 days
!    and a year boundary at day 30, which day 31's step performs, with a disturbance and
!    patch fusion: every site's stand is the one-site stand, and the polygon totals are
!    the one site's to round-off. Two sites on different soils each evolve as a one-site
!    polygon on its own soil, and their soil columns differ.
! 5. THREADS (ed_stepper: advance_one_step's patch loop): the same polygon on one thread
!    and on four.
! 6. RESTART (ed_init_restart): the run of 1 restarts from its state file with the
!    config's sites (region_open), and the two-soil stand of 4 comes back from a state
!    file (read_restart_file) with each patch's site and each site's patches, plant carbon
!    and soil water.
! 7. WHAT edmore REFUSES (ed_config_check: check_sites, check_key_use): shares that
!    do not sum to 1, a table without area_fraction, a key no table takes, a census run
!    without census_site_id, [[sites]] in a region, and a site's soil keys beside the soil
!    maps.
!
!   The configs of 1-3, 6 and 7 are derived from example 04's ed_config_eval.toml (the
! test runs from its directory) by derive_config: a case's own keys, then a run block that
! cuts the run to a few days without the fast loop or forcing, then the template, each key
! written once, at its first value. Everything the test writes goes under the work
! directory, the first argument; edmore is the second.
!---~---
program test_sites
   use iso_c_binding,           only : c_int, c_size_t, c_double
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : rho_h2o
   use ed_config,               only : ed_config_t
   use ed_config_io,            only : load_ed_config
   use ed_polygon_state_vars,   only : polygon_t
   use ed_init_bare_ground,     only : init_bare_sites, add_cohort, finalize_init
   use ed_init,                 only : init_fast_state
   use ed_init_restart,         only : restart_file_t, restart_file_open, restart_file_write, restart_file_close, &
                                       read_restart_file
   use ed_stepper,              only : advance_one_step, run_ledgers_t
   use ed_output_reduce,        only : total_agb, total_nplant, total_area
   use ed_region,               only : ed_region_t, region_open, region_step, region_done, region_finalize, &
                                       region_free, RUN_OK
   use ed_time,                 only : ed_time_t
   use ed_netcdf_c
   use ed_test_support,         only : banner, build_test_config, check, check_close, test_report
   use ed_test_run,             only : derive_config, main_refuses
   implicit none

   character(len=*), parameter :: TEMPLATE = 'ed_config_eval.toml'
   integer(ik), parameter :: N_DAY = 31_ik            ! two 15-day months, and the step that turns the year
   real(wp),    parameter :: RTOL  = 1.0e-12_wp       ! relative: check_close scales it by the expected value
   character(len=512) :: work, edmore
   character(len=:), allocatable :: cfg_path, nl
   type(ed_config_t)   :: cfg, cfg_shallow, sites_cfg
   type(ed_region_t)   :: reg
   type(polygon_t)     :: one, two, two_mt, one_shallow, mixed, back
   type(ed_time_t)     :: t_read
   logical     :: ok, found
   integer(ik) :: status, i, is
   real(wp)    :: biggest

   call get_command_argument(1, work)
   call get_command_argument(2, edmore)
   if (len_trim(work) == 0 .or. len_trim(edmore) == 0)                                           &
      error stop 'test_sites: pass the work directory and the edmore executable'
   call execute_command_line('rm -rf "'//trim(work)//'" && mkdir -p "'//trim(work)//'"')
   nl = new_line('a')
   cfg_path = trim(work)//'/sites.toml'

   !--- 1. Two sites from [[sites]]: the config, the polygon, a three-day run.
   call banner('the sites from [[sites]]')
   call derive(cfg_path, run_block('2024-07-01', '2024-07-04', 0_ik, '')//                      &
                         '[output]'//nl//'state_write = true'//nl//two_sites())
   call load_ed_config(cfg_path, sites_cfg)
   call check(allocated(sites_cfg%sites), 'the config holds the [[sites]]')
   if (allocated(sites_cfg%sites)) then
      call check(size(sites_cfg%sites) == 2, 'two [[sites]] tables are two sites')
      call check(sites_cfg%sites(1)%own_ground .and. .not. sites_cfg%sites(2)%own_ground,         &
                 'a table with [soil] keys has its own ground, one without takes the run''s')
      call check_close(sites_cfg%sites(1)%soil%depth, 0.6_wp, 1.0e-12_wp, 'the first site is 0.6 m deep')
      call check_close(sites_cfg%sites(1)%soil%theta_sat, sites_cfg%soil%theta_sat, 0.0_wp,        &
                       'a key the table leaves out takes [soil]''s value')
      call check(sites_cfg%sites(1)%has_elevation .and. .not. sites_cfg%sites(2)%has_elevation,   &
                 'the first site has its own elevation, the second the polygon''s')
   end if

   call region_open(cfg_path, reg, ok, verbose=.false.)
   call check(ok, 'the run of two sites opens')
   associate (poly => reg%poly(1))
      call check(poly%site%n == 2_ik, 'the polygon has the two sites')
      call check_close(poly%site%area(1), 0.25_wp, 1.0e-15_wp, 'the first site is a quarter of it')
      call check_close(poly%site%area(2), 0.75_wp, 1.0e-15_wp, 'the second the rest')
      call check(poly%site%surface(1)%soil_profile%n_active < poly%site%surface(2)%soil_profile%n_active,          &
                 'the shallow site has fewer soil layers')
      call check_close(poly%site%surface(1)%elevation_m, 400.0_wp, 0.0_wp, 'the first site''s elevation')
      call check_close(poly%site%surface(2)%elevation_m, poly%elevation_m, 0.0_wp,                &
                       'the second site sits at the polygon''s elevation')
      call check(poly%site%surface(1)%slope_deg == 12.0_wp .and. poly%site%surface(1)%aspect_deg == 135.0_wp, &
                 'the first site''s slope and aspect')
      call check(all(poly%site%patch_count(1:2) > 0_ik), 'each site has bare-ground patches')
      call check(all(poly%patch%owner_site(1:poly%patch%n - 1_ik) <= poly%patch%owner_site(2:poly%patch%n)), &
                 'the patches are sorted by site')
   end associate
   do while (.not. region_done(reg))
      call region_step(reg, status)
      if (status /= RUN_OK) exit
   end do
   call check(status == RUN_OK, 'the two sites run three days')
   call region_finalize(reg, status)
   call check(status == RUN_OK, 'the run finishes with its area conserved')
   call region_free(reg)

   !--- 2. The site axis in the run's daily file (a cohort tier's files are monthly).
   call banner('the site axis in the files')
   call check_output(trim(work)//'/sites-D-202407.nc')

   !--- 3. A census read by census_site_id.
   call banner('a census by census_site_id')
   call write_census(trim(work)//'/sites_census.csv')
   call derive(cfg_path, run_block('2024-07-01', '2024-07-02', 1_ik, '')//                      &
                         '[init]'//nl//'census_file = "'//trim(work)//'/sites_census.csv"'//nl//     &
                         '[[sites]]'//nl//'area_fraction = 0.4'//nl//'census_site_id = 9'//nl//       &
                         '[[sites]]'//nl//'area_fraction = 0.6'//nl//'census_site_id = 7'//nl)
   call region_open(cfg_path, reg, ok, verbose=.false.)
   call check(ok, 'a census run of two sites opens')
   if (ok) then
      associate (poly => reg%poly(1))
         call check(poly%site%n == 2_ik .and. all(poly%site%patch_count(1:2) > 0_ik),            &
                    'each site has its census patches')
         biggest = 0.0_wp
         do i = 1_ik, poly%cohort%n
            if (poly%patch%owner_site(poly%cohort%owner_patch(i)) == 1_ik) biggest = max(biggest, poly%cohort%dbh(i))
         end do
         call check(biggest > 40.0_wp, 'the first site holds census site 9''s big tree')
      end associate
   end if
   call region_free(reg)

   !--- 4. The stands evolve as one-site polygons.
   call banner('the stands evolve as one-site polygons')
   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   cfg_shallow = cfg
   cfg_shallow%soil%depth = 0.6_wp
   !---~---
   !   Two identical sites, 0.3 and 0.7 of the polygon, against one; the same on four threads
   ! for 5.
   !---~---
   call stand(one, cfg, [cfg], [1.0_wp])
   call stand(two, cfg, [cfg, cfg], [0.3_wp, 0.7_wp])
   call stand(two_mt, cfg, [cfg, cfg], [0.3_wp, 0.7_wp])
   call run(one, 1_ik)
   call run(two, 1_ik)
   call run(two_mt, 4_ik)
   call check(one%patch%n > 2_ik, 'the year boundary added a disturbance gap to the one site')
   do is = 1_ik, 2_ik
      call check(two%site%patch_count(is) == one%patch%n,                                            &
                 'each identical site has the one site''s patches after the disturbance and fusion')
      call check_close(site_carbon(two, is), site_carbon(one, 1_ik), RTOL,                             &
                       'each identical site holds the one site''s plant carbon')
      call check_close(site_water(two, is), site_water(one, 1_ik), RTOL,                              &
                       'each identical site holds the one site''s soil water')
   end do
   call check(all(two%patch%owner_site(1:two%patch%n) == [(1_ik, is = 1_ik, two%site%patch_count(1)),  &
                                                            (2_ik, is = 1_ik, two%site%patch_count(2))]), &
              'the patch block stays sorted by site through the disturbance and fusion')
   call check_close(total_agb(two), total_agb(one), RTOL,                                           &
                    'the polygon''s AGB is the one site''s to round-off')
   call check_close(total_nplant(two), total_nplant(one), RTOL,                                     &
                    'the polygon''s stem density is the one site''s to round-off')
   call check_close(total_area(two), 1.0_wp, 1.0e-12_wp, 'the patches still cover the polygon')
   do is = 1_ik, 2_ik
      associate (o => two%site%patch_offset(is), c => two%site%patch_count(is))
         call check_close(sum(two%patch%area(o:o + c - 1_ik)), 1.0_wp, 1.0e-12_wp,                  &
                          'each site''s patch areas still sum to 1')
      end associate
   end do
   !--- Two soils: each site evolves on its own.
   call stand(one_shallow, cfg_shallow, [cfg_shallow], [1.0_wp])
   call stand(mixed, cfg, [cfg, cfg_shallow], [0.5_wp, 0.5_wp])
   call run(one_shallow, 1_ik)
   call run(mixed, 1_ik)
   call check(mixed%site%surface(1)%soil_profile%n_active /= mixed%site%surface(2)%soil_profile%n_active,           &
              'the two sites have different soil columns')
   call check(abs(site_water(mixed, 1_ik) - site_water(mixed, 2_ik)) > 1.0_wp,                      &
              'the two soils hold different water')
   call check_close(site_carbon(mixed, 1_ik), site_carbon(one, 1_ik), RTOL,                         &
                    'the deep site grows the one-site stand on the deep soil')
   call check_close(site_carbon(mixed, 2_ik), site_carbon(one_shallow, 1_ik), RTOL,                 &
                    'the shallow site grows the one-site stand on its soil')
   call check_close(site_water(mixed, 2_ik), site_water(one_shallow, 1_ik), RTOL,                   &
                    'the shallow site keeps the one-site water on its soil')

   !--- 5. One thread or four: the same polygon.
   call banner('one thread or four')
   call check(all(two_mt%patch%area(1:two%patch%n) == two%patch%area(1:two%patch%n)) .and.        &
              site_carbon(two_mt, 1_ik) == site_carbon(two, 1_ik) .and.                             &
              site_carbon(two_mt, 2_ik) == site_carbon(two, 2_ik) .and.                             &
              site_water(two_mt, 2_ik) == site_water(two, 2_ik),                                    &
              'four patch threads give the polygon one thread gives, bit for bit')

   !--- 6. Restart.
   call banner('restart')
   !--- The run of 1 from its state file, with the config's [[sites]].
   call derive(cfg_path, run_block('2024-07-04', '2024-07-05', 2_ik,                            &
                                   trim(work)//'/sites-S-20240704000000.nc')//two_sites())
   call region_open(cfg_path, reg, ok, verbose=.false.)
   call check(ok, 'the run restarts from its state file')
   if (ok) then
      call check(reg%poly(1)%site%n == 2_ik .and. reg%poly(1)%site%area(1) == 0.25_wp,            &
                 'the restart reads the two sites and their shares')
      call check(reg%poly(1)%site%surface(1)%soil_profile%n_active < reg%poly(1)%site%surface(2)%soil_profile%n_active, &
                 'each site keeps its soil column')
   end if
   call region_free(reg)
   !--- The two-soil stand of 4, grown with cohorts, through a state file.
   call write_state(mixed, ed_time_t(year=2000, month=1, day=31))
   call read_restart_file(back, cfg, trim(work)//'/stands-S-20000131000000.nc', t_read, found,      &
                          polygon=1_ik, site_cfg=[cfg, cfg_shallow])
   call check(found .and. back%site%n == 2_ik, 'the state file holds the two sites')
   call check(all(back%patch%owner_site(1:back%patch%n) == mixed%patch%owner_site(1:mixed%patch%n)), &
              'each patch keeps its site')
   do is = 1_ik, 2_ik
      call check(back%site%patch_count(is) == mixed%site%patch_count(is), 'each site keeps its patches')
      call check_close(site_carbon(back, is), site_carbon(mixed, is), RTOL,                         &
                       'each site keeps its plant carbon')
      call check(site_water(back, is) == site_water(mixed, is), 'each site keeps its soil water')
   end do

   !--- 7. What edmore refuses.
   call banner('what edmore refuses')
   call refused('shares that do not sum to 1',                                                  &
                '[[sites]]'//nl//'area_fraction = 0.5'//nl//'[[sites]]'//nl//'area_fraction = 0.4'//nl, &
                'must sum to 1')
   call refused('a table without area_fraction',                                                &
                '[[sites]]'//nl//'slope = 5.0'//nl, 'sites.area_fraction')
   call refused('a key no [[sites]] table takes',                                               &
                '[[sites]]'//nl//'area_fraction = 1.0'//nl//'w_pond_max = 1.0'//nl,              &
                'sites.1.w_pond_max')
   call refused('a census run without census_site_id',                                          &
                '[init]'//nl//'init_mode = 1'//nl//'[[sites]]'//nl//'area_fraction = 1.0'//nl,    &
                'census_site_id')
   !---~---
   !   A region takes one site per cell, and with the soil maps a site takes its cell's soil.
   ! The maps are never opened: the config is refused first.
   !---~---
   call refused('[[sites]] in a region', '[run]'//nl//'mode = "region"'//nl//'[region]'//nl//       &
                'box_nwse = [10.0, -80.0, 9.0, -79.0]'//nl//'[[sites]]'//nl//'area_fraction = 1.0'//nl, &
                'a region refuses [[sites]]')
   call refused('a [[sites]] soil key beside the soil maps', '[maps]'//nl//                     &
                'soil_depth_map = "'//trim(work)//'/soil_depth.nc"'//nl//                         &
                'soil_texture_map = "'//trim(work)//'/soil_texture.nc"'//nl//                     &
                'soil_chemistry_map = "'//trim(work)//'/soil_chemistry.nc"'//nl//                 &
                '[[sites]]'//nl//'area_fraction = 1.0'//nl//'depth = 1.0'//nl,                    &
                'every site takes its cell''s soil')

   call test_report('test_sites')

contains

   !---~---
   !   The run block: a short run from `mode` (0 bare ground, 1 census, 2 restart), no fast
   ! loop or forcing, daily water and structure output into the work directory.
   !---~---
   function run_block(start, end, mode, restart) result(s)
      character(len=*), intent(in) :: start, end, restart
      integer(ik),      intent(in) :: mode
      character(len=:), allocatable :: s
      character(len=8) :: m
      write(m, '(i0)') mode
      s = '[run]'//nl//'start_time = "'//start//'"'//nl//'end_time = "'//end//'"'//nl//            &
          'dt_slow = "1d"'//nl//'[fast]'//nl//'fast_biophysics_on = false'//nl//                   &
          'gpp_ref = 0.3'//nl//                                                                    &
          '[forcing]'//nl//'forcing_on = false'//nl//                                              &
          '[init]'//nl//'init_mode = '//trim(m)//nl//'restart_file = "'//restart//'"'//nl//         &
          '[output]'//nl//'state_file_dir = "'//trim(work)//'"'//nl//'state_file_prefix = "sites"'//nl// &
          'dir = "'//trim(work)//'"'//nl//'prefix = "sites"'//nl//                                  &
          'io_config = ""'//nl//                                                                    &
          '[output.fast]'//nl//'enabled = false'//nl//'[output.daily]'//nl//'enabled = true'//nl//  &
          'file_chunk = "run"'//nl
   end function run_block

   !---~---
   !   Two sites: a quarter of the polygon on 0.6 m of soil, at 400 m on a 12-degree slope
   ! facing south-east; the rest on [soil].
   !---~---
   function two_sites() result(s)
      character(len=:), allocatable :: s
      s = '[[sites]]'//nl//'area_fraction = 0.25'//nl//'depth = 0.6'//nl//'elevation = 400.0'//nl//  &
          'slope = 12.0'//nl//'aspect = 135.0'//nl//'[[sites]]'//nl//'area_fraction = 0.75'//nl
   end function two_sites

   !--- A config from these blocks and the template, each key once (ed_test_run).
   subroutine derive(path, overrides)
      character(len=*), intent(in) :: path, overrides
      call derive_config(path, overrides, TEMPLATE)
   end subroutine derive

   !---~---
   !   edmore must stop on the config, saying `message`. The case's keys come first, so they
   ! win over the run block's.
   !---~---
   subroutine refused(what, overrides, message)
      character(len=*), intent(in) :: what, overrides, message
      call derive(trim(work)//'/refused.toml', overrides//run_block('2024-07-01', '2024-07-02', 0_ik, ''))
      call check(main_refuses(trim(edmore), trim(work)//'/refused.toml', trim(work)//'/refused.log', message), &
                 'edmore refuses '//what)
   end subroutine refused

   !--- Census sites 7 and 9, two patches each; site 9 holds the one tree over 40 cm.
   subroutine write_census(path)
      character(len=*), intent(in) :: path
      integer :: u
      open(newunit=u, file=path, status='replace', action='write')
      write(u, '(a)') 'site_id,patch_id,dbh,pft,nplant'
      write(u, '(a)') '7,1,12.0,1,0.02'
      write(u, '(a)') '7,2,8.0,1,0.05'
      write(u, '(a)') '9,1,45.0,1,0.005'
      write(u, '(a)') '9,2,5.0,1,0.10'
      close(u)
   end subroutine write_census

   !---~---
   !   The daily file: a site axis of two, each site's share, and the shallow site's soil
   ! coordinate filled below its bottom layer.
   !---~---
   subroutine check_output(path)
      character(len=*), intent(in) :: path
      integer(c_int)    :: ncid, vid
      integer(c_size_t) :: n
      real(c_double)    :: area(2), z(1)
      call check(nc_open_f(path, NC_NOWRITE, ncid) == NC_NOERR, 'the daily file is written')
      call check(nc_inq_dimlen_f(ncid, 'site', n) == NC_NOERR .and. n == 2_c_size_t,              &
                 'its site axis holds the two sites')
      call nc_check(nc_inq_varid_f(ncid, 'area_si', vid), 'area_si')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [1_c_size_t, 2_c_size_t], area), &
                    'get site_area')
      call check(area(1) == 0.25_c_double .and. area(2) == 0.75_c_double, 'it holds each site''s share')
      call nc_check(nc_inq_varid_f(ncid, 'soil_z', vid), 'soil_z')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t, 6_c_size_t],            &
                    [1_c_size_t, 1_c_size_t, 1_c_size_t], z), 'get soil_z')
      call check(z(1) > 1.0e30_c_double, 'the shallow site''s soil coordinate stops at its bottom layer')
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_output

   !---~---
   !   A two-patch stand in each site: a canopy tree and an understorey in one patch, a
   ! sapling in the other.
   !---~---
   subroutine stand(poly, run_cfg, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(ed_config_t),   intent(in)  :: run_cfg, site_cfg(:)
      real(wp),            intent(in)  :: site_area(:)
      integer(ik) :: k, ns
      ns = int(size(site_cfg), ik)
      call init_bare_sites(poly, run_cfg, [(2_ik, k = 1_ik, ns)], site_cfg, site_area)
      do k = 0_ik, ns - 1_ik
         call add_cohort(poly, run_cfg, 2_ik*k + 1_ik, 1_ik, 0.05_wp, 30.0_wp)
         call add_cohort(poly, run_cfg, 2_ik*k + 1_ik, 2_ik, 0.20_wp,  4.0_wp)
         call add_cohort(poly, run_cfg, 2_ik*k + 2_ik, 3_ik, 0.50_wp,  2.0_wp)
      end do
      call finalize_init(poly)
      call init_fast_state(poly, run_cfg)
   end subroutine stand

   !---~---
   !   N_DAY daily steps, a month's restructuring owed every 15 days and a year's every 30,
   ! which the next step performs first, as a run does.
   !---~---
   subroutine run(poly, patch_threads)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: patch_threads
      type(run_ledgers_t) :: ledgers
      integer(ik) :: day
      do day = 1_ik, N_DAY
         call advance_one_step(poly, cfg, ledgers, patch_threads=patch_threads)
         poly%restructure_pending  = mod(day, 15_ik) == 0_ik
         poly%restructure_new_year = mod(day, 30_ik) == 0_ik
      end do
   end subroutine run

   !--- The polygon as the one polygon of a state file in the work directory (id 0, at 0 N 0 E).
   subroutine write_state(poly, t)
      type(polygon_t),   intent(in) :: poly
      type(ed_time_t),   intent(in) :: t
      type(restart_file_t) :: sf
      call restart_file_open(sf, trim(work), 'stands', t, poly%n_pft,                               &
                           int(size(poly%cohort%growth_hist, 1), ik), poly%cohort%n, poly%patch%n,    &
                           poly%site%n, polygon_id=[0_ik], polygon_lat=[0.0_wp], polygon_lon=[0.0_wp])
      call restart_file_write(sf, poly, t, 1_ik)
      call restart_file_close(sf)
   end subroutine write_state

   !--- Plant carbon [kgC/m2 of the site]: its patches' cohorts, weighted by patch area.
   real(wp) function site_carbon(poly, is) result(c)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      integer(ik) :: ip, i
      c = 0.0_wp
      do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
         do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik
            c = c + poly%patch%area(ip) * poly%cohort%nplant(i) * (poly%cohort%leaf_carbon(i)           &
                  + poly%cohort%fineroot_carbon(i) + poly%cohort%wood_carbon(i)                       &
                  + poly%cohort%nonstructural_carbon(i))
         end do
      end do
   end function site_carbon

   !--- Soil water [kg/m2 of the site] over its patches' active layers.
   real(wp) function site_water(poly, is) result(w)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      integer(ik) :: ip, k
      w = 0.0_wp
      associate (soil => poly%site%surface(is)%soil_profile)
         do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
            do k = 1_ik, soil%n_active
               w = w + poly%patch%area(ip) * poly%patch%soil_water(ip)%theta(k) * soil%dz(k) * rho_h2o
            end do
         end do
      end associate
   end function site_water

end program test_sites
