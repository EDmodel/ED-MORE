! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_init -- how a run's polygon starts (src/io/init; docs/configuration.md, How a run     !
! starts).                                                                                  !
!   1. THE CENSUS READER (meds_init_census): a site is built from a cohort CSV -- the site     !
!      filter, the patch map, the allometry from dbh, comment and blank lines, a header-less    !
!      file read positionally, columns by name in any order, patch_area normalized to the site, !
!      and restructure_census_stand conserving the stand's stems and biomass and honouring its  !
!      switches. A missing file or a site with no row is reported, not a crash.                  !
!   2. THE CENSUS REFUSALS: every reason the reader gives for a file it cannot use, each from a !
!      file written to give exactly that reason. The driver stops the run with the reason.       !
!   3. THE INITIAL SOIL STATE (meds_init_bare_ground, meds_init): [init].soil_temp and           !
!      soil_theta reach every layer of every patch through build_soil_profile, with the         !
!      defaults pinned (288 K, theta = 0.30); a start below freezing holds the liquid its        !
!      freezing curve allows; a texture per layer starts each layer at its own field capacity;  !
!      and bare ground is empty unless a PFT's nbg_init_density seeds it at min_cohort_height.  !
!   4. THE INITIAL SNOWPACK (meds_init: init_fast_state; docs/science/snow_biophysics.md): no   !
!      [fast].snow_init_swe seeds no pack; one seeds every patch with one frozen layer of that   !
!      water equivalent at [fast].snow_init_temp, as deep as its water over the snow density,   !
!      at a configured density too.                                                            !
!   5. THE CANOPY AIR'S DEPTH (meds_init: init_polygon, which ends every start with the depth  !
!      the stand sets): a bare-ground start, which has no stand, opens each patch's canopy air !
!      at [aerodynamics].min_canopy_depth.                                                     !
!==========================================================================================!
program test_init
   use meds_kinds,              only : wp, ik
   use meds_config,             only : meds_config_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_init_census,        only : init_from_census, restructure_census_stand
   use meds_output_reduce,      only : total_nplant
   use meds_soil_derive,        only : soil_layers, build_soil_profile
   use meds_soil_params,        only : soil_profile_t
   use meds_init_bare_ground,   only : init_bare_ground, seed_bare_ground, finalize_init
   use meds_constants,          only : n_soil_layer_max, rho_h2o
   use meds_soil_heat,          only : soil_layer_temp
   use meds_soil_hydr_lib,      only : soil_liquid_capacity
   use meds_init,               only : init_fast_state, init_polygon
   use meds_met_types,          only : met_forcing_t
   use meds_time,               only : meds_time_t
   use meds_test_support,       only : banner, build_test_config, check, check_close
   implicit none

   call test_census()
   call test_census_refusals()
   call test_soil_state()
   call test_initial_snowpack()
   call test_canopy_air_floor()
   write(*,'(a)') 'test_init: ALL PASSED'

contains

   !----- 1. The census reader. --------------------------------------------------------------!
   subroutine test_census()
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      integer(ik)         :: u, i, ip
      logical             :: found
      real(wp)            :: n0, agb0
      character(len=:), allocatable :: why


      call banner('the census reader')
      cfg = build_test_config()

      !----- Write a small census file (with a comment + blank line to exercise skipping). ----!
      open(newunit=u, file='test_census_tmp.csv', status='replace', action='write')
      write(u,'(a)') '# pseudo cohort census (two patches in site 1, one row for site 2)'
      write(u,'(a)') 'site_id,patch_id,cohort_id,dbh,height,pft,nplant'
      write(u,'(a)') ''
      write(u,'(a)') '1,1,1,30.0,21.0,3,0.05'
      write(u,'(a)') '1,1,2,10.0,12.0,2,0.20'
      write(u,'(a)') '1,2,3, 5.0, 8.0,1,0.40'
      write(u,'(a)') '2,1,1,40.0,27.0,3,0.02'
      close(u)

      !=== Default: first site_id (1) selected; site 2 row excluded. ==========================!
      call init_from_census(poly, cfg, 'test_census_tmp.csv', found)
      call check(found, 'census file should be usable')
      call check(poly%patch%n  == 2_ik, 'site 1 should yield 2 patches')
      call check(poly%cohort%n == 3_ik, 'site 1 should yield 3 cohorts (site 2 row excluded)')
      call check(sum(poly%patch%cohort_count(1:poly%patch%n)) == poly%cohort%n, 'CSR count mismatch')
      !----- dbh drove the allometry: every cohort has a positive derived height/agb/leaf area. !
      call check(all(poly%cohort%height(1:poly%cohort%n)    > 0.0_wp), 'derived height not positive')
      call check(all(poly%cohort%agb(1:poly%cohort%n)       > 0.0_wp), 'derived agb not positive')
      call check(all(poly%cohort%leaf_area(1:poly%cohort%n) > 0.0_wp), 'derived leaf area not positive')

      !=== Explicit site selection: site 2 -> one patch, one (climax, pft 3) cohort. ==========!
      call init_from_census(poly, cfg, 'test_census_tmp.csv', found, use_site_id = 2_ik)
      call check(found, 'site 2 should be usable')
      call check(poly%patch%n  == 1_ik, 'site 2 should yield 1 patch')
      call check(poly%cohort%n == 1_ik, 'site 2 should yield 1 cohort')
      call check(poly%cohort%pft(1) == 3_ik, 'site 2 cohort should be pft 3')

      !=== Header-LESS CSV: the first data line must be KEPT, not silently dropped as a header. ==!
      open(newunit=u, file='test_census_nohdr.csv', status='replace', action='write')
      write(u,'(a)') '1,1,1,30.0,21.0,3,0.05'      ! no header row -> the first line is DATA
      write(u,'(a)') '1,1,2,10.0,12.0,2,0.20'
      write(u,'(a)') '1,2,3, 5.0, 8.0,1,0.40'
      close(u)
      call init_from_census(poly, cfg, 'test_census_nohdr.csv', found)
      call check(found, 'header-less census should be usable')
      call check(poly%cohort%n == 3_ik, 'header-less census must keep all 3 cohorts (first not dropped)')
      open(newunit=u, file='test_census_nohdr.csv', status='old', action='read')
      close(u, status='delete')

      !=== Columns by NAME, in any order, with patch_area and without site_id. =================!
      !    Patch 7 has three times the area of patch 9, so the site fractions are 0.75 and 0.25.  !
      open(newunit=u, file='test_census_named.csv', status='replace', action='write')
      write(u,'(a)') '# columns reordered, areas in m2'
      write(u,'(a)') 'nplant,pft,dbh,patch_area,patch_id'
      write(u,'(a)') '0.05,1,30.0,300.0,7'
      write(u,'(a)') '0.20,1,10.0,300.0,7'
      write(u,'(a)') '0.40,1,5.0,100.0,9'
      close(u)
      call init_from_census(poly, cfg, 'test_census_named.csv', found)
      call check(found, 'named-column census should be usable')
      call check(poly%patch%n  == 2_ik, 'named-column census should yield 2 patches')
      call check(poly%cohort%n == 3_ik, 'named-column census should yield 3 cohorts')
      call check_close(sum(poly%patch%area(1:2)), 1.0_wp, 1.0e-12_wp, 'patch areas must sum to the site')
      call check_close(poly%patch%area(1), 0.75_wp, 1.0e-12_wp, 'patch_area not normalized to the site (7)')
      call check_close(poly%patch%area(2), 0.25_wp, 1.0e-12_wp, 'patch_area not normalized to the site (9)')
      do i = 1_ik, poly%cohort%n
         if (abs(poly%cohort%dbh(i) - 5.0_wp) < 1.0e-9_wp) then
            call check(poly%cohort%owner_patch(i) == 2_ik, 'the 5 cm cohort must sit in patch 9')
            call check_close(poly%cohort%nplant(i), 0.40_wp, 1.0e-12_wp, 'nplant not read by name')
         end if
      end do
      open(newunit=u, file='test_census_named.csv', status='old', action='read')
      close(u, status='delete')
      !=== Restructuring a census stand: 20 identical cells, each with two rows of the same size,   !
      !    fuse to one patch and one cohort per size, conserving the site's stems and biomass.     !
      open(newunit=u, file='test_census_cells.csv', status='replace', action='write')
      write(u,'(a)') 'patch_id,patch_area,dbh,pft,nplant'
      do ip = 1_ik, 20_ik
         write(u,'(i0,a)') ip, ',400.0,25.0,1,0.0025'
         write(u,'(i0,a)') ip, ',400.0,25.0,1,0.0025'
         write(u,'(i0,a)') ip, ',400.0,4.0,1,0.0500'
      end do
      close(u)
      call init_from_census(poly, cfg, 'test_census_cells.csv', found)
      call check(found .and. poly%patch%n == 20_ik .and. poly%cohort%n == 60_ik, 'the cell census did not load')
      n0   = total_nplant(poly)
      agb0 = sum([(poly%patch%area(poly%cohort%owner_patch(i)) * poly%cohort%nplant(i) * poly%cohort%agb(i), &
                   i = 1_ik, poly%cohort%n)])
      call restructure_census_stand(poly, cfg)
      call check(poly%patch%n == 1_ik,  'identical census cells must fuse to one patch')
      call check(poly%cohort%n == 2_ik, 'identical sizes must fuse to one cohort each')
      call check_close(total_nplant(poly), n0, 1.0e-12_wp, 'restructuring must conserve the site''s stems')
      call check_close(sum([(poly%patch%area(poly%cohort%owner_patch(i)) * poly%cohort%nplant(i)             &
                             * poly%cohort%agb(i), i = 1_ik, poly%cohort%n)]), agb0, 1.0e-12_wp,         &
                       'restructuring must conserve the site''s biomass')
      !----- With demography off it leaves the stand alone. ----------------------------------!
      call init_from_census(poly, cfg, 'test_census_cells.csv', found)
      cfg%demography%demography_on = .false.
      call restructure_census_stand(poly, cfg)
      call check(poly%patch%n == 20_ik .and. poly%cohort%n == 60_ik, 'restructuring must honour the switches')
      cfg%demography%demography_on = .true.
      open(newunit=u, file='test_census_cells.csv', status='old', action='read')
      close(u, status='delete')
      !=== A missing file, or a site the file has no row for, is reported (found=.false.) with  !
      !    the reason, not a crash: the driver stops the run with it. ========================!
      call init_from_census(poly, cfg, 'no_such_census_file.csv', found, why=why)
      call check(.not. found, 'missing census file should report found=.false.')
      call check(why == 'cannot be opened', 'a missing census file should say it cannot be opened')
      call init_from_census(poly, cfg, 'test_census_tmp.csv', found, use_site_id=7_ik, why=why)
      call check(.not. found, 'a site with no row should report found=.false.')
      call check(why == 'holds no cohort for the site', 'a site with no row should say the file holds none')

      open(newunit=u, file='test_census_tmp.csv', status='old', action='read')
      close(u, status='delete')

   end subroutine test_census

   !----- 2. Every reason the census reader gives for a file it cannot use. Each case is a header !
   !      and up to two rows written to give exactly that reason, and the reason is matched in  !
   !      full: a file the reader accepted, or refused for another reason, fails. -------------!
   subroutine test_census_refusals()
      integer, parameter  :: N_CASE = 12
      character(len=64)   :: header(N_CASE), row1(N_CASE), row2(N_CASE)
      character(len=72)   :: expect(N_CASE)
      character(len=256)  :: wide
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      logical             :: found
      character(len=:), allocatable :: why
      integer             :: i, j, u
      call banner('the census refusals')
      cfg = build_test_config()
      header = 'patch_id,dbh,pft,nplant'
      row1 = ''
      row2 = ''
      header(1) = 'patch_id,basal,dbh,pft,nplant'      ; row1(1) = '1,0.5,10.0,1,0.01'
      expect(1) = 'has an unknown column "basal"'
      header(2) = 'patch_id,patch_area,dbh,pft,nplant' ; row1(2) = '1,400.0,10.0,1,0.01' ; row2(2) = '1,300.0,20.0,1,0.01'
      expect(2) = 'has two different patch_area values for patch_id 1'
      header(3) = 'patch_id,dbh,dbh,pft,nplant'        ; row1(3) = '1,10.0,10.0,1,0.01'
      expect(3) = 'has the column "dbh" twice'
      header(4) = 'dbh,pft,nplant'                     ; row1(4) = '10.0,1,0.01'
      expect(4) = 'has no patch_id column'
      header(5) = 'patch_id,pft,nplant'                ; row1(5) = '1,1,0.01'
      expect(5) = 'has no dbh column'
      header(6) = 'patch_id,dbh,nplant'                ; row1(6) = '1,10.0,0.01'
      expect(6) = 'has no pft column'
      header(7) = 'patch_id,dbh,pft'                   ; row1(7) = '1,10.0,1'
      expect(7) = 'has no nplant column'
      row1(8)   = '1,10.0,1'
      expect(8) = 'has a row with too few columns'
      row1(9)   = '1,ten,1,0.01'
      expect(9) = 'has a malformed row'
      row1(10)  = '1,-10.0,1,0.01'
      expect(10) = 'has a row with a non-positive dbh, height, nplant or patch_area'
      row1(11)  = '1,10.0,99,0.01'                                     ! the fixture has 3 PFTs
      expect(11) = 'has a pft index out of range'
      expect(12) = 'has a header with too many columns'                ! 33 names, one past the 32 it takes
      do i = 1, N_CASE
         open(newunit=u, file='test_init_refused_tmp.csv', status='replace', action='write')
         if (i == N_CASE) then
            write(wide, '(*(a,i0,:,","))') ('c', j, j = 1, 33)
            write(u, '(a)') trim(wide)
         else
            write(u, '(a)') trim(header(i))
            write(u, '(a)') trim(row1(i))
            if (len_trim(row2(i)) > 0) write(u, '(a)') trim(row2(i))
         end if
         close(u)
         call init_from_census(poly, cfg, 'test_init_refused_tmp.csv', found, why=why)
         call check(.not. found, 'the census was accepted, but '//trim(expect(i))//' should refuse it')
         call check(why == trim(expect(i)), 'the census was refused as "'//why//'", not as "'//trim(expect(i))//'"')
      end do
      open(newunit=u, file='test_init_refused_tmp.csv', status='old', action='read')
      close(u, status='delete')
   end subroutine test_census_refusals

   !----- 3. The initial soil state. ---------------------------------------------------------!
   subroutine test_soil_state()
      type(meds_config_t)  :: cfg
      type(polygon_t)      :: poly
      type(soil_profile_t) :: prof
      integer(ik)          :: ip, nsl, k, n
      real(wp)             :: bottoms(n_soil_layer_max), tdiag, fdiag, theta_liq, dtheta
      type(met_forcing_t)  :: ref   !< the reference climate a start's canopy air is seeded at


      call banner('the initial soil state')
      cfg = build_test_config()

      !=== 1. The defaults: 288 K and theta = 0.30. ===========================================!
      call check_close(cfg%init%soil_temp,  288.0_wp, 1.0e-14_wp, 'default soil_temp changed (was 288 K)')
      call check_close(cfg%init%soil_theta, 0.30_wp,  1.0e-14_wp, 'default soil_theta changed (was 0.30)')

      !=== 2. Set values reach the site, then every layer of every patch. ======================!
      cfg%init%soil_temp  = 298.65_wp
      cfg%init%soil_theta = 0.25_wp
      call build_soil_profile(cfg%soil, prof)
      call init_bare_ground(poly, cfg, 3_ik)
      call check_close(poly%site%surface(1)%theta_start(1), 0.25_wp, 1.0e-14_wp, 'soil_theta did not reach the site')
      call finalize_init(poly)
      call init_fast_state(poly, cfg)
      nsl = prof%n_active
      do ip = 1_ik, poly%patch%n
         call check(all(abs(poly%patch%soil_energy(ip)%soil_temp(1:nsl) - 298.65_wp) < 1.0e-12_wp),         &
                    'a soil layer did not start at [init].soil_temp')
         call check(all(abs(poly%patch%soil_water(ip)%theta(1:nsl) - 0.25_wp) < 1.0e-12_wp),              &
                    'a soil layer did not start at [init].soil_theta')
         call check(poly%patch%cas(ip)%can_temp == ref%air_temp .and. poly%patch%cas(ip)%can_shv == ref%air_shv, &
                    'a laid-out patch''s canopy air is not at the reference climate')
      end do

      !=== 2b. A start below freezing: the water starts as the liquid its freezing curve holds and  !
      !        the rest ice, and the stored energy gives back the start temperature. ==================!
      cfg%init%soil_temp = 260.0_wp
      call init_bare_ground(poly, cfg, 3_ik)
      associate (se => poly%patch%soil_energy(2), surf => poly%site%surface(1))
         do k = 1_ik, nsl
            call soil_layer_temp(se%soil_energy(k), surf%theta_start(k), surf%soil_profile%dry_heat_capacity(k),   &
                                 surf%soil_profile%layer_curve(k), tdiag, fdiag)
            call check_close(tdiag, 260.0_wp, 1.0e-10_wp, 'a frozen start does not give back its temperature')
            call soil_liquid_capacity(surf%soil_profile%layer_curve(k), 260.0_wp, theta_liq, dtheta)
            call check_close(se%soil_fliq(k), theta_liq / surf%theta_start(k), 1.0e-10_wp,                        &
                             'a frozen start does not hold its freezing curve''s liquid')
            call check(fdiag == se%soil_fliq(k) .and. fdiag < 0.5_wp, 'a frozen start is not mostly ice')
         end do
      end associate
      cfg%init%soil_temp = 298.65_wp

      !=== 3. A texture per layer: each layer's own curve, each starting at its field capacity. ======!
      call soil_layers(cfg%soil, bottoms, n)
      associate (sc => cfg%soil)
         sc%n_layer_values = n
         do k = 1_ik, n
            sc%layer_theta_sat(k)   = 0.40_wp + 0.01_wp * real(k, wp)
            sc%layer_theta_res(k)   = 0.05_wp
            sc%layer_ksat(k)        = 1.0e-6_wp / real(k, wp)
            sc%layer_curve_par_a(k) = 2.0_wp + 0.1_wp * real(k, wp)
            sc%layer_curve_par_n(k) = 1.3_wp + 0.02_wp * real(k, wp)
            sc%layer_solid_conductivity(k) = 3.0_wp ; sc%layer_dry_conductivity(k) = 0.2_wp
            sc%layer_dry_heat_capacity(k)  = 1.5e6_wp + 1.0e4_wp * real(k, wp)
         end do
      end associate
      call build_soil_profile(cfg%soil, prof)
      associate (soil => prof)
         call check(soil%n_active == n .and. all(soil%theta_sat(1:n) == cfg%soil%layer_theta_sat(1:n))         &
                    .and. all(soil%curve_par_n(1:n) == cfg%soil%layer_curve_par_n(1:n)),                              &
                    'a layer did not take its own curve')
         call check(all(prof%dry_heat_capacity(1:n) ==                                                   &
                        cfg%soil%layer_dry_heat_capacity(1:n)), 'a layer did not take its own heat capacity')
      end associate
      call init_bare_ground(poly, cfg, 3_ik)
      associate (surf => poly%site%surface(1))
         call check(all(surf%theta_start(1:n) == surf%soil_profile%theta_fc(1:n)) .and.                     &
                    surf%soil_profile%theta_fc(1) /= surf%soil_profile%theta_fc(n),                          &
                    'a layer did not start at its own field capacity')
      end associate
      call finalize_init(poly)
      call init_fast_state(poly, cfg)
      call check(all(poly%patch%soil_water(2)%theta(1:n) == prof%theta_fc(1:n)),                                          &
                 'a patch''s soil did not start at field capacity')

      !=== 4. Bare ground is empty by default; nbg_init_density seeds each patch with each PFT. ======!
      call init_bare_ground(poly, cfg, 3_ik)
      call seed_bare_ground(poly, cfg)
      call check(poly%cohort%n == 0_ik, 'bare ground is not empty with nbg_init_density = 0')
      cfg%pft%nbg_init_density = [0.5_wp, 0.0_wp, 0.2_wp]
      call init_bare_ground(poly, cfg, 3_ik)
      call seed_bare_ground(poly, cfg)
      call check(poly%cohort%n == 6_ik, 'not one cohort per seeded PFT per patch')
      call check(all(abs(poly%cohort%height(1:6) - cfg%pft%min_cohort_height) < 1.0e-9_wp),                    &
                 'a seeded cohort is not at min_cohort_height')
      call check(all(poly%cohort%nplant(1:6) == 0.5_wp .or. poly%cohort%nplant(1:6) == 0.2_wp) .and.             &
                 count(poly%cohort%pft(1:6) == 1_ik) == 3 .and. count(poly%cohort%pft(1:6) == 3_ik) == 3,       &
                 'a seeded cohort has the wrong PFT or density')
   end subroutine test_soil_state

   !----- 4. The initial snowpack: none by default; [fast].snow_init_swe = 40 kg/m2 at 265 K lays one   !
   !      frozen layer on every patch, whose energy gives back 265 K, as deep as the pack's water over  !
   !      the snow density (depth = swe / rho_snow, the relation the pack keeps every step). ---------!
   subroutine test_initial_snowpack()
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      integer(ik)         :: ip
      real(wp), parameter :: SWE = 40.0_wp, T_PACK = 265.0_wp   ! [kg/m2], [K]
      call banner('the initial snowpack')
      cfg = build_test_config()
      call init_bare_ground(poly, cfg, 3_ik)
      call finalize_init(poly)
      call init_fast_state(poly, cfg)
      do ip = 1_ik, poly%patch%n
         call check(poly%patch%snow(ip)%nlayer == 0_ik .and. poly%patch%snow(ip)%swe(1) == 0.0_wp,               &
                    'a patch holds snow with no [fast].snow_init_swe')
      end do
      cfg%snow_init_swe  = SWE
      cfg%snow_init_temp = T_PACK
      call init_bare_ground(poly, cfg, 3_ik)
      call finalize_init(poly)
      call init_fast_state(poly, cfg)
      do ip = 1_ik, poly%patch%n
         associate (snow => poly%patch%snow(ip))
            call check(snow%nlayer == 1_ik, 'a seeded patch does not hold one snow layer')
            call check_close(snow%swe(1), SWE, 1.0e-15_wp, 'a seeded pack does not hold [fast].snow_init_swe')
            call check_close(snow%snow_temp(1), T_PACK, 1.0e-10_wp, 'a seeded pack''s energy does not give back its temperature')
            call check(snow%snow_fliq(1) == 0.0_wp, 'a seeded pack below freezing holds liquid')
            call check_close(snow%snow_depth(1), SWE / cfg%snow%rho_snow, 1.0e-12_wp,                               &
                             'a seeded pack is not as deep as its water over the snow density')
         end associate
      end do
      !----- At another density: a seed at the default density's depth would disagree with its water   !
      !      until the pack's first exchange, and the first window's canopy burial and snow cover     !
      !      would read the wrong depth. -----------------------------------------------------------!
      cfg%snow%rho_snow = 350.0_wp
      call init_bare_ground(poly, cfg, 3_ik)
      call finalize_init(poly)
      call init_fast_state(poly, cfg)
      call check_close(poly%patch%snow(1)%snow_depth(1), SWE / cfg%snow%rho_snow, 1.0e-12_wp,                       &
                       'a pack seeded at [snow].rho_snow = 350 is as deep as its water over 350 kg/m3')
   end subroutine test_initial_snowpack

   !----- 5. The canopy air's depth at a bare-ground start: the floor, not the type's default.     !
   !      The floor is set above the freeboard over the tallest cohort (5 m by default), so a     !
   !      depth taken from the freeboard on the empty stand would not pass. ----------------------!
   subroutine test_canopy_air_floor()
      real(wp), parameter :: FLOOR = 7.0_wp   ! [m]
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      type(meds_time_t)   :: now
      logical             :: ok
      call banner('the canopy air''s depth at a bare-ground start')
      cfg = build_test_config()
      cfg%aero%min_canopy_depth = FLOOR
      now = cfg%start_time
      call init_polygon(poly, cfg, 'bare ground', .false., 1_ik, now, ok, 0_ik, 0.0_wp, 0.0_wp, 0.0_wp)
      call check(ok .and. poly%cohort%n == 0_ik, 'the bare-ground start did not open, or not empty')
      call check(all(poly%patch%cas(1:poly%patch%n)%can_depth == FLOOR),                                  &
                 'a bare-ground stand does not open at the canopy-air floor')
   end subroutine test_canopy_air_floor

end program test_init
