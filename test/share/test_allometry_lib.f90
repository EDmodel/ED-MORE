! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_allometry_lib -- meds_allometry_lib, the size allometry (docs/science/plant_carbon_allocation.md §2). !
!   1. THE LAWS: dbh <-> height in both regimes (below and at the cap), AGB and leaf area     !
!      increasing with size and wood density, agb_to_dbh and wood_to_dbh inverting their      !
!      forward laws, the carbon targets as re-expressions of the size allometry, Chave et al.  !
!      (2014) eq. 6a against BIOMASS's retrieveH, a local log-quadratic fit, and one record    !
!      per PFT.                                                                                !
!   2. THE SHIPPED DEFAULTS: meds_config_pft.toml's [allometry] block equals the library's    !
!      initializers, and those are Chave et al. (2014) eq. 4 in carbon and ED2's BAAD leaf fit !
!      over C2B. The initializers exist for the Python library, which reaches the allometry    !
!      without a config; a model run takes every value from the PFT config, so the two copies  !
!      can drift apart silently. Both fits are dry mass or its equivalent and MEDS keeps        !
!      carbon: taken without the division by C2B, AGB is dry mass labelled as carbon and leaf  !
!      area twice ED2's.                                                                       !
! The path of meds_config_pft.toml is the first argument (ctest passes the source tree's).  !
!==========================================================================================!
program test_allometry_lib
   use meds_kinds,         only : wp, ik
   use meds_toml,          only : toml_table_t, toml_parse_file, toml_has, toml_real, toml_string, toml_real_array
   use meds_allometry_lib, only : dbh_to_height, height_to_dbh, dbh_to_agb, agb_to_dbh, dbh_to_leaf_area,  &
                                  allometry_t, size_to_leaf_carbon, size_to_wood_carbon, wood_to_dbh,      &
                                  HEIGHT_CHAVE, HEIGHT_POWER
   use meds_test_assert,   only : banner, check, check_close
   implicit none

   real(wp), parameter :: C2B              = 2.0_wp          ! ED2 C2B, carbon -> biomass
   real(wp), parameter :: CHAVE14_SCALE    = 0.0673_wp       ! Chave et al. (2014) eq. 4 [kg dry mass]
   real(wp), parameter :: CHAVE14_EXPONENT = 0.976_wp        ! Chave et al. (2014) eq. 4
   real(wp), parameter :: ED2_BL_XX_SCALE  = 0.46769540_wp   ! ED2 c14f15_bl_xx(1)
   real(wp), parameter :: RTOL             = 1.0e-12_wp
   type(toml_table_t)  :: t          !< the shipped PFT file

   call test_allometry_laws()
   call test_shipped_defaults()
   write(*,'(a)') 'test_allometry_lib: ALL PASSED'

contains

   !----- 1. The laws. ----------------------------------------------------------------------!
   subroutine test_allometry_laws()
      real(wp), parameter :: RHO  = 0.60_wp
      real(wp), parameter :: HCAP = 46.0_wp        ! the record's height cap (its initializer)
      real(wp), parameter :: SLA  = 13.0_wp        ! [m2/kgC] specific leaf area
      real(wp), parameter :: AGF  = 0.70_wp        ! aboveground fraction of woody carbon
      real(wp) :: dbh, h, agb, dbh_back, wood_c
      integer(ik) :: i
      type(allometry_t) :: allom, two(2)
      real(wp)          :: hh(2)
      real(wp), dimension(3) :: dvec = [ 5.0_wp, 30.0_wp, 80.0_wp ]


      call banner('the allometry laws')
      !----- The shipped pan-tropical coefficients: the record's initializers. ----------------!
      allom = allometry_t()
      call check(allom%hgt_max == HCAP .and. allom%aboveground_frac == AGF, 'the record''s cap and aboveground share')

      !=== Height <-> diameter round-trip below the height cap. ==============================!
      do i = 1_ik, 3_ik
         dbh = dvec(i)
         h   = dbh_to_height(allom, dbh)
         call check(h < HCAP, 'test sizes must stay below the height cap')
         call check_close(height_to_dbh(allom, h), dbh, 1.0e-9_wp, 'dbh<->height round-trip failed')
      end do

      !=== Height saturates at the cap for very large diameters. ============================!
      call check_close(dbh_to_height(allom, 500.0_wp), HCAP, 1.0e-12_wp, 'height did not cap at hgt_max')

      !=== AGB and leaf area increase monotonically with diameter. ==========================!
      call check(dbh_to_agb(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp), RHO) >                              &
                 dbh_to_agb(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp), RHO), 'AGB not increasing with dbh')
      call check(dbh_to_leaf_area(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp)) >                             &
                 dbh_to_leaf_area(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp)), 'leaf area not increasing with dbh')
      call check(dbh_to_agb(allom, 30.0_wp, dbh_to_height(allom, 30.0_wp), RHO) > 0.0_wp, 'AGB must be positive')

      !=== AGB heavier for denser wood at the same size. ====================================!
      call check(dbh_to_agb(allom, 30.0_wp, dbh_to_height(allom, 30.0_wp), 0.90_wp) >                          &
                 dbh_to_agb(allom, 30.0_wp, dbh_to_height(allom, 30.0_wp), 0.40_wp), 'denser wood should weigh more')

      !=== agb_to_dbh inverts dbh_to_agb in BOTH regimes. ===================================!
      !----- uncapped (small): height(30) < hcap. ------------------------------------------------!
      dbh      = 30.0_wp
      agb      = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), RHO)
      dbh_back = agb_to_dbh(allom, agb, RHO)
      call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (uncapped)')
      !----- capped (large): height(130) saturates at hcap. --------------------------------------!
      dbh      = 130.0_wp
      call check_close(dbh_to_height(allom, dbh), HCAP, 1.0e-12_wp, 'large stem should be at the cap')
      agb      = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), RHO)
      dbh_back = agb_to_dbh(allom, agb, RHO)
      call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (capped)')

      !=== Carbon-pool targets: EXACT re-expressions of the size allometry. =====================!
      do i = 1_ik, 3_ik
         dbh = dvec(i)
         h   = dbh_to_height(allom, dbh)
         !----- leaf_carbon * SLA reproduces the leaf area exactly. -------------------------------!
         call check_close(size_to_leaf_carbon(allom, dbh, h, SLA) * SLA, dbh_to_leaf_area(allom, dbh, h), 1.0e-12_wp, &
                          'size_to_leaf_carbon * sla /= leaf_area')
         !----- aboveground_frac * wood_carbon reproduces the Chave AGB exactly. ------------------!
         call check_close(size_to_wood_carbon(allom, dbh, h, RHO) * AGF, dbh_to_agb(allom, dbh, h, RHO), 1.0e-12_wp, &
                          'aboveground_frac * wood_carbon /= agb')
      end do

      !=== Carbon targets increase with diameter. ==========================================!
      call check(size_to_leaf_carbon(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp), SLA) >                         &
                 size_to_leaf_carbon(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp), SLA), 'leaf carbon not increasing')
      call check(size_to_wood_carbon(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp), RHO) >                         &
                 size_to_wood_carbon(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp), RHO), 'wood carbon not increasing')

      !=== wood_to_dbh inverts size_to_wood_carbon in BOTH regimes (the carbon anchor). =============!
      dbh    = 30.0_wp                                        ! uncapped
      wood_c = size_to_wood_carbon(allom, dbh, dbh_to_height(allom, dbh), RHO)
      call check_close(wood_to_dbh(allom, wood_c, RHO), dbh, 1.0e-8_wp, 'wood_to_dbh inverse failed (uncapped)')
      dbh    = 130.0_wp                                       ! capped
      wood_c = size_to_wood_carbon(allom, dbh, dbh_to_height(allom, dbh), RHO)
      call check_close(wood_to_dbh(allom, wood_c, RHO), dbh, 1.0e-8_wp, 'wood_to_dbh inverse failed (capped)')

      !=== Chave et al. (2014) eq. 6a at E = 0.2 and its two inverses. The expected heights are     !
      !    BIOMASS's retrieveH: exp(0.893 - E + 0.760 ln D - 0.0340 (ln D)^2 + 0.243^2/2). =========!
      allom = allometry_t(height_allometry = HEIGHT_CHAVE, chave_e = 0.2_wp, hgt_max = 100.0_wp)
      call check_close(dbh_to_height(allom, 1.0_wp), 2.0596261794_wp, 1.0e-9_wp, 'Chave height at 1 cm')
      call check_close(dbh_to_height(allom, 10.0_wp), 9.8969302660_wp, 1.0e-9_wp, 'Chave height at 10 cm')
      call check_close(dbh_to_height(allom, 30.0_wp), 18.4325408427_wp, 1.0e-9_wp, 'Chave height at 30 cm')
      call check_close(dbh_to_height(allom, 60.0_wp), 26.1607976122_wp, 1.0e-9_wp, 'Chave height at 60 cm')
      call check_close(dbh_to_height(capped(allom, 30.0_wp), 100.0_wp), 30.0_wp, 1.0e-12_wp, 'Chave height caps at hgt_max')
      do i = 1_ik, 3_ik
         dbh = dvec(i)
         call check_close(height_to_dbh(allom, dbh_to_height(allom, dbh)), dbh, 1.0e-9_wp, 'Chave dbh<->height round trip')
         agb = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), RHO)
         call check_close(agb_to_dbh(allom, agb, RHO), dbh, 1.0e-9_wp, 'Chave agb_to_dbh inverse (uncapped)')
      end do
      dbh = 120.0_wp
      call check_close(dbh_to_height(capped(allom, 30.0_wp), dbh), 30.0_wp, 1.0e-12_wp, 'a 120 cm stem is at a 30 m cap')
      agb = dbh_to_agb(allom, dbh, dbh_to_height(capped(allom, 30.0_wp), dbh), RHO)
      call check_close(agb_to_dbh(capped(allom, 30.0_wp), agb, RHO), dbh, 1.0e-9_wp, 'Chave agb_to_dbh inverse (capped)')
      wood_c = size_to_wood_carbon(allom, 80.0_wp, dbh_to_height(allom, 80.0_wp), RHO)
      call check_close(wood_to_dbh(allom, wood_c, RHO), 80.0_wp, 1.0e-9_wp, 'Chave wood_to_dbh inverse')
      !----- A PFT's own log-quadratic coefficients replace Chave's: a flatter curve is shorter. --!
      call check(dbh_to_height(allometry_t(height_allometry = HEIGHT_CHAVE, chave_b = 0.6_wp), 30.0_wp) <              &
                 dbh_to_height(allometry_t(height_allometry = HEIGHT_CHAVE), 30.0_wp), 'a smaller chave_b gives shorter trees')

      !=== A local log-quadratic (Barro Colorado's: chave_a, chave_b, chave_c, no E, no RSE term)   !
      !    follows its own coefficients, and both inverses hold on it. =======================!
      allom = allometry_t(height_allometry = HEIGHT_CHAVE, chave_a = 0.8898_wp, chave_b = 0.7908_wp,          &
                          chave_c = -0.04923_wp, chave_rse = 0.0_wp, chave_e = 0.0_wp, hgt_max = 100.0_wp)
      call check_close(dbh_to_height(allom, 1.0_wp), exp(0.8898_wp), 1.0e-12_wp, 'local fit: height at 1 cm is exp(chave_a)')
      call check_close(dbh_to_height(allom, 100.0_wp),                                                          &
                       exp(0.8898_wp + 0.7908_wp * log(100.0_wp) - 0.04923_wp * log(100.0_wp)**2), 1.0e-12_wp, &
                       'local fit: height at 100 cm follows its coefficients')
      do i = 1_ik, 3_ik
         dbh = dvec(i)
         call check_close(height_to_dbh(allom, dbh_to_height(allom, dbh)), dbh, 1.0e-9_wp, 'local fit: dbh<->height')
         agb = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), RHO)
         call check_close(agb_to_dbh(allom, agb, RHO), dbh, 1.0e-9_wp, 'local fit: agb_to_dbh inverse')
      end do

      !=== Each PFT has its own record: a short log-quadratic understorey PFT and a power-law PFT     !
      !    keep their own heights and caps side by side (elemental over the records). ============!
      two(1) = allometry_t(height_allometry = HEIGHT_CHAVE, chave_a = 0.6_wp, chave_b = 0.6_wp, chave_c = -0.03_wp,  &
                           chave_rse = 0.0_wp, hgt_max = 12.0_wp)
      two(2) = allometry_t()
      hh = dbh_to_height(two, 40.0_wp)
      call check_close(hh(1), min(exp(0.6_wp + 0.6_wp * log(40.0_wp) - 0.03_wp * log(40.0_wp)**2), 12.0_wp), 1.0e-12_wp, &
                       'PFT 1 follows its own log-quadratic and cap')
      call check_close(hh(2), dbh_to_height(allometry_t(), 40.0_wp), 0.0_wp, 'PFT 2 follows the power law')
      call check(hh(1) < hh(2), 'the understorey PFT is shorter at the same diameter')

   end subroutine test_allometry_laws

   !----- 2. The shipped [allometry] block and the initializers behind it. ----------------!
   subroutine test_shipped_defaults()
      character(len=512) :: arg
      logical            :: ok
      type(allometry_t)  :: d          ! a record at its initializers

      call banner('shipped allometry defaults')

      !----- CMake passes the SOURCE-tree path; ctest runs from the build tree. -------------!
      call get_command_argument(1, arg)
      call check(len_trim(arg) > 0, 'the meds_config_pft.toml path was passed on the command line')
      call toml_parse_file(trim(arg), t, ok)
      call check(ok, 'the shipped meds_config_pft.toml parses')

      !----- The Fortran initializers: Chave 2014 in carbon, and ED2's leaf fit over C2B. ----!
      call check_close(d%agb_c1, CHAVE14_SCALE / C2B, RTOL, 'agb_c1 initializer = Chave 2014 eq. 4 scale / C2B')
      call check_close(d%agb_c2, CHAVE14_EXPONENT, RTOL, 'agb_c2 initializer = Chave 2014 eq. 4 exponent')
      call check_close(d%lai_b1, ED2_BL_XX_SCALE / C2B, RTOL, 'lai_b1 initializer = ED2 c14f15_bl_xx(1)/C2B')

      !----- The shipped config equals the initializers, key by key. -----------------------!
      call same('b1Ht', d%b1Ht)
      call same('b2Ht', d%b2Ht)
      call same('agb_c1', d%agb_c1)
      call same('agb_c2', d%agb_c2)
      call same('lai_b1', d%lai_b1)
      call same('lai_b2', d%lai_b2)
      !----- The keys the shipped config lists per PFT: every entry equals the initializer. --!
      call same_per_pft('hgt_max',          d%hgt_max)
      call same_per_pft('wai_b1',           d%wai_b1)
      call same_per_pft('wai_b2',           d%wai_b2)
      call same_per_pft('sapwood_area_b1',  d%sapwood_area_b1)
      call same_per_pft('sapwood_area_b2',  d%sapwood_area_b2)
      call same_per_pft('aboveground_frac', d%aboveground_frac)
      !----- The log-quadratic initializers are Chave et al. (2014) eq. 6a. -----------------!
      call check_close(d%chave_a, 0.893_wp, RTOL, 'chave_a initializer = Chave 2014 eq. 6a')
      call check_close(d%chave_b, 0.760_wp, RTOL, 'chave_b initializer = Chave 2014 eq. 6a')
      call check_close(d%chave_c, -0.0340_wp, RTOL, 'chave_c initializer = Chave 2014 eq. 6a')
      call check_close(d%chave_rse, 0.243_wp, RTOL, 'chave_rse initializer = Chave 2014 eq. 6a')

      !----- The height curve: the shipped config and the initializer are both the power law. ---!
      call check(trim(toml_string(t, 'allometry.height_allometry', '')) == 'power',                    &
                 'meds_config_pft.toml sets [allometry].height_allometry = "power"')
      call check(d%height_allometry == HEIGHT_POWER, 'the height-curve initializer is the power law')

   end subroutine test_shipped_defaults

   subroutine same(key, initializer)
      character(len=*), intent(in) :: key
      real(wp),         intent(in) :: initializer
      call check(toml_has(t, 'allometry.'//key), 'meds_config_pft.toml sets [allometry].'//key)
      call check_close(toml_real(t, 'allometry.'//key, -huge(1.0_wp)), initializer, RTOL,         &
                       'meds_config_pft.toml [allometry].'//key//' equals the meds_allometry_lib initializer')
   end subroutine same

   subroutine same_per_pft(key, initializer)
      character(len=*), intent(in) :: key
      real(wp),         intent(in) :: initializer
      real(wp)    :: buf(64)
      integer(ik) :: n
      call toml_real_array(t, 'allometry.'//key, buf, n)
      call check(n > 0, 'meds_config_pft.toml lists [allometry].'//key//' per PFT')
      call check(all(abs(buf(1:n) - initializer) <= RTOL * abs(initializer)),                         &
                 'meds_config_pft.toml [allometry].'//key//' equals the meds_allometry_lib initializer for every PFT')
   end subroutine same_per_pft

   !----- A copy of a record with another height cap. ------------------------------------!
   pure function capped(a, hgt_max) result(b)
      type(allometry_t), intent(in) :: a
      real(wp),          intent(in) :: hgt_max
      type(allometry_t)             :: b
      b = a
      b%hgt_max = hgt_max
   end function capped

end program test_allometry_lib
