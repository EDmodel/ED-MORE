! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_c_api_demography -- COVERAGE FOR THE DEMOGRAPHY C-API SHIM.                         !
!                                                                                          !
! WHY THIS TEST EXISTS: outside this target the shim is built only by the optional         !
! `-DMEDS_BUILD_PYLIB=ON` library, so a signature change under it (in `apply_recruitment`, !
! say) leaves the C API unable to compile while every other test passes. Every C-API shim  !
! is compiled into a ctest target so that an interface mismatch fails the build            !
! (docs/building.md §Running the tests).                                                   !
!                                                                                          !
! Two things are therefore load-bearing here and neither is an assertion:                  !
!   * the shim is COMPILED by a mandatory target, so a signature or struct change is a     !
!     build break;                                                                         !
!   * the test runs from the SOURCE directory against the SHIPPED example config, so that  !
!     config has a build-time consumer. A config nothing loads falls behind the config     !
!     reader unnoticed, exactly as a shim nothing compiles falls behind its interfaces.    !
!                                                                                          !
! The assertions cover what neither of those can: that the opaque-handle registry hands out         !
! distinct handles, releases them, and that a site actually advances behind one; and that the       !
! caller-supplied rates of example 03's path (meds_apply_rates) act on the stand as rates: growth    !
! moves each cohort's dbh by g dt, mortality thins each cohort and the stand's biomass by            !
! exp(-m dt), and recruitment waits in its pool until the month's boundary founds it, r dt a step.   !
!==========================================================================================!
program test_c_api_demography
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char, c_long
   use meds_kinds,            only : wp
   use meds_c_api_demography,  only : meds_config_load, meds_config_n_pft, meds_config_dt_years,   &
                                     meds_site_create, meds_site_init_bare, meds_site_free,       &
                                     meds_site_init_census,                                       &
                                     meds_advance_slow, meds_site_n_patch, meds_site_n_cohort,    &
                                     meds_site_total_agb, meds_site_generation, meds_apply_rates, &
                                     meds_site_total_nplant, meds_site_get_real, meds_site_get_int
   use meds_test_support, only : banner, check, check_close
   implicit none

   character(len=*), parameter :: CFG = 'examples/example03_demography/example_config_main.toml'
   !----- The root config names data/census_example.csv, a small pseudo census, as its census.  !
   character(len=*), parameter :: CFG_CENSUS = 'meds_config_main.toml'
   integer(c_int) :: ch, ch2, sh, sh2, npft, npatch
   real(c_double) :: dt_yr, agb0, agb1
   integer(c_long) :: gen0, gen1
   integer :: k
   character(kind=c_char), allocatable :: cpath(:)

   call banner('demography C-API shim')

   !=== 1. The SHIPPED example config loads. ================================================!
   !      This is the assertion that catches the example config falling behind the reader (a !
   !      required key it lacks): it makes this test a consumer of                           !
   !      examples/example03_demography/.                                                    !
   cpath = to_c(CFG)
   ch = meds_config_load(cpath, int(len(CFG), c_int))
   call check(ch > 0_c_int, 'the shipped example config loads through the C-API')
   npft  = meds_config_n_pft(ch)
   dt_yr = meds_config_dt_years(ch)
   call check(npft >= 1_c_int, 'config reports a sane PFT count')
   call check(dt_yr > 0.0_wp .and. dt_yr <= 1.0_wp, 'config reports a sane dt_years')

   !=== 2. The opaque-handle registry hands out DISTINCT handles and releases them. ==========!
   !      This is the shim's own state -- the one thing no kernel test can cover -- and the    !
   !      reason this file is not split further (see its header).                              !
   sh  = meds_site_create()
   sh2 = meds_site_create()
   call check(sh > 0_c_int .and. sh2 > 0_c_int, 'site_create returns valid handles')
   call check(sh /= sh2, 'two live sites get DISTINCT handles')
   call meds_site_free(sh2)
   sh2 = meds_site_create()
   call check(sh2 /= sh, 'a handle freed and re-taken does not collide with the live one')
   call meds_site_free(sh2)

   !=== 3. A site behind a handle initialises and ADVANCES. =================================!
   call meds_site_init_bare(sh, ch, 2_c_int)
   npatch = meds_site_n_patch(sh)
   call check(npatch == 2_c_int, 'init_bare honours the requested patch count')
   call check(meds_site_n_cohort(sh) >= 0_c_int, 'a bare site reports a cohort count')

   gen0 = meds_site_generation(sh)
   agb0 = meds_site_total_agb(sh)
   !----- A month of daily steps: enough to cross the monthly fuse/fission cadence once, which  !
   !      is where the recruitment path this shim drives actually does something.  -------------!
   do k = 1, 31
      call meds_advance_slow(sh, ch, merge(1_c_int, 0_c_int, k == 1), merge(1_c_int, 0_c_int, k == 1))
   end do
   gen1 = meds_site_generation(sh)
   agb1 = meds_site_total_agb(sh)
   call check(gen1 > gen0, 'advance_slow bumps the site generation counter')
   call check(agb1 >= 0.0_wp, 'total_agb is readable after advancing')

   call meds_site_free(sh)

   !=== 4. A site STARTS FROM A CENSUS: the file the config names is read and restructured. ==!
   cpath = to_c(CFG_CENSUS)
   ch2 = meds_config_load(cpath, int(len(CFG_CENSUS), c_int))
   call check(ch2 > 0_c_int, 'the root config loads through the C-API')
   sh2 = meds_site_create()
   call check(meds_site_init_census(sh2, ch2) == 1_c_int, 'init_census reads the configured census')
   call check(meds_site_n_patch(sh2) >= 1_c_int .and. meds_site_n_cohort(sh2) >= 1_c_int,           &
              'the census stand has patches and cohorts')
   call check(meds_site_total_agb(sh2) > 0.0_wp, 'the census stand carries biomass')
   call meds_site_free(sh2)

   !=== 5. CALLER-SUPPLIED RATES (meds_apply_rates), example 03's path, on the census stand. ==!
   call test_apply_rates(ch2)
   print '(a)', 'test_c_api_demography: ALL PASSED'

contains

   !----- The census stand driven by rates the caller supplies, one step of dt at a time, off the   !
   !      calendar's boundaries unless a step says otherwise. Each rate acts over its step: growth    !
   !      g [cm/yr] moves every cohort's dbh by g dt (through the on-allometry carbon the shim       !
   !      derives for the new size and the size read back from the wood); mortality m [1/yr] thins   !
   !      every cohort by exp(-m dt), and the stand's biomass with it; recruitment r [plant/m2/yr]   !
   !      fills each patch's pool by r dt and founds cohorts from it only at the month's boundary,    !
   !      so after a month the stand holds r times the month's length more stems than a twin given   !
   !      none. Cohorts are matched by their persistent id: the shim sorts every step. ---------------!
   subroutine test_apply_rates(ch)
      integer(c_int), intent(in) :: ch
      real(c_double), parameter :: G_RATE = 0.8_c_double, M_RATE = 0.2_c_double, R_RATE = 1.0_c_double
      integer, parameter :: MONTH = 30
      integer(c_int) :: sa, sb, n, np, npft
      integer(c_int), allocatable :: id0(:), id1(:)
      real(c_double), allocatable :: grow(:), mort(:), recr(:), none(:), dbh0(:), dbh1(:), npl0(:), npl1(:)
      real(c_double) :: dt, nplant0, agb0
      integer(c_long) :: gen0
      integer :: i, j, k
      dt   = meds_config_dt_years(ch)
      npft = meds_config_n_pft(ch)
      sa = meds_site_create()
      call check(meds_site_init_census(sa, ch) == 1_c_int, 'apply_rates: the census stand starts')
      n  = meds_site_n_cohort(sa)
      np = meds_site_n_patch(sa)
      allocate(grow(n), mort(n), id0(n), id1(n), dbh0(n), dbh1(n), npl0(n), npl1(n), recr(npft * np), none(npft * np))
      none = 0.0_c_double

      !----- (a) Growth alone: each cohort's dbh moves by g dt, its density stays. ------------------!
      call cohort_state(sa, id0, dbh0, npl0)
      grow = G_RATE
      mort = 0.0_c_double
      gen0 = meds_site_generation(sa)
      call meds_apply_rates(sa, ch, grow, mort, none, 0_c_int, 0_c_int)
      call check(meds_site_generation(sa) == gen0 + 1_c_long, 'apply_rates bumps the site generation once')
      call check(meds_site_n_cohort(sa) == n, 'apply_rates: off the boundaries the cohorts stay as they are')
      call cohort_state(sa, id1, dbh1, npl1)
      do i = 1, n
         j = findloc(id1, id0(i), dim=1)
         call check(j > 0, 'apply_rates: every cohort keeps its id')
         call check_close(dbh1(j) - dbh0(i), G_RATE * dt, 1.0e-8_wp, 'apply_rates: growth moves a cohort''s dbh by g dt')
         call check_close(npl1(j), npl0(i), 1.0e-14_wp, 'apply_rates: growth alone keeps a cohort''s density')
      end do

      !----- (b) Mortality alone: every cohort, so the stand's stems and biomass, thin by exp(-m dt). -!
      grow = 0.0_c_double
      mort = M_RATE
      nplant0 = meds_site_total_nplant(sa)
      agb0    = meds_site_total_agb(sa)
      call cohort_state(sa, id0, dbh0, npl0)
      call meds_apply_rates(sa, ch, grow, mort, none, 0_c_int, 0_c_int)
      call cohort_state(sa, id1, dbh1, npl1)
      do i = 1, n
         j = findloc(id1, id0(i), dim=1)
         call check_close(npl1(j), npl0(i) * exp(-M_RATE * dt), 1.0e-12_wp,                                       &
                          'apply_rates: mortality thins a cohort by exp(-m dt)')
         call check_close(dbh1(j), dbh0(i), 1.0e-12_wp, 'apply_rates: mortality alone keeps a cohort''s size')
      end do
      call check_close(meds_site_total_nplant(sa), nplant0 * exp(-M_RATE * dt), 1.0e-12_wp,                       &
                       'apply_rates: mortality thins the stand''s stems by exp(-m dt)')
      call check_close(meds_site_total_agb(sa), agb0 * exp(-M_RATE * dt), 1.0e-12_wp,                             &
                       'apply_rates: and its biomass with them, the survivors'' unchanged')

      !----- (c) Recruitment: r for the first PFT in every patch of one twin, none in the other, over  !
      !      a month whose last step is the boundary. -----------------------------------------------!
      sb = meds_site_create()
      call check(meds_site_init_census(sb, ch) == 1_c_int, 'apply_rates: the twin census stand starts')
      call meds_site_free(sa)
      sa = meds_site_create()
      call check(meds_site_init_census(sa, ch) == 1_c_int, 'apply_rates: the census stand starts again')
      recr = 0.0_c_double
      do j = 1, np
         recr((j - 1) * npft + 1) = R_RATE
      end do
      grow = 0.0_c_double
      mort = 0.0_c_double
      do k = 1, MONTH
         call meds_apply_rates(sa, ch, grow, mort, recr, merge(1_c_int, 0_c_int, k == MONTH), 0_c_int)
         call meds_apply_rates(sb, ch, grow, mort, none, merge(1_c_int, 0_c_int, k == MONTH), 0_c_int)
         if (k == MONTH - 1) call check_close(meds_site_total_nplant(sa), meds_site_total_nplant(sb), 1.0e-14_wp,  &
                                              'apply_rates: recruits wait in the pool until the month''s boundary')
      end do
      call check(meds_site_n_cohort(sa) > meds_site_n_cohort(sb), 'apply_rates: the boundary founds the recruits')
      call check_close(meds_site_total_nplant(sa) - meds_site_total_nplant(sb), R_RATE * MONTH * dt, 1.0e-9_wp,    &
                       'apply_rates: a month of recruitment founds r times the month''s length in stems')
      call meds_site_free(sa)
      call meds_site_free(sb)
   end subroutine test_apply_rates

   !----- A site's cohorts: their persistent ids, dbh [cm] and density [plant/m2]. --------------!
   subroutine cohort_state(sh, id, dbh, nplant)
      integer(c_int), intent(in)  :: sh
      integer(c_int), intent(out) :: id(:)
      real(c_double), intent(out) :: dbh(:), nplant(:)
      call meds_site_get_int(sh, 2_c_int, id)
      call meds_site_get_real(sh, 0_c_int, dbh)
      call meds_site_get_real(sh, 2_c_int, nplant)
   end subroutine cohort_state

   !----- A NUL-free c_char array; the shim takes (buffer, length), not a C string. ----------!
   function to_c(s) result(a)
      character(len=*), intent(in)        :: s
      character(kind=c_char), allocatable :: a(:)
      integer :: i
      allocate(a(len(s)))
      do i = 1, len(s)
         a(i) = s(i:i)
      end do
   end function to_c

end program test_c_api_demography
