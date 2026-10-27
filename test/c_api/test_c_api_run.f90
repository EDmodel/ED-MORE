! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_c_api_run -- COVERAGE FOR THE FULL-COUPLED-MODEL C-API SHIM (`ed.model`).
!
!   Two things here are load-bearing and neither is an assertion:
!
!   * `src/c_api/ed_c_api_run.f90` is COMPILED by a MANDATORY target, so a change to
!     `ed_region`'s open/step/finalize seam is a build break rather than a silent rot in
!     an optional library. Every C-API shim is compiled into a ctest target for this
!     reason (docs/building.md §Running the tests).
!
!   * this test runs from `examples/example04_column_biophysics/` and LOADS BOTH of that
!     example's shipped configs, the default and the calibrated run. The example runs them with
!     scripts the build knows nothing about, and a config nothing loads rots like a shim nothing
!     compiles.
!
!   The assertions cover the shim's own logic: the invalid-handle guard (every accessor has to
! return a sentinel rather than index a stale registry slot), distinct handles, a real
! open -> step -> finalize -> free cycle through the C boundary, and a failed open that leaves its
! slot free for the next one.
!
!   The cycle runs on a config DERIVED from the shipped default one: fast biophysics and live met
! forcing off (the forcing netCDF is not tracked in git), bare ground, three days. The derivation
! writes an override block, then the original less the keys the block sets (derive_config), since
! `ed_toml` refuses a key set twice. Deriving it, rather than writing a config from scratch,
! means this test also proves the shipped file is complete enough to INITIALISE, not merely to
! parse. The derived configs go into the work directory (the first argument), the only place the
! test writes.
!---~---
program test_c_api_run
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char
   use ed_kinds,          only : wp
   use ed_config,         only : ed_config_t
   use ed_config_io,      only : load_ed_config
   use ed_c_api_run,     only : ed_run_open, ed_run_step, ed_run_is_done,                 &
                                 ed_run_finalize, ed_run_free, ed_run_year, ed_run_day,          &
                                 ed_run_istep, ed_run_n_patch, ed_run_n_cohort,                  &
                                 ed_run_total_agb, ed_run_total_lai, ed_run_soil_carbon,         &
                                 ed_run_n_site, ed_run_get_patch_int, ed_run_get_site_real
   use ed_test_support,   only : banner, check
   use ed_test_run,       only : derive_config
   implicit none

   character(len=*), parameter :: CFG_DEFAULT    = 'ed_config_eval.toml'
   character(len=*), parameter :: CFG_CALIBRATED = 'calibration/ed_config_calibrated.toml'
   character(len=512)            :: work
   character(len=:), allocatable :: cfg_derived, cfg_no_census, pft_seeded, pft_file
   integer                       :: npft
   integer(c_int)      :: h, h2, status, d0, d1, nsteps
   real(c_double)      :: agb, lai, soilc
   type(ed_config_t)   :: cfg

   call banner('full-model C-API shim (ed.model)')
   call get_command_argument(1, work)
   if (len_trim(work) == 0) error stop 'test_c_api_run: pass the work directory'
   call execute_command_line('mkdir -p "'//trim(work)//'"')
   cfg_derived   = trim(work)//'/derived.toml'
   cfg_no_census = trim(work)//'/no_census.toml'
   pft_seeded    = trim(work)//'/pft_seeded.toml'

   !---~---
   !   1. Both SHIPPED example04_column_biophysics configs load.
   ! The schema-drift guard: the example's scripts are outside the build, so nothing else
   ! here would notice its configs falling behind the config reader.
   !---~---
   call load_ed_config(CFG_DEFAULT, cfg)
   call check(cfg%pft%n >= 1_c_int, 'the shipped default config loads and carries PFTs')
   npft     = int(cfg%pft%n)
   pft_file = trim(cfg%init%pft_config)
   call check(cfg%fast_biophysics_on, 'the default config requests fast biophysics')
   call check(cfg%soil_carbon%soil_carbon_on, 'the default config requests SOIL CARBON')
   call load_ed_config(CFG_CALIBRATED, cfg)
   call check(cfg%pft%n >= 1_c_int, 'the shipped calibrated config loads and carries PFTs')
   call check(cfg%fast_biophysics_on, 'the calibrated config requests fast biophysics')
   call check(cfg%soil_carbon%soil_carbon_on, 'the calibrated config requests SOIL CARBON')

   !---~---
   !   2. The invalid-handle guard.
   ! Every accessor takes a bare integer from the caller. One that was never opened, or
   ! already freed, must return a sentinel -- not index a stale registry slot.
   !---~---
   call check(ed_run_step(99_c_int)        == -1_c_int, 'step on a bogus handle returns -1')
   call check(ed_run_is_done(99_c_int)     ==  1_c_int, 'is_done on a bogus handle reports done')
   call check(ed_run_finalize(0_c_int)     == -1_c_int, 'finalize on handle 0 returns -1')
   call check(ed_run_year(-3_c_int)        == -1_c_int, 'a negative handle yields the year sentinel')
   call check(ed_run_total_agb(99_c_int)   == 0.0_wp,   'aggregates on a bogus handle read 0')

   !--- 3. A real open -> step -> finalize -> free cycle across the C boundary.
   call write_seeded_pft(pft_seeded, pft_file, npft)
   call write_derived_config(cfg_derived, pft=pft_seeded)
   h = ed_run_open(to_c(cfg_derived), int(len(cfg_derived), c_int), 0_c_int)
   call check(h > 0_c_int, 'a run opens from a config through the C-API')
   call check(ed_run_is_done(h) == 0_c_int, 'a freshly opened run is not already done')
   call check(ed_run_n_patch(h) >= 1_c_int, 'the opened run reports its patch count')
   block
      integer(c_int) :: owner(64)
      real(c_double) :: site_area(1)
      call check(ed_run_n_site(h) == 1_c_int, 'the opened run reports its one site')
      call ed_run_get_site_real(h, 0_c_int, site_area)
      call check(site_area(1) == 1.0_c_double, 'its one site is the whole polygon')
      call ed_run_get_patch_int(h, 0_c_int, owner)
      call check(all(owner(1:ed_run_n_patch(h)) == 1_c_int), 'every patch is in the one site')
   end block

   h2 = ed_run_open(to_c(cfg_derived), int(len(cfg_derived), c_int), 0_c_int)
   call check(h2 > 0_c_int .and. h2 /= h, 'two live runs get DISTINCT handles')
   call ed_run_free(h2)
   call check(ed_run_is_done(h2) == 1_c_int, 'a freed handle reads as done, not as a live run')

   !--- Step to the end. The CALLER owns the loop -- that is the whole point of this shim.
   d0 = ed_run_day(h)
   nsteps = 0_c_int
   do while (ed_run_is_done(h) == 0_c_int)
      status = ed_run_step(h)
      call check(status == 0_c_int, 'each step reports RUN_OK')
      nsteps = nsteps + 1_c_int
      if (nsteps > 10_c_int) exit                     ! the derived config is a 3-day window
   end do
   d1 = ed_run_day(h)
   call check(nsteps == 3_c_int, 'a 3-day window takes exactly 3 daily steps')
   call check(d1 /= d0, 'the calendar advanced')
   call check(ed_run_istep(h) == nsteps, 'the driver istep counter matches the steps taken')

   !--- Stepping past the end is a STATUS, not a crash: the Python loop relies on it.
   call check(ed_run_step(h) == 1_c_int, 'a step past end_time returns RUN_FINISHED')

   !---~---
   !   3b. A REUSED registry slot must not inherit the previous run's counters.
   ! The registry is module-`save`, so a freed handle is handed out again. `ed_region_t`'s
   ! default initialisers apply to a fresh variable, not to a slot being reopened: what makes
   ! the slot fresh is region_free, which leaves it as a run never opened.
   !
   !   A Python process that runs one config after another (a calibration's trials) reuses
   ! slots this way. The step counter, which the C API reads, is what is checked: a slot
   ! that kept it would have kept the rest of the old run too.
   !---~---
   call ed_run_free(h)
   h = ed_run_open(to_c(cfg_derived), int(len(cfg_derived), c_int), 0_c_int)
   call check(h > 0_c_int, 'a second run opens after the first is freed')
   call check(ed_run_istep(h) == 0_c_int,                                                        &
              'a reopened registry slot starts at step 0, not the previous run''s count')
   call check(ed_run_is_done(h) == 0_c_int, 'and is not already done')
   status = ed_run_step(h)
   call check(ed_run_istep(h) == 1_c_int, 'the second run counts its own steps from 1')

   !---~---
   !   3c. A run that FAILS to open holds nothing.
   ! A census start whose file is missing stops in the polygons' start, after the run has
   ! allocated them. The slot it failed in is the next one handed out, so the open after
   ! it must start clean rather than find the failed attempt's arrays still there.
   !---~---
   call write_derived_config(cfg_no_census, census='no_such_census.csv')
   h2 = ed_run_open(to_c(cfg_no_census), int(len(cfg_no_census), c_int), 0_c_int)
   call check(h2 == -2_c_int, 'a run whose census file is missing does not open')
   h2 = ed_run_open(to_c(cfg_derived), int(len(cfg_derived), c_int), 0_c_int)
   call check(h2 > 0_c_int .and. h2 /= h, 'a run opens in the slot a failed open left')
   call check(ed_run_step(h2) == 0_c_int, 'and steps from its own start')
   call ed_run_free(h2)

   !--- 4. State is readable, and readable AFTER finalize (that is why free is separate).
   agb   = ed_run_total_agb(h)
   lai   = ed_run_total_lai(h)
   soilc = ed_run_soil_carbon(h)
   !---~---
   !   The derived run's bare ground is seeded (write_seeded_pft), so the stepped stand has
   ! cohorts, biomass, leaves and soil carbon for the shim to read back.
   !---~---
   call check(agb > 0.0_wp .and. lai > 0.0_wp, 'the seeded stand reads back its biomass and leaves')
   call check(soilc > 0.0_wp, 'and its soil carbon')
   call check(ed_run_n_cohort(h) > 0_c_int, 'and its cohorts')

   status = ed_run_finalize(h)
   call check(status == 0_c_int, 'finalize reports no area-conservation failure')
   call check(ed_run_total_agb(h) == agb, 'the state is still readable after finalize')
   call ed_run_free(h)
   call check(ed_run_total_agb(h) == 0.0_wp, 'a freed handle no longer reads state')

   print '(a)', 'test_c_api_run: ALL PASSED'

contains

   !---~---
   !   The shipped default config under an override block, its own values of the block's keys
   ! left out (derive_config). What is turned off: the fast loop and live met forcing (the
   ! forcing netCDF is untracked), the output streams and the state checkpoint, and the window
   ! is cut to three days from bare ground, or from `census`, a census file. The shipped
   ! config's own paths stay relative to the example's directory, the test's working
   ! directory.
   !---~---
   subroutine write_derived_config(path, census, pft)
      character(len=*),           intent(in) :: path
      character(len=*), optional, intent(in) :: census, pft
      character(len=1), parameter :: NL = new_line('a')
      character(len=:), allocatable :: init
      if (present(census)) then
         init = 'init_mode = 1'//NL//'census_file = "'//census//'"'//NL
      else
         init = 'init_mode = 0'//NL
      end if
      if (present(pft)) init = init//'pft_config = "'//pft//'"'//NL
      call derive_config(path, '# GENERATED by test_c_api_run.'//NL//                                   &
                         '[run]'//NL//'start_time = "2024-07-01"'//NL//'end_time   = "2024-07-04"'//NL//      &
                         'dt_slow    = "1d"'//NL//'[fast]'//NL//'fast_biophysics_on = false'//NL//            &
                         'gpp_ref = 0.3'//NL//'[forcing]'//NL//'forcing_on = false'//NL//'[init]'//NL//init// &
                         '[output]'//NL//'state_write = false'//NL//'enabled = false', CFG_DEFAULT)
   end subroutine write_derived_config

   !---~---
   !   The example's PFT file with every PFT seeding a bare-ground start (pft.nbg_init_density,
   ! 0.1 plants/m2 per patch), so the derived run's stand has cohorts.
   !---~---
   subroutine write_seeded_pft(path, source, npft)
      character(len=*), intent(in) :: path, source
      integer,          intent(in) :: npft
      integer :: uin, uout, ios
      character(len=2048) :: line
      open(newunit=uout, file=path, status='replace', action='write')
      open(newunit=uin, file=source, status='old', action='read')
      do
         read(uin,'(a)',iostat=ios) line
         if (ios /= 0) exit
         write(uout,'(a)') trim(line)
         if (trim(adjustl(line)) == '[pft]')                                                      &
            write(uout,'(3a)') 'nbg_init_density = [', repeat('0.1, ', npft - 1), '0.1]'
      end do
      close(uin)
      close(uout)
   end subroutine write_seeded_pft

   !--- A NUL-free c_char array; the shim takes (buffer, length), not a C string.
   function to_c(s) result(a)
      character(len=*), intent(in)        :: s
      character(kind=c_char), allocatable :: a(:)
      integer :: i
      allocate(a(len(s)))
      do i = 1, len(s)
         a(i) = s(i:i)
      end do
   end function to_c

end program test_c_api_run
