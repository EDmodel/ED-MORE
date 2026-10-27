! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_main -- the program behind the edmore executable: parse the command line, drive the
! run, and turn a run status into an exit code.
!
!   The model lives in ed_region (open / step / finalize), and this program is one of its
! two callers; the other is the C-API shim ed_c_api_run, which lets Python drive the
! coupled model (edmore.model). Anything added here is unreachable from Python, so the
! program stays this thin.
!
!   Usage:  edmore [main.toml]      (default ./ed_config_main.toml; a missing file or
!                                     required key is an error: docs/configuration.md)
!         edmore --dump-io-config [main.toml]
!                                     write ed_io_config.toml, listing every output
!                                     variable with a comment, and exit
!---~---
program ed_main
   use ed_kinds,              only : ik
   use ed_config,             only : ed_config_t
   use ed_config_io,          only : load_ed_config
   use ed_output_types,       only : output_registry_t
   use ed_output_registry,    only : build_output_registry, dump_io_config
   use ed_region,             only : ed_region_t, region_open, region_step_month, region_done,  &
                                     region_finalize, region_free, RUN_OK, RUN_ERR_NAN, RUN_ERR_AREA, &
                                     RUN_ERR_BOUNDS
   implicit none

   type(ed_region_t)   :: reg
   integer(ik)         :: status, step_status
   logical             :: ok
   character(len=256)  :: path

   path = 'ed_config_main.toml'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)

   !---~---
   !   `--dump-io-config [main.toml]` is the discoverability half of per-variable output
   ! control: it lists every name an override can switch. It builds the registry and exits
   ! without running anything, so it stays here in the CLI rather than in the run.
   !---~---
   if (trim(path) == '--dump-io-config') then
      block
         type(ed_config_t)       :: cfg
         type(output_registry_t) :: reg_dump
         path = 'ed_config_main.toml'
         if (command_argument_count() >= 2_ik) call get_command_argument(2, path)
         call load_ed_config(trim(path), cfg)
         write(*,'(2a)') ' config: ', trim(path)
         call build_output_registry(reg_dump, cfg)
         call dump_io_config(reg_dump, 'ed_io_config.toml')
      end block
      stop
   end if

   call region_open(trim(path), reg, ok)
   if (.not. ok) error stop 'ed_main: could not open the run'
   step_status = RUN_OK
   do while (.not. region_done(reg))
      call region_step_month(reg, step_status)
      if (step_status /= RUN_OK) exit
   end do
   !---~---
   !   The run returns a failed step as a status, so that a library caller (Python) survives it;
   ! in the executable it is fatal. The run is finalized first either way: after a failure it
   ! writes nothing more, and it closes the files, so the failed month's records stay readable.
   !---~---
   call region_finalize(reg, status)
   call region_free(reg)
   if (step_status == RUN_ERR_NAN)    error stop 'ed_main: NaN in the state or the energy/water ledger'
   if (step_status == RUN_ERR_BOUNDS) error stop 'ed_main: a cohort or soil-carbon state out of bounds'
   if (status == RUN_ERR_AREA) error stop 'ed_main: a polygon''s area was not conserved'
   write(*,'(a)') ' OK: simulation completed, area conserved, no NaNs.'

end program ed_main
