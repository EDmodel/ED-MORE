! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_main -- the executable: parse the command line, drive the run, and turn a run       !
! status into an exit code.                                                                !
!                                                                                          !
! The model lives in meds_region (open / step / finalize), and this program is one of its  !
! two callers; the other is the C-API shim meds_c_api_run, which lets Python drive the     !
! coupled model (meds.model). Anything added here is unreachable from Python, so the       !
! program stays this thin.                                                                 !
!                                                                                          !
! Usage:  meds_main [main.toml]      (default ./meds_config_main.toml; a missing file or   !
!                                     required key is an error: docs/configuration.md)     !
!         meds_main --dump-io-config [main.toml]                                           !
!                                     write meds_io_config.toml, listing every output      !
!                                     variable with a comment, and exit                    !
!==========================================================================================!
program meds_main
   use meds_kinds,            only : ik
   use meds_config,           only : meds_config_t
   use meds_config_io,        only : load_meds_config
   use meds_output_types,     only : output_registry_t
   use meds_output_registry,  only : build_output_registry, dump_io_config
   use meds_region,           only : meds_region_t, region_open, region_step_month, region_done,  &
                                     region_finalize, region_free, RUN_OK, RUN_ERR_NAN, RUN_ERR_AREA, &
                                     RUN_ERR_BOUNDS
   implicit none

   type(meds_region_t) :: reg
   integer(ik)         :: status, step_status
   logical             :: ok
   character(len=256)  :: path

   path = 'meds_config_main.toml'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)

   !----- `--dump-io-config [main.toml]` is the discoverability half of per-variable output        !
   !      control: it lists every name an override can switch. It builds the registry and exits    !
   !      without running anything, so it stays here in the CLI rather than in the run. -----------!
   if (trim(path) == '--dump-io-config') then
      block
         type(meds_config_t)     :: cfg
         type(output_registry_t) :: reg_dump
         path = 'meds_config_main.toml'
         if (command_argument_count() >= 2_ik) call get_command_argument(2, path)
         call load_meds_config(trim(path), cfg)
         write(*,'(2a)') ' config: ', trim(path)
         call build_output_registry(reg_dump, cfg)
         call dump_io_config(reg_dump, 'meds_io_config.toml')
      end block
      stop
   end if

   call region_open(trim(path), reg, ok)
   if (.not. ok) error stop 'meds_main: could not open the run'
   step_status = RUN_OK
   do while (.not. region_done(reg))
      call region_step_month(reg, step_status)
      if (step_status /= RUN_OK) exit
   end do
   !----- The run returns a failed step as a status, so that a library caller (Python) survives it;  !
   !      in the executable it is fatal. The run is finalized first either way: after a failure it    !
   !      writes nothing more, and it closes the files, so the failed month's records stay readable.  !
   call region_finalize(reg, status)
   call region_free(reg)
   if (step_status == RUN_ERR_NAN)    error stop 'meds_main: NaN in the state or the energy/water ledger'
   if (step_status == RUN_ERR_BOUNDS) error stop 'meds_main: a cohort or soil-carbon state out of bounds'
   if (status == RUN_ERR_AREA) error stop 'meds_main: a polygon''s area was not conserved'
   write(*,'(a)') ' OK: simulation completed, area conserved, no NaNs.'

end program meds_main
