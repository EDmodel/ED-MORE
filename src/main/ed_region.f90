! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_region -- a run: one polygon or many, in one process                                  !
! (docs/science/order_of_processes.md §1; the configuration is docs/configuration.md).     !
! Every run goes through here: ed_main steps it a month at a time, and the C-API shim      !
! ed_c_api_run (Python's ed.model) one slow step at a time.                                !
!                                                                                          !
! [run].mode chooses the polygons: "site" (the default) is one polygon at the [site]       !
! location, on either forcing format; "region" is every selected cell of a box of the      !
! ED_ERA5land archive. The forcing loads a month for every cell at once; then each polygon !
! runs through it, reading its own cell; then the I/O phase writes every polygon's queued  !
! records. Nothing inside the polygon loop touches a file, so the polygons of a month run  !
! side by side on [run].n_threads threads, each with a one-thread patch loop; a run of one !
! polygon gives its patch loop the threads instead.                                        !
!                                                                                          !
! Output and state files carry a polygon axis, the cohort and patch axes padded to the     !
! largest polygon's. A failed step is returned as a status, not an `error stop`: this      !
! module is in the shared library Python loads, where an `error stop` ends the             !
! interpreter. ed_main makes the status fatal.                                             !
!==========================================================================================!
module ed_region
   use ed_kinds,                    only : wp, ik
   use ed_constants,                only : day_sec, yr_day
   use ed_numerics,                 only : ascending_order
   use ed_config,                   only : ed_config_t, INIT_RESTART,                                &
                                           RUN_MODE_REGION, soil_maps_on, site_configs, apply_map_inputs
   use ed_config_check,             only : validate_config
   use ed_config_io,                only : load_ed_config, write_pft_params_csv
   use ed_time,                     only : ed_time_t, calendar_step_t, calendar_step, time_lt,       &
                                           time_to_string, years_between
   use ed_polygon_state_vars,       only : polygon_t, free_polygon
   use ed_init_restart,             only : restart_file_t, restart_file_open, restart_file_write, restart_file_close, &
                                           restart_file_slots
   use ed_init,                     only : init_polygon
   use ed_met_types,                only : met_source_t, met_grid_t, met_cells_t
   use ed_forcing_params,           only : MET_BACKEND_ERA5LAND
   use ed_met,                      only : met_open, met_open_co2, met_close, met_prefetch
   use ed_met_share,                only : MET_OK
   use ed_maps,                     only : apply_cell_mask, map_status_message, MAP_OK, read_map_inputs
   use ed_maps_params,              only : map_inputs_t
   use ed_met_grid,                 only : select_box_cells, select_site_cell, GRID_OK, GRID_ERR_NO_CELL
   use ed_met_era5land,             only : era5land_grid
   use ed_output_reduce,            only : print_summary, total_area, total_agb, total_lai
   use ed_polygon,                  only : polygon_run_t, polygon_prepare, polygon_step, polygon_report,  &
                                           RUN_OK, RUN_FINISHED, RUN_ERR_NAN, RUN_ERR_AREA, RUN_ERR_BOUNDS
   use ed_output_types,             only : output_files_t, output_buffers_t
   use ed_output_params,            only : output_params_t, STATE_DAYS, STATE_MONTHS, STATE_YEARS
   use ed_output_stream,            only : output_serialize, output_close
   use ed_output_setup,             only : open_output_files, attach_output
   use ed_toml,                     only : toml_write_record
   implicit none
   private

   public :: ed_region_t, region_open, region_step, region_step_month, region_done, region_finalize
   public :: region_free
   !----- The status codes live with the step (ed_polygon); re-exported for the run's              !
   !      callers. --------------------------------------------------------------------------------!
   public :: RUN_OK, RUN_FINISHED, RUN_ERR_NAN, RUN_ERR_AREA, RUN_ERR_BOUNDS

   !----- How far a polygon's patch areas may sum from one at the end of a run                     !
   !      (RUN_ERR_AREA). ----------------------------------------------------------------------!
   real(wp), parameter :: AREA_TOL = 1.0e-5_wp

   !----- Everything a run needs between steps. No state hides in module scope, so two runs can be !
   !      open at once (the C API's handles). -----------------------------------------------------!
   type :: ed_region_t
      type(ed_config_t)     :: cfg
      type(met_source_t)    :: met_src          !< one reader for every polygon (opened only if forcing_on)
      type(output_files_t)  :: out_files           !< the run's files (built only if output.enabled)
      !----- Each polygon as its state (the stand) and what the run keeps to step and check it     !
      !      (ed_polygon). ----------------------------------------------------------------------  !
      type(polygon_t),     allocatable :: poly(:)
      type(polygon_run_t), allocatable :: run(:)
      !----- Each polygon's share of the files, one per polygon. A contiguous array of its own, so !
      !      the I/O phase passes it whole: a section run(:)%bufs of a type with allocatable       !
      !      components makes gfortran copy it through a temporary whose copy-out dangles          !
      !      them. --------------------------------------------------------------------------------!
      type(output_buffers_t), allocatable :: out_bufs(:)
      type(ed_time_t)       :: now
      integer(ik)           :: step_days = 1_ik
      integer(ik)           :: istep = 0_ik       !< slow steps taken
      integer(ik)           :: iyear = 0_ik       !< year boundaries stepped across, for the state files
      integer(ik)           :: imonth = 0_ik      !< month boundaries stepped across, likewise
      !----- The forcing's load time (met_prefetch) over the run and this year, and the archive    !
      !      month files read before this year, for the summaries. --------------------------------!
      real(wp)              :: read_seconds = 0.0_wp, year_read_seconds = 0.0_wp
      integer(ik)           :: year_loads0 = 0_ik
      !----- Last month's wall time per polygon: the next month starts the costliest first. -------!
      real(wp),    allocatable :: cost(:)
      !----- The restructuring the calendar owes after the last step, the same for every polygon   !
      !      still running: a state file takes it from here, not from a polygon that has           !
      !      failed. ------------------------------------------------------------------------------!
      logical               :: restructure_pending = .false., restructure_new_year = .false.
      logical               :: verbose = .true.
   end type ed_region_t

contains

   pure logical function region_done(reg)
      type(ed_region_t), intent(in) :: reg
      region_done = .not. time_lt(reg%now, reg%cfg%end_time)
   end function region_done

   !---------------------------------------------------------------------------------------!
   ! region_open -- everything before the first step: the configuration, the polygons (one !
   ! at the [site], or the box's cells), the forcing, each polygon's start and             !
   ! preparation, and the files. A run that cannot open is freed, so the variable (a C-API !
   ! handle's slot) can open another.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine region_open(path, reg, ok, verbose)
      character(len=*),    intent(in)    :: path
      type(ed_region_t),   intent(inout) :: reg
      logical,             intent(out)   :: ok
      logical, optional,   intent(in)    :: verbose
      if (present(verbose)) reg%verbose = verbose
      call open_run(path, reg, ok)
      if (.not. ok) call region_free(reg)
   end subroutine region_open

   !----- Everything region_open does, in order; ok = .false. (after saying why) when the run      !
   !      cannot start. ---------------------------------------------------------------------------!
   subroutine open_run(path, reg, ok)
      character(len=*),    intent(in)    :: path
      type(ed_region_t),   intent(inout) :: reg
      logical,             intent(out)   :: ok
      type(met_cells_t)          :: cells, places
      type(map_inputs_t)         :: inputs
      integer(ik), allocatable   :: ids(:), slot(:)
      real(wp)                   :: distance_km
      character(len=:), allocatable :: msg
      logical                    :: have_maps

      ok = .false.
      call load_ed_config(trim(path), reg%cfg)
      if (reg%verbose) write(*,'(2a)') ' config: ', trim(path)
      associate (cfg => reg%cfg)
      !----- A region steps each polygon on an OpenMP worker thread, whose stack is OMP_STACKSIZE  !
      !      (4 MB by default with ifx): a grown stand's automatic arrays overflow it within a few !
      !      decades, a SIGSEGV. Say so when neither variable sets it. ----------------------------!
      if (cfg%mode == RUN_MODE_REGION .and. cfg%n_threads > 1_ik .and.                              &
          .not. (env_set('OMP_STACKSIZE') .or. env_set('KMP_STACKSIZE')))                           &
         write(*,'(2a)') ' WARNING: region: OMP_STACKSIZE is not set, and a grown stand overflows a ', &
                         'worker thread''s default stack (4 MB with ifx); set it, e.g. OMP_STACKSIZE=512M'

      if (.not. choose_polygons(cfg, reg%verbose, cells, places, distance_km, inputs, have_maps)) return
      allocate(ids(places%ncell), source=0_ik)
      if (cfg%mode == RUN_MODE_REGION) ids = places%row * places%nlon + places%col   ! the cells' row-major index

      !----- A restart (init_mode = 2): each polygon's row of the state file, by its id. A file    !
      !      may hold more polygons than the run: a region whose cell failed resumes with that     !
      !      cell masked out. ---------------------------------------------------------------------!
      allocate(slot(places%ncell), source=0_ik)
      if (cfg%init%init_mode == INIT_RESTART) then
         call restart_file_slots(trim(cfg%init%restart_file), ids, slot, msg)
         if (len(msg) > 0) then
            write(*,'(4a)') ' init  : restart (mode 2): the state file ', trim(cfg%init%restart_file), ' ', msg
            return
         end if
      end if

      if (.not. open_forcing(reg, cells, distance_km)) return
      if (.not. start_polygons(reg, places, ids, slot, inputs, have_maps)) return
      reg%step_days = max(1_ik, nint(cfg%dt_slow / day_sec, ik))
      allocate(reg%cost(places%ncell), source=0.0_wp)

      if (reg%verbose) then
         write(*,'(a)') '==================== ED-MORE run ===================='
         if (cfg%mode == RUN_MODE_REGION)                                                          &
            write(*,'(a,4(1x,f0.3),a,i0,a,f0.2,a)') ' region: box_nwse =', cfg%region%box_nwse, ', ',    &
               places%ncell, ' polygons (land fraction >= ', cfg%region%land_fraction_min, ')'
         write(*,'(5a)') ' run   : ', time_to_string(cfg%start_time), ' -> ', time_to_string(cfg%end_time)
         write(*,'(a,f0.2,a,i0)') ' span  : ', years_between(reg%now, cfg%end_time), ' yr  steps/yr~',    &
               max(1_ik, nint(yr_day / real(reg%step_days, wp), ik))
      end if
      call open_files(reg, ids, places)
      if (reg%verbose .and. places%ncell == 1_ik) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         call print_summary(reg%poly(1), 'start')
      end if
      end associate
      ok = .true.
   end subroutine open_run

   !----- The polygons, and the archive cells they read. A region's are its box's cells, each      !
   !      polygon at its cell's centre; a site run's is one polygon at the [site] location,        !
   !      reading the archive cell nearest it when the run reads the archive's forcing or maps     !
   !      (distance_km away), and none on a point forcing file. With [maps], each cell's map       !
   !      inputs. .false. (after saying why) when the cells or the maps cannot be had. ------------!
   logical function choose_polygons(cfg, verbose, cells, places, distance_km, inputs, have_maps) result(ok)
      type(ed_config_t),   intent(in)  :: cfg
      logical,             intent(in)  :: verbose
      type(met_cells_t),   intent(out) :: cells        !< the archive cells the run reads (none: ncell = 0)
      type(met_cells_t),   intent(out) :: places       !< each polygon's location and grid indices
      real(wp),            intent(out) :: distance_km  !< a site's distance to its cell
      type(map_inputs_t),  intent(out) :: inputs
      logical,             intent(out) :: have_maps
      type(met_grid_t)              :: grid
      character(len=:), allocatable :: msg
      integer(ik)                   :: st
      logical                       :: region, archive
      ok = .false.
      distance_km = 0.0_wp
      region    = cfg%mode == RUN_MODE_REGION
      have_maps = soil_maps_on(cfg%maps) .or. len_trim(cfg%maps%chave_e_map) > 0
      archive   = cfg%fast_biophysics_on .and. cfg%forcing%forcing_on .and.                          &
                  cfg%forcing%format == MET_BACKEND_ERA5LAND
      if (region .or. have_maps .or. archive) grid = era5land_grid(cfg%forcing)
      if (region) then
         if (.not. select_cells(cfg, grid, cells, verbose)) return
         places = cells
      else
         if (have_maps .or. archive) then
            call select_site_cell(grid, cfg%forcing%latitude, cfg%forcing%longitude,                 &
                                  cfg%forcing%max_distance_km, cells, distance_km, st)
            if (st == GRID_ERR_NO_CELL) then
               write(*,'(a,f0.1,a,f0.4,a,f0.4,a)') ' force : no valid ED_ERA5land cell within ',        &
                     cfg%forcing%max_distance_km, ' km of the site (', cfg%forcing%latitude, ', ',        &
                     cfg%forcing%longitude, ')'
               return
            else if (st /= GRID_OK) then
               write(*,'(2a)') ' force : cannot read the ED_ERA5land static file ', trim(grid%static_file)
               return
            end if
         end if
         call site_place(cfg, cells, places)
      end if
      if (have_maps) then
         call read_map_inputs(cfg%maps, grid, cells, inputs, msg)
         if (len(msg) > 0) then
            write(*,'(2a)') merge(' region: ', ' maps  : ', region), msg
            return
         end if
         if (verbose .and. .not. region) write(*,'(a)') ' maps  : soil and site index from the [maps] files'
      end if
      ok = .true.
   end function choose_polygons

   !----- A one-polygon run's place: the [site] location, which the polygon runs at, with the grid !
   !      indices of the archive cell it reads (0 when it reads none). ----------------------------!
   subroutine site_place(cfg, cells, place)
      type(ed_config_t),   intent(in)  :: cfg
      type(met_cells_t),   intent(in)  :: cells
      type(met_cells_t),   intent(out) :: place
      place%ncell = 1_ik
      allocate(place%row(1), place%col(1), source=0_ik)
      allocate(place%lat(1), source=cfg%forcing%latitude)
      allocate(place%lon(1), source=cfg%forcing%longitude)
      allocate(place%elevation(1), source=cfg%forcing%elevation)
      if (cells%ncell == 1_ik) then
         place%nlat = cells%nlat
         place%nlon = cells%nlon
         place%row  = cells%row
         place%col  = cells%col
      end if
   end subroutine site_place

   !----- The forcing: the archive at the run's cells, or a point forcing file; without forcing,   !
   !      the CO2 series alone, at whose CO2 the fast loop reads the reference climate             !
   !      (docs/science/forcing.md §12). .false. when it cannot be opened (met_open says why). ----!
   logical function open_forcing(reg, cells, distance_km) result(ok)
      type(ed_region_t),   intent(inout) :: reg
      type(met_cells_t),   intent(in)    :: cells
      real(wp),            intent(in)    :: distance_km
      integer(ik) :: st
      ok = .true.
      associate (cfg => reg%cfg)
      if (.not. cfg%fast_biophysics_on) return
      if (cfg%forcing%forcing_on) then
         if (cells%ncell > 0_ik) then
            call met_open(reg%met_src, cfg%forcing, stat=st, run_start=cfg%start_time,              &
                          run_end=cfg%end_time, cells=cells)
         else
            call met_open(reg%met_src, cfg%forcing, stat=st, run_start=cfg%start_time,              &
                          run_end=cfg%end_time)
         end if
         ok = st == MET_OK
         if (.not. ok) return
         if (cells%ncell > 0_ik) then
            if (cfg%mode == RUN_MODE_REGION) then
               write(*,'(a,i0,a)') ' force : ED_ERA5land region, ', cells%ncell, ' cells'
            else
               write(*,'(a,f0.3,a,f0.3,a,f7.2,a,f0.1,a)') ' force : ED_ERA5land cell (', cells%lat(1), ', ',  &
                     cells%lon(1), '), ', distance_km, ' km from the site; orography ', cells%elevation(1), ' m'
            end if
         end if
         if (reg%verbose) then
            if (cfg%forcing%format == MET_BACKEND_ERA5LAND) then
               write(*,'(3a)') ' force : met forcing ON (ED_ERA5land archive ', trim(cfg%forcing%data_path), ')'
            else
               write(*,'(3a)') ' force : met forcing ON (', trim(cfg%forcing%path), ')'
            end if
         end if
      else
         call met_open_co2(reg%met_src, cfg%forcing, cfg%start_time, cfg%end_time, stat=st)
         ok = st == MET_OK
      end if
      end associate
   end function open_forcing

   !----- Each polygon at its place, from its start, on its own configuration: the run's, with its !
   !      cell's map inputs. Then what follows the start (polygon_prepare); one polygon gives its  !
   !      patch loop the run's threads, a region gives them to its polygons. .false. when a        !
   !      start's file cannot serve it. -----------------------------------------------------------!
   logical function start_polygons(reg, places, ids, slot, inputs, have_maps) result(ok)
      type(ed_region_t),   intent(inout) :: reg
      type(met_cells_t),   intent(in)    :: places
      integer(ik),         intent(in)    :: ids(:), slot(:)
      type(map_inputs_t),  intent(in)    :: inputs
      logical,             intent(in)    :: have_maps
      type(ed_config_t)                :: pcfg
      type(ed_config_t),   allocatable :: site_cfg(:)     !< [[sites]]: each site's configuration ...
      real(wp),            allocatable :: site_area(:)    !< ... and share of the polygon
      integer(ik) :: p, n
      ok = .false.
      n = places%ncell
      associate (cfg => reg%cfg)
      allocate(reg%poly(n), reg%run(n), reg%out_bufs(n))
      reg%now = cfg%start_time
      do p = 1_ik, n
         associate (poly => reg%poly(p), run => reg%run(p))
            pcfg = cfg
            if (cfg%mode == RUN_MODE_REGION)                                                        &
               write(run%label,'(a,i0,a,f0.2,a,f0.2,a)') 'polygon ', ids(p), ' (', places%lat(p), ', ',  &
                                                       places%lon(p), ')'
            if (have_maps) then
               call apply_map_inputs(inputs, p, pcfg)
               call validate_config(pcfg)
            end if
            if (allocated(site_cfg)) deallocate(site_cfg, site_area)
            if (allocated(pcfg%sites)) call site_configs(pcfg, site_cfg, site_area)
            call init_polygon(poly, pcfg, polygon_name(run), reg%verbose .and. n == 1_ik, slot(p), reg%now,     &
                              ok, ids(p), places%lat(p), places%lon(p), places%elevation(p),                &
                              site_cfg=site_cfg, site_area=site_area)
            if (.not. ok) return
            call polygon_prepare(cfg, reg%met_src, p, poly, run,                                       &
                                 patch_threads=merge(cfg%n_threads, 1_ik, n == 1_ik),               &
                                 verbose=reg%verbose .and. n == 1_ik)
         end associate
      end do
      !----- A restart's polygons all read the restructuring their file owes; a new stand owes     !
      !      none. --------------------------------------------------------------------------------!
      reg%restructure_pending  = reg%poly(1)%restructure_pending
      reg%restructure_new_year = reg%poly(1)%restructure_new_year
      end associate
   end function start_polygons

   !----- The run's files: the state files' directory and the run's PFT table; the diagnostic      !
   !      output ([output].enabled), behind the polygon axis, each polygon writing its own layers, !
   !      cohorts and patches; and the parameter record beside the output and the state, every     !
   !      setting the loader read, from which file, set or defaulted, and its value, so a caller   !
   !      can tell whether a key was read. --------------------------------------------------------!
   subroutine open_files(reg, ids, places)
      type(ed_region_t),   intent(inout) :: reg
      integer(ik),         intent(in)    :: ids(:)
      type(met_cells_t),   intent(in)    :: places
      integer(ik) :: p
      associate (cfg => reg%cfg, out => reg%cfg%output)
      if (out%state_write) then
         call ensure_output_dir(trim(out%state_file_dir))
         call write_pft_params_csv(cfg, trim(out%state_file_dir)//'/'//trim(out%state_file_prefix)//     &
                                   '_pft_parameters.csv')
      end if
      if (out%enabled) then
         call ensure_output_dir(trim(out%dir))
         call open_output_files(cfg, reg%out_files, reg%verbose, ids, reg%poly(:)%site%n, places%lat,    &
                                places%lon, places%row, places%col)
         do p = 1_ik, size(reg%poly, kind=ik)
            call attach_output(cfg, reg%out_files, reg%out_bufs(p), reg%poly(p), p)
         end do
         if (reg%verbose) write(*,'(a)') ' output: diagnostic aggregation ON ([output])'
      end if
      if (out%enabled) call toml_write_record(trim(out%dir)//'/'//trim(out%prefix)//'_parameters.csv')
      if (out%state_write)                                                                         &
         call toml_write_record(trim(out%state_file_dir)//'/'//trim(out%state_file_prefix)//'_parameters.csv')
      end associate
   end subroutine open_files

   !----- Create an output directory that does not exist yet. The model's libraries never touch    !
   !      the file system; the run does. ----------------------------------------------------------!
   subroutine ensure_output_dir(dir)
      character(len=*), intent(in) :: dir
      integer :: stat
      if (len_trim(dir) == 0 .or. trim(dir) == '.') return
      call execute_command_line('mkdir -p "'//trim(dir)//'"', wait=.true., exitstat=stat)
      if (stat /= 0) write(*,'(3a)') ' warning: could not create output dir "', trim(dir), '"'
   end subroutine ensure_output_dir

   !----- A region's cells: the box's valid cells with enough land, in row-major order, less those   !
   !      the cell mask drops. .false. (after saying why) when none is left or a file cannot be read. !
   logical function select_cells(cfg, grid, cells, verbose) result(ok)
      type(ed_config_t),   intent(in)  :: cfg
      type(met_grid_t),    intent(in)  :: grid
      type(met_cells_t),   intent(out) :: cells
      logical,             intent(in)  :: verbose
      integer(ik) :: st, n, j
      ok = .false.
      call select_box_cells(grid, cfg%region%box_nwse, cells, st, cfg%region%land_fraction_min)
      if (st == GRID_ERR_NO_CELL) then
         write(*,'(a,4(f0.3,1x),a,f0.2)') ' region: no valid cell in box_nwse = ', cfg%region%box_nwse, &
               'with land fraction >= ', cfg%region%land_fraction_min
         return
      else if (st /= GRID_OK) then
         write(*,'(2a)') ' region: cannot read valid, elevation and land_fraction from ', trim(grid%static_file)
         return
      end if
      if (len_trim(cfg%region%cell_mask) > 0) then
         n = cells%ncell
         call apply_cell_mask(cells, trim(cfg%region%cell_mask), trim(cfg%region%cell_mask_var),       &
                              cfg%region%cell_mask_min, grid, j, st)
         if (st /= MAP_OK) then
            write(*,'(4a)') ' region: the cell mask ', trim(cfg%region%cell_mask), ' ', map_status_message(st)
            return
         end if
         if (verbose) write(*,'(a,i0,a,i0,3a)') ' region: the cell mask keeps ', cells%ncell, ' of ', n, &
               ' cells (', trim(cfg%region%cell_mask_var), ')'
         if (cells%ncell == 0_ik) then
            write(*,'(2a)') ' region: the cell mask keeps no cell of the box: ', trim(cfg%region%cell_mask)
            return
         end if
      end if
      ok = .true.
   end function select_cells

   !----- How messages name a polygon: its label, or "the site" for a one-polygon run. ------------!
   function polygon_name(run) result(name)
      type(polygon_run_t), intent(in) :: run
      character(len=:), allocatable :: name
      if (len_trim(run%label) > 0) then
         name = trim(run%label)
      else
         name = 'the site'
      end if
   end function polygon_name

   !----- Is the environment variable set? --------------------------------------------------------!
   logical function env_set(name) result(yes)
      character(len=*), intent(in) :: name
      integer :: n, st
      call get_environment_variable(name, length=n, status=st)
      yes = st == 0 .and. n > 0
   end function env_set

   !---------------------------------------------------------------------------------------!
   ! region_step -- one slow step of every polygon, as the C API steps a run: the step's   !
   ! forcing (loaded when a month starts), the polygons, and after a month's last step the !
   ! I/O phase.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine region_step(reg, status)
      type(ed_region_t),   intent(inout) :: reg
      integer(ik),         intent(out)   :: status
      if (region_done(reg)) then
         status = RUN_FINISHED
         return
      end if
      call step_polygons(reg, list_steps(reg, 1_ik), status)
   end subroutine region_step

   !---------------------------------------------------------------------------------------!
   ! region_step_month -- every polygon through a month: to the next month boundary, to    !
   ! end_time, or to a state file that falls inside the month.                             !
   !---------------------------------------------------------------------------------------!
   subroutine region_step_month(reg, status)
      type(ed_region_t),   intent(inout) :: reg
      integer(ik),         intent(out)   :: status
      if (region_done(reg)) then
         status = RUN_FINISHED
         return
      end if
      call step_polygons(reg, list_steps(reg, huge(1_ik)), status)
   end subroutine region_step_month

   !----- The steps from the run's date, at most max_steps of them, ending at the first month      !
   !      boundary, at end_time, or where a state file every N days falls due. A step's forcing    !
   !      lies in the month its start date belongs to (daily steps from midnight, and              !
   !      validate_config puts the recycle seam on a month boundary), so one prefetch at the first !
   !      step's start loads every step's forcing. ------------------------------------------------!
   function list_steps(reg, max_steps) result(steps)
      type(ed_region_t),   intent(in) :: reg
      integer(ik),         intent(in) :: max_steps
      type(calendar_step_t), allocatable :: steps(:)
      type(calendar_step_t) :: step
      integer(ik)           :: n
      allocate(steps(0))
      step%end = reg%now
      do
         step  = calendar_step(step%end, reg%step_days)
         steps = [steps, step]
         n = size(steps, kind=ik)
         if (n >= max_steps .or. step%new_month) exit
         if (.not. time_lt(step%end, reg%cfg%end_time)) exit
         if (state_file_due(reg%cfg%output, reg%istep + n, step, reg%imonth, reg%iyear)) exit
      end do
   end function list_steps

   !---------------------------------------------------------------------------------------!
   ! step_polygons -- the listed steps of every polygon, then the I/O phase. The polygons  !
   ! run side by side on [run].n_threads threads, costliest first; they share only what    !
   ! they read, so the results do not depend on the thread count. One polygon runs on the  !
   ! calling thread and gives its patch loop the threads. A polygon that fails stops; the  !
   ! others finish, the output is written, and the status returned is the first failure's. !
   !---------------------------------------------------------------------------------------!
   subroutine step_polygons(reg, steps, status)
      type(ed_region_t),     intent(inout) :: reg
      type(calendar_step_t), intent(in)    :: steps(:)
      integer(ik),           intent(out)   :: status
      logical,     allocatable :: was_running(:)
      integer(ik), allocatable :: order(:), failed_step(:)
      integer(ik) :: nstep, s, p, k, st, npoly
      integer(ik) :: tick0, tick1, tick_rate

      status = RUN_OK
      associate (cfg => reg%cfg, last => steps(size(steps)))
      npoly = size(reg%poly, kind=ik)
      nstep = size(steps, kind=ik)
      if (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on) then
         call system_clock(tick0, tick_rate)
         call met_prefetch(reg%met_src, steps(1)%start)
         call system_clock(tick1)
         reg%year_read_seconds = reg%year_read_seconds + real(tick1 - tick0, wp) / real(max(tick_rate, 1_ik), wp)
      end if

      !----- The compute phase: no file is touched here; each polygon writes only its own state    !
      !      and buffers. -------------------------------------------------------------------------!
      allocate(order(npoly), failed_step(npoly), was_running(npoly))
      call ascending_order(-reg%cost, npoly, order)
      was_running = reg%run(:)%status == RUN_OK
      failed_step = 0_ik
      !$omp parallel do schedule(dynamic, 1) num_threads(cfg%n_threads) default(shared) if(npoly > 1)   &
      !$omp    private(k, p, s, st, tick0, tick1, tick_rate)
      do k = 1_ik, npoly
         p = order(k)
         if (.not. was_running(p)) cycle
         call system_clock(tick0, tick_rate)
         do s = 1_ik, nstep
            call polygon_step(cfg, reg%met_src, reg%out_files, reg%out_bufs(p), reg%poly(p), reg%run(p),  &
                              steps(s), st)
            if (st /= RUN_OK) then
               failed_step(p) = s
               exit
            end if
         end do
         call system_clock(tick1)
         reg%cost(p) = real(tick1 - tick0, wp) / real(max(tick_rate, 1_ik), wp)
      end do
      !$omp end parallel do

      do p = 1_ik, npoly
         if (failed_step(p) == 0_ik) cycle
         if (npoly > 1_ik) write(*,'(4a)') ' region: ', trim(reg%run(p)%label), ' failed on ',       &
                                           time_to_string(steps(failed_step(p))%end)
         if (status == RUN_OK) status = reg%run(p)%status
      end do
      reg%now   = last%end
      reg%istep = reg%istep + nstep
      if (last%new_month) reg%imonth = reg%imonth + 1_ik
      reg%restructure_pending  = last%new_month
      reg%restructure_new_year = last%new_year

      !----- A year's end: the summary line(s), and the forcing-read accounting. ------------------!
      if (last%new_year) then
         reg%iyear = reg%iyear + 1_ik
         if (reg%verbose) call year_summary(reg)
         reg%read_seconds = reg%read_seconds + reg%year_read_seconds
         reg%year_read_seconds = 0.0_wp
         reg%year_loads0 = reg%met_src%n_loads
      end if

      !----- The I/O phase: the records the steps closed, after a month's last step, a failure or  !
      !      the run's last step, and the state file when one is due. -----------------------------!
      if (cfg%output%enabled .and. (last%new_month .or. status /= RUN_OK .or. region_done(reg)))  &
         call output_serialize(reg%out_files, reg%out_bufs)
      if (state_file_due(cfg%output, reg%istep, last, reg%imonth, reg%iyear)) call write_state_file(reg)
      end associate
   end subroutine step_polygons

   !----- Is a state file due after `step`, the step that brings the run to `istep` steps, the run !
   !      having crossed imonth months and iyear years by then? Every state_every days from the    !
   !      start, or at every state_every-th month or year boundary ([output].state_interval). -----!
   pure logical function state_file_due(out, istep, step, imonth, iyear) result(due)
      type(output_params_t), intent(in) :: out
      integer(ik),           intent(in) :: istep, imonth, iyear
      type(calendar_step_t), intent(in) :: step
      due = .false.
      if (.not. out%state_write) return
      select case (out%state_unit)
      case (STATE_DAYS)
         due = mod(istep * step%days, out%state_every) == 0_ik
      case (STATE_MONTHS)
         due = step%new_month .and. mod(imonth, out%state_every) == 0_ik
      case (STATE_YEARS)
         due = step%new_year .and. mod(iyear, out%state_every) == 0_ik
      end select
   end function state_file_due

   !----- The state file: every polygon's row, at the current date (ed_init_restart). The          !
   !      cohort, patch and site axes are the largest polygon's. The restructuring it owes is the  !
   !      calendar's (the run's), which a polygon that has failed no longer tracks. ---------------!
   subroutine write_state_file(reg)
      type(ed_region_t), intent(in) :: reg
      type(restart_file_t)       :: sf
      integer(ik), allocatable :: ids(:)
      real(wp),    allocatable :: lat(:), lon(:)
      integer(ik)              :: p, n_cohort, n_patch, n_site
      n_cohort = 0_ik
      n_patch = 0_ik
      n_site = 0_ik
      allocate(ids(size(reg%poly)), lat(size(reg%poly)), lon(size(reg%poly)))
      do p = 1_ik, size(reg%poly, kind=ik)
         ids(p)   = reg%poly(p)%id
         lat(p)   = reg%poly(p)%latitude_deg
         lon(p)   = reg%poly(p)%longitude_deg
         n_cohort = max(n_cohort, reg%poly(p)%cohort%n)
         n_patch  = max(n_patch, reg%poly(p)%patch%n)
         n_site   = max(n_site, reg%poly(p)%site%n)
      end do
      call restart_file_open(sf, trim(reg%cfg%output%state_file_dir), trim(reg%cfg%output%state_file_prefix), reg%now, &
                           reg%cfg%pft%n, int(size(reg%poly(1)%cohort%growth_hist, 1), ik), n_cohort,       &
                           n_patch, n_site, polygon_id=ids, polygon_lat=lat, polygon_lon=lon,              &
                           restructure_pending=reg%restructure_pending,                                    &
                           restructure_new_year=reg%restructure_new_year)
      do p = 1_ik, size(reg%poly, kind=ik)
         call restart_file_write(sf, reg%poly(p), reg%now, p)
      end do
      call restart_file_close(sf)
   end subroutine write_state_file

   !----- A year's summary: one polygon's stand every 5 years, or a line for the polygons' mean    !
   !      and range of AGB and LAI every year, with what the forcing reads cost. ------------------!
   subroutine year_summary(reg)
      type(ed_region_t),   intent(in) :: reg
      character(len=19) :: when
      real(wp)    :: agb(size(reg%poly)), lai(size(reg%poly))
      integer(ik) :: p
      when = time_to_string(reg%now)
      if (size(reg%poly) == 1) then
         if (mod(reg%iyear, 5_ik) == 0_ik .or. reg%iyear == 1_ik) call print_summary(reg%poly(1), when(1:10))
         return
      end if
      do p = 1_ik, size(reg%poly, kind=ik)
         agb(p) = total_agb(reg%poly(p))
         lai(p) = total_lai(reg%poly(p))
      end do
      write(*,'(a,a10,a,3(1x,f0.3),a,3(1x,f0.3))') ' ', when(1:10), '  AGB kgC/m2 min/mean/max:',     &
            minval(agb), sum(agb) / size(agb), maxval(agb), '   LAI:', minval(lai),                  &
            sum(lai) / size(lai), maxval(lai)
      write(*,'(a,f0.2,a,i0,a)') '             forcing: ', reg%year_read_seconds, ' s loading, ',        &
            reg%met_src%n_loads - reg%year_loads0, ' archive month(s) read from file'
   end subroutine year_summary

   !---------------------------------------------------------------------------------------!
   ! region_finalize -- the terminal state file, the end-of-run report, and every file and !
   ! the forcing closed; RUN_ERR_AREA if a polygon's patch areas do not sum to 1. After !
   ! a failure the run writes nothing more, so a closed chunk file is never made again     !
   ! from scratch.                                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine region_finalize(reg, status)
      type(ed_region_t),     intent(inout) :: reg
      integer(ik), optional, intent(out)   :: status
      integer(ik) :: p, st, n_bad
      real(wp)    :: a1
      logical     :: failed

      st = RUN_OK
      n_bad = 0_ik
      failed = any(reg%run(:)%status /= RUN_OK)
      !----- The true terminal state, so a restart resumes exactly here. --------------------------!
      if (reg%cfg%output%state_write .and. .not. failed) call write_state_file(reg)
      if (size(reg%poly) == 1) then
         associate (poly => reg%poly(1), run => reg%run(1))
            if (reg%verbose) call print_summary(poly, 'final')
            a1 = total_area(poly)
            if (reg%verbose) then
               write(*,'(a)') '-----------------------------------------------------------------------------'
               write(*,'(a,f12.9,a,f12.9)') ' site area start=', run%area_start, '  end=', a1
            end if
            if (abs(a1 - 1.0_wp) > AREA_TOL) st = RUN_ERR_AREA
            if (reg%verbose) call polygon_report(reg%cfg, run)
         end associate
      else
         if (reg%verbose) then
            write(*,'(a)') '-----------------------------------------------------------------------------'
            write(*,'(a)') '  polygon      lat      lon  AGB kgC/m2    LAI  cohorts  area error  energy/water fails'
         end if
         do p = 1_ik, size(reg%poly, kind=ik)
            associate (poly => reg%poly(p), run => reg%run(p))
               a1 = total_area(poly)
               if (abs(a1 - 1.0_wp) > AREA_TOL) then
                  st = RUN_ERR_AREA
                  n_bad = n_bad + 1_ik
               end if
               if (reg%verbose)                                                                    &
                  write(*,'(i9,2f9.2,f12.4,f7.3,i9,es12.2,2(1x,i0))') poly%id,                      &
                        poly%latitude_deg, poly%longitude_deg, total_agb(poly),                     &
                        total_lai(poly), poly%cohort%n, a1 - 1.0_wp,                          &
                        run%ledgers%energy%n_fail, run%ledgers%water%n_fail
            end associate
         end do
         if (n_bad > 0_ik) write(*,'(a,i0,a)') ' ERROR: ', n_bad, ' polygon(s) did not conserve area'
         if (reg%verbose) write(*,'(a,f0.2,a,i0,a)') ' forcing: ',                                  &
               reg%read_seconds + reg%year_read_seconds, ' s loading in all, ', reg%met_src%n_loads,    &
               ' archive month(s) read from file'
      end if

      if (reg%cfg%output%enabled) call output_close(reg%out_files, reg%out_bufs, flush_partial=.not. failed)
      if (reg%cfg%fast_biophysics_on) call met_close(reg%met_src)
      if (present(status)) status = st
   end subroutine region_finalize

   !---------------------------------------------------------------------------------------!
   ! region_free -- close a run's files (an open period is not written) and its forcing,   !
   ! free its polygons, and leave the variable as a run never opened. Apart from finalize, !
   ! so a caller can read the final state after the files are closed.                      !
   !---------------------------------------------------------------------------------------!
   subroutine region_free(reg)
      type(ed_region_t), intent(inout) :: reg
      type(ed_region_t) :: never_opened     ! default-initialised, not built by a constructor (nvfortran)
      integer(ik) :: p
      if (allocated(reg%out_bufs)) call output_close(reg%out_files, reg%out_bufs, flush_partial=.false.)
      call met_close(reg%met_src)
      if (allocated(reg%poly)) then
         do p = 1_ik, size(reg%poly, kind=ik)
            call free_polygon(reg%poly(p))
         end do
      end if
      reg = never_opened
   end subroutine region_free

end module ed_region
