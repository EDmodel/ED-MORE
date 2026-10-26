! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_met_types -- the runtime meteorological-forcing types (docs/science/forcing.md §2):   !
! met_forcing_t, the instantaneous atmospheric state the fast loop reads; met_record_t,    !
! one ingested file record; met_source_t, the reader state a run shares, and               !
! met_cursor_t, each polygon's own; met_grid_t and met_cells_t, a gridded product's grid   !
! and the cells a run reads; met_month_t, the archive's month buffer; co2_series_t, the    !
! prescribed CO2 series. Uses only ed_kinds, ed_time and ed_forcing_params; the            !
! reader (ed_met) owns the only mutable forcing state.                                     !
!==========================================================================================!
module ed_met_types
   use ed_kinds,            only : wp, sp, ik
   use ed_time,             only : ed_time_t
   use ed_forcing_params,   only : forcing_params_t, MET_BACKEND_DEFAULT, MET_PATH_LEN
   implicit none
   private

   public :: met_forcing_t, met_record_t, met_source_t, met_cursor_t, met_grid_t, met_cells_t, met_month_t
   public :: co2_series_t
   public :: HUMIDITY_QAIR, HUMIDITY_RHAIR, HUMIDITY_TDEW
   public :: N_ED_FIELD, ED_FIELD, FLD_TAIR, FLD_QAIR, FLD_RHAIR, FLD_TDEW, FLD_PSURF, FLD_WIND
   public :: FLD_U10, FLD_V10, FLD_RAINF, FLD_LWDOWN, FLD_SWDOWN, FLD_PAR_BEAM, FLD_PAR_DIFFUSE
   public :: FLD_NIR_BEAM, FLD_NIR_DIFFUSE

   !----- The fields an ED-MORE forcing file may carry, listed once. met_open finds each one's      !
   !      column in the series it reads (met_source_t%series_col); read_record looks them up by      !
   !      index. ------------------------------------------------------------------------------------!
   integer(ik), parameter :: N_ED_FIELD = 15_ik
   integer(ik), parameter :: FLD_TAIR = 1_ik, FLD_QAIR = 2_ik, FLD_RHAIR = 3_ik, FLD_TDEW = 4_ik,     &
                             FLD_PSURF = 5_ik, FLD_WIND = 6_ik, FLD_U10 = 7_ik, FLD_V10 = 8_ik,        &
                             FLD_RAINF = 9_ik, FLD_LWDOWN = 10_ik, FLD_SWDOWN = 11_ik,                 &
                             FLD_PAR_BEAM = 12_ik, FLD_PAR_DIFFUSE = 13_ik, FLD_NIR_BEAM = 14_ik,      &
                             FLD_NIR_DIFFUSE = 15_ik
   character(len=24), parameter :: ED_FIELD(N_ED_FIELD) = [character(len=24) ::                   &
      'Tair', 'Qair', 'RHair', 'Tdew', 'PSurf', 'Wind', 'u10', 'v10', 'Rainf', 'LWdown', 'SWdown',     &
      'SWdown_par_beam', 'SWdown_par_diffuse', 'SWdown_nir_beam', 'SWdown_nir_diffuse']

   !----- The humidity a source carries, as its source measured it, converted with the model's own !
   !      saturation curve (docs/science/forcing.md §7). Qair is for a model-made source. ----------!
   integer(ik), parameter :: HUMIDITY_QAIR  = 0_ik          !< specific humidity [kg/kg], used as is
   integer(ik), parameter :: HUMIDITY_RHAIR = 1_ik          !< relative humidity over liquid water [1]
   integer(ik), parameter :: HUMIDITY_TDEW  = 2_ik          !< dewpoint temperature [K]

   !==========================================================================================!
   !  met_forcing_t -- the instantaneous atmospheric state the fast loop reads: a polygon's       !
   !  sample, or a site's or a patch's after the vertical moves (docs/science/forcing.md §8).     !
   !  Plain scalars, so it copies cheaply. The defaults are a reference climate, the air a run   !
   !  without forcing reads at every sub-step (the four streams sum to 400 W/m2).                 !
   !==========================================================================================!
   type :: met_forcing_t
      real(wp) :: air_temp     = 288.0_wp     !< [K]        air temperature at reference height
      real(wp) :: air_shv      = 0.008_wp     !< [kg/kg]    specific humidity
      real(wp) :: air_pressure = 101325.0_wp  !< [Pa]       surface pressure
      !----- Total precipitation, split into rain and snow per patch (docs/science/forcing.md §7). ---!
      real(wp) :: precip       = 0.0_wp       !< [kg/m2/s]  total precipitation rate
      real(wp) :: wind         = 2.0_wp       !< [m/s]      wind speed at reference height
      !----- The wind vector, when the source has components (docs/science/forcing.md §1). ----------!
      real(wp) :: wind_u       = 2.0_wp       !< [m/s]      eastward component at reference height
      real(wp) :: wind_v       = 0.0_wp       !< [m/s]      northward component at reference height
      logical  :: has_wind_vector = .false.   !< .true. when the source supplied components, not just speed
      real(wp) :: lw_down      = 380.0_wp     !< [W/m2]     downwelling longwave (positive down)
      real(wp) :: par_beam     = 180.0_wp     !< [W/m2]     direct-beam PAR at canopy top
      real(wp) :: par_diffuse  = 40.0_wp      !< [W/m2]     diffuse PAR
      real(wp) :: nir_beam     = 150.0_wp     !< [W/m2]     direct-beam NIR
      real(wp) :: nir_diffuse  = 30.0_wp      !< [W/m2]     diffuse NIR   (Sigma = 400 W/m2)
      real(wp) :: co2          = 420.0_wp     !< [umol/mol] free-atmosphere CO2 ([forcing].co2_source)
      real(wp) :: cosz         = 0.5_wp       !< [-]        cos(solar zenith); derived each sub-step
      real(wp) :: air_density  = 1.2_wp       !< [kg/m3]    derived from air_temp, air_pressure and air_shv
   contains
      procedure :: sw_down        => met_sw_down         !< total downward SW = sum of the four streams
   end type met_forcing_t

   !----- One record as ingested (read_record), before interpolation: converted, lapsed when the run !
   !      lapses, its shortwave in four streams. No CO2 (docs/science/forcing.md §12). ------------!
   type :: met_record_t
      type(ed_time_t) :: when
      real(wp) :: air_temp = 288.0_wp, air_shv = 0.008_wp, air_pressure = 101325.0_wp
      real(wp) :: rain = 0.0_wp, wind = 2.0_wp, lw_down = 380.0_wp
      real(wp) :: wind_u = 2.0_wp, wind_v = 0.0_wp          !< wind vector (only when the source has it)
      real(wp) :: par_beam = 180.0_wp, par_diffuse = 40.0_wp, nir_beam = 150.0_wp, nir_diffuse = 30.0_wp
   end type met_record_t

   !==========================================================================================!
   !  met_grid_t -- a gridded product's grid (ed_met_grid): its static file, which holds the       !
   !  regular lat/lon axes and each cell's valid mask, orography and land fraction under these     !
   !  names, and the edge of the square chunks its data files are stored in.                      !
   !==========================================================================================!
   type :: met_grid_t
      character(len=MET_PATH_LEN) :: static_file = ''
      character(len=32) :: lat_var = 'lat', lon_var = 'lon'
      character(len=32) :: valid_var = 'valid', elevation_var = 'elevation', land_fraction_var = 'land_fraction'
      integer(ik)       :: chunk = 1_ik                     !< [cells] the data files' chunk edge
   end type met_grid_t

   !==========================================================================================!
   !  met_cells_t -- the cells of a grid a run reads: one for a site, a box's valid cells for a   !
   !  region. Grid indices are 0-based, as the netCDF C API counts. The cells are grouped by the  !
   !  grid's chunk, so a month load reads each touched chunk column once.                        !
   !==========================================================================================!
   type :: met_cells_t
      integer(ik) :: ncell = 0_ik
      integer(ik) :: nlat = 0_ik, nlon = 0_ik               !< grid size
      integer(ik) :: chunk = 1_ik                           !< [cells] the chunk edge they are grouped by
      integer(ik), allocatable :: row(:), col(:)            !< 0-based (lat, lon) index of each cell
      real(wp),    allocatable :: lat(:), lon(:)            !< [deg] cell centre
      real(wp),    allocatable :: elevation(:)              !< [m] static-file orography
      real(wp),    allocatable :: land_fraction(:)          !< [-] static-file land fraction (box selection only)
      integer(ik) :: nchunk = 0_ik
      integer(ik), allocatable :: chunk_row(:), chunk_col(:)!< 0-based first row/col of the box read in each chunk
      integer(ik), allocatable :: chunk_nrow(:), chunk_ncol(:) !< that box's extent (1 x 1 for a site)
      integer(ik), allocatable :: chunk_first(:)            !< by_chunk(chunk_first(k):chunk_first(k+1)-1)
      integer(ik), allocatable :: by_chunk(:)               !< cell numbers grouped by chunk
   end type met_cells_t

   !----- One archive month for every cell: values(hour, cell, variable), float32 as stored, from  !
   !      01:00 on the 1st to 00:00 on the 1st of the next month (end-stamped). -------------------!
   type :: met_month_t
      integer(ik) :: year = 0_ik, month = 0_ik, nt = 0_ik
      real(sp), allocatable :: values(:,:,:)
   end type met_month_t

   !==========================================================================================!
   !  co2_series_t -- an ED-MORE CO2 file as read (co2_series_read): each period's mean, placed at  !
   !  the period's middle (docs/science/forcing.md eq. (18)). Times are seconds after the first    !
   !  period's start.                                                                              !
   !==========================================================================================!
   type :: co2_series_t
      integer(ik)       :: n = 0_ik                         !< rows (periods)
      type(ed_time_t)   :: first_start                      !< start of the first period: the time origin
      real(wp)          :: span_end_sec = 0.0_wp            !< [s] end of the last period, after first_start
      real(wp), allocatable :: mid_sec(:)                   !< [s] each period's middle, after first_start
      real(wp), allocatable :: co2(:)                       !< [umol/mol] each period's mean
   end type co2_series_t

   !==========================================================================================!
   !  met_source_t -- the reader state a run shares, read by every polygon (docs/science/forcing.md  !
   !  §2): the config, the time axis and recycle window, the cells, and the records in memory. Only  !
   !  met_open, met_prefetch and met_close change it, all outside the compute phase, so polygons     !
   !  can share it across threads.                                                                    !
   !==========================================================================================!
   type :: met_source_t
      type(forcing_params_t)  :: fcfg                       !< [forcing]/[site] config
      integer(ik) :: backend    = MET_BACKEND_DEFAULT    !< MET_BACKEND_*: ED_DEFAULT | ED_ERA5LAND
      integer(ik) :: humidity   = HUMIDITY_TDEW             !< the humidity the source carries (HUMIDITY_*)
      integer(ik) :: ncid       = -1_ik                     !< NetCDF handle (via ed_netcdf_c), -1 if closed
      integer(ik) :: grid_index = 1_ik                      !< ED-MORE forcing file: the location (1..ngrid) read
      integer(ik) :: ngrid      = 1_ik                      !< total locations in the file (cells for the archive)
      integer(ik) :: nrec       = 0_ik                      !< total time records on the axis
      integer(ik) :: rec_first  = 1_ik                      !< the record a cursor's first bracket starts at
      real(wp)    :: dt_forcing  = 3600.0_wp                !< [s] native forcing interval
      type(ed_time_t)  :: base_time                         !< time-coordinate anchor ("seconds since base_time")
      real(wp), allocatable :: time_sec(:)                  !< [s] cached time coordinate (seconds since base_time)
      !----- The declared recycle window, checked against the file at open (docs/science/forcing.md §9). !
      type(ed_time_t) :: cycle_anchor                       !< = fcfg%recycle_start; the cycle's first instant
      integer(ik) :: n_cycle_years   = 0_ik                 !< whole calendar years the declared window spans
      integer(ik) :: irec_cycle_first = 0_ik                !< record index of cycle_anchor (exact match, 1-based)
      integer(ik) :: irec_cycle_last  = 0_ik                !< record index of the last record strictly inside the window
      logical     :: has_wind_vector = .false.              !< the source supplies u10/v10
      !----- The archive: its months laid end to end on time_sec (seconds since 1970-01-01); month k !
      !      holds records month_rec0(k)+1 .. month_rec0(k)+hours, and one month sits in `buffer`. ---!
      character(len=MET_PATH_LEN) :: file_template = ''     !< resolved file template
      type(met_grid_t)  :: grid                             !< the archive's grid
      type(met_cells_t) :: cells                            !< the cells read: one for a site, many for a region
      integer(ik), allocatable :: month_year(:), month_month(:), month_rec0(:)
      type(met_month_t) :: buffer                           !< the loaded month (year = 0 before the first)
      !----- What a step reads beyond the loaded month, in memory before the step                  !
      !      (docs/science/forcing.md §2 and §9). -----------------------------------------------!
      real(sp),    allocatable :: carry(:,:)                !< (cell, variable) the record before the loaded month
      integer(ik) :: carry_rec = 0_ik                       !< its record index on the axis (0 = none)
      !----- The recycle window's first day, which a step across the seam reads after the window's !
      !      last record, a month the buffer does not hold. -------------------------------------!
      real(sp),    allocatable :: head(:,:,:)               !< (record, cell, variable) from irec_cycle_first on
      integer(ik) :: n_head    = 0_ik                       !< records in head (0 = none)
      integer(ik) :: n_loads   = 0_ik                       !< archive month loads so far (tests: none inside a step)
      real(wp),    allocatable :: series(:,:)               !< (record, column) ED-MORE forcing file at grid_index
      integer(ik) :: series_col(N_ED_FIELD) = 0_ik          !< each ED_FIELD's column in series (0 = absent)
      !----- The prescribed CO2 (co2_source = "file"): one global series, read at open. ----------!
      type(co2_series_t) :: co2
   end type met_source_t

   !==========================================================================================!
   !  met_cursor_t -- one polygon's reader state: its cell and location, and the two records        !
   !  bracketing the model time. Plain scalars and records, so a region holds one per polygon       !
   !  cheaply. The polygon's location ([site] in a site run, its cell centre in a region) drives    !
   !  the solar geometry and the elevation lapse.                                                   !
   !==========================================================================================!
   type :: met_cursor_t
      integer(ik) :: cell = 1_ik                            !< this polygon's cell in the source (1 for a site run)
      real(wp)    :: latitude_deg     = 0.0_wp              !< [deg +N] solar geometry
      real(wp)    :: longitude_deg    = 0.0_wp              !< [deg +E] local solar time from the UTC clock
      real(wp)    :: elevation_m      = 0.0_wp              !< [m] the polygon's elevation (lapse target)
      real(wp)    :: grid_elevation_m = 0.0_wp              !< [m] forcing-cell elevation (lapse origin)
      integer(ik) :: irec_prev  = 0_ik                      !< cursor: index of the "previous" bracketing record
      type(met_record_t) :: rec_prev, rec_next              !< the two records bracketing the model time
      logical     :: at_wrap_seam    = .false.              !< current bracket is the cycle-boundary seam (last -> first)
      !----- The longwave synthesis holds the last daytime clearness index through the night, from   !
      !      clear sky at the start (docs/science/forcing.md §11). --------------------------------!
      real(wp)    :: kt_last_day    = 1.0_wp                !< [-] last daytime clearness index
   end type met_cursor_t

contains

   pure real(wp) function met_sw_down(met) result(sw)
      class(met_forcing_t), intent(in) :: met
      sw = met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse
   end function met_sw_down


end module ed_met_types
