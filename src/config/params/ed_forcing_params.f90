! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_forcing_params -- the [forcing] and [site] blocks (docs/science/forcing.md) and the
! selector codes the config and the reader share; the shortwave partition's and the
! clear-sky emissivity's are ed_radiation_lib's. A config leaf (share only), so
! ed_config carries forcing_params_t with no edge back to the forcing layer. The
! defaults are an Ithaca NY / ERA5-Land reference site.
!---~---
module ed_forcing_params
   use ed_kinds, only : wp, ik
   use ed_radiation_lib, only : SWPART_CLEARIDX, LW_CLEAR_BRUTSAERT
   use ed_time,  only : ed_time_t
   implicit none
   private

   public :: forcing_params_t
   public :: MET_BACKEND_DEFAULT, MET_BACKEND_ERA5LAND
   public :: MET_PATH_LEN, ARCHIVE_DT_SEC, MAX_RECYCLE_YEARS
   public :: METAVG_END, METAVG_BEGIN
   public :: LW_FILE, LW_SYNTHESIZE
   public :: CLAMP_ERROR, CLAMP_HOLD
   public :: GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   public :: CO2_SOURCE_CONST, CO2_SOURCE_FILE
   public :: HEIGHT_ABOVE_ZERO_PLANE, HEIGHT_ABOVE_GROUND, WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL

   !--- The ED_ERA5land archive's record spacing: hourly, end-stamped.
   real(wp), parameter :: ARCHIVE_DT_SEC = 3600.0_wp        !--- [s]
   !--- Search bound of the whole-year test on the recycle window [calendar years].
   integer(ik), parameter :: MAX_RECYCLE_YEARS = 200_ik
   !--- [forcing].format: the forcing file, the ED_ERA5land archive, or no file.
   integer(ik), parameter :: MET_BACKEND_DEFAULT  = 1_ik  !--- the multi-grid forcing NetCDF (format = "default")
   integer(ik), parameter :: MET_BACKEND_ERA5LAND = 2_ik  !--- the monthly ED_ERA5land archive (format = "ERA5land")

   integer, parameter :: MET_PATH_LEN = 1024                !--- length of every forcing path field

   !--- Timestamp of a forcing record, avg_convention (docs/science/forcing.md §5).
   integer(ik), parameter :: METAVG_END     = 1_ik   !--- flux mean over the interval ending at the stamp (ERA5-Land)
   integer(ik), parameter :: METAVG_BEGIN   = 2_ik   !--- mean over the interval beginning at the stamp

   !--- Downwelling longwave source (docs/science/forcing.md §11).
   integer(ik), parameter :: LW_FILE       = 0_ik   !--- read from file (ERA5-Land has strd)
   integer(ik), parameter :: LW_SYNTHESIZE = 1_ik   !--- Brutsaert/Idso clear-sky synthesis (source lacking LW)

   !--- Model-start-before-first-record policy.
   integer(ik), parameter :: CLAMP_ERROR = 0_ik   !--- hard stop
   integer(ik), parameter :: CLAMP_HOLD  = 1_ik   !--- clamp to the first record (w_next = 0)


   !--- How a polygon binds to the file's grid (docs/science/forcing.md §1).
   integer(ik), parameter :: GRIDMATCH_EXPLICIT = 0_ik  !--- use grid_index verbatim (default)
   integer(ik), parameter :: GRIDMATCH_NEAREST  = 1_ik  !--- pick the grid cell nearest [site] lat/lon (great-circle)

   !---~---
   !   Where the free-atmosphere CO2 comes from (docs/science/forcing.md §12): one
   ! source for every backend, read on model time, never from the met file.
   !---~---
   integer(ik), parameter :: CO2_SOURCE_CONST = 0_ik   !--- co2_const, held for the whole run
   integer(ik), parameter :: CO2_SOURCE_FILE  = 1_ik   !--- an ED-MORE CO2 file (format: docs/science/forcing.md §12)

   !---~---
   !   The forcing's own vertical frame: the forcing is moved from its heights to the
   ! top of each patch's canopy air space (met_to_cas_top; docs/science/forcing.md
   ! §8 Step 2).
   !---~---
   integer(ik), parameter :: HEIGHT_ABOVE_ZERO_PLANE = 0_ik  !--- above the patch's displacement height (a reanalysis has no d)
   integer(ik), parameter :: HEIGHT_ABOVE_GROUND     = 1_ik  !--- heights above the ground (a flux tower)
   integer(ik), parameter :: WIND_EXPOSURE_OPEN_TERRAIN = 0_ik  !--- an open-terrain diagnostic, as ERA5's 10 m wind
   integer(ik), parameter :: WIND_EXPOSURE_LOCAL        = 1_ik  !--- the wind was measured over this canopy

   !---~---
   !   The [forcing] and [site] blocks, plain scalars. ED-MORE never gap-fills: a missing
   ! required value stops the run (docs/science/forcing.md §10).
   !---~---
   type :: forcing_params_t
      logical            :: forcing_on   = .false.               !--- master gate (indep. of fast_biophysics_on)
      integer(ik)        :: format       = MET_BACKEND_DEFAULT !--- MET_BACKEND_*: "default" | "ERA5land"
      character(len=MET_PATH_LEN) :: path = ''                   !--- forcing NetCDF path (format = "default")
      !--- The archive (docs/science/forcing.md §1); '' derives from data_path.
      character(len=MET_PATH_LEN) :: data_path     = ''          !--- archive folder
      character(len=MET_PATH_LEN) :: file_template = ''          !--- '' -> {data_path}/ED_ERA5land_{var}_{yyyy}{mm}.nc
      character(len=MET_PATH_LEN) :: static_file   = ''          !--- '' -> {data_path}/ED_ERA5land_static.nc
      real(wp)           :: max_distance_km = 15.0_wp            !--- [km] reach of the search for the nearest valid cell
      integer(ik)        :: grid_index   = 1_ik                  !--- which (time,grid) location this polygon reads
      real(wp)           :: timestep     = 3600.0_wp             !--- [s] native interval (hourly ERA5-Land)
      integer(ik)        :: avg_convention = METAVG_END          !--- METAVG_*; ERA5-Land's flux means end at the stamp
      integer(ik)        :: sw_partition = SWPART_CLEARIDX        !--- SWPART_*; ERA5-Land carries total SW only
      integer(ik)        :: lwdown_source = LW_FILE              !--- file (ERA5-Land strd) | synthesize
      !--- LW synthesis (docs/science/forcing.md §11, eq. (17)).
      integer(ik)        :: lw_clear_form = LW_CLEAR_BRUTSAERT    !--- LW_CLEAR_BRUTSAERT | LW_CLEAR_IDSO
      real(wp)           :: lw_cloud_a    = 0.22_wp               !--- [-] a of the cloud term (1 + a(1-kt)); 0: clear sky only
      integer(ik)        :: co2_source   = CO2_SOURCE_CONST      !--- "const" | "file"
      real(wp)           :: co2_const    = 420.0_wp              !--- [umol/mol] co2_source = "const"
      character(len=MET_PATH_LEN) :: co2_file = ''               !--- the ED-MORE CO2 file (co2_source = "file")
      !--- Recycling, with its declared window (docs/science/forcing.md §9).
      logical            :: recycle      = .false.               !--- cycle the record when the run outruns the file
      type(ed_time_t)    :: recycle_start                        !--- first instant of the cycle (inclusive)
      type(ed_time_t)    :: recycle_end                          !--- first instant after the cycle (exclusive)
      integer(ik)        :: start_clamp  = CLAMP_ERROR           !--- model start < base_time: error | hold the first record
      integer(ik)        :: grid_match   = GRIDMATCH_EXPLICIT    !--- explicit grid_index | nearest [site] lat/lon
      !--- [site]; every forcing clock is UTC (docs/science/forcing.md eq. (2)).
      real(wp)           :: latitude      = 42.44_wp             !--- [deg +N] Ithaca NY
      real(wp)           :: longitude     = -76.50_wp            !--- [deg +E] Ithaca NY (solar time)
      real(wp)           :: elevation     = 320.0_wp             !--- [m] site elevation (Ithaca ~320 m; elevation lapse)
      !--- Forcing heights, ERA5-Land's by default (docs/science/forcing.md §8).
      real(wp)           :: tq_height            = 2.0_wp       !--- [m] height of air temperature and humidity
      real(wp)           :: wind_height          = 10.0_wp      !--- [m] height of the wind
      integer(ik)        :: height_above         = HEIGHT_ABOVE_ZERO_PLANE   !--- what those heights are above
      integer(ik)        :: wind_exposure        = WIND_EXPOSURE_OPEN_TERRAIN !--- open-terrain diagnostic | local
      real(wp)           :: wind_exposure_z0     = 0.03_wp      !--- [m] roughness the open-terrain wind was made for
      real(wp)           :: wind_blending_height = 40.0_wp      !--- [m] height it was brought down from
      !--- Terrain lapse (docs/science/forcing.md eqs. (13)-(14)).
      logical            :: apply_elevation_lapse = .false.      !--- lapse T, P, q, LW from the grid-cell elevation to the site
      real(wp)           :: lapse_rate_tair(12) = 0.0065_wp      !--- [K/m] monthly environmental lapse (> 0: cooling upward)
      real(wp)           :: grid_elevation    = 320.0_wp         !--- [m] elevation of the forcing source grid cell
   end type forcing_params_t

end module ed_forcing_params
