! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_params -- the [output] block (docs/science/diagnostics.md) and the codes it
! uses: temporal tiers (FREQ_*), variable groups (GRP_*), file chunks (FC_*) and sync
! (SYNC_*). A config leaf (ed_kinds only), so ed_config carries
! output_params_t with no edge back to the io layer; the io-only codes (AGG_*, DIM_*) are
! in ed_output_types.
!---~---
module ed_output_params
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: output_params_t
   public :: FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE, N_FREQ
   public :: TIER_FAST, TIER_DAILY, TIER_MONTHLY, TIER_ANNUAL
   public :: GRP_STRUCTURE, GRP_CARBON, GRP_WATER, GRP_ENERGY, GRP_RADIATION, GRP_ECOPHYS,       &
             GRP_BIOGEOCHEM, GRP_NUMERICS, GRP_FORCING, N_GRP
   public :: FC_DAY, FC_MONTH, FC_YEAR, FC_ALL
   public :: SYNC_NEVER, SYNC_FLUSH
   public :: STATE_DAYS, STATE_MONTHS, STATE_YEARS
   public :: MAX_DBH_EDGE
   public :: freq_tier_index, freq_letter

   !---~---
   !   Temporal tiers, as bits of a per-variable membership mask (ior/iand); a FAST
   ! record spans fast_interval_steps sub-steps (docs/science/diagnostics.md §4).
   !---~---
   integer(ik), parameter :: FREQ_NONE    = 0_ik
   integer(ik), parameter :: FREQ_FAST    = 1_ik
   integer(ik), parameter :: FREQ_DAILY   = 2_ik
   integer(ik), parameter :: FREQ_MONTHLY = 4_ik
   integer(ik), parameter :: FREQ_ANNUAL  = 8_ik
   integer(ik), parameter :: N_FREQ       = 4_ik   ! number of temporal tiers (FAST/DAILY/MONTHLY/ANNUAL)
   !--- Each tier's index into the per-tier arrays (bit FREQ_* = 2**(index - 1)).
   integer(ik), parameter :: TIER_FAST = 1_ik, TIER_DAILY = 2_ik, TIER_MONTHLY = 3_ik, TIER_ANNUAL = 4_ik

   !---~---
   !   Variable groups, the headings of the --dump-io-config list (docs/science/diagnostics.md
   ! §5).
   !---~---
   integer(ik), parameter :: GRP_STRUCTURE = 1_ik  ! demography/size: nplant, dbh, agb, patch geometry, ids
   integer(ik), parameter :: GRP_CARBON    = 2_ik  ! GPP, NPP, Ra, allocation, vegetation carbon stocks
   integer(ik), parameter :: GRP_WATER     = 3_ik  ! ET, soil moisture and potential, snow, surface water
   integer(ik), parameter :: GRP_ENERGY    = 4_ik  ! LE, H, net radiation, temperatures, canopy-air state
   integer(ik), parameter :: GRP_NUMERICS  = 5_ik  ! integrator work counters and the energy and water budget residuals
   integer(ik), parameter :: GRP_RADIATION  = 6_ik  ! albedo, absorbed SW, PAR and LW, ground radiation
   integer(ik), parameter :: GRP_ECOPHYS    = 7_ik  ! per-cohort leaf gas exchange and hydraulics (the highest volume)
   integer(ik), parameter :: GRP_BIOGEOCHEM = 8_ik  ! soil-carbon pools, Rh, litterfall
   integer(ik), parameter :: GRP_FORCING    = 9_ik  ! the forcing the run used, after the reader's partitions
   integer(ik), parameter :: N_GRP          = 9_ik

   !--- DBH size-class edge ceiling (n_class <= MAX_DBH_EDGE-1).
   integer(ik), parameter :: MAX_DBH_EDGE = 33_ik

   !---~---
   !   File-chunk buckets: the span of one file of a tier ([output.<tier>].file_chunk), at least
   ! the tier's own period.
   !---~---
   integer(ik), parameter :: FC_DAY   = 1_ik  ! one file per sim-day
   integer(ik), parameter :: FC_MONTH = 2_ik  ! one file per sim-month
   integer(ik), parameter :: FC_YEAR  = 3_ik  ! one file per sim-year
   integer(ik), parameter :: FC_ALL   = 4_ik  ! one file for all of the run's records (annual default)

   !---~---
   !   The unit of [output].state_interval: a state file every N days, N month boundaries or N
   ! year boundaries.
   !---~---
   integer(ik), parameter :: STATE_DAYS   = 1_ik
   integer(ik), parameter :: STATE_MONTHS = 2_ik
   integer(ik), parameter :: STATE_YEARS  = 3_ik

   !--- nc_sync durability policy: never, or after every flush.
   integer(ik), parameter :: SYNC_NEVER = 0_ik  ! no nc_sync (fastest; crash forfeits the open chunk)
   integer(ik), parameter :: SYNC_FLUSH = 1_ik  ! nc_sync after each flush (durable; default)

   !---~---
   !   The [output] block: plain scalars and fixed arrays. What is written is the registry's
   ! defaults and the lines of [output].io_config (docs/science/diagnostics.md §5).
   !---~---
   type :: output_params_t
      logical            :: enabled     = .false.               ! master switch
      character(len=256) :: dir         = '.'                   ! output directory
      character(len=256) :: prefix      = 'ED'                  ! filename stem
      character(len=256) :: io_config   = ''                    ! path to ed_io_config.toml ('' -> the defaults)
      integer(ik)        :: cohort_max  = 4096_ik               ! cohort slab size: a run stops if its live count exceeds it
      integer(ik)        :: patch_max   = 256_ik                ! patch slab size: likewise
      integer(ik)        :: sync_every  = SYNC_FLUSH            ! nc_sync policy: SYNC_NEVER | SYNC_FLUSH
      logical            :: overwrite   = .false.               ! replace existing diagnostic files instead of stopping
      integer(ik)        :: fast_interval_steps = 4_ik          ! fast tier flushes every N*dt_fast
      integer(ik)        :: diel_interval_hours = 1_ik          ! a diel_ variable's bins [h]: 1, 2, 3 or 4
      !--- Per-tier enable and file chunk, indexed FAST, DAILY, MONTHLY, ANNUAL.
      logical            :: freq_on(N_FREQ)    = [.false., .true., .true., .true.]
      integer(ik)        :: file_chunk(N_FREQ) = [FC_DAY, FC_MONTH, FC_YEAR, FC_ALL]
      !---~---
      !   The state files (<state_file_prefix>-S-<date>.nc), the model's whole state to restart
      ! from: every state_every days, months or years (state_unit) and at the run's end.
      !---~---
      logical            :: state_write       = .false.
      integer(ik)        :: state_every       = 1_ik
      integer(ik)        :: state_unit        = STATE_YEARS
      character(len=256) :: state_file_dir    = '.'            ! their directory, and the PFT table's
      character(len=256) :: state_file_prefix = 'ED'           ! their file-name stem
      !--- DBH class edges [cm]; empty here, setup_output_files installs the default.
      integer(ik)        :: n_dbh_class = 0_ik
      real(wp)           :: dbh_class_edges(MAX_DBH_EDGE) = 0.0_wp
   end type output_params_t

   !---~---
   !   Default DBH class edges [cm], the common tropical forest-inventory set: 0-10,
   ! 10-20, 20-30, 30-50, 50-70 and 70+, the last open at the top (dbh_class_index;
   ! docs/science/diagnostics.md §3).
   !---~---
   integer(ik), parameter, public :: N_DBH_CLASS_DEFAULT = 6_ik
   real(wp),    parameter, public :: DBH_EDGES_DEFAULT(7) =                                      &
      [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp, 50.0_wp, 70.0_wp, 100.0_wp]

contains

   !--- Tier index (1..N_FREQ) of a single FREQ_* bit; 0 if not a single tier bit.
   pure integer(ik) function freq_tier_index(freq_bit) result(idx)
      integer(ik), intent(in) :: freq_bit
      select case (freq_bit)
      case (FREQ_FAST)
         idx = 1_ik
      case (FREQ_DAILY)
         idx = 2_ik
      case (FREQ_MONTHLY)
         idx = 3_ik
      case (FREQ_ANNUAL)
         idx = 4_ik
      case default
         idx = 0_ik
      end select
   end function freq_tier_index

   !--- Single-letter stream tag F/D/M/Y for a FREQ_* bit (file naming).
   pure character(len=1) function freq_letter(freq_bit) result(c)
      integer(ik), intent(in) :: freq_bit
      select case (freq_bit)
      case (FREQ_FAST)
         c = 'F'
      case (FREQ_DAILY)
         c = 'D'
      case (FREQ_MONTHLY)
         c = 'M'
      case (FREQ_ANNUAL)
         c = 'Y'
      case default
         c = '?'
      end select
   end function freq_letter

end module ed_output_params
