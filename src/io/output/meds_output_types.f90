! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_output_types -- the plain-data types of the diagnostic output: a variable's descriptor  !
! (var_desc_t), the running reduction of one (variable, tier) pair (integ_buffer_t), the      !
! registry, a closed record and its queue, a stream file's handles, and the run's file set    !
! and one polygon's buffers (output_files_t, output_buffers_t); also the codes AGG_* and      !
! DIM_* (FREQ_*, GRP_* and FC_* are in meds_output_params). netCDF-free                       !
! (docs/science/diagnostics.md §1). Design record: MEDS_IO_DESIGN.md §3.                      !
!==========================================================================================!
module meds_output_types
   use meds_kinds,         only : wp, ik
   use meds_time,          only : meds_time_t
   use meds_output_params, only : N_FREQ
   implicit none
   private

   public :: var_desc_t, integ_buffer_t, output_registry_t, diag_setup_t
   public :: pending_record_t, record_queue_t, stream_file_t
   public :: output_files_t, output_buffers_t
   public :: AGG_MEAN, AGG_RATE, AGG_LAST, AGG_TMEAN
   public :: DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE
   public :: XTYPE_DOUBLE, XTYPE_INT
   public :: MISSING_VALUE, MISSING_INT, MAX_OUTPUT_VARS, MAX_DBH_CLASS
   public :: ragged_dim, slab_col

   !----- Temporal reduction operators (the `agg` of a variable): a mean over the period, or its   !
   !      value at the period's end (docs/science/diagnostics.md §4). ----------------------------!
   integer(ik), parameter :: AGG_MEAN  = 1_ik  !< equal-weight arithmetic mean of the samples
   integer(ik), parameter :: AGG_RATE  = 2_ik  !< mean rate of a per-step amount: Sum(amount) / Sum(dt)
   integer(ik), parameter :: AGG_LAST  = 5_ik  !< end-of-period snapshot (ids / CSR / instantaneous)
   integer(ik), parameter :: AGG_TMEAN = 6_ik  !< dt-weighted state mean (the physical-stock default)

   !----- The trailing axis of a variable. The cohort and patch slots are fixed within a window     !
   !      (docs/science/diagnostics.md §4), so every axis folds by direct slot index. -------------!
   integer(ik), parameter :: DIM_SCALAR = 0_ik
   integer(ik), parameter :: DIM_COHORT = 1_ik
   integer(ik), parameter :: DIM_PATCH  = 2_ik
   !----- (soil layer x site): each site's column on its own layers, flattened as DIM_SOIL_PATCH    !
   !      is below, at (is-1)*n_soil_layer_max + k, and written (time, polygon, site, soil). -------!
   integer(ik), parameter :: DIM_SOIL   = 3_ik
   integer(ik), parameter :: DIM_PFT    = 4_ik
   integer(ik), parameter :: DIM_SIZE   = 5_ik   !< DBH size class (edges from [output].dbh_class_edges)
   !----- 2-D (soil layer x patch), flattened into a 1-D slab at (ip-1)*n_soil_layer_max + k, so the  !
   !      slab routines serve it unchanged. The stride is the layer ceiling, never a live count      !
   !      (docs/science/diagnostics.md §3); the serializer writes it (time, polygon, patch, soil). --!
   integer(ik), parameter :: DIM_SOIL_PATCH = 6_ik
   !----- One value per site of the polygon (time, polygon, site): what differs between sites,     !
   !      their forcing after their terrain and their growth temperature. -------------------------!
   integer(ik), parameter :: DIM_SITE = 7_ik

   !----- On-disk numeric type, declared here to keep the module netCDF-free; the serializer maps   !
   !      XTYPE_DOUBLE to NC_DOUBLE and XTYPE_INT to NC_INT. -------------------------------------!
   integer(ik), parameter :: XTYPE_DOUBLE = 1_ik
   integer(ik), parameter :: XTYPE_INT    = 2_ik

   !----- Fill values: those of netCDF's NC_FILL_DOUBLE and NC_FILL_INT, so a reader recognizes     !
   !      them, declared here to keep the module netCDF-free. -------------------------------------!
   real(wp),    parameter :: MISSING_VALUE = 9.9692099683868690e+36_wp
   integer(ik), parameter :: MISSING_INT   = -2147483647_ik

   !----- Registry capacity: the registered set (docs/science/diagnostics.md §7) plus headroom; the !
   !      list is short enough for a linear find_var_index. --------------------------------------!
   integer(ik), parameter :: MAX_OUTPUT_VARS = 512_ik
   !----- DBH size-class ceiling (DIM_SIZE). Edges come from [output].dbh_class_edges. ---------!
   integer(ik), parameter :: MAX_DBH_CLASS   = 32_ik

   !==========================================================================================!
   ! The run-dependent parameters the derived diagnostics need (the PFT count, the DBH class  !
   ! edges, whether the soil has a curve to read), resolved in setup_output_files and handed  !
   ! to extract_variable as plain data, so the per-step tick does not need the configuration. !
   !==========================================================================================!
   type :: diag_setup_t
      integer(ik) :: n_pft       = 0_ik      !< run-time PFT count (the DIM_PFT axis length)
      !----- The fast loop runs, so the soil water has a curve to read (attach_output). ----------!
      logical     :: soil_ready  = .false.   !< .false. => psi/wetness diagnostics emit _FillValue
      !----- DBH size classes (DIM_SIZE): n_class+1 ascending edges [cm]. ----------------------!
      integer(ik) :: n_dbh_class = 0_ik
      real(wp)    :: dbh_edges(MAX_DBH_CLASS + 1_ik) = 0.0_wp
   end type diag_setup_t

   !==========================================================================================!
   ! One diagnostic variable: plain data, no pointers into the state.                         !
   !==========================================================================================!
   type :: var_desc_t
      character(len=32) :: name      = ''          !< netCDF variable name (scale-suffixed, e.g. 'mean_agb_co')
      character(len=96) :: long_name = ''          !< CF long_name attribute
      character(len=24) :: units     = ''          !< CF units attribute
      integer(ik) :: dim             = DIM_SCALAR   !< DIM_* : the trailing axis
      integer(ik) :: agg             = AGG_MEAN     !< AGG_* : temporal reduction operator
      integer(ik) :: group           = 1_ik         !< GRP_* : the group a high-level toggle switches
      integer(ik) :: streams         = 0_ik         !< ior(FREQ_*) effective membership (after config resolution)
      integer(ik) :: streams_default = 0_ik         !< the registry default membership (restored by a `true` override)
      logical     :: enabled         = .true.       !< master per-variable on/off
      integer(ik) :: xtype           = XTYPE_DOUBLE !< XTYPE_DOUBLE | XTYPE_INT
      integer(ik) :: source_id       = 0_ik         !< the field the D, M and Y tiers read (0: none)
      integer(ik) :: source_fast     = 0_ik         !< the staged sub-step value the F tier reads (0: none)
      logical     :: flux            = .false.      !< a flux or a rate: written as its mean only, no snapshot
      logical     :: coord           = .false.      !< the axis's bookkeeping (ids, CSR): written with the axis,
                                                    !  not selected (resolve_streams)
      logical     :: snap_staged     = .false.      !< a snapshot the D, M and Y tiers read from the last staged
                                                    !  sub-step (a fast-loop diagnostic's end-of-period value)
      !----- How the field reduces to this variable's axis, so one field serves every scale. The   !
      !      weight kind and mean together state extensive vs intensive (docs/science/diagnostics.md §2). !
      integer(ik) :: weight          = 0_ik         !< W_* per-cohort weight kind (0 = W_NONE)
      logical     :: mean            = .false.      !< .true. weighted mean (intensive); .false. sum (extensive)
      real(wp)    :: scale           = 1.0_wp       !< unit conversion of source_id, applied where it is read
   end type var_desc_t

   !==========================================================================================!
   ! One running reduction of one (variable, tier) pair; an axis folds by direct slot index.  !
   !==========================================================================================!
   type :: integ_buffer_t
      integer(ik) :: var_id = 0_ik          !< index into registry%var(:)
      integer(ik) :: freq   = 0_ik          !< the single FREQ_* bit this buffer serves
      integer(ik) :: agg    = AGG_MEAN
      integer(ik) :: dim    = DIM_SCALAR
      logical     :: active = .false.       !< allocated + participating (this var is on in this tier)
      !----- scalar accumulators (DIM_SCALAR). -------------------------------------------------!
      real(wp)    :: scal  = 0.0_wp         !< running reduction
      real(wp)    :: wsum  = 0.0_wp         !< Sum(dt) weight (AGG_TMEAN, AGG_RATE)
      integer(ik) :: nsamp = 0_ik           !< equal-weight sample count (AGG_MEAN, AGG_LAST)
      !----- per-slot accumulators of an axis variable, sized to the axis's capacity. ------------!
      real(wp),    allocatable :: slab(:)
      real(wp),    allocatable :: wsum_slab(:)
      integer(ik), allocatable :: hits(:)
      integer(ik) :: n_slab = 0_ik          !< live index count this window (the [1:n] cohort/patch count)
   end type integ_buffer_t

   !==========================================================================================!
   ! The registration list and, per tier, the indices of its live variables.                  !
   !==========================================================================================!
   type :: output_registry_t
      type(var_desc_t), allocatable :: var(:)      !< the registration list [1:nvar]
      integer(ik)                   :: nvar = 0_ik
      integer(ik), allocatable      :: idx_freq(:,:) !< (MAX_OUTPUT_VARS, N_FREQ): var indices live in each tier
      integer(ik)                   :: nidx(N_FREQ) = 0_ik
      logical                       :: tier_on(N_FREQ) = .true.  !< [output.<tier>].enabled: off writes nothing
   end type output_registry_t

   !==========================================================================================!
   ! A closed period. close_tier normalizes a tier into its scratch record and copies that to    !
   ! the tier's queue, which output_serialize drains in the I/O phase                            !
   ! (docs/science/order_of_processes.md §6). Indexed by registry variable; only the tier's live !
   ! variables are filled.                                                                       !
   !==========================================================================================!
   type :: pending_record_t
      logical           :: used     = .false.
      integer(ik)       :: freq     = 0_ik      !< the FREQ_* bit of the tier this record closes
      type(meds_time_t) :: t_open              !< period-start stamp (calendar companions)
      integer(ik)       :: n_cohort = 0_ik, n_patch = 0_ik  !< live slab lengths this record
      real(wp),    allocatable :: sval(:)       !< (nvar) normalized scalar value
      logical,     allocatable :: svalid(:)     !< (nvar)
      real(wp),    allocatable :: slab(:,:)     !< (rows, columns) normalized slab values
      logical,     allocatable :: slabvalid(:,:)!< (rows, columns)
      integer(ik), allocatable :: nslab(:)      !< (nvar) slab length (0 for scalar vars)
      !----- A record stores only its tier's live slab variables: slab(:, col(k)) is registry      !
      !      variable k's slab (slab_col), in the scratch record and the queued ones alike. ----------!
      integer(ik), allocatable :: col(:)        !< (nvar) slab column of each variable (0 = none)
   end type pending_record_t

   !----- The closed records of one tier awaiting the I/O phase, in closing order. A queued record  !
   !      keeps only the slab rows the writer reads, so a month of fast-tier records costs the live !
   !      cohort count, not cohort_max. Elements are reused across months; `n` counts the live ones. !
   type :: record_queue_t
      integer(ik) :: n = 0_ik
      type(pending_record_t), allocatable :: rec(:)
   end type record_queue_t

   !==========================================================================================!
   ! One open netCDF stream file, the serializer's per-tier handle. Its ids are plain integers  !
   ! so the type stays netCDF-free; the serializer casts them to c_int.                         !
   !==========================================================================================!
   type :: stream_file_t
      integer(ik) :: freq         = 0_ik
      integer(ik) :: ncid         = -1_ik   !< open handle (-1 = none)
      integer(ik) :: nrec         = 0_ik    !< records written to the current file
      integer(ik) :: chunk_bucket = -1_ik   !< current time-chunk bucket key (-1 = no file open)
      logical     :: has_cohort   = .false. !< this tier defines the cohort dim
      logical     :: has_patch    = .false. !< this tier defines the patch dim
      logical     :: has_soil     = .false. !< this tier defines the soil dim
      logical     :: has_pft      = .false. !< this tier defines the pft dim
      logical     :: has_size     = .false. !< this tier defines the dbh-class dim
      integer(ik) :: d_time = -1_ik, d_cohort = -1_ik, d_patch = -1_ik, d_soil = -1_ik
      integer(ik) :: d_pft = -1_ik, d_size = -1_ik
      integer(ik) :: d_polygon = -1_ik       !< a region file's polygon dimension
      integer(ik) :: d_site = -1_ik          !< the site dimension, behind the polygon's
      integer(ik) :: v_pft = -1_ik, v_dbh_lower = -1_ik, v_dbh_upper = -1_ik  !< self-describing axis coords
      integer(ik) :: v_soil_z = -1_ik                                        !< soil layer node depths [m]
      integer(ik) :: cohort_dim = 0_ik, patch_dim = 0_ik   !< the file's trimmed cohort/patch axis lengths
      integer(ik) :: v_time = -1_ik, v_year = -1_ik, v_month = -1_ik, v_day = -1_ik
      integer(ik) :: v_hour = -1_ik, v_minute = -1_ik, v_second = -1_ik   !< FAST-tier sub-daily companions
      integer(ik) :: v_ncohort = -1_ik, v_npatch = -1_ik
      integer(ik), allocatable :: vid(:)    !< (nvar) netCDF varid per registry var (-1 if not in tier)
   end type stream_file_t

   !==========================================================================================!
   ! The run's output, as two plain-data types (only output_serialize touches netCDF):         !
   !   output_files_t   -- the run's one file set: a row per polygon, the registry, the        !
   !                       diagnostic parameters, file settings and stream handles; read-only  !
   !                       while a step runs.                                                  !
   !   output_buffers_t -- one polygon's buffers (one set per polygon): its running reductions, !
   !                       open windows, closed records and fast-tier staging, the only output !
   !                       state a step writes.                                                !
   !==========================================================================================!
   type :: output_files_t
      logical                 :: enabled = .false.
      !----- The run has a fast loop. Without one the fast-loop diagnostic blocks stay off, so their  !
      !      variables read as missing rather than as a 0 nothing computed. -------------------------!
      logical                 :: fast_loop_on = .true.
      type(output_registry_t) :: reg
      type(diag_setup_t)      :: diag       !< run-dependent params the derived diagnostics need
      integer(ik)       :: cohort_max = 0_ik, patch_max = 0_ik, max_slab = 0_ik
      logical           :: cohort_axis = .false., patch_axis = .false.  !< a live variable has that axis
      type(stream_file_t)    :: stream(N_FREQ)
      character(len=256)     :: dir = '.', prefix = 'meds'
      !----- Forcing provenance, a global attribute on every file: from the ED_ERA5land archive the  !
      !      model computes qair itself, so the formula is recorded (docs/science/forcing.md §7). --!
      character(len=512)     :: forcing_qair = ''
      integer(ik)           :: file_chunk(N_FREQ) = 0_ik
      integer(ik)           :: sync_every = 1_ik
      integer(ik)           :: fast_interval_steps = 4_ik   !< fast tier closes every N*dt_fast sub-steps
      !----- The polygon axis every file carries: one entry per polygon's buffers, with these        !
      !      coordinates (a site run is one polygon, id 0, at row 0 col 0 unless its forcing file   !
      !      places it). ---------------------------------------------------------------------------!
      integer(ik)              :: n_polygon = 0_ik
      integer(ik), allocatable :: polygon_id(:)             !< the cell's row-major index on the forcing grid
      integer(ik), allocatable :: polygon_row(:), polygon_col(:)   !< 0-based grid indices
      real(wp),    allocatable :: polygon_lat(:), polygon_lon(:)   !< [deg] cell centre
      !----- Each polygon's sites, as it runs them, recorded by attach_output: each one's area,     !
      !      terrain and soil column, the files' soil coordinate and static variables. The site     !
      !      axis is the largest polygon's; a polygon's other entries hold the fill value. ---------!
      integer(ik)              :: site_max = 1_ik
      integer(ik), allocatable :: polygon_n_site(:)          !< the polygon's sites
      real(wp),    allocatable :: site_area(:,:)             !< (site, polygon) [-] fraction of the polygon
      real(wp),    allocatable :: site_elevation(:,:)        !< (site, polygon) [m]
      real(wp),    allocatable :: site_slope(:,:)            !< (site, polygon) [degrees]
      real(wp),    allocatable :: site_aspect(:,:)           !< (site, polygon) [degrees from north]
      integer(ik), allocatable :: site_n_soil(:,:)           !< (site, polygon) active soil layers
      real(wp),    allocatable :: site_soil_z(:,:,:)         !< (layer, site, polygon) [m] node depths, <= 0
      real(wp),    allocatable :: site_soil_face(:,:,:)      !< (layer + 1, site, polygon) [m] layer faces, <= 0
      real(wp),    allocatable :: site_soil_depth(:,:)       !< (site, polygon) [m] soil column depth
      integer(ik), allocatable :: site_bottom_bc(:,:)        !< (site, polygon) the bottom water boundary (SOIL_BC_*)
      real(wp),    allocatable :: site_drain_factor(:,:)     !< (site, polygon) [-] drainage's f; MISSING_VALUE on an aquifer
      real(wp),    allocatable :: site_deep_temp(:,:)        !< (site, polygon) [K] deep soil temperature
      real(wp),    allocatable :: polygon_chave_e(:)         !< [-] Chave's E; only with that allometry
   end type output_files_t

   type :: output_buffers_t
      type(integ_buffer_t), allocatable :: buf(:,:)   !< (nvar, N_FREQ) running reductions
      logical           :: has_data(N_FREQ) = .false. !< tier's current window has >=1 sample
      type(meds_time_t) :: t_open(N_FREQ)             !< period-start of each tier's current window
      type(pending_record_t) :: pending(N_FREQ)          !< per-tier scratch that close_tier normalizes into
      type(record_queue_t)   :: queue(N_FREQ)            !< closed records awaiting the I/O phase
      !----- FAST-tier staging, filled per (patch, sub-step) by the fast loop and replayed into        !
      !      buf(:, TIER_FAST) by output_integrate_fast (docs/science/diagnostics.md §4): each patch's !
      !      row of the patch block (PA_*) and soil column, their polygon mean, the forcing (PY_*)    !
      !      and the per-cohort rows by the polygon's cohort slot. They follow the live counts, which  !
      !      the output layer checks against output.patch_max and cohort_max. ---------------------!
      logical              :: fast_on = .false.                !< output on, and the F tier or a snapshot reads the staging
      type(meds_time_t),   allocatable :: fast_time(:)         !< (n_fast_sub) each sub-step's start
      real(wp),            allocatable :: fast_forcing(:,:)     !< (PY_*, n_fast_sub) the sub-step's forcing
      real(wp),            allocatable :: fast_patch(:,:,:)     !< (PA_*, patch, n_fast_sub) each patch's row
      real(wp),            allocatable :: fast_site(:,:)        !< (PA_*, n_fast_sub) the polygon mean of the rows
      !----- The soil columns, per patch in the soil-by-patch slab's layout, (ip-1)*n_soil_layer_max  !
      !      + k (DIM_SOIL_PATCH), and per site, each site's patches area-weighted, at                !
      !      (is-1)*n_soil_layer_max + k (DIM_SOIL). The masks mark the layers each column has. ------!
      real(wp),            allocatable :: fast_soil_temp_patch(:,:)   !< (layer x patch, n_fast_sub) [K]
      real(wp),            allocatable :: fast_soil_water_patch(:,:)  !< (layer x patch, n_fast_sub) [m3/m3]
      real(wp),            allocatable :: fast_soil_temp(:,:)   !< (layer x site, n_fast_sub) [K]
      real(wp),            allocatable :: fast_soil_water(:,:)  !< (layer x site, n_fast_sub) [m3/m3]
      logical,             allocatable :: fast_soil_live(:)        !< (layer x site) a layer of its site
      logical,             allocatable :: fast_soil_patch_live(:)  !< (layer x patch) a layer of its patch's site
      real(wp),            allocatable :: fast_cohort_gpp(:,:)     !< (cohort, n_fast_sub) per-cohort GPP [umol/plant/s]
      !< (cohort, n_fast_sub) per-cohort height [m] (tallest post-proc)
      real(wp),            allocatable :: fast_cohort_height(:,:)
      !< (cohort, CO_*, n_fast_sub) the per-cohort diagnostic block's rows, staged while it is on
      real(wp),            allocatable :: fast_cohort_diag(:,:,:)
      integer(ik)          :: n_fast_sub   = 0_ik              !< sub-steps staged this slow step
      integer(ik)          :: fast_n_site  = 0_ik              !< sites in the fast soil slabs
      integer(ik)          :: fast_n_cohort = 0_ik             !< the polygon's live cohorts in the fast cohort slabs
      integer(ik)          :: fast_n_patch  = 0_ik             !< live patches in the fast patch rows
      logical              :: fast_ready   = .false.           !< the staging holds this slow step's sub-steps
      integer(ik)          :: fast_step_total = 0_ik           !< sub-steps replayed over the run (output_tick)
      !----- The diagnostic parameters this polygon's extraction reads: its file set's, with         !
      !      soil_ready set when its fast loop runs (attach_output). --------------------------------!
      type(diag_setup_t)   :: diag
   end type output_buffers_t

contains

   !----- The slab column of registry variable k in a record (identity for the scratch record). --!
   pure integer(ik) function slab_col(pr, k) result(c)
      type(pending_record_t), intent(in) :: pr
      integer(ik),            intent(in) :: k
      if (allocated(pr%col)) then
         c = pr%col(k)
      else
         c = k
      end if
   end function slab_col


   !----- The axes whose length and slots follow the stand: cohort, patch, and soil by patch. Fixed  !
   !      only within a month (docs/science/diagnostics.md §3), so a file of them spans at most a   !
   !      month, the annual stream refuses them, and a file pads them to its largest live count. ---!
   pure logical function ragged_dim(dim) result(yes)
      integer(ik), intent(in) :: dim
      yes = dim == DIM_COHORT .or. dim == DIM_PATCH .or. dim == DIM_SOIL_PATCH
   end function ragged_dim

end module meds_output_types
