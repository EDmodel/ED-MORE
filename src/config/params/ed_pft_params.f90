! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_pft_params -- the plant functional type (PFT) trait table, filled from the PFT file
! (docs/configuration.md, "The PFT file").
!
!   Structure-of-arrays: one allocatable array per trait, indexed by PFT. Each PFT's size
! allometry (allom, ed_allometry_lib), hydraulics (hydro, plant_hydraulics_t) and
! phenology (pheno) are records. Every trait is the PFT file's: a trait that follows a
! fitted relation (the sapwood conductivity or the mortality hazards from wood density,
! the plasticity slopes from Vcmax) is computed outside the model, by scripts/prepare_input/pft_traits.
!
!   The table also holds three tables derived from the traits once at load
! (ed_pft_derive): each PFT's leaf-photosynthesis parameters, its hydraulics with its
! Kirchhoff lookup, and the two-stream's optics. A per-PFT table is <process>_table_t and
! its entries <process>_pft_t.
!---~---
module ed_pft_params
   use ed_kinds,            only : wp, ik
   use ed_allometry_lib,    only : allometry_t
   use ed_plant_hydr_lib,   only : water_curve_t
   implicit none
   private

   public :: pft_table_t, phenology_traits_t
   public :: CUE_NONE, CUE_TEMP, CUE_LIGHT, CUE_WATER, CUE_ALL
   public :: leaf_photo_pft_t, leaf_photo_table_t, hydraulics_pft_t, hydraulics_table_t, rad_pft_optics_t
   public :: plant_hydraulics_t, kirchhoff_table_t, KIRCHHOFF_TABLE_NTAB, KIRCHHOFF_TABLE_RMAX
   public :: HYD_UNSET
   public :: PATH_C3, PATH_C4
   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT

   !--- Photosynthetic-pathway flags (per-PFT trait values).
   integer(ik), parameter :: PATH_C3 = 1_ik   ! C3 (Farquhar-von Caemmerer-Berry)
   integer(ik), parameter :: PATH_C4 = 2_ik   ! C4 (Collatz et al. 1992)

   !--- The optics table's radiation bands, indices into its band dimension.
   integer(ik), parameter :: RAD_VIS = 1_ik   ! visible / PAR (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_NIR = 2_ik   ! near-infrared (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_LW  = 3_ik   ! thermal / longwave (diffuse + emission, no beam)
   integer(ik), parameter :: N_RAD_BAND_DEFAULT = 3_ik

   !---~---
   !   A PFT's leaf-photosynthesis parameters (docs/science/leaf_gas_exchange.md), flat, so
   ! the leaf kernel needs no configuration: its traits and [leaf_physiology]'s constants.
   !---~---
   type :: leaf_photo_pft_t
      integer(ik) :: pathway        ! PATH_C3 | PATH_C4
      !--- Per-PFT capacities at 25 degC and stomatal/water-stress traits.
      real(wp) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(wp) :: g0, g1, d0, quantum_yield
      !---~---
      !   theta_j is the electron-transport curvature; the theta_*_c3/c4 pairs are the co-limitation
      ! curvatures, one pair per pathway (docs/science/leaf_gas_exchange.md §2).
      !---~---
      real(wp) :: theta_j
      real(wp) :: theta_cj_c3, theta_ip_c3
      real(wp) :: theta_cj_c4, theta_ic_c4
      !---~---
      !   psi_open/psi_close [MPa]: the leaf potential at which beta_nonstomata is 1 and 0, per PFT
      ! (pft_wstress_psi_open/_close: the turgor-loss point and twice it unless the PFT file sets them).
      !---~---
      real(wp) :: lambda25, psi_open, psi_close, lambda_psi_exp, sref_stomata
      !---~---
      !   Predawn potential [MPa] below which the stomatal stress begins, eq. (10); a run sets it per
      ! PFT (pft_stomata_psi_onset); the 0 here is for a record built by hand.
      !---~---
      real(wp) :: psi_onset = 0.0_wp
      !--- Shared biochemistry constants at 25 degC + activation/deactivation terms.
      real(wp) :: kc25, ko25, gstar25
      real(wp) :: ea_kc, ea_ko, ea_gstar, ea_vcmax, ea_jmax, ea_rd
      real(wp) :: hd_vcmax, hd_jmax, hd_rd, ds_vcmax, ds_jmax, ds_rd
      real(wp) :: o2_mol_frac, absorptance, phi_psii
      real(wp) :: vpd_min = 0.05_wp   ! [kPa] the VPD the Medlyn and Katul models use at least ([leaf_physiology])
   end type leaf_photo_pft_t

   !---~---
   !   Every PFT's leaf parameters. A cohort's plastic capacities override the table's, Jmax25 and
   ! TPU25 scaling with Vcmax25 by the two ratios; the solver selectors ride along.
   !---~---
   type :: leaf_photo_table_t
      integer(ik) :: n_pft = 0_ik
      type(leaf_photo_pft_t),    allocatable :: pft(:)              ! per-PFT parameters (table capacities)
      real(wp),                  allocatable :: jmax_vcmax_ratio(:) ! [--] Jmax25 / Vcmax25 per PFT
      real(wp),                  allocatable :: tpu_vcmax_ratio(:)  ! [--] TPU25 / Vcmax25 per PFT
      integer(ik) :: stomatal_model     = 0_ik    ! SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form = 0_ik    ! TEMP_RESPONSE_ARRHENIUS | TEMP_RESPONSE_PEAKED
      integer(ik) :: colimitation       = 0_ik    ! COLIM_MIN | COLIM_QUADRATIC
      logical     :: use_boundary_layer = .false. ! couple through the leaf boundary layer (gb)
   end type leaf_photo_table_t

   !---~---
   !   One PFT's hydraulic traits, from its PFT file's [hydraulics] block
   ! (docs/science/plant_hydraulics.md, Parameters). The traits have no defaults: a
   ! PFT file gives every one, for every PFT.
   !---~---
   type :: plant_hydraulics_t
      type(water_curve_t) :: leaf_curve, wood_curve      ! the tissues' pressure-volume curves
      real(wp) :: wood_psi50           = 0.0_wp  ! [MPa, < 0] xylem potential at 50 % loss of conductance
      real(wp) :: wood_kexp            = 0.0_wp  ! [-] vulnerability-curve shape
      real(wp) :: wood_kmax            = 0.0_wp  ! [kg/m/s/MPa] sapwood specific conductivity
      real(wp) :: vessel_curl          = 0.0_wp  ! [-] xylem path length per unit height
      real(wp) :: root_beta            = 0.0_wp  ! [-] the untruncated profile's share below the rooting depth (eq. (13))
      real(wp) :: root_depth_b1        = 0.0_wp  ! [m] rooting depth = b1 * height[m]^b2
      real(wp) :: root_depth_b2        = 0.0_wp  ! [-]
      real(wp) :: specific_root_length = 0.0_wp  ! [m/kgC] fine-root length per unit fine-root carbon
      real(wp) :: fine_root_radius     = 0.0_wp  ! [m] fine-root radius
   end type plant_hydraulics_t

   !---~---
   !   The Kirchhoff lookup table for a general vulnerability exponent: G(r) on a
   ! uniform grid (docs/science/plant_hydraulics.md §3), built by
   ! build_kirchhoff_table (ed_pft_derive). A fixed-size record, so a PFT's
   ! parameters copy as plain data.
   !---~---
   integer,  parameter :: KIRCHHOFF_TABLE_NTAB = 512        ! intervals over [0, KIRCHHOFF_TABLE_RMAX]
   real(wp), parameter :: KIRCHHOFF_TABLE_RMAX = 8.0_wp     ! [-] table domain in r = psi/psi50
   type :: kirchhoff_table_t
      real(wp)    :: kexp  = 0.0_wp                           ! [-] vulnerability exponent the table holds; its tail uses it too
      real(wp)    :: r_max = KIRCHHOFF_TABLE_RMAX                 ! [-] table domain [0, r_max]
      real(wp)    :: inv_h = 0.0_wp                           ! [-] 1 / grid spacing
      integer(ik) :: n     = int(KIRCHHOFF_TABLE_NTAB, ik)        ! number of intervals
      real(wp)    :: g(0:KIRCHHOFF_TABLE_NTAB) = 0.0_wp           ! g(j) = integral_0^{j*h} du/(1+u^kexp)
   end type kirchhoff_table_t

   !--- A PFT's hydraulic traits ([hydraulics]) and the Kirchhoff lookup its wood_kexp sets.
   type, extends(plant_hydraulics_t) :: hydraulics_pft_t
      !---~---
      !   Kirchhoff lookup from wood_kexp (build_kirchhoff_table), read only for wood_kexp outside {1,2},
      ! where no closed form exists (docs/science/plant_hydraulics.md §3).
      !---~---
      type(kirchhoff_table_t) :: vuln_table
      real(wp) :: rwc_close   ! [-] leaf relative water content at the PFT's stomatal closing potential
   end type hydraulics_pft_t

   !---~---
   !   Every PFT's hydraulics; each entry carries its own Kirchhoff lookup, so each PFT keeps
   ! its vulnerability shape.
   !---~---
   type :: hydraulics_table_t
      integer(ik) :: n_pft = 0_ik
      type(hydraulics_pft_t), allocatable :: pft(:)
   end type hydraulics_table_t

   !---~---
   !   Per-PFT optics (docs/science/canopy_radiation_transfer.md eqs. (4)-(6)), leaf and
   ! wood apart, blended per cohort by clumped area; lidf is kept to evaluate G(mu).
   !---~---
   type :: rad_pft_optics_t
      integer(ik) :: n_pft  = 0_ik
      integer(ik) :: n_band = 0_ik
      real(wp), allocatable :: omega_leaf(:,:)     ! (band,pft) single-scatter albedo, leaves
      real(wp), allocatable :: omega_wood(:,:)     ! (band,pft) single-scatter albedo, wood
      real(wp), allocatable :: g_leaf(:,:)         ! (band,pft) bf*(rho-tau)/omega, leaves
      real(wp), allocatable :: g_wood(:,:)         ! (band,pft) bf*(rho-tau)/omega, wood
      real(wp), allocatable :: clumping_leaf(:)    ! (pft) leaf clumping factor (0,1]
      real(wp), allocatable :: clumping_wood(:)    ! (pft) wood clumping factor (0,1]
      real(wp), allocatable :: lidf(:,:)           ! (class,pft) leaf-angle distribution weights
      real(wp), allocatable :: bf(:)               ! (pft) <cos^2(theta_leaf)>
      logical,  allocatable :: has_beam(:)         ! (band) band has a collimated (solar) beam
      logical,  allocatable :: has_emission(:)     ! (band) band has thermal emission (LW)
   end type rad_pft_optics_t

   !---~---
   !   Cue bits: flush_cue_mask and shed_cue_mask are sums of them, selecting the cues of
   ! each side (docs/science/plant_phenology.md §2-§3).
   !---~---
   integer(ik), parameter :: CUE_NONE  = 0_ik   ! no cues (always flushing / no senescence)
   integer(ik), parameter :: CUE_TEMP  = 1_ik   ! temperature: warmth sum (flush), cold sum (shed)
   integer(ik), parameter :: CUE_LIGHT = 2_ik   ! light: running-mean hours a day of PAR above par_min
   integer(ik), parameter :: CUE_WATER = 4_ik   ! water: predawn leaf psi summed against the TLP
   integer(ik), parameter :: CUE_ALL   = CUE_TEMP + CUE_LIGHT + CUE_WATER

   !---~---
   !   A PFT's [phenology] block (docs/science/plant_phenology.md), the record the
   ! phenology kernel reads. The defaults are a PFT with no [phenology] section: no cues,
   ! so only the rates and timescales act. leaf_psi_tlp is derived at load.
   !---~---
   type :: phenology_traits_t
      integer(ik) :: flush_cue_mask        = CUE_NONE    ! cues of the flush side (a sum of CUE_*)
      integer(ik) :: shed_cue_mask         = CUE_NONE    ! cues of the shed side
      real(wp)    :: flush_cue_timescale   = 5.0_wp      ! [day] smoothing of the flush tendency
      real(wp)    :: shed_cue_timescale    = 5.0_wp      ! [day] smoothing of the shed tendency
      real(wp)    :: flush_rate_max        = 0.06667_wp  ! [1/day] leaf growth at full tendency
      real(wp)    :: shed_rate_max         = 0.05_wp     ! [1/day] senescence at full tendency
      real(wp)    :: min_leaf_cover        = 0.0_wp      ! [-] leaf cover senescence stops at
      real(wp)    :: bare_leaf_cover       = 0.02_wp     ! [-] a dormant canopy below this goes bare
      real(wp)    :: flush_base_temp       = 278.15_wp   ! [K] warmth base
      real(wp)    :: flush_degree_days     = 100.0_wp    ! [K day] warmth requirement
      real(wp)    :: flush_temp_sharpness  = 0.04_wp     ! [1/(K day)]
      real(wp)    :: shed_base_temp        = 290.15_wp   ! [K] cold base
      real(wp)    :: shed_degree_days      = 50.0_wp     ! [K day] cold requirement after midsummer
      real(wp)    :: shed_temp_sharpness   = 0.1_wp      ! [1/(K day)]
      real(wp)    :: par_min               = 5.0_wp      ! [umol/m2/s] PAR that counts as an hour of light
      real(wp)    :: flush_light_hours     = 12.0_wp     ! [h/day] hours of light that permit flushing
      real(wp)    :: flush_light_sharpness = 1.0_wp      ! [1/h]
      real(wp)    :: shed_light_hours      = 11.0_wp     ! [h/day] hours of light that trigger shedding
      real(wp)    :: shed_light_sharpness  = -1.0_wp     ! [1/h]; < 0: short or dim days trigger
      real(wp)    :: light_window          = 10.0_wp     ! [day] running mean of the hours of light
      real(wp)    :: flush_water_sum       = 10.0_wp     ! [MPa day] wet sum that permits flushing
      real(wp)    :: flush_water_sharpness = 0.5_wp      ! [1/(MPa day)]
      real(wp)    :: shed_water_sum        = 10.0_wp     ! [MPa day] dry sum that triggers shedding
      real(wp)    :: shed_water_sharpness  = 0.5_wp      ! [1/(MPa day)]
      real(wp)    :: leaf_psi_tlp          = -2.0_wp     ! [MPa] the leaf's turgor-loss point (derive_pft_tables)
   end type phenology_traits_t

   !--- PFT trait table (SoA).  Units in brackets.
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !--- The size allometry and the wood-density axis.
      type(allometry_t),        allocatable :: allom(:) ! each PFT's allometry ([allometry])
      type(plant_hydraulics_t), allocatable :: hydro(:) ! each PFT's hydraulic traits ([hydraulics])
      real(wp), allocatable :: wood_density(:)   ! [g/cm3] rho: the wood's dry mass per volume
      !--- Reproduction and recruitment.
      real(wp), allocatable :: reproduction_investment_fraction(:) ! [--] growth fraction diverted to reproduction
      real(wp), allocatable :: min_repro_height(:)                 ! [m]  the height a cohort reproduces from
      real(wp), allocatable :: recruit_carbon_efficiency(:)        ! [--] share of reproduction carbon that becomes recruits
      real(wp), allocatable :: recruit_shade_decay(:)              ! [1/(m2/m2)] that share falls as exp(-k * patch LAI)
      real(wp), allocatable :: nbg_init_density(:)                 ! [plant/m2] seedlings of a bare-ground start, per patch
      !--- Mortality hazard gamma + alpha*exp(-beta*growth) (Camac et al. 2018).
      real(wp), allocatable :: mort_gamma(:)     ! [1/yr]  growth-independent baseline hazard
      real(wp), allocatable :: mort_alpha(:)     ! [1/yr]  low-growth hazard magnitude
      real(wp), allocatable :: mort_beta(:)      ! [yr/cm] growth sensitivity of the hazard
      !--- Recruitment.
      real(wp),    allocatable :: seed_rain_recruits(:) ! [plant/m2/yr] baseline external seed rain
      integer(ik), allocatable :: include_pft(:)        ! 1 = PFT may recruit, 0 = excluded
      !--- Shared height thresholds.
      real(wp) :: min_cohort_height          ! [m] smallest tracked cohort; recruits born here
      !--- Leaf photosynthesis and stomata (docs/science/leaf_gas_exchange.md).
      integer(ik), allocatable :: photosynthetic_pathway(:) ! PATH_C3 | PATH_C4
      real(wp), allocatable :: vcmax25(:)            ! [umol/m2/s] max carboxylation at 25 degC
      real(wp), allocatable :: jmax_vcmax_ratio(:)   ! [--]    Jmax25 / Vcmax25
      real(wp), allocatable :: tpu_vcmax_ratio(:)    ! [--]    TPU25 / Vcmax25 (C3 product limit)
      real(wp), allocatable :: rd_vcmax_ratio(:)     ! [--]    Rd25 / Vcmax25 (leaf respiration)
      real(wp), allocatable :: kp25(:)               ! [mol/m2/s] C4 PEPcase initial slope at 25 degC
      real(wp), allocatable :: stomatal_g0(:)        ! [mol/m2/s] cuticular/residual conductance
      real(wp), allocatable :: stomatal_g1(:)        ! Leuning [--] / Medlyn [kPa^0.5] slope
      real(wp), allocatable :: stomatal_d0(:)        ! [Pa]    Leuning humidity sensitivity
      real(wp), allocatable :: quantum_yield_c4(:)   ! [mol CO2/mol photon] C4 light-limited slope
      real(wp), allocatable :: theta_j(:)            ! [--]    non-rectangular hyperbola curvature (C3 J)
      !--- Co-limitation, not theta_j (docs/science/leaf_gas_exchange.md §2).
      real(wp), allocatable :: theta_cj_c3(:)        ! [--]    C3 Ac/Aj co-limitation curvature
      real(wp), allocatable :: theta_ip_c3(:)        ! [--]    C3 (Ac,Aj)/Ap co-limitation curvature
      real(wp), allocatable :: theta_cj_c4(:)        ! [--]    C4 co-limitation curvature 1
      real(wp), allocatable :: theta_ic_c4(:)        ! [--]    C4 co-limitation curvature 2
      real(wp), allocatable :: katul_lambda25(:)     ! [umol CO2/mol H2O] Katul marginal water-use efficiency
      !--- Water stress (docs/science/leaf_gas_exchange.md §4, eqs. (9)-(10)).
      real(wp), allocatable :: wstress_psi_open(:)   ! [MPa] beta_nonstomata starts to fall here; unset: the turgor-loss point
      real(wp), allocatable :: wstress_psi_close(:)  ! [MPa] ... and reaches 0 here; unset: twice the turgor-loss point
      real(wp), allocatable :: wstress_lambda_exp(:) ! [--]    Katul lambda water-stress exponent
      real(wp), allocatable :: wstress_sref_stomata(:) ! [1/MPa] Sabot slope, beta_stomata = exp(sref*(psi_pd - onset))
      real(wp), allocatable :: stomata_psi_onset(:)  ! [MPa] predawn onset of beta_stomata; unset: half the turgor-loss point
      !--- Canopy-film capacities (rain and dew), 0.1 by default (CLM5's dewmx).
      real(wp), allocatable :: leaf_surf_water_max(:)   ! [kg/m2 leaf]
      real(wp), allocatable :: wood_surf_water_max(:)   ! [kg/m2 wood]
      !--- Leaf capacities derived from vcmax25 (derive_leaf_params).
      real(wp), allocatable :: jmax25(:)             ! [umol/m2/s] DERIVED = jmax_vcmax_ratio * vcmax25
      real(wp), allocatable :: tpu25(:)              ! [umol/m2/s] DERIVED = tpu_vcmax_ratio  * vcmax25
      real(wp), allocatable :: rd25(:)               ! [umol/m2/s] DERIVED = rd_vcmax_ratio   * vcmax25
      !--- Carbon allocation, in carbon (docs/science/plant_carbon_allocation.md).
      real(wp),    allocatable :: sla(:)                    ! [m2/kgC] specific leaf area
      !--- Leaf and branch sizes (docs/science/canopy_aerodynamics.md §4).
      real(wp),    allocatable :: leaf_width(:)             ! [m]  characteristic leaf width (leaf boundary layer)
      real(wp),    allocatable :: branch_diameter(:)        ! [m]  characteristic branch diameter (wood boundary layer)
      !--- Allocation targets (docs/science/plant_carbon_allocation.md §2).
      real(wp),    allocatable :: root_to_leaf_ratio(:)     ! [--]     fine-root:leaf target ratio
      real(wp),    allocatable :: storage_cushion(:)        ! [--]     storage target as a multiple of the leaf + fine-root target
      !--- Resorption (docs/science/plant_carbon_allocation.md §1).
      real(wp),    allocatable :: retained_carbon_fraction(:)  ! [-] share of leaf and fine-root loss resorbed
      !--- Sink limits, 0 = off (docs/science/plant_carbon_allocation.md eq. (3)).
      real(wp),    allocatable :: max_relative_growth_rate(:) ! [1/yr] maximum relative growth rate of wood
      real(wp),    allocatable :: max_absolute_growth_rate(:) ! [cm/yr] maximum diameter growth of a 1 cm stem
      real(wp),    allocatable :: max_absolute_growth_exponent(:) ! [--] its scaling with dbh [cm]
      real(wp),    allocatable :: growth_resp_factor(:)     ! [--]     construction cost (fraction of metabolic NPP)
      !--- Stem and root respiration (docs/science/plant_respiration.md §2-§3).
      logical,     allocatable :: is_woody(:)               ! .false. (e.g. grass) => stem respiration is 0
      real(wp),    allocatable :: stem_resp_factor25(:)     ! [umol CO2/m2 stem/s @25C] stem rate per surface
      real(wp),    allocatable :: stem_resp_sapwood25(:)    ! [umol CO2/m3 sapwood/s @25C] stem rate per sapwood volume
      real(wp),    allocatable :: root_resp_factor25(:)     ! [umol CO2/kgC fine root/s @25C] baseline root rate
      !--- Optics, LW as emissivity (docs/science/canopy_radiation_transfer.md §2).
      real(wp),    allocatable :: leaf_reflect_vis(:), leaf_transmit_vis(:)   ! [--] leaf VIS
      real(wp),    allocatable :: leaf_reflect_nir(:), leaf_transmit_nir(:)   ! [--] leaf NIR
      real(wp),    allocatable :: leaf_emissivity(:)                          ! [--] leaf thermal emissivity
      real(wp),    allocatable :: wood_reflect_vis(:), wood_transmit_vis(:)   ! [--] wood/bark VIS
      real(wp),    allocatable :: wood_reflect_nir(:), wood_transmit_nir(:)   ! [--] wood/bark NIR
      real(wp),    allocatable :: wood_emissivity(:)                          ! [--] wood thermal emissivity
      real(wp),    allocatable :: leaf_clumping(:), wood_clumping(:)          ! [--] clumping factor in (0,1]
      real(wp),    allocatable :: leaf_angle_mean(:)                          ! [deg] mean leaf inclination
      real(wp),    allocatable :: leaf_angle_std(:)                           ! [deg] its standard deviation
      !--- Tissue turnover.
      real(wp),    allocatable :: leaf_turnover_rate_toc(:) ! [1/yr]   top-of-canopy baseline leaf turnover
      real(wp),    allocatable :: fineroot_turnover_rate(:) ! [1/yr]   baseline fine-root turnover
      !--- Necromass split into litter pools (docs/science/soil_carbon.md §5).
      real(wp),    allocatable :: f_labile_leaf(:)          ! [--] labile fraction of leaf, storage and fine-root necromass
      real(wp),    allocatable :: f_labile_stem(:)          ! [--] labile (vs structural) fraction of wood/CWD necromass
      real(wp),    allocatable :: struct_lignin_frac(:)     ! [--] lignin fraction of the structural-litter stream
      !--- Light-plasticity slopes (docs/science/plant_traits.md eq. (1)).
      real(wp),    allocatable :: kplastic_sla(:)           ! [1/(m2/m2)] SLA light-response slope (>0)
      real(wp),    allocatable :: kplastic_vcmax(:)         ! [1/(m2/m2)] Vcmax light-response slope (<0)
      real(wp),    allocatable :: kplastic_rd(:)            ! [1/(m2/m2)] Rd light-response slope
      real(wp),    allocatable :: kplastic_leaf_turnover(:) ! [1/(m2/m2)] leaf-turnover light-response slope (<=0)
      type(phenology_traits_t), allocatable :: pheno(:) ! each PFT's phenology ([phenology])
      !--- The tables the fast loop reads, derived from the traits above (derive_pft_tables).
      type(leaf_photo_table_t) :: leaf_photo   ! DERIVED each PFT's leaf-photosynthesis parameters
      type(hydraulics_table_t) :: hydraulics   ! DERIVED each PFT's hydraulics and Kirchhoff lookup
      type(rad_pft_optics_t)   :: optics       ! DERIVED the two-stream's per-PFT optics
      !---~---
      !   A recruit's carbon [kgC/plant] at min_cohort_height on its PFT's allometry, the
      ! value of one plant of the recruit pool; derived on a polygon's own copy of the table,
      ! which carries the cell's allometry (derive_recruit_carbon).
      !---~---
      real(wp), allocatable :: recruit_carbon(:)   ! DERIVED (pft)
   end type pft_table_t

   !---~---
   !   Sentinel of an optional stomatal stress point the PFT file leaves out: large and
   ! negative, so no leaf potential is mistaken for it.
   !---~---
   real(wp), parameter :: HYD_UNSET = -1.0e30_wp

end module ed_pft_params
