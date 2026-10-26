! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_params -- the [soil] block, the ground a site is set on (soil_params_t: its    !
! column, its surface and its two bottom boundaries), and the profile the soil kernels     !
! read, soil_profile_t (the layers, their hydraulic and thermal texture, the field         !
! capacity and the water network's layer curves), which meds_soil_derive builds from it    !
! once per site (build_soil_profile). A config leaf (src/share only): the soil kernels     !
! import it without meds_config.                                                           !
!==========================================================================================!
module meds_soil_params
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : n_soil_layer_max
   use meds_soil_hydr_lib,    only : SOIL_RETENTION_VG, soil_layer_curve_t
   implicit none
   private

   public :: soil_params_t, soil_texture_layers_t, soil_profile_t
   public :: SOIL_BC_DRAIN, SOIL_BC_AQUIFER, ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET
   public :: ANNUAL_DAMPING_DEPTH

   !----- Soil-water bottom boundary: drainage at a fraction f of the deepest layer's      !
   !      conductivity, q = f K(theta_N), or an aquifer (docs/science/soil_biophysics.md,  !
   !      "Boundary conditions"). ---------------------------------------------------------!
   integer(ik), parameter :: SOIL_BC_DRAIN   = 1_ik     !< drainage q = f K(theta_N)  (default)
   integer(ik), parameter :: SOIL_BC_AQUIFER = 2_ik     !< water table at the column base (wet sites only)

   !----- Soil-energy bottom boundary (docs/science/soil_biophysics.md). ------------------!
   integer(ik), parameter :: ENERGY_BC_GEOTHERMAL   = 1_ik   !< bottom: prescribed flux, held at zero (Neumann)
   integer(ik), parameter :: ENERGY_BC_DIRICHLET     = 2_ik   !< bottom: conduction to a fixed deep temperature
   !----- Annual damping depth d = sqrt(2*alpha/omega) of a mid-latitude mineral soil at   !
   !      theta = 0.3. The Dirichlet anchor reflects the annual wave least d/sqrt(2) below !
   !      the bottom node (docs/science/soil_biophysics.md eq. (13)): deep_depth's         !
   !      default, and a soil map's column. -----------------------------------------------!
   real(wp), parameter :: ANNUAL_DAMPING_DEPTH = 1.973_wp   !< [m]

   !----- [soil] (docs/science/soil_biophysics.md), each key optional at the default below: !
   !      the ground a site is set on. Each site holds its own (site_surface_t%soil), with  !
   !      its [[sites]] table's values or its cell's soil maps in it. ----------------------!
   type :: soil_params_t
      !----- Layers: layer_bottoms or SOIL_LAYER_TEMPLATE, cut at depth. ------------------!
      real(wp)    :: depth           = 2.0_wp    !< [m]  total column depth (positive; z is negative down)
      integer(ik) :: n_layer_bottoms = 0_ik      !< [-]  entries of layer_bottoms; 0 = the template
      real(wp)    :: layer_bottoms(n_soil_layer_max) = 0.0_wp   !< [m] layer bottoms, positive down
      !----- Hydraulic texture, one for the column (or per layer, below). -----------------!
      integer(ik) :: retention   = SOIL_RETENTION_VG  !< van_genuchten | campbell
      real(wp)    :: theta_sat   = 0.43_wp    !< [m3/m3] porosity
      real(wp)    :: theta_res   = 0.078_wp   !< [m3/m3] residual water content
      real(wp)    :: ksat        = 2.89e-6_wp !< [m/s]   saturated hydraulic conductivity
      real(wp)    :: curve_par_a = 3.6_wp     !< [1/m] van Genuchten alpha, or [m] Campbell psi_sat (< 0)
      real(wp)    :: curve_par_n = 1.56_wp    !< [-]   van Genuchten n (>1), or [-] Campbell b
      real(wp)    :: psi_fc      = -3.37_wp   !< [m]   field-capacity matric head (derives theta_fc)
      !----- Thermal texture, one for the column (or per layer, below). -------------------!
      real(wp)    :: solid_conductivity = 3.0_wp    !< [W/m/K]  mineral-solid conductivity
      real(wp)    :: dry_conductivity   = 0.15_wp   !< [W/m/K]  dry-matrix conductivity
      real(wp)    :: dry_heat_capacity  = 2.0e6_wp  !< [J/m3/K] dry-matrix volumetric heat capacity
      !----- The surface: the pond, the dry surface layer (soil evaporation), and the bare-soil !
      !      optics, blended with snow where there is snow. ---------------------------------!
      real(wp)    :: w_pond_max        = 5.0_wp     !< [kg/m2] ponding capacity before surface runoff
      real(wp)    :: dsl_dmax          = 0.015_wp   !< [m] DSL max thickness (soil evap)
      real(wp)    :: dsl_theta_init    = 0.8_wp     !< DSL initiation theta_1/theta_sat
      real(wp)    :: ground_albedo_vis = 0.15_wp    !< [-] soil albedo, visible
      real(wp)    :: ground_albedo_nir = 0.30_wp    !< [-] soil albedo, near infrared
      real(wp)    :: ground_emissivity = 0.95_wp    !< [-] soil longwave emissivity
      !----- The water bottom boundary, with drainage's factor f of q = f K. ---------------!
      integer(ik) :: water_bottom_bc = SOIL_BC_DRAIN    !< SOIL_BC_DRAIN | SOIL_BC_AQUIFER
      real(wp)    :: drain_factor    = 1.0_wp           !< [-] f in [0, 1]: free_drain 1, bedrock 0
      !----- The thermal bottom boundary (docs/science/soil_biophysics.md eqs. (12)-(13)). The  !
      !      backward-Euler soil heat solve is unconditionally stable and the outer march      !
      !      controls its accuracy, so there are no sub-step controls. ----------------------!
      integer(ik) :: thermal_bottom_bc = ENERGY_BC_GEOTHERMAL   !< ENERGY_BC_GEOTHERMAL | ENERGY_BC_DIRICHLET
      real(wp)    :: deep_temp   = -1.0_wp        !< [K] anchor temperature; < 0: not set (dirichlet requires it)
      real(wp)    :: deep_depth  = 3.145_wp       !< [m] anchor plane depth: the default column's optimum, eq. (13)
      !----- Per-layer texture from the soil maps (meds_maps); no config key. -------------!
      integer(ik) :: n_layer_values = 0_ik      !< [-] the column's layer count; 0: the one texture above
      real(wp)    :: layer_theta_sat(n_soil_layer_max)          = 0.0_wp
      real(wp)    :: layer_theta_res(n_soil_layer_max)          = 0.0_wp
      real(wp)    :: layer_ksat(n_soil_layer_max)               = 0.0_wp
      real(wp)    :: layer_curve_par_a(n_soil_layer_max)        = 0.0_wp
      real(wp)    :: layer_curve_par_n(n_soil_layer_max)        = 0.0_wp
      real(wp)    :: layer_solid_conductivity(n_soil_layer_max) = 0.0_wp
      real(wp)    :: layer_dry_conductivity(n_soil_layer_max)   = 0.0_wp
      real(wp)    :: layer_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp
   end type soil_params_t

   !----- A profile's texture, one value per layer: the per-layer arrays, or the one       !
   !      texture repeated (soil_layer_textures). -----------------------------------------!
   type :: soil_texture_layers_t
      real(wp) :: theta_sat(n_soil_layer_max), theta_res(n_soil_layer_max), ksat(n_soil_layer_max)
      real(wp) :: par_a(n_soil_layer_max), par_n(n_soil_layer_max)
      real(wp) :: k_solid(n_soil_layer_max), k_dry(n_soil_layer_max), dry_cvol(n_soil_layer_max)
   end type soil_texture_layers_t

   !----- A site's soil profile, derived once from its [soil] (or its cell's soil maps): the  !
   !      layers' geometry, elevations z <= 0 below ground, dz and dz_node > 0, and each     !
   !      layer's hydraulic and thermal texture, the soil water and heat solvers' input.     !
   !      Fixed size (n_soil_layer_max), so the soil kernels stay allocatable-free and       !
   !      GPU-eligible. -------------------------------------------------------------------!
   type :: soil_profile_t
      integer(ik) :: n_active  = n_soil_layer_max            !< active layer count (<= n_soil_layer_max)
      integer(ik) :: retention = SOIL_RETENTION_VG           !< curve family
      real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp  !< [m] interface elevations (<= 0)
      real(wp) :: z_node(n_soil_layer_max)  = 0.0_wp         !< [m] node (mid) elevations (<= 0)
      real(wp) :: dz(n_soil_layer_max)      = 0.0_wp         !< [m] layer thickness (> 0)
      real(wp) :: dz_node(n_soil_layer_max) = 0.0_wp         !< [m] internode spacing (> 0)
      real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp       !< [m3/m3] porosity
      real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp       !< [m3/m3] residual water content
      real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp       !< [m/s] saturated conductivity
      !----- One pair for both curve families (docs/science/soil_biophysics.md). ----------!
      real(wp) :: curve_par_a(n_soil_layer_max) = 0.0_wp     !< [1/m] VG alpha, or [m] Campbell psi_sat
      real(wp) :: curve_par_n(n_soil_layer_max) = 0.0_wp     !< [-] VG n, or Campbell b
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
      !----- What the soil-plant water network reads (DERIVED): each layer's curve in the   !
      !      network's unknown, and the lower layer's weight in each face's conductivity. --!
      type(soil_layer_curve_t) :: layer_curve(n_soil_layer_max)
      real(wp) :: face_weight(n_soil_layer_max) = 0.0_wp     !< [-] the face below layer k
      !----- Thermal texture. ---------------------------------------------------------------!
      real(wp) :: solid_conductivity(n_soil_layer_max) = 0.0_wp   !< [W/m/K] kappa_solid
      real(wp) :: dry_conductivity(n_soil_layer_max)   = 0.0_wp   !< [W/m/K] kappa_dry
      real(wp) :: dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   !< [J/m3/K] dry-matrix vol. heat cap
   end type soil_profile_t

end module meds_soil_params
