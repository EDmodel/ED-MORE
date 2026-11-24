! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_ground_surface_params -- the [ground_surface_layer] block: the layers on the soil, top
! down the water layer (snow, ice and ponded water), the litter and the moss, each a material
! with its own table (docs/science/snow_biophysics.md, docs/science/surface_organic_layer.md).
! Types and codes only: a config leaf the kernels import without ed_config.
!---~---
module ed_ground_surface_params
   use ed_kinds,         only : wp, ik
   use ed_soil_hydr_lib, only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   implicit none
   private

   public :: water_material_t, porous_material_t, ground_surface_layer_params_t
   public :: POROUS_LITTER, POROUS_MOSS, N_POROUS
   public :: MASS_FROM_NONE, MASS_FROM_METABOLIC_LITTER, MASS_FROM_SITE_DEPTH

   !--- The porous layers' slots, top down: litter over moss.
   integer(ik), parameter :: POROUS_LITTER = 1_ik, POROUS_MOSS = 2_ik, N_POROUS = 2_ik

   !---~---
   !   Where a porous layer's dry matter comes from (its mass_source): nowhere (no layer), the
   ! above-ground metabolic litter pool, or a site's [soil].moss_depth at the material's bulk
   ! density.
   !---~---
   integer(ik), parameter :: MASS_FROM_NONE = 0_ik, MASS_FROM_METABOLIC_LITTER = 1_ik, MASS_FROM_SITE_DEPTH = 2_ik

   !---~---
   !   [ground_surface_layer.water]: the water layer, snow, ice and the liquid they hold, one
   ! layer. Its ice matrix conducts at conductivity_at_fresh at the fresh density, scaled to its own
   ! density by Calonne et al. (2011)'s law; its albedo blends snow's, solid ice's and standing
   ! water's (CLM5's lake ice and open water). Its free water beyond what its ice holds is the
   ! pond, which runs off above pond_capacity.
   !---~---
   type :: water_material_t
      real(wp) :: fresh_density          = 250.0_wp    ! [kg/m3] snowfall's bulk density
      real(wp) :: conductivity_at_fresh  = 0.15_wp     ! [W/m/K] the matrix's thermal conductivity at fresh_density
      real(wp) :: albedo_vis_fresh       = 0.85_wp     ! [-] fresh snow, visible
      real(wp) :: albedo_nir_fresh       = 0.65_wp     ! [-] fresh snow, near infrared
      real(wp) :: albedo_vis_aged        = 0.518_wp    ! [-] wet (aged) snow, visible
      real(wp) :: albedo_nir_aged        = 0.435_wp    ! [-] wet (aged) snow, near infrared
      real(wp) :: albedo_vis_ice         = 0.60_wp     ! [-] solid ice, visible
      real(wp) :: albedo_nir_ice         = 0.40_wp     ! [-] solid ice, near infrared
      real(wp) :: albedo_vis_liquid      = 0.10_wp     ! [-] standing water, visible
      real(wp) :: albedo_nir_liquid      = 0.10_wp     ! [-] standing water, near infrared
      real(wp) :: emissivity             = 0.97_wp     ! [-] longwave emissivity
      real(wp) :: liquid_holding_frac    = 0.10_wp     ! [-] the layer's mass share its ice holds as liquid
      real(wp) :: cover_scale_factor     = 2.5_wp      ! [-] Niu & Yang (2007)'s cover scale over the roughness
      real(wp) :: cover_density_exponent = 1.0_wp      ! [-] its growth with the layer's density
      real(wp) :: pond_capacity          = 5.0_wp      ! [kg/m2] the free water above which it runs off
      real(wp) :: initial_mass           = 0.0_wp      ! [kg/m2] the layer at a fresh start
      real(wp) :: initial_temp           = 270.0_wp    ! [K]     its temperature
   end type water_material_t

   !---~---
   !   [ground_surface_layer.litter], [ground_surface_layer.moss]: a porous layer's material, its
   ! dry mass per volume, its thermal conductivity dry and saturated, and its water retention and
   ! conductivity as a soil layer's (ed_soil_hydr_lib): the curve's form, the water at
   ! saturation and at residual, its two parameters and the saturated conductivity; and where its
   ! matter comes from.
   !---~---
   type :: porous_material_t
      integer(ik) :: mass_source  = MASS_FROM_NONE     ! MASS_FROM_*
      real(wp)    :: bulk_density = 0.0_wp             ! [kg/m3] dry mass per volume
      real(wp)    :: k_dry        = 0.0_wp             ! [W/m/K] thermal conductivity, dry
      real(wp)    :: k_wet        = 0.0_wp             ! [W/m/K] thermal conductivity, saturated
      integer(ik) :: retention    = SOIL_RETENTION_VG  ! SOIL_RETENTION_VG or SOIL_RETENTION_CAMPBELL
      real(wp)    :: theta_sat    = 0.0_wp             ! [m3/m3] water at saturation
      real(wp)    :: theta_res    = 0.0_wp             ! [m3/m3] residual water
      real(wp)    :: curve_par_a  = 0.0_wp             ! vG alpha [1/m], or Campbell's air-entry head [m, < 0]
      real(wp)    :: curve_par_n  = 0.0_wp             ! vG n [-], or Campbell's b [-]
      real(wp)    :: ksat         = 0.0_wp             ! [m/s] saturated conductivity
   end type porous_material_t

   !---~---
   !   The [ground_surface_layer] block. A layer is a layer of the ground from min_depth: thinner
   ! matter stays where it came from. The litter is the above-ground metabolic litter pool's dry
   ! matter (Napoly et al. 2017's forest-floor values), retaining and conducting water as Letts et
   ! al. (2000)'s surface moss, the curve CLASS gives its litter and moss (Wu et al. 2016); the
   ! moss is a site's [soil].moss_depth, at the MAT tundra's top organic 5 cm (ARC LTER
   ! knb-lter-arc.20129.1), Chadburn et al. (2015)'s conductivities, Zhang et al. (2010)'s
   ! retention of surface moss and O'Connor et al. (2020)'s saturated conductivity of the Arctic
   ! acrotelm.
   !---~---
   type :: ground_surface_layer_params_t
      real(wp) :: min_depth = 1.0e-4_wp                ! [m]
      type(water_material_t)  :: water
      type(porous_material_t) :: porous(N_POROUS) = [                                                      &
         porous_material_t(mass_source=MASS_FROM_METABOLIC_LITTER, bulk_density=45.0_wp, k_dry=0.10_wp,        &
                           k_wet=0.1294_wp, retention=SOIL_RETENTION_CAMPBELL, theta_sat=0.98_wp,              &
                           theta_res=0.01_wp, curve_par_a=-0.0103_wp, curve_par_n=2.3_wp, ksat=1.83e-3_wp),    &
         porous_material_t(mass_source=MASS_FROM_SITE_DEPTH, bulk_density=43.0_wp, k_dry=0.06_wp,             &
                           k_wet=0.50_wp, retention=SOIL_RETENTION_VG, theta_sat=0.85_wp, theta_res=0.18_wp,   &
                           curve_par_a=50.0_wp, curve_par_n=1.45_wp, ksat=1.3e-3_wp)]
   end type ground_surface_layer_params_t

end module ed_ground_surface_params
