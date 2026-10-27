! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_site_state_vars -- the sites of a polygon: each site's area and its slice of the
! patch block.
!
!   A site is the part of a polygon with one soil and one terrain. Its patches are one
! contiguous slice of the polygon's patch block, which is sorted by site and then by age,
! and the site block maps each site to its slice with an offset and a count, as the patch
! block maps each patch to its cohorts. A patch's area is a fraction of its site, and a
! site's area a fraction of the polygon. Each site has its own surface: its soil column and
! the parameters derived from it, the column's bottom, and its terrain; its own climate,
! the polygon's forcing after its terrain (docs/science/forcing.md §8); and its own leaf
! table, acclimated to that climate.
!---~---
module ed_site_state_vars
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : n_soil_layer_max
   use ed_state_defaults,       only : UNSET_POSITIVE
   use ed_diag_vars,            only : site_diag_block
   use ed_config,               only : ed_config_t
   use ed_soil_derive,          only : build_soil_profile, soil_start_theta
   use ed_soil_params,          only : soil_params_t, soil_profile_t
   use ed_pft_params,           only : leaf_photo_table_t
   implicit none
   private

   public :: site_block, site_surface_t, alloc_site_block, free_site_block, make_site_surface

   !---~---
   !   What makes a site: its soil and its terrain. It is configuration held on the site:
   ! never integrated, blended or reordered, and its parameters are derived once, when the
   ! site is made (make_site_surface).
   !---~---
   type :: site_surface_t
      !---~---
      !   The ground: its [soil], with its [[sites]] table's values or its cell's soil maps
      ! in it, and the soil profile derived from it.
      !---~---
      type(soil_params_t)  :: soil
      type(soil_profile_t) :: soil_profile
      real(wp)    :: theta_start(n_soil_layer_max) = 0.0_wp   ! [m3/m3] each layer's starting water
      !--- The terrain: what corrects the polygon's forcing to the site.
      real(wp)    :: elevation_m = 0.0_wp          ! [m]
      real(wp)    :: slope_deg   = 0.0_wp          ! [degrees]
      real(wp)    :: aspect_deg  = 0.0_wp          ! [degrees from north, clockwise]
   end type site_surface_t

   !--- All sites of a polygon: an array per number.
   type :: site_block
      !--- 0. BLOCK RECORD
      integer(ik) :: n   = 0_ik       ! sites in use
      integer(ik) :: cap = 0_ik       ! allocated capacity
      !--- 1. AREA
      real(wp),    allocatable :: area(:)           ! fraction of the polygon, sum = 1
      !--- 2. PATCH MAP: a CSR map into the patch block (rebuild_patch_map derives it)
      integer(ik), allocatable :: patch_offset(:)   ! first patch of the site (1-based)
      integer(ik), allocatable :: patch_count(:)    ! patches in the site
      !--- 3. SURFACE
      type(site_surface_t), allocatable :: surface(:)
      !--- 4. CLIMATE MEMORY
      !---~---
      !   GROWTH TEMPERATURE for thermal acclimation (docs/science/leaf_gas_exchange.md eq.
      ! (4)): an exponential running mean of the site's daily-mean air temperature, window
      ! [leaf_physiology].acclim_window_days (~30 d), which is what Kattge & Knorr (2007)
      ! fitted dS against. PROGNOSTIC and slow, so it is written to the state file:
      ! rebuilding a month of memory after every restart would make the acclimated optimum
      ! jump. Until it is seeded, the first slow step adopts that day's mean outright
      ! rather than relaxing from an arbitrary origin.
      !---~---
      real(wp),    allocatable :: t_growth_avg(:)   ! [K] running-mean growth temperature
      !---~---
      !   Each site's leaf table: the PFT table's, its entropy terms and Jmax:Vcmax ratio
      ! acclimated to the site's growth temperature (acclimate_site_leaf_tables). Derived
      ! from t_growth_avg, so never restarted.
      !---~---
      type(leaf_photo_table_t), allocatable :: leaf_photo(:)
      !--- 5. FAST-TO-SLOW ACCUMULATORS (the model reads them)
      !---~---
      !   The daily-mean air-temperature accumulator for phenology: the fast loop sums the
      ! site's air temperature over the slow step (reset each step), and the slow-loop
      ! phenology and growth temperature read the mean (sum/n) then. Transient (never
      ! restarted).
      !---~---
      real(wp),    allocatable :: pheno_tair_sum(:)
      integer(ik), allocatable :: pheno_tair_n(:)
      !--- 6. OUTPUT ACCUMULATORS (only output reads them; never restarted)
      type(site_diag_block) :: diag                 ! the forcing each site receives (SI_*)
   end type site_block

contains

   !--- THE SITE LOCKSTEP: the routines that list every site field, in section order

   subroutine alloc_site_block(site, cap, leaf_photo)
      type(site_block),         intent(inout) :: site
      integer(ik),              intent(in)    :: cap
      type(leaf_photo_table_t), intent(in)    :: leaf_photo   ! the PFT table's, each site's until it acclimates
      site%cap = cap
      site%n = 0_ik
      !--- 1. Area
      allocate(site%area(cap))
      site%area = 0.0_wp
      !--- 2. Patch map
      allocate(site%patch_offset(cap), site%patch_count(cap))
      site%patch_offset = 0_ik
      site%patch_count = 0_ik
      !--- 3. Surface
      allocate(site%surface(cap))                                ! default-initialised by its type
      !--- 4. Climate memory
      allocate(site%t_growth_avg(cap))
      site%t_growth_avg = UNSET_POSITIVE
      allocate(site%leaf_photo(cap), source=leaf_photo)
      !--- 5. Fast-to-slow accumulators
      allocate(site%pheno_tair_sum(cap), site%pheno_tair_n(cap))
      site%pheno_tair_sum = 0.0_wp
      site%pheno_tair_n = 0_ik
   end subroutine alloc_site_block

   !--- Free every per-site array.
   subroutine free_site_block(site)
      type(site_block), intent(inout) :: site
      site%n = 0_ik
      site%cap = 0_ik
      if (allocated(site%area)) deallocate(site%area, site%patch_offset, site%patch_count, site%surface,     &
                                           site%t_growth_avg, site%leaf_photo, site%pheno_tair_sum,          &
                                           site%pheno_tair_n)
      if (allocated(site%diag%v)) deallocate(site%diag%v, site%diag%w)
      site%diag%active = .false.
   end subroutine free_site_block

   !--- A SITE'S SURFACE

   !---~---
   !   A site's surface from a configuration: its ground, from [soil] or the soil maps of
   ! its cell (apply_map_inputs), and the profile derived from it. The terrain is the
   ! caller's: it knows where the site is.
   !---~---
   subroutine make_site_surface(surf, cfg)
      type(site_surface_t), intent(inout) :: surf
      type(ed_config_t),  intent(in)      :: cfg
      surf%soil = cfg%soil
      call build_soil_profile(cfg%soil, surf%soil_profile)
      call soil_start_theta(cfg%soil, surf%soil_profile, cfg%init%soil_theta, surf%theta_start)
   end subroutine make_site_surface

end module ed_site_state_vars
