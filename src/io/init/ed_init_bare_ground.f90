! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_init_bare_ground -- init_mode 0, a near-bare-ground start, and the builders every     !
! start that makes cohorts shares (docs/science/order_of_processes.md §1, setup step 2;    !
! docs/configuration.md, How a run starts). init_bare_ground (one site) and                !
! init_bare_sites (several) lay out empty patches and seed their stores; seed_bare_ground  !
! plants each PFT's nbg_init_density; add_cohort appends a cohort and finalize_init maps   !
! and sorts them.                                                                          !
!==========================================================================================!
module ed_init_bare_ground
   use ed_kinds,        only : wp, ik
   use ed_config,       only : ed_config_t, growth_window_steps
   use ed_therm_lib,    only : cas_temp_to_enthalpy
   use ed_soil_heat,    only : soil_layer_energy, soil_layer_temp
   use ed_met_types,    only : met_forcing_t
   use ed_allometry_lib, only : height_to_dbh
   use ed_cohort_state_vars,       only : expand_cohort_block, init_cohort
   use ed_site_state_vars,         only : make_site_surface
   use ed_polygon_state_vars,      only : polygon_t, alloc_polygon, rebuild_patch_map, rebuild_cohort_map,      &
                                          assign_cohort_id, assign_patch_id
   use ed_patch_state_vars,        only : patch_physical_state, DIST_PRIMARY
   use ed_slow_update, only : sort_cohorts
   implicit none
   private

   public :: init_bare_ground, init_bare_sites, seed_bare_ground, add_cohort, finalize_init, seed_patch_stores
   public :: N_PATCH_INIT

   integer(ik), parameter :: N_PATCH_INIT = 6_ik     !< a bare-ground start's patches per site

contains

   !---------------------------------------------------------------------------------------!
   ! A bare polygon of one site: n_patch identical empty patches sharing the site's area.  !
   !---------------------------------------------------------------------------------------!
   subroutine init_bare_ground(poly, cfg, n_patch)
      type(polygon_t),     intent(out) :: poly
      type(ed_config_t),   intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch
      call init_bare_sites(poly, cfg, [n_patch], [cfg], [1.0_wp])
   end subroutine init_bare_ground

   !-------------------------------------------------------------------------------------------!
   ! A near-bare-ground polygon of several sites: site `is` takes its share `site_area(is)` of !
   ! the polygon, its surface from `site_cfg(is)` (its soil, from [soil] or its own maps) and  !
   ! n_patch(is) identical empty patches sharing its area. The patches are numbered site by    !
   ! site, so the block is sorted by site, and their stores are seeded (seed_patch_stores), so  !
   ! every cohort set up on them is born in its patch's state. Everything else comes from      !
   ! `cfg`, the run's.                                                                         !
   !-------------------------------------------------------------------------------------------!
   subroutine init_bare_sites(poly, cfg, n_patch, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(ed_config_t),   intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch(:)     !< patches per site
      type(ed_config_t),   intent(in)  :: site_cfg(:)
      real(wp),            intent(in)  :: site_area(:)   !< [-] fractions of the polygon
      integer(ik) :: ip, is, k, ns
      ns = int(size(n_patch), ik)
      call alloc_polygon(poly, cfg%pft%n, coh_cap = 64_ik, pat_cap = max(sum(n_patch), 1_ik),  &
                      n_growth_window = growth_window_steps(cfg), pft = cfg%pft,                    &
                      n_site = ns)
      poly%patch%n = sum(n_patch)
      ip = 0_ik
      do is = 1_ik, ns
         call make_site_surface(poly%site%surface(is), site_cfg(is))
         poly%site%area(is) = site_area(is)
         do k = 1_ik, n_patch(is)
            ip = ip + 1_ik
            poly%patch%area(ip)       = 1.0_wp / real(n_patch(is), wp)
            poly%patch%age(ip)        = 0.0_wp
            poly%patch%dist_type(ip)  = DIST_PRIMARY
            poly%patch%owner_site(ip) = is
            call assign_patch_id(poly, ip)
         end do
      end do
      poly%patch%recruit_pool = 0.0_wp
      poly%cohort%n = 0_ik
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
      call seed_patch_stores(poly, cfg)
   end subroutine init_bare_sites

   !----- A start's patch stores, the same on every patch of a site: the soil at its site's        !
   !      theta_start and [init].soil_temp, its water liquid and ice as the layer's freezing curve !
   !      holds them there (soil_layer_energy); the pond dry; and the canopy air at the reference  !
   !      climate (met_forcing_t's air, at the run's CO2), which the fast loop's first window      !
   !      replaces with the first forcing sample (docs/science/canopy_air_space_biophysics.md, The !
   !      start). ---------------------------------------------------------------------------------!
   subroutine seed_patch_stores(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      type(met_forcing_t) :: ref   !< the reference climate
      integer(ik) :: ip, k, nsl
      ref = met_forcing_t()
      do ip = 1_ik, poly%patch%n
         associate (cas => poly%patch%cas(ip), se => poly%patch%soil_energy(ip), sw => poly%patch%soil_water(ip), &
                    surf => poly%site%surface(poly%patch%owner_site(ip)))
            nsl = surf%soil_profile%n_active
            sw%theta(1:nsl)     = surf%theta_start(1:nsl)
            sw%w_surface        = 0.0_wp
            sw%w_surface_energy = 0.0_wp        ! dry pond -> zero enthalpy
            do k = 1_ik, nsl
               se%soil_energy(k) = soil_layer_energy(cfg%init%soil_temp, surf%theta_start(k),                  &
                                                     surf%soil_profile%dry_heat_capacity(k),                    &
                                                     surf%soil_profile%layer_curve(k))
               call soil_layer_temp(se%soil_energy(k), surf%theta_start(k), surf%soil_profile%dry_heat_capacity(k), &
                                    surf%soil_profile%layer_curve(k), se%soil_temp(k), se%soil_fliq(k))
            end do
            cas%can_shv      = ref%air_shv
            cas%can_co2      = cfg%forcing%co2_const
            cas%can_enthalpy = cas_temp_to_enthalpy(ref%air_temp, ref%air_shv)
            cas%can_temp     = ref%air_temp
         end associate
      end do
   end subroutine seed_patch_stores

   !---------------------------------------------------------------------------------------!
   ! A near-bare-ground start's seedlings: every patch gets one cohort of each PFT with    !
   ! nbg_init_density > 0, that many plants per m2, at min_cohort_height on the PFT's own  !
   ! height curve, where recruits are born. All densities 0 leave the patches empty.       !
   !---------------------------------------------------------------------------------------!
   subroutine seed_bare_ground(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      real(wp)    :: dbh
      integer(ik) :: ip, ipft
      if (.not. any(cfg%pft%nbg_init_density > 0.0_wp)) return
      do ip = 1_ik, poly%patch%n
         do ipft = 1_ik, cfg%pft%n
            if (cfg%pft%nbg_init_density(ipft) <= 0.0_wp) cycle
            dbh = height_to_dbh(poly%cohort%pft_table%allom(ipft), cfg%pft%min_cohort_height)
            call add_cohort(poly, cfg, ip, ipft, cfg%pft%nbg_init_density(ipft), dbh)
         end do
      end do
      call finalize_init(poly)
   end subroutine seed_bare_ground

   !---------------------------------------------------------------------------------------!
   ! Append one cohort to patch ip (caller invokes finalize_init afterwards).               !
   !---------------------------------------------------------------------------------------!
   subroutine add_cohort(poly, cfg, ip, ipft, nplant, dbh)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      integer(ik),         intent(in)    :: ip, ipft
      real(wp),            intent(in)    :: nplant, dbh
      integer(ik) :: m
      real(wp)    :: t_birth, psi_birth
      call expand_cohort_block(poly%cohort, poly%cohort%n + 1_ik)
      m = poly%cohort%n + 1_ik
      !----- Born in its patch's state (patch_physical_state), which init_bare_sites seeded. ----------!
      call patch_physical_state(poly%patch%cas(ip), poly%patch%soil_water(ip),                         &
                                poly%site%surface(poly%patch%owner_site(ip)), t_birth, psi_birth)
      call init_cohort(poly%cohort, m, cfg%pft, ipft, ip, nplant, dbh, t_birth, psi_birth)
      call assign_cohort_id(poly, m)
      poly%cohort%n = m
   end subroutine add_cohort

   !----- The cohorts' map onto their patches, then each patch's cohorts sorted by height. --!
   subroutine finalize_init(poly)
      type(polygon_t), intent(inout) :: poly
      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine finalize_init

end module ed_init_bare_ground
