! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_init -- the start of a polygon: init_polygon runs the start [init].init_mode names    !
! (bare ground, a census or a restart file; docs/configuration.md, How a run starts) and   !
! that start's own finishing: a census stand is restructured and its plastic traits        !
! acclimated at once, a restart's traits re-acclimated when [init].reacclimate_traits asks !
! (docs/science/order_of_processes.md §1, setup step 3). Bare ground and a census then     !
! settle their cohorts in their patches (init_fast_state); a restart read them from its    !
! file. Every start ends with its sites' terrain, the canopy-air depth, and the            !
! soil-carbon spin-up a new stand asks for, so the polygon's state is complete when        !
! init_polygon returns.                                                                    !
!==========================================================================================!
module ed_init
   use ed_kinds,                 only : wp, ik
   use ed_time,                  only : ed_time_t, time_to_string
   use ed_config,                only : ed_config_t, INIT_CENSUS, INIT_RESTART
   use ed_soil_derive,           only : top_soil_psi
   use ed_therm_lib,             only : temp_to_internal_energy, internal_energy_to_temp
   use ed_cohort_state_vars,     only : set_cohort_tissue_water, hold_cohort_tissue_temps
   use ed_polygon_state_vars,    only : polygon_t
   use ed_slow_update,           only : update_overtopping_lai, max_cohort_count, commit_trait_carbon,  &
                                        refresh_canopy_depth
   use ed_slow_ledger,           only : slow_ledger_t
   use ed_soil_carbon,           only : soil_carbon_steady_pools
   use ed_patch_state_vars,      only : soil_carbon_state_t
   use ed_slow_state,            only : advance_plant_traits, reacclimate_plant_traits
   use ed_init_bare_ground,      only : init_bare_ground, init_bare_sites, seed_bare_ground, seed_patch_stores, &
                                         N_PATCH_INIT
   use ed_init_census,           only : init_from_census, restructure_census_stand
   use ed_init_restart,          only : read_restart_file
   implicit none
   private

   public :: init_polygon
   public :: init_fast_state   ! used only by the tests

contains

   !---------------------------------------------------------------------------------------!
   ! Start `poly` as cfg%init%init_mode says, at its place (polygon_id and location, set   !
   ! once the start has allocated it). A restart reads row `restart_row` of the restart    !
   ! file and sets `now` to the file's time; a census reads the census file's rows for the !
   ! polygon's site(s); bare ground lays out N_PATCH_INIT patches per site and seeds them. !
   ! Then each site's terrain, the canopy-air depth and the soil-carbon spin-up. site_cfg  !
   ! and site_area give each site of a polygon of several. ok = .false. (with the reason   !
   ! printed, naming `label`) when the start's file cannot serve it.                       !
   !---------------------------------------------------------------------------------------!
   subroutine init_polygon(poly, cfg, label, verbose, restart_row, now, ok, polygon_id, latitude_deg,     &
                           longitude_deg, elevation_m, site_cfg, site_area)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      character(len=*),    intent(in)    :: label        !< the polygon's name, for messages
      integer(ik),         intent(in)    :: polygon_id   !< the polygon's id: its cell's index, 0 for a site
      real(wp),            intent(in)    :: latitude_deg, longitude_deg, elevation_m   !< its location
      logical,             intent(in)    :: verbose      !< report the start
      integer(ik),         intent(in)    :: restart_row  !< the polygon's row of the restart file
      type(ed_time_t),     intent(inout) :: now          !< the run's start; a restart's time on return
      logical,             intent(out)   :: ok
      type(ed_config_t),   intent(in), optional :: site_cfg(:)   !< each site's configuration ...
      real(wp),            intent(in), optional :: site_area(:)  !< ... and share of the polygon
      character(len=:), allocatable :: msg
      type(ed_time_t)   :: restart_time
      logical           :: found
      real(wp), allocatable :: leaf_carbon_new(:), storage_new(:)   !< [kgC/plant] the traits' leaf and storage
      type(slow_ledger_t)       :: no_ledger   !< inactive: the start books nothing
      type(soil_carbon_state_t) :: pools

      ok = .false.
      select case (cfg%init%init_mode)
      case (INIT_RESTART)
         call read_restart_file(poly, cfg, trim(cfg%init%restart_file), restart_time, found, restart_row,   &
                                site_cfg=site_cfg, why=msg)
         if (.not. found) then
            write(*,'(7a)') ' init  : restart (mode 2): cannot read ', trim(label),                       &
                            ' from the state file ', trim(cfg%init%restart_file), ', which ', msg
            return
         end if
         now = restart_time
         if (verbose) call report_restart(cfg, now, poly%restructure_pending)
      case (INIT_CENSUS)
         if (present(site_cfg)) then
            call init_from_census(poly, cfg, trim(cfg%init%census_file), found, why=msg,                   &
                                  site_ids=cfg%sites%census_site_id, site_cfg=site_cfg, site_area=site_area)
         else
            call init_from_census(poly, cfg, trim(cfg%init%census_file), found, why=msg)
         end if
         if (.not. found) then
            write(*,'(4a)') ' init  : census (mode 1): the census file ', trim(cfg%init%census_file), ' ', msg
            return
         end if
         if (verbose) then
            write(*,'(2a)') ' init  : census (mode 1) ', trim(cfg%init%census_file)
            write(*,'(a,i0,a,i0,a)') '         read ', poly%patch%n, ' patches, ', poly%cohort%n, ' cohorts'
         end if
         !----- A census stand is restructured by the slow step's own operators before the first  !
         !      step: it arrives with a cohort per measured size and a patch per plot cell.  -------!
         call restructure_census_stand(poly, cfg)
         if (verbose) write(*,'(a,i0,a,i0,a,i0,a)') '         restructured to ', poly%patch%n,          &
                                                     ' patches, ', poly%cohort%n, ' cohorts (at most ',   &
                                                     max_cohort_count(poly), ' in a patch)'
      case default
         if (present(site_cfg)) then
            call init_bare_sites(poly, cfg, spread(N_PATCH_INIT, 1, size(site_cfg)), site_cfg, site_area)
         else
            call init_bare_ground(poly, cfg, N_PATCH_INIT)
         end if
         call seed_bare_ground(poly, cfg)
         if (verbose) write(*,'(a)') ' init  : bare ground (mode 0)'
      end select
      !----- The polygon's place, on the stand the start built (allocating it reset the polygon). -----!
      poly%id            = polygon_id
      poly%latitude_deg  = latitude_deg
      poly%longitude_deg = longitude_deg
      poly%elevation_m   = elevation_m

      !----- Census cohorts sit in an established stand but carry no trait history, so with      !
      !      plasticity on their leaf traits acclimate to their light at once. Bare ground starts !
      !      at top of canopy; a restart keeps the traits it read, unless                         !
      !      [init].reacclimate_traits asks for this run's PFT file instead.                      !
      if (cfg%trait_dynamics%trait_plasticity_on .and. cfg%init%init_mode == INIT_CENSUS) then
         call update_overtopping_lai(poly)
         call advance_plant_traits(poly, cfg%dt_years, leaf_carbon_new, storage_new, instantaneous=.true.)
         call commit_trait_carbon(poly, leaf_carbon_new, storage_new)
      end if
      if (cfg%init%reacclimate_traits .and. cfg%init%init_mode == INIT_RESTART) then
         call reacclimate_plant_traits(poly, cfg, leaf_carbon_new, storage_new)
         call commit_trait_carbon(poly, leaf_carbon_new, storage_new)
         if (verbose) write(*,'(a)') ' init  : plastic traits re-acclimated to this run''s PFT file'
      end if
      !----- The final stand's cohorts in their patches, and the snowpack, with the fast loop on or   !
      !      off, and the canopy air owed the first forcing sample; a restart read them instead. -----!
      if (cfg%init%init_mode /= INIT_RESTART) then
         call init_fast_state(poly, cfg)
         poly%seed_canopy_air = .true.
         if (verbose .and. cfg%snow_init_swe > 0.0_wp)                                             &
            write(*,'(a,f6.1,a)') ' snow  : seeded initial pack SWE = ', cfg%snow_init_swe, ' kg/m2'
      end if
      !----- What every start ends with: each site's terrain, the canopy-air depth the stand sets (a   !
      !      plain geometry update, so no ledger books it), and, unless a restart restored them, the  !
      !      soil-carbon pools at their spin-up steady state when [soil_carbon].spinup_steady asks. --!
      call set_site_terrain(poly, cfg)
      call refresh_canopy_depth(poly, cfg, no_ledger)
      if (cfg%soil_carbon%soil_carbon_on .and. cfg%soil_carbon%spinup_steady .and.                  &
          cfg%init%init_mode /= INIT_RESTART) then
         pools = soil_carbon_steady_pools(cfg%soil_carbon)
         poly%patch%soil_carbon(1:poly%patch%n) = pools
         if (verbose) write(*,'(a,f8.3,a)') ' soilc : steady-state spin-up, total = ',            &
            sum(pools%carbon), ' kgC/m2'
      end if
      ok = .true.
   end subroutine init_polygon

   !----- Each site's terrain: the polygon's elevation, flat, unless its [[sites]] table gives its own. !
   subroutine set_site_terrain(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      integer(ik) :: is
      poly%site%surface(1:poly%site%n)%elevation_m = poly%elevation_m
      if (.not. allocated(cfg%sites)) return
      do is = 1_ik, poly%site%n
         associate (spec => cfg%sites(is), surf => poly%site%surface(is))
            if (spec%has_elevation) surf%elevation_m = spec%elevation
            surf%slope_deg  = spec%slope
            surf%aspect_deg = spec%aspect
         end associate
      end do
   end subroutine set_site_terrain

   !----- The final stand's fast stores: each patch's soil and canopy air seeded again             !
   !      (seed_patch_stores, as init_bare_sites did when it laid the patches out, so a census       !
   !      stand's fusion leaves no blend behind); each cohort settled in its patch as a recruit is  !
   !      born, its tissues at the canopy-air temperature holding the water they hold at the top   !
   !      soil layer's potential (the restructuring and the traits' acclimation change the tissues !
   !      the setup cohorts were born with); and the snowpack [fast].snow_init_swe asks for. ---!
   subroutine init_fast_state(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      integer(ik) :: ip, i
      real(wp)    :: psi
      call seed_patch_stores(poly, cfg)
      do ip = 1_ik, poly%patch%n
         associate (cas => poly%patch%cas(ip), sw => poly%patch%soil_water(ip),                       &
                    surf => poly%site%surface(poly%patch%owner_site(ip)))
            psi = top_soil_psi(surf%soil_profile, sw%theta(1))
            do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik
               poly%cohort%leaf_temp(i) = cas%can_temp
               poly%cohort%wood_temp(i) = cas%can_temp
               associate (hyd => cfg%pft%hydro(poly%cohort%pft(i)))
                  call set_cohort_tissue_water(poly%cohort, i, psi, hyd%leaf_curve, hyd%wood_curve)
               end associate
               call hold_cohort_tissue_temps(poly%cohort, i)
            end do
         end associate
      end do
      if (cfg%snow_init_swe > 0.0_wp) call seed_snow(poly, cfg)
   end subroutine init_fast_state

   !----- The initial snowpack, [fast].snow_init_swe at [fast].snow_init_temp, one layer per patch. -!
   subroutine seed_snow(poly, cfg)
      type(polygon_t),      intent(inout) :: poly
      type(ed_config_t),    intent(in)    :: cfg
      integer(ik) :: ipp
      do ipp = 1_ik, poly%patch%n
         poly%patch%snow(ipp)%swe(1)         = cfg%snow_init_swe
         poly%patch%snow(ipp)%snow_energy(1) =                                                   &
            temp_to_internal_energy(0.0_wp, cfg%snow_init_swe, cfg%snow_init_temp, 0.0_wp)
         poly%patch%snow(ipp)%snow_depth(1)  = cfg%snow_init_swe / cfg%snow%rho_snow
         poly%patch%snow(ipp)%nlayer         = 1_ik
         call internal_energy_to_temp(poly%patch%snow(ipp)%snow_energy(1),                       &
                                      cfg%snow_init_swe, 0.0_wp,                                 &
                                      poly%patch%snow(ipp)%snow_temp(1),                         &
                                      poly%patch%snow(ipp)%snow_fliq(1))
      end do
   end subroutine seed_snow

   !----- What a restart read: its file and time, and a restructuring the stand still owes. ------!
   subroutine report_restart(cfg, now, pending)
      type(ed_config_t),   intent(in) :: cfg
      type(ed_time_t),     intent(in) :: now
      logical,             intent(in) :: pending
      write(*,'(4a)') ' init  : restart (mode 2) ', trim(cfg%init%restart_file), ' @ ', time_to_string(now)
      if (pending) write(*,'(a)') '         (the boundary''s restructuring is pending: it runs before the first step)'
   end subroutine report_restart

end module ed_init
