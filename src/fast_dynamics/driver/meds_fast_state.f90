! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_state -- moves the state between the polygon and a fast step's records.        !
! gather_patch fills a patch's static state (its cohorts' structure and traits, what that  !
! structure fixes -- roots, water paths, dry heat capacities (derive_cohort_terms) -- the  !
! canopy geometry and canopy-air top, the slow loop's daily inputs, the site's leaf table) !
! and its dynamic state (the patch's stores and its cohorts' tissue state); scatter_patch  !
! writes the dynamic state back. reset_fast_to_slow and accumulate_fast_to_slow keep the   !
! window's sums the slow step reads: GPP, maintenance respiration, light hours, the        !
! predawn leaf potential and the decomposition scalars. seed_start_canopy_air puts a       !
! start's canopy air at the first forcing sample. static_state_fixture builds a test's     !
! static state through the birth path.                                                     !
!==========================================================================================!
module meds_fast_state
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : cm2_to_m2, c2b_wood, t_3ple, umol_2_kgc, day_sec, par_w_2_umol, pi,      &
                                        tiny_num, tissue_c_leaf, tissue_c_wood_dry
   use meds_numerics,            only : ascending_order
   use meds_config,              only : meds_config_t
   use meds_pft_params,          only : pft_table_t, hydraulics_table_t
   use meds_soil_params,         only : soil_profile_t
   use meds_plant_hydraulics,    only : cohort_root_profile, crown_gravity_head
   use meds_plant_hydr_lib,      only : xylem_max_conductance, rhizosphere_conductance
   use meds_therm_lib,           only : tissue_dry_heat_capacity, cas_temp_to_enthalpy
   use meds_cohort_state_vars,   only : cohort_block, alloc_cohort_block, init_cohort, cohort_tissue_temps,       &
                                        set_cohort_tissue_temps, hold_cohort_tissue_temps
   use meds_met_types,           only : met_forcing_t
   use meds_patch_state_vars,    only : decomp_xi_t
   use meds_polygon_state_vars,  only : polygon_t
   use meds_state_defaults,      only : DMAX_PSI_LEAF_ACCUM_RESET
   use meds_canopy_aerodynamics, only : canopy_roughness, open_canopy_fraction, canopy_burial_fraction
   use meds_soil_heat,           only : soil_layer_temp
   use meds_fast_integ_types,          only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,  &
                                        ensure_static_capacity
   implicit none
   private

   public :: reset_fast_to_slow, gather_patch, scatter_patch, accumulate_fast_to_slow, seed_start_canopy_air
   public :: derive_cohort_terms
   public :: static_state_fixture   ! used only by the tests

   !----- [K] round-off allowed between a cohort's stored tissue temperatures and those its energy gives. !
   real(wp), parameter :: TISSUE_TEMP_TOL    = 1.0e-6_wp

contains

   !----- A start's first window (docs/science/canopy_air_space_biophysics.md, the start): each     !
   !      patch's canopy air takes its site's first forcing sample -- the air's temperature,        !
   !      humidity and CO2 -- and each cohort's tissues are held at that temperature, their water   !
   !      as it is. A restart keeps the canopy air its state file holds. -----------------------!
   subroutine seed_start_canopy_air(poly, met_first)
      type(polygon_t),     intent(inout) :: poly
      type(met_forcing_t), intent(in)    :: met_first(:)   !< (site) each site's first sub-step sample
      integer(ik) :: ip, i
      do ip = 1_ik, poly%patch%n
         associate (cas => poly%patch%cas(ip), m => met_first(poly%patch%owner_site(ip)))
            cas%can_temp     = m%air_temp
            cas%can_shv      = m%air_shv
            cas%can_co2      = m%co2
            cas%can_enthalpy = cas_temp_to_enthalpy(m%air_temp, m%air_shv)
            do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik
               poly%cohort%leaf_temp(i) = m%air_temp
               poly%cohort%wood_temp(i) = m%air_temp
               call hold_cohort_tissue_temps(poly%cohort, i)
            end do
         end associate
      end do
   end subroutine seed_start_canopy_air

   !----- Open the window: zero the fast->slow sums the slow step reads after it, and roll the        !
   !      predawn water status over. Today's maximum leaf potential becomes the value the slow       !
   !      stomatal control uses tomorrow (docs/science/leaf_gas_exchange.md §4.2); a cohort with no  !
   !      sample today keeps its value (psi_leaf <= 0, so the guard tests the reset value). ---------!
   subroutine reset_fast_to_slow(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      poly%cohort%gpp_accum(1:poly%cohort%n)       = 0.0_wp
      poly%cohort%leaf_resp_accum(1:poly%cohort%n) = 0.0_wp
      poly%cohort%stem_resp_accum(1:poly%cohort%n) = 0.0_wp
      poly%cohort%root_resp_accum(1:poly%cohort%n) = 0.0_wp
      poly%cohort%light_hours_accum(1:poly%cohort%n) = 0.0_wp
      where (poly%cohort%dmax_psi_leaf_accum(1:poly%cohort%n) > 0.5_wp * DMAX_PSI_LEAF_ACCUM_RESET)          &
         poly%cohort%dmax_psi_leaf(1:poly%cohort%n) = poly%cohort%dmax_psi_leaf_accum(1:poly%cohort%n)
      poly%cohort%dmax_psi_leaf_accum(1:poly%cohort%n) = DMAX_PSI_LEAF_ACCUM_RESET
      !----- The day's decomposition scalars ([soil_carbon].soil_carbon_on). ----------------------!
      if (cfg%soil_carbon%soil_carbon_on) poly%patch%decomp_xi(1:poly%patch%n) = decomp_xi_t()
   end subroutine reset_fast_to_slow

   !---------------------------------------------------------------------------------------!
   ! Gather patch ip for the window. The static state: the cohorts and what they fix on the  !
   ! site's soil layers, the canopy geometry, the canopy-air top with its roughness and      !
   ! displacement, the slow loop's daily inputs and the site's leaf table. The dynamic       !
   ! state: the patch's stores and each cohort's tissue energy, water and films. Every field !
   ! is assigned, so neither needs seeding.                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine gather_patch(poly, cfg, ip, static, dynamic)
      type(polygon_t),            intent(in)    :: poly
      type(meds_config_t),        intent(in)    :: cfg
      integer(ik),                intent(in)    :: ip
      type(fast_static_state_t),  intent(inout) :: static    !< sized to the largest patch (ensure_static_capacity)
      type(fast_dynamic_state_t), intent(inout) :: dynamic   !< sized to the largest patch (ensure_dynamic_capacity)
      real(wp)    :: sum_lai, sum_wai, t_leaf_e, t_wood_e
      integer(ik) :: i0, ncoh, i, j

      ncoh = poly%patch%cohort_count(ip)
      i0   = poly%patch%cohort_offset(ip)
      call gather_cohorts(static, poly%cohort, i0, ncoh)
      call derive_cohort_terms(static, cfg%pft%hydraulics, poly%site%surface(poly%patch%owner_site(ip))%soil_profile)

      !----- Canopy geometry. The open fraction exp(-(LAI + WAI)) weighs bare ground in the ground  !
      !      conductance (docs/science/canopy_aerodynamics.md eq. (7)), so a patch without cohorts  !
      !      couples as bare soil; the pack at the window's start buries the canopy's roughness,    !
      !      eq. (1). The canopy is at least the height whose roughness is the bare soil's, so a    !
      !      patch without cohorts, or with only shorter ones, has the soil's roughness. -----------!
      sum_lai = 0.0_wp
      sum_wai = 0.0_wp
      do j = 1_ik, ncoh
         sum_lai = sum_lai + static%lai(j)
         sum_wai = sum_wai + static%wai(j)
      end do
      static%geom%veg_height = cfg%aero%soil_rough / cfg%aero%z0m_ratio
      do j = 1_ik, ncoh
         static%geom%veg_height = max(static%geom%veg_height, static%height(j))
      end do
      static%geom%opencan_frac = open_canopy_fraction(sum_lai, sum_wai)
      static%geom%snowfac      = canopy_burial_fraction(poly%patch%snow(ip)%snow_depth(1), static%geom%veg_height)
      !----- The canopy-air top and roughness, fixed for the day. Each forcing sample is moved to   !
      !      that top along the patch's log profile (docs/science/forcing.md §8). ----------------!
      static%z_top = poly%patch%cas(ip)%can_depth
      call canopy_roughness(cfg%aero, static%geom, static%rough, static%displace)
      !----- The slow loop's daily inputs, held for the day and never written back: the soil-carbon !
      !      pools (soil_carbon_step is their only writer; 0 with soil carbon off), the turnover    !
      !      shed water and the slow-step CO2 rate. ---------------------------------------------!
      static%soil_carbon     = poly%patch%soil_carbon(ip)
      static%shed_water_rate = poly%patch%shed_water_rate(ip)
      static%shed_water_enth = poly%patch%shed_water_enth(ip)
      static%slow_co2_rate   = poly%patch%slow_co2_rate(ip)
      !----- The site's leaf table, acclimated to its growth temperature by the slow step. ---------!
      static%leaf_photo = poly%site%leaf_photo(poly%patch%owner_site(ip))

      !----- The patch's stores. -----------------------------------------------------------------!
      dynamic%cas           = poly%patch%cas(ip)
      dynamic%soil_energy   = poly%patch%soil_energy(ip)
      dynamic%soil_water    = poly%patch%soil_water(ip)
      !----- The soil's temperatures and liquid fractions are what its energies give, whatever wrote    !
      !      them last (a restart, a patch fusion's blend): the same function the commit diagnoses      !
      !      with, so a consistent state reads back unchanged. -----------------------------------------!
      associate (prof => poly%site%surface(poly%patch%owner_site(ip))%soil_profile)
         call soil_layer_temp(dynamic%soil_energy%soil_energy(1:prof%n_active),                           &
                              dynamic%soil_water%theta(1:prof%n_active), prof%dry_heat_capacity(1:prof%n_active), &
                              prof%layer_curve(1:prof%n_active), dynamic%soil_energy%soil_temp(1:prof%n_active), &
                              dynamic%soil_energy%soil_fliq(1:prof%n_active))
      end associate
      dynamic%snow          = poly%patch%snow(ip)
      dynamic%adapt_dt_last = poly%patch%adapt_dt_last(ip)   ! per-patch, not loop-carried
      do j = 1_ik, ncoh
         i = i0 + j - 1_ik
         !----- Tissue energy is the state the fast loop advances (per m2 of ground here). A slow    !
         !      operator that changed a tissue's mass without hold_cohort_tissue_temps leaves the    !
         !      stored temperatures off the ones the energy gives, and stops the run here. ---------!
         call cohort_tissue_temps(poly%cohort, i, t_leaf_e, t_wood_e)
         if (abs(t_leaf_e - poly%cohort%leaf_temp(i)) > TISSUE_TEMP_TOL .or.                          &
             abs(t_wood_e - poly%cohort%wood_temp(i)) > TISSUE_TEMP_TOL) then
            write(*, '(a,i0,4(a,es12.5))') ' fast_dynamics: cohort ', i, ' leaf ', poly%cohort%leaf_temp(i),    &
                  ' K stored vs ', t_leaf_e, ' K from its energy; wood ', poly%cohort%wood_temp(i), ' vs ', t_wood_e
            error stop 'fast_dynamics: a cohort''s tissue temperatures are not the ones its energy gives'
         end if
         dynamic%leaf_energy(j) = poly%cohort%leaf_energy(i) * poly%cohort%nplant(i)
         dynamic%wood_energy(j) = poly%cohort%wood_energy(i) * poly%cohort%nplant(i)
         dynamic%leaf_temp(j)   = poly%cohort%leaf_temp(i)
         dynamic%wood_temp(j)   = poly%cohort%wood_temp(i)
         !----- Tissue water is adopted as it is: fitting it to the capacity the biomass allows is   !
         !      the slow step's (docs/science/soil_plant_water.md §2). ------------------------------!
         dynamic%leaf_water_mass(j) = poly%cohort%leaf_water_mass(i)
         dynamic%wood_water_mass(j) = poly%cohort%wood_water_mass(i)
         !----- The canopy films, per plant on the cohort and per m2 of ground here (0 is dry). ------!
         dynamic%leaf_surf_water(j) = poly%cohort%leaf_surf_water(i) * poly%cohort%nplant(i)
         dynamic%wood_surf_water(j) = poly%cohort%wood_surf_water(i) * poly%cohort%nplant(i)
      end do
   end subroutine gather_patch

   !----- Write patch ip's evolved state back: its stores and each cohort's tissue state. ----------!
   subroutine scatter_patch(dynamic, poly, ip)
      type(fast_dynamic_state_t), intent(in)    :: dynamic
      type(polygon_t),            intent(inout) :: poly
      integer(ik),                intent(in)    :: ip
      integer(ik) :: i0, i, j
      i0 = poly%patch%cohort_offset(ip)
      poly%patch%cas(ip)           = dynamic%cas
      poly%patch%soil_energy(ip)   = dynamic%soil_energy
      poly%patch%soil_water(ip)    = dynamic%soil_water
      poly%patch%snow(ip)          = dynamic%snow
      poly%patch%adapt_dt_last(ip) = dynamic%adapt_dt_last
      do j = 1_ik, poly%patch%cohort_count(ip)
         i = i0 + j - 1_ik
         poly%cohort%leaf_energy(i)     = dynamic%leaf_energy(j) / poly%cohort%nplant(i)
         poly%cohort%wood_energy(i)     = dynamic%wood_energy(j) / poly%cohort%nplant(i)
         poly%cohort%leaf_water_mass(i) = dynamic%leaf_water_mass(j)
         poly%cohort%wood_water_mass(i) = dynamic%wood_water_mass(j)
         poly%cohort%leaf_surf_water(i) = dynamic%leaf_surf_water(j) / poly%cohort%nplant(i)
         poly%cohort%wood_surf_water(i) = dynamic%wood_surf_water(j) / poly%cohort%nplant(i)
         call set_cohort_tissue_temps(poly%cohort, i)
      end do
   end subroutine scatter_patch

   !---------------------------------------------------------------------------------------!
   ! Add one fast step of patch ip to the sums the slow step reads. The cohorts' slots are   !
   ! disjoint across patches, so the threads need no reduction.                              !
   !---------------------------------------------------------------------------------------!
   subroutine accumulate_fast_to_slow(poly, cfg, ip, forcing, out)
      type(polygon_t),      intent(inout) :: poly
      type(meds_config_t),  intent(in)    :: cfg
      integer(ik),          intent(in)    :: ip
      type(fast_forcing_t), intent(in)    :: forcing
      type(fast_output_t),  intent(in)    :: out
      real(wp)    :: dt_fast_days
      integer(ik) :: i0, i, j
      !----- The day's integrals of the per-pool environment scalar and the matrix Rh, which the     !
      !      daily soil_carbon_step consumes (docs/science/soil_carbon.md §6). ----------------------!
      if (cfg%soil_carbon%soil_carbon_on) then
         dt_fast_days = cfg%dt_fast / day_sec
         associate (dx => poly%patch%decomp_xi(ip))
            dx%xi_days = dx%xi_days + out%report%exports%xi_step * dt_fast_days
            dx%rh_fast = dx%rh_fast + out%report%exports%rh_matrix_step * dt_fast_days
         end associate
      end if
      i0 = poly%patch%cohort_offset(ip)
      do j = 1_ik, poly%patch%cohort_count(ip)
         i = i0 + j - 1_ik
         !----- Gross GPP and the maintenance losses [umol/plant/s -> kgC/plant], kept apart: the      !
         !      slow step's compute_carbon_allocation nets them. -----------------------------------!
         poly%cohort%gpp_accum(i)       = poly%cohort%gpp_accum(i)       + out%gpp_coh(j)       * cfg%dt_fast * umol_2_kgc
         poly%cohort%leaf_resp_accum(i) = poly%cohort%leaf_resp_accum(i) + out%leaf_resp_coh(j) * cfg%dt_fast * umol_2_kgc
         poly%cohort%stem_resp_accum(i) = poly%cohort%stem_resp_accum(i) + out%stem_resp_coh(j) * cfg%dt_fast * umol_2_kgc
         poly%cohort%root_resp_accum(i) = poly%cohort%root_resp_accum(i) + out%root_resp_coh(j) * cfg%dt_fast * umol_2_kgc
         !----- Hours of light at the cohort's top, the phenology light cue: a step counts when the    !
         !      PAR there exceeds the PFT's par_min. -----------------------------------------------!
         if (forcing%par_top(j) * par_w_2_umol > cfg%pft%pheno(poly%cohort%pft(i))%par_min)               &
            poly%cohort%light_hours_accum(i) = poly%cohort%light_hours_accum(i) + cfg%dt_fast / 3600.0_wp
         !----- Running daily maximum of psi_leaf, the predawn value that drives tomorrow's stomatal   !
         !      control (docs/science/leaf_gas_exchange.md §4.2). -----------------------------------!
         poly%cohort%dmax_psi_leaf_accum(i) = max(poly%cohort%dmax_psi_leaf_accum(i), out%psi_leaf_coh(j))
      end do
   end subroutine accumulate_fast_to_slow

   !---------------------------------------------------------------------------------------!
   ! Copy one patch's contiguous cohort section [i0, i0+ncoh-1] into the static state. Each  !
   ! line is a copy or a per-ground index; the geometry is cached on the cohort block        !
   ! (set_cohort_size_from_carbon), and the indices are formed here, so none carries a stale !
   ! plant density; then the cohorts' order from the ground up. A new per-cohort input to    !
   ! the fast step is one field on the cohort block, one on fast_static_state_t and one line !
   ! here.                                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine gather_cohorts(static, cohort, i0, ncoh)
      type(fast_static_state_t), intent(inout) :: static
      type(cohort_block),        intent(in)    :: cohort
      integer(ik),               intent(in)    :: i0      !< first cohort of the patch (CSR offset)
      integer(ik),               intent(in)    :: ncoh    !< cohorts in the patch
      integer(ik) :: i, j, pf
      call ensure_static_capacity(static, ncoh)
      associate (pft => cohort%pft_table)
      do j = 1_ik, ncoh
         i = i0 + j - 1_ik
         pf = cohort%pft(i)
         !----- Demographic inputs. ---------------------------------------------------------!
         static%pft(j)           = cohort%pft(i)
         static%nplant(j)        = cohort%nplant(i)
         static%dbh(j)           = cohort%dbh(i)
         static%height(j)        = cohort%height(i)
         static%leaf_area(j)     = cohort%leaf_area(i)
         !----- Per-ground indices: nplant * per-plant area, formed here, never stored. -------!
         static%lai(j)           = cohort%nplant(i) * cohort%leaf_area(i)
         static%wai(j)           = cohort%nplant(i) * cohort%wood_area(i)
         !----- Carbon pools the fast loop reads (capacitance + thermal stores). --------------!
         static%bleaf(j)         = cohort%leaf_carbon(i)
         static%broot(j)         = cohort%fineroot_carbon(i)
         static%bsap(j)          = cohort%sapwood_carbon(i)   ! sapwood ring -> hydraulics
         static%bwood(j)         = cohort%wood_carbon(i)      ! all wood     -> thermal store
         static%sap_area(j)      = cohort%sapwood_area(i) * cm2_to_m2   ! [cm2] -> [m2] (hydraulics)
         !----- Plastic leaf capacities + yesterday's daily-max leaf potential. ---------------!
         static%vcmax25(j)       = cohort%vcmax25(i)
         static%rd25(j)          = cohort%rd25(i)
         static%dmax_psi_leaf(j) = cohort%dmax_psi_leaf(i)
         !----- Canopy-element geometry, emissivities and respiration traits (the PFT's, from the  !
         !      block's table). --------------------------------------------------------------!
         static%leaf_width(j)          = pft%leaf_width(pf)
         static%branch_diam(j)         = pft%branch_diameter(pf)
         static%aboveground_frac(j)    = pft%allom(pf)%aboveground_frac
         static%is_woody(j)            = pft%is_woody(pf)
         static%stem_resp_factor25(j)  = pft%stem_resp_factor25(pf)
         static%stem_resp_sapwood25(j) = pft%stem_resp_sapwood25(pf)
         !----- Sapwood volume [m3/plant]: its carbon as dry mass over the PFT's wood density       !
         !      [g/cm3 = 1000 kg/m3]; biomass grows as density^0.976, so it barely depends on it. ---!
         static%sap_volume(j)          = cohort%sapwood_carbon(i) * c2b_wood / (1.0e3_wp * pft%wood_density(pf))
         static%root_resp_factor25(j)  = pft%root_resp_factor25(pf)
         static%leaf_emissivity(j)     = pft%leaf_emissivity(pf)
         static%wood_emissivity(j)     = pft%wood_emissivity(pf)
      end do
      end associate
      call ascending_order(static%height, ncoh, static%order)
   end subroutine gather_cohorts

   !---------------------------------------------------------------------------------------!
   ! What the cohorts' structure fixes for the slow step on the site's soil layers: each     !
   ! cohort's fine roots by layer, within its own rooting depth, and the patch's profile,    !
   ! the cohorts' weighted by their fine-root carbon (a patch with no roots puts its weight  !
   ! in the top layer: it has no root respiration to weight); the water network's paths --   !
   ! the intact xylem from the sapwood, the head to the crown, each layer's rhizosphere from !
   ! the roots' length there, half the distance between the roots of every cohort in it     !
   ! (docs/science/plant_hydraulics.md eq. (14)), and the carbon of the wood store (sapwood  !
   ! and fine roots); and the dry heat capacity of the leaves and of all the wood            !
   ! (docs/science/vegetation_energy_dynamics.md eq. (8)).                                   !
   !---------------------------------------------------------------------------------------!
   subroutine derive_cohort_terms(static, hyd, profile)
      type(fast_static_state_t), intent(inout) :: static
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(soil_profile_t),      intent(in)    :: profile
      real(wp)    :: z_top(profile%n_active), z_bot(profile%n_active), root_carbon, half_spacing
      real(wp)    :: root_len(profile%n_active, static%n)   !< [m/plant] each cohort's fine roots by layer
      integer(ik) :: i, k, n, nsl
      n   = static%n
      nsl = profile%n_active
      z_top(1:nsl) = -profile%soil_layer_z(1:nsl)
      z_bot(1:nsl) = -profile%soil_layer_z(2:nsl+1)
      static%root_frac_patch = 0.0_wp
      do i = 1_ik, n
         associate (hp => hyd%pft(static%pft(i)))
            call cohort_root_profile(static%height(i), hp, z_top(1:nsl), z_bot(1:nsl), static%root_frac(1:nsl, i))
            static%root_frac_patch(1:nsl) = static%root_frac_patch(1:nsl)                                    &
                                            + static%nplant(i)*static%broot(i)*static%root_frac(1:nsl, i)
            root_len(1:nsl, i)   = static%broot(i) * static%root_frac(1:nsl, i) * hp%specific_root_length
            static%k_xylem(i)    = static%nplant(i) * xylem_max_conductance(hp%wood_kmax, static%sap_area(i),   &
                                                                            static%height(i), hp%vessel_curl)
            static%crown_head(i) = crown_gravity_head(static%height(i))
         end associate
         static%bwood_water(i)   = static%bsap(i) + static%broot(i)
         static%leaf_dry_hcap(i) = tissue_dry_heat_capacity(static%nplant(i), static%bleaf(i), tissue_c_leaf)
         static%wood_dry_hcap(i) = tissue_dry_heat_capacity(static%nplant(i), static%bwood(i), tissue_c_wood_dry)
      end do
      root_carbon = sum(static%root_frac_patch(1:nsl))
      if (root_carbon > tiny_num) then
         static%root_frac_patch(1:nsl) = static%root_frac_patch(1:nsl) / root_carbon
      else
         static%root_frac_patch(1:nsl) = 0.0_wp
         static%root_frac_patch(1) = 1.0_wp
      end if
      do k = 1_ik, nsl
         half_spacing = 1.0_wp / sqrt(pi * max(sum(static%nplant(1:n) * root_len(k, 1:n)) / profile%dz(k),    &
                                               tiny_num))
         do i = 1_ik, n
            static%g_rhizo(k, i) = rhizosphere_conductance(static%nplant(i), root_len(k, i), half_spacing,      &
                                                           hyd%pft(static%pft(i))%fine_root_radius)
         end do
      end do
   end subroutine derive_cohort_terms

   !----- A fixture's static state: n cohorts born through init_cohort as a real recruit is,      !
   !      gathered by gather_cohorts with what they fix on the soil profile, and the PFT table's   !
   !      leaf table. The caller gives only each cohort's PFT, diameter and density and the soil;  !
   !      the canopy geometry and the daily inputs keep their defaults. A test that edits a       !
   !      cohort's structure afterwards calls derive_cohort_terms again. ------------------------!
   subroutine static_state_fixture(static, pft, ipft, dbh, nplant, profile)
      type(fast_static_state_t), intent(inout) :: static
      type(pft_table_t),         intent(in)    :: pft
      integer(ik),               intent(in)    :: ipft(:)    !< (n) PFT index per cohort
      real(wp),                  intent(in)    :: dbh(:)     !< (n) [cm]       diameter
      real(wp),                  intent(in)    :: nplant(:)  !< (n) [plant/m2] density
      type(soil_profile_t),      intent(in)    :: profile
      type(cohort_block) :: block
      integer(ik)        :: j, n
      n = int(size(ipft), ik)
      call alloc_cohort_block(block, n, 1_ik)   ! growth ring buffer irrelevant to a fixture
      block%pft_table = pft                     ! as alloc_polygon sets it
      block%n = n
      do j = 1_ik, n
         !----- The static state carries no tissue water or temperature, so the birth state reaches  !
         !      nothing: saturated (0 MPa) at the triple point. ----------------------------------!
         call init_cohort(block, j, pft, ipft(j), 1_ik, nplant(j), dbh(j), t_3ple, 0.0_wp)
      end do
      call gather_cohorts(static, block, 1_ik, n)
      call derive_cohort_terms(static, pft%hydraulics, profile)
      static%leaf_photo = pft%leaf_photo
   end subroutine static_state_fixture

end module meds_fast_state
