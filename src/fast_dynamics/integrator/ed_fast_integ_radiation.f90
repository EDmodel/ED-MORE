! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_radiation -- the step's radiation: the sky (fast_forcing_t's) through the
! two-stream canopy onto the ground as the ground's stage left it, so the canopy, the ground's
! optics and its surface's balance see one ground (docs/science/canopy_radiation_transfer.md
! §4, docs/science/ground_surface_layers.md, "The ground's radiation").
!---~---
module ed_fast_integ_radiation
   use ed_kinds,             only : wp, ik
   use ed_constants,         only : tiny_num, par_w_2_umol
   use ed_config,            only : ed_config_t
   use ed_pft_params,        only : N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use ed_canopy_radiation,  only : rad_env_t, rad_out_t, alloc_rad_env, ground_optics_env_t, canopy_radiation,  &
                                    ground_optics
   use ed_surface,           only : snow_blend_ground_optics, ground_tile_radiation
   use ed_patch_state_vars,  only : water_layer_state_t, ground_top_temp
   use ed_fast_integ_types,  only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t
   implicit none
   private

   public :: canopy_ground_radiation

contains

   !---~---
   !   Run the patch's two-stream canopy radiation and fill each cohort's absorbed SW,
   ! incident-equivalent PAR and net LW, the ground's SW and LW, each ground tile's own, and the
   ! top-of-canopy fluxes. The ground's optics are the bare soil's blended toward the pack's by
   ! the cover its stage set, and it emits at its surface's temperature: the porous layers'
   ! where there are any, else the top soil node's. Cohorts arrive tallest first and the solver
   ! takes them bottom to top, in static%order, so outputs go back through it.
   !---~---
   subroutine canopy_ground_radiation(cfg, static, dynamic, pack, cover, forcing)
      type(ed_config_t),          intent(in)    :: cfg
      type(fast_static_state_t),  intent(in)    :: static     ! the cohorts (gather order, top first)
      type(fast_dynamic_state_t), intent(in)    :: dynamic    ! the tissues' and the ground's temperatures
      type(water_layer_state_t),  intent(in)    :: pack       ! the pack after the ground's stage
      real(wp),                   intent(in)    :: cover      ! [-] its cover of the ground
      type(fast_forcing_t),       intent(inout) :: forcing    ! the sky in; the absorbed radiation out
      integer(ik) :: j, i, ig, ncoh
      integer(ik) :: pft_bt(static%n)
      real(wp)    :: lai_bt(static%n), wai_bt(static%n), hgt_bt(static%n)
      real(wp)    :: tleaf_bt(static%n), twood_bt(static%n)
      real(wp)    :: soil_albedo(N_RAD_BAND_DEFAULT), soil_temp_top
      type(rad_env_t)           :: rf
      type(rad_out_t)           :: flux
      type(ground_optics_env_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)

      ncoh = static%n
      !--- The uncovered ground's albedo per band from [soil] (no longwave reflection) and its emission temperature.
      soil_albedo   = [cfg%soil%ground_albedo_vis, cfg%soil%ground_albedo_nir, 0.0_wp]
      soil_temp_top = ground_top_temp(dynamic%ground_surface, dynamic%soil_energy%soil_temp(1))

      !--- The cohorts bottom -> top; the tissues emit at their start-of-step temperatures.
      do j = 1_ik, ncoh
         i = static%order(j)
         pft_bt(j)   = static%pft(i)
         lai_bt(j)   = static%lai(i)
         hgt_bt(j)   = static%height(i)
         wai_bt(j)   = static%wai(i)
         tleaf_bt(j) = dynamic%leaf_temp(i)
         twood_bt(j) = dynamic%wood_temp(i)
      end do

      call alloc_rad_env(rf, N_RAD_BAND_DEFAULT)
      rf%cosz = forcing%cosz
      rf%incid_beam(RAD_VIS) = forcing%par_beam
      rf%incid_diff(RAD_VIS) = forcing%par_diffuse
      rf%incid_beam(RAD_NIR) = forcing%nir_beam
      rf%incid_diff(RAD_NIR) = forcing%nir_diffuse
      rf%incid_beam(RAD_LW)  = 0.0_wp
      rf%incid_diff(RAD_LW)  = forcing%lw_down

      !--- The ground's optics: the bare soil's, blended toward the pack's by its cover.
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = soil_albedo
      surf%soil_emiss  = cfg%soil%ground_emissivity
      surf%soil_temp   = soil_temp_top
      call snow_blend_ground_optics(pack, cover, cfg%ground_surface%water, surf%soil_albedo(RAD_VIS),            &
                                    surf%soil_albedo(RAD_NIR), surf%soil_emiss, surf%soil_temp)
      he = [.false., .false., .true.]
      call ground_optics(surf, N_RAD_BAND_DEFAULT, he, rf%grnd_refl, rf%grnd_emiss)

      call canopy_radiation(cfg%pft%optics, rf, ncoh, pft_bt, hgt_bt, lai_bt, wai_bt, tleaf_bt, twood_bt, flux)

      !--- inverse-scatter: RT index j (bottom->top) maps to gather index static%order(j).
      forcing%apar = sum(flux%abs_leaf(RAD_VIS, 1:ncoh))      ! the leaves' absorbed PAR, for the diagnostics
      do j = 1_ik, ncoh
         ig = static%order(j)
         forcing%abs_sw(ig)      = flux%abs_leaf(RAD_VIS, j) + flux%abs_leaf(RAD_NIR, j)   ! absorbed leaf SW
         forcing%abs_par(ig)     = flux%abs_leaf(RAD_VIS, j)                                          &
                                   / max(cfg%leaf_physiology%leaf_absorptance, tiny_num)   ! incident PAR
         forcing%abs_lw(ig)      = flux%abs_leaf(RAD_LW, j)                  ! net leaf LW, its own emission included
         forcing%abs_sw_wood(ig) = flux%abs_wood(RAD_VIS, j) + flux%abs_wood(RAD_NIR, j)   ! absorbed wood SW
         forcing%abs_lw_wood(ig) = flux%abs_wood(RAD_LW, j)                  ! net wood LW, its own emission included
         forcing%par_top(ig)     = flux%incid_top(RAD_VIS, j)                ! PAR at the cohort's top
      end do
      forcing%abs_sw_ground = (flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS))                      &
                            + (flux%dn_ground(RAD_NIR) - flux%up_ground(RAD_NIR))
      forcing%abs_lw_ground = flux%dn_ground(RAD_LW) - flux%up_ground(RAD_LW)   ! net ground LW, emission included
      !--- What each tile takes of the radiation reaching the ground, at its own optics.
      call ground_tile_radiation(flux%dn_ground(RAD_VIS), flux%dn_ground(RAD_NIR), flux%dn_ground(RAD_LW),     &
                                 soil_albedo(RAD_VIS), soil_albedo(RAD_NIR), cfg%soil%ground_emissivity,       &
                                 soil_temp_top, pack, cfg%ground_surface%water, forcing%rnet_soil, forcing%sw_snow, &
                                 forcing%sw_ground)
      forcing%lw_down_ground = flux%dn_ground(RAD_LW)
      forcing%par_per_w      = par_w_2_umol                    ! true VIS absorbed -> photon flux
      !---~---
      !   Top-of-canopy fluxes: each band's albedo times the incident it was formed from, kept as
      ! fluxes (docs/science/canopy_radiation_transfer.md, Surface radiative fluxes on the
      ! diagnostic path).
      !---~---
      forcing%sw_in_vis = rf%incid_beam(RAD_VIS) + rf%incid_diff(RAD_VIS)
      forcing%sw_in_nir = rf%incid_beam(RAD_NIR) + rf%incid_diff(RAD_NIR)
      forcing%sw_up_vis = flux%albedo(RAD_VIS) * forcing%sw_in_vis
      forcing%sw_up_nir = flux%albedo(RAD_NIR) * forcing%sw_in_nir
      forcing%lw_up     = flux%albedo(RAD_LW)  * (rf%incid_beam(RAD_LW) + rf%incid_diff(RAD_LW))
   end subroutine canopy_ground_radiation

end module ed_fast_integ_radiation
