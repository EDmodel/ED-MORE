! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_fast_diag -- what the output reads from the fast loop. size_staging sizes the         !
! diagnostic accumulators and the FAST tier's staging to the stand and stamps the          !
! sub-steps; record_substep forms a (patch, sub-step)'s patch row and stages and           !
! accumulates it with the cohorts' rows; fold_window, after the patch loop, folds the FAST !
! means, the forcing echo and each site's forcing (docs/science/diagnostics.md).           !
!==========================================================================================!
module ed_fast_diag
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : latent_heat_vap, stefan_boltzmann, n_soil_layer_max
   use ed_therm_lib,            only : specific_humidity_to_vpd
   use ed_time,                 only : ed_time_t, time_advance_seconds
   use ed_output_types,         only : output_buffers_t
   use ed_diag_vars,            only : N_CO_FAST_ROWS, N_PA_ROWS, patch_diag_block, polygon_diag_block, PA_LE, PA_H, &
                                       PA_RNET, PA_SW_GROUND,                                                    &
                                       PA_LW_GROUND, PA_SW_IN_VIS, PA_SW_IN_NIR, PA_SW_UP_VIS, PA_SW_UP_NIR,     &
                                       PA_LW_UP, PA_USTAR, PA_GGNET, PA_ROUGH, PA_DISPLACE, PA_CAS_TEMP,         &
                                       PA_CAS_SHV, PA_CAS_CO2, PA_GPP, PA_NEE, PA_TRANSP, PA_NPP, PA_RECO,       &
                                       PA_SW_UP, PA_SKIN_TEMP, PA_SOIL_TEMP_TOP, PA_CAS_VPD, PA_W_SURFACE,       &
                                       PA_RESID_ENERGY, PA_RESID_WATER, PA_WIND_CAS_TOP, PA_TAIR_CAS_TOP,        &
                                       PA_Z_CAS_TOP, PY_SW_IN, PY_PRECIP, PY_TAIR, PY_QAIR, PY_PSURF, PY_WIND,   &
                                       PY_LWDOWN, PY_PAR_BEAM, PY_PAR_DIFFUSE, PY_NIR_BEAM, PY_NIR_DIFFUSE,      &
                                       PY_SNOWFALL, PY_CO2, PY_COSZ, PY_RHO_AIR, PY_RAINF, N_PY_ROWS,            &
                                       expand_cohort_diag, expand_patch_diag, site_diag_block, SI_SW_IN,         &
                                       SI_PAR_BEAM, SI_NIR_BEAM, SI_LWDOWN, SI_TAIR, SI_QAIR, SI_PSURF,          &
                                       SI_RHO_AIR, N_SI_ROWS
   use ed_met_types,            only : met_forcing_t
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area
   use ed_fast_integ_types,         only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t
   implicit none
   private

   public :: size_staging, record_substep, fold_window

contains

   !---------------------------------------------------------------------------------------!
   ! Size the diagnostic accumulators to the stand, and, when the FAST tier or a snapshot of  !
   ! a fast-loop value reads them (do_fast), the FAST staging, with the sub-steps' stamps.    !
   ! The accumulators are not zeroed: they run from one output tick to the next, and          !
   ! reset_step_diagnostics zeroes them (docs/science/order_of_processes.md §6).              !
   !---------------------------------------------------------------------------------------!
   subroutine size_staging(poly, nsub, dt_fast, do_fast, step_start, out_bufs)
      type(polygon_t),        intent(inout)           :: poly
      integer(ik),            intent(in)              :: nsub
      real(wp),               intent(in)              :: dt_fast
      logical,                intent(in)              :: do_fast
      type(ed_time_t),      intent(in),      optional :: step_start   !< the slow step's start (the stamps)
      type(output_buffers_t), intent(inout), optional :: out_bufs     !< present when do_fast
      integer(ik) :: ip, is, isub, ks, nl, npatch
      npatch = poly%patch%n
      if (poly%cohort%diag%active) then
         call expand_cohort_diag(poly%cohort%diag, max(poly%cohort%n, 1_ik))
         poly%cohort%diag%n = poly%cohort%n
      end if
      if (poly%patch%diag%active) then
         call expand_patch_diag(poly%patch%diag, max(npatch, 1_ik))
         poly%patch%diag%n = npatch
      end if
      if (.not. do_fast) return
      nl = n_soil_layer_max
      call size_fast_staging(out_bufs, nsub, npatch, poly%site%n, poly%cohort%n, poly%cohort%diag%active)
      out_bufs%n_fast_sub    = nsub
      !----- Only the active layers are data: a site's column and each of its patches' have the site's !
      !      layers, and the rest of the slab reads back as fill. ---------------------------------!
      out_bufs%fast_n_site = poly%site%n
      out_bufs%fast_soil_live = .false.
      out_bufs%fast_soil_patch_live = .false.
      do is = 1_ik, poly%site%n
         ks = (is - 1_ik) * nl
         out_bufs%fast_soil_live(ks + 1_ik:ks + poly%site%surface(is)%soil_profile%n_active) = .true.
      end do
      do ip = 1_ik, npatch
         ks = (ip - 1_ik) * nl
         out_bufs%fast_soil_patch_live(ks + 1_ik:ks + poly%site%surface(poly%patch%owner_site(ip))%soil_profile%n_active) &
            = .true.
      end do
      out_bufs%fast_n_cohort = poly%cohort%n
      out_bufs%fast_n_patch  = npatch
      out_bufs%fast_site = 0.0_wp
      out_bufs%fast_soil_temp = 0.0_wp
      out_bufs%fast_soil_water = 0.0_wp
      out_bufs%fast_cohort_gpp = 0.0_wp
      out_bufs%fast_cohort_height = 0.0_wp
      if (poly%cohort%diag%active) out_bufs%fast_cohort_diag = 0.0_wp
      !----- A FAST record is stamped by its sub-step's start, as the slower tiers stamp theirs. -----!
      if (present(step_start)) then
         do isub = 1_ik, nsub
            out_bufs%fast_time(isub) = time_advance_seconds(step_start, real(isub - 1_ik, wp) * dt_fast)
         end do
      end if
   end subroutine size_staging

   !---------------------------------------------------------------------------------------!
   ! One (patch, sub-step): the patch row, staged and accumulated as it is, and the cohorts'  !
   ! rows. The patch's and its cohorts' slots are disjoint across patches, so the threads     !
   ! need no reduction. NPP is GPP net of maintenance respiration only: the slow allocator    !
   ! charges growth respiration.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine record_substep(poly, ip, isub, dt_fast, static, dynamic, forcing, out, do_fast, out_bufs)
      type(polygon_t),            intent(inout)           :: poly
      integer(ik),                intent(in)              :: ip, isub
      real(wp),                   intent(in)              :: dt_fast
      type(fast_static_state_t),  intent(in)              :: static
      type(fast_dynamic_state_t), intent(in)              :: dynamic
      type(fast_forcing_t),       intent(in)              :: forcing
      type(fast_output_t),        intent(in)              :: out
      logical,                    intent(in)              :: do_fast
      type(output_buffers_t),     intent(inout), optional :: out_bufs   !< present when do_fast
      real(wp)    :: prow(N_PA_ROWS)                !< the patch block's row for this (patch, sub-step)
      real(wp)    :: gpp_patch, npp_patch, rnet, rnet_ts, lw_up
      integer(ik) :: ncoh, i0, i, j, kb, nl
      logical     :: do_cdiag, do_pdiag
      do_cdiag = poly%cohort%diag%active
      do_pdiag = poly%patch%diag%active
      ncoh = static%n
      i0   = poly%patch%cohort_offset(ip)
      nl   = n_soil_layer_max
      if (do_pdiag .or. do_fast) then
         gpp_patch = sum(out%gpp_coh(1:ncoh) * static%nplant(1:ncoh))    ! [umol/m2/s]
         !----- The net radiation the ground, the leaves and the wood absorbed, as their balances      !
         !      booked it over the step. The two-stream's, at the start-of-step temperatures, differs  !
         !      by the change in their emission within the step, which leaves to the sky: the          !
         !      upwelling longwave carries it, so the reported budget closes on the booked net         !
         !      radiation (docs/science/vegetation_energy_dynamics.md, "The emission base"). ----------!
         rnet_ts   = forcing%abs_sw_ground + forcing%abs_lw_ground                                    &
                     + sum(forcing%abs_sw(1:ncoh)) + sum(forcing%abs_lw(1:ncoh))                      &
                     + sum(forcing%abs_sw_wood(1:ncoh)) + sum(forcing%abs_lw_wood(1:ncoh))
         rnet      = out%rnet
         lw_up     = forcing%lw_up + (rnet_ts - rnet)
         npp_patch = gpp_patch - sum((out%leaf_resp_coh(1:ncoh) + out%stem_resp_coh(1:ncoh)              &
                                      + out%root_resp_coh(1:ncoh)) * static%nplant(1:ncoh))
         call patch_diag_row(prow, dt_fast, out%le_flux, out%h_flux, rnet,                              &
                             forcing%abs_sw_ground, forcing%abs_lw_ground,                              &
                             out%aero%ustar, out%aero%ggnet, out%aero%rough, out%aero%displace,         &
                             dynamic%cas%can_temp, dynamic%cas%can_shv, dynamic%cas%can_co2,            &
                             gpp_patch, out%report%exports%nee_last, npp_patch,                         &
                             dynamic%soil_energy%soil_temp(1), out%report%ledgers%whole_energy%resid,   &
                             out%report%ledgers%whole_water%resid, forcing%wind, forcing%air_temp,      &
                             static%z_top, forcing%sw_in_vis, forcing%sw_in_nir, forcing%sw_up_vis,     &
                             forcing%sw_up_nir, lw_up, dynamic%soil_water%w_surface,                    &
                             forcing%press)
      end if
      !----- FAST-tier staging in the patch's and its cohorts' own slots. ---------------------------!
      if (do_fast) then
         out_bufs%fast_patch(:, ip, isub) = prow
         kb = (ip - 1_ik) * nl
         out_bufs%fast_soil_temp_patch(kb+1_ik:kb+nl, isub)  = dynamic%soil_energy%soil_temp(1:nl)
         out_bufs%fast_soil_water_patch(kb+1_ik:kb+nl, isub) = dynamic%soil_water%theta(1:nl)
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            out_bufs%fast_cohort_gpp(i,isub)    = out%gpp_coh(j)
            out_bufs%fast_cohort_height(i,isub) = static%height(j)
            if (do_cdiag) out_bufs%fast_cohort_diag(i, :, isub) = out%cdiag(:, j)
         end do
      end if
      !----- The dt-weighted accumulators (ed_diag_vars). ------------------------------------------  !
      if (do_cdiag) then
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            poly%cohort%diag%v(:, i) = poly%cohort%diag%v(:, i) + out%cdiag(:, j) * dt_fast
            poly%cohort%diag%w(i)    = poly%cohort%diag%w(i)    + dt_fast
         end do
      end if
      if (do_pdiag) call accumulate_patch_diag(poly%patch%diag, ip, dt_fast, prow)
   end subroutine record_substep

   !---------------------------------------------------------------------------------------!
   ! After the patch loop, in patch order, so the sums do not move with the thread count      !
   ! (docs/science/numerical_scheme.md §6a): the FAST means, the forcing echo per sub-step    !
   ! (the polygon's sample with the precipitation split as the patches received it,           !
   ! docs/science/forcing.md §7) and each site's forcing for the coarse tiers.                !
   !---------------------------------------------------------------------------------------!
   subroutine fold_window(poly, dt_fast, met_sample, met_site, rain_sub, snow_sub, do_fast, out_bufs)
      type(polygon_t),        intent(inout)           :: poly
      real(wp),               intent(in)              :: dt_fast
      type(met_forcing_t),    intent(in)              :: met_sample(:)        !< (sub-step) the polygon's samples
      type(met_forcing_t),    intent(in)              :: met_site(:,:)        !< (sub-step, site)
      real(wp),               intent(in)              :: rain_sub(:), snow_sub(:)   !< (sub-step) [kg/m2/s]
      logical,                intent(in)              :: do_fast
      type(output_buffers_t), intent(inout), optional :: out_bufs   !< present when do_fast
      real(wp)    :: echo(N_PY_ROWS)
      integer(ik) :: ip, is, isub, kb, ks, nl, nsub
      nsub = size(met_sample, kind=ik)
      nl   = n_soil_layer_max
      !----- FAST-tier means for each sub-step: the polygon's patch rows, each patch weighted by its    !
      !      share of the polygon, and each site's soil column, each of its patches weighted by its    !
      !      area in the site. -------------------------------------------------------------------------!
      if (do_fast) then
         do isub = 1_ik, nsub
            do ip = 1_ik, poly%patch%n
               kb = (ip - 1_ik) * nl
               ks = (poly%patch%owner_site(ip) - 1_ik) * nl
               out_bufs%fast_site(:, isub) = out_bufs%fast_site(:, isub)                                  &
                                             + patch_polygon_area(poly, ip) * out_bufs%fast_patch(:, ip, isub)
               out_bufs%fast_soil_temp(ks+1_ik:ks+nl, isub)  = out_bufs%fast_soil_temp(ks+1_ik:ks+nl, isub)  &
                                             + poly%patch%area(ip) * out_bufs%fast_soil_temp_patch(kb+1_ik:kb+nl, isub)
               out_bufs%fast_soil_water(ks+1_ik:ks+nl, isub) = out_bufs%fast_soil_water(ks+1_ik:ks+nl, isub) &
                                             + poly%patch%area(ip) * out_bufs%fast_soil_water_patch(kb+1_ik:kb+nl, isub)
            end do
         end do
      end if
      if (do_fast .or. poly%diag%active) then
         do isub = 1_ik, nsub
            call forcing_echo(met_sample(isub), rain_sub(isub), snow_sub(isub), echo)
            if (do_fast) out_bufs%fast_forcing(:, isub) = echo
            if (poly%diag%active) call accumulate_polygon_diag(poly%diag, dt_fast, echo)
         end do
      end if
      if (poly%site%diag%active) then
         do isub = 1_ik, nsub
            do is = 1_ik, poly%site%n
               call accumulate_site_diag(poly%site%diag, is, dt_fast, met_site(isub, is))
            end do
         end do
      end if
      if (do_fast) out_bufs%fast_ready = .true.   ! signal main: the staging holds this slow step
   end subroutine fold_window

   !=======================================================================================!
   !  The patch block's row for one (patch, sub-step): the numbers the physics used, none  !
   !  recomputed (docs/science/diagnostics.md §4, The patch row). The slow rows stay 0.    !
   !=======================================================================================!
   pure subroutine patch_diag_row(row, dt, le_flux, h_flux, rnet, sw_ground, lw_ground,             &
                                  ustar, ggnet, rough, displace, cas_temp, cas_shv, cas_co2, gpp, nee, &
                                  npp, soil_temp_top, resid_energy, resid_water, wind_top, tair_top,  &
                                  z_top, sw_in_vis, sw_in_nir, sw_up_vis, sw_up_nir, lw_up, w_surface, &
                                  press)
      real(wp), intent(out) :: row(N_PA_ROWS)
      real(wp), intent(in)  :: dt                        !< [s]        the sub-step
      real(wp), intent(in)  :: le_flux, h_flux           !< [W/m2]     CAS -> atmosphere
      real(wp), intent(in)  :: rnet                      !< [W/m2]     net all-wave radiation absorbed
      real(wp), intent(in)  :: sw_ground, lw_ground      !< [W/m2]     ground SW / net LW
      real(wp), intent(in)  :: ustar, ggnet              !< [m/s]      friction velocity, ground conductance
      real(wp), intent(in)  :: rough, displace           !< [m]        roughness, displacement height
      real(wp), intent(in)  :: cas_temp, cas_shv, cas_co2 !< [K],[kg/kg],[umol/mol] canopy air
      real(wp), intent(in)  :: gpp, nee, npp             !< [umol/m2/s] gross uptake, net exchange (+ to atm), NPP
      real(wp), intent(in)  :: soil_temp_top             !< [K]        top soil-node temperature
      real(wp), intent(in)  :: resid_energy, resid_water !< [J/m2],[kg/m2] this step's signed ledger residuals
      real(wp), intent(in)  :: wind_top, tair_top        !< [m/s],[K] the forcing at the canopy-air top
      real(wp), intent(in)  :: z_top                     !< [m]   that top (the canopy-air depth)
      !----- Top-of-canopy radiative fluxes per band. -----------------------------------------!
      real(wp), intent(in)  :: sw_in_vis, sw_in_nir, sw_up_vis, sw_up_nir, lw_up
      real(wp), intent(in)  :: w_surface                 !< [kg/m2]   ponded surface water
      real(wp), intent(in)  :: press                     !< [Pa]      air pressure
      row = 0.0_wp
      row(PA_LE)           = le_flux
      row(PA_H)            = h_flux
      row(PA_RNET)         = rnet
      row(PA_SW_GROUND)    = sw_ground
      row(PA_LW_GROUND)    = lw_ground
      row(PA_SW_IN_VIS)    = sw_in_vis
      row(PA_SW_IN_NIR)    = sw_in_nir
      row(PA_SW_UP_VIS)    = sw_up_vis
      row(PA_SW_UP_NIR)    = sw_up_nir
      row(PA_SW_UP)        = sw_up_vis + sw_up_nir
      row(PA_LW_UP)        = lw_up
      row(PA_SKIN_TEMP)    = sqrt(sqrt(max(lw_up, 0.0_wp) / stefan_boltzmann))
      row(PA_USTAR)        = ustar
      row(PA_GGNET)        = ggnet
      row(PA_ROUGH)        = rough
      row(PA_DISPLACE)     = displace
      row(PA_CAS_TEMP)     = cas_temp
      row(PA_CAS_SHV)      = cas_shv
      row(PA_CAS_CO2)      = cas_co2
      !----- VPD is nonlinear in temperature, so it is formed per sub-step, not from mean states. --!
      row(PA_CAS_VPD)      = specific_humidity_to_vpd(cas_temp, cas_shv, press)
      row(PA_W_SURFACE)    = w_surface
      !----- The patch's own forcing, at its canopy-air top (docs/science/forcing.md §8). ---------!
      row(PA_WIND_CAS_TOP) = wind_top
      row(PA_TAIR_CAS_TOP) = tair_top
      row(PA_Z_CAS_TOP)    = z_top
      !----- nee drives the CAS CO2 box (positive to the atmosphere); Reco = NEE + GPP. -----------!
      row(PA_GPP)          = gpp
      row(PA_NEE)          = nee
      row(PA_NPP)          = npp
      row(PA_RECO)         = nee + gpp
      !----- The canopy air's whole evaporative flux to the atmosphere [kg/m2/s], not transpiration  !
      !      alone; each cohort's own transpiration is CO_TRANSP. -------------------------------!
      row(PA_TRANSP)       = le_flux / latent_heat_vap
      row(PA_SOIL_TEMP_TOP) = soil_temp_top
      !----- The whole-column ledgers' signed residuals as leak rates [W/m2, kg/m2/s], positive when  !
      !      store appears from nowhere. --------------------------------------------------------!
      row(PA_RESID_ENERGY) = resid_energy / dt
      row(PA_RESID_WATER)  = resid_water / dt
   end subroutine patch_diag_row

   !----- Fold one (patch, sub-step) row into the patch accumulator, dt-weighted; read as v / w. ---!
   subroutine accumulate_patch_diag(pd, ip, dt, row)
      type(patch_diag_block), intent(inout) :: pd
      integer(ik),            intent(in)    :: ip
      real(wp),               intent(in)    :: dt        !< [s] sample weight
      real(wp),               intent(in)    :: row(N_PA_ROWS)
      pd%v(:, ip) = pd%v(:, ip) + row * dt
      pd%w(ip)    = pd%w(ip)    + dt
   end subroutine accumulate_patch_diag

   !----- Size the FAST tier's staging. Sub-step and site counts are run-constant; the patch and     !
   !      cohort buffers grow to at least twice their size when outgrown. The output layer checks   !
   !      the counts against output.patch_max and cohort_max. ---------------------------------------!
   subroutine size_fast_staging(b, nsub, npatch, nsite, ncohort, with_cdiag)
      type(output_buffers_t), intent(inout) :: b
      integer(ik),            intent(in)    :: nsub, npatch, nsite, ncohort
      logical,                intent(in)    :: with_cdiag   !< the per-cohort diagnostic block is on
      integer(ik) :: nl, np, nc, ns
      nl = n_soil_layer_max
      ns = max(nsite, 1_ik)
      if (.not. allocated(b%fast_time))                                                            &
         allocate(b%fast_time(nsub), b%fast_forcing(N_PY_ROWS, nsub), b%fast_site(N_PA_ROWS, nsub),    &
                  b%fast_soil_temp(nl * ns, nsub), b%fast_soil_water(nl * ns, nsub), b%fast_soil_live(nl * ns))
      np = max(npatch, 1_ik)
      if (allocated(b%fast_patch)) then
         if (size(b%fast_patch, 2, kind=ik) < np) then
            np = max(np, 2_ik * size(b%fast_patch, 2, kind=ik))
            deallocate(b%fast_patch, b%fast_soil_temp_patch, b%fast_soil_water_patch, b%fast_soil_patch_live)
         end if
      end if
      if (.not. allocated(b%fast_patch))                                                           &
         allocate(b%fast_patch(N_PA_ROWS, np, nsub), b%fast_soil_temp_patch(nl * np, nsub),           &
                  b%fast_soil_water_patch(nl * np, nsub), b%fast_soil_patch_live(nl * np))
      nc = max(ncohort, 1_ik)
      if (allocated(b%fast_cohort_gpp)) then
         if (size(b%fast_cohort_gpp, 1, kind=ik) < nc) then
            nc = max(nc, 2_ik * size(b%fast_cohort_gpp, 1, kind=ik))
            deallocate(b%fast_cohort_gpp, b%fast_cohort_height)
         end if
      end if
      if (.not. allocated(b%fast_cohort_gpp))                                                      &
         allocate(b%fast_cohort_gpp(nc, nsub), b%fast_cohort_height(nc, nsub))
      !----- The per-cohort diagnostic rows, staged only while their FAST-tier variables are on. ---!
      if (with_cdiag) then
         nc = size(b%fast_cohort_gpp, 1, kind=ik)
         if (allocated(b%fast_cohort_diag)) then
            if (size(b%fast_cohort_diag, 1, kind=ik) < nc) deallocate(b%fast_cohort_diag)
         end if
         if (.not. allocated(b%fast_cohort_diag)) allocate(b%fast_cohort_diag(nc, N_CO_FAST_ROWS, nsub))
      end if
   end subroutine size_fast_staging

   !----- The polygon's forcing for one sub-step as the fast loop used it (docs/science/forcing.md). -!
   subroutine accumulate_polygon_diag(d, dt, f)
      type(polygon_diag_block), intent(inout) :: d
      real(wp),                 intent(in)    :: dt     !< [s] sample weight
      real(wp),                 intent(in)    :: f(N_PY_ROWS)   !< the sub-step's forcing echo
      d%v = d%v + f * dt
      d%w = d%w + dt
   end subroutine accumulate_polygon_diag

   !----- One sub-step of a site's forcing, in the site block's table (SI_*). ---------------------!
   subroutine accumulate_site_diag(d, is, dt, met)
      type(site_diag_block), intent(inout) :: d
      integer(ik),           intent(in)    :: is
      real(wp),              intent(in)    :: dt     !< [s] sample weight
      type(met_forcing_t),   intent(in)    :: met    !< the site's forcing
      real(wp) :: f(N_SI_ROWS)
      f(SI_SW_IN)    = met%sw_down()
      f(SI_PAR_BEAM) = met%par_beam
      f(SI_NIR_BEAM) = met%nir_beam
      f(SI_LWDOWN)   = met%lw_down
      f(SI_TAIR)     = met%air_temp
      f(SI_QAIR)     = met%air_shv
      f(SI_PSURF)    = met%air_pressure
      f(SI_RHO_AIR)  = met%air_density
      d%v(:, is) = d%v(:, is) + f * dt
      d%w(is)    = d%w(is) + dt
   end subroutine accumulate_site_diag

   !----- The forcing one sub-step used, as the polygon block's row (PY_*): the one list of the     !
   !      forcing echo, accumulated by the coarse tiers and staged by the FAST tier. --------------!
   pure subroutine forcing_echo(met, rainf, snowfall, f)
      type(met_forcing_t), intent(in)  :: met
      real(wp),            intent(in)  :: rainf, snowfall   !< [kg/m2/s] the split the patches used, area-weighted
      real(wp),            intent(out) :: f(N_PY_ROWS)
      f(PY_SW_IN)       = met%sw_down()
      f(PY_PRECIP)      = met%precip
      f(PY_TAIR)        = met%air_temp
      f(PY_QAIR)        = met%air_shv
      f(PY_PSURF)       = met%air_pressure
      f(PY_WIND)        = met%wind
      f(PY_LWDOWN)      = met%lw_down
      f(PY_PAR_BEAM)    = met%par_beam
      f(PY_PAR_DIFFUSE) = met%par_diffuse
      f(PY_NIR_BEAM)    = met%nir_beam
      f(PY_NIR_DIFFUSE) = met%nir_diffuse
      f(PY_SNOWFALL)    = snowfall
      f(PY_CO2)         = met%co2
      f(PY_COSZ)        = met%cosz
      f(PY_RHO_AIR)     = met%air_density
      f(PY_RAINF)       = rainf
   end subroutine forcing_echo

end module ed_fast_diag
