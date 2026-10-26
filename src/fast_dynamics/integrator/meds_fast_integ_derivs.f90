! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_integ_derivs -- the right-hand sides the fast loop's schemes evaluate, free of !
! side effects: no flux hidden in a backward-Euler denominator and no store mutated, so    !
! any tableau can combine the evaluations (docs/science/numerical_scheme.md §3).           !
!   * surface_derivs -- leaf and wood energy balance, ground skin, the three canopy-air    !
!     twins;                                                                               !
!   * march_derivs -- surface_derivs plus the soil-heat column and the soil and plant      !
!     water, which move at the water network's frozen rates (§2);                          !
!   * cas_conductances -- the canopy-air <-> atmosphere conductances at a live canopy-air  !
!     state;                                                                               !
!   * water_derivs -- the soil-plant water network's flows, inflows and Jacobian at one    !
!     state, and set_water_state its state from the model's                                !
!     (docs/science/soil_plant_water.md §1-§2).                                            !
!==========================================================================================!
module meds_fast_integ_derivs
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o, grav_head
   use meds_therm_lib, only : cas_enthalpy_to_temp, sat_specific_humidity, sat_specific_humidity_temp_deriv
   use meds_soil_heat, only : soil_heat_env_t, soil_energy_time_deriv, set_soil_heat_env, soil_layer_temp
   use meds_surface, only : surface_env_t, ground_surface_exchange
   use meds_patch_state_vars,   only : soil_energy_state_t
   use meds_canopy_air_space, only : cas_env_t, cas_box_time_deriv, cas_condensation
   use meds_canopy_aerodynamics, only : canopy_air_conductances
   use meds_plant_energy, only : plant_energy_env_t, canopy_tissue_energy
   use meds_canopy_water, only : canopy_water_out_t
   use meds_fast_integ_types, only : surface_state_t, surface_tend_t, cas_boundary_t,                                  &
                               march_state_t, frozen_t, march_tend_t, fast_static_state_t,                    &
                               water_drive_t, water_tend_t, water_work_t, ensure_water_capacity, POND, WOOD, LEAF
   use meds_site_state_vars, only : site_surface_t
   use meds_soil_hydr_lib, only : soil_water_relations, soil_theta_to_unknown, pond_water_relations,              &
                                  pond_water_to_head, soil_ice_t
   use meds_plant_hydr_lib, only : tissue_water_relations, tissue_water_to_rwc, tissue_rwc_to_unknown,            &
                                   tissue_psi_to_rwc
   use meds_plant_hydraulics, only : supply_limited_transpiration, root_water_uptake_profile, xylem_sapflow
   use meds_soil_water, only : soil_water_flux, ground_evaporation_rate, bottom_flux
   use meds_linear_algebra, only : arrowhead_system_t
   use meds_soil_params, only : soil_profile_t, soil_params_t
   use meds_pft_params, only : hydraulics_table_t
   implicit none
   private

   public :: surface_derivs, march_derivs, cas_conductances
   public :: set_water_state, water_derivs

contains

   !---------------------------------------------------------------------------------------!
   ! refresh_cas_conductances -- re-solve only the bulk Monin-Obukhov surface layer at a live canopy-  !
   ! air state for the three CAS <-> atmosphere conductances (docs/science/canopy_aerodynamics.md      !
   ! eq. (6)). The per-cohort boundary layers stay frozen (docs/science/numerical_scheme.md §2), as do !
   ! the geometry and the reference-level state: only (T_cas, q_cas) is live. Vapour and CO2 reuse the !
   ! heat transfer factor (z0q = z0h), CO2 on the molar capacity at the live humidity.                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp) :: tcas
      tcas = cas_enthalpy_to_temp(cas_enthalpy, cas_shv)
      !----- The pre-pass solved at tcas too, so only the state differs from its solve. ----------!
      call canopy_air_conductances(cas%aero_cfg, cas%env, cas%geom, cas%rho, tcas, cas_shv,              &
                                   g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine refresh_cas_conductances

   !---------------------------------------------------------------------------------------!
   ! cas_conductances -- the one place that decides whether to re-solve the surface layer: a fresh    !
   ! solve when the bundle carries live MO inputs, else the bundle's own values. One rule keeps a      !
   ! tendency and the boundary flux charged against it in agreement (explicit_stage_bflux, RK45).   !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      if (cas%mo_live) then
         call refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      else
         g_atm_heat = cas%box%atm_conductance_enthalpy
         g_atm_vapour = cas%box%atm_conductance_vapor
         g_atm_co2 = cas%box%atm_conductance_co2
      end if
   end subroutine cas_conductances

   !---------------------------------------------------------------------------------------!
   ! surface_derivs -- the CAS surface-block RHS (leaf/wood energy balance + ground skin + the three !
   ! CAS twins), side-effect-free. Every tissue flux is linearized about the current tcas, which is  !
   ! also the longwave emission base (docs/science/vegetation_energy_dynamics.md eqs. (1)-(2)); the  !
   ! twins' tendencies are docs/science/canopy_air_space_biophysics.md eq. (1).                      !
   !---------------------------------------------------------------------------------------!
   pure subroutine surface_derivs(y, cas, tissue, film, ground, t_ground, n, f)
      type(surface_state_t),        intent(in)  :: y
      type(cas_boundary_t),         intent(in)  :: cas      !< CAS <-> atmosphere boundary (a stage may pass its own conductances)
      type(plant_energy_env_t),     intent(in)  :: tissue   !< per-cohort leaf/wood energy coefficients
      type(canopy_water_out_t),     intent(in)  :: film     !< interception film
      type(surface_env_t),          intent(in)  :: ground   !< what the ground surface reads, the snow pack among it
      real(wp),               intent(in)  :: t_ground   !< [K] soil-top temperature at this evaluation (live, not frozen)
      integer(ik),            intent(in)  :: n
      !----- Filled in place; the per-cohort arrays are allocated only when missing or resized, so a  !
      !      Newton solve that evaluates this many times on one record allocates them once. ---------!
      type(surface_tend_t),   intent(inout) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt
      real(wp)    :: coh_h, coh_qw, coh_transp, coh_rnet, coh_film_evap   !< the canopy's, into the canopy air
      type(surface_env_t) :: senv      !< the ground surface's surroundings at this evaluation
      real(wp)    :: gah_l, gaw_l, gac_l      !< the conductances this evaluation uses
      type(cas_env_t) :: cas_env

      call ensure_surface_tend(f, n)

      tcas   = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, cas%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, cas%press)

      !----- Every cohort's wood and leaf at this canopy air (docs/science/vegetation_energy_dynamics.md);  !
      !      with canopy_water_on off, f_wet_c = 0 and the films' pathway vanishes. ------------------------!
      call canopy_tissue_energy(tissue, film, tcas, qcas, qsat_c, dqdt, cas%rho, n, f%leaf_temp, f%wood_temp,     &
                                f%transp_c, f%film_evap_leaf, f%film_evap_wood, coh_h, coh_qw, coh_transp,        &
                                coh_film_evap, coh_rnet)

      !----- The ground surface at this evaluation's canopy air and soil top, as the tissues' is at       !
      !      theirs: the snowfac-weighted pack and the bare soil (meds_surface). ----------------------!
      senv            = ground
      senv%can_temp   = tcas
      senv%can_shv    = qcas
      senv%rho_air    = cas%rho
      senv%press      = cas%press
      senv%t_soil_top = t_ground
      call ground_surface_exchange(senv, f%h_ground, f%le_ground, f%g_top, f%ground_rad, f%snow_enth_in,         &
                                   f%snow_vap_out)

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + coh_film_evap + ground%soil_evap + f%snow_vap_out
      f%coh_rnet   = coh_rnet

      !----- The condensation sink (dew, fog; meds_canopy_air_space), taken from the canopy air's vapour !
      !      and, as the condensate's liquid enthalpy, from its enthalpy. ------------------------------!
      call cas_condensation(cas%box%air_mass_capacity, y%cas_shv, qsat_c, tcas, f%cond, f%cond_enth)
      f%src_vap  = f%src_vap  - f%cond
      f%src_enth = f%src_enth - f%cond_enth

      !----- CAS box tendencies (shared kernel; the condensation adjustment above is folded into    !
      !      src_enth/src_vap, so the box math is scheme-agnostic). ---------------------------------!
      cas_env = cas%box
      cas_env%surface_enthalpy_source = f%src_enth
      cas_env%surface_vapor_source    = f%src_vap
      cas_env%biotic_co2_source       = cas%nee_biotic
      !----- The CAS <-> atmosphere conductances follow the live canopy-air state: RK45 re-solves them   !
      !      here. ESDIRK2 arrives with mo_live cleared on its stage copy, having solved once for the    !
      !      stage, so its Newton's residual evaluations reuse that value. ------------------------------!
      call cas_conductances(cas, y%cas_enthalpy, y%cas_shv, gah_l, gaw_l, gac_l)
      cas_env%atm_conductance_enthalpy = gah_l
      cas_env%atm_conductance_vapor    = gaw_l
      cas_env%atm_conductance_co2      = gac_l
      call cas_box_time_deriv(y%cas_enthalpy, y%cas_shv, y%cas_co2, cas_env,                 &
                                 f%d_cas_enthalpy, f%d_cas_shv, f%d_cas_co2)
   end subroutine surface_derivs

   !----- Allocate a surface-tendency record's per-cohort arrays for n cohorts, unless they are. -!
   pure subroutine ensure_surface_tend(f, n)
      type(surface_tend_t), intent(inout) :: f
      integer(ik),          intent(in)    :: n
      if (allocated(f%leaf_temp) .and. allocated(f%wood_temp) .and. allocated(f%transp_c) .and.   &
          allocated(f%film_evap_leaf) .and. allocated(f%film_evap_wood)) then
         if (size(f%leaf_temp) == n .and. size(f%wood_temp) == n .and. size(f%transp_c) == n .and. &
             size(f%film_evap_leaf) == n .and. size(f%film_evap_wood) == n) return
      end if
      if (allocated(f%leaf_temp))      deallocate(f%leaf_temp)
      if (allocated(f%wood_temp))      deallocate(f%wood_temp)
      if (allocated(f%transp_c))       deallocate(f%transp_c)
      if (allocated(f%film_evap_leaf)) deallocate(f%film_evap_leaf)
      if (allocated(f%film_evap_wood)) deallocate(f%film_evap_wood)
      allocate(f%leaf_temp(n), f%wood_temp(n), f%transp_c(n), f%film_evap_leaf(n), f%film_evap_wood(n))
   end subroutine ensure_surface_tend

   !---------------------------------------------------------------------------------------!
   ! march_derivs -- the whole-column RHS an RK45 stage evaluates; commits nothing. Diagnoses the    !
   ! soil-top temperature from the state, runs the surface block, and assembles the soil-heat, soil-  !
   ! water and plant-water-mass tendencies. Water moves at the soil-plant water network's time-mean   !
   ! rates and the leaves transpire what the network gave them, so every stage sees the same water    !
   ! (docs/science/soil_plant_water.md §4, docs/science/numerical_scheme.md §5b).                     !
   !---------------------------------------------------------------------------------------!
   pure subroutine march_derivs(y, frozen, site, n, nsl, f, sf_out)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      integer(ik),           intent(in)  :: n      !< number of cohorts
      integer(ik),           intent(in)  :: nsl    !< number of active soil layers
      type(march_tend_t),   intent(out) :: f
      !----- This stage's own surface tendencies, for RK45's whole-column boundary-flux ledger (as   !
      !      backward_euler_stage's surf_tend serves ESDIRK2's). -------------------------------------!
      type(surface_tend_t),  intent(out), optional :: sf_out

      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_state_t)  :: soil_energy
      type(soil_heat_env_t)     :: eforc
      real(wp)                   :: t_ground, fliq1, transp_i, e_drain, water_heat(nsl)
      integer(ik)                :: i

      allocate(f%d_leaf_water_mass(n), f%d_wood_water_mass(n), f%leaf_temp(n))
      allocate(f%d_leaf_surf_water(n), f%d_wood_surf_water(n))

      !----- The soil-top temperature from the current state, for the ground skin. -----------------!
      call soil_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                &
                           site%soil_profile%layer_curve(1), t_ground, fliq1)

      !----- 1. Surface block (leaf + ground + CAS twins). ------------------------------------!
      y_stage%cas_enthalpy = y%cas_enthalpy
      y_stage%cas_shv = y%cas_shv
      y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground,                 &
                          t_ground, n, surf_tend)
      f%d_cas_enthalpy = surf_tend%d_cas_enthalpy
      f%d_cas_shv = surf_tend%d_cas_shv
      f%d_cas_co2 = surf_tend%d_cas_co2
      f%g_top = surf_tend%g_top
      f%leaf_temp(1:n) = surf_tend%leaf_temp(1:n)

      !----- 2. Soil water at the network's steady rate, which any consistent tableau integrates      !
      !      exactly, so RK45's theta ends where the network's does; so do drainage and uptake. ----!
      f%dtheta_dt(1:nsl) = frozen%hydrology%theta_rate(1:nsl)
      f%drainage_rate    = frozen%hydrology%drainage
      f%uptake_rate      = frozen%roots%uptake

      !----- 3. Soil heat: g_top from the surface, and the enthalpy of the network's flows, the numbers !
      !      the ESDIRK2 stage uses (docs/science/numerical_scheme.md §5b, soil_plant_water.md §5). The !
      !      leaf pays its transpiration's full vapour enthalpy in surface_derivs. ------------------!
      soil_energy%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      call set_soil_heat_env(nsl, surf_tend%g_top, y%theta, frozen%roots%root_heat, frozen%hydrology%exfil_enth,  &
                             frozen%hydrology%w_flux_frozen, frozen%hydrology%infiltration,                       &
                             frozen%hydrology%t_infil, f%drainage_rate, frozen%hydrology%t_bot, eforc,            &
                             water_heat, e_drain)
      call soil_energy_time_deriv(soil_energy, eforc, site%soil_profile, site%soil,   &
                                  f%dedt, bottom_face=f%g_bottom)

      !----- 4. Plant water mass per plant: the network's sapflow and uptake in, and out the          !
      !      transpiration this stage charged the canopy air, the network's own when the stages hold it. !
      do i = 1_ik, n
         transp_i = surf_tend%transp_c(i) / max(frozen%plant%nplant(i), tiny_num)
         f%d_leaf_water_mass(i) = frozen%plant%sapflow_frozen(i) - transp_i
         f%d_wood_water_mass(i) = frozen%plant%uptake_frozen(i)  - frozen%plant%sapflow_frozen(i)
      end do

      !----- 5. Canopy film: frozen interception in, per-stage film evaporation out, both per m2 of   !
      !      ground. Zero with canopy_water_on off. -------------------------------------------------!
      f%d_leaf_surf_water(1:n) = frozen%film%intercept_leaf(1:n) - surf_tend%film_evap_leaf(1:n)
      f%d_wood_surf_water(1:n) = frozen%film%intercept_wood(1:n) - surf_tend%film_evap_wood(1:n)
      if (present(sf_out)) sf_out = surf_tend
   end subroutine march_derivs

   !=======================================================================================!
   !  The soil-plant water network (docs/science/soil_plant_water.md §1-§2): the pond, the   !
   !  soil layers, and each cohort's wood and leaves. The stores' curves and the paths'      !
   !  flows are the process libraries'; these join them.                                    !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! set_water_state -- the network's state from the model's: each soil layer's water      !
   ! content and ice, the pond's water and enthalpy, and each cohort's leaf and wood water !
   ! [kg per plant]; each store's unknown is read from its water. A leafless cohort's      !
   ! leaves hold none, so their unknown starts at the wood's potential less the head to    !
   ! the crown, where an empty leaf comes to rest.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine set_water_state(static, profile, hyd, theta, ice, w_pond, e_pond, leaf_water, wood_water, water)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(hydraulics_table_t),  intent(in)    :: hyd
      real(wp),                  intent(in)    :: theta(:)                       !< (nsl) [m3/m3] liquid and ice
      type(soil_ice_t),          intent(in)    :: ice(:)                         !< (nsl) each layer's ice
      real(wp),                  intent(in)    :: w_pond, e_pond                 !< [kg/m2], [J/m2]
      real(wp),                  intent(in)    :: leaf_water(:), wood_water(:)   !< (ncoh) [kg/plant]
      type(water_work_t),        intent(inout) :: water
      integer(ik) :: k, i
      real(wp)    :: psi, dpsi, w, dw, kk, dk, head1
      call ensure_water_capacity(water, profile%n_active, static%n)
      associate (st => water%state)
         do k = 1_ik, profile%n_active
            st%tau_soil(1_ik + k)   = soil_theta_to_unknown(profile%layer_curve(k), ice(k), theta(k))
            st%water_soil(1_ik + k) = rho_h2o * profile%dz(k) * theta(k)
         end do
         !----- The pond starts at its water's head, but no higher than the head at which it neither   !
         !      gives nor takes water from the top layer: over dry soil a trace pond empties within the  !
         !      step, and its own head would start Newton at a flow thousands of times what it holds. --!
         st%water_soil(POND) = w_pond
         st%e_pond = e_pond
         call soil_water_relations(profile%layer_curve(1), ice(1), st%tau_soil(2), psi, dpsi, w, dw, kk, dk)
         head1 = psi / grav_head
         st%tau_soil(POND) = min(pond_water_to_head(w_pond), head1 + profile%z_node(1))
         do i = 1_ik, static%n
            associate (hp => hyd%pft(static%pft(i)))
               st%tau_coh(WOOD, i)   = tissue_rwc_to_unknown(                                                  &
                                          tissue_water_to_rwc(hp%wood_curve, static%bwood_water(i), wood_water(i)))
               st%water_coh(WOOD, i) = wood_water(i) * static%nplant(i)
               if (hp%leaf_curve%water_sat * static%bleaf(i) > tiny_num) then
                  st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(                                               &
                                             tissue_water_to_rwc(hp%leaf_curve, static%bleaf(i), leaf_water(i)))
                  st%water_coh(LEAF, i) = leaf_water(i) * static%nplant(i)
               else
                  call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i),            &
                                              st%tau_coh(WOOD, i), psi, dpsi, w, dw)
                  st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(tissue_psi_to_rwc(psi - static%crown_head(i),   &
                                             hp%leaf_curve%pi0, hp%leaf_curve%elastic_mod))
                  st%water_coh(LEAF, i) = 0.0_wp
               end if
            end associate
         end do
      end associate
   end subroutine set_water_state

   !---------------------------------------------------------------------------------------!
   ! water_derivs -- the network at one state: each store's potential and water from its     !
   ! curve; the flow along every path (pond to top layer, between layers, evaporation,       !
   ! drainage, each cohort's uptake, sapflow and transpiration); the net inflow into each    !
   ! store; and, when asked, the flows' Jacobian in the stores' unknowns (sys).              !
   !---------------------------------------------------------------------------------------!
   subroutine water_derivs(static, profile, soil, hyd, drive, tau_soil, tau_coh, ev, sys, with_jacobian)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(in)    :: drive
      real(wp),                  intent(in)    :: tau_soil(:), tau_coh(:,:)
      type(water_tend_t),        intent(inout) :: ev
      type(arrowhead_system_t),  intent(inout) :: sys   !< the flows' Jacobian, written when with_jacobian
      logical,                   intent(in)    :: with_jacobian
      integer(ik) :: k, i, ca, cb, c, nsl, nsoil
      real(wp) :: q, dq_pa, dq_pb, dq_ka, dq_kb, dqa, dqb, u, s, ds_w, ds_l
      real(wp) :: e, de_dtheta, de_dpsi, de, qbot, dqbot_dpsi, dqbot_dk, theta1
      real(wp) :: psi_w, dpsi_w, u_sum, dqb_sum, du_dpsi(profile%n_active), du_dk(profile%n_active)

      nsl   = profile%n_active
      nsoil = 1_ik + nsl

      ev%inflow_soil = 0.0_wp
      ev%inflow_coh = 0.0_wp
      !----- The Jacobian's entries are written afresh: the soil block's band, which the flows add  !
      !      into, is cleared; every other entry an evaluation touches is set outright below, and the  !
      !      rest stay zero from arrowhead_allocate. ---------------------------------------------------!
      if (with_jacobian) then
         do c = 1_ik, nsoil
            sys%a(c, c) = 0.0_wp
            if (c < nsoil) then
               sys%a(c, c + 1_ik) = 0.0_wp
               sys%a(c + 1_ik, c) = 0.0_wp
            end if
         end do
      end if

      !----- Storage. ------------------------------------------------------------------!
      call pond_water_relations(tau_soil(POND), ev%psi_soil(1), ev%dpsi_soil(1), ev%water_soil(1),             &
                                ev%dwater_soil(1))
      !----- The pond conducts as saturated soil; its path lies in the top layer, so it takes that  !
      !      layer's ice impedance in full. --------------------------------------------------------!
      ev%k_soil(1) = profile%ksat(1) * drive%ice(1)%impedance
      ev%dk_soil(1) = 0.0_wp
      do k = 1_ik, nsl
         c = 1_ik + k
         call soil_water_relations(profile%layer_curve(k), drive%ice(k), tau_soil(c), ev%psi_soil(c),             &
                                   ev%dpsi_soil(c), ev%water_soil(c), ev%dwater_soil(c), ev%k_soil(c),            &
                                   ev%dk_soil(c))
      end do
      do i = 1_ik, static%n
         associate (hp => hyd%pft(static%pft(i)))
            call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i), tau_coh(WOOD, i), &
                                        ev%psi_coh(WOOD, i), ev%dpsi_coh(WOOD, i), ev%water_coh(WOOD, i),        &
                                        ev%dwater_coh(WOOD, i))
            call tissue_water_relations(hp%leaf_curve, static%bleaf(i), static%nplant(i), tau_coh(LEAF, i),       &
                                        ev%psi_coh(LEAF, i), ev%dpsi_coh(LEAF, i), ev%water_coh(LEAF, i),        &
                                        ev%dwater_coh(LEAF, i))
         end associate
      end do

      !----- The pond and the top layer: Darcy, the pond a saturated layer as thick as the top one. !
      ca = POND
      cb = 2_ik
      call soil_water_flux(ev%k_soil(ca), ev%k_soil(cb), 0.5_wp, ev%psi_soil(ca), ev%psi_soil(cb),                &
                           -profile%z_node(1), -profile%z_node(1), q, dq_pa, dq_pb, dq_ka, dq_kb)
      ev%pond_flux = q
      call add_soil_flow(ca, cb, q, dq_pa * ev%dpsi_soil(ca), dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb))

      !----- Between soil layers. --------------------------------------------------------!
      do k = 1_ik, nsl - 1_ik
         ca = 1_ik + k
         cb = 2_ik + k
         call soil_water_flux(ev%k_soil(ca), ev%k_soil(cb), profile%face_weight(k), ev%psi_soil(ca),              &
                              ev%psi_soil(cb), profile%dz_node(k), profile%dz_node(k), q, dq_pa, dq_pb, dq_ka, dq_kb)
         ev%face_flux(k) = q
         dqa = dq_pa * ev%dpsi_soil(ca) + dq_ka * ev%dk_soil(ca)
         dqb = dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb)
         call add_soil_flow(ca, cb, q, dqa, dqb)
      end do

      !----- Rain into the pond; evaporation from the top layer; drainage from the bottom one. The  !
      !      dry surface layer reads the liquid and dsl_theta_init of the ice, so its onset is      !
      !      dsl_theta_init of the pores the ice leaves (CLM5's effective porosity). ---------------!
      ev%inflow_soil(POND) = ev%inflow_soil(POND) + drive%rain
      c = 2_ik
      theta1 = ev%water_soil(c) / (rho_h2o * profile%dz(1)) - (1.0_wp - soil%dsl_theta_init) * drive%ice(1)%theta_ice
      call ground_evaporation_rate(theta1, ev%psi_soil(c), drive%ground, soil, e, de_dtheta, de_dpsi)
      ev%evap = e
      ev%inflow_soil(c) = ev%inflow_soil(c) - e
      if (with_jacobian) sys%a(c, c) = sys%a(c, c)                                                              &
                         + de_dtheta * ev%dwater_soil(c) / (rho_h2o * profile%dz(1))                               &
                         + de_dpsi * ev%dpsi_soil(c)
      c = 1_ik + nsl
      call bottom_flux(profile, soil, nsl, ev%psi_soil(c), ev%k_soil(c), drive%ice(nsl)%impedance, qbot,           &
                       dqbot_dpsi, dqbot_dk)
      ev%drain = rho_h2o * qbot
      ev%inflow_soil(c) = ev%inflow_soil(c) - ev%drain
      if (with_jacobian) sys%a(c, c) = sys%a(c, c)                                                              &
                         + rho_h2o * (dqbot_dpsi * ev%dpsi_soil(c) + dqbot_dk * ev%dk_soil(c))

      !----- Each cohort: uptake from its layers, sapflow to its leaves, and transpiration. The      !
      !      cohort's uptake comes from its whole root profile at once; its own sums are kept in       !
      !      scalars over its layers and written once. ------------------------------------------!
      do i = 1_ik, static%n
         psi_w = ev%psi_coh(WOOD, i)
         dpsi_w = ev%dpsi_coh(WOOD, i)
         call root_water_uptake_profile(nsl, static%g_rhizo(1:nsl, i), ev%k_soil(2:nsoil), ev%psi_soil(2:nsoil), &
                                        profile%z_node(1:nsl), psi_w, ev%uptake(1:nsl, i), du_dpsi, du_dk)
         u_sum = 0.0_wp
         dqb_sum = 0.0_wp
         do k = 1_ik, nsl
            c = 1_ik + k
            if (static%g_rhizo(k, i) <= 0.0_wp) then
               if (with_jacobian) then
                  sys%b(c, i) = 0.0_wp
                  sys%c(c, i) = 0.0_wp
               end if
               cycle
            end if
            u = ev%uptake(k, i)
            ev%inflow_soil(c) = ev%inflow_soil(c) - u
            u_sum = u_sum + u
            if (with_jacobian) then
               dqa = du_dpsi(k) * ev%dpsi_soil(c) + du_dk(k) * ev%dk_soil(c)
               dqb = -du_dpsi(k) * dpsi_w
               sys%a(c, c) = sys%a(c, c) + dqa
               sys%b(c, i) = dqb
               sys%c(c, i) = -dqa
               dqb_sum = dqb_sum - dqb
            end if
         end do
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) + u_sum
         associate (hp => hyd%pft(static%pft(i)))
            call xylem_sapflow(static%k_xylem(i), ev%psi_coh(WOOD, i), ev%psi_coh(LEAF, i), static%crown_head(i),  &
                               hp%wood_psi50, hp%wood_kexp, hp%vuln_table, s, ds_w, ds_l)
         end associate
         ev%sapflow(i) = s
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) - s
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) + s
         call supply_limited_transpiration(drive%e_demand(i), tau_coh(LEAF, i), hyd%pft(static%pft(i))%rwc_close,   &
                                           e, de)
         ev%transp(i) = e
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) - e
         if (with_jacobian) then
            dqa = ds_w * ev%dpsi_coh(WOOD, i)
            dqb = ds_l * ev%dpsi_coh(LEAF, i)
            sys%d(WOOD, WOOD, i) = dqb_sum + dqa
            sys%d(WOOD, LEAF, i) = dqb
            sys%d(LEAF, WOOD, i) = -dqa
            sys%d(LEAF, LEAF, i) = -dqb + de
         end if
      end do

   contains

      !----- A flow q from soil node a to soil node b, with dq/dtau_a and dq/dtau_b: the inflows   !
      !      and, when wanted, the Jacobian of the residual (W - W0)/h - inflow. ----------------!
      subroutine add_soil_flow(a, b, q, dqa, dqb)
         integer(ik), intent(in) :: a, b
         real(wp),    intent(in) :: q, dqa, dqb
         ev%inflow_soil(a) = ev%inflow_soil(a) - q
         ev%inflow_soil(b) = ev%inflow_soil(b) + q
         if (.not. with_jacobian) return
         sys%a(a, a) = sys%a(a, a) + dqa
         sys%a(a, b) = sys%a(a, b) + dqb
         sys%a(b, a) = sys%a(b, a) - dqa
         sys%a(b, b) = sys%a(b, b) - dqb
      end subroutine add_soil_flow

   end subroutine water_derivs

end module meds_fast_integ_derivs
