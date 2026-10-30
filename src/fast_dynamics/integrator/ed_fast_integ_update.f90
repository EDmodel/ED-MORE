! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_update -- the end of a fast step, the same for both schemes.
! update_dynamic_state commits the march's end state, the water network's soil water and
! pond, the tissues' energy, the condensate and the snow pack: the one place a step writes
! the dynamic state. close_ledgers closes the seven conservation ledgers on the march's
! boundary-flux amounts (docs/science/numerical_scheme.md §5b). Also the commit's parts
! and the store totals the ledgers read.
!---~---
module ed_fast_integ_update
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : rho_h2o, tiny_num, mmh2o, latent_heat_vap
   use ed_config,             only : ed_config_t
   use ed_site_state_vars,    only : site_surface_t
   use ed_fast_integ_types,       only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,     &
                                     step_report_t, march_state_t, frozen_t, boundary_flux_t,              &
                                     diagnose_soil_temps
   use ed_plant_energy,       only : commit_tissue_energy
   use ed_canopy_water,       only : clamp_surface_water
   use ed_patch_state_vars,   only : soil_water_store, soil_energy_store
   use ed_diag_vars,          only : CO_TRANSP
   use ed_surface, only : snow_take_exchange
   use ed_therm_lib,          only : cas_enthalpy_to_temp
   use ed_budget_check,       only : budget_check, BUDGET_ENERGY_RATE_FLOOR, BUDGET_WATER_RATE_FLOOR,           &
                                     BUDGET_CO2_RATE_FLOOR
   implicit none
   private

   public :: commit_amounts_t, update_dynamic_state, close_ledgers
   public :: unpack_column_state, transp_per_leaf_area   ! used only by the tests

   !---~---
   !   What the update moved besides the march, which the ledgers book: the films' overflow at
   ! their capacity and deficit below zero, and their net enthalpy at the tissues'
   ! start-of-step temperatures, the condensate deposited in the top soil layer with its
   ! enthalpy, the tissue store (films included) before and after, and the snow pack after
   ! its exchange.
   !---~---
   type :: commit_amounts_t
      real(wp) :: surf_overflow = 0.0_wp, surf_deficit = 0.0_wp   ! [kg/m2]
      real(wp) :: surf_net_enth = 0.0_wp                          ! [J/m2] (overflow - deficit) x u_liq
      real(wp) :: cond_mass     = 0.0_wp                          ! [kg/m2]
      real(wp) :: cond_enth     = 0.0_wp                          ! [J/m2]
      real(wp) :: tissue_store0 = 0.0_wp, tissue_store1 = 0.0_wp  ! [J/m2]
      real(wp) :: snow_swe1     = 0.0_wp                          ! [kg/m2] the pack's mass at the step's end
      real(wp) :: snow_enth1    = 0.0_wp                          ! [J/m2]  its internal energy then
   end type commit_amounts_t

contains

   !---~---
   !   update_dynamic_state -- commit one fast step. Soil water and the pond are the water
   ! network's: the march moved theta at the network's steady rate, and its end-of-step theta
   ! and pond are committed once, here. A film evaporating at its frozen conductance can
   ! overdraw below zero, which would create water in the next step's interception, so the
   ! films are floored at zero and capped at capacity, and the update reports both. Each
   ! tissue's energy is committed from its time-average temperature over the march, the
   ! complement of what the canopy air received (docs/science/vegetation_energy_dynamics.md
   ! §One formulation). The condensate, dew and fog, lands inside the column, so it goes,
   ! mass and enthalpy paired, into soil layer 1 (docs/science/numerical_scheme.md §7, item 3)
   ! before the soil's temperatures are diagnosed.
   !---~---
   subroutine update_dynamic_state(dt_fast, cfg, site, static, frozen, acc, dt_warm, y_out, dynamic, out, moved)
      real(wp),                   intent(in)    :: dt_fast
      type(ed_config_t),          intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(frozen_t),      intent(in)    :: frozen
      type(boundary_flux_t),       intent(in)    :: acc       ! the march's boundary-flux amounts
      real(wp),                   intent(in)    :: dt_warm   ! [s] the step the next march starts from
      type(march_state_t),       intent(inout) :: y_out     ! the march's end state, committed as stored
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      type(fast_output_t),        intent(inout) :: out
      type(commit_amounts_t),     intent(out)   :: moved
      real(wp)    :: t_leaf_avg(static%n), t_wood_avg(static%n), leaf_surf_water0(static%n), wood_surf_water0(static%n)
      integer(ik) :: n, nsl

      n   = static%n
      leaf_surf_water0(1:n) = dynamic%leaf_surf_water(1:n)              ! the films at the step's start
      wood_surf_water0(1:n) = dynamic%wood_surf_water(1:n)
      nsl = site%soil_profile%n_active
      y_out%theta(1:nsl)     = frozen%hydrology%theta1(1:nsl)
      y_out%w_surface        = frozen%hydrology%w_surface1
      y_out%w_surface_energy = frozen%hydrology%w_surface_energy1
      if (cfg%canopy_water_on) then
         call clamp_surface_water(y_out%leaf_surf_water, y_out%wood_surf_water, static%lai, static%wai,             &
                                cfg%pft%leaf_surf_water_max(static%pft(1:n)),                                  &
                                cfg%pft%wood_surf_water_max(static%pft(1:n)), frozen%tissue%t_leaf0,           &
                                frozen%tissue%t_wood0, n, moved%surf_overflow, moved%surf_deficit,             &
                                moved%surf_net_enth)
      end if
      !---~---
      !   The condensate lands before the soil's temperatures are diagnosed, so the top layer's
      ! temperature includes its water and enthalpy.
      !---~---
      if (acc%whole_cond > 0.0_wp) then
         moved%cond_mass = acc%whole_cond
         moved%cond_enth = acc%whole_cond_enth       ! b-weighted at the stage CAS temperatures (boundary_flux_t)
         call deposit_condensate(y_out, site%soil_profile%dz(1), moved%cond_mass, moved%cond_enth)
      end if

      call unpack_column_state(y_out, n, nsl, dynamic)
      dynamic%soil_water%w_surface        = y_out%w_surface
      dynamic%soil_water%w_surface_energy = y_out%w_surface_energy
      call diagnose_soil_temps(y_out, site%soil_profile, nsl, dynamic%soil_energy%soil_temp,                    &
                               dynamic%soil_energy%soil_fliq)
      t_leaf_avg(1:n) = acc%tissue_leaf_int(1:n) / dt_fast
      t_wood_avg(1:n) = acc%tissue_wood_int(1:n) / dt_fast
      call commit_tissue_energy(frozen%tissue, static%nplant, dt_fast, t_leaf_avg, t_wood_avg, n,                 &
                                dynamic%leaf_water_mass, dynamic%wood_water_mass, leaf_surf_water0, wood_surf_water0,         &
                                dynamic%leaf_surf_water, dynamic%wood_surf_water, dynamic%leaf_energy,            &
                                dynamic%wood_energy, dynamic%leaf_temp, dynamic%wood_temp, moved%tissue_store0,   &
                                moved%tissue_store1)
      !--- Each cohort's transpiration over the step, as the leaf water was debited.
      if (out%capture .and. allocated(acc%transp_coh))                                                       &
         out%cdiag(CO_TRANSP, 1:n) = transp_per_leaf_area(acc%transp_coh(1:n), dt_fast, static%nplant(1:n),     &
                                                          static%leaf_area(1:n))
      !---~---
      !   The pack takes its exchange over the march, the energy and vapour the canopy air and the
      ! soil received from it (docs/science/snow_biophysics.md §The snow store); and the warm start
      ! of the next march.
      !---~---
      dynamic%snow = frozen%ground%pack
      call snow_take_exchange(dynamic%snow, acc%snow_enth_in, acc%snow_vap_out, frozen%ground%params)
      moved%snow_swe1       = dynamic%snow%swe(1)
      moved%snow_enth1      = dynamic%snow%snow_energy(1)
      dynamic%adapt_dt_last = dt_warm
   end subroutine update_dynamic_state

   !---~---
   !   close_ledgers -- the seven conservation ledgers, from the march's b-weighted boundary-flux
   ! amounts, which close to round-off (docs/science/numerical_scheme.md §5b): the canopy air's
   ! energy, water and CO2, the soil's energy and water, and the whole column's water and
   ! energy; and the step's canopy-air -> atmosphere exports. y is the march's start state and
   ! y_out the committed end. With halt ([fast].debug_error) a budget that does not close stops
   ! the run; the tolerance scales with the boundary flux, not the store.
   !---~---
   subroutine close_ledgers(dt_fast, site, forcing, static, frozen, y, y_out, acc, moved, halt, report)
      real(wp),                  intent(in)    :: dt_fast
      type(site_surface_t),      intent(in)    :: site
      type(fast_forcing_t),      intent(in)    :: forcing
      type(fast_static_state_t), intent(in)    :: static
      type(frozen_t),     intent(in)    :: frozen
      type(march_state_t),      intent(in)    :: y, y_out
      type(boundary_flux_t),      intent(in)    :: acc
      type(commit_amounts_t),    intent(in)    :: moved
      logical,                   intent(in)    :: halt
      type(step_report_t),       intent(inout) :: report
      real(wp)    :: cas_mass_capacity, cas_molar_capacity, enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, w_plant0, w_plant1, w_surface0, e_pond0, e_pond1
      real(wp)    :: surf_water0, surf_water1, intercept_total
      integer(ik) :: n, nsl

      n   = static%n
      nsl = site%soil_profile%n_active
      cas_mass_capacity  = frozen%cas%box%air_mass_capacity
      cas_molar_capacity = frozen%cas%box%air_molar_capacity
      enth0 = y%cas_enthalpy
      shv0  = y%cas_shv
      co20  = y%cas_co2
      enth1 = y_out%cas_enthalpy
      shv1  = y_out%cas_shv
      co21  = y_out%cas_co2
      !--- The pond at the start, as the march started from it, and at the end, the network's.
      w_surface0 = y%w_surface
      e_pond0    = y%w_surface_energy
      e_pond1    = frozen%hydrology%w_surface_energy1
      e_soil0 = soil_energy_store(y%soil_energy,     site%soil_profile%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, site%soil_profile%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     site%soil_profile%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, site%soil_profile%dz, nsl)
      !---~---
      !   Plant water is a store: it absorbs the step's mismatch between transpiration and uptake
      ! (docs/science/plant_hydraulics.md eqs. (1)-(2)).
      !---~---
      w_plant0 = plant_water_store(static%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(static%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !---~---
      !   The films, per m2 of ground (0 with canopy_water_on off). Their energy is their
      ! tissues', in the tissue store.
      !---~---
      surf_water0 = surface_water_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = surface_water_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      intercept_total = sum(frozen%veg_surface%intercept_leaf(1:n) + frozen%veg_surface%intercept_wood(1:n))

      report%exports%atm_heat_export = acc%atm_heat_out
      report%exports%atm_vap_export  = acc%atm_vap_out
      !---~---
      !   The vapour sources the canopy air was credited (surface_derivs' src_vap): transpiration,
      ! the films', the ground's at the network's rate held over the step, and the snow pack's;
      ! and the condensate it was debited.
      !---~---
      report%exports%transp_export      = 0.0_wp
      if (allocated(acc%transp_coh)) report%exports%transp_export = sum(acc%transp_coh(1:n))
      report%exports%veg_surface_evap_export   = acc%veg_surface_evap
      report%exports%ground_evap_export = frozen%ground%soil_evap * dt_fast + acc%snow_vap_out
      report%exports%cond_export        = acc%whole_cond
      !---~---
      !   The budgets' other terms, so both close from the output (docs/science/diagnostics.md §4,
      ! The budgets): the same amounts the two whole-column ledgers below set against their
      ! stores. The water's enthalpy is split by direction: the precipitation's in (the rain the
      ! surfaces catch, the snowfall, the rain reaching the ground), and the outflows' out (the
      ! drainage, the runoff, the surface water's overflow, and the exported vapour's beyond LE).
      !---~---
      report%exports%intercept_export    = intercept_total * dt_fast
      report%exports%infiltration_export = frozen%hydrology%infiltration * dt_fast
      report%exports%melt_export         = frozen%snow%melt_rate * dt_fast
      report%exports%sublimation_export  = acc%snow_vap_out
      report%exports%shed_water_export   = static%shed_water_rate * dt_fast
      report%exports%veg_surface_overflow_export  = moved%surf_overflow - moved%surf_deficit
      report%exports%water_store_gain    = (w_soil1 + cas_mass_capacity*shv1 + frozen%hydrology%w_surface1 + w_plant1 &
                                           + surf_water1 + moved%snow_swe1)                                       &
                                           - (w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0 &
                                           + frozen%snow%swe0)
      report%exports%veg_surface_water_store          = surf_water1
      report%exports%soil_water_store    = w_soil1
      report%exports%plant_water_store   = w_plant1
      report%exports%heat_bottom_export  = acc%whole_enth_in - acc%rad_net
      report%exports%precip_enth_export  = intercept_total*dt_fast*frozen%veg_surface%rain_liquid_enthalpy                  &
                                           + frozen%snow%acc_enth + frozen%hydrology%ground_inflow_enth*dt_fast
      report%exports%outflow_enth_export = acc%whole_enth_out - acc%atm_heat_out - latent_heat_vap*acc%atm_vap_out  &
                                           + frozen%hydrology%runoff_enth*dt_fast + moved%surf_net_enth
      report%exports%heat_gain_soil      = e_soil1 - e_soil0
      report%exports%heat_gain_canopy    = cas_mass_capacity*(enth1 - enth0) + moved%tissue_store1 - moved%tissue_store0
      report%exports%heat_gain_surface   = moved%snow_enth1 - frozen%snow%enth0 + e_pond1 - e_pond0
      report%exports%h_canopy_export     = acc%h_canopy
      report%exports%h_ground_export     = acc%h_ground
      call budget_check(report%ledgers%cas_energy,                                                           &
                        cas_mass_capacity*enth0, cas_mass_capacity*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'cas_energy', halt)
      call budget_check(report%ledgers%cas_water,                                                            &
                        cas_mass_capacity*shv0,  cas_mass_capacity*shv1,  acc%cas_vap_in,  acc%cas_vap_out,  &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'cas_water', halt)
      call budget_check(report%ledgers%cas_co2,                                                              &
                        cas_molar_capacity*co20, cas_molar_capacity*co21, acc%cas_co2_in,  acc%cas_co2_out,  &
                        dt_fast, BUDGET_CO2_RATE_FLOOR, 'cas_co2', halt)
      !---~---
      !   The condensate deposit is an inflow to the soil budgets, though the whole-column ledger
      ! sees it as a transfer between two of its stores.
      !---~---
      call budget_check(report%ledgers%soil_energy, e_soil0, e_soil1, acc%soil_enth_in + moved%cond_enth,          &
                        acc%soil_enth_out, dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'soil_energy', halt)
      !---~---
      !   Soil water: theta^n -> theta1, in through the top face net of evaporation and
      ! exfiltration (q_top), out through the bottom and to the roots, plus the condensate.
      !---~---
      call budget_check(report%ledgers%soil_water,  w_soil0, w_soil1,                                              &
                        frozen%hydrology%q_top*rho_h2o*dt_fast + moved%cond_mass,                          &
                        (frozen%hydrology%drainage + frozen%roots%uptake)*dt_fast,                         &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'soil_water', halt)
      !---~---
      !   Whole water: rainfall, snowfall and shed water in; the CAS vapour export, drainage and
      ! runoff out; soil, pond, plant water, films, snow pack and canopy-air vapour stored.
      ! Sublimation leaves through the vapour export and meltwater moves pack -> pond, so
      ! neither needs a term.
      !---~---
      call budget_check(report%ledgers%whole_water,                                                                &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0 + frozen%snow%swe0,  &
                        w_soil1 + cas_mass_capacity*shv1 + frozen%hydrology%w_surface1 + w_plant1 + surf_water1 &
                        + moved%snow_swe1,                                                                 &
                        acc%whole_wat_in + (forcing%rainfall + forcing%snowfall + static%shed_water_rate)*dt_fast, &
                        acc%whole_wat_out + (frozen%hydrology%runoff_surf + frozen%hydrology%drainage)*dt_fast  &
                        + moved%surf_overflow - moved%surf_deficit,                                        &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'whole_water', halt)
      !---~---
      !   Whole energy: meltwater passes pack -> pond, both stores, so it telescopes. The
      ! films' capture comes in at the rain's liquid enthalpy, and the update's cut leaves at
      ! the tissues' (moved%surf_net_enth).
      !---~---
      call budget_check(report%ledgers%whole_energy,                                                               &
                        e_soil0 + cas_mass_capacity*enth0 + frozen%snow%enth0 + e_pond0 + moved%tissue_store0,      &
                        e_soil1 + cas_mass_capacity*enth1 + moved%snow_enth1 + e_pond1 + moved%tissue_store1,       &
                        acc%whole_enth_in + intercept_total*dt_fast*frozen%veg_surface%rain_liquid_enthalpy               &
                                          + frozen%snow%acc_enth + frozen%hydrology%ground_inflow_enth*dt_fast, &
                        acc%whole_enth_out + moved%surf_net_enth + frozen%hydrology%runoff_enth*dt_fast,         &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'whole_energy', halt)
   end subroutine close_ledgers

   !--- The commit's parts and the store totals the ledgers read.

   !---~---
   !   A cohort's transpiration over a step [kg/m2 ground] as the rate per leaf area the
   ! diagnostics report [mol H2O/m2 leaf/s]; zero for a leafless cohort.
   !---~---
   elemental pure function transp_per_leaf_area(amount, dt, nplant, leaf_area) result(rate)
      real(wp), intent(in) :: amount      ! [kg/m2 ground] over dt
      real(wp), intent(in) :: dt          ! [s]
      real(wp), intent(in) :: nplant      ! [plant/m2]
      real(wp), intent(in) :: leaf_area   ! [m2 leaf/plant]
      real(wp) :: rate
      rate = 0.0_wp
      if (nplant * leaf_area > tiny_num) rate = amount / dt / (nplant * leaf_area) / mmh2o
   end function transp_per_leaf_area

   !---~---
   !   The plant water's and the canopy films' totals for the conservation ledgers [per m2
   ! ground]; the soil column's are the patch state's (soil_water_store,
   ! soil_energy_store). Shared by both schemes, so a ledger built from these is the same
   ! arithmetic on both.
   !---~---
   pure function plant_water_store(nplant, leaf_mass, wood_mass, n) result(w)   ! [kg/m2] from [kg/plant]
      real(wp),    intent(in) :: nplant(:), leaf_mass(:), wood_mass(:)
      integer(ik), intent(in) :: n
      real(wp) :: w
      w = sum(nplant(1:n) * (leaf_mass(1:n) + wood_mass(1:n)))
   end function plant_water_store

   pure function surface_water_store(leaf_surf, wood_surf, n) result(w)   ! [kg/m2 ground], already ground-referenced
      real(wp),    intent(in) :: leaf_surf(:), wood_surf(:)
      integer(ik), intent(in) :: n
      real(wp) :: w
      w = sum(leaf_surf(1:n) + wood_surf(1:n))
   end function surface_water_store

   !---~---
   !   Deposit condensed canopy-air vapour into soil layer 1 as paired mass [kg/m2] + liquid
   ! enthalpy [J/m2] (an internal CAS -> soil transfer; the ledger needs no boundary term).
   !---~---
   pure subroutine deposit_condensate(y_out, dz1, mass, enth)
      type(march_state_t), intent(inout) :: y_out
      real(wp),             intent(in)    :: dz1, mass, enth
      y_out%theta(1)       = y_out%theta(1) + mass / (rho_h2o * dz1)
      y_out%soil_energy(1) = y_out%soil_energy(1) + enth / dz1
   end subroutine deposit_condensate

   !---~---
   !   Unpack the committed column state into the dynamic state: CAS twins (and their
   ! temperature), soil energy and moisture, plant internal water, canopy films; not the pond.
   !---~---
   pure subroutine unpack_column_state(y_out, n, nsl, dynamic)
      type(march_state_t),       intent(in)    :: y_out
      integer(ik),                intent(in)    :: n, nsl
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      dynamic%cas%can_enthalpy = y_out%cas_enthalpy
      dynamic%cas%can_shv = y_out%cas_shv
      dynamic%cas%can_co2 = y_out%cas_co2
      dynamic%cas%can_temp = cas_enthalpy_to_temp(y_out%cas_enthalpy, y_out%cas_shv)
      dynamic%soil_energy%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      dynamic%soil_water%theta(1:nsl)       = y_out%theta(1:nsl)
      dynamic%leaf_water_mass(1:n) = y_out%leaf_water_mass(1:n)
      dynamic%wood_water_mass(1:n) = y_out%wood_water_mass(1:n)
      dynamic%leaf_surf_water(1:n) = y_out%leaf_surf_water(1:n)
      dynamic%wood_surf_water(1:n) = y_out%wood_surf_water(1:n)
      !---~---
      !   The pond (w_surface, w_surface_energy) is left out on purpose: each scheme commits it from
      ! the water network, so writing it here would commit it twice (test_fast_integ_kernels).
      !---~---
   end subroutine unpack_column_state

end module ed_fast_integ_update
