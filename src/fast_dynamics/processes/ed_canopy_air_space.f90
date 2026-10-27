! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_canopy_air_space -- the canopy air space (CAS)
! (docs/science/canopy_air_space_biophysics.md): one well-mixed layer of air whose
! specific enthalpy, specific humidity and CO2 mixing ratio each change by their summed
! surface sources and their exchange with the atmosphere above. cas_box_time_deriv gives
! their rates of change; cas_box_step_implicit advances them by one backward-Euler step,
! implicit in the atmospheric exchange; cas_condensation is the condensation sink of
! supersaturated air. Both forms of the update read one record, cas_env_t.
!---~---
module ed_canopy_air_space
   use ed_kinds,              only : wp, ik
   use ed_therm_lib,          only : internal_energy_liquid
   implicit none
   private

   public :: cas_env_t
   public :: cas_box_time_deriv, cas_box_step_implicit, cas_condensation

   !---~---
   !   [s] e-folding time of the condensation sink, dq/dt = -(q - qsat)/TAU_COND: short against
   ! dt_fast, so supersaturated canopy air condenses within a step, and a smooth rate, not a
   ! clamp, so the step's error estimate stays small. Any short value pins q near qsat.
   !---~---
   real(wp), parameter :: TAU_COND = 300.0_wp

   !---~---
   !   What the canopy air reads, per patch: its capacities, the canopy-air <-> atmosphere
   ! conductances and the atmosphere's state, held fixed over the step, and the summed
   ! sources into it, condensation included, so both forms of the update see it.
   !---~---
   type :: cas_env_t
      real(wp) :: air_mass_capacity  = 0.0_wp   ! [kg/m2]  CAS air mass per ground area (for enthalpy and humidity)
      real(wp) :: air_molar_capacity = 0.0_wp   ! [mol/m2] CAS dry-air molar capacity (for CO2)
      real(wp) :: atm_conductance_enthalpy = 0.0_wp  ! [kg/m2/s]  atm<->CAS enthalpy conductance (rho*ustar*heat_transfer_factor)
      real(wp) :: atm_conductance_vapor    = 0.0_wp  ! [kg/m2/s]  atm<->CAS vapour conductance
      real(wp) :: atm_conductance_co2      = 0.0_wp  ! [mol/m2/s] atm<->CAS CO2 conductance
      real(wp) :: atm_enthalpy         = 0.0_wp ! [J/kg]      reference-level specific enthalpy
      real(wp) :: atm_specific_humidity= 0.0_wp ! [kg/kg]     reference-level specific humidity
      real(wp) :: atm_co2              = 400.0_wp! [umol/mol]  free-atmosphere CO2
      real(wp) :: surface_enthalpy_source = 0.0_wp   ! [W/m2]      sensible + latent into the CAS
      real(wp) :: surface_vapor_source    = 0.0_wp   ! [kg/m2/s]   vapour into the CAS
      real(wp) :: biotic_co2_source       = 0.0_wp   ! [umol/m2/s] respiration - GPP (Reco - GPP)
   end type cas_env_t

contains

   !---~---
   !   The rates of change of the canopy air's enthalpy, humidity and CO2 at the current state
   ! (docs/science/canopy_air_space_biophysics.md eq. (1)): the sources plus the exchange
   ! with the atmosphere, over the capacity.
   !---~---
   pure subroutine cas_box_time_deriv(cas_enthalpy, cas_shv, cas_co2, env,                &
                                         d_enthalpy, d_shv, d_co2)
      real(wp),           intent(in)  :: cas_enthalpy, cas_shv, cas_co2   ! [J/kg], [kg/kg], [umol/mol]
      type(cas_env_t),    intent(in)  :: env   ! capacities, conductances, the atmosphere and the sources
      real(wp),           intent(out) :: d_enthalpy, d_shv, d_co2   ! [J/kg/s], [kg/kg/s], [umol/mol/s]
      d_enthalpy = (env%surface_enthalpy_source                                             &
                    + env%atm_conductance_enthalpy * (env%atm_enthalpy - cas_enthalpy))  &
                   / env%air_mass_capacity
      d_shv      = (env%surface_vapor_source                                                &
                    + env%atm_conductance_vapor * (env%atm_specific_humidity - cas_shv))  &
                   / env%air_mass_capacity
      d_co2      = (env%biotic_co2_source                                                    &
                    + env%atm_conductance_co2 * (env%atm_co2 - cas_co2))                  &
                   / env%air_molar_capacity
   end subroutine cas_box_time_deriv

   !---~---
   !   The canopy air's enthalpy, humidity and CO2 after one backward-Euler step of dt
   ! (docs/science/canopy_air_space_biophysics.md eq. (2)): implicit in the atmospheric
   ! exchange, the stiff term, so the update is L-stable; the sources stay fixed over dt.
   !---~---
   pure subroutine cas_box_step_implicit(cas_enthalpy, cas_shv, cas_co2, env, dt,         &
                                            enthalpy_new, shv_new, co2_new)
      real(wp),           intent(in)  :: cas_enthalpy, cas_shv, cas_co2   ! [J/kg], [kg/kg], [umol/mol] at the start
      type(cas_env_t),    intent(in)  :: env   ! capacities, conductances, the atmosphere and the sources
      real(wp),           intent(in)  :: dt    ! [s] step length
      real(wp),           intent(out) :: enthalpy_new, shv_new, co2_new   ! [J/kg], [kg/kg], [umol/mol] at the end
      enthalpy_new = (env%air_mass_capacity * cas_enthalpy                                  &
                      + dt * (env%surface_enthalpy_source                                   &
                              + env%atm_conductance_enthalpy * env%atm_enthalpy))        &
                     / (env%air_mass_capacity + dt * env%atm_conductance_enthalpy)
      shv_new      = (env%air_mass_capacity * cas_shv                                       &
                      + dt * (env%surface_vapor_source                                      &
                              + env%atm_conductance_vapor * env%atm_specific_humidity))  &
                     / (env%air_mass_capacity + dt * env%atm_conductance_vapor)
      co2_new      = (env%air_molar_capacity * cas_co2                                      &
                      + dt * (env%biotic_co2_source                                         &
                              + env%atm_conductance_co2 * env%atm_co2))                  &
                     / (env%air_molar_capacity + dt * env%atm_conductance_co2)
   end subroutine cas_box_step_implicit

   !---~---
   !   The smooth condensation sink (dew, fog): supersaturated canopy air relaxes to saturation
   ! over TAU_COND. The condensate leaves as liquid at the air's temperature, taking only its
   ! liquid enthalpy; the latent heat stays in the air, as in real condensation.
   !---~---
   elemental subroutine cas_condensation(air_mass_capacity, shv, qsat, temp, cond, cond_enth)
      real(wp), intent(in)  :: air_mass_capacity   ! [kg/m2]   the canopy air's mass
      real(wp), intent(in)  :: shv                 ! [kg/kg]   its specific humidity
      real(wp), intent(in)  :: qsat                ! [kg/kg]   saturation at its temperature
      real(wp), intent(in)  :: temp                ! [K]       its temperature
      real(wp), intent(out) :: cond                ! [kg/m2/s] condensation rate
      real(wp), intent(out) :: cond_enth           ! [W/m2]    the liquid enthalpy the condensate carries away
      cond      = (air_mass_capacity / TAU_COND) * max(0.0_wp, shv - qsat)
      cond_enth = cond * internal_energy_liquid(temp)
   end subroutine cas_condensation

end module ed_canopy_air_space
