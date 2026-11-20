! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_organic_lib -- the organic layer on the soil, leaf litter over moss
! (docs/science/surface_organic_layer.md): the [surface_layer] block's materials, the layer's
! per-patch record, and its laws -- its sums as one thermal node, its temperature, its energy
! shared among the slots at that temperature, the water it takes and the exchange it takes. A
! law library, as ed_allometry_lib is, because the configuration, the state, the fast loop
! and the slow operators all read it.
!---~---
module ed_organic_lib
   use ed_kinds,     only : wp, ik
   use ed_constants, only : n_organic_slot, cp_organic, rho_h2o, t_3ple
   use ed_therm_lib, only : internal_energy_to_temp, temp_to_internal_energy
   implicit none
   private

   public :: organic_material_t, organic_params_t, organic_state_t, SLOT_LITTER, SLOT_MOSS
   public :: organic_present, organic_slot_material, organic_layer_sums, organic_diagnose, organic_split_energy
   public :: organic_take_water, organic_take_exchange, organic_spill

   !---~---
   !   One organic material (docs/science/surface_organic_layer.md): its dry mass per volume,
   ! the water it holds per volume, and its conductivity dry and at that capacity, linear in
   ! between.
   !---~---
   type :: organic_material_t
      real(wp) :: bulk_density   = 45.0_wp      ! [kg/m3] dry mass per volume
      real(wp) :: water_capacity = 0.12_wp      ! [m3/m3] water held per volume
      real(wp) :: k_dry          = 0.10_wp      ! [W/m/K] thermal conductivity, dry
      real(wp) :: k_wet          = 0.1036_wp    ! [W/m/K] thermal conductivity at its water capacity
   end type organic_material_t

   !---~---
   !   The [surface_layer] block: the ground's organic layer, leaf litter over moss. The litter
   ! is the above-ground metabolic litter pool's dry matter (Napoly et al. 2017's forest-floor values);
   ! the moss is a site's [soil].moss_depth, at the MAT tundra's top organic 5 cm (ARC LTER
   ! knb-lter-arc.20129.1) and Chadburn et al. (2015)'s porosity and conductivities.
   !---~---
   type :: organic_params_t
      logical :: litter_on = .true.                                                ! the litter slot
      type(organic_material_t) :: litter = organic_material_t(45.0_wp, 0.12_wp, 0.10_wp, 0.1036_wp)
      type(organic_material_t) :: moss   = organic_material_t(43.0_wp, 0.90_wp, 0.06_wp, 0.50_wp)
   end type organic_params_t

   !---~---
   !   The ground's organic layer above the soil, its slots top down (leaf litter over moss),
   ! each with its dry matter, its water (liquid and ice), its internal energy and its bulk
   ! density, which makes the matter a depth. The slots share one temperature, read off their
   ! summed energy, water and dry heat capacity; a slot with no matter is absent
   ! (docs/science/surface_organic_layer.md).
   !---~---
   integer(ik), parameter :: SLOT_LITTER = 1_ik, SLOT_MOSS = 2_ik
   type :: organic_state_t
      real(wp) :: mass(n_organic_slot)    = 0.0_wp   ! [kg/m2] dry matter
      real(wp) :: water(n_organic_slot)   = 0.0_wp   ! [kg/m2] water, liquid and ice (prognostic)
      real(wp) :: energy(n_organic_slot)  = 0.0_wp   ! [J/m2]  internal energy (prognostic)
      real(wp) :: density(n_organic_slot) = 0.0_wp   ! [kg/m3] bulk density, dry matter over depth
      real(wp) :: temp = t_3ple                      ! [K]     the layer's, diagnosed
      real(wp) :: fliq = 1.0_wp                      ! [-]     its water's liquid fraction, diagnosed
   end type organic_state_t

contains

   !--- Whether a patch has an organic layer: dry matter in any slot.
   pure logical function organic_present(org)
      type(organic_state_t), intent(in) :: org
      organic_present = any(org%mass > 0.0_wp)
   end function organic_present

   !--- Slot k's material: the litter's on top, the moss's below it.
   pure function organic_slot_material(params, k) result(m)
      type(organic_params_t), intent(in) :: params
      integer(ik),            intent(in) :: k
      type(organic_material_t)           :: m
      if (k == SLOT_LITTER) then
         m = params%litter
      else
         m = params%moss
      end if
   end function organic_slot_material

   !---~---
   !   The organic layer as one thermal node (docs/science/surface_organic_layer.md): its dry heat
   ! capacity, water, energy and water capacity summed over the slots, and its thermal
   ! resistance, each slot's depth over its conductivity in series, the conductivity linear in
   ! the slot's fill between its dry and its wet value.
   !---~---
   pure subroutine organic_layer_sums(org, params, dry_cap, water, energy, w_max, resistance)
      type(organic_state_t),  intent(in)  :: org
      type(organic_params_t), intent(in)  :: params
      real(wp),               intent(out) :: dry_cap      ! [J/m2/K]
      real(wp),               intent(out) :: water        ! [kg/m2]
      real(wp),               intent(out) :: energy       ! [J/m2]
      real(wp),               intent(out) :: w_max        ! [kg/m2] the water the slots hold at capacity
      real(wp),               intent(out) :: resistance   ! [m2 K/W]
      type(organic_material_t) :: m
      real(wp)    :: depth, cap, fill
      integer(ik) :: k
      dry_cap    = 0.0_wp
      water      = 0.0_wp
      energy     = 0.0_wp
      w_max      = 0.0_wp
      resistance = 0.0_wp
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp) cycle
         m      = organic_slot_material(params, k)
         depth  = org%mass(k) / org%density(k)
         cap    = m%water_capacity * rho_h2o * depth
         fill   = min(1.0_wp, org%water(k) / cap)
         dry_cap    = dry_cap + cp_organic * org%mass(k)
         water      = water + org%water(k)
         energy     = energy + org%energy(k)
         w_max      = w_max + cap
         resistance = resistance + depth / (m%k_dry + (m%k_wet - m%k_dry) * fill)
      end do
   end subroutine organic_layer_sums

   !--- The layer's temperature and liquid fraction, from its summed energy, water and dry heat capacity.
   pure subroutine organic_diagnose(org)
      type(organic_state_t), intent(inout) :: org
      if (.not. organic_present(org)) then
         org%temp = t_3ple
         org%fliq = 1.0_wp
         return
      end if
      call internal_energy_to_temp(sum(org%energy), sum(org%water), cp_organic * sum(org%mass), org%temp, org%fliq)
   end subroutine organic_diagnose

   !---~---
   !   The slots share one temperature: a total energy split among them, each its matter's and
   ! its water's at the layer's temperature and liquid fraction, the last present slot taking
   ! the remainder so the sum is the total exactly.
   !---~---
   pure subroutine organic_split_energy(org, total)
      type(organic_state_t), intent(inout) :: org
      real(wp),              intent(in)    :: total     ! [J/m2]
      real(wp)    :: rest
      integer(ik) :: k, last
      org%energy = 0.0_wp
      if (.not. organic_present(org)) return
      call internal_energy_to_temp(total, sum(org%water), cp_organic * sum(org%mass), org%temp, org%fliq)
      last = 0_ik
      do k = 1_ik, n_organic_slot
         if (org%mass(k) > 0.0_wp) last = k
      end do
      rest = total
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp .or. k == last) cycle
         org%energy(k) = temp_to_internal_energy(cp_organic * org%mass(k), org%water(k), org%temp, org%fliq)
         rest = rest - org%energy(k)
      end do
      org%energy(last) = rest
   end subroutine organic_split_energy

   !---~---
   !   Liquid water reaching the layer fills its slots top down, each to its capacity, and brings
   ! its energy; taken returns what the slots took, the rest passing on to the pond.
   !---~---
   pure subroutine organic_take_water(org, params, mass, energy, taken)
      type(organic_state_t),  intent(inout) :: org
      type(organic_params_t), intent(in)    :: params
      real(wp),               intent(in)    :: mass      ! [kg/m2] the water reaching the layer
      real(wp),               intent(in)    :: energy    ! [J/m2]  its energy
      real(wp),               intent(out)   :: taken     ! [kg/m2]
      type(organic_material_t) :: m
      real(wp)    :: room, add
      integer(ik) :: k
      taken = 0.0_wp
      if (mass <= 0.0_wp .or. .not. organic_present(org)) return
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp) cycle
         m    = organic_slot_material(params, k)
         room = max(0.0_wp, m%water_capacity * rho_h2o * org%mass(k) / org%density(k) - org%water(k))
         add  = min(room, mass - taken)
         org%water(k) = org%water(k) + add
         taken        = taken + add
      end do
      if (taken > 0.0_wp) call organic_split_energy(org, sum(org%energy) + energy * (taken / mass))
   end subroutine organic_take_water

   !---~---
   !   The layer takes an exchange's amounts: the energy in and the vapour out, the vapour from
   ! the slots top down, as far as each holds water; dew (vapour in) joins the top slot.
   !---~---
   pure subroutine organic_take_exchange(org, energy_amount, vapour_amount)
      type(organic_state_t), intent(inout) :: org
      real(wp),              intent(in)    :: energy_amount   ! [J/m2] in
      real(wp),              intent(in)    :: vapour_amount   ! [kg/m2] out
      real(wp)    :: left, out
      integer(ik) :: k
      if (.not. organic_present(org)) return
      left = vapour_amount
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp) cycle
         if (left < 0.0_wp) then                               ! dew: the top slot takes it
            org%water(k) = org%water(k) - left
            left = 0.0_wp
         end if
         out = min(max(left, 0.0_wp), org%water(k))
         org%water(k) = org%water(k) - out
         left = left - out
      end do
      call organic_split_energy(org, sum(org%energy) + energy_amount)
   end subroutine organic_take_exchange

   !---~---
   !   The water a slot holds beyond its capacity (dew on a full layer) spills: its water and energy,
   ! at the layer's temperature and liquid fraction, returned for the pond.
   !---~---
   pure subroutine organic_spill(org, params, w_out, e_out)
      type(organic_state_t),  intent(inout) :: org
      type(organic_params_t), intent(in)    :: params
      real(wp),               intent(out)   :: w_out     ! [kg/m2]
      real(wp),               intent(out)   :: e_out     ! [J/m2]
      type(organic_material_t) :: m
      real(wp)    :: cap, excess, total
      integer(ik) :: k
      w_out = 0.0_wp
      e_out = 0.0_wp
      if (.not. organic_present(org)) return
      total = sum(org%energy)
      call organic_diagnose(org)
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp) cycle
         m      = organic_slot_material(params, k)
         cap    = m%water_capacity * rho_h2o * org%mass(k) / org%density(k)
         excess = max(0.0_wp, org%water(k) - cap)
         if (excess <= 0.0_wp) cycle
         org%water(k) = cap
         w_out = w_out + excess
         e_out = e_out + temp_to_internal_energy(0.0_wp, excess, org%temp, org%fliq)
      end do
      if (w_out > 0.0_wp) call organic_split_energy(org, total - e_out)
   end subroutine organic_spill

end module ed_organic_lib
