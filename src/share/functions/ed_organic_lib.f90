! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_organic_lib -- the organic layer on the soil, leaf litter over moss
! (docs/science/surface_organic_layer.md): the [surface_layer] block's materials, the layer's
! per-patch record, and its laws -- its sums as one thermal node, its temperature, its energy
! shared among the slots at that temperature, each slot's depth and its water retention as a
! node of the soil-plant water network, the exchange it takes, and the ground top's temperature
! it sets. A law library, as ed_allometry_lib is, because the configuration, the state, the fast
! loop and the slow operators all read it.
!---~---
module ed_organic_lib
   use ed_kinds,     only : wp, ik
   use ed_constants, only : n_organic_slot, cp_organic, rho_h2o, t_3ple
   use ed_therm_lib, only : internal_energy_to_temp, temp_to_internal_energy
   use ed_soil_hydr_lib, only : soil_layer_curve_t, build_soil_layer_curve, SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   use ed_numerics,      only : within_roundoff
   implicit none
   private

   public :: organic_material_t, organic_params_t, organic_state_t, SLOT_LITTER, SLOT_MOSS
   public :: organic_present, organic_slot_material, organic_layer_sums, organic_diagnose, organic_split_energy
   public :: organic_take_exchange, organic_slot_curve, organic_slot_capacity, organic_spill
   public :: ORGANIC_MIN_DEPTH, organic_registered_mass, organic_slot_depth, ground_top_temp

   !---~---
   !   One organic material (docs/science/surface_organic_layer.md): its dry mass per volume, its
   ! thermal conductivity dry and saturated, and its water retention and conductivity as a soil
   ! layer's (ed_soil_hydr_lib): the curve's form, the water at saturation and at residual, its
   ! two parameters and the saturated conductivity. Its values are its slot's, in organic_params_t.
   !---~---
   type :: organic_material_t
      real(wp)    :: bulk_density   ! [kg/m3] dry mass per volume
      real(wp)    :: k_dry          ! [W/m/K] thermal conductivity, dry
      real(wp)    :: k_wet          ! [W/m/K] thermal conductivity, saturated
      integer(ik) :: retention      ! SOIL_RETENTION_VG or SOIL_RETENTION_CAMPBELL
      real(wp)    :: theta_sat      ! [m3/m3] water at saturation
      real(wp)    :: theta_res      ! [m3/m3] residual water
      real(wp)    :: par_a          ! vG alpha [1/m], or Campbell's air-entry head [m, < 0]
      real(wp)    :: par_n          ! vG n [-], or Campbell's b [-]
      real(wp)    :: ksat           ! [m/s] saturated conductivity
   end type organic_material_t

   !---~---
   !   The [surface_layer] block: the ground's organic layer, leaf litter over moss. The litter
   ! is the above-ground metabolic litter pool's dry matter (Napoly et al. 2017's forest-floor values),
   ! retaining and conducting water as Letts et al. (2000)'s surface moss, the curve CLASS gives its
   ! litter and moss (Wu et al. 2016); the moss is a site's [soil].moss_depth, at the MAT tundra's
   ! top organic 5 cm (ARC LTER knb-lter-arc.20129.1), Chadburn et al. (2015)'s porosity and
   ! conductivities, Zhang et al. (2010)'s retention of surface moss and O'Connor et al. (2020)'s
   ! saturated conductivity of the Arctic acrotelm.
   !---~---
   type :: organic_params_t
      logical :: litter_on = .true.                                                ! the litter slot
      type(organic_material_t) :: litter = organic_material_t(bulk_density=45.0_wp, k_dry=0.10_wp,               &
                                              k_wet=0.1294_wp, retention=SOIL_RETENTION_CAMPBELL,                 &
                                              theta_sat=0.98_wp, theta_res=0.01_wp, par_a=-0.0103_wp,             &
                                              par_n=2.3_wp, ksat=1.83e-3_wp)
      type(organic_material_t) :: moss   = organic_material_t(bulk_density=43.0_wp, k_dry=0.06_wp,               &
                                              k_wet=0.50_wp, retention=SOIL_RETENTION_VG,                         &
                                              theta_sat=0.85_wp, theta_res=0.18_wp, par_a=50.0_wp,                &
                                              par_n=1.45_wp, ksat=1.3e-3_wp)
   end type organic_params_t

   !---~---
   !   The ground's organic layer above the soil, its slots top down (leaf litter over moss),
   ! each with its dry matter, its water (liquid and ice) and its internal energy; its depth is
   ! its matter over its material's bulk density. The slots share one temperature, read off their
   ! summed energy, water and dry heat capacity; a slot with no matter is absent
   ! (docs/science/surface_organic_layer.md).
   !---~---
   integer(ik), parameter :: SLOT_LITTER = 1_ik, SLOT_MOSS = 2_ik
   !---~---
   !   [m] The depth from which a slot's matter is a layer of the ground. Thinner, a film of a few
   ! micrometres, it holds a fraction of a joule per kelvin and a trace of water, which no balance
   ! can carry; its matter stays in its pool, and the ground's surface is what lies below it.
   !---~---
   real(wp), parameter :: ORGANIC_MIN_DEPTH = 1.0e-4_wp

   type :: organic_state_t
      real(wp) :: mass(n_organic_slot)    = 0.0_wp   ! [kg/m2] dry matter
      real(wp) :: water(n_organic_slot)   = 0.0_wp   ! [kg/m2] water, liquid and ice (prognostic)
      real(wp) :: energy(n_organic_slot)  = 0.0_wp   ! [J/m2]  internal energy (prognostic)
      real(wp) :: temp = t_3ple                      ! [K]     the layer's, diagnosed
      real(wp) :: fliq = 1.0_wp                      ! [-]     its water's liquid fraction, diagnosed
   end type organic_state_t

contains

   !--- The dry matter [kg/m2] a slot of mass [kg/m2] at bulk density [kg/m3] holds as a layer: none below ORGANIC_MIN_DEPTH.
   elemental real(wp) function organic_registered_mass(mass, density) result(registered)
      real(wp), intent(in) :: mass, density
      registered = 0.0_wp
      if (density > 0.0_wp .and. mass >= ORGANIC_MIN_DEPTH * density) registered = mass
   end function organic_registered_mass

   !--- Whether a patch has an organic layer: dry matter in any slot.
   pure logical function organic_present(org)
      type(organic_state_t), intent(in) :: org
      organic_present = any(org%mass > 0.0_wp)
   end function organic_present

   !--- [m] Slot k's depth: its matter over its material's bulk density.
   pure real(wp) function organic_slot_depth(params, org, k) result(depth)
      type(organic_params_t), intent(in) :: params
      type(organic_state_t),  intent(in) :: org
      integer(ik),            intent(in) :: k
      type(organic_material_t) :: m
      m     = organic_slot_material(params, k)
      depth = org%mass(k) / m%bulk_density
   end function organic_slot_depth

   !---~---
   !   [K] The temperature of the ground's top, which the canopy air and the radiation see: the
   ! organic layer's where there is one, else the top soil node's.
   !---~---
   pure real(wp) function ground_top_temp(org, soil_temp_top) result(t)
      type(organic_state_t), intent(in) :: org
      real(wp),              intent(in) :: soil_temp_top   ! [K]
      t = soil_temp_top
      if (organic_present(org)) t = org%temp
   end function ground_top_temp

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
   ! capacity, water, energy and the water it holds saturated summed over the slots, and its
   ! thermal resistance. Each slot holds its water in its lower part and dries from its top: its
   ! dried top, its depth times the share of saturation it does not hold, is what the layer's own
   ! vapour crosses, and the heat crosses it too, at the slot's dry conductivity, in series with
   ! its wet part at its saturated conductivity. The slots are in series.
   !---~---
   pure subroutine organic_layer_sums(org, params, dry_cap, water, energy, w_max, resistance, dry_depth)
      type(organic_state_t),  intent(in)  :: org
      type(organic_params_t), intent(in)  :: params
      real(wp),               intent(out) :: dry_cap      ! [J/m2/K]
      real(wp),               intent(out) :: water        ! [kg/m2]
      real(wp),               intent(out) :: energy       ! [J/m2]
      real(wp),               intent(out) :: w_max        ! [kg/m2] the water the slots hold saturated
      real(wp),               intent(out) :: resistance   ! [m2 K/W]
      real(wp),               intent(out) :: dry_depth    ! [m] the dried top of the slots, summed
      type(organic_material_t) :: m
      real(wp)    :: depth, cap, fill
      integer(ik) :: k
      dry_cap    = 0.0_wp
      water      = 0.0_wp
      energy     = 0.0_wp
      w_max      = 0.0_wp
      resistance = 0.0_wp
      dry_depth  = 0.0_wp
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp) cycle
         m      = organic_slot_material(params, k)
         depth  = org%mass(k) / m%bulk_density
         cap    = organic_slot_capacity(params, org, k)
         fill   = min(1.0_wp, org%water(k) / cap)
         dry_cap    = dry_cap + cp_organic * org%mass(k)
         water      = water + org%water(k)
         energy     = energy + org%energy(k)
         w_max      = w_max + cap
         resistance = resistance + depth * ((1.0_wp - fill) / m%k_dry + fill / m%k_wet)
         dry_depth  = dry_depth + depth * (1.0_wp - fill)
      end do
   end subroutine organic_layer_sums

   !--- [kg/m2] The water slot k holds saturated: its material's theta_sat over its depth, its matter over its bulk density.
   pure real(wp) function organic_slot_capacity(params, org, k) result(cap)
      type(organic_params_t), intent(in) :: params
      type(organic_state_t),  intent(in) :: org
      integer(ik),            intent(in) :: k
      type(organic_material_t) :: m
      m   = organic_slot_material(params, k)
      cap = m%theta_sat * rho_h2o * org%mass(k) / m%bulk_density
   end function organic_slot_capacity

   !---~---
   !   Slot k's water retention and conductivity, its material's as a soil layer's
   ! (build_soil_layer_curve), at its depth.
   !---~---
   pure function organic_slot_curve(params, org, k) result(c)
      type(organic_params_t), intent(in) :: params
      type(organic_state_t),  intent(in) :: org
      integer(ik),            intent(in) :: k
      type(soil_layer_curve_t)           :: c
      type(organic_material_t) :: m
      m = organic_slot_material(params, k)
      call build_soil_layer_curve(m%retention, m%theta_sat, m%theta_res, m%par_a, m%par_n, m%ksat,                &
                                  org%mass(k) / m%bulk_density, c)
   end function organic_slot_curve

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
   !   The layer takes an exchange's amounts after the step's march: the energy in, and the vapour
   ! each slot gave the canopy air as the water network took it (dew, negative, joined the top
   ! slot), so the water that left carried its share of that energy.
   !---~---
   pure subroutine organic_take_exchange(org, energy_amount, vapour_amount)
      type(organic_state_t), intent(inout) :: org
      real(wp),              intent(in)    :: energy_amount                    ! [J/m2] in
      real(wp),              intent(in)    :: vapour_amount(n_organic_slot)    ! [kg/m2] out of each slot
      if (.not. organic_present(org)) return
      if (.not. all(org%water - vapour_amount >= 0.0_wp                                                     &
                    .or. within_roundoff(org%water - vapour_amount, org%water + abs(vapour_amount))))         &
         error stop 'organic_take_exchange: a slot gave the canopy air more vapour than it held'
      org%water = max(0.0_wp, org%water - vapour_amount)   ! clamp-ok: the network's own amounts, to round-off (asserted)
      call organic_split_energy(org, sum(org%energy) + energy_amount)
   end subroutine organic_take_exchange

   !---~---
   !   The water a slot holds above saturation spills: dew or frost on a saturated or frozen slot,
   ! which its pores cannot take and its ice keeps from passing on. Its water and energy, at the
   ! layer's temperature and liquid fraction, are returned for the pond, which takes them at a
   ! step's start.
   !---~---
   pure subroutine organic_spill(org, params, w_out, e_out)
      type(organic_state_t),  intent(inout) :: org
      type(organic_params_t), intent(in)    :: params
      real(wp),               intent(out)   :: w_out     ! [kg/m2]
      real(wp),               intent(out)   :: e_out     ! [J/m2]
      real(wp)    :: excess, total
      integer(ik) :: k
      w_out = 0.0_wp
      e_out = 0.0_wp
      if (.not. organic_present(org)) return
      total = sum(org%energy)
      call organic_diagnose(org)
      do k = 1_ik, n_organic_slot
         if (org%mass(k) <= 0.0_wp) cycle
         excess = max(0.0_wp, org%water(k) - organic_slot_capacity(params, org, k))
         if (excess <= 0.0_wp) cycle
         org%water(k) = org%water(k) - excess
         w_out = w_out + excess
         e_out = e_out + temp_to_internal_energy(0.0_wp, excess, org%temp, org%fliq)
      end do
      if (w_out > 0.0_wp) call organic_split_energy(org, total - e_out)
   end subroutine organic_spill

end module ed_organic_lib
