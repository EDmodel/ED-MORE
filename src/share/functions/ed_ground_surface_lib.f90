! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_ground_surface_lib -- the laws of the ground's surface layers and of the porous media
! below and among them (docs/science/ground_surface_layers.md,
! docs/science/soil_biophysics.md):
!   * a porous medium's phase (soil, litter, moss): its temperature and liquid share from its
!     energy on its freezing curve, and back;
!   * a porous layer's depth, the water it holds saturated, the matter that makes a layer, and
!     its thermal resistance as it dries from its top;
!   * the water layer (snow, ice and ponded water): the liquid its ice holds and the free water
!     beyond, its depth over its ice matrix, its matrix's depth as snow falls, water freezes and
!     ice leaves, its cover of the ground, its albedo, and its thermal resistance.
! Laws only: every routine takes its parameters as numbers, so the configuration, the state,
! the fast loop and the slow operators can all call them.
!---~---
module ed_ground_surface_lib
   use ed_kinds,         only : wp, ik
   use ed_constants,     only : rho_h2o, rho_ice, cp_liq, cp_ice, tsupercool_liq, t_3ple, tiny_num, k_water,     &
                               latent_heat_fusion
   use ed_therm_lib,     only : internal_energy_liquid
   use ed_soil_hydr_lib, only : soil_layer_curve_t, soil_liquid_capacity
   implicit none
   private

   !--- A porous medium's phase.
   public :: porous_layer_temp, porous_energy_state, porous_layer_energy, porous_energy_at_temp
   public :: film_on_porous_layer
   !--- A porous layer of the ground's surface.
   public :: porous_layer_depth, porous_layer_capacity, porous_layer_fill, porous_layer_resistance
   public :: porous_layer_matter
   !--- The water layer.
   public :: water_layer_free_water, water_layer_depth, matrix_depth_update
   public :: water_layer_bulk_density
   public :: water_layer_cover, water_layer_albedo, water_layer_optics_shares, water_layer_resistance
   public :: snow_matrix_conductivity

   !---~---
   !   The inversion of a porous medium's energy on its freezing curve: the most Newton steps,
   ! and the temperature step [K] below which it has converged.
   !---~---
   integer(ik), parameter :: PHASE_TEMP_ITER_MAX = 60_ik
   real(wp),    parameter :: PHASE_TEMP_TOL      = 1.0e-10_wp
   !--- [kg/m3] The density Niu & Yang (2007)'s cover scale is referenced to.
   real(wp),    parameter :: COVER_DENSITY_REF = 100.0_wp
   !---~---
   !   Calonne et al. (2011)'s regression of snow's effective thermal conductivity [W/m/K] on its
   ! density rho [kg/m3], k = A2 rho^2 + A1 rho + A0, fitted at 271 K to conductivities computed on
   ! 3-D images of seasonal snow from about 100 to 550 kg/m3.
   !---~---
   real(wp),    parameter :: CALONNE_A2 = 2.5e-6_wp     ! [W m5/kg2/K]
   real(wp),    parameter :: CALONNE_A1 = -1.23e-4_wp   ! [W m2/kg/K]
   real(wp),    parameter :: CALONNE_A0 = 0.024_wp      ! [W/m/K]

contains

   !--- A porous medium's phase.

   !---~---
   !   One porous layer's temperature and liquid fraction from its internal energy and water
   ! content, on its freezing curve (soil_liquid_capacity; docs/science/soil_biophysics.md,
   ! "Frozen soil", eq. (15)): all liquid at and above t_3ple, or below it while the curve
   ! holds the layer's water; below that the liquid the curve holds at the temperature, the rest
   ! ice. Water a layer stores beyond saturation freezes at t_3ple itself, where the curve holds
   ! theta_sat, so an energy between the curve's there and the all-liquid one sits at t_3ple.
   ! The energy rises with the temperature, and a partly frozen layer is warmer than its
   ! all-liquid temperature and colder than its all-ice one and t_3ple, so Newton's iterates are
   ! kept inside that bracket.
   !---~---
   elemental pure subroutine porous_layer_temp(energy, theta, dry_heat_capacity, curve, temp, fliq)
      real(wp),                 intent(in)  :: energy              ! [J/m3]   internal energy
      real(wp),                 intent(in)  :: theta               ! [m3/m3]  water content
      real(wp),                 intent(in)  :: dry_heat_capacity   ! [J/m3/K] heat capacity of the dry medium
      type(soil_layer_curve_t), intent(in)  :: curve               ! the layer's retention curve
      real(wp),                 intent(out) :: temp                ! [K]      temperature
      real(wp),                 intent(out) :: fliq                ! [-]      liquid fraction of the water
      real(wp) :: de
      call porous_energy_state(energy, theta, dry_heat_capacity, curve, temp, fliq, de)
   end subroutine porous_layer_temp

   !---~---
   !   porous_layer_temp, with the apparent heat capacity de [J/m3/K] at the temperature found
   ! (porous_energy_at_temp's), which a Newton on the energy reads.
   !---~---
   elemental pure subroutine porous_energy_state(energy, theta, dry_heat_capacity, curve, temp, fliq, de)
      real(wp),                 intent(in)  :: energy, theta, dry_heat_capacity   ! as porous_layer_temp's
      type(soil_layer_curve_t), intent(in)  :: curve
      real(wp),                 intent(out) :: temp, fliq, de      ! [K], [-], [J/m3/K]
      real(wp)    :: wmass, theta_liq, dtheta, lo, hi, e, step
      integer(ik) :: it
      wmass = theta * rho_h2o
      temp  = (energy + wmass * cp_liq * tsupercool_liq) / (dry_heat_capacity + wmass * cp_liq)
      fliq  = 1.0_wp
      de    = dry_heat_capacity + wmass * cp_liq
      if (wmass <= 0.0_wp .or. temp >= t_3ple) return        ! dry, or liquid at and above t_3ple
      call soil_liquid_capacity(curve, temp, theta_liq, dtheta)
      if (theta_liq >= theta) return                         ! the curve holds all of it: all liquid
      !--- At t_3ple, between the curve's energy there and the all-liquid one.
      call porous_energy_at_temp(t_3ple, theta, dry_heat_capacity, curve, e, de, fliq)
      if (energy >= e) then
         temp = t_3ple
         fliq = (energy - (dry_heat_capacity + wmass * cp_ice) * t_3ple)                                   &
                / (wmass * (internal_energy_liquid(t_3ple) - cp_ice * t_3ple))
         return
      end if
      lo = temp
      hi = min(t_3ple, energy / (dry_heat_capacity + wmass * cp_ice))
      temp = hi
      do it = 1_ik, PHASE_TEMP_ITER_MAX
         call porous_energy_at_temp(temp, theta, dry_heat_capacity, curve, e, de, fliq)
         if (e > energy) then
            hi = temp
         else
            lo = temp
         end if
         step = (energy - e) / de
         if (temp + step <= lo .or. temp + step >= hi) step = 0.5_wp * (lo + hi) - temp   ! bisect
         temp = temp + step
         if (abs(step) < PHASE_TEMP_TOL) exit
      end do
      call porous_energy_at_temp(temp, theta, dry_heat_capacity, curve, e, de, fliq)
   end subroutine porous_energy_state

   !---~---
   !   The inverse: one porous layer's internal energy [J/m3] at a temperature, with the liquid
   ! its freezing curve holds there.
   !---~---
   elemental pure function porous_layer_energy(temp, theta, dry_heat_capacity, curve) result(energy)
      real(wp),                 intent(in) :: temp, theta, dry_heat_capacity   ! [K], [m3/m3], [J/m3/K]
      type(soil_layer_curve_t), intent(in) :: curve
      real(wp) :: energy, de, fliq
      call porous_energy_at_temp(temp, theta, dry_heat_capacity, curve, energy, de, fliq)
   end function porous_layer_energy

   !---~---
   !   A porous layer's internal energy [J/m3] at a temperature, its slope, the apparent heat
   ! capacity [J/m3/K] (the sensible heat and the latent heat of the liquid the curve frees per
   ! kelvin), and its liquid fraction: all liquid above t_3ple, and at t_3ple the bottom of the
   ! plateau on which water stored beyond saturation freezes.
   !---~---
   elemental pure subroutine porous_energy_at_temp(temp, theta, dry_heat_capacity, curve, e, de, fliq)
      real(wp),                 intent(in)  :: temp, theta, dry_heat_capacity   ! [K], [m3/m3], [J/m3/K]
      type(soil_layer_curve_t), intent(in)  :: curve
      real(wp),                 intent(out) :: e, de, fliq                      ! [J/m3], [J/m3/K], [-]
      real(wp) :: wmass, theta_liq, dtheta, dfliq
      wmass = theta * rho_h2o
      fliq  = 1.0_wp
      dfliq = 0.0_wp
      if (wmass > 0.0_wp .and. temp <= t_3ple) then
         call soil_liquid_capacity(curve, temp, theta_liq, dtheta)
         if (theta_liq < theta) then
            fliq  = theta_liq / theta
            dfliq = dtheta / theta
         end if
      end if
      e  = dry_heat_capacity * temp + wmass * (fliq * internal_energy_liquid(temp)                        &
                                             + (1.0_wp - fliq) * cp_ice * temp)
      de = dry_heat_capacity + wmass * (fliq * cp_liq + (1.0_wp - fliq) * cp_ice                          &
                                       + dfliq * (internal_energy_liquid(temp) - cp_ice * temp))
   end subroutine porous_energy_at_temp

   !---~---
   !   A film of water w_film [kg/m2] on a porous layer, the two at one temperature (the water layer's
   ! film rule on a porous layer): their total energy e_total [J/m2] gives the temperature, the film's
   ! liquid fraction and its energy, the film's water on the pure-water law and the layer's on its
   ! freezing curve (porous_layer_temp, with the layer's dry heat capacity dry_cap [J/m2/K], its water
   ! w_layer [kg/m2], its depth [m] and its curve). Below t_3ple the film is ice, a heat capacity
   ! beside the layer's matter; above it both are liquid; between, the two sit at t_3ple and the film
   ! melts first, the layer's water stored beyond saturation after it.
   !---~---
   pure subroutine film_on_porous_layer(e_total, w_film, w_layer, dry_cap, depth, curve, temp, fliq_film, e_film)
      real(wp),                 intent(in)  :: e_total, w_film, w_layer, dry_cap, depth   ! [J/m2], [kg/m2], [kg/m2], [J/m2/K], [m]
      type(soil_layer_curve_t), intent(in)  :: curve
      real(wp),                 intent(out) :: temp, fliq_film, e_film   ! [K], [-], [J/m2]
      real(wp) :: theta, e_ice3, e_liq3, e_bottom3, fl, de
      theta = w_layer / (rho_h2o * depth)
      !--- At t_3ple: the film all ice with the layer at the bottom of its plateau, and both all liquid.
      call porous_energy_at_temp(t_3ple, theta, dry_cap / depth, curve, e_bottom3, de, fl)
      e_bottom3 = e_bottom3 * depth
      e_ice3 = e_bottom3 + w_film * cp_ice * t_3ple
      e_liq3 = dry_cap * t_3ple + (w_layer + w_film) * internal_energy_liquid(t_3ple)
      if (e_total < e_ice3) then
         call porous_layer_temp(e_total / depth, theta, (dry_cap + w_film * cp_ice) / depth, curve, temp, fl)
         fliq_film = 0.0_wp
         e_film    = w_film * cp_ice * temp
      else if (e_total > e_liq3) then
         temp      = (e_total + (w_layer + w_film) * cp_liq * tsupercool_liq) / (dry_cap + (w_layer + w_film) * cp_liq)
         fliq_film = 1.0_wp
         e_film    = w_film * internal_energy_liquid(temp)
      else
         temp   = t_3ple
         e_film = min(w_film * internal_energy_liquid(t_3ple), max(w_film * cp_ice * t_3ple, e_total - e_bottom3))
         fliq_film = 0.0_wp
         if (w_film > 0.0_wp) fliq_film = (e_film - w_film * cp_ice * t_3ple) / (w_film * latent_heat_fusion)
      end if
   end subroutine film_on_porous_layer

   !--- A porous layer of the ground's surface (docs/science/ground_surface_layers.md §1).

   !--- [m] A porous layer's depth: its dry matter [kg/m2] over its material's bulk density [kg/m3].
   elemental real(wp) function porous_layer_depth(matter, bulk_density) result(depth)
      real(wp), intent(in) :: matter, bulk_density
      depth = matter / bulk_density
   end function porous_layer_depth

   !---~---
   !   [kg/m2] The water a porous layer holds saturated: its material's water at saturation
   ! [m3/m3] over its depth, its dry matter [kg/m2] over its bulk density [kg/m3].
   !---~---
   elemental real(wp) function porous_layer_capacity(matter, bulk_density, theta_sat) result(capacity)
      real(wp), intent(in) :: matter, bulk_density, theta_sat
      capacity = theta_sat * rho_h2o * matter / bulk_density
   end function porous_layer_capacity

   !--- [-] The share of its saturation a porous layer's water [kg/m2] fills, at most 1, of its capacity [kg/m2].
   elemental real(wp) function porous_layer_fill(water, capacity) result(fill)
      real(wp), intent(in) :: water, capacity
      fill = min(1.0_wp, water / capacity)
   end function porous_layer_fill

   !---~---
   !   [m2 K/W] A porous layer's thermal resistance top to bottom. It holds its water in its lower
   ! part and dries from its top: its dried top, its depth times the share of saturation its
   ! water does not fill, conducts at the material's dry conductivity, in series with its wet
   ! part at the saturated one [W/m/K].
   !---~---
   elemental real(wp) function porous_layer_resistance(depth, fill, k_dry, k_wet) result(resistance)
      real(wp), intent(in) :: depth, fill, k_dry, k_wet   ! [m], [-], [W/m/K], [W/m/K]
      resistance = depth * ((1.0_wp - fill) / k_dry + fill / k_wet)
   end function porous_layer_resistance

   !---~---
   !   [kg/m2] The dry matter [kg/m2] a porous layer at its bulk density [kg/m3] holds as a layer:
   ! all of it from min_depth [m], none below. Thinner, a film of a few micrometres holds a
   ! fraction of a joule per kelvin and a trace of water, which no balance can carry; its matter
   ! stays where it came from.
   !---~---
   elemental real(wp) function porous_layer_matter(matter, bulk_density, min_depth) result(registered)
      real(wp), intent(in) :: matter, bulk_density, min_depth
      registered = 0.0_wp
      if (bulk_density > 0.0_wp .and. matter >= min_depth * bulk_density) registered = matter
   end function porous_layer_matter

   !--- The water layer: snow, ice and ponded water (docs/science/ground_surface_layers.md §1).

   !---~---
   !   [kg/m2] The free water of a water layer of water [kg/m2] at a liquid fraction [-]: the
   ! liquid beyond what its ice holds, where the ice holds holding_frac of the layer's mass
   ! (LEAF-3's 1:9 at 0.1). That is the liquid L less h/(1-h) of the ice I; all the liquid when
   ! there is no ice; none while the ice holds it all.
   !---~---
   elemental real(wp) function water_layer_free_water(water, fliq, holding_frac) result(free)
      real(wp), intent(in) :: water, fliq, holding_frac   ! [kg/m2], [-], [-]
      free = max(0.0_wp, water * (fliq - holding_frac) / max(1.0_wp - holding_frac, tiny_num))
      free = min(free, water)
   end function water_layer_free_water

   !---~---
   !   [m] A water layer's depth: its ice matrix's, or the volume of its ice and liquid where
   ! that is more. Liquid fills the matrix's pores first, so rain on snow adds no depth until the
   ! pores are full, and a pond (no matrix) is as deep as its water.
   !---~---
   elemental real(wp) function water_layer_depth(matrix_depth, ice, liquid) result(depth)
      real(wp), intent(in) :: matrix_depth, ice, liquid   ! [m], [kg/m2], [kg/m2]
      depth = max(matrix_depth, ice / rho_ice + liquid / rho_h2o)
   end function water_layer_depth

   !---~---
   !   [m] The ice matrix's depth after a step in which snowfall [kg/m2] fell at its fresh density
   ! [kg/m3] and the layer's ice went from ice_old to ice_new [kg/m2] besides that snowfall. Snow
   ! adds its own volume. Ice that forms fills the matrix's pores first, so the matrix grows only
   ! once they are full; ice lost to melt or vapour leaves at the matrix's density, so the matrix
   ! shrinks in proportion.
   !---~---
   elemental real(wp) function matrix_depth_update(matrix_depth, ice_old, ice_new, snowfall, fresh_density) &
      result(depth)
      real(wp), intent(in) :: matrix_depth, ice_old, ice_new   ! [m], [kg/m2], [kg/m2]
      real(wp), intent(in) :: snowfall, fresh_density         ! [kg/m2], [kg/m3]
      real(wp) :: ice_before
      depth      = matrix_depth + snowfall / fresh_density
      ice_before = ice_old + snowfall
      if (ice_new >= ice_before) then
         depth = max(depth, ice_new / rho_ice)
      else if (ice_before > 0.0_wp) then
         depth = depth * (ice_new / ice_before)
      end if
   end function matrix_depth_update

   !--- [kg/m3] A water layer's bulk density: its water [kg/m2] over its depth [m]; 0 with no depth.
   elemental real(wp) function water_layer_bulk_density(water, depth) result(rho)
      real(wp), intent(in) :: water, depth
      rho = 0.0_wp
      if (depth > 0.0_wp) rho = water / depth
   end function water_layer_bulk_density

   !---~---
   !   [-] The share of the ground a water layer of a depth [m] and bulk density [kg/m3] covers,
   ! tanh(depth / scale) (Niu & Yang 2007): the scale is the ground's roughness [m], the
   ! micro-relief a thin layer has to fill, times scale_factor, and grows with the density by its
   ! exponent, so a dense layer covers less of the ground at the same depth.
   !---~---
   elemental real(wp) function water_layer_cover(depth, bulk_density, scale_factor, density_exponent, ground_rough) &
      result(cover)
      real(wp), intent(in) :: depth, bulk_density, scale_factor, density_exponent, ground_rough
      real(wp) :: scale
      scale = max(scale_factor * ground_rough * (bulk_density / COVER_DENSITY_REF) ** density_exponent, tiny_num)
      cover = tanh(depth / scale)
   end function water_layer_cover

   !---~---
   !   [-] A water layer's albedo in one band: the water standing above its ice matrix (its share
   ! standing) at standing water's, the rest at its matrix's: snow's, fresh to aged by the layer's
   ! liquid fraction, over the matrix's porous share, solid ice's over the rest.
   !---~---
   elemental real(wp) function water_layer_albedo(standing, porous_share, fliq, albedo_fresh, albedo_aged,       &
                                                  albedo_ice, albedo_liquid) result(albedo)
      real(wp), intent(in) :: standing, porous_share, fliq           ! [-] (water_layer_optics_shares), [-]
      real(wp), intent(in) :: albedo_fresh, albedo_aged, albedo_ice, albedo_liquid
      real(wp) :: snow
      snow   = (1.0_wp - fliq) * albedo_fresh + fliq * albedo_aged
      albedo = standing * albedo_liquid + (1.0_wp - standing) * (porous_share * snow + (1.0_wp - porous_share) * albedo_ice)
   end function water_layer_albedo

   !---~---
   !   The shares a water layer's optics weigh (water_layer_albedo), from its ice matrix's depth
   ! [m], its ice and its liquid [kg/m2]: standing, the share of its depth that its ice and liquid
   ! fill beyond the matrix, water standing over it (none in snow whose pores are not full, all of
   ! a pond); and porous_share, how far the matrix's density lies from solid ice's towards the
   ! fresh density [kg/m3] (one for fresh snow, none for ice).
   !---~---
   elemental subroutine water_layer_optics_shares(matrix_depth, ice, liquid, fresh_density, standing, porous_share)
      real(wp), intent(in)  :: matrix_depth, ice, liquid, fresh_density   ! [m], [kg/m2], [kg/m2], [kg/m3]
      real(wp), intent(out) :: standing, porous_share                     ! [-]
      real(wp) :: depth, rho_matrix
      depth = water_layer_depth(matrix_depth, ice, liquid)
      standing     = 0.0_wp
      porous_share = 0.0_wp
      if (.not. (depth > 0.0_wp)) return
      standing = min(1.0_wp, max(0.0_wp, ice / rho_ice + liquid / rho_h2o - matrix_depth) / depth)
      if (.not. (matrix_depth > 0.0_wp)) return
      rho_matrix   = ice / matrix_depth
      porous_share = min(1.0_wp, max(0.0_wp, (rho_ice - rho_matrix) / (rho_ice - fresh_density)))
   end subroutine water_layer_optics_shares

   !---~---
   !   [W/m/K] The thermal conductivity of a water layer's ice matrix at its density [kg/m3]:
   ! Calonne et al. (2011)'s density law as a shape, scaled so the matrix conducts k_fresh at the
   ! fresh density. The law is fitted from about 100 to 550 kg/m3; a denser matrix, up to solid
   ! ice, takes the law's shape beyond that range, still scaled by k_fresh, so from 0.15 W/m/K in
   ! fresh snow it reaches about 2 at solid ice's.
   !---~---
   elemental real(wp) function snow_matrix_conductivity(rho_matrix, k_fresh, fresh_density) result(k)
      real(wp), intent(in) :: rho_matrix, k_fresh, fresh_density   ! [kg/m3], [W/m/K], [kg/m3]
      k = k_fresh * calonne(rho_matrix) / calonne(fresh_density)
   contains
      elemental real(wp) function calonne(rho) result(kc)
         real(wp), intent(in) :: rho
         kc = (CALONNE_A2 * rho + CALONNE_A1) * rho + CALONNE_A0
      end function calonne
   end function snow_matrix_conductivity

   !---~---
   !   [m2 K/W] A water layer's thermal resistance top to bottom: its ice matrix's depth at the
   ! matrix's conductivity [W/m/K], in series with the water standing above the matrix at
   ! liquid water's.
   !---~---
   elemental real(wp) function water_layer_resistance(depth, matrix_depth, k_matrix) result(resistance)
      real(wp), intent(in) :: depth, matrix_depth, k_matrix   ! [m], [m], [W/m/K]
      resistance = max(matrix_depth, tiny_num) / max(k_matrix, tiny_num) + max(0.0_wp, depth - matrix_depth) / k_water
   end function water_layer_resistance

end module ed_ground_surface_lib
