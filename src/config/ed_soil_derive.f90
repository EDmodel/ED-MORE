! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_soil_derive -- a site's soil profile (soil_profile_t) from its [soil]
! (soil_params_t), built once per site and per map cell: the layers, cut from
! layer_bottoms or SOIL_LAYER_TEMPLATE at the depth; each layer's texture, the soil
! maps' or the one texture; each layer's retention curve, Ks and thermal texture; and
! the field capacity. Also the water a layer starts the run with, the top layer's
! potential, and a map cell's soil over a layer.
!---~---
module ed_soil_derive
   use ed_kinds,           only : wp, ik
   use ed_constants,       only : n_soil_layer_max
   use ed_soil_hydr_lib,   only : soil_head_to_theta, soil_matric_potential, soil_texture_t,  &
                                  build_soil_layer_curve
   use ed_soil_params,     only : soil_params_t, soil_texture_layers_t, soil_profile_t
   use ed_maps_params,     only : map_inputs_t
   implicit none
   private

   public :: SOIL_LAYER_TEMPLATE, soil_layer_grid, soil_layers, soil_layer_textures
   public :: build_soil_profile, soil_start_theta, top_soil_psi, map_layer_texture
   public :: build_profile_hydraulics, build_profile_thermal

   !---~---
   !   The default layer bottoms [m, down], cut at each site's depth: SoilGrids 2.0's
   ! depths to 1 m, then thickening to 27 m (docs/science/soil_biophysics.md, "The
   ! soil layers").
   !---~---
   real(wp), parameter :: SOIL_LAYER_TEMPLATE(15) = [0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp,     &
                          1.50_wp, 2.00_wp, 3.00_wp, 4.50_wp, 6.50_wp, 9.00_wp, 12.0_wp, 16.0_wp, 21.0_wp, 27.0_wp]
   !---~---
   !   A remainder below the cut thinner than this share of the list layer it falls in
   ! joins the layer above, so no profile ends in a sliver.
   !---~---
   real(wp), parameter :: LAYER_MERGE_SHARE = 0.5_wp
   real(wp), parameter :: GRID_TOL = 1.0e-6_wp     ! [m] a depth this close to a list bottom ends on it

contains

   !---~---
   !   The profile: its layers, its texture per layer, and the profile the soil kernels
   ! read, built once per site and per map cell.
   !---~---

   !---~---
   !   The layers: the bottoms [m, positive down] of `layer_bottoms` when [soil] sets it,
   ! else SOIL_LAYER_TEMPLATE, cut at `depth` (soil_layer_grid).
   !---~---
   pure subroutine soil_layers(soil, bottoms, n_layer)
      type(soil_params_t), intent(in)  :: soil
      real(wp),            intent(out) :: bottoms(n_soil_layer_max)
      integer(ik),         intent(out) :: n_layer
      if (soil%n_layer_bottoms > 0_ik) then
         call soil_layer_grid(soil%depth, soil%layer_bottoms, soil%n_layer_bottoms, bottoms, n_layer)
      else
         call soil_layer_grid(soil%depth, SOIL_LAYER_TEMPLATE, size(SOIL_LAYER_TEMPLATE, kind=ik),   &
                              bottoms, n_layer)
      end if
   end subroutine soil_layers

   !--- The texture per layer: the maps' per-layer values, or the one texture.
   pure subroutine soil_layer_textures(soil, tex)
      type(soil_params_t),         intent(in)  :: soil
      type(soil_texture_layers_t), intent(out) :: tex
      if (soil%n_layer_values > 0_ik) then
         tex%theta_sat = soil%layer_theta_sat
         tex%theta_res = soil%layer_theta_res
         tex%ksat = soil%layer_ksat
         tex%par_a     = soil%layer_curve_par_a
         tex%par_n   = soil%layer_curve_par_n
         tex%k_solid   = soil%layer_solid_conductivity
         tex%k_dry = soil%layer_dry_conductivity
         tex%dry_cvol  = soil%layer_dry_heat_capacity
      else
         tex%theta_sat = soil%theta_sat
         tex%theta_res = soil%theta_res
         tex%ksat = soil%ksat
         tex%par_a     = soil%curve_par_a
         tex%par_n   = soil%curve_par_n
         tex%k_solid   = soil%solid_conductivity
         tex%k_dry = soil%dry_conductivity
         tex%dry_cvol  = soil%dry_heat_capacity
      end if
   end subroutine soil_layer_textures

   !---~---
   !   The soil profile [soil] describes: its layers (soil_layers), each layer's retention
   ! curve, Ks and thermal texture. The sites, the slow ledger and validate_config all
   ! build the profile here.
   !---~---
   pure subroutine build_soil_profile(soil, profile)
      type(soil_params_t),  intent(in)  :: soil
      type(soil_profile_t), intent(out) :: profile
      type(soil_texture_layers_t) :: tex
      real(wp)    :: bottoms(n_soil_layer_max)
      integer(ik) :: n_layer
      call soil_layers(soil, bottoms, n_layer)
      call soil_layer_textures(soil, tex)
      call build_profile_hydraulics_layers(n_layer, soil%retention, bottoms, tex%theta_sat, tex%theta_res, &
                                           tex%ksat, tex%par_a, tex%par_n, soil%psi_fc, profile)
      call build_profile_thermal_layers(n_layer, tex%k_solid, tex%k_dry, tex%dry_cvol, profile)
   end subroutine build_soil_profile

   !---~---
   !   The water each layer starts the run with [m3/m3]: [init].soil_theta, or each layer's
   ! field capacity for a profile with a texture per layer (the soil maps).
   !---~---
   pure subroutine soil_start_theta(soil, profile, init_soil_theta, theta)
      type(soil_params_t),  intent(in)  :: soil
      type(soil_profile_t), intent(in)  :: profile
      real(wp),             intent(in)  :: init_soil_theta   ! [m3/m3] [init].soil_theta
      real(wp),             intent(out) :: theta(n_soil_layer_max)
      if (soil%n_layer_values > 0_ik) then
         theta = profile%theta_fc
      else
         theta = init_soil_theta
      end if
   end subroutine soil_start_theta

   !---~---
   !   The water potential [MPa] of the profile's top layer at water content theta_top: the
   ! predawn potential of a plant that has relaxed to the soil overnight. It seeds the
   ! stomatal predawn potential and a new cohort's tissue water.
   !---~---
   pure real(wp) function top_soil_psi(profile, theta_top) result(psi)
      type(soil_profile_t), intent(in) :: profile
      real(wp),             intent(in) :: theta_top   ! [m3/m3]
      psi = soil_matric_potential(profile%retention, theta_top, profile%theta_sat(1), profile%theta_res(1),   &
                                  profile%curve_par_a(1), profile%curve_par_n(1))
   end function top_soil_psi

   !---~---
   !   Map cell k's soil over the layer [top, bottom]: the map layers it overlaps, each
   ! weighted by the overlap. The deepest map layer reaches down without end.
   !---~---
   pure function map_layer_texture(inp, k, top, bottom) result(tex)
      type(map_inputs_t), intent(in) :: inp
      integer(ik),        intent(in) :: k
      real(wp),           intent(in) :: top, bottom
      type(soil_texture_t) :: tex
      real(wp)    :: w, wsum, hi
      integer(ik) :: m, nm
      nm = size(inp%layer_top, kind=ik)
      tex%sand = 0.0_wp
      tex%silt = 0.0_wp
      tex%clay = 0.0_wp
      tex%org_carbon = 0.0_wp
      tex%cec  = 0.0_wp
      tex%ph   = 0.0_wp
      tex%bulk_density = 0.0_wp
      tex%coarse_frag = 0.0_wp
      wsum = 0.0_wp
      do m = 1_ik, nm
         hi = bottom
         if (m < nm) hi = min(bottom, inp%layer_bottom(m))
         w = hi - max(top, inp%layer_top(m))
         if (w <= 0.0_wp) cycle
         wsum = wsum + w
         associate (t => inp%texture(m, k))
            tex%sand         = tex%sand         + w * t%sand
            tex%silt         = tex%silt         + w * t%silt
            tex%clay         = tex%clay         + w * t%clay
            tex%org_carbon   = tex%org_carbon   + w * t%org_carbon
            tex%cec          = tex%cec          + w * t%cec
            tex%ph           = tex%ph           + w * t%ph
            tex%bulk_density = tex%bulk_density + w * t%bulk_density
            tex%coarse_frag  = tex%coarse_frag  + w * t%coarse_frag
         end associate
      end do
      tex%sand = tex%sand / wsum
      tex%silt = tex%silt / wsum
      tex%clay = tex%clay / wsum
      tex%org_carbon = tex%org_carbon / wsum
      tex%cec = tex%cec / wsum
      tex%ph = tex%ph / wsum
      tex%bulk_density = tex%bulk_density / wsum
      tex%coarse_frag = tex%coarse_frag / wsum
   end function map_layer_texture

   !---~---
   !   The profile's builders: pure, state-free constructors that fill a soil profile from
   ! geometry and texture; the tests build their profiles with them too.
   !---~---

   !---~---
   !   A profile's layer bottoms [m, positive down]: `list` (SOIL_LAYER_TEMPLATE or a
   ! config list) cut at `depth`. The last layer ends at the depth; a remainder thinner
   ! than LAYER_MERGE_SHARE of the list layer it falls in joins the layer above. A list
   ! that ends above the depth gets one more layer down to it (the config refuses one).
   !---~---
   pure subroutine soil_layer_grid(depth, list, n_list, bottoms, n_layer)
      real(wp),    intent(in)  :: depth, list(:)
      integer(ik), intent(in)  :: n_list
      real(wp),    intent(out) :: bottoms(n_soil_layer_max)
      integer(ik), intent(out) :: n_layer
      integer(ik) :: k
      real(wp)    :: top
      bottoms = 0.0_wp
      n_layer = 0_ik
      top     = 0.0_wp
      do k = 1_ik, n_list
         if (list(k) >= depth - GRID_TOL) exit                 ! list layer k holds the cut
         n_layer = n_layer + 1_ik
         bottoms(n_layer) = list(k)
         top = list(k)
      end do
      if (k <= n_list .and. n_layer > 0_ik) then
         if (depth - top < LAYER_MERGE_SHARE * (list(k) - top)) n_layer = n_layer - 1_ik
      end if
      n_layer = n_layer + 1_ik
      bottoms(n_layer) = depth
   end subroutine soil_layer_grid

   !---~---
   !   A soil profile's geometry and hydraulics: the negative-z geometry from the layer
   ! bottoms (soil_layer_grid), each layer's own texture, and what derives from them: the
   ! field capacity, and the water network's layer curves and face weights. The arrays
   ! hold one value per active layer: from a soil map's layers, or one value broadcast
   ! (build_profile_hydraulics). It starts the profile afresh; build_profile_thermal_layers
   ! adds the thermal texture. The roots are the cohorts' own
   ! (ed_plant_hydraulics::cohort_root_profile), not the profile's.
   !---~---
   pure subroutine build_profile_hydraulics_layers(n_active, retention, layer_bottoms, theta_sat,       &
                                                   theta_res, ksat, par_a, par_n, psi_fc_m, profile)
      integer(ik),          intent(in)  :: n_active, retention
      real(wp),             intent(in)  :: layer_bottoms(:)    ! [m] bottoms of the active layers, positive down
      real(wp),             intent(in)  :: theta_sat(:), theta_res(:), ksat(:), par_a(:), par_n(:)
      real(wp),             intent(in)  :: psi_fc_m
      type(soil_profile_t), intent(out) :: profile
      integer(ik) :: k

      profile%n_active  = n_active
      profile%retention = retention

      !--- Interface elevations soil_layer_z(k) <= 0, k = 1 at the surface.
      profile%soil_layer_z(1)              = 0.0_wp
      profile%soil_layer_z(2:n_active + 1) = -layer_bottoms(1:n_active)

      !--- Thicknesses, node elevations, internode spacings (dz, dz_node > 0).
      do k = 1_ik, n_active
         profile%dz(k)     = profile%soil_layer_z(k) - profile%soil_layer_z(k+1)
         profile%z_node(k) = 0.5_wp * (profile%soil_layer_z(k) + profile%soil_layer_z(k+1))
      end do
      do k = 1_ik, n_active - 1_ik
         profile%dz_node(k)     = profile%z_node(k) - profile%z_node(k+1)
         profile%face_weight(k) = profile%dz(k) / (profile%dz(k) + profile%dz(k+1))
      end do
      profile%dz_node(n_active) = profile%dz(n_active)      ! unused at the bottom face; kept > 0

      !--- Each layer's texture.
      profile%theta_sat(1:n_active) = theta_sat(1:n_active)
      profile%theta_res(1:n_active) = theta_res(1:n_active)
      profile%ksat(1:n_active)      = ksat(1:n_active)
      profile%curve_par_a(1:n_active) = par_a(1:n_active)
      profile%curve_par_n(1:n_active) = par_n(1:n_active)

      !--- Derived threshold (from each layer's retention curve).
      profile%theta_fc(1:n_active) = soil_head_to_theta(retention, psi_fc_m, theta_sat(1:n_active),    &
                                        theta_res(1:n_active), par_a(1:n_active), par_n(1:n_active))

      !--- Each layer's curve in the water network's unknown.
      do k = 1_ik, n_active
         call build_soil_layer_curve(retention, theta_sat(k), theta_res(k), par_a(k), par_n(k), ksat(k),    &
                                     profile%dz(k), profile%layer_curve(k))
      end do
   end subroutine build_profile_hydraulics_layers

   !--- One texture for every layer: the per-layer build with each value broadcast.
   pure subroutine build_profile_hydraulics(n_active, retention, layer_bottoms, theta_sat, theta_res,  &
                                            ksat, par_a, par_n, psi_fc_m, profile)
      integer(ik),          intent(in)  :: n_active, retention
      real(wp),             intent(in)  :: layer_bottoms(:)    ! [m] bottoms of the active layers, positive down
      real(wp),             intent(in)  :: theta_sat, theta_res, ksat, par_a, par_n, psi_fc_m
      type(soil_profile_t), intent(out) :: profile
      real(wp) :: ts(n_active), tr(n_active), ks(n_active), pa(n_active), pn(n_active)
      ts = theta_sat
      tr = theta_res
      ks = ksat
      pa = par_a
      pn = par_n
      call build_profile_hydraulics_layers(n_active, retention, layer_bottoms, ts, tr, ks, pa, pn,        &
                                           psi_fc_m, profile)
   end subroutine build_profile_hydraulics

   !---~---
   !   A soil profile's thermal texture, one value per active layer, added to the
   ! geometry and hydraulics build_profile_hydraulics_layers set.
   !---~---
   pure subroutine build_profile_thermal_layers(n_active, k_solid, k_dry, dry_cvol, profile)
      integer(ik),          intent(in)    :: n_active
      real(wp),             intent(in)    :: k_solid(:), k_dry(:), dry_cvol(:)
      type(soil_profile_t), intent(inout) :: profile
      profile%solid_conductivity(1:n_active) = k_solid(1:n_active)
      profile%dry_conductivity(1:n_active)   = k_dry(1:n_active)
      profile%dry_heat_capacity(1:n_active)  = dry_cvol(1:n_active)
   end subroutine build_profile_thermal_layers

   !--- One thermal texture for every layer.
   pure subroutine build_profile_thermal(n_active, k_solid, k_dry, dry_cvol, profile)
      integer(ik),          intent(in)    :: n_active
      real(wp),             intent(in)    :: k_solid, k_dry, dry_cvol
      type(soil_profile_t), intent(inout) :: profile
      real(wp) :: ks(n_active), kd(n_active), cv(n_active)
      ks = k_solid
      kd = k_dry
      cv = dry_cvol
      call build_profile_thermal_layers(n_active, ks, kd, cv, profile)
   end subroutine build_profile_thermal

end module ed_soil_derive
