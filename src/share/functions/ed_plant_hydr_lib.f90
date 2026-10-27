! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_plant_hydr_lib -- the water relations of a plant's tissues and xylem: how
! much water a tissue holds at a given potential, its pressure-volume curve
! (docs/science/plant_hydraulics.md §2), the vulnerability curve of its xylem (§3),
! and the conductances its sapwood and roots have from their size (§3-§4).
!
!   The read-offs (tissue_psi_to_rwc, tissue_water_to_psi, ...) carry floors that keep
! a diagnostic finite. The solver's relations (tissue_water_relations) have none and
! come with their derivatives in the solver's unknown (soil_plant_water.md §2). The
! flows that use these laws, sapflow and root uptake, are ed_plant_hydraulics'; the
! soil's side is ed_soil_hydr_lib.
!---~---
module ed_plant_hydr_lib
   use ed_kinds,       only : wp
   use ed_constants,   only : grav_head, rho_h2o, pi, tiny_num
   implicit none
   private

   !--- Tissue pressure-volume curve: read-offs and the solver's relations.
   public :: water_curve_t
   public :: tissue_psi_tlp, tissue_psi_to_rwc, tissue_rwc_to_psi
   public :: tissue_psi_to_water, tissue_capacitance, tissue_water_to_psi, tissue_capped_water
   public :: tissue_water_relations, tissue_rwc_to_unknown, tissue_water_to_rwc
   !--- Vulnerability curve; the sapwood's and the roots' conductances.
   public :: conductance_retained_fraction, KEXP_TOL
   public :: xylem_max_conductance, root_fraction_profile, rhizosphere_conductance

   !---~---
   !   One tissue's pressure-volume curve: the four traits every tissue-water routine
   ! reads, passed together with the tissue's biomass.
   !---~---
   type :: water_curve_t
      real(wp) :: pi0           = 0.0_wp   ! [MPa] osmotic potential at full turgor (< 0)
      real(wp) :: elastic_mod   = 0.0_wp   ! [MPa] bulk elastic modulus
      real(wp) :: apoplast_frac = 0.0_wp   ! [-]   share of the saturated water outside the living cells
      real(wp) :: water_sat     = 0.0_wp   ! [kg H2O / kgC] water per unit biomass at saturation
   end type water_curve_t

   !---~---
   !   [-] Floor on the symplastic relative water content in the read-offs; the
   ! solver's relations have none.
   !---~---
   real(wp), parameter :: RWC_FLOOR = 1.0e-4_wp

   !---~---
   !   [kg/s/MPa per plant] Intact conductance of a plant without sapwood, so its leaves stay
   ! joined to its wood.
   !---~---
   real(wp), parameter :: XYLEM_COND_MIN = 1.0e-15_wp

   !---~---
   !   [-] How close the vulnerability exponent kexp must be to 1 or 2 to take that exponent's
   ! closed forms.
   !---~---
   real(wp), parameter :: KEXP_TOL = 1.0e-9_wp

   !---~---
   !   [-] Floor of half-spacing / root radius in the single-root (Gardner)
   ! conductance: closer, the roots all but touch and the form fails.
   !---~---
   real(wp), parameter :: MIN_SPACING_RATIO = 2.0_wp

contains

   !---~---
   !   Tissue pressure-volume curve (Bartlett / Tyree-Hammel;
   ! docs/science/plant_hydraulics.md §2, eqs. (3)-(6)): the read-offs diagnostics and
   ! a store's set-up use, floored at RWC_FLOOR so they stay finite.
   !---~---

   !---~---
   !   Turgor loss point [MPa] (Bartlett et al. 2012 eq. 1;
   ! docs/science/plant_hydraulics.md eq. (3)).
   !---~---
   elemental real(wp) function tissue_psi_tlp(pi0, elastic_mod) result(psi_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
   end function tissue_psi_tlp

   !---~---
   !   Symplastic relative water content at turgor loss (Bartlett et al. 2012 eq. 2;
   ! docs/science/plant_hydraulics.md eq. (3)).
   !---~---
   elemental real(wp) function tissue_rwc_tlp(pi0, elastic_mod) result(rwc_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
   end function tissue_rwc_tlp

   !---~---
   !   Water potential [MPa] at symplastic relative water content rwc, floored at
   ! RWC_FLOOR (docs/science/plant_hydraulics.md eq. (4)).
   !---~---
   elemental real(wp) function tissue_rwc_to_psi(rwc, pi0, elastic_mod) result(psi)
      real(wp), intent(in) :: rwc, pi0, elastic_mod
      real(wp) :: rwc_tlp, r
      rwc_tlp = tissue_rwc_tlp(pi0, elastic_mod)
      r       = max(rwc, RWC_FLOOR)
      if (r >= rwc_tlp) then
         psi = elastic_mod*(r - rwc_tlp) + pi0/r   ! turgor + osmotic
      else
         psi = pi0/r                               ! turgor lost
      end if
   end function tissue_rwc_to_psi

   !---~---
   !   Symplastic relative water content at water potential psi, floored at RWC_FLOOR
   ! (docs/science/plant_hydraulics.md eq. (4)).
   !---~---
   elemental real(wp) function tissue_psi_to_rwc(psi, pi0, elastic_mod) result(rwc)
      real(wp), intent(in) :: psi, pi0, elastic_mod
      real(wp) :: psi_tlp, b, disc
      psi_tlp = tissue_psi_tlp(pi0, elastic_mod)
      if (psi >= psi_tlp) then                ! turgid
         !--- The larger root of elastic_mod R^2 - b R + pi0 = 0.
         b    = psi + elastic_mod + pi0
         disc = max(b*b - 4.0_wp*elastic_mod*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*elastic_mod)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, RWC_FLOOR)
   end function tissue_psi_to_rwc

   !---~---
   !   Total tissue water [kg per plant] at potential psi: the symplast's saturated
   ! water times R, plus the apoplast's, a constant reservoir
   ! (docs/science/plant_hydraulics.md eq. (5)).
   !---~---
   elemental real(wp) function tissue_psi_to_water(psi, curve, biomass) result(w)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, r
      w_sat = curve%water_sat*biomass
      r     = tissue_psi_to_rwc(psi, curve%pi0, curve%elastic_mod)
      w     = (1.0_wp - curve%apoplast_frac)*w_sat*r + curve%apoplast_frac*w_sat
   end function tissue_psi_to_water

   !---~---
   !   Capacitance dW/dpsi [kg/MPa per plant] at potential psi
   ! (docs/science/plant_hydraulics.md eqs. (5)-(6)).
   !---~---
   elemental real(wp) function tissue_capacitance(psi, curve, biomass) result(c)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      associate (pi0 => curve%pi0, elastic_mod => curve%elastic_mod)
         psi_tlp   = tissue_psi_tlp(pi0, elastic_mod)
         w_sat_sym = (1.0_wp - curve%apoplast_frac)*curve%water_sat*biomass
         r         = tissue_psi_to_rwc(psi, pi0, elastic_mod)
         if (psi >= psi_tlp) then
            cr = 1.0_wp / (elastic_mod - pi0/(r*r))   ! dR/dpsi, turgid
         else
            cr = -(r*r)/pi0                           ! dR/dpsi, flaccid
         end if
      end associate
      c = w_sat_sym*cr
   end function tissue_capacitance

   !---~---
   !   Water potential [MPa] of a tissue holding w [kg per plant], the exact inverse of
   ! tissue_psi_to_water (docs/science/plant_hydraulics.md eq. (5)).
   !---~---
   elemental real(wp) function tissue_water_to_psi(w, curve, biomass) result(psi)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, rwc
      w_sat = curve%water_sat*biomass
      rwc   = (w - curve%apoplast_frac*w_sat) / max((1.0_wp - curve%apoplast_frac)*w_sat, tiny_num)
      psi   = tissue_rwc_to_psi(rwc, curve%pi0, curve%elastic_mod)
   end function tissue_water_to_psi

   !---~---
   !   Tissue water [kg per plant] capped at saturation, water_sat * biomass, for a
   ! store that lost biomass without its water; the caller sheds the excess (w
   ! result) to the ground (docs/science/soil_plant_water.md §2).
   !---~---
   elemental real(wp) function tissue_capped_water(w, curve, biomass) result(w_capped)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      w_capped = min(w, curve%water_sat*biomass)
   end function tissue_capped_water

   !---~---
   !   Tissue water relations for the soil-plant water solver, without the read-offs'
   ! floors.
   !---~---

   !---~---
   !   A tissue at symplastic relative water content r: potential [MPa] and water
   ! [kg/m2 of ground] for nplant plants of `biomass` each, with their derivatives in
   ! r. The pressure-volume curve without tissue_rwc_to_psi's floor: the solver keeps
   ! r above zero (docs/science/plant_hydraulics.md eqs. (4)-(5)).
   !---~---
   elemental subroutine tissue_water_relations(curve, biomass, nplant, r, psi, dpsi, w, dw)
      type(water_curve_t), intent(in)  :: curve
      real(wp),            intent(in)  :: biomass   ! [kgC per plant]
      real(wp),            intent(in)  :: nplant    ! [plants/m2]
      real(wp),            intent(in)  :: r         ! [-] symplastic relative water content
      real(wp),            intent(out) :: psi, dpsi, w, dw
      real(wp) :: r_tlp, w_sat
      associate (pi0 => curve%pi0, emod => curve%elastic_mod, apo => curve%apoplast_frac)
         r_tlp = tissue_rwc_tlp(pi0, emod)
         if (r >= r_tlp) then
            psi  = emod * (r - r_tlp) + pi0 / r
            dpsi = emod - pi0 / (r * r)
         else
            psi  = pi0 / r
            dpsi = -pi0 / (r * r)
         end if
         w_sat = nplant * curve%water_sat * biomass
         w     = w_sat * ((1.0_wp - apo) * r + apo)
         dw    = w_sat * (1.0_wp - apo)
      end associate
   end subroutine tissue_water_relations

   !---~---
   !   Symplastic relative water content of a tissue holding w [kg per plant], the
   ! inverse of tissue_water_relations's water. Undefined with no biomass: the caller
   ! then starts the tissue from a potential instead (tissue_psi_to_rwc).
   !---~---
   elemental function tissue_water_to_rwc(curve, biomass, w) result(r)
      type(water_curve_t), intent(in) :: curve
      real(wp),            intent(in) :: biomass, w
      real(wp) :: r, w_sat
      w_sat = curve%water_sat * biomass
      r = (w / w_sat - curve%apoplast_frac) / (1.0_wp - curve%apoplast_frac)
   end function tissue_water_to_rwc

   !---~---
   !   A tissue's unknown is its relative water content, kept strictly positive: zero
   ! is an infinitely dry tissue.
   !---~---
   elemental function tissue_rwc_to_unknown(r) result(tau)
      real(wp), intent(in) :: r
      real(wp) :: tau
      real(wp), parameter :: R_START = 1.0e-6_wp
      tau = max(r, R_START)   ! clamp-ok: an emptied store starts just inside the curve
   end function tissue_rwc_to_unknown

   !--- Vulnerability curve (docs/science/plant_hydraulics.md §3, eq. (7)).

   !---~---
   !   Retained fraction of xylem conductance, 1 - PLC = 1 / (1 + r^kexp), r
   ! psi/psi50, and one at or above zero potential (docs/science/plant_hydraulics.md
   ! eq. (7)). r = 0 is kept off the general power, which a compiler may take through
   ! log(0) under floating-point trapping; exponent 2 uses a product, cheaper than
   ! the general power.
   !---~---
   elemental real(wp) function conductance_retained_fraction(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      if (r <= 0.0_wp) then
         f = 1.0_wp
      else if (abs(kexp - 2.0_wp) < KEXP_TOL) then
         f = 1.0_wp / (1.0_wp + r * r)
      else
         f = 1.0_wp / (1.0_wp + r**kexp)
      end if
   end function conductance_retained_fraction

   !--- Xylem (docs/science/plant_hydraulics.md eqs. (9)-(10)).

   !---~---
   !   One plant's intact-xylem conductance [kg/s/MPa]: the sapwood's specific
   ! conductivity over its cross-section and its path, height x vessel_curl, floored
   ! for a plant with no sapwood yet (docs/science/plant_hydraulics.md eq. (10)).
   !---~---
   elemental function xylem_max_conductance(wood_kmax, sap_area, height, vessel_curl) result(k)
      real(wp), intent(in) :: wood_kmax    ! [kg/m/s/MPa] sapwood specific conductivity
      real(wp), intent(in) :: sap_area     ! [m2] sapwood cross-section
      real(wp), intent(in) :: height       ! [m]
      real(wp), intent(in) :: vessel_curl  ! [-] path length per unit height
      real(wp) :: k
      k = max(wood_kmax * sap_area / max(height * vessel_curl, tiny_num), XYLEM_COND_MIN)
   end function xylem_max_conductance

   !--- Roots (docs/science/plant_hydraulics.md eqs. (13)-(15)).

   !---~---
   !   The share of a cohort's roots in the soil layer between depths z_top and z_bot,
   ! beta^(z_top/D) - beta^(z_bot/D) with D = root_depth and the depths clamped to
   ! [0, D], so a layer below the rooting depth gets none. Over [0, D] the shares sum
   ! to 1 - beta, which the caller divides by (docs/science/plant_hydraulics.md eq.
   ! (13)).
   !---~---
   elemental real(wp) function root_fraction_profile(root_beta, root_depth, z_top, z_bot) result(frac)
      real(wp), intent(in) :: root_beta   ! [-]  root-profile decay (0,1); smaller => shallower
      real(wp), intent(in) :: root_depth  ! [m]  maximum rooting depth (> 0)
      real(wp), intent(in) :: z_top       ! [m]  depth of the layer top    (>= 0)
      real(wp), intent(in) :: z_bot       ! [m]  depth of the layer bottom (> z_top)
      real(wp) :: inv_d, a, b
      inv_d = 1.0_wp / max(root_depth, tiny_num)
      a = min(max(z_top, 0.0_wp), root_depth) * inv_d
      b = min(max(z_bot, 0.0_wp), root_depth) * inv_d
      frac = root_beta**a - root_beta**b
   end function root_fraction_profile

   !---~---
   !   A cohort's rhizosphere conductance in one soil layer per unit soil conductivity
   ! [kg/s/MPa per m2 of ground, per m/s]: steady radial flow to a single root
   ! (Gardner 1960), 2 pi L / ln(r_half / r_root) for each plant, made a potential
   ! gradient by rho/grav_head. It depends on the roots alone, so it is set once a
   ! step (docs/science/plant_hydraulics.md eq. (14)).
   !---~---
   elemental function rhizosphere_conductance(nplant, root_length, half_spacing, root_radius)      &
                      result(g_rhizo)
      real(wp), intent(in) :: nplant        ! [plants/m2]
      real(wp), intent(in) :: root_length   ! [m per plant] fine-root length in the layer
      real(wp), intent(in) :: half_spacing  ! [m] half the distance between the layer's roots
      real(wp), intent(in) :: root_radius   ! [m] absorbing-root radius
      real(wp) :: g_rhizo
      g_rhizo = nplant * 2.0_wp * pi * root_length * rho_h2o                                      &
               / (grav_head * log(max(half_spacing / root_radius, MIN_SPACING_RATIO)))
   end function rhizosphere_conductance

end module ed_plant_hydr_lib
