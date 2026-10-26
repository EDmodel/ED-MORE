! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_hydr_lib -- the water relations of the soil: a layer's van Genuchten curve and !
! saturated conductivity from what the soil maps report (docs/science/soil_biophysics.md   !
! §Hydraulic parameters from soil texture, eqs. (6)-(8)), its retention curve and          !
! conductivity (soil_biophysics.md eq. (2)), and the pond's depression store.              !
!                                                                                          !
! The read-offs (soil_head_to_theta, soil_theta_to_head) take a state and return a         !
! potential or a water content, with floors that keep a diagnostic finite. The water       !
! relations (soil_water_relations, pond_water_relations) give the soil-plant water network !
! (water_derivs) a store's potential, water and conductivity with their derivatives in     !
! the variable it solves for, with no floors on the curve itself: only a dry layer's start !
! (SE_MIN), a very dry head (LOG_HEAD_MAX) and an underflowing slope are guarded           !
! (docs/science/soil_plant_water.md §2). Flow                                              !
! between soil layers is in meds_soil_water; the plant's side is meds_plant_hydr_lib.      !
!==========================================================================================!
module meds_soil_hydr_lib
   use meds_kinds,     only : wp, ik
   use meds_constants, only : grav_head, rho_h2o, day_sec, tiny_num, grav, t_3ple, latent_heat_fusion
   implicit none
   private

   !----- A layer's curve from the soil maps. ---------------------------------------------!
   public :: soil_texture_t, soil_texture_to_hydraulics
   !----- Retention curve: theta(head), head(theta) and potential(theta) read-offs. -------!
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, VG_PORE_CONNECTIVITY
   public :: soil_head_to_theta, soil_theta_to_head, soil_matric_potential, soil_theta_to_saturation
   public :: soil_liquid_capacity
   !----- The water relations the soil-plant water solver uses, in the variable it solves  !
   !      for: a soil layer's saturation (dry) or potential (wet), the pond's head. -------!
   public :: SOIL_SPECIFIC_STORAGE, POND_DEPRESSION_DEPTH
   public :: soil_layer_curve_t, build_soil_layer_curve, soil_water_relations, soil_theta_to_unknown
   public :: soil_ice_t, soil_ice_state
   public :: pond_water_relations, pond_water_to_head

   !----- One soil layer as the soil maps describe it (SoilGrids 2.0's quantities). -------!
   type :: soil_texture_t
      real(wp) :: sand         = 0.0_wp   !< [%]           sand share of the fine earth
      real(wp) :: silt         = 0.0_wp   !< [%]           silt share of the fine earth
      real(wp) :: clay         = 0.0_wp   !< [%]           clay share of the fine earth
      real(wp) :: org_carbon   = 0.0_wp   !< [g/kg]        organic carbon of the fine earth
      real(wp) :: cec          = 0.0_wp   !< [cmol(c)/kg]  cation exchange capacity
      real(wp) :: ph           = 0.0_wp   !< [-]           pH in water
      real(wp) :: bulk_density = 0.0_wp   !< [g/cm3]       bulk density of the fine earth
      real(wp) :: coarse_frag  = 0.0_wp   !< [% volume]    coarse fragments (> 2 mm) in the layer
   end type soil_texture_t

   !----- Hodnett & Tomasella (2002) give alpha per kPa. ----------------------------------!
   real(wp), parameter :: KPA_PER_MPA      = 1000.0_wp                !< [kPa/MPa]
   !----- Ottoni et al. (2019): Ks = a (theta_sat - theta(-33 kPa))^b, Ks in m/day. -------!
   real(wp), parameter :: HEAD_33KPA       = -0.033_wp / grav_head   !< [m] -33 kPa as water head
   real(wp), parameter :: OTTONI_KS_COEF   = 19.31_wp                !< [m/day]
   real(wp), parameter :: OTTONI_KS_EXP    = 1.948_wp                !< [-]
   !----- Density of the mineral grains and of the fragments (Saxton & Rawls 2006). -------!
   real(wp), parameter :: PARTICLE_DENSITY = 2.65_wp                 !< [g/cm3]

   !----- Retention-curve selector codes (soil_profile_t%retention). -----------------------!
   integer(ik), parameter :: SOIL_RETENTION_VG       = 1_ik  !< van Genuchten-Mualem (default)
   integer(ik), parameter :: SOIL_RETENTION_CAMPBELL = 2_ik  !< Campbell / Clapp-Hornberger (option)

   !----- Mualem's pore-connectivity exponent L in the van Genuchten conductivity          !
   !      (docs/science/soil_biophysics.md eq. (2)), set here only: -1 after Schaap & Leij !
   !      (2000, SSSA J. 64:843), which keeps a soil of n near 1.4 conducting down to      !
   !      about -2 MPa. -------------------------------------------------------------------!
   real(wp), parameter :: VG_PORE_CONNECTIVITY = -1.0_wp

   !----- [-] Floor on the effective saturation in the read-offs; the solver's curves have !
   !      none. ---------------------------------------------------------------------------!
   real(wp), parameter :: SE_MIN = 1.0e-9_wp
   !----- The largest logarithm the dry-soil head is evaluated at: e^700 m of head is far  !
   !      beyond any physical soil and still finite in double precision. ------------------!
   real(wp), parameter :: LOG_HEAD_MAX = 700.0_wp

   !----- [1/m] Water a saturated soil layer gains per metre of pressure head, from the    !
   !      compressibility of the soil matrix and of water: a typical value for clays and   !
   !      loams (Freeze and Cherry 1979, table 2.5). It acts only above saturation         !
   !      (docs/science/soil_plant_water.md §2). ------------------------------------------!
   real(wp), parameter :: SOIL_SPECIFIC_STORAGE = 1.0e-4_wp

   !----- One soil layer's curve, set up once per layer: the retention parameters and the  !
   !      point where the solver's unknown changes from saturation to potential. ----------!
   type :: soil_layer_curve_t
      integer(ik) :: retention    = 1_ik    !< SOIL_RETENTION_VG or SOIL_RETENTION_CAMPBELL
      real(wp)    :: theta_sat    = 0.0_wp  !< [m3/m3]
      real(wp)    :: theta_res    = 0.0_wp  !< [m3/m3]
      real(wp)    :: par_a        = 0.0_wp  !< vG alpha [1/m], or Campbell's air-entry head [m, < 0]
      real(wp)    :: par_n        = 0.0_wp  !< vG n [-], or Campbell's b [-]
      real(wp)    :: ksat         = 0.0_wp  !< [m/s] saturated conductivity
      real(wp)    :: dz           = 0.0_wp  !< [m] layer thickness
      real(wp)    :: se_switch    = 1.0_wp  !< [-] saturation where the unknown changes branch
      real(wp)    :: head_switch  = 0.0_wp  !< [m] the matric head there
      real(wp)    :: slope_switch = 1.0_wp  !< [1/m] dSe/dhead there, the unknown's scale on the wet side
   end type soil_layer_curve_t

   !----- A layer's ice, held through a water step (docs/science/soil_biophysics.md,      !
   !      "Frozen soil"). The liquid takes the curve's potential at its own content, as   !
   !      if the ice were air, so in a frozen layer it sits at the freezing curve's       !
   !      potential; the ice slows the flow by the impedance 10^(-ICE_IMPEDANCE_EXP       !
   !      theta_ice/theta_sat) (CLM5, Swenson et al. 2012). Once the liquid fills the     !
   !      pores the ice leaves, the liquid's head rises from the head they filled at with !
   !      no more water, and above zero through SOIL_SPECIFIC_STORAGE, as a thawed        !
   !      layer's does. A record at its defaults is a layer with no ice. -----------------!
   type :: soil_ice_t
      real(wp) :: theta_ice  = 0.0_wp          !< [m3/m3] the ice
      real(wp) :: impedance  = 1.0_wp          !< [-] its factor on the conductivity
      real(wp) :: tau_full   = huge(1.0_wp)    !< the unknown at which the liquid fills the pores
      real(wp) :: head_full  = 0.0_wp          !< [m] the head there
      real(wp) :: dhead_full = 0.0_wp          !< [m per unit tau] the head's slope there, kept beyond
      real(wp) :: k_full     = 0.0_wp          !< [m/s] the liquid's conductivity there, unimpeded
   end type soil_ice_t

   !----- [-] CLM5's ice impedance exponent (Swenson et al. 2012, e_ice). --------------------!
   real(wp), parameter :: ICE_IMPEDANCE_EXP = 6.0_wp

   !----- [m] Depth of the surface's small hollows, which the pond fills before a head     !
   !      builds, so its storage stays smooth as it fills or empties                       !
   !      (docs/science/soil_plant_water.md §2). ------------------------------------------!
   real(wp), parameter :: POND_DEPRESSION_DEPTH = 1.0e-3_wp

contains

   !=======================================================================================!
   !     A layer's curve from the soil maps (docs/science/soil_biophysics.md §Hydraulic    !
   !     parameters from soil texture): the fine earth's van Genuchten curve after Hodnett !
   !     & Tomasella (2002), its saturated conductivity after Ottoni et al. (2019), and    !
   !     coarse fragments that hold no water and slow the flow (Saxton & Rawls 2006).      !
   !=======================================================================================!

   !----- The van Genuchten curve and Ks of a soil layer (docs/science/soil_biophysics.md  !
   !      eqs. (6)-(8)). Alpha comes out per kPa and is turned into per metre of water     !
   !      head, the unit of the model's soil potential. -----------------------------------!
   elemental subroutine soil_texture_to_hydraulics(texture, theta_sat, theta_res, vg_alpha, vg_n, ksat)
      type(soil_texture_t), intent(in)  :: texture
      real(wp),             intent(out) :: theta_sat   !< [m3/m3] water content at saturation, per layer volume
      real(wp),             intent(out) :: theta_res   !< [m3/m3] residual water content, per layer volume
      real(wp),             intent(out) :: vg_alpha    !< [1/m]   van Genuchten alpha
      real(wp),             intent(out) :: vg_n        !< [-]     van Genuchten n
      real(wp),             intent(out) :: ksat        !< [m/s]   saturated hydraulic conductivity of the layer
      real(wp) :: oc_percent, theta_33, drainable, rock_volume, rock_mass, matrix_density_ratio

      associate(sand => texture%sand, silt => texture%silt, clay => texture%clay, cec => texture%cec,  &
                ph => texture%ph, bulk_density => texture%bulk_density)
         oc_percent = texture%org_carbon / 10.0_wp

         !----- The fine earth's curve (Hodnett & Tomasella 2002, Table 8). ---------------!
         vg_alpha  = exp((-2.294_wp - 3.526_wp*silt + 2.440_wp*oc_percent - 0.076_wp*cec - 11.331_wp*ph  &
                          + 0.019_wp*silt**2) / 100.0_wp) * grav_head * KPA_PER_MPA
         vg_n      = exp((62.986_wp - 0.833_wp*clay - 0.529_wp*oc_percent + 0.593_wp*ph                 &
                          + 0.007_wp*clay**2 - 0.014_wp*sand*silt) / 100.0_wp)
         theta_sat = (81.799_wp + 0.099_wp*clay - 31.42_wp*bulk_density + 0.018_wp*cec + 0.451_wp*ph     &
                      - 0.0005_wp*sand*clay) / 100.0_wp
         !----- The regression dips below zero on sandy, alkaline, low-CEC soils, which    !
         !      hold no residual water. ---------------------------------------------------!
         theta_res = max((22.733_wp - 0.164_wp*sand + 0.235_wp*cec - 0.831_wp*ph + 0.0018_wp*clay**2  &
                          + 0.0026_wp*sand*clay) / 100.0_wp, 0.0_wp)

         !----- Ks of the fine earth from the water that drains between saturation, at     !
         !      zero head on this curve, and -33 kPa (Ottoni et al. 2019). ----------------!
         theta_33  = soil_head_to_theta(SOIL_RETENTION_VG, HEAD_33KPA, theta_sat, theta_res, vg_alpha, vg_n)
         drainable = max(theta_sat - theta_33, tiny_num)
         ksat      = OTTONI_KS_COEF / day_sec * drainable**OTTONI_KS_EXP

         !----- Coarse fragments hold no water and slow the flow: Saxton & Rawls (2006)    !
         !      eq. 22, on the fragments' share of the mass, which eq. 19 gives from their !
         !      share of the volume. ------------------------------------------------------!
         rock_volume          = texture%coarse_frag / 100.0_wp
         matrix_density_ratio = bulk_density / PARTICLE_DENSITY
         rock_mass = rock_volume / max(matrix_density_ratio*(1.0_wp - rock_volume) + rock_volume, tiny_num)
         theta_sat = theta_sat * (1.0_wp - rock_volume)
         theta_res = theta_res * (1.0_wp - rock_volume)
         ksat      = ksat * (1.0_wp - rock_mass) / (1.0_wp - rock_mass*(1.0_wp - 1.5_wp*matrix_density_ratio))
      end associate
   end subroutine soil_texture_to_hydraulics

   !=======================================================================================!
   !     Soil retention curves: van Genuchten-Mualem (default) or Campbell                 !
   !     (docs/science/soil_biophysics.md eq. (2)). par_a and par_n are alpha [1/m] and n  !
   !     for van Genuchten, the air-entry head [m, < 0] and b for Campbell; heads are in   !
   !     metres of water, below zero when unsaturated. Every soil-water routine evaluates  !
   !     the curve through the four kernels below, so a diagnostic potential is the        !
   !     solver's own.                                                                     !
   !=======================================================================================!

   !----- Volumetric water content [m3/m3] at matric head [m]: theta_sat at or above zero, !
   !      and the saturation floored at SE_MIN. -------------------------------------------!
   elemental function soil_head_to_theta(retention, head, theta_sat, theta_res, par_a, par_n)  &
                      result(theta)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: head, theta_sat, theta_res, par_a, par_n
      real(wp)                :: theta, se, dse_dhead
      if (head >= 0.0_wp) then
         theta = theta_sat
         return
      end if
      call head_to_saturation(retention, par_a, par_n, head, se, dse_dhead)
      se    = min(max(se, SE_MIN), 1.0_wp)
      theta = theta_res + (theta_sat - theta_res) * se
   end function soil_head_to_theta

   !----- The soil freezing curve: the liquid water a layer can hold at temperature temp  !
   !      [m3/m3], and its slope [m3/m3/K]. Below t_3ple the layer's own retention curve  !
   !      at the freezing-point-depression head L_f (T - t_3ple) / (g T), the head at     !
   !      which liquid water and ice are in equilibrium (Clausius-Clapeyron; Niu & Yang   !
   !      2006, CLM5's supercooled liquid); saturation at and above it                    !
   !      (docs/science/soil_biophysics.md, "Frozen soil"). ------------------------------!
   elemental subroutine soil_liquid_capacity(curve, temp, theta_liq, dtheta_dtemp)
      type(soil_layer_curve_t), intent(in)  :: curve          !< the layer's curve
      real(wp),                 intent(in)  :: temp           !< [K]
      real(wp),                 intent(out) :: theta_liq      !< [m3/m3]   the liquid it can hold
      real(wp),                 intent(out) :: dtheta_dtemp   !< [m3/m3/K] its slope
      real(wp) :: head, se, dse_dhead
      if (temp >= t_3ple) then
         theta_liq    = curve%theta_sat
         dtheta_dtemp = 0.0_wp
         return
      end if
      head = latent_heat_fusion * (temp - t_3ple) / (grav * temp)
      call head_to_saturation(curve%retention, curve%par_a, curve%par_n, head, se, dse_dhead)
      if (se <= SE_MIN) dse_dhead = 0.0_wp
      se = min(max(se, SE_MIN), 1.0_wp)
      theta_liq    = curve%theta_res + (curve%theta_sat - curve%theta_res) * se
      dtheta_dtemp = (curve%theta_sat - curve%theta_res) * dse_dhead * latent_heat_fusion * t_3ple       &
                     / (grav * temp * temp)
   end subroutine soil_liquid_capacity

   !----- Relative saturation [-] of a layer holding theta [m3/m3]: where it sits between    !
   !      residual and saturation, in [0, 1], so a solver excursion outside the curve's domain  !
   !      reads as exactly dry or saturated. ---------------------------------------------------!
   elemental function soil_theta_to_saturation(theta, theta_sat, theta_res) result(se)
      real(wp), intent(in) :: theta, theta_sat, theta_res
      real(wp)             :: se
      se = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      se = min(max(se, 0.0_wp), 1.0_wp)
   end function soil_theta_to_saturation

   !----- Matric head [m, <= 0] of a layer holding theta [m3/m3], its saturation floored   !
   !      at SE_MIN; zero at saturation. --------------------------------------------------!
   elemental function soil_theta_to_head(retention, theta, theta_sat, theta_res, par_a, par_n) &
                      result(head)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n
      real(wp)                :: head, se, dse_dhead
      se = max(soil_theta_to_saturation(theta, theta_sat, theta_res), SE_MIN)
      if (se >= 1.0_wp) then
         head = 0.0_wp
         return
      end if
      call saturation_to_head(retention, par_a, par_n, se, head, dse_dhead)
   end function soil_theta_to_head

   !----- Matric potential [MPa] of a layer holding theta [m3/m3]: its head through        !
   !      grav_head, the constant plant hydraulics uses, so a soil potential and a root's  !
   !      are on one scale. ---------------------------------------------------------------!
   elemental function soil_matric_potential(retention, theta, theta_sat, theta_res, par_a, par_n) &
                      result(psi_mpa)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n
      real(wp)                :: psi_mpa
      psi_mpa = grav_head * soil_theta_to_head(retention, theta, theta_sat, theta_res, par_a, par_n)
   end function soil_matric_potential

   !=======================================================================================!
   !     Water relations for the soil-plant water solver: soil layers.                     !
   !=======================================================================================!

   !----- Set up one layer's curve and its branch point, where the curve holds the most    !
   !      water per unit head: (alpha |head|)^n = m, Se = (1+m)^-m, for van Genuchten; the !
   !      air-entry point for Campbell (docs/science/soil_plant_water.md §2). -------------!
   pure subroutine build_soil_layer_curve(retention, theta_sat, theta_res, par_a, par_n, ksat, dz, c)
      integer(ik),              intent(in)  :: retention
      real(wp),                 intent(in)  :: theta_sat, theta_res, par_a, par_n, ksat, dz
      type(soil_layer_curve_t), intent(out) :: c
      real(wp) :: m, x
      c%retention = retention ; c%theta_sat = theta_sat ; c%theta_res = theta_res
      c%par_a = par_a ; c%par_n = par_n ; c%ksat = ksat ; c%dz = dz
      if (retention == SOIL_RETENTION_CAMPBELL) then
         c%se_switch    = 1.0_wp
         c%head_switch   = par_a
         c%slope_switch = -1.0_wp / (par_n * par_a)
      else
         m = 1.0_wp - 1.0_wp / par_n
         x = m ** (1.0_wp / par_n)
         c%se_switch    = (1.0_wp + m) ** (-m)
         c%head_switch   = -x / par_a
         c%slope_switch = par_a * m * par_n * c%se_switch * (m / (1.0_wp + m)) / x
      end if
   end subroutine build_soil_layer_curve

   !----- A soil layer's potential [MPa], water [kg/m2] and conductivity [m/s], with their !
   !      derivatives in the solver's unknown tau (docs/science/soil_plant_water.md §2).   !
   !      Dry of the branch point tau is the effective saturation Se; wet of it tau runs   !
   !      linearly with the potential, at the slope that joins the two branches smoothly.  !
   !      Above saturation the layer stores water through SOIL_SPECIFIC_STORAGE instead of !
   !      being clipped. With ice, tau measures the liquid alone: the layer's water is the !
   !      liquid's and the ice's, and its conductivity the liquid's times the ice's        !
   !      impedance (soil_ice_t). ---------------------------------------------------------!
   elemental subroutine soil_water_relations(c, ice, tau, psi, dpsi, w, dw, kcond, dkcond)
      type(soil_layer_curve_t), intent(in)  :: c
      type(soil_ice_t),         intent(in)  :: ice
      real(wp),                 intent(in)  :: tau
      real(wp),                 intent(out) :: psi, dpsi    !< [MPa], [MPa per unit tau]
      real(wp),                 intent(out) :: w, dw        !< [kg/m2], [kg/m2 per unit tau]
      real(wp),                 intent(out) :: kcond, dkcond !< [m/s], [m/s per unit tau]
      real(wp) :: head
      if (tau > ice%tau_full) then
         !----- The pores are full: the liquid's head rises from where they filled with no more  !
         !      water, and above zero the layer stores through the storage term, as a thawed one !
         !      does. ---------------------------------------------------------------------------!
         head   = ice%head_full + (tau - ice%tau_full) * ice%dhead_full
         psi    = grav_head * head
         dpsi   = grav_head * ice%dhead_full
         w      = rho_h2o * c%dz * (c%theta_sat + SOIL_SPECIFIC_STORAGE * max(head, 0.0_wp))
         dw     = 0.0_wp
         if (head > 0.0_wp) dw = rho_h2o * c%dz * SOIL_SPECIFIC_STORAGE * ice%dhead_full
         kcond  = ice%k_full * ice%impedance
         dkcond = 0.0_wp
         return
      end if
      call liquid_water_relations(c, tau, psi, dpsi, w, dw, kcond, dkcond)
      w      = w + rho_h2o * c%dz * ice%theta_ice
      kcond  = kcond * ice%impedance
      dkcond = dkcond * ice%impedance
   end subroutine soil_water_relations

   !----- soil_water_relations of a layer with no ice. -----------------------------------!
   elemental subroutine liquid_water_relations(c, tau, psi, dpsi, w, dw, kcond, dkcond)
      type(soil_layer_curve_t), intent(in)  :: c
      real(wp),                 intent(in)  :: tau
      real(wp),                 intent(out) :: psi, dpsi, w, dw, kcond, dkcond   !< as soil_water_relations's
      real(wp) :: head, dhead, se, dse_dhead, dk_dse, dk_dhead, dtheta
      dtheta = c%theta_sat - c%theta_res
      if (tau <= c%se_switch) then
         !----- Dry branch: tau is the saturation itself. ---------------------------------!
         se = tau
         call saturation_to_head(c%retention, c%par_a, c%par_n, se, head, dse_dhead)
         call saturation_to_conductivity(c%retention, c%par_n, c%ksat, se, kcond, dk_dse)
         psi    = grav_head * head
         dpsi   = grav_head / max(dse_dhead, tiny_num)   ! clamp-ok: a slope underflowed in extreme dry soil
         w      = rho_h2o * c%dz * (c%theta_res + dtheta * se)
         dw     = rho_h2o * c%dz * dtheta
         dkcond = dk_dse
      else
         !----- Wet branch: tau measures the potential. -----------------------------------!
         head  = c%head_switch + (tau - c%se_switch) / c%slope_switch
         dhead = 1.0_wp / c%slope_switch
         psi   = grav_head * head
         dpsi  = grav_head * dhead
         if (head < 0.0_wp) then
            call head_to_saturation(c%retention, c%par_a, c%par_n, head, se, dse_dhead)
            call head_to_conductivity(c%retention, c%par_a, c%par_n, c%ksat, head, kcond, dk_dhead)
            w      = rho_h2o * c%dz * (c%theta_res + dtheta * se)
            dw     = rho_h2o * c%dz * dtheta * dse_dhead * dhead
            dkcond = dk_dhead * dhead
         else
            w      = rho_h2o * c%dz * (c%theta_sat + SOIL_SPECIFIC_STORAGE * head)
            dw     = rho_h2o * c%dz * SOIL_SPECIFIC_STORAGE * dhead
            kcond  = c%ksat
            dkcond = 0.0_wp
         end if
      end if
   end subroutine liquid_water_relations

   !----- A layer's ice record (soil_ice_t) for theta_ice [m3/m3]: its impedance, and     !
   !      where on the curve the liquid fills the pores the ice leaves. A layer that is   !
   !      nearly all ice fills at the curve's driest start, SE_MIN. ----------------------!
   elemental function soil_ice_state(c, theta_ice) result(ice)
      type(soil_layer_curve_t), intent(in) :: c
      real(wp),                 intent(in) :: theta_ice
      type(soil_ice_t) :: ice
      real(wp) :: psi, dpsi, w, dw, dk
      if (theta_ice <= 0.0_wp) return
      ice%theta_ice = theta_ice
      ice%impedance = 10.0_wp ** (-ICE_IMPEDANCE_EXP * theta_ice / c%theta_sat)
      ice%tau_full  = liquid_theta_to_unknown(c, c%theta_sat - theta_ice)
      call liquid_water_relations(c, ice%tau_full, psi, dpsi, w, dw, ice%k_full, dk)
      ice%head_full  = psi / grav_head
      ice%dhead_full = dpsi / grav_head
   end function soil_ice_state

   !----- The unknown tau of a layer holding theta [m3/m3] of water, liquid and ice, the      !
   !      inverse of soil_water_relations's water. -------------------------------------------!
   elemental function soil_theta_to_unknown(c, ice, theta) result(tau)
      type(soil_layer_curve_t), intent(in) :: c
      type(soil_ice_t),         intent(in) :: ice
      real(wp),                 intent(in) :: theta
      real(wp) :: tau
      if (ice%theta_ice > 0.0_wp .and. theta > c%theta_sat) then
         tau = ice%tau_full + ((theta - c%theta_sat) / SOIL_SPECIFIC_STORAGE - ice%head_full) / ice%dhead_full
      else
         tau = liquid_theta_to_unknown(c, theta - ice%theta_ice)
      end if
   end function soil_theta_to_unknown

   !----- The unknown of a layer holding theta [m3/m3] of liquid and no ice. A layer at   !
   !      or below residual starts just inside the curve, at SE_MIN; a Campbell layer at  !
   !      saturation sits at the branch point, its air-entry potential. ------------------!
   elemental function liquid_theta_to_unknown(c, theta) result(tau)
      type(soil_layer_curve_t), intent(in) :: c
      real(wp),                 intent(in) :: theta
      real(wp) :: tau
      real(wp) :: se, head, dse_dhead
      se = (theta - c%theta_res) / max(c%theta_sat - c%theta_res, tiny_num)
      if (se <= c%se_switch) then
         tau = max(se, SE_MIN)     ! clamp-ok: a state at residual starts just inside the curve
         return
      end if
      if (theta > c%theta_sat) then
         head = (theta - c%theta_sat) / SOIL_SPECIFIC_STORAGE
      else if (se >= 1.0_wp) then
         head = 0.0_wp
      else
         call saturation_to_head(c%retention, c%par_a, c%par_n, se, head, dse_dhead)
      end if
      tau = c%se_switch + c%slope_switch * (head - c%head_switch)
   end function liquid_theta_to_unknown

   !=======================================================================================!
   !     The curve's kernels.                                                              !
   !=======================================================================================!

   !----- Matric head [m] at saturation se in (0, 1), and dSe/dhead [1/m]. The van         !
   !      Genuchten head is evaluated through its logarithm, so a very dry layer gives a   !
   !      large but finite number, and its slope Se (1 - Se^(1/m)) / x is written without  !
   !      a cancellation. -----------------------------------------------------------------!
   elemental subroutine saturation_to_head(retention, par_a, par_n, se, head, dse_dhead)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_a, par_n, se
      real(wp),    intent(out) :: head, dse_dhead
      real(wp) :: m, u, logx, x
      if (retention == SOIL_RETENTION_CAMPBELL) then
         head      = par_a * se ** (-par_n)
         dse_dhead = -se ** (1.0_wp + par_n) / (par_n * par_a)
      else
         m    = 1.0_wp - 1.0_wp / par_n
         u    = 1.0_wp - se ** (1.0_wp / m)
         logx = (-log(se) / m + log(u)) / par_n
         x    = exp(min(logx, LOG_HEAD_MAX))   ! clamp-ok: keeps an extreme dry head finite
         head = -x / par_a
         dse_dhead = par_a * m * par_n * se * u / x
      end if
   end subroutine saturation_to_head

   !----- Saturation and dSe/dhead [1/m] at a head below zero, written in x = alpha |head| !
   !      so that 1 - Se^(1/m) = x^n / (1 + x^n) stays exact near saturation. A Campbell   !
   !      layer between its air-entry point and zero is saturated. ------------------------!
   elemental subroutine head_to_saturation(retention, par_a, par_n, head, se, dse_dhead)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_a, par_n, head
      real(wp),    intent(out) :: se, dse_dhead
      real(wp) :: m, x, xn, u
      if (retention == SOIL_RETENTION_CAMPBELL) then
         if (head < par_a) then
            se        = (head / par_a) ** (-1.0_wp / par_n)
            dse_dhead = -se ** (1.0_wp + par_n) / (par_n * par_a)
         else
            se = 1.0_wp ; dse_dhead = 0.0_wp
         end if
         return
      end if
      m  = 1.0_wp - 1.0_wp / par_n
      x  = par_a * abs(head)
      xn = x ** par_n
      se = (1.0_wp + xn) ** (-m)
      u  = xn / (1.0_wp + xn)                      ! 1 - Se^(1/m)
      dse_dhead = par_a * m * par_n * se * u / max(x, tiny_num)
   end subroutine head_to_saturation

   !----- Conductivity [m/s] and dK/dSe at saturation se: van Genuchten-Mualem with L =    !
   !      VG_PORE_CONNECTIVITY (docs/science/soil_biophysics.md eq. (2)), or Campbell's Ks !
   !      Se^(2b+3). No floor: K goes to zero with Se. The van Genuchten slope is          !
   !      unbounded at Se = 1; the solver calls this only dry of the branch point. The     !
   !      slope's powers come from the ones K takes, a division each: Se^(L-1) = Se^L/Se,  !
   !      Se^(1/m-1) = Se^(1/m)/Se and u^(m-1) = u^m/u, with u = 1 - Se^(1/m). ------------!
   elemental subroutine saturation_to_conductivity(retention, par_n, ksat, se, kcond, dk_dse)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_n, ksat, se
      real(wp),    intent(out) :: kcond, dk_dse
      real(wp) :: m, se_1m, se_l, u, u_m, t
      if (retention == SOIL_RETENTION_CAMPBELL) then
         kcond  = ksat * se ** (2.0_wp * par_n + 3.0_wp)
         dk_dse = (2.0_wp * par_n + 3.0_wp) * kcond / se
      else
         m     = 1.0_wp - 1.0_wp / par_n
         se_1m = se ** (1.0_wp / m)
         se_l  = se ** VG_PORE_CONNECTIVITY
         u     = 1.0_wp - se_1m
         u_m   = u ** m
         t     = 1.0_wp - u_m
         kcond  = ksat * se_l * t * t
         dk_dse = ksat * (VG_PORE_CONNECTIVITY * (se_l / se) * t * t + 2.0_wp * se_l * t * (u_m / u) * (se_1m / se))
      end if
   end subroutine saturation_to_conductivity

   !----- Conductivity [m/s] and dK/dhead [1/s] at a head below zero, in x = alpha |head|  !
   !      like head_to_saturation. The slope is taken in the head directly: near           !
   !      saturation dK/dSe grows without bound while dSe/dhead goes to zero. -------------!
   elemental subroutine head_to_conductivity(retention, par_a, par_n, ksat, head, kcond, dk_dhead)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_a, par_n, ksat, head
      real(wp),    intent(out) :: kcond, dk_dhead
      real(wp) :: m, x, xn, u, t, se, dse_dhead, dk_dse
      if (retention == SOIL_RETENTION_CAMPBELL) then
         call head_to_saturation(retention, par_a, par_n, head, se, dse_dhead)
         if (head < par_a) then
            call saturation_to_conductivity(retention, par_n, ksat, se, kcond, dk_dse)
            dk_dhead = dk_dse * dse_dhead
         else
            kcond = ksat ; dk_dhead = 0.0_wp
         end if
         return
      end if
      m  = 1.0_wp - 1.0_wp / par_n
      x  = par_a * abs(head)
      xn = x ** par_n
      se = (1.0_wp + xn) ** (-m)
      u  = xn / (1.0_wp + xn)                      ! 1 - Se^(1/m)
      dse_dhead = par_a * m * par_n * se * u / max(x, tiny_num)
      t     = 1.0_wp - u ** m
      kcond = ksat * se ** VG_PORE_CONNECTIVITY * t * t
      !----- dK/dhead = dK/dSe dSe/dhead, the Mualem term rewritten in u and x. -----------!
      dk_dhead = ksat * (VG_PORE_CONNECTIVITY * t * t * se ** (VG_PORE_CONNECTIVITY - 1.0_wp) * dse_dhead          &
                 + 2.0_wp * par_a * m * par_n * t * se ** (VG_PORE_CONNECTIVITY + 1.0_wp / m) * u ** m             &
                   / max(x, tiny_num))
   end subroutine head_to_conductivity

   !=======================================================================================!
   !     The pond.                                                                         !
   !=======================================================================================!

   !----- The pond at head h [m]: potential [MPa], water [kg/m2] and their derivatives     !
   !      (docs/science/soil_plant_water.md §2). At or below zero head it holds nothing    !
   !      and its potential follows the soil it drains into. ------------------------------!
   elemental subroutine pond_water_relations(h, psi, dpsi, w, dw)
      real(wp), intent(in)  :: h
      real(wp), intent(out) :: psi, dpsi, w, dw
      psi  = grav_head * h
      dpsi = grav_head
      if (h <= 0.0_wp) then
         w = 0.0_wp ; dw = 0.0_wp
      else if (h < POND_DEPRESSION_DEPTH) then
         w  = rho_h2o * h * h / (2.0_wp * POND_DEPRESSION_DEPTH)
         dw = rho_h2o * h / POND_DEPRESSION_DEPTH
      else
         w  = rho_h2o * (h - 0.5_wp * POND_DEPRESSION_DEPTH)
         dw = rho_h2o
      end if
   end subroutine pond_water_relations

   !----- The head [m] of a pond holding w [kg/m2]. An empty pond has no head to name; it  !
   !      returns zero, and the solver finds the head from the soil below. ----------------!
   elemental function pond_water_to_head(w) result(h)
      real(wp), intent(in) :: w
      real(wp) :: h, depth
      depth = w / rho_h2o
      if (depth <= 0.0_wp) then
         h = 0.0_wp
      else if (depth < 0.5_wp * POND_DEPRESSION_DEPTH) then
         h = sqrt(2.0_wp * POND_DEPRESSION_DEPTH * depth)
      else
         h = depth + 0.5_wp * POND_DEPRESSION_DEPTH
      end if
   end function pond_water_to_head

end module meds_soil_hydr_lib
