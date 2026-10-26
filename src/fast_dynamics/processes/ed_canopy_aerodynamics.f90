! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_canopy_aerodynamics -- the canopy's turbulent exchange                                !
! (docs/science/canopy_aerodynamics.md), a function of the canopy-air and vegetation state !
! with no state of its own: the CLM5 Monin-Obukhov surface layer and the canopy-air <->    !
! atmosphere conductances for the canopy air's enthalpy, humidity and CO2 (§2), the ground !
! <-> canopy-air conductance (§3), and the in-canopy wind and the Nusselt leaf and wood    !
! boundary layers (§4). Main routines: canopy_aerodynamics, mo_surface_layer,              !
! cas_atm_conductances, canopy_air_conductances, canopy_roughness. All pure; per-cohort    !
! outputs go into caller-owned arrays (CLAUDE.md, compiler traps).                         !
!==========================================================================================!
module ed_canopy_aerodynamics
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : grav, pi, tiny_num, von_karman, safe_exp, virtual_temp_factor
   use ed_biophysics_params, only : aero_params_t
   use ed_therm_lib,          only : cas_molar_density
   use ed_numerics,           only : equal_height_layers
   implicit none
   private

   !----- The sweeps stop when the canopy air's coupling to the atmosphere, u* times the scalar    !
   !      profile factor, moves by less than this share of itself from one sweep to the next. -----!
   real(wp), parameter :: MO_TOL = 1.0e-4_wp       !< [-]
   !----- Bounds of the Monin-Obukhov stability, from CLM5 (Zeng et al. 1998, J. Climate 11,       !
   !      2628-2644): the first guess uses Ri_b <= 0.19, where 1 - 5 Ri_b stays positive, and      !
   !      starts at least 0.01 from neutral; the guess and every sweep stay at least -100. --------!
   real(wp), parameter :: RIB_GUESS_MAX      = 0.19_wp      !< [-]
   real(wp), parameter :: ZETA_GUESS_NEUTRAL = 0.01_wp      !< [-]
   real(wp), parameter :: ZETA_UNSTABLE_MIN  = -100.0_wp    !< [-]
   !----- Each Monin-Obukhov sweep's convective velocity (CLM5): w_c = beta w_*, with w_* the      !
   !      convective velocity scale of a mixed layer MIXED_LAYER_DEPTH deep, from the sweep's own  !
   !      buoyancy flux (docs/science/canopy_aerodynamics.md eq. (4)). ----------------------------!
   real(wp), parameter :: CONVECTIVE_BETA   = 1.0_wp      !< [-]
   real(wp), parameter :: MIXED_LAYER_DEPTH = 1000.0_wp   !< [m]
   !----- Cap on the under-canopy stability parameter stab of CLM4's dense-canopy ground           !
   !      conductance (whose cs_dense and gamma_g it uses): at 10, its factor 1/(1+gamma_g*stab)   !
   !      is 1/6. ---------------------------------------------------------------------------------!
   real(wp), parameter :: GROUND_STAB_MAX = 10.0_wp               !< [-]
   !----- The in-canopy wind's extinction per unit leaf area of a layer of closed crowns, through  !
   !      the whole layer and to its mid-depth (ED2's canopy_struct_dynamics, crown area 1). ------!
   real(wp), parameter :: WIND_EXTINCTION_LAYER = 0.50_wp   !< [-]
   real(wp), parameter :: WIND_EXTINCTION_MID   = 0.25_wp   !< [-]
   !----- The least ratio of the height above the displacement to the roughness length in the wind !
   !      profile, so its logarithm stays positive. -----------------------------------------------!
   real(wp), parameter :: WIND_HEIGHT_ROUGH_MIN = 1.001_wp  !< [-]
   !----- ED2's cap on the canopy's snow burial fraction (canopy_struct_dynamics), so a buried     !
   !      canopy keeps a trace of its own roughness. ----------------------------------------------!
   real(wp), parameter :: BURIAL_MAX = 0.99_wp              !< [-]
   !----- The bare soil's transfer coefficient from the ground to its roughness height, a (z0 u*/nu)^b  !
   !      (Zeng & Dickinson 1998, J. Climate 11, 2628-2644; CLM5's soil coefficient). ---------------!
   real(wp), parameter :: SOIL_TRANSFER_A   = 0.13_wp       !< [-]
   real(wp), parameter :: SOIL_TRANSFER_EXP = 0.45_wp       !< [-]

   public :: canopy_aerodynamics             !< the per-patch entry point: winds, conductances, boundary layers
   public :: canopy_roughness                !< roughness length + displacement height of a patch
   public :: mo_surface_layer                !< the Monin-Obukhov surface-layer solve (CLM5)
   public :: cas_atm_conductances            !< CAS<->atmosphere bulk conductances from ustar + profile factors
   public :: canopy_air_conductances         !< the same, the surface layer re-solved at a canopy-air state
   public :: open_canopy_fraction            !< a patch's open fraction W from its leaf and wood area
   public :: canopy_burial_fraction          !< the share of a patch's canopy height a snow pack buries
   public :: aero_env_t, aero_args_t, aero_out_t, ensure_aero_out_capacity
   public :: set_aero_env_atm, set_aero_env_canopy
   public :: alloc_aero_out                  ! used only by the tests

   !----- One patch's reference-level forcing and its canopy-air and ground state. ----------------!
   type :: aero_env_t
      real(wp) :: u_ref     = 2.0_wp                !< [m/s]      wind at reference height
      real(wp) :: zref      = 30.0_wp               !< [m]        reference height: the canopy-air top, where the forcing is sampled
      real(wp) :: theta_atm = 298.15_wp             !< [K]        potential temp, referenced to zref (= the air temperature there)
      real(wp) :: shv_atm   = 0.010_wp              !< [kg/kg]    specific humidity at zref
      real(wp) :: can_temp  = 298.15_wp             !< [K]        CAS temperature, which at zref is also its potential temperature
      real(wp) :: can_shv   = 0.010_wp              !< [kg/kg]    CAS specific humidity
      real(wp) :: t_ground  = 298.15_wp             !< [K]        ground skin temp (ground-conductance stability)
   end type aero_env_t

   !----- Per-patch canopy geometry. -------------------------------------------------------!
   type :: aero_args_t
      real(wp) :: veg_height   = 20.0_wp            !< [m]  canopy top height
      real(wp) :: opencan_frac = 0.0_wp             !< [-]  open fraction W = exp(-(LAI+WAI)): 1 bare, 0 dense
      real(wp) :: snowfac      = 0.0_wp             !< [-]  snow burial fraction of the canopy
   end type aero_args_t

   !----- Outputs: per-patch scalars + per-cohort arrays (caller-owned; written in place). ---!
   type :: aero_out_t
      integer(ik) :: n_coh = 0_ik                   !< number of cohorts
      real(wp) :: ustar = 0.0_wp, tstar = 0.0_wp   !< [m/s],[K] friction velocity, temperature scale
      real(wp) :: heat_transfer_factor   = 0.0_wp   !< [-] scalar profile factor for heat (CLM's temp1): g_atm_heat = rho ustar this
      real(wp) :: vapour_transfer_factor = 0.0_wp   !< [-] the same for vapour and CO2 (CLM's temp2)
      real(wp) :: zeta = 0.0_wp, rib = 0.0_wp, obu = 0.0_wp   !< [-],[-],[m] stability, bulk Richardson number, Obukhov length
      real(wp) :: ggbare = 0.0_wp, ggveg = 0.0_wp, ggnet = 0.0_wp   !< [m/s] ground conductances (r_aero = 1/ggnet)
      real(wp) :: rough = 0.0_wp, displace = 0.0_wp   !< [m]
      real(wp) :: uh = 0.0_wp                        !< [m/s] canopy-top wind (diagnostic)
      real(wp), allocatable :: wind(:)               !< [m/s] per-cohort in-canopy wind
      real(wp), allocatable :: leaf_gbh(:), leaf_gbw(:)   !< [m/s] leaf boundary-layer heat/vapour conductance
      real(wp), allocatable :: wood_gbh(:), wood_gbw(:)   !< [m/s] wood boundary-layer heat/vapour conductance
   end type aero_out_t

contains

   !---------------------------------------------------------------------------------------!
   ! The canopy-air <-> atmosphere bulk conductances (docs/science/canopy_aerodynamics.md  !
   ! eq. (6)), for enthalpy and vapour [kg/m2/s] and for CO2 [mol/m2/s].                   !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_atm_conductances(rho, can_dmol, ustar, temp_heat, temp_vapour, g_atm_heat, g_atm_vapour, g_atm_co2)
      real(wp), intent(in)  :: rho          !< [kg/m3]  canopy-air density
      real(wp), intent(in)  :: can_dmol     !< [mol/m3] canopy-air dry-air molar density
      real(wp), intent(in)  :: ustar        !< [m/s]    friction velocity
      real(wp), intent(in)  :: temp_heat    !< [-]      scalar profile factor for heat (heat_transfer_factor)
      real(wp), intent(in)  :: temp_vapour  !< [-]      scalar profile factor for vapour / CO2 (vapour_transfer_factor)
      real(wp), intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2   !< [kg/m2/s], [kg/m2/s], [mol/m2/s]
      g_atm_heat = rho      * ustar * temp_heat
      g_atm_vapour = rho      * ustar * temp_vapour
      g_atm_co2 = can_dmol * ustar * temp_vapour
   end subroutine cas_atm_conductances

   !----- A patch's open fraction W = exp(-(LAI + WAI)) [-], the weight of bare ground in the      !
   !      ground conductance (docs/science/canopy_aerodynamics.md eq. (7)). -----------------------!
   elemental real(wp) function open_canopy_fraction(lai, wai) result(w)
      real(wp), intent(in) :: lai, wai   !< [m2/m2] the patch's total leaf and wood area indices
      w = safe_exp(-(lai + wai))
   end function open_canopy_fraction

   !----- The canopy's snow burial fraction: the pack's depth over the canopy height, capped at    !
   !      BURIAL_MAX (ED2's snowfac_can; docs/science/canopy_aerodynamics.md eq. (1)). ------------!
   elemental real(wp) function canopy_burial_fraction(snow_depth, veg_height) result(f)
      real(wp), intent(in) :: snow_depth   !< [m] the pack's depth
      real(wp), intent(in) :: veg_height   !< [m] the canopy height, positive
      f = min(BURIAL_MAX, snow_depth / veg_height)
   end function canopy_burial_fraction

   !----- Roughness length and displacement height of a patch from its canopy height               !
   !      (docs/science/canopy_aerodynamics.md eq. (1)). ------------------------------------------!
   pure subroutine canopy_roughness(cfg, geom, rough, displace)
      type(aero_params_t), intent(in)  :: cfg        !< aerodynamic parameters
      type(aero_args_t),   intent(in)  :: geom       !< canopy height and snow burial
      real(wp),            intent(out) :: rough      !< [m] roughness length (z0m)
      real(wp),            intent(out) :: displace   !< [m] zero-plane displacement height
      rough    = (1.0_wp - geom%snowfac) * cfg%z0m_ratio * geom%veg_height                     &
                 + geom%snowfac * cfg%snow_rough
      rough    = max(rough, cfg%snow_rough)
      displace = cfg%d_ratio * geom%veg_height * (1.0_wp - geom%snowfac)
   end subroutine canopy_roughness

   !---------------------------------------------------------------------------------------!
   ! Fill an aero_out_t for one patch (docs/science/canopy_aerodynamics.md). The cohort    !
   ! arrays are in any order; order(k) is the k-th cohort from the ground up. Only the     !
   ! wind cascade and the boundary layers depend on the order, so the cohorts are put      !
   ! bottom to top, solved, and their outputs put back. The reordering is an identity for  !
   ! n <= 1.                                                                               !
   !---------------------------------------------------------------------------------------!
   pure subroutine canopy_aerodynamics(cfg, env, geom, n, order, height, lai, leaf_temp, wood_temp, &
                                       leaf_width, branch_diam, out)
      type(aero_params_t), intent(in)    :: cfg               !< aerodynamic parameters
      type(aero_env_t),    intent(in)    :: env               !< reference-level forcing, canopy-air and ground state
      type(aero_args_t),   intent(in)    :: geom              !< canopy height, open fraction, snow burial
      integer(ik),         intent(in)    :: n                 !< number of cohorts
      integer(ik),         intent(in)    :: order(:)          !< the cohorts from the ground up
      real(wp),            intent(in)    :: height(:)         !< [m]  cohort height
      real(wp),            intent(in)    :: lai(:)            !< [m2/m2] leaf area index
      real(wp),            intent(in)    :: leaf_temp(:)      !< [K]  leaf temperature
      !----- Wood has its own temperature, which lags the leaves', for its free convection. -------!
      real(wp),            intent(in)    :: wood_temp(:)      !< [K]  wood temperature
      real(wp),            intent(in)    :: leaf_width(:)     !< [m]  leaf width (boundary layer)
      real(wp),            intent(in)    :: branch_diam(:)    !< [m]  branch diameter (wood boundary layer)
      type(aero_out_t),    intent(inout) :: out               !< the patch's winds, conductances and boundary layers
      integer(ik) :: k, i
      real(wp)    :: h_bt(n), lai_bt(n), lt_bt(n), wt_bt(n), lw_bt(n), bd_bt(n)
      real(wp)    :: wind_bt(n), lgbh_bt(n), lgbw_bt(n), wgbh_bt(n), wgbw_bt(n)

      do k = 1_ik, n
         i = order(k)
         h_bt(k)   = height(i)
         lai_bt(k) = lai(i)
         lt_bt(k)  = leaf_temp(i)
         wt_bt(k)  = wood_temp(i)
         lw_bt(k)  = leaf_width(i)
         bd_bt(k)  = branch_diam(i)
      end do

      call aerodynamics_bottom_to_top(cfg, env, geom, n, h_bt, lai_bt, lt_bt, wt_bt, lw_bt, bd_bt, out)

      !----- out%*(k) runs bottom to top; copy out, then put back in the callers' order. ----------!
      do k = 1_ik, n
         wind_bt(k) = out%wind(k)
         lgbh_bt(k) = out%leaf_gbh(k)
         lgbw_bt(k) = out%leaf_gbw(k)
         wgbh_bt(k) = out%wood_gbh(k)
         wgbw_bt(k) = out%wood_gbw(k)
      end do
      do k = 1_ik, n
         out%wind(order(k))     = wind_bt(k)
         out%leaf_gbh(order(k)) = lgbh_bt(k)
         out%leaf_gbw(order(k)) = lgbw_bt(k)
         out%wood_gbh(order(k)) = wgbh_bt(k)
         out%wood_gbw(order(k)) = wgbw_bt(k)
      end do
   end subroutine canopy_aerodynamics

   !---------------------------------------------------------------------------------------!
   ! The patch's aerodynamics with its cohorts ordered bottom (1) to top (n), in the order !
   ! of docs/science/canopy_aerodynamics.md §1-§4: roughness, the surface layer, the       !
   ! ground conductance, the in-canopy wind and the leaf and wood boundary layers.         !
   !---------------------------------------------------------------------------------------!
   pure subroutine aerodynamics_bottom_to_top(cfg, env, geom, n, height, lai, leaf_temp, wood_temp,  &
                                              leaf_width, branch_diam, out)
      type(aero_params_t), intent(in)    :: cfg     !< aerodynamic parameters
      type(aero_env_t),    intent(in)    :: env     !< reference-level forcing, canopy-air and ground state
      type(aero_args_t),   intent(in)    :: geom    !< canopy height, open fraction, snow burial
      integer(ik),         intent(in)    :: n       !< number of cohorts
      real(wp),            intent(in)    :: height(n), lai(n)   !< [m], [m2/m2]
      real(wp),            intent(in)    :: leaf_temp(n), wood_temp(n), leaf_width(n), branch_diam(n)   !< [K], [K], [m], [m]
      type(aero_out_t),    intent(inout) :: out     !< the patch's winds, conductances and boundary layers

      real(wp)    :: rough, displace, zldis, uh, stab, denom, lai_layer
      integer(ik) :: ico, ico0, il, n_layer, layer(n)

      !----- 1. Roughness + displacement from canopy height (snow blend). ------------------!
      call canopy_roughness(cfg, geom, rough, displace)
      zldis    = max(env%zref - displace, 2.0_wp * rough)
      out%rough = rough
      out%displace = displace
      out%n_coh = n

      !----- 2. CLM Monin-Obukhov surface layer (canopy z0h = z0m). ------------------------!
      call mo_surface_layer(cfg, env%u_ref, env%zref, displace, rough, env%theta_atm,          &
                            env%shv_atm, env%can_temp, env%can_shv, out%ustar, out%heat_transfer_factor,          &
                            out%zeta, out%rib, out%obu)
      out%vapour_transfer_factor = out%heat_transfer_factor   ! z0q = z0h => scalar factors equal
      out%tstar = out%heat_transfer_factor * (env%theta_atm - env%can_temp)

      !----- 3. Ground <-> canopy-air conductance: the bare-ground value and CLM4's dense-canopy   !
      !      one, blended by the open fraction W. The under-canopy stability stab = g h (T_cas -   !
      !      T_ground) / (T_mean ustar^2), kept in [0, GROUND_STAB_MAX], lowers the dense-canopy   !
      !      value when the ground is colder than the canopy air. ---------------------------------!
      out%ggbare = soil_conductance(cfg, out%ustar, env%can_temp)                 ! [m/s] bare-ground scalar conductance
      denom      = (env%can_temp + env%t_ground) * out%ustar ** 2
      stab       = 0.0_wp
      if (denom > tiny_num) then
         stab = max(0.0_wp, min(GROUND_STAB_MAX, 2.0_wp * grav * geom%veg_height               &
                    * (env%can_temp - env%t_ground) / denom))
      end if
      out%ggveg = cfg%cs_dense * out%ustar / (1.0_wp + cfg%gamma_g * stab)
      out%ggnet = geom%opencan_frac * out%ggbare + (1.0_wp - geom%opencan_frac) * out%ggveg

      !----- 4. In-canopy wind, top -> bottom through each layer's closed crowns, a layer being   !
      !      the cohorts of one height (equal_height_layers). ------------------------------------!
      uh     = reduced_wind(cfg, out%ustar, out%zeta, geom%veg_height, displace, rough, zldis)
      out%uh = uh
      call equal_height_layers(height, n, layer, n_layer)
      ico = n
      do il = n_layer, 1_ik, -1_ik
         ico0 = ico
         lai_layer = lai(ico)
         do while (ico0 > 1_ik)
            if (layer(ico0 - 1_ik) /= il) exit
            ico0 = ico0 - 1_ik
            lai_layer = lai_layer + lai(ico0)
         end do
         out%wind(ico0:ico) = max(cfg%ugbmin, uh * safe_exp(-WIND_EXTINCTION_MID * lai_layer))   ! wind at crown mid-depth
         uh                 = uh * safe_exp(-WIND_EXTINCTION_LAYER * lai_layer)                    ! attenuate for the layer below
         ico = ico0 - 1_ik
      end do

      !----- 5. Per-cohort leaf (flat plate) + wood (cylinder) boundary layers [m/s]. ------!
      do ico = 1_ik, n
         out%leaf_gbh(ico) = boundary_layer_conductance(cfg, out%wind(ico), leaf_temp(ico), env%can_temp, &
                                              leaf_width(ico), .true.)
         out%leaf_gbw(ico) = cfg%gbh_2_gbw * out%leaf_gbh(ico)
         out%wood_gbh(ico) = boundary_layer_conductance(cfg, out%wind(ico), wood_temp(ico), env%can_temp, &
                                              branch_diam(ico), .false.)
         out%wood_gbw(ico) = cfg%gbh_2_gbw * out%wood_gbh(ico)
      end do
   end subroutine aerodynamics_bottom_to_top

   !---------------------------------------------------------------------------------------!
   ! CLM5 Monin-Obukhov surface layer (docs/science/canopy_aerodynamics.md eqs. (2)-(5)):  !
   ! sweeps of the four-range denominators from the bulk-Richardson guess until the        !
   ! coupling settles to MO_TOL, at most n_iter_mo of them.                                !
   ! Returns ustar, the scalar profile factor heat_transfer_factor = vonk/D_h (CLM's       !
   ! temp1) and the stability diagnostics. Every sweep's zeta is held in                   !
   ! [-100, zeta_max_stable], as CLM5 holds it. zeta divides only in the very-unstable     !
   ! and very-stable branches, so the solve is safe under -fpe0.                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine mo_surface_layer(cfg, u_ref, zref, displace, rough, theta_atm, shv_atm,     &
                                    can_theta, can_shv, ustar, heat_transfer_factor, zeta, rib, obu)
      type(aero_params_t), intent(in)  :: cfg   !< aerodynamic parameters
      real(wp),            intent(in)  :: u_ref, zref, displace, rough, theta_atm, shv_atm   !< [m/s],[m],[m],[m],[K],[kg/kg]
      real(wp),            intent(in)  :: can_theta, can_shv   !< [K], [kg/kg] canopy-air temperature and humidity
      real(wp),            intent(out) :: ustar, heat_transfer_factor, zeta, rib, obu   !< [m/s], [-], [-], [-], [m]
      real(wp)    :: uref, zldis, z0m, thv_atm, thv_can, dthv, dth, dqh, um, thvstar, wc, coupling_prev
      integer(ik) :: it

      uref    = max(u_ref, cfg%ubmin)
      z0m     = rough
      zldis   = max(zref - displace, 2.0_wp * z0m)
      thv_atm = theta_atm * (1.0_wp + virtual_temp_factor * shv_atm)
      thv_can = can_theta * (1.0_wp + virtual_temp_factor * can_shv)
      dthv    = thv_atm - thv_can
      dth     = theta_atm - can_theta
      dqh     = shv_atm   - can_shv

      !----- MoninObukIni: a first convective velocity if unstable, then the analytic Rib -> zeta  !
      !      guess, bounded. ----------------------------------------------------------------------!
      if (dthv >= 0.0_wp) then
         um = uref                                    ! already at least cfg%ubmin
      else
         um = sqrt(uref * uref + cfg%wc * cfg%wc)
      end if
      rib = grav * zldis * dthv / (thv_atm * um * um)
      if (rib >= 0.0_wp) then
         zeta = rib * log(zldis / z0m) / (1.0_wp - 5.0_wp * min(rib, RIB_GUESS_MAX))
         zeta = min(cfg%zeta_max_stable, max(zeta, ZETA_GUESS_NEUTRAL))
      else
         zeta = rib * log(zldis / z0m)
         zeta = max(ZETA_UNSTABLE_MIN, min(zeta, -ZETA_GUESS_NEUTRAL))
      end if

      !----- Sweep until the coupling u* x heat_transfer_factor settles to MO_TOL of itself, at    !
      !      most n_iter_mo sweeps. Under hot canopy air and a weak wind the sweeps swing about    !
      !      the root, and a few of them can miss it by a fifth; near neutral two suffice. --------!
      ustar = cfg%ustmin
      heat_transfer_factor = von_karman / log(zldis / z0m)
      coupling_prev = -1.0_wp
      do it = 1_ik, cfg%n_iter_mo
         ustar   = max(cfg%ustmin, von_karman * um / d_mom(zeta, zldis, z0m, cfg))
         heat_transfer_factor   = von_karman / d_heat(zeta, zldis, z0m, cfg)
         if (abs(ustar * heat_transfer_factor - coupling_prev) <= MO_TOL * ustar * heat_transfer_factor) exit
         coupling_prev = ustar * heat_transfer_factor
         thvstar = (heat_transfer_factor * dth) * (1.0_wp + virtual_temp_factor * shv_atm)                         &
                 + virtual_temp_factor * theta_atm * (heat_transfer_factor * dqh)
         zeta    = zldis * von_karman * grav * thvstar / (ustar * ustar * thv_atm)
         !----- CLM5's bounds on each sweep. Unbounded, a stable night runs past zeta = 1, where   !
         !      u* falls to ustmin and the canopy air loses its coupling to the air above,         !
         !      whatever the wind. CLM5 also holds |zeta| >= 0.01; that bound is left out, since   !
         !      nothing here divides by zeta near neutral and it would put a jump in u* there. ----!
         zeta    = min(cfg%zeta_max_stable, max(zeta, ZETA_UNSTABLE_MIN))
         !----- The convective velocity follows this sweep's buoyancy flux, so it fades to zero as !
         !      the air turns neutral and the wind has no jump there. -----------------------------!
         if (zeta >= 0.0_wp) then
            um = uref
         else
            wc = CONVECTIVE_BETA * (-grav * ustar * thvstar * MIXED_LAYER_DEPTH / thv_atm) ** (1.0_wp / 3.0_wp)
            um = sqrt(uref * uref + wc * wc)
         end if
      end do
      !----- ustar and heat_transfer_factor at the last sweep's zeta. ------------------------------!
      ustar = max(cfg%ustmin, von_karman * um / d_mom(zeta, zldis, z0m, cfg))
      heat_transfer_factor = von_karman / d_heat(zeta, zldis, z0m, cfg)
      obu   = zldis / sign(max(abs(zeta), 1.0e-6_wp), zeta)   ! clamp-ok: a finite Obukhov length near neutral
   end subroutine mo_surface_layer

   !----- Momentum profile denominator D_m (ustar = vonk*um/D_m), four CLM5 ranges. ---------------!
   pure function d_mom(zeta, zldis, z0m, cfg) result(dm)
      real(wp),            intent(in) :: zeta, zldis, z0m   !< [-], [m], [m] stability, height above displacement, roughness
      type(aero_params_t), intent(in) :: cfg                !< aerodynamic parameters
      real(wp) :: dm, zeta0m                                !< [-] D_m (the result), z0m / L
      zeta0m = z0m * zeta / zldis                                  ! = z0m / L  (no division by zeta)
      if (zeta < -cfg%zeta_m) then                                 ! very unstable
         dm = log(-cfg%zeta_m * zldis / (zeta * z0m)) - stabfunc1(-cfg%zeta_m) + stabfunc1(zeta0m) &
              + 1.14_wp * ((-zeta) ** 0.333_wp - cfg%zeta_m ** 0.333_wp)
      else if (zeta < 0.0_wp) then                                 ! unstable
         dm = log(zldis / z0m) - stabfunc1(zeta) + stabfunc1(zeta0m)
      else if (zeta <= 1.0_wp) then                                ! stable (Businger-Dyer)
         dm = log(zldis / z0m) + 5.0_wp * zeta - 5.0_wp * zeta0m
      else                                                         ! very stable (log cap)
         dm = log(zldis / (zeta * z0m)) + 5.0_wp - 5.0_wp * zeta0m + (5.0_wp * log(zeta) + zeta - 1.0_wp)
      end if
   end function d_mom

   !----- Scalar profile denominator D_h (heat_transfer_factor = vonk/D_h), four CLM5 ranges. -----!
   pure function d_heat(zeta, zldis, z0h, cfg) result(dh)
      real(wp),            intent(in) :: zeta, zldis, z0h   !< [-], [m], [m] stability, height above displacement, roughness
      type(aero_params_t), intent(in) :: cfg                !< aerodynamic parameters
      real(wp) :: dh, zeta0h                                !< [-] D_h (the result), z0h / L
      zeta0h = z0h * zeta / zldis
      if (zeta < -cfg%zeta_t) then                                 ! very unstable
         dh = log(-cfg%zeta_t * zldis / (zeta * z0h)) - stabfunc2(-cfg%zeta_t) + stabfunc2(zeta0h) &
              + 0.8_wp * (cfg%zeta_t ** (-0.333_wp) - (-zeta) ** (-0.333_wp))
      else if (zeta < 0.0_wp) then                                 ! unstable
         dh = log(zldis / z0h) - stabfunc2(zeta) + stabfunc2(zeta0h)
      else if (zeta <= 1.0_wp) then                                ! stable
         dh = log(zldis / z0h) + 5.0_wp * zeta - 5.0_wp * zeta0h
      else                                                         ! very stable
         dh = log(zldis / (zeta * z0h)) + 5.0_wp - 5.0_wp * zeta0h + (5.0_wp * log(zeta) + zeta - 1.0_wp)
      end if
   end function d_heat

   !----- CLM StabilityFunc1 (psi_m of eq. (5), unstable branch: call only with zeta < 0) ---!
   pure function stabfunc1(zeta) result(psi)
      real(wp), intent(in) :: zeta     !< [-] stability, negative
      real(wp) :: psi, x               !< [-] psi_m (the result), (1 - 16 zeta)^(1/4)
      x   = (1.0_wp - 16.0_wp * zeta) ** 0.25_wp
      psi = 2.0_wp * log((1.0_wp + x) * 0.5_wp) + log((1.0_wp + x * x) * 0.5_wp)                &
            - 2.0_wp * atan(x) + pi * 0.5_wp
   end function stabfunc1

   !----- CLM StabilityFunc2 (psi_h of eq. (5), unstable branch: call only with zeta < 0) ---!
   pure function stabfunc2(zeta) result(psi)
      real(wp), intent(in) :: zeta     !< [-] stability, negative
      real(wp) :: psi, x2              !< [-] psi_h (the result), (1 - 16 zeta)^(1/2)
      x2  = sqrt(1.0_wp - 16.0_wp * zeta)
      psi = 2.0_wp * log((1.0_wp + x2) * 0.5_wp)
   end function stabfunc2

   !----- The log-profile wind at `height`, floored at ugbmin: the canopy-top wind                 !
   !      (docs/science/canopy_aerodynamics.md eq. (8)). It is u*/kappa times the momentum         !
   !      denominator d_mom from the roughness up to height - displace, at the surface layer's     !
   !      Obukhov length, so the profile and u* share CLM5's four ranges. height - displace is     !
   !      held just above the roughness length, so the log stays positive. ------------------------!
   pure function reduced_wind(cfg, ustar, zeta, height, displace, rough, zldis) result(u)
      type(aero_params_t), intent(in) :: cfg   !< aerodynamic parameters
      real(wp),            intent(in) :: ustar, zeta, height, displace, rough, zldis   !< [m/s], [-], [m], [m], [m], [m]
      real(wp) :: u, hd                        !< [m/s] the wind (the result); [m] height - displace
      hd = max(height - displace, rough * WIND_HEIGHT_ROUGH_MIN)
      u  = (ustar / von_karman) * d_mom(zeta * hd / zldis, hd, rough, cfg)
      u  = max(u, cfg%ugbmin)
   end function reduced_wind

   !---------------------------------------------------------------------------------------!
   ! Leaf (flat-plate) or wood (cylinder) boundary-layer heat conductance gbh_mos [m/s],   !
   ! the sum of forced (Reynolds) and free (Grashof) convection Nusselt numbers            !
   ! (docs/science/canopy_aerodynamics.md eqs. (10)-(11)).                                 !
   !---------------------------------------------------------------------------------------!
   pure function boundary_layer_conductance(cfg, wind, elem_temp, can_temp, lchar, is_leaf) result(gbh_mos)
      type(aero_params_t), intent(in) :: cfg     !< aerodynamic parameters
      real(wp),            intent(in) :: wind, elem_temp, can_temp, lchar   !< [m/s], [K], [K], [m]
      logical,             intent(in) :: is_leaf !< true: a leaf (flat plate); false: wood (cylinder)
      real(wp) :: gbh_mos, nu_v, al_th, re, gr, nu_forced, nu_free   !< [m/s] the result; [m2/s] nu_v, al_th
      nu_v  = air_kinematic_viscosity(cfg, can_temp)
      al_th = cfg%th_diff0  * (1.0_wp + cfg%dth_diff  * (can_temp - cfg%t_ref_air))
      re    = max(wind * lchar / nu_v, tiny_num)
      gr    = max(grav / (can_temp * nu_v * nu_v) * abs(elem_temp - can_temp) * lchar ** 3, 0.0_wp)
      if (is_leaf) then                                            ! flat plate
         nu_forced = max(cfg%aflat_lami * re ** cfg%nflat_lami, cfg%aflat_turb * re ** cfg%nflat_turb)
         nu_free   = max(cfg%bflat_lami * gr ** cfg%mflat_lami, cfg%bflat_turb * gr ** cfg%mflat_turb)
      else                                                         ! cylinder
         nu_forced = max(cfg%ocyli_lami + cfg%acyli_lami * re ** cfg%ncyli_lami,                &
                         cfg%ocyli_turb + cfg%acyli_turb * re ** cfg%ncyli_turb)
         nu_free   = max(cfg%bcyli_lami * gr ** cfg%mcyli_lami, cfg%bcyli_turb * gr ** cfg%mcyli_turb)
      end if
      gbh_mos = max(cfg%gbhmos_min, (nu_forced + nu_free) * al_th / lchar)
   end function boundary_layer_conductance

   !----- The air's kinematic viscosity [m2/s], linear in temperature about t_ref_air. -----------!
   elemental real(wp) function air_kinematic_viscosity(cfg, temp) result(nu)
      type(aero_params_t), intent(in) :: cfg    !< aerodynamic parameters
      real(wp),            intent(in) :: temp   !< [K] the air's temperature
      nu = cfg%kin_visc0 * (1.0_wp + cfg%dkin_visc * (temp - cfg%t_ref_air))
   end function air_kinematic_viscosity

   !----- The bare soil's conductance to the canopy air [m/s] (docs/science/canopy_aerodynamics.md !
   !      eq. (7)): u* times CLM5's soil coefficient, the thin layer between the soil and its      !
   !      roughness height only, kappa / (a (z0 u*/nu)^b). The turbulent layer above it is the     !
   !      canopy air's own exchange with the atmosphere, so the two in series are CLM5's bare-soil !
   !      resistance from the soil to zref. -------------------------------------------------------!
   pure real(wp) function soil_conductance(cfg, ustar, air_temp) result(g)
      type(aero_params_t), intent(in) :: cfg        !< aerodynamic parameters
      real(wp),            intent(in) :: ustar      !< [m/s] friction velocity
      real(wp),            intent(in) :: air_temp   !< [K] canopy-air temperature, for the viscosity
      real(wp) :: roughness_re   !< [-] the soil's roughness Reynolds number
      roughness_re = cfg%soil_rough * ustar / air_kinematic_viscosity(cfg, air_temp)
      g = von_karman * ustar / (SOIL_TRANSFER_A * roughness_re ** SOIL_TRANSFER_EXP)
   end function soil_conductance

   !----- Allocate the per-cohort output arrays of an aero_out_t (the kernel writes in place). !
   subroutine alloc_aero_out(out, n_coh)
      type(aero_out_t), intent(out) :: out     !< the record to allocate
      integer(ik),      intent(in)  :: n_coh   !< number of cohorts
      out%n_coh = n_coh
      allocate(out%wind(n_coh), out%leaf_gbh(n_coh), out%leaf_gbw(n_coh),                      &
               out%wood_gbh(n_coh), out%wood_gbw(n_coh))
      out%wind = 0.0_wp
      out%leaf_gbh = 0.0_wp
      out%leaf_gbw = 0.0_wp
      out%wood_gbh = 0.0_wp
      out%wood_gbw = 0.0_wp
   end subroutine alloc_aero_out

   !----- Grow-only capacity of aero_out_t: canopy_aerodynamics overwrites every active element. --!
   subroutine ensure_aero_out_capacity(out, n_coh)
      type(aero_out_t), intent(inout) :: out     !< the record to grow
      integer(ik),       intent(in)    :: n_coh  !< number of cohorts
      if (.not. allocated(out%wind)) then
         call alloc_aero_out(out, n_coh)
      else if (size(out%wind) < n_coh) then
         call alloc_aero_out(out, n_coh)
      else
         out%n_coh = n_coh
      end if
   end subroutine ensure_aero_out_capacity

   !---------------------------------------------------------------------------------------!
   ! Fill aero_env_t's reference level from the forcing at zref, the canopy-air top, where !
   ! potential temperature equals the actual air temperature                               !
   ! (docs/science/canopy_aerodynamics.md §2, "The reference level"). Set all four fields  !
   ! here: a theta_atm left at its 298.15 K default gives the Monin-Obukhov solve a false  !
   ! stratification.                                                                       !
   !---------------------------------------------------------------------------------------!
   pure subroutine set_aero_env_atm(aenv, u_ref, zref, air_temp, shv_atm)
      type(aero_env_t), intent(inout) :: aenv       !< the environment to fill
      real(wp),         intent(in)    :: u_ref      !< [m/s]   wind at zref
      real(wp),         intent(in)    :: zref       !< [m]     the reference height, the canopy-air top
      real(wp),         intent(in)    :: air_temp   !< [K]     actual air temperature at zref
      real(wp),         intent(in)    :: shv_atm    !< [kg/kg] specific humidity at zref
      aenv%u_ref     = u_ref
      aenv%zref      = zref
      aenv%theta_atm = air_temp
      aenv%shv_atm   = shv_atm
   end subroutine set_aero_env_atm

   !----- Fill aero_env_t's canopy-air and ground fields from the model state; a stale canopy-air  !
   !      temperature biases the stability solve as a stale theta_atm does. -----------------------!
   pure subroutine set_aero_env_canopy(aenv, can_temp, can_shv, t_ground)
      type(aero_env_t), intent(inout) :: aenv   !< the environment to fill
      real(wp),         intent(in)    :: can_temp, can_shv, t_ground   !< [K], [kg/kg], [K]
      aenv%can_temp  = can_temp
      aenv%can_shv   = can_shv
      aenv%t_ground  = t_ground
   end subroutine set_aero_env_canopy

   !----- The canopy-air <-> atmosphere conductances at a given canopy-air temperature and         !
   !      humidity (docs/science/canopy_aerodynamics.md eq. (6)): the surface layer solved again   !
   !      against that state, with the reference level in env and the canopy geometry in geom, so  !
   !      the conductances follow the canopy air's own stability. Vapour and CO2 take the heat     !
   !      factor (z0q = z0h); CO2 uses the dry-air molar density at that humidity. ----------------!
   pure subroutine canopy_air_conductances(cfg, env, geom, rho, can_temp, can_shv, g_atm_heat,       &
                                           g_atm_vapour, g_atm_co2)
      type(aero_params_t), intent(in)  :: cfg        !< aerodynamic parameters
      type(aero_env_t),    intent(in)  :: env        !< the reference level (its canopy-air fields are not read)
      type(aero_args_t),   intent(in)  :: geom       !< canopy height and snow burial
      real(wp),            intent(in)  :: rho        !< [kg/m3]  canopy-air density
      real(wp),            intent(in)  :: can_temp   !< [K]      canopy-air temperature
      real(wp),            intent(in)  :: can_shv    !< [kg/kg]  canopy-air specific humidity
      real(wp),            intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2   !< [kg/m2/s], [kg/m2/s], [mol/m2/s]
      real(wp) :: rough, displace, ustar, heat_factor, zeta, rib, obu, can_dmol
      call canopy_roughness(cfg, geom, rough, displace)
      call mo_surface_layer(cfg, env%u_ref, env%zref, displace, rough, env%theta_atm, env%shv_atm,      &
                            can_temp, can_shv, ustar, heat_factor, zeta, rib, obu)
      can_dmol = cas_molar_density(rho, can_shv)
      call cas_atm_conductances(rho, can_dmol, ustar, heat_factor, heat_factor, g_atm_heat, g_atm_vapour,  &
                                g_atm_co2)
   end subroutine canopy_air_conductances

end module ed_canopy_aerodynamics
