! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_fast_integ_rk45 -- the explicit fast-loop scheme a run can choose ("rk45"):
! rk45_step is one Cash-Karp 5(4) step over the whole march state (canopy air, soil energy
! and water, plant water, films) through the shared right-hand side march_derivs, from the
! same frozen record as ESDIRK2, with its embedded error and b-weighted boundary-flux
! amounts; rk45_state_railed checks the march's end. It is also the explicit reference the
! tests compare ESDIRK2 against. The stepper marches it, and marches the step again with
! ESDIRK2 when it cannot resolve it (docs/science/numerical_scheme.md §3 `rk45`, §5b).
!---~---
module ed_fast_integ_rk45
   use ed_kinds,              only : wp, ik
   use ed_constants,          only : cp_air
   use ed_therm_lib,          only : cas_enthalpy_to_temp, internal_energy_liquid
   use ed_fast_integ_derivs, only : march_derivs, cas_conductances
   use ed_site_state_vars,  only : site_surface_t
   use ed_fast_integ_types,       only : march_state_t, frozen_t, march_tend_t, surface_tend_t, stage_bflux_t, &
                                     boundary_flux_t, state_init, state_accum, state_sub, clamp_cas,               &
                                     diagnose_soil_temps, at_temp_bound
   implicit none
   private

   public :: rk45_step, rk45_state_railed
   public :: RK45_P_ORDER, RK45_WORK_CAP

   !---~---
   !   Cash-Karp 5(4) tableau (Cash & Karp 1990, ACM TOMS 16:201; Numerical Recipes' rkck). The
   ! frozen record has no explicit time dependence, so the stage times never enter.
   !---~---
   real(wp), parameter :: A21 = 0.2_wp
   real(wp), parameter :: A31 = 3.0_wp/40.0_wp,      A32 = 9.0_wp/40.0_wp
   real(wp), parameter :: A41 = 0.3_wp,               A42 = -0.9_wp,             A43 = 1.2_wp
   real(wp), parameter :: A51 = -11.0_wp/54.0_wp,     A52 = 2.5_wp,              A53 = -70.0_wp/27.0_wp, &
                          A54 = 35.0_wp/27.0_wp
   real(wp), parameter :: A61 = 1631.0_wp/55296.0_wp, A62 = 175.0_wp/512.0_wp,   A63 = 575.0_wp/13824.0_wp, &
                          A64 = 44275.0_wp/110592.0_wp, A65 = 253.0_wp/4096.0_wp
   !--- 5th-order (committed) and embedded 4th-order weights; b2 = 0 in both, b5 = 0 at 5th order.
   real(wp), parameter :: B1 = 37.0_wp/378.0_wp,    B3 = 250.0_wp/621.0_wp, B4 = 125.0_wp/594.0_wp, &
                          B6 = 512.0_wp/1771.0_wp
   real(wp), parameter :: BS1 = 2825.0_wp/27648.0_wp, BS3 = 18575.0_wp/48384.0_wp,                &
                          BS4 = 13525.0_wp/55296.0_wp, BS5 = 277.0_wp/14336.0_wp, BS6 = 0.25_wp
   !--- The embedded estimate's order, p in docs/science/numerical_scheme.md eq. (4).
   integer(ik), parameter :: RK45_P_ORDER = 4_ik

   !---~---
   !   Work cap per dt_fast. A healthy step takes a few sub-steps (docs/science/numerical_scheme.md
   ! §3 `rk45`); on a dense, cold canopy the explicit march instead thrashes and still rails, and
   ! the stepper marches the step again with ESDIRK2, so stopping at the cap bounds the wasted
   ! work.
   !---~---
   integer(ik), parameter :: RK45_WORK_CAP = 64_ik

contains

   !---~---
   !   explicit_stage_bflux -- one explicit stage's boundary-flux rates, the counterpart of
   ! backward_euler_stage's: the state is the stage input y_stage, and k and surf_tend are
   ! march_derivs' tendency and surface diagnostic there, condensation included at the
   ! stage's own canopy-air temperature. The canopy air is charged at the conductances the
   ! tendency used, re-solved at the stage's state (docs/science/numerical_scheme.md §5b);
   ! cas_conductances is pure, so this matches march_derivs, and carrying them on
   ! surface_tend_t instead would perturb ifx's inlining and FMA choices in surface_derivs.
   ! The soil's and the column's bottom face is k%g_bottom, positive up.
   !---~---
   pure subroutine explicit_stage_bflux(y_stage, frozen, k, surf_tend, bf)
      type(march_state_t),  intent(in)  :: y_stage
      type(frozen_t), intent(in)  :: frozen
      type(march_tend_t),   intent(in)  :: k
      type(surface_tend_t),  intent(in)  :: surf_tend
      type(stage_bflux_t),   intent(out) :: bf
      real(wp) :: gah, gaw, gac, e_infil, e_drain
      call cas_conductances(frozen%cas, y_stage%cas_enthalpy, y_stage%cas_shv, gah, gaw, gac)
      e_infil = frozen%hydrology%infiltration * internal_energy_liquid(frozen%hydrology%t_infil)
      e_drain = k%drainage_rate * internal_energy_liquid(frozen%hydrology%t_bot)
      bf%cas_enth_in   = surf_tend%src_enth + gah*frozen%cas%box%atm_enthalpy
      bf%cas_enth_out  = gah*y_stage%cas_enthalpy
      bf%cas_vap_in    = surf_tend%src_vap  + gaw*frozen%cas%box%atm_specific_humidity
      bf%cas_vap_out   = gaw*y_stage%cas_shv
      bf%cas_co2_in    = frozen%cas%nee_biotic + gac*frozen%cas%box%atm_co2
      bf%cas_co2_out   = gac*y_stage%cas_co2
      bf%soil_enth_in  = k%g_top + k%g_bottom + e_infil
      bf%soil_enth_out = sum(frozen%roots%root_heat) + frozen%hydrology%exfil_enth + e_drain
      !---~---
      !   Soil water's storage change and flows are the network's, booked once per step by
      ! close_ledgers: the per-stage soil-water amounts are zero, and the whole-column water
      ! here is only the canopy-air vapour exchange.
      !---~---
      bf%soil_wat_in   = 0.0_wp
      bf%soil_wat_out  = 0.0_wp
      bf%whole_enth_in = surf_tend%coh_rnet + surf_tend%ground_rad + k%g_bottom
      bf%whole_enth_out= gah*(y_stage%cas_enthalpy - frozen%cas%box%atm_enthalpy) + e_drain
      bf%whole_wat_in  = 0.0_wp
      bf%whole_wat_out = gaw*(y_stage%cas_shv - frozen%cas%box%atm_specific_humidity)
      bf%whole_cond    = surf_tend%cond
      bf%whole_cond_enth = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
      bf%snow_enth_in  = surf_tend%snow_enth_in
      bf%snow_vap_out  = surf_tend%snow_vap_out
      bf%atm_heat_out  = gah*cp_air*(cas_enthalpy_to_temp(y_stage%cas_enthalpy, y_stage%cas_shv)          &
                                     - frozen%cas%env%theta_atm)
      bf%atm_vap_out   = gaw*(y_stage%cas_shv - frozen%cas%box%atm_specific_humidity)
      bf%rad_net       = surf_tend%coh_rnet + surf_tend%ground_rad
   end subroutine explicit_stage_bflux

   !---~---
   !   rk45_step -- one Cash-Karp 5(4) step over march_derivs, committing the 5th-order
   ! solution; y_err is the 5th - 4th difference, and bf the step's boundary-flux amounts,
   ! b-weighted like the commit, with the tissue-temperature time integrals and each cohort's
   ! transpiration. Each stage input's canopy air is clamped (clamp_cas) so that an oversized
   ! trial gives a finite derivative, which the embedded error then rejects. Soil energy needs
   ! no clamp: any finite energy has a temperature, and a stage reads it only through terms
   ! linear in it. The committed state is not clamped: at the floor a step is committed with
   ! the controller out of moves, and a clamp there would move mass or energy with no ledger
   ! entry; a railed commit is caught by rk45_state_railed.
   !---~---
   pure subroutine rk45_step(y, frozen, site, n, nsl, dt, y_out, y_err, bf, clamp_stage_n)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(march_state_t),  intent(out) :: y_out, y_err
      type(boundary_flux_t), optional, intent(out) :: bf   ! b-weighted boundary-flux amounts over dt (ledger)
      !--- Stage-clamp count (the report's WORK_CLAMP_STAGE), expected on an oversized trial.
      integer(ik), optional, intent(inout) :: clamp_stage_n

      type(march_tend_t)  :: k1, k2, k3, k4, k5, k6
      type(march_state_t) :: y_stage, y_4th
      type(surface_tend_t) :: surf_tend
      type(stage_bflux_t)  :: s(6)
      !---~---
      !   Per-stage tissue temperatures and transpiration, for the time integrals the tissues'
      ! energy and the diagnostics are set from.
      !---~---
      real(wp) :: tleaf_s(n,6), twood_s(n,6), transp_s(n,6)

      call march_derivs(y, frozen, site, n, nsl, k1, sf_out=surf_tend)
      call explicit_stage_bflux(y, frozen, k1, surf_tend, s(1))
      tleaf_s(1:n,1) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,1) = surf_tend%wood_temp(1:n)
      transp_s(1:n,1) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A21, k1, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call march_derivs(y_stage, frozen, site, n, nsl, k2, sf_out=surf_tend)
      call explicit_stage_bflux(y_stage, frozen, k2, surf_tend, s(2))
      tleaf_s(1:n,2) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,2) = surf_tend%wood_temp(1:n)
      transp_s(1:n,2) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A31, k1, n, nsl)
      call state_accum(y_stage, dt*A32, k2, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call march_derivs(y_stage, frozen, site, n, nsl, k3, sf_out=surf_tend)
      call explicit_stage_bflux(y_stage, frozen, k3, surf_tend, s(3))
      tleaf_s(1:n,3) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,3) = surf_tend%wood_temp(1:n)
      transp_s(1:n,3) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A41, k1, n, nsl)
      call state_accum(y_stage, dt*A42, k2, n, nsl)
      call state_accum(y_stage, dt*A43, k3, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call march_derivs(y_stage, frozen, site, n, nsl, k4, sf_out=surf_tend)
      call explicit_stage_bflux(y_stage, frozen, k4, surf_tend, s(4))
      tleaf_s(1:n,4) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,4) = surf_tend%wood_temp(1:n)
      transp_s(1:n,4) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A51, k1, n, nsl)
      call state_accum(y_stage, dt*A52, k2, n, nsl)
      call state_accum(y_stage, dt*A53, k3, n, nsl)
      call state_accum(y_stage, dt*A54, k4, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call march_derivs(y_stage, frozen, site, n, nsl, k5, sf_out=surf_tend)
      call explicit_stage_bflux(y_stage, frozen, k5, surf_tend, s(5))
      tleaf_s(1:n,5) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,5) = surf_tend%wood_temp(1:n)
      transp_s(1:n,5) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A61, k1, n, nsl)
      call state_accum(y_stage, dt*A62, k2, n, nsl)
      call state_accum(y_stage, dt*A63, k3, n, nsl)
      call state_accum(y_stage, dt*A64, k4, n, nsl)
      call state_accum(y_stage, dt*A65, k5, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call march_derivs(y_stage, frozen, site, n, nsl, k6, sf_out=surf_tend)
      call explicit_stage_bflux(y_stage, frozen, k6, surf_tend, s(6))
      tleaf_s(1:n,6) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,6) = surf_tend%wood_temp(1:n)
      transp_s(1:n,6) = surf_tend%transp_c(1:n)

      !--- y_out = y + dt*(B1*k1 + B3*k3 + B4*k4 + B6*k6), the 5th-order solution.
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt*B1, k1, n, nsl)
      call state_accum(y_out, dt*B3, k3, n, nsl)
      call state_accum(y_out, dt*B4, k4, n, nsl)
      call state_accum(y_out, dt*B6, k6, n, nsl)
      !---~---
      !   y_4th, the embedded 4th-order solution, is its own temporary: passing y_err as both
      ! state_sub's b and out would alias two dummies, one modified, which Fortran forbids.
      !---~---
      call state_init(y, n, nsl, y_4th)
      call state_accum(y_4th, dt*BS1, k1, n, nsl)
      call state_accum(y_4th, dt*BS3, k3, n, nsl)
      call state_accum(y_4th, dt*BS4, k4, n, nsl)
      call state_accum(y_4th, dt*BS5, k5, n, nsl)
      call state_accum(y_4th, dt*BS6, k6, n, nsl)
      call state_sub(y_out, y_4th, n, nsl, y_err)   ! y_err := y_5th - y_4th (the embedded estimate)

      if (present(bf)) then
         call bflux_cash_karp(s, dt, bf)
         allocate(bf%tissue_leaf_int(n), bf%tissue_wood_int(n), bf%transp_coh(n))
         bf%tissue_leaf_int(1:n) = dt * (B1*tleaf_s(1:n,1) + B3*tleaf_s(1:n,3)                       &
                                         + B4*tleaf_s(1:n,4) + B6*tleaf_s(1:n,6))
         bf%tissue_wood_int(1:n) = dt * (B1*twood_s(1:n,1) + B3*twood_s(1:n,3)                       &
                                         + B4*twood_s(1:n,4) + B6*twood_s(1:n,6))
         bf%transp_coh(1:n)      = dt * (B1*transp_s(1:n,1) + B3*transp_s(1:n,3)                     &
                                         + B4*transp_s(1:n,4) + B6*transp_s(1:n,6))
      end if
   end subroutine rk45_step

   !---~---
   !   The step's boundary-flux amounts over dt: the stage rates b-weighted like the commit
   ! (b2 = b5 = 0), summed over the four stages in order.
   !---~---
   pure subroutine bflux_cash_karp(s, dt, acc)
      type(stage_bflux_t),  intent(in)    :: s(6)
      real(wp),             intent(in)    :: dt
      type(boundary_flux_t), intent(inout) :: acc
      acc%cas_enth_in    = (B1*s(1)%cas_enth_in    + B3*s(3)%cas_enth_in    + B4*s(4)%cas_enth_in              &
                            + B6*s(6)%cas_enth_in) * dt
      acc%cas_enth_out   = (B1*s(1)%cas_enth_out   + B3*s(3)%cas_enth_out   + B4*s(4)%cas_enth_out             &
                            + B6*s(6)%cas_enth_out) * dt
      acc%cas_vap_in     = (B1*s(1)%cas_vap_in     + B3*s(3)%cas_vap_in     + B4*s(4)%cas_vap_in               &
                            + B6*s(6)%cas_vap_in) * dt
      acc%cas_vap_out    = (B1*s(1)%cas_vap_out    + B3*s(3)%cas_vap_out    + B4*s(4)%cas_vap_out              &
                            + B6*s(6)%cas_vap_out) * dt
      acc%cas_co2_in     = (B1*s(1)%cas_co2_in     + B3*s(3)%cas_co2_in     + B4*s(4)%cas_co2_in               &
                            + B6*s(6)%cas_co2_in) * dt
      acc%cas_co2_out    = (B1*s(1)%cas_co2_out    + B3*s(3)%cas_co2_out    + B4*s(4)%cas_co2_out              &
                            + B6*s(6)%cas_co2_out) * dt
      acc%soil_enth_in   = (B1*s(1)%soil_enth_in   + B3*s(3)%soil_enth_in   + B4*s(4)%soil_enth_in             &
                            + B6*s(6)%soil_enth_in) * dt
      acc%soil_enth_out  = (B1*s(1)%soil_enth_out  + B3*s(3)%soil_enth_out  + B4*s(4)%soil_enth_out            &
                            + B6*s(6)%soil_enth_out) * dt
      acc%soil_wat_in    = (B1*s(1)%soil_wat_in    + B3*s(3)%soil_wat_in    + B4*s(4)%soil_wat_in              &
                            + B6*s(6)%soil_wat_in) * dt
      acc%soil_wat_out   = (B1*s(1)%soil_wat_out   + B3*s(3)%soil_wat_out   + B4*s(4)%soil_wat_out             &
                            + B6*s(6)%soil_wat_out) * dt
      acc%whole_enth_in  = (B1*s(1)%whole_enth_in  + B3*s(3)%whole_enth_in  + B4*s(4)%whole_enth_in            &
                            + B6*s(6)%whole_enth_in) * dt
      acc%whole_enth_out = (B1*s(1)%whole_enth_out + B3*s(3)%whole_enth_out + B4*s(4)%whole_enth_out           &
                            + B6*s(6)%whole_enth_out) * dt
      acc%whole_wat_in   = (B1*s(1)%whole_wat_in   + B3*s(3)%whole_wat_in   + B4*s(4)%whole_wat_in             &
                            + B6*s(6)%whole_wat_in) * dt
      acc%whole_wat_out  = (B1*s(1)%whole_wat_out  + B3*s(3)%whole_wat_out  + B4*s(4)%whole_wat_out            &
                            + B6*s(6)%whole_wat_out) * dt
      acc%whole_cond     = (B1*s(1)%whole_cond     + B3*s(3)%whole_cond     + B4*s(4)%whole_cond               &
                            + B6*s(6)%whole_cond) * dt
      acc%whole_cond_enth = (B1*s(1)%whole_cond_enth + B3*s(3)%whole_cond_enth + B4*s(4)%whole_cond_enth       &
                             + B6*s(6)%whole_cond_enth) * dt
      acc%snow_enth_in   = (B1*s(1)%snow_enth_in   + B3*s(3)%snow_enth_in   + B4*s(4)%snow_enth_in             &
                            + B6*s(6)%snow_enth_in) * dt
      acc%snow_vap_out   = (B1*s(1)%snow_vap_out   + B3*s(3)%snow_vap_out   + B4*s(4)%snow_vap_out             &
                            + B6*s(6)%snow_vap_out) * dt
      acc%atm_heat_out   = (B1*s(1)%atm_heat_out   + B3*s(3)%atm_heat_out   + B4*s(4)%atm_heat_out             &
                            + B6*s(6)%atm_heat_out) * dt
      acc%atm_vap_out    = (B1*s(1)%atm_vap_out    + B3*s(3)%atm_vap_out    + B4*s(4)%atm_vap_out              &
                            + B6*s(6)%atm_vap_out) * dt
      acc%rad_net        = (B1*s(1)%rad_net        + B3*s(3)%rad_net        + B4*s(4)%rad_net                  &
                            + B6*s(6)%rad_net) * dt
   end subroutine bflux_cash_karp

   !---~---
   !   Did the explicit march end on a railed state? Its 5th- and 4th-order solutions can rail
   ! together, so the controller sees err ~ 0 and accepts an impossible temperature; the
   ! stepper then marches the step again with ESDIRK2. The canopy air's and each soil layer's
   ! temperature at the march's end, at or beyond the bounds the stage clamp holds.
   !---~---
   pure function rk45_state_railed(y_out, site, nsl) result(railed)
      type(march_state_t), intent(in) :: y_out
      type(site_surface_t), intent(in) :: site
      integer(ik),          intent(in) :: nsl
      logical :: railed
      real(wp) :: t_cas, t_soil(nsl), f_liq(nsl)
      t_cas = cas_enthalpy_to_temp(y_out%cas_enthalpy, y_out%cas_shv)
      call diagnose_soil_temps(y_out, site%soil_profile, nsl, t_soil, f_liq)
      railed = at_temp_bound(t_cas) .or. any(at_temp_bound(t_soil))
   end function rk45_state_railed

end module ed_fast_integ_rk45
