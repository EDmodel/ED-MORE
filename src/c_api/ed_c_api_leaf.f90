! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_c_api_leaf -- the C-API shim for LEAF GAS EXCHANGE (`edmore.plant.leaf`).                !
!                                                                                          !
! One shim per subsystem, mirroring the Fortran tree (src/README.md), so a change to one       !
! subsystem's ABI rebuilds and re-tests only that shim.                                        !
!                                                                                          !
! EVERY bind(c) struct here is an ABI CONTRACT with `python/edmore/plant/_ffi.py`: the field   !
! ORDER must match, member for member. It is compiled by a mandatory ctest target              !
! (`test_c_api_leaf`) precisely so a mismatch is a BUILD failure, not a silent wrong answer: a !
! component inserted mid-type in `leaf_photo_pft_t` breaks the C API while the rest of the    !
! suite stays green, and only a target that compiles this file sees it.                        !
!==========================================================================================!
module ed_c_api_leaf
   use iso_c_binding,          only : c_double, c_int
   use ed_kinds,               only : wp, ik
   use ed_leaf_gas_exchange,   only : leaf_env_t, leaf_out_t
   use ed_pft_params,          only : leaf_photo_pft_t
   use ed_leaf_gas_exchange,   only : solve_leaf_gas_exchange
   use ed_leaf_gas_exchange,   only : assimilation_demand_c3, electron_transport_j, rubisco_kinetics
   use ed_therm_lib,           only : peaked_arrhenius_scale, arrhenius_scale
   implicit none
   private

   public :: leaf_env_c, leaf_args_c, leaf_flux_c, leaf_c3_demand_c, leaf_kinetics_c
   public :: ed_leaf_solve, ed_leaf_solve_batch, ed_assimilation_demand_c3, ed_electron_transport_j
   public :: ed_rubisco_kinetics
   public :: ed_peaked_arrhenius, ed_arrhenius

   !----- C-interoperable mirror of leaf_env_t (8 doubles). --------------------------------!
   type, bind(c) :: leaf_env_c
      real(c_double) :: par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi
   end type leaf_env_c

   !----- C-interoperable mirror of leaf_out_t (7 doubles + 2 ints; converged 0/1). --------!
   type, bind(c) :: leaf_flux_c
      real(c_double) :: A_net, A_gross, gs, ci, cs, transpiration, rd
      integer(c_int) :: limitation
      integer(c_int) :: converged
   end type leaf_flux_c

   !----- C-interoperable mirror of leaf_photo_pft_t (1 int + 37 doubles, same order), less    !
   !      psi_onset and vpd_min (see to_params). ---------------------------------------------!
   type, bind(c) :: leaf_args_c
      integer(c_int) :: pathway
      real(c_double) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(c_double) :: g0, g1, d0, quantum_yield, theta_j
      real(c_double) :: theta_cj_c3, theta_ip_c3, theta_cj_c4, theta_ic_c4
      real(c_double) :: lambda25, psi_open, psi_close, lambda_psi_exp, sref_stomata
      real(c_double) :: kc25, ko25, gstar25
      real(c_double) :: ea_kc, ea_ko, ea_gstar, ea_vcmax, ea_jmax, ea_rd
      real(c_double) :: hd_vcmax, hd_jmax, hd_rd, ds_vcmax, ds_jmax, ds_rd
      real(c_double) :: o2_mol_frac, absorptance, phi_psii
   end type leaf_args_c

   !----- C-interoperable C3 demand rates at a prescribed Ci (4 doubles). -------------------!
   type, bind(c) :: leaf_c3_demand_c
      real(c_double) :: A_gross, Ac, Aj, Ap
   end type leaf_c3_demand_c

   !----- C-interoperable C3 Rubisco kinetics as mole fractions [umol/mol] (4 doubles). -----!
   type, bind(c) :: leaf_kinetics_c
      real(c_double) :: kc, ko, gstar, o2
   end type leaf_kinetics_c

contains

   !---------------------------------------------------------------------------------------!
   ! Coupled leaf gas-exchange: unpack the C structs into the model types, solve, pack back. !
   ! sm, temp_response_form and colim are the SM_*, TEMP_RESPONSE_* and COLIM_* integer      !
   ! codes; use_boundary_layer is 0/1.                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine ed_leaf_solve(env_c, p_c, sm, temp_response_form, colim, use_boundary_layer, flux_c) bind(c, name="ed_leaf_solve")
      type(leaf_env_c),    intent(in)  :: env_c
      type(leaf_args_c),   intent(in)  :: p_c
      integer(c_int), value, intent(in) :: sm, temp_response_form, colim, use_boundary_layer
      type(leaf_flux_c),   intent(out) :: flux_c
      type(leaf_env_t)          :: env
      type(leaf_photo_pft_t)    :: p
      type(leaf_out_t)         :: flux

      env = to_env(env_c)
      p = to_params(p_c)

      call solve_leaf_gas_exchange(env, p, int(sm), int(temp_response_form), int(colim), use_boundary_layer /= 0_c_int, flux)

      flux_c%A_net = flux%A_net
      flux_c%A_gross = flux%A_gross
      flux_c%gs = flux%gs
      flux_c%ci = flux%ci
      flux_c%cs = flux%cs
      flux_c%transpiration = flux%transpiration
      flux_c%rd = flux%rd
      flux_c%limitation = int(flux%limitation, c_int)
      flux_c%converged = merge(1_c_int, 0_c_int, flux%converged)
   end subroutine ed_leaf_solve

   !---------------------------------------------------------------------------------------!
   ! The same solve for n leaves sharing one parameter set: one call instead of n, for a caller !
   ! that sweeps drivers (a light response, a canopy). env_c and flux_c are arrays of the       !
   ! mirror structs; each leaf is exactly ed_leaf_solve's.                                       !
   !---------------------------------------------------------------------------------------!
   subroutine ed_leaf_solve_batch(n, env_c, p_c, sm, temp_response_form, colim, use_boundary_layer, flux_c) &
                                    bind(c, name="ed_leaf_solve_batch")
      integer(c_int), value, intent(in)  :: n
      type(leaf_env_c),      intent(in)  :: env_c(n)
      type(leaf_args_c),     intent(in)  :: p_c
      integer(c_int), value, intent(in)  :: sm, temp_response_form, colim, use_boundary_layer
      type(leaf_flux_c),     intent(out) :: flux_c(n)
      integer(c_int) :: i
      do i = 1_c_int, n
         call ed_leaf_solve(env_c(i), p_c, sm, temp_response_form, colim, use_boundary_layer, flux_c(i))
      end do
   end subroutine ed_leaf_solve_batch

   !---------------------------------------------------------------------------------------!
   ! Raw C3 FvCB demand at a PRESCRIBED intercellular CO2 (assimilation_demand_c3), stomata bypassed !
   ! and NO temperature scaling: the caller passes already-in-situ values -- vcmax, j (the      !
   ! electron-transport RATE, from ed_electron_transport_j), tpu, and the mole-fraction          !
   ! kinetics gstar/kc/ko/o2 [umol/mol]. Returns gross A and the Ac/Aj/Ap limitation rates      !
   ! (net = gross - Rd is the caller's business). colim is a COLIM_* code; theta_cj and theta_ip  !
   ! are the TWO C3 co-limitation curvatures -- Ac/Aj and then (Ac,Aj)/Ap. Neither is            !
   ! theta_j, the electron-transport hyperbola's curvature (docs/science/leaf_gas_exchange.md §2, !
   ! eq. (2)). Use COLIM_MIN for a sharp min().                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine ed_assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colim, theta_cj,      &
                                   theta_ip, dem_c) bind(c, name="ed_assimilation_demand_c3")
      real(c_double), value, intent(in) :: ci, vcmax, j, tpu, gstar, kc, ko, o2, theta_cj, theta_ip
      integer(c_int), value, intent(in) :: colim
      type(leaf_c3_demand_c), intent(out) :: dem_c
      real(wp) :: a_gross, ac, aj, ap
      call assimilation_demand_c3(real(ci, wp), real(vcmax, wp), real(j, wp), real(tpu, wp),           &
                           real(gstar, wp), real(kc, wp), real(ko, wp), real(o2, wp),           &
                           int(colim), real(theta_cj, wp), real(theta_ip, wp), a_gross, ac, aj, ap)
      dem_c%A_gross = a_gross
      dem_c%Ac = ac
      dem_c%Aj = aj
      dem_c%Ap = ap
   end subroutine ed_assimilation_demand_c3

   !---------------------------------------------------------------------------------------!
   ! The C3 Rubisco kinetics Kc, Ko, Gamma* and O2 as mole fractions [umol/mol] at a leaf    !
   ! temperature [K] and pressure [Pa] (rubisco_kinetics; docs/science/leaf_gas_exchange.md  !
   ! eq. (3)): the solver's own values, for callers of ed_assimilation_demand_c3.            !
   !---------------------------------------------------------------------------------------!
   subroutine ed_rubisco_kinetics(params_c, t_leaf, pressure, kin_c) bind(c, name="ed_rubisco_kinetics")
      type(leaf_args_c),     intent(in)  :: params_c
      real(c_double), value, intent(in)  :: t_leaf, pressure
      type(leaf_kinetics_c), intent(out) :: kin_c
      type(leaf_photo_pft_t)   :: p
      real(wp) :: kc, ko, gstar, o2
      p = to_params(params_c)
      call rubisco_kinetics(p, real(t_leaf, wp), real(pressure, wp), kc, ko, gstar, o2)
      kin_c%kc = kc
      kin_c%ko = ko
      kin_c%gstar = gstar
      kin_c%o2 = o2
   end subroutine ed_rubisco_kinetics

   !---------------------------------------------------------------------------------------!
   ! Electron-transport rate J from Jmax and incident PAR (non-rectangular hyperbola).      !
   !---------------------------------------------------------------------------------------!
   function ed_electron_transport_j(par, absorptance, phi_psii, jmax, theta) result(j)       &
                                      bind(c, name="ed_electron_transport_j")
      real(c_double), value, intent(in) :: par, absorptance, phi_psii, jmax, theta
      real(c_double)                    :: j
      j = electron_transport_j(par, absorptance, phi_psii, jmax, theta)
   end function ed_electron_transport_j

   !---------------------------------------------------------------------------------------!
   ! Peaked-Arrhenius temperature response (k25 anchored; Ea/Hd/dS in J/mol, T in K).       !
   !---------------------------------------------------------------------------------------!
   function ed_peaked_arrhenius(k25, ea, hd, ds, t_leaf) result(y) bind(c, name="ed_peaked_arrhenius")
      real(c_double), value, intent(in) :: k25, ea, hd, ds, t_leaf
      real(c_double)                    :: y
      y = peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf)
   end function ed_peaked_arrhenius

   !---------------------------------------------------------------------------------------!
   ! Plain Arrhenius temperature response.                                                 !
   !---------------------------------------------------------------------------------------!
   function ed_arrhenius(k25, ea, t_leaf) result(y) bind(c, name="ed_arrhenius")
      real(c_double), value, intent(in) :: k25, ea, t_leaf
      real(c_double)                    :: y
      y = arrhenius_scale(k25, ea, t_leaf)
   end function ed_arrhenius

   !---------------------------------------------------------------------------------------!
   ! Unpack the C mirror structs into the model's derived types (used by ed_leaf_solve).     !
   !---------------------------------------------------------------------------------------!
   pure function to_env(env_c) result(env)
      type(leaf_env_c), intent(in) :: env_c
      type(leaf_env_t)             :: env
      !----- KEYWORD form, deliberately. A positional structure constructor silently re-binds every !
      !      value when a component is inserted mid-type, and components WITH defaults cannot be    !
      !      skipped positionally at all; leaf_photo_pft_t has such components (see to_params      !
      !      below). -------------------------------------------------------------------------------!
      env = leaf_env_t(par=env_c%par, leaf_temp=env_c%leaf_temp, vpd=env_c%vpd, ca=env_c%ca,      &
                       pressure=env_c%pressure, psi_leaf=env_c%psi_leaf, gb=env_c%gb,             &
                       psi=env_c%psi)
   end function to_env

   pure function to_params(p_c) result(p)
      type(leaf_args_c),   intent(in) :: p_c
      type(leaf_photo_pft_t)          :: p
      !----- KEYWORD form (see to_env). A positional constructor here would shift every trailing value !
      !      by one each time leaf_photo_pft_t gains a component mid-type, and could not omit         !
      !      psi_onset and vpd_min, which carry defaults.                                              !
      !                                                                                          !
      !      psi_onset and vpd_min are absent from leaf_args_c and take their type defaults            !
      !      (0 MPa, 0.05 kPa; edmore.canopy builds its parameters from the run's config instead).     !
      !      Adding a field to a bind(c) struct breaks the ABI for every existing caller, so they stay !
      !      out. psi_open/psi_close do cross: the caller sets the ramp of the fast capacity limb      !
      !      (docs/science/leaf_gas_exchange.md §4.1, eq. (9)), which is always applied; sref_stomata  !
      !      with the default psi_onset sets the slow stomatal factor (§4.2, eq. (10)). ---------------!
      p = leaf_photo_pft_t(pathway=p_c%pathway, vcmax25=p_c%vcmax25, jmax25=p_c%jmax25,          &
             tpu25=p_c%tpu25, rd25=p_c%rd25, kp25=p_c%kp25, g0=p_c%g0, g1=p_c%g1, d0=p_c%d0,      &
             quantum_yield=p_c%quantum_yield, theta_j=p_c%theta_j,                                &
             theta_cj_c3=p_c%theta_cj_c3, theta_ip_c3=p_c%theta_ip_c3,                            &
             theta_cj_c4=p_c%theta_cj_c4, theta_ic_c4=p_c%theta_ic_c4,                            &
             lambda25=p_c%lambda25, psi_open=p_c%psi_open,                                        &
             psi_close=p_c%psi_close, lambda_psi_exp=p_c%lambda_psi_exp,                          &
             sref_stomata=p_c%sref_stomata,                                                       &
             kc25=p_c%kc25, ko25=p_c%ko25, gstar25=p_c%gstar25, ea_kc=p_c%ea_kc,                  &
             ea_ko=p_c%ea_ko, ea_gstar=p_c%ea_gstar, ea_vcmax=p_c%ea_vcmax,                       &
             ea_jmax=p_c%ea_jmax, ea_rd=p_c%ea_rd, hd_vcmax=p_c%hd_vcmax, hd_jmax=p_c%hd_jmax,    &
             hd_rd=p_c%hd_rd, ds_vcmax=p_c%ds_vcmax, ds_jmax=p_c%ds_jmax, ds_rd=p_c%ds_rd,        &
             o2_mol_frac=p_c%o2_mol_frac, absorptance=p_c%absorptance, phi_psii=p_c%phi_psii)
   end function to_params

end module ed_c_api_leaf
