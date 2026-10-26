! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_c_api_leaf -- COVERAGE FOR THE LEAF C-API SHIM.                                     !
!                                                                                          !
! WHY THIS TEST EXISTS -- and why it is MANDATORY, not gated on the pylib option. Outside  !
! this target the shim is compiled only by the optional shared library                     !
! (`-DMEDS_BUILD_PYLIB=ON`), which is in neither the default build nor ctest. So a change  !
! to `leaf_env_t` or `leaf_photo_pft_t` -- a component inserted mid-type shifts every     !
! trailing value of a POSITIONAL structure constructor -- could pass the whole suite on    !
! both back ends while leaving the C API unable to compile. Every C-API shim is therefore  !
! compiled into a ctest target, so an interface mismatch fails the build (docs/building.md !
! §Running the tests).                                                                     !
!                                                                                          !
! The constructors in `to_env` and `to_params` are keyword-form, so the remaining failure  !
! mode is a compile error, and this target makes it impossible to miss. The assertions     !
! below then cover the part a compile cannot: that the bind(c) structs actually map onto   !
! the Fortran types field-for-field.                                                       !
!                                                                                          !
! Deliberately does NOT need the .so, Python, or MEDS_BUILD_PYLIB -- it calls the bind(c) procedures !
! directly as Fortran.                                                                              !
!                                                                                          !
! It calls every entry point of the leaf shim (meds_c_api_leaf) -- the solve and its batch, !
! the scalar temperature responses, the C3 demand, the electron transport and the kinetics -- !
! and of the canopy shim (meds_c_api_canopy) on the reference configuration, which it opens !
! from the source directory: its PFT count, the plastic capacities, the leaf solve over     !
! many leaves and the two-stream over a stand, whose shortwave must close.                  !
!==========================================================================================!
program test_c_api_leaf
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use meds_kinds,        only : wp, ik
   use meds_c_api_leaf,    only : leaf_env_c, leaf_args_c, leaf_flux_c, meds_leaf_solve,          &
                                 meds_arrhenius, meds_peaked_arrhenius, leaf_kinetics_c,          &
                                 meds_rubisco_kinetics, meds_leaf_solve_batch, leaf_c3_demand_c,  &
                                 meds_assimilation_demand_c3, meds_electron_transport_j
   use meds_c_api_canopy,  only : meds_canopy_open, meds_canopy_n_pft, meds_canopy_plastic, meds_canopy_leaf,  &
                                 meds_canopy_two_stream
   use meds_config,        only : meds_config_t
   use meds_config_io,     only : load_meds_config
   use meds_pft_params, only : PATH_C3
   use meds_leaf_params,  only : SM_MEDLYN, COLIM_QUADRATIC, COLIM_MIN
   use meds_therm_lib,     only : TEMP_RESPONSE_PEAKED
   use meds_test_support, only : banner, check, check_close
   implicit none

   type(leaf_env_c)    :: env
   type(leaf_args_c)   :: p
   type(leaf_flux_c)   :: f0, f1, f2, fb
   type(leaf_kinetics_c) :: kin
   real(wp)            :: k25, arr

   call banner('leaf C-API shim')

   p = c3_defaults()

   !=== 1. The shim solves at all -- i.e. to_env/to_params bound every component. ===========!
   env = make_env(0.0_c_double)
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TEMP_RESPONSE_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 0_c_int, f0)
   call check(f0%converged /= 0_c_int, 'C-API solve converges on well-watered C3 defaults')
   call check(f0%A_net > 0.0_wp,  'C-API A_net > 0 in full light')
   call check(f0%gs    > 0.0_wp,  'C-API gs > 0 when unstressed')
   call check(f0%ci > 0.0_wp .and. f0%ci < 400.0_wp, 'C-API ci is between 0 and ca')

   !=== 2. `psi` actually reaches beta_stomata (docs/science/leaf_gas_exchange.md eq. (10)). =!
   !      This is the ABI assertion a compile cannot make: if leaf_env_c's components were         !
   !      mis-bound, `psi` would land on some other field and gs would not respond.                !
   env = make_env(-1.0_c_double)
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TEMP_RESPONSE_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 0_c_int, f1)
   env = make_env(-3.0_c_double)
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TEMP_RESPONSE_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 0_c_int, f2)
   call check(f1%gs < f0%gs, 'C-API gs falls when psi drops 0 -> -1 MPa (beta_stomata is live)')
   call check(f2%gs < f1%gs, 'C-API gs falls further at psi = -3 MPa (monotone in psi)')
   call check(f2%A_net < f0%A_net, 'C-API A_net falls with water stress')

   !=== 3. A scalar entry point, to catch a by-value/by-reference regression. ===============!
   !      meds_arrhenius takes its arguments BY VALUE; a binding slip here shows up as a wrong      !
   !      number rather than a link error, so assert against the closed form at the reference        !
   !      temperature, where the response must be exactly 1x k25. ---------------------------------!
   k25 = 60.0_wp
   arr = real(meds_arrhenius(real(k25, c_double), 65330.0_c_double, 298.15_c_double), wp)
   call check_close(arr, k25, 1.0e-9_wp, 'C-API meds_arrhenius returns k25 at the 25 C reference')

   !=== 4. The boundary-layer switch is an int flag, not a logical -- assert it is honoured. =!
   env = make_env(0.0_c_double) ; env%gb = 0.5_c_double
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TEMP_RESPONSE_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 1_c_int, fb)
   call check(fb%converged /= 0_c_int, 'C-API solve converges with the boundary layer ON')
   call check(fb%cs < 400.0_wp, 'C-API boundary layer draws cs below ca (flag reached the kernel)')

   !=== 5. The kinetics entry reads the parameters it is given: at 25 C and 21 % O2 Gamma* is    !
   !      gstar25 at 100 kPa as a mole fraction (42.75 umol/mol), Kc and Ko are their partial     !
   !      pressures over the air's pressure, and O2 is o2_mol_frac (leaf_gas_exchange.md eq. (3)). !
   call meds_rubisco_kinetics(p, 298.15_c_double, 101325.0_c_double, kin)
   call check_close(real(kin%gstar, wp), 42.75_wp, 1.0e-12_wp, 'C-API Gamma* is 42.75 umol/mol at 25 C')
   call check_close(real(kin%kc, wp), 40.49_wp / 101325.0_wp * 1.0e6_wp, 1.0e-12_wp,            &
                    'C-API Kc is kc25 over the pressure at 25 C')
   call check_close(real(kin%ko, wp), 27840.0_wp / 101325.0_wp * 1.0e6_wp, 1.0e-12_wp,          &
                    'C-API Ko is ko25 over the pressure at 25 C')
   call check_close(real(kin%o2, wp), 209000.0_wp, 1.0e-12_wp, 'C-API O2 is o2_mol_frac in umol/mol')

   !=== 6. The batch entry is the single solve leaf by leaf, bit for bit: the three environments of !
   !      checks 1 and 2 in one call. ------------------------------------------------------------!
   block
      type(leaf_env_c)  :: envs(3)
      type(leaf_flux_c) :: fs(3)
      envs = [make_env(0.0_c_double), make_env(-1.0_c_double), make_env(-3.0_c_double)]
      call meds_leaf_solve_batch(3_c_int, envs, p, int(SM_MEDLYN, c_int), int(TEMP_RESPONSE_PEAKED, c_int),     &
                                 int(COLIM_QUADRATIC, c_int), 0_c_int, fs)
      call check(fs(1)%A_net == f0%A_net .and. fs(2)%A_net == f1%A_net .and. fs(3)%A_net == f2%A_net,         &
                 'C-API batch: A_net is the single solve''s, leaf by leaf')
      call check(fs(1)%gs == f0%gs .and. fs(2)%gs == f1%gs .and. fs(3)%gs == f2%gs .and.                       &
                 fs(1)%ci == f0%ci .and. fs(3)%ci == f2%ci .and. all(fs%converged == 1_c_int),                 &
                 'C-API batch: gs, ci and convergence are the single solve''s')
   end block

   !=== 7. The C3 demand at a prescribed Ci, by value (docs/science/leaf_gas_exchange.md eq. (2)): the    !
   !      hand-computed Ac = 100 x 260 / (300 + 400 (1 + 210000/275000)) = 25.8593, Aj = 120 x 260 /     !
   !      (4 x 300 + 8 x 40) = 20.5263, Ap = 3 x 15 = 45, and the sharp minimum picks Aj. --------------!
   block
      type(leaf_c3_demand_c) :: dem
      call meds_assimilation_demand_c3(300.0_c_double, 100.0_c_double, 120.0_c_double, 15.0_c_double,          &
                                       40.0_c_double, 400.0_c_double, 275000.0_c_double, 210000.0_c_double,   &
                                       int(COLIM_MIN, c_int), 0.98_c_double, 0.95_c_double, dem)
      call check_close(real(dem%Ac, wp), 25.85930_wp, 1.0e-4_wp, 'C-API C3 demand: Ac known-answer')
      call check_close(real(dem%Aj, wp), 20.52632_wp, 1.0e-4_wp, 'C-API C3 demand: Aj known-answer')
      call check_close(real(dem%Ap, wp), 45.0_wp, 1.0e-6_wp, 'C-API C3 demand: Ap known-answer')
      call check_close(real(dem%A_gross, wp), real(dem%Aj, wp), 1.0e-12_wp, 'C-API C3 demand: min picks Aj')
   end block

   !=== 8. The electron transport, by value: none without light; in full light it lies below both the  !
   !      light PSII uses, I2 = phi absorptance PAR / 2, and Jmax, and is the smaller root of the      !
   !      non-rectangular hyperbola theta J^2 - (I2 + Jmax) J + I2 Jmax = 0. --------------------------!
   block
      real(wp), parameter :: PAR = 1500.0_wp, ABSORPTANCE = 0.85_wp, PHI = 0.85_wp, JMAX = 108.0_wp, THETA = 0.85_wp
      real(wp) :: j, i2
      j = real(meds_electron_transport_j(0.0_c_double, real(ABSORPTANCE, c_double), real(PHI, c_double),      &
                                         real(JMAX, c_double), real(THETA, c_double)), wp)
      call check(j == 0.0_wp, 'C-API electron transport: none without light')
      j  = real(meds_electron_transport_j(real(PAR, c_double), real(ABSORPTANCE, c_double), real(PHI, c_double), &
                                          real(JMAX, c_double), real(THETA, c_double)), wp)
      i2 = 0.5_wp * PHI * ABSORPTANCE * PAR
      call check(j > 0.0_wp .and. j < min(i2, JMAX), 'C-API electron transport: below the light and Jmax')
      call check(abs(THETA * j * j - (i2 + JMAX) * j + i2 * JMAX) <= 1.0e-9_wp * i2 * JMAX,                   &
                 'C-API electron transport: the non-rectangular hyperbola''s root')
   end block

   !=== 9. The peaked response, by value: exactly k25 at the 25 C reference. ======================!
   arr = real(meds_peaked_arrhenius(real(k25, c_double), 65330.0_c_double, 200000.0_c_double, 650.0_c_double,   &
                                    298.15_c_double), wp)
   call check_close(arr, k25, 1.0e-9_wp, 'C-API meds_peaked_arrhenius returns k25 at the 25 C reference')

   !=== 10. The canopy shim on the reference configuration, read from the source directory. ==========!
   call test_canopy_shim()

   write(*,'(a)') '   PASS'

contains

   !----- The canopy shim (meds_c_api_canopy) on meds_config_main.toml and the PFT file it names: it   !
   !      opens and reports the configuration's PFTs; the plastic capacities are positive and do not  !
   !      rise with the leaf area above, and are the PFT's own with plasticity off; a lit leaf          !
   !      assimilates through open stomata with Ci below Ca and a dark one respires; and the            !
   !      two-stream over two cohorts reflects, absorbs or passes to the ground every photon of VIS,   !
   !      absorbs nothing at night, and gives the taller cohort more light per unit of leaf area. -----!
   subroutine test_canopy_shim()
      character(len=*), parameter :: CONFIG = 'meds_config_main.toml'
      type(meds_config_t) :: cfg
      integer(c_int)      :: status, pft2(2), pft3(3)
      real(c_double)      :: lai_above(3), vc(3), rd(3)
      real(c_double)      :: a_gross(2), a_net(2), gs(2), ci(2)
      real(c_double), parameter :: COSZ(3) = [0.0_c_double, 0.3_c_double, 0.9_c_double]
      real(c_double), parameter :: VB(3) = [0.0_c_double, 100.0_c_double, 300.0_c_double]
      real(c_double), parameter :: VD(3) = [0.0_c_double, 80.0_c_double, 100.0_c_double]
      real(c_double), parameter :: HEIGHT(2) = [25.0_c_double, 5.0_c_double], LAI(2) = [3.0_c_double, 1.5_c_double]
      real(c_double), parameter :: WAI(2) = [0.4_c_double, 0.1_c_double]
      real(c_double) :: up_vis(3), up_nir(3), leaf_vis(3), wood_vis(3), ground_vis(3), coh_leaf_vis(3, 2)
      real(wp)       :: albedo, budget
      integer        :: h
      call load_meds_config(CONFIG, cfg)
      status = meds_canopy_open(CONFIG, int(len(CONFIG), c_int))
      call check(status == 0_c_int, 'C-API canopy: the reference configuration opens')
      call check(meds_canopy_n_pft() == int(cfg%pft%n, c_int), 'C-API canopy: it reports the configuration''s PFTs')

      !----- Plastic capacities under 0, 2 and 5 m2/m2 of leaf area above. ---------------------------!
      pft3      = 1_c_int
      lai_above = [0.0_c_double, 2.0_c_double, 5.0_c_double]
      call meds_canopy_plastic(3_c_int, pft3, lai_above, vc, rd)
      call check(all(vc > 0.0_c_double) .and. all(rd > 0.0_c_double), 'C-API canopy: plastic capacities are positive')
      call check(vc(2) <= vc(1) .and. vc(3) <= vc(2), 'C-API canopy: Vcmax25 does not rise with the leaf area above')
      if (.not. cfg%trait_dynamics%trait_plasticity_on) then
         call check(all(vc == cfg%pft%vcmax25(1)) .and. all(rd == cfg%pft%rd25(1)),                           &
                    'C-API canopy: with plasticity off, the PFT''s own capacities')
      end if

      !----- Two leaves of PFT 1 at its top-of-canopy capacities, 25 C, 1.5 kPa: lit and dark. --------!
      pft2 = 1_c_int
      call meds_canopy_leaf(2_c_int, pft2, [vc(1), vc(1)], [rd(1), rd(1)], [1500.0_c_double, 0.0_c_double],         &
                            [298.15_c_double, 298.15_c_double], [1500.0_c_double, 1500.0_c_double],                 &
                            [400.0_c_double, 400.0_c_double], [101325.0_c_double, 101325.0_c_double],              &
                            [0.0_c_double, 0.0_c_double], [0.0_c_double, 0.0_c_double], [0.0_c_double, 0.0_c_double], &
                            a_gross, a_net, gs, ci)
      call check(a_net(1) > 0.0_c_double .and. a_gross(1) > a_net(1), 'C-API canopy: a lit leaf assimilates')
      call check(gs(1) > 0.0_c_double .and. ci(1) > 0.0_c_double .and. ci(1) < 400.0_c_double,                 &
                 'C-API canopy: a lit leaf has open stomata and Ci below Ca')
      call check(a_gross(2) == 0.0_c_double .and. a_net(2) < 0.0_c_double, 'C-API canopy: a dark leaf respires')

      !----- The two-stream over three hours: night, a low sun and a high one; NIR at 1.1 x VIS. ------!
      call meds_canopy_two_stream(3_c_int, COSZ, VB, VD, 1.1_c_double * VB, 1.1_c_double * VD, 2_c_int, pft2,    &
                                  HEIGHT, LAI, WAI, up_vis, up_nir, leaf_vis, wood_vis, ground_vis, coh_leaf_vis)
      albedo = cfg%soil%ground_albedo_vis
      do h = 2, 3
         budget = real(up_vis(h) + leaf_vis(h) + wood_vis(h), wp) + (1.0_wp - albedo) * real(ground_vis(h), wp)
         call check(abs(budget - real(VB(h) + VD(h), wp)) <= 1.0e-9_wp * real(VB(h) + VD(h), wp),             &
                    'C-API canopy: every VIS photon is reflected, absorbed by the stand or by the ground')
         call check(abs(sum(coh_leaf_vis(h, :)) - leaf_vis(h)) <= 1.0e-12_wp * leaf_vis(h),                  &
                    'C-API canopy: the cohorts'' leaf VIS sums to the stand''s')
      end do
      call check(leaf_vis(1) == 0.0_c_double .and. leaf_vis(3) > leaf_vis(2) .and. leaf_vis(2) > 0.0_c_double,  &
                 'C-API canopy: no light absorbed at night, more under the higher sun')
      call check(coh_leaf_vis(3, 1) / LAI(1) > coh_leaf_vis(3, 2) / LAI(2),                                  &
                 'C-API canopy: the taller cohort absorbs more per unit of leaf area')
   end subroutine test_canopy_shim

   !----- Env at full light / 25 C, parameterised only by the stress potential under test. ---!
   pure function make_env(psi) result(e)
      real(c_double), intent(in) :: psi
      type(leaf_env_c) :: e
      e%par = 1500.0_c_double ; e%leaf_temp = 298.0_c_double ; e%vpd = 1500.0_c_double
      e%ca  = 400.0_c_double  ; e%pressure  = 101325.0_c_double
      e%psi_leaf = 0.0_c_double          ! demand-side limb off, so only `psi` moves gs
      e%gb  = 0.0_c_double
      e%psi = psi                         ! supply-side: what beta_stomata keys on
   end function make_env

   !----- A fixed C3 trait set for the shim's pass-through checks. The Python wrapper ships its !
   !      own defaults (meds/plant/leaf.py), which differ in theta_j, psi_open, psi_close and   !
   !      phi_psii. ------------------------------------------------------------------------------!
   pure function c3_defaults() result(q)
      type(leaf_args_c)   :: q
      q%pathway = int(PATH_C3, c_int)
      q%vcmax25 = 60.0_c_double  ; q%jmax25 = 108.0_c_double ; q%tpu25 = 1.0e6_c_double
      q%rd25    = 0.9_c_double   ; q%kp25   = 0.0_c_double
      q%g0      = 0.01_c_double  ; q%g1     = 4.0_c_double   ; q%d0    = 1500.0_c_double
      q%quantum_yield = 0.0_c_double
      q%theta_j     = 0.85_c_double
      q%theta_cj_c3 = 0.98_c_double ; q%theta_ip_c3 = 0.95_c_double
      q%theta_cj_c4 = 0.80_c_double ; q%theta_ic_c4 = 0.95_c_double
      q%lambda25 = 600.0_c_double
      q%psi_open = -0.5_c_double ; q%psi_close = -2.5_c_double ; q%lambda_psi_exp = 1.0_c_double
      q%sref_stomata = 2.0_c_double
      q%kc25 = 40.49_c_double ; q%ko25 = 27840.0_c_double ; q%gstar25 = 4.275_c_double
      q%ea_kc = 79430.0_c_double ; q%ea_ko = 36380.0_c_double ; q%ea_gstar = 37830.0_c_double
      q%ea_vcmax = 65330.0_c_double ; q%ea_jmax = 43540.0_c_double ; q%ea_rd = 46390.0_c_double
      q%hd_vcmax = 200000.0_c_double ; q%hd_jmax = 200000.0_c_double ; q%hd_rd = 1.0e9_c_double
      q%ds_vcmax = 650.0_c_double ; q%ds_jmax = 640.0_c_double ; q%ds_rd = 490.0_c_double
      q%o2_mol_frac = 0.209_c_double ; q%absorptance = 0.85_c_double ; q%phi_psii = 0.85_c_double
   end function c3_defaults

end program test_c_api_leaf
