! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_c_api_canopy -- the C-API shim for the CANOPY kernels (`edmore.canopy`): the leaf solve   !
! over many leaves and the two-stream over a frozen stand, with every parameter taken from a   !
! run's own configuration: the leaf table and the optics its loader derives on the PFT table    !
! (ed_pft_derive), so a caller outside the model computes what the fast loop computes.          !
!                                                                                          !
! The fast calibration (scripts/calibrate_fast) uses it for its optics and photosynthesis       !
! stages: the radiation solver alone against the tower's albedo, and a canopy of leaf solves     !
! against its GPP, without running the coupled model.                                           !
!                                                                                          !
! ONE configuration is open at a time (ed_canopy_open replaces it), and the calls are not        !
! thread-safe. A configuration with thermal acclimation on is refused: its leaf table depends on !
! a growth temperature the model carries as state, which a caller here does not have.            !
!==========================================================================================!
module ed_c_api_canopy
   use iso_c_binding
   use ed_kinds,                  only : wp, ik
   use ed_config,                 only : ed_config_t
   use ed_config_io,              only : load_ed_config
   use ed_leaf_gas_exchange, only : leaf_env_t, leaf_out_t
   use ed_pft_params,             only : leaf_photo_pft_t, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use ed_leaf_gas_exchange,      only : solve_leaf_gas_exchange, leaf_params_at_capacity
   use ed_plant_traits,           only : light_plastic_traits
   use ed_canopy_radiation, only : rad_env_t, rad_out_t, ground_optics_env_t, alloc_rad_env
   use ed_canopy_radiation,       only : canopy_radiation, ground_optics
   use ed_numerics,               only : ascending_order
   implicit none
   private

   public :: ed_canopy_open, ed_canopy_n_pft, ed_canopy_plastic, ed_canopy_leaf
   public :: ed_canopy_two_stream

   type(ed_config_t),   save :: g_cfg
   logical,             save :: g_open = .false.

contains

   pure function f_string(cbuf, n) result(s)
      character(kind=c_char), intent(in) :: cbuf(*)
      integer(c_int),         intent(in) :: n
      character(len=n) :: s
      integer :: i
      do i = 1, n
         s(i:i) = cbuf(i)
      end do
   end function f_string

   !----- Load a run's main configuration (and the PFT file it names), with the leaf table and the  !
   !      radiation optics its loader derives. 0 on success; -1 when thermal acclimation is on. ----!
   function ed_canopy_open(path, path_len) result(status) bind(c, name="ed_canopy_open")
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value,  intent(in) :: path_len
      integer(c_int)                     :: status
      g_open = .false.
      call load_ed_config(f_string(path, path_len), g_cfg)
      if (g_cfg%leaf_physiology%thermal_acclimation) then
         status = -1_c_int
         return
      end if
      g_open = .true.
      status = 0_c_int
   end function ed_canopy_open

   function ed_canopy_n_pft() result(n) bind(c, name="ed_canopy_n_pft")
      integer(c_int) :: n
      n = 0_c_int
      if (g_open) n = int(g_cfg%pft%n, c_int)
   end function ed_canopy_n_pft

   !----- Each cohort's plastic Vcmax25 and Rd25 from its PFT and the leaf area above it, as a      !
   !      restart re-acclimates them (reacclimate_plant_traits): the top-of-canopy value times the   !
   !      light gradient with plasticity on, the PFT's value with it off. ------------------------!
   subroutine ed_canopy_plastic(n, pft, lai_above, vcmax25, rd25) bind(c, name="ed_canopy_plastic")
      integer(c_int), value, intent(in)  :: n
      integer(c_int),        intent(in)  :: pft(n)
      real(c_double),        intent(in)  :: lai_above(n)
      real(c_double),        intent(out) :: vcmax25(n), rd25(n)
      integer(ik) :: i, pf
      real(wp)    :: sla_t, vc_t, rd_t, turnover_t
      associate (t => g_cfg%pft)
         do i = 1_ik, int(n, ik)
            pf = int(pft(i), ik)
            if (g_cfg%trait_dynamics%trait_plasticity_on) then
               call light_plastic_traits(real(lai_above(i), wp), t%sla(pf), t%vcmax25(pf), t%rd25(pf),  &
                                         t%leaf_turnover_rate_toc(pf), t%kplastic_sla(pf), t%kplastic_vcmax(pf), &
                                         t%kplastic_rd(pf), t%kplastic_leaf_turnover(pf), sla_t, vc_t, rd_t,   &
                                         turnover_t)
            else
               vc_t = t%vcmax25(pf)
               rd_t = t%rd25(pf)
            end if
            vcmax25(i) = vc_t
            rd25(i) = rd_t
         end do
      end associate
   end subroutine ed_canopy_plastic

   !----- The leaf solve for n leaves, each with its own drivers and plastic capacities: the loop   !
   !      of leaf_gas_exchange_batch, with canopy-air CO2 and pressure per leaf rather than shared. --!
   subroutine ed_canopy_leaf(n, pft, vcmax25, rd25, par, leaf_temp, vpd, ca, pressure, psi_leaf,   &
                               gb, psi, a_gross, a_net, gs, ci) bind(c, name="ed_canopy_leaf")
      integer(c_int), value, intent(in)  :: n
      integer(c_int),        intent(in)  :: pft(n)
      real(c_double),        intent(in)  :: vcmax25(n), rd25(n)
      real(c_double),        intent(in)  :: par(n), leaf_temp(n), vpd(n), ca(n), pressure(n)
      real(c_double),        intent(in)  :: psi_leaf(n), gb(n), psi(n)
      real(c_double),        intent(out) :: a_gross(n), a_net(n), gs(n), ci(n)
      type(leaf_env_t)          :: env
      type(leaf_out_t)         :: flux
      type(leaf_photo_pft_t)    :: p
      integer(ik) :: i
      do i = 1_ik, int(n, ik)
         env%par = par(i)
         env%leaf_temp = leaf_temp(i)
         env%vpd = vpd(i)
         env%ca = ca(i)
         env%pressure = pressure(i)
         env%psi_leaf = psi_leaf(i)
         env%gb = gb(i)
         env%psi = psi(i)
         associate (leaf_photo => g_cfg%pft%leaf_photo)
            p = leaf_params_at_capacity(leaf_photo, int(pft(i), ik), real(vcmax25(i), wp), real(rd25(i), wp))
            call solve_leaf_gas_exchange(env, p, leaf_photo%stomatal_model, leaf_photo%temp_response_form,   &
                                         leaf_photo%colimitation, leaf_photo%use_boundary_layer, flux)
         end associate
         a_gross(i) = flux%A_gross
         a_net(i) = flux%A_net
         gs(i) = flux%gs
         ci(i) = flux%ci
      end do
   end subroutine ed_canopy_leaf

   !----- The two-stream for ONE patch over nh hours of shortwave, the stand fixed: per hour the     !
   !      upwelling VIS and NIR at the canopy top, the VIS absorbed by leaves and by wood, the VIS     !
   !      reaching the ground, and each cohort's leaf-absorbed VIS (nh x ncoh, hour fastest). Cohorts  !
   !      come in any order and are sorted by height as the fast loop does (apply_rt_forcing). The     !
   !      ground is the configuration's bare soil; no snow. ------------------------------------------!
   subroutine ed_canopy_two_stream(nh, cosz, vis_beam, vis_diff, nir_beam, nir_diff, ncoh, pft,    &
                                    height, lai, wai, up_vis, up_nir, leaf_vis, wood_vis, ground_vis, &
                                    coh_leaf_vis) bind(c, name="ed_canopy_two_stream")
      integer(c_int), value, intent(in)  :: nh, ncoh
      real(c_double),        intent(in)  :: cosz(nh), vis_beam(nh), vis_diff(nh), nir_beam(nh), nir_diff(nh)
      integer(c_int),        intent(in)  :: pft(ncoh)
      real(c_double),        intent(in)  :: height(ncoh), lai(ncoh), wai(ncoh)
      real(c_double),        intent(out) :: up_vis(nh), up_nir(nh), leaf_vis(nh), wood_vis(nh), ground_vis(nh)
      real(c_double),        intent(out) :: coh_leaf_vis(nh, ncoh)
      integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT
      integer(ik) :: perm(max(ncoh,1)), pft_bt(max(ncoh,1))
      real(wp)    :: hgt_bt(max(ncoh,1)), lai_bt(max(ncoh,1)), wai_bt(max(ncoh,1)), tcan_bt(max(ncoh,1))
      type(rad_env_t)         :: rf
      type(rad_out_t)            :: flux
      type(ground_optics_env_t) :: surf
      logical     :: he(NB)
      integer(ik) :: h, j, nc
      nc = int(ncoh, ik)
      call ascending_order(real(height, wp), nc, perm(1:nc))
      do j = 1_ik, nc
         pft_bt(j) = int(pft(perm(j)), ik)
         hgt_bt(j) = height(perm(j))
         lai_bt(j) = lai(perm(j))
         wai_bt(j) = wai(perm(j))
      end do
      tcan_bt = 298.15_wp                                     ! longwave only; its outputs are unused
      call alloc_rad_env(rf, NB)
      surf%n_band = NB
      allocate(surf%soil_albedo(NB))
      surf%soil_albedo = [g_cfg%soil%ground_albedo_vis, g_cfg%soil%ground_albedo_nir, 0.0_wp]
      surf%soil_emiss  = g_cfg%soil%ground_emissivity
      surf%soil_temp   = 298.15_wp
      he = [.false., .false., .true.]
      call ground_optics(surf, NB, he, rf%grnd_refl, rf%grnd_emiss)
      do h = 1_ik, int(nh, ik)
         rf%cosz = cosz(h)
         rf%incid_beam(RAD_VIS) = vis_beam(h)
         rf%incid_diff(RAD_VIS) = vis_diff(h)
         rf%incid_beam(RAD_NIR) = nir_beam(h)
         rf%incid_diff(RAD_NIR) = nir_diff(h)
         rf%incid_beam(RAD_LW)  = 0.0_wp
         rf%incid_diff(RAD_LW)  = 0.0_wp
         call canopy_radiation(g_cfg%pft%optics, rf, nc, pft_bt(1:nc), hgt_bt(1:nc), lai_bt(1:nc),       &
                               wai_bt(1:nc), tcan_bt(1:nc), tcan_bt(1:nc), flux)
         up_vis(h)     = flux%albedo(RAD_VIS) * (vis_beam(h) + vis_diff(h))
         up_nir(h)     = flux%albedo(RAD_NIR) * (nir_beam(h) + nir_diff(h))
         leaf_vis(h)   = sum(flux%abs_leaf(RAD_VIS, 1:nc))
         wood_vis(h)   = sum(flux%abs_wood(RAD_VIS, 1:nc))
         ground_vis(h) = flux%dn_ground(RAD_VIS)
         do j = 1_ik, nc
            coh_leaf_vis(h, perm(j)) = flux%abs_leaf(RAD_VIS, j)
         end do
      end do
   end subroutine ed_canopy_two_stream

end module ed_c_api_canopy
