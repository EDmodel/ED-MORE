! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_pft_derive -- what the PFT file's traits imply, derived once at load                 !
! (derive_parameters): the leaf capacities from Vcmax25, the stomatal stress points a PFT  !
! file may leave out, and the three tables the fast loop reads, kept on the PFT table: each !
! PFT's leaf-photosynthesis parameters, its hydraulics with its Kirchhoff lookup (built     !
! here, build_kirchhoff_table), and the                                                    !
! two-stream's optics. The laws these evaluate are share/functions'.                       !
!==========================================================================================!
module ed_pft_derive
   use ed_kinds,            only : wp, ik
   use ed_pft_params,       only : pft_table_t, HYD_UNSET, leaf_photo_pft_t, leaf_photo_table_t,           &
                                   hydraulics_table_t, rad_pft_optics_t, N_RAD_BAND_DEFAULT, RAD_VIS,       &
                                   RAD_NIR, RAD_LW, kirchhoff_table_t, KIRCHHOFF_TABLE_NTAB,                &
                                   KIRCHHOFF_TABLE_RMAX
   use ed_plant_hydr_lib,   only : tissue_psi_tlp, tissue_psi_to_rwc
   use ed_allometry_lib,    only : min_cohort_carbon
   use ed_optics_lib,       only : N_LEAF_CLASS, beta_lidf, leaf_bf, scatter_pair
   implicit none
   private

   public :: derive_leaf_params, derive_pft_tables, derive_recruit_carbon
   public :: pft_leaf_psi_tlp, pft_stomata_psi_onset, pft_wstress_psi_open, pft_wstress_psi_close
   public :: derive_rad_optics, moments_to_beta_params
   public :: build_kirchhoff_table

contains

   !---------------------------------------------------------------------------------------!
   ! Derive the per-PFT leaf-photosynthesis capacities (Jmax25, TPU25, Rd25) as fixed      !
   ! ratios of Vcmax25. Call after vcmax25 and the ratio arrays are set from the PFT file. !
   !---------------------------------------------------------------------------------------!
   subroutine derive_leaf_params(pft)
      type(pft_table_t), intent(inout) :: pft
      pft%jmax25 = pft%jmax_vcmax_ratio * pft%vcmax25
      pft%tpu25  = pft%tpu_vcmax_ratio  * pft%vcmax25
      pft%rd25   = pft%rd_vcmax_ratio   * pft%vcmax25
   end subroutine derive_leaf_params

   !----- The three tables the fast loop reads, and the phenology's turgor-loss point, onto the  !
   !      PFT table. The leaf table takes                                                         !
   !      [leaf_physiology]: its constants (`shared`, the part every PFT's parameters share) and    !
   !      its solver selectors. Call after derive_leaf_params. ------------------------------------!
   subroutine derive_pft_tables(pft, shared, stomatal_model, temp_response_form, colimitation,          &
                                use_boundary_layer)
      type(pft_table_t),      intent(inout) :: pft
      type(leaf_photo_pft_t), intent(in)    :: shared
      integer(ik),            intent(in)    :: stomatal_model, temp_response_form, colimitation
      logical,                intent(in)    :: use_boundary_layer
      type(leaf_photo_table_t) :: leaf_photo
      type(hydraulics_table_t) :: hydraulics
      type(rad_pft_optics_t)   :: optics
      integer(ik)              :: ipft
      call build_leaf_photo_table(pft, shared, leaf_photo)
      leaf_photo%stomatal_model     = stomatal_model
      leaf_photo%temp_response_form = temp_response_form
      leaf_photo%colimitation       = colimitation
      leaf_photo%use_boundary_layer = use_boundary_layer
      call build_hydraulics_table(pft, hydraulics)
      call build_rad_optics(pft, optics)
      pft%leaf_photo = leaf_photo
      pft%hydraulics = hydraulics
      pft%optics     = optics
      !----- The turgor-loss point the phenology's drought counter reads. -----------------!
      do ipft = 1_ik, pft%n
         pft%pheno(ipft)%leaf_psi_tlp = pft_leaf_psi_tlp(pft, ipft)
      end do
   end subroutine derive_pft_tables

   !----- Each PFT's recruit carbon, from the table's own allometry: called on a polygon's copy, !
   !      once the cell's allometry is in it (docs/science/demography.md §3). -------------------!
   subroutine derive_recruit_carbon(pft)
      type(pft_table_t), intent(inout) :: pft
      integer(ik) :: ipft
      if (allocated(pft%recruit_carbon)) deallocate(pft%recruit_carbon)
      allocate(pft%recruit_carbon(pft%n))
      do ipft = 1_ik, pft%n
         pft%recruit_carbon(ipft) = min_cohort_carbon(pft%allom(ipft), pft%min_cohort_height,           &
                                                      pft%wood_density(ipft))
      end do
   end subroutine derive_recruit_carbon

   !=======================================================================================!
   !  The stomatal stress points.                                                          !
   !=======================================================================================!

   !----- A PFT's leaf turgor-loss point [MPa], from its own leaf pressure-volume curve, so the      !
   !      stomatal closure and the phenology's drought counter see the curve the solver stores leaf  !
   !      water on. -----------------------------------------------------------------------------------!
   pure real(wp) function pft_leaf_psi_tlp(pft, ipft) result(psi_tlp)
      type(pft_table_t), intent(in) :: pft
      integer(ik),       intent(in) :: ipft
      associate (leaf_curve => pft%hydro(ipft)%leaf_curve)
         psi_tlp = tissue_psi_tlp(leaf_curve%pi0, leaf_curve%elastic_mod)
      end associate
   end function pft_leaf_psi_tlp

   !----- Where a PFT's stomatal water stress begins [MPa of predawn leaf potential]: the [pft]       !
   !      stomata_psi_onset where given, else half the PFT's turgor-loss point, so the gravity head  !
   !      of a tall tree in wet soil does not read as drought (docs/science/leaf_gas_exchange.md     !
   !      eq. (10)). ----------------------------------------------------------------------------------!
   pure real(wp) function pft_stomata_psi_onset(pft, ipft) result(psi_onset)
      type(pft_table_t), intent(in) :: pft
      integer(ik),       intent(in) :: ipft
      psi_onset = 0.5_wp * pft_leaf_psi_tlp(pft, ipft)
      if (allocated(pft%stomata_psi_onset)) then
         if (pft%stomata_psi_onset(ipft) > HYD_UNSET) psi_onset = pft%stomata_psi_onset(ipft)
      end if
   end function pft_stomata_psi_onset

   !----- The leaf potential [MPa] at which the fast water-stress limb on Vcmax/Jmax/TPU starts and   !
   !      completes: the [pft] wstress_psi_open / wstress_psi_close where given, else the PFT's      !
   !      turgor-loss point and twice it (docs/science/leaf_gas_exchange.md eq. (9)). ----------------!
   pure real(wp) function pft_wstress_psi_open(pft, ipft) result(psi_open)
      type(pft_table_t), intent(in) :: pft
      integer(ik),       intent(in) :: ipft
      psi_open = pft_leaf_psi_tlp(pft, ipft)
      if (allocated(pft%wstress_psi_open)) then
         if (pft%wstress_psi_open(ipft) > HYD_UNSET) psi_open = pft%wstress_psi_open(ipft)
      end if
   end function pft_wstress_psi_open

   pure real(wp) function pft_wstress_psi_close(pft, ipft) result(psi_close)
      type(pft_table_t), intent(in) :: pft
      integer(ik),       intent(in) :: ipft
      psi_close = 2.0_wp * pft_leaf_psi_tlp(pft, ipft)
      if (allocated(pft%wstress_psi_close)) then
         if (pft%wstress_psi_close(ipft) > HYD_UNSET) psi_close = pft%wstress_psi_close(ipft)
      end if
   end function pft_wstress_psi_close

   !=======================================================================================!
   !  The leaf table.                                                                      !
   !=======================================================================================!

   !----- Every PFT's leaf parameters: the shared constants, each PFT's traits and capacities, and   !
   !      the two ratios a cohort's plastic Vcmax25 scales Jmax25 and TPU25 by. -------------------!
   subroutine build_leaf_photo_table(pft, shared, table)
      type(pft_table_t),        intent(in)  :: pft
      type(leaf_photo_pft_t),   intent(in)  :: shared
      type(leaf_photo_table_t), intent(out) :: table
      integer(ik) :: ipft
      table%n_pft = pft%n
      allocate(table%pft(table%n_pft), table%jmax_vcmax_ratio(table%n_pft), table%tpu_vcmax_ratio(table%n_pft))
      do ipft = 1_ik, table%n_pft
         table%pft(ipft) = leaf_photo_for_pft(pft, ipft, shared)
         table%jmax_vcmax_ratio(ipft) = pft%jmax_vcmax_ratio(ipft)
         table%tpu_vcmax_ratio(ipft)  = pft%tpu_vcmax_ratio(ipft)
      end do
   end subroutine build_leaf_photo_table

   !----- PFT ipft's leaf-photosynthesis parameters: the shared constants with the PFT's traits and   !
   !      its table Vcmax25/Jmax25/TPU25/Rd25. leaf_params_at_capacity applies a cohort's plastic      !
   !      overrides on top of it. ----------------------------------------------------------------------!
   pure function leaf_photo_for_pft(pft, ipft, shared) result(p)
      type(pft_table_t),      intent(in) :: pft
      integer(ik),            intent(in) :: ipft
      type(leaf_photo_pft_t), intent(in) :: shared
      type(leaf_photo_pft_t)             :: p
      p = shared
      p%pathway        = pft%photosynthetic_pathway(ipft)
      p%vcmax25        = pft%vcmax25(ipft)
      p%jmax25         = pft%jmax25(ipft)
      p%tpu25          = pft%tpu25(ipft)
      p%rd25           = pft%rd25(ipft)
      p%kp25           = pft%kp25(ipft)
      p%g0             = pft%stomatal_g0(ipft)
      p%g1             = pft%stomatal_g1(ipft)
      p%d0             = pft%stomatal_d0(ipft)
      p%quantum_yield  = pft%quantum_yield_c4(ipft)
      p%theta_j        = pft%theta_j(ipft)
      p%theta_cj_c3    = pft%theta_cj_c3(ipft)
      p%theta_ip_c3    = pft%theta_ip_c3(ipft)
      p%theta_cj_c4    = pft%theta_cj_c4(ipft)
      p%theta_ic_c4    = pft%theta_ic_c4(ipft)
      p%lambda25       = pft%katul_lambda25(ipft)
      p%psi_open       = pft_wstress_psi_open(pft, ipft)    ! turgor loss unless the PFT file sets it
      p%psi_close      = pft_wstress_psi_close(pft, ipft)   ! twice it unless the PFT file sets it
      p%lambda_psi_exp = pft%wstress_lambda_exp(ipft)
      p%sref_stomata   = pft%wstress_sref_stomata(ipft)
      p%psi_onset      = pft_stomata_psi_onset(pft, ipft)
   end function leaf_photo_for_pft

   !=======================================================================================!
   !  The hydraulics table.                                                                !
   !=======================================================================================!

   !----- Each PFT's [hydraulics] traits, the Kirchhoff lookup built from its own wood_kexp         !
   !      (docs/science/plant_hydraulics.md eq. (8)), which makes a non-integer vulnerability       !
   !      exponent affordable on the hot path, and the leaf water content at which its stomata      !
   !      close, the water network's supply limit (soil_plant_water.md eq. (3)). -----------------!
   subroutine build_hydraulics_table(pft, table)
      type(pft_table_t),        intent(in)  :: pft
      type(hydraulics_table_t), intent(out) :: table
      integer(ik) :: i
      table%n_pft = pft%n
      allocate(table%pft(max(pft%n, 1_ik)))
      do i = 1_ik, pft%n
         table%pft(i)%plant_hydraulics_t = pft%hydro(i)
         call build_kirchhoff_table(table%pft(i)%vuln_table, table%pft(i)%wood_kexp)
         table%pft(i)%rwc_close = tissue_psi_to_rwc(pft_wstress_psi_close(pft, i), table%pft(i)%leaf_curve%pi0,   &
                                                    table%pft(i)%leaf_curve%elastic_mod)
      end do
   end subroutine build_hydraulics_table

   !----- Fill a PFT's Kirchhoff table (docs/science/plant_hydraulics.md eq. (8)), one     !
   !      7-point panel per interval added to the node below, so it rises monotonically.   !
   !      A subroutine, not an array-valued function (docs/building.md, Compilers). ------!
   pure subroutine build_kirchhoff_table(tab, kexp)
      type(kirchhoff_table_t), intent(out) :: tab
      real(wp),                intent(in)  :: kexp
      real(wp) :: h
      integer  :: j
      tab%kexp  = kexp
      tab%r_max = KIRCHHOFF_TABLE_RMAX
      tab%n     = int(KIRCHHOFF_TABLE_NTAB, ik)
      h         = KIRCHHOFF_TABLE_RMAX / real(KIRCHHOFF_TABLE_NTAB, wp)
      tab%inv_h = 1.0_wp / h
      tab%g(0) = 0.0_wp
      do j = 1, KIRCHHOFF_TABLE_NTAB
         tab%g(j) = tab%g(j-1) + kirchhoff_panel(real(j-1, wp)*h, real(j, wp)*h, kexp)
      end do
   end subroutine build_kirchhoff_table

   !----- integral_a^b du/(1+u^kexp) on one panel (0 <= a <= b) by 7-point Gauss-Legendre, !
   !      exact for polynomials up to degree 13. The sum is written out rather than handed !
   !      a function (docs/building.md, no procedure arguments); u = 0 is guarded as in    !
   !      conductance_retained_fraction. --------------------------------------------------!
   pure real(wp) function kirchhoff_panel(a, b, kexp) result(g)
      real(wp), intent(in) :: a, b, kexp
      real(wp), parameter :: NODE(7)   = [ -0.9491079123427585_wp, -0.7415311855993945_wp,     &
                                           -0.4058451513773972_wp,  0.0000000000000000_wp,     &
                                            0.4058451513773972_wp,  0.7415311855993945_wp,     &
                                            0.9491079123427585_wp ]
      real(wp), parameter :: WEIGHT(7) = [  0.1294849661688697_wp,  0.2797053914892766_wp,     &
                                            0.3818300505051189_wp,  0.4179591836734694_wp,     &
                                            0.3818300505051189_wp,  0.2797053914892766_wp,     &
                                            0.1294849661688697_wp ]
      real(wp)    :: mid, half, u, acc
      integer(ik) :: k
      mid  = 0.5_wp * (a + b)                    ! the nodes on [-1,1], mapped onto [a,b]
      half = 0.5_wp * (b - a)
      acc  = 0.0_wp
      do k = 1_ik, 7_ik
         u = mid + half * NODE(k)
         if (u <= 0.0_wp) then
            acc = acc + WEIGHT(k)
         else
            acc = acc + WEIGHT(k) / (1.0_wp + u**kexp)
         end if
      end do
      g = half * acc
   end function kirchhoff_panel

   !=======================================================================================!
   !  The optics table (docs/science/canopy_radiation_transfer.md §1-§2).                  !
   !=======================================================================================!

   !----- The two-stream's per-PFT optics from the [pft] traits; the canopy C API reads the same. ---!
   subroutine build_rad_optics(pft, optics)
      type(pft_table_t),      intent(in)  :: pft
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT
      integer(ik) :: np, ipf
      real(wp), allocatable :: rl(:,:), tl(:,:), rw(:,:), tw(:,:), cl(:), cw(:), bp(:), bq(:)
      logical  :: hb(NB), he(NB)
      real(wp) :: bpp, bqq
      np = pft%n
      allocate(rl(NB,np), tl(NB,np), rw(NB,np), tw(NB,np), cl(np), cw(np), bp(np), bq(np))
      !----- Shortwave comes as reflectance and transmittance; longwave as emissivity, whose band   !
      !      reflectance is 1 - emissivity with zero transmittance (a leaf is opaque there). ---!
      rl(RAD_VIS,1:np) = pft%leaf_reflect_vis(1:np)
      tl(RAD_VIS,1:np) = pft%leaf_transmit_vis(1:np)
      rl(RAD_NIR,1:np) = pft%leaf_reflect_nir(1:np)
      tl(RAD_NIR,1:np) = pft%leaf_transmit_nir(1:np)
      rl(RAD_LW ,1:np) = 1.0_wp - pft%leaf_emissivity(1:np)
      tl(RAD_LW,1:np) = 0.0_wp
      rw(RAD_VIS,1:np) = pft%wood_reflect_vis(1:np)
      tw(RAD_VIS,1:np) = pft%wood_transmit_vis(1:np)
      rw(RAD_NIR,1:np) = pft%wood_reflect_nir(1:np)
      tw(RAD_NIR,1:np) = pft%wood_transmit_nir(1:np)
      rw(RAD_LW ,1:np) = 1.0_wp - pft%wood_emissivity(1:np)
      tw(RAD_LW,1:np) = 0.0_wp
      cl(1:np) = pft%leaf_clumping(1:np)
      cw(1:np) = pft%wood_clumping(1:np)
      !----- The Beta leaf-angle shape parameters, from the angle's measured mean and std. ---!
      do ipf = 1_ik, np
         call moments_to_beta_params(pft%leaf_angle_mean(ipf), pft%leaf_angle_std(ipf), bpp, bqq)
         bp(ipf) = bpp
         bq(ipf) = bqq
      end do
      hb = [.true.,  .true.,  .false.]                         ! VIS/NIR have a beam; LW does not
      he = [.false., .false., .true. ]                         ! only LW emits
      call derive_rad_optics(NB, np, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, optics)
   end subroutine build_rad_optics

   !----- Beta (p, q) from the mean and standard deviation of the leaf angle [deg], by     !
   !      moment matching on t = theta/90 (Goel & Strebel 1984). --------------------------!
   pure subroutine moments_to_beta_params(mean_deg, std_deg, p, q)
      real(wp), intent(in)  :: mean_deg, std_deg
      real(wp), intent(out) :: p, q
      real(wp) :: mt, vt, kappa
      real(wp), parameter :: VT_MIN = 1.0e-6_wp   ! variance floor (std_deg floor ~0.09 deg)
      mt = min(max(mean_deg / 90.0_wp, 1.0e-3_wp), 1.0_wp - 1.0e-3_wp)   ! clamp-ok: the beta mean inside (0, 1)
      vt = (std_deg / 90.0_wp) ** 2
      vt = min(vt, mt * (1.0_wp - mt) * (1.0_wp - 1.0e-6_wp))             ! a Beta needs vt < mt (1 - mt)
      vt = max(vt, VT_MIN)                                                ! keep kappa finite at std_deg = 0
      kappa = mt * (1.0_wp - mt) / vt - 1.0_wp                            ! concentration, p + q
      p = mt * kappa
      q = (1.0_wp - mt) * kappa
   end subroutine moments_to_beta_params

   !---------------------------------------------------------------------------------------!
   ! Fill the per-PFT optics table (docs/science/canopy_radiation_transfer.md eqs. (1),    !
   ! (2) and (4)) from (band, pft) spectra: reflectance and transmittance for a shortwave  !
   ! band, 1 - emissivity and 0 for the thermal band; beta_p, beta_q shape the Beta LIDF.  !
   !---------------------------------------------------------------------------------------!
   subroutine derive_rad_optics(n_band, n_pft, reflect_leaf, trans_leaf, reflect_wood,       &
                                trans_wood, clumping_leaf, clumping_wood, beta_p, beta_q,     &
                                has_beam, has_emission, optics)
      integer(ik), intent(in) :: n_band, n_pft
      real(wp),    intent(in) :: reflect_leaf(n_band, n_pft), trans_leaf(n_band, n_pft)
      real(wp),    intent(in) :: reflect_wood(n_band, n_pft), trans_wood(n_band, n_pft)
      real(wp),    intent(in) :: clumping_leaf(n_pft), clumping_wood(n_pft)
      real(wp),    intent(in) :: beta_p(n_pft), beta_q(n_pft)
      logical,     intent(in) :: has_beam(n_band), has_emission(n_band)
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik) :: b, p

      call alloc_rad_pft_optics(optics, n_band, n_pft, N_LEAF_CLASS)
      optics%clumping_leaf = clumping_leaf
      optics%clumping_wood = clumping_wood
      optics%has_beam      = has_beam
      optics%has_emission  = has_emission

      do p = 1_ik, n_pft
         optics%lidf(:,p) = beta_lidf(beta_p(p), beta_q(p))
         optics%bf(p)     = leaf_bf(optics%lidf(:,p))
         do b = 1_ik, n_band
            call scatter_pair(reflect_leaf(b,p), trans_leaf(b,p), optics%bf(p),               &
                              optics%omega_leaf(b,p), optics%g_leaf(b,p))
            call scatter_pair(reflect_wood(b,p), trans_wood(b,p), optics%bf(p),               &
                              optics%omega_wood(b,p), optics%g_wood(b,p))
         end do
      end do
   end subroutine derive_rad_optics

   subroutine alloc_rad_pft_optics(optics, n_band, n_pft, n_class)
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik),            intent(in)  :: n_band, n_pft, n_class
      optics%n_band = n_band
      optics%n_pft  = n_pft
      allocate(optics%omega_leaf(n_band, n_pft), optics%omega_wood(n_band, n_pft))
      allocate(optics%g_leaf(n_band, n_pft),     optics%g_wood(n_band, n_pft))
      allocate(optics%clumping_leaf(n_pft),      optics%clumping_wood(n_pft))
      allocate(optics%lidf(n_class, n_pft),      optics%bf(n_pft))
      allocate(optics%has_beam(n_band),          optics%has_emission(n_band))
      optics%omega_leaf = 0.0_wp
      optics%omega_wood = 0.0_wp
      optics%g_leaf = 0.0_wp
      optics%g_wood = 0.0_wp
      optics%clumping_leaf = 1.0_wp
      optics%clumping_wood = 1.0_wp
      optics%lidf = 0.0_wp
      optics%bf = 0.0_wp
      optics%has_beam = .false.
      optics%has_emission = .false.
   end subroutine alloc_rad_pft_optics

end module ed_pft_derive
