! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_canopy_radiation -- the canopy two-stream radiative transfer                          !
! (docs/science/canopy_radiation_transfer.md), on plain arrays with no state of its own:   !
! cohorts run bottom (1) to top (ncoh) and every flux is absolute W/m2. ground_optics sets !
! the ground's reflectance and emission, the lower boundary (§4). canopy_radiation blends  !
! each cohort's leaf and wood optics (§2), solves each band by the adding method           !
! (solve_band, layer_rt; §3) and splits each cohort's absorbed radiation into leaves and   !
! wood, each tissue emitting at its own temperature (§4). The per-PFT optics table is      !
! built at load (ed_pft_derive).                                                           !
!==========================================================================================!
module ed_canopy_radiation
   use ed_kinds,          only : wp, ik
   use ed_constants,      only : stefan_boltzmann, tiny_num, lnexp_max
   use ed_pft_params,     only : rad_pft_optics_t
   use ed_optics_lib,     only : gfun_direct
   use ed_numerics,       only : equal_height_layers
   implicit none
   private

   public :: ground_optics
   public :: canopy_radiation
   public :: rad_env_t, rad_out_t, ground_optics_env_t, alloc_rad_env
   public :: layer_beam_sources   ! used only by the tests

   !----- The largest single-scatter albedo a layer takes, off conservative scattering (omega =    !
   !      1), where the two-stream's diffuse rate m vanishes. -------------------------------------!
   real(wp), parameter :: OMEGA_MAX = 1.0_wp - 1.0e-6_wp     !< [-]
   !----- The relative distance the beam's extinction is moved off the diffuse rate m, where the   !
   !      in-layer beam solution has a removable singularity. -------------------------------------!
   real(wp), parameter :: RESONANCE_REL = 1.0e-6_wp          !< [-]

   !----- Per-band incident forcing and the ground boundary. In the thermal band the diffuse is    !
   !      the sky's longwave and grnd_emiss the ground's emission; in a shortwave band grnd_emiss  !
   !      is 0. -----------------------------------------------------------------------------------!
   type :: rad_env_t
      integer(ik) :: n_band = 0_ik                 !< number of bands
      real(wp)    :: cosz   = 1.0_wp               !< [-] cosine of the solar zenith angle (the beam extinction guards cosz -> 0)
      real(wp), allocatable :: incid_beam(:)       !< (band) [W/m2] direct-beam incident at canopy top
      real(wp), allocatable :: incid_diff(:)       !< (band) [W/m2] diffuse incident at canopy top
      real(wp), allocatable :: grnd_refl(:)        !< (band) [-] ground reflectance (albedo, or 1 - emissivity)
      real(wp), allocatable :: grnd_emiss(:)       !< (band) [W/m2] ground thermal emission (0 for SW)
   end type rad_env_t

   !----- Returned fluxes: per band and cohort, the radiation absorbed by leaves and by wood [W/m2 !
   !      of ground]; per band, the albedo and the fluxes below the canopy. -----------------------!
   type :: rad_out_t
      integer(ik) :: n_band = 0_ik, n_coh = 0_ik   !< number of bands, number of cohorts
      real(wp), allocatable :: abs_leaf(:,:)       !< (band,cohort) [W/m2] absorbed by leaves
      real(wp), allocatable :: abs_wood(:,:)       !< (band,cohort) [W/m2] absorbed by wood
      real(wp), allocatable :: incid_top(:,:)      !< (band,cohort) [W/m2] beam + downward diffuse at its layer's top
      real(wp), allocatable :: albedo(:)           !< (band) [-] upward over incident at top; LW's includes emission (can be > 1)
      real(wp), allocatable :: dn_ground(:)        !< (band) [W/m2] downwelling below canopy (to ground)
      real(wp), allocatable :: up_ground(:)        !< (band) [W/m2] upwelling from ground into canopy
   end type rad_out_t

   !----- The ground's optical state, the two-stream's lower boundary: bare soil, blended toward   !
   !      the snowpack's optics (docs/science/canopy_radiation_transfer.md §4). -------------------!
   type :: ground_optics_env_t
      integer(ik)           :: n_band = 0_ik     !< number of bands
      real(wp), allocatable :: soil_albedo(:)    !< (band) [-] shortwave ground albedo; unused for emission bands
      real(wp)              :: soil_emiss = 0.95_wp   !< [-] thermal emissivity of the ground
      real(wp)              :: soil_temp  = 298.0_wp  !< [K] ground (skin) temperature
   end type ground_optics_env_t

contains

   !---------------------------------------------------------------------------------------!
   ! Blend leaf and wood optics for every cohort, for one band, and add the solar geometry !
   ! (docs/science/canopy_radiation_transfer.md eqs. (5)-(7)): omega, beta, beta0, k and   !
   ! leaf_frac, from the clumping-corrected leaf and wood area indices elai, ewai. A       !
   ! cohort with no leaf or wood area takes the leaf optics, and one that absorbs nothing  !
   ! has leaf_frac = 1.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine blend_cohort_optics(optics, band, ncoh, pft, elai, ewai, kdir_in,              &
                                  omega, beta, beta0, kdir, leaf_frac)
      type(rad_pft_optics_t), intent(in)  :: optics     !< the per-PFT optics table
      integer(ik),            intent(in)  :: band, ncoh !< band index, number of cohorts
      integer(ik),            intent(in)  :: pft(ncoh)  !< each cohort's PFT
      real(wp),               intent(in)  :: elai(ncoh), ewai(ncoh)   !< [m2/m2] clumping-corrected leaf, wood area index
      real(wp),               intent(in)  :: kdir_in(ncoh)   !< [-] beam extinction per unit area index, the same in every band
      real(wp),               intent(out) :: omega(ncoh), beta(ncoh), beta0(ncoh)   !< [-] the cohort's blended optics
      real(wp),               intent(out) :: kdir(ncoh), leaf_frac(ncoh)   !< [-] beam extinction (0 with no beam), leaf share
      integer(ik) :: i, ip
      real(wp)    :: etai, wl, ww, gc, aleaf, awood

      do i = 1_ik, ncoh
         ip   = pft(i)
         etai = elai(i) + ewai(i)
         if (etai > tiny_num) then
            wl = elai(i) / etai
            ww = ewai(i) / etai
         else
            wl = 1.0_wp
            ww = 0.0_wp
         end if
         !----- Area-weighted single-scatter albedo and leaf-angle asymmetry. --------------!
         omega(i) = wl * optics%omega_leaf(band,ip) + ww * optics%omega_wood(band,ip)
         gc       = wl * optics%g_leaf(band,ip)     + ww * optics%g_wood(band,ip)
         beta(i)  = 0.5_wp * (1.0_wp + gc)
         !----- Direct-beam extinction and upscatter, for beam bands only. -----------------!
         if (optics%has_beam(band)) then
            kdir(i)  = kdir_in(i)                              ! band-independent Ross-G extinction
            beta0(i) = 0.5_wp * (1.0_wp + gc / max(kdir(i), tiny_num))
         else
            kdir(i)  = 0.0_wp
            beta0(i) = 0.0_wp
         end if
         !----- Leaf share of absorbed radiation (absorptivity 1-omega, clumping-weighted). -!
         aleaf = (1.0_wp - optics%omega_leaf(band,ip)) * elai(i)
         awood = (1.0_wp - optics%omega_wood(band,ip)) * ewai(i)
         if (aleaf + awood > tiny_num) then
            leaf_frac(i) = aleaf / (aleaf + awood)
         else
            leaf_frac(i) = 1.0_wp
         end if
      end do
   end subroutine blend_cohort_optics

   !=======================================================================================!
   ! Ground optics: the two-stream's lower boundary, bare soil blended toward the snowpack !
   ! (docs/science/canopy_radiation_transfer.md §4).                                       !
   !=======================================================================================!

   !----- Per-band ground reflectance and emission: the ground albedo and 0 for a shortwave band,  !
   !      1 - emissivity and emissivity * sigma * T_soil^4 for a thermal band. --------------------!
   subroutine ground_optics(surf, n_band, has_emission, grnd_refl, grnd_emiss)
      type(ground_optics_env_t), intent(in)  :: surf   !< the ground's albedo, emissivity and temperature
      integer(ik),           intent(in)  :: n_band     !< number of bands
      logical,               intent(in)  :: has_emission(n_band)   !< true for a thermal (emitting) band
      real(wp),              intent(out) :: grnd_refl(n_band)      !< [-] ground reflectance
      real(wp),              intent(out) :: grnd_emiss(n_band)     !< [W/m2] ground emission
      integer(ik) :: b
      do b = 1_ik, n_band
         if (has_emission(b)) then
            grnd_refl(b)  = 1.0_wp - surf%soil_emiss
            grnd_emiss(b) = surf%soil_emiss * stefan_boltzmann * surf%soil_temp ** 4
         else
            grnd_refl(b)  = surf%soil_albedo(b)
            grnd_emiss(b) = 0.0_wp
         end if
      end do
   end subroutine ground_optics

   !=======================================================================================!
   !  Two-stream solver: the single-band adding-method solve over the layer stack.         !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Solve one band over the layers, bottom to top, by the adding method                   !
   ! (docs/science/canopy_radiation_transfer.md §3). Returns each layer's absorbed flux    !
   ! and the beam + downward diffuse at its top, the ground's downward and upward fluxes,  !
   ! and the albedo (upward over incident at the top; 0 when nothing is incident).         !
   !---------------------------------------------------------------------------------------!
   subroutine solve_band(ncoh, etai, omega, beta, beta0, kdir, emission, incid_beam,         &
                         incid_diff, grnd_refl, grnd_emiss, has_beam, has_emission,           &
                         absorbed, incident, dn_ground, up_ground, albedo)
      integer(ik), intent(in)  :: ncoh   !< number of layers (one cohort, or a group of equally tall cohorts)
      real(wp),    intent(in)  :: etai(ncoh), omega(ncoh), beta(ncoh), beta0(ncoh), kdir(ncoh) !< [m2/m2] etai; rest [-]
      real(wp),    intent(in)  :: emission(ncoh)   !< [W/m2] blackbody emission sigma T^4 (thermal band)
      real(wp),    intent(in)  :: incid_beam, incid_diff, grnd_refl, grnd_emiss   !< [W/m2], [W/m2], [-], [W/m2]
      logical,     intent(in)  :: has_beam, has_emission   !< the band carries a direct beam; it carries emission
      real(wp),    intent(out) :: absorbed(ncoh), incident(ncoh)   !< [W/m2] absorbed; beam + down diffuse at the top
      real(wp),    intent(out) :: dn_ground, up_ground, albedo     !< [W/m2], [W/m2], [-]

      real(wp), dimension(ncoh)   :: rdd, tdd, sup, sdn          ! layer R/T + up/down sources
      real(wp), dimension(ncoh+1) :: down0                       ! direct-beam profile (interfaces)
      real(wp), dimension(0:ncoh) :: rbel, sbel                  ! reflectance & up-source of stack below
      real(wp), dimension(ncoh)   :: din, dbot, uup, utop        ! diffuse fluxes per layer
      real(wp) :: beam_up, beam_dn, emit_face, denom, incid_tot
      integer(ik) :: i

      !----- Direct-beam profile: Beer's law from the top down (each layer covers the ground). -!
      if (has_beam) then
         down0(ncoh+1) = incid_beam
         do i = ncoh, 1_ik, -1_ik
            down0(i) = down0(i+1) * exp(-min(kdir(i) * etai(i), lnexp_max))
         end do
      else
         down0(1:ncoh+1) = 0.0_wp                    ! no beam: zero the whole profile, incl. top face
      end if

      !----- Per-layer diffuse reflectance/transmittance and the sources at its faces: the beam it !
      !      scatters, as it leaves the layer, and its emission. ----------------------------------!
      do i = 1_ik, ncoh
         call layer_rt(omega(i), beta(i), etai(i), rdd(i), tdd(i))
         beam_up = 0.0_wp
         beam_dn = 0.0_wp
         if (has_beam) call layer_beam_sources(omega(i), beta(i), beta0(i), kdir(i), etai(i), down0(i+1),  &
                                                beam_up, beam_dn)
         emit_face = 0.0_wp
         if (has_emission) emit_face = (1.0_wp - rdd(i) - tdd(i)) * emission(i)
         sup(i) = beam_up + emit_face                             ! up-source (top of layer)
         sdn(i) = beam_dn + emit_face                             ! down-source (bottom of layer)
      end do

      !----- Bottom-up sweep: reflectance and upward source of {ground + layers 1..i}. The         !
      !      ground's source is its emission plus its reflected direct beam; rbel(0) reflects the  !
      !      diffuse. -----------------------------------------------------------------------------!
      rbel(0) = grnd_refl
      sbel(0) = grnd_emiss + grnd_refl * down0(1)
      do i = 1_ik, ncoh
         denom   = 1.0_wp - rbel(i-1) * rdd(i)
         rbel(i) = rdd(i) + tdd(i) * tdd(i) * rbel(i-1) / denom
         sbel(i) = sup(i) + tdd(i) * (sbel(i-1) + rbel(i-1) * sdn(i)) / denom
      end do

      !----- Top-down sweep: recover diffuse fluxes at every interface. --------------------!
      din(ncoh) = incid_diff
      do i = ncoh, 1_ik, -1_ik
         denom   = 1.0_wp - rdd(i) * rbel(i-1)
         dbot(i) = (tdd(i) * din(i) + sdn(i) + rdd(i) * sbel(i-1)) / denom
         uup(i)  = sbel(i-1) + rbel(i-1) * dbot(i)                ! up flux at bottom of layer i
         utop(i) = rdd(i) * din(i) + tdd(i) * uup(i) + sup(i)     ! up flux at top of layer i
         if (i > 1_ik) din(i-1) = dbot(i)
      end do

      !----- Absorbed per layer = beam divergence + diffuse divergence (in - out); incident = the  !
      !      beam and downward diffuse reaching the layer's top (the light a cohort sees). --------!
      do i = 1_ik, ncoh
         absorbed(i) = (down0(i+1) - down0(i)) + (din(i) + uup(i)) - (utop(i) + dbot(i))
         incident(i) = down0(i+1) + din(i)
      end do

      !----- Ground and top-of-canopy diagnostics. ----------------------------------------!
      dn_ground = dbot(1) + down0(1)                              ! diffuse + beam reaching the ground
      up_ground = uup(1)                                          ! upwelling leaving the ground
      albedo    = 0.0_wp
      incid_tot = incid_beam + incid_diff
      if (incid_tot > tiny_num) albedo = utop(ncoh) / incid_tot
   end subroutine solve_band

   !----- Diffuse reflectance rdd and transmittance tdd of a homogeneous layer, the Stokes         !
   !      formulas (docs/science/canopy_radiation_transfer.md eqs. (8)-(9)). ----------------------!
   pure subroutine layer_rt(omega, beta, etai, rdd, tdd)
      real(wp), intent(in)  :: omega, beta, etai   !< [-] single-scatter albedo, [-] backscatter, [m2/m2] area index
      real(wp), intent(out) :: rdd, tdd            !< [-] diffuse reflectance and transmittance
      real(wp) :: w, sigb, sigf, att, m, rinf, e, e2, denom
      w = min(omega, OMEGA_MAX)                    ! guard the conservative-scattering singularity
      sigb = w * beta                              ! backward diffuse scattering
      sigf = w * (1.0_wp - beta)                   ! forward diffuse scattering
      att  = 1.0_wp - sigf
      if (sigb < tiny_num) then                    ! no backscatter -> no diffuse reflectance
         rdd = 0.0_wp
         tdd = exp(-min(att * etai, lnexp_max))
         return
      end if
      m    = sqrt(max(att * att - sigb * sigb, 0.0_wp))
      rinf = (att - m) / sigb
      e    = exp(-min(m * etai, lnexp_max))
      e2   = e * e
      denom = 1.0_wp - rinf * rinf * e2
      tdd  = (1.0_wp - rinf * rinf) * e / denom
      rdd  = rinf * (1.0_wp - e2) / denom
   end subroutine layer_rt

   !----- The diffuse light a homogeneous layer scatters out of the beam entering its top, as it   !
   !      leaves the layer's two faces (docs/science/canopy_radiation_transfer.md eqs.             !
   !      (10)-(11)): the two-stream equations with the beam's source omega k S0 exp(-k x)         !
   !      through the layer's depth, split up by beta0 and down by 1 - beta0, solved exactly with  !
   !      no diffuse entering either face. Their solution is a particular part in exp(-k x) plus   !
   !      the two diffuse modes of layer_rt. The adding method then combines the layers exactly,   !
   !      so the light does not depend on how a canopy is cut into cohorts. -----------------------!
   pure subroutine layer_beam_sources(omega, beta, beta0, kdir, etai, beam_top, src_up, src_dn)
      real(wp), intent(in)  :: omega, beta, beta0   !< [-] single-scatter albedo, diffuse backscatter, beam upscatter
      real(wp), intent(in)  :: kdir, etai            !< [-] beam extinction per unit area; [m2/m2] area index
      real(wp), intent(in)  :: beam_top              !< [W/m2] beam entering the layer's top
      real(wp), intent(out) :: src_up, src_dn        !< [W/m2] scattered beam leaving the top and the bottom
      real(wp) :: w, sigb, att, m, rinf, k, det, a_up, b_dn, e, ek, c1, c2
      w    = min(omega, OMEGA_MAX)
      sigb = w * beta                                ! backward diffuse scattering
      att  = 1.0_wp - w * (1.0_wp - beta)            ! diffuse attenuation
      m    = sqrt(max(att * att - sigb * sigb, 0.0_wp))
      rinf = 0.0_wp
      if (sigb >= tiny_num) rinf = (att - m) / sigb
      k = kdir
      if (abs(m - k) < RESONANCE_REL * m) k = m * (1.0_wp + RESONANCE_REL)
      !----- The particular solution, up a_up e^(-k x) and down b_dn e^(-k x). --------------------!
      det  = m * m - k * k
      a_up = w * k * beam_top * (beta0 * (att - k) + sigb * (1.0_wp - beta0)) / det
      b_dn = w * k * beam_top * ((1.0_wp - beta0) * (att + k) + sigb * beta0) / det
      !----- The diffuse modes' weights, c1 for e^(-m x) and c2 for e^(m (x - etai)), so that no   !
      !      diffuse enters the top (down = 0 at x = 0) or the bottom (up = 0 at x = etai). -------!
      e  = exp(-min(m * etai, lnexp_max))
      ek = exp(-min(k * etai, lnexp_max))
      c1 = (rinf * e * a_up * ek - b_dn) / (1.0_wp - rinf * rinf * e * e)
      c2 = -c1 * rinf * e - a_up * ek
      src_up = c1 * rinf + c2 * e + a_up
      src_dn = c1 * e + c2 * rinf + b_dn * ek
   end subroutine layer_beam_sources

   !=======================================================================================!
   ! The patch's radiation: every band for one patch, the module's entry point.            !
   !=======================================================================================!

   !----- Solve every band for one patch (docs/science/canopy_radiation_transfer.md §4), cohorts   !
   !      bottom (1) to top (ncoh), and split each cohort's absorbed radiation into leaves and     !
   !      wood. lai and wai are unclumped (clumping is applied here); leaf_temp and wood_temp      !
   !      drive the thermal emission. -------------------------------------------------------------!
   subroutine canopy_radiation(optics, forcing, ncoh, pft, height, lai, wai, leaf_temp, wood_temp, flux)
      type(rad_pft_optics_t), intent(in)  :: optics     !< the per-PFT optics table
      type(rad_env_t),    intent(in)  :: forcing        !< per-band incident radiation and ground boundary
      integer(ik),            intent(in)  :: ncoh       !< number of cohorts
      integer(ik),            intent(in)  :: pft(ncoh)  !< each cohort's PFT
      real(wp),               intent(in)  :: height(ncoh)   !< [m] ascending; exact ties are co-dominant
      real(wp),               intent(in)  :: lai(ncoh), wai(ncoh)               !< [m2/m2]
      real(wp),               intent(in)  :: leaf_temp(ncoh), wood_temp(ncoh)   !< [K]
      type(rad_out_t),       intent(out) :: flux        !< absorbed, incident and ground fluxes, and the albedo

      real(wp), dimension(max(ncoh,1)) :: elai, ewai, etai, omega, beta, beta0, kdir, leaf_frac
      real(wp), dimension(max(ncoh,1)) :: emission, absorbed, kdir_coh
      !----- Each tissue's thermal emission from both faces, and their sum. -----------------------!
      real(wp), dimension(max(ncoh,1)) :: emit_leaf, emit_wood, emit_coh
      !----- Layer arrays: one entry per layer, a group of exactly tied cohorts. ------------------!
      real(wp), dimension(max(ncoh,1)) :: etai_l, omega_l, beta_l, beta0_l, kdir_l, emis_l, abs_l
      real(wp), dimension(max(ncoh,1)) :: inc_l, emit_l
      real(wp), dimension(max(ncoh,1)) :: wabs, wsum_l
      integer(ik) :: lyr(max(ncoh,1)), nmemb(max(ncoh,1))
      integer(ik) :: nlayer, il
      !----- Each PFT's beam extinction coefficient, computed once for its cohorts. ---------------!
      real(wp)    :: kdir_pft(size(optics%lidf, 2))
      logical     :: kdir_done(size(optics%lidf, 2))
      real(wp)    :: gc_l, dn_ground, up_ground, albedo, incid_tot
      integer(ik) :: b, i, ip

      call alloc_rad_out(flux, forcing%n_band, ncoh)

      !----- No resolvable cohorts: all radiation reaches the ground. ---------------------!
      if (ncoh < 1_ik) then
         do b = 1_ik, forcing%n_band
            incid_tot = forcing%incid_beam(b) + forcing%incid_diff(b)
            flux%dn_ground(b) = incid_tot
            flux%up_ground(b) = forcing%grnd_emiss(b) + forcing%grnd_refl(b) * incid_tot
            if (incid_tot > tiny_num) flux%albedo(b) = flux%up_ground(b) / incid_tot
         end do
         return
      end if

      !----- Clumping-corrected leaf & wood area indices (per cohort); total etai. ---------!
      do i = 1_ik, ncoh
         ip      = pft(i)
         elai(i) = optics%clumping_leaf(ip) * lai(i)
         ewai(i) = optics%clumping_wood(ip) * wai(i)
         etai(i) = elai(i) + ewai(i)
      end do

      !----- Cohorts of exactly equal height share one layer (equal_height_layers), as the        !
      !      demography's update_overtopping_lai treats them; a one-member layer is its cohort,   !
      !      bit for bit. ------------------------------------------------------------------------!
      call equal_height_layers(height, ncoh, lyr(1:ncoh), nlayer)
      nmemb(1:nlayer) = 0_ik
      do i = 1_ik, ncoh
         nmemb(lyr(i)) = nmemb(lyr(i)) + 1_ik
      end do

      !----- Ross-G beam extinction depends only on the PFT's leaf angles and the sun, not the     !
      !      band: once per PFT present, shared by its cohorts and every band. --------------------!
      kdir_coh(1:ncoh) = 0.0_wp
      if (any(optics%has_beam(1:forcing%n_band))) then
         kdir_done = .false.
         do i = 1_ik, ncoh
            ip = pft(i)
            if (.not. kdir_done(ip)) then
               kdir_pft(ip)  = gfun_direct(optics%lidf(:, ip), forcing%cosz) / max(forcing%cosz, tiny_num)
               kdir_done(ip) = .true.
            end if
            kdir_coh(i) = kdir_pft(ip)
         end do
      end if

      !----- Loop the configured bands. A band with no light and no emission (VIS and NIR at       !
      !      night) absorbs, reflects and passes nothing: its fluxes keep the zeros alloc_rad_out  !
      !      set, which are what the solve would return. ------------------------------------------!
      do b = 1_ik, forcing%n_band
         if (.not. optics%has_emission(b) .and. forcing%incid_beam(b) == 0.0_wp .and.                    &
             forcing%incid_diff(b) == 0.0_wp .and. forcing%grnd_emiss(b) == 0.0_wp) cycle
         call blend_cohort_optics(optics, b, ncoh, pft, elai(1:ncoh), ewai(1:ncoh), kdir_coh(1:ncoh), &
                                  omega(1:ncoh), beta(1:ncoh), beta0(1:ncoh), kdir(1:ncoh),         &
                                  leaf_frac(1:ncoh))

         !----- A cohort's absorptivity-weighted area (1-omega)*etai, the weight leaf_frac is     !
         !      built from and the layer's members share its absorption by (eq. (14)). -----------!
         wsum_l(1:nlayer) = 0.0_wp
         do i = 1_ik, ncoh
            wabs(i) = (1.0_wp - omega(i)) * etai(i)
            wsum_l(lyr(i)) = wsum_l(lyr(i)) + wabs(i)
         end do

         !----- Thermal emission (§3): each tissue emits from both faces of its area at its       !
         !      absorptivity (Kirchhoff), 2*(1-omega)*AI*sigma*T^4, and the cohort emits as one   !
         !      blackbody at the mean of its tissues' sigma*T^4, weighted as their emission is    !
         !      (leaf_frac). ---------------------------------------------------------------------!
         if (optics%has_emission(b)) then
            do i = 1_ik, ncoh
               ip = pft(i)
               emit_leaf(i) = 2.0_wp * (1.0_wp - optics%omega_leaf(b,ip)) * elai(i)                     &
                              * stefan_boltzmann * leaf_temp(i) ** 4
               emit_wood(i) = 2.0_wp * (1.0_wp - optics%omega_wood(b,ip)) * ewai(i)                     &
                              * stefan_boltzmann * wood_temp(i) ** 4
               emission(i)  = stefan_boltzmann * (leaf_frac(i) * leaf_temp(i) ** 4                         &
                                                  + (1.0_wp - leaf_frac(i)) * wood_temp(i) ** 4)
            end do
         else
            emit_leaf(1:ncoh) = 0.0_wp
            emit_wood(1:ncoh) = 0.0_wp
            emission(1:ncoh)  = 0.0_wp
         end if
         emit_coh(1:ncoh) = emit_leaf(1:ncoh) + emit_wood(1:ncoh)

         !----- Aggregate cohorts onto their layer by area index etai, the weight                  !
         !      blend_cohort_optics mixes leaf and wood by; the emission, like a cohort's, by      !
         !      absorptivity-weighted area. gc is recovered from beta = 0.5*(1+gc) and averaged,   !
         !      which keeps beta0 (gc over kdir) consistent with beta. ----------------------------!
         etai_l(1:nlayer) = 0.0_wp
         omega_l(1:nlayer) = 0.0_wp
         beta_l(1:nlayer) = 0.0_wp
         beta0_l(1:nlayer) = 0.0_wp
         kdir_l(1:nlayer) = 0.0_wp
         emis_l(1:nlayer) = 0.0_wp
         emit_l(1:nlayer) = 0.0_wp
         do i = 1_ik, ncoh
            il = lyr(i)
            etai_l(il)  = etai_l(il)  + etai(i)
            omega_l(il) = omega_l(il) + omega(i) * etai(i)
            beta_l(il)  = beta_l(il)  + (2.0_wp*beta(i) - 1.0_wp) * etai(i)    ! accumulate gc*etai
            kdir_l(il)  = kdir_l(il)  + kdir(i)  * etai(i)
            emis_l(il)  = emis_l(il)  + emission(i) * wabs(i)
            emit_l(il)  = emit_l(il)  + emit_coh(i)
         end do
         do il = 1_ik, nlayer
            if (nmemb(il) == 1_ik) cycle       ! single member: filled verbatim below, bit-identical
            if (etai_l(il) > tiny_num) then
               omega_l(il) = omega_l(il) / etai_l(il)
               gc_l        = beta_l(il)  / etai_l(il)
               kdir_l(il)  = kdir_l(il)  / etai_l(il)
            else
               omega_l(il) = 0.0_wp
               gc_l = 0.0_wp
               kdir_l(il) = 0.0_wp
            end if
            if (wsum_l(il) > tiny_num) then
               emis_l(il) = emis_l(il) / wsum_l(il)
            else
               emis_l(il) = 0.0_wp
            end if
            beta_l(il) = 0.5_wp * (1.0_wp + gc_l)
            if (optics%has_beam(b)) then
               beta0_l(il) = 0.5_wp * (1.0_wp + gc_l / max(kdir_l(il), tiny_num))
            else
               beta0_l(il) = 0.0_wp
            end if
         end do
         !----- One-member layers take the cohort's own values, so a tie-free canopy is unchanged. -!
         do i = 1_ik, ncoh
            if (nmemb(lyr(i)) /= 1_ik) cycle
            il = lyr(i)
            etai_l(il) = etai(i)
            omega_l(il) = omega(i)
            beta_l(il) = beta(i)
            beta0_l(il) = beta0(i)
            kdir_l(il) = kdir(i)
            emis_l(il) = emission(i)
         end do

         call solve_band(nlayer, etai_l(1:nlayer), omega_l(1:nlayer),                          &
                         beta_l(1:nlayer), beta0_l(1:nlayer), kdir_l(1:nlayer), emis_l(1:nlayer), &
                         forcing%incid_beam(b), forcing%incid_diff(b), forcing%grnd_refl(b),    &
                         forcing%grnd_emiss(b), optics%has_beam(b), optics%has_emission(b),      &
                         abs_l(1:nlayer), inc_l(1:nlayer), dn_ground, up_ground, albedo)

         !----- The solve's absorbed flux is net of the layer's emission. Every element of a layer  !
         !      sees the same radiation, so what the layer takes in, that net plus all its members  !
         !      emit, goes to its members by absorptivity-weighted area (eq. (14)), and each is     !
         !      charged its own emission; members that absorb nothing get nothing. -----------------!
         do i = 1_ik, ncoh
            il = lyr(i)
            if (nmemb(il) == 1_ik) then
               absorbed(i) = abs_l(il)
            else if (wsum_l(il) > tiny_num) then
               absorbed(i) = (abs_l(il) + emit_l(il)) * wabs(i) / wsum_l(il) - emit_coh(i)
            else
               absorbed(i) = 0.0_wp
            end if
         end do

         !----- The same within a cohort (eq. (15)): leaves and wood share what it takes in by      !
         !      leaf_frac, and each is charged its own emission. Every member of a layer sees the   !
         !      light reaching the layer's top. ----------------------------------------------------!
         do i = 1_ik, ncoh
            flux%abs_leaf(b,i)  = (absorbed(i) + emit_coh(i)) * leaf_frac(i) - emit_leaf(i)
            flux%abs_wood(b,i)  = (absorbed(i) + emit_coh(i)) * (1.0_wp - leaf_frac(i)) - emit_wood(i)
            flux%incid_top(b,i) = inc_l(lyr(i))
         end do
         flux%albedo(b)    = albedo
         flux%dn_ground(b) = dn_ground
         flux%up_ground(b) = up_ground
      end do
   end subroutine canopy_radiation

   !----- Allocate the per-band forcing for n_band bands, all zero. -------------------------------!
   subroutine alloc_rad_env(f, n_band)
      type(rad_env_t), intent(out) :: f            !< the forcing record
      integer(ik),         intent(in)  :: n_band   !< number of bands
      f%n_band = n_band
      allocate(f%incid_beam(n_band), f%incid_diff(n_band), f%grnd_refl(n_band), f%grnd_emiss(n_band))
      f%incid_beam = 0.0_wp
      f%incid_diff = 0.0_wp
      f%grnd_refl  = 0.0_wp
      f%grnd_emiss = 0.0_wp
   end subroutine alloc_rad_env

   !----- Allocate the returned fluxes for n_band bands and n_coh cohorts, all zero. --------------!
   subroutine alloc_rad_out(flux, n_band, n_coh)
      type(rad_out_t), intent(out) :: flux             !< the returned fluxes
      integer(ik),      intent(in)  :: n_band, n_coh   !< number of bands, number of cohorts
      flux%n_band = n_band
      flux%n_coh  = n_coh
      allocate(flux%abs_leaf(n_band, n_coh), flux%abs_wood(n_band, n_coh), flux%incid_top(n_band, n_coh))
      allocate(flux%albedo(n_band), flux%dn_ground(n_band), flux%up_ground(n_band))
      flux%abs_leaf = 0.0_wp
      flux%abs_wood = 0.0_wp
      flux%incid_top = 0.0_wp
      flux%albedo = 0.0_wp
      flux%dn_ground = 0.0_wp
      flux%up_ground = 0.0_wp
   end subroutine alloc_rad_out

end module ed_canopy_radiation
