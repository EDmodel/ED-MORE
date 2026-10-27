! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_radiation_lib -- the radiation above the canopy: the clearness index, the total
! shortwave split into beam and diffuse PAR and NIR (Erbs et al. 1982; Weiss & Norman
! 1985), the clear-sky emissivity, the longwave synthesis and its lapse, and the beam on
! a slope (docs/science/forcing.md §6, §8 and §11). Stateless laws over plain scalars,
! with their selector codes (SWPART_*, LW_CLEAR_*), which [forcing] takes.
!---~---
module ed_radiation_lib
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : t_3ple, tiny_num, p_std, stefan_boltzmann
   use ed_therm_lib,   only : specific_humidity_to_vapor_pressure
   implicit none
   private

   public :: clearness_index, partition_shortwave, clear_sky_emissivity, synthesize_lwdown
   public :: lapse_longwave, slope_beam_factor
   public :: SWPART_PASSTHROUGH, SWPART_CLEARIDX, SWPART_WEISS_NORMAN, LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO
   public :: COSZ_MIN

   !--- sw_partition: total SW into four streams (docs/science/forcing.md §6).
   integer(ik), parameter :: SWPART_PASSTHROUGH  = 0_ik   !--- file already carries the four streams
   integer(ik), parameter :: SWPART_CLEARIDX     = 1_ik   !--- clearness-index (Erbs) split, eq. (7)
   integer(ik), parameter :: SWPART_WEISS_NORMAN = 2_ik   !--- Weiss-Norman (1985) band-specific split, eqs. (8)-(10)
   !--- The clear-sky emissivity's forms ([forcing].lw_clear_form).
   integer(ik), parameter :: LW_CLEAR_BRUTSAERT = 0_ik   !--- Brutsaert (1975), vapour-pressure power law
   integer(ik), parameter :: LW_CLEAR_IDSO      = 1_ik   !--- Idso & Jackson (1969), temperature only

   real(wp), parameter :: SOLAR_CONSTANT = 1361.0_wp   !--- [W/m2] TOA normal irradiance (mean; eccentricity ~+-3% ignored)
   real(wp), parameter :: COSZ_MIN       = 1.0e-3_wp   !--- [-] cosz floor: below this it is night (SW = 0)
   real(wp), parameter :: FPAR_BEAM      = 0.43_wp     !--- [-] PAR fraction of direct-beam SW energy
   real(wp), parameter :: FPAR_DIFFUSE   = 0.52_wp     !--- [-] PAR fraction of diffuse SW energy (diffuse is PAR-enriched)

   !---~---
   !   The Weiss & Norman (1985) partition's constants, as WN85 publish them. With them the
   ! scheme's PPFD is within 4 % of three towers' quantum sensors (docs/science/forcing.md §6).
   !---~---
   real(wp), parameter :: WN_VIS_TOA = 600.0_wp    !--- [W/m2] visible beam above the atmosphere (WN85 eq. 1)
   real(wp), parameter :: WN_NIR_TOA = 720.0_wp    !--- [W/m2] near-infrared beam above the atmosphere (WN85 eq. 4)
   real(wp), parameter :: WN_W10_TOA = 1320.0_wp   !--- [W/m2] scale of the NIR water-vapour absorption (WN85 eq. 6)
   real(wp), parameter :: WN_PAR_BEAM_EXPEXT = -0.185_wp   !--- visible-beam atmospheric extinction (WN85 eq. 1)
   real(wp), parameter :: WN_NIR_BEAM_EXPEXT = -0.060_wp   !--- NIR-beam atmospheric extinction (WN85 eq. 4)
   real(wp), parameter :: WN_PAR2DIFF_SUN    = 0.400_wp    !--- visible beam->diffuse potential fraction (WN85 eq. 3)
   real(wp), parameter :: WN_NIR2DIFF_SUN    = 0.600_wp    !--- NIR beam->diffuse potential fraction (WN85 eq. 5)
   real(wp), parameter :: WN_W10_A = -1.1950_wp, WN_W10_B = 0.4459_wp, WN_W10_C = -0.0345_wp  !--- NIR water absorption (WN85 eq. 6)
   real(wp), parameter :: WN_PAR_ACT_A = 0.90_wp, WN_PAR_ACT_B = 0.70_wp   !--- visible actual-beam fraction (WN85 eq. 11)
   real(wp), parameter :: WN_NIR_ACT_A = 0.88_wp, WN_NIR_ACT_B = 0.68_wp   !--- NIR actual-beam fraction (WN85 eq. 12)
   real(wp), parameter :: WN_FVIS_DIFF = 0.52_wp, WN_FNIR_DIFF = 0.48_wp   !--- dawn all-diffuse visible/NIR fractions
   real(wp), parameter :: WN_TWOTHIRDS = 2.0_wp / 3.0_wp
   real(wp), parameter :: WN_COSZ_MIN  = 0.017452406437_wp  !--- cos(89 deg): below this the sun is grazing -> all diffuse

contains

   !---~---
   !   SHORTWAVE PARTITION: total SWdown -> the four (beam/diffuse) x (PAR/NIR) streams
   ! (docs/science/forcing.md §6). SWPART_CLEARIDX (default): the Erbs et al. (1982) diffuse
   ! fraction of the clearness index kt = SW/(S0*cosz) (eq. (7)), then fixed PAR fractions.
   ! SWPART_WEISS_NORMAN: the band-specific, pressure-dependent Weiss & Norman (1985) split
   ! (eqs. (8)-(10)). Both conserve energy exactly: the four streams sum to sw.
   !---~---
   pure subroutine partition_shortwave(swdown_total, cosz, psurf_pa, scheme,                   &
                                       par_beam, par_diffuse, nir_beam, nir_diffuse)
      real(wp),    intent(in)  :: swdown_total, cosz, psurf_pa
      integer(ik), intent(in)  :: scheme
      real(wp),    intent(out) :: par_beam, par_diffuse, nir_beam, nir_diffuse
      real(wp) :: sw, i_toa, kt, fdiff, beam, diff

      sw = max(swdown_total, 0.0_wp)
      if (sw <= tiny_num .or. cosz <= COSZ_MIN) then          ! night / no light
         par_beam = 0.0_wp
         par_diffuse = 0.0_wp
         nir_beam = 0.0_wp
         nir_diffuse = 0.0_wp
         return
      end if

      select case (scheme)
      case (SWPART_WEISS_NORMAN)
         call weiss_norman_partition(sw, cosz, psurf_pa, par_beam, par_diffuse, nir_beam, nir_diffuse)
      case default   ! SWPART_CLEARIDX (and anything else -- the reader rejects unknown spellings)
         i_toa = SOLAR_CONSTANT * cosz                        ! TOA irradiance on a horizontal surface
         kt    = min(1.0_wp, max(0.0_wp, sw / max(i_toa, tiny_num)))
         fdiff = erbs_diffuse_fraction(kt)
         diff = fdiff * sw
         beam = sw - diff
         par_beam    = FPAR_BEAM    * beam
         nir_beam    = beam - par_beam
         par_diffuse = FPAR_DIFFUSE * diff
         nir_diffuse = diff - par_diffuse
      end select
   end subroutine partition_shortwave

   !---~---
   !   WEISS & NORMAN (1985) band-specific partition (docs/science/forcing.md eqs. (8)-(10)):
   ! clear-sky beam and diffuse per band from the air-mass optical depth, scaled to the observed
   ! total, whose ratio to the clear-sky total sets each band's beam fraction (cloudier -> more
   ! diffuse). The four streams sum to sw exactly. A grazing sun (cosz <= WN_COSZ_MIN) is all
   ! diffuse.
   !---~---
   pure subroutine weiss_norman_partition(sw, cosz, psurf_pa, par_beam, par_diffuse, nir_beam, nir_diffuse)
      real(wp), intent(in)  :: sw, cosz, psurf_pa
      real(wp), intent(out) :: par_beam, par_diffuse, nir_beam, nir_diffuse
      real(wp) :: secz, log10secz, prat, w10
      real(wp) :: par_beam_pot, par_diff_pot, par_full_pot
      real(wp) :: nir_beam_pot, nir_diff_pot, nir_full_pot, ratio, par_full, nir_full
      real(wp) :: aux_par, aux_nir, fvis_beam, fnir_beam

      if (cosz <= WN_COSZ_MIN) then                           ! grazing sun: secant unstable -> all diffuse
         par_beam = 0.0_wp
         nir_beam = 0.0_wp
         par_diffuse = WN_FVIS_DIFF * sw
         nir_diffuse = WN_FNIR_DIFF * sw   ! 0.52+0.48 = 1 (conserves sw)
         return
      end if

      secz      = 1.0_wp / cosz
      log10secz = log10(secz)
      prat      = psurf_pa / p_std

      !---~---
      !   Clear-sky visible beam and diffuse (WN85 eqs. 1, 3, 9; docs/science/forcing.md
      ! eq. (8)). The diffuse term scales the beam it subtracts by cosz once more than other
      ! implementations of WN85 eq. 3 do (docs/science/forcing.md §6).
      !---~---
      par_beam_pot = WN_VIS_TOA * exp(WN_PAR_BEAM_EXPEXT * prat * secz) * cosz
      par_diff_pot = WN_PAR2DIFF_SUN * (WN_VIS_TOA - par_beam_pot) * cosz
      par_full_pot = par_beam_pot + par_diff_pot

      !---~---
      !   Clear-sky NIR: the water-vapour absorption w10, then beam and diffuse (WN85 eqs. 6,
      ! 4, 5, 10; docs/science/forcing.md eq. (9)).
      !---~---
      w10          = WN_W10_TOA * 10.0_wp ** (WN_W10_A + log10secz * (WN_W10_B + WN_W10_C * log10secz))
      nir_beam_pot = (WN_NIR_TOA * exp(WN_NIR_BEAM_EXPEXT * prat * secz) - w10) * cosz
      nir_diff_pot = WN_NIR2DIFF_SUN * (WN_NIR_TOA - nir_beam_pot - w10) * cosz
      nir_full_pot = nir_beam_pot + nir_diff_pot

      !---~---
      !   Scale the clear-sky streams to the observed total (WN85 eqs. 7-8). Just above the
      ! grazing cutoff (cosz ~ 0.0175-0.029) w10 exceeds the attenuated NIR beam, so nir_diffuse
      ! can go slightly negative, and PAR exceed sw, while the four streams still sum to sw
      ! exactly; the shortwave is tiny there (docs/science/forcing.md §6).
      !---~---
      ratio    = sw / max(par_full_pot + nir_full_pot, tiny_num)
      par_full = ratio * par_full_pot
      nir_full = ratio * nir_full_pot

      !---~---
      !   The beam fraction of each band, clipped to [0,1] (WN85 eqs. 11-12;
      ! docs/science/forcing.md eq. (10)).
      !---~---
      aux_par   = min(WN_PAR_ACT_A, max(0.0_wp, ratio))
      aux_nir   = min(WN_NIR_ACT_A, max(0.0_wp, ratio))
      fvis_beam = min(1.0_wp, max(0.0_wp,                                                       &
                  par_beam_pot * (1.0_wp - ((WN_PAR_ACT_A - aux_par) / WN_PAR_ACT_B) ** WN_TWOTHIRDS) &
                  / max(par_full_pot, tiny_num)))
      fnir_beam = min(1.0_wp, max(0.0_wp,                                                       &
                  nir_beam_pot * (1.0_wp - ((WN_NIR_ACT_A - aux_nir) / WN_NIR_ACT_B) ** WN_TWOTHIRDS) &
                  / max(nir_full_pot, tiny_num)))

      par_beam    = fvis_beam * par_full
      par_diffuse = (1.0_wp - fvis_beam) * par_full
      nir_beam    = fnir_beam * nir_full
      nir_diffuse = (1.0_wp - fnir_beam) * nir_full
   end subroutine weiss_norman_partition

   !---~---
   !   LONGWAVE SYNTHESIS, for a source without downwelling longwave (docs/science/forcing.md
   ! §11): LWdown = eps_clear * sigma * Tair^4 * (1 + a*(1 - kt)), eq. (17).
   !---~---

   !---~---
   !   Clearness index kt = SW / (S0 cosz) in [0,1], the one index the shortwave partition and
   ! the longwave synthesis share. Negative at night, where kt is undefined: a zero would
   ! read as overcast (docs/science/forcing.md §11).
   !---~---
   pure function clearness_index(swdown_total, cosz) result(kt)
      real(wp), intent(in) :: swdown_total, cosz
      real(wp) :: kt, i_toa
      if (cosz <= COSZ_MIN) then
         kt = -1.0_wp                                   ! night: undefined, not zero
         return
      end if
      i_toa = SOLAR_CONSTANT * cosz
      kt    = min(1.0_wp, max(0.0_wp, max(swdown_total, 0.0_wp) / max(i_toa, tiny_num)))
   end function clearness_index

   !---~---
   !   Clear-sky emissivity: Brutsaert (1975), a power law in the screen-level vapour
   ! pressure over temperature, or Idso & Jackson (1969), temperature alone. The vapour
   ! pressure comes from q; a saturation curve at Tair would assume saturation.
   !---~---
   pure function clear_sky_emissivity(form, tair_k, qair, psurf_pa) result(eps)
      integer(ik), intent(in) :: form                   !--- LW_CLEAR_BRUTSAERT | LW_CLEAR_IDSO
      real(wp),    intent(in) :: tair_k, qair, psurf_pa
      real(wp) :: eps, e_hpa, dt
      select case (form)
      case (LW_CLEAR_IDSO)
         dt  = t_3ple - tair_k
         eps = 1.0_wp - 0.261_wp * exp(-7.77e-4_wp * dt * dt)
      case default                                      ! LW_CLEAR_BRUTSAERT
         e_hpa = 0.01_wp * specific_humidity_to_vapor_pressure(max(qair, 0.0_wp), psurf_pa)   ! [hPa]
         eps   = 1.24_wp * (max(e_hpa, tiny_num) / max(tair_k, tiny_num)) ** (1.0_wp / 7.0_wp)
      end select
      eps = min(1.0_wp, max(0.5_wp, eps))   ! clamp-ok: the formulas give >~0.6; 0.5 only guards malformed input
   end function clear_sky_emissivity

   !---~---
   !   The synthesis (docs/science/forcing.md eq. (17)); at night the caller passes the last
   ! daytime kt.
   !---~---
   pure function synthesize_lwdown(form, tair_k, qair, psurf_pa, kt, cloud_a) result(lw)
      integer(ik), intent(in) :: form
      real(wp),    intent(in) :: tair_k, qair, psurf_pa
      real(wp),    intent(in) :: kt                     !--- [-] clearness index in [0,1]
      real(wp),    intent(in) :: cloud_a                !--- [-] cloud-correction coefficient
      real(wp) :: lw, eps, kt_use
      eps    = clear_sky_emissivity(form, tair_k, qair, psurf_pa)
      kt_use = min(1.0_wp, max(0.0_wp, kt))
      lw     = eps * stefan_boltzmann * tair_k ** 4 * (1.0_wp + max(cloud_a, 0.0_wp) * (1.0_wp - kt_use))
   end function synthesize_lwdown

   !--- Erbs et al. (1982) diffuse fraction of hourly global radiation vs clearness index.
   pure function erbs_diffuse_fraction(kt) result(fd)
      real(wp), intent(in) :: kt
      real(wp) :: fd
      if (kt <= 0.22_wp) then
         fd = 1.0_wp - 0.09_wp * kt
      else if (kt <= 0.80_wp) then
         fd = 0.9511_wp - 0.1604_wp*kt + 4.388_wp*kt**2 - 16.638_wp*kt**3 + 12.336_wp*kt**4
      else
         fd = 0.165_wp
      end if
      fd = min(1.0_wp, max(0.0_wp, fd))
   end function erbs_diffuse_fraction

   !---~---
   !   File longwave at the site: scaled by the clear-sky emission eps*T^4 at the site over the
   ! cell (NLDAS; docs/science/forcing.md eq. (14)), with the same clear-sky emissivity the
   ! longwave synthesis uses.
   !---~---
   pure function lapse_longwave(lw_grid, form, t_grid, q_grid, p_grid, t_site, q_site, p_site) result(lw_site)
      real(wp),    intent(in) :: lw_grid                   !--- [W/m2] downward longwave at the cell
      integer(ik), intent(in) :: form                      !--- LW_CLEAR_BRUTSAERT | LW_CLEAR_IDSO
      real(wp),    intent(in) :: t_grid, q_grid, p_grid, t_site, q_site, p_site
      real(wp) :: lw_site
      lw_site = lw_grid * (clear_sky_emissivity(form, t_site, q_site, p_site) * t_site**4)      &
                        / (clear_sky_emissivity(form, t_grid, q_grid, p_grid) * t_grid**4)
   end function lapse_longwave

   !---~---
   !   The factor the direct beam takes on a slope: the beam normal to the sun, beam / cos z,
   ! times cos i, over the beam on level ground. The beam normal is held at most the solar
   ! constant, which keeps a low sun's small cos z from inflating the beam a slope facing it
   ! receives. 1 with no beam to scale.
   !---~---
   pure real(wp) function slope_beam_factor(beam, cosz, cosi) result(f)
      real(wp), intent(in) :: beam   !--- [W/m2] direct beam on level ground (all bands)
      real(wp), intent(in) :: cosz   !--- [-] cosine of the solar zenith angle
      real(wp), intent(in) :: cosi   !--- [-] cosine of the incidence on the slope
      f = 1.0_wp
      if (beam <= 0.0_wp .or. cosz <= COSZ_MIN) return
      f = min(cosi / cosz, SOLAR_CONSTANT * cosi / beam)
   end function slope_beam_factor

end module ed_radiation_lib
