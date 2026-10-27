! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_allometry_lib -- a stem's size from its diameter: height, aboveground biomass
! (AGB), leaf area, wood area and sapwood; the way back from carbon to diameter; and the
! carbon-pool targets the allocator fills. The curves, their shipped coefficients and the
! inversions are docs/science/plant_carbon_allocation.md §2 (The size allometry).
!
!   Each PFT has its own record, allometry_t, the first argument of every routine: its
! height curve and cap, its biomass, leaf-area, wood-area and sapwood allometries, and the
! share of its wood above ground. The records are the PFT table's (pft%allom(pft)); each
! polygon holds its own copy, whose chave_e is its cell's (docs/configuration.md §Regional
! runs). Biomass and leaf area are in carbon.
!---~---
module ed_allometry_lib
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : tiny_num, pio4
   implicit none
   private

   public :: allometry_t, HEIGHT_POWER, HEIGHT_CHAVE
   public :: dbh_to_height, height_to_dbh
   public :: dbh_to_agb, agb_to_dbh, min_cohort_carbon
   public :: dbh_to_leaf_area, dbh_to_wai, sapwood_fraction
   public :: size_to_leaf_carbon, size_to_wood_carbon, wood_to_dbh, carbon_to_structure, storage_carbon_target

   !---~---
   !   Height curves: a power law, or the log-quadratic of Chave et al. (2014) eq. 6a,
   ! whose coefficients a PFT may fit to its own site.
   !---~---
   integer(ik), parameter :: HEIGHT_POWER = 1_ik, HEIGHT_CHAVE = 2_ik

   !---~---
   !   One PFT's allometry, from the PFT file's [allometry] block. The defaults are the
   ! shipped pan-tropical values, so a record made without a config (a unit test, a
   ! Python call) is usable.
   !---~---
   type :: allometry_t
      !--- Height: the curve, and the cap on every curve.
      integer(ik) :: height_allometry = HEIGHT_POWER   !--- the height curve: power law or log-quadratic
      real(wp)    :: hgt_max     = 46.0_wp        !--- [m]   the largest height
      real(wp)    :: b1Ht        = 1.139963_wp    !--- [--]  power law: intercept of ln H on ln dbh
      real(wp)    :: b2Ht        = 0.564899_wp    !--- [--]  power law: slope of ln H on ln dbh
      real(wp)    :: chave_a     = 0.893_wp       !--- [--]  log-quadratic: intercept (Chave et al. 2014 eq. 6a)
      real(wp)    :: chave_b     = 0.760_wp       !--- [--]  log-quadratic: slope on ln dbh
      real(wp)    :: chave_c     = -0.0340_wp     !--- [--]  log-quadratic: coefficient of (ln dbh)^2
      real(wp)    :: chave_rse   = 0.243_wp       !--- [--]  log-quadratic: residual standard error of the fit
      real(wp)    :: chave_e     = 0.0_wp         !--- [--]  the site's stress index E: the larger, the shorter
      !--- Biomass, leaf area, wood area and sapwood.
      real(wp)    :: agb_c1      = 0.03365_wp     !--- [kgC] AGB scale: Chave et al. (2014) eq. 4, halved to carbon
      real(wp)    :: agb_c2      = 0.976_wp       !--- [--]  AGB exponent on rho and on dbh^2*h
      real(wp)    :: lai_b1      = 0.23384770_wp  !--- [--]  per-stem leaf-area scale
      real(wp)    :: lai_b2      = 0.6410495_wp   !--- [--]  per-stem leaf-area exponent
      real(wp)    :: wai_b1      = 0.0096_wp      !--- [--]  per-stem wood-area scale
      real(wp)    :: wai_b2      = 2.0947_wp      !--- [--]  per-stem wood-area exponent
      real(wp)    :: sapwood_area_b1 = 1.582_wp   !--- [cm2/cm^b2] sapwood-area scale
      real(wp)    :: sapwood_area_b2 = 1.764_wp   !--- [--]  sapwood-area exponent
      real(wp)    :: aboveground_frac = 0.7_wp    !--- [--]  share of the wood carbon above ground
   end type allometry_t

contains

   !--- Height: the PFT's curve, capped at hgt_max, and its inverse.

   !--- Diameter [cm] -> height [m] on the PFT's curve, capped at hgt_max.
   elemental pure function dbh_to_height(allom, dbh) result(h)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh
      real(wp)                      :: h, x
      select case (allom%height_allometry)
      case (HEIGHT_CHAVE)
         x  = log(max(dbh, tiny_num))
         h  = exp(chave_intercept(allom) + allom%chave_b * x + allom%chave_c * x * x)
      case default
         h  = exp(allom%b1Ht + allom%b2Ht * log(max(dbh, tiny_num)))
      end select
      h = min(h, allom%hgt_max)
   end function dbh_to_height

   !--- Height [m] -> diameter [cm], the inverse of the uncapped curve.
   elemental pure function height_to_dbh(allom, h) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: h
      real(wp)                      :: dbh
      select case (allom%height_allometry)
      case (HEIGHT_CHAVE)
         dbh = exp(chave_small_root(allom%chave_c, allom%chave_b, chave_intercept(allom) - log(max(h, tiny_num))))
      case default
         dbh = exp((log(max(h, tiny_num)) - allom%b1Ht) / allom%b2Ht)
      end select
   end function height_to_dbh

   !---~---
   !   The intercept of the log-quadratic curve, ln H = intercept + chave_b x + chave_c
   ! x^2 with x = ln(dbh): the fit's, less the site's E, plus RSE^2/2, which turns
   ! the fit's median into its mean.
   !---~---
   pure real(wp) function chave_intercept(allom)
      type(allometry_t), intent(in) :: allom
      chave_intercept = allom%chave_a - allom%chave_e + 0.5_wp * allom%chave_rse * allom%chave_rse
   end function chave_intercept

   !---~---
   !   The smaller root x of a x^2 + b x + c = 0: the curve's rising branch, which for
   ! Chave's coefficients holds up to dbh = exp(-b / (2a)), about 70,000 cm. Written
   ! as 2(-c)/(b + sqrt(.)) so it keeps its digits when c is small.
   !---~---
   pure real(wp) function chave_small_root(a, b, c) result(x)
      real(wp), intent(in) :: a, b, c
      x = -2.0_wp * c / (b + sqrt(max(b * b - 4.0_wp * a * c, 0.0_wp)))
   end function chave_small_root

   !--- Aboveground biomass (Chave et al. 2014 eq. 4, in carbon) and its inverse.

   !--- Aboveground biomass [kgC/plant], Chave et al. (2014) eq. 4 in carbon.
   elemental pure function dbh_to_agb(allom, dbh, h, rho) result(agb)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, rho
      real(wp)                      :: agb
      agb = allom%agb_c1 * rho ** allom%agb_c2 * (dbh * dbh * h) ** allom%agb_c2
   end function dbh_to_agb

   !---~---
   !   AGB -> diameter [cm], which recovers a diameter from the conserved carbon. Under
   ! the power law agb = K dbh^P on both sides of the height cap: solve the uncapped
   ! form, and redo it in the capped form if that stem is taller than hgt_max.
   !---~---
   elemental pure function agb_to_dbh(allom, agb, rho) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: agb, rho
      real(wp)                      :: dbh, k_un, p_un, k_cap, p_cap, h
      if (allom%height_allometry == HEIGHT_CHAVE) then
         dbh = agb_to_dbh_chave(allom, agb, rho)
         return
      end if
      k_un = allom%agb_c1 * rho ** allom%agb_c2 * exp(allom%agb_c2 * allom%b1Ht)
      p_un = allom%agb_c2 * (2.0_wp + allom%b2Ht)
      dbh  = (max(agb, tiny_num) / k_un) ** (1.0_wp / p_un)
      h    = exp(allom%b1Ht + allom%b2Ht * log(max(dbh, tiny_num)))
      if (h > allom%hgt_max) then
         k_cap = allom%agb_c1 * rho ** allom%agb_c2 * allom%hgt_max ** allom%agb_c2
         p_cap = 2.0_wp * allom%agb_c2
         dbh   = (max(agb, tiny_num) / k_cap) ** (1.0_wp / p_cap)
      end if
   end function agb_to_dbh

   !---~---
   !   AGB -> diameter under Chave's height, in closed form: below the cap, the biomass
   ! law 2x + ln h(x) = target is quadratic in x = ln(dbh); x is its smaller root. A
   ! stem that root makes taller than hgt_max has 2x + ln(hgt_max) = target.
   !---~---
   elemental pure function agb_to_dbh_chave(allom, agb, rho) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: agb, rho
      real(wp)                      :: dbh, target, x
      target = log(max(agb, tiny_num) / (allom%agb_c1 * rho ** allom%agb_c2)) / allom%agb_c2
      x = chave_small_root(allom%chave_c, 2.0_wp + allom%chave_b, chave_intercept(allom) - target)
      if (chave_intercept(allom) + allom%chave_b * x + allom%chave_c * x * x > log(allom%hgt_max))     &
         x = 0.5_wp * (target - log(allom%hgt_max))
      dbh = exp(x)
   end function agb_to_dbh_chave

   !---~---
   !   The AGB [kgC/plant] of a recruit at min_cohort_height: the carbon per recruit
   ! that turns reproduction carbon into recruit numbers.
   !---~---
   elemental pure function min_cohort_carbon(allom, min_cohort_height, rho) result(carbon_min)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: min_cohort_height, rho
      real(wp)                      :: carbon_min
      carbon_min = dbh_to_agb(allom, height_to_dbh(allom, min_cohort_height), min_cohort_height, rho)
   end function min_cohort_carbon

   !--- Leaf area, wood area and sapwood.

   !--- One-sided leaf area of one stem [m2]; the cohort's LAI is nplant times it.
   elemental pure function dbh_to_leaf_area(allom, dbh, h) result(leaf_area)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h
      real(wp)                      :: leaf_area
      leaf_area = allom%lai_b1 * (dbh * dbh * h) ** allom%lai_b2
   end function dbh_to_leaf_area

   !---~---
   !   Wood area index [m2/m2], nplant wai_b1 dbh^wai_b2. It follows the stem, not the
   ! leaf area, because the wood's heat store and exchanges scale with the stem
   ! (docs/science/vegetation_energy_dynamics.md, below eq. (8)).
   !---~---
   elemental pure function dbh_to_wai(allom, dbh, nplant) result(wai)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh      !--- [cm]    diameter at breast height
      real(wp),          intent(in) :: nplant   !--- [pl/m2] plant density
      real(wp)                      :: wai      !--- [m2/m2] wood area index
      wai = nplant * allom%wai_b1 * max(dbh, tiny_num) ** allom%wai_b2
   end function dbh_to_wai

   !---~---
   !   Sapwood share of basal area, sapwood_area_b1 dbh^sapwood_area_b2/(pi/4 dbh^2)
   ! (Xu 2018), capped at one since a small stem is all sapwood
   ! (docs/science/vegetation_energy_dynamics.md, below eq. (8)).
   !---~---
   elemental pure function sapwood_fraction(allom, dbh) result(f_sap)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh      !--- [cm] diameter at breast height
      real(wp)                      :: f_sap    !--- [--] sapwood area / basal area, in (0,1]
      real(wp) :: d, sapw_area, basal_area
      d          = max(dbh, tiny_num)
      sapw_area  = allom%sapwood_area_b1 * d ** allom%sapwood_area_b2
      basal_area = pio4 * d * d
      f_sap      = min(1.0_wp, sapw_area / max(basal_area, tiny_num))
   end function sapwood_fraction

   !--- The carbon pools: their targets, and the geometry they set.

   !---~---
   !   The leaf carbon target [kgC/plant]: the allometric leaf area over SLA, so a
   ! cohort at its target shows exactly that leaf area
   ! (docs/science/plant_carbon_allocation.md §2).
   !---~---
   elemental pure function size_to_leaf_carbon(allom, dbh, h, sla) result(leaf_carbon)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, sla        !--- [cm],[m],[m2/kgC]
      real(wp)                      :: leaf_carbon         !--- [kgC/plant]
      leaf_carbon = dbh_to_leaf_area(allom, dbh, h) / max(sla, tiny_num)
   end function size_to_leaf_carbon

   !---~---
   !   A stem's total wood carbon [kgC/plant], coarse roots included: its AGB over the
   ! aboveground fraction, so the aboveground share is exactly the AGB.
   !---~---
   elemental pure function size_to_wood_carbon(allom, dbh, h, rho) result(wood_carbon)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, rho
      real(wp)                      :: wood_carbon         !--- [kgC/plant] total woody (above + belowground)
      wood_carbon = dbh_to_agb(allom, dbh, h, rho) / max(allom%aboveground_frac, tiny_num)
   end function size_to_wood_carbon

   !---~---
   !   The inverse of size_to_wood_carbon: wood carbon is the prognostic size, and the
   ! diameter follows from it.
   !---~---
   elemental pure function wood_to_dbh(allom, wood_carbon, rho) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: wood_carbon, rho
      real(wp)                      :: dbh                 !--- [cm]
      dbh = agb_to_dbh(allom, wood_carbon * allom%aboveground_frac, rho)
   end function wood_to_dbh

   !---~---
   !   A cohort's geometry from its carbon pools: the diameter from wood carbon, then
   ! the height, basal area and AGB, and the leaf area from leaf carbon. It is the
   ! body of set_cohort_size_from_carbon, the one writer of a cohort's geometry.
   !---~---
   pure subroutine carbon_to_structure(allom, wood_carbon, leaf_carbon, rho, sla, dbh, height, basal_area,  &
                                       agb, leaf_area)
      type(allometry_t), intent(in)  :: allom
      real(wp),          intent(in)  :: wood_carbon, leaf_carbon    !--- [kgC/plant] prognostic pools
      real(wp),          intent(in)  :: rho, sla                    !--- [g/cm3], [m2/kgC]
      real(wp),          intent(out) :: dbh, height, basal_area, agb, leaf_area
      dbh        = wood_to_dbh(allom, wood_carbon, rho)
      height     = dbh_to_height(allom, dbh)
      basal_area = pio4 * dbh * dbh
      agb        = allom%aboveground_frac * wood_carbon
      leaf_area  = leaf_carbon * sla
   end subroutine carbon_to_structure

   !---~---
   !   A cohort's storage target: storage_cushion times the carbon one full flush
   ! rebuilds, its leaves plus its fine roots. It is also the starting storage of a
   ! cohort set on allometry.
   !---~---
   elemental pure function storage_carbon_target(storage_cushion, leaf_carbon, fineroot_carbon) result(storage)
      real(wp), intent(in) :: storage_cushion              !--- [--]
      real(wp), intent(in) :: leaf_carbon, fineroot_carbon !--- [kgC/plant]
      real(wp)             :: storage                      !--- [kgC/plant]
      storage = storage_cushion * (leaf_carbon + fineroot_carbon)
   end function storage_carbon_target

end module ed_allometry_lib
