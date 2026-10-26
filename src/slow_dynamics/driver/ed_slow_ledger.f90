! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_slow_ledger -- the site conservation ledger across a slow step                        !
! (docs/science/conservation.md). At each operator boundary it snapshots the site's        !
! carbon, water and energy store and attributes what the operator did not declare to its   !
! phase (§3, eq. (1)); a breach of eq. (2) is reported at the end of the run (§5). It      !
! measures; it corrects nothing. Routines: slow_ledger_open, _rebase, _declare, _mark,     !
! _report; slow_ledger_heat_open/_close, an operator's thermal-mass exchange (§4);         !
! slow_site_store (§2); the declared amounts slow_fast_carbon_handover, slow_co2_handoff   !
! and litter_carbon_total (§3); check_plausible (§6). slow_ledger_soil_check keeps the     !
! run's worst soil-carbon checks.                                                          !
!==========================================================================================!
module ed_slow_ledger
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_cohort_state_vars,    only : cohort_film_water, cohort_carbon_supply
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area
   use ed_patch_state_vars,     only : cas_capacities, soil_water_store, soil_energy_store, litter_input_t
   use ed_budget_check,         only : BUDGET_RTOL_FLUX
   use ed_soil_carbon,          only : soil_carbon_bad_pool, soil_carbon_out_t
   use ed_constants,            only : rho_air_ref, umol_2_kgc
   use ed_time,                 only : ed_time_t,   time_to_string
   use ed_biogeochem_params,    only : N_SOIL_POOL, SOIL_POOL_NAME
   implicit none
   private

   public :: slow_store_t, slow_ledger_t, slow_site_store, slow_fast_carbon_handover
   public :: slow_tissue_energy, slow_co2_handoff, litter_carbon_total
   public :: slow_ledger_open, slow_ledger_rebase, slow_ledger_declare, slow_ledger_mark, slow_ledger_report
   public :: slow_ledger_heat_open, slow_ledger_heat_close, slow_ledger_soil_check, slow_ledger_failures
   public :: SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW, SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT,          &
             SLOW_PHASE_DISTURB, SLOW_PHASE_PATCH, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC,              &
             SLOW_PHASE_TISSUE_WATER

   !----- The phases of a slow step (conservation.md §3), each marked at the end of the    !
   !      operators that own it, so a residual names the operator that produced it. -------!
   integer(ik), parameter :: SLOW_PHASE_ALLOCATE = 1_ik  !< phenology, traits, carbon allocation, turnover shed
   integer(ik), parameter :: SLOW_PHASE_GROW     = 2_ik  !< tendency build + state commit (growth and mortality)
   integer(ik), parameter :: SLOW_PHASE_RECRUIT  = 3_ik  !< apply_recruitment
   integer(ik), parameter :: SLOW_PHASE_COHORT   = 4_ik  !< cohort fuse / terminate / split
   integer(ik), parameter :: SLOW_PHASE_DISTURB  = 5_ik  !< apply_patch_disturbance
   integer(ik), parameter :: SLOW_PHASE_PATCH    = 6_ik  !< patch sort / fuse / terminate (+ trailing cohort ops)
   integer(ik), parameter :: SLOW_PHASE_CANOPY   = 7_ik  !< refresh_canopy_depth (the open CAS control volume)
   integer(ik), parameter :: SLOW_PHASE_SOILC    = 8_ik  !< the CENTURY daily step
   integer(ik), parameter :: SLOW_PHASE_TISSUE_WATER = 9_ik  !< tissue water above capacity shed, before the fast loop
   integer(ik), parameter :: N_SLOW_PHASE        = 9_ik

   !----- The tolerance, eq. (2) of conservation.md §5: an absolute floor per mark, plus   !
   !      the relative test on the declared flux that the fast budgets use. ---------------!
   real(wp), parameter :: SLOW_ATOL_CARBON = 1.0e-12_wp   !< [kgC/m2] per mark
   real(wp), parameter :: SLOW_ATOL_WATER  = 1.0e-9_wp    !< [kg/m2]  per mark
   real(wp), parameter :: SLOW_ATOL_ENERGY = 1.0e-3_wp    !< [J/m2]   per mark
   !----- ...plus a round-off allowance on the store (§5): a phase that declares nothing   !
   !      still rearranges a large store. -------------------------------------------------!
   real(wp), parameter :: SLOW_ROUNDOFF_FRAC = 1.0e-11_wp !< [-] of the store, per mark

   !----- One total per currency, per m2 of site, area-weighted over the patches (§2). ----!
   type :: slow_store_t
      real(wp) :: carbon = 0.0_wp   !< [kgC/m2] live pools + necromass + canopy air + recruit pool
      real(wp) :: water  = 0.0_wp   !< [kg/m2]  soil + pond + snow + canopy vapour + tissue + film
      real(wp) :: energy = 0.0_wp   !< [J/m2]   soil + pond + snow + canopy air + tissue heat
   end type slow_store_t

   !----- The run's accumulator, held in the polygon's run record (polygon_run_t). --------!
   type :: slow_ledger_t
      logical             :: active = .false.       !< off => every entry point is a no-op
      !----- Rolling state within the current step. ---------------------------------------!
      type(slow_store_t)  :: at_mark                !< store when the current phase opened
      type(slow_store_t)  :: decl_in, decl_out      !< declared boundary terms for the current phase
      !----- Per-phase sums over the run; the signed sum shows a one-signed drift that a   !
      !      running maximum misses (§5). -------------------------------------------------!
      type(slow_store_t)  :: resid_sum(N_SLOW_PHASE)
      type(slow_store_t)  :: resid_abs(N_SLOW_PHASE)
      type(slow_store_t)  :: resid_worst(N_SLOW_PHASE)
      type(slow_store_t)  :: decl_gross(N_SLOW_PHASE)
      integer(ik)         :: n_mark(N_SLOW_PHASE) = 0_ik
      integer(ik)         :: n_step = 0_ik
      !----- Whole-run endpoints, for the drift line of the report. -----------------------!
      type(slow_store_t)  :: store_first, store_last
      logical             :: has_first = .false.
      !----- Plausibility (§6): the first phase that left a CENTURY pool that cannot       !
      !      exist, and where; everything after it is downstream of that state. -----------!
      integer(ik)         :: bad_phase = 0_ik      !< phase that first left an impossible store
      integer(ik)         :: bad_patch = 0_ik      !< which patch
      integer(ik)         :: bad_pool  = 0_ik      !< which CENTURY pool (index)
      real(wp)            :: bad_value = 0.0_wp    !< its value
      integer(ik)         :: n_bad     = 0_ik      !< how many marks saw an impossible store
      !----- A thermal-mass bracket's opening tissue energy [J/m2] (slow_ledger_heat_open). --!
      real(wp)            :: heat_at_open = 0.0_wp
      !----- The run's worst soil-carbon checks (soil_carbon.md §7, §9), kept whether or not  !
      !      the ledger is active, with the step and patch count where the Rh gap peaked. --!
      integer(ik)         :: n_seam       = 0_ik     !< soil-carbon steps checked
      real(wp)            :: worst_rh_gap = 0.0_wp   !< [kgC/m2/day] max |the step's Rh - the fast loop's|
      character(len=19)   :: worst_rh_gap_when   = ''   !< the step it peaked on (its start)
      integer(ik)         :: worst_rh_gap_npatch = 0_ik !< and the patch count then
      real(wp)            :: worst_lignin = 0.0_wp   !< [kgC/m2/day] max |lignin passive-tracer residual|
      real(wp)            :: worst_lambda = 0.0_wp   !< [-] max fraction of any pool withdrawn in one step
      integer(ik)         :: lambda_pool  = 0_ik     !< the pool that attained it
      character(len=19)   :: step_label   = ''       !< the current step's start, for the dates above
   end type slow_ledger_t

contains

   !---------------------------------------------------------------------------------------!
   ! Open a slow step: snapshot the store as the first phase's baseline and clear the      !
   ! declared terms. The step's start labels the soil-carbon checks, active or not.        !
   !---------------------------------------------------------------------------------------!
   subroutine slow_ledger_open(ledger, poly, cfg, step_start)
      type(slow_ledger_t), intent(inout) :: ledger
      type(polygon_t),     intent(in)    :: poly
      type(ed_config_t),   intent(in)    :: cfg
      type(ed_time_t),     intent(in), optional :: step_start
      ledger%step_label = ''
      if (present(step_start)) ledger%step_label = time_to_string(step_start)
      if (.not. ledger%active) return
      ledger%at_mark  = slow_site_store(poly, cfg)
      ledger%decl_in  = slow_store_t()
      ledger%decl_out = slow_store_t()
      ledger%n_step   = ledger%n_step + 1_ik
      if (.not. ledger%has_first) then
         ledger%store_first = ledger%at_mark
         ledger%has_first   = .true.
      end if
   end subroutine slow_ledger_open

   !----- Take the site's store as the next mark's reference without opening a step, for   !
   !      what runs before any step has opened the ledger (a resumed run's owed            !
   !      restructuring). It declares nothing and counts no step. -------------------------!
   subroutine slow_ledger_rebase(ledger, poly, cfg)
      type(slow_ledger_t), intent(inout) :: ledger
      type(polygon_t),     intent(in)    :: poly
      type(ed_config_t),   intent(in)    :: cfg
      if (.not. ledger%active) return
      ledger%at_mark  = slow_site_store(poly, cfg)
      ledger%decl_in  = slow_store_t()
      ledger%decl_out = slow_store_t()
      if (.not. ledger%has_first) then
         ledger%store_first = ledger%at_mark
         ledger%has_first   = .true.
      end if
   end subroutine slow_ledger_rebase

   !---------------------------------------------------------------------------------------!
   ! Declare a boundary term for the phase now open: a quantity that genuinely crosses the !
   ! ledger's boundary and that the model means to (conservation.md §3). It is subtracted  !
   ! from the phase's residual, so declaring a defect would hide it.                       !
   !---------------------------------------------------------------------------------------!
   pure subroutine slow_ledger_declare(ledger, carbon_in, carbon_out, water_in, water_out,         &
                                       energy_in, energy_out)
      type(slow_ledger_t), intent(inout) :: ledger
      real(wp), optional,  intent(in)    :: carbon_in, carbon_out, water_in, water_out
      real(wp), optional,  intent(in)    :: energy_in, energy_out
      if (.not. ledger%active) return
      if (present(carbon_in))  ledger%decl_in%carbon  = ledger%decl_in%carbon  + carbon_in
      if (present(carbon_out)) ledger%decl_out%carbon = ledger%decl_out%carbon + carbon_out
      if (present(water_in))   ledger%decl_in%water   = ledger%decl_in%water   + water_in
      if (present(water_out))  ledger%decl_out%water  = ledger%decl_out%water  + water_out
      if (present(energy_in))  ledger%decl_in%energy  = ledger%decl_in%energy  + energy_in
      if (present(energy_out)) ledger%decl_out%energy = ledger%decl_out%energy + energy_out
   end subroutine slow_ledger_declare

   !----- A thermal-mass bracket: the operators between open and close change the tissues'    !
   !      thermal mass at unchanged temperatures, an exchange the close declares             !
   !      (conservation.md §4). ------------------------------------------------------------!
   subroutine slow_ledger_heat_open(ledger, poly)
      type(slow_ledger_t), intent(inout) :: ledger
      type(polygon_t),     intent(in)    :: poly
      if (.not. ledger%active) return
      ledger%heat_at_open = slow_tissue_energy(poly)
   end subroutine slow_ledger_heat_open

   subroutine slow_ledger_heat_close(ledger, poly)
      type(slow_ledger_t), intent(inout) :: ledger
      type(polygon_t),     intent(in)    :: poly
      if (.not. ledger%active) return
      call slow_ledger_declare(ledger, energy_in = slow_tissue_energy(poly) - ledger%heat_at_open)
   end subroutine slow_ledger_heat_close

   !----- Keep the run's worst of one patch's soil-carbon checks (soil_carbon.md §7, §9),    !
   !      active or not, and the step and patch count where the Rh gap peaked: the gap is    !
   !      zero by construction, so a nonzero peak names the step and patch count to look     !
   !      at. -------------------------------------------------------------------------------!
   subroutine slow_ledger_soil_check(ledger, out, n_patch)
      type(slow_ledger_t),     intent(inout) :: ledger
      type(soil_carbon_out_t), intent(in)    :: out
      integer(ik),             intent(in)    :: n_patch
      ledger%n_seam = ledger%n_seam + 1_ik
      if (abs(out%rh_seam_gap) > ledger%worst_rh_gap) then
         ledger%worst_rh_gap        = abs(out%rh_seam_gap)
         ledger%worst_rh_gap_when   = ledger%step_label
         ledger%worst_rh_gap_npatch = n_patch
      end if
      ledger%worst_lignin = max(ledger%worst_lignin, abs(out%lignin_resid))
      if (out%lambda_max > ledger%worst_lambda) then
         ledger%worst_lambda = out%lambda_max
         ledger%lambda_pool  = out%lambda_pool
      end if
   end subroutine slow_ledger_soil_check

   !---------------------------------------------------------------------------------------!
   ! Close the phase now open: add its residual (conservation.md eq. (1)) to the phase's   !
   ! sums and reopen at the new store. A mark with no operator since the last adds zero,   !
   ! so a caller may mark around a conditional block.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine slow_ledger_mark(ledger, poly, cfg, phase)
      type(slow_ledger_t), intent(inout) :: ledger
      type(polygon_t),     intent(in)    :: poly
      type(ed_config_t),   intent(in)    :: cfg
      integer(ik),         intent(in)    :: phase
      type(slow_store_t) :: now, r
      if (.not. ledger%active) return
      if (phase < 1_ik .or. phase > N_SLOW_PHASE) return
      now = slow_site_store(poly, cfg)
      r%carbon = (now%carbon - ledger%at_mark%carbon) - (ledger%decl_in%carbon - ledger%decl_out%carbon)
      r%water  = (now%water  - ledger%at_mark%water ) - (ledger%decl_in%water  - ledger%decl_out%water )
      r%energy = (now%energy - ledger%at_mark%energy) - (ledger%decl_in%energy - ledger%decl_out%energy)

      associate (s => ledger%resid_sum(phase), a => ledger%resid_abs(phase),                       &
                 w => ledger%resid_worst(phase), g => ledger%decl_gross(phase))
         s%carbon = s%carbon + r%carbon
         s%water = s%water + r%water
         s%energy = s%energy + r%energy
         a%carbon = a%carbon + abs(r%carbon)
         a%water  = a%water  + abs(r%water)
         a%energy = a%energy + abs(r%energy)
         w%carbon = max(w%carbon, abs(r%carbon))
         w%water  = max(w%water,  abs(r%water))
         w%energy = max(w%energy, abs(r%energy))
         g%carbon = g%carbon + abs(ledger%decl_in%carbon) + abs(ledger%decl_out%carbon)
         g%water  = g%water  + abs(ledger%decl_in%water)  + abs(ledger%decl_out%water)
         g%energy = g%energy + abs(ledger%decl_in%energy) + abs(ledger%decl_out%energy)
      end associate
      ledger%n_mark(phase) = ledger%n_mark(phase) + 1_ik

      ledger%at_mark    = now
      ledger%store_last = now
      ledger%decl_in    = slow_store_t()
      ledger%decl_out   = slow_store_t()

      !----- ...then the plausibility test (§6): can what this phase left exist? ----------!
      call check_plausible(ledger, poly, phase)
   end subroutine slow_ledger_mark

   !---------------------------------------------------------------------------------------!
   ! check_plausible -- does any patch now hold a CENTURY pool that cannot exist           !
   ! (soil_carbon_bad_pool, soil_carbon.md §9)? Counts the marks that see one and records  !
   ! the first, whose phase is the one at fault (conservation.md §6).                      !
   !---------------------------------------------------------------------------------------!
   subroutine check_plausible(ledger, poly, phase)
      type(slow_ledger_t), intent(inout) :: ledger
      type(polygon_t),     intent(in)    :: poly
      integer(ik),         intent(in)    :: phase
      real(wp)    :: x(N_SOIL_POOL)
      integer(ik) :: ip, k
      do ip = 1_ik, poly%patch%n
         k = soil_carbon_bad_pool(poly%patch%soil_carbon(ip))
         if (k == 0_ik) cycle
         ledger%n_bad = ledger%n_bad + 1_ik
         if (ledger%bad_phase == 0_ik) then
            ledger%bad_phase = phase
            ledger%bad_patch = ip
            ledger%bad_pool  = k
            x = poly%patch%soil_carbon(ip)%carbon
            ledger%bad_value = x(k)
         end if
         return
      end do
   end subroutine check_plausible

   !=======================================================================================!
   !  The store: every conserved quantity the site holds, per m2 of site (conservation.md  !
   !  §2). Water and energy use the fast budgets' own store functions, so the two tiers    !
   !  agree on what exists; carbon is added here. The recruit pool is valued at c_min per  !
   !  plant; the shed water is a handoff, declared, not stored (§3).                       !
   !=======================================================================================!
   pure function slow_site_store(poly, cfg) result(store)
      type(polygon_t),     intent(in) :: poly
      type(ed_config_t),   intent(in) :: cfg
      type(slow_store_t)              :: store
      real(wp)    :: a, dmol, cap_mass, cap_mol, e_soil, w_soil
      integer(ik) :: ip, i, i0, i1, pf

      do ip = 1_ik, poly%patch%n
         a = patch_polygon_area(poly, ip)
         if (a <= 0.0_wp) cycle

         !----- Soil column: the fast budgets' totals over the patch's site's layers. -----!
         associate (soil => poly%site%surface(poly%patch%owner_site(ip))%soil_profile)
            w_soil = soil_water_store(poly%patch%soil_water(ip)%theta, soil%dz, soil%n_active)
            e_soil = soil_energy_store(poly%patch%soil_energy(ip)%soil_energy, soil%dz, soil%n_active)
         end associate
         store%water  = store%water  + a * w_soil
         store%energy = store%energy + a * e_soil

         !----- Pond and snow: both already extensive per m2 of patch ground. -------------!
         store%water  = store%water  + a * poly%patch%soil_water(ip)%w_surface
         store%energy = store%energy + a * poly%patch%soil_water(ip)%w_surface_energy
         store%water  = store%water  + a * sum(poly%patch%snow(ip)%swe)
         store%energy = store%energy + a * sum(poly%patch%snow(ip)%snow_energy)

         !----- Canopy air at the fast step's capacities; a depth change is the            !
         !      open-volume exchange refresh_canopy_depth declares. -----------------------!
         associate (cas => poly%patch%cas(ip))
            call cas_capacities(rho_air_ref, cas%can_shv, cas%can_depth, dmol, cap_mass, cap_mol)
            store%water  = store%water  + a * cap_mass * cas%can_shv
            store%energy = store%energy + a * cap_mass * cas%can_enthalpy
            store%carbon = store%carbon + a * cap_mol  * cas%can_co2 * umol_2_kgc
         end associate

         !----- The seven CENTURY pools; the lignin is part of the structural pools. ------!
         associate (sc => poly%patch%soil_carbon(ip))
            store%carbon = store%carbon + a * sum(sc%carbon)
         end associate

         !----- Carbon in transit: the recruit carry-forward pool. ------------------------!
         do pf = 1_ik, poly%cohort%pft_table%n
            store%carbon = store%carbon + a * poly%patch%recruit_pool(pf, ip)                   &
                                        * poly%cohort%pft_table%recruit_carbon(pf)
         end do

         !----- Cohorts: nplant is per m2 of patch; the area weight converts to site. -----!
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         associate (c => poly%cohort)
            do i = i0, i1
               store%carbon = store%carbon + a * c%nplant(i) * (c%leaf_carbon(i)                   &
                            + c%fineroot_carbon(i) + c%wood_carbon(i) + c%nonstructural_carbon(i))
               store%water  = store%water  + a * c%nplant(i) * (c%leaf_water_mass(i) + c%wood_water_mass(i))
               store%water  = store%water  + a * cohort_film_water(c, i)
               !----- Tissue energy [J/plant], the cohort's state. ------------------------!
               store%energy = store%energy + a * c%nplant(i) * (c%leaf_energy(i) + c%wood_energy(i))
            end do
         end associate
      end do

   end function slow_site_store

   !---------------------------------------------------------------------------------------!
   ! The site's tissue energy [J/m2 site], the slice of slow_site_store's energy a driver  !
   ! needs to declare a thermal-mass change (conservation.md §4).                          !
   !---------------------------------------------------------------------------------------!
   pure function slow_tissue_energy(poly) result(e)
      type(polygon_t), intent(in) :: poly
      real(wp)    :: e
      integer(ik) :: ip, i, i0, i1
      e = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            e = e + patch_polygon_area(poly, ip) * poly%cohort%nplant(i)                         &
                  * (poly%cohort%leaf_energy(i) + poly%cohort%wood_energy(i))
         end do
      end do
   end function slow_tissue_energy

   !---------------------------------------------------------------------------------------!
   ! The day's net carbon from the fast loop [kgC/m2 site]: GPP less the maintenance       !
   ! respiration already paid back to the canopy air. It sits in the cohorts'              !
   ! accumulators, not a store, so it is declared (conservation.md §3). It takes           !
   ! cohort_carbon_demand's branch, the fast-off stub included, to declare what the        !
   ! allocator consumes.                                                                   !
   !---------------------------------------------------------------------------------------!
   pure function slow_fast_carbon_handover(poly, cfg) result(net_c)
      type(polygon_t),     intent(in) :: poly
      type(ed_config_t),   intent(in) :: cfg
      real(wp)    :: net_c, gross_gpp, resp_maint, a
      integer(ik) :: ip, i, i0, i1
      net_c = 0.0_wp
      do ip = 1_ik, poly%patch%n
         a  = patch_polygon_area(poly, ip)
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         associate (ch => poly%cohort)
            do i = i0, i1
               call cohort_carbon_supply(ch, i, cfg%fast_biophysics_on, cfg%gpp_ref, cfg%dt_years,     &
                                         gross_gpp, resp_maint)
               net_c = net_c + a * ch%nplant(i) * (gross_gpp - resp_maint)
            end do
         end associate
      end do
   end function slow_fast_carbon_handover

   !----- The CO2 this step owes the canopy air [kgC/m2 site]: growth respiration less the          !
   !       starvation deficit. A declared handoff, not a store: it reaches the canopy air          !
   !       in the next fast window (docs/science/conservation.md §3). -------------------------------!
   pure function slow_co2_handoff(poly, cfg) result(c)
      type(polygon_t),     intent(in) :: poly
      type(ed_config_t),   intent(in) :: cfg
      real(wp)    :: c
      integer(ik) :: ip
      c = 0.0_wp
      do ip = 1_ik, poly%patch%n
         c = c + patch_polygon_area(poly, ip) * poly%patch%slow_co2_rate(ip) * umol_2_kgc * cfg%dt_slow
      end do
   end function slow_co2_handoff

   !----- Site-total carbon in this step's litter accumulator [kgC/m2 site], declared leaving at    !
   !       SLOW_PHASE_GROW and arriving at SLOW_PHASE_SOILC; with soil carbon off nothing arrives   !
   !       (docs/science/conservation.md §3). ------------------------------------------------------!
   pure function litter_carbon_total(poly, lit) result(c)
      type(polygon_t),      intent(in) :: poly
      type(litter_input_t), intent(in) :: lit(:)
      real(wp)    :: c
      integer(ik) :: ip
      c = 0.0_wp
      do ip = 1_ik, poly%patch%n
         c = c + patch_polygon_area(poly, ip) * (lit(ip)%labile_grnd + lit(ip)%labile_soil          &
                                      + lit(ip)%struct_grnd + lit(ip)%struct_soil)
      end do
   end function litter_carbon_total

   !----- Phase labels for the report. ----------------------------------------------------!
   pure function slow_phase_name(phase) result(nm)
      integer(ik), intent(in) :: phase
      character(len=16)       :: nm
      select case (phase)
      case (SLOW_PHASE_ALLOCATE)
         nm = 'allocate'
      case (SLOW_PHASE_GROW)
         nm = 'grow+mortality'
      case (SLOW_PHASE_RECRUIT)
         nm = 'recruit'
      case (SLOW_PHASE_COHORT)
         nm = 'cohort fuse/fiss'
      case (SLOW_PHASE_DISTURB)
         nm = 'disturbance'
      case (SLOW_PHASE_PATCH)
         nm = 'patch fuse/term'
      case (SLOW_PHASE_CANOPY)
         nm = 'canopy depth'
      case (SLOW_PHASE_SOILC)
         nm = 'soil carbon'
      case (SLOW_PHASE_TISSUE_WATER)
         nm = 'tissue water'
      case default
         nm = '?'
      end select
   end function slow_phase_name

   !---------------------------------------------------------------------------------------!
   ! End-of-run report (conservation.md §5): a row per phase and currency that moved, the  !
   ! signed sum first, since a leak shows in it and not in an unsigned maximum.            !
   !---------------------------------------------------------------------------------------!
   subroutine slow_ledger_report(ledger)
      type(slow_ledger_t), intent(in) :: ledger
      integer(ik) :: p, nbad
      !----- The soil-carbon seam: a number that should be machine-zero, reported whether or not  !
      !      it ever breached. Lambda, the fraction of a pool one slow step withdraws, says whether !
      !      freezing the pool across the step is sound; it is expected small but nonzero, and a    !
      !      WARNING as it approaches 1, long before a pool could go negative. ---------------------!
      if (ledger%n_seam > 0_ik) then
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_rh]  worst |pool debit - fast Rh| = ',        &
               ledger%worst_rh_gap, ' kgC/m2'
         if (ledger%worst_rh_gap > 0.0_wp)                                                        &
            write(*,'(3a,i0,a)') '                       worst at ', ledger%worst_rh_gap_when,    &
                  ' with ', ledger%worst_rh_gap_npatch, ' patches'
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_lignin]  worst passive-tracer resid = ',      &
               ledger%worst_lignin, ' kgC/m2/day'
         write(*,'(a,es12.3,2a)') ' freeze[soil_carbon]   worst lambda = dt.K.xi / pool  = ',     &
               ledger%worst_lambda, ' [-]  pool ', trim(SOIL_POOL_NAME(ledger%lambda_pool))
         if (ledger%worst_lambda > 0.1_wp)                                                        &
            write(*,'(a)') '   WARNING: a slow step withdraws >10% of a soil-carbon pool --'      &
                        // ' the frozen-pool approximation is degrading; shorten dt_slow.'
      end if
      if (.not. ledger%active .or. ledger%n_step == 0_ik) return
      print '(a)', ''
      print '(a)', '=== slow-loop conservation ledger (site totals, per m2) ==================='
      print '(a,i0,a)', 'slow steps: ', ledger%n_step,                                             &
                        '   [residual = store change - declared boundary terms]'
      !----- `declared`, the phase's gross declared flux, is the scale its residual is     !
      !      read against (§5). -----------------------------------------------------------!
      print '(a)', 'phase              currency    signed sum      sum |r|       worst     declared   marks'
      do p = 1_ik, N_SLOW_PHASE
         if (ledger%n_mark(p) == 0_ik) cycle
         call row(p, 'carbon [kgC]', ledger%resid_sum(p)%carbon, ledger%resid_abs(p)%carbon,       &
                  ledger%resid_worst(p)%carbon, ledger%decl_gross(p)%carbon)
         call row(p, 'water   [kg]', ledger%resid_sum(p)%water,  ledger%resid_abs(p)%water,        &
                  ledger%resid_worst(p)%water,  ledger%decl_gross(p)%water)
         call row(p, 'energy   [J]', ledger%resid_sum(p)%energy, ledger%resid_abs(p)%energy,       &
                  ledger%resid_worst(p)%energy, ledger%decl_gross(p)%energy)
      end do
      print '(a)', '--------------------------------------------------------------------------'
      print '(a,es12.4,a,es12.4,a,es12.4)', 'store drift   carbon ',                               &
            ledger%store_last%carbon - ledger%store_first%carbon, '   water ',                     &
            ledger%store_last%water  - ledger%store_first%water,  '   energy ',                    &
            ledger%store_last%energy - ledger%store_first%energy
      !----- The verdict, eq. (2): each phase and currency against its declared flux and   !
      !      the per-mark floors; a breach is reported, not fatal (§5). -------------------!
      nbad = 0_ik
      do p = 1_ik, N_SLOW_PHASE
         nbad = nbad + slow_ledger_failures(ledger, p)
      end do
      if (nbad == 0_ik) then
         print '(a)', 'VERDICT: every phase closes within tolerance on all three currencies.'
      else
         print '(a,i0,a)', 'WARNING: ', nbad, ' phase/currency residual(s) exceed tolerance --'
         print '(a)', '         a slow-step operator is moving something it does not declare.'
      end if
      !----- ...and the plausibility verdict (§6), reported apart: a run can close on all  !
      !      three currencies while holding a state that cannot exist. --------------------!
      if (ledger%bad_phase /= 0_ik) then
         print '(a)', '--------------------------------------------------------------------------'
         print '(a)', 'IMPOSSIBLE STORE: a phase left a CENTURY pool that cannot exist.'
         print '(5a,i0,a,es12.5)', '   first at phase "', trim(slow_phase_name(ledger%bad_phase)),      &
               '", pool ', trim(SOIL_POOL_NAME(ledger%bad_pool)),                                  &
               ', patch ', ledger%bad_patch, ', value ', ledger%bad_value
         print '(a,i0,a)', '   seen at ', ledger%n_bad, ' operator boundaries in all.'
         print '(a)', '   Conservation says the carbon was accounted for; this says the state is'
         print '(a)', '   not physical. The phase named above is where to look.'
      end if
      print '(a)', '=========================================================================='
   contains
      subroutine row(pp, label, s, aa, w, g)
         integer(ik),      intent(in) :: pp
         character(len=*), intent(in) :: label
         real(wp),         intent(in) :: s, aa, w, g
         if (aa <= 0.0_wp .and. g <= 0.0_wp) return
         print '(a18,1x,a12,4(1x,es12.4),1x,i7)', slow_phase_name(pp), label, s, aa, w, g,         &
               ledger%n_mark(pp)
      end subroutine row
   end subroutine slow_ledger_report

   !----- How many of a phase's three currencies breach eq. (2) of conservation.md §5: the     !
   !      declared flux's share plus the per-mark floors. 0 for a phase that never ran. ---------!
   pure integer(ik) function slow_ledger_failures(ledger, p) result(nbad)
      type(slow_ledger_t), intent(in) :: ledger
      integer(ik),         intent(in) :: p   !< the phase (SLOW_PHASE_*)
      nbad = 0_ik
      if (ledger%n_mark(p) == 0_ik) return
      nbad = breaches(ledger%resid_sum(p)%carbon, ledger%decl_gross(p)%carbon, SLOW_ATOL_CARBON,      &
                      ledger%store_last%carbon, ledger%n_mark(p))                                      &
           + breaches(ledger%resid_sum(p)%water,  ledger%decl_gross(p)%water,  SLOW_ATOL_WATER,       &
                      ledger%store_last%water,  ledger%n_mark(p))                                      &
           + breaches(ledger%resid_sum(p)%energy, ledger%decl_gross(p)%energy, SLOW_ATOL_ENERGY,      &
                      ledger%store_last%energy, ledger%n_mark(p))
   end function slow_ledger_failures

   !----- 1 if this residual breaches the mixed relative/absolute tolerance, else 0. ---------!
   pure integer(ik) function breaches(resid, gross, atol_per_mark, store, nmark)
      real(wp),    intent(in) :: resid, gross, atol_per_mark, store
      integer(ik), intent(in) :: nmark
      real(wp) :: tol
      tol = BUDGET_RTOL_FLUX * abs(gross)                                                          &
          + real(nmark, wp) * (atol_per_mark + SLOW_ROUNDOFF_FRAC * abs(store))
      breaches = 0_ik
      if (abs(resid) > tol) breaches = 1_ik
   end function breaches

end module ed_slow_ledger
