! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_carbon_allocation -- how a cohort spends, or fails to pay, its daily carbon   !
! (docs/science/plant_carbon_allocation.md).                                               !
!   * cohort_carbon_demand    -- the step's turnover, taken first, and the                  !
!                                flush-capped growth demands (§1; plant_phenology.md §5).  !
!   * growth_window_so_far    -- the growth window the sink limit sees (eq. (3)).          !
!   * growth_sink_limitation  -- the most wood a cohort can build in a step, whatever its  !
!                                carbon (§2, "Sink limitation of wood growth", eq. (3)).   !
!   * plant_carbon_allocation -- the priority ladder: growth by pool, the net storage      !
!                                change, growth respiration, root exudate and the unpaid   !
!                                maintenance (§2-§3).                                      !
!   * growth_respiration      -- g x the carbon built (eq. (4)).                           !
! All amounts are carbon [kgC/plant]. Turnover is taken before the ladder, so the ladder   !
! only adds carbon; it closes eq. (5).                                                     !
!==========================================================================================!
module meds_plant_carbon_allocation
   use meds_kinds,           only : wp, ik
   use meds_allometry_lib,   only : allometry_t, size_to_wood_carbon, dbh_to_height
   use meds_constants,       only : yr_day
   use meds_plant_phenology, only : leaf_turnover_step, turnover_gate
   implicit none
   private

   public :: cohort_carbon_demand, growth_window_so_far, growth_sink_limitation
   public :: plant_carbon_allocation, growth_respiration

   real(wp), parameter :: DBH_FLOOR = 1.0e-3_wp   !< [cm] keeps dbh^exponent defined for any exponent

contains

   !---------------------------------------------------------------------------------------!
   ! One cohort's carbon demand: apply turnover first (leaf_turnover_step for              !
   ! the leaves; the fine roots at their own rate under the leaves' gate, not dropping with      !
   ! them), and form the flush-capped growth demands from the post-shed pools                    !
   ! (docs/science/plant_phenology.md eq. (9) and §5).                                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine cohort_carbon_demand(                                                          &
            leaf_flush_tendency, leaf_shed_tendency, flush_rate_max, shed_rate_max, leaf_turn,     &
            fineroot_turnover_rate, min_leaf_cover, bare_leaf_cover,                               &
            leaf_carbon, fineroot_carbon, leaf_target, dt_day, r2l,                                &
            leaf_shed_c, fineroot_shed_c, leaf_demand, fineroot_demand,                            &
            retained_frac, leaf_retained_c, fineroot_retained_c)
      real(wp), intent(in)  :: leaf_flush_tendency, leaf_shed_tendency  !< [-]
      real(wp), intent(in)  :: flush_rate_max, shed_rate_max            !< [1/day]
      real(wp), intent(in)  :: leaf_turn, fineroot_turnover_rate        !< [1/yr]
      real(wp), intent(in)  :: min_leaf_cover, bare_leaf_cover          !< [-]
      real(wp), intent(in)  :: leaf_carbon, fineroot_carbon, leaf_target, dt_day, r2l
      real(wp), intent(out) :: leaf_shed_c, fineroot_shed_c
      real(wp), intent(out) :: leaf_demand, fineroot_demand
      real(wp), intent(in)  :: retained_frac        !< [-] share of the loss resorbed; the shed stays the full loss
      real(wp), intent(out) :: leaf_retained_c      !< [kgC/plant] leaf carbon resorbed this step
      real(wp), intent(out) :: fineroot_retained_c  !< [kgC/plant] fine-root carbon resorbed this step
      real(wp) :: senescence, background, flush_cap, leaf_post, root_post, root_pool, retained

      !----- Turnover first: the leaf and fine-root losses, then the post-shed pools. ---------!
      call leaf_turnover_step(leaf_carbon, leaf_target, leaf_flush_tendency, leaf_shed_tendency,    &
                              flush_rate_max, shed_rate_max, leaf_turn, min_leaf_cover,              &
                              bare_leaf_cover, dt_day, senescence, background, flush_cap)
      leaf_shed_c     = senescence + background
      root_pool       = max(fineroot_carbon, 0.0_wp)
      fineroot_shed_c = min(max(fineroot_turnover_rate, 0.0_wp) / yr_day                            &
                            * turnover_gate(leaf_flush_tendency, leaf_shed_tendency) * root_pool * dt_day, &
                            root_pool)
      leaf_post = leaf_carbon     - leaf_shed_c
      root_post = fineroot_carbon - fineroot_shed_c
      retained            = min(max(retained_frac, 0.0_wp), 1.0_wp)
      leaf_retained_c     = retained * leaf_shed_c
      fineroot_retained_c = retained * fineroot_shed_c
      !----- Flush-capped growth demands toward the target, from the post-shed pools. --------!
      leaf_demand     = min(max(0.0_wp, leaf_target       - leaf_post), flush_cap)
      fineroot_demand = min(max(0.0_wp, r2l * leaf_target - root_post), flush_cap * r2l)
   end subroutine cohort_carbon_demand

   !---------------------------------------------------------------------------------------!
   ! A cohort's growth window as the sink limit sees it this step: the diameter already grown   !
   ! in it [cm] and its length [yr] with this step in it. Once the history is full, this step's !
   ! sample replaces the oldest; a new cohort's window is as long as its life                   !
   ! (docs/science/plant_carbon_allocation.md eq. (3); demography.md eq. (2)).                  !
   !---------------------------------------------------------------------------------------!
   pure subroutine growth_window_so_far(count, accum, oldest, n_window, dt_yr, window_yr, grown)
      integer(ik), intent(in)  :: count       !< samples in the ring buffer
      real(wp),    intent(in)  :: accum       !< [cm/yr] their sum
      real(wp),    intent(in)  :: oldest      !< [cm/yr] the sample in the slot this step overwrites
      integer(ik), intent(in)  :: n_window    !< the buffer's length [steps]
      real(wp),    intent(in)  :: dt_yr       !< [yr] step length
      real(wp),    intent(out) :: window_yr   !< [yr] the window this step closes
      real(wp),    intent(out) :: grown       !< [cm] diameter grown in it before this step
      if (count < n_window) then
         window_yr = real(count + 1_ik, wp) * dt_yr
         grown     = accum * dt_yr
      else
         window_yr = real(n_window, wp) * dt_yr
         grown     = (accum - oldest) * dt_yr
      end if
   end subroutine growth_window_so_far

   !---------------------------------------------------------------------------------------!
   ! The sink-limited wood growth of one step [kgC/plant]. The absolute limit is the wood  !
   ! carbon, on the model's own allometry, between today's diameter and today's plus the   !
   ! window's budget less grown_in_window, so a step it binds brings the window's growth   !
   ! to exactly the cap; window_yr = dt_yr with grown_in_window = 0 is a per-step cap.     !
   ! With both limits off the result is huge(1.0_wp): no limit.                            !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_sink_limitation(allom, wood_carbon, dbh, wood_density,                &
                                                  max_relative_growth_rate,                          &
                                                  max_absolute_growth_rate, max_absolute_growth_exponent, &
                                                  dt_yr, window_yr, grown_in_window) result(growth_max)
      type(allometry_t), intent(in) :: allom            !< the cohort's PFT's allometry
      real(wp), intent(in) :: wood_carbon               !< [kgC/plant] wood carbon at the step's start
      real(wp), intent(in) :: dbh                       !< [cm] stem diameter at the step's start
      real(wp), intent(in) :: wood_density              !< [g/cm3] the cohort's wood density
      real(wp), intent(in) :: max_relative_growth_rate  !< [1/yr] the PFT's maximum relative growth rate of wood
      real(wp), intent(in) :: max_absolute_growth_rate  !< [cm/yr] the PFT's maximum diameter growth of a 1 cm stem
      real(wp), intent(in) :: max_absolute_growth_exponent !< [--] how that maximum scales with dbh [cm]
      real(wp), intent(in) :: dt_yr                     !< [yr] step length
      real(wp), intent(in) :: window_yr                 !< [yr] the growth window this step closes (>= dt_yr)
      real(wp), intent(in) :: grown_in_window           !< [cm] diameter already grown in it, before this step
      real(wp)             :: growth_max                !< [kgC/plant] the most wood the step can build
      real(wp) :: dt, dbh_room, dbh_max, wood_now, wood_max

      dt = max(dt_yr, 0.0_wp)
      growth_max = huge(1.0_wp)
      if (max_relative_growth_rate > 0.0_wp)                                                        &
         growth_max = min(growth_max, max_relative_growth_rate * max(wood_carbon, 0.0_wp) * dt)
      if (max_absolute_growth_rate > 0.0_wp) then
         dbh_room   = max_absolute_growth_rate * max(dbh, DBH_FLOOR) ** max_absolute_growth_exponent       &
                    * max(window_yr, dt) - max(grown_in_window, 0.0_wp)
         dbh_max    = dbh + max(dbh_room, 0.0_wp)
         wood_now   = size_to_wood_carbon(allom, dbh,     dbh_to_height(allom, dbh),     wood_density)
         wood_max   = size_to_wood_carbon(allom, dbh_max, dbh_to_height(allom, dbh_max), wood_density)
         growth_max = min(growth_max, max(wood_max - wood_now, 0.0_wp))
      end if
   end function growth_sink_limitation

   !---------------------------------------------------------------------------------------!
   ! Daily allocation for one cohort down the ladder of plant_carbon_allocation.md §2:     !
   ! maintenance debt (from storage when net = gpp - resp_maint < 0) > leaf and fine-root  !
   ! growth (NPP, then storage) > storage refill (NPP only) > reproduction > wood (the     !
   ! residual), at (1+g) per unit of tissue and 1:1 for storage. The leaf and fine-root    !
   ! demands arrive flush-capped and net of the step's shed. Maintenance no pool can pay   !
   ! is `deficit`, which the driver nets out of the cohort's CO2 (plant_respiration.md §6). !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine plant_carbon_allocation(gpp, resp_maint, growth_resp_frac, storage, &
                                    leaf_demand, fineroot_demand, storage_demand, repro_frac,     &
                                    growth_leaf, growth_fineroot, growth_wood, npp_store,         &
                                    growth_repro, growth_resp, deficit, wood_growth_max, exudate)
      real(wp), intent(in)  :: gpp              !< [kgC/plant] gross primary production this step (>= 0)
      real(wp), intent(in)  :: resp_maint       !< [kgC/plant] maintenance respiration this step (>= 0)
      real(wp), intent(in)  :: growth_resp_frac !< [--] construction-cost fraction g (PFT trait)
      real(wp), intent(in)  :: storage          !< [kgC/plant] current nonstructural carbon (available to draw)
      real(wp), intent(in)  :: leaf_demand      !< [kgC/plant] flush-capped leaf growth demand (post-shed; >= 0)
      real(wp), intent(in)  :: fineroot_demand  !< [kgC/plant] flush-capped fine-root growth demand (post-shed; >= 0)
      real(wp), intent(in)  :: storage_demand   !< [kgC/plant] deficit toward the storage target (>= 0)
      real(wp), intent(in)  :: repro_frac       !< [--] fraction of the post-storage residual -> reproduction
      real(wp), intent(out) :: growth_leaf      !< [kgC/plant] carbon grown into leaf (>= 0)
      real(wp), intent(out) :: growth_fineroot  !< [kgC/plant] carbon grown into fine root (>= 0)
      real(wp), intent(out) :: growth_wood      !< [kgC/plant] carbon grown into wood (residual; >= 0)
      real(wp), intent(out) :: npp_store        !< [kgC/plant] net storage change (refill - drawdown; signed)
      real(wp), intent(out) :: growth_repro     !< [kgC/plant] carbon -> reproduction (>= 0; -> recruits)
      real(wp), intent(out) :: growth_resp      !< [kgC/plant] growth respiration on realized growth (>= 0)
      real(wp), intent(out) :: deficit          !< [kgC/plant] unpaid maintenance after storage exhausted (>= 0)
      real(wp), intent(in)  :: wood_growth_max  !< [kgC/plant] sink-limited wood growth (growth_sink_limitation)
      real(wp), intent(out) :: exudate          !< [kgC/plant] residual the sink cannot use (>= 0)
      real(wp) :: net, g, cost, avail, store_left, draw, debt, storage_refill, c_repro, c_wood, c_exudate

      growth_leaf = 0.0_wp
      growth_fineroot = 0.0_wp
      growth_wood = 0.0_wp
      growth_repro = 0.0_wp
      storage_refill = 0.0_wp
      draw = 0.0_wp
      growth_resp = 0.0_wp
      deficit = 0.0_wp
      store_left = max(storage, 0.0_wp)
      net = gpp - resp_maint
      g   = max(growth_resp_frac, 0.0_wp)
      cost = 1.0_wp + g

      if (net >= 0.0_wp) then
         avail = net
      else
         !----- Pay the maintenance debt from storage first; leftover reserves may still grow. ---!
         debt = -net
         draw = min(debt, store_left)
         store_left = store_left - draw
         if (draw < debt) deficit = debt - draw
         avail = 0.0_wp
      end if

      !----- P1: leaf and fine-root growth toward the flush-capped demand, from NPP and then !
      !          storage, so a leafless canopy can flush on its reserves. -----------------!
      call fill_carbon_demand(leaf_demand,     g, cost, avail, store_left, draw, growth_leaf,     growth_resp)
      call fill_carbon_demand(fineroot_demand, g, cost, avail, store_left, draw, growth_fineroot, growth_resp)
      !----- P2: refill storage from remaining NPP only (no construction cost). -----------------!
      storage_refill = min(max(storage_demand, 0.0_wp), max(avail, 0.0_wp))
      avail   = avail - storage_refill
      !----- P3: reproduction -- a fraction of the post-storage residual, construction-charged. --!
      c_repro      = max(repro_frac, 0.0_wp) * max(avail, 0.0_wp)
      growth_repro = c_repro / cost
      growth_resp  = growth_resp + growth_respiration(growth_repro, g)
      avail        = avail - c_repro
      !----- P4: wood takes the residual up to the sink limit; the rest is root exudate at the !
      !          same (1+g) cost, so the limit leaves growth respiration unchanged. -------!
      c_wood      = max(avail, 0.0_wp)
      growth_wood = c_wood / cost
      c_exudate   = 0.0_wp
      if (growth_wood > wood_growth_max) then
         c_exudate   = growth_wood - max(wood_growth_max, 0.0_wp)
         growth_wood = max(wood_growth_max, 0.0_wp)
      end if
      growth_resp = growth_resp + growth_respiration(growth_wood + c_exudate, g)
      exudate     = c_exudate

      npp_store = storage_refill - draw
   end subroutine plant_carbon_allocation

   !---------------------------------------------------------------------------------------!
   ! Fund a tissue demand from NPP first, then from storage, at (1+g) carbon per unit built; !
   ! add the tissue to `built`, the storage drawn to `draw` and the respired carbon to     !
   ! `growth_resp` (docs/science/plant_carbon_allocation.md §3).                           !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine fill_carbon_demand(demand, g, cost, avail, store_left, draw, built, growth_resp)
      real(wp), intent(in)    :: demand, g, cost
      real(wp), intent(inout) :: avail, store_left, draw, built, growth_resp
      real(wp) :: want, build
      want = max(demand, 0.0_wp)
      if (want <= 0.0_wp) return
      !----- from NPP first (carbon-limited: at most avail/cost units of tissue). --------------!
      build       = min(want, max(avail, 0.0_wp) / cost)
      avail       = avail - build * cost
      built       = built + build
      growth_resp = growth_resp + growth_respiration(build, g)
      want        = want - build
      !----- then from storage, if still needed and available. --------------------------------!
      if (want > 0.0_wp .and. store_left > 0.0_wp) then
         build       = min(want, store_left / cost)
         store_left  = store_left - build * cost
         draw        = draw + build * cost
         built       = built + build
         growth_resp = growth_resp + growth_respiration(build, g)
      end if
   end subroutine fill_carbon_demand

   !---------------------------------------------------------------------------------------!
   ! Growth (construction) respiration [kgC/plant] = g x the carbon actually built this    !
   ! step (docs/science/plant_carbon_allocation.md eq. (4)). A tissue- or                  !
   ! temperature-dependent construction cost would replace this one form.                  !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_respiration(npp_growth, growth_resp_factor) result(rg)
      real(wp), intent(in) :: npp_growth          !< [kgC/plant] carbon allocated to new growth this step
      real(wp), intent(in) :: growth_resp_factor  !< [--] per-PFT construction-cost fraction
      real(wp)             :: rg                   !< [kgC/plant] growth respiration
      rg = growth_resp_factor * max(0.0_wp, npp_growth)
   end function growth_respiration

end module meds_plant_carbon_allocation
