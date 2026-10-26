! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_carbon_allocation -- unit tests for the elemental, GROWTH-ONLY daily carbon      !
! allocation kernel + growth_respiration. (Leaf loss is applied upstream by the driver from   !
! ed_plant_phenology::leaf_turnover_step -- tested in test_plant_phenology.) The equations are     !
! docs/science/plant_carbon_allocation.md eqs. (4)-(5).                                    !
!                                                                                          !
!   1. GROWTH RESP     : growth_respiration = g * max(0, npp_growth).                            !
!   2. CLOSURE         : (growth pools + npp_store) - deficit = (gpp - resp) - growth_resp.       !
!   3. GROWTH-RESP CHG : charged (1+g) on REALIZED growth only; storage refill is exempt.         !
!   4. WOOD RESIDUAL   : a budget that only covers the leaf demand => wood 0.                     !
!   5. STORAGE GROWTH  : storage funds leaf growth EVEN when net < 0 (spring leaf-out).            !
!   6. STARVING        : maintenance debt beyond storage => a deficit, no growth.                   !
!   7. SINK LIMIT      : relative = rate x wood x dt; absolute = the wood that grows the diameter    !
!                        by rate x dbh^exponent x dt, under and over the height cap; the tighter     !
!                        applies; 0 = off.                                                          !
!   8. SINK CAP        : wood stops at the limit, the rest is exudate at cost (1+g), carbon closes, !
!                        growth respiration is unchanged; a limit that does not bind changes       !
!                        nothing.                                                                   !
!   9. MATURE COHORT   : with reproduction and growth respiration on, the books still close,  !
!                        and a surplus grows wood.                                         !
!==========================================================================================!
program test_plant_carbon_allocation
   use ed_test_assert, only : check_close, check_true, test_report
   use ed_kinds,           only : wp, ik
   use ed_plant_carbon_allocation,   only : plant_carbon_allocation, growth_respiration,                 &
                                            growth_sink_limitation
   use ed_allometry_lib,             only : allometry_t, size_to_wood_carbon, wood_to_dbh, dbh_to_height
   implicit none

   !----- A stem for the sink-limit tests: wood density, height cap and aboveground share. -----!
   real(wp), parameter :: RHO = 0.6_wp, HMAX = 35.0_wp, AGF = 0.7_wp
   type(allometry_t), parameter :: ALLOM = allometry_t(hgt_max = HMAX, aboveground_frac = AGF)   ! shipped coefficients

   call test_growth_respiration()
   call test_closure()
   call test_growth_resp_charge()
   call test_wood_residual()
   call test_storage_growth()
   call test_starving()
   call test_sink_limitation()
   call test_sink_cap()
   call test_mature_cohort()

   call test_report('test_plant_carbon_allocation')

contains



   !----- The growth-side carbon-closure identity the kernel guarantees on every call. -----!
   subroutine check_closure(name, gpp, resp, gl, gf, gw, gr, gs, gresp, def)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: gpp, resp, gl, gf, gw, gr, gs, gresp, def
      call check_close(name, (gl + gf + gw + gr + gs) - def, (gpp - resp) - gresp)
   end subroutine check_closure

   !----- 1. growth_respiration = g * max(0, npp_growth). -----------------------------------!
   subroutine test_growth_respiration()
      call check_close('growth_resp: 0.3 * 10 = 3', growth_respiration(10.0_wp, 0.3_wp), 3.0_wp)
      call check_close('growth_resp: npp<=0 => 0',  growth_respiration(-5.0_wp, 0.3_wp), 0.0_wp)
   end subroutine test_growth_respiration

   !----- 2. Carbon closes (g = 0 for the simplest arithmetic): surplus -> wood. ------------!
   subroutine test_closure()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.3_wp, 0.2_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
      call check_close('closure: growth_leaf = demand', gl, 0.3_wp)
      call check_close('closure: growth_wood = residual', gw, 0.5_wp)
      call check_true ('closure: no deficit', def == 0.0_wp)
      call check_closure('closure: g=0', 1.0_wp, 0.0_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_closure

   !----- 3. Growth respiration charged (1+g) on realized growth; storage refill exempt. ----!
   subroutine test_growth_resp_charge()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      ! g = 0.5: leaf demand 0.3 built at cost 1.5 (gresp += 0.15); storage 0.2 at 1:1 (exempt);
      ! residual 0.35 -> wood 0.35/1.5.
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.2_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
      call check_close('growth resp: on realized growth', gresp, 0.15_wp + 0.5_wp*(0.35_wp/1.5_wp))
      call check_close('growth resp: storage exempt (npp_store = 0.2)', gs, 0.2_wp)
      call check_close('growth resp: wood tissue = residual/(1+g)', gw, 0.35_wp/1.5_wp)
      call check_closure('growth resp: closes', 1.0_wp, 0.0_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_growth_resp_charge

   !----- 4. Wood is the residual sink: a budget that only covers the leaf demand => wood 0. !
   subroutine test_wood_residual()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      call plant_carbon_allocation(0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
      call check_close('residual: growth_leaf = demand', gl, 0.5_wp)
      call check_close('residual: growth_wood = 0',      gw, 0.0_wp)
   end subroutine test_wood_residual

   !----- 5. Storage funds leaf growth even when net < 0 (spring leaf-out from reserves). ---!
   subroutine test_storage_growth()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      ! net = -0.05 (paid from storage), then 0.1 leaf built from the remaining reserves (g=0).
      call plant_carbon_allocation(0.0_wp, 0.05_wp, 0.0_wp, 1.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
      call check_true ('storage growth: leaf grows despite net<0', gl > 0.0_wp)
      call check_close('storage growth: growth_leaf = demand',     gl, 0.1_wp)
      call check_close('storage growth: storage draw = debt+growth', gs, -0.15_wp)
      call check_true ('storage growth: no deficit',               def == 0.0_wp)
      call check_closure('storage growth: closes', 0.0_wp, 0.05_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_storage_growth

   !----- 6. Maintenance debt beyond storage => a deficit; no growth. ----------------------!
   subroutine test_starving()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      call plant_carbon_allocation(0.0_wp, 1.0_wp, 0.0_wp, 0.3_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
      call check_true ('starving: deficit left',       def > 0.0_wp)
      call check_close('starving: deficit = debt - storage', def, 0.7_wp)
      call check_close('starving: no leaf growth',      gl, 0.0_wp)
      call check_close('starving: storage drained',     gs, -0.3_wp)
      call check_closure('starving: closes', 0.0_wp, 1.0_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_starving

   !----- 7. The sink limits: relative = rate x wood x dt; absolute = the wood that grows the  !
   !          diameter by rate x the window, less what the window has grown, exactly, on either  !
   !          side of the height cap; the tighter applies; 0 = off. --------------------------!
   subroutine test_sink_limitation()
      real(wp) :: d, w, cap
      integer  :: k
      call check_close('sink limit: 0.2/yr x 5 kgC x 0.5 yr', limit(5.0_wp, 10.0_wp, 0.2_wp, 0.0_wp, 0.5_wp), 0.5_wp)
      call check_true ('sink limit: both rates 0 => no limit', limit(5.0_wp, 10.0_wp, 0.0_wp, 0.0_wp, 0.5_wp) > 1.0e30_wp)
      call check_close('sink limit: no wood => no growth',  limit(-1.0_wp, 10.0_wp, 0.2_wp, 0.0_wp, 0.5_wp), 0.0_wp)
      !----- A sapling below the height cap and a big tree above it (35 m is reached near 60 cm). -!
      do k = 1, 2
         d   = merge(3.0_wp, 90.0_wp, k == 1)
         w   = size_to_wood_carbon(ALLOM, d, dbh_to_height(ALLOM, d), RHO)
         cap = limit(w, d, 0.0_wp, 0.4_wp, 0.5_wp)
         call check_close('sink limit: absolute => dbh grows by rate x dt',                          &
                          wood_to_dbh(ALLOM, w + cap, RHO) - d, 0.2_wp, 1.0e-9_wp)
         !----- With exponent 0.5 the cap is rate x sqrt(dbh): 0.4 x sqrt(d) x 0.5 yr. ----------!
         cap = limit(w, d, 0.0_wp, 0.4_wp, 0.5_wp, 0.5_wp)
         call check_close('sink limit: absolute x dbh^0.5 => dbh grows by rate x sqrt(dbh) x dt',     &
                          wood_to_dbh(ALLOM, w + cap, RHO) - d, 0.2_wp * sqrt(d), 1.0e-9_wp)
      end do
      !----- Both set: the tighter one. At 3 cm a 0.4 cm/yr cap is ~30 % of wood a year. -------!
      d = 3.0_wp
      w = size_to_wood_carbon(ALLOM, d, dbh_to_height(ALLOM, d), RHO)
      call check_close('sink limit: both => the tighter (relative)', limit(w, d, 0.05_wp, 0.4_wp, 1.0_wp), 0.05_wp * w)
      call check_close('sink limit: both => the tighter (absolute)', limit(w, d, 5.0_wp, 0.4_wp, 1.0_wp),  &
                       limit(w, d, 0.0_wp, 0.4_wp, 1.0_wp))
      !----- The absolute cap is a budget over the growth window. A daily step of a stem that grew  !
      !      nothing in its one-year window may grow the whole year's 0.4 cm; one that has grown    !
      !      0.1 cm may grow 0.3; one that has grown the 0.4 may grow nothing. ---------------------!
      cap = limit(w, d, 0.0_wp, 0.4_wp, 1.0_wp / 365.0_wp, window=1.0_wp, grown=0.0_wp)
      call check_close('sink limit: window, nothing grown => the window''s growth in one step',       &
                       wood_to_dbh(ALLOM, w + cap, RHO) - d, 0.4_wp, 1.0e-9_wp)
      cap = limit(w, d, 0.0_wp, 0.4_wp, 1.0_wp / 365.0_wp, window=1.0_wp, grown=0.1_wp)
      call check_close('sink limit: window, 0.1 cm grown => 0.3 cm left',                             &
                       wood_to_dbh(ALLOM, w + cap, RHO) - d, 0.3_wp, 1.0e-9_wp)
      call check_close('sink limit: window used up => no wood',                                       &
                       limit(w, d, 0.0_wp, 0.4_wp, 1.0_wp / 365.0_wp, window=1.0_wp, grown=0.5_wp), 0.0_wp)
   end subroutine test_sink_limitation

   !----- growth_sink_limitation on the test stem; the absolute cap's dbh exponent defaults to 0, and   !
   !      its window to this one step with nothing grown in it (a per-step cap). ---------------------!
   real(wp) function limit(wood, dbh, rgr, agr, dt, expo, window, grown)
      real(wp), intent(in) :: wood, dbh, rgr, agr, dt
      real(wp), intent(in), optional :: expo, window, grown
      real(wp) :: e, win, gro
      e = 0.0_wp
      win = dt
      gro = 0.0_wp
      if (present(expo))   e   = expo
      if (present(window)) win = window
      if (present(grown))  gro = grown
      limit = growth_sink_limitation(ALLOM, wood_carbon=wood, dbh=dbh, wood_density=RHO,                       &
                                     max_relative_growth_rate=rgr,                                        &
                                     max_absolute_growth_rate=agr, max_absolute_growth_exponent=e,        &
                                     dt_yr=dt, window_yr=win, grown_in_window=gro)
   end function limit

   !----- 8. Wood stops at the sink limit; the rest is exuded, construction-charged. -------!
   subroutine test_sink_cap()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      real(wp) :: gl0, gf0, gw0, gs0, gr0, gresp0, def0
      ! g = 0.5: leaf 0.3 at cost 1.5 leaves 0.55, which would build 0.55/1.5 of wood. The sink
      ! takes 0.1 (cost 0.15); the other 0.40 makes 0.40/1.5 of exudate at the same cost, so growth
      ! respiration is what it would be without the limit.
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def,                                                &
           wood_growth_max=limit(10.0_wp, 10.0_wp, 0.01_wp, 0.0_wp, 1.0_wp), exudate=ex)
      call check_close('sink cap: wood = the limit',                 gw, 0.1_wp)
      call check_close('sink cap: exudate = (residual - wood cost)/(1+g)', ex, 0.4_wp/1.5_wp)
      call check_close('sink cap: growth resp on wood + exudate', gresp, 0.15_wp + 0.5_wp*(0.55_wp/1.5_wp))
      call check_close('sink cap: closes with the exudate', (gl + gf + gw + gr + gs + ex) - def,  &
                       (1.0_wp - 0.0_wp) - gresp)
      ! A limit above the residual changes nothing, to the bit.
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl0, gf0, gw0, gs0, gr0, gresp0, def0, huge(1.0_wp), ex)
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def,                                                &
           wood_growth_max=limit(10.0_wp, 10.0_wp, 0.0_wp, 0.0_wp, 1.0_wp), exudate=ex)
      call check_true('sink cap: no limit => bit-identical wood', gw == gw0 .and. gresp == gresp0)
      call check_close('sink cap: no limit => no exudate', ex, 0.0_wp)
   end subroutine test_sink_cap

   !----- 9. A mature cohort with a surplus: the books close with reproduction and growth        !
   !      respiration in them, and the surplus grows wood. ---------------------------------!
   subroutine test_mature_cohort()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def, ex
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.05_wp, 0.05_wp, 0.02_wp, 0.3_wp, &
           gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
      call check_close('mature: (growth + store) - deficit = net - growth resp', (gl+gf+gw+gr+gs) - def,     &
                       1.0_wp - gresp, 1.0e-9_wp)
      call check_true ('mature: a mature cohort allocates to reproduction', gr > 0.0_wp)
      call check_true ('mature: the surplus grows wood', gw > 0.0_wp)
      call check_true ('mature: growth respiration is charged on realized growth', gresp > 0.0_wp)
   end subroutine test_mature_cohort

end program test_plant_carbon_allocation
