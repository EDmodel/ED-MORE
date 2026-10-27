! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_plant_traits -- unit tests for ed_plant_traits: the light-plasticity kernel
! (docs/science/plant_traits.md eqs. (1)-(2)) and the thermal acclimation of a leaf table
! (docs/science/leaf_gas_exchange.md §2, eq. (4)).
!
!   1. NO SHADE      : cum_lai_above = 0 => every target equals its top-of-canopy value.
!   2. SHADE SIGNS   : shaded => SLA up, Vcmax/Rd and leaf turnover down (ED2 directions).
!   3. GRADIENT      : target = trait_toc * exp(kplastic * cum_lai_above); the turnover
!                      rate's target is the reciprocal of a lifespan that lengthens as
!                      exp(-k L).
!   4. RELAX RATE    : update_plastic_trait moves by f = 1 - exp(-rate dt); dt->0 barely
!                      moves, a slow-turnover (shaded) canopy acclimates slowly,
!                      dt >> 1/rate reaches target, and the instantaneous (census) mode
!                      reaches it at once.
!   5. THERMAL       : acclimate_leaf_photo_table (Kattge & Knorr 2007): a warmer growth
!      ACCLIMATION     temperature lowers both entropy terms and the Jmax:Vcmax ratio and
!                      moves the Vcmax and Jmax optima up, the Vcmax optimum of stands
!                      grown at 10 C and 30 C being the page's 28.2 C and 38.2 C; at the
!                      growth temperature where the fit gives the table's own values the
!                      table is unchanged; no other trait moves; the table depends on the
!                      growth temperature alone.
!---~---
program test_plant_traits
   use ed_test_assert,    only : check_close, check_true, test_report
   use ed_kinds,          only : wp, ik
   use ed_constants,      only : t_kelvin
   use ed_therm_lib,      only : peaked_arrhenius_scale
   use ed_pft_params,     only : leaf_photo_table_t, PATH_C3
   use ed_plant_traits,   only : light_plastic_traits, update_plastic_trait, acclimate_leaf_photo_table
   implicit none

   real(wp)    :: sla_t, vc_t, rd_t, turnover_t

   !--- 1. No shade: targets == top-of-canopy.
   call light_plastic_traits(0.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 0.5_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, -0.05_wp, sla_t, vc_t, rd_t, turnover_t)
   call check_close('no shade: sla = toc', sla_t, 15.0_wp)
   call check_close('no shade: vcmax = toc', vc_t, 50.0_wp)
   call check_close('no shade: rd = toc', rd_t, 1.0_wp)
   call check_close('no shade: leaf turnover rate = toc', turnover_t, 0.5_wp)

   !--- 2 + 3. Shaded (cum_lai_above = 3): signs + exact exponential gradient.
   call light_plastic_traits(3.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 0.5_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, -0.05_wp, sla_t, vc_t, rd_t, turnover_t)
   call check_true('shade: SLA increases', sla_t > 15.0_wp)
   call check_true('shade: Vcmax decreases', vc_t  < 50.0_wp)
   call check_true('shade: Rd decreases', rd_t  < 1.0_wp)
   call check_true('shade: leaf turnover slows', turnover_t < 0.5_wp)
   call check_close('shade: SLA gradient', sla_t, 15.0_wp * exp( 0.1_wp * 3.0_wp))
   call check_close('shade: Vcmax gradient', vc_t, 50.0_wp * exp(-0.1_wp * 3.0_wp))
   call check_close('shade: turnover gradient', turnover_t, 0.5_wp * exp(-0.05_wp * 3.0_wp))
   !--- A lifespan of 2 yr lengthening as exp(0.05 L) in shade is this rate, exactly.
   call check_close('shade: turnover = 1 / the lifespan form', turnover_t, 1.0_wp / (2.0_wp * exp(0.05_wp * 3.0_wp)))

   !--- 4. Replacement-weighted update (gradual) + the instantaneous (census) mode.
   call check_close('update: rate dt = 1 => f = 1-exp(-1)',                                                           &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.5_wp, 2.0_wp, .false.),                                           &
           15.0_wp + (1.0_wp - exp(-1.0_wp)) * 5.0_wp)
   call check_true('update: tiny dt barely moves',                                                                    &
           abs(update_plastic_trait(15.0_wp, 20.0_wp, 0.5_wp, 1.0e-3_wp, .false.) - 15.0_wp) < 0.01_wp)
   call check_true('update: slow turnover => slow',                                                                   &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.01_wp, 1.0_wp, .false.) < 15.1_wp)
   call check_true('update: dt >> 1/rate => target',                                                                  &
           abs(update_plastic_trait(15.0_wp, 20.0_wp, 10.0_wp, 10.0_wp, .false.) - 20.0_wp) < 1.0e-6_wp)
   call check_close('update: instant => target regardless of rate dt',                                                &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.01_wp, 1.0e-6_wp, .true.), 20.0_wp)
   call check_close('update: a rate of 0 replaces no leaves',                                                         &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.0_wp, 1.0_wp, .false.), 15.0_wp)

   call test_thermal_acclimation()

   call test_report('test_plant_traits')

contains

   !---~---
   !   5. Thermal acclimation (acclimate_leaf_photo_table; docs/science/leaf_gas_exchange.md
   ! §2, eq. (4)), at Kattge & Knorr (2007)'s coefficients, the [leaf_physiology] defaults,
   ! on a two-PFT table with the reference configuration's Vcmax and Jmax responses.
   !---~---
   subroutine test_thermal_acclimation()
      real(wp), parameter :: DSV_A = 668.39_wp, DSV_B = 1.07_wp   ! [J/mol/K], [J/mol/K2] dS of Vcmax
      real(wp), parameter :: DSJ_A = 659.70_wp, DSJ_B = 0.75_wp   ! [J/mol/K], [J/mol/K2] dS of Jmax
      real(wp), parameter :: JV_A  = 2.59_wp,   JV_B  = 0.035_wp  ! [-], [1/K] Jmax25:Vcmax25
      real(wp), parameter :: T_COLD = t_kelvin + 10.0_wp          ! [K] growth temperatures
      real(wp), parameter :: T_WARM = t_kelvin + 30.0_wp
      real(wp), parameter :: T_REF  = t_kelvin + 25.0_wp
      type(leaf_photo_table_t) :: own, cold, warm, table
      real(wp)    :: topt_cold, topt_warm, dsv_a_ref, dsj_a_ref, jv_a_ref
      integer(ik) :: ip

      own  = photo_table()
      cold = own
      warm = own
      call acclimate_leaf_photo_table(T_COLD, DSV_A, DSV_B, DSJ_A, DSJ_B, JV_A, JV_B, cold)
      call acclimate_leaf_photo_table(T_WARM, DSV_A, DSV_B, DSJ_A, DSJ_B, JV_A, JV_B, warm)

      !---~---
      !   A warm-grown stand has lower entropy terms, and invests relatively less in electron
      ! transport, every PFT.
      !---~---
      do ip = 1_ik, own%n_pft
         call check_true('warmer growth: lower dS of Vcmax', warm%pft(ip)%ds_vcmax < cold%pft(ip)%ds_vcmax,        &
                         cold%pft(ip)%ds_vcmax - warm%pft(ip)%ds_vcmax)
         call check_true('warmer growth: lower dS of Jmax', warm%pft(ip)%ds_jmax < cold%pft(ip)%ds_jmax,           &
                         cold%pft(ip)%ds_jmax - warm%pft(ip)%ds_jmax)
         call check_true('warmer growth: lower Jmax:Vcmax', warm%jmax_vcmax_ratio(ip) < cold%jmax_vcmax_ratio(ip),  &
                         cold%jmax_vcmax_ratio(ip) - warm%jmax_vcmax_ratio(ip))
      end do

      !---~---
      !   What the shift is for: the optimum, measured on the table's own response, sits where
      ! the page puts it for stands grown at 10 C and 30 C.
      !---~---
      topt_cold = peak_temperature(cold%pft(1)%ea_vcmax, cold%pft(1)%hd_vcmax, cold%pft(1)%ds_vcmax)
      topt_warm = peak_temperature(warm%pft(1)%ea_vcmax, warm%pft(1)%hd_vcmax, warm%pft(1)%ds_vcmax)
      call check_close('Vcmax optimum of a stand grown at 10 C [C]', topt_cold - t_kelvin, 28.2_wp, 0.05_wp)
      call check_close('Vcmax optimum of a stand grown at 30 C [C]', topt_warm - t_kelvin, 38.2_wp, 0.05_wp)
      call check_true('and a higher Jmax optimum',                                                                &
                      peak_temperature(warm%pft(1)%ea_jmax, warm%pft(1)%hd_jmax, warm%pft(1)%ds_jmax) >           &
                      peak_temperature(cold%pft(1)%ea_jmax, cold%pft(1)%hd_jmax, cold%pft(1)%ds_jmax))

      !--- Only the entropy terms and the ratio move: every other trait stays the table's.
      do ip = 1_ik, own%n_pft
         associate (a => warm%pft(ip), b => own%pft(ip))
            call check_true('acclimation moves no other trait',                                                   &
                            a%vcmax25 == b%vcmax25 .and. a%jmax25 == b%jmax25 .and. a%ds_rd == b%ds_rd .and.      &
                            a%ea_vcmax == b%ea_vcmax .and. a%hd_vcmax == b%hd_vcmax .and.                         &
                            a%ea_jmax == b%ea_jmax .and. a%hd_jmax == b%hd_jmax .and.                             &
                            warm%tpu_vcmax_ratio(ip) == own%tpu_vcmax_ratio(ip))
         end associate
      end do

      !---~---
      !   At the growth temperature where the fit returns PFT 1's own entropy terms and ratio
      ! (coefficients whose a - b T_C passes through them at 25 C), PFT 1's table is its own.
      !---~---
      dsv_a_ref = own%pft(1)%ds_vcmax + DSV_B * (T_REF - t_kelvin)
      dsj_a_ref = own%pft(1)%ds_jmax  + DSJ_B * (T_REF - t_kelvin)
      jv_a_ref  = own%jmax_vcmax_ratio(1) + JV_B * (T_REF - t_kelvin)
      table = own
      call acclimate_leaf_photo_table(T_REF, dsv_a_ref, DSV_B, dsj_a_ref, DSJ_B, jv_a_ref, JV_B, table)
      call check_close('at the fit''s reference: dS of Vcmax is the table''s', table%pft(1)%ds_vcmax,               &
                       own%pft(1)%ds_vcmax, 1.0e-12_wp * own%pft(1)%ds_vcmax)
      call check_close('at the fit''s reference: dS of Jmax is the table''s', table%pft(1)%ds_jmax,                 &
                       own%pft(1)%ds_jmax, 1.0e-12_wp * own%pft(1)%ds_jmax)
      call check_close('at the fit''s reference: Jmax:Vcmax is the table''s', table%jmax_vcmax_ratio(1),            &
                       own%jmax_vcmax_ratio(1), 1.0e-12_wp)

      !---~---
      !   The table depends on the growth temperature alone: acclimating the cold table to 30 C
      ! gives the warm one, and a repeat call changes nothing.
      !---~---
      table = cold
      call acclimate_leaf_photo_table(T_WARM, DSV_A, DSV_B, DSJ_A, DSJ_B, JV_A, JV_B, table)
      call check_true('an earlier acclimation leaves no trace', all(table%pft%ds_vcmax == warm%pft%ds_vcmax) .and.  &
                      all(table%pft%ds_jmax == warm%pft%ds_jmax) .and.                                            &
                      all(table%jmax_vcmax_ratio == warm%jmax_vcmax_ratio))
      call acclimate_leaf_photo_table(T_WARM, DSV_A, DSV_B, DSJ_A, DSJ_B, JV_A, JV_B, table)
      call check_true('a repeat call changes nothing', all(table%pft%ds_vcmax == warm%pft%ds_vcmax) .and.           &
                      all(table%pft%ds_jmax == warm%pft%ds_jmax) .and.                                            &
                      all(table%jmax_vcmax_ratio == warm%jmax_vcmax_ratio))
   end subroutine test_thermal_acclimation

   !---~---
   !   The temperature [K] where a peaked response is highest, by a scan at 0.01 K from 0 to
   ! 50 C, so the optimum is measured from the response rather than derived.
   !---~---
   pure real(wp) function peak_temperature(ea, hd, ds) result(t_opt)
      real(wp), intent(in) :: ea, hd, ds   ! [J/mol], [J/mol], [J/mol/K]
      real(wp)    :: t, k, k_best
      integer(ik) :: j
      k_best = -1.0_wp
      t_opt  = 0.0_wp
      do j = 0_ik, 5000_ik
         t = t_kelvin + 0.01_wp * real(j, wp)
         k = peaked_arrhenius_scale(1.0_wp, ea, hd, ds, t)
         if (k > k_best) then
            k_best = k
            t_opt  = t
         end if
      end do
   end function peak_temperature

   !---~---
   !   Two C3 PFTs' leaf tables, with the reference configuration's temperature responses
   ! ([leaf_physiology]) and each its own capacities and ratios.
   !---~---
   function photo_table() result(table)
      type(leaf_photo_table_t) :: table
      integer(ik) :: ip
      table%n_pft = 2_ik
      allocate(table%pft(2), table%jmax_vcmax_ratio(2), table%tpu_vcmax_ratio(2))
      table%jmax_vcmax_ratio = [1.8_wp, 1.7_wp]
      table%tpu_vcmax_ratio  = [0.167_wp, 0.160_wp]
      table%pft(1)%vcmax25 = 60.0_wp
      table%pft(2)%vcmax25 = 45.0_wp
      do ip = 1_ik, 2_ik
         table%pft(ip)%pathway  = PATH_C3
         table%pft(ip)%jmax25   = table%jmax_vcmax_ratio(ip) * table%pft(ip)%vcmax25
         table%pft(ip)%ea_vcmax = 65330.0_wp
         table%pft(ip)%ea_jmax  = 43540.0_wp
         table%pft(ip)%hd_vcmax = 200000.0_wp
         table%pft(ip)%hd_jmax  = 200000.0_wp
         table%pft(ip)%ds_vcmax = 650.0_wp
         table%pft(ip)%ds_jmax  = 640.0_wp
         table%pft(ip)%ds_rd    = 490.0_wp
      end do
   end function photo_table

end program test_plant_traits
