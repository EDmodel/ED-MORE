! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_slow_diag -- the slow step's own diagnostic rows (docs/science/diagnostics.md §4,     !
! "The slow rows"). open_slow_cohort_rows sizes and clears the slow cohort block for the   !
! step; write_slow_patch_rows adds the step's litterfall and recruitment to the patch      !
! rows. A routine that computes an amount the output reports writes its row where it       !
! computes it: the allocation split, the tendencies, the mortality carbon by pathway.      !
!==========================================================================================!
module ed_slow_diag
   use ed_kinds,                only : wp, ik
   use ed_config,               only : ed_config_t
   use ed_polygon_state_vars,   only : polygon_t
   use ed_diag_vars,            only : expand_cohort_diag, reset_cohort_diag, PA_LITTER_LEAF,            &
                                       PA_LITTER_FINEROOT, PA_LITTER_WOOD, PA_LITTER_REPRO,              &
                                       PA_LITTER_STRUCT, PA_RECRUIT_NPLANT
   implicit none
   private

   public :: open_slow_cohort_rows, write_slow_patch_rows

contains

   !----- Size and zero the slow cohort block (one sample per step) before anything writes to it. !
   subroutine open_slow_cohort_rows(poly)
      type(polygon_t), intent(inout) :: poly
      if (.not. poly%cohort%sdiag%active) return
      call expand_cohort_diag(poly%cohort%sdiag, max(poly%cohort%n, 1_ik))
      call reset_cohort_diag(poly%cohort%sdiag)
      poly%cohort%sdiag%n = poly%cohort%n
   end subroutine open_slow_cohort_rows

   !---------------------------------------------------------------------------------------!
   ! The slow patch rows: litterfall and recruitment, written before the restructuring     !
   ! permutes the patch axis. Each adds its per-year rate times dt_slow in seconds, a      !
   ! litter amount made a rate first; the litter rows and the soil's litter input are one  !
   ! number (docs/science/diagnostics.md §4, "The slow rows").                             !
   !---------------------------------------------------------------------------------------!
   subroutine write_slow_patch_rows(poly, cfg, recruitment)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t),   intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)   !< [plant/m2/yr] (pft, patch)
      integer(ik) :: ip

      if (.not. poly%patch%diag%active) return
      !----- Without the fast loop, the step sets the block's patch count itself. ------------!
      if (.not. cfg%fast_biophysics_on) poly%patch%diag%n = poly%patch%n
      do ip = 1_ik, poly%patch%n
         !----- Litter by the tissue it came from, not by the pool it reaches. ---------------!
         poly%patch%diag%v(PA_LITTER_LEAF,     ip) = poly%patch%diag%v(PA_LITTER_LEAF,     ip)    &
                                                   + (poly%patch%litter_in(ip)%leaf_c                &
                                                      / cfg%dt_years) * cfg%dt_slow
         poly%patch%diag%v(PA_LITTER_FINEROOT, ip) = poly%patch%diag%v(PA_LITTER_FINEROOT, ip)    &
                                                   + (poly%patch%litter_in(ip)%fineroot_c            &
                                                      / cfg%dt_years) * cfg%dt_slow
         poly%patch%diag%v(PA_LITTER_WOOD,     ip) = poly%patch%diag%v(PA_LITTER_WOOD,     ip)    &
                                                   + ((poly%patch%litter_in(ip)%wood_c               &
                                                     + poly%patch%litter_in(ip)%storage_c)           &
                                                      / cfg%dt_years) * cfg%dt_slow
         poly%patch%diag%v(PA_LITTER_REPRO,    ip) = poly%patch%diag%v(PA_LITTER_REPRO,    ip)    &
                                                   + (poly%patch%litter_in(ip)%repro_c               &
                                                      / cfg%dt_years) * cfg%dt_slow
         poly%patch%diag%v(PA_LITTER_STRUCT,   ip) = poly%patch%diag%v(PA_LITTER_STRUCT,   ip)    &
                                                   + ((poly%patch%litter_in(ip)%struct_grnd          &
                                                     + poly%patch%litter_in(ip)%struct_soil)         &
                                                      / cfg%dt_years) * cfg%dt_slow
         poly%patch%diag%v(PA_RECRUIT_NPLANT,  ip) = poly%patch%diag%v(PA_RECRUIT_NPLANT,  ip)    &
                                                   + sum(recruitment(:, ip)) * cfg%dt_slow
         !----- Without the fast loop, the step adds the block's weight itself, or every row reads 0. !
         if (.not. cfg%fast_biophysics_on) poly%patch%diag%w(ip) = poly%patch%diag%w(ip) + cfg%dt_slow
      end do
   end subroutine write_slow_patch_rows

end module ed_slow_diag
