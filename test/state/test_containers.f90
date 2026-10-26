! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_containers -- the state layer's containers (src/state):                             !
!   1. THE COHORT CSR (ed_polygon_state_vars, ed_init_bare_ground: finalize_init;           !
!      ed_slow_update: terminate_cohorts): counts, offsets and owners, the height sort,     !
!      persistent unique ids, and termination's compaction, the diagnostic blocks' counts   !
!      included.                                                                            !
!   2. THE FAST DIAGNOSTIC BLOCK (ed_diag_vars: reorder_cohort_diag, copy_cohort_diag,      !
!      clear_cohort_diag, fuse_cohort_diag, cohort_diag_value): it moves in lockstep with   !
!      the cohorts, and a fused cohort's diagnostics blend by each field's own kind.        !
!==========================================================================================!
program test_containers
   use ed_kinds,             only : wp, ik
   use ed_config,            only : ed_config_t
   use ed_polygon_state_vars,         only : polygon_t
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_slow_update, only : terminate_cohorts
   use ed_diag_vars,         only : cohort_diag_block, alloc_cohort_diag, reorder_cohort_diag, copy_cohort_diag,   &
                                    fuse_cohort_diag, clear_cohort_diag, cohort_diag_value, CO_LEAF_TEMP,         &
                                    CO_SAPFLOW, CO_ABS_SW, CO_WOOD_TEMP, CO_GSW
   use ed_test_support, only : banner, build_test_config, check, check_close
   use ed_slow_state,  only : kill_report_t
   implicit none
   type(kill_report_t)  :: kills   !< what the cull and the treefall killed
   real(wp)            :: shed_water   !< an operator's shed water, unread here

   type(ed_config_t) :: cfg
   type(polygon_t)  :: poly
   integer(ik)         :: ip, k, i0, i1, nbefore, a, b, gid_dropped
   integer(ik)         :: nread
   real(wp)            :: xbuf(8), xshort(2)

   call banner('containers: CSR + sort + termination')
   cfg = build_test_config()
   call init_bare_ground(poly, cfg, 2_ik)

   !=== 1. THE COHORT CSR. =================================================================!
   !----- Add cohorts (same PFT so height order == DBH order), unsorted. ------------------!
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.2_wp,  5.0_wp)
   call add_cohort(poly, cfg, 2_ik, 1_ik, 0.2_wp, 20.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 15.0_wp)
   call add_cohort(poly, cfg, 2_ik, 1_ik, 0.15_wp, 3.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.1_wp,  9.0_wp)
   call finalize_init(poly)

   !----- CSR: counts sum to n, offsets contiguous, owners match their slice. -------------!
   call check(sum(poly%patch%cohort_count(1:poly%patch%n)) == poly%cohort%n, 'cohort_count sum /= n')
   call check(poly%patch%cohort_offset(1) == 1_ik, 'cohort_offset(1) /= 1')
   call check(poly%patch%cohort_count(1) == 3_ik, 'patch 1 should have 3 cohorts')
   call check(poly%patch%cohort_count(2) == 2_ik, 'patch 2 should have 2 cohorts')
   do ip = 1_ik, poly%patch%n
      i0 = poly%patch%cohort_offset(ip)
      i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
      if (ip < poly%patch%n) call check(poly%patch%cohort_offset(ip+1) == i1 + 1_ik, 'cohort_offset not contiguous')
      do k = i0, i1
         call check(poly%cohort%owner_patch(k) == ip, 'owner_patch mismatch in slice')
      end do
      !----- Sorted height-descending within the slice. ---------------------------------!
      do k = i0, i1 - 1_ik
         call check(poly%cohort%height(k) >= poly%cohort%height(k+1), 'cohorts not height-descending')
      end do
   end do

   !----- Persistent global ids: assigned (>0) and unique across the live set, even after   !
   !       the sort permuted the cohort order. --------------------------------------------!
   call check(all(poly%cohort%global_id(1:poly%cohort%n) > 0_ik), 'global_id not assigned')
   do a = 1_ik, poly%cohort%n
      do b = a + 1_ik, poly%cohort%n
         call check(poly%cohort%global_id(a) /= poly%cohort%global_id(b), 'global_id not unique')
      end do
   end do

   !----- Termination removes a sub-threshold cohort and keeps CSR consistent. ------------!
   !                                                                                          !
   !      THE DIAGNOSTIC BLOCKS RIDE THE SAME COMPACTION, COUNT INCLUDED. They are armed   !
   !      here because a cull is the operator that SHRINKS the cohort array, and a block   !
   !      whose `n` stays at the pre-cull value claims more live slots than the array has. !
   !      extract_variable sizes its scratch buffer from cohort%n and cohort_diag_value    !
   !      fills 1..d%n into it, so a stale count is an out-of-bounds WRITE -- silent in    !
   !      Release, a SIGSEGV some steps later. A stand that never shrinks cannot show it   !
   !      (the Ithaca reference goes 114 -> 125 cohorts over three years). ----------------!
   call alloc_cohort_diag(poly%cohort%diag,  poly%cohort%cap, .true.)
   call alloc_cohort_diag(poly%cohort%sdiag, poly%cohort%cap, .true., slow_rows=.true.)
   poly%cohort%diag%n  = poly%cohort%n
   poly%cohort%sdiag%n = poly%cohort%n

   nbefore = poly%cohort%n
   gid_dropped = poly%cohort%global_id(poly%patch%cohort_offset(1))
   poly%cohort%nplant(poly%patch%cohort_offset(1)) = 1.0e-12_wp     ! below negligible_nplant
   call terminate_cohorts(poly, cfg, shed_water, kills)
   call check(poly%cohort%n == nbefore - 1_ik, 'termination did not drop exactly one cohort')
   call check(sum(poly%patch%cohort_count(1:poly%patch%n)) == poly%cohort%n, 'CSR inconsistent after terminate')
   !----- The culled cohort's id never reappears (ids are not reused). --------------------!
   call check(all(poly%cohort%global_id(1:poly%cohort%n) /= gid_dropped), 'dropped global_id reused')

   call check(poly%cohort%diag%n  == poly%cohort%n, 'fast diag block kept a stale n across the cull')
   call check(poly%cohort%sdiag%n == poly%cohort%n, 'slow diag block kept a stale n across the cull')

   !----- ...and independently, the extractor must never write past the caller's array, however  !
   !      the counts got out of step. `xshort` is deliberately smaller than the live cohort       !
   !      count: the reader has to report what it actually wrote, not what the block holds. ------!
   call cohort_diag_value(poly%cohort%diag, 1_ik, xbuf, nread, 0.0_wp)
   call check(nread == poly%cohort%n, 'extractor did not report the live cohort count')
   call cohort_diag_value(poly%cohort%diag, 1_ik, xshort, nread, 0.0_wp)
   call check(nread <= int(size(xshort), ik), 'extractor reported more than the caller''s array holds')

   call test_diag_lockstep()

   write(*,'(a)') '   PASS'

contains

   !=======================================================================================!
   !  2. THE LOCKSTEP + FUSION CONTRACT of the fast diagnostic block (ed_diag_vars). This is    !
   !  the highest-risk part of the output subsystem: the block is filled by the fast loop and    !
   !  read at the output tick, with restructuring in between, so slot i must keep meaning        !
   !  cohort i across every permutation -- and a fused cohort's diagnostics must blend by each   !
   !  field's own kind.                                                                          !
   !                                                                                          !
   !  The fields exercised below are deliberately one of EACH kind, so a blend that used one      !
   !  rule for all of them cannot pass: leaf_temp and wood_temp fuse on their tissue's heat         !
   !  capacity, gsw is INTENSIVE (leaf-area-weighted), sapflow is EXTENSIVE (nplant-weighted),      !
   !  abs_sw is GROUND-referenced (summed).                                                         !
   !=======================================================================================!
   subroutine test_diag_lockstep()
      type(cohort_diag_block) :: d
      real(wp)    :: x(8)
      integer(ik) :: perm(4), n

      call alloc_cohort_diag(d, 8_ik, .true.)
      d%n = 4_ik
      d%v(CO_LEAF_TEMP, 1:4) = [290.0_wp, 300.0_wp, 310.0_wp, 320.0_wp]
      d%v(CO_WOOD_TEMP, 1:4) = [280.0_wp, 285.0_wp, 295.0_wp, 305.0_wp]
      d%v(CO_GSW,       1:4) = [  0.1_wp,   0.2_wp,   0.3_wp,   0.4_wp]
      d%v(CO_SAPFLOW,   1:4) = [  1.0_wp,   2.0_wp,   3.0_wp,   4.0_wp]
      d%v(CO_ABS_SW,    1:4) = [ 10.0_wp,  20.0_wp,  30.0_wp,  40.0_wp]
      d%w(1:4)               = 1.0_wp

      !----- REORDER: one statement must permute EVERY field, not the ones someone remembered. !
      perm = [4_ik, 3_ik, 2_ik, 1_ik]
      call reorder_cohort_diag(d, perm, 4_ik)
      call cohort_diag_value(d, CO_LEAF_TEMP, x, n, 0.0_wp)
      call check(n == 4_ik, 'reorder keeps the live count')
      call check_close(x(1), 320.0_wp, 1.0e-12_wp, 'reorder permuted leaf_temp')
      call cohort_diag_value(d, CO_SAPFLOW, x, n, 0.0_wp)
      call check_close(x(1), 4.0_wp, 1.0e-12_wp, 'reorder permuted sapflow in lockstep')
      call cohort_diag_value(d, CO_ABS_SW, x, n, 0.0_wp)
      call check_close(x(1), 40.0_wp, 1.0e-12_wp, 'reorder permuted abs_sw in lockstep')

      !----- COPY SLOT (a cohort split): the daughter is a full copy. -----!
      call copy_cohort_diag(d, 5_ik, 1_ik)
      call check_close(d%v(CO_SAPFLOW, 5_ik), d%v(CO_SAPFLOW, 1_ik), 1.0e-12_wp, 'copy_slot copies every field')
      call check_close(d%w(5_ik), d%w(1_ik), 1.0e-12_wp, 'copy_slot carries the dt weight')

      !----- CLEAR SLOT (a recruit into a reused, stale cull slot). -----!
      call clear_cohort_diag(d, 5_ik)
      call check(d%v(CO_SAPFLOW, 5_ik) == 0.0_wp, 'clear_slot zeroes a reused slot')
      call check(d%w(5_ik) == 0.0_wp,             'clear_slot zeroes the weight')
      !----- ...and reads as the caller's "no samples" value, not as 0 K. -----!
      d%n = 5_ik
      call cohort_diag_value(d, CO_LEAF_TEMP, x, n, -1.0_wp)
      call check(x(5) == -1.0_wp .and. x(1) == 320.0_wp, 'a slot with no samples reads the empty value')
      d%n = 4_ik

      !----- FUSE slot 2 into slot 1 with leaf-area weights (3, 1), nplant weights (2, 6), leaf    !
      !      heat capacities (1, 4) and wood heat capacities (5, 1). After the reorder above:     !
      !      slot1 = (leaf 320, wood 305, gsw 0.4, sapflow 4, abs_sw 40), slot2 = (310, 295, 0.3,   !
      !      3, 30).                                                                             !
      !        LEAF HEAT  leaf_temp = (1*320 + 4*310)/5      = 312.0  (317.5 on leaf area)         !
      !        WOOD HEAT  wood_temp = (5*305 + 1*295)/6      = 303.33 (302.5 on leaf area)         !
      !        INTENSIVE  gsw       = (3*0.4 + 1*0.3)/4      = 0.375                               !
      !        EXTENSIVE  sapflow   = (2*4   + 6*3  )/8      = 3.25                                !
      !        GROUND     abs_sw    =  40 + 30                = 70                                  !
      call fuse_cohort_diag(d, 1_ik, 2_ik, 3.0_wp, 1.0_wp, 2.0_wp, 6.0_wp, 1.0_wp, 4.0_wp, 5.0_wp, 1.0_wp)
      call cohort_diag_value(d, CO_LEAF_TEMP, x, n, 0.0_wp)
      call check_close(x(1), 312.0_wp, 1.0e-12_wp, 'leaf temperature fuses on leaf heat capacity')
      call cohort_diag_value(d, CO_WOOD_TEMP, x, n, 0.0_wp)
      call check_close(x(1), 1820.0_wp / 6.0_wp, 1.0e-12_wp, 'wood temperature fuses on wood heat capacity')
      call cohort_diag_value(d, CO_GSW, x, n, 0.0_wp)
      call check_close(x(1), 0.375_wp, 1.0e-12_wp, 'INTENSIVE field fuses leaf-area-weighted')
      call cohort_diag_value(d, CO_SAPFLOW, x, n, 0.0_wp)
      call check_close(x(1), 3.25_wp, 1.0e-12_wp, 'EXTENSIVE field fuses nplant-weighted')
      call cohort_diag_value(d, CO_ABS_SW, x, n, 0.0_wp)
      call check_close(x(1), 70.0_wp, 1.0e-12_wp, 'GROUND-referenced field fuses by SUM')
      !----- The five answers are mutually distinct, so a blend that applied one rule to all     !
      !      of them would fail at least four of the assertions above. -----------------------!

      !----- An INACTIVE block is a total no-op: nothing allocated, every entry point safe. ---!
      block
         type(cohort_diag_block) :: dead
         call alloc_cohort_diag(dead, 8_ik, .false.)
         call reorder_cohort_diag(dead, perm, 4_ik)      ! must not touch anything
         call cohort_diag_value(dead, CO_SAPFLOW, x, n, 0.0_wp)
         call check(n == 0_ik, 'inactive block reports no data and does not crash')
      end block
   end subroutine test_diag_lockstep

end program test_containers
