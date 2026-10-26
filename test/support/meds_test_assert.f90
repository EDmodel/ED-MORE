! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_test_assert -- the shared assertion helpers for the CTest programs.                 !
!                                                                                          !
! SEPARATE from meds_test_support on purpose. `build_test_config` needs the config, PFT and !
! demography layers; assertions need nothing but a kind. Many tests deliberately link one     !
! narrow library each (`meds_share`, `meds_config`, `meds_fast_kernels`, `meds_met`) to    !
! keep those layers standalone-buildable, and they must be able to assert without dragging     !
! the whole model in. That is why this module exists and why it uses meds_kinds only.           !
!                                                                                          !
! TWO assertion families, deliberately. They differ in FAILURE BEHAVIOUR, not just signature:  !
!                                                                                          !
!   * FATAL, condition-first -- check(cond, msg) / check_close(got, expect, rtol, msg).         !
!     Stops at the first failure. Right when later assertions are meaningless once an early      !
!     one fails (a roundtrip that did not round-trip).                                            !
!   * ACCUMULATING, name-first -- check(name, got, expect, atol), check_true, check_close,         !
!     check_int. Counts failures and keeps going, so one run reports EVERY failure. Right for a     !
!     sweep of independent cases. The program ends with test_report(), which error-stops on a       !
!     nonzero count.                                                                                 !
!                                                                                          !
! `check` and `check_close` are GENERIC over the two families: the first dummy is logical/real   !
! in the fatal form and character in the accumulating form, which is what makes them              !
! distinguishable.                                                                         !
!==========================================================================================!
module meds_test_assert
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: check, check_close, check_true, check_int, banner
   public :: test_report, test_reset, test_failures

   !----- Failure counter for the ACCUMULATING family. Module state on purpose: it is what      !
   !      lets a test call the helpers from any internal procedure without threading a counter   !
   !      through every one. --------------------------------------------------------------!
   integer(ik) :: n_fail = 0_ik

   interface check
      module procedure check_cond          !< (cond, msg)                -- fatal
      module procedure check_named         !< (name, got, expect, atol)  -- accumulating
   end interface check

   interface check_close
      module procedure check_close_rtol    !< (got, expect, rtol, msg)     -- fatal
      module procedure check_close_named   !< (name, got, expect [, atol]) -- accumulating
   end interface check_close

contains

   !=======================================================================================!
   !  FATAL family -- stops at the first failure.                                            !
   !=======================================================================================!
   subroutine check_cond(cond, msg)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) then
         write(*,'(2a)') 'FAIL: ', msg
         error stop 1
      end if
   end subroutine check_cond

   subroutine check_close_rtol(a, b, rtol, msg)
      real(wp),         intent(in) :: a, b, rtol
      character(len=*), intent(in) :: msg
      real(wp) :: tol
      tol = rtol * max(abs(b), 1.0e-30_wp) + 1.0e-12_wp
      if (abs(a - b) > tol) then
         write(*,'(2a)') 'FAIL: ', msg
         write(*,'(a,es16.8,a,es16.8,a,es10.2)') '   got=', a, ' expected=', b, ' rtol=', rtol
         error stop 1
      end if
   end subroutine check_close_rtol

   !=======================================================================================!
   !  ACCUMULATING family -- counts failures, reports every one, error-stops in test_report. !
   !  The ok line carries the value; the FAIL line carries got, expected and the tolerance !
   !  that was breached.                                                                   !
   !=======================================================================================!
   subroutine check_named(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5,a)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         n_fail = n_fail + 1_ik
         print '(a,a,a,es13.5,a,es13.5,a,es10.2)', '  FAIL : ', name, '  got ', got,          &
               ' expected ', expect, '  |diff|>', atol
      end if
   end subroutine check_named

   !----- atol ABSENT: 1e-9 relative, with a 1e-9 absolute floor where |expect| < 1. ------!
   subroutine check_close_named(name, got, expect, atol)
      character(len=*),   intent(in) :: name
      real(wp),           intent(in) :: got, expect
      real(wp), optional, intent(in) :: atol
      real(wp) :: tol
      tol = 1.0e-9_wp * max(1.0_wp, abs(expect)) ; if (present(atol)) tol = atol
      call check_named(name, got, expect, tol)
   end subroutine check_close_named

   subroutine check_true(name, cond, val)
      character(len=*),   intent(in) :: name
      logical,            intent(in) :: cond
      real(wp), optional, intent(in) :: val
      character(len=:), allocatable  :: tag
      tag = '  ok   : ' ; if (.not. cond) tag = '  FAIL : '
      if (.not. cond) n_fail = n_fail + 1_ik
      if (present(val)) then
         print '(a,a,a,es13.5,a)', tag, name, '  (', val, ')'
      else
         print '(a,a)', tag, name
      end if
   end subroutine check_true

   subroutine check_int(name, got, expect)
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: got, expect
      if (got == expect) then
         print '(a,a)', '  ok   : ', name
      else
         n_fail = n_fail + 1_ik
         print '(a,a,i0,a,i0)', '  FAIL : ', name, got, ' expected ', expect
      end if
   end subroutine check_int

   !----- End-of-program verdict for the accumulating family. error stop here and NOT at each   !
   !      failure is the whole point: one run lists every broken case. --------------------------!
   subroutine test_report(name)
      character(len=*), intent(in) :: name
      if (n_fail == 0_ik) then
         print '(2a)', trim(name), ': ALL PASSED'
      else
         print '(2a,i0,a)', trim(name), ': ', n_fail, ' FAILED'
         error stop 1
      end if
   end subroutine test_report

   subroutine test_reset()
      n_fail = 0_ik
   end subroutine test_reset

   integer(ik) function test_failures()
      test_failures = n_fail
   end function test_failures

   subroutine banner(name)
      character(len=*), intent(in) :: name
      write(*,'(2a)') '[test] ', name
   end subroutine banner
end module meds_test_assert
