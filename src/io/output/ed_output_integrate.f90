! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_integrate -- stage [4] of the diagnostic output, the temporal reduction: fold a
! variable's value (ed_output_extract) into its (variable, tier) running reduction
! (integrate_scalar, integrate_slab), close a period (normalize_*,
! close_tier) and re-seed (reset_buffer). output_integrate is the per-step tick and
! output_integrate_fast the FAST tier's replay. netCDF-free. Each tier integrates raw state
! on its own, which for these operators equals chaining (docs/science/diagnostics.md §1, §4).
!---~---
module ed_output_integrate
   use ed_kinds,            only : wp, ik
   use ed_time,             only : ed_time_t, calendar_step_t
   use ed_output_params,      only : N_FREQ, FREQ_FAST, TIER_FAST, TIER_DAILY, TIER_MONTHLY, TIER_ANNUAL
   use ed_output_types,     only : var_desc_t, integ_buffer_t, output_files_t, output_buffers_t,                  &
                                   pending_record_t, record_queue_t, AGG_MEAN, AGG_RATE, AGG_LAST, AGG_TMEAN,     &
                                   DIM_SCALAR, DIM_COHORT, DIM_PATCH,                                             &
                                   DIM_SOIL, DIM_SOIL_PATCH, MISSING_VALUE, N_WEIGHT, slab_col, slab_var
   use ed_polygon_state_vars, only : polygon_t
   use ed_constants,     only : n_soil_layer_max
   use ed_diag_vars,          only : N_CO_FAST_ROWS
   use ed_output_extract,   only : SRC_F_CO0, SRC_F_COH_GPP, SRC_F_COH_HEIGHT, SRC_F_COH_NPP, SRC_F_PA0,           &
                                   SRC_F_SOIL_LAYER, SRC_F_SOIL_WATER, SRC_NONE, extract_fast_scalar,             &
                                   extract_variable
   implicit none
   private

   public :: alloc_integ_buffer, reset_buffer, integrate_scalar, integrate_slab
   public :: normalize_scalar, normalize_slab
   public :: output_tick, output_integrate, output_integrate_fast, close_tier
   public :: write_order, find_ids

contains

   !---~---
   !   A polygon's output tick after a slow step: replay the FAST samples its fast loop
   ! staged, closing that tier every fast_interval_steps of them over the run, and fold each
   ! into its bin of the mean diel cycles; then fold the step into the slower tiers
   ! (output_integrate), whose snapshots of fast-loop values read the last staged sub-step.
   ! Closed records are queued for the I/O phase.
   !---~---
   subroutine output_tick(files, bufs, poly, step, dt_slow, dt_fast)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(in)    :: poly
      type(calendar_step_t),  intent(in)    :: step
      real(wp),               intent(in)    :: dt_slow, dt_fast   ! [s] the slow step, a fast sub-step
      integer(ik) :: isub
      if (bufs%fast_ready .and. files%reg%nidx(TIER_FAST) > 0_ik) then
         do isub = 1_ik, bufs%n_fast_sub
            call output_integrate_fast(files, bufs, isub, dt_fast, poly)
            bufs%fast_step_total = bufs%fast_step_total + 1_ik
            if (mod(bufs%fast_step_total, max(files%fast_interval_steps, 1_ik)) == 0_ik)          &
               call close_tier(files, bufs, TIER_FAST)
         end do
      end if
      !---~---
      !   The sub-steps belong to the step output_integrate folds next, so they go into the
      ! windows it folds into, before it closes the periods the step's end has left.
      !---~---
      if (bufs%fast_ready .and. files%diel_live) then
         do isub = 1_ik, bufs%n_fast_sub
            call output_integrate_diel(files, bufs, poly, isub, dt_fast)
         end do
      end if
      call output_integrate(files, bufs, poly, step%start, dt_slow, step%new_day, step%new_month,     &
                            step%new_year)
      bufs%fast_ready = .false.
   end subroutine output_tick

   !--- Buffer lifecycle

   !---~---
   !   Allocate + seed one (variable, tier) buffer. slab arrays sized to `cap` for a
   ! non-scalar dim; scalar dims allocate none.
   !---~---
   subroutine alloc_integ_buffer(buf, v, freq, cap)
      type(integ_buffer_t), intent(out) :: buf
      type(var_desc_t),     intent(in)  :: v
      integer(ik),          intent(in)  :: freq, cap
      buf%var_id = 0_ik      ! set by the caller (registry index)
      buf%freq   = freq
      buf%agg    = v%agg
      !---~---
      !   The F tier folds staged sub-step samples, rates and states alike: a mean variable their
      ! dt-weighted mean, a snapshot the last.
      !---~---
      if ((freq == FREQ_FAST .or. v%diel) .and. v%agg /= AGG_LAST) buf%agg = AGG_TMEAN
      buf%dim    = v%dim
      buf%active = .true.
      if (slab_var(v)) then
         allocate(buf%slab(cap), buf%wsum_slab(cap), buf%hits(cap))
      end if
      call reset_buffer(buf)
   end subroutine alloc_integ_buffer

   !--- Re-seed a buffer for the next period.
   subroutine reset_buffer(buf)
      type(integ_buffer_t), intent(inout) :: buf
      buf%scal  = 0.0_wp
      buf%wsum  = 0.0_wp
      buf%nsamp = 0_ik
      buf%n_slab = 0_ik
      if (allocated(buf%slab)) then
         buf%slab(:)      = 0.0_wp
         buf%wsum_slab(:) = 0.0_wp
         buf%hits(:)      = 0_ik
      end if
   end subroutine reset_buffer


   !--- Integrate (fold one step)

   !--- Fold one scalar sample x over a step of length dt [s] into a scalar buffer.
   elemental subroutine integrate_scalar(buf, x, dt)
      type(integ_buffer_t), intent(inout) :: buf
      real(wp),             intent(in)    :: x, dt
      if (x == MISSING_VALUE) return        ! nothing to report this step (see not_simulated)
      select case (buf%agg)
      case (AGG_MEAN)
         buf%scal = buf%scal + x
         buf%nsamp = buf%nsamp + 1_ik
      case (AGG_RATE)
         buf%scal = buf%scal + x
         buf%wsum  = buf%wsum + dt
      case (AGG_LAST)
         buf%scal = x
         buf%nsamp = buf%nsamp + 1_ik
      case (AGG_TMEAN)
         buf%scal = buf%scal + x*dt
         buf%wsum  = buf%wsum + dt
      end select
   end subroutine integrate_scalar

   !---~---
   !   Fold a slab of n samples by direct slot index. `valid` (optional) marks the slots with
   ! something to report this step; the rest (an empty patch, PFT or class, a layer below the
   ! column) are skipped and close as _FillValue, not as a 0 a reader would take for a
   ! measurement (docs/science/diagnostics.md §2). Without `valid` every slot is folded.
   !---~---
   subroutine integrate_slab(buf, x, n, dt, valid)
      type(integ_buffer_t), intent(inout) :: buf
      real(wp),             intent(in)    :: x(:)
      integer(ik),          intent(in)    :: n
      real(wp),             intent(in)    :: dt
      logical, optional,    intent(in)    :: valid(:)
      integer(ik) :: i
      logical     :: use_v
      buf%n_slab = max(buf%n_slab, n)
      use_v = present(valid)
      do i = 1_ik, n
         if (use_v) then
            if (.not. valid(i)) cycle
         end if
         select case (buf%agg)
         case (AGG_MEAN)
            buf%slab(i) = buf%slab(i) + x(i)
            buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_RATE)
            buf%slab(i) = buf%slab(i) + x(i)
            buf%wsum_slab(i) = buf%wsum_slab(i) + dt
            buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_LAST)
            buf%slab(i) = x(i)
            buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_TMEAN)
            buf%slab(i) = buf%slab(i) + x(i)*dt
            buf%wsum_slab(i) = buf%wsum_slab(i) + dt
            buf%hits(i) = buf%hits(i) + 1_ik
         end select
      end do
   end subroutine integrate_slab

   !---~---
   !   Fold a slab of n samples into bin `bin` of each slot of a mean diel cycle, dt-weighted:
   ! slot i's nb bins are entries (i-1)*nb+1 .. i*nb, the bin varying fastest. Invalid slots
   ! are skipped, as in integrate_slab.
   !---~---
   subroutine integrate_diel(buf, x, valid, n, bin, nb, dt)
      type(integ_buffer_t), intent(inout) :: buf
      real(wp),             intent(in)    :: x(:)
      logical,              intent(in)    :: valid(:)
      integer(ik),          intent(in)    :: n, bin, nb
      real(wp),             intent(in)    :: dt
      integer(ik) :: i, m
      buf%n_slab = max(buf%n_slab, n * nb)
      do i = 1_ik, n
         if (.not. valid(i)) cycle
         m = (i - 1_ik) * nb + bin
         buf%slab(m)      = buf%slab(m) + x(i) * dt
         buf%wsum_slab(m) = buf%wsum_slab(m) + dt
         buf%hits(m)      = buf%hits(m) + 1_ik
      end do
   end subroutine integrate_diel

   !--- Normalize (close a period). valid=.false. -> emit MISSING_VALUE (no NaN, no ±huge).
   subroutine normalize_scalar(buf, out, valid)
      type(integ_buffer_t), intent(in)  :: buf
      real(wp),             intent(out) :: out
      logical,              intent(out) :: valid
      valid = .true.
      out = MISSING_VALUE
      select case (buf%agg)
      case (AGG_MEAN)
         if (buf%nsamp > 0_ik) then
            out = buf%scal / real(buf%nsamp, wp)
         else
            valid = .false.
         end if
      case (AGG_TMEAN, AGG_RATE)
         if (buf%wsum > 0.0_wp) then
            out = buf%scal / buf%wsum
         else
            valid = .false.
         end if
      case (AGG_LAST)
         if (buf%nsamp > 0_ik) then
            out = buf%scal
         else
            valid = .false.
         end if
      end select
   end subroutine normalize_scalar

   !--- Slab normalize: one value + validity per live slot [1:n_out].
   subroutine normalize_slab(buf, out, valid, n_out)
      type(integ_buffer_t), intent(in)  :: buf
      real(wp),             intent(out) :: out(:)
      logical,              intent(out) :: valid(:)
      integer(ik),          intent(out) :: n_out
      integer(ik) :: i
      n_out = buf%n_slab
      !---~---
      !   A record row per live slot. A shorter record means the slab was sized before the
      ! registry was final, and the loop below would write past this column into the next.
      !---~---
      if (n_out > size(out, kind=ik))                                                            &
         error stop 'ed_output_integrate: a slab is longer than its record (max_slab sized too early)'
      out(:)   = MISSING_VALUE
      valid(:) = .false.
      do i = 1_ik, n_out
         select case (buf%agg)
         case (AGG_MEAN)
            if (buf%hits(i) > 0_ik) then
               out(i) = buf%slab(i) / real(buf%hits(i), wp)
               valid(i) = .true.
            end if
         case (AGG_TMEAN, AGG_RATE)
            if (buf%wsum_slab(i) > 0.0_wp) then
               out(i) = buf%slab(i) / buf%wsum_slab(i)
               valid(i) = .true.
            end if
         case (AGG_LAST)
            if (buf%hits(i) > 0_ik) then
               out(i) = buf%slab(i)
               valid(i) = .true.
            end if
         end select
      end do
   end subroutine normalize_slab

   !---~---
   !   The per-step tick. The slow step just taken covers [prev, now): it is folded into the
   ! window that holds prev, and only then is each window that now has left closed, so a
   ! step belongs to the period it starts in (docs/science/diagnostics.md §4, Which steps
   ! a record holds). Closed records go to bufs%queue. The calendar restructuring runs
   ! after this tick, at the next step's start (ed_stepper).
   !---~---
   subroutine output_integrate(files, bufs, poly, prev, dt, is_new_day, is_new_month, is_new_year)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(in)    :: poly
      type(ed_time_t),        intent(in)    :: prev          ! the start of the step just taken
      real(wp),               intent(in)    :: dt
      logical,                intent(in)    :: is_new_day, is_new_month, is_new_year  ! its end left prev's period
      integer(ik) :: t, j, k, n_out, nc
      real(wp)    :: scal
      real(wp)    :: slab(max(files%max_slab, 1_ik))
      logical     :: vslab(max(files%max_slab, 1_ik))
      integer(ik) :: at_live(max(poly%cohort%n, 1_ik))    ! each live cohort's slot in the month's order
      integer(ik) :: at_staged(max(bufs%fast_n_cohort, 1_ik))
      if (.not. files%enabled) return
      !---~---
      !   The cohort and patch slabs are output.cohort_max and patch_max long. Stop at the step the
      ! stand outgrows one, not when its record is written, up to a month later.
      !---~---
      if (files%cohort_axis .and. poly%cohort%n > files%cohort_max) call cap_exceeded('cohort', poly%cohort%n, files%cohort_max)
      if (files%patch_axis  .and. poly%patch%n  > files%patch_max)  call cap_exceeded('patch',  poly%patch%n,  files%patch_max)
      !---~---
      !   A cohort value is folded into its cohort's slot of the month (month_slots), from the
      ! live stand's order or the staging's, and the live order is kept as the order the
      ! month's records are written in if the month ends here (write_order).
      !---~---
      if (files%cohort_axis) then
         nc = poly%cohort%n
         call month_slots(files, bufs, poly%cohort%global_id(1:nc), at_live)
         if (bufs%fast_ready) call month_slots(files, bufs, bufs%fast_cohort_id(1:bufs%fast_n_cohort), at_staged)
         bufs%end_id(1:nc) = poly%cohort%global_id(1:nc)
         bufs%n_end = nc
      end if
      !---~---
      !   Fold the step into the daily, monthly and annual tiers. The FAST tier is fed from the
      ! staged sub-steps (output_integrate_fast).
      !---~---
      do t = TIER_DAILY, N_FREQ
         if (files%reg%nidx(t) == 0_ik) cycle
         if (.not. bufs%has_data(t)) bufs%t_open(t) = prev
         do j = 1_ik, files%reg%nidx(t)
            k = files%reg%idx_freq(j, t)
            if (files%reg%var(k)%diel) cycle               ! folded from the sub-steps (output_integrate_diel)
            if (files%reg%var(k)%snap_staged) then
               call staged_snapshot(bufs, poly, files%reg%var(k), scal, slab, vslab, n_out)
            else
               call extract_variable(poly, bufs%diag, files%reg%var(k), scal, slab, vslab, n_out)
            end if
            if (files%reg%var(k)%dim == DIM_SCALAR) then
               call integrate_scalar(bufs%buf(k,t), scal, dt)
            else if (files%reg%var(k)%dim == DIM_COHORT) then
               if (files%reg%var(k)%snap_staged .and. files%reg%var(k)%source_fast /= SRC_NONE) then
                  call integrate_cohort_slab(bufs, bufs%buf(k,t), slab, vslab, n_out, at_staged, dt)
               else
                  call integrate_cohort_slab(bufs, bufs%buf(k,t), slab, vslab, n_out, at_live, dt)
               end if
            else
               call integrate_slab(bufs%buf(k,t), slab, n_out, dt, vslab)
            end if
         end do
         bufs%has_data(t) = .true.
      end do
      !--- then close the periods the step's end has left (tiers independent).
      if (is_new_year  .and. bufs%has_data(TIER_ANNUAL))  call close_tier(files, bufs, TIER_ANNUAL)
      if (is_new_month .and. bufs%has_data(TIER_MONTHLY)) call close_tier(files, bufs, TIER_MONTHLY)
      if (is_new_day   .and. bufs%has_data(TIER_DAILY))   call close_tier(files, bufs, TIER_DAILY)
      !---~---
      !   The month's records are closed: the calendar restructuring that follows may recruit,
      ! fuse, split or cull cohorts, so the next month opens its own order.
      !---~---
      if (is_new_month) bufs%n_month = 0_ik
   end subroutine output_integrate

   !---~---
   !   Fold a cohort slab, in the order `at` maps it from (live or staged), into the slots of
   ! the month's order: slot at(i) takes entry i, and a slot nothing maps to is not folded
   ! this step.
   !---~---
   subroutine integrate_cohort_slab(bufs, buf, x, valid, n, at, dt)
      type(output_buffers_t), intent(in)    :: bufs
      type(integ_buffer_t),   intent(inout) :: buf
      real(wp),               intent(in)    :: x(:)
      logical,                intent(in)    :: valid(:)
      integer(ik),            intent(in)    :: n, at(:)
      real(wp),               intent(in)    :: dt
      real(wp)    :: xm(max(bufs%n_month, 1_ik))
      logical     :: vm(max(bufs%n_month, 1_ik))
      integer(ik) :: i
      xm = MISSING_VALUE
      vm = .false.
      do i = 1_ik, n
         xm(at(i)) = x(i)
         vm(at(i)) = valid(i)
      end do
      call integrate_slab(buf, xm, bufs%n_month, dt, vm)
   end subroutine integrate_cohort_slab

   !---~---
   !   Each of `ids`' slot in the month's cohort order, at(i). The month opens with the first
   ! ids it sees; a cohort it has not seen (none in a month the calendar restructures only at
   ! its turn) takes the next slot.
   !---~---
   subroutine month_slots(files, bufs, ids, at)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: ids(:)
      integer(ik),            intent(out)   :: at(:)
      integer(ik) :: i
      call find_ids(bufs%month_id(1:bufs%n_month), ids, at)
      do i = 1_ik, size(ids, kind=ik)
         if (at(i) > 0_ik) cycle
         bufs%n_month = bufs%n_month + 1_ik
         if (bufs%n_month > files%cohort_max) call cap_exceeded('cohort', bufs%n_month, files%cohort_max)
         bufs%month_id(bufs%n_month) = ids(i)
         at(i) = bufs%n_month
      end do
   end subroutine month_slots

   !---~---
   !   Where each of `ids` sits in `ref`: at(i) = j with ref(j) == ids(i), or 0 when it is not
   ! there. A re-sort moves a cohort only within its patch, so the search starts at i and
   ! widens outward: a stand in the same order costs one comparison per cohort.
   !---~---
   pure subroutine find_ids(ref, ids, at)
      integer(ik), intent(in)  :: ref(:), ids(:)
      integer(ik), intent(out) :: at(:)
      integer(ik) :: i, d, nr
      nr = size(ref, kind=ik)
      do i = 1_ik, size(ids, kind=ik)
         at(i) = 0_ik
         do d = 0_ik, max(nr, i)
            if (i + d <= nr) then
               if (ref(i + d) == ids(i)) then
                  at(i) = i + d
                  exit
               end if
            end if
            if (d > 0_ik .and. i - d >= 1_ik .and. i - d <= nr) then
               if (ref(i - d) == ids(i)) then
                  at(i) = i - d
                  exit
               end if
            end if
            if (i + d > nr .and. i - d < 1_ik) exit
         end do
      end do
   end subroutine find_ids

   !--- A live count past its output cap: say which, by how much, and what to raise.
   subroutine cap_exceeded(axis, n, cap)
      character(len=*), intent(in) :: axis
      integer(ik),      intent(in) :: n, cap
      write(*,'(5a,i0,a,i0,a)') ' output: the stand has more ', axis, 's than output.', axis, '_max (', n, ' > ', &
                                cap, '); raise it in [output]'
      error stop 'output: a live count exceeds its output cap (see the message above)'
   end subroutine cap_exceeded

   !--- Normalize a tier's buffers into its pending record + reset them (staging).
   subroutine close_tier(files, bufs, t)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: t
      integer(ik) :: j, k, c, ns, nsl, nb
      integer(ik) :: slab_k(files%reg%nidx(t))
      integer(ik) :: patch_rows(files%reg%nidx(t))
      bufs%pending(t)%used   = .true.
      bufs%pending(t)%freq   = ishft(1_ik, t - 1_ik)
      bufs%pending(t)%t_open = bufs%t_open(t)
      bufs%pending(t)%n_cohort = 0_ik
      bufs%pending(t)%n_patch = 0_ik
      bufs%pending(t)%sval(:)   = MISSING_VALUE
      bufs%pending(t)%svalid(:) = .false.
      bufs%pending(t)%nslab(:)  = 0_ik
      do j = 1_ik, files%reg%nidx(t)
         k = files%reg%idx_freq(j, t)
         if (.not. slab_var(files%reg%var(k))) then
            call normalize_scalar(bufs%buf(k,t), bufs%pending(t)%sval(k), bufs%pending(t)%svalid(k))
         else
            c = slab_col(bufs%pending(t), k)
            call normalize_slab(bufs%buf(k,t), bufs%pending(t)%slab(:,c),                          &
                                bufs%pending(t)%slabvalid(:,c), ns)
            bufs%pending(t)%nslab(k) = ns
            !---~---
            !   The record's axis lengths. A soil-by-patch slab is patch-major, strided by the
            ! layer ceiling, so it holds ns / n_soil_layer_max patches; a mean diel cycle holds
            ! n_diel entries a slot.
            !---~---
            nb = merge(files%n_diel, 1_ik, files%reg%var(k)%diel)
            select case (files%reg%var(k)%dim)
            case (DIM_COHORT)
               bufs%pending(t)%n_cohort = max(bufs%pending(t)%n_cohort, ns / nb)
            case (DIM_PATCH)
               bufs%pending(t)%n_patch  = max(bufs%pending(t)%n_patch,  ns / nb)
            case (DIM_SOIL_PATCH)
               bufs%pending(t)%n_patch  = max(bufs%pending(t)%n_patch,  ns / (n_soil_layer_max * nb))
            end select
         end if
         call reset_buffer(bufs%buf(k,t))
      end do
      bufs%has_data(t) = .false.
      !--- The record's cohort rows are in the month's order: keep their ids (write_order).
      if (bufs%pending(t)%n_cohort > 0_ik) then
         if (.not. allocated(bufs%pending(t)%cohort_id)) allocate(bufs%pending(t)%cohort_id(files%cohort_max))
         bufs%pending(t)%cohort_id(1:bufs%pending(t)%n_cohort) = bufs%month_id(1:bufs%pending(t)%n_cohort)
      end if
      nsl = 0_ik                                          ! this tier's slab variables, for the queue
      do j = 1_ik, files%reg%nidx(t)
         k = files%reg%idx_freq(j, t)
         if (.not. slab_var(files%reg%var(k))) cycle
         nsl = nsl + 1_ik
         slab_k(nsl) = k
         patch_rows(nsl) = 0_ik
         if (files%reg%var(k)%dim == DIM_SOIL_PATCH)                                              &
            patch_rows(nsl) = n_soil_layer_max * merge(files%n_diel, 1_ik, files%reg%var(k)%diel)
      end do
      call enqueue_record(bufs%queue(t), bufs%pending(t), slab_k(1:nsl), patch_rows(1:nsl))
      bufs%pending(t)%used = .false.
   end subroutine close_tier

   !---~---
   !   Append a copy of the scratch record to the tier's queue. Only what the writer reads is
   ! copied: the tier's slab variables (`slab_k`), each to its nslab rows, or the whole patch x
   ! layer block (n_patch * patch_rows: the layers, times the bins of a mean diel cycle) for a
   ! soil-patch variable. The copy is exact, so a queued record writes exactly what an
   ! immediate write from the scratch would have.
   !---~---
   subroutine enqueue_record(q, src, slab_k, patch_rows)
      type(record_queue_t),   intent(inout) :: q
      type(pending_record_t), intent(in)    :: src
      integer(ik),            intent(in)    :: slab_k(:)
      integer(ik),            intent(in)    :: patch_rows(:)   ! rows a patch takes; 0 off the soil-patch axis
      type(pending_record_t), allocatable :: grown(:)
      integer(ik) :: rows, nv, ncol, c
      if (.not. allocated(q%rec)) allocate(q%rec(8))
      if (q%n == size(q%rec, kind=ik)) then                   ! double the capacity (rare: once per run)
         allocate(grown(2 * size(q%rec)))
         grown(1:q%n) = q%rec(1:q%n)
         call move_alloc(grown, q%rec)
      end if
      q%n = q%n + 1_ik
      nv = size(src%sval, kind=ik)
      ncol = size(slab_k, kind=ik)
      rows = 1_ik
      do c = 1_ik, ncol
         rows = max(rows, src%nslab(slab_k(c)))
         rows = max(rows, src%n_patch * patch_rows(c))
      end do
      rows = min(rows, size(src%slab, 1, kind=ik))
      associate (r => q%rec(q%n))
         r%used = .true.
         r%freq = src%freq
         r%t_open = src%t_open
         r%n_cohort = src%n_cohort
         r%n_patch = src%n_patch
         if (allocated(r%sval)) then
            if (size(r%sval) /= nv .or. size(r%slab, 1) < rows .or. size(r%slab, 2) < max(ncol, 1_ik)) &
               deallocate(r%sval, r%svalid, r%nslab, r%col, r%slab, r%slabvalid)
         end if
         if (.not. allocated(r%sval))                                                             &
            allocate(r%sval(nv), r%svalid(nv), r%nslab(nv), r%col(nv),                            &
                     r%slab(rows, max(ncol, 1_ik)), r%slabvalid(rows, max(ncol, 1_ik)))
         if (src%n_cohort > 0_ik) then
            if (allocated(r%cohort_id)) then
               if (size(r%cohort_id) < src%n_cohort) deallocate(r%cohort_id)
            end if
            if (.not. allocated(r%cohort_id)) allocate(r%cohort_id(src%n_cohort))
            r%cohort_id(1:src%n_cohort) = src%cohort_id(1:src%n_cohort)
         end if
         r%sval = src%sval
         r%svalid = src%svalid
         r%nslab = src%nslab
         r%col = 0_ik
         do c = 1_ik, ncol
            r%col(slab_k(c)) = c
            r%slab(1:rows, c) = src%slab(1:rows, slab_col(src, slab_k(c)))
            r%slabvalid(1:rows, c) = src%slabvalid(1:rows, slab_col(src, slab_k(c)))
         end do
      end associate
   end subroutine enqueue_record

   !---~---
   !   Fold one staged sub-step (isub) into the FAST tier's buffers. output_tick calls it
   ! once per staged sub-step and closes the tier every fast_interval_steps folds. It reads
   ! the buffers' staging, never the polygon's live state (the end of the slow step by now),
   ! but for the axes' bookkeeping, which the slow step leaves as it was but for the order,
   ! and which is folded by cohort id like every cohort value (`poly`; without it, not read).
   !---~---
   subroutine output_integrate_fast(files, bufs, isub, dt, poly)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: isub
      real(wp),               intent(in)    :: dt
      type(polygon_t), optional, intent(in) :: poly
      integer(ik) :: j, k, n
      real(wp)    :: scal
      real(wp)    :: slab(max(files%max_slab, 1_ik))
      logical     :: valid(max(files%max_slab, 1_ik))
      integer(ik) :: at_staged(max(bufs%fast_n_cohort, 1_ik))   ! each staged cohort's slot in the month's order
      integer(ik), allocatable :: at_live(:)                    ! each live cohort's (the bookkeeping)
      if (.not. files%enabled) return
      if (files%reg%nidx(TIER_FAST) == 0_ik) return
      !---~---
      !   The staging follows the stand; the slabs are output.cohort_max and patch_max long. The
      ! same check as the coarse tiers' (output_integrate), before this tier writes into them.
      !---~---
      if (files%cohort_axis .and. bufs%fast_n_cohort > files%cohort_max)                           &
         call cap_exceeded('cohort', bufs%fast_n_cohort, files%cohort_max)
      if (files%patch_axis .and. bufs%fast_n_patch > files%patch_max)                              &
         call cap_exceeded('patch', bufs%fast_n_patch, files%patch_max)
      if (files%cohort_axis) then
         call month_slots(files, bufs, bufs%fast_cohort_id(1:bufs%fast_n_cohort), at_staged)
         if (present(poly)) then
            allocate(at_live(max(poly%cohort%n, 1_ik)))
            call month_slots(files, bufs, poly%cohort%global_id(1:poly%cohort%n), at_live)
         end if
      end if
      if (.not. bufs%has_data(TIER_FAST)) bufs%t_open(TIER_FAST) = bufs%fast_time(isub)
      do j = 1_ik, files%reg%nidx(TIER_FAST)
         k = files%reg%idx_freq(j, TIER_FAST)
         if (files%reg%var(k)%coord) then
            if (.not. present(poly)) cycle
            call extract_variable(poly, bufs%diag, files%reg%var(k), scal, slab, valid, n)
            if (files%reg%var(k)%dim == DIM_COHORT) then
               call integrate_cohort_slab(bufs, bufs%buf(k,TIER_FAST), slab, valid, n, at_live, dt)
            else
               call integrate_slab(bufs%buf(k,TIER_FAST), slab, n, dt, valid)
            end if
            cycle
         end if
         call staged_sample(bufs, files%reg%var(k), isub, scal, slab, valid, n, poly)
         if (files%reg%var(k)%dim == DIM_SCALAR) then
            call integrate_scalar(bufs%buf(k,TIER_FAST), scal, dt)
         else if (files%reg%var(k)%dim == DIM_COHORT) then
            call integrate_cohort_slab(bufs, bufs%buf(k,TIER_FAST), slab, valid, n, at_staged, dt)
         else
            call integrate_slab(bufs%buf(k,TIER_FAST), slab, n, dt, valid)
         end if
      end do
      bufs%has_data(TIER_FAST) = .true.
   end subroutine output_integrate_fast

   !---~---
   !   Fold one staged sub-step (isub) into its bin of the coarse tiers' mean diel cycles (the
   ! bin of the hour the sub-step starts in), each read from the staging as the F tier reads
   ! it, a cohort value into its cohort's slot of the month.
   !---~---
   subroutine output_integrate_diel(files, bufs, poly, isub, dt)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(in)    :: poly
      integer(ik),            intent(in)    :: isub
      real(wp),               intent(in)    :: dt
      integer(ik) :: t, j, k, n, i, bin
      real(wp)    :: scal
      real(wp)    :: slab(max(files%max_slab, 1_ik))
      logical     :: valid(max(files%max_slab, 1_ik))
      real(wp)    :: xm(max(files%cohort_max, 1_ik))
      logical     :: vm(max(files%cohort_max, 1_ik))
      integer(ik) :: at_staged(max(bufs%fast_n_cohort, 1_ik))   ! each staged cohort's slot in the month's order
      if (.not. files%enabled) return
      if (files%cohort_axis) call month_slots(files, bufs, bufs%fast_cohort_id(1:bufs%fast_n_cohort), at_staged)
      bin = bufs%fast_time(isub)%hour / (24_ik / files%n_diel) + 1_ik
      do t = TIER_DAILY, N_FREQ
         do j = 1_ik, files%reg%nidx(t)
            k = files%reg%idx_freq(j, t)
            if (.not. files%reg%var(k)%diel) cycle
            call staged_sample(bufs, files%reg%var(k), isub, scal, slab, valid, n, poly)
            select case (files%reg%var(k)%dim)
            case (DIM_SCALAR)
               n = 1_ik
               slab(1)  = scal
               valid(1) = scal /= MISSING_VALUE
            case (DIM_COHORT)
               xm(1:bufs%n_month) = MISSING_VALUE
               vm(1:bufs%n_month) = .false.
               do i = 1_ik, n
                  xm(at_staged(i)) = slab(i)
                  vm(at_staged(i)) = valid(i)
               end do
               n = bufs%n_month
               slab(1:n)  = xm(1:n)
               valid(1:n) = vm(1:n)
            end select
            call integrate_diel(bufs%buf(k,t), slab, valid, n, bin, files%n_diel, dt)
         end do
      end do
   end subroutine output_integrate_diel

   !---~---
   !   One variable's staged value of sub-step isub on its own axis, as the F tier and the mean
   ! diel cycles read it. A per-cohort row reduced to a coarser axis and a soil-layer field are
   ! formed against the stand (`poly`); without it they have nothing to report.
   !---~---
   subroutine staged_sample(bufs, v, isub, scal, slab, valid, n, poly)
      type(output_buffers_t),    intent(in)  :: bufs
      type(var_desc_t),          intent(in)  :: v
      integer(ik),               intent(in)  :: isub
      real(wp),                  intent(out) :: scal, slab(:)
      logical,                   intent(out) :: valid(:)
      integer(ik),               intent(out) :: n
      type(polygon_t), optional, intent(in)  :: poly
      scal = MISSING_VALUE
      n = 0_ik
      if (reduced_cohort_row(v)) then
         if (present(poly)) call staged_cohort_reduced(bufs, poly, v, isub, scal, slab, valid, n)
      else if (v%source_fast == SRC_F_SOIL_LAYER) then
         if (present(poly)) call staged_soil_layer(bufs, poly, v, isub, slab, valid, n)
      else
         call staged_value(bufs, v%source_fast, v%dim, isub, scal, slab, valid, n)
      end if
   end subroutine staged_sample

   !---~---
   !   One staged source's value in sub-step isub, on the variable's axis: the polygon mean of
   ! the patch rows or a forcing row, each patch's row, a soil column, or the per-cohort
   ! values. Every slot is valid but a layer below its column. A per-cohort diagnostic row
   ! is staged only while the block is on (a fast loop that ran); else it has no slots.
   !---~---
   subroutine staged_value(bufs, src, dm, isub, scal, slab, valid, n)
      type(output_buffers_t), intent(in)  :: bufs
      integer(ik),            intent(in)  :: src, dm, isub
      real(wp),               intent(out) :: scal, slab(:)
      logical,                intent(out) :: valid(:)
      integer(ik),            intent(out) :: n
      scal = MISSING_VALUE
      n = 0_ik
      select case (dm)
      case (DIM_SCALAR)
         scal = extract_fast_scalar(src, bufs%fast_site(:, isub), bufs%fast_forcing(:, isub))
      case (DIM_PATCH)
         n = bufs%fast_n_patch
         slab(1:n)  = bufs%fast_patch(src - SRC_F_PA0, 1:n, isub)
         valid(1:n) = .true.
      case (DIM_SOIL)
         n = bufs%fast_n_site * n_soil_layer_max
         if (src == SRC_F_SOIL_WATER) then
            slab(1:n) = bufs%fast_soil_water(1:n, isub)
         else
            slab(1:n) = bufs%fast_soil_temp(1:n, isub)
         end if
         valid(1:n) = bufs%fast_soil_live(1:n)
      case (DIM_SOIL_PATCH)
         n = bufs%fast_n_patch * n_soil_layer_max
         if (src == SRC_F_SOIL_WATER) then
            slab(1:n) = bufs%fast_soil_water_patch(1:n, isub)
         else
            slab(1:n) = bufs%fast_soil_temp_patch(1:n, isub)
         end if
         valid(1:n) = bufs%fast_soil_patch_live(1:n)
      case (DIM_COHORT)
         n = bufs%fast_n_cohort
         select case (src)
         case (SRC_F_COH_GPP)
            slab(1:n) = bufs%fast_cohort_gpp(1:n, isub)
         case (SRC_F_COH_NPP)
            slab(1:n) = bufs%fast_cohort_npp(1:n, isub)
         case (SRC_F_COH_HEIGHT)
            slab(1:n) = bufs%fast_cohort_height(1:n, isub)
         case default   ! SRC_F_CO0 + CO_*: a row of the per-cohort diagnostic block
            if (allocated(bufs%fast_cohort_diag)) then
               slab(1:n) = bufs%fast_cohort_diag(1:n, src - SRC_F_CO0, isub)
            else
               n = 0_ik
            end if
         end select
         valid(1:n) = .true.
      end select
   end subroutine staged_value

   !---~---
   !   A snapshot of a fast-loop diagnostic: its value over the last fast step of the slow step
   ! just taken, the last staged sub-step (docs/science/diagnostics.md §3). A per-cohort row
   ! reduced to a coarser axis is that sub-step's row, reduced as its mean is. With nothing
   ! staged this step (no fast loop) it is missing, and the period keeps its last value.
   !---~---
   subroutine staged_snapshot(bufs, poly, v, scal, slab, valid, n)
      type(output_buffers_t), intent(in)  :: bufs
      type(polygon_t),        intent(in)  :: poly
      type(var_desc_t),       intent(in)  :: v
      real(wp),               intent(out) :: scal, slab(:)
      logical,                intent(out) :: valid(:)
      integer(ik),            intent(out) :: n
      scal = MISSING_VALUE
      n = 0_ik
      if (.not. bufs%fast_ready) return
      if (reduced_cohort_row(v)) then
         call staged_cohort_reduced(bufs, poly, v, bufs%n_fast_sub, scal, slab, valid, n)
      else if (v%source_fast /= SRC_NONE) then
         call staged_value(bufs, v%source_fast, v%dim, bufs%n_fast_sub, scal, slab, valid, n)
      end if
   end subroutine staged_snapshot

   !---~---
   !   A soil-layer field of sub-step isub (the matric potential, the relative saturation, the
   ! liquid fraction), formed per patch from the staged water content and liquid fraction as the
   ! coarse tiers form it from the live ones, then reduced to v's axis.
   !---~---
   subroutine staged_soil_layer(bufs, poly, v, isub, slab, valid, n)
      type(output_buffers_t), intent(in)  :: bufs
      type(polygon_t),        intent(in)  :: poly
      type(var_desc_t),       intent(in)  :: v
      integer(ik),            intent(in)  :: isub
      real(wp),               intent(out) :: slab(:)
      logical,                intent(out) :: valid(:)
      integer(ik),            intent(out) :: n
      real(wp)    :: theta(n_soil_layer_max, max(poly%patch%n, 1_ik)), fliq(n_soil_layer_max, max(poly%patch%n, 1_ik))
      real(wp)    :: scal
      integer(ik) :: ip, np
      n = 0_ik
      np = bufs%fast_n_patch
      if (np /= poly%patch%n .or. np < 1_ik) return
      do ip = 1_ik, np
         theta(:, ip) = bufs%fast_soil_water_patch((ip - 1_ik)*n_soil_layer_max + 1_ik:ip*n_soil_layer_max, isub)
         fliq(:, ip)  = bufs%fast_soil_fliq_patch((ip - 1_ik)*n_soil_layer_max + 1_ik:ip*n_soil_layer_max, isub)
      end do
      call extract_variable(poly, bufs%diag, v, scal, slab, valid, n, soil_theta=theta(:, 1:np), soil_fliq=fliq(:, 1:np))
   end subroutine staged_soil_layer

   !--- A per-cohort fast-loop row on a coarser axis: the patch, the polygon (staged_row).
   pure logical function reduced_cohort_row(v) result(yes)
      type(var_desc_t), intent(in) :: v
      yes = v%dim /= DIM_COHORT .and. v%source_fast > SRC_F_CO0 .and. v%source_fast <= SRC_F_CO0 + N_CO_FAST_ROWS
   end function reduced_cohort_row

   !---~---
   !   A staged per-cohort row of sub-step isub, reduced to v's axis as its mean is. The staged
   ! row is in the fast loop's order: each value goes with its own cohort in the live order,
   ! and is reduced with the weights the fast loop used. With no staged cohort block, or a
   ! stand that does not match the staging, nothing is reported.
   !---~---
   subroutine staged_cohort_reduced(bufs, poly, v, isub, scal, slab, valid, n)
      type(output_buffers_t), intent(in)  :: bufs
      type(polygon_t),        intent(in)  :: poly
      type(var_desc_t),       intent(in)  :: v
      integer(ik),            intent(in)  :: isub
      real(wp),               intent(out) :: scal, slab(:)
      logical,                intent(out) :: valid(:)
      integer(ik),            intent(out) :: n
      integer(ik) :: nc
      integer(ik) :: at(max(poly%cohort%n, 1_ik))          ! each live cohort's slot in the staging
      real(wp)    :: x(max(poly%cohort%n, 1_ik)), w(max(poly%cohort%n, 1_ik), N_WEIGHT)
      scal = MISSING_VALUE
      n = 0_ik
      if (.not. allocated(bufs%fast_cohort_diag)) return
      nc = bufs%fast_n_cohort
      if (nc /= poly%cohort%n) return
      call find_ids(bufs%fast_cohort_id(1:nc), poly%cohort%global_id(1:nc), at(1:nc))
      if (.not. all(at(1:nc) > 0_ik)) return
      x(1:nc)    = bufs%fast_cohort_diag(at(1:nc), v%source_fast - SRC_F_CO0, isub)
      w(1:nc, :) = bufs%fast_cohort_w(at(1:nc), :)
      call extract_variable(poly, bufs%diag, v, scal, slab, valid, n, cohort_x=x(1:nc), cohort_w=w(1:nc, :))
   end subroutine staged_cohort_reduced

   !---~---
   !   Put every queued cohort record of each polygon in its month's end order, the live order
   ! of the latest tick (end_id), before the I/O phase writes it (docs/science/diagnostics.md
   ! §3). A record's cohort ids name its rows; a cohort the end order lacks keeps its place
   ! after those it holds.
   !---~---
   subroutine write_order(files, bufs)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      integer(ik) :: p, t, i, j, k, n, c
      integer(ik), allocatable :: at(:), key(:), perm(:)
      do p = 1_ik, size(bufs, kind=ik)
         do t = 1_ik, N_FREQ
            do i = 1_ik, bufs(p)%queue(t)%n
               associate (r => bufs(p)%queue(t)%rec(i))
                  n = r%n_cohort
                  if (n <= 0_ik .or. .not. allocated(r%cohort_id)) cycle
                  if (n == bufs(p)%n_end) then
                     if (all(r%cohort_id(1:n) == bufs(p)%end_id(1:n))) cycle
                  end if
                  allocate(at(n), key(n), perm(n))
                  call find_ids(bufs(p)%end_id(1:bufs(p)%n_end), r%cohort_id(1:n), at)
                  do j = 1_ik, n
                     key(j) = merge(at(j), bufs(p)%n_end + j, at(j) > 0_ik)
                     perm(j) = j
                  end do
                  call sort_by_key(key, perm)
                  r%cohort_id(1:n) = r%cohort_id(perm)
                  do k = 1_ik, files%reg%nvar
                     if (files%reg%var(k)%dim /= DIM_COHORT) cycle
                     c = slab_col(r, k)
                     if (c <= 0_ik .or. r%nslab(k) <= 0_ik) cycle
                     if (files%reg%var(k)%diel) then
                        call permute_blocks(r%slab(:, c), r%slabvalid(:, c), perm, files%n_diel)
                     else
                        r%slab(1:n, c)      = r%slab(perm, c)
                        r%slabvalid(1:n, c) = r%slabvalid(perm, c)
                     end if
                  end do
                  deallocate(at, key, perm)
               end associate
            end do
         end do
      end do
   end subroutine write_order

   !--- Reorder a slab of blocks of nb entries, one block a cohort: block j takes block perm(j).
   pure subroutine permute_blocks(x, valid, perm, nb)
      real(wp),    intent(inout) :: x(:)
      logical,     intent(inout) :: valid(:)
      integer(ik), intent(in)    :: perm(:), nb
      real(wp)    :: y(size(perm) * nb)
      logical     :: w(size(perm) * nb)
      integer(ik) :: j, n
      n = size(perm, kind=ik)
      do j = 1_ik, n
         y((j - 1_ik)*nb + 1_ik:j*nb) = x((perm(j) - 1_ik)*nb + 1_ik:perm(j)*nb)
         w((j - 1_ik)*nb + 1_ik:j*nb) = valid((perm(j) - 1_ik)*nb + 1_ik:perm(j)*nb)
      end do
      x(1:n*nb)     = y
      valid(1:n*nb) = w
   end subroutine permute_blocks

   !--- Order perm by ascending key (insertion sort: a record is nearly in order already).
   pure subroutine sort_by_key(key, perm)
      integer(ik), intent(inout) :: key(:), perm(:)
      integer(ik) :: i, j, kk, pp
      do i = 2_ik, size(key, kind=ik)
         kk = key(i)
         pp = perm(i)
         j = i - 1_ik
         do while (j >= 1_ik)
            if (key(j) <= kk) exit
            key(j + 1_ik) = key(j)
            perm(j + 1_ik) = perm(j)
            j = j - 1_ik
         end do
         key(j + 1_ik) = kk
         perm(j + 1_ik) = pp
      end do
   end subroutine sort_by_key

end module ed_output_integrate
