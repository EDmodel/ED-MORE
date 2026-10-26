! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_output_integrate -- stage [4] of the diagnostic output, the temporal reduction: fold a    !
! variable's value (ed_output_extract) into its (variable, tier) running reduction             !
! (integrate_scalar, integrate_slab), close a period (normalize_*,                             !
! close_tier) and re-seed (reset_buffer). output_integrate is the per-step tick and            !
! output_integrate_fast the FAST tier's replay. netCDF-free. Each tier integrates raw state    !
! on its own, which for these operators equals chaining (docs/science/diagnostics.md §1, §4).  !
!==========================================================================================!
module ed_output_integrate
   use ed_kinds,            only : wp, ik
   use ed_time,             only : ed_time_t, calendar_step_t
   use ed_output_params,      only : N_FREQ,   FREQ_FAST, TIER_FAST, TIER_DAILY, TIER_MONTHLY, TIER_ANNUAL
   use ed_output_types,   only :   var_desc_t, integ_buffer_t, output_files_t, output_buffers_t,                  &
                                   pending_record_t, record_queue_t, AGG_MEAN, AGG_RATE, AGG_LAST, AGG_TMEAN,     &
                                   DIM_SCALAR, DIM_COHORT, DIM_PATCH,                                             &
                                   DIM_SOIL, DIM_SOIL_PATCH, MISSING_VALUE, slab_col
   use ed_polygon_state_vars,   only : polygon_t
   use ed_constants,     only   : n_soil_layer_max
   use ed_output_extract, only :   SRC_F_CO0, SRC_F_COH_GPP, SRC_F_COH_HEIGHT, SRC_F_PA0,                         &
                                   SRC_F_SOIL_WATER, SRC_NONE, FLD_C_DIAG0, extract_fast_scalar, extract_variable
   implicit none
   private

   public :: alloc_integ_buffer, reset_buffer, integrate_scalar, integrate_slab
   public :: normalize_scalar, normalize_slab
   public :: output_tick, output_integrate, output_integrate_fast, close_tier

contains

   !=======================================================================================!
   !  A polygon's output tick after a slow step: replay the FAST samples its fast loop      !
   !  staged, closing that tier every fast_interval_steps of them over the run, then fold   !
   !  the step into the slower tiers (output_integrate), whose snapshots of fast-loop values !
   !  read the last staged sub-step. Closed records are queued for the I/O phase.          !
   !=======================================================================================!
   subroutine output_tick(files, bufs, poly, step, dt_slow, dt_fast)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(in)    :: poly
      type(calendar_step_t),  intent(in)    :: step
      real(wp),               intent(in)    :: dt_slow, dt_fast   !< [s] the slow step, a fast sub-step
      integer(ik) :: isub
      if (bufs%fast_ready .and. files%reg%nidx(TIER_FAST) > 0_ik) then
         do isub = 1_ik, bufs%n_fast_sub
            call output_integrate_fast(files, bufs, isub, dt_fast)
            bufs%fast_step_total = bufs%fast_step_total + 1_ik
            if (mod(bufs%fast_step_total, max(files%fast_interval_steps, 1_ik)) == 0_ik)          &
               call close_tier(files, bufs, TIER_FAST)
         end do
      end if
      call output_integrate(files, bufs, poly, step%start, dt_slow, step%new_day, step%new_month,     &
                            step%new_year)
      bufs%fast_ready = .false.
   end subroutine output_tick

   !=======================================================================================!
   !  Buffer lifecycle                                                                      !
   !=======================================================================================!

   !----- Allocate + seed one (variable, tier) buffer. slab arrays sized to `cap` for a       !
   !      non-scalar dim; scalar dims allocate none.                                            !
   subroutine alloc_integ_buffer(buf, v, freq, cap)
      type(integ_buffer_t), intent(out) :: buf
      type(var_desc_t),     intent(in)  :: v
      integer(ik),          intent(in)  :: freq, cap
      buf%var_id = 0_ik      ! set by the caller (registry index)
      buf%freq   = freq
      buf%agg    = v%agg
      !----- The F tier folds staged sub-step samples, rates and states alike: a mean variable their  !
      !      dt-weighted mean, a snapshot the last. ------------------------------------------------!
      if (freq == FREQ_FAST .and. v%agg /= AGG_LAST) buf%agg = AGG_TMEAN
      buf%dim    = v%dim
      buf%active = .true.
      if (v%dim /= DIM_SCALAR) then
         allocate(buf%slab(cap), buf%wsum_slab(cap), buf%hits(cap))
      end if
      call reset_buffer(buf)
   end subroutine alloc_integ_buffer

   !----- Re-seed a buffer for the next period. ---------------------------------------------!
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


   !=======================================================================================!
   !  Integrate (fold one step)                                                             !
   !=======================================================================================!

   !----- Fold one scalar sample x over a step of length dt [s] into a scalar buffer. ------!
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

   !----- Fold a slab of n samples by direct slot index. `valid` (optional) marks the slots with  !
   !      something to report this step; the rest (an empty patch, PFT or class, a layer below the !
   !      column) are skipped and close as _FillValue, not as a 0 a reader would take for a        !
   !      measurement (docs/science/diagnostics.md §2). Without `valid` every slot is folded. -----!
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

   !=======================================================================================!
   !  Normalize (close a period). valid=.false. -> emit MISSING_VALUE (no NaN, no ±huge).     !
   !=======================================================================================!
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

   !----- Slab normalize: one value + validity per live slot [1:n_out]. --------------------!
   subroutine normalize_slab(buf, out, valid, n_out)
      type(integ_buffer_t), intent(in)  :: buf
      real(wp),             intent(out) :: out(:)
      logical,              intent(out) :: valid(:)
      integer(ik),          intent(out) :: n_out
      integer(ik) :: i
      n_out = buf%n_slab
      !----- A record row per live slot. A shorter record means the slab was sized before the      !
      !      registry was final, and the loop below would write past this column into the next. --!
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

   !=======================================================================================!
   !  The per-step tick. The slow step just taken covers [prev, now): it is folded into the  !
   !  window that holds prev, and only then is each window that now has left closed, so a   !
   !  step belongs to the period it starts in (docs/science/diagnostics.md §4, Which steps  !
   !  a record holds). Closed records go to bufs%queue. The calendar restructuring runs     !
   !  after this tick, at the next step's start (ed_stepper).                               !
   !=======================================================================================!
   subroutine output_integrate(files, bufs, poly, prev, dt, is_new_day, is_new_month, is_new_year)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(in)    :: poly
      type(ed_time_t),        intent(in)    :: prev          !< the start of the step just taken
      real(wp),               intent(in)    :: dt
      logical,                intent(in)    :: is_new_day, is_new_month, is_new_year  !< its end left prev's period
      integer(ik) :: t, j, k, n_out
      real(wp)    :: scal
      real(wp)    :: slab(max(files%max_slab, 1_ik))
      logical     :: vslab(max(files%max_slab, 1_ik))
      if (.not. files%enabled) return
      !----- The cohort and patch slabs are output.cohort_max and patch_max long. Stop at the step the  !
      !      stand outgrows one, not when its record is written, up to a month later. ------------------!
      if (files%cohort_axis .and. poly%cohort%n > files%cohort_max) call cap_exceeded('cohort', poly%cohort%n, files%cohort_max)
      if (files%patch_axis  .and. poly%patch%n  > files%patch_max)  call cap_exceeded('patch',  poly%patch%n,  files%patch_max)
      !----- Fold the step into the daily, monthly and annual tiers. The FAST tier is fed from the !
      !      staged sub-steps (output_integrate_fast). -------------------------------------------!
      do t = TIER_DAILY, N_FREQ
         if (files%reg%nidx(t) == 0_ik) cycle
         if (.not. bufs%has_data(t)) bufs%t_open(t) = prev
         do j = 1_ik, files%reg%nidx(t)
            k = files%reg%idx_freq(j, t)
            if (files%reg%var(k)%snap_staged) then
               call staged_snapshot(bufs, poly, files%reg%var(k), scal, slab, vslab, n_out)
            else
               call extract_variable(poly, bufs%diag, files%reg%var(k), scal, slab, vslab, n_out)
            end if
            if (files%reg%var(k)%dim == DIM_SCALAR) then
               call integrate_scalar(bufs%buf(k,t), scal, dt)
            else
               call integrate_slab(bufs%buf(k,t), slab, n_out, dt, vslab)
            end if
         end do
         bufs%has_data(t) = .true.
      end do
      !----- then close the periods the step's end has left (tiers independent). ----------------!
      if (is_new_year  .and. bufs%has_data(TIER_ANNUAL))  call close_tier(files, bufs, TIER_ANNUAL)
      if (is_new_month .and. bufs%has_data(TIER_MONTHLY)) call close_tier(files, bufs, TIER_MONTHLY)
      if (is_new_day   .and. bufs%has_data(TIER_DAILY))   call close_tier(files, bufs, TIER_DAILY)
   end subroutine output_integrate

   !----- A live count past its output cap: say which, by how much, and what to raise. -----------!
   subroutine cap_exceeded(axis, n, cap)
      character(len=*), intent(in) :: axis
      integer(ik),      intent(in) :: n, cap
      write(*,'(5a,i0,a,i0,a)') ' output: the stand has more ', axis, 's than output.', axis, '_max (', n, ' > ', &
                                cap, '); raise it in [output]'
      error stop 'output: a live count exceeds its output cap (see the message above)'
   end subroutine cap_exceeded

   !----- Normalize a tier's buffers into its pending record + reset them (staging). ----------!
   subroutine close_tier(files, bufs, t)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: t
      integer(ik) :: j, k, c, ns, nsl
      integer(ik) :: slab_k(files%reg%nidx(t))
      logical     :: soil_patch(files%reg%nidx(t))
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
         if (files%reg%var(k)%dim == DIM_SCALAR) then
            call normalize_scalar(bufs%buf(k,t), bufs%pending(t)%sval(k), bufs%pending(t)%svalid(k))
         else
            c = slab_col(bufs%pending(t), k)
            call normalize_slab(bufs%buf(k,t), bufs%pending(t)%slab(:,c),                          &
                                bufs%pending(t)%slabvalid(:,c), ns)
            bufs%pending(t)%nslab(k) = ns
            !----- The record's axis lengths. A soil-by-patch slab is patch-major, strided by the  !
            !      layer ceiling, so it holds ns / n_soil_layer_max patches. ---------------------!
            select case (files%reg%var(k)%dim)
            case (DIM_COHORT)
               bufs%pending(t)%n_cohort = max(bufs%pending(t)%n_cohort, ns)
            case (DIM_PATCH)
               bufs%pending(t)%n_patch  = max(bufs%pending(t)%n_patch,  ns)
            case (DIM_SOIL_PATCH)
               bufs%pending(t)%n_patch  = max(bufs%pending(t)%n_patch,  ns / n_soil_layer_max)
            end select
         end if
         call reset_buffer(bufs%buf(k,t))
      end do
      bufs%has_data(t) = .false.
      nsl = 0_ik                                          ! this tier's slab variables, for the queue
      do j = 1_ik, files%reg%nidx(t)
         k = files%reg%idx_freq(j, t)
         if (files%reg%var(k)%dim == DIM_SCALAR) cycle
         nsl = nsl + 1_ik
         slab_k(nsl) = k
         soil_patch(nsl) = files%reg%var(k)%dim == DIM_SOIL_PATCH
      end do
      call enqueue_record(bufs%queue(t), bufs%pending(t), slab_k(1:nsl), soil_patch(1:nsl))
      bufs%pending(t)%used = .false.
   end subroutine close_tier

   !----- Append a copy of the scratch record to the tier's queue. Only what the writer reads is   !
   !      copied: the tier's slab variables (`slab_k`), each to its nslab rows, or the whole patch x !
   !      layer block (n_patch * n_soil_layer_max) for a soil-patch variable. The copy is exact, so  !
   !      a queued record writes exactly what an immediate write from the scratch would have. -----!
   subroutine enqueue_record(q, src, slab_k, soil_patch)
      type(record_queue_t),   intent(inout) :: q
      type(pending_record_t), intent(in)    :: src
      integer(ik),            intent(in)    :: slab_k(:)
      logical,                intent(in)    :: soil_patch(:)
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
         if (soil_patch(c)) rows = max(rows, src%n_patch * n_soil_layer_max)
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

   !=======================================================================================!
   !  Fold one staged sub-step (isub) into the FAST tier's buffers. output_tick calls it   !
   !  once per staged sub-step and closes the tier every fast_interval_steps folds. It     !
   !  reads only the buffers' staging, never the polygon's live state (the end of the slow !
   !  step by now).                                                                        !
   !=======================================================================================!
   subroutine output_integrate_fast(files, bufs, isub, dt)
      type(output_files_t),   intent(in)    :: files
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: isub
      real(wp),               intent(in)    :: dt
      integer(ik) :: j, k, n
      real(wp)    :: scal
      real(wp)    :: slab(max(files%max_slab, 1_ik))
      logical     :: valid(max(files%max_slab, 1_ik))
      if (.not. files%enabled) return
      if (files%reg%nidx(TIER_FAST) == 0_ik) return
      !----- The staging follows the stand; the slabs are output.cohort_max and patch_max long. The  !
      !      same check as the coarse tiers' (output_integrate), before this tier writes into them. --!
      if (files%cohort_axis .and. bufs%fast_n_cohort > files%cohort_max)                           &
         call cap_exceeded('cohort', bufs%fast_n_cohort, files%cohort_max)
      if (files%patch_axis .and. bufs%fast_n_patch > files%patch_max)                              &
         call cap_exceeded('patch', bufs%fast_n_patch, files%patch_max)
      if (.not. bufs%has_data(TIER_FAST)) bufs%t_open(TIER_FAST) = bufs%fast_time(isub)
      do j = 1_ik, files%reg%nidx(TIER_FAST)
         k = files%reg%idx_freq(j, TIER_FAST)
         call staged_value(bufs, files%reg%var(k)%source_fast, files%reg%var(k)%dim, isub, scal, slab, valid, n)
         if (files%reg%var(k)%dim == DIM_SCALAR) then
            call integrate_scalar(bufs%buf(k,TIER_FAST), scal, dt)
         else
            call integrate_slab(bufs%buf(k,TIER_FAST), slab, n, dt, valid)
         end if
      end do
      bufs%has_data(TIER_FAST) = .true.
   end subroutine output_integrate_fast

   !----- One staged source's value in sub-step isub, on the variable's axis: the polygon mean of  !
   !      the patch rows or a forcing row, each patch's row, a soil column, or the per-cohort      !
   !      values. Every slot is valid but a layer below its column. A per-cohort diagnostic row    !
   !      is staged only while the block is on (a fast loop that ran); else it has no slots. ------!
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

   !----- A snapshot of a fast-loop diagnostic: its value over the last fast step of the slow step  !
   !      just taken, the last staged sub-step (docs/science/diagnostics.md §3). A per-cohort row   !
   !      reduced to a coarser axis is that sub-step's row, reduced as its mean is. With nothing    !
   !      staged this step (no fast loop) it is missing, and the period keeps its last value. -----!
   subroutine staged_snapshot(bufs, poly, v, scal, slab, valid, n)
      type(output_buffers_t), intent(in)  :: bufs
      type(polygon_t),        intent(in)  :: poly
      type(var_desc_t),       intent(in)  :: v
      real(wp),               intent(out) :: scal, slab(:)
      logical,                intent(out) :: valid(:)
      integer(ik),            intent(out) :: n
      integer(ik) :: isub, nc
      scal = MISSING_VALUE
      n = 0_ik
      if (.not. bufs%fast_ready) return
      isub = bufs%n_fast_sub
      if (v%source_fast /= SRC_NONE) then
         call staged_value(bufs, v%source_fast, v%dim, isub, scal, slab, valid, n)
      else if (allocated(bufs%fast_cohort_diag)) then
         nc = bufs%fast_n_cohort
         if (nc == poly%cohort%n)                                                                    &
            call extract_variable(poly, bufs%diag, v, scal, slab, valid, n,                          &
                                  cohort_x=bufs%fast_cohort_diag(1:nc, v%source_id - FLD_C_DIAG0, isub))
      end if
   end subroutine staged_snapshot

end module ed_output_integrate
