! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_output_reduce -- stage [3] of the diagnostic output (docs/science/diagnostics.md
! §1-§2): one weighted reduction of a per-cohort field to each patch, the polygon, each PFT
! or each DBH class, and of a per-patch field to the polygon or to each site's soil column.
! The weight kind and the sum or mean come from the variable's descriptor. An empty set is
! not valid (the caller writes _FillValue), except a sum over a PFT or class, a true 0.
! Also the polygon totals the driver and the C API read. Keep each sum(array-expression) as
! written: a loop may associate differently and move every reported total in its last bit.
!---~---
module ed_output_reduce
   use ed_kinds,       only : wp, ik
   use ed_constants,   only : tiny_num, cm2_to_m2, n_soil_layer_max
   use ed_polygon_state_vars,   only : polygon_t, patch_polygon_area, patch_polygon_areas
   use ed_output_types,         only : MISSING_VALUE, W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA, W_AGB
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   implicit none
   private

   !--- The generic API (stage 3).
   public :: W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA, W_AGB
   public :: reduce_cohort_to_polygon, reduce_cohort_to_patch
   public :: reduce_cohort_to_pft, reduce_cohort_to_size
   public :: reduce_patch_to_polygon, reduce_patch_columns_to_sites, gather_patch_columns
   !---~---
   !   Polygon totals for the callers that want one number (the driver's summaries, the C
   ! API, the tests).
   !---~---
   public :: total_nplant, total_basal_area, total_agb, total_lai, total_area
   public :: print_summary
   public :: dbh_class_index

contains

   !--- Weights

   !---~---
   !   Fill w(1:n) with the per-cohort weight of `wkind` over the polygon's whole cohort
   ! array: the live stand's, or `w_given`'s column of that kind (the weights the fast loop
   ! used, for a staged sub-step's row). A cohort whose x is missing (one recruited after the
   ! fast loop has no samples yet) carries no weight, so it is left out of a sum and of both
   ! halves of a mean.
   !---~---
   pure subroutine cohort_weight(poly, wkind, x, w, n, w_given)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in)  :: wkind
      real(wp),     intent(in)  :: x(:)
      real(wp),     intent(out) :: w(:)
      integer(ik),  intent(out) :: n
      real(wp), optional, intent(in) :: w_given(:,:)   ! (cohort, W_NPLANT:W_AGB)
      n = poly%cohort%n
      if (n <= 0_ik) return
      if (present(w_given) .and. wkind /= W_NONE) then
         w(1:n) = w_given(1:n, wkind)
         where (x(1:n) == MISSING_VALUE) w(1:n) = 0.0_wp
         return
      end if
      select case (wkind)
      case (W_NPLANT)
         w(1:n) = poly%cohort%nplant(1:n)
      case (W_LEAF_AREA)
         w(1:n) = poly%cohort%nplant(1:n) * poly%cohort%leaf_area(1:n)
      case (W_BASAL_AREA)
         w(1:n) = poly%cohort%nplant(1:n) * poly%cohort%basal_area(1:n)
      case (W_AGB)
         w(1:n) = poly%cohort%nplant(1:n) * poly%cohort%agb(1:n)
      case default
         w(1:n) = 1.0_wp                       ! W_NONE
      end select
      where (x(1:n) == MISSING_VALUE) w(1:n) = 0.0_wp
   end subroutine cohort_weight

   !--- cohort -> polygon

   !---~---
   !   Reduce a per-cohort field x(1:ncohort) to one polygon scalar, each patch weighted by
   ! its share of the polygon. `scale` multiplies each patch's contribution before it is
   ! summed; applying it once at the end would change the rounding of every total.
   !---~---
   pure subroutine reduce_cohort_to_polygon(poly, x, wkind, mean, val, valid, scale, w_given)
      type(polygon_t),    intent(in)  :: poly
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind
      logical,            intent(in)  :: mean
      real(wp),           intent(out) :: val
      logical,            intent(out) :: valid
      real(wp), optional, intent(in)  :: scale
      real(wp), optional, intent(in)  :: w_given(:,:)   ! the fast loop's weights (cohort_weight)
      real(wp)    :: w(max(poly%cohort%n, 1_ik)), share(max(poly%patch%n, 1_ik)), num, den, sc
      integer(ik) :: ip, i0, i1, nw
      sc = 1.0_wp
      if (present(scale)) sc = scale
      call cohort_weight(poly, wkind, x, w, nw, w_given)
      call patch_polygon_areas(poly, share)
      num = 0.0_wp
      den = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         if (i1 < i0) cycle
         num = num + share(ip) * sum(w(i0:i1) * x(i0:i1)) * sc
         if (mean) den = den + share(ip) * sum(w(i0:i1))
      end do
      if (mean) then
         valid = (abs(den) > tiny_num)
         val   = 0.0_wp
         if (valid) val = num / den
      else
         valid = .true.
         val = num
      end if
   end subroutine reduce_cohort_to_polygon

   !--- cohort -> patch

   !---~---
   !   Reduce a per-cohort field to one value per patch, with no area factor: a patch value
   ! is per unit of its own ground. An empty patch is not valid (_FillValue).
   !---~---
   pure subroutine reduce_cohort_to_patch(poly, x, wkind, mean, out, valid, n_out, scale, w_given)
      type(polygon_t),    intent(in)  :: poly
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind
      logical,            intent(in)  :: mean
      real(wp),           intent(out) :: out(:)
      logical,            intent(out) :: valid(:)
      integer(ik),        intent(out) :: n_out
      real(wp), optional, intent(in)  :: scale
      real(wp), optional, intent(in)  :: w_given(:,:)   ! the fast loop's weights (cohort_weight)
      real(wp)    :: w(max(poly%cohort%n, 1_ik)), num, den, sc
      integer(ik) :: ip, i0, i1, nw
      sc = 1.0_wp
      if (present(scale)) sc = scale
      call cohort_weight(poly, wkind, x, w, nw, w_given)
      n_out = poly%patch%n
      out(1:max(n_out,1_ik)) = 0.0_wp
      valid(1:max(n_out,1_ik)) = .false.
      do ip = 1_ik, n_out
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         if (i1 < i0) cycle                                   ! empty patch -> _FillValue
         num = sum(w(i0:i1) * x(i0:i1)) * sc
         if (mean) then
            den = sum(w(i0:i1))
            if (abs(den) > tiny_num) then
               out(ip) = num / den
               valid(ip) = .true.
            end if
         else
            out(ip) = num
            valid(ip) = .true.
         end if
      end do
   end subroutine reduce_cohort_to_patch

   !--- cohort -> PFT   (DIM_PFT)

   !---~---
   !   Reduce a per-cohort field onto the PFT axis, each patch weighted by its share of the
   ! polygon, so a value is per m2 of polygon ground and the PFTs sum to the polygon total.
   !---~---
   pure subroutine reduce_cohort_to_pft(poly, x, wkind, mean, n_pft, out, valid, n_out, scale, w_given)
      type(polygon_t),    intent(in)  :: poly
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind, n_pft
      logical,            intent(in)  :: mean
      real(wp),           intent(out) :: out(:)
      logical,            intent(out) :: valid(:)
      integer(ik),        intent(out) :: n_out
      real(wp), optional, intent(in)  :: scale
      real(wp), optional, intent(in)  :: w_given(:,:)   ! the fast loop's weights (cohort_weight)
      real(wp)    :: w(max(poly%cohort%n, 1_ik)), sc, aw
      real(wp)    :: num(max(n_pft,1_ik)), den(max(n_pft,1_ik))
      integer(ik) :: ip, i0, i1, i, k, nw
      sc = 1.0_wp
      if (present(scale)) sc = scale
      call cohort_weight(poly, wkind, x, w, nw, w_given)
      n_out = n_pft
      num = 0.0_wp
      den = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         aw = patch_polygon_area(poly, ip)
         do i = i0, i1
            k = poly%cohort%pft(i)
            if (k < 1_ik .or. k > n_pft) cycle
            num(k) = num(k) + aw * w(i) * x(i) * sc
            den(k) = den(k) + aw * w(i)
         end do
      end do
      do k = 1_ik, n_pft
         if (mean) then
            valid(k) = (abs(den(k)) > tiny_num)
            out(k)   = 0.0_wp
            if (valid(k)) out(k) = num(k) / den(k)
         else
            !---~---
            !   A sum over a PFT with no members is a true 0, so the closure of
            ! docs/science/diagnostics.md eq. (3) holds when a PFT goes locally extinct.
            !---~---
            valid(k) = .true.
            out(k) = num(k)
         end if
      end do
   end subroutine reduce_cohort_to_pft

   !--- cohort -> DBH size class   (DIM_SIZE)

   !---~---
   !   Reduce a per-cohort field onto the DBH-class axis: each cohort whole in the class of
   ! its mean dbh, weighted and closed over empty classes as on the PFT axis above.
   !---~---
   pure subroutine reduce_cohort_to_size(poly, x, wkind, mean, edges, n_class, out, valid,     &
                                         n_out, scale, w_given)
      type(polygon_t),    intent(in)  :: poly
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind, n_class
      logical,            intent(in)  :: mean
      real(wp),           intent(in)  :: edges(:)
      real(wp),           intent(out) :: out(:)
      logical,            intent(out) :: valid(:)
      integer(ik),        intent(out) :: n_out
      real(wp), optional, intent(in)  :: scale
      real(wp), optional, intent(in)  :: w_given(:,:)   ! the fast loop's weights (cohort_weight)
      real(wp)    :: w(max(poly%cohort%n, 1_ik)), sc, aw
      real(wp)    :: num(max(n_class,1_ik)), den(max(n_class,1_ik))
      integer(ik) :: ip, i0, i1, i, k, nw
      sc = 1.0_wp
      if (present(scale)) sc = scale
      call cohort_weight(poly, wkind, x, w, nw, w_given)
      n_out = n_class
      num = 0.0_wp
      den = 0.0_wp
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         aw = patch_polygon_area(poly, ip)
         do i = i0, i1
            k = dbh_class_index(poly%cohort%dbh(i), edges, n_class)
            num(k) = num(k) + aw * w(i) * x(i) * sc
            den(k) = den(k) + aw * w(i)
         end do
      end do
      do k = 1_ik, n_class
         if (mean) then
            valid(k) = (abs(den(k)) > tiny_num)
            out(k)   = 0.0_wp
            if (valid(k)) out(k) = num(k) / den(k)
         else
            valid(k) = .true.
            out(k) = num(k)
         end if
      end do
   end subroutine reduce_cohort_to_size

   !--- patch -> polygon

   !---~---
   !   Reduce a per-patch field to one polygon scalar, each patch weighted by its share of
   ! the polygon. mean = .false. gives the plain weighted sum, already the mean of a
   ! per-ground-area quantity since the shares sum to 1 (the driver checks it at the end of
   ! a run); mean = .true. divides by their sum too, which holds even where they do not.
   !---~---
   pure subroutine reduce_patch_to_polygon(poly, x, mean, val, valid)
      type(polygon_t), intent(in) :: poly
      real(wp),     intent(in)  :: x(:)
      logical,      intent(in)  :: mean
      real(wp),     intent(out) :: val
      logical,      intent(out) :: valid
      real(wp)    :: share(max(poly%patch%n, 1_ik)), num, den
      integer(ik) :: ip
      call patch_polygon_areas(poly, share)
      num = 0.0_wp
      den = 0.0_wp
      do ip = 1_ik, poly%patch%n
         num = num + share(ip) * x(ip)
         den = den + share(ip)
      end do
      if (mean) then
         valid = (den > tiny_num)
         val   = 0.0_wp
         if (valid) val = num / den
      else
         valid = (poly%patch%n > 0_ik)
         val = num
      end if
   end subroutine reduce_patch_to_polygon

   !---~---
   !   Area-weight each site's per-patch soil columns into the site's column, flattened at
   ! (is-1)*n_soil_layer_max + k (DIM_SOIL). A site's layers are its own column's; the slots
   ! below its bottom layer are not valid.
   !---~---
   pure subroutine reduce_patch_columns_to_sites(poly, col, out, valid, n_out)
      type(polygon_t), intent(in) :: poly
      real(wp),     intent(in)  :: col(:,:)   ! (layer, patch)
      real(wp),     intent(out) :: out(:)
      logical,      intent(out) :: valid(:)
      integer(ik),  intent(out) :: n_out
      real(wp)    :: share(max(poly%patch%n, 1_ik)), wsum, s(n_soil_layer_max)
      integer(ik) :: is, ip, nl, base
      call patch_polygon_areas(poly, share)
      n_out = poly%site%n * n_soil_layer_max
      out(1:n_out)   = 0.0_wp
      valid(1:n_out) = .false.
      do is = 1_ik, poly%site%n
         nl   = poly%site%surface(is)%soil_profile%n_active
         base = (is - 1_ik) * n_soil_layer_max
         s(1:nl) = 0.0_wp
         wsum = 0.0_wp
         do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
            s(1:nl) = s(1:nl) + share(ip) * col(1:nl, ip)
            wsum = wsum + share(ip)
         end do
         if (wsum > tiny_num) s(1:nl) = s(1:nl) / wsum
         out(base + 1_ik:base + nl)   = s(1:nl)
         valid(base + 1_ik:base + nl) = .true.
      end do
   end subroutine reduce_patch_columns_to_sites

   !---~---
   !   Flatten per-patch columns into the (layer, patch) slab at (ip-1)*n_soil_layer_max + k,
   ! with no reduction (DIM_SOIL_PATCH). Each patch's layers are its site's; the stride is
   ! the layer ceiling, never a live count (docs/science/diagnostics.md §3).
   !---~---
   pure subroutine gather_patch_columns(poly, col, out, valid, n_out)
      type(polygon_t), intent(in) :: poly
      real(wp),     intent(in)  :: col(:,:)   ! (layer, patch)
      real(wp),     intent(out) :: out(:)
      logical,      intent(out) :: valid(:)
      integer(ik),  intent(out) :: n_out
      integer(ik) :: ip, k, base
      n_out = poly%patch%n * n_soil_layer_max
      out(1:max(n_out,1_ik))   = 0.0_wp
      valid(1:max(n_out,1_ik)) = .false.
      do ip = 1_ik, poly%patch%n
         base = (ip - 1_ik) * n_soil_layer_max
         do k = 1_ik, poly%site%surface(poly%patch%owner_site(ip))%soil_profile%n_active
            out(base + k)   = col(k, ip)
            valid(base + k) = .true.
         end do
      end do
   end subroutine gather_patch_columns

   !--- Polygon totals (thin wrappers over the reduction).

   !--- Polygon plant number [plant m-2 ground].
   pure real(wp) function total_nplant(poly) result(tot)
      type(polygon_t), intent(in) :: poly
      logical :: ok
      call reduce_cohort_to_polygon(poly, poly%cohort%nplant, W_NONE, .false., tot, ok)
   end function total_nplant

   !--- Polygon basal area [m2 m-2]: per-plant cm2 -> per-ground m2.
   pure real(wp) function total_basal_area(poly) result(tot)
      type(polygon_t), intent(in) :: poly
      logical :: ok
      call reduce_cohort_to_polygon(poly, poly%cohort%basal_area, W_NPLANT, .false., tot, ok,    &
                                 scale=cm2_to_m2)
   end function total_basal_area

   !--- Polygon aboveground biomass [kgC m-2 ground].
   pure real(wp) function total_agb(poly) result(tot)
      type(polygon_t), intent(in) :: poly
      logical :: ok
      call reduce_cohort_to_polygon(poly, poly%cohort%agb, W_NPLANT, .false., tot, ok)
   end function total_agb

   !--- Polygon leaf area index [m2 m-2 ground].
   pure real(wp) function total_lai(poly) result(tot)
      type(polygon_t), intent(in) :: poly
      logical :: ok
      call reduce_cohort_to_polygon(poly, poly%cohort%leaf_area, W_NPLANT, .false., tot, ok)
   end function total_lai

   !--- The patches' shares of the polygon, summed: 1 while no patch area is lost.
   pure real(wp) function total_area(poly) result(tot)
      type(polygon_t), intent(in) :: poly
      real(wp) :: share(max(poly%patch%n, 1_ik))
      call patch_polygon_areas(poly, share)
      tot = sum(share(1:poly%patch%n))
   end function total_area

   !--- Basal-area-weighted mean diameter [cm] (a simple structural index).
   pure real(wp) function mean_dbh(poly) result(dm)
      type(polygon_t), intent(in) :: poly
      logical :: ok
      call reduce_cohort_to_polygon(poly, poly%cohort%dbh, W_BASAL_AREA, .true., dm, ok)
   end function mean_dbh


   !--- Health checks / reporting


   !--- Human-readable one-line summary.
   subroutine print_summary(poly, label)
      type(polygon_t), intent(in) :: poly
      character(len=*), intent(in) :: label
      write(*,'(a,t18,a,i6,a,i4,a,f9.4,a,f8.4,a,f8.3,a,f7.2)')                             &
         trim(label), 'cohorts=', poly%cohort%n, '  patches=', poly%patch%n,                    &
         '  N=', total_nplant(poly), '  LAI=', total_lai(poly),                            &
         '  AGB=', total_agb(poly), '  Dmean=', mean_dbh(poly)
   end subroutine print_summary

   !---~---
   !   1-based DBH class of a dbh (docs/science/diagnostics.md §3). Over n_class+1 ascending
   ! edges, bin k is [edges(k), edges(k+1)); the last bin also takes every dbh above its
   ! top edge, and a dbh below edges(1) falls in bin 1, so no plant is dropped and
   ! sum(mean_nplant_size_py) == mean_nplant_py holds.
   !---~---
   pure integer(ik) function dbh_class_index(dbh, edges, n_class) result(k)
      real(wp),    intent(in) :: dbh
      real(wp),    intent(in) :: edges(:)   ! (n_class+1) ascending
      integer(ik), intent(in) :: n_class
      integer(ik) :: j
      k = 1_ik
      if (n_class <= 0_ik) return
      do j = n_class, 1_ik, -1_ik
         if (dbh >= edges(j)) then
            k = j
            return
         end if
      end do
   end function dbh_class_index

end module ed_output_reduce
