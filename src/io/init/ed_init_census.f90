! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_init_census -- init_mode 1, a start from a cohort census CSV (docs/configuration.md,   !
! How a run starts). read_census_file reads the file into a census_table_t, one row per    !
! cohort; init_from_census builds the polygon's patches and cohorts from it, on the        !
! bare-ground builders; restructure_census_stand then fuses and splits the stand with the  !
! slow step's own operators before the first step.                                         !
!==========================================================================================!
module ed_init_census
   use ed_kinds,        only : wp, ik
   use ed_config,       only : ed_config_t
   use ed_polygon_state_vars,   only : polygon_t
   use ed_init_bare_ground,     only : init_bare_ground, init_bare_sites, add_cohort, finalize_init
   use ed_slow_state,    only : kill_report_t, kill_litter
   use ed_slow_update,   only : sort_cohorts, fuse_cohorts, terminate_cohorts, split_cohorts,        &
                                sort_patches, fuse_patches, terminate_patches
   implicit none
   private

   public :: census_table_t, read_census_file, init_from_census, restructure_census_stand

   !----- The census file's limits: columns in the header, characters in a field. ---------!
   integer, parameter :: CENSUS_MAX_COLS = 32, CENSUS_FIELD_LEN = 64
   !----- Two rows of one patch agree on its area within this fraction of it. -------------!
   real(wp), parameter :: PATCH_AREA_RTOL = 1.0e-9_wp
   !----- The census file's columns: each one's position in a row, 0 if absent. -----------!
   type :: census_columns_t
      logical :: header = .false.
      integer :: n = 0
      integer :: site = 0, patch = 0, cohort = 0, dbh = 0, height = 0, pft = 0, nplant = 0, area = 0
   end type census_columns_t

   !----- A census file's rows, in file order: one cohort each. ------------------------------!
   type :: census_table_t
      integer(ik) :: n = 0_ik                           !< rows
      logical     :: has_area = .false.                 !< the file has a patch_area column
      integer(ik), allocatable :: site_id(:)            !< 1 where the file has no site_id column
      integer(ik), allocatable :: patch_id(:)
      integer(ik), allocatable :: pft(:)
      real(wp),    allocatable :: dbh(:)                !< [cm]
      real(wp),    allocatable :: nplant(:)             !< [plants per m2 of the patch]
      real(wp),    allocatable :: patch_area(:)         !< [any unit]; 1 without the column
   end type census_table_t

contains

   !---------------------------------------------------------------------------------------!
   ! Read a census CSV into `table` (docs/configuration.md, How a run starts). Columns are  !
   ! matched by header name: patch_id, dbh [cm], pft and nplant required; patch_area,      !
   ! site_id, cohort_id and height [m] optional. A first data line of seven numbers is     !
   ! read positionally as site_id, patch_id, cohort_id, dbh, height, pft, nplant. Blank    !
   ! and '#' lines are skipped. why says what is wrong with the file, '' when it reads: it  !
   ! cannot be opened, or has an unknown, repeated or missing column, a malformed row or a !
   ! non-positive value.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine read_census_file(path, table, why)
      character(len=*),              intent(in)  :: path
      type(census_table_t),          intent(out) :: table
      character(len=:), allocatable, intent(out) :: why
      integer,     parameter :: LINELEN = 1024
      character(len=LINELEN) :: line
      type(census_columns_t) :: cols
      integer(ik)            :: sid, pid, ipft
      real(wp)               :: dbh, height, nplant, area
      logical                :: have_layout
      integer                :: u, ios

      allocate(table%site_id(0), table%patch_id(0), table%pft(0), table%dbh(0), table%nplant(0),      &
               table%patch_area(0))
      why = 'cannot be opened'
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      why = ''
      have_layout = .false.
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (.not. is_data_line(line)) cycle
         if (.not. have_layout) then
            call census_layout(line, cols, why)
            if (len(why) > 0) exit
            have_layout = .true.
            if (cols%header) cycle
         end if
         call parse_census_line(line, cols, sid, pid, dbh, height, ipft, nplant, area, why)
         if (len(why) > 0) exit
         table%site_id    = [table%site_id, sid]
         table%patch_id   = [table%patch_id, pid]
         table%pft        = [table%pft, ipft]
         table%dbh        = [table%dbh, dbh]
         table%nplant     = [table%nplant, nplant]
         table%patch_area = [table%patch_area, area]
      end do
      close(u)
      table%n = int(size(table%site_id), ik)
      table%has_area = cols%area > 0
   end subroutine read_census_file

   !---------------------------------------------------------------------------------------!
   ! Build a polygon from a cohort census CSV (read_census_file). Each (site_id, patch_id)  !
   ! is a primary patch of age 0, first seen first within its site; with patch_area a      !
   ! site's patches take the file's areas normalized to the site, else equal shares. One   !
   ! site is built from use_site_id (else the first site_id seen); with site_ids, site k   !
   ! from the rows of site_ids(k), with site_cfg(k) and site_area(k). found = .false. and  !
   ! `why` say what is wrong with the file (read_census_file), that it has no row for the  !
   ! site, a patch whose rows disagree on patch_area, or a pft out of range.               !
   !---------------------------------------------------------------------------------------!
   subroutine init_from_census(poly, cfg, path, found, use_site_id, why, site_ids, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(ed_config_t),   intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      logical,             intent(out) :: found
      integer(ik), intent(in), optional :: use_site_id
      character(len=:), allocatable, intent(out), optional :: why   !< why found = .false.; '' when found
      integer(ik),         intent(in), optional :: site_ids(:)      !< several sites: each one's census site_id,
      type(ed_config_t),   intent(in), optional :: site_cfg(:)      !< configuration
      real(wp),            intent(in), optional :: site_area(:)     !< and share of the polygon
      type(census_table_t)     :: table
      integer(ik)              :: npatch, ip, k, i, j, r
      integer(ik), allocatable :: sel(:), site_of(:), patch_id_of(:), slot_of(:), n_patch(:)
      real(wp),    allocatable :: patch_area_of(:), site_patch_area(:)
      character(len=:), allocatable :: msg
      character(len=24)        :: pid_str

      found = .false.
      call read_census_file(path, table, msg)
      if (len(msg) > 0) then
         if (present(why)) why = msg
         return
      end if
      if (present(why)) why = 'holds no cohort for the site'

      !----- The census site of each site built: the given ones, or the first one the file names. -!
      if (present(site_ids)) then
         sel = site_ids
      else if (present(use_site_id)) then
         sel = [use_site_id]
      else if (table%n > 0_ik) then
         sel = [table%site_id(1)]
      else
         allocate(sel(0))
      end if

      !----- The distinct patches of the selected sites, first seen first, and their areas. ----!
      allocate(site_of(0), patch_id_of(0), patch_area_of(0))
      do r = 1_ik, table%n
         k = int(findloc(sel, table%site_id(r), dim=1), ik)
         if (k == 0_ik) cycle
         ip = census_patch(k, table%patch_id(r))
         if (ip == 0_ik) then
            site_of       = [site_of, k]
            patch_id_of   = [patch_id_of, table%patch_id(r)]
            patch_area_of = [patch_area_of, table%patch_area(r)]
         else if (table%has_area) then
            if (abs(table%patch_area(r) - patch_area_of(ip)) > PATCH_AREA_RTOL * patch_area_of(ip)) then
               write(pid_str, '(i0)') table%patch_id(r)
               if (present(why)) why = 'has two different patch_area values for patch_id '//trim(pid_str)
               return
            end if
         end if
      end do
      npatch = int(size(patch_id_of), ik)
      if (npatch == 0_ik) return                           ! no rows for the selected sites

      !----- The patches, site by site, first seen first within a site: equal areas, or the file's   !
      !      areas normalized to their site. -------------------------------------------------------!
      allocate(slot_of(npatch), n_patch(size(sel)))
      j = 0_ik
      do k = 1_ik, int(size(sel), ik)
         n_patch(k) = int(count(site_of == k), ik)
         do i = 1_ik, npatch
            if (site_of(i) /= k) cycle
            j = j + 1_ik
            slot_of(i) = j
         end do
      end do
      if (present(site_ids)) then
         call init_bare_sites(poly, cfg, n_patch, site_cfg, site_area)
      else
         call init_bare_ground(poly, cfg, npatch)
      end if
      if (table%has_area) then
         do k = 1_ik, int(size(sel), ik)
            site_patch_area = pack(patch_area_of, site_of == k)
            site_patch_area = site_patch_area / sum(site_patch_area)
            j = 0_ik
            do i = 1_ik, npatch
               if (site_of(i) /= k) cycle
               j = j + 1_ik
               poly%patch%area(slot_of(i)) = site_patch_area(j)
            end do
         end do
      end if

      !----- Every cohort of the selected sites, into its patch, in file order. ------------------!
      do r = 1_ik, table%n
         k = int(findloc(sel, table%site_id(r), dim=1), ik)
         if (k == 0_ik) cycle
         ip = census_patch(k, table%patch_id(r))
         if (table%pft(r) < 1_ik .or. table%pft(r) > cfg%pft%n) then
            if (present(why)) why = 'has a pft index out of range'
            return
         end if
         call add_cohort(poly, cfg, slot_of(ip), table%pft(r), table%nplant(r), table%dbh(r))
      end do

      call finalize_init(poly)
      found = .true.
      if (present(why)) why = ''

   contains

      !----- The patch of site k with census patch_id pid, in first-seen order (0 if new). ----!
      integer(ik) function census_patch(k, pid) result(ip)
         integer(ik), intent(in) :: k, pid
         integer(ik) :: i
         ip = 0_ik
         do i = 1_ik, int(size(patch_id_of), ik)
            if (site_of(i) == k .and. patch_id_of(i) == pid) then
               ip = i
               return
            end if
         end do
      end function census_patch
   end subroutine init_from_census

   !---------------------------------------------------------------------------------------!
   ! Restructure a census stand before the first step with the slow step's own operators   !
   ! and switches: its monthly cohort block, then its yearly patch block, less recruitment !
   ! and disturbance. A census carries a cohort per measured size and a patch per plot     !
   ! cell, far more than a step should carry (docs/configuration.md, How a run starts).    !
   !---------------------------------------------------------------------------------------!
   subroutine restructure_census_stand(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(ed_config_t), intent(in)      :: cfg
      real(wp) :: water_shed, kill_carbon   !< what the cull handed on; a start books no ledger
      type(kill_report_t) :: kills
      if (.not. cfg%demography%demography_on) return
      if (cfg%demography%cohort_fusion_on) call fuse_cohorts(poly, cfg)
      call terminate_cohorts(poly, cfg, water_shed, kills)
      call kill_litter(poly, kills, kill_carbon)
      if (cfg%demography%cohort_fission_on) call split_cohorts(poly, cfg)
      call sort_cohorts(poly)
      call sort_patches(poly)
      if (cfg%demography%patch_fusion_on) call fuse_patches(poly, cfg)
      call terminate_patches(poly, cfg)
      if (cfg%demography%cohort_fusion_on) call fuse_cohorts(poly, cfg)
      call terminate_cohorts(poly, cfg, water_shed, kills)
      call kill_litter(poly, kills, kill_carbon)
      call sort_cohorts(poly)
   end subroutine restructure_census_stand

   !----- A line carries data unless it is blank or a '#' comment. ------------------------!
   logical function is_data_line(line)
      character(len=*), intent(in) :: line
      character(len=len(line))     :: s
      s = adjustl(line)
      is_data_line = (len_trim(s) > 0) .and. (s(1:1) /= '#')
   end function is_data_line

   !---------------------------------------------------------------------------------------!
   ! The file's layout, from its first data line: seven numbers are a header-less row read   !
   ! positionally; anything else is the header, and each name is mapped to its column.       !
   !---------------------------------------------------------------------------------------!
   subroutine census_layout(line, cols, why)
      character(len=*),              intent(in)    :: line
      type(census_columns_t),        intent(out)   :: cols
      character(len=:), allocatable, intent(inout) :: why   !< what is wrong with the header ('' if nothing)
      character(len=CENSUS_FIELD_LEN) :: names(CENSUS_MAX_COLS)
      integer(ik) :: sid, pid, cid, ipft
      real(wp)    :: dbh, height, nplant
      logical     :: numeric
      integer     :: k

      call parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, numeric)
      if (numeric) then
         cols%header = .false.
         cols%n = 7
         cols%site = 1
         cols%patch = 2
         cols%cohort = 3
         cols%dbh = 4
         cols%height = 5
         cols%pft = 6
         cols%nplant = 7
         cols%area = 0
         return
      end if
      cols%header = .true.
      call split_csv(line, names, cols%n)
      if (cols%n > CENSUS_MAX_COLS) then
         why = 'has a header with too many columns'
         return
      end if
      do k = 1, cols%n
         select case (lower(trim(names(k))))
         case ('site_id')
            call claim(cols%site,   k, 'site_id')
         case ('patch_id')
            call claim(cols%patch,  k, 'patch_id')
         case ('cohort_id')
            call claim(cols%cohort, k, 'cohort_id')
         case ('dbh')
            call claim(cols%dbh,    k, 'dbh')
         case ('height')
            call claim(cols%height, k, 'height')
         case ('pft')
            call claim(cols%pft,    k, 'pft')
         case ('nplant')
            call claim(cols%nplant, k, 'nplant')
         case ('patch_area')
            call claim(cols%area,   k, 'patch_area')
         case default
            why = 'has an unknown column "'//trim(names(k))//'"'
         end select
         if (len(why) > 0) return
      end do
      if (cols%patch  == 0) why = 'has no patch_id column'
      if (cols%dbh    == 0) why = 'has no dbh column'
      if (cols%pft    == 0) why = 'has no pft column'
      if (cols%nplant == 0) why = 'has no nplant column'
   contains
      subroutine claim(slot, k, name)
         integer,          intent(inout) :: slot
         integer,          intent(in)    :: k
         character(len=*), intent(in)    :: name
         if (slot /= 0) why = 'has the column "'//name//'" twice'
         slot = k
      end subroutine claim
   end subroutine census_layout

   !----- One data row under the layout; an absent site, area or height is 1. -------------!
   subroutine parse_census_line(line, cols, sid, pid, dbh, height, ipft, nplant, area, why)
      character(len=*),              intent(in)    :: line
      type(census_columns_t),        intent(in)    :: cols
      integer(ik),                   intent(out)   :: sid, pid, ipft
      real(wp),                      intent(out)   :: dbh, height, nplant, area
      character(len=:), allocatable, intent(inout) :: why   !< what is wrong with the row ('' if nothing)
      character(len=CENSUS_FIELD_LEN) :: f(CENSUS_MAX_COLS)
      integer(ik) :: cid
      integer     :: nf, ios
      logical     :: ok

      area = 1.0_wp
      if (.not. cols%header) then
         call parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, ok)
         if (.not. ok) then
            why = 'has a malformed row'
            return
         end if
      else
         call split_csv(line, f, nf)
         if (nf < cols%n) then
            why = 'has a row with too few columns'
            return
         end if
         sid = 1_ik
         height = 1.0_wp
         ios = 0
         if (cols%site   > 0 .and. ios == 0) read(f(cols%site),   *, iostat=ios) sid
         if (ios == 0) read(f(cols%patch),  *, iostat=ios) pid
         if (ios == 0) read(f(cols%dbh),    *, iostat=ios) dbh
         if (ios == 0) read(f(cols%pft),    *, iostat=ios) ipft
         if (ios == 0) read(f(cols%nplant), *, iostat=ios) nplant
         if (cols%height > 0 .and. ios == 0) read(f(cols%height), *, iostat=ios) height
         if (cols%area   > 0 .and. ios == 0) read(f(cols%area),   *, iostat=ios) area
         if (ios /= 0) then
            why = 'has a malformed row'
            return
         end if
      end if
      if (dbh <= 0.0_wp .or. height <= 0.0_wp .or. nplant <= 0.0_wp .or. area <= 0.0_wp)         &
         why = 'has a row with a non-positive dbh, height, nplant or patch_area'
   end subroutine parse_census_line

   !----- Parse one comma-separated census row (commas -> blanks, then list-directed read). !
   subroutine parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, ok)
      character(len=*), intent(in)  :: line
      integer(ik),      intent(out) :: sid, pid, cid, ipft
      real(wp),         intent(out) :: dbh, height, nplant
      logical,          intent(out) :: ok
      character(len=len(line)) :: buf
      integer :: k, ios
      buf = line
      do k = 1, len(buf)
         if (buf(k:k) == ',') buf(k:k) = ' '
      end do
      read(buf, *, iostat=ios) sid, pid, cid, dbh, height, ipft, nplant
      ok = (ios == 0)
   end subroutine parse_census_row

   !----- Split a comma-separated line into its fields, each left-adjusted. nf counts them  !
   !      all, so a caller can tell a line with more fields than `fields` holds.  ----------!
   subroutine split_csv(line, fields, nf)
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: fields(:)
      integer,          intent(out) :: nf
      integer :: k, i0, n
      n = len_trim(line)
      nf = 0
      i0 = 1
      fields = ' '
      do k = 1, n + 1
         if (k > n) then
            nf = nf + 1
            if (nf <= size(fields)) fields(nf) = adjustl(line(i0:n))
         else if (line(k:k) == ',') then
            nf = nf + 1
            if (nf <= size(fields)) fields(nf) = adjustl(line(i0:k-1))
            i0 = k + 1
         end if
      end do
   end subroutine split_csv

   !----- Lower-case copy of an ASCII string (column names are matched case-blind). -------!
   pure function lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: t
      integer :: k, c
      t = s
      do k = 1, len(s)
         c = iachar(s(k:k))
         if (c >= iachar('A') .and. c <= iachar('Z')) t(k:k) = achar(c + 32)
      end do
   end function lower

end module ed_init_census
