! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_test_run -- the configs a test derives from a template config, and the edmore
! runs it expects to be refused.
!
!   derive_config writes a config from override blocks and a template. The blocks come
! first and are layered: the first value a block gives a key is the one written, and a
! later block's value for it is left out. Then the template follows, less every key the
! blocks set and any key the caller drops. Each key is written once, as ed_toml requires:
! it refuses a key set twice. Keys are named as ed_toml names them: `section.key`, and
! `name.<i>.key` in the i-th `[[name]]` table, so two tables may each set a key.
!
!   main_refuses runs edmore on a config and tells whether it stopped and its report
! holds a message.
!---~---
module ed_test_run
   implicit none
   private

   public :: derive_config, main_refuses

   integer, parameter :: LINE_LEN = 2048    ! longest config line read
   integer, parameter :: KEY_LEN  = 128     ! longest qualified key
   integer, parameter :: MAX_ARRAYS = 8     ! [[name]] arrays in one config

   !--- The keys written so far and the [[name]] tables met, as ed_toml counts them.
   type :: key_set_t
      integer :: n = 0
      character(len=KEY_LEN), allocatable :: key(:)
      character(len=KEY_LEN) :: section = ''
      integer :: n_array = 0
      character(len=KEY_LEN) :: array_name(MAX_ARRAYS) = ''
      integer :: array_count(MAX_ARRAYS) = 0
   end type key_set_t

contains

   !---~---
   !   Write `path`: the lines of `overrides` (one string, its lines split by new_line), each
   ! key at its first value, then the lines of the file `template` whose keys the overrides
   ! do not set and `drop` (keys as section.key) does not name. Comments and headers pass
   ! through.
   !---~---
   subroutine derive_config(path, overrides, template, drop)
      character(len=*),           intent(in) :: path, overrides, template
      character(len=*), optional, intent(in) :: drop(:)
      character(len=LINE_LEN) :: line
      character(len=KEY_LEN)  :: key
      type(key_set_t) :: set
      integer :: uin, uout, ios, i0, i1
      allocate(set%key(64))
      open(newunit=uout, file=path, status='replace', action='write')
      !--- The override blocks, line by line.
      i0 = 1
      do while (i0 <= len(overrides))
         i1 = index(overrides(i0:), new_line('a'))
         if (i1 == 0) then
            i1 = len(overrides) + 1
         else
            i1 = i0 + i1 - 1
         end if
         line = overrides(i0:i1 - 1)
         i0 = i1 + 1
         key = line_key(set, line)
         if (len_trim(key) > 0) then
            if (holds(set, key)) cycle
            call add_key(set, key)
         end if
         write(uout, '(a)') trim(line)
      end do
      !---~---
      !   The template, less the keys set above and those dropped. Its keys are not added to
      ! the set: a key the template itself sets twice is written twice, for ed_toml to
      ! refuse, as it would refuse the template.
      !---~---
      open(newunit=uin, file=template, status='old', action='read')
      do
         read(uin, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == LINE_LEN) error stop 'ed_test_run: a template line fills the buffer'
         key = line_key(set, line)
         if (len_trim(key) > 0) then
            if (holds(set, key)) cycle
            if (present(drop)) then
               if (any(drop == key)) cycle
            end if
         end if
         write(uout, '(a)') trim(line)
      end do
      close(uin)
      close(uout)
   end subroutine derive_config

   !---~---
   !   Run edmore on `config`, its output into `log`: true when it stopped with an error
   ! and a line of the log holds `message`.
   !---~---
   logical function main_refuses(edmore, config, log, message) result(refused)
      character(len=*), intent(in) :: edmore, config, log, message
      character(len=LINE_LEN) :: line
      integer :: es, u, ios
      logical :: said
      call execute_command_line('"'//edmore//'" "'//config//'" > "'//log//'" 2>&1', exitstat=es)
      said = .false.
      open(newunit=u, file=log, status='old', action='read')
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, message) > 0) said = .true.
      end do
      close(u)
      refused = es /= 0 .and. said
   end function main_refuses

   !---~---
   !   The qualified key a line sets, or '' for a comment, a blank or a header. A header
   ! moves the section; a [[name]] header opens that array's next table.
   !---~---
   function line_key(set, raw) result(key)
      type(key_set_t),  intent(inout) :: set
      character(len=*), intent(in)    :: raw
      character(len=KEY_LEN) :: key
      character(len=len(raw)) :: line
      integer :: nb, ieq, ia
      key = ''
      line = adjustl(strip_comment(raw))
      if (len_trim(line) == 0) return
      if (line(1:2) == '[[') then
         nb = index(line, ']]')
         if (nb > 3) then
            ia = findloc(set%array_name(1:set%n_array), adjustl(line(3:nb - 1)), dim=1)
            if (ia == 0) then
               if (set%n_array >= MAX_ARRAYS) error stop 'ed_test_run: too many [[array]] names'
               set%n_array = set%n_array + 1
               ia = set%n_array
               set%array_name(ia) = adjustl(line(3:nb - 1))
            end if
            set%array_count(ia) = set%array_count(ia) + 1
            write(set%section, '(a,".",i0)') trim(set%array_name(ia)), set%array_count(ia)
         end if
         return
      end if
      if (line(1:1) == '[') then
         nb = index(line, ']')
         if (nb > 2) set%section = adjustl(line(2:nb - 1))
         return
      end if
      ieq = index(line, '=')
      if (ieq < 2) return
      if (len_trim(set%section) > 0) then
         key = trim(set%section)//'.'//trim(adjustl(line(1:ieq - 1)))
      else
         key = trim(adjustl(line(1:ieq - 1)))
      end if
   end function line_key

   !--- A line without its trailing `# comment`, a `#` inside a quoted string kept.
   pure function strip_comment(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, nq
      out = s
      nq = 0
      do i = 1, len(s)
         if (s(i:i) == '"') nq = nq + 1
         if (s(i:i) == '#' .and. mod(nq, 2) == 0) then
            out = s(1:i - 1)
            return
         end if
      end do
   end function strip_comment

   pure logical function holds(set, key)
      type(key_set_t),  intent(in) :: set
      character(len=*), intent(in) :: key
      holds = any(set%key(1:set%n) == key)
   end function holds

   subroutine add_key(set, key)
      type(key_set_t),  intent(inout) :: set
      character(len=*), intent(in)    :: key
      character(len=KEY_LEN), allocatable :: grown(:)
      if (set%n == size(set%key)) then
         allocate(grown(2 * set%n))
         grown(1:set%n) = set%key(1:set%n)
         call move_alloc(grown, set%key)
      end if
      set%n = set%n + 1
      set%key(set%n) = key
   end subroutine add_key

end module ed_test_run
