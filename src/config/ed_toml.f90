! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_toml -- a minimal TOML reader for the ED-MORE configuration: `# comments`,            !
! `[section]` headers, `[[name]]` tables (keys `name.<i>.key`), scalar                     !
! int/real/bool/string values and `[a, b, c]` arrays, each key stored as `section.key`.    !
! toml-f is not used: its compiled .mod is compiler-specific, as netCDF-Fortran's is       !
! (docs/building.md). A getter returns its default only for an absent key; a present value !
! it cannot parse stops the run, so a typo in a required value cannot pass as zero; and a  !
! key set twice stops the parse, so a second value cannot be ignored in silence.           !
!==========================================================================================!
module ed_toml
   use ed_kinds, only : wp, ik
   implicit none
   private

   public :: toml_table_t, toml_parse_file
   public :: toml_has, toml_has_section, toml_int, toml_real, toml_logical, toml_string, toml_real_array
   public :: toml_is_array, toml_string_array
   public :: toml_table_count
   public :: toml_write_record, toml_record_reset

   !----- VALLEN is also the line buffer, so it bounds a whole `key = value` line. 1024 holds     !
   !      the forcing path fields (MET_PATH_LEN) with room for the key. --------------------------!
   integer, parameter :: KEYLEN = 64, VALLEN = 1024, MAXKEYS = 512, MAXARRAYS = 8

   type :: toml_table_t
      integer(ik)                     :: n = 0_ik
      character(len=KEYLEN), allocatable :: key(:)
      character(len=VALLEN), allocatable :: val(:)
      character(len=256)              :: source = ''   !< the file read into it (the record's label)
      !----- The arrays of tables, `[[name]]`, and how many tables each holds. ------------------!
      integer(ik)           :: n_array = 0_ik
      character(len=KEYLEN) :: array_name(MAXARRAYS) = ''
      integer(ik)           :: array_count(MAXARRAYS) = 0_ik
   end type toml_table_t

   !----- The parameter record (docs/configuration.md, "The parameter philosophy"): every         !
   !      value a toml_* getter returns, with its file and whether the file sets the key or       !
   !      the default stands. A key never read has no row, so a caller can check that a key it    !
   !      set is read. A re-read replaces its row; toml_record_reset starts a new record for      !
   !      each config a process loads. -----------------------------------------------------------!
   integer, parameter :: MAXREC = 4096, RECVAL = 512
   integer(ik), save :: rec_n = 0_ik
   character(len=256),    save :: rec_source(MAXREC)
   character(len=KEYLEN), save :: rec_key(MAXREC)
   integer(ik),           save :: rec_index(MAXREC)      !< 0 for a scalar, else the array element
   logical,               save :: rec_present(MAXREC)    !< .true. = set in the file, .false. = default
   character(len=RECVAL), save :: rec_value(MAXREC)

contains

   !---------------------------------------------------------------------------------------!
   ! Parse a TOML file into the (section.key -> raw-value) table; ok = .false. if it       !
   ! cannot be opened.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine toml_parse_file(path, t, ok)
      character(len=*),    intent(in)  :: path
      type(toml_table_t),  intent(out) :: t
      logical,             intent(out) :: ok
      integer                :: unit, ios, ieq, nb, ia
      character(len=VALLEN)  :: line
      character(len=KEYLEN)  :: section
      character(len=VALLEN)  :: raw

      ok = .false.
      allocate(t%key(MAXKEYS), t%val(MAXKEYS))
      t%n = 0_ik
      t%source = path
      section = ''

      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      ok = .true.

      do
         read(unit, '(a)', iostat=ios) raw
         if (ios /= 0) exit
         line = strip_comment(raw)
         line = adjustl(line)
         if (len_trim(line) == 0) cycle

         if (line(1:2) == '[[') then                      ! [[name]]: the array's next table
            nb = index(line, ']]')
            if (nb > 3) then
               section = adjustl(line(3:nb-1))
               ia = int(findloc(t%array_name(1:t%n_array), section, dim=1))
               if (ia == 0) then
                  if (t%n_array >= MAXARRAYS) error stop 'ed_toml: too many [[array]] names in '//trim(path)
                  t%n_array = t%n_array + 1_ik
                  ia = int(t%n_array)
                  t%array_name(ia) = section
               end if
               t%array_count(ia) = t%array_count(ia) + 1_ik
               write(section, '(a,".",i0)') trim(t%array_name(ia)), t%array_count(ia)
            end if
            cycle
         end if
         if (line(1:1) == '[') then                       ! [section] header
            nb = index(line, ']')
            if (nb > 2) section = adjustl(line(2:nb-1))
            cycle
         end if

         ieq = index(line, '=')
         if (ieq < 2) cycle                               ! not a key = value line
         if (t%n >= MAXKEYS) exit
         t%n = t%n + 1_ik
         if (len_trim(section) > 0) then
            t%key(t%n) = trim(section) // '.' // trim(adjustl(line(1:ieq-1)))
         else
            t%key(t%n) = trim(adjustl(line(1:ieq-1)))
         end if
         t%val(t%n) = adjustl(line(ieq+1:))
         !----- TOML forbids setting a key twice, and a getter would read only the first value, so a !
         !      second one would be ignored in silence. A [[name]] table's keys carry its index, so  !
         !      two tables may each set a key. ------------------------------------------------------!
         if (find_key(t, t%key(t%n)) < t%n) then
            write(*,'(4a)') ' ed_toml:   key "', trim(t%key(t%n)), '" is set twice in ', trim(path)
            error stop 'ed_toml: a key   set twice'
         end if
      end do
      close(unit)
   end subroutine toml_parse_file

   !----- Strip a trailing `# comment` that is not inside a quoted string. ----------------!
   pure function strip_comment(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: out
      integer :: i, nq
      out = s
      nq = 0
      do i = 1, len(s)
         if (s(i:i) == '"') nq = nq + 1
         if (s(i:i) == '#' .and. mod(nq, 2) == 0) then
            out = s(1:i-1)
            return
         end if
      end do
   end function strip_comment

   !----- Index of a fully-qualified key, or 0 if absent. ---------------------------------!
   pure integer(ik) function find_key(t, key) result(idx)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      integer(ik) :: i
      idx = 0_ik
      do i = 1_ik, t%n
         if (trim(t%key(i)) == trim(key)) then
            idx = i
            return
         end if
      end do
   end function find_key

   logical function toml_has(t, key) result(yes)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      yes = find_key(t, key) > 0_ik
   end function toml_has

   !----- How many tables the array `[[name]]` holds (0 when the file has none). -----------------!
   pure integer(ik) function toml_table_count(t, name) result(n)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: name
      integer(ik) :: ia
      n = 0_ik
      do ia = 1_ik, t%n_array
         if (trim(t%array_name(ia)) == trim(name)) n = t%array_count(ia)
      end do
   end function toml_table_count

   !----- Is there any key under `[section]`? Keys are stored flattened (section.key), so a       !
   !      header leaves no entry for toml_has to find. The gate for an optional block whose       !
   !      keys are required once it is present: gating on one key inside it would skip, in        !
   !      silence, a block that misspells that key. ----------------------------------------------!
   pure logical function toml_has_section(t, section) result(yes)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: section
      character(len=:), allocatable  :: pre
      integer(ik) :: i, n
      pre = trim(section)//'.'
      n   = len(pre)
      yes = .false.
      do i = 1_ik, t%n
         if (len_trim(t%key(i)) > n) then
            if (t%key(i)(1:n) == pre) then
               yes = .true.
               return
            end if
         end if
      end do
   end function toml_has_section

   !----- A present value that cannot be parsed: print its key and raw text, and stop. ---------!
   subroutine toml_parse_error(key, raw, what)
      character(len=*), intent(in) :: key, raw, what
      write(*,'(5a)') ' ed_toml:   cannot parse ', trim(what), ' for key "', trim(key), '"'
      write(*,'(3a)') '   raw value: "', trim(raw), '"'
      error stop 'ed_toml: malformed configuration value'
   end subroutine toml_parse_error

   integer(ik) function toml_int(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      integer(ik),        intent(in) :: default
      integer(ik) :: idx, ios
      v = default
      idx = find_key(t, key)
      if (idx /= 0_ik) then
         read(t%val(idx), *, iostat=ios) v
         if (ios /= 0) call toml_parse_error(key, t%val(idx), 'integer')
      end if
      call record_value(t, key, 0_ik, idx /= 0_ik, int_text(v))
   end function toml_int

   real(wp) function toml_real(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      real(wp),           intent(in) :: default
      integer(ik) :: idx, ios
      v = default
      idx = find_key(t, key)
      if (idx /= 0_ik) then
         read(t%val(idx), *, iostat=ios) v
         if (ios /= 0) call toml_parse_error(key, t%val(idx), 'real')
      end if
      call record_value(t, key, 0_ik, idx /= 0_ik, real_text(v))
   end function toml_real

   logical function toml_logical(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      logical,            intent(in) :: default
      integer(ik) :: idx
      character(len=VALLEN) :: s
      v = default
      idx = find_key(t, key)
      if (idx == 0_ik) then                               ! absent -> default
         call record_value(t, key, 0_ik, .false., merge('true ', 'false', v))
         return
      end if
      s = adjustl(t%val(idx))
      !----- Exact match (documented lowercase true/false + common Fortran/case variants),       !
      !      not a prefix test: 'trueish' is an error, and '.true.' is .true. -------------------!
      select case (trim(s))
      case ('true', '.true.', 'True', 'TRUE')
         v = .true.
      case ('false', '.false.', 'False', 'FALSE')
         v = .false.
      case default
         call toml_parse_error(key, t%val(idx), 'logical')
      end select
      call record_value(t, key, 0_ik, .true., merge('true ', 'false', v))
   end function toml_logical

   !----- String value with surrounding double quotes stripped. ---------------------------!
   function toml_string(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key, default
      character(len=VALLEN)          :: v
      integer(ik) :: idx
      character(len=VALLEN) :: s
      v = default
      idx = find_key(t, key)
      if (idx /= 0_ik) then
         s = adjustl(t%val(idx))
         if (s(1:1) == '"') then
            if (index(s(2:), '"') == 0) call toml_parse_error(key, s, 'string (no closing quote)')
            v = s(2:index(s(2:), '"'))
         else
            v = trim(s)
         end if
      end if
      call record_value(t, key, 0_ik, idx /= 0_ik, v)
   end function toml_string

   !----- True when the key's value is an array literal, [...]. -----------------------------!
   pure logical function toml_is_array(t, key) result(yes)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      integer(ik) :: idx
      yes = .false.
      idx = find_key(t, key)
      if (idx /= 0_ik) yes = index(adjustl(t%val(idx)), '[') == 1
   end function toml_is_array

   !----- A string array, ["a", "b", ...], or one string, which counts as one entry. ------!
   subroutine toml_string_array(t, key, out, nout)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      character(len=*),   intent(out) :: out(:)
      integer(ik),        intent(out) :: nout
      integer(ik) :: idx, lb, rb, q1, q2
      character(len=VALLEN) :: s
      nout = 0_ik
      out  = ''
      idx = find_key(t, key)
      if (idx == 0_ik) return
      if (.not. toml_is_array(t, key)) then
         nout = 1_ik
         out(1) = toml_string(t, key, '')
         return
      end if
      s  = adjustl(t%val(idx))
      lb = index(s, '[')
      rb = index(s, ']', back=.true.)
      if (rb <= lb) call toml_parse_error(key, t%val(idx), 'string array (expected [..] literal)')
      s = s(lb+1:rb-1)
      do
         q1 = index(s, '"')
         if (q1 == 0) exit
         q2 = index(s(q1+1:), '"')
         if (q2 == 0) call toml_parse_error(key, t%val(idx), 'string array (unclosed quote)')
         nout = nout + 1_ik
         if (nout > size(out, kind=ik)) call toml_parse_error(key, t%val(idx), 'string array (too many entries)')
         out(nout) = s(q1+1:q1+q2-1)
         call record_value(t, key, nout, .true., out(nout))
         s = s(q1+q2+1:)
      end do
   end subroutine toml_string_array

   !----- Parse `[a, b, c]` into out(:) reals; nout = the count read (0 if absent). --------!
   subroutine toml_real_array(t, key, out, nout)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      real(wp),           intent(out) :: out(:)
      integer(ik),        intent(out) :: nout
      integer(ik) :: idx, i, lb, rb, ios
      character(len=VALLEN) :: s
      nout = 0_ik
      idx = find_key(t, key)
      if (idx == 0_ik) return                       ! absent -> nout = 0 (caller records missing key)
      s  = t%val(idx)
      lb = index(s, '[')
      rb = index(s, ']', back=.true.)
      if (lb == 0 .or. rb <= lb) call toml_parse_error(key, t%val(idx), 'real array (expected [..] literal)')
      s = s(lb+1:rb-1)
      do i = 1, len_trim(s)                         ! commas -> spaces for list-directed read
         if (s(i:i) == ',') s(i:i) = ' '
      end do
      block
         real(wp)    :: tmp(size(out))
         integer(ik) :: k, ntok
         ntok = count_tokens(s)
         if (ntok == 0_ik) return                   ! empty [] -> nout = 0
         !----- Bound by the receiving buffer before reading, so an over-long list stops         !
         !      instead of overrunning it. ------------------------------------------------------!
         if (ntok > size(out, kind=ik)) call toml_parse_error(key, t%val(idx), 'real array (too many entries)')
         tmp = 0.0_wp
         read(s, *, iostat=ios) (tmp(k), k = 1_ik, ntok)
         if (ios /= 0) call toml_parse_error(key, t%val(idx), 'real array')   ! any unparseable token -> hard error
         out(1:ntok) = tmp(1:ntok)
         nout = ntok                                ! the parsed count, never more
         do k = 1_ik, ntok
            call record_value(t, key, k, .true., real_text(out(k)))
         end do
      end block
   end subroutine toml_real_array

   !----- The parameter record (see the module-level note). ------------------------------------!
   subroutine record_value(t, key, index, present, text)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key, text
      integer(ik),        intent(in) :: index
      logical,            intent(in) :: present
      integer(ik) :: i, slot
      slot = 0_ik
      do i = 1_ik, rec_n
         if (rec_index(i) == index .and. rec_key(i) == key .and. rec_source(i) == t%source) then
            slot = i
            exit
         end if
      end do
      if (slot == 0_ik) then
         if (rec_n >= MAXREC) return
         rec_n = rec_n + 1_ik
         slot = rec_n
      end if
      rec_source(slot) = t%source
      rec_key(slot) = key
      rec_index(slot) = index
      rec_present(slot) = present
      rec_value(slot) = adjustl(text)
   end subroutine record_value

   !----- Empty the record: what a run reads from here on is the whole of its record. -------!
   subroutine toml_record_reset()
      rec_n = 0_ik
   end subroutine toml_record_reset

   !----- 17 significant digits: the value the run used, to the last bit of a double. -------!
   function real_text(x) result(s)
      real(wp), intent(in) :: x
      character(len=RECVAL) :: s
      write(s, '(es24.16e3)') x
      s = adjustl(s)
   end function real_text

   function int_text(i) result(s)
      integer(ik), intent(in) :: i
      character(len=RECVAL) :: s
      write(s, '(i0)') i
   end function int_text

   !----- Write the record as CSV: source,key,index,present,value, one row per scalar key and one  !
   !      per array element. Text fields are quoted. A failed open warns and returns. ------------!
   subroutine toml_write_record(path)
      character(len=*), intent(in) :: path
      integer     :: u, ios
      integer(ik) :: i
      open(newunit=u, file=path, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write(*,'(3a)') ' warning: could not write the parameter record "', trim(path), '"'
         return
      end if
      write(u, '(a)') 'source,key,index,present,value'
      do i = 1_ik, rec_n
         write(u, '(a,",",a,",",i0,",",a,",",a)') csv_quote(rec_source(i)), trim(rec_key(i)), rec_index(i), &
              trim(merge('true ', 'false', rec_present(i))), csv_quote(rec_value(i))
      end do
      close(u)
   end subroutine toml_write_record

   function csv_quote(text) result(q)
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: q
      integer :: i
      q = '"'
      do i = 1, len_trim(text)
         if (text(i:i) == '"') then
            q = q // '""'
         else
            q = q // text(i:i)
         end if
      end do
      q = q // '"'
   end function csv_quote

   !----- Number of whitespace-separated tokens in a string. ------------------------------!
   pure integer(ik) function count_tokens(s) result(nt)
      character(len=*), intent(in) :: s
      integer(ik) :: i
      logical     :: in_tok
      nt = 0_ik
      in_tok = .false.
      do i = 1_ik, len_trim(s)
         if (s(i:i) /= ' ') then
            if (.not. in_tok) nt = nt + 1_ik
            in_tok = .true.
         else
            in_tok = .false.
         end if
      end do
   end function count_tokens

end module ed_toml
