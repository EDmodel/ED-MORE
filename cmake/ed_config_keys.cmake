# SPDX-License-Identifier: Apache-2.0
# The keys of the two reference configs, ed_config_main.toml and ed_config_pft.toml.
#
# Every key ED-MORE reads is listed in one of them, set or commented out at its default. That list is
# what a run checks its own config against: a key it does not hold stops the run (ed_config_check).
# The build turns the two files into a Fortran include (ed_write_config_keys), and the test
# check_config_keys.cmake holds them equal to the keys the loader reads.

# ed_toml_keys(<file> <out_var>): every key <file> lists, as "section.key". A key counts whether
# it is set ("key = value") or commented out ("# key = value", at most one space after the #), and
# a section header whether it is live ("[name]") or commented out alone on its line ("# [name]"). A
# table of an array, "[[name]]", is the section "name".
function(ed_toml_keys file out_var)
   file(READ "${file}" text)
   # A CMake list splits on ';' and groups on '[ ]': park all three before splitting into lines.
   string(REPLACE ";" "<semi>" text "${text}")
   string(REPLACE "[" "<lb>" text "${text}")
   string(REPLACE "]" "<rb>" text "${text}")
   string(REPLACE "\n" ";" lines "${text}")
   set(section "")
   set(keys "")
   foreach(line IN LISTS lines)
      if(line MATCHES "^[ \t]*#?[ \t]*<lb><lb>([A-Za-z0-9_.]+)<rb><rb>[ \t]*$")
         set(section "${CMAKE_MATCH_1}")
      elseif(line MATCHES "^[ \t]*#?[ \t]*<lb>([A-Za-z0-9_.]+)<rb>[ \t]*$")
         set(section "${CMAKE_MATCH_1}")
      elseif(line MATCHES "^[ \t]*#[ ]?([A-Za-z_][A-Za-z0-9_.]*)[ \t]*=")
         list(APPEND keys "${section}.${CMAKE_MATCH_1}")
      elseif(line MATCHES "^[ \t]*([A-Za-z_][A-Za-z0-9_.]*)[ \t]*=")
         list(APPEND keys "${section}.${CMAKE_MATCH_1}")
      endif()
   endforeach()
   list(REMOVE_DUPLICATES keys)
   set(${out_var} "${keys}" PARENT_SCOPE)
endfunction()

# A Fortran parameter array of character(len=KEY_LEN) holding <keys>, several to a line.
function(ed_fortran_key_array name keys out_var)
   set(body "")
   set(line "")
   list(LENGTH keys n)
   set(i 0)
   foreach(key IN LISTS keys)
      string(LENGTH "${key}" klen)
      if(klen GREATER 64)
         message(FATAL_ERROR "config key ${key} is longer than KEY_LEN = 64")
      endif()
      math(EXPR i "${i} + 1")
      if(i LESS n)
         set(piece "'${key}',")
      else()
         set(piece "'${key}'")
      endif()
      string(LENGTH "${line} ${piece}" llen)
      if(llen GREATER 100 AND NOT line STREQUAL "")
         string(APPEND body "      ${line} &\n")
         set(line "")
      endif()
      if(line STREQUAL "")
         set(line "${piece}")
      else()
         set(line "${line} ${piece}")
      endif()
   endforeach()
   string(APPEND body "      ${line}]\n")
   set(${out_var} "   character(len=KEY_LEN), parameter :: ${name}(*) = [character(len=KEY_LEN) :: &\n${body}"
       PARENT_SCOPE)
endfunction()

# ed_write_config_keys(<main toml> <pft toml> <include file>): write MAIN_KEYS and PFT_KEYS. The
# file is rewritten only when its text changes, so a re-configure does not force a rebuild.
function(ed_write_config_keys main_toml pft_toml out_inc)
   ed_toml_keys("${main_toml}" main_keys)
   ed_toml_keys("${pft_toml}" pft_keys)
   ed_fortran_key_array(MAIN_KEYS "${main_keys}" main_decl)
   ed_fortran_key_array(PFT_KEYS "${pft_keys}" pft_decl)
   set(text "   !----- Written by CMake (cmake/ed_config_keys.cmake) from ed_config_main.toml and   !\n")
   string(APPEND text "   !      ed_config_pft.toml. Do not edit: list a new key in the reference instead. ---!\n")
   string(APPEND text "${main_decl}${pft_decl}")
   file(WRITE "${out_inc}.new" "${text}")
   configure_file("${out_inc}.new" "${out_inc}" COPYONLY)
endfunction()
