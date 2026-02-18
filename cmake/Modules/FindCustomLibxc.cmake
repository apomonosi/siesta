#
# Copyright (C) 2023 Siesta developers group.
#
#
#  Search for libxc library.
#  This is a custom module, as it needs to support several ways of discovery,
#[=======================================================================[.rst:
FindCustomLibxc
----------------

This will find and locate the XC library using pkg-config.

On entry the following variables may be used:

LIBXC_Fortran_INTERFACE :
  a cmake-list (; separated) with the modules that may be
  used in Siesta. It also holds the order of their usage.
  For instance LIBXC_Fortran_INTERFACE=f90;f03 will search for both
  f90 and f03, however, if f90 is found, it will _not_ search for the f03
  module.

Notes:
  for now this find-package is a custom wrapper to allow future libxc versions
  to ship their own cmake package.

The target:
   Libxc::xc
   Libxc::xc_Fortran
will be usable upon return.
While it might have been more obvious to use libxc::xc we need to use the
same upstream name. Libxc project uses Libxc::xc|xcf03

#]=======================================================================]
# The searched fortran interfaces may be controlled with
#   LIBXC_Fortran_INTERFACE = f90;f03
# which will prefer libxcf90, and if not found, then libxcf03.
#
find_package(PkgConfig REQUIRED QUIET)

message(CHECK_START "Searching for libXC library")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# This will be a required module
pkg_check_modules(LIBXC_C libxc)

message(DEBUG "LIBXC_C_LIBDIR: ${LIBXC_C_LIBDIR}")
message(DEBUG "LIBXC_C_LIBRARIES: ${LIBXC_C_LIBRARIES}")
message(DEBUG "LIBXC_C_LINK_LIBRARIES: ${LIBXC_C_LINK_LIBRARIES}")
message(DEBUG "LIBXC_C_INCLUDEDIR: ${LIBXC_C_INCLUDEDIR}")
message(DEBUG "LIBXC_C_INCLUDE_DIRS: ${LIBXC_C_INCLUDE_DIRS}")

if(NOT TARGET Libxc::xc)
  # we need to add the target if it does not exist
  add_library(Libxc::xc INTERFACE IMPORTED)

  target_link_libraries(
    Libxc::xc
    INTERFACE "${LIBXC_C_LINK_LIBRARIES}"
  )
  target_include_directories(
    Libxc::xc
    INTERFACE "${LIBXC_C_INCLUDE_DIRS}" "${LIBXC_C_INCLUDEDIR}"
  )
endif()


# Now we need to search for the fortran libraries
set(LIBXC_Fortran_INTERFACE "f03;f90" CACHE STRING
  "Which fortran libxc interfaces to use (in order of priority) (LIBXC_Fortran_INTERFACE)")
mark_as_advanced(LIBXC_Fortran_INTERFACE)
foreach(xcv IN LISTS LIBXC_Fortran_INTERFACE)
  message(CHECK_START "Using pkg-config to search for libxc${xcv}")
  string(TOUPPER "${xcv}" xcV)

  # search for the fortran interface
  pkg_check_modules(LIBXC_${xcV} libxc${xcv})

  # Pass upwards without knowing interface name
  set(LIBXC_Fortran_FOUND "${LIBXC_${xcV}_FOUND}")

  if( "${LIBXC_${xcV}_FOUND}" )
    message(CHECK_PASS "found")
    
    message(DEBUG "LIBXC_${xcV}_LIBDIR: ${LIBXC_${xcV}_LIBDIR}")
    message(DEBUG "LIBXC_${xcV}_LIBRARIES: ${LIBXC_${xcV}_LIBRARIES}")
    message(DEBUG "LIBXC_${xcV}_LINK_LIBRARIES: ${LIBXC_${xcV}_LINK_LIBRARIES}")
    message(DEBUG "LIBXC_${xcV}_INCLUDEDIR: ${LIBXC_${xcV}_INCLUDEDIR}")
    message(DEBUG "LIBXC_${xcV}_INCLUDE_DIRS: ${LIBXC_${xcV}_INCLUDE_DIRS}")

    add_library(Libxc::xc${xcv} INTERFACE IMPORTED)

    target_link_libraries(
      Libxc::xc${xcv}
      INTERFACE "${LIBXC_${xcV}_LINK_LIBRARIES}"
    )

    target_include_directories(
      Libxc::xc${xcv}
      INTERFACE "${LIBXC_${xcV}_INCLUDE_DIRS}" "${LIBXC_${xcV}_INCLUDEDIR}"
    )
  
    add_library(Libxc::xc_Fortran INTERFACE IMPORTED)
    target_link_libraries(Libxc::xc_Fortran INTERFACE Libxc::xc${xcv} Libxc::xc)

    break()
  else()
    message(CHECK_FAIL "not found")
  endif()

endforeach()


if(NOT TARGET Libxc::xc_Fortran)
  set(LIBXC_Fortran_FOUND FALSE)

  message(STATUS "Could not find any libxc fortran libraries.")
  message(STATUS "Searched (in priority order):")
  foreach(xcv IN LISTS LIBXC_Fortran_INTERFACE)
    string(TOUPPER "${xcv}" xcV)
    message(STATUS "  - libxc${xcv} (FOUND = ${LIBXC_${xcV}_FOUND})")
  endforeach()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  message(CHECK_FAIL "not found")

else()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  message(CHECK_PASS "found")

endif()

