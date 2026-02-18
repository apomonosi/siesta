# Distributed under the OSI-approved BSD 2-Clause License.
#
# Copyright (C) 2022  DFTB+ developers group (derivation of FindCustomLapack.cmake)
#

#[=======================================================================[.rst:
FindCustomBLAS
----------------

Finds the BLAS library

This is a wrapper around CMakes FindBLAS module with the additional
possibility to customize the library name manually. In latter case the module will
check the existence of those libraries and stop if they are not found.

Note: The module is named FindCustomBLAS (and not FindBLAS) to avoid name
collision with CMakes built-in FindBLAS module.


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``BLAS::BLAS``
  The BLAS library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``BLAS_FOUND``
  True if the system has the BLAS library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``BLAS_DETECTION``
  Whether BLAS libraries should be detected (default: True). If set to False,
  the settings in ``BLAS_LIBRARY`` will be used without any further checks.

``BLAS_LIBRARY``
  Customized BLAS library/libraries to use (instead of autodetected ones).  If
  no BLAS library is required (e.g. the linker automatically links it) set
  ``BLAS_LIBRARY="NONE"``. If not set or empty, the built-in BLAS finder
  (the findBLAS module) will be invoked. Otherwise, the listed libraries
  will be checked for existence (unless disabled in ``BLAS_DETECTION``) and
  the variable is overwritten to contain the libraries with their with full
  path.

``BLAS_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

``BLAS_LINKER_FLAG``
  Flags to use when linking BLAS

Additionally, the cache variables of the built-in FindBLAS modules may used to
influence the BLAS detection if the built-in module is invoked.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(SiestaCustomLibraryFinder)

# indent for cleaner output
message(STATUS "Parsing BLAS options")
list(APPEND CMAKE_MESSAGE_INDENT "  ")


if(TARGET BLAS::BLAS)

  set(CUSTOMBLAS_FOUND True)
  set(CustomBlas_FOUND True)
  set(BLAS_FOUND True)
  set(Blas_FOUND True)

else()

  option(BLAS_DETECTION "Whether BLAS library should be detected" TRUE)
  message(CHECK_START "Locating BLAS library")

  if(BLAS_DETECTION)
    
    if(NOT "${BLAS_LIBRARY_DIR}" STREQUAL "")
      message(STATUS "Searching in: ${BLAS_LIBRARY_DIR}")
    endif()
    
    # BLAS has either not been found yet or it was found by an older built-in findBLAS module.
    # which does not provide the imported target BLAS::BLAS

    if("${BLAS_LIBRARY}" STREQUAL "")

      # No user customized BLAS library, try built-in finder
      if(NOT BLAS_FOUND)
        find_package(BLAS)
        if( BLAS_FOUND AND "${BLAS_LIBRARIES}" STREQUAL "")
          set(BLAS_LIBRARIES "NONE")
        endif()
      endif()
      set(BLAS_LIBRARY "${BLAS_LIBRARIES}" CACHE STRING "BLAS library to link" FORCE)
      set(BLAS_LINKER_FLAG "${BLAS_LINKER_FLAGS}" CACHE STRING
        "Linker flags to use when linking BLAS" FORCE)

    elseif(NOT "${BLAS_LIBRARY}" STREQUAL "NONE")

      # BLAS explicitely set by the user, search for those libraries
      # ON = find_quietly
      siesta_find_custom_libraries("${BLAS_LIBRARY}" "${BLAS_LIBRARY_DIR}" ON _libs)
      set(BLAS_LIBRARY "${_libs}" CACHE STRING "List of BLAS libraries to link" FORCE)
      unset(_libs)

    endif()

    set(BLAS_DETECTION False CACHE BOOL "Whether BLAS libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomBlas REQUIRED_VARS BLAS_LIBRARY)

  set(CUSTOMBLAS_FOUND ${CustomBlas_FOUND})
  set(BLAS_FOUND ${CustomBlas_FOUND})
  set(Blas_FOUND ${CustomBlas_FOUND})

  if( BLAS_FOUND )
    message(CHECK_PASS "found")
    message(STATUS "BLAS library: ${BLAS_LIBRARY}")
    message(STATUS "BLAS link flags: ${BLAS_LINKER_FLAG}")
  else()
    message(CHECK_FAIL "not found")
  endif()

  if (BLAS_FOUND AND NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    if(NOT "${BLAS_LIBRARY}" STREQUAL "NONE")
      target_link_libraries(BLAS::BLAS INTERFACE "${BLAS_LIBRARY}")
    endif()
    if(NOT "${BLAS_LIBRARY_DIR}" STREQUAL "")
      target_link_directories(BLAS::BLAS BEFORE INTERFACE "${BLAS_LIBRARY_DIR}")
    endif()
    if(NOT "${BLAS_LINKER_FLAG}" STREQUAL "")
      target_link_options(BLAS::BLAS INTERFACE "${BLAS_LINKER_FLAG}")
    endif()
  endif()

  mark_as_advanced(BLAS_DETECTION BLAS_LIBRARY BLAS_LIBRARY_DIR BLAS_LINKER_FLAG)

endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
