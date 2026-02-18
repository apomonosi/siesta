# Distributed under the OSI-approved BSD 2-Clause License.
#
# Copyright (C) 2020 DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomScalapack
-------------------

Finds the ScaLAPACK library.

This is a simple auto-detection module for the ScaLAPACK library (looks
basically for a library with the name 'scalapack'). It also assumes that
ScaLAPACK is MPI-based and defines a respective dependency on
``MPI::MPI_Fortran.``


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``scalapack``
  The ScaLAPACK library. The name was chosen to be compatible with the name used by the CMake
  export file (which is unfortunately broken in several binary packages).

``Scalapack::Scalapack``
  Name space (and CMake-style capitalized) variant of the scalapack target.


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``SCALAPACK_FOUND``
  True if the system has the SCALAPACK library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``SCALAPACK_DETECTION``
  Whether ScaLAPACK libraries should be detected (default: True). If set to False,
  the settings in ``SCALAPACK_LIBRARY`` will be used without any further checks.

``SCALAPACK_LIBRARY``
  Customized ScaLAPACK library/libraries to use.  If no SCALAPACK library is
  required (e.g. the linker automatically links it) set
  ``SCALAPACK_LIBRARY="NONE"``. If not set or empty, it will use
  ``find_package(scalapack)`` to find the scalapack library (relying on the
  existence of a CMake export file for scalapack). Otherwise, the listed
  libraries will be checked for existence (unless disabled in
  ``SCALAPACK_DETECTION``) and the variable is overwritten to contain the
  libraries with their with full path.

  Note: On Ubuntu Focal (20.4 LTS) system, the installed CMake export file is broken.
  You should not realy on the CMake-autodetection on those systems.

``SCALAPACK_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(SiestaCustomLibraryFinder)

set(_target "SCALAPACK::SCALAPACK")
# indent for cleaner output
message(STATUS "Parsing ScaLAPACK options")
list(APPEND CMAKE_MESSAGE_INDENT "  ")


if(TARGET scalapack)
  set(CUSTOMSCALAPACK_FOUND True)
  set(CustomScalapack_FOUND True)
  set(SCALAPACK_FOUND True)
  set(Scalapack_FOUND True)

  message(STATUS "ScaLAPACK already defined")

else()
  # We are already locating MPI at the top-level, or at least we should do that!
  find_package(MPI ${Find_Scalapack_REQUIRED} QUIET)

  option(SCALAPACK_DETECTION "Whether ScaLAPACK library should be detected" TRUE)
  message(CHECK_START "Locating ScaLAPACK library")
  if(SCALAPACK_DETECTION)
    if(NOT "${SCALAPACK_LIBRARY_DIR}" STREQUAL "")
      message(STATUS "Searching in: ${SCALAPACK_LIBRARY_DIR}")
    endif()

    if("${SCALAPACK_LIBRARY}" STREQUAL "")

      # Try Scalapack via CMake export file
      find_package(scalapack)
      if(scalapack_FOUND)
      	message(STATUS "Found intrinsic package")
        get_target_property(_scalapack_library scalapack INTERFACE_LINK_LIBRARIES)

        if("${_scalapack_library}" STREQUAL "")
          set(_scalapack_library "NONE")
        endif()
        set(SCALAPACK_LIBRARY "${_scalapack_library}" CACHE STRING "ScaLAPACK library to link" FORCE)
        unset(_scalapack_library)
      else()

      	message(STATUS "Trying to use pkg-config")

        # Very simple ScaLAPACK auto-detection: looking for a library called scalapack
        # The following logic comes from SIRIUS, and we add our own SCALAPACK_LIBRARY_DIR
      	# to the list of hints
      	find_package(PkgConfig REQUIRED QUIET)

        # If found, this command will set _SCALAPACK_LIBRARY_DIRS
        pkg_search_module(_SCALAPACK scalapack QUIET)
        find_library(SCALAPACK_LIBRARY
         NAMES scalapack scalapack-openmpi
         HINTS
         ${SCALAPACK_LIBRARY_DIR}
         ${_SCALAPACK_LIBRARY_DIRS}
         ENV SCALAPACK_ROOT
         /usr
         PATH_SUFFIXES lib
         DOC "scalapack library path")

      endif()

    elseif(NOT "${SCALAPACK_LIBRARY}" STREQUAL "NONE")

      message(STATUS "Using user-defined variables")

      # ON = find_quietly
      siesta_find_custom_libraries("${SCALAPACK_LIBRARY}" "${SCALAPACK_LIBRARY_DIR}" ON _libs)
      set(SCALAPACK_LIBRARY "${_libs}" CACHE STRING "List of ScaLAPACK libraries to link" FORCE)
      unset(_libs)

    endif()

    set(SCALAPACK_DETECTION False CACHE BOOL "Whether ScaLAPACK libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomScalapack 
    REQUIRED_VARS SCALAPACK_LIBRARY MPI_Fortran_FOUND)

  set(CUSTOMSCALAPACK_FOUND ${CustomScalapack_FOUND})
  set(SCALAPACK_FOUND ${CustomScalapack_FOUND})
  set(Scalapack_FOUND ${CustomScalapack_FOUND})

  if( SCALAPACK_FOUND )
    message(CHECK_PASS "found")
    message(STATUS "ScaLAPACK library: ${SCALAPACK_LIBRARY}")
    message(STATUS "ScaLAPACK link flags: ${SCALAPACK_LINKER_FLAG}")

    if(NOT TARGET scalapack)
      add_library(scalapack INTERFACE IMPORTED)
      if(NOT "${SCALAPACK_LIBRARY}" STREQUAL "NONE")
        target_link_libraries(scalapack INTERFACE "${SCALAPACK_LIBRARY}")
      endif()

    endif()

  endif()

  mark_as_advanced(SCALAPACK_DETECTION SCALAPACK_LIBRARY SCALAPACK_LIBRARY_DIR)

endif()

# Add namespaced library name variant
if(TARGET scalapack AND NOT TARGET ${_target})

  add_library(${_target} INTERFACE IMPORTED)
  target_link_libraries(${_target}
    INTERFACE
      scalapack
      MPI::MPI_Fortran
  )
  if(TARGET LAPACK::LAPACK)
    # lapack should have the logic for adding BLAS::BLAS, if needed
    target_link_libraries(${_target} INTERFACE LAPACK::LAPACK)
  elseif(TARGET BLAS::BLAS)
    target_link_libraries(${_target} INTERFACE BLAS::BLAS)
  endif()
endif()

if(CustomScalapack_FIND_REQUIRED AND NOT TARGET ${_target})
  message(FATAL_ERROR "Required package SCALAPACK cannot be found")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
